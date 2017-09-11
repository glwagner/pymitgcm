from __future__ import division
import os, subprocess, shutil, time
import numpy as np
import matplotlib.pyplot as plt
import f90nml

import gentools

# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Prepare set up
gcmdir  = '/Users/glwagner/Software/MITgcm'
optfile = '{}/tools/build_options/neve'.format(gcmdir)

# Domain
Lx = 13.3e3
Ly = 5.0e3
Lz = 200.0

nx = 100
ny = 1
nz = 100

# Run and make parallelism
nprun  = 1
npmake = 2

dx, dy, dz = Lx/nx, Ly/ny, Lz/nz

# Physical parameters
g     = 9.81
alpha = 2.0e-4
N     = 1.0e-3

# Cell centers
x = np.arange(dx, Lx+dx, dx)
z = np.arange(-0.5*dz, -Lz, -dz)

# Stratification
Tz = N**2.0 / (g*alpha)
Tref = Tz*z - (Tz*z).mean()

# Initial temperature profile
init_T = np.ones((nx,))[:, np.newaxis] * Tref[np.newaxis, :]

# Topography
hbump = 0.6*Lz
lbump = 0.1*Lx
topo = -Lz + hbump * np.exp( -(x-0.5*Lx)**2.0 / (2.0*lbump**2.0) )
topo[-1] = 0.0

# Open boundaries
U0 = 0.02
obc_U = U0 * np.ones((nz,))
obc_T = init_T[0, :].squeeze()


# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Directories
workdir  = os.getcwd()

setupdirs = {
    'builddir' : '{}/build'.format(workdir),
    'codedir'  : '{}/code'.format(workdir),
    'inputdir' : '{}/input'.format(workdir),
    'rundir'   : '{}/run'.format(workdir),
}

templdir = '/Users/glwagner/Numerics/gcmprocess/templates/leewave'
nametempldir = '{}/namelists'.format(templdir)
codetempldir = '{}/code'.format(templdir)

# Remove and remake paths.
for path in setupdirs.values():
    if os.path.exists(path):
        shutil.rmtree(path)

for path in setupdirs.values():
    os.makedirs(path)

# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Manage compiled 'code' files

# Copy
for filename in os.listdir(codetempldir):
    shutil.copy('{}/{}'.format(codetempldir, filename), 
        '{}/{}'.format(setupdirs['codedir'], filename)) 

# Modify
if nx % nprun == 0:
    sizevars = {'sNx': int(nx/nprun), 'Nr': nz, 'nPx': nprun}
else:
    raise ValueError("Number of grid points must be a mulitple of the number "
        "of run processes.")

for var, value in sizevars.items():
    gentools.change_size_var(var, value, codedir=setupdirs['codedir'])

# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Custom namelist params
patches = { 
    'data': {
        'parm01': {
            'Tref'    : list(np.round(Tref*100)/100),
            #'Sref'    : [ 35.0 for i in range(nz) ],
            'viscAz'  : 2.0e-3,
            'viscAh'  : 1.0e-2,
            'diffKhT' : 1.0e-2,
            'diffKzT' : 1.0e-3,
        },
        'parm03': {
            'nTimeSteps' : 10000,
            'deltaT'     : 100.0,
            'pChkptFreq' : 100000.0,
            'chkptFreq'  : 100000.0,
        },
        'parm04': {
            'delZ'       : [ round(dz, 3) for i in range(nz) ],
        }
    },
}

# Copy namelists
#for filename in os.listdir(nametempldir):
#    if filename not in patches.keys():
#        shutil.copy('{}/{}'.format(nametempldir, filename), 
#            '{}/{}'.format(setupdirs['inputdir'], filename)) 

# Patch namelist
#for namefile in patches.keys():
#    f90nml.patch('{}/{}'.format(nametempldir, namefile), patches[namefile],
#        out_path='{}/{}'.format(setupdirs['inputdir'], namefile))


# Read namelists
namelists = {}
for filename in os.listdir(nametempldir):
    namelists[filename] = f90nml.read('{}/{}'.format(nametempldir, filename))

# Merge namelists, giving priority to namelists in the dict 'patches'.
for nmlfile in patches.keys():
    for nml in patches[nmlfile].keys():
        for var in patches[nmlfile][nml].keys():
            nmlfileL, nmlL, varL = nmlfile.lower(), nml.lower(), var.lower()
            namelists[nmlfileL][nmlL][varL] = patches[nmlfile][nml][var]
                

# Save namelists
for filename in namelists.keys():
    savename = '{}/{}'.format(setupdirs['inputdir'], filename)
    with open(savename, 'w') as namefile:
        namelists[filename].write(savename, force=True)

viscAz = gentools.sift_nested_dict('viscAz', namelists)
print('Vertical viscosity is {}'.format(viscAz))

# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Save binary input files
gentools.convert_and_save(
    {'topo'  : topo, 
     'obc_U' : obc_U, 
     'obc_T' : obc_T, 
     'init_T': init_T, 
     'dx'    : dx*np.ones((nx,))}, 
    savedir=setupdirs['inputdir'])


# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Compile
if nprun > 1:
    compilecmd = [
        "{}/tools/genmake2".format(gcmdir),
            "-optfile={}".format(optfile),
        "-enable=mnc",
        "-mpi",
        "-mods=../code/",
        "-rootdir={}".format(gcmdir),
    ]
else:
    compilecmd = [ 
        "{}/tools/genmake2".format(gcmdir),
            "-optfile={}".format(optfile),
        "-enable=mnc",
        "-mods=../code/",
        "-rootdir={}".format(gcmdir),
    ]

# Genmake and make
os.chdir(setupdirs['builddir'])

starttime = time.time()
with open('out_genmake2.txt', 'w') as genmakeout:
    with open('err_genmake2.txt', 'w') as genmakeerr:
        process = subprocess.call(compilecmd, stdout=genmakeout, stderr=genmakeerr)
print('Genmake time: {:3f} s'.format(time.time()-starttime))

if process is not 0:
    with open('err_genmake2.txt', 'r') as errfile:
        for line in errfile:
            print(line.rstrip('\n'))
    raise RuntimeError('Genmake2 failed.')

starttime = time.time()
with open('out_makedepend.txt', 'w') as makedependout:
    with open('err_makedepend.txt', 'w') as makedependerr:
        process = subprocess.call(['make', 'depend'], 
            stdout=makedependout, stderr=makedependerr)
print('Make depend time: {:3f} s'.format(time.time()-starttime))

if process is not 0:
    with open('err_makedepend.txt', 'r') as errfile:
        for line in errfile:
            print(line.rstrip('\n'))
    raise RuntimeError('make depend failed.')

starttime = time.time()
with open('out_make.txt', 'w') as makeout:
    with open('err_make.txt', 'w') as makeerr:
        process = subprocess.call(['make', '-j{}'.format(npmake)], 
            stdout=makeout, stderr=makeerr)
print('Make time: {:3f} s'.format(time.time()-starttime))

if process is not 0:
    with open('err_make.txt', 'r') as errfile:
        for line in errfile:
            print(line.rstrip('\n'))
    raise RuntimeError('make failed.')

# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Run
for filename in os.listdir(setupdirs['rundir']):
    os.remove('{}/{}'.format(setupdirs['rundir'], filename))

for filename in os.listdir(setupdirs['inputdir']):
    os.symlink(
        '{}/{}'.format(setupdirs['inputdir'], filename),
        '{}/{}'.format(setupdirs['rundir'], filename) 
    )
os.chdir(setupdirs['rundir'])

starttime = time.time()
with open('./output.txt', 'w') as output:
    _, stopmsg = subprocess.Popen('{}/mitgcmuv'.format(setupdirs['builddir']), 
        stdout=output, stderr=subprocess.PIPE).communicate()
print('Run time: {:.3f}'.format(time.time() - starttime))

if 'STOP NORMAL END' not in stopmsg:
    print(stopmsg)
