from __future__ import division
import os, subprocess, shutil, time
import numpy as np
import matplotlib.pyplot as plt
import f90nml

import gentools

workdir  = os.getcwd()
inputdir = '{}/input'.format(workdir)
rundir   = '{}/run'.format(workdir)
builddir = '{}/build'.format(workdir)
templdir = '{}/template'.format(workdir)

plotting = False

# Domain
Lx = 13.3e3
Ly = 5.0e3
Lz = 200.0

nx = 60
ny = 1
nz = 20

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

# Plot
if plotting:
    fig1 = plt.figure()
    plt.plot(x*1.0e-3, topo)

    fig2 = plt.figure()
    plt.plot(init_T[0, :], z)
    
    plt.show()

# Manage input directory
for filename in os.listdir(inputdir):
    os.remove('{}/{}'.format(inputdir, filename))

# Custom params
patches = { 
    'data': {
        'parm01': {
            'viscAz'  : 1.0e-3,
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
    },
}

# Copy namelists
for filename in os.listdir(templdir):
    if filename not in patches.keys():
        shutil.copy('{}/{}'.format(templdir, filename), 
            '{}/{}'.format(inputdir, filename)) 

# Patch namelist
for namefile in patches.keys():
    f90nml.patch('{}/{}'.format(templdir, namefile), patches[namefile],
        out_path='{}/{}'.format(inputdir, namefile))

# Read namelists
namelists = {}
for filename in os.listdir(inputdir):
    namelists[filename] = f90nml.read('{}/{}'.format(templdir, filename))


gentools.convert_and_save(
    {'topo'  : topo, 
     'obc_U' : obc_U, 
     'obc_T' : obc_T, 
     'init_T': init_T, 
     'dx'    : dx*np.ones((nx,))}, 
    savedir=inputdir)

viscAz = gentools.sift_nested_dict('viscAz', namelists)
print('Vertical viscosity is {}'.format(viscAz))


# Run
for filename in os.listdir(rundir):
    os.remove('{}/{}'.format(rundir, filename))

for filename in os.listdir(inputdir):
    os.symlink(
        '{}/{}'.format(inputdir, filename),
        '{}/{}'.format(rundir, filename) 
    )
    
os.chdir(rundir)

start = time.time()
with open('./output.txt', 'w') as output:
    _, stopmsg = subprocess.Popen('{}/mitgcmuv'.format(builddir), 
        stdout=output, stderr=subprocess.PIPE).communicate()
elapsed = time.time() - start

if 'STOP NORMAL END' not in stopmsg:
    print(stopmsg)

print("Run time: {:.3f} s".format(elapsed))
