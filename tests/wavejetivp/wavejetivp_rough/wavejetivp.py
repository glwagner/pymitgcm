""" Test creating a model, changing the resolution, and 
compiling, running, and looking at the output. """

import os, shutil, sys; sys.path.append('../../..')
import numpy as np; from numpy import pi
import matplotlib.pyplot as plt
import pymitgcm

gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'
templatepath = os.path.join(os.getcwd(), '../../../templates/wavejetivp')
optfile = 'neve'


# Parameters
nx, ny, nz = 200,   1,  100
Lx, Ly, Lz = 1e6, 1e3, 4000

f, N     = 1e-4, 0e-4
g, alpha = 9.81, 2e-4

dt = 0.01 * 2*pi/f
nt = 1 * int(np.round(2*pi/(f*dt)))

m = pymitgcm.RectangularModel(nx, ny, nz, Lx, Ly, Lz, nprun=1,
    templatepath=templatepath, gcmpath=gcmpath)

params = {
    'f0'             : f,
    'beta'           : 0.0,
    'gravity'        : g,
    'talpha'         : alpha,
    'viscaz'         : 0.001,
    'viscah'         : 1e3,
    'diffkzt'        : 1e-5,
    'diffkht'        : 1e3,
    'diffkzs'        : 1e-5,
    'diffkhs'        : 1e3,
    'nonhydrostatic' : True,
    'ntimesteps'     : nt,
    'deltat'         : dt, 
    'dumpfreq'       : dt #2*pi/f
}

# Initial condition: constant stratification, near-inertial wave
m.init_ic()
m.ic.set_constant_stratification(N=N, g=g, alpha=alpha, T0=2.0)


# Wave and jet parameters
U0, du = 0.01, 400
Ro, dv, x0 = 0.1, Lx/20, Lx/2

u0 = U0*np.exp(m.Z/du)

v0 = Ro*f*dv*np.exp(-(m.X-x0)**2/(2*dv**2))
v0 -= v0.mean()

m.ic.add_uvel(u0)
m.ic.add_vvel(v0)

# Add rough topography
dh, lh, x0 = 1.0*Lz/100, Lx/20, Lx/2
m.topo = -Lz + dh*np.exp(-(m.X-x0)**2/(2*lh**2))

# Shenanigans
m.gensetup(cleansetup=False, params=params)
m.run(overwrite=True, 
    compile=False, npmake=2, optfile='neve', remove=['*.nc'])
