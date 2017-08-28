""" Test creating a model, changing the resolution, and 
compiling, running, and looking at the output. """

import os, shutil, sys; sys.path.append('../../')
import numpy as np; from numpy import pi
import matplotlib.pyplot as plt
import pymitgcm

gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'
templatepath = os.path.join(os.getcwd(), '../../templates/wavejetivp')
optfile = 'neve'


# Parameters
nx, ny, nz = 200,   1,  100
Lx, Ly, Lz = 1e6, 1e3, 4000

f, N     = 1e-4, 2e-3
g, alpha = 9.81, 2e-4

m = pymitgcm.RectangularModel(nx, ny, nz, Lx, Ly, Lz, nprun=1,
    templatepath=templatepath, gcmpath=gcmpath)

params = {
    'f0'         : f,
    'beta'       : 0.0,
    'gravity'    : g,
    'talpha'     : alpha,
    'ntimesteps' : 20*100, 
    'deltat'     : 0.01 * 2*pi/f,    # [seconds]
    'basetime'   : 1000,             # ?
    'dumpfreq'   : 2*pi/f,           # Save frequency
}

# Initial condition: constant stratification, near-inertial wave
m.init_ic()
m.ic.set_constant_stratification(N=N, g=g, alpha=alpha, T0=2.0)

# NIW parameters: initial amplitude and decay scale
u0, du = 0.01, 400
m.ic.add_uvel(u0*np.exp(m.Z/du))

# Jet parameters
Ro, dv, x0 = 0.1, Lx/20, Lx/2

v0 = Ro*f*dv*np.exp(-(m.X-x0)**2/(2*dv**2))
v0 -= v0.mean()
m.ic.add_vvel(v0)

# Shenanigans
m.gensetup(cleansetup=False, params=params)
m.run(overwrite=True, 
    compile=True, npmake=2, optfile='neve', remove=['*.nc'])
