""" Test creating a model, changing the resolution, and 
compiling, running, and looking at the output. """

import os, shutil, sys; sys.path.append('../../')
import numpy as np; from numpy import pi
import matplotlib.pyplot as plt
import pymitgcm

gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'
templatepath = os.path.join(os.getcwd(), '../../templates/leewave')
optfile = 'neve'


# Parameters
nx, nz = 400, 100
Lx, Lz = 1e6, 2000
dt, nt = 100, 100
compile = True

f, N     = 1e-4, 1e-2
g, alpha = 9.81, 2e-4

params = {
    'ntimesteps'     : nt, 
    'deltat'         : dt, 
    'dumpfreq'       : 100*dt, 
    'viscah'         : 1e-1,
    'diffkht'        : 1e-1,
    'viscaz'         : 1e-3,
    'diffkzt'        : 0.0,
    'visca4'         : 1e8
}

m = pymitgcm.RectangularModel(nx, 1, nz, Lx, 1e3, Lz, nprun=1,
    templatepath=templatepath, gcmpath=gcmpath)

m.eqns['nonhydrostatic'] = False
m.eqns['rigidlid']       = False
m.eqns['saltstepping']   = False


# Initial condition
u0 = 0.04
m.init_ic()
m.ic.set_constant_stratification(N=N, T0=2.0)
m.ic.add_uvel(u0)


# OBCS
m.init_obcs(['east', 'west'], copyic=m.ic, orlanski=True)
m.init_obcs_sponge(50, urelaxbound=dt/10, vrelaxbound=0.0)


# Topography
lbump, hbump = Lx/20, 0.2*Lz
x0, y0 = m.X.mean(), m.Y.mean()
m.topo = -Lz + hbump*np.exp( -((m.X-x0)**2.0+(m.Y-y0)**2.0)/(2*lbump**2.0) )


# Shenanigans
m.gensetup(cleansetup=False, params=params)
m.run(overwrite=True, 
    compile=compile, npmake=2, optfile='neve', remove=['*.nc'])
