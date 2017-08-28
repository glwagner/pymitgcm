""" Test creating a model, changing the resolution, and 
compiling, running, and looking at the output. """

import os, shutil, sys; sys.path.append('../../')
import numpy as np; from numpy import pi
import matplotlib.pyplot as plt
import pymitgcm

gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'
templatepath = os.path.join(os.getcwd(), '../../templates/verif_exp4')
optfile = 'neve'


# Parameters
nx, ny, nz = 160,    80,   16
Lx, Ly, Lz = 4e5,   2e5, 4500

f, N     = 1e-4, 1e-3
g, alpha = 9.81, 2e-4

m = pymitgcm.RectangularModel(nx, ny, nz, Lx, Ly, Lz, nprun=1,
    templatepath=templatepath, gcmpath=gcmpath)


# Packages
m.pkgs['useobcs']     = True
m.pkgs['useptracers'] = False
m.pkgs['usemnc']      = True
m.pkgs['userbcs']     = False

params = {
    'ntimesteps'     : 1000, 
    'deltat'         : 100, 
    'basetime'       : 1000, 
    'dumpfreq'       : 1000, 
}


# Topography
lbump, hbump = 25e3, 0.6*Lz
x0, y0 = m.X.mean(), m.Y.mean()
m.topo = -Lz + hbump*np.exp( -((m.X-x0)**2.0+(m.Y-y0)**2.0)/(2*lbump**2.0) )


# OBCS
m.init_obcs(['north', 'south', 'east', 'west'])

u0, w0 = 0.25, 0.001
du = 0.00*u0

w1 = np.zeros((ny, nz, 1)) + w0*np.sin(m.z*pi/m.Lz)
s1 = 36*np.ones((ny, nz, 1))

umerid = np.concatenate((u0*np.ones((nx, nz, 1)), 
    u0*np.ones((nx, nz, 1))), axis=2)
uwest  = np.concatenate(((u0+du)*np.ones((ny, nz, 1)), 
    (u0-du)*np.ones((ny, nz, 1))), axis=2)
ueast  = np.concatenate(((u0+du)*np.ones((ny, nz, 1)), 
    (u0-du)*np.ones((ny, nz, 1))), axis=2)
wwest  = np.concatenate((w1, w1), axis=2)
swest  = np.concatenate((s1, s1), axis=2)


m.obcs['north'].fields['U'] = umerid
m.obcs['south'].fields['U'] = umerid
m.obcs['east'].fields['U'] = ueast
m.obcs['west'].fields['U'] = uwest
m.obcs['west'].fields['W'] = wwest
m.obcs['west'].fields['S'] = swest


# Shenanigans
m.gensetup(cleansetup=False, params=params)
m.run(overwrite=True, 
    compile=True, npmake=2, optfile='neve', remove=['*.nc'])

# Initial condition


