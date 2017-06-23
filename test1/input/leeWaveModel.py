import os
import numpy as np
import f90nml 

# Numerical parameters
dt = 600    # Timestep in secs
nt = 100    # Number of timesteps

nx = 256
nz = 256

# Grid, centered around x=0
Lx, Lz = 100.0e3, 4000.0
dx, dz =   Lx/nx,   Lz/nz

x = np.arange(-Lx/2.0, Lx/2.0, dx)
z = np.arange(-Lz, 0, dz)
X, Z = np.meshgrid(x, z)

# Rotation, stratification, gravity, thermal expansion coeff.
f0 = 1.0e-4
N0 = 20.0*f0
g  = 9.81
alpha = 2.0e-4

# Viscosity
nuz = 1.0e-2        # Vertical Laplacian viscosity 
nuh = 1.0e4         # Horizontal Laplacian viscosity
nu4 = 1.0e8         # Small Horizontal hyperviscosity

kapz = 1.0e-5       # Vertical diffusivity
kaph = 1.0e3        # Horizonztal diffusivity 

# Constant temperature gradient corresponding to N0
Tz = N0**2.0 / (g*alpha)
Tref = Tz*z - (Tz*z).mean()

# Gaussian bump
hbump, lbump = 1000.0, 10.0e3
topo = -Lz + hbump*np.exp( -X**2.0 / (2.0*lbump**2.0) )

# Open boundaries
OB_Ieast = 1
OB_Iwest = -1
