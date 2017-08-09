from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

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
N     = 2.0e-3

# Cell centers
x = np.arange(dx, Lx+dx, dx)
z = np.arange(-0.5*dz, -Lz, -dz)

# Stratification
Tz = N**2.0 / (g*alpha)
Tref = Tz*z - (Tz*z).mean()

# Initial temperature profile
T0 = np.ones((nx,))[:, np.newaxis] * Tref[np.newaxis, :]

# Topography
hbump = 0.5*Lz
lbump = 0.05*Lx
topo = -Lz + hbump*np.exp( (x-0.5*Lx)**2.0 / (2.0*lbump**2.0) )

# Plot
fig1 = plt.figure()
plt.plot(x, topo)

fig2 = plt.figure()
plt.plot(T0, z)

plt.show()
