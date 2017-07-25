from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

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
hbump = 0.2*Lz
lbump = 0.1*Lx
topo = -Lz + hbump * np.exp( -(x-0.5*Lx)**2.0 / (2.0*lbump**2.0) )
topo[-1] = 0.0

# Open boundaries
U0 = 0.05
obc_U = U0 * np.ones((nz,))
obc_T = init_T[0, :].squeeze()

# dx
del_X = dx * np.ones((nx,))

# Plot
if plotting:
    fig1 = plt.figure()
    plt.plot(x*1.0e-3, topo)

    fig2 = plt.figure()
    plt.plot(init_T[0, :], z)
    #
    plt.show()


# Convert to proper format
savetype = '>f8'

topo = topo.astype(savetype)
del_X = del_X.astype(savetype)
obc_U = obc_U.astype(savetype)
obc_T = obc_T.astype(savetype)
init_T = init_T.astype(savetype)

# Save
with open('topo.bin', 'wb') as file:
    topo.tofile(file)

with open('init_T.bin', 'wb') as file:
    for k in range(nz):
        init_T[:, k].tofile(file)

with open('obc_U.bin', 'wb') as file:
    obc_U.tofile(file)

with open('obc_T.bin', 'wb') as file:
    obc_T.tofile(file)

with open('dx.bin', 'wb') as file:
    del_X.tofile(file)
