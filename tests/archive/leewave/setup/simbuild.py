from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

inputdir = '../input'
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
U0 = 0.01
obc_U = U0 * np.ones((nz,))
obc_T = init_T[0, :].squeeze()

# Plot
if plotting:
    fig1 = plt.figure()
    plt.plot(x*1.0e-3, topo)

    fig2 = plt.figure()
    plt.plot(init_T[0, :], z)
    #
    plt.show()

def convert_and_save(vars, savedir):

    # Convert to proper format
    savetype = '>f8'

    for var in vars.keys():
        vars[var] = vars[var].astype(savetype) 
        varshape = vars[var].shape

        with open('{}/{}.bin'.format(savedir, var), 'wb') as file:
            if len(varshape) is 1:
                vars[var].tofile(file)            
            elif len(varshape) is 2:
                for k in range(varshape[1]):
                    vars[var][:, k].tofile(file)

vars = {}
vars['topo']   = topo
vars['obc_U']  = obc_U
vars['obc_T']  = obc_T
vars['init_T'] = init_T
vars['dx']     = dx * np.ones((nx,))

convert_and_save(vars, inputdir)
