import os, sys; sys.path.append('../../')

import numpy as np
import matplotlib.pyplot as plt
import pymitgcm

#gcmdir = '/data5/glwagner/Numerics/pymitgcm/MITgcm'
gcmdir = '/Users/glwagner/Software/MITgcm'

# Parameters
nx, nz = 100, 100
Lx, Ly, Lz = 13.3e3, 5e3, 200.0
hbump, lbump = 0.5*Lz, 0.1*Lx
U, N = 0.01, 1e-3

def add_gaussian_topo(model, H, L):
    model.topo = -model.Lz + H*np.exp( 
        -(model.x-0.5*model.Lx)**2.0 / (2.0*L**2.0) )




# Make model
m = pymitgcm.RectangularModel(nx, 1, nz, Lx, Ly, Lz, gcmdir=gcmdir, 
        templatedir=os.path.abspath('../../templates/openboundary2d'))
    
m.eqns['gravity'] = 9.81
m.eqns['talpha'] = 2e-4
m.eqns['sbeta'] = 0.0
m.eqns['viscaz'] = 0.002

m.init_ic()
m.init_obcs(['east', 'west'])

m.ic.set_constant_stratification(N=N, 
    g=m.eqns['gravity'], alpha=m.eqns['talpha'])    

m.obcs['east'].continue_ic(m.ic)
m.obcs['west'].continue_ic(m.ic)
m.obcs['east'].add_uvel(U)
m.obcs['west'].add_uvel(U)

add_gaussian_topo(m, hbump, lbump)

# Check sanity
fig, axs = plt.subplots(nrows=3, ncols=2)

print(m.ic.fields['T'].shape)

axs[0, 0].imshow(m.ic.fields['T'].squeeze())
axs[0, 1].plot(m.x, m.topo)
axs[1, 0].plot(m.obcs['east'].fields['T'].squeeze(), m.z)
axs[1, 1].plot(m.obcs['west'].fields['T'].squeeze(), m.z)
axs[2, 0].plot(m.obcs['east'].fields['U'].squeeze(), m.z)
axs[2, 1].plot(m.obcs['west'].fields['U'].squeeze(), m.z)

plt.show()

m.gensetup(cleansetup=True)

m.compile(optfilename='neve', npmake=8)
m.run()
