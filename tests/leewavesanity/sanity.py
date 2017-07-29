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
U, N = 0.02, 1e-3

def add_gaussian_topo(model, H, L):
    model.topo = -model.Lz + H*np.exp( 
        -(model.x-0.5*model.Lx)**2.0 / (2.0*L**2.0) )

# Make model
m = pymitgcm.RectangularModel(nx, 1, nz, Lx, Ly, Lz, gcmdir=gcmdir, 
    templatedir='../../templates/leewavesanity')

# Timestepper
m.tstepper['deltat']     = 1.0
m.tstepper['ntimesteps'] = 100
m.tstepper['pchkptfreq'] = 100
m.tstepper['chkptfreq']  = 100
m.tstepper['dumpfreq']   = ntimesteps//10

m.init_ic()
m.ic.set_constant_stratification(N=N, 
    g=9.81, alpha=0.0002)

m.init_obcs(['east', 'west'])
m.obcs['east'].copy_ic(m.ic)
m.obcs['west'].copy_ic(m.ic)
m.obcs['east'].add_uvel(U)
m.obcs['west'].add_uvel(U)

fig, axs = plt.subplots(ncols=2)
axs[0].plot(m.obcs['east'].fields['T'])
axs[1].plot(m.ic.fields['T'][0, 0, :])
plt.show()

m.gensetup(cleansetup=True, nprun=1)
#m.compile(optfilename='neve', npmake=2)

m.run()

    
#m.eqns['gravity'] = 9.81
#m.eqns['talpha'] = 2e-4
#m.eqns['sbeta'] = 0.0
#m.eqns['viscaz'] = 0.002
#
#
#
#
#add_gaussian_topo(m, hbump, lbump)


#plt.show()
