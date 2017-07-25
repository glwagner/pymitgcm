import sys; sys.path.append('../../')

import pymitgcm
from pymitgcm.setups import leewaves2d

# Parameters
nx, nz = 100, 100
Lx, Ly, Lz = 20e3, 5e3, 200.0
U0 = 0.02
hbump, lbump = 0.6*Lz, 0.1*Lx


# Make model
m = leewaves2d.TwoDimLeeWaveModel(nx, nz, Lx, Lz, Ly=Ly, U0=U0, 
    gcmdir='/Users/glwagner/Numerics/MITgcm')

leewaves2d.add_gaussian_topo(m, H=hbump, L=lbump)

m.gensetup(cleansetup=True)
m.compile(optfilename='neve', npmake=2)
m.run()
