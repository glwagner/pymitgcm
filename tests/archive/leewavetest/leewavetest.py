import sys; sys.path.append('../../')

import matplotlib.pyplot as plt

import pymitgcm
from pymitgcm.setups import leewaves2d

gcmdir = '/data5/glwagner/Numerics/pymitgcm/MITgcm'


# Parameters
nx, nz = 200, 200
Lx, Ly, Lz = 13.3e3, 5e3, 200.0
hbump, lbump = 0.5*Lz, 0.1*Lx

# Make model
m = leewaves2d.TwoDimLeeWaveModel(nx, nz, Lx, Lz, Ly=Ly,
    U0=0.01, N=1e-3, gcmdir=gcmdir)
    
leewaves2d.add_gaussian_topo(m, H=hbump, L=lbump)

m.gensetup(cleansetup=True, namepatch={'data': {'parm01': {'viscaz': 0.002}}})

#fig = plt.figure()
#plt.plot(m.ic.fields['T'].mean(axis=(0, 1)), m.z)
#
#fig = plt.figure()
#plt.plot(m.x, m.topo.squeeze())
#
#plt.show()

m.compile(optfilename='sverdrup_glwagner', npmake=8)
m.run()
