import os, shutil, sys; sys.path.append('../../')
import numpy as np
import matplotlib.pyplot as plt
import pymitgcm

setup = pymitgcm.Setup()
nt = pymitgcm.getoutputiters(setup.run)

x, y, z, t, xg, yg = pymitgcm.getgrid(setup.run)

# Kilometers
x, y, xg, yg  = 1e-3*x, 1e-3*y, 1e-3*xg, 1e-3*yg


# x,z slice
figxz, axxz = plt.subplots(nrows=2, ncols=1, sharex=True)
imxz = [None, None]
for i in range(nt):
    vars = pymitgcm.getoutputstate(setup.run, iter=i)

    axxz[0].cla()
    axxz[1].cla()

    imxz[0] = axxz[0].pcolormesh(xg, z, vars['u'][:, setup.ny//2, :].squeeze(), 
        cmap='RdBu_r', vmin=0.0, vmax=0.5)

    imxz[1] = axxz[1].pcolormesh(xg, z, vars['w'][:, setup.ny//2, :].squeeze(),
        cmap='RdBu_r', vmin=-0.02, vmax=0.02)

    if i is 0:
        plt.colorbar(mappable=imxz[0], ax=axxz[0])
        plt.colorbar(mappable=imxz[1], ax=axxz[1])

    plt.pause(0.1)

plt.show()
