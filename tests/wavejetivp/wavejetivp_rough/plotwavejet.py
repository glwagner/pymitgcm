import os, sys; sys.path.append('../../../')
import numpy as np; from numpy import pi
import matplotlib.pyplot as plt
import pymitgcm

# Create setup object in current directory
setup = pymitgcm.Setup(os.getcwd())
f0 = setup.getparam('f0')

nt = pymitgcm.getoutputiters(setup.run)

# Kilometers
x, y, z, t, xg, yg = pymitgcm.getgrid(setup.run)
x, y, xg, yg  = 1e-3*x, 1e-3*y, 1e-3*xg, 1e-3*yg

# Get initial condition
t0, vars0 = pymitgcm.getoutputstate(setup.run, iter=0)

# Plot limits
umax, wmax, Tmax = 0.01, 4e-5, 1e-6
nlevs = 20
ulevs = np.linspace(-umax, umax, nlevs+1)
wlevs = np.linspace(-wmax, wmax, nlevs+1)
Tlevs = np.linspace(-Tmax, Tmax, nlevs+1)

# x,z slice
fig, axs = plt.subplots(nrows=2, ncols=1)
im, cb = [None]*2, [None]*2
for i in range(nt):

    axs[0].cla(); axs[1].cla()

    t, vars = pymitgcm.getoutputstate(setup.run, iter=i)
    w = vars['w'][:, setup.ny//2, :].squeeze()
    u = vars['u'][:, setup.ny//2, :].squeeze()
    Tp = (vars['T'][:, setup.ny//2, :].squeeze() 
        - vars0['T'][:, setup.ny//2, :].squeeze())

    im[0] = axs[0].contourf(xg, z, u, 
        cmap='RdBu_r', levels=ulevs)
    #im[0] = axs[0].contourf(x, z, Tp, 
    #    cmap='RdBu_r', levels=Tlevs)
    im[1] = axs[1].contourf(x, z, w, 
        cmap='RdBu_r', levels=wlevs)

    titles = ["$u$ (m/s)", "$w$ (m/s)"]
    message = "$t = {:.4f}$ inertial periods".format(t*f0/(2*pi))
    axs[0].text(0.0, 1.03, message, transform=axs[0].transAxes)
    axs[0].text(1.0, 1.03, titles[0], transform=axs[0].transAxes, 
        HorizontalAlignment='right')
    axs[1].text(1.0, 1.03, titles[1], transform=axs[1].transAxes,
        HorizontalAlignment='right')


    axs[0].set_xlim(x[[0, -1]])
    axs[1].set_xlim(x[[0, -1]])
    
    axs[0].set_xticklabels([])
    axs[0].tick_params(length=0)
    axs[1].tick_params(length=0)
    axs[1].set_xlabel("$x$ (km)")

    axs[0].set_ylabel("$z$ (m)")
    axs[1].set_ylabel("$z$ (m)")

    if i is 0:
        
        uticks, wticks = ulevs[::4], wlevs[::4]
        cb[0] = plt.colorbar(mappable=im[0], ax=axs[0], ticks=uticks)
        cb[1] = plt.colorbar(mappable=im[1], ax=axs[1], ticks=wticks)

        uticklabels = ["{:.1e}".format(tick) for tick in uticks]
        wticklabels = ["{:.1e}".format(tick) for tick in wticks]

        cb[0].ax.set_yticklabels(uticklabels)
        cb[1].ax.set_yticklabels(wticklabels)

        cb[0].ax.tick_params(length=0)
        cb[1].ax.tick_params(length=0)

    plt.tight_layout()
    plt.pause(0.1)

plt.show()
