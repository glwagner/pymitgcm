import os, sys; sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import pymitgcm


iters = range(0, 1000, 18)
zlev = 0

def loaditer(field, iter):
    U = 0.0
    for i in range(1, 3):
        for j in range(1, 3):
            U += pymitgcm.rdmds(
                './run/{}.{:010d}.{:03d}.{:03d}'.format(field, iter, i, j))

    return U

# Load and plot
fig = plt.figure()
for i in iters:
    U = loaditer('U', i)
    plt.pcolormesh(U[zlev])
    plt.pause(0.1)
