import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

np = 4
rundir = './run'

# ----------------------------------------------------------------------------- 
state, ptracers, grid = [], [], []
for i in range(np):
    grid.append(
        Dataset("{}/grid.t{:03d}.nc".format(rundir, i+1), 'r'))
    ptracers.append(    
        Dataset("{}/ptracers.{:010d}.t{:03d}.nc".format(rundir, 0, i+1), 'r'))
    state.append(
        Dataset("{}/state.{:010d}.t{:03d}.nc".format(rundir, 0, i+1), 'r'))


print(ptracers[0].variables.keys())
print(state[0].variables.keys())

T = []
for i in range(np):
    T.append(state[i].variables['Temp'][-1, 0, :, :])


fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
for i in range(np):

    j = int(i/2)
    k = int(i - j*(np/2))
    print(j, k)

    axs[j, k].imshow(T[i])

#plt.colorbar(origin='lower')

test = state[0].variables['Temp'][-1, :, -1, :]
fig2 = plt.figure()
plt.imshow(test)

plt.show()
