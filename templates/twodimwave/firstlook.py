import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

np = 1
rundir = './run'

i = 0
grid  = Dataset("{}/grid.t{:03d}.nc".format(rundir, i+1), 'r')
state = Dataset("{}/state.{:010d}.t{:03d}.nc".format(rundir, 0, i+1), 'r')
print(state.variables.keys())
print(state.variables['Temp'])

nt = state.variables['Temp'][:, 0, 0, 0].size
T = state.variables['Temp'][0, :, 0, :]

fig = plt.figure()
plt.imshow(T)
plt.colorbar()
    
for i in range(nt):
    T = state.variables['Temp'][i, :, 0, :]
    T = T.squeeze()

    plt.clf()
    plt.imshow(T)
    plt.pause(0.1)
