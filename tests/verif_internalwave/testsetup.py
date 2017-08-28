""" Test compiling and running a setup. """

import os, sys; sys.path.append('../../')
import numpy as np
import matplotlib.pyplot as plt
import pymitgcm

# Number of time steps
ntimesteps = 1000

# Paths and parameters custom to Greg's system
gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'
setuppath = os.path.join(os.getcwd(), '../../templates/verif_internalwave')
optfile = 'neve'

# Create and compile setup
setup = pymitgcm.Setup(setuppath)
setup.compilesetup(gcmpath, optfile=optfile)

# Get current ntimesteps
oldntimesteps = setup.getparam('ntimesteps')
print('Old ntimesteps: {}'.format(oldntimesteps))

# Set new ntimesteps
setup.setparam('ntimesteps', ntimesteps)
newntimesteps = setup.getparam('ntimesteps')
print('New ntimesteps: {}'.format(newntimesteps))

# Run the setup and return ntimesteps to its former value
setup.runsetup(overwrite=True)
setup.setparam('ntimesteps', 100)
