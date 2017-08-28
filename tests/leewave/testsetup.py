""" Test compiling and running MITgcm's exp4 verification experiment. """

import os, sys; sys.path.append('../../')
import numpy as np
import matplotlib.pyplot as plt
import pymitgcm

ntimesteps = 100

gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'
templatepath = os.path.join(os.getcwd(), '../../templates/verif_exp4')
optfile = 'neve'

# Create and compile the setup
setup = pymitgcm.Setup(templatepath)
setup.compilesetup(gcmpath, optfile=optfile)


# ----------------------------------------------------------------------------- 
# Test extrating  and changing the setup's ntimesteps parameter
oldntimesteps = setup.getparam('ntimesteps')
print("Original ntimesteps: {}".format(oldntimesteps))

setup.setparam('ntimesteps', ntimesteps)
newntimesteps = setup.getparam('ntimesteps')
print("New ntimesteps: {}".format(newntimesteps))

setup.runsetup(overwrite=True)

# Return ntimesteps to its former value
setup.setparam('ntimesteps', oldntimesteps)
