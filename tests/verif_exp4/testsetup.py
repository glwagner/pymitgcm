""" Test compiling and running a setup. """

import os, sys; sys.path.append('../../')
import numpy as np
import matplotlib.pyplot as plt
import pymitgcm

templatepath = os.path.join(os.getcwd(), '../../templates/verif_exp4')
optfilename = 'neve'

if os.path.exists(os.path.join(os.getcwd(), '../../MITgcm')):
    gcmpath = os.path.join(os.getcwd(), '../../MITgcm')
else:
    gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'



setup = pymitgcm.Setup(templatepath)

if not setup.hasexecutable():
    setup.compilesetup(gcmpath, optfilename=optfilename)

ntimesteps = setup.getparam('ntimesteps')

setup.setparam('ntimesteps', 1000)
setup.runsetup(overwrite=True)

# Return ntimesteps to its former value
setup.getparam('ntimesteps', ntimesteps)

