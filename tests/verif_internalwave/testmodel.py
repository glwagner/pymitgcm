""" Test creating a model, changing the resolution, and 
compiling, running, and looking at the output. """

import sys; sys.path.append('../../')
import numpy as np
import matplotlib.pyplot as plt
import pymitgcm

templatepath = '/Users/glwagner/Numerics/pymitgcm/templates/internalwaveverif'
gcmpath = '/Users/glwagner/Numerics/pymitgcm/MITgcm'
optfile = 'neve'

m = pymitgcm.Model(templatepath=templatepath, gcmpath=gcmpath)

print(m.template.nx)
print(m.template.nz)
print(m.template.nprun)
