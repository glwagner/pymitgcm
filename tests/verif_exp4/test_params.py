""" Test changing parameter values. """

import os, sys; sys.path.append('../../')
sys.path.append('.')
import numpy as np
import matplotlib.pyplot as plt
import pymitgcm

ntimesteps = 100

def test_ntimesteps():
  # ----------------------------------------------------------------------------- 
  # Test extrating  and changing the setup's ntimesteps parameter
  oldntimesteps = setup.getparam('ntimesteps')
  print("Original ntimesteps: {}".format(oldntimesteps))

  # if they are the same, then the test doesn't actually test anything
  assert oldntimesteps != ntimesteps

  setup.setparam('ntimesteps', ntimesteps)
  newntimesteps = setup.getparam('ntimesteps')
  print("New ntimesteps: {}".format(newntimesteps))

  assert newntimesteps == ntimesteps

  # Return ntimesteps to its former value
  setup.setparam('ntimesteps', oldntimesteps)
