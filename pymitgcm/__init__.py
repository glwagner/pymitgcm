# ----------------------------------------------------------------------------- 
# pymitgcm
#
# A python code for interacting with, building, running, and analyzing
# MITgcm models.
#
# ----------------------------------------------------------------------------- 

from . import initialconditions
from . import openboundaryconditions


# models.py contains the core Setup and Model classes that interface with
# MITgcm input files and code
from .models                 import Setup, Model, RectangularModel


# initialconditions.py and openboundaryconditions.py contain classes that 
# help set and save initial and open boundary conditions
from .initialconditions      import InitialCondition
from .openboundaryconditions import OpenBoundaryCondition


# gcmutils.py defines miscellaneous helper functions
# that are used to save input, call MITgcm command-line utilities,
# interact with MITgcm namelists and source files, etc.
from .gcmutils import (
    rungcm, compile, genmake, makedepend, make, 
    readsizevars, changesizevar, changesizevars, siftdict, savegcminput, 
    truncate, getcompacttiling
)


# analysis.py defines simple functions for plotting and analyzing the 
# output of simple MITgcm simulations.
from .analysis import (
    quicklook, getoutputiters, getoutputstate,
    getgrid, rdmds, readmeta
)
