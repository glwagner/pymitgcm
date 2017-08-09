from . import initialconditions
from . import openboundaryconditions

from .models                 import Model, Setup
from .initialconditions      import InitialCondition
from .openboundaryconditions import OpenBoundaryCondition

from .gcmutils import (
    rungcm, compile, genmake, makedepend, make, 
    readsizevars, changesizevar, changesizevars, siftdict, savegcminput, 
    quicklook, truncate,
)
