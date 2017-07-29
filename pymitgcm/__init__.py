from . import initialconditions
from . import openboundaryconditions

from .models                 import RectangularModel
from .initialconditions      import InitialCondition
from .openboundaryconditions import OpenBoundaryCondition

from .gcmutils import (
    rungcm, compile, genmake, makedepend, make, 
    changesizevar, changesizevars, siftdict, savegcminput, truncate,
)
