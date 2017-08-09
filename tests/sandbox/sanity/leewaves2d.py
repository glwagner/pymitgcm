import os
import numpy as np

from ..models                   import RectangularModel
from ..initialconditions        import InitialCondition
from ..openboundaryconditions   import OpenBoundaryCondition


class TwoDimLeeWaveModel(RectangularModel):
    def __init__(self, nx, nz, Lx, Lz, Ly=1.0, N=1e-3, U0=0.0, gcmdir=None):
        """ Instantiate a two-dimensional MITgcm model for lee waves. """

        RectangularModel.__init__(self, nx, 1, nz, Lx, Ly, Lz, gcmdir=gcmdir)

        self.templatedir = os.path.abspath('../../templates/openboundary2d')

        self.phys['gravity'] = 9.81
        self.phys['talpha'] = 2e-4
        self.phys['sbeta'] = 0.0

        self.ic = InitialCondition(self)
        self.ic.set_constant_stratification(N=N, 
            g=self.phys['gravity'], alpha=self.phys['talpha'])    

        self.obcs = {}
        self.obcs['east'] = OpenBoundaryCondition(self, 'east')
        self.obcs['west'] = OpenBoundaryCondition(self, 'west')

        self.obcs['east'].continue_ic(self.ic)
        self.obcs['west'].continue_ic(self.ic)
        self.obcs['east'].add_uvel(U=U0)
        self.obcs['west'].add_uvel(U=U0)



def constant_stratification(model, N=1.0e-3, g=9.81, alpha=2e-4):

    model.N = N
    model.g = g
    model.alpha = alpha

    model.Tz = N**2.0 / (g*alpha)
    model.Tref = model.Tz*model.z - (model.Tz*model.z).mean()
    model.init_T = (
        np.ones((model.nx,))[:, np.newaxis] * model.Tref[np.newaxis, :] )




def add_gaussian_topo(model, H=None, L=None):
    if H is None: H=0.6*model.Lz
    if L is None: L=0.1*model.Lx

    model.topo = -model.Lz + H*np.exp( 
        -(model.x-0.5*model.Lx)**2.0 / (2.0*L**2.0) )




def uniform_flow(model, U0=0.02):

    model.obcs = {
        'east': OpenBoundaryCondition(model, 'east'), 
        'west': OpenBoundaryCondition(model, 'west'),
    }

    model.obcs['east'].set_U(U0 * np.ones((model.nz,)))
    model.obcs['east'].set_T(init_T[0, :].squeeze())

    model.obcs['west'].set_U(U0 * np.ones((model.nz,)))
    model.obcs['west'].set_T(init_T[0, :].squeeze())
