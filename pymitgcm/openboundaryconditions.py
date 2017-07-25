import numpy as np

from . import gcmutils




class OpenBoundaryCondition:
    def __init__(self, model, edge):
        """ Instantiate an open boundary condition for an MITgcm model. """


        filenameform = 'obc_{}_{}.bin'
        namelistform = 'ob{}{}file'
        varlist = ['U', 'V', 'W', 'T', 'S']

        self.namelistnames, self.filenames = {}, {}
        for var in varlist:
            self.namelistnames[var] = namelistform.format(edge[0], var)
            self.filenames[var] = filenameform.format(edge, var)

        if edge is 'south':
            self.n = model.x
            self.nn = model.nx
            self.I = np.arange(0, model.nx)
            self.J = np.zeros((model.nx,), dtype=np.int64)
        elif edge is 'north':
            self.n = model.x
            self.nn = model.nx
            self.I = np.arange(0, model.nx)
            self.J = (model.ny-1) * np.ones((model.nx,), dtype=np.int64)
        elif edge is 'east':
            self.n = model.y
            self.nn = model.ny
            self.I = np.zeros((model.ny,), dtype=np.int64)
            self.J = np.arange(0, model.ny)
        elif edge is 'west':
            self.n = model.y
            self.nn = model.ny
            self.I = (model.nx-1) * np.ones((model.ny,), dtype=np.int64)
            self.J = np.arange(0, model.ny)


        self.model = model
        self.edge = edge
        self.nz = model.nz
        self.fields = {}


    def save(self, savedir):
        """ Save open boundary condition as binary file in MITgcm format. """
    
        # Build savevars dictionary
        savevars = {}

        for varname, var in self.fields.items():
            savevars[self.filenames[varname]] = var

        gcmutils.savegcminput(savevars, savedir=savedir)


    def set_uvel(U):
        """ Set the eastward velocity of the open boundary condition. """
        self.fields['U'] = U


    def set_vvel(V):
        """ Set the northward velocity of the open boundary condition. """
        self.fields['V'] = V


    def set_theta(T):
        """ Set the temperature distribution the open boundary condition. """
        self.fields['T'] = T


    def set_salt(S):
        """ Set the salinity distribution the open boundary condition. """
        self.fields['S'] = S


    def continue_ic(ic):
        """ Continue an initial condition into the open boundary """

        for fieldname, field in ic.fields.items():
            self.fields[fieldname] = field[self.I, self.J, :].squeeze()
