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
        elif edge is 'west':
            self.n = model.y
            self.nn = model.ny
            self.I = np.zeros((model.ny,), dtype=np.int64)
            self.J = np.arange(0, model.ny)
        elif edge is 'east':
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


    def add_uvel(self, U):
        """ Set the eastward velocity of the open boundary condition. """
        if 'U' not in self.fields.keys():
            self.fields['U'] = np.zeros((self.nn, self.nz)) + U
        else:
            self.fields['U'] = self.fields['U'] + U
        self.fields['U'] = self.fields['U'].squeeze()


    def add_vvel(self, V):
        """ Set the northward velocity of the open boundary condition. """
        if 'V' not in self.fields.keys():
            self.fields['V'] = np.zeros((self.nn, self.nz)) + V
        else:
            self.fields['V'] = self.fields['V'] + V
        self.fields['V'] = self.fields['V'].squeeze()


    def add_theta(self, T):
        """ Set the temperature distribution the open boundary condition. """
        if 'T' not in self.fields.keys():
            self.fields['T'] = np.zeros((self.nn, self.nz)) + T
        else:
            self.fields['T'] = self.fields['T'] + T
        self.fields['T'] = self.fields['T'].squeeze()


    def add_salt(self, S):
        """ Set the salinity distribution the open boundary condition. """
        if 'S' not in self.fields.keys():
            self.fields['S'] = np.zeros((self.nn, self.nz)) + S
        else:
            self.fields['S'] = self.fields['S'] + S
        self.fields['S'] = self.fields['S'].squeeze()
        

    def copy_ic(self, ic):
        """ Continue an initial condition into the open boundary """

        for fieldname, field in ic.fields.items():
            self.fields[fieldname] = field[self.I, self.J, :].squeeze()

