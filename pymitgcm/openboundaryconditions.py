import numpy as np

from . import gcmutils




class OpenBoundaryCondition:
    def __init__(self, model, edge, nt=1, dt=None, copyic=None, 
        orlanski=False):
        """ Instantiate an open boundary condition for an MITgcm model. 

        Args:
            model: The pymitgcm.Model to which the obc's belong.
    
            edge (str): The edge of the obcs (either south, north, west, 
                or east).

            nt (int): The number of time-points on the obc.

            dt (float): The time-interval between open boundary conditions.
            
            ic: A pymitgcm initial condition to be continued into the open 
                boundary.
        """

        filenameform = 'obc_{}_{}.bin'
        namelistform = 'ob{}{}file'
        varlist = ['U', 'V', 'W', 'T', 'S']

        self.nt = nt
        if nt > 1: 
            if dt is None:
                raise ValueError("Open boundary condition dt must be set if "
                    "nt > 1. dt is the interval between open boundary time "
                    "slices.")
            else:
                self.dt = dt
                

        self.namelistnames, self.filenames = {}, {}
        for var in varlist:
            self.namelistnames[var] = namelistform.format(edge[0], var)
            self.filenames[var] = filenameform.format(edge, var)

        if edge is 'south':
            self.n = model.x
            self.nn = model.nx
            # MITgcm indices
            self.I = np.arange(1, model.nx+1)
            self.J = np.ones((model.nx,), dtype=np.int64)
        elif edge is 'north':
            self.n = model.x
            self.nn = model.nx
            # MITgcm indices
            self.I = np.arange(1, model.nx+1)
            self.J = model.ny * np.ones((model.nx,), dtype=np.int64)
        elif edge is 'west':
            self.n = model.y
            self.nn = model.ny
            # MITgcm indices
            self.I = np.ones((model.ny,), dtype=np.int64)
            self.J = np.arange(1, model.ny+1)
        elif edge is 'east':
            self.n = model.y
            self.nn = model.ny
            # MITgcm indices
            self.I = model.nx * np.ones((model.ny,), dtype=np.int64)
            self.J = np.arange(1, model.ny+1)

        self.model = model
        self.edge = edge
        self.nz = model.nz
        self.fields = {}
        self.orlanski = orlanski

        if copyic is not None:
            self.copy_ic(copyic)


    def save(self, savepath):
        """ Save open boundary condition as binary file in MITgcm format. """
    
        # Build savevars dictionary
        savevars = {}

        for varname, var in self.fields.items():
            savevars[self.filenames[varname]] = var

        gcmutils.savegcminput(savevars, savepath=savepath)


    def add_uvel(self, U):
        """ Set the eastward velocity of the open boundary condition. """
        if 'U' not in self.fields.keys():
            self.fields['U'] = np.zeros((self.nn, self.nz, self.nt)) + U
        else:
            self.fields['U'] = self.fields['U'] + U
        self.fields['U'] = self.fields['U'].squeeze()


    def add_vvel(self, V):
        """ Set the northward velocity of the open boundary condition. """
        if 'V' not in self.fields.keys():
            self.fields['V'] = np.zeros((self.nn, self.nz, self.nt)) + V
        else:
            self.fields['V'] = self.fields['V'] + V
        self.fields['V'] = self.fields['V'].squeeze()


    def add_wvel(self, W):
        """ Set the northward velocity of the open boundary condition. """
        if 'W' not in self.fields.keys():
            self.fields['W'] = np.zeros((self.nn, self.nz, self.nt)) + W
        else:
            self.fields['W'] = self.fields['W'] + W
        self.fields['W'] = self.fields['W'].squeeze()


    def add_theta(self, T):
        """ Set the temperature distribution the open boundary condition. """
        if 'T' not in self.fields.keys():
            self.fields['T'] = np.zeros((self.nn, self.nz, self.nt)) + T
        else:
            self.fields['T'] = self.fields['T'] + T
        self.fields['T'] = self.fields['T'].squeeze()


    def add_salt(self, S):
        """ Set the salinity distribution the open boundary condition. """
        if 'S' not in self.fields.keys():
            self.fields['S'] = np.zeros((self.nn, self.nz, self.nt)) + S
        else:
            self.fields['S'] = self.fields['S'] + S
        self.fields['S'] = self.fields['S'].squeeze()
        

    def copy_ic(self, ic):
        """ Continue an initial condition into the open boundary """

        for fieldname, field in ic.fields.items():
            self.fields[fieldname] = field[self.I-1, self.J-1, :]

            #self.fields[fieldname] = (self.fields[fieldname]*
            #    np.ones((self.nn, self.nz, self.nt)))
