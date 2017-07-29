import numpy as np

from . import gcmutils



class InitialCondition:
    def __init__(self, model):
        """ Instantiate an initial condition. """

        filenameform = 'init_{}.bin'

        self.namelistnames = {
            'U': 'uvelinitfile', 
            'V': 'vvelinitfile', 
            'W': 'wvelinitfile', 
            'T': 'hydrogthetafile', 
            'S': 'hydrogsaltfile', 
        }

        self.filenames = {}
        for var in self.namelistnames.keys():
            self.filenames[var] = filenameform.format(var)

        self.model = model
        self.nx, self.ny, self.nz = model.nx, model.ny, model.nz 
        self.fields = {}


    def save(self, savedir):
        """ Save initial condition as binary file in MITgcm format. """
    
        # Build savevars dictionary
        savevars = {}

        for varname, var in self.fields.items():
            savevars[self.filenames[varname]] = var

        gcmutils.savegcminput(savevars, savedir=savedir)
            

    def set_theta(self, T):
        """ Set potential temperature profile of the initial condition. """ 
        self.fields['T'] = T


    def set_salt(self, S):
        """ Set the salinity profile of the initial condition. """ 
        self.fields['S'] = S


    def add_uvel(self, U):
        """ Add an eastward velocity to the initial condition. """ 

        U = U*np.ones((self.nx, self.ny, self.nz))  # Broadcast to proper size

        if 'U' not in self.fields.keys():
            self.fields['U'] = U
        else:
            self.fields['U'] += U


    def add_vvel(self, V):
        """ Add northward velocity of the initial condition. """ 

        V = V*np.ones((self.nx, self.ny, self.nz))  # Broadcast to proper size

        if 'V' not in self.fields.keys():
            self.fields['V'] = V
        else:
            self.fields['V'] += V



    def add_theta(self, T):
        """ Add potential temperature profile of the initial condition. """ 

        T = T*np.ones((self.nx, self.ny, self.nz))  # Broadcast to proper size

        if 'T' not in self.fields.keys():
            self.fields['T'] = T
        else:
            self.fields['T'] += T


    def add_salt(self, S):
        """ Add salinity profile of the initial condition. """

        S = S*np.ones((self.nx, self.ny, self.nz))  # Broadcast to proper size

        if 'S' not in self.fields.keys():
            self.fields['S'] = S
        else:
            self.fields['S'] += S


    def set_constant_stratification(self, N=1e-3, g=9.81, alpha=2e-4, T0=0.0):
        """ Set initial condition to a constant temperature stratification. """

        Tz = N**2.0 / (g*alpha) 
        T = Tz*self.model.z
        T += T0 - T.min()

        self.set_theta(
              np.ones((self.nx, self.ny))[:, :, np.newaxis]
            * T[np.newaxis, np.newaxis, :]
        )


    def add_barotropic_flow(self, U=0.0, V=0.0):
        """ Add a barotropic flow to the initial condition. """
        
        if U != 0.0:
            self.add_uvel(U*np.ones((self.nx, self.ny, self.nz)))

        if V != 0.0:
            self.add_vvel(U*np.ones((self.nx, self.ny, self.nz)))
