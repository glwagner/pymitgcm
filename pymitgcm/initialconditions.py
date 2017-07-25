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
            

    def set_theta(T):
        """ Set potential temperature profile of the initial condition. """ 
        self.fields['T'] = T


    def set_salt(S):
        """ Set the salinity profile of the initial condition. """ 
        self.fields['S'] = S


    def add_uvel(U):
        """ Add an eastward velocity to the initial condition. """ 

        if 'U' not in self.fields.keys():
            self.fields['U'] = U
        else:
            self.fields['U'] += U


    def add_vvel(V):
        """ Add northward velocity of the initial condition. """ 

        if 'V' not in self.fields.keys():
            self.fields['V'] = V
        else:
            self.fields['V'] += V



    def add_theta(T):
        """ Add potential temperature profile of the initial condition. """ 

        if 'T' not in self.fields.keys():
            self.fields['T'] = T
        else:
            self.fields['T'] += T


    def add_salt(S):
        """ Add salinity profile of the initial condition. """

        if 'S' not in self.fields.keys():
            self.fields['S'] = S
        else:
            self.fields['S'] += S


    def set_constant_stratification(N=1e-3, g=9.81, alpha=2e-3, T0=2.0):
        """ Set initial condition to a constant temperature stratification. """

        Tz = N**2.0 / (g*alpha) 
        Tref = Tz*model.z - (Tz*model.z).mean()

        self.set_theta(
              np.ones((self.nx, self.ny))[:, :, np.newaxis]
            * Tref[np.newaxis, np.newaxis, :]
        )

        model.N = N
        model.g = g
        model.alpha = alpha


    def add_barotropic_flow(U=0.0, V=0.0):
        """ Add a barotropic flow to the initial condition. """
        
        if U != 0.0:
            self.add_uvel(U*np.ones((self.nx, self.ny, self.nz)))

        if V != 0.0:
            self.add_vvel(U*np.ones((self.nx, self.ny, self.nz)))
