import os, shutil
import numpy as np
import time
import f90nml

from . import gcmutils




class RectangularModel:
    def __init__(self, nx, ny, nz, Lx, Ly, Lz, gcmdir=None):
        """ Instantiate a rectangular MITgcm model on a Cartesian grid. 
        
        Args:
            nx, ny, nz: Number of grid points of the rectangular domain in 
                x (east), y (north), and z (vertical).

            Lx, Ly, Lz: Extent of the model domain in x, y, and z.

            gcmdir: Path to local version of MITgcm. If unset, gcmdir defaults
                to the path '../MITgcm' if it exists.
        """

        self.Lx, self.Ly, self.Lz = Lx, Ly, Lz
        self.nx, self.ny, self.nz = nx, ny, nz

        self.topofilename = 'topo.bin'

        if gcmdir is None and os.path.exists('../MITgcm'):
            self.gcmdir = os.path.abspath('../MITgcm')
        else:
            self.gcmdir = gcmdir

        self.initgrid()

        # Initialize namelist patch
        self.namelists = {
            'data': {
                'parm01': {},       # Continuous equation
                'parm03': {},       # Time-stepping
                'parm04': {},       # Gridding
                'parm05': {},       # Bathymetry initial condition files
            },
        }

        # Some default names for bathymetry files
        self.namelists['data']['parm04']['usingcartesiangrid'] = True

        # Convenience abbreviations
        self.phys = self.namelists['data']['parm01']


    def initgrid(self):
        """ Initialize the grid for a two-dimensional model in x, z """

        self.dx = gcmutils.truncate(self.Lx/self.nx) * np.ones((self.nx,))
        self.dy = gcmutils.truncate(self.Ly/self.ny) * np.ones((self.ny,))
        self.dz = gcmutils.truncate(self.Lz/self.nz) * np.ones((self.nz,))

        self.x = self.dx.cumsum()
        self.y = self.dy.cumsum()

        zf = np.concatenate(([0.0], -gcmutils.truncate(self.dz.cumsum())))

        self.z = 0.5*(zf[0:-1] + zf[1:])

        self.Y, self.X = np.meshgrid(self.y, self.x)


    def set_thetaref(self, Tref=None):
        """ Set the reference temperature profile. """

        if Tref is not None:
            self.Tref = Tref
        elif hasattr(self, 'ic') and 'T' in self.ic.fields.keys():
            self.Tref = self.ic.fields['T'].mean(axis=(0, 1))
        elif hasattr(self, 'obcs'):
            for obc in self.obcs.keys():
                if 'T' in self.obcs[obc].fields.keys():
                    if not hasattr(self, 'Tref'):
                        self.Tref = self.obcs[obc].fields['T'].mean(axis=0)
                    else:
                        self.Tref = np.concatenate(
                            (self.Tref, self.obcs[obc].fields['T'].mean(axis=0)))

        if not hasattr(self, 'Tref'):
            self.Tref = np.zeros((self.nz,))

        self.Tref = gcmutils.truncate(self.Tref, digits=5)


    def set_saltref(self, Sref=None):
        """ Set the reference salinity profile. """

        if Sref is not None:
            self.Sref = Sref
        elif hasattr(self, 'ic') and 'S' in self.ic.fields.keys():
            self.Sref = self.ic.fields['S'].mean(axis=(0, 1))
        elif hasattr(self, 'obcs'):
            for obc in self.obcs.keys():
                if 'S' in self.obcs[obc].fields.keys():
                    if not hasattr(self, 'Sref'):
                        self.Sref = self.obcs[obc].fields['S'].mean(axis=0)
                    else:
                        self.Sref = np.concatenate(
                            (self.Sref, self.obcs[obc].fields['S'].mean(axis=0)))

        if not hasattr(self, 'Sref'):
            self.Sref = 35.0*np.ones((self.nz,))

        self.Sref = gcmutils.truncate(self.Sref, digits=5)

    def updatenamelists(self):
        """ Push essential fields into the model's namelists attribute. """

        # Clean house before updating the namelists
        if not hasattr(self, 'Tref'): self.set_thetaref()
        if not hasattr(self, 'Sref'): self.set_saltref()

        # Update namepatch with required items...
        self.namelists['data']['parm01']['Tref'] = list(self.Tref)
        self.namelists['data']['parm01']['Sref'] = list(self.Sref)

        self.namelists['data']['parm04']['delx'] = list(self.dx)
        self.namelists['data']['parm04']['dely'] = list(self.dy)
        self.namelists['data']['parm04']['delz'] = list(self.dz)

        # Names of the topography file
        if hasattr(self, 'topo'):
            self.namelists['data']['parm05']['bathyfile'] = self.topofilename

        # Names of the binary files that specify initial conditions
        if hasattr(self, 'ic'):
            for var in self.ic.fields.keys():
                self.namelists['data']['parm05'][
                    self.ic.namelistnames[var]] = (self.ic.filenames[var])

        if hasattr(self, 'obcs'):
            self.namelists['data.pkg'] = { 'packages': { 'useobcs': True } }
            self.namelists['data.obcs'] = { 'obcs_parm01' : {} }

            for obc in self.obcs.keys():

                if obc is 'east' or obc is 'west':      idx = 'I'
                elif obc is 'south' or obc is 'north':  idx = 'J'

                self.namelists['data.obcs']['obcs_parm01'][
                    'ob_'+idx+obc] = list(getattr(self.obcs[obc], idx))

                for var in self.obcs[obc].fields.keys():
                    self.namelists['data.obcs']['obcs_parm01'][
                        self.obcs[obc].namelistnames[var]] = (
                        self.obcs[obc].filenames[var] )    


    def saveinput(self):
        """ Save available grid files, initial conditions, and boundary 
        conditions to disk. """


        # Topography
        if hasattr(self, 'topo'):
            topovar = { self.topofilename : self.topo }
            gcmutils.savegcminput(topovar, self.setupdirs['inputdir'])

        # Initial condition and boundary conditions
        if hasattr(self, 'ic'):
            self.ic.save(self.setupdirs['inputdir'])

        if hasattr(self, 'obcs'):
            for obc in self.obcs.keys():
                self.obcs[obc].save(self.setupdirs['inputdir'])




    def gensetup(self, namepatch=None, nprun=1, templatedir=None, 
        workdir=None, cleansetup=False):
        """ Generate the MITgcm setup for the process model using a template
        setup.

        Args:
            namepatch (dict): A namelist patch dictionary with the structure
                namepatch[namelistfilename][namelist][variable] = value.
                If specified, it will be merged with the existing default.

            nprun (int): Number of processors for the run.

            templatedir (str): Path to a pymitgcm template directory with
                /namelist and /code subdirectories. If 'None', defaults to 
                the attribute self.templatedir, which must exist.

            workdir (str): Path to the working directory in which to put
                the setup.

            cleansetup (bool): Boolean indicating whether to erase an
                existing setup.
        """

        if self.nx % nprun != 0:
            raise ValueError("Number of grid points must be a mulitple of the "
                "number of run processes.")

        # Initialization: save inputs and patch the model's namelist
        self.initsetupdirs(workdir=workdir, cleansetup=cleansetup)
        self.saveinput()
        self.updatenamelists()
        if namepatch is not None: 
            self.namelists = {**namepatch, **self.namelists}

        # Identify templatedir and daughter code and namelist dirs
        if templatedir is None: templatedir = self.templatedir
        nametempldir = '{}/namelists'.format(templatedir)
        codetempldir = '{}/code'.format(templatedir)

        # Copy code
        for filename in os.listdir(codetempldir):
            shutil.copy('{}/{}'.format(codetempldir, filename),
                '{}/{}'.format(self.setupdirs['codedir'], filename))

        # Modify SIZE.h
        gcmutils.changesizevars(
            {'sNx': int(self.nx/nprun), 
             'Nr' : self.nz,
             'nPx': nprun }, sizedir=self.setupdirs['codedir'])

        #for var, value in sizevars.items():
        #gcmutils.changesizevar(var, value, codedir=self.setupdirs['codedir'])

        # Direct copy namelists not in the model namelist dictionary
        for filename in os.listdir(nametempldir):
            if filename not in self.namelists.keys():
                shutil.copy('{}/{}'.format(nametempldir, filename), 
                    '{}/{}'.format(self.setupdirs['inputdir'], filename)) 

        # Merge remaining template namelists with model namelist
        fullnamelists = {}
        for nmlfile in self.namelists.keys():

            fullnamelists[nmlfile] = f90nml.read('{}/{}'.format(
                nametempldir, nmlfile))

            for nml in self.namelists[nmlfile].keys():
                for var in self.namelists[nmlfile][nml].keys():
                    fullnamelists[nmlfile][nml][var] = (
                        self.namelists[nmlfile][nml][var])

        # Save namelists
        for filename in fullnamelists.keys():
            savename = '{}/{}'.format(self.setupdirs['inputdir'], filename)
            with open(savename, 'w') as namefile:
                fullnamelists[filename].write(savename, force=True)



    def initsetupdirs(self, workdir=None, cleansetup=False):
        """ Initialize the directory structure of a pymitgcm setup. """

        if workdir is None:
            if hasattr(self, 'setupdirs'):
                workdir = self.setupdirs['workdir']
            else:
                workdir = os.path.abspath('.')

        setupdirs = { 
            'workdir'  : workdir,
            'builddir' : '{}/build'.format(workdir),
            'codedir'  : '{}/code'.format(workdir),
            'inputdir' : '{}/input'.format(workdir),
            'rundir'   : '{}/run'.format(workdir),
        }

        if cleansetup:
            # Remove and remake paths.
            for dir, path in setupdirs.items():
                if dir is not 'workdir' and os.path.exists(path):
                    shutil.rmtree(path)

        # Make directories if they don't exist
        for path in setupdirs.values():
            if not os.path.exists(path):
                os.makedirs(path)

        self.setupdirs = setupdirs


    def run(self):
        """ Run the model. """

        starttime = time.time()
        msg = gcmutils.rungcm(self.setupdirs['rundir'], 
            self.setupdirs['builddir'], self.setupdirs['inputdir'])
        print('Run time: {:.3f}'.format(time.time() - starttime))

        print(msg.decode('utf-8'))    


    def compile(self, optfilename=None, npmake=1, nprun=1):
        """ Compile the model. """

        if nprun > 1:   mpi = True
        else:           mpi = False

        os.chdir(self.setupdirs['builddir'])
            
        starttime = time.time()
        gcmutils.genmake(self.gcmdir, optfilename=optfilename, mpi=mpi)
        print('Genmake time: {:3f} s'.format(time.time()-starttime))
            
        starttime = time.time()
        gcmutils.makedepend()
        print('Make depend time: {:3f} s'.format(time.time()-starttime))

        starttime = time.time()
        gcmutils.make(npmake=npmake)
        print('Make time: {:3f} s'.format(time.time()-starttime))
