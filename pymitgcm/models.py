import os, shutil
import os, time
import numpy as np
import f90nml

from . import gcmutils
from .initialconditions      import InitialCondition
from .openboundaryconditions import OpenBoundaryCondition



# Globals
topofilename = 'topo.bin'
validtemplatedirs = ['code', 'input', 'run']


class Model:
    def __init__(templatedir=None, gcmdir=None):

        try:    self.gcmdir = os.path.abspath(gcmdir)
        except: self.gcmdir = os.path.abspath('../MITgcm')

        try: self.templatedir = os.path.abspath(templatedir)
        except: pass

        # Initialize namelist patch
        self.namelists = {
            'data': {
                'parm01': {},       # Continuous equation
                'parm03': {},       # Time-stepping
                'parm04': {},       # Gridding
                'parm05': {},       # Bathymetry initial condition files
            },
        }

        # Convenience names for common namelists.
        self.eqns  = self.namelists['data']['parm01']
        self.grid  = self.namelists['data']['parm04']
        self.files = self.namelists['data']['parm05']
        self.timestepper = self.namelists['data']['parm03']



class RectangularModel:
    def __init__(self, nx, ny, nz, Lx, Ly, Lz, gcmdir=None, 
        templatedir=None):
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

        try:    self.gcmdir = os.path.abspath(gcmdir)
        except: self.gcmdir = os.path.abspath('../MITgcm')

        try: self.templatedir = os.path.abspath(templatedir)
        except: pass
            
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

        # Convenience names for common namelists.
        self.eqns  = self.namelists['data']['parm01']
        self.timestepper = self.namelists['data']['parm03']
        self.grid  = self.namelists['data']['parm04']
        self.files = self.namelists['data']['parm05']

        self.grid['usingcartesiangrid'] = True


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


    def init_ic(self):
        """ Initialize the model's initial condition. """
        self.ic = InitialCondition(self)


    def init_obcs(self, edges):
        """ Initialize open boundary conditions. """
        self.obcs = {}
        for edge in edges:
            self.obcs[edge] = OpenBoundaryCondition(self, edge)


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

        self.Tref = gcmutils.truncate(self.Tref, digits=4)


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


    def updatesizevars(self, nprun=None):
        """ Push essential fields into the model's sizevars attribute. """

        # Initialize
        if nprun is None: nprun = 1
        if not hasattr(self, 'sizevars'): self.sizevars = {}

        # Vertical grid
        self.sizevars['Nr'] = self.nz

        # Parallelism
        if self.nx % nprun == 0.0:
            self.sizevars['sNx'] = int(self.nx/nprun)
            self.sizevars['nSx'] = int(nprun)
            self.sizevars['nPx'] = int(nprun)


    def updatenamelists(self):
        """ Push essential fields into the model's namelists attribute. """

        # Clean house before updating the namelists
        if not hasattr(self, 'Tref'): self.set_thetaref()
        if not hasattr(self, 'Sref'): self.set_saltref()

        # Update namepatch with required items...
        if hasattr(self, 'Tref'): self.eqns['Tref'] = list(self.Tref)
        if hasattr(self, 'Sref'): self.eqns['Sref'] = list(self.Sref)

        # TODO: detect properties of template setup and choose setup accordingly.
        #self.grid['delx'] = list(self.dx)
        #self.grid['dely'] = list(self.dy)
        #self.grid['delz'] = list(self.dz)

        # Names of the topography file
        if hasattr(self, 'topo'):
            self.files['bathyfile'] = topofilename

        # Names of the binary files that specify initial conditions
        if hasattr(self, 'ic'):
            for var in self.ic.fields.keys():
                self.files[self.ic.namelistnames[var]] = self.ic.filenames[var]
                    
        if hasattr(self, 'obcs'):
            self.namelists['data.pkg'] = { 'packages': { 'useobcs': True } }
            self.namelists['data.obcs'] = { 'obcs_parm01' : {} }

            for obc in self.obcs.keys():

                if obc is 'east' or obc is 'west':      idx = 'I'
                elif obc is 'south' or obc is 'north':  idx = 'J'

                #self.namelists['data.obcs']['obcs_parm01'][
                #    'ob_'+idx+obc] = list(getattr(self.obcs[obc], idx))

                for var in self.obcs[obc].fields.keys():
                    self.namelists['data.obcs']['obcs_parm01'][
                        self.obcs[obc].namelistnames[var]] = (
                        self.obcs[obc].filenames[var] )    


    def saveinput(self):
        """ Save available grid files, initial conditions, and boundary 
        conditions to disk. """

        # Topography
        if hasattr(self, 'topo'):
            topovar = { topofilename : self.topo }
            gcmutils.savegcminput(topovar, self.setupdirs['input'])

        # Initial condition and boundary conditions
        if hasattr(self, 'ic'):
            self.ic.save(self.setupdirs['input'])

        if hasattr(self, 'obcs'):
            for obc in self.obcs.keys():
                self.obcs[obc].save(self.setupdirs['input'])


    def gensetup(self, namepatch=None, nprun=1, templatedir=None, 
        setupdir=None, cleansetup=False):
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

            setupdir (str): Path to the working directory in which to put
                the setup.

            cleansetup (bool): Boolean indicating whether to erase an
                existing setup.
        """

        if namepatch is not None: 
            self.namelists = {**namepatch, **self.namelists}

        if templatedir is not None: 
            self.templatedir = os.path.abspath(templatedir)
        elif not hasattr(self, 'templatedir'):
            self.templatedir = os.getcwd()

        if self.nx % nprun != 0:
            raise ValueError("Number of grid points must be a mulitple of the "
                "number of run processes.")
        else:
            self.nprun = nprun

        # Initialization: save inputs and patch the model's namelist
        template = {}
        for dir in validtemplatedirs:
            daughterdir = os.path.abspath('{}/{}'.format(self.templatedir, dir))
            if os.path.exists(daughterdir):
                template[dir] = daughterdir

        if len(template.keys()) is 0:
            raise RuntimeError("Template is either empty or does not exist.")

        self.initsetupdirs(setupdir=setupdir, cleansetup=cleansetup)
        self.saveinput()
        self.updatenamelists()
        self.updatesizevars(nprun=nprun)
            
        # Copy template
        for dir in template.keys():
            for filename in os.listdir(template[dir]):
                shutil.copy('{}/{}'.format(template[dir], filename),
                    '{}/{}'.format(self.setupdirs[dir], filename))

        # Modify SIZE.h
        gcmutils.changesizevars(self.sizevars, sizedir=self.setupdirs['code'])

        # Merge template namelists with model namelist
        fullnamelists = {}
        for nmlfile in self.namelists.keys():

            # Load namelist files from template directory, not newly-created
            # setup
            fullnamelists[nmlfile] = f90nml.read('{}/{}'.format(
                template['input'], nmlfile))

            for nml in self.namelists[nmlfile].keys():
                for var in self.namelists[nmlfile][nml].keys():
                    fullnamelists[nmlfile][nml][var] = (
                        self.namelists[nmlfile][nml][var])

        # Save namelists
        for filename in fullnamelists.keys():
            savename = '{}/{}'.format(self.setupdirs['input'], filename)
            with open(savename, 'w') as namefile:
                fullnamelists[filename].write(savename, force=True)



    def initsetupdirs(self, setupdir=None, cleansetup=False):
        """ Initialize the directory structure of a pymitgcm setup. """

        if setupdir is None: setupdir = os.getcwd()

        setupdirs = { 
            'build' : os.path.join(setupdir, 'build'),
            'code'  : os.path.join(setupdir, 'code'),
            'input' : os.path.join(setupdir, 'input'),
            'run'   : os.path.join(setupdir, 'run'),
        }

        if cleansetup:
            # Remove and remake paths.
            for dir, path in setupdirs.items():
                if dir is not 'base' and os.path.exists(path):
                    shutil.rmtree(path)

        # Make directories if they don't exist
        for path in setupdirs.values():
            if not os.path.exists(path):
                os.makedirs(path)

        self.setupdirs = setupdirs


    def run(self, clean=False):
        """ Run the model. """

        starttime = time.time()

        msg = gcmutils.rungcm(
                       self.setupdirs['run'], 
            inputdir = self.setupdirs['input'],
            clean=clean)

        print('Run time: {:.3f}'.format(time.time() - starttime))
        print(msg.decode('utf-8'))    


    def compile(self, optfilename=None, npmake=1):
        """ Compile the model. 
    
        Args:
            optfilename (str): Name of the optfile to use with genmake2
            npmake (int): Number of processors to use during compilation
            nprun (int): Number of processors the model will be run with. Used
                to determine whether or not to compile the model with
                mpi enabled.
            mpi (bool): An flag alternative to specifying nprun that directs
                genmake2 to compile with mpi enabled.
        """

        # The setup filestructure and setupdirs attribute must be 
        # initialized to compile the model.
        os.chdir(self.setupdirs['build'])

        if self.nprun > 1: mpi=True
        else:              mpi=False
            
        starttime = time.time()
        gcmutils.genmake(self.gcmdir, optfilename=optfilename, mpi=mpi)
        print('Genmake time: {:3f} s'.format(time.time()-starttime))
            
        starttime = time.time()
        gcmutils.makedepend()
        print('Make depend time: {:3f} s'.format(time.time()-starttime))

        starttime = time.time()
        gcmutils.make(npmake=npmake)
        print('Make time: {:3f} s'.format(time.time()-starttime))

        # Copy executable to run directory
        shutil.copy(os.path.join(self.setupdirs['build'], 'mitgcmuv'),
            os.path.join(self.setupdirs['run'], ''))
