import os, shutil, time, glob, logging
import numpy as np
import f90nml

from . import gcmutils
from .initialconditions      import InitialCondition
from .openboundaryconditions import OpenBoundaryCondition



# Globals
topofilename = 'topo.bin'
gridfilenames = { 'dx': 'dx.bin', 'dy': 'dy.bin', 'dz': 'dz.bin' }
validnamelists = ['data', 'data.pkg', 'data.obcs', 'data.mnc', 'data.exf', 
    'eedata', 'data.exch2']

logger = logging.getLogger(name='models')
logger.setLevel(0)



# ----------------------------------------------------------------------------- 
# Notes
# -----
#
# MITgcm model file structure:
#
#   base/
#       build/
#       code/
#           SIZE.h
#           ... etc
#       input/
#           namelists [data, data.pkg, ...]
#           binaries [ics, obcs, ...]
#       run/
#           executable
#
# Essential properties readable from template
#   - Resolution (from SIZE.h)
#   - Namelist properties (from input namelists)
#   - Model extent (but only by parsing namelists or reading namelist binaries)
#       Can implement delxfile and delx for Cartesian geometries.
#
# Sources of fragility:
#   - Changing model size and partially-copying binaries 
#       (example: change initial condition but not obcs.)
#       Solution is to only allow changes in initial condition or obcs if
#       grid is known?
#
#
# ----------------------------------------------------------------------------- 
class Setup:
    def __init__(self, *args, init=False, clean=False, build=None, code=None, 
        input=None, run=None):
        """ Instantiate a setup directory tree. Keyword arguments that
        correspond to elements of validsetupdirs defined below will be 
        assigned as attributes to the directory tree object. 

        Args:
            path (str): Path to the setup. If not provided, the current
                working directory is used.

        Keyword args:
            init (bool): Whether or not to initalize the setup.

            clean (bool): Whether or not to remove existing directories
                and files in the setup path.

            build (str): Path to the setup's build directory.

            code (str): Path to the setup's code directory.

            input (str): Path to the setup's input directory.

            run (str): Path to the setup's run directory.
        """
    
        if len(args) > 1:
            raise TypeError("Expected at most 1 argument, got %d" % len(args))

        try:    self.path = os.path.abspath(args[0])
        except: self.path = os.getcwd()

        validsetupdirs = ['build', 'code', 'input', 'run']
        self.execname = 'mitgcmuv'
        self.dirs = []

        # Initialize the directory tree
        for dir in validsetupdirs:

            try: # to use non-default keywords:
                dirpath = os.path.join(self.path, locals()[dir])
            except: # use default names:
                dirpath = os.path.join(self.path, dir) 

            if os.path.exists(dirpath):
                if clean: self.cleandir(dirpath)
                setattr(self, dir, dirpath)
            elif init:
                os.makedirs(dirpath)
                setattr(self, dir, dirpath)
            else:
                setattr(self, dir, None)

        for dir in validsetupdirs:
            if getattr(self, dir) is not None:
                self.dirs.append(dir)

        self.getinputfilenames()
        try: self.getsize()
        except: pass


    def cleandir(self, dirpath):
        """ Remove files and folders from a directory. """

        for file in os.listdir(dirpath):
            filepath = os.path.join(dirpath, file)
            try:
                if os.path.isfile(filepath):  os.remove(filepath)
                elif os.path.isdir(filepath): shutil.rmtree(filepath)
            except:
                pass


    def getsize(self):
        """ Parse SIZE.h and assign nx, ny, nz, and nprun to setup object. """
        self.size = gcmutils.readsizevars(self.code)

        self.nx = self.size['nSx']*self.size['sNx']
        self.ny = self.size['nSy']*self.size['sNy']
        self.nz = self.size['Nr']
        self.nprun = self.size['nPx']*self.size['nPy']

        return self.nx, self.ny, self.nz, self.nprun

    
    def buildparamdict(self):
        """ Build a dictionary containing all the setup parameters. """
 
        paramdict = {}
        for filename in os.listdir(self.input):
            filepath = os.path.join(self.input, filename)
            try: 
                nml = f90nml.read(filepath)
                if len(nml) > 0: paramdict[filename] = nml
            except: 
                pass

        try: paramdict['size'] = gcmutils.readsizevars(self.code)
        except: pass

        return paramdict
       

    def getinputfilenames(self):
        """ Extract key namelist parameters that correspond to filenames. """

        filenameparams = ['delxfile', 'delyfile', 'delrfile', 
            'delzfile', 'bathyfile']
        paramdict = self.buildparamdict()

        self.inputfilenames = {}
        for param in filenameparams:
            filename = gcmutils.siftdict(param, paramdict)[0]
            if filename is not None:
                self.inputfilenames[param] = filename


    def hascorrectinput(self):
        """ Dummy placeholder. Want to return whether or not setup inputs are 
        correct --- which means they must exist, correspond to filenames in 
        namelists, and that input grids and files have a size that corresponds 
        to 'SIZE.h' """

        return False


    def iscompiled(self):
        """ Return whether setup is compiled. At the moment it only tests 
        to see if the setup has an executable. """

        if not self.hasexecutable(): 
            return False
        elif not self.hascorrectinput():
            return False
        else:
            return True


    def getparam(self, param):
        """ Look through input namelists and SIZE.h in search of a parameter
        and return its value. """

        return gcmutils.siftdict(param, self.buildparamdict())[0]


    def setparam(self, param, value, checksize=True):
        """ Look through input namelists and SIZE.h in search of a parameter,
        find its namelist file and namelist, and change its value. """
        
        val, dpath, level = gcmutils.siftdict(param, self.buildparamdict())

        if dpath[0] is 'size' and checksize:
            gcmutils.changesizevars({param: value}, sizepath=self.input)
        else:
            self.editnamelist(dpath[0], dpath[1], param, value)


    def editnamelist(self, namelistfile, namelist, param, value):
        """ Edit a parameter within one of the setup's namelists. """

        patch = {namelist: {param: value}}
        namelistpath = os.path.join(self.input, namelistfile)

        oldnamelistpath = namelistpath + '_temp'
        shutil.copy(namelistpath, oldnamelistpath)
        os.remove(namelistpath)

        f90nml.patch(oldnamelistpath, patch, namelistpath)
        os.remove(oldnamelistpath)
     

    def hasexecutable(self):
        """ Return true if an mitgcmuv executable is found in the Setup's 
        build or run dirs. """
        try:
            if os.path.exists(os.path.join(self.run, self.execname)):   
                hasexecutable = True
            elif os.path.exists(os.path.join(self.build, self.execname)): 
                hasexecutable = True
            else: 
                hasexecutable = False
        except:
            hasexecutable = False

        return hasexecutable


    def runsetup(self, overwrite=False):
        """ Run the setup by invoking its executable. 
        
        Args:
            overwrite (bool): Overwrite existing input/output. If output exists,
                runsetup(overwrite=False) will fail.
        """

        logger.info("Running a setup with properties:\n"
          + "   nx, ny : {}, {}\n".format(self.nx, self.ny)
          + "   nz     : {}\n".format(self.nz)
          + "   nprun  : {}\n".format(self.nprun))

        starttime = time.time()
        msg = gcmutils.rungcm(self.run, inputpath=self.input,
            overwrite=overwrite)

        logger.info("Run time: {:.3f}".format(time.time() - starttime))
        logger.info(msg.decode('utf-8'))


    def compilesetup(self, gcmpath, optfile=None, npmake=1, mnc=True):
        """ Compile a setup.
    
        Args:
            optfile (str): Name of the optfile to use with genmake2

            npmake (int): Number of processors to use during compilation
    
            mnc (bool): Whether or not to compile with NetCDF saving enabled.
        """

        self.getsize()

        if self.nprun > 1: mpi=True
        else:              mpi=False

        if self.build is None:
            buildpath = os.path.join(self.path, 'build')
            os.makedirs(buildpath)
            self.build = buildpath

        if mnc: mncmsg = 'Enabled'
        else:   mncmsg = 'Disabled'

        logger.info("Compiling MITgcm setup in {} with\n".format(self.path)
          + "   nx, ny : {}, {}\n".format(self.nx, self.ny)
          + "   nz     : {}\n".format(self.nz)
          + "   nprun  : {}\n".format(self.nprun)
          + "   npmake : {}\n".format(npmake)
          + "      mnc : {}\n".format(mncmsg)
          + "...\n\n")

        os.chdir(self.build)

        starttime = time.time()
        gcmutils.genmake(gcmpath, optfile=optfile, mpi=mpi, mnc=mnc)
        logger.info("   Genmake completed in     {.3f} s".format(
            time.time()-starttime))
            
        starttime = time.time()
        gcmutils.makedepend()
        logger.info("   Make depend completed in {.3f} s".format(
            time.time()-starttime))

        starttime = time.time()
        gcmutils.make(npmake=npmake)
        logger.info("   Make completed in        {.3f} s".format(
            time.time()-starttime))

        # Copy executable to run directory
        shutil.copy(os.path.join(self.build, self.execname),
            os.path.join(self.run, ''))


    


    
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
# M O D E L C L A S S 
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
class Model:
    def __init__(self, templatepath, gcmpath=None, mnc=True):
        """ Initialize an pymitgcm model. 

        Args:
            templatepath (str): Path to the model's template setup.
        
            gcmpath (str): Path to MITgcm to use to compile the model.
            
            mnc (bool): Whether or not to use NetCDF to save output.
        """

        self.template = Setup(os.path.abspath(templatepath))

        try:              self.gcmpath = os.path.abspath(gcmpath)
        except TypeError: self.gcmpath = os.environ['MITGCMPATH']
        except KeyError:  pass

        # Initialize namelist patch
        self.namelists = {
            'data': {
                'parm01': {},       # Continuous equation
                'parm03': {},       # Time-stepping
                'parm04': {},       # Gridding
                'parm05': {},       # Bathymetry initial condition files
            },
            'data.pkg': {
                'packages': {},     # Packages to load
            }, 
        }

        # Convenience names for common namelists.
        self.eqns  = self.namelists['data']['parm01']
        self.grid  = self.namelists['data']['parm04']
        self.files = self.namelists['data']['parm05']
        self.time  = self.namelists['data']['parm03']
        self.pkgs  = self.namelists['data.pkg']['packages']
        self.size  = {}

        self.mnc = mnc


    def run(self, overwrite=False, remove=None, compile=False, optfile=None, 
        npmake=1):
        """ Run the model. 

        Args:
            overwrite (bool): Whether or not to overwrite input
                in the run directory. Included for safety so that files
                are not erased accidentally.

            remove (list): List of strings and string patterns specifying files
                to remove from run directory before running the model.

            compile (bool): Whether or not to compile the executable prior
                to running.

            optfile (str): Name of the optfile to be used when
                compiling the model. Only used if compile=True.

            npmake (int): Number of processors to use during compilation.
                Only used if compile=True.
        """

        output = 'run.info'
        defaultremove = ['grid*', 'monitor*', 'phiHyd*', output]

        if remove is None:
            remove = defaultremove
        else:
            try:
                remove = remove + defaultremove
            except TypeError: # because 'remove' is not a list
                remove = [remove] + defaultremove

        for pattern in remove:
            for filename in glob.glob(os.path.join(self.setup.run, pattern)):
                try:    os.remove(filename)
                except: pass

        if compile:
            self.compile(optfile=optfile, npmake=npmake)

        logger.setLevel(20)
        logger.error("\nExecuting MITgcm setup:\n"
          + "   Setup path : {}\n".format(self.setup.path)
          + "   nx, ny, nz : {}, {}, {}\n".format(self.nx, self.ny, self.nz)
          + "        nprun : {}\n".format(self.nprun)
          + "   ntimesteps : {}\n".format(self.setup.getparam('ntimesteps'))
          + "       deltaT : {}\n\n".format(self.setup.getparam('deltat'))
          + "Running...")
     
        starttime = time.time()
        msg = gcmutils.rungcm(self.setup.run, inputpath=self.setup.input, 
            overwrite=overwrite, outputname=output)

        logger.error("  Run completed in {:.3f} s.\n".format(
            time.time()-starttime))
        logger.error("MITgcm messages: \n{}".format(msg.decode('utf-8')))


    def compile(self, optfile=None, npmake=1):
        """ Compile the model. 
    
        Args:
            optfile (str): Name of the optfile to use with genmake2

            npmake (int): Number of processors to use during compilation
        """

        # The setup filestructure and setup attribute must be 
        # initialized to compile the model.
        shutil.rmtree(self.setup.build)
        os.makedirs(self.setup.build)
        os.chdir(self.setup.build)

        mpi = True if self.nprun > 1 else False
        withobcs = True if self.pkgs['useobcs'] else False

        logger.setLevel(20)
        logger.error("\nInitiating compilation of MITgcm setup:\n"
          + "   Setup path : {}\n".format(self.setup.path)
          + "  MITgcm path : {}\n".format(self.gcmpath)
          + "   nx, ny, nz : {}, {}, {}\n".format(self.nx, self.ny, self.nz)
          + "        nprun : {}\n".format(self.nprun)
          + "       npmake : {}\n".format(npmake)
          + "          mnc : {}\n".format('Yes' if self.mnc else 'No')
          + "         obcs : {}\n".format('Yes' if withobcs else 'No')
          + "\nCompiling...")
            
        totaltime, starttime = time.time(), time.time()
        gcmutils.genmake(self.gcmpath, optfile=optfile, mpi=mpi, mnc=self.mnc, 
            obcs=withobcs)
        logger.error("           Genmake2 completed in {:7.3f} s...".format(
            time.time()-starttime))
            
        starttime = time.time()
        gcmutils.makedepend()
        logger.error("        Make depend completed in {:7.3f} s...".format(
            time.time()-starttime))

        starttime = time.time()
        gcmutils.make(npmake=npmake)
        logger.error("               Make completed in {:7.3f} s...".format(
            time.time()-starttime))

        # Copy executable to run directory
        shutil.copy(os.path.join(self.setup.build, 'mitgcmuv'),
            os.path.join(self.setup.run, ''))


        logger.error("Build complete. Total build time {:7.3f} s.\n".format(
            time.time()-totaltime))


    def init_ic(self):
        """ Initialize the model's initial condition. """
        self.ic = InitialCondition(self)


    def init_obcs_namelists(self):
        """ Initialize the obcs namelists. """
        if 'data.obcs' not in self.namelists.keys():
            self.namelists['data.obcs'] = { 'obcs_parm01' : {},
                                            'obcs_parm02' : {}  }

        try:
            self.obparams = self.namelists['data.obcs']['obcs_parm01']
        except KeyError:
            self.namelists['data.obcs']['obcs_parm01'] = {}
            self.obparams = self.namelists['data.obcs']['obcs_parm01']

        try:
            self.orlanski = self.namelists['data.obcs']['obcs_parm02']
        except KeyError:
            self.namelists['data.obcs']['obcs_parm02'] = {}
            self.orlanski = self.namelists['data.obcs']['obcs_parm02']


    def init_obcs(self, edges, nt=1, dt=None, copyic=None, orlanski=False, 
        balanced=True):
        """ Initialize open boundary conditions. """

        self.init_obcs_namelists()
        self.obparams['useobcsbalance'] = True if balanced else False

        self.obcs = {}
        for edge in edges:
            self.obcs[edge] = OpenBoundaryCondition(self, edge, nt=nt, dt=dt,
                copyic=copyic, orlanski=orlanski)
        

    def init_obcs_sponge(self, thickness, urelaxbound=0.0, vrelaxbound=0.0,
        urelaxinner=0.0, vrelaxinner=0.0):
        """ Initialize an obcs sponge layer. Sponge layers damp velocities
        normal to openboundaries; thus a 'urelax' parameter controls 
        east-west boundaries and a 'vrelax' parameter control north-south.

        Args:
            thickness (int): Thickness of the sponge layer in grid points.
                Sets the 'spongethickness' obcs parameter.

            urelaxinner (float): relaxation time scale at the innermost sponge
                layer point of an east-west open boundary. Sets the
                'urelaxobcsinner' obcs parameter.

            vrelaxinner (float): relaxation time scale at the innermost sponge
                layer point of a north-south open boundary. Sets the
                'vrelaxobcsinner' obcs parameter.

            urelaxbound (float): relaxation time scale at the outermost sponge
                layer point of an east-west open boundary. Sets the
                'urelaxobcsbound' obcs parameter.

            vrelaxbound (float): relaxation time scale at the outermost sponge
                layer point of a north-south open boundary. Sets the
                'vrelaxobcsbound' obcs parameter.
        """

        self.namelists['data.obcs'] = { 'obcs_parm03': {} }
        self.sponge = self.namelists['data.obcs']['obcs_parm03']

        self.sponge['spongethickness'] = thickness

        self.sponge['urelaxobcsinner'] = urelaxinner
        self.sponge['vrelaxobcsinner'] = vrelaxinner
        self.sponge['urelaxobcsbound'] = urelaxbound
        self.sponge['vrelaxobcsbound'] = vrelaxbound


    def set_thetaref(self, Tref=None):
        """ Set the reference temperature profile. """

        if Tref is not None:
            self.Tref = Tref
        elif hasattr(self, 'ic') and 'T' in self.ic.fields.keys():
            self.Tref = self.ic.fields['T'].mean(axis=(0, 1))
       
        if not hasattr(self, 'Tref'): # set default Tref:
            self.Tref = [0.0]*self.nz

        self.Tref = gcmutils.truncate(self.Tref, digits=4)


    def set_saltref(self, Sref=None):
        """ Set the reference salinity profile. """

        if Sref is not None:
            self.Sref = Sref
        elif hasattr(self, 'ic') and 'S' in self.ic.fields.keys():
            self.Sref = self.ic.fields['S'].mean(axis=(0, 1))
        
        if not hasattr(self, 'Sref'): # set default Sref:
            self.Sref = [35.0]*self.nz

        self.Sref = gcmutils.truncate(self.Sref, digits=5)


    def updatesize(self, nprun=None):
        """ Push essential fields into the model's size attribute. """

        if nprun is not None:
            self.nprun = nprun

        # Horizontal tiling
        try:
            self.tiling = gcmutils.getcompacttiling(
                self.nx, self.ny, self.nprun)
            
            self.size['sNx'] = int(self.nx/self.tiling[0])
            self.size['sNy'] = int(self.ny/self.tiling[1])
            self.size['nSx'] = self.tiling[0]
            self.size['nSy'] = self.tiling[1]
            self.size['nPx'] = self.tiling[0]
            self.size['nPy'] = self.tiling[1]

        except TypeError:    
            raise ValueError(
                "The grid nx, ny = {}, {}".format(self.nx, self.ny)
              + " has no valid tiling for np = {}.".format(nprun))

        # Vertical grid
        self.size['Nr'] = self.nz
                        

    def updatenamelists(self):
        """ Push essential fields into the model's namelists attribute. """

        # Clean house before updating the namelists
        if not hasattr(self, 'Tref'): self.set_thetaref()
        if not hasattr(self, 'Sref'): self.set_saltref()

        # Update namepatch with required items...
        if hasattr(self, 'Tref'): self.eqns['Tref'] = list(self.Tref)
        if hasattr(self, 'Sref'): self.eqns['Sref'] = list(self.Sref)

        self.grid['delxfile'] = gridfilenames['dx']
        self.grid['delyfile'] = gridfilenames['dy']
        self.grid['delz']     = list(self.dz)

        # Names of the topography file
        if hasattr(self, 'topo'):
            self.files['bathyfile'] = topofilename

        # Names of the binary files that specify initial conditions
        if hasattr(self, 'ic'):
            for var in self.ic.fields.keys():
                self.files[self.ic.namelistnames[var]] = self.ic.filenames[var]
                    
        if hasattr(self, 'obcs'):

            self.pkgs['useobcs'] = True
            self.init_obcs_namelists()

            # Ensure property specification of obcs time-dependence        
            obc0 = [obc for obc in self.obcs.values()][0]
            if obc0.nt > 1: 
                self.time['periodicexternalforcing'] = True 
                self.time['externforcingperiod'] = obc0.dt
                self.time['externforcingcycle'] = obc0.dt*obc0.nt
            else:
                self.time['periodicexternalforcing'] = False
           
            # Set obcs params
            for obc in self.obcs.keys():

                if obc is 'east' or obc is 'west':      idx = 'I'
                elif obc is 'south' or obc is 'north':  idx = 'J'

                self.obparams['ob_'+idx+obc] = list(getattr(self.obcs[obc], idx))
                if self.obcs[obc].orlanski:
                    self.obparams['useorlanski'+obc] = True
                    # Default cmax and cveltimescale
                    if 'cmax' not in self.orlanski.keys(): 
                        self.orlanski['cmax'] = 0.45
                    if 'cveltimescale' not in self.orlanski.keys(): 
                        self.orlanski['cveltimescale'] = 1000
                else:
                    self.obparams['useorlanski'+obc] = False

                for var in self.obcs[obc].fields.keys():
                    self.namelists['data.obcs']['obcs_parm01'][
                        self.obcs[obc].namelistnames[var]] = (
                        self.obcs[obc].filenames[var] )    


        if self.mnc:
            self.pkgs['usemnc'] = True
            self.namelists['data.mnc'] = { 'mnc_01': {} }
            self.namelists['data.mnc']['mnc_01'] = {'mnc_use_outdir'  : False}
            self.namelists['data.mnc']['mnc_01'] = {'monitor_mnc'     : False}
            self.namelists['data.mnc']['mnc_01'] = {'pickup_read_mnc' : False}
            self.namelists['data.mnc']['mnc_01'] = {'pickup_write_mnc': False}


    def saveinput(self):
        """ Save available grid files, initial conditions, and boundary 
        conditions to disk. """

        # Grid
        gridvars = { gridfilenames['dx']: self.dx,
                     gridfilenames['dy']: self.dy }
        gcmutils.savegcminput(gridvars, self.setup.input)

        # Topography
        if hasattr(self, 'topo'):
            topovar = { topofilename : self.topo }
            gcmutils.savegcminput(topovar, self.setup.input)

        # Initial condition and boundary conditions
        if hasattr(self, 'ic'):
            self.ic.save(self.setup.input)

        if hasattr(self, 'obcs'):
            for obc in self.obcs.keys():
                self.obcs[obc].save(self.setup.input)


    def gensetup(self, namepatch=None, params=None, templatepath=None, 
        setuppath=None, cleansetup=False):
        """ Generate the MITgcm setup for the process model using a template
        setup.

        Args:
            namepatch (dict): A namelist patch dictionary with the structure
                namepatch[namelistfilename][namelist][variable] = value.
                If specified, it will be merged with the existing default.

            params (dict): A dictionary of parameters with the form 
                {paramname: value}. Each paramname is looked up in the 
                template namelists and, if found, changed to value.
                Each param must have a pre-existing template value for any 
                change to take place.

            setuppath (str): Path to the working directory in which to put
                the setup.

            cleansetup (bool): Boolean indicating whether to erase an
                existing setup.
        """

        if namepatch is not None: 
            self.namelists = {**namepatch, **self.namelists}

        
        self.setup = Setup(setuppath, init=True, clean=cleansetup)
        self.saveinput()
        self.updatenamelists()
        self.updatesize()

        # Copy template
        dontcopylist = ['.swp', 'mitgcmuv']
        for dir in ['code', 'input']:
            for filename in os.listdir(getattr(self.template, dir)):
                if not any(bad in filename for bad in dontcopylist):
                    shutil.copy(
                        os.path.join(getattr(self.template, dir), filename),
                        os.path.join(getattr(self.setup, dir), filename))

        # Change size vars
        gcmutils.changesizevars(self.size, sizepath=self.setup.code)

        # Merge template namelists with model namelist
        fullnamelists = {}
        for nmlfile in self.namelists.keys():

            # Load namelist files from template directory, not newly-created
            # setup
            try:
                fullnamelists[nmlfile] = f90nml.read(
                    os.path.join(self.template.input, nmlfile))
            except FileNotFoundError:
                fullnamelists[nmlfile] = f90nml.Namelist(self.namelists[nmlfile])

            # Paste-in new vars
            for nml in self.namelists[nmlfile].keys():
                try:
                    for var in self.namelists[nmlfile][nml].keys():
                        fullnamelists[nmlfile][nml][var] = (
                            self.namelists[nmlfile][nml][var])
                except KeyError:
                    fullnamelists[nmlfile][nml] = f90nml.Namelist(
                        self.namelists[nmlfile][nml])
                    
        # Homogenize references to 'r' or 'z' in data namelist
        fullnamelists['data'] = gcmutils.correctdatanamelist(
            fullnamelists['data'], correctto='z')

        fullnamelists = self.trimnamelists(fullnamelists)

        # Save namelists
        for filename in fullnamelists.keys():
            savename = os.path.join(self.setup.input, filename)
            with open(savename, 'w') as namefile:
                fullnamelists[filename].write(savename, force=True)

        # Change parameters
        if params is not None:
            for param, value in params.items():
                self.setup.setparam(param, value, checksize=False)


    def trimnamelists(self, namelists):
        """ Trim spurious parameters from namelists. """    

        # Shortcuts
        eqns = namelists['data']['parm01']


        # Compatabilities in the continuous equation params.
        try:
            if eqns['rigidlid']: # No free surface evolution
                eqns['implicitfreesurface'] = False
                eqns['exactconserv'] = False
        except: 
            pass

        # This code uses only delxfile, delyfile, and delz
        for gridvar in ['delx', 'dely', 'delr']:
            try:    namelists['data']['parm04'].pop(gridvar)
            except: pass
        
        # Initial condition params
        icnames = { 
            'U': 'uvelinitfile', 
            'V': 'vvelinitfile', 
            'W': 'wvelinitfile', 
            'T': 'hydrogthetafile', 
            'S': 'hydrogsaltfile', 
        } 

        if hasattr(self, 'ic'):
            for var in icnames.keys():
                if var not in self.ic.fields.keys():
                    try:    namelists['data']['parm05'].pop(icnames[var])
                    except: pass

        # OBCS params
        edges = ['north', 'south', 'east', 'west']
        vars = ['U', 'V', 'W', 'T', 'S']

        if hasattr(self, 'obcs') and 'data.obcs' in namelists.keys():
            obparams = namelists['data.obcs']['obcs_parm01']
            for edge in edges:

                if   edge is 'north' or edge is 'south': idx = 'J'
                elif edge is 'west'  or edge is 'east' : idx = 'I'

                if edge not in self.obcs.keys():
                    param = 'ob_' + idx + edge
                    try:    obparams.pop(param)
                    except: pass

                for v in vars:
                    try: 
                        _ = self.obcs[edge].fields[v].shape
                    except:
                        param = 'ob' + edge[0] + v + 'file'
                        try:    obparams.pop(param)
                        except: pass

        return namelists


class RectangularModel(Model):
    def __init__(self, nx, ny, nz, Lx, Ly, Lz, templatepath, 
        nprun=1, gcmpath=None, mnc=True):
        """ Instantiate a rectangular MITgcm model on a Cartesian grid. 
        
        Args:
            nx, ny, nz: Number of grid points of the rectangular domain in 
                x (east), y (north), and z (vertical).

            Lx, Ly, Lz: Extent of the model domain in x, y, and z.

            templatepath (str): Path to the model's template setup.

            gcmpath (str): Path to local version of MITgcm. If unset, gcmpath 
                defaults to the path '../MITgcm' if it exists.
        """

        Model.__init__(self, templatepath, gcmpath=gcmpath, mnc=mnc)

        self.Lx, self.Ly, self.Lz = Lx, Ly, Lz
        self.nx, self.ny, self.nz = nx, ny, nz

        self.updatesize(nprun=nprun)
            
        self.grid['usingcartesiangrid'] = True
        self.initgrid()


    def initgrid(self):
        """ Initialize the grid for a model in Cartesian x,y,z-space. """

        self.dx = gcmutils.truncate(self.Lx/self.nx) * np.ones((self.nx,))
        self.dy = gcmutils.truncate(self.Ly/self.ny) * np.ones((self.ny,))
        self.dz = gcmutils.truncate(self.Lz/self.nz) * np.ones((self.nz,))

        # Horizontal grid variables
        self.xg = np.concatenate(([0.0], self.dx.cumsum()))
        self.yg = np.concatenate(([0.0], self.dy.cumsum()))
            
        self.x = 0.5*(self.xg[1:]+self.xg[0:-1])
        self.y = 0.5*(self.yg[1:]+self.yg[0:-1])

        # For reference
        self.xc = self.x
        self.yc = self.y
        self.xu = self.xg
        self.yu = self.yc
        self.xv = self.xc
        self.yv = self.yg

        # Vertical grid
        self.zf = np.concatenate(([0.0], -self.dz.cumsum()))
        self.z = 0.5*(self.zf[0:-1] + self.zf[1:])
        self.zc = self.z

        self.Y, self.X, self.Z = np.meshgrid(self.y, self.x, self.z)


    def stretchgridedge(self, edge, width=1, factor=2, 
        stretching='linear'): 
        """ Coarsen the grid at boundaries and refine in center. """
