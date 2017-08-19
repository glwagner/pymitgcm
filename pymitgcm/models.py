import os, shutil
import os, time
import numpy as np
import f90nml

from . import gcmutils
from .initialconditions      import InitialCondition
from .openboundaryconditions import OpenBoundaryCondition



# Globals
topofilename = 'topo.bin'
validnamelists = ['data', 'data.pkg', 'data.obcs', 'data.mnc', 'data.exf', 
    'eedata', 'data.exch2']



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
        correspond to elements of ._validdirs defined below will be 
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
        elif len(args) is 0:
            self.path = os.getcwd()
        elif len(args) is 1:
            self.path = os.path.abspath(args[0])

        self.execname = 'mitgcmuv'
        self._validdirs = ['build', 'code', 'input', 'run']
        self.dirs = []

        # Initialize the directory tree
        for dir in self._validdirs:

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

        for dir in self._validdirs:
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

        filenameparams = ['delxfile', 'delyfile', 'delrfile']
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


    def setparam(self, param, value):
        """ Look through input namelists and SIZE.h in search of a parameter,
        find its namelist file and namelist, and change its value. """
        
        val, dpath, level = gcmutils.siftdict(param, self.buildparamdict())

        if dpath[0] is 'size':
            gcmutils.changesizevars({param: value}, sizepath=self.input)
        else:
            self.editnamelist(dpath[0], dpath[1], param, value)


    def editnamelist(self, namelistfile, namelist, param, value):
        """ Edit a parameter within one of the setup's namelists. """

        patch = {namelist: {param: value}}
        namelistpath = os.path.join(self.input, namelistfile)
        oldnamelistpath = namelistpath + '_old'

        shutil.copy(namelistpath, oldnamelistpath)
        os.remove(namelistpath)
        f90nml.patch(oldnamelistpath, patch, namelistpath)
     

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


    def runsetup(self, overwrite=False, params=None):
        """ Run the setup by invoking its executable. 
        
        Args:
            overwrite (bool): Overwrite existing input/output. If output exists,
                runsetup(overwrite=False) will fail.

            params (dict): Dictionary of {param: value} pairs to be changed
                before running the setup.
        """

        if params is not None:
            pass

        starttime = time.time()
        msg = gcmutils.rungcm(self.run, inputpath=self.input,
            overwrite=overwrite)

        print('Run time: {:.3f}'.format(time.time() - starttime))
        print(msg.decode('utf-8'))    


    def compilesetup(self, gcmpath, optfilename=None, npmake=1):
        """ Compile a setup.
    
        Args:
            optfilename (str): Name of the optfile to use with genmake2
            npmake (int): Number of processors to use during compilation
            nprun (int): Number of processors the model will be run with. Used
                to determine whether or not to compile the model with
                mpi enabled.
            mpi (bool): An flag alternative to specifying nprun that directs
                genmake2 to compile with mpi enabled.
        """

        self.getsize()

        if self.nprun > 1: mpi=True
        else:              mpi=False

        if self.build is None:
            buildpath = os.path.join(self.path, 'build')
            os.makedirs(buildpath)
            self.build = buildpath

        os.chdir(self.build)

        starttime = time.time()
        gcmutils.genmake(gcmpath, optfilename=optfilename, mpi=mpi)
        print('Genmake time: {:3f} s'.format(time.time()-starttime))
            
        starttime = time.time()
        gcmutils.makedepend()
        print('Make depend time: {:3f} s'.format(time.time()-starttime))

        starttime = time.time()
        gcmutils.make(npmake=npmake)
        print('Make time: {:3f} s'.format(time.time()-starttime))

        # Copy executable to run directory
        shutil.copy(os.path.join(self.build, self.execname),
            os.path.join(self.run, ''))


    


    
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
# M O D E L C L A S S 
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
class Model:
    def __init__(self, templatepath=None, gcmpath=None):
        """ Initialize an pymitgcm model. """

        try:              self.gcmpath = os.path.abspath(gcmpath)
        except TypeError: self.gcmpath = os.environ['MITGCMPATH']
        except KeyError:  pass

        try:              templatepath = os.path.abspath(templatepath)
        except TypeError: templatepath = os.getcwd()

        self.template = Setup(templatepath)

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


    def run(self, overwrite=False, optfilename=None, npmake=1):
        """ Run the model. """

        if not self.setup.iscompiled(): 
            self.compile(optfilename=optfilename, npmake=npmake)

        starttime = time.time()
        msg = gcmutils.rungcm(self.setup.run, inputpath=self.setup.input, 
            overwrite=overwrite)

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

        # The setup filestructure and setup attribute must be 
        # initialized to compile the model.
        os.chdir(self.setup.build)

        if self.nprun > 1: mpi=True
        else:              mpi=False
            
        starttime = time.time()
        gcmutils.genmake(self.gcmpath, optfilename=optfilename, mpi=mpi)
        print('Genmake time: {:3f} s'.format(time.time()-starttime))
            
        starttime = time.time()
        gcmutils.makedepend()
        print('Make depend time: {:3f} s'.format(time.time()-starttime))

        starttime = time.time()
        gcmutils.make(npmake=npmake)
        print('Make time: {:3f} s'.format(time.time()-starttime))

        # Copy executable to run directory
        shutil.copy(os.path.join(self.setup.build, 'mitgcmuv'),
            os.path.join(self.setup.run, ''))


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

        # Assume we are using MNC and set default vars --- for now.
        self.namelists['data.mnc'] = { 'mnc_01': {} }
        self.namelists['data.mnc']['mnc_01'] = {'mnc_use_outdir'  : False}
        self.namelists['data.mnc']['mnc_01'] = {'monitor_mnc'     : False}
        self.namelists['data.mnc']['mnc_01'] = {'pickup_read_mnc' : False}
        self.namelists['data.mnc']['mnc_01'] = {'pickup_write_mnc': False}


    def saveinput(self):
        """ Save available grid files, initial conditions, and boundary 
        conditions to disk. """

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


    def gensetup(self, namepatch=None, nprun=1, templatepath=None, 
        setuppath=None, cleansetup=False):
        """ Generate the MITgcm setup for the process model using a template
        setup.

        Args:
            namepatch (dict): A namelist patch dictionary with the structure
                namepatch[namelistfilename][namelist][variable] = value.
                If specified, it will be merged with the existing default.

            nprun (int): Number of processors for the run.

            templatepath (str): Path to a pymitgcm template directory with
                /namelist and /code subdirectories. If 'None', defaults to 
                the attribute self.templatepath, which must exist.

            setuppath (str): Path to the working directory in which to put
                the setup.

            cleansetup (bool): Boolean indicating whether to erase an
                existing setup.
        """

        if namepatch is not None: 
            self.namelists = {**namepatch, **self.namelists}

        if templatepath is not None:    
            self.template = Setup(templatepath)
        elif not hasattr(self, 'template'):
            raise RunTimeError("No template has been specified for the model.")
        
        if self.nx % nprun != 0:
            raise ValueError("Number of grid points must be a mulitple of the "
                "number of run processes.")
        else:
            self.nprun = nprun

        self.setup = Setup(setuppath, init=True, clean=cleansetup)
        self.saveinput()
        self.updatenamelists()
        self.updatesizevars(nprun=nprun)
            
        # Copy template
        for dir in self.template.dirs:
            for filename in os.listdir(getattr(self.template, dir)):
                shutil.copy(
                    os.path.join(getattr(self.template, dir), filename),
                    os.path.join(getattr(self.setup, dir), filename))

        gcmutils.changesizevars(self.sizevars, sizepath=self.setup.code)

        # Merge template namelists with model namelist
        fullnamelists = {}
        for nmlfile in self.namelists.keys():

            # Load namelist files from template directory, not newly-created
            # setup
            fullnamelists[nmlfile] = f90nml.read(
                os.path.join(template['input'], nmlfile))

            for nml in self.namelists[nmlfile].keys():
                for var in self.namelists[nmlfile][nml].keys():
                    fullnamelists[nmlfile][nml][var] = (
                        self.namelists[nmlfile][nml][var])

        # Save namelists
        for filename in fullnamelists.keys():
            savename = os.path.join(self.setup.input, filename)
            with open(savename, 'w') as namefile:
                fullnamelists[filename].write(savename, force=True)






class RectangularModel(Model):
    def __init__(self, nx, ny, nz, Lx, Ly, Lz, gcmpath=None, 
        templatepath=None):
        """ Instantiate a rectangular MITgcm model on a Cartesian grid. 
        
        Args:
            nx, ny, nz: Number of grid points of the rectangular domain in 
                x (east), y (north), and z (vertical).

            Lx, Ly, Lz: Extent of the model domain in x, y, and z.

            gcmpath: Path to local version of MITgcm. If unset, gcmpath defaults
                to the path '../MITgcm' if it exists.
        """

        Model.__init__(self, templatepath=templatepath, gcmpath=gcmpath)

        self.Lx, self.Ly, self.Lz = Lx, Ly, Lz
        self.nx, self.ny, self.nz = nx, ny, nz
            
        self.grid['usingcartesiangrid'] = True
        self.initgrid()


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


