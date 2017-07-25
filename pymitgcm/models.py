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

        self.topofile = 'topo.bin'

        if gcmdir is None and os.path.exists('../MITgcm'):
            self.gcmdir = os.path.abspath('../MITgcm')
        else:
            self.gcmdir = gcmdir

        self.initgrid()

        # Initialize namelist patch
        self.namepatch = {
            'data': {
                'parm01': {},       # Continuous equation
                'parm04': {},       # Gridding
                'parm05': {},       # Bathymetry initial condition files
            },
        }

        # Some default names for bathymetry files
        self.namepatch['data']['parm04']['usingcartesiangrid'] = True


    def initgrid(self):
        """ Initialize the grid for a two-dimensional model in x, z """

        self.dx = gcmutils.truncate(self.Lx/self.nx) * np.ones((self.nx,))
        self.dy = gcmutils.truncate(self.Ly/self.ny) * np.ones((self.ny,))
        self.dz = gcmutils.truncate(self.Lz/self.nz) * np.ones((self.nz,))

        self.x = self.dx.cumsum()
        self.y = self.dy.cumsum()

        zf = np.concatenate(([0.0], -gcmutils.truncate(self.dz.cumsum())))
        self.z = 0.5*(zf[0:-2] + zf[1:-1])

        self.Y, self.X = np.meshgrid(self.y, self.x)

    
    def updatenamepatch(self):
        """ Push essential fields into the namelist patch. """

        # Update namepatch with required namepatch items...

        self.namepatch['data']['parm04']['delx'] = list(
            gcmutils.truncate(self.dx)) 
        self.namepatch['data']['parm04']['dely'] = list(
            gcmutils.truncate(self.dy)) 
        self.namepatch['data']['parm04']['delz'] = list(
            gcmutils.truncate(self.dz)) 

        # Tref and Sref: search desperately for req. items
        if hasattr(self, 'Tref'):
            self.namepatch['data']['parm01']['Tref'] = list(self.Tref)
        elif hasattr(self, 'ic') and 'T' in self.ic.fields.keys():
            self.namepatch['data']['parm01']['Tref'] = list(
                self.ic.fields['T'].mean(axis=(0, 1)))
        else:
            self.namepatch['data']['parm01']['Tref'] = list(
                np.zeros((self.nz,)))

        if hasattr(self, 'Sref'):
            self.namepatch['data']['parm01']['Sref'] = list(self.Sref)
        elif hasattr(self, 'ic') and 'S' in self.ic.fields.keys():
            self.namepatch['data']['parm01']['Sref'] = list(
                self.ic.fields['S'].mean(axis=(0, 1)))
        else:
            self.namepatch['data']['parm01']['Sref'] = list(
                np.zeros((self.nz,)))

        # Names of the topography file
        if hasattr(self, 'topo'):
            self.namepatch['data']['parm05']['bathyfile'] = self.topofile

        # Names of the binary files that specify initial conditions
        if hasattr(self, 'ic'):
            for var in self.ic.fields.keys():
                self.namepatch['data']['parm05'][
                    self.ic.namelistnames[var]] = (self.ic.filenames[var])

        
        if hasattr(self, 'obcs'):
            self.namepatch['data.pkg'] = { 'packages': { 'useobcs': True } }
            self.namepatch['data.obcs'] = { 'obcs_parm01' : {} }
            for obc in self.obcs.keys():
                for var in self.obcs[obc].fields.keys():
                    self.namepatch['data.obcs']['obcs_parm01'][
                        self.obcs[obc].namelistnames[var]] = (
                        self.obcs[obc].filenames[var] )    


    def saveinput(self):
        """ Save available grid files, initial conditions, and boundary 
        conditions to disk. """


        # Topography
        if hasattr(self, 'topo'):
            topovar = { self.topofile : self.topo }
            gcmutils.savegcminput(topovar, self.setupdirs['inputdir'])

        # Initial condition and boundary conditions
        if hasattr(self, 'ic'):
            self.ic.save(self.setupdirs['inputdir'])

        if hasattr(self, 'obcs'):
            for obc in self.obcs.keys():
                self.obcs[obc].save(self.setupdirs['inputdir'])
        
                
            


    def gensetup(self, namepatch=None, nprun=1, templatedir=None, 
        workdir='.', cleansetup=False):
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

        if not hasattr(self, 'setupdirs'): 
            self.initsetupdirs(workdir=workdir, cleansetup=cleansetup)
    

        # Identify templatedir and daughter code and namelist dirs
        if templatedir is None: templatedir = self.templatedir
        nametempldir = '{}/namelists'.format(templatedir)
        codetempldir = '{}/code'.format(templatedir)


        # Update namepatch and merge with keyword-specific patch
        self.updatenamepatch()
        if namepatch is not None:
            self.namepatch = {**namepatch, **self.namepatch}

        self.saveinput()
        
        # Copy code
        for filename in os.listdir(codetempldir):
            shutil.copy('{}/{}'.format(codetempldir, filename),
                '{}/{}'.format(self.setupdirs['codedir'], filename))


        # Modify SIZE.h
        if self.nx % nprun == 0:
            sizevars = { 'sNx': int(self.nx/nprun), 
                         'Nr' : self.nz,        
                         'nPx': nprun }
        else:
            raise ValueError("Number of grid points must be a mulitple of the "
                "number of run processes.")

        for var, value in sizevars.items():
            gcmutils.changesize(var, value, codedir=self.setupdirs['codedir'])


        # Copy namelists
        for filename in os.listdir(nametempldir):
            if filename not in self.namepatch.keys():
                shutil.copy('{}/{}'.format(nametempldir, filename), 
                    '{}/{}'.format(self.setupdirs['inputdir'], filename)) 


        # Merge namelists, giving priority to namelists in the 
        # dict 'namepatches'.
        namelists = {}
        for nmlfile in self.namepatch.keys():

            namelists[nmlfile] = f90nml.read('{}/{}'.format(
                nametempldir, nmlfile))

            for nml in self.namepatch[nmlfile].keys():
                for var in self.namepatch[nmlfile][nml].keys():
                    namelists[nmlfile][nml][var] = (
                        self.namepatch[nmlfile][nml][var])


        # Save namelists
        for filename in namelists.keys():
            savename = '{}/{}'.format(self.setupdirs['inputdir'], filename)
            with open(savename, 'w') as namefile:
                namelists[filename].write(savename, force=True)

        


    def initsetupdirs(self, workdir='.', cleansetup=False):
        """ Initialize the directory structure of a pymitgcm setup. """

        workdir = os.path.abspath(workdir)

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
        gcmutils.rungcm(self.setupdirs['rundir'], 
            self.setupdirs['builddir'], self.setupdirs['inputdir'])
        print('Run time: {:.3f}'.format(time.time() - starttime))
            


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
