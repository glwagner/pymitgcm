import os
import numpy as np
import time, subprocess



def rungcm(rundir, builddir, inputdir):
    """ Run a compiled MITgcm model.
        
    Args:
        rundir (str): Path the directory in which to run the model.

        builddir (str): Path to the build directory of the model that 
            contains the executable 'mitgcmuv'.

        inputdir (str): Path to the directory that contains input 
            namelists and binary files.
    """

    for filename in os.listdir(rundir):
        os.remove('{}/{}'.format(rundir, filename))

    for filename in os.listdir(inputdir):
        os.symlink(
            '{}/{}'.format(inputdir, filename),
            '{}/{}'.format(rundir, filename) 
        )

    os.chdir(rundir)

    with open('./output.txt', 'w') as output:
        _, stopmsg = subprocess.Popen('{}/mitgcmuv'.format(
            builddir), stdout=output, 
            stderr=subprocess.PIPE).communicate()
            

    if b'STOP NORMAL END' not in stopmsg:
        print(stopmsg)




def genmake(gcmdir, optfilename=None, mpi=False, mnc=True, builddir=None):
    """ Execute MITgcm's 'genmake2' with appropriate compile options in the 
    build directory.
    
    Args:
        gcmdir (str): Path to MITgcm.

        optfilename (str): Name of the optfile. Must be an optfile in 
            gcmdir/tools/build_options.

        mpi (bool): Whether or not to compile with mpi enabled.

        mnc (bool): Whether or not to compile with NetCDF enabled.

        builddir (str): Directory in which to initiate the build.
    """

    # Change directories if specified
    if builddir is not None: os.chdir(builddir)

    # Generate compile command
    compilecmd = ["{}/tools/genmake2".format(gcmdir)]
    options = {}

    options['mods'] = "-mods=../code/"
    options['rootdir'] = "-rootdir={}".format(gcmdir)

    if mnc: options['mnc'] = "-enable=mnc"
    if mpi: options['mpi'] = "-mpi"
        
    if optfilename is not None:
        optfile = '{}/tools/build_options/{}'.format(gcmdir, optfilename)
        if os.path.exists(optfile):
            options['optfile'] = "-optfile={}".format(optfile)

    for opt in options.values():
        compilecmd.append(opt)

    with open('out_genmake2.txt', 'w') as genmakeout:
        with open('err_genmake2.txt', 'w') as genmakeerr:
            process = subprocess.call(compilecmd, stdout=genmakeout, 
                stderr=genmakeerr)

    if process is not 0:
        with open('err_genmake2.txt', 'r') as errfile:
            for line in errfile:
                print(line.rstrip('\n'))
        raise RuntimeError('Genmake2 failed.')


def makedepend(builddir=None):
    """ Execute 'make depend' in a build directory. 'genmake2' must be
    successfully executed in the directory prior to running this command.

    Args:
        builddir (str): Build directory where genmake2 was run. If unset
            the current directory is assumed.
    """ 

    if builddir is not None: os.chdir(builddir)

    with open('out_makedepend.txt', 'w') as makedependout:
        with open('err_makedepend.txt', 'w') as makedependerr:
            process = subprocess.call(['make', 'depend'],
                stdout=makedependout, stderr=makedependerr)

    if process is not 0:
        with open('err_makedepend.txt', 'r') as errfile:
            for line in errfile:
                print(line.rstrip('\n'))
        raise RuntimeError('make depend failed.')


def make(npmake=1, builddir=None):
    """ 
        npmake (int): Number of processors to use for compilation.
    """

    if builddir is not None: os.chdir(builddir)

    with open('out_make.txt', 'w') as makeout:
        with open('err_make.txt', 'w') as makeerr:
            process = subprocess.call(['make', '-j{}'.format(npmake)],
                stdout=makeout, stderr=makeerr)

    if process is not 0:
        with open('err_make.txt', 'r') as errfile:
            for line in errfile:
                print(line.rstrip('\n'))
        raise RuntimeError('make failed.')


def compile(gcmdir, optfilename=None, npmake=1,
    mpi=False, mnc=True, builddir=None):
    """ Compile an MITgcm setup using genmake2, make depend, and make. 
    
    Args:
        gcmdir (str): Path to MITgcm.

        optfilename (str): Name of the optfile. Must be an optfile in 
            gcmdir/tools/build_options.

        npmake (int): Number of processors to use for compilation.

        mpi (bool): Whether or not to compile with mpi enabled.

        mnc (bool): Whether or not to compile with NetCDF enabled.

        builddir (str): Directory in which to initiate the build.
    """

    if builddir is not None: os.chdir(builddir)

    genmake(self.gcmdir, optfilename=optfilename, mpi=mpi, mnc=mnc)
    makedepend()
    make(npmake=npmake)



def changesize(var, value, codedir='.', origdir=None):
    """Change a variable in the header file SIZE.h.
        
    Args:
        var (str): Name of the variable to be changed. This variable must be
            contained in a small list of 'validvars'.
        
        value (int): Value of the variable to be changed. Because it 
            is always associated with some aspect of the grid size, it must 
            be an integer.

        codedir (str): Path to the 'code' directory where the new SIZE.h will 
            be saved.

        origdir (str): Path to the template directory with an original copy of 
            SIZE.h. If this path is not specified, origdir is set to codedir.
            
            
    """

    sizename = 'SIZE.h'
    prespace = '     &           '
    validvars = ['sNx', 'sNy', 'OLx', 'OLy', 'nSx', 'nSy', 
        'nPx', 'nPy', 'Nr']

    if origdir is None:
        origdir = codedir

    if var not in validvars:
        raise ValueError("The parameter 'var' {} must be one of {}!".format(
            var, validvars))
    elif type(value) is not int:
        raise ValueError("The parameter 'value' must be an integer.")

    # Make output 'code' directory if it does not exist
    if not os.path.exists(codedir):
        os.makedirs(codedir)

    # Read file
    with open('{}/{}'.format(origdir, sizename), 'r') as sizefile:
        sizelines = sizefile.readlines()

    linenum = 0
    for line in sizelines:
        if (len(line) > 1 and 
            line.lstrip()[0] is not 'C' and 
            '{:<3s} ='.format(var) in line):
            if var is 'Nr':
                sizelines[linenum] = '{}{:<3s} ={:>4d})\n'.format(
                    prespace, var, value)
            else:
                sizelines[linenum] = '{}{:<3s} ={:>4d},\n'.format(
                    prespace, var, value)
        linenum += 1

    with open('{}/{}'.format(codedir, sizename), 'w') as sizefile:
        for line in sizelines:
            sizefile.write(line)
            
    

def get_namelist_param(param, nmldir='.'):
    """Returns the value of a parameter extracted from within the namelists 
    in nmldir."""

    # Read namelists
    namefile = {}
    for filename in os.listdir(nmldir):
        namefiles[filename] = f90nml.read('{}/{}'.format(nmldir, filename))

    # Sift namelists to find param
    return sift_nested_dict(param, namefiles)
    

def siftdict(param, d, value=None):
    """Sift recurively through a nested dictionary in search of 
    the key 'param' and either return its value or throw an error."""

    if param.lower() in d.keys():
        value = d[param]
    else: 
        for k, v in d.items():
            if isinstance(v, dict):
                value = sift_nested_dict(param, v, value=value)
        
    return value


    

def savegcminput(savevars, savedir='.'):
    """Convert and save variables in proper format for input into
    MITgcm.

    Args:
        savevars (dict): A dictionary of variables in the form savename: var.
        
        savedir (str):  A string specifying the directory in which to save
            the binary files.
    """

    # Convert variables to big endian, double precision for saving.
    savetype = '>f8'

    for savename in savevars.keys():
        savevars[savename] = savevars[savename].astype(savetype) 
        varshape = savevars[savename].shape

        with open('{}/{}'.format(savedir, savename), 'wb') as file:

            if len(varshape) is 1:
                savevars[savename].tofile(file)            

            elif len(varshape) is 2:
                for k in range(varshape[1]):
                    savevars[savename][:, k].tofile(file)

            elif len(varshape) is 3:
                for k in range(varshape[2]):
                    for l in range(varshape[1]):
                        savevars[savename][:, l, k].tofile(file)


def get_diag_list_header():
    """Return the standard comment header to the &diagnostics_list namelist
    in data.diagnostics as a string."""

    header = (
        "# Diagnostic Package Choices\n"
        "#-------------------------------------------------------------------\n"
        "# for each output-stream:\n"
        "#\n"
        "#  filename(n) : prefix of the output file name (only 8.c long) \n"
        "#      for outp.stream n\n"
        "#\n"
        "#  frequency(n):< 0 : write snap-shot output every multiple of \n"
        "#      |frequency| (iter)\n"
        "#               > 0 : write time-average output every multiple of \n"
        "#      frequency (iter)\n"
        "#\n"
        "#  levels(:,n) : list of levels to write to file (Notes: declared \n"
        "#      as REAL) when this entry is missing, select all common levels\n" 
        "#      of this list\n"
        "#\n"
        "#  fields(:,n) : list of diagnostics fields (8.c) \n"
        "#      (see 'available_diagnostics' file for the list of all\n"
        "#      available diag. in this particular config)\n"
        "#\n"                 
        "#-------------------------------------------------------------------"
    )

    return header


def get_diag_stat_header():
    """Return the standard comment header to the &diag_statis_param namelist
    in data.diagnostics as a string."""

    header = (
        "# Parameter for Diagnostics of per level statistics:\n"
        "#-------------------------------------------------------------------\n"
        "# for each output-stream:\n"
        "#\n"
        "#  stat_fname(n) : prefix of the output file name (only 8.c long) \n"
        "#      for outp.stream n\n"
        "#\n"
        "#  stat_freq(n):< 0 : write snap-shot output every |stat_freq| \n"
        "#                     seconds\n"
        "#               > 0 : write t-average output every stat_freq seconds\n"
        "#\n"
        "#  stat_phase(n)    : write at time = stat_phase + multiple of \n"
        "#                     |stat_freq|\n"
        "#\n"
        "#  stat_region(:,n) : list of 'regions' (deft: 1 region only=global)\n"
        "#\n"
        "#  stat_fields(:,n) : list of diagnostics fields (8.c) (see \n"
        "#      'available_diagnostics.log' file for the list of all \n"
        "#      available diag. in this particular config)\n"
        "#--------------------------------------------------------------------"
    )

    return header 


 
def replace_namelist_param(namepath, params, savepath='.'):
    """ Find and replace a parameter in the MITgcm namelist
    specified by datapath. 

    Args:
        namepath (str): Path to the namelist to be changed.
    
        params (dict): Dictionary containing the parameter:value pairs to be
            replaced.

        saveorig (bool): A boolean specifying whether or not to save the
            original file.
    """

    # Characters that identify the end of a parameter specification
    endchars = ['#', '&', '/'] 
    validtypes = [bool, str, int, float, list]
    linewidth = 8

    # Check for errors
    if not os.path.isfile(namepath):
        raise ValueError("The file {} does not exist!".format(namepath))
    elif type(paramval) not in validtypes:
        raise ValueError("The parameter type is not valid")

    # Determine whether original will be replaced or not
    if savepath is os.path.dirname(namepath):
        replaceorig = True
    else:
        replaceorig = False

    # Read the file into a list of strings
    with open(namepath, 'r') as namefile:
        namelist = namefile.readlines()

    # Loop over items in iteritems
    oldnamelist = list(namelist)
    for paramname, paramval in params.items():

        # Convert numpy array to list 
        if type(paramval) is np.ndarray:
            if len(paramval.shape) > 1:
                raise ValueError("The numpy array that specifies the "
                        "parameter value can only have one dimension.")

            paramval = paramval.tolist()

        # Find the first line where the parameter is specified
        paramline, linenum = None, 0
        for line in namelist:
            line = line.lstrip()
            if len(line) > 0:
                if (line[0] not in endchars and 
                    paramname in line and 
                    '=' in line):
                    if paramline is not None:
                        raise RuntimeError("Two lines have been found that \n"
                            "contain the parameter name, the equals sign \n"
                            "'=', and do not appear to be a comment. Check \n"
                            "the input data file and remove the duplicate or "
                            "otherwise offensive line.")
                    else:
                        paramline = linenum

            linenum += 1 

        # For multiline parameters, find last line where parameter is specified
        (linenum, endfound) = (paramline, False)
        while not endfound:
            linenum += 1
            # Determine if the end has been found
            line = namelist[linenum].lstrip()
            if len(line) == 0: 
                endfound = True
            elif line[0] in endchars or '=' in line:
                endfound = True
            elif linenum == len(namelist):
                raise RuntimeError("The end of the apparently multi-line \n"
                        "parameter specification for {} \n".format(paramname) +
                        "could not be found.")

        paramendline = linenum

        # Generate the new namelist line by line:
        newnamelist = []

        #   1. Write unchanged lines 
        for line in range(paramline):
            newnamelist.append(namelist[line])

        #   2. Generate the new parameter lines
        if type(paramval) is str:
            newtext = " {} = '{}',\n".format(paramname, paramval)

        elif type(paramval) is bool:
            if paramval is False: paramtext = '.FALSE.'
            else: paramtext = '.TRUE.'
            newtext = " {} = {},".format(paramname, paramtext)

        elif type(paramval) in [int, float]:
            newtext = " {} = {},\n".format(paramname, paramval)

        elif type(paramval) is list:
            lineindent = 4 + len(paramname)
            newtext = " {} = ".format(paramname)
            for i in range(len(paramval)):
                newtext = newtext + "{:6}, ".format(paramval[i])
                if (i+1) % linewidth == 0 and (i+1) != len(paramval):
                    newtext = newtext + "\n{}".format(' '.ljust(lineindent))

            newtext = newtext + "\n"    

        newnamelist.append(newtext)

        #   3. Write the rest of the file
        for line in range(paramendline, len(namelist)):
            newnamelist.append(namelist[line])

        # Save new namelist for next iteration.
        namelist = list(newnamelist)

    # If replacing the original, remove original and save in renamed file
    if replaceorig: 
        os.remove(namepath)
        with (namepath + "_orig", 'w') as namefile:
            for line in oldnamelist:
                namefile.write(line)

    # Save the file
    newnamepath = '{}/{}'.format(savepath, os.path.basename(namepath))
    with open(newnamepath, 'w') as newnamefile:
        for line in newnamelist:
           newnamefile.write(line) 


def write_data_diagnostics(fields, freqs, levels, 
    savepath='.', overwrite=False):
    """ Write a diagnostics file for MITgcm diagnostics.
    Each diagnostic is defined by a list of field names, frequency
    of output, vertical levels in the case of 3D diagnostics, and a 
    filename. 

    Args:
        fields[diagname] : List of fields in the diagnostic.

        freqs[diagname]  : Number (float) giving the frequency at which 
            the diagnostic is outputted.

        levels[diagname] : For 3D diagnostics, the vertical levels from which
            the diag will be extracted.

        savepath (str)   : Path in which to save output.

        overwrite (bool) : Flag to indicate whether or not to overwrite 
            an existing file at savepath.
    """
                
    # Parameters
    levelLinewidth = 10
    fieldLinewidth = 3
    diagFilename = "{}/data.diagnostics".format(savedir)
    ndiags = len(fields.keys())

    # Check to see if file exists
    if os.path.isfile(diagFilename):
        if overwrite:
            os.remove(diagFilename)
        else:
            raise ValueError("File {} exists! "
                    "Either delete it or set overwrite=True.".format(
                    diagFilename))

    # Header to the diagnostics_list namelist
    listheader = get_diag_list_header()

    # Construct body
    listbody = (
        " &diagnostics_list\n"
        "\n"
        " dumpatlast = .TRUE.,\n"
        "\n"
    )

    i = 0
    for diag in fields.keys():

        i += 1
        nfields = len(fields[diag])

        filetext  = " filename({:d}) = 'diags/{}',\n".format(i, diag)
        freqtext  = " frequency({:d}) = {:.2f},\n".format(i, freqs[diag])

        # Construct text block for field data
        fieldtext = " fields(1:{:d},{:d}) = ".format(nfields, i)
        for j in range(nfields):
            fieldtext = fieldtext + "'{}', ".format(fields[diag][j])
            if (j+1) % fieldLinewidth == 0 and (j+1) != nfields:
                fieldtext = fieldtext + "\n{:16}".format(' ')

        fieldtext = fieldtext + "\n"

        # Construct text block for level data
        if levels[diag] is None: 
            leveltext = ''
        else:
            nlevels = len(levels[diag])
            leveltext = " levels(1:{:d},{:d}) = ".format(nlevels, i)
            for j in range(nlevels):
                leveltext = leveltext + "{:3d}., ".format(levels[diag][j])
                if (j+1) % levelLinewidth == 0 and (j+1) != nlevels:
                    leveltext = leveltext + "\n{:18}".format(' ')
            leveltext = leveltext + "\n"

        listbody = listbody + filetext + freqtext + fieldtext + leveltext

        if i != ndiags-1:
            listbody = listbody + "\n\n"
        else:
            listbody = listbody + " /\n"


    # Header to the diag_statis_parms namelist.
    statheader = get_diag_stat_header()

    # This function does not support a definition of parameters in the 
    # diag_statis_parms namelist
    statbody = ''

    # Construct full file text for data.diagnostics
    fulldiagtext = (listheader + '\n' + listbody + '\n' 
        + statheader + statbody)

    # Write the file
    with open(savepath, 'w') as diagfile:
        diagfile.write(fulldiagtext)


def truncate(a, digits=3):
    """ Truncate the array "a" to an integer number of digits. """
    atrunc = (np.round(a
        / 10**(np.floor(np.log10(a))-digits))
        * 10**(np.floor(np.log10(a))-digits))
    return atrunc

