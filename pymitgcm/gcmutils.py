import os
import time, subprocess
import numpy as np
import matplotlib.pyplot as plt
import netCDF4



def rungcm(runpath, buildpath=None, inputpath=None, overwrite=False, 
    outputname='run.info'):
    """ Run a compiled MITgcm model.
        
    Args:
        runpath (str): Path the directory in which to run the model.

        buildpath (str): Path to directory that contains the executable 
            'mitgcmuv'. Assumes runpath if unset.

        inputpath (str): Path to the directory that contains input 
            namelists and binary files. 
    """

    curdir = os.getcwd()

    runpath = os.path.abspath(runpath)
    os.chdir(runpath)

    if buildpath is None: buildpath=runpath

    # Symlink input into run directory
    if inputpath is not None: 
        for filename in os.listdir(inputpath):
            inputfilepath = os.path.join(inputpath, filename)
            runfilepath = os.path.join(runpath, filename)
            if os.path.isfile(inputfilepath):
                if overwrite:
                    try: os.remove(runfilepath)
                    except: pass

                os.symlink(inputfilepath, runfilepath)
                

    with open(os.path.join(runpath, outputname), 'w') as output:
        _, stopmsg = subprocess.Popen(os.path.join(buildpath, 'mitgcmuv'), 
            stdout=output, stderr=subprocess.PIPE).communicate()
            
    os.chdir(curdir)

    return stopmsg





def genmake(gcmpath, optfilename=None, mpi=False, mnc=True, buildpath=None):
    """ Execute MITgcm's 'genmake2' with appropriate compile options in the 
    build directory.
    
    Args:
        gcmpath (str): Path to MITgcm.

        optfilename (str): Name of the optfile. Must be an optfile in 
            gcmpath/tools/build_options.

        mpi (bool): Whether or not to compile with mpi enabled.

        mnc (bool): Whether or not to compile with NetCDF enabled.

        buildpath (str): Directory in which to initiate the build.
    """

    # Change directories if specified
    if buildpath is not None: os.chdir(buildpath)

    # Generate compile command
    compilecmd = ["{}/tools/genmake2".format(gcmpath)]
    options = {}

    options['mods']    = "-mods=../code/"
    options['rootdir'] = "-rootdir={}".format(gcmpath)

    if mnc: options['mnc'] = "-enable=mnc"
    if mpi: options['mpi'] = "-mpi"
        
    if optfilename is not None:
        optfile = os.path.join(gcmpath, 'tools', 'build_options', optfilename)
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


def makedepend(buildpath=None):
    """ Execute 'make depend' in a build directory. 'genmake2' must be
    successfully executed in the directory prior to running this command.

    Args:
        buildpath (str): Build directory where genmake2 was run. If unset
            the current directory is assumed.
    """ 

    if buildpath is not None: os.chdir(buildpath)

    with open('out_makedepend.txt', 'w') as makedependout:
        with open('err_makedepend.txt', 'w') as makedependerr:
            process = subprocess.call(['make', 'depend'],
                stdout=makedependout, stderr=makedependerr)

    if process is not 0:
        with open('err_makedepend.txt', 'r') as errfile:
            for line in errfile:
                print(line.rstrip('\n'))
        raise RuntimeError('make depend failed.')


def make(npmake=1, buildpath=None):
    """ 
        npmake (int): Number of processors to use for compilation.
    """

    if buildpath is not None: os.chdir(buildpath)

    with open('out_make.txt', 'w') as makeout:
        with open('err_make.txt', 'w') as makeerr:
            process = subprocess.call(['make', '-j{}'.format(npmake)],
                stdout=makeout, stderr=makeerr)

    if process is not 0:
        with open('err_make.txt', 'r') as errfile:
            for line in errfile:
                print(line.rstrip('\n'))
        raise RuntimeError('make failed.')


def compile(gcmpath, optfilename=None, npmake=1,
    mpi=False, mnc=True, buildpath=None):
    """ Compile an MITgcm setup using genmake2, make depend, and make. 
    
    Args:
        gcmpath (str): Path to MITgcm.

        optfilename (str): Name of the optfile. Must be an optfile in 
            gcmpath/tools/build_options.

        npmake (int): Number of processors to use for compilation.

        mpi (bool): Whether or not to compile with mpi enabled.

        mnc (bool): Whether or not to compile with NetCDF enabled.

        buildpath (str): Directory in which to initiate the build.
    """

    if buildpath is not None: os.chdir(buildpath)

    genmake(self.gcmpath, optfilename=optfilename, mpi=mpi, mnc=mnc)
    makedepend()
    make(npmake=npmake)




def readsizevars(*args):
    """Change a variable in the header file SIZE.h.
        
    Args:
        sizepath (str): Path to the 'code' directory where the new SIZE.h will 
            be parsed.

    Returns: A dictionary with parsed variables.
    """

    if len(args) > 1:
        raise TypeError("Expected at most 1 argument, got %d" % len(args))
    elif len(args) is 0:
        sizepath = os.getcwd()
    else:
        sizepath = os.path.abspath(args[0])

    sizename = 'SIZE.h'
    varnames = ['sNx', 'sNy', 'OLx', 'OLy', 'nSx', 'nSy', 
        'nPx', 'nPy', 'Nr']

    # Read file
    with open(os.path.join(sizepath, sizename), 'r') as sizefile:
        sizelines = sizefile.readlines()

    sizevars = {}
    for line in sizelines:
        for var in varnames:
            if (len(line) > 1 and 
                line.lstrip()[0] is not 'C' and 
                'integer' not in line.lower() and
                'max_olx' not in line.lower() and
                'max_oly' not in line.lower() and
                var.lower() in line.lower()):

                # Find text between '=' and ',' or ')' and assign to sizevars.
                sep = '|'
                line = line.replace(' ', '').lower()
                for char in ['=', ',', ')']: line=line.replace(char, sep)

                try:    sizevars[var] = int(line.split(sep)[1])
                except: pass
                    
    return sizevars



def changesizevars(vars, sizepath=None):
    """Change a variable in the header file SIZE.h.
        
    Args:
        vars (dict): Dictionary of the form {varname: value}. The varname
            must be contained in a small list of 'validvars'.
        
        sizepath (str): Path to the 'code' directory where the new SIZE.h will 
            be saved.
    """

    sizename = 'SIZE.h'
    prespace = '     &           '
    validvarnames = ['sNx', 'sNy', 'OLx', 'OLy', 'nSx', 'nSy', 
        'nPx', 'nPy', 'Nr']

    if sizepath is None: sizepath = os.getcwd()

    # Read file
    with open(os.path.join(sizepath, sizename), 'r') as sizefile:
        sizelines = sizefile.readlines()

           

    for varname, value in vars.items():

        if varname not in validvarnames:
            raise ValueError("The parameter 'varname' {} must be one of "
                "{}!".format(
                varname, validvarnames))
        elif type(value) is not int:
            raise ValueError("The parameter 'value' must be an integer.")

        # Overwrite the line in sizelines containing a non-commented
        # and valid reference to 'var'
        linenum = 0
        for line in sizelines:
            if (len(line) > 1 and 
                line.lstrip()[0] is not 'C' and 
                '{:<3s} ='.format(varname) in line):
                if varname is 'Nr':
                    sizelines[linenum] = '{}{:<3s} ={:>4d})\n'.format(
                        prespace, varname, value)
                else:
                    sizelines[linenum] = '{}{:<3s} ={:>4d},\n'.format(
                        prespace, varname, value)
            linenum += 1

    with open(os.path.join(sizepath, sizename), 'w') as sizefile:
        for line in sizelines:
            sizefile.write(line)



def changesizevar(var, value, codepath='.', origpath=None):
    """Change a variable in the header file SIZE.h.
        
    Args:
        var (str): Name of the variable to be changed. This variable must be
            contained in a small list of 'validvars'.
        
        value (int): Value of the variable to be changed. Because it 
            is always associated with some aspect of the grid size, it must 
            be an integer.

        codepath (str): Path to the 'code' directory where the new SIZE.h will 
            be saved.

        origpath (str): Path to the template directory with an original copy of 
            SIZE.h. If this path is not specified, origpath is set to codepath.
            
            
    """

    sizename = 'SIZE.h'
    prespace = '     &           '
    validvars = ['sNx', 'sNy', 'OLx', 'OLy', 'nSx', 'nSy', 
        'nPx', 'nPy', 'Nr']

    if origpath is None:
        origpath = codepath

    if var not in validvars:
        raise ValueError("The parameter 'var' {} must be one of {}!".format(
            var, validvars))
    elif type(value) is not int:
        raise ValueError("The parameter 'value' must be an integer.")

    # Make output 'code' directory if it does not exist
    if not os.path.exists(codepath):
        os.makedirs(codepath)

    # Read file
    with open(os.path.join(origpath, sizename), 'r') as sizefile: 
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

    with open(os.path.join(codepath, sizename), 'w') as sizefile:
        for line in sizelines:
            sizefile.write(line)
            
    

def get_namelist_param(param, nmlpath='.'):
    """Returns the value of a parameter extracted from within the namelists 
    in nmlpath."""

    # Read namelists
    namefile = {}
    for filename in os.listdir(nmlpath):
        namefiles[filename] = f90nml.read(os.path.join(nmlpath, filename))

    # Sift namelists to find param
    return sift_nested_dict(param, namefiles)
    


def siftdict(param, d, value=None, dpath=None, level=1):
    """Sift recursively and case-inesensitively through a nested 
    dictionary in search of the key 'param' and either return its 
    value or None if not found. """

    if dpath is None: dpath=[]

    lowered = {k.lower():v for k, v in d.items()}
    if param.lower() in lowered.keys():
        value = lowered[param.lower()]
        dpath.append(param.lower())
    else: 
        for k, v in d.items():
            if isinstance(v, dict) and value is None:
                dpath.append(k)
                value, dpath, level = siftdict(param, v, value=value,
                    level=level+1, dpath=dpath)

    if value is None:
        return value, dpath[:-1], level-1
    else:
        return value, dpath, level


def savegcminput(savevars, savepath='.'):
    """Convert and save variables in proper format for input into
    MITgcm.

    Args:
        savevars (dict): A dictionary of variables in the form savename: var.
        
        savepath (str):  A string specifying the directory in which to save
            the binary files.
    """

    # Convert variables to big endian, double precision for saving.
    savetype = '>f8'

    for savename in savevars.keys():

        savevars[savename] = savevars[savename].astype(savetype) 
        varshape = savevars[savename].shape

        print(savename)
        print(varshape)
        print(len(varshape))

        with open(os.path.join(savepath, savename), 'wb') as file: 
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
    newnamepath = os.path.join(savepath, os.path.basename(namepath))
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
    diagfilename = os.path.join(savepath, 'data.diagnostics')
    ndiags = len(fields.keys())

    # Check to see if file exists
    if os.path.isfile(diagfilename):
        if overwrite:
            os.remove(diagfilename)
        else:
            raise ValueError("File {} exists! "
                    "Either delete it or set overwrite=True.".format(
                    diagfilename))

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
    """ Return the array a, truncated to a specified number of 
    digits. """
    return np.round(np.round(a
        / 10**np.floor(np.log10(np.abs(a)+1e-15)-digits))
        * 10**np.floor(np.log10(np.abs(a)+1e-15)-digits), 14)



def quicklook(runpath, var='Temp'):
    """ Show rudimentary plots of the model state. """

    grid = netCDF4.Dataset(os.path.join(runpath,
        'grid.t{:03d}.nc'.format(1)), 'r') 
    
    state = netCDF4.Dataset(os.path.join(runpath,
        'state.{:010d}.t{:03d}.nc'.format(0, 1)), 'r')

    nt  = state.variables[var][:, 0, 0, 0].size
    phi = state.variables[var][0, :, 0, :].squeeze()

    fig = plt.figure()
    plt.imshow(phi)
    plt.colorbar()
        
    for i in range(nt):
        print("Savetime %d".format(i))

        phi = state.variables[var][i, :, 0, :].squeeze()

        plt.clf()
        plt.imshow(phi)
        plt.colorbar()

        plt.pause(0.1)
