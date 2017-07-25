import os
import numpy as np
import f90nml



class Namelist:
    def __init__(self, namepath):
        """An object for storing, reading, and writing namelist 
        objects associated with MITgcm simulations."""

        self.base  = f90nml.read(namepath + '/data')
        self.cal   = f90nml.read(namepath + '/data.cal')
        self.exch2 = f90nml.read(namepath + '/data.exch2')
        self.exf   = f90nml.read(namepath + '/data.exf')
        self.kpp   = f90nml.read(namepath + '/data.kpp')
        self.obcs  = f90nml.read(namepath + '/data.obcs')
        self.pkg   = f90nml.read(namepath + '/data.pkg')
        self.ee    = f90nml.read(namepath + '/eedata')

        self.namepath = namepath

        
        


def edit_namelist(namelistobj, namelist, params):
    """Edits a namelist inside an f90nml namelist object to either change
    an existing namelist parameter or add a new one if it does not exist.

    Args:
        namelistobj: The f90nml namelist object to be edited.

        namelist (str): The namelist to be edited. The namelist must be found
            inside the namelistobj.

        params (dict): A dictionary of parameters to be changed with the specified
            namelist.

    Returns:
        An edited f90nml namelist object.
    """

    
def replace_namelist_param(namepath, param, saveorig=True):
    
    """ Find and replace a parameter in the MITgcm namelist
    specified by datapath. The input param is a dictionary whose
    key gives the parameter and whose value gives the value.
    This program automatically determines and formats the value
    according to MITgcm's namelist requirements."""

    # Characters that identify the end of a parameter specification
    endchars = ['#', '&', '/'] 
    validtypes = [bool, str, int, float, list]
    linewidth = 8

    paramname = param.keys()[0]
    paramval = param.values()[0]

    # Convert numpy array to list 
    if type(paramval) is np.ndarray:
        if len(paramval.shape) > 1:
            raise ValueError("The numpy array that specifies the "
                    "parameter value can only have one dimension.")
        paramval = paramval.tolist()

    if not os.path.isfile(namepath):
        raise ValueError("The file {} does not exist!".format(namepath))

    if type(paramval) not in validtypes:
        raise ValueError("The parameter type is not valid")

    # Read the file into a list of strings
    with open(namepath, 'r') as namefile:
        namelist = namefile.readlines()

    # Find the first line where the parameter is specified
    linenum = 0
    for line in namelist:
        line = line.lstrip()
        if len(line) > 0:
            if line[0] not in endchars and paramname in line and '=' in line:
                if 'paramline' in locals():
                    raise RuntimeError( 
                        "Two lines have been found that contain the parameter name, \n"
                        "the equals sign '=', and do not appear to be a comment. \n"
                        "Check the input data file and remove the duplicate or "
                        "otherwise offensive line.")
                else:
                    paramline = linenum

        linenum += 1 

    # For multiline parameters, find the last line where the parameter is specified
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
            raise RuntimeError("The end of the apparently multi-line parameter\n"
                    "specification could not be found.")
    paramendline = linenum

    # Save original file
    with open(namepath + '_orig', 'w') as namefile:
        for line in range(len(namelist)):
            namefile.write(namelist[line])

    # Write new file
    with open(namepath, 'w') as namefile:

        # Write the unmodified beginning lines
        for line in range(paramline):
            namefile.write(namelist[line])

        # Generate the new parameter lines
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

        # Write the new text into file
        print("The text replacing the parameter {} is".format(paramname)
            + "\n\n{}".format(newtext)
        )
        namefile.write(newtext)

        # Write the rest of the file
        for line in range(paramendline, len(namelist)):
            namefile.write(namelist[line])


def write_data_diagnostics(fields, freqs, levels, savepath='.', overwrite=False):
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

        overwrite (bool) : Flag to indicate whether or not to overwrite an existing
            file at savepath.
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
                    "Either delete it or set overwrite=True.".format(diagFilename))

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



