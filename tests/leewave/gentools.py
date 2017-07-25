import os
import numpy as np


def change_size_var(var, value, codedir='.', origdir=None):
    """Change a variable in the header file SIZE.h.
        
    Args:
        var (str): Name of the variable to be changed. This variable must be
            contained in a small list of 'validvars'.
        
        value (int): Value of the variable to be changed. Because it 
            is always associated with some aspect of the grid size, it must 
            be an integer.

        codedir (str): Path to the 'code' directory where SIZE.h is to be kept.

        origdir (str): Path to the template directory where an example SIZE.h
            is kept.
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
    with open('{}/{}'.format(templdir, sizename), 'r') as sizefile:
        sizelines = sizefile.readlines()

    linenum = 0
    for line in sizelines:
        if line.lstrip()[0] is not 'C' and var in line:
            sizelines[linenum] = '{}{:<3s} ={:>4d},\n'.format(
                prespace, var, value)
        linenum += 1

    with open('{}/{}'.format(codedir, sizename), 'w') as sizefile:
        for line in sizelines:
            sizefile.write(line)
            
    

def get_param(param, namelists):
    """Returns the value of a parameter extracted from within the namelists 
    in nmldir."""

    # Read namelists
    namefile = {}
    for filename in os.listdir(nmldir):
        namefiles[filename] = f90nml.read('{}/{}'.format(nmldir, filename))

    # Sift namelists to find param
    return sift_nested_dict(param, namefiles)
    

def sift_nested_dict(param, d, value=None):
    """Sift recurively through a nested dictionary in search of 
    the key 'param' and either return its value or throw an error."""

    if param.lower() in d.keys():
        value = d[param]
    else: 
        for k, v in d.iteritems():
            if isinstance(v, dict):
                value = sift_nested_dict(param, v, value=value)
        
    return value


    

def convert_and_save(vars, savedir='.'):
    """Convert and save variables in proper format for input into
    MITgcm.

    Args:
        vars (dict):    A dictionary of variables.
        
        savedir (str):  A string specifying the directory in which to save
            the binary files.
    """

    # Convert variables to big endian, double precision for saving.
    savetype = '>f8'

    for var in vars.keys():
        vars[var] = vars[var].astype(savetype) 
        varshape = vars[var].shape

        with open('{}/{}.bin'.format(savedir, var), 'wb') as file:
            if len(varshape) is 1:
                vars[var].tofile(file)            
            elif len(varshape) is 2:
                for k in range(varshape[1]):
                    vars[var][:, k].tofile(file)


