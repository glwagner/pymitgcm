import os, glob, re
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

from functools import reduce
from operator import mul

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
# Plotting
def getoutputiters(runpath):
    """ Retrieve the number of time slices saved in a NetCDF output file. """
    state = netCDF4.Dataset(os.path.join(runpath,
        'state.{:010d}.t{:03d}.nc'.format(0, 1)), 'r')
    return len(state.variables['T'])


def getgrid(runpath):

    state = netCDF4.Dataset(os.path.join(runpath,
        'state.{:010d}.t{:03d}.nc'.format(0, 1)), 'r')

    x = state.variables['X'][:]
    y = state.variables['Y'][:]
    z = state.variables['Z'][:]

    xg = state.variables['Xp1'][:]
    yg = state.variables['Yp1'][:]

    t = state.variables['T'][:]

    return x, y, z, t, xg, yg


def getoutputstate(runpath, iter=0, fields='all'):

    # Somehow reconstruct full fields when multiple processors are used.
    #nprun = len(glob.glob('state.{:010d}.t*.nc'.format(0)))

    fieldkey = {
        't': 'T',
        'u': 'U',
        'v': 'V',
        'w': 'W',
        'T': 'Temp',
        'S': 'S',
        'h': 'Eta',
    }

    if fields is 'all':
        fields = ['u', 'v', 'w', 'T', 'S', 'h']

    state = netCDF4.Dataset(os.path.join(runpath,
        'state.{:010d}.t{:03d}.nc'.format(0, 1)), 'r')

    vars = {}
    for fld in fields:
        vars[fld] = state.variables[fieldkey[fld]][iter]

    t = state.variables['T'][iter]

    return t, vars


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
        logger.info("Savetime {:d}".format(i))

        phi = state.variables[var][i, :, 0, :].squeeze()

        plt.clf()
        plt.imshow(phi)
        plt.colorbar()

        plt.pause(0.1)





# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
# Copied from mds.py
debug = False

_typeprefixes = {   
    'ieee-be': '>',
    'b'      : '>',
    '>'      : '>',
    'ieee-le': '<',
    'l'      : '<',
    '<'      : '<',
}

_typesuffixes = {
    'float32': 'f4',
    'float64': 'f8',
}

_comment_pattern = re.compile(
    r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
    re.DOTALL | re.MULTILINE
)

_closing = {'[':']',
            '{':'}',
}

_string_pattern = re.compile(r"'(.*)'$")

# These deal with comments in the metafile:
def _comment_replacer(match):
    s = match.group(0)
    if s.startswith('/'): return ''
    else:                 return s
        

def strip_comments(text):
    """ Strips C and C++ style comments from text """
    return re.sub(_comment_pattern, _comment_replacer, text)


def parse1(s):
    """ Convert one item to appropriate type """
    m = _string_pattern.match(s)
    if m:
        s = m.group(1)
        # unquote quotes
        s = re.sub(r"''","'",s)
        return s

    if '.' in s or 'e' in s.lower():
        return float(s)
    else:
        try:
            return int(s)
        except ValueError:
            raise ParseError("Cannot parse value: " + s)


def aslist(listlikeitem):
    """ If listlikeitem is iterable, convert to list. If not, form new 
    single-element list with item as first value. """
    try:              listitem = list(listlikeitem)
    except TypeError: listitem = [listlikeitem]
    return listitem


def fromfileshape(filename,dtype, shape=None, **kwargs):
    return np.fromfile(filename, dtype, **kwargs).reshape(shape)


def rdmds(fnamearg, itrs=-1, machineformat='b', rec=None, fill_value=0,
          returnmeta=False, astype=float, region=None, lev=(),
          usememmap=False, mm=False, squeeze=True, verbose=False):
    """ 
    a = rdmds(fname,...)
    a = rdmds(fname,itrs,...)
    a, its, meta = rdmds(fname,...,returnmeta=True)

    Read meta-data files as written by MITgcm.

    Without itrs, will try to read

      fname.meta or fname.001.001.meta, ...

    If itrs is a list of integers or an integer, it will read the corresponding

      fname.000000iter.meta, ...

    If itrs is NaN, it will read all iterations for which files are found.
    If itrs is Inf, it will read the highest iteration found.

    fname may contain shell wildcards, which is useful for tile files organized
    into directories, e.g.,

      T = rdmds('prefix*/T', 2880)

    will read prefix0000/T.0000002880.*, prefix0001/T.0000002880.*, ...
    (and any others that match the wildcard, so be careful how you name things!)

    Returns:

        a    :: numpy array of the data read
        its  :: list of iteration numbers read (only if itrs=NaN or Inf)
        meta :: dictionary of metadata (only if returnmeta=True)

    Keyword arguments:

        machineformat :: endianness ('b' or 'l', default 'b')
        rec           :: list of records to read (default all)
                         useful for pickups and multi-field diagnostics files
        fill_value    :: fill value for missing (blank) tiles (default 0)
        astype        :: data type to return (default: double precision)
                         None: keep data type/precision of file
        region        :: read only (x0,x1,y0,y1). Default (0,nx,0,ny)
        lev           :: list of levels to read, or, for multiple dimensions
                         (excluding x,y), tuple(!) of lists (see examples below)
        usememmap     :: if True, use memory map for reading data (default False)
                         recommended when using lev, or region with global files
                         to save memory and, possibly, time

    Examples:

        XC = rdmds('XC')
        XC = rdmds('res_*/XC')
        T = rdmds('T.0000002880')
        T = rdmds('T', 2880)
        T2 = rdmds('T', [2880,5760])
        T, its = rdmds('T', numpy.Inf)
        VVEL = rdmds('pickup', 2880, rec=range(50,100))
        a5 = rdmds('diags', 2880, rec=0, lev=[5])
        a = rdmds('diags', 2880, rec=0, lev=([0], [0,1,5,6,7]))

        from numpy import r_
        a = rdmds('diags', 2880, rec=0, lev=([0], r_[:2, 5:8])) 
        a = rdmds('diags', 2880, rec=0)[0, [0, 1, 5, 6, 7], ...]
        a = rdmds('diags', 2880)[0, 0, [0, 1, 5, 6, 7], ...]  
    """
    usememmap = usememmap or mm
    if usememmap:
        readdata = np.memmap
    else:
        readdata = fromfileshape

    # add iteration number to file name unless itrs is -1
    additrs = itrs != -1
    if itrs is np.nan:
        # all iterations
        itrs = scanforfiles(fnamearg)
        if verbose: warning('Reading {0} time levels: '.format(len(itrs)), *itrs)
        returnits = True
        itrsislist = True
    elif itrs is np.inf:
        # last iteration
        itrs = scanforfiles(fnamearg)
        if len(itrs):
            if verbose: warning(
                'Found {0} time levels, reading'.format(len(itrs)), itrs[-1])
        else:
            if verbose: warning('Found 0 time levels for {}'.format(fnamearg))
        itrs = itrs[-1:]
        returnits = True
        itrsislist = False
    else:
        returnits = False
        itrsislist = np.iterable(itrs)

    # always make itrs a list
    itrs = aslist(itrs)

    allrec = rec is None
    reclist = aslist(rec)
    if not isinstance(lev, tuple):
        lev = (lev, )
    levs = tuple( aslist(l) for l in lev )
    levdims = tuple(len(l) for l in levs)
    levinds = np.ix_(*levs)
    nlev = len(levdims)

    if usememmap:
        recsatonce = True
        readdata = np.memmap
    else:
        recsatonce = allrec
        readdata = fromfileshape

    try:
        typepre = _typeprefixes[machineformat]
    except KeyError:
        raise ValueError('Allowed machineformats: ' + ' '.join(_typeprefixes))

    arr = None
    metaref = {}
    timeStepNumbers = []
    timeIntervals = []
    for iit,it in enumerate(itrs):
        if additrs:
            fname = fnamearg + '.{0:010d}'.format(int(it))
        else:
            fname = fnamearg

        metafiles = glob.glob(
            fname + 2*('.'+3*'[0-9]') + '.meta') or glob.glob(fname+'.meta')
        if len(metafiles) == 0:
            raise IOError('No files found for ' + fname + '.meta')

        if verbose: warning(metafiles[0])
        if debug: warning('Found', len(metafiles), 'metafiles for iteration', it)

        for metafile in metafiles:
            gdims, i0s, ies, timestep, timeinterval, map2gl, meta = (
                readmeta(metafile))
            if arr is None:
                # initialize, allocate
                try:
                    dataprec, = meta['dataprec']
                except KeyError:
                    dataprec, = meta['format']
                tp = typepre + _typesuffixes[dataprec]
                size = np.dtype(tp).itemsize
                if astype is None: astype = tp
                recshape = tuple( ie-i0 for i0,ie in zip(i0s,ies) )
                count = reduce(mul, recshape)
                nrecords, = meta['nrecords']
                tileshape = (nrecords,) + recshape
                if allrec:
                    reclist = range(nrecords)
                    recinds = np.s_[:,] + levinds
                else:
                    recinds = np.ix_(reclist, *levs)

                if region is None:
                    ri0,rie,rj0,rje = 0,gdims[-1],0,gdims[-2]
                else:
                    ri0,rie,rj0,rje = region
                    if ri0 < 0: ri0 += gdims[-1]
                    if rie < 0: rie += gdims[-1]
                    if rj0 < 0: rj0 += gdims[-2]
                    if rje < 0: rje += gdims[-2]

                assert nlev+2 <= len(gdims)
                rdims = levdims + gdims[len(levdims):-2] + (rje-rj0,rie-ri0)
                # always include itrs and rec dimensions and squeeze later
                arr = np.empty((len(itrs),len(reclist))+rdims, astype)
                arr[...] = fill_value
                metaref = meta
            else:
                if meta != metaref:
                    raise ValueError('Meta files not compatible')

            datafile = metafile[:-4] + 'data'

            if region is not None:
                if map2gl is None:
                    # overlap of tile with region:
                    i0 = min(rie, max(ri0, i0s[-1]))
                    ie = min(rie, max(ri0, ies[-1]))
                    j0 = min(rje, max(rj0, i0s[-2]))
                    je = min(rje, max(rj0, ies[-2]))
                    # source indices
                    I0 = i0 - i0s[-1]
                    Ie = ie - i0s[-1]
                    J0 = j0 - i0s[-2]
                    Je = je - i0s[-2]
                    # target indices
                    i0s[-1] = i0 - ri0
                    ies[-1] = ie - ri0
                    i0s[-2] = j0 - rj0
                    ies[-2] = je - rj0
                else:
                    raise NotImplementedError(
                        "Region selection is not implemented for "
                      + "map2glob != [0,1]")

            sl = tuple( slice(i0,ie) for i0,ie in zip(i0s,ies) )
            if map2gl is None:
                # part of arr that will receive tile (all records)
                arrtile = arr[(iit,slice(None))+sl]
            else:
                ny,nx = arr.shape[-2:]
                i0 = i0s[-1]
                j0 = i0s[-2]
                ie = ies[-1]
                je = ies[-2]
                # "flat" stride for j
                jstride = map2gl[1]*nx + map2gl[0]
                n = (je-j0)*jstride
                # start of a jstride by je-j0 block that contains this tile
                ii0 = min(i0+nx*j0, nx*ny-n)
                # tile starts at ioff+i0
                ioff = nx*j0 - ii0
                # flatten x,y dimensions
                arrflat = arr.reshape(arr.shape[:-2]+(nx*ny,))
                # extract tile
                arrmap = arrflat[...,ii0:ii0+n].reshape(
                    arr.shape[:-2]+(je-j0,jstride))[...,:,ioff+i0:ioff+ie]
                # slice non-x,y dimensions (except records)
                arrtile = arrmap[(iit,slice(None))+sl[:-2]]
                del arrflat,arrmap

            if recsatonce:
                if region is None:
                    arrtile[...] = readdata(datafile, tp, 
                        shape=tileshape)[recinds]
                else:
                    if Ie > I0 and Je > J0:
                        if debug: message(datafile, I0,Ie,J0,Je)
                        arrtile[...] = readdata(datafile, tp, 
                            shape=tileshape)[recinds + np.s_[..., J0:Je, I0:Ie]]
            else:
                f = open(datafile)
                for irec,recnum in enumerate(reclist):
                    if recnum < 0: recnum += nrecords
                    f.seek(recnum*count*size)
                    if region is None:
                        arrtile[irec] = np.fromfile(f, tp, count=count).reshape(
                            recshape)[levinds]
                    else:
                        if Ie > I0 and Je > J0:
                            if debug: message(datafile, I0,Ie,J0,Je)
                            tilerec = np.fromfile(f, tp, count=count).reshape(
                                recshape)
                            arrtile[irec] = (
                                tilerec[levinds + np.s_[..., J0:Je, I0:Ie]])
                f.close()

        if timestep is not None:
            timeStepNumbers.extend(timestep)

        if timeinterval is not None:
            timeIntervals.append(timeinterval)

    # put list of iteration numbers back into metadata dictionary
    if len(timeStepNumbers):
        metaref['timeStepNumber'] = timeStepNumbers

    if len(timeIntervals):
        metaref['timeInterval'] = timeIntervals

    if arr is None:
        arr = np.array([])
    else:
        # Squeeze singleton iteration, record and level dimensions 
        # like matlab version
        dims = (len(itrs),len(reclist)) + levdims
        if squeeze:
            # squeeze all singleton dimensions
            squeezed = tuple( d for d in dims if d > 1 )
        else:
            # squeeze all that came from scalar arguments
            keepers = ([itrsislist, np.iterable(rec)] 
                + [np.iterable(l) for l in lev])
            squeezed = tuple(d for d, keep in zip(dims, keepers) if keep)

        arr = arr.reshape(squeezed+arr.shape[2+nlev:])

    if returnmeta:
        meta = dict((k.lower(), v) for k, v in metaref.items())
        return arr, itrs, meta
    elif returnits:
        return arr, itrs
    else:
        return arr



def parsemeta(metafile):
    """ Parses metafile (file object or filename) into a dictionary of lists
        of floats, ints or strings
    """
    global _currentline

    try:
        lines = open(metafile)
    except TypeError:
        lines = iter(metafile)

    d = {}
    for line in lines:
        line = strip_comments(line)
        # skip empty lines
        if re.match(r'\s*$', line):
            continue

        m = re.match(r' *(\w*) *= *(.*?) *$', line)
        if m:
            key,line = m.groups()
        else:
            raise ParseError(metafile,line)

        # look for the opening delimiter ('[' or '{')
        opening = line[0]
        try:
            closing = _closing[opening]
        except KeyError:
            raise ParseError(metafile, line, 
                "Values must be enclosed in [] or {}.")

        # read more lines until a matching closing delimiter is found
        while closing not in line:
            try:
                nextline = next(lines)
            except StopIteration:
                raise ParseError(metafile,line,'No closing ' + closing 
                    + ' found.')

            line += ' ' + strip_comments(nextline).rstrip()

        if line[-2:] != closing + ';':
            raise ParseError(metafile,line,
                             'Values must be enclosed in "[ ];" or "{ };".')

        # remove delimiters
        line = line[1:-2].strip(" ,")
        _currentline = line

        if opening == '[':
            # [] can contain any type of values, separated by commas
            val = [ parse1(s) for s in re.split(r',? *',line) ]
        else:
            # {} can only contain single quote-delimited strings 
            # separated by space
            val = [ s.rstrip() for s in re.split(r"'  *'", line.strip("'")) ]

        d[key] = val

    return d



def readmeta(f):
    """ Read meta file and extract tile/timestep-specific parameters. """

    meta = parsemeta(f)
    dimList = meta.pop('dimList')

    gdims = tuple(dimList[-3::-3])
    i0s   = [ i-1 for i in dimList[-2::-3] ]
    ies   = dimList[-1::-3]

    # Remove file-specific parameters
    timeInterval = meta.pop('timeInterval', None)
    timeStepNumber = meta.pop('timeStepNumber', None)
    map2gl = meta.pop('map2glob', None)

    # Put back only global dimensions
    meta['dimList'] = list(gdims[::-1])
    return gdims, i0s, ies, timeStepNumber, timeInterval, map2gl, meta
