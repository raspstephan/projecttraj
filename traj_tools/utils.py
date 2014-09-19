"""
utils module
------------

Contains helpful functions outside of core functionality.

"""

import glob
import os
import numpy as np
import netCDF4 as nc
import cosmo_utils.pywgrib as pwg
import fortran.futils as futils
import cPickle
import scipy.ndimage as ndi



####################
### Helper Functions
####################


def create_startfile(lonmin, lonmax, dlon, 
                     latmin, latmax, dlat, 
                     zmin, zmax, dz, outname,
                     zfilter = True, cfile = None,
                     refdate = None):
    """
    Creates a trajectory start file in output directory
    
    Parameters
    ----------

    lonmin : float
      Lower boundary of longitude
    lonmax : float
      Upper boundary of longitude
    dlon : float
      longitude increment
    latmin : float
      Lower boundary of latitude
    latmax : float
      Upper boundary of latitude
    dlat : float
      latitude increment
    zmin : float
      Lower height boundary [m]
    zmax : float
      Upper height boundary [m]
    dz : float
      Height increment [m]
    outname : str
      Path to output directory, plus prefix
    zfilter : bool
      If True, trajectories which are "under ground" are filtered out. 
      cfile has to be specified! 
    cfile : str
      Location of COSMO EU constants file
    refdate : str
      Reference date to be written as a header
      
    Examples
    --------
    
    >>> import traj_tools as trj
    >>> cfile = '/home/scratch/users/stephan.rasp/Case1_20070720/d2euout/lfff00000000c'
    >>> trj.utils.create_startfile(-5, 5, 0.25, -7, 7, 0.25, 1000, 5000, 250, 
    >>>                            './case1_coarse.sf', zfilter = True,
    >>>                            cfile, '20070720')
    
    """
    
    suff = '.sf'
    
    f = open(outname + suff, 'w+')
    if refdate == None:
        f.write('Reference Date N/A\n')
    else:
        f.write('Reference Date ' + refdate + '\n')
    f.write('lon lat z\n')
    f.write('-----------------\n')
    
    # Create lon and lat arrays
    lonlist = list(np.arange(lonmin, lonmax + dlon, dlon))
    latlist = list(np.arange(latmin, latmax + dlat, dlat))
    zlist = list(np.arange(zmin, zmax + dz, dz))
    
    if zfilter:
        if cfile == None:
            raise Exception('No cfile give for filter')
        hh = pwg.getfieldobj(cfile, 'HH')
        
    print 'Output file is:', outname + suff
    print 'Unfiltered number of trajectories:', \
          len(lonlist) * len(latlist) * len(zlist)
    
    trjcount = 0
    
    for lon in lonlist:
        for lat in latlist:
            if zfilter:
                # Interpolate to closest x and y indices
                xind = (np.abs(hh.lons - lon)).argmin()
                yind = (np.abs(hh.lats - lat)).argmin()
                
                # Get nearest z value
                hhint = hh.data[-1, yind, xind]
                
                # Filter out values in zlist which are below hhint
                zlisttmp = list(np.array(zlist)[np.array(zlist) > hhint])
                
            else:
                zlisttmp = zlist
               
            for z in zlisttmp:
                lontmp = '%3.3f' % (lon)
                lontmp = (8 - len(lontmp)) * ' ' + lontmp
                lattmp = '%2.3f' % (lat)
                lattmp = (7 - len(lattmp)) * ' ' + lattmp
                ztmp = '%4.3f' % (z)
                ztmp = (8 - len(ztmp)) * ' ' + ztmp
                line = (' ' + lontmp + ' ' + lattmp + ' ' + ztmp)
                assert (len(line) == 26), \
                        'Start file line does not have correct length'
                f.write(line + '\n')
                trjcount += 1
    
    print 'Total number of trajectories:', trjcount
    print 'Filtered out:', len(lonlist) * len(latlist) * len(zlist) - trjcount
    
    f.close() 
                    
                

def calc_theta(files):
    """
    Adds Potential Temperature as a Variable to given netCDF files. 
    The new variable is called 'THETA'.
    
    Parameters
    ----------
    
    files : List or string
      Either a list of all netCDF files to be considered, 
      or the path to the directory containing .nc files. 
    
    """
    
    # Define constants
    P0 = 1.e5   # reference pressure [Pa]
    R = 287.    # specific gas constant dry air [J K-1 kg-1]
    CP = 1004.  # specific heat at constant pressure [J K-1 kg-1]
    
    
    # Checking if filelist needs to be created
    if type(files) == list:
        filelist = files
    elif type(files) == str:
        files = os.path.normpath(files) + '/'
        filelist = sorted(glob.glob(files + '*.nc'))
    else:
        raise Exception('Wrong type for files.')
    
    assert (len(filelist) > 0), 'No files selected.'
    
    # Create new files
    thetalist = []
    for f in filelist:
        thetaf = f.rstrip('.nc') + '_theta.nc'
        thetalist.append(thetaf)
        os.system('cp ' + f + ' ' + thetaf)
    
    # Iterate over files in filelist
    for f in thetalist:
        print 'Open file:', f 
        rootgrp = nc.Dataset(f, 'a')
        
        # Read file arrays needed for calculation
        pmat = rootgrp.variables['P'][:, :] * 100.   # Convert to SI Units
        tmat = rootgrp.variables['T'][:, :]

        # Add new array to netCDF file
        theta = rootgrp.createVariable('THETA', 'f4', ('time', 'id'))
        theta[:, :] = tmat * ((P0 / pmat) ** (R / CP))
        rootgrp.close()
        
        
    return thetalist
    

def convert_pickle2netcdf(indir, outdir):
    """
    Converts pickled trajectory files to NetCDF files, so they can be used with
    this version of traj_tools
    
    Parameters
    ----------
    indir : string
      Path to directory containing GRIB files
    outdir : string
      Path to directory where netcdf files are created
    """
    
    # Specifics for Case0
    pnum = 7   # ptracer
    trjstart = 360 * 60   # Trj start after model start in secs
    dt = 300   # Save time interval in secs
    tsteps = 1369   # Max time steps in pickled files
    
    # Create filelist
    inprefix = '/test_case_'
    filelist = sorted(glob.glob(indir + inprefix + '*'))
    
    # Output prefix
    savebase = outdir + '/traj_t' 
    
    # Set counter
    fcount = 1
    end = False
    
    # loop over all files
    for fn in filelist:
        print 'Converting file:', fn
        
        # Load pickled files and extract data
        f = open(fn, 'r')
        spec = cPickle.load(f)
        mat = cPickle.load(f)
        f.close()
        
        # Check for starting time
        tstart = np.nonzero(mat[0, pnum, :])[0][0]   # Start index of first traj
        
        # Check if mat needs to be split
        if np.where(mat[:, pnum, tstart] == 0)[0].shape[0] != 0:
            splitind = np.where(mat[:, pnum, tstart] == 0)[0][0]
            if np.nonzero(mat[splitind, pnum, :])[0].shape[0] != 0:
                tstart2 = np.nonzero(mat[splitind, pnum, :])[0][0]
            else:
                end = True

            
            # Swap axes
            conmat = np.swapaxes(mat[:splitind, :, tstart:], 0, 2)
            conmat2 = np.swapaxes(mat[splitind:, :, tstart2:], 0, 2)
            
            # First file
            _write2netcdf(conmat, tstart, savebase, fcount)
            fcount = 1
            # Second file
            if end:
                print 'Arrays contain only zeros! END!'
                break
            else:
                _write2netcdf(conmat2, tstart2, savebase, fcount)
            fcount += 1
        
        else:
            # Swap axes
            conmat = np.swapaxes(mat[:, :, tstart:], 0, 2)  
            
            # Just one file
            _write2netcdf(conmat, tstart, savebase, fcount)
            fcount += 1
        
        
def _write2netcdf(conmat, tstart, savebase, fcount):
    """
    To be used inside this function
    
    Parameters
    ----------
    conmat : np.array
        Converted and trimmed matrix
    savename : string
        Path of new netCDF file
    
    """
    # Specifics for Case0, see above (repeated here for convenience)
    trjstart = 360 * 60   # Trj start after model start in secs
    dt = 300   # Save time interval in secs
    tsteps = 1369   # Max time steps in pickled files
    
    # Create savename
    savename = (savebase + str(360 + tstart * 5).zfill(6) + 
                '_p' + str(fcount).zfill(3) + '.nc')
    
    # Create time array
    tarray = np.arange(tstart * dt + trjstart, 
                       tsteps* dt + trjstart, dt)
    
    # Create NetCDF file
    print 'Creating netCDF file:', savename
    rootgrp = nc.Dataset(savename, 'w')
    # Write global attributes
    rootgrp.ref_year = 2012
    rootgrp.ref_month = 10
    rootgrp.ref_day = 13
    rootgrp.ref_hour = 0
    rootgrp.ref_min = 0
    rootgrp.ref_sec = 0
    rootgrp.duration_in_sec = tsteps* dt + trjstart
    rootgrp.pollon = -165
    rootgrp.pollat = 35
    rootgrp.output_timestep_in_sec = 300
    # Create dimensions
    time = rootgrp.createDimension('time', tarray.shape[0])
    id = rootgrp.createDimension('id', conmat.shape[2])
    # Create variables
    times = rootgrp.createVariable('time', 'f4', ('time', ))
    lon = rootgrp.createVariable('longitude', 'f4', ('time', 'id', ))
    lat = rootgrp.createVariable('latitude', 'f4', ('time', 'id', ))
    u = rootgrp.createVariable('u', 'f4', ('time', 'id', ))
    v = rootgrp.createVariable('v', 'f4', ('time', 'id', ))
    w = rootgrp.createVariable('w', 'f4', ('time', 'id', ))
    t = rootgrp.createVariable('T', 'f4', ('time', 'id', ))
    z = rootgrp.createVariable('z', 'f4', ('time', 'id', ))
    p = rootgrp.createVariable('P', 'f4', ('time', 'id', ))
    # Fill variables
    times[:] = tarray
    lon[:, :] = conmat[:, 0, :]
    lat[:, :] = conmat[:, 1, :]
    z[:, :] = conmat[:, 2, :]
    u[:, :] = conmat[:, 3, :]
    v[:, :] = conmat[:, 4, :]
    w[:, :] = conmat[:, 5, :]
    t[:, :] = conmat[:, 6, :]
    p[:, :] = conmat[:, 7, :]
    
    rootgrp.close()
        
    


####################################################
# Functions used in core
####################################################

def _delta(filelist, tracer):
    """
    Calculate difference of given tracer between min und max 'P' values
    
    Parameters
    ----------
    filelist : list
      List of saved NetDCF file locations
    yspan : float
      Ascent criterion in y-direction
    tracer : str 
      COSMO name of y-axis variable
      
    Returns
    -------
    deltaarray : np.array
    
    """
    # Initialize lists, convert to np.array later
    deltaarray = []
    
    for f in filelist:
        print 'Opening file:', f
        pmat = nc.Dataset(f, 'r').variables['P'][:, :]
        trcmat = nc.Dataset(f, 'r').variables[tracer][:, :]
        for j in range(pmat.shape[1]):
            parray = pmat[:, j][np.isfinite(pmat[:, j])]
            trcarray = trcmat[:, j][np.isfinite(trcmat[:, j])]
            minind = parray.argmin()
            maxind = parray.argmax()
            if minind < maxind:
                delta = trcarray[maxind] - trcarray[minind]
            else:
                delta = trcarray[minind] - trcarray[maxind]
            if tracer == 'POT_VORTIC' and delta > 0.00005:
                delta = np.nan
                print 'Detected irregulat PV value, set to NaN!'
            deltaarray.append(delta)

    return np.array(deltaarray)
            
    
    


def _minasct(filelist, yspan, tracer, dtrj, interpolate = False):
    """
    Calculate minimum ascent time for all trajectories from NetCDF files.
    
    Parameters
    ----------
    filelist : list
      List of saved NetDCF file locations
    yspan : float
      Ascent criterion in y-direction
    tracer : str 
      COSMO name of y-axis variable
    dtrj : float
      Timestep between saved values
    interpolate : bool
      If True, ascent time will be interpolated
     
    Returns
    -------
    ascdata : tuple
      Tuple containing three np.arrays:
      * Minimum ascent time
      * Index of ascent start
      * Index of ascent stop
      * Value of tracer at start index
      * Value of tracer at stop index
    
    """
    
    # Initialize lists, convert to np.array later
    asct = []
    ascstart = []
    ascstop = []
    ascstartval = []
    ascstopval = []
    
    if tracer == 'P':
        flip = True
    else:
        flip = False
    
    for f in filelist:
        print 'Opening file:', f
        rootgrp = nc.Dataset(f, 'r')
        mat = rootgrp.variables[tracer][:, :]
        trjstart = int(rootgrp.variables['time'][0] / 60) / dtrj
        for j in range(mat.shape[1]):
            asctup = _minxspan(mat[:, j], yspan, flip)
            asct.append(asctup[0])
            ascstart.append(asctup[1] + trjstart)
            ascstop.append(asctup[2] + trjstart)
            ascstartval.append(asctup[3])
            ascstopval.append(asctup[4])
    assert (len(asct) == len(ascstart) == len(ascstop) == len(ascstartval) == 
            len (ascstopval)), \
            'Array lengths do not match'
    if interpolate:
        if tracer == 'P':
            diff = np.array(ascstartval) - np.array(ascstopval)
        else: 
            diff = np.array(ascstopval) - np.array(ascstartval)
        asct = (np.array(asct) * yspan) / diff
        
    ascdata = (np.array(asct) * dtrj, np.array(ascstart) * dtrj, 
               np.array(ascstop) * dtrj, np.array(ascstartval), 
               np.array(ascstopval))
    return ascdata


def _allasct(filelist, yspan, xmax, tracer, dtrj):
    """
    Finds and returns info about all areas in all given trajectories where 
    yspan is covered in xmax.
    
    Parameters
    ----------
    filelist : list
      List of saved NetDCF file locations
    yspan : float
      Ascent criterion in y-direction
    xmax : float
      Maximum x-span for yspan
    tracer : str 
      COSMO name of y-axis variable
    dtrj : float
      Timestep between saved values
    
     
    Returns
    -------
    alllist : np.array with type object
      Contains following information for all seperate areas for all trjs
      [file index in trjlist, index in file, lonstart, lonstop, latstart, 
      latstop, startt, stopt in mins after modelstart, startval, stopval]
    """
    
    # Initialize list
    alllist = []
    
    if tracer == 'P':
        flip = True
    else:
        flip = False
    fi = 0   # File Index
    for f in filelist:
        print 'Opening file:', f
        rootgrp = nc.Dataset(f, 'r')
        mat = rootgrp.variables[tracer][:, :]
        lon = rootgrp.variables['longitude'][:, :]
        lat = rootgrp.variables['latitude'][:, :]
        trjstart = int(rootgrp.variables['time'][0] / 60)
        for j in range(mat.shape[1]):
            alllist.append([])
            tuplist = _allxspan(mat[:, j], yspan, xmax, flip)
            for tup in tuplist:
                xstart = lon[tup[0], j]
                xstop = lon[tup[1] -1, j]
                ystart = lat[tup[0], j]
                ystop = lat[tup[1] - 1, j]
                alllist[-1].append( (fi, j, xstart, xstop, ystart, ystop) + 
                                  (tup[0] * dtrj + trjstart, 
                                   tup[1] * dtrj + trjstart, tup[2], tup[3]) )
        fi +=1
        
    return np.array(alllist)


def _allasct_cd(filelist, diff, sigma, tracer, dtrj):
    """
    TODO
    """
    
    # Initialize list
    alllist = []
    
    if tracer == 'P':
        flip = True
    else:
        flip = False
    fi = 0
    
    for f in filelist:
        print 'Opening file:', f
        rootgrp = nc.Dataset(f, 'r')
        mat = rootgrp.variables[tracer][:, :]
        pmat = rootgrp.variables['P'][:, :]
        lon = rootgrp.variables['longitude'][:, :]
        lat = rootgrp.variables['latitude'][:, :]            
        trjstart = int(rootgrp.variables['time'][0] / 60)
        
        for j in range(mat.shape[1]):
            alllist.append([])
            
            # Smooth, filter and get CD
            array = mat[:, j]
            array = array[np.isfinite(array)]
            array = array[array != 0]
            array = np.gradient(array)
            array = ndi.filters.gaussian_filter(array, sigma)

            # Get slices 
            mask = np.ma.masked_where(array > (diff * 5), array)
            slices = np.ma.notmasked_contiguous(mask)
            
            # Create lists
            tuplist = []
            
            # Loop
            if slices == None:
                return ()
            elif type(slices) != list:
                istart = slices.start
                istop = slices.stop
                tuplist.append( (istart, istop, pmat[istart, j], 
                                pmat[istop - 1, j]) )
            else:
                for s in slices:
                    istart = s.start
                    istop = s.stop
                    tuplist.append( (istart, istop, pmat[istart, j], 
                                    pmat[istop - 1, j]) )
            for tup in tuplist:
                xstart = lon[tup[0], j]
                xstop = lon[tup[1] -1, j]
                ystart = lat[tup[0], j]
                ystop = lat[tup[1] - 1, j]
                alllist[-1].append( (fi, j, xstart, xstop, ystart, ystop) + 
                                (tup[0] * dtrj + trjstart, 
                                tup[1] * dtrj + trjstart, tup[2], tup[3]) )
        fi +=1
        
    return np.array(alllist)
           

def _minxspan(array, yspan, flip = False):
    """
    Returns the minimum time steps needed to conver given criterion.
    Automatically filters out zero and nan values. If yspan is not fulfilled, 
    returns np.nan. 
    
    Parameters
    ----------
    array : np.array
      Input array
    yspan : float
      Ascent criterion in y-direction
    flip : bool
      If True array will be flipped (e.g. use for 'P')
    
    Returns
    -------
    xspan : int
      Time steps for criterion
    istart : int
      Start index of ascent
    istop : int
      End index of ascent
    startval : float 
      Value of tracer at start index
    stopval : float
      Value of tracer at stop index
    """
    
    # Filter out nans and zeros
    array = array[np.isfinite(array)]
    array = array[array != 0]
    # Flip array if needed
    if flip:
        array = -array 

    # Check if criterion is met
    #print array
    if array.shape[0] == 0:
        xspan = np.nan
        istart = np.nan
        istop = np.nan
        startval = np.nan
        stopval = np.nan
    elif np.amax(array) - np.amin(array) < yspan:
        xspan = np.nan
        istart = np.nan
        istop = np.nan
        startval = np.nan
        stopval = np.nan
    else:
        # Use Fortran implementation, use 0 as error values
        xspan, istart, istop, startval, stopval  = futils.futils.minxspan(array, 
                                                   yspan, len(array) + 1, 0, 0, 0, 0)
        
        # Check if output is correct. NOTE: THIS IS A POTENTIAL BUG!!!
        if (istart < 0) and (istop < 0):
            xspan = np.nan
            istart = np.nan
            istop = np.nan
            startval = np.nan
            stopval = np.nan
            
        #assert ((xspan > 0) and (xspan <= len(array)) and (istart >= 0) 
                #and (istop >= 0)), \
                #'One of the minxspan outputs is zero or negative.'
    if flip:
        startval = -startval
        stopval = -stopval
    return [xspan, istart, istop, startval, stopval]
 
 
 
def _allxspan(array, yspan, xmax, flip = False):
    """
    Finds all areas in array where yspan is covered in xmax and returns 
    information.
    Automatically filters out zero and nan values. If yspan is not fulfilled, 
    returns empty tuple. 
    
    Parameters
    ----------
    array : np.array
      Input array
    yspan : float
      Ascent criterion in y-direction
    xmax : float
      Maximum x-span for yspan
    flip : bool
      If True array will be flipped (e.g. use for 'P')
    
    Returns
    -------
    tuplist : list of tuples
      If criterion is not met for array, returns empty tuple.
      Otherwise returns list of tuples, e.g. [(...), (...)],
      where each tuple contains the following information:
      (startindex, stopindex, startvalue, stopvalue)
      
    """
    
    # Filter out nans and zeros
    array = array[np.isfinite(array)]
    array = array[array != 0]
    
    # Flip array if needed
    if flip:
        array = -array 

    # Check if criterion is met
    if array.shape[0] == 0:
        xspan = np.nan

    elif np.amax(array) - np.amin(array) < yspan:
        xspan = np.array([len(array) + 1] * array.shape[0])
    else:
        # Use Fortran implementation, use 0 arrays as error values
        
        xspan, startval, stopval  = futils.futils.allxspan(
            array, yspan, np.array([len(array) + 1] * array.shape[0]),
            np.zeros(array.shape[0]), np.zeros(array.shape[0]))

        # Check if output is correct. NOTE: THIS IS A POTENTIAL BUG!!!
        #if (istart < 0) and (istop < 0):
            #xspan = np.nan
            #startval = np.nan
            #stopval = np.nan
            
        #assert ((xspan > 0) and (xspan <= len(array)) and (istart >= 0) 
                #and (istop >= 0)), \
                #'One of the minxspan outputs is zero or negative.'
    if flip:
        array = -array   # Flip array back

    
    # Get slices
    mask = np.ma.masked_where(xspan > xmax, xspan)
    slices = np.ma.notmasked_contiguous(mask)
    
    # Create lists
    tuplist = []
    
    # Loop
    if slices == None:
        return ()
    else:
        for s in slices:
            istart = s.start
            istop = s.stop
            tuplist.append( (istart, istop, array[istart], 
                            array[istop - 1 + xspan[istop -1]]) )
        return tuplist   
    
    
# Testing
if __name__ == '__main__':   
    #from random import randrange
    #fn = '/home/scratch/users/stephan.rasp/Case1_20070720/d4deout/traj_t001800_p001.nc'
    #rootgrp = nc.Dataset(fn, 'r')
    #array = rootgrp.variables['P'][:, 182]
    #print _allxspan(array, 100, 10, True)
    
    fl = ['/home/scratch/users/stephan.rasp/Case1_20070720/d4deout/traj_t001800_p001.nc',
          '/home/scratch/users/stephan.rasp/Case1_20070720/d4deout/traj_t001800_p002.nc']
    tout = _allasct(fl, 100, 10, 'P', 5)
    
    
    
    
    
    
    
    

