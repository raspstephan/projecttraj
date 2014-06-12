"""
utils module
------------

Contains helpful functions outside of core functionality.

"""

import glob
import os.path
import numpy as np
import netCDF4 as nc
import cosmo_utils.pywgrib as pwg



####################
### Helper Functions
####################


def create_startfile(lonmin, lonmax, dlon, 
                     latmin, latmax, dlat, 
                     zmin, zmax, dz, outdir,
                     zfilter = False, cfile = None):
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
    outdir : str
      Path to output directory
    zfilter : bool
      If True, trajectories which are "under ground" are filtered out. 
      cfile has to be specified! 
    cfile : str
      Location of COSMO EU constants file
    
    """
    
    suff = ''
    
    f = open(outdir + suff, 'w+')
    
    f.write('Reference Date somedate\n')  # Write first line, slo use CaseSpecs
    f.write('lon lat z\n')
    f.write('-----------------\n')
    
    # Create lon and lat arrays
    
    lonlist = list(np.arange(lonmin, lonmax, dlon))
    latlist = list(np.arange(latmin, latmax, dlat))
    zlist = list(np.arange(zmin, zmax, dz))
    
    if zfilter:
        hh = pwg.getfieldobj(cfile, 'HH')
    
    print('Total number of trajectories:', 
          len(lonlist) * len(latlist) * len(zlist))
    
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
    
    # Iterate over files in filelist
    for f in filelist:
        print('Open file:', f)
        rootgrp = nc.Dataset(f, 'a')
        
        # Read file arrays needed for calculation
        pmat = rootgrp.variables['P'][:, :] * 100.   # Convert to SI Units
        tmat = rootgrp.variables['T'][:, :]

        # Add new array to netCDF file
        theta = rootgrp.createVariable('THETA', 'f4', ('time', 'id'))
        theta[:, :] = tmat * ((P0 / pmat) ** (R / CP))
        rootgrp.close()
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

