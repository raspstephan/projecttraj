"""
utils module
------------

Contains helpful functions outside of core functionality.

"""

import glob
import numpy as np
import netCDF4 as nc



####################
### Helper Functions
####################


def create_startfile(lonmin, lonmax, dlon, 
                     latmin, latmax, dlat, 
                     zmin, zmax, dz, outdir):
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
      
    """

    suff = 'test.txt'  # Later: get name from CaseSpecs
    
    f = open(outdir + suff, 'w+')
    
    f.write('Reference Date somedate\n')  # Write first line, slo use CaseSpecs
    f.write('lon lat z\n')
    f.write('-----------------\n')
    
    # Create lon and lat arrays
    
    lonlist = list(np.arange(lonmin, lonmax, dlon))
    latlist = list(np.arange(latmin, latmax, dlat))
    zlist = list(np.arange(zmin, zmax, dz))
    
    print('Total number of trajectories:', 
          len(lonlist) * len(latlist) * len(zlist))
    
    for j in range(len(lonlist)):
        for n in range(len(latlist)):
            for m in range(len(zlist)):
                lontmp = '%3.3f' % (lonlist[j])
                lontmp = (8 - len(lontmp)) * ' ' + lontmp
                lattmp = '%2.3f' % (latlist[n])
                lattmp = (7 - len(lattmp)) * ' ' + lattmp
                ztmp = '%4.3f' % (zlist[m])
                ztmp = (8 - len(ztmp)) * ' ' + ztmp
                line = (' ' + lontmp + ' ' + lattmp + ' ' + ztmp)
                assert (len(line) == 26), \
                        'Start file line does not have correct length'
                f.write(line + '\n')
    
    f.close() 
                    
                

def calc_totlh(files):
    """
    Adds Total Latent Heat as a Variable to given netCDF files.
    
    Parameters
    ----------
    
    files : List or string
      Either a list of all netCDF files to be considered, 
      or the path to the directory containing .nc files. 
    
    """
    
    # Checking if filelist needs to be created
    if type(files) == list:
        filelist = files
    elif type(files) == str:
        filelist = sorted(glob.glob(files + '*.nc'))
    else:
        raise Exception('Wrong type for files.')
    
    # Iterate over files in filelist
    for f in filelist:
        rootgrp = nc.Dataset(f, 'r')
        
        # Read file arrays needed for calculation: QV, QI, etc. eg:
        m1 = rootgrp.variables['QV'][:, :]
        # ...
        
        # Do the calculation
        # Should return [:, :] array
        
        # Add new array to netCDF file
        totlh = rootgrp.createVariable('TOT_LH', 'f4', ('time', 'id'))
        totlh = outarray
        
        rootgrp.close()
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

