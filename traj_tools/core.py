"""
core module
-----------

This module contains the heart of the package: the TrajPack class.

"""

import cPickle
import pickle
import numpy as np
import netCDF4 as nc
import glob
import plots
import utils
import fortran.futils as futils


class TrajPack(object):
    """
    Class containing relevant information for 
    a set (or pack) of trajectories.
    
    This class is the heart of this module. 
    It links the raw data and the evaluation and plotting
    functions.
    
    All diagnostics are also performed in a class object.
    
    Inside the class, arrays are numpy.arrays to make masking 
    easier. The output type is as specified under attributes.
    
    
    Parameters
    ----------
    datadir : string
      Path to directory containing output of both, COSMO and trajectory files.
      For now, all files have to be in the same directory.
    xlim : tuple
      Domain dimensions in x-direction in non-rotated coordinates 
      Has to be specified for xy-plotting!
    ylim : tuple
      Domain dimensions in y-direction in non-rotated coordinates. 
      Has to be specified for xy-plotting!
    pollon : float
      Rotated pole longitude. Default value does not change input coordinates.
    pollat : float
      Rotated pole latitude. Default value does not change input coordinates.
      
    
    Attributes
    ----------
    CaseSpecs : Class object
      Class object containing constant specifications of case
    inddict : dictionary
      Dictionary connecting variable name to index in data list
    data : list
      List of numpy.arrays with information about each trajectory

    """
    
    def __init__(self, datadir, xlim, ylim, 
                 pollon = 180., pollat = 90.):

        self.datadir = datadir
        self._init_prop()
        self.pollon = pollon   # Default: 180
        self.pollat = pollat   # Default: 90
        self.xlim = self._nrot2rot(xlim, 'lon')
        self.ylim = self._nrot2rot(ylim, 'lat')
        
        
        
    def _init_prop(self):
        """
        Initialize dictionary and list of numpy arrays with fixed properties:
        
        * filename: location of file containing respective trajectory
        * trjind: index of trajectory in file
        * startt: starting time of trajectory after model start
        
        """
        
        # Extracting data from output directory
        rsuff = "_5m"
        psuff = "p_5m"
        rfiles = glob.glob(self.datadir + "*" + rsuff)
        rfiles = sorted(rfiles)
        pfiles = glob.glob(self.datadir + "*" + psuff)
        self.pfiles = sorted(pfiles)
        self.rfiles = [x for x in rfiles if x not in pfiles]
        self.cfile = glob.glob(self.datadir + '*00c*')[0]
        assert len(glob.glob(self.datadir + '*00c*')) == 1, \
               'More than one cfile detected'
        trjfiles = glob.glob(self.datadir + 'traj_*')
        self.trjfiles = sorted(trjfiles)
        

        # Setting up lists and dictionaries
        self.inddict = dict(startt = 0)
        self.filename = []
        self.trajid = []
        self.data = [[]]
        self.filtdict = dict()
        self.filtlist = []
        
        
        # Looping over all files, initializing lists instead of np.arrays
        # in order to dynamically add files
        
        ntot = 0   # Counter for assert check below
        
        for i in range(len(self.trjfiles)):
            rootgrp = nc.Dataset(self.trjfiles[i], 'r')
            ntmp = len(rootgrp.dimensions['id'])
            tmpt = self.trjfiles[i].split('/')[-1]
            tmpt = int(tmpt.split('_')[1].lstrip('t'))
            self.filename.extend([self.trjfiles[i]] * ntmp)
            self.trajid.extend(range(ntmp))
            self.data[0].extend([tmpt] * ntmp)
            ntot += ntmp
            
        # Convert lists to np.array
        self.ntrj = ntot
        self.filename = np.array(self.filename)
        self.trajid = np.array(self.trajid)
        self.data[0] = np.array(self.data[0])

        assert (self.data[0].shape[0] == self.filename.shape[0] == 
                self.trajid.shape[0] == ntot), \
                "Error while initializing properties for class: Attribute arrays do not have same shape!"
                               
  
    def new_prop_array(self, name, array):
        """
        Add one or more arrays to data. 
        If more than one array is added, names and data arrays must
        be in tuples in the right order.
        
        Parameters
        ----------
        name : string or tuple
          name(s) of the array(s) to be added
        array : np.array of tuple
          data array(s) to be added         
        
        """
        
        if type(name) == str:
            assert (array.shape[0] == self.ntrj), \
                    'Array shaped do not match.'
            assert (name not in self.inddict), 'Name already exists.' 
            self.inddict[name] = len(self.data)
            self.data.append(array)
            print(name, 'has been added.')
            
            
        elif type(name) == tuple or type(name) == list:
            for i in range(len(name)):
                assert (name[i] not in self.inddict), \
                        'One of the names already exists.'
                assert (array[i].shape[0] == ntrj), \
                        'One or all of the array shapes do not match.'
            self.inddict.update(dict(zip(name), 
                                     range(len(self.data), 
                                           len(self.data) + len(name))))
            self.data.extend(array)
            print(name, 'have been added.')
        
    
    
    def new_prop_asc(self, yspan, tracer = 'P'):
        """
        Adds a new ascent filter to data.
        
        Parameters
        ----------
        yspan : float
          Ascent criterion in y-direction
        tracer : str (default = 'P')
          COSMO name of y-axis variable
        
        """
        
        # Update dictionary
        code = tracer + str(yspan)
        self.inddict[code] = len(self.data)
        self.inddict[code + '_start'] = len(self.data) + 1
        self.inddict[code + '_stop'] = len(self.data) + 2
        
        ascdata = minasct(self.trjfiles, yspan, tracer)
        assert (ascdata[0].shape[0] == self.ntrj), \
                'Array shapes do not match. Look for error in source code.'
        self.data.extend(ascdata)
    
       
    def create_filter(self, name, filters):
        """
        Adds a filter to filtlist. Adds name to filtdict.
        
        Parameters
        ----------
        name : string
          Name of given filter
        filters : tuple
          List of filter tuples, see example.
        
        Examples
        --------
        >>> filters = (('P600', 2880, 0), ('P300', 1000, 200))
        
        """
        
        assert (name not in self.filtdict), 'Filter name already exists.'
        
        mask = np.array([True] * self.ntrj)   # Initialize mask 
        
        for i in range(len(filters)):

            ind = self.inddict[filters[i][0]]   # Get index in self.data
            
            if len(filters[i]) == 3:
                mx = filters[i][1]
                mn = filters[i][2]
                
                mask &= self.data[ind] <= mx
                mask &= self.data[ind] >= mn
                
            elif len(filters[i]) == 2:
                crit = filters[i][1]
                
                mask &= self.data[ind] == crit
                
            else:
                raise Exception('Wrong input for filter tuple.')
       
        self.filtdict[name] = len(self.filtlist)
        self.filtlist.append(mask)


    def _nrot2rot(self, nrcoord, mode):
        """
        Converts non-rotated coordinates to rotated coordinates.
        pollon and pollat have to be specified before.
        mode is either 'lon' or 'lat'!
        
        Parameters
        ----------
        nrcoord : float or tuple
          Non-rotated coordinates.
          Either a single value or two values in a tuple
        mode : string
          Must be either 'lat' or 'lon'
        
        Returns
        -------
        rcoord : float or tuple
          Rotated coordinates. In same format as imput.
        
        """
        if mode not in ['lon', 'lat']:
            raise Exception('Invalid input for mode.')
        
        if type(nrcoord) in [float, int]:
            if mode == 'lon':
                rcoord = nrcoord + 180 - self.pollon
            if mode == 'lat':
                rcoord = nrcoord + 90 - self.pollat
                
        elif type(nrcoord) in [tuple, list]:
            if mode == 'lon':
                tmp1 = nrcoord[0] + 180 - self.pollon
                tmp2 = nrcoord[1] + 180 - self.pollon
            if mode == 'lat':
                tmp1 = nrcoord[0] + 90 - self.pollat
                tmp2 = nrcoord[1] + 90 - self.pollat
            rcoord = (tmp1, tmp2)
        else:
            raise Exception('Wrong type for nrcoord.')
        return rcoord
    
    
        
    def _mask_iter(self, filtername):
        """
        Converts mask to iterable lists. To be used in plotting functions, etc.
        
        Parameters
        ----------
        filtername : np.array
          Mask 
        
        Returns
        -------
        uniqueloc : list
          List of all unique save locations.
          Use as first iterator.
        idlist : list
          List of list. Trajectory IDs for every element in 
          uniquelist.yspan
          Use as second iterator.
          
        """
        maskid = self.filtdict[filtername]
        mask = self.filtlist[maskid]
        uniqueloc = np.unique(self.filename[mask])
        idlist = []
        for i in range(len(uniqueloc)):
            locmask = (self.filename == uniqueloc[i])
            idlist.append(list((self.trajid[locmask & mask])))

        return list(uniqueloc), idlist
    
    
    def _mask_array(self, filtername, dataname):
        """
        Returns one masked array specified by given mask and dataname.
        
        Parameters
        ----------
        filtername : np.array
          Mask 
        dataname : string
          Identifier of data array
          
        Returns
        -------
        marray         : np.array
          Masked data array
        
        """
        
        maskid = self.filtdict[filtername]
        mask = self.filtlist[maskid]
        
        dataid = self.inddict[dataname]
        data = self.data[dataid]
        
        return data[mask]
        
    
    def saveme(self, savename):
        """
        Saves class object under savename.
        As pickled file.
        
        Parameters
        ----------
        savename : string
          Save path
          
        """
        
        f = open(savename, 'w+')
        pickle.dump(self, f, 2)   # Using Protocol 2
        f.close()
        
    def calc_theta(self):
        """
        Class method to add Potential Temperature as a Variable to the netCDF files.
        
        """
        
        utils.calc_theta(self.trjfiles)
        
        
    ######################
    # Plotting functions
    ######################
    
    def draw_hist(self, dataname, filtername = None, savebase = None, 
                  starts = False):
        """
        Draws a Histogram of data specified by dataname.
        If filtername is given, plots only filtered trajectories.
        If starts is True, plots one histogram seperately for each start time.
        
        Parameters
        ----------
        dataname : string
          Identifier of wanted data array
        filtername : string
          Identifier of wanted filter
        savebase : string
          Path to output directory
        starts : bool
          If False, plots histogram for entire array given by filter.
          If True, plots histograms for each start time.
        
        """
        if filtername != None:
            array = self._mask_array(filtername, dataname)
        else:
            array = self.data[self.inddict[dataname]]
        if savebase != None:    
            savename = savebase + 'hist_' + dataname + '_' + filtername + '.png'
        else:
            savename = savebase
        plots.draw_hist(array, savename = savename)
        
     
    def draw_xy(self, varlist, filtername = None, savebase = None, 
                starts = False):
        """
        Draws XY Plot of trajectories with color as a function of 'P'.
        If filtername is given, plots only filetered trajectories.
        If starts is True, plots one xy plot seperately for each start time. 
        
        Parameters
        ----------
        varlist : list
          List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
        filtername : string
          Identiefier of wanted filter
        savebase : string
          Path to output directory
        starts : bool
          If False, plots xy plot for entire filtered array (NOT recommended,
          if start time is not used in filter!)
          If True, plots xy plot seperately for each start time.
        
        """
        assert (self.xlim != None and self.ylim != None), \
                'xlim and/or ylim are not specified!'
        
        if savebase != None:    
            savename = savebase + 'xy_' + filtername + '.png'
        else:
            savename = savebase
            
        loclist, idlist = self._mask_iter(filtername)
        
        plots.draw_xy(varlist, loclist, idlist, self.cfile,
                      self.rfiles, self.pfiles, savename = savename, 
                      pollon = self.pollon, pollat = self.pollat, 
                      xlim = self.xlim, ylim = self.ylim)
         

def loadme(savename):
    """
    Unpickles file saved under savename.
    
    Parameters
    ----------
    savename : string
      Save path
          
    Returns
    -------
    obj : object
      Saved object
      
    """
    
    f = open(savename, 'r')
    obj = pickle.load(f)
    f.close()
    return obj 


def minasct(filelist, yspan, tracer):
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
     
    Returns
    -------
    ascdata : tuple
      Tuple containing three np.arrays:
      * Minimum ascent time
      * Index of ascent start
      * Index of ascent stop
    
    """
    
    # Initialize lists, convert to np.array later
    asct = []
    ascstart = []
    ascstop = []
     
    for i in range(len(filelist)):
        mat = nc.Dataset(filelist[i], 'r').variables[tracer][:, :]
        for j in range(mat.shape[1]):
            asctup = minxspan(mat[:, j], yspan)
            asct.append(asctup[0])
            ascstart.append(asctup[1])
            ascstop.append(asctup[2])
    assert (len(asct) == len(ascstart) == len(ascstop)), \
            'Array lengths do not match'
    ascdata = (np.array(asct), np.array(ascstart), np.array(ascstop))
    return ascdata
            

def minxspan(array, yspan, flip = False):
    """
    Returns the minimum time steps needed to conver given criterion.
    Automatically filters out zero values. If yspan is not fulfilled, returns
    np.nan. 
    
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
    
    """
    
    # Filter out zeros, adjust for offset
    offset = np.where(array != 0)[0][0]
    array = array[array != 0]
    
    # Flip array if needed
    if flip:
        array = -array
    
    # Check if criterion is met
    if np.amax(array) - np.amin(array) < crit:
        xspan = np.nan
        istart = np.nan
        istop = np.nan
    else:
        # Use Fortran implementation, use 0 as error values
        xspan, istart, istop  = futils.futils.minxspan(array, crit, 
                                                        len(array) + 1, 0, 0)

        assert ((xspan > 0) and (xspan <= len(array)) and (istart >= 0) 
                and (istop >= 0)), \
                'One of the minxspan outputs is zero or negative.'
            
    return (xspan, istart + offset, istop + offset)


if __name__ == '__main__':
    """
    Test routine for minxspan algorithm.
    
    """
    
    print('************************************')
    print('Running algorithm check')
    
    alen = 1000
    crit = 200
    import timeit
    import matplotlib.pyplot as plt
    
    def _wrapper(func, *args, **kwargs):
        """
        Function wrapper for use in timeit module
        """
        def _wrapped():
            return func(*args, **kwargs)
        return _wrapped
    

    print('********* Test array 1 *************')
    a1 = np.array([x for x in range(alen)])
    wrap = _wrapper(minxspan, a1, crit)
    print('Time taken for array 1:', timeit.timeit(wrap, number = 1000))
    print('Results for array 1:', minxspan(a1, crit))
    plt.plot(a1)
    span, a, b = minxspan(a1, crit)
    plt.scatter([a, b], [a1[a], a1[b]])
    
    print('********* Test array 2 *************')
    a2 = np.array([(x**2 / 1000) for x in range(alen)])
    wrap = _wrapper(minxspan, a2, crit)
    print('Time taken for array 2:', timeit.timeit(wrap, number = 1000))
    print('Results for array 2:', minxspan(a2, crit))
    plt.plot(a2)
    span, a, b = minxspan(a2, crit)
    plt.scatter([a, b], [a2[a], a2[b]])
    
    print('********* Test array 3 *************')
    a3 = np.array([(500 *(np.sin(0.1 * (x - 500.1))) / 
                    (0.1 * (x - 500.1)) + 300) for x in range(alen)])
    wrap = _wrapper(minxspan, a3, crit)
    print('Time taken for array 3:', timeit.timeit(wrap, number = 1000))
    print('Results for array 3:', minxspan(a3, crit))
    plt.plot(a3)
    span, a, b = minxspan(a3, crit)
    plt.scatter([a, b], [a3[a], a3[b]])
    
    print('********* Test array 4 *************')
    a4 = np.array([500 for x in range(alen)])
    wrap = _wrapper(minxspan, a4, crit)
    print('Time taken for array 4:', timeit.timeit(wrap, number = 1000))
    print('Results for array 4:', minxspan(a4, crit))
       
    print('********* Test array 5 *************')
    a5 = np.array([(500 *(np.sin(0.1 * (x - 500.1))) / 
                    (0.1 * (x - 500.1)) - 100) for x in range(alen)])
    wrap = _wrapper(minxspan, a5, crit, flip = True)
    print('Time taken for array 5:', timeit.timeit(wrap, number = 1000))
    print('Results for array 5:', minxspan(a5, crit, flip = True))
    plt.plot(a5)
    span, a, b = minxspan(a5, crit, flip = True)
    plt.scatter([a, b], [a5[a], a5[b]])
    
    print('********* Test array 6 *************')
    a6 = np.array([(500 *(np.sin(0.05 * x)) - 500) for x in range(alen)])
    a6[:300] = 0
    a6[700:] = 0
    wrap = _wrapper(minxspan, a6, crit)
    print('Time taken for array 6:', timeit.timeit(wrap, number = 1000))
    print('Results for array 6:', minxspan(a6, crit))
    plt.plot(a6)
    span, a, b = minxspan(a6, crit)
    plt.scatter([a, b], [a6[a], a6[b]])
    
    plt.show()
    
    
    
    
    
    
    

