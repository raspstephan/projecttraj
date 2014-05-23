"""
filters module!!
Submodule of traj_tools
Opening and Filtering saved trajectory files
Naming conventions:
Files = all trajectories
Matrix = one trajectory file (i.e. 1000 trajectories)
Array = one trajectory
"""
import common
import cPickle
import pickle
import numpy as np
import loadbin
import netCDF4 as nc
import glob
from . import plots


class TrajProp(object):
    """
    Class containing relevant information for 
    a set of trajectories.
    
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
      
    
    Attributes
    ----------
    CaseSpecs : Class object
      Class object containing constant specifications of case
    inddict : dictionary
      Dictionary connecting variable name to index in data list
    data : list
      List of numpy.arrays with information about each trajectory

    """
    
    def __init__(self,
                 datadir):

        self.datadir = datadir
        self._init_prop()
        
        
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
        
        
    def _mask_iter(self, filtername):
        """
        Converts mask to iterable lists. To be used in plotting functions, etc.
        
        Paremeters
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
        marray : np.array
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
        if savebase != None:    
            savename = savebase + 'xy_' + filtername + '.png'
        else:
            savename = savebase
        
        plots.draw_xy(varlist, self._mask_iter(filtername), self.cfile,
                      self.rfiles, self.pfiles, savename)
         


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


kate-swp

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
            

def minxspan(Array, Criterion, mode = 2):
    """
    Returns the min x span for required criterion
    Automatically converts min to max
    
    ALGORITHM HAS TO BE CHECKED, IMPROVED!!!
    """
    a = np.nan
    b = np.nan
    off = np.where(Array != 0)[0][0]
    Array = Array[Array != 0] # Removes Zero values
    if np.amax(Array) - np.amin(Array) < Criterion:
        asc_span = np.nan
        a = np.nan
        b = np.nan
    else:
        #import pdb; pdb.set_trace()
        min_tot = Array.argmin()
        max_tot = Array.argmax()
        if min_tot > max_tot:
            Array = Array[::-1]
            min_tot = Array.argmin()
            max_tot = Array.argmax()
        asc_span = max_tot - min_tot
        for i in range(max_tot - min_tot):
            where = np.where(Array > (Array[min_tot+i] + Criterion))[0]
            where = where[where > min_tot]
            if where.shape[0] == 0:
                
                break
            asc_ind = where[0]
            if (asc_ind - (min_tot+i)) < asc_span:
                asc_span = asc_ind - (min_tot+i)
                a = asc_ind
                b = min_tot+i
    #assert asc_span > 0 or np.isnan(asc_span), "asc_span is 0 or negative"
    if mode == 1:
        return asc_span
    elif mode == 2:
        #print a, Array.shape[0] ,Array.shape[0] - a
        return (asc_span, Array.shape[0] - a + off, Array.shape[0] - b + off)





