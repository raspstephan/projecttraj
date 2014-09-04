"""
core module
-----------

This module contains the heart of the package: the TrjObj class.

"""

import cPickle
import pickle
import numpy as np
import netCDF4 as nc
import glob
import plots
import utils
import cosmo_utils.pywgrib as pwg
import datetime as dt
import re


class TrjObj(object):
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
    pollon : float
      Rotated pole longitude. Default value does not change input coordinates.
    pollat : float
      Rotated pole latitude. Default value does not change input coordinates.
      
    
    Attributes
    ----------
    rfiles : list
      List of regular (meaning non-pressure levels) COSMO output files
    pfiles : list
      List of pressure-level COSMO output files
    cfile : string
      Location of COSMO constants file
    trjfiles : list
      List of trajectory netCDF files
    ntrj : int
      Total number of trajectories
    filename : np.array
      Sorted array containing the file location for each trajectory
    trjid : np.array
      Sorted array containing the unique ID for each trajectory
    data : list
      List containing data in np.arrays
    datadict : dictionary
      Dictionary containing data names and matching IDs in data
    filtlist : list
      List containing filters/masks in np.arrays
    filtdict : dictionary
      Dictionary containing filter names and matching filter IDs in filtlist
    
    
    Examples
    --------
    >>> import traj_tools as trj
    >>> # Define directory containing trajectory and COSMO files
    >>> datadir = '/home/scratch/users/stephan.rasp/Case1_20070720/d4deout_eight'
    >>> # Define lon and lat of rotated pole (not necessary, but strongly recommended!)
    >>> rotlon = -170
    >>> rotlat = 40
    >>> # Allocate class object
    >>> case1 = trj.TrajPack(datadir, rotlon, rotlat)
    
    Get an overview:
    
    >>> case1
    >>> # print case1 gives the same output
    ---------------
    TrajPack object
    ---------------
    File directory =        /home/scratch/users/stephan.rasp/Case1_20070720/d4deout_eight///
    Contains the following:
    Data arrays =   ['P400', 'startt', 'P600', 'P300']
    Filters =       ['WCB_360', 'WCB']
    ---------------

    

    """
    
    def __init__(self, datadir, pollon = 180., pollat = 90.):

        self.datadir = datadir + '/'
        self.pollon = pollon   # Default: 180
        self.pollat = pollat   # Default: 90
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
        if len(glob.glob(self.datadir + '*00c*')) != 1:
            print 'More than one c file detected, take first one!'
            self.cfile = sorted(glob.glob(self.datadir + '*00c*'))[0]
        trjfiles = glob.glob(self.datadir + 'traj_*')
        if 'theta' in ''.join(trjfiles):
            for i in range(len(trjfiles)):
                trjfiles = [x for x in trjfiles if 'theta' in x]
        self.trjfiles = sorted(trjfiles)
        
        # Getting xlim, ylim from COSMO data
        tmpfobj = pwg.getfobj(self.cfile, 'FR_LAND_S')
        xlim = (tmpfobj.lons[0], tmpfobj.lons[-1])
        ylim = (tmpfobj.lats[0], tmpfobj.lats[-1])
        self.xlim = self._nrot2rot(xlim, 'lon')
        self.ylim = self._nrot2rot(ylim, 'lat')     
        
        # Set COSMO and trajectory output interval
        # Open temporary trj file to read information
        tmprg = nc.Dataset(trjfiles[0])
        self.date = dt.datetime(tmprg.ref_year, tmprg.ref_month, 
                                tmprg.ref_day, tmprg.ref_hour)
        self.dtrj = int(tmprg.output_timestep_in_sec / 60)   # Minutes
        self.dcosmo = 5

        # Setting up lists and dictionaries
        self.datadict = dict(trjid = 0, startt = 1)
        self.filename = []
        self.data = [[], []]
        self.filtdict = dict()
        self.filtlist = []
        
        self.maxmins = 6120
        
        
        # Looping over all files, initializing lists instead of np.arrays
        # in order to dynamically add files
        
        ntot = 0   # Counter for assert check below
        
        for i in range(len(self.trjfiles)):
            rootgrp = nc.Dataset(self.trjfiles[i], 'r')
            ntmp = len(rootgrp.dimensions['id'])
            tmpt = self.trjfiles[i].split('/')[-1]
            tmpt = int(tmpt.split('_')[1].lstrip('t'))
            self.filename.extend([self.trjfiles[i]] * ntmp)
            self.data[0].extend(range(ntmp))
            self.data[1].extend([tmpt] * ntmp)
            ntot += ntmp
            
        # Convert lists to np.array
        self.ntrj = ntot
        self.filename = np.array(self.filename)
        self.data[0] = np.array(self.data[0])
        self.data[1] = np.array(self.data[1])

        assert (self.data[0].shape[0] == self.filename.shape[0] == 
                ntot), \
                "Error while initializing properties for class: Attribute arrays do not have same shape!"
                               
  
    def new_array(self, name, array):
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
          
        Examples
        --------
        >>> case1.new_prop_array('id', np.array(range(case1.ntrj)))
        id has been added.
        
        """
        
        if type(name) == str:
            assert (array.shape[0] == self.ntrj), \
                    'Array shaped do not match.'
            assert (name not in self.datadict), 'Name already exists.' 
            self.datadict[name] = len(self.data)
            self.data.append(array)
            print name, 'has been added.'
            
            
        elif type(name) == tuple or type(name) == list:
            for i in range(len(name)):
                assert (name[i] not in self.datadict), \
                        'One of the names already exists.'
                assert (array[i].shape[0] == ntrj), \
                        'One or all of the array shapes do not match.'
            self.datadict.update(dict(zip(name), 
                                     range(len(self.data), 
                                           len(self.data) + len(name))))
            self.data.extend(array)
            print name, 'have been added.'
        
    
    
    def new_asc(self, yspan, tracer = 'P', interpolate = True):
        """
        Adds a new ascent filter to data.
        
        Parameters
        ----------
        yspan : float
          Ascent criterion in y-direction
        tracer : str (default = 'P')
          COSMO name of y-axis variable
        
        Examples
        --------
        >>> case1.new_prop_asc(600, 'P')
        P600 has been added
        
        """

        ascdata = utils._minasct(self.trjfiles, yspan, tracer, self.dtrj, 
                                 interpolate)
        
        # Update dictionary
        code = tracer + str(yspan)
        self.datadict[code] = len(self.data)
        self.datadict[code + '_start'] = len(self.data) + 1
        self.datadict[code + '_stop'] = len(self.data) + 2
        self.datadict[code + '_start_val'] = len(self.data) + 3
        self.datadict[code + '_stop_val'] = len(self.data) + 4
        
        assert (ascdata[0].shape[0] == self.ntrj), \
                'Array shapes do not match. Look for error in source code.'
        self.data.extend(ascdata)
        print code, 'has been added.'
    
    
    def new_delta(self, tracer):
        """
        Adds a new delta array to data.
        
        Parameters
        ----------
        tracer : str (default = 'P')
          COSMO name of y-axis variable
          
        """

        deltaarray = utils._delta(self.trjfiles, tracer)
        
        code = 'delta' + tracer
        self.datadict[code] = len(self.data)
        
        assert (deltaarray.shape[0] == self.ntrj), \
                'Array shapes do not match. Look for error in source code.'
        self.data.append(deltaarray)
        print code, 'has been added.'
        
    
    
       
    def create_filter(self, name, filters):
        """
        Adds a filter to filtlist. Adds name to filtdict.
        
        Parameters
        ----------
        name : string
          Name of given filter
        filters : list
          List of filter tuples, see example.
        
        Examples
        --------
        >>> case1.create_filter('WCB_start360', [('P600', 0, 2880), ('startt', 360)])
        WCB_startt360 has been added.
        
        If only one filter is given, it still has to be a tuple in a list, like
        
        >>> case1.create_filter('WCB', [('P600', 0, 2880)])
        WCB has been added.
        
        """
        
        assert (name not in self.filtdict), 'Filter name already exists.'
        
        mask = np.array([True] * self.ntrj)   # Initialize mask 
        
        for i in range(len(filters)):

            ind = self.datadict[filters[i][0]]   # Get index in self.data
            
            if len(filters[i]) == 3:
                mn = filters[i][1]
                mx = filters[i][2]
                
                mask &= self.data[ind] <= mx
                mask &= self.data[ind] >= mn
                
            elif len(filters[i]) == 2:
                crit = filters[i][1]
                
                mask &= self.data[ind] == crit
                
            else:
                raise Exception('Wrong input for filter tuple.')
       
        self.filtdict[name] = len(self.filtlist)
        self.filtlist.append(mask)
        print name, 'has been added.'


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
    
    
        
    def _mask_iter(self, filtername, addmask = None):
        """
        Converts mask to iterable lists. To be used in plotting functions, etc.
        
        Parameters
        ----------
        filtername : np.array
          Mask 
        addmask : logical np.array
          Additional mask
        
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
        mask = np.array(self.filtlist[maskid], copy = True)
        
        if type(addmask) == np.ndarray:
            mask &= addmask
        elif type(addmask) in [tuple, list]:
            for m in addmask:
                mask &= m
        
        uniqueloc = np.unique(self.filename[mask])
        idlist = []
        for i in range(len(uniqueloc)):
            locmask = (self.filename == uniqueloc[i])
            idlist.append(list((self.data[0][locmask & mask])))

        return list(uniqueloc), idlist
    
    
    def _mask_array(self, filtername, dataname, addmask = None):
        """
        Returns one masked array specified by given mask and dataname.
        
        Parameters
        ----------
        filtername : np.array
          Name of filter in dictionary
        dataname : string
          Identifier of data array
        addmask : logical np.mask
          Additional mask 
          
        Returns
        -------
        marray         : np.array
          Masked data array
        
        """
        
        maskid = self.filtdict[filtername]
        mask = np.array(self.filtlist[maskid], copy = True)
        
        if addmask != None:
            mask &= addmask
        
        dataid = self.datadict[dataname]
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
        
        Examples
        --------
        >>> case1.saveme('./case1.tp')
        Saved as ./case1.tp
        """
        
        f = open(savename, 'w+')
        # 1. Save datadir
        pickle.dump(self.datadir, f, 2)
        # 2. Save pollon, pollat
        pickle.dump((self.pollon, self.pollat), f, 2)
        # 3. Save xlim, ylim
        pickle.dump((self.xlim, self.ylim), f, 2)
        # 4. Save data
        pickle.dump(self.data, f, 2)
        # 5. Save datadict
        pickle.dump(self.datadict, f, 2)
        # 6. Save filtlist
        pickle.dump(self.filtlist, f, 2)
        # 7. Save filtdict
        pickle.dump(self.filtdict, f, 2)
        f.close()
        print 'Saved as', savename
        
    def calc_theta(self):
        """
        Class method to add Potential Temperature as a Variable to the 
        netCDF files.
        
        Examples
        --------
        >>> case1.calc_theta()
        """
        
        thetalist = utils.calc_theta(self.trjfiles)
        self.trjfiles = thetalist
        
        
    ######################
    # Plotting functions
    ######################
    
    
    def draw_vs_t(self, dataname, totind, savename = None):
        """
        """
        plots.draw_vs_t(dataname, self.filename[totind], self.data[0][totind],
                        savename = savename)
        
    def draw_scatter(self, dataname1, dataname2, factor1 = 1, factor2 = 1, 
                     carray = None, filtername = None, idtext = '', 
                     savebase = None):
        """
        Make a scatter plot of two data arrays, multiplied by factors.
        
        Parameters
        ----------
        
        dataname1 : string
          Name of parameter on x-axis
        dataname2 : string
          Name of parameter on y-axis
        factor1 : float
          Multiplication factory for first parameter
        factor2 : float
          Multiplication factory for second parameter
        carray : string
          Name of criterion for color coding
        filtername : string
          Name of filter to be applied 
        idtext : string
          Text to be displayed in plot
        savebase : string
          Path to output directory
          
        """
        
        # Retrieve the two arrays
        if filtername != None:
            array1 = self._mask_array(filtername, dataname1)
            array2 = self._mask_array(filtername, dataname2)
        else:
            array1 = self.data[self.datadict[dataname1]]
            array2 = self.data[self.datadict[dataname2]]
        
        # Multiply by factor if given
        array1 = array1 * factor1
        array2 = array2 * factor2
        
        # Get carray
        if carray != None:
            if filtername != None:
                startval = self._mask_array(filtername, carray + '_start_val')
                stopval = self._mask_array(filtername, carray + '_stop_val')
            else:
                startval = self.data[self.datadict[carray + '_start_val']]
                stopval = self.data[self.datadict[carray + '_stop_val']]
            carray = (startval + stopval) / 2
        
        # Create savename and label names
        savename = (savebase + '/scatter_' + dataname1 + '_' + dataname2 + '_' +
                    str(filtername) + '_' + str(idtext))
        xlabel = 'Time x ' + str(factor1) +' for ' + str(dataname1) + ' [hrs]'
        ylabel = 'Time x ' + str(factor2) +' for ' + str(dataname2) + ' [hrs]'
        
        # Pass parameters to plots function
        plots.draw_scatter(array1, array2, carray, idtext, xlabel, ylabel, 
                           savename)
        
        
        
    def draw_hist(self, data, filtername = None, idtext = '', savebase = None):
        """
        Make a histogram of the ratio of two parameters.
        dataname1 / dataname 2 adjusted by multiplication factor.
        
        Parameters
        ----------
        
        data : string or list
          Name of data. If '_val' is at the end, takes the average value.
          Or list with [dataname1, dataname2, factor 1, factor2]
        filtername : string
          Name of filter to be applied 
        idtext : string
          Text to be displayed in plot
        savebase : string
          Path to output directory
          
        """
        
        # Check data type
        if (type(data) == str):
            
            # Histogram of ascent location
            if '_val' in data:
                dataname = data[:4]   # NOTE: Lazy, could lead to problems
                if filtername != None:
                    startval = self._mask_array(filtername, 
                                                dataname + '_start_val')
                    stopval = self._mask_array(filtername, 
                                               dataname + '_stop_val')
                else:
                    startval = self.data[self.datadict[dataname + '_start_val']]
                    stopval = self.data[self.datadict[dataname + '_stop_val']]
                array = (startval + stopval) / 2
                savename = (savebase + '/hist_' + dataname + '_loc' 
                            + '_' + str(filtername) + '_' + str(idtext))
                xlabel = ('Location [hPa] of fastest ' + dataname + ' ascent')
            
            # Regular histogram of given array
            else:
                if filtername != None:
                    array = self._mask_array(filtername, data)
                else:
                    array = self.data[self.datadict[dataname]]
                # Create savename and label names
                savename = (savebase + '/hist_' + data 
                            + '_' + str(filtername) + '_' + str(idtext))
                xlabel = data
                    
        # Histogram of ratio of two arrays            
        elif len(data) == 4:
            dataname1 = data[0]
            dataname2 = data[1]
            factor1 = data[2]
            factor2 = data[3]
            
            # Retrieve the two arrays
            if filtername != None:
                array1 = self._mask_array(filtername, dataname1)
                array2 = self._mask_array(filtername, dataname2)
            else:
                array1 = self.data[self.datadict[dataname1]]
                array2 = self.data[self.datadict[dataname2]]
            # Multiply by factor if given and get ratio
            array1 = array1 * factor1
            array2 = array2 * factor2
            array = array1 / array2
            # Create savename and label names
            savename = (savebase + '/hist_' + dataname1 + 'vs' + dataname2 + '_' 
                        + str(filtername) + '_' + str(idtext))
            xlabel = ('Ratio ' + dataname1 + ' x ' + str(factor1) + ' / ' + 
                      dataname2 + ' x ' + str(factor2))
            
        # Unknown input
        else:
            raise Exception('Wrong input for data')

        # Pass parameters to plots function
        plots.draw_hist(array, idtext, xlabel, savename)

        
     
    def draw_trj_all(self, varlist, filtername = None, savebase = None, 
                     starts = False, onlyasc = None, trjstart = None, 
                     idtext = '', linewidth = 0.7):
        """
        Draws XY Plot of trajectories with color as a function of 'P'.
        If filtername is given, plots only filetered trajectories.
        If starts is True, plots one xy plot seperately for each start time. 
        
        Parameters
        ----------
        varlist : list
          List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
          'CUM_PREC' for cumulative precipitation
        filtername : string
          Identiefier of wanted filter
        savebase : string
          Path to output directory
        starts : bool
          If False, plots xy plot for entire filtered array (NOT recommended,
          if start time is not used in filter!)
          If True, plots xy plot seperately for each start time.
        onlyasc : str
          If name of ascent property given, plots trajectories only during
          ascen time. E.g. 'P600'
        
        """
        assert (self.xlim != None and self.ylim != None), \
                'xlim and/or ylim are not specified!'
        
        #if starts:
        #Not possible until mem leak is resolved!!!
            
        
        if savebase != None:    
            savename = savebase + 'xy_' + filtername + '.png'
        else:
            savename = savebase
            
        loclist, idlist = self._mask_iter(filtername)
        
        if onlyasc != None:
            starttarray = self._mask_array(filtername, 'startt')
            startarray = (self._mask_array(filtername, onlyasc + '_start') - 
                          starttarray) / self.dtrj
            stoparray = (self._mask_array(filtername, onlyasc + '_stop') - 
                          starttarray) / self.dtrj

            onlybool = True
        else:
            startarray = None
            stoparray = None
            onlybool = False
        
        plots.draw_trj(varlist, loclist, idlist, self.cfile,
                        self.rfiles, self.pfiles, savename = savename, 
                        pollon = self.pollon, pollat = self.pollat, 
                        xlim = self.xlim, ylim = self.ylim, onlybool = onlybool,
                        startarray = startarray, stoparray = stoparray, 
                        trjstart = trjstart, idtext = idtext, 
                        linewidth = linewidth)
        
    def draw_trj_evo(self, varlist, tafter = None, interval = None, 
                     filtername = None, savebase = None, trjstart = None):
        """
        Draws trajectories at certain times after trajectory start 
        with correct background plots
        """
        
        # if trjstart is None, check if all trajectories start at same time
        # TODO
        
        if interval == None:
            tlist = [tafter]
        else:
            tlist = range(0, self.maxmins, interval)
            
        # if trjstart is given, create temporary mask
        tmpmask = None
        if not trjstart == None:
            tmpmask = self.data[self.datadict['startt']] == trjstart
            
        loclist, idlist = self._mask_iter(filtername, addmask = tmpmask)
        
        for t in tlist:
            
            if savebase != None:  
                # TODO Change name
                savename = (savebase + 'xy_' + filtername + '_' + 
                            str(t).zfill(4) + '.png')
            else:
                savename = savebase

            plots.draw_trj_evo(varlist, loclist, idlist, t, self.cfile, 
                            self.rfiles, self.pfiles, savename = savename, 
                            pollon = self.pollon, pollat = self.pollat, 
                            xlim = self.xlim, ylim = self.ylim, 
                            dtrj = self.dtrj, dcosmo = self.dcosmo)
        
    
    def draw_trj_dot(self, varlist, tplus = None, interval = None, 
                     filtername = None, savebase = None, trjstart = None,
                     onlyasc = None):
        """
        Draws trajectoriy position as dots with correct background plots.
        Tplus is now time after model start!
        """
        
        # if trjstart is None, check if all trajectories start at same time
        # TODO
        
        if interval == None:
            tlist = [tplus]
        else:
            tlist = range(0, self.maxmins, interval)
            
        # if trjstart is given, create temporary mask
        tmpmask1 = None
        if not trjstart == None:
            tmpmask1 = self.data[self.datadict['startt']] == trjstart
        masklist = [tmpmask1]
        ascind = self.datadict[onlyasc]
        for t in tlist:
           
            tmpmask2 = (t > self.data[ascind+1]) & (t < self.data[ascind+2])
            #print (t > self.data[ascind+1]) & (t < self.data[ascind+2])
            masklist = [x for x in masklist if not x == None]
            #print masklist
            
            loclist, idlist = self._mask_iter(filtername, 
                                              addmask = tmpmask2)

            
            
            if savebase != None:  
                # TODO Change name
                savename = (savebase + 'xy_' + filtername + '_' + 
                            str(t).zfill(4) + '.png')
            else:
                savename = savebase

            plots.draw_trj_dot(self,varlist, loclist, idlist, t, 
                               savename = savename)
    
    
        
    def draw_contour(self, varlist, time, savebase = None, interval = None,
                     trjstart = None):
        """
        Draws a countourplot of given variables.
        
        Parameters
        ----------
        varlist : list
          List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
        time : integer
          In steps after model start
        savebase : string
          Path to output directory
        Interval : integer
          NOT IMPLEMENTED
        """
        if savebase != None:    
            savename = savebase + 'contour_' + str(time) + '.png'
        else:
            savename = savebase
        
        if interval == None:
            timelist = [time]
        else:
            timelist = range(time, self.maxmins, interval)
        for time in timelist:
            if savebase != None:    
                savename = savebase + 'contour_' + str(time).zfill(4) + '.png'
            else:
                savename = savebase
            print 'Plotting for time:', time
            plots.draw_contour(self, varlist, time, savename = savename)
            
            
            
            
    def __repr__(self):
        
        return self.__str__()
    
    def __str__(self):
        """
        What you see if you type the object name or print(obj)
        """
        
        o =  '---------------\n'
        o += 'TrajPack object\n'
        o += '---------------\n'
        
        o += 'File directory = \t' + self.datadir + '\n'
        o += 'Contains the following:\n'
        tmpdata = self.datadict.keys()
        tmpdata = [x for x in tmpdata if not '_stop' in x]
        tmpdata = [x for x in tmpdata if not '_start' in x]
        o += 'Data arrays =\t' + str(tmpdata) + '\n'
        o += 'Filters = \t' +str(self.filtdict.keys()) + '\n'
        
        o += '---------------'
        
        return o
        
    

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
    # 1. Load datadir
    datadir = pickle.load(f)
    # 2. Load pollon, pollat
    pollon, pollat = pickle.load(f)
    # 3. Load xlim, ylim
    xlim, ylim = pickle.load(f)
    # Allocate object
    obj = TrjObj(datadir, pollon, pollat)
    # 4. Load data
    obj.data = pickle.load(f)
    # 5. Load datadict
    obj.datadict = pickle.load(f)
    # 6. Load filtlist
    obj.filtlist = pickle.load(f)
    # 7. Load filtdict
    obj.filtdict = pickle.load(f)
    f.close()
    
    return obj 





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
    
    
    
    
    
    
    

