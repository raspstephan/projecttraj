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
import os
import scipy.ndimage as ndi
try:
    from scipy.interpolate import RegularGridInterpolator
except ImportError:
    print "Interpolate_2d not available"


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
      List containing data in np.arrays. data[0] contains trajectory id within
      savefile. data[1] contains starttime[mins].
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
    File directory =/home/scratch/users/stephan.rasp/Case1_20070720/d4deout_eight///
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
        self.dprcosmo = 5   # Cosmo output interval for p and r files
        halfsuff = '_30m'
        self.dhalfcosmo = 30   # Interval for half hour files
        rfiles = glob.glob(self.datadir + "*" + rsuff)
        rfiles = sorted(rfiles)
        pfiles = glob.glob(self.datadir + "*" + psuff)
        self.pfiles = sorted(pfiles)
        self.rfiles = [x for x in rfiles if x not in pfiles]
        self.cfile = glob.glob(self.datadir + '*00c*')[0]
        halffiles = glob.glob(self.datadir + "*" + halfsuff)
        self.halffiles = sorted(halffiles)
        tmpcfiles = glob.glob(self.datadir + '*00c*')
        afiles = sorted(glob.glob(self.datadir + "lfff*"))
        self.afiles = [x for x in afiles 
                       if x not in (self.rfiles + self.pfiles + tmpcfiles + 
                                    halffiles)]
        self.dacosmo = 60   # Interval for a(ll) files
        if len(glob.glob(self.datadir + '*00c*')) != 1:
            print 'More than one c file detected, take first one!'
            self.cfile = sorted(glob.glob(self.datadir + '*00c*'))[0]
        trjfiles = glob.glob(self.datadir + 'traj_*')
        # Take longest files (most added variables)
        maxlen = len(max(trjfiles, key=len))  # Get longest file
        trjfiles = [x for x in trjfiles if len(x) == maxlen]
        self.trjfiles = sorted(trjfiles)
        
        # Getting xlim, ylim from COSMO data
        tmpfobj = pwg.getfobj(self.cfile, 'FR_LAND_S')
        xlim = (tmpfobj.rlons[0][0], tmpfobj.rlons[0][-1])
        ylim = (tmpfobj.rlats[0][0], tmpfobj.rlats[-1][0])
        self.xlim = self._nrot2rot(xlim, 'lon')
        self.ylim = self._nrot2rot(ylim, 'lat')     
        
        # Set COSMO and trajectory output interval
        # Open temporary trj file to read information, if trjs present
        try:
            tmprg = nc.Dataset(trjfiles[0])
            self.date = dt.datetime(tmprg.ref_year, tmprg.ref_month, 
                                    tmprg.ref_day, tmprg.ref_hour)
            self.dtrj = int(tmprg.output_timestep_in_sec / 60)   # Minutes
        except IndexError:
            print 'No trajectory files found. Set date to default'
            self.date = dt.datetime(2000, 01, 01, 00)

        # Setting up lists and dictionaries
        self.datadict = dict(trjid = 0, startt = 1)
        self.filename = []
        self.data = [[], []]
        self.filtdict = dict()
        self.filtlist = []
        
        self.maxmins = (len(self.pfiles) - 1) * self.dprcosmo
        
        
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
    
    
    def new_max_diff(self, tracer, flip = False):
        """
        TODO
        """
        
        diffarray = utils._max_diff(self, tracer, flip)
        
        code = tracer + '_max_diff'
        self.datadict[code] = len(self.data)
        
        self.data.append(diffarray)
        print code, 'has been added.'
    
    
    def new_allasct(self, yspan, xmax, tracer = 'P'):
        """
        Creats new array in data, with informations about trajectory ascent 
        positions. 
        
        Parameters
        ----------
        yspan : float
          Criterion to be covered on y-axis
        xmax : float
          Maximum span on x-axis
        tracer : string
          COSMO name of y-axis variable
        """
        xsteps = int(xmax / self.dtrj)
        allmat = utils._allasct(self.trjfiles, yspan, xsteps, tracer, self.dtrj)
        print allmat.shape
        # Update dictionary
        code = tracer + str(yspan) + 'in' + str(xmax)
        self.datadict[code] = len(self.data)
        
        self.data.append(allmat)
        print code, 'has been added.'

    
    def new_allasct_cd(self, diff, sigma, tracer = 'z'):
        """
        Creats new array in data, with informations about trajectory ascent 
        positions. Filtered by velocities after Gaussian smoothing
        
        Parameters
        ----------
        diff : float
          threshold difference in one time 
        sigma : float
          sigma for gaussian smoothing
        tracer : string
          COSMO name of y-axis variable
          
        """
        
        allmat = utils._allasct_cd(self.trjfiles, diff, sigma, tracer, self.dtrj)
        print allmat.shape
        
        # Update dictionary
        code = tracer + str(diff) + 'per_min'
        self.datadict[code] = len(self.data)
        
        self.data.append(allmat)
        print code, 'has been added.'
        
        
    
    def new_loc_filter(self, xmin, xmax, ymin, ymax):
        """
        Adds a boolian array indicating if respective trajectories passed 
        through the given rectangle. 
        
        NOTE: For now in rotated coords.
        
        xmin : float
          Lower x boundary
        xmax : float
          Upper x boundary
        ymin : float
          Lower y boundary
        ymax : float
          Upper y boundary
        """
        
        boolarr = utils._loc_filter(self.trjfiles, xmin, xmax, ymin, ymax)
        
        code = str(xmin) + str(xmax) + str(ymin) + str(ymax)
        self.datadict[code] = len(self.data)
        
        self.data.append(boolarr)
        print code, 'has been added.'
        
    
    def new_delta(self, tracer, mode = 'minmax'):
        """
        Adds a new delta array to data.
        
        Parameters
        ----------
        tracer : str (default = 'P')
          COSMO name of y-axis variable
        mode : string
          Type of delta calculation
          'minmax' = diff between minimum and maximum
          'climb' = cumulative climbed value
          'climb_r' = reverse of climb
          'lifespan' = life span of the individual trajectories (tracer does not 
                       matter)
          
        """

        deltaarray = utils._delta(self.trjfiles, tracer, mode = mode)
        
        code = 'delta' + tracer + mode
        self.datadict[code] = len(self.data)
        
        assert (deltaarray.shape[0] == self.ntrj), \
                'Array shapes do not match. Look for error in source code.'
        if mode == 'lifespan':
            deltaarray *= self.dtrj
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
        if filtername == None:
            mask = np.array([True] * self.ntrj)  # New mask
        else:
            maskid = self.filtdict[filtername]
            mask = np.array(self.filtlist[maskid], copy = True)
        
        if type(addmask) == np.ndarray:
            mask &= addmask
        elif type(addmask) in [tuple, list]:
            for m in addmask:
                mask &= m
        
        uniqueloc= np.unique(self.filename[mask])
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
        
        if filtername == None:
            mask = np.array([True] * self.ntrj)  # New mask
        else:
            maskid = self.filtdict[filtername]
            mask = np.array(self.filtlist[maskid], copy = True)
        
        if addmask != None:
            mask &= addmask
        
        dataid = self.datadict[dataname]
        data = self.data[dataid]
        
        return data[mask]
    
    
    def _get_index(self, varname, mins):
        """
        Returns the index and filelist of Cosmofile containing the variable at 
        time: mins after model start. If variable appears in more than one 
        file type, take the one with smallest output interval. 
        
        Parameters
        ----------
        varname : string
          Name of variable
        mins : float
          Time in mins after model start
          
        Returns
        -------
        cosmoind : integer
          Index for correct cosmofile
        filelist : lists
          Filelist containing variable
          
        """
        
        if varname in pwg.get_fieldtable(self.rfiles[0]).fieldnames:
            cosmoind = int(mins / self.dprcosmo)
            filelist = self.rfiles
        elif varname in pwg.get_fieldtable(self.pfiles[0]).fieldnames:
            cosmoind = int(mins / self.dprcosmo)
            filelist = self.pfiles
        elif varname in pwg.get_fieldtable(self.halffiles[0]).fieldnames:
            cosmoind = int(mins / self.dhalfcosmo)
            filelist = self.halffiles
        elif varname in pwg.get_fieldtable(self.afiles[0]).fieldnames:
            cosmoind = int(mins / self.dacosmo)
            filelist = self.afiles
        else:
            raise Exception (varname, 'not found!')
        return cosmoind, filelist
        
        
    
    
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
    
    
    def interpolate_2d(self, varname):
        """
        Adds a surface variable as a tracer to trajectory files.
        
        NOTE: For now only variables in rfiles!!!
        
        Parameters
        ----------
        varname : string
          Name of surface variable
          
        """
        
        # Get rotated lon/lat and height fields for interpolation
        hhobj = pwg.getfobj(self.cfile, 'HH')
        lons = hhobj.rlons[0, :]
        lats = hhobj.rlats[:, 0]
        # hhfield = hhobj.data
        del hhobj
        
        # Create new filelist
        newlist = []
        for trjfn in self.trjfiles:
            newfn = trjfn.rstrip('.nc') + '_' + varname + '.nc'
            newlist.append(newfn)
            os.system('cp ' + trjfn + ' ' + newfn)
            
        # Start iteration over trajectory files
        for fn in newlist:
            print 'Opening file:', fn
            
            # Open trajectory file
            rootgrp = nc.Dataset(fn, 'a')
            
            # Read position matrices
            lon = rootgrp.variables['longitude'][:, :]
            lat = rootgrp.variables['latitude'][:, :]
            #z = rootgrp.variables['z'][:, :]
            
            # Allocate new netCDF array
            newvar = rootgrp.createVariable(varname, 'f4', ('time', 'id'))
            
            # Retrieve trajectory start time
            trjstart = rootgrp.variables['time'][0] / 60   # In mins       
            
            # Iterate over trj times
            for itrj in range(lon.shape[0]):
                if itrj % 100 == 0:
                    print 'Interpolating time step', itrj, 'of', lon.shape[0]
                
                # Get cosmo index
                icosmo = int((itrj * self.dtrj + trjstart) / self.dprcosmo)
                field = pwg.getfield(self.rfiles[icosmo], varname)
                field = np.transpose(field)   # flip x and y values
                
                # Get 2D position array
                lonlattrj = np.array([lon[itrj, :], lat[itrj, :]]) 
                lonlattrj = np.transpose(lonlattrj)
                
                # Interpolate with scipy.interpolate
                intobj = RegularGridInterpolator((lons, lats), field, 
                                                 fill_value = None)
                newvar[itrj, :] = intobj(lonlattrj)
            
            # Clean up 
            rootgrp.close()
        
        self.trjfiles = newlist
       
    def interpolate_3d(self, varname, outint = 60):
        """
        TODO
        """
        
        newlist = utils._interpolate_3d(self, varname, outint)
        self.trjfiles = newlist
        
                
    def get_val_start(self, ascname, tracer):
        """
        TODO
        """
        
        starttarray = self._mask_array(None, 'startt')
        startarray = (self._mask_array(None, ascname + '_start') - 
                        starttarray) / self.dtrj
        
        
        valarray = utils._get_val_start(self, startarray, tracer)
        
        code = tracer + '_at_' + ascname
        self.datadict[code] = len(self.data)
        
        self.data.append(valarray)
        print code, 'has been added.'
        
    
    
    def interpolate_value(self, totind, time, trjtracer, cosmotracer):
        """
        TODO
        """
        
        # Get correct time indices
        aind = int(time / self.dafiles)
        trjstart = self.data[1][totind]
        trjind = int((time - trjstart) / self.dtrj)
        
        
        # Get values from trajectory
        rootgrp = nc.Dataset(self.filename[totind], 'r')
        trjval = rootgrp.variables[trjtracer][trjind, self.data[0][totind]]
        lon = rootgrp.variables['longitude'][trjind, self.data[0][totind]]
        lat = rootgrp.variables['latitude'][trjind, self.data[0][totind]]
        z = rootgrp.variables['z'][trjind, self.data[0][totind]]
        
        print 'Trajectory value is: ', trjval
        
        # Get COSMO indices (closest to trj pos, no interpolation)
        fobj = pwg.getfobj(self.cfile, 'HH_S')
        lons = fobj.rlons[0, :]
        lats = fobj.rlats[:, 0]
        xind = np.argmin(np.abs(lons - lon))
        yind = np.argmin(np.abs(lats - lat))
        hh = pwg.getfield(self.cfile, "HH")
        zind = np.argmin(np.abs(hh[:, xind, yind] - z))
        
        # Get COSMO value
        cval = pwg.getfield(self.afiles[aind], cosmotracer)[zind, xind, yind]
        
        print 'COSMO value is: ', cval
        
        
    ######################
    # Plotting functions
    ######################
    
    
    def draw_vs_t(self, tracer, filtername, savename = None, sigma = None):
        """
        Draws tracer of trajectories in filter against time.
        
        """
        
        # Get iterable lists
        loclist, idlist = self._mask_iter(filtername)
        
        
        plots.draw_vs_t(self, tracer, loclist, idlist,
                        savename = savename, sigma = sigma)
  
  
    def draw_centered_vs_t(self, tracer, filtername, carray, 
                           savename = None, sigma = None):
        """
        TODO
        """
        loclist, idlist = self._mask_iter(filtername)
        startval = self._mask_array(filtername, carray + '_start')
        stopval = self._mask_array(filtername, carray + '_stop')
        carray = (startval + stopval) / 2
        plots.draw_centered_vs_t(self, loclist, idlist, tracer, carray, 
                                 savename)
 
    
    
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
        if savebase != None:
            savename = (savebase + '/scatter_' + dataname1 + '_' + dataname2 
                        + '_' + str(filtername) + '_' + str(idtext))
        else: 
            savename = savebase
        #xlabel = 'Time x ' + str(factor1) +' for ' + str(dataname1) + ' [hrs]'
        #ylabel = 'Time x ' + str(factor2) +' for ' + str(dataname2) + ' [hrs]'
        
        xlabel = dataname1
        ylabel = dataname2
        
        # Pass parameters to plots function
        plots.draw_scatter(array1, array2, carray, idtext, xlabel, ylabel, 
                           savename)
        
    
    def draw_avg(self, dataname, filtername, idtext = '', centdiff = False, 
                 savebase = None):
        """
        Draws the average of a certain paramters over all arrays in filter.
        NOTE: For now only works with arrays of same start time. Also, no 
        zero filter.
        
        Parameters
        ----------
        dataname : string
          Name of parameter. 
        filtername : string
          Name of filter to be applied
        idtext : string
          Text to be displayed in plot
        centdiff : bool
          If true centered difference of tracer will be plotted
        savebase : string
          Path to output directory
        """
        
        loclist, idlist = self._mask_iter(filtername)
        
        # Create savename
        if centdiff:
            cdtext = '_cd'
        else:
            cdtext = ''
        savename = (savebase + 'avg_' + dataname + cdtext + '_' + filtername + 
                    '_' + idtext)
        
        plots.draw_avg(dataname, loclist, idlist, idtext, centdiff, savename)
    
    
        
    def draw_hist(self, data, filtername = None, idtext = '', savebase = None,
                  log = False, **kwargs):
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
        log : bool
          If true x axis is logarithmic NOTE: range is hardcoded!
          
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
                if savebase != None:
                    savename = (savebase + '/hist_' + dataname + '_loc' 
                                + '_' + str(filtername) + '_' + str(idtext))
                else:
                    savename = savebase
                xlabel = ('Location [hPa] of fastest ' + dataname + ' ascent')
            
            # Regular histogram of given array
            else:
                if filtername != None:
                    array = self._mask_array(filtername, data)
                else:
                    array = self.data[self.datadict[data]]
                # Create savename and label names
                if savebase != None:
                    savename = (savebase + '/hist_' + data 
                                + '_' + str(filtername) + '_' + str(idtext))
                else:
                    savename = savebase
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
        plots.draw_hist(array, idtext, xlabel, savename, log = log, **kwargs)

    
    def draw_mult_hist(self, dataname, objs, filtername = None, idtext = '', 
                       savebase = None, xlabel = '', **kwargs):
        """
        TODO
        """
        
        # Apply filters to arrays and package them
        arraylist = [self._mask_array(filtername, dataname)]
        
        if type(objs) in [list, tuple]:
            for obj in objs:
                arraylist.append(obj._mask_iter(filtername, dataname))
        else:
            arraylist.append(objs._mask_iter(filtername, dataname))
        
        if savebase != None:    
            savename = savebase + 'multi_hist_' + dataname + '.png'
        else:
            savename = savebase

        plots.draw_mult_hist(arraylist, idtext, xlabel, savename, **kwargs)
        


    def draw_hist_2d(self, varlist, time, filtername = None, tracerange = None, 
                     idtext = '', savebase = None):
        """
        Draws a 2D histogram of trajectories over a contour map.
        
        Parameters
        ----------
        varlist : list
          List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
          'CUM_PREC' for cumulative precipitation
        time : float
          Time in minutes after model start
        filtername : string
          Identiefier of wanted filter
        tracerange : tuple
          Tuple eg ('P', 700, 1000), only consider trajectories within this 
          range.
        idtext : string
          Text to be displayed in plot
        savebase : string
          Path to output directory
          
        """
        
        tmask = self.data[1] < time
        filelist, idlist = self._mask_iter(filtername, addmask = tmask)
        
        
        if savebase != None:    
            savename = savebase + 'hist2d_' + str(time).zfill(4) + '.png'
        else:
            savename = savebase
            
        plots.draw_hist_2d(self, varlist, filelist, idlist, time, tracerange, 
                           idtext, savename)

    def draw_intersect_hor(self, level, leveltype = 'P', idtext = '', 
                           filtername = None, savebase = None):
        """
        TODO
        """
        
        filelist, idlist = self._mask_iter(filtername)
        
        if savebase != None:    
            savename = (savebase + 'intersect_hor_' + str(level) + 
                        '_' + idtext + '.png')
        else:
            savename = savebase
        
        velarray = plots.draw_intersect_hor(self, filelist, idlist, level, 
                                            leveltype, idtext, savename)
        return velarray
        
    
    
     
    def draw_trj_all(self, varlist, filtername = None, savebase = None, 
                     starts = False, onlyasc = None, trjstart = None, 
                     idtext = '', linewidth = 0.7, carray = 'P',
                     centdiff = False, sigma = None):
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
        
        plots.draw_trj(self, varlist, loclist, idlist, self.cfile,
                        self.rfiles, self.pfiles, savename = savename, 
                        pollon = self.pollon, pollat = self.pollat, 
                        xlim = self.xlim, ylim = self.ylim, onlybool = onlybool,
                        startarray = startarray, stoparray = stoparray, 
                        trjstart = trjstart, idtext = idtext, 
                        linewidth = linewidth, carray = carray,
                        centdiff = centdiff, sigma = sigma)
        
    def draw_trj_evo(self, varlist, filtername = None, tafter = None, 
                     interval = None, idtext = '', onlyasc = None, 
                     savebase = None):
        """
        Draws trajectories at certain times after trajectory start 
        with correct background plots. If several starting times are given in 
        filter, trajecoties are drawn at correct times.
        
        Parameters
        ----------
        varlist : list
          List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
          'CUM_PREC' for cumulative precipitation
        filtername : string
          Identiefier of wanted filter
        tafter : int
          Plotting time after simulation start (if interval: first plot time)
        interval : int
          If given, plots several figures, seperated by interval time
        idtext : string
          Text to be displayed in plot
        onlyasc : str
          If name of ascent property given, plots trajectories only during
          ascen time. E.g. 'P600'
        savebase : string
          Path to output directory
        
        """

        # Get plotting times [mins after simulation start]
        if interval == None:
            tlist = [tafter]
        else:
            if tafter == None:
                tlist = range(0, self.maxmins, interval)
            else: 
                tlist = range(tafter, self.maxmins, interval)
        
        # Get indeces of savefiles and trjids to be plotted
        loclist, idlist= self._mask_iter(filtername)  
        
        # Get stop and start arrays if onlyasc
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
        
        
        # Plot for each time
        for tplot in tlist: 
            # Create savename
            if savebase != None:  
                savename = (savebase + 'xy_' + filtername + '_' + 
                            str(tplot).zfill(4) + '.png')
            else:
                savename = savebase

            plots.draw_trj_evo(self, varlist, loclist, idlist, tplot, 
                               idtext = idtext, onlybool = onlybool, 
                               startarray = startarray, stoparray = stoparray, 
                               savename = savename)
    
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
    
    def draw_asc_loc(self, dataname, varlist, filtername, tplot, tspan, 
                     idtext = '', savebase = None):
        """
        Draws ascent locations as dots for all trjs in filter where ascent is 
        happening within tspan of tplot.
        
        Parameters
        ----------
        dataname : string
          Name of ascent locations to be plotted. (Created with new_allasct)
        varlist : list
          List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
          'CUM_PREC' for cumulative precipitation
        filtername : string
          Identiefier of wanted filter
        tplot : float
          Time [in mins after model start] of plot
        tspan : float
          Tolerance interval. All asc plotted which are in 
          [tplot - tspan, tplot + tspan]
        idtext : string
          Text to be displayed in plot
        savebase : string
          Path to output directory
          
        """
        
        # Filter by filter
        mat = self._mask_array(filtername, dataname)
        
        # Rearrange data
        newlist = []
        for i in list(mat):
            for j in i:
                newlist.append(j)
        mat = np.array(newlist)
        
        # Filter by time
        tavg = (mat[:, 6] + mat[:, 7] ) / 2 
        mask = tavg >= (tplot - tspan)
        mask &= tavg <= (tplot + tspan)
        
        # Filter
        
        # Get avg arrays
        lonavg = (mat[:, 2] + mat[:, 3]) / 2
        latavg = (mat[:, 4] + mat[:, 5]) / 2
        pavg = (mat[:, 8] + mat[:, 9]) / 2
        lonavg = lonavg[mask]
        latavg = latavg[mask]
        pavg = pavg[mask]
        
        # Create Savename
        if savebase != None:  
            # TODO Change name
            savename = (savebase + 'xy_' + filtername + '_' + dataname + '_' +
                        str(tplot).zfill(4) + '.png')
        else:
            savename = savebase
        
        plots.draw_asc_loc(self, lonavg, latavg, pavg, varlist, tplot, idtext, 
                           savename)
    
    
    
        
    def draw_contour(self, varlist, time, savebase = None, interval = None, 
                     idtext = ''):
        """
        Draws a countourplot of given variables.
        
        Parameters
        ----------
        varlist : list
          List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
        time : integer
          In minutes after model start
        savebase : string
          Path to output directory
        Interval : integer
          If specified, plots for all in range(time, self.maxmins, interval)
        idtext : string
          Text to be diplayed
          
        """
        if savebase != None:    
            savename = savebase + 'contour_' + str(time) 
        else:
            savename = savebase
        
        if interval == None:
            timelist = [time]
        else:
            timelist = range(time, self.maxmins, interval)
        for time in timelist:
            if savebase != None:    
                savename = savebase + 'contour_' + str(time).zfill(5)
            else:
                savename = savebase
            print 'Plotting for time:', time
            plots.draw_contour(self, varlist, time, savename = savename, 
                               idtext = idtext)
            
            
    ##################
    # String functions
    ##################
            
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
    
    
    
    
    
    
    

