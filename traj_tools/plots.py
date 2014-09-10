"""
plots module
------------

This module contains all plotting functions.
To make movie use from command line.

ffmpeg -r 2 -pattern_type glob -i '*.png' -c:v libx264 movie.mkv

"""


import numpy as np
import matplotlib
import os
# taken from cosmo_utils:
try:
    os.environ["DISPLAY"]
    print "X-Server detected, using tkagg backend for plotting"
except KeyError:
    if matplotlib.get_backend() in matplotlib.rcsetup.interactive_bk:
        matplotlib.use("Agg") 
        # interface "Agg" can plot without x-server connection
    print "No X-Server detected, using agg backend for plotting"
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.collections as col
from mpl_toolkits.basemap import Basemap
import cosmo_utils.pywgrib as pwg
import scipy.ndimage as ndi
from datetime import datetime, timedelta
import netCDF4 as nc


def draw_vs_t(dataname, fileloc, fileid, savename = None):
    """
    """
    rootgrp = nc.Dataset(fileloc, 'r')
    array = rootgrp.variables[dataname][:, fileid]
    plt.plot(array)
    if savename != None:
        plt.savefig(savename)
        plt.close('all')

def draw_scatter(array1, array2, carray = None, idtext = '', xlabel = None, 
                 ylabel = None, savename = None):
    """
    Returns/Saves a scatter plot of two arrays.
    
    Parameters
    ----------
    
    array1 : np.array
      Data on x-axis
    array2 : np.array
      Data on y-axis
    carray : np.array
      Array for color coding
    idtext : string
      Text to be displayed in plot
    xlabel : string
      Text for xlabel
    ylabel : string
      Text for ylabel
    savename : string
      Full path of file to be saved
      
    """
    
    # Set up figure size
    fig = plt.figure(figsize = (10, 10))
    ax = plt.gca()   
    #ax.set_aspect('equal')
    if carray == None:
        carray = 'b'   # Set to default
        
    # Convert to hours
    array1 = array1 / 60
    array2 = array2 / 60
    
    # Plot scatter plot, add labels
    sca = ax.scatter(array1, array2, c = carray, s = 8,
                     cmap=plt.get_cmap('Spectral'), 
                     norm=plt.Normalize(100, 1000), linewidths = 0)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Set limits, plot diagnonal line, set ticks
    xmx = np.amax(array1)
    ymx = np.amax(array2)
    mx = max(xmx, ymx)
    plt.xlim(0, xmx)
    plt.ylim(0, ymx+1)
    plt.plot([0, mx], [0, mx])
    plt.xticks(np.arange(0, xmx+3, 3))
    plt.yticks(np.arange(0, ymx+3, 3))
    
    # Add idtext
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    cb = fig.colorbar(sca, shrink = 0.7)
    cb.set_label('p')
    cb.ax.invert_yaxis()
    ax.set_axis_bgcolor('grey')
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename)
        plt.close('all')


def draw_hist(array, idtext = '', xlabel =  None, savename = None, log = False):
    """
    Returns/Saves a histogram of given array
    
    array : np.array
      Array to be used for histogram
    idtext : string
      Text to be displayed in plot
    xlabel : string
      Text for xlabel
    savename : string
      Full path of file to be saved
      
    """
    # Set up figure
    fig = plt.figure()
    
    # Plot histogram, remove nans
    if log:
        plt.hist(array[np.isfinite(array)], bins = np.logspace(-1, 3, 144))
        plt.gca().set_xscale('log')
    else:
        plt.hist(array[np.isfinite(array)], bins = 144)
        
    
    # Add labels and text
    plt.ylabel("Number of trajectories")
    plt.xlabel(xlabel)
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename)
        plt.close('all')
        plt.clf()
        

def draw_contour(obj, varlist, time, idtext, savename = None):
    """
    Plots a contour plot of the given variables at the specified time.
    
    Parameters
    ----------
    obj : TrjObj object
      Object
    varlist : list
      List of variables to be plotted
    time : int
      Time in minutes after simulation start
    savename : string
      Full path to file to be saved
      
    """
    
    # Getting index for COSMO files
    cosmoind = int(time / obj.dcosmo)
    
    # Setting up figure
    fig = plt.figure(figsize = (12,8))
    ax = plt.gca()   
    ax.set_aspect('equal')
    
    # Draw basemap
    basemap(obj.cfile, obj.xlim, obj.ylim)
    
    # Plotting all contour fields
    for i in range(len(varlist)):   
        # NOTE: 'CUM_PREC' not implemented right now
        #if varlist[i] == 'CUM_PREC':
            #contour(rfiles, varlist[i], cosmoind, xlim, ylim, 
                    #trjstart = trjstart)
        # Check if variable is in rfiles
        if varlist[i] in pwg.get_fieldtable(obj.rfiles[cosmoind]).fieldnames:
            contour(obj.rfiles, varlist[i], cosmoind, obj.xlim, obj.ylim)
        # Check if variable is in pfiles
        elif varlist[i] in pwg.get_fieldtable(obj.pfiles[cosmoind]).fieldnames:
            contour(obj.pfiles, varlist[i], cosmoind, obj.xlim, obj.ylim)
        else:
            raise Exception('Variable' + varlist[i] + 'not available!')
        
    # Set plot properties
    plt.xlim(obj.xlim)
    plt.ylim(obj.ylim)
    
    # Set labels and title
    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.title(obj.date + timedelta(minutes = time))
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    
        
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400)
        plt.close('all')
        plt.clf()


def draw_trj(obj, varlist, filelist, idlist, cfile, rfiles, pfiles, 
             savename = False,pollon = None, pollat = None, xlim = None, 
             ylim = None, onlybool = False, startarray = None, 
             stoparray = None, trjstart = None, idtext = '', linewidth = 0.7):
    """
    Plots one xy plot with trajectories.
    
    Parameters
    ----------
    obj : TrjObj object
      self object
    filelist : list
      List of unique file locations
    idlist : list
      List of list of trajectory IDs 
    cfile : string
      Constants COSMO file
    rfiles : list
      List of regular COSMO output files
    pfiles : list
      List of pressure lvl COSMO output files
    savename : string
      Name of save location
    pollon : float
      Rotated pole longitude
    pollat : float
      Rotated pole latitude
    xlim : tuple
      Tuple with x axis bounds
    ylim : tuple
      Tuple with y axis bounds
    onlybool : bool
      If True trajectories will be plotted during the ascent time only.
    startarray : np.array
      
    """
    
    # Get trj_time after model start 
    tmpt = filelist[0].split('/')[-1]
    tmpt = int(tmpt.split('_')[1].lstrip('t'))
    
    # Plot contours
    draw_contour(obj, varlist, tmpt, idtext = idtext)
    
    # Plot trajectories
    cnt = 0   # continuous counter for startarray and stoparray
    for i in range(len(filelist)):
        print 'Plotting file', i+1, 'of', len(filelist)
        rootgrp = nc.Dataset(filelist[i], 'r')
        lonmat = rootgrp.variables['longitude'][:, :]
        latmat = rootgrp.variables['latitude'][:, :]
        pmat = rootgrp.variables['P'][:, :]
    
        lonmat[:, :] += (180 - pollon)   # Convert to real coordinates
        latmat[:, :] += (90 - pollat)
        
        for j in idlist[i]:
            # Filter out zero values!
            if onlybool:
                parray = pmat[startarray[cnt]:stoparray[cnt], j]
                lonarray = lonmat[startarray[cnt]:stoparray[cnt], j]
                latarray = latmat[startarray[cnt]:stoparray[cnt], j]
                cnt += 1
            else:
                parray = pmat[:, j][pmat[:, j] != 0]
                lonarray = lonmat[:, j][pmat[:, j] != 0]
                latarray = latmat[:, j][pmat[:, j] != 0]
            
            single_trj(lonarray, latarray, parray, linewidth = linewidth)

    cb = plt.colorbar(lc, shrink = 0.7)
    cb.set_label('p')
    cb.ax.invert_yaxis()
    plt.tight_layout()
    
    
    # Save Plot
    if savename != False:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400)
        plt.close('all')
        plt.clf()


def draw_trj_evo(varlist, loclist, idlist, tplus, cfile, rfiles, 
                 pfiles, savename = None, pollon = None, pollat = None, 
                 xlim = None, ylim = None, dtrj = None, dcosmo = None):
    """
    TODO
    """
    
    # Setting up indices
    rootgrp = nc.Dataset(loclist[0], 'r')
    off = int(rootgrp.variables['time'][0] / 60) # also trjstart time
    cosmoind = (off + tplus) / dcosmo
    trjind = tplus / dtrj   
    
    # Set up figure
    fig = plt.figure(figsize = (12,8))
    ax = plt.gca()   
    ax.set_aspect('equal')
    basemap(cfile, xlim, ylim)
    
    
    # Plotting all contour fields
    for i in range(len(varlist)):   # Plotting all contour fields
        if varlist[i] == 'CUM_PREC':
            contour(rfiles, varlist[i], cosmoind, xlim, ylim, 
                    trjstart = trjstart)
        elif varlist[i] in pwg.get_fieldtable(rfiles[cosmoind]).fieldnames:
            contour(rfiles, varlist[i], cosmoind, xlim, ylim)
        elif varlist[i] in pwg.get_fieldtable(pfiles[cosmoind]).fieldnames:
            contour(pfiles, varlist[i], cosmoind, xlim, ylim)
        else:
            raise Exception('Variable' + varlist[i] + 'not available!')
    
    # Plot trajectories
    
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        lonmat = rootgrp.variables['longitude'][:trjind, :]
        latmat = rootgrp.variables['latitude'][:trjind, :]
        pmat = rootgrp.variables['P'][:trjind, :]
    
        lonmat[:, :] += (180 - pollon)   # Convert to real coordinates
        latmat[:, :] += (90 - pollat)
        
        for j in idlist[i]:
            # Filter out zero values!
            
            parray = pmat[:, j][pmat[:, j] != 0]
            lonarray = lonmat[:, j][pmat[:, j] != 0]
            latarray = latmat[:, j][pmat[:, j] != 0]
            
            single_trj(lonarray, latarray, parray)
    
    # Set plot properties
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)

    cb = fig.colorbar(lc, shrink = 0.7)
    cb.set_label('p')
    cb.ax.invert_yaxis()
    plt.tight_layout()
    
    # Save Plot
    if savename != False:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400)
        plt.close('all')
        plt.clf()

def draw_trj_dot(obj, varlist, loclist, idlist, tplus, 
                 savename = None):
    """
    tplus = time after MODEL start
    """
    
    # Setting up indices
    cosmoind = tplus / obj.dcosmo 
    
    # Set up figure
    fig = plt.figure(figsize = (12,8))
    ax = plt.gca()   
    ax.set_aspect('equal')
    basemap(obj.cfile, obj.xlim, obj.ylim)
    
    
    # Plotting all contour fields
    for i in range(len(varlist)):   # Plotting all contour fields
        if varlist[i] == 'CUM_PREC':
            contour(obj.rfiles, varlist[i], cosmoind, obj.xlim, obj.ylim, 
                    trjstart = trjstart)
        elif varlist[i] in pwg.get_fieldtable(obj.rfiles[cosmoind]).fieldnames:
            contour(obj.rfiles, varlist[i], cosmoind, obj.xlim, obj.ylim)
        elif varlist[i] in pwg.get_fieldtable(obj.pfiles[cosmoind]).fieldnames:
            contour(obj.pfiles, varlist[i], cosmoind, obj.xlim, obj.ylim)
        else:
            raise Exception('Variable' + varlist[i] + 'not available!')
    
    # Plot trajectories
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        trjstart = int(rootgrp.variables['time'][0] / 60)
        trjind = (tplus - trjstart) / obj.dtrj
        lonarray = rootgrp.variables['longitude'][trjind, idlist[i]]
        latarray = rootgrp.variables['latitude'][trjind, idlist[i]]
        parray = rootgrp.variables['P'][trjind, idlist[i]]
    
        lonarray += (180 - obj.pollon)   # Convert to real coordinates
        latarray += (90 - obj.pollat)
        
        plt.scatter(lonarray, latarray, c = parray, s = 10,
                    cmap = plt.get_cmap('Spectral'), linewidth = 0.1,
                    norm=plt.Normalize(100, 1000))
        
    # Set plot properties
    plt.xlim(obj.xlim)
    plt.ylim(obj.ylim)
    plt.title(obj.date + timedelta(minutes = tplus))
        
    # Save Plot
    if savename != False:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400)
        plt.close('all')
        plt.clf()

def draw_asc_loc(obj, lon, lat, varlist, tplot, idtext = '', savename = None):
    """
    TODO
    """

    
    # Plot contours
    draw_contour(obj, varlist, tplot, idtext = idtext, savename = None)
    
    # Convert to Real coords
    lon += (180 - obj.pollon)   
    lat += (90 - obj.pollat)

    print lon, lat
    
    # Plot trajectory dots
    plt.scatter(lon, lat)
    
    # Save Plot
    if savename != False:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400)
        plt.close('all')
        plt.clf()


#########################################
### Back End Plotting Functions
#########################################


def basemap(cfile, xlim, ylim):
    """
    Draws Land-Sea Map
    
    Parameters
    ----------
    cfile : string
      COSMO constants file location
    xlim : tuple
      Dimensions in x-direction in rotated coordinates
    ylim : tuple
      Dimensions in y-direction in rotated coordinates

    """
    dx = 0.025   # Assumes COSMO 2.2 resolution
    field = pwg.getfield(cfile, "FR_LAND_S", invfile = False)
    
    # Setting up grid
    ny, nx = field.shape
    x = np.linspace(xlim[0], xlim[1], nx)
    y = np.linspace(ylim[0], ylim[1], ny)

    X, Y = np.meshgrid(x, y)
    
    # Plot and filter field
    field[field != 0] = 1.
    plt.contourf(X, Y, field, levels = [0.5, 2], colors = ["0.85"], zorder=0.4)
    #plt.contour(X, Y, field, 2, colors = "0.5")
    
    del field
    
    
def contour(filelist, variable, cosmoind, xlim, ylim, trjstart = None):
    """
    Draws contour plot of one variable.
    
    Parameters
    ----------
    filelist : list
      List of COSMO output files
    variable : string
      COSMO identifier of variable
    cosmoind : int
      index in filelist
    xlim : tuple
      Dimensions in x-direction in rotated coordinates
    ylim : tuple
      Dimensions in y-direction in rotated coordinates
      
    """
    print 'Plotting:', variable
    dt = 5. * 60
    dx = 0.025
    
    # Retrieve or evaluate fields
    
    # NOTE: 'CUM_PREC' note implemented right now!
    #if variable == 'CUM_PREC':
        #field1 = pwg.getfield(filelist[trjstart/5], "TOT_PREC_S", 
                              #invfile = False)
        #field2 = pwg.getfield(filelist[-1], "TOT_PREC_S", 
                              #invfile = False)
        #diff = len(filelist)-1 - trjstart/5
        #field = (field2 - field1)/(dt*diff)*3600
    
    # Get precipitation difference
    if variable == "TOT_PREC_S":
        try:
            field1 = pwg.getfield(filelist[cosmoind - 1], "TOT_PREC_S", 
                                invfile = False)
            field2 = pwg.getfield(filelist[cosmoind + 1], "TOT_PREC_S", 
                                invfile = False)
            field = (field2 - field1)/(dt*2)*3600
        except IndexError:
            # Reached end of array, use precip from previous time step
            field1 = pwg.getfield(filelist[cosmoind - 2], "TOT_PREC_S", 
                                invfile = False)
            field2 = pwg.getfield(filelist[cosmoind], "TOT_PREC_S", 
                                invfile = False)
            field = (field2 - field1)/(dt*2)*3600
    # Retrieve regular fields
    else:
        field = pwg.getfield(filelist[cosmoind], variable, invfile = False)
     
    # Setting up grid
    ny, nx = field.shape
    x = np.linspace(xlim[0], xlim[1], nx)
    y = np.linspace(ylim[0], ylim[1], ny)
    X, Y = np.meshgrid(x, y)
    
    # Plot fields with special format
    if variable == "FI":   # Geopotential field
        field = smoothfield(field, 8)
        levels = list(np.arange(400,600,8)) 
        plt.contour(X, Y, field/100, levels = levels, colors = "k", 
                    linewidths = 2)
    elif variable == "T":   # Temperature field
        field = smoothfield(field, 8)
        plt.contourf(X, Y, field, alpha = 0.5)
        plt.colorbar()
        # plt.contour(X, Y, field, levels = list(np.arange(150, 350, 4)), 
                     # colors = "grey", linewidths = 2)
    elif variable in ["TOT_PREC_S", 'CUM_PREC']:   # Precipitation fields
        cmPrec =( (0    , 0.627 , 1    ),
                  (0.137, 0.235 , 0.98 ),
                  (0.1  , 0.1   , 0.784),
                  (0.392, 0     , 0.627),
                  (0.784, 0     , 0.627),
                  (1    , 0.3   , 0.9  ) )   # Tobi's colormap
        levels = [0.1, 0.3, 1, 3, 10, 100]
        plt.contourf(X, Y, field, levels, colors=cmPrec, extend='max', 
                     alpha = 0.8, zorder = 1)
        cbar = plt.colorbar(shrink = 0.7)
        cbar.set_label('Precipitation [cm/h]', rotation = 90)
    elif variable == "PMSL":   # Surface pressure
        field = smoothfield(field, 8)/100
        levels = list(np.arange(900, 1100, 5))
        CS = plt.contour(X, Y, field, levels = levels, colors = 'k', 
                         linewidths = 1, zorder = 0.5, alpha = 0.5)
        plt.clabel(CS, fontsize = 7, inline = 1, fmt = "%.0f")
    elif variable == 'var145_S':   # CAPE
        field = smoothfield(field, 8)
        levels = list(np.arange(100, 3000, 100))
        plt.contourf(X, Y, field, cmap = plt.get_cmap('hot_r'), 
                     extend = 'max', levels = levels, alpha = 0.8, zorder = 0.8)
        cbar = plt.colorbar(shrink = 0.7)
        cbar.set_label('CAPE [J/kg]', rotation = 90)
        
    else:   # All other fields, unformatted
        plt.contour(X, Y, field)
        #plt.colorbar()
    del field
    
    
def single_trj(lonarray, latarray, parray, linewidth = 0.7):
    """
    Plots XY Plot of one trajectory, with color as a function of p
    Helper Function for DrawXYTraj
    """
    global lc
    x = lonarray
    y = latarray
    p = parray

    points = np.array([x,y]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
    lc = col.LineCollection(segments, cmap=plt.get_cmap('Spectral'), 
                            norm=plt.Normalize(100, 1000), alpha = 0.5)
    lc.set_array(p)
    lc.set_linewidth(linewidth)
    plt.gca().add_collection(lc)
    
    
def smoothfield(field, sigma = 8):
    """
    Smoothes given field
    """
    newfield = ndi.filters.gaussian_filter(field, sigma)
    return newfield
        
    
    
    
    
    
    
    
