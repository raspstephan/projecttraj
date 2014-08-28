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

def draw_scatter(array1, array2, idtext, xlabel, ylabel, savename = None):
    """
    TODO
    """
    fig = plt.figure(figsize = (10, 10))
    plt.scatter(array1, array2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    mx = max(np.amax(array1), np.amax(array2))
    plt.xlim(0, mx)
    plt.ylim(0, mx)
    plt.plot([0, mx], [0, mx])
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    if savename != None:
        plt.savefig(savename)
        plt.close('all')



def draw_hist(array, savename = None, xlim = None):
    """
    Plots a histogram of given array.
    
    Parameters
    ----------
    array : np.array
      Array to be used for histogram
    savename : string
      Full path of wanted save names
      
    """
    
    fig = plt.figure()
    
    # Remove nan values, they cause an error in histogram!
    plt.hist(array[np.isfinite(array)], bins = 144, range = xlim)  
    plt.ylabel("Number of trajectories")
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename)
        plt.close('all')
        plt.clf()


def draw_contour(varlist, time, cfile, rfiles, pfiles, savename = False, 
                 pollon = None, pollat = None, xlim = None, ylim = None, 
                 trjstart = None):
    """
    Plots one contour plot. TEST FOR NOW!
    """
    
    # Getting index for COSMO files
    outint = 5   # COSMO output interval
    cosmoind = time / outint
    
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
       
     
    # Set plot properties
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
        
    if savename != False:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400)
        plt.close('all')
        plt.clf()


def draw_trj(varlist, filelist, idlist, cfile, rfiles, pfiles, 
              savename = False,pollon = None, pollat = None, xlim = None, 
              ylim = None, onlybool = False, startarray = None, 
              stoparray = None, trjstart = None):
    """
    Plots one xy plot with trajectories.
    
    Parameters
    ----------
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
         
    # Getting index for COSMO files
    outint = 5   # COSMO output interval
    tmpt = filelist[0].split('/')[-1]
    tmpt = int(tmpt.split('_')[1].lstrip('t'))
    cosmoind = int(tmpt/outint)
    
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
    
    #Get field
    if variable == 'CUM_PREC':
        field1 = pwg.getfield(filelist[trjstart/5], "TOT_PREC_S", 
                              invfile = False)
        field2 = pwg.getfield(filelist[-1], "TOT_PREC_S", 
                              invfile = False)
        diff = len(filelist)-1 - trjstart/5
        field = (field2 - field1)/(dt*diff)*3600

    elif variable == "TOT_PREC_S":
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
    else:
        field = pwg.getfield(filelist[cosmoind], variable, invfile = False)
     
    # Setting up grid
    ny, nx = field.shape
    x = np.linspace(xlim[0], xlim[1], nx)
    y = np.linspace(ylim[0], ylim[1], ny)

    X, Y = np.meshgrid(x, y)
    
    
    if variable == "FI":
        field = smoothfield(field, 8)
        levels = list(np.arange(400,600,8))   # Needs smoothing?
        plt.contour(X, Y, field/100, levels = levels, colors = "k", 
                    linewidths = 2)
    elif variable == "T":
        field = smoothfield(field, 8)
        plt.contourf(X, Y, field, alpha = 0.5)
        plt.colorbar()
        # plt.contour(X, Y, field, levels = list(np.arange(150, 350, 4)), 
                     # colors = "grey", linewidths = 2)
    elif variable in ["TOT_PREC_S", 'CUM_PREC']:
        cmPrec =( (0    , 0.627 , 1    ),
                  (0.137, 0.235 , 0.98 ),
                  (0.1  , 0.1   , 0.784),
                  (0.392, 0     , 0.627),
                  (0.784, 0     , 0.627),
                  (1    , 0.3   , 0.9  ) )   # Tobi's colormap
        levels = [0.1, 0.3, 1, 3, 10, 100]
        plt.contourf(X, Y, field, levels, colors=cmPrec, extend='max', 
                     alpha = 0.8, zorder = 1)
        plt.colorbar(shrink = 0.7)
    elif variable == "PMSL":
        field = smoothfield(field, 8)/100
        levels = list(np.arange(900, 1100, 5))
        CS = plt.contour(X, Y, field, levels = levels, colors = 'k', 
                         linewidths = 1, zorder = 0.5, alpha = 0.5)
        plt.clabel(CS, fontsize = 7, inline = 1, fmt = "%.0f")
    else:
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
        
    
    
    
    
    
    
    
