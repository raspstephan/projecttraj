"""
plots module
------------

This module contains all plotting functions.

"""


import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.collections as col
from mpl_toolkits.basemap import Basemap
import cosmo_utils.pywgrib as pwg
import scipy.ndimage as ndi
from datetime import datetime, timedelta
import netCDF4 as nc


def draw_hist(array, savename = None):
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
    plt.hist(array[np.isfinite(array)], bins = 150)  
    
    if savename != None:
        plt.savefig(savename)
        plt.xlabel('hmmm???')
        plt.ylabel("Number of trajectories")
        plt.close('all')
        plt.clf()


def draw_xy(varlist, filelist, idlist, cfile, rfiles, pfiles, savename = False, 
            pollon = None, pollat = None, xlim = None, ylim = None):
    """
    Plots one xy plot.
    
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
        try:
            contour(pfiles, varlist[i], cosmoind, xlim, ylim)
        except:
            pass
        else:
            contour(rfiles, varlist[i], cosmoind, xlim, ylim)
    
    # Plot trajectories
    for i in range(len(filelist)):
        rootgrp = nc.Dataset(filelist[i], 'r')
        lonmat = rootgrp.variables['longitude'][:, :]
        latmat = rootgrp.variables['latitude'][:, :]
        pmat = rootgrp.variables['P'][:, :]
    
        lonmat[:, :] += (180 - pollon)   # Convert to real coordinates
        latmat[:, :] += (90 - pollat)
        
        for j in idlist[i]:
            # Filter out zero values!
            parray = pmat[:, j][pmat[:, j] != 0]
            lonarray = lonmat[:, j][pmat[:, j] != 0]
            latarray = latmat[:, j][pmat[:, j] != 0]
            
            single_traj(lonarray, latarray, parray)
    
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
    plt.contourf(X, Y, field, levels = [0.5, 2], colors = ["0.85"])
    #plt.contour(X, Y, field, 2, colors = "0.5")
    
    del field
    
    
def contour(filelist, variable, cosmoind, xlim, ylim):
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
    dt = 5. * 60
    dx = 0.025
    
    #Get field
    if variable == "TOT_PREC_S":
        field1 = pwg.getfield(filelist[cosmoind - 1], "TOT_PREC_S", 
                              invfile = False)
        field2 = pwg.getfield(filelist[cosmoind + 1], "TOT_PREC_S", 
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
    elif variable == "TOT_PREC_S":
        cmPrec =( (0    , 0.627 , 1    ),
                  (0.137, 0.235 , 0.98 ),
                  (0.1  , 0.1   , 0.784),
                  (0.392, 0     , 0.627),
                  (0.784, 0     , 0.627),
                  (1    , 0.3   , 0.9  ) )   # Tobi's colormap
        levels = [0.1, 0.3, 1, 3, 10, 100]
        plt.contourf(X, Y, field, levels, colors=cmPrec, extend='max', 
                     alpha = 0.8, zorder = 10)
        plt.colorbar(shrink = 0.7)
    elif variable == "PMSL":
        field = smoothfield(field, 8)/100
        levels = list(np.arange(900, 1100, 5))
        CS = plt.contour(X, Y, field, levels = levels, colors = 'k', 
                         linewidths = 1, zorder = 9, alpha = 0.5)
        plt.clabel(CS, fontsize = 7, inline = 1, fmt = "%.0f")
    else:
        plt.contour(X, Y, field)
        #plt.colorbar()
    del field
    
    
def single_traj(lonarray, latarray, parray, linewidth = 0.7):
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
                            norm=plt.Normalize(100, 1000), alpha = 0.8)
    lc.set_array(p)
    lc.set_linewidth(linewidth)
    plt.gca().add_collection(lc)
    
    
def smoothfield(field, sigma = 8):
    """
    Smoothes given field
    """
    newfield = ndi.filters.gaussian_filter(field, sigma)
    return newfield
        
    
    
    
    
    
    
    
