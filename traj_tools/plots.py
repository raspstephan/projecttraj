"""
plots module
------------

This module contains all plotting functions.
To make movie use from command line.

ffmpeg -r 2 -pattern_type glob -i '*.png' -c:v libx264 movie.mkv

mencoder mf://*.png -mf w=1000:fps=2:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output_fast.avi
"""


import numpy as np
import matplotlib
import os
# taken from cosmo_utils, does not work
try:
    os.environ["DISPLAY"]
    matplotlib.use('Qt4Agg')
    print "X-Server detected, using Qt4Agg backend for plotting"
except KeyError:
    if matplotlib.get_backend() in matplotlib.rcsetup.interactive_bk:
        matplotlib.use("Agg") 
        # interface "Agg" can plot without x-server connection
    print "No X-Server detected, using agg backend for plotting"
# matplotlib.use('Qt4Agg')   # Use at your own risk
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.collections as col
from mpl_toolkits.basemap import Basemap
import cosmo_utils.pywgrib as pwg
import scipy.ndimage as ndi
from scipy.stats import nanmean
from datetime import datetime, timedelta
import netCDF4 as nc




def draw_vs_t(obj, tracer, loclist, idlist, savename = None, sigma = None):
    """
    TODO
    """
    
    
    # 
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        
    #if dataname == 'z_CD':
        #array = rootgrp.variables['z'][:, fileid]
        
        #array = np.gradient(array)
        #if sigma != None:
            #array = ndi.filters.gaussian_filter(array, sigma)
    #else:
        #array = rootgrp.variables[dataname][:, fileid]
        #if sigma != None:
            #array = ndi.filters.gaussian_filter(array, sigma)
            
        tracerarray = rootgrp.variables[tracer][:, :]
        tarray = rootgrp.variables['time'][:] / 60   # Convert to minutes
        
        # Convert zeros to nans
        tracerarray[tracerarray == 0] = np.nan
        
        # Plot trajectories in idlist
        plt.plot(tarray, tracerarray[:, idlist[i]])
            
            
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')


def draw_vs_p(obj, tracer, loclist, idlist, startarray, stoparray, xlim, 
              savename = None, sigma = 1, idtext = '', ylim = None, 
              binwidth = 5.):
    """
    TODO
    """
    counter = 0
    
    # Initialize p array
    parray = np.arange(0, 1050., binwidth) 
    nbins = parray.shape[0]
    tracerlist = [[] for i in range(nbins)]     
    
    for i in range(len(loclist)):
        
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        tracemat = rootgrp.variables[tracer][:, :]
        if tracer in ['POT_VORTIC', 'var4']:
            tracemat = tracemat * 1.e6
        pmat = rootgrp.variables['P'][:, :]
        
        # Get P data and set zeros to nan
        zmask = np.ma.mask_or(pmat == 0, np.isnan(pmat))
        tracemat[zmask] = np.nan
        tracemat[tracemat == 0] = np.nan
        
        # Convert pmat to indices
        pindmat = np.around((pmat - xlim[0]) / binwidth)
        
        for j in idlist[i]:
            
            tracearray = tracemat[startarray[counter]:stoparray[counter], j]
            pindarray = pindmat[startarray[counter]:stoparray[counter], j]
            
            for k in range(tracearray.shape[0]):
                #print pindarray[k], nbins
                tracerlist[int(pindarray[k])].append(tracearray[k])
                
            counter += 1
    
    countlist = []
    meanlist = []
    per5 = []
    per25 = []
    per50 = []
    per75 = []
    per95 = []

    # Set lower density limit TODO
    denslim = 100

    for i in range(nbins):
        tmparray = np.array(tracerlist[i])
        tmparray = tmparray[np.isfinite(tmparray)]
        countlist.append(tmparray.shape[0])
        if tmparray.shape[0] >= denslim:
            meanlist.append(np.mean(tmparray))
            per5.append(np.percentile(tmparray, 5))
            per25.append(np.percentile(tmparray, 25))
            per50.append(np.percentile(tmparray, 50))
            per75.append(np.percentile(tmparray, 75))
            per95.append(np.percentile(tmparray, 95))
        else: 
            meanlist.append(np.nan)
            per5.append(np.nan)
            per25.append(np.nan)
            per50.append(np.nan)
            per75.append(np.nan)
            per95.append(np.nan)
    
    fig = plt.figure(figsize = (10, 8))
    ax = plt.gca()
    ax.grid(color = 'dimgrey', linestyle = '-')
    ax.set_frame_on(False)
    plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    plt.plot(parray, meanlist)
    ax.fill_between(parray, per5, per95, 
                        facecolor = 'lightgrey', edgecolor = 'lightgrey',
                        label = '90%')
    ax.fill_between(parray, per25, per75, 
                        facecolor = 'darkgrey', edgecolor = 'darkgrey',
                        label = '50%')  
    l2, = plt.plot(parray, per50, 'black')
    l1, = plt.plot(parray, meanlist, 'firebrick', linewidth = 2)
    # Get filled colors as legends
    r1 = plt.Rectangle((0, 0), 1, 1, fc="lightgrey")
    r2 = plt.Rectangle((0, 0), 1, 1, fc="darkgrey")
    plt.legend([l1, l2, r1, r2], ['mean', 'median', '50%', '90%'], loc = 2)
    
    # Plot second axis
    ax2 = ax.twinx()
    ax2.bar(parray, countlist, linewidth = 0, color = 'sage', width = 5, alpha = 0.8)
    ax2.set_yticks(np.arange(0, np.max(countlist) + 5000, 5000))
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax2.set_xlim(xlim)
    ax2.set_ylim((0, np.max(countlist) * 4))
    ax.invert_xaxis()

    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight', dpi = 300)
        plt.close('all')


def draw_centered_vs_t(obj, loclist, idlist, tracer, carray, savename = None,
                       plottype = 'Smooth', idtext = '', ylim = None, 
                       xlim = None, sigma = 1):
    """
    Draws evolution of a tracer of all trajectories given by filter,
    centered around midpoint of ascent, as given by carray.
    Now also allows tracer "CD_w" = Centered Difference vertical velocity.
    
    Parameters
    ----------
    obj : TrjObj object
      self 
    loclist : list
      List of trjfiles
    idlist : list of lists
      Indices for each trjfile
    tracer : string
      NetCDF name of tracer
    carray : string
      Ascent array to be used for centering
    savename : string
      Name of file to be saved
    plottype : string
      Type of plot. Only 'Smooth' is up to date!
    idtext : string
      Id string to be displayed
    ylim : tuple, list
      Tuple or list of y-axis limits
    xlim : tuple, list
      Tuple or list of x-axis limits
    sigma : float
      Sigma value for smoothing

    NOTE: Can use a lot of RAM!
    """
    istart = 0
    matlist = []

    # Round carray to closes value divisable by dtrj
    carray = obj.dtrj * np.around(carray / obj.dtrj)
    exttarray = np.arange(-obj.maxmins, obj.maxmins + obj.dtrj, obj.dtrj,
                          dtype = 'float')
    
    for i in range(len(loclist)):
        istop = len(idlist[i]) + istart
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        tarray = rootgrp.variables['time'][:] / 60   # Convert to minutes
        
        # Create extended matrix and tarray
        exttracemat = np.empty((exttarray.shape[0], len(idlist[i])))
        exttracemat.fill(np.nan)
        
        
        # Get data and relative times
        if tracer == 'CD_w':
            zmat = rootgrp.variables['z'][:, idlist[i]]
            if zmat.shape[1] == 1:
                tracemat = zmat
                grad = np.gradient(zmat[:, 0])
                tracemat[:, 0] = grad
            else:
                tracemat = np.gradient(zmat)[0]
                #print tracemat
            tracemat = tracemat / 60. / obj.dtrj   # m/s
            #print tracemat
        else:
            tracemat = rootgrp.variables[tracer][:, idlist[i]]
        
        # Get P data and set zeros to nan
        pmat = rootgrp.variables['P'][:, idlist[i]]
        zmask = np.ma.mask_or(pmat == 0, np.isnan(pmat))
        tracemat[zmask] = np.nan
        tracemat[tracemat == 0] = np.nan
        
        if tracer in ['var4', 'POT_VORTIC']:
            tracemat = tracemat * 1.e6   # PVU
        tmat = np.array([tarray] * len(idlist[i])).transpose()
        reltmat = tmat - carray[istart:istop]
        istart = istop
        
        for j in range(len(idlist[i])):
            ind = np.where(exttarray == reltmat[0, j])[0][0]
            exttracemat[ind:(ind + tracemat[:, j].shape[0]), j] = tracemat[:, j]
                    
        matlist.append(exttracemat)
    
    totmat = np.hstack(matlist)    
    meanarray = np.nanmean(totmat, axis = 1)
    
    mask = np.isfinite(meanarray)
    
    # Convert time array to hours
    exttarray = exttarray / 60.
    
    # Set up figure
    fig = plt.figure(figsize = (10, 8))
    ax = plt.gca()
    
    # Set plot properties
    if ylim == None:
        ylim = (1.5 * np.nanmin(meanarray), 1.5 * np.nanmax(meanarray))
    if xlim == None:
        xlim = [exttarray[mask].min(), exttarray[mask].max()]
    
    plt.plot([0, 0], ylim, color = 'dimgrey', linewidth = 2)
    plt.plot(xlim, [0, 0], color = 'dimgrey', linewidth = 2)
    
    
    if plottype == 'All':
        plt.plot(exttarray, totmat, 'grey')
        plt.plot(exttarray[mask], meanarray[mask], 'r')
    elif plottype == 'Std':
        stdarray = np.nanstd(totmat, axis = 1)
        per5 = nanpercentile(totmat, 5)
        per95 = nanpercentile(totmat, 95)
        plt.plot(exttarray[mask], (meanarray - stdarray)[mask], 'grey')
        plt.plot(exttarray[mask], (meanarray + stdarray)[mask], 'grey')
        plt.plot(exttarray[mask], meanarray[mask], 'r')
    elif plottype == 'Range':
        maxarray = np.nanmax(totmat, axis = 1)
        minarray = np.nanmin(totmat, axis = 1)
        plt.plot(exttarray[mask], maxarray[mask], 'grey')
        plt.plot(exttarray[mask], minarray[mask], 'grey')
        plt.plot(exttarray[mask], meanarray[mask], 'r')
    elif plottype == 'Smooth':
        per5 = nanpercentile(totmat, 5)
        per95 = nanpercentile(totmat, 95)
        smoothper5 = ndi.filters.gaussian_filter(per5[mask], sigma)
        smoothper95 = ndi.filters.gaussian_filter(per95[mask], sigma)
        plt.plot(exttarray[mask], smoothper5, 'lightgrey')
        plt.plot(exttarray[mask], smoothper95, 'lightgrey')
        per25 = nanpercentile(totmat, 25)
        per75 = nanpercentile(totmat, 75)
        smoothper25 = ndi.filters.gaussian_filter(per25[mask], sigma)
        smoothper75 = ndi.filters.gaussian_filter(per75[mask], sigma)
        plt.plot(exttarray[mask], smoothper25, 'darkgrey')
        plt.plot(exttarray[mask], smoothper75, 'darkgrey')
        smootharray = ndi.filters.gaussian_filter(meanarray[mask], sigma)
        plt.plot(exttarray[mask], smootharray, 'r')
    elif plottype == 'Fill':
        per5 = nanpercentile(totmat, 5)
        per95 = nanpercentile(totmat, 95)
        smoothper5 = ndi.filters.gaussian_filter(per5[mask], sigma)
        smoothper95 = ndi.filters.gaussian_filter(per95[mask], sigma)
        ax.fill_between(exttarray[mask], smoothper5, smoothper95, 
                        facecolor = 'lightgrey', edgecolor = 'lightgrey',
                        label = '90%')
        
        per25 = nanpercentile(totmat, 25)
        per75 = nanpercentile(totmat, 75)
        smoothper25 = ndi.filters.gaussian_filter(per25[mask], sigma)
        smoothper75 = ndi.filters.gaussian_filter(per75[mask], sigma)
        ax.fill_between(exttarray[mask], smoothper25, smoothper75, 
                        facecolor = 'darkgrey', edgecolor = 'darkgrey',
                        label = '50%')
        
        per50 = nanpercentile(totmat, 50)
        smoothper50 = ndi.filters.gaussian_filter(per50[mask], sigma)
        l2, = plt.plot(exttarray[mask], smoothper50, 'black')
        
        smootharray = ndi.filters.gaussian_filter(meanarray[mask], sigma)
        l1, = plt.plot(exttarray[mask], smootharray, 'firebrick', linewidth = 2)
        # Get filled colors as legends
        r1 = plt.Rectangle((0, 0), 1, 1, fc="lightgrey")
        r2 = plt.Rectangle((0, 0), 1, 1, fc="darkgrey")
        plt.legend([l1, l2, r1, r2], ['mean', 'median', '50%', '90%'])
        
    del totmat
    
    ax.set_ylim(ylim)
    if tracer == 'P':
        ax.invert_yaxis()
    if (xlim[1] - xlim[0]) < 24:
        dtick = 1
    elif (xlim[1] - xlim[0]) < 48:
        dtick = 6
    else:
        dtick = 12
    ax.xaxis.set_ticks(np.arange(-120, 120, dtick))
    ax.set_xlim(xlim)
    ax.grid(color = 'dimgrey', linestyle = '-')
    ax.set_frame_on(False)
    plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    if tracer == 'var4':
        tracer == 'PVU'
    plt.ylabel(tracer)
    plt.xlabel('Time [hrs] relative to center') 
    
    
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight', dpi = 300)
        plt.close('all')
    
def nanpercentile(a, per):
    """
    Percentile function ignoring nans
    
    Parameters
    ----------
    a : np.array
      Input array
    per : float
      Percentile value
    
    Returns
    -------
    out : np.array
      Output values
      
    """
    a = np.array(a)
    if len(a.shape) == 1:
        out = np.percentile(a[np.isfinite(a)], per)
    out = np.empty(a.shape[0])
    for i in range(a.shape[0]):
        b = a[i][np.isfinite(a[i])]
        if b.shape[0] == 0:
            out[i] = np.nan
        else:
            out[i] = np.percentile(b[np.isfinite(b)], per)
    return out


def draw_hist_2d(obj, varname1, varname2, loclist, idlist, carray, dplus):
    """
    Draws a 2D histogram of two variables in filter.
    'plus' mins after 'after' ascent.
    
    Parameters
    ----------
    obj : TrjObj object
      self
    varname1 : string
      x-axis variable
    varname2 : string
      y-axis variable
    loclist : list
      List of trjfiles
    idlist : list of lists
      Indices for each trjfile
    carray : string
      Name of ascent after which to be plotted
    dplus : int
      Time [mins] after ascent to be plotted
    
    TODO
    - only round carray when necessary!
    """
    
    varlist1 = np.array([])
    varlist2 = np.array([])
    acosmobytrj = obj.dacosmo / obj.dtrj
    carray = acosmobytrj * np.around(carray / acosmobytrj)   # Round to nearest 60
    
    cnt = 0   # continuous counter for startarray and stoparray
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        
        var1 = rootgrp.variables[varname1][:, :]
        var2 = rootgrp.variables[varname2][:, :]
        p = rootgrp.variables['P'][:, :]
        
        for j in idlist[i]:
            #mask = (np.isfinite(p[:, j])) & (p[:, j] != 0)
            #mask &= (np.isfinite(var1[:, j])) & (np.isfinite(var2[:, j]))
            #if not stoparray == None:
                #mask[:stoparray[cnt]] = False
            #varlist1 = np.append(varlist1, var1[:, j][mask])
            #varlist2 = np.append(varlist2, var2[:, j][mask])
            ind = carray[cnt] + dplus
            if ind < var1[:, j].shape[0]:
                if (np.isfinite(p[ind, j])) & (p[ind, j] != 0):
                    varlist1 = np.append(varlist1, var1[ind, j])
                    varlist2 = np.append(varlist2, var2[ind, j])
            cnt += 1
    
    varlist2 *= 1.e6
    #plt.scatter(varlist1, varlist2)
    x_edges = np.arange(0, 1000 + 20, 20)
    y_edges = np.arange(-25, 25 + 1., 1.)
    plt.hist2d(varlist1, varlist2, cmin = 1, 
               norm = clr.LogNorm(1, 1000),
               bins = [x_edges, y_edges], zorder = 10)
    plt.colorbar(use_gridspec = True)
    plt.tight_layout()
    plt.gcf().subplots_adjust(bottom = 0.1, left = 0.1)
    #plt.gca().set_ylim(-10.e-6, 10.e-6)


def draw_scatter_3(obj, varname, loclist, idlist, carray, dplus):
    """
    TODO
    - only round carray when necessary!
    """
    draw_contour(obj, [], 1440, idtext = '')
    
    varlist = np.array([])
    lonlist = np.array([])
    latlist = np.array([])
    acosmobytrj = obj.dacosmo / obj.dtrj
    carray = acosmobytrj * np.around(carray / acosmobytrj)   # Round to nearest 60
    
    cnt = 0   # continuous counter for startarray and stoparray
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        
        var = rootgrp.variables[varname][:, :]
        lon = rootgrp.variables['longitude'][:, :]
        lat = rootgrp.variables['latitude'][:, :]
        p = rootgrp.variables['P'][:, :]
        
        for j in idlist[i]:
            #mask = (np.isfinite(p[:, j])) & (p[:, j] != 0)
            #mask &= (np.isfinite(var1[:, j])) & (np.isfinite(var2[:, j]))
            #if not stoparray == None:
                #mask[:stoparray[cnt]] = False
            #varlist1 = np.append(varlist1, var1[:, j][mask])
            #varlist2 = np.append(varlist2, var2[:, j][mask])
            if cnt % 3 == 0: 
                ind = carray[cnt] + dplus
                if ind < var[:, j].shape[0]:
                    if (np.isfinite(p[ind, j])) & (p[ind, j] != 0):
                        varlist = np.append(varlist, var[ind, j])
                        lonlist = np.append(lonlist, lon[ind, j])
                        latlist = np.append(latlist, lat[ind, j])
            cnt += 1
    
    varlist *= 1.e6
    plt.scatter(lonlist, latlist, c = varlist, s = 10, 
                norm = plt.Normalize(-15, 15), linewidth = 0.1)
    plt.colorbar(use_gridspec = True)
    plt.tight_layout()

    

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
        carray = 'blue'   # Set to default
        
    # Convert to hours
    array1 = array1 / 60
    #array2 = array2 / 60
    #plt.scatter(array1, array2)
    # Plot scatter plot, add labels
    sca = ax.scatter(array1, array2) #c = carray, s = 8,
                     #cmap=plt.get_cmap('Spectral'), 
                     #norm=plt.Normalize(100, 1000), linewidths = 0)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Set limits, plot diagnonal line, set ticks
    xmx = np.amax(array1[np.isfinite(array1)])
    ymx = np.amax(array2[np.isfinite(array2)])
    mx = max(xmx, ymx)
    plt.xlim(0, xmx+0.1*xmx)
    plt.ylim(0, ymx+0.1*ymx)
    #plt.plot([0, mx], [0, mx])
    #plt.xticks(np.arange(0, xmx+3, 3))
    #plt.yticks(np.arange(0, ymx+3, 3))
    
    # Add idtext
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    ##cb = fig.colorbar(sca, shrink = 0.7)
    ##cb.set_label('p')
    ##cb.ax.invert_yaxis()
    ##ax.set_axis_bgcolor('grey')
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')


def draw_avg(dataname, loclist, idlist, idtext = '', centdiff = False, 
             savename = None):
    """
    Draws average of data over all given trajectories.
    
    Parameters
    ---------
    dataname : string
      Tracer to be averaged
    loclist : list
      List of unique file locations
    idlist : list
      List of list of trajectory IDs 
    idtext : string
      Text to be displayed on top right of picture
    centdiff : bool
      If true centered difference of tracer will be plotted
    savename : string
      Full path of file to be saved 
    
    """
    fig = plt.figure()
    newmat = []
    
    # Add idtext
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        mat = rootgrp.variables[dataname][:, :]
        time = rootgrp.variables['time'][:] / 60 / 60   # in hrs
        for j in idlist[i]:
            array = mat[:, j]
            # Convert zeros to nan
            array[array == 0] = np.nan
            if centdiff:
                array = np.gradient(array)
                array = ndi.filters.gaussian_filter(array, 5)
                array = array / 5   # Convert to Minutes
            newmat.append(array)
            plt.rc('axes', color_cycle = ['cyan', 'sage', 
                                          'salmon', 'gold', 'violet', 'brown'])
            plt.gca().plot(time, array, linewidth = 0.6, alpha = 0.5)
    newmat = np.array(newmat)
    avg = nanmean(newmat, axis = 0)
    plt.plot(time, avg, 'black', linewidth = 2)
    plt.gca().set_axis_bgcolor('grey')
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')
        


def draw_hist(array, idtext = '', xlabel =  None, savename = None, log = False,
              mintohrs = False, **kwargs):
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
    
    # Convert to np array
    array = np.array(array)
    
    if mintohrs:
        array = array / 60.
    
    # Plot histogram, remove nans
    if log:
        plt.hist(array[np.isfinite(array)], bins = np.logspace(-1, 3, 144), 
                 **kwargs)
        plt.gca().set_xscale('log')
    else:
        plt.hist(array[np.isfinite(array)], bins = 144, **kwargs)
        
    
    # Add labels and text
    plt.title('Total Number of trajectories: ' + str( array.shape[0]))
    plt.ylabel("Number of trajectories")
    plt.xlabel(xlabel)
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()



def draw_mult_hist(arrays, idtext = '', xlabel = '', savename = None, **kwargs):
    """
    TODO
    Does not work!
    """
    
    # Set up figure
    fig = plt.figure()
    print arrays
    # Loop through arrays and plot
    for array in arrays:
        array = array[np.isfinite(array)]
        plt.hist(array, histtype = 'step', bins = 100, normed  = True, 
                 **kwargs)
    
    # Add labels and text
    plt.ylabel("Number of trajectories")
    plt.xlabel(xlabel)
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()
    


def draw_trj_dens(obj, varlist, filelist, idlist, tplot, tracerange = None, 
                 idtext = '', savename = None):
    """
    Draws a 2D histogram of trajectories over a contour map.
        
    Parameters
    ----------
    varlist : list
      List of variables to be plotted. E.g. ["PMSL", "TOT_PREC_S"]
      'CUM_PREC' for cumulative precipitation
    time : float
      Time in minutes after model start
    filelist : list
      List of unique file locations
    idlist : list
      List of list of trajectory IDs 
    tracerange : tuple
      Tuple eg ('P', 700, 1000), only consider trajectories within this 
      range.
    idtext : string
      Text to be displayed in plot
    savename : string
      Name of file to be saved
        
    """
    
    # Draw contour maps
    draw_contour(obj, varlist, tplot, idtext = idtext)
    
    # Get arrays for histogram
    lonlist = []
    latlist = []
    tracelist = []
    for i in range(len(filelist)):
        rootgrp = nc.Dataset(filelist[i], 'r')
        startt = rootgrp.variables['time'][0] / 60   # Convert to minutes
        trjind = int((tplot - startt) / obj.dtrj)
        lonlist.append(rootgrp.variables['longitude'][trjind, idlist[i]])
        latlist.append(rootgrp.variables['latitude'][trjind, idlist[i]])
        if not tracerange == None:
            tracelist.append(rootgrp.variables[tracerange[0]][trjind, 
                                                              idlist[i]])
        rootgrp.close()
    
    if len(lonlist) != 0:    
        lonlist = np.concatenate(lonlist)
        latlist = np.concatenate(latlist)
    
    if not tracerange == None:
        tracelist = np.concatenate(tracelist)
        mask = (tracelist > tracerange[1]) & (tracelist < tracerange[2])
        lonlist = lonlist[mask]
        latlist = latlist[mask]

    # Plot 2D Histogram
    norm = plt.Normalize(0, 500)
    cmap = clr.ListedColormap(['khaki'] + ['yellow'] * 2 + 
                              ['greenyellow'] * 3 + ['darksage'] * 4 + 
                              ['darkolivegreen'])
    dbin = 1   # degrees
    x_edges = np.arange(obj.xlim[0], obj.xlim[1] + dbin, dbin)
    y_edges = np.arange(obj.ylim[0], obj.ylim[1] + dbin, dbin)
    plt.hist2d(lonlist, latlist, alpha = 0.5, zorder = 10, cmin = 1, 
               range = [obj.xlim, obj.ylim], norm = norm, cmap = cmap, 
               bins = [x_edges, y_edges])
    plt.colorbar(shrink = 0.7, extend = 'max')
    plt.tight_layout()
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()


def draw_intersect_hor(obj, tracer, filelist, idlist, level, leveltype = 'P', 
                       idtext = '', savename = None):
    """
    Draws value of tracer at intersection with level over bg map
    
    Parameters
    ----------
    tracer : string
      NetCDF Name
    filelist : list
      List of unique file locations
    idlist : list
      List of list of trajectory IDs 
    level : float
      value of level
    leveltype : string
      Variable name of level
    idtext : string
      id text
    savename : string
      path to output directory
      
    Returns
    -------
    traceplot : array
      values at intersection
    """
    
    # Draw contour maps
    draw_contour(obj, [], 0, idtext = idtext)
    
    # Initialize plotlists
    lonplot = []
    latplot = []
    traceplot = []
    # Go through trajectory files
    for i in range(len(filelist)):
        print 'Calculating intersections for', filelist[i]
        rootgrp = nc.Dataset(filelist[i], 'r')
        p = rootgrp.variables['P'][:, :]
        lons = rootgrp.variables['longitude'][:, :]
        lats = rootgrp.variables['latitude'][:, :]
        z = rootgrp.variables['z'][:, :]
        if tracer == 'CD_w':
            tracemat = np.gradient(z)[0] / obj.dtrj / 60.   # m/s
            positive = True
        else:
            tracemat = rootgrp.variables[tracer][:, :]
            if tracer in ['var4', 'POT_VORTIC']:
                tracemat = tracemat * 1.e6   # PVU
            positive = False
        
        for j in idlist[i]:
            if leveltype == 'P':
                array = p[:, j]
            elif leveltype == 'z':
                array = z[:, j]
            else:
                raise Exception('Wrong leveltype')
            # Filter for zeros and nans
            
            array = array[np.isfinite(tracemat[:, j])]
            tracearray = tracemat[:, j][np.isfinite(tracemat[:, j])]
            lonarray = lons[:, j][np.isfinite(tracemat[:, j])]
            latarray = lats[:, j][np.isfinite(tracemat[:, j])]
            array = array[np.isfinite(array)]
            array = array[array != 0]
            
            mask = np.ma.masked_where(array > level, array)
            slices = np.ma.notmasked_contiguous(mask)

            if (np.any(mask.mask) == False) or (np.all(mask.mask) == True):
                pass
            elif type(slices) != list:
                if slices.start != 0:
                    start = slices.start - 1
                    stop = slices.start
                    traceplot, lonplot, latplot = interp_vert_vel(j,
                            array, tracearray, lonarray, latarray, start, stop, 
                            traceplot, 
                            lonplot, latplot, level, positive)
                if slices.stop != array.shape[0]:
                    start = slices.stop - 1
                    stop = slices.stop
                    traceplot, lonplot, latplot = interp_vert_vel(j,
                            array, tracearray, lonarray, latarray, start, stop, 
                            traceplot, 
                            lonplot, latplot, level, positive)
                    
            else:
                for s in slices:
                    if s.start != 0:
                        start = s.start - 1
                        stop = s.start
                        traceplot, lonplot, latplot = interp_vert_vel(j,
                            array, tracearray, lonarray, latarray, start, stop, 
                            traceplot, 
                            lonplot, latplot, level, positive)
                    if s.stop != array.shape[0]:
                        start = s.stop - 1
                        stop = s.stop
                        traceplot, lonplot, latplot = interp_vert_vel(j,
                            array, tracearray, lonarray, latarray, start, stop, 
                            traceplot, 
                            lonplot, latplot, level, positive)
            
            traceplot = np.array(traceplot) / obj.dtrj / 60 # Convert to seconds
            plt.scatter(lonplot, latplot, c = traceplot,
                        cmap = plt.get_cmap('spectral_r'),
                        norm = clr.LogNorm(0.01, 10),
                        linewidth = 0.1, s = 5)
            plt.title('Velocity at intersection with ' + leveltype + str(level))
            cb = plt.colorbar(extend = 'max')
            cb.set_label('Vertical velocity [m/s]')
                    
    # Save Plot
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()
        
    return traceplot
        
def interp_vert_vel(j, array, tracearray, lonarray, latarray, start, stop, 
                    traceplot, lonplot, latplot, level, positive = True):
    """
    Function for use in draw_intersect_hor
    TODO: docstring
    
    """
    w1 = np.abs(array[stop] - level) / np.abs(array[stop] - array[start])
    w2 = 1. - w1
    intvel = w1 * tracearray[start] + w2 * tracearray[stop]
    if positive:
        if intvel > 0:
            lonplot.append(w1 * lonarray[start] + w2 * lonarray[stop])
            latplot.append(w1 * latarray[start] + w2 * latarray[stop])
            traceplot.append(w1 * tracearray[start] + w2 * tracearray[stop])
    else:
        lonplot.append(w1 * lonarray[start] + w2 * lonarray[stop])
        latplot.append(w1 * latarray[start] + w2 * latarray[stop])
        traceplot.append(w1 * tracearray[start] + w2 * tracearray[stop])
        
    #traceplot.append((z[stop, j] - z[start, j]))
    return traceplot, lonplot, latplot

                
            
def draw_field(obj, field, fieldname, savename = None):
    """
    TODO
    """
    # Setting up figure
    fig = plt.figure(figsize = (12,8))
    ax = plt.gca()   
    ax.set_aspect('equal')
    
    ny, nx = field.shape
    x = np.linspace(obj.xlim[0], obj.xlim[1], nx)
    y = np.linspace(obj.ylim[0], obj.ylim[1], ny)
    X, Y = np.meshgrid(x, y)
    
    if fieldname == 'var4':
        sigma = 0
        field = ndi.filters.gaussian_filter(field * 1.e6, sigma)
        
        cdict2 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 1.0),
                   (1.0, 0.1, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.1),
                   (0.5, 1.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
        clrlist = (3 * ['blue'] + ['lightblue', 'khaki'] +  
                  ['orange', 'salmon', 'red', 'darkred', 'darkgrey', 'black'])
        blue_red2 = clr.LinearSegmentedColormap('BlueRed2', cdict2)
        
        plt.contourf(X, Y, field, levels = range(-8, 16, 2), 
                     colors = clrlist)
        plt.colorbar()
        plt.contour(X, Y, field, levels = [2.], colors = 'darkgreen', 
                    extend = 'max')
    
    
    
    
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
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
    idtext : string
          Text to be displayed in plot
    savename : string
      Full path to file to be saved
      
    """
    
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
                    
        # Get index and filelist
        if type(varlist[i]) in [list, tuple]:
            var = varlist[i][0]
            zlevel = varlist[i][1]
            cosmoind, filelist = obj._get_index(var, time)
            contour(filelist, var, cosmoind, obj.xlim, obj.ylim, zlevel)
        else: 
            cosmoind, filelist = obj._get_index(varlist[i], time)
            contour(filelist, varlist[i], cosmoind, obj.xlim, obj.ylim)
        
    # Set plot properties
    plt.xlim(obj.xlim)
    plt.ylim(obj.ylim)
    
    # Set labels and title
    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.title(obj.date + timedelta(minutes = time))
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    plt.tight_layout()
    
        
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()


def draw_trj(obj, varlist, filelist, idlist, cfile, rfiles, pfiles, 
             savename = False,pollon = None, pollat = None, xlim = None, 
             ylim = None, onlybool = False, startarray = None, 
             stoparray = None, trjstart = None, idtext = '', linewidth = 0.7,
             carray = 'P', centdiff = False, sigma = None):
    """
    Plots one xy plot with trajectories.
    
    Parameters
    ----------
    obj : TrjObj object
      self object
    varlist : list
      List of contours to be plotted
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
    
    carray : string
      NetCDF code of array used for color coding
    centdiff : bool
      If true apply centered differencing on carray
    sigma : float
      If value given, apply filter on carray
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
        pmat = rootgrp.variables[carray][:, :]
    
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
            
            if centdiff:
                if sigma != None:
                    parray = ndi.filters.gaussian_filter(parray, sigma)
                parray = np.gradient(parray)
            
            single_trj(lonarray, latarray, parray, linewidth = linewidth, 
                       carray = carray)

    cb = plt.colorbar(lc, shrink = 0.7) 
    cb.set_label(carray)
    if carray == 'P':
        cb.ax.invert_yaxis()
    plt.tight_layout()
    
    
    # Save Plot
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()


def draw_trj_evo(obj, varlist, loclist, idlist, tplot, 
                 idtext = '', onlybool = False, startarray = None, 
                 stoparray  = None, savename = None, inrange = None,
                 linewidth = 0.7, limited = None):
    """
    obj : TrjObj object
      self object
    varlist : list
      List of contours to be plotted
    filelist : list
      List of unique file locations
    idlist : list
      List of list of trajectory IDs 
    tplot : int 
      Time to be plotted in mins after simulation start
    idtext : string
      Text to be displayed in plot
    savename : string
      Name of save location
      
    """
    
    # Plot contours
    draw_contour(obj, varlist, tplot, idtext = idtext)
    
    cnt = 0   # continuous counter for startarray and stoparray
    # Plot trajectories
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        # Retrieve starting time
        rootgrp = nc.Dataset(loclist[i], 'r')
        startt = rootgrp.variables['time'][0] / 60   # Convert to minutes

        # Only plot if trajectories start before tplot
        if startt <= tplot:
            
            # Get trjind for end of arrays
            trjind = int((tplot - startt) / obj.dtrj)

            # Retrieve arrays
            lonmat = rootgrp.variables['longitude'][:trjind, :]
            latmat = rootgrp.variables['latitude'][:trjind, :]
            pmat = rootgrp.variables['P'][:trjind, :]
            
            lonmat[:, :] += (180 - obj.pollon)   # Convert to real coordinates
            latmat[:, :] += (90 - obj.pollat)
            
            for j in idlist[i]:
                # Filter out zero values!
                if onlybool:
                    stop = min(stoparray[cnt], pmat.shape[0] - 1)
                    parray = pmat[startarray[cnt]:stop, j]
                    lonarray = lonmat[startarray[cnt]:stop, j]
                    latarray = latmat[startarray[cnt]:stop, j]
                    cnt += 1
                else:
                    parray = pmat[:, j][pmat[:, j] != 0]
                    lonarray = lonmat[:, j][pmat[:, j] != 0]
                    latarray = latmat[:, j][pmat[:, j] != 0]
                
                
                if inrange != None:
                    lonarray = lonarray[(parray > inrange[0]) & 
                                        (parray < inrange[1])]
                    latarray = latarray[(parray > inrange[0]) & 
                                        (parray < inrange[1])]
                    parray = parray[(parray > inrange[0]) & 
                                        (parray < inrange[1])]
                if limited != None:
                    lim = limited / obj.dtrj

                    lonarray = lonarray[-lim:]
                    latarray = latarray[-lim:]
                    parray = parray[-lim:]
                    
                single_trj(lonarray, latarray, parray, linewidth = linewidth)
    
    cb = plt.colorbar(lc, shrink = 0.7)
    cb.set_label('p')
    cb.ax.invert_yaxis()
    plt.tight_layout()
    
    # Save Plot
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()

def draw_trj_dot(obj, varlist, loclist, idlist, tplus, 
                 savename = None, idtext = '', inrange = None, cafter = None):
    """
    tplus = time after MODEL start
    """
    
    draw_contour(obj, varlist, tplus, idtext = idtext)
    cnt = 0   # continuous counter for startarray and stoparray
    
    # Plot trajectories
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
        trjstart = int(rootgrp.variables['time'][0] / 60)
        trjind = (tplus - trjstart) / obj.dtrj
        if trjind <= 0:
            break   # Break out of loop
        lonarray = rootgrp.variables['longitude'][trjind, idlist[i]]
        latarray = rootgrp.variables['latitude'][trjind, idlist[i]]
        parray = rootgrp.variables['P'][trjind, idlist[i]]
        
        if not inrange == None:
            norm = plt.Normalize(inrange[0], inrange[1])
            carray = parray[(parray > inrange[0]) & (parray < inrange[1])]
            if cafter != None:
                carray = cafter[cnt:cnt+parray.shape[0]]
                carray = carray[(parray > inrange[0]) & (parray < inrange[1])]
                carray = tplus - carray
                norm = plt.Normalize(-2880, 2880)
                
            lonarray = lonarray[(parray > inrange[0]) & (parray < inrange[1])]
            latarray = latarray[(parray > inrange[0]) & (parray < inrange[1])]
        else:
            norm = plt.Normalize(100, 1000)
        
        cmap = clr.ListedColormap(['darkorange', 'orange', 'khaki', 'beige',
                                   'greenyellow', 'lawngreen', 'green', 
                                   'darkolivegreen'])
        
        plt.scatter(lonarray, latarray, c = carray, s = 10,
                    cmap = cmap, linewidth = 0.1,
                    norm = norm)
        
        
    # Set plot properties
    plt.xlim(obj.xlim)
    plt.ylim(obj.ylim)
    plt.title(obj.date + timedelta(minutes = tplus))
    cb = plt.colorbar(shrink = 0.7) 
    #cb.set_label('P')
    #cb.ax.invert_yaxis()
    plt.tight_layout()
    
    # Save Plot
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()
        
        

def draw_asc_loc(obj, lon, lat, p, varlist, tplot, idtext = '', savename = None):
    """
    Plots map at tplot with colorcoded (key: p) dots (lon, lat). 
    
    obj : TrjObj object
      self object
    lon : np.array
      Array of x positions
    lat : np.array
      Array of y positions
    p : np.array
      Array for color coding
    varlist : list
      List of contours to be plotted
    tplot : float
      Time of contour plot in mins after model start
    idtext : string
      text diplayed on plot
    savename : string
      Name of save location
    """
    
    
    # Plot contours
    draw_contour(obj, varlist, tplot, idtext = idtext, savename = None)
    
    # Convert to Real coords
    lon += (180 - obj.pollon)   
    lat += (90 - obj.pollat)

    
    # Plot trajectory dots
    plt.scatter(lon, lat, c = p, cmap = plt.get_cmap('Spectral'), 
                norm=plt.Normalize(100, 1000), linewidth = 0.1)
    
    # Draw colorbar
    if lon.shape[0] != 0:
        cb = plt.colorbar(shrink = 0.7)
        cb.set_label('p')
        cb.ax.invert_yaxis()
    plt.tight_layout()
    
    # Save Plot
    if savename != False:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
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
    
    
def contour(filelist, variable, cosmoind, xlim, ylim, zlevel = None):
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
    elif zlevel != None:
        print zlevel
        field = pwg.getfield(filelist[cosmoind], variable, levs = zlevel)
    else:
        field = pwg.getfield(filelist[cosmoind], variable)
     
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
        #plt.contourf(X, Y, field, alpha = 0.5, zorder = 0.45)
        #plt.colorbar(shrink = 0.7)
        plt.contour(X, Y, field, levels = list(np.arange(150, 350, 4)), 
                     colors = "r", linewidths = 0.5, zorder = 0.45)
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
        cbar = plt.colorbar(shrink = 0.7, pad = 0)
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
                     extend = 'max', levels = levels, alpha = 0.8, 
                     zorder = 0.45)
        cbar = plt.colorbar(shrink = 0.7)
        cbar.set_label('CAPE [J/kg]', rotation = 90)
    elif variable == 'var4':   # PV
        #field = smoothfield(field, 2)
        field = field * 1.e6   # PVU
        levels = list(np.arange(-3, 8, 1))
        plt.contourf(X, Y, field, cmap = plt.get_cmap('coolwarm'), 
                     extend = 'both', levels = levels, alpha = 0.5,
                     zorder = 0.45)
        cbar = plt.colorbar(shrink = 0.7, pad = 0)
        cbar.set_label('PV [PVU]', rotation = 90)
        plt.contour(X, Y, field, color = 'black', levels = [2.], zorder = 0.5, 
                    linewidth = 4.)
        
    else:   # All other fields, unformatted
        plt.contourf(X, Y, field)
        #plt.colorbar()
    del field
    
    
def single_trj(lonarray, latarray, parray, linewidth = 0.7, carray = 'P'):
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
    
    # get colormap and normalize according to carray
    if carray == 'P':
        #cmap = clr.ListedColormap(['linen', 'yellow', 'greenyellow', 
                                   #'darkgreen', 'grey', 'grey', 'grey', 'grey',
                                   #'grey'])
        cmap = plt.get_cmap('Spectral')
        norm = plt.Normalize(100, 1000)
    elif carray == 'w':
        cmap = plt.get_cmap('Reds')
        norm = plt.Normalize(0, 2)
    elif carray == 'z':   #NOTE: Temporary solution
        cmap = clr.ListedColormap(['lime', 'aqua', 'darkblue', 
                                       'lightblue', 'lightsalmon', 
                                       'r', 'violet', 'purple'])
        norm = plt.Normalize(-200, 200)
    elif carray == 'POT_VORTIC':
        cmap = plt.get_cmap('Spectral')
        norm = plt.Normalize(-5e-6, 10e-6)
    else:
        #cmap = clr.ListedColormap(['white', 'linen', 'yellow', 'greenyellow', 
                                   #'darkgreen', 'grey', 'grey', 'grey', 'grey',
                                   #'grey'])
        cmap = plt.get_cmap('Spectral')
        norm = plt.Normalize(0, 1000)
    lc = col.LineCollection(segments, cmap = cmap, norm = norm, alpha = 0.5)
    lc.set_array(p)
    lc.set_linewidth(linewidth)
    plt.gca().add_collection(lc)
    
    
def smoothfield(field, sigma = 8):
    """
    Smoothes given field
    """
    newfield = ndi.filters.gaussian_filter(field, sigma)
    return newfield
        
    
    
    
    
    
    
    
