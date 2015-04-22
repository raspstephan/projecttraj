"""
plots module
------------

This module contains all plotting functions.
To make movie use from command line.

ffmpeg -r 3 -pattern_type glob -i '*.png' -vcodec mpeg4 output.mp4
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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.collections as col
from matplotlib.colorbar import make_axes, ColorbarBase
from mpl_toolkits.basemap import Basemap
import cosmo_utils.pywgrib as pwg
import scipy.ndimage as ndi
from scipy.stats import nanmean, pearsonr
from datetime import datetime, timedelta
import netCDF4 as nc
import utils
from mpl_toolkits.mplot3d import Axes3D
import scipy as sp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import RectBivariateSpline
import random
mpl.rcParams['font.size'] = 8
mpl.rcParams['font.family'] = 'sans-serif'




def draw_vs_t(obj, tracer, loclist, idlist, savename = None, sigma = None):
    """
    TODO
    """
    sigma = 1
    ylim = (0, 1000)
    xlim = None
    idtext = ''
    matlist = []
    exttarray = np.arange(0, obj.maxmins + obj.dtrj, obj.dtrj, dtype = 'float')
    exttarray = exttarray / 60.
    for i in range(len(loclist)):
        print 'Plotting file', i+1, 'of', len(loclist)
        rootgrp = nc.Dataset(loclist[i], 'r')
            
        tracemat = rootgrp.variables[tracer][:, :]
        tstart = rootgrp.variables['time'][0] / 60. / obj.dtrj   # Convert to minutes
        print tstart
        tstart = int(tstart)
        exttracemat = np.empty((exttarray.shape[0], len(idlist[i])))
        exttracemat.fill(np.nan)
        
        # Convert zeros to nans
        tracemat[tracemat == 0] = np.nan
        

            
        
        # Write in file
        exttracemat[tstart:] = tracemat[:, idlist[i]]
        matlist.append(exttracemat)
        
    totmat = np.hstack(matlist)
    meanarray = np.nanmean(totmat, axis = 1)
    countlist = np.sum(np.isfinite(totmat), axis = 1)
    mask = np.isfinite(meanarray)
    
    # Set up figure
    fig = plt.figure(figsize = (10, 8))
    ax = plt.gca()
    
    # Set plot properties
    if ylim == None:
        ylim = (1.5 * np.nanmin(meanarray), 1.5 * np.nanmax(meanarray))
    if xlim == None:
        xlim = [exttarray[mask].min(), exttarray[mask].max()]
        
    per5 = nanpercentile(totmat, 5)
    per95 = nanpercentile(totmat, 95)
    smoothper5 = ndi.filters.gaussian_filter(per5[mask], sigma)
    smoothper95 = ndi.filters.gaussian_filter(per95[mask], sigma)
    ax.fill_between(exttarray[mask], smoothper5, smoothper95, 
                    facecolor = 'lightgrey', edgecolor = 'lightgrey',
                    label = '90%', alpha = 0.5)
    
    per25 = nanpercentile(totmat, 25)
    per75 = nanpercentile(totmat, 75)
    smoothper25 = ndi.filters.gaussian_filter(per25[mask], sigma)
    smoothper75 = ndi.filters.gaussian_filter(per75[mask], sigma)
    ax.fill_between(exttarray[mask], smoothper25, smoothper75, 
                    facecolor = 'darkgrey', edgecolor = 'darkgrey',
                    label = '50%', alpha = 0.5)
    
    per50 = nanpercentile(totmat, 50)
    smoothper50 = ndi.filters.gaussian_filter(per50[mask], sigma)
    l2, = plt.plot(exttarray[mask], smoothper50, 'black')
    
    smootharray = ndi.filters.gaussian_filter(meanarray[mask], sigma)
    l1, = plt.plot(exttarray[mask], smootharray, 'firebrick', linewidth = 2)
    # Get filled colors as legends
    r1 = plt.Rectangle((0, 0), 1, 1, fc="lightgrey")
    r2 = plt.Rectangle((0, 0), 1, 1, fc="darkgrey")
    plt.legend([l1, l2, r2, r1], ['mean', 'median', '50%', '90%'])
    
    del totmat
    
    ax.set_ylim(ylim)
    if tracer == 'P':
        ax.invert_yaxis()
    dtick = 3
    #ax.xaxis.set_ticks(np.arange(0, 120, dtick))
    ax.set_xlim(xlim)
    ax.grid(color = 'dimgrey', linestyle = '-')
    ax.set_frame_on(False)
    plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    maxbin = np.max(countlist)
    
    ax2 = ax.twinx()
    ax2.bar(exttarray[mask], countlist[mask], linewidth = 0, color = 'darkgreen', 
            width = 5, alpha = 0.8)
    ax.set_zorder(2)
    ax2.set_zorder(1)
    inc = 500
    ax2.set_yticks(np.arange(inc, np.max(countlist) + inc, inc))
    ax2.set_ylabel('Number of Trajectories', position = (0.1, 0.175))
    ax2.set_xlim(xlim)
    ax2.set_ylim((0, np.max(countlist) * 4))
    
    ax.set_ylabel(tracer)
    ax.set_xlabel('Time [hrs] relative to center') 
    
    
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')


def draw_spaghetti(obj, tracer, loclist, idlist, savename = None, 
                   limit = (0, 50), locind = False, xlim = None):
    """
    TODO
    """
    if not locind == False:
        locind = locind
    else:
        locind = random.randint(0, len(loclist)-1)
    print locind
    rootgrp = nc.Dataset(loclist[locind], 'r')
    
    tmat = rootgrp.variables['time'][:]
    pmat = rootgrp.variables['P'][:, idlist[locind]]
    tracemat = rootgrp.variables[tracer][:, idlist[locind]]
    tracemat[pmat == 0] = np.nan
    tracemat[np.isnan(pmat)] = np.nan
    print tracemat.shape[1]
    tmat = tmat / 3600.
    
    
    fig = plt.figure(figsize = (5, 4))
    plt.gca().set_color_cycle(['red', 'maroon', 'darkorange', 'yellow', 
                        'lawngreen', 'darkgreen', 'cyan', 'navy', 'fuchsia'])
    for j in range(limit[1]):
        ind = random.randint(0, tracemat.shape[1]-1)
        plt.plot(tmat, tracemat[:, ind], linewidth = 2.5)
    if tracer == 'P':
        plt.gca().invert_yaxis()
    plt.xlabel('Time after simulation start [hrs]')
    plt.ylabel('Pressure [hPa]')
    plt.tight_layout()
    if xlim != None:
        plt.xlim(xlim)
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, dpi = 300, bbox_inches = 'tight')
        plt.close('all')
    

def draw_slice_hist(obj, tracer, loclist, idlist, threshold, startarray, 
                    stoparray, savename = None, mult = -1):
    """
    TODO
    """
    
    startlist = []
    stoplist = []
    cnt = 0
    #print loclist
    for i in range(len(loclist)):
        print 'Opening file:', loclist[i]
        rootgrp = nc.Dataset(loclist[i], 'r')
        tracemat = rootgrp.variables[tracer][:, idlist[i]]
        pmat = rootgrp.variables['P'][:, idlist[i]]
        
        for j in range(tracemat.shape[1]):
            #print j
            tracearray = tracemat[startarray[cnt]:stoparray[cnt], j] * mult
            parray = pmat[startarray[cnt]:stoparray[cnt]+1, j]
            
            mask = np.ma.masked_where(tracearray < threshold, tracearray)
            slices = np.ma.notmasked_contiguous(mask)
            
            if type(slices) == list:
                for s in slices:
                    istart = s.start
                    istop = s.stop - 1
                    diff = istop - istart
                    if diff > 1:
                        startlist.append(parray[istart])
                        stoplist.append(parray[istop])
            elif slices != None:
                istart = slices.start
                istop = slices.stop - 1
                diff = istart- istop
                if diff > 1:
                    startlist.append(parray[istart])
                    stoplist.append(parray[istop])
            cnt += 1  
        rootgrp.close()
    #print 'Plotting'        
    fig = plt.figure(figsize = (3.15, 3.15))
    ax = plt.gca()
    ax.set_aspect('equal')
    H, _, _ = np.histogram2d(startlist, stoplist, bins = 45, 
                         range = [[150, 1050], [150, 1050]])
    normmax = int(np.max(H))
    plt.hist2d(startlist, stoplist, bins = 45, range = [[150, 1050], [150, 1050]],
              norm = plt.Normalize(1, normmax), cmin = 1)
    plt.xlabel('Start of ascent [hPa]')
    plt.ylabel('End of ascent [hPa]')
    #ax.set_xlim(300, 1050)
    #ax.set_ylim(300, 1050)
    ax.invert_xaxis()
    ax.invert_yaxis()
    cb = plt.colorbar(shrink = 0.78)
    cb.set_label('Frequency')
    if savename != None:
        print 'Save figure as', savename
        fig.subplots_adjust(left = 0.17, right = 0.95, top = 0.95, bottom = 0.15)
        fig.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
    
    
    fig2 = plt.figure(figsize = (3.15, 2.3))
    diff = np.array(startlist) - np.array(stoplist)
    plt.hist(diff, bins = 50, range = [0, 1000], color = '#CCCCCC', 
             linewidth = 0.3)
    plt.xlabel('Ascent range [hPa]')
    plt.ylabel('Frequency')
    if savename != None:
        savename = savename + 'diffhist'
        print 'Save figure as', savename
        plt.subplots_adjust(left = 0.15, right = 0.925, top = 0.95, bottom = 0.19)
        fig2.savefig(savename, dpi = 400)
        plt.close('all')
    
     


        

def draw_vs_p(obj, tracer, loclist, idlist, startarray, stoparray, xlim, 
              savename = None, sigma = 1, idtext = '', ylim = None, 
              binwidth = 5., ylabel = '', log = False, legnames = None,
              legpos = 1, ax2upper = 300, mult = 1, letter = ''):
    """
    TODO
    """
    print tracer
    
    returnlist = []
    color = []
    if len(loclist) == 1:
        color.append(['lightgray', 'darkgray', 'black', 'black'])
        hatch = ['X']
    else:
        color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
        #color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
        color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
        hatch = ['/', '\\']
    
    fig = plt.figure(figsize = (3.15, 2.5))
    ax = plt.gca()
    ax.grid(color = 'dimgrey', linestyle = '-', zorder = 0.01)
    #ax.set_frame_on(False)
    plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes,
            fontsize = 6)
    plt.text(0.93, 0.93, letter, transform = plt.gca().transAxes,
            fontsize = 8)
    if not ax2upper == None:
        ax2 = ax.twinx()
    
    # Initialize p array
    parray = np.arange(0, 1050., binwidth) 
    returnlist.append(parray)
    nbins = parray.shape[0]
    leglist = []
    for n in range(len(loclist)):
        counter = 0
        total = 0
        tracerlist = [[] for j in range(nbins)]     
        for i in range(len(loclist[n])):
            total += len(idlist[n][i])
            print 'Plotting file', i+1, 'of', len(loclist[n])
            rootgrp = nc.Dataset(loclist[n][i], 'r')
            tracemat = rootgrp.variables[tracer[n]][:, :]
            tracemat = tracemat * mult
            if tracer[n] in ['POT_VORTIC', 'var4']:
                tracemat = tracemat * 1.e6
            pmat = rootgrp.variables['P'][:, :]

            # Get P data and set zeros to nan
            zmask = np.ma.mask_or(pmat == 0, np.isnan(pmat))
            tracemat[zmask] = np.nan
            
            # Convert pmat to indices
            pindmat = np.around((pmat) / binwidth)
            
            for j in idlist[n][i]:
                
                tracearray = tracemat[startarray[n][counter]:
                                          stoparray[n][counter], j]
                pindarray = pindmat[startarray[n][counter]:
                                        stoparray[n][counter], j]
                
                for k in range(tracearray.shape[0]):
                    #print pindarray[k], nbins
                    tracerlist[int(pindarray[k])].append(tracearray[k])
                    #print tracearray[k]
                    
                counter += 1
        
        countlist = []
        meanlist = []
        per5 = []
        per25 = []
        per50 = []
        per75 = []
        per95 = []
        
        # Set lower density limit TODO 
        denslim = 0.001
        
        for j in range(nbins):
            tmparray = np.array(tracerlist[j])
            tmparray[tmparray > 1e20] = np.nan
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
        
        per5 = np.array(per5)
        per25 = np.array(per25)
        per50 = np.array(per50)
        per75 = np.array(per75)
        per95 = np.array(per95)
        meanlist = np.array(meanlist)
        #mask = np.isfinite(meanlist)
        
        dashes1 = [1, 1]
        dashes2 = [3, 1]
        ax.fill_between(parray, per5, per95,
                        facecolor = color[n][0], edgecolor = 'none',
                        label = '90%', alpha = 0.5, zorder = 0.1)
        ax.fill_between(parray, per25, per75,
                        facecolor = color[n][1], edgecolor = 'none',
                        label = '50%', alpha = 0.5, zorder = 0.1)
        l2, = ax.plot(parray, per50, color[n][2], linestyle = '--', 
                      linewidth = 1.5)
        l1, = ax.plot(parray, meanlist, color[n][2], linewidth = 1.5)
        leglist.append(l1)
        l3, = ax.plot(parray, per5, 
                color = color[n][3], linestyle = ':', linewidth = 1)
        l4, = ax.plot(parray, per95, color = color[n][3], linestyle = ':', 
                      linewidth = 1)
        l5, = ax.plot(parray, per25, color = color[n][3], linestyle = '-.', 
                      linewidth = 1)
        l6, = ax.plot(parray, per75, color = color[n][3], linestyle = '-.', 
                      linewidth = 1)
        
        l3.set_dashes(dashes1)
        l4.set_dashes(dashes1)
        l5.set_dashes(dashes2)
        l6.set_dashes(dashes2)
        # Get filled colors as legends
        #r1 = plt.Rectangle((0, 0), 1, 1, fc="lightgrey")
        #r2 = plt.Rectangle((0, 0), 1, 1, fc="darkgrey")
        #plt.legend([l1, l2, r2, r1], ['mean', 'median', '50%', '90%'], loc = 2)
        # Plot second axis
        zero = np.zeros(len(countlist))
        relcountlist = np.array(countlist) / float(total) * 100.
        returnlist.append(relcountlist)
        if not ax2upper == None:
            fill_between_steps(parray, np.array(relcountlist), zero, ax = ax2, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
    if legnames != None:
        plt.legend(leglist, legnames, loc = legpos, fontsize = 8)  
    maxbin = np.max(countlist)
    print maxbin
    if not ax2upper == None:
        maxbin = ax2upper
        inc = 50
        ax2.set_yticks(np.arange(inc, maxbin + inc, inc))
        ax2.set_ylim((0, maxbin * 4))
        ax2.set_ylabel('Percent of trajectories in bin', position = (0.1, 0.175))
            
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('p [hPa]')
    ax.set_ylabel(ylabel)
    if not ax2upper == None:
        ax2.set_xlim(xlim)
    ax.invert_xaxis()
    if log:
        ax.set_yscale('log')


    if savename != None:
        print 'Save figure as', savename
        plt.subplots_adjust(left = 0.17, right = 0.95, top = 0.95, bottom = 0.15)
        plt.savefig(savename, dpi = 400)
        plt.close('all')
    
    return returnlist

def fill_between_steps(x, y1, y2=0, h_align='mid', ax=None, **kwargs):
    ''' Fills a hole in matplotlib: fill_between for step plots.
    Parameters :
    ------------
    x : array-like
    Array/vector of index values. These are assumed to be equally-spaced.
    If not, the result will probably look weird...
    y1 : array-like
    Array/vector of values to be filled under.
    y2 : array-Like
    Array/vector or bottom values for filled area. Default is 0.
    **kwargs will be passed to the matplotlib fill_between() function.
    '''
    # If no Axes opject given, grab the current one:
    if ax is None:
        ax = plt.gca()
    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = sp.repeat((x[1:] - x[:-1]), 2)
    xstep = sp.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    # Now: add one step at end of row.
    xx = sp.append(xx, xx.max() + xstep[-1])
    
    # Make it possible to chenge step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep
    
    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)#[:-1]
    if type(y2) == sp.ndarray:
        y2 = y2.repeat(2)#[:-1]
    
    # now to the plotting part:
    ax.fill_between(xx, y1, y2=y2, **kwargs)
    
    return ax


def draw_centered_vs_t(obj, loclist, idlist, tracer, carray, savename = None,
                       plottype = 'Smooth', idtext = '', ylim = None, 
                       xlim = None, sigma = 1, select = False,
                       legnames = None, legpos = 2, ax2 = True, mult = 1, 
                       letter = '', ylabel = ''):
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
    returnlist = []
    color = []
    leglist = []
    if len(loclist) == 1:
        color.append(['lightgray', 'darkgray', 'black', 'black'])
        hatch = ['X']
    else:
        color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
        #color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
        color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
        hatch = ['/', '\\']
    
    # Set up figure
    fig = plt.figure(figsize = (3.15, 2.5))
    ax = plt.gca()
    if ax2:
        ax2 = ax.twinx()
    # Set plot properties
    
    ax.plot([0, 0], ylim, color = 'dimgrey', linewidth = 2)
    ax.plot(xlim, [0, 0], color = 'dimgrey', linewidth = 2)
    
    for i in range(len(loclist)):
    
        totmat, exttarray = utils._centered_mat(obj[i], loclist[i], idlist[i], tracer[i], 
                                                    carray[i], mult)
        #print exttarray.shape
        if select != False:
            # Get start and end indices
            startind = np.argmin(np.abs(exttarray - (select[0] * 60)))
            stopind = np.argmin(np.abs(exttarray - (select[1] * 60)))
            #print totmat.shape
            #print startind, stopind
            #totmat = totmat[startind:stopind, :]
            #exttarray = exttarray[startind:stopind]
            #print totmat.shape
            selectarray = np.where(np.isfinite(np.sum(totmat[startind:stopind, :],
                                                    axis = 0)))
            print selectarray
            #print selectarray, np.sum(totmat, axis = 0)
            totmat = totmat[:, selectarray[0]]
            
            #totmat[totmat == 0] = np.nan
            #print totmat.shape
            
        meanarray = np.nanmean(totmat, axis = 1)
        countlist = np.sum(np.isfinite(totmat), axis = 1)
        
        print "Mean at t = 0:", meanarray[np.round(meanarray.shape[0]/2)]

        mask = np.isfinite(meanarray)
        
        # Convert time array to hours
        exttarray = exttarray / 60.

        if plottype == 'Fill':
            
            dashes1 = [1, 1]
            dashes2 = [3, 1]
            per5 = nanpercentile(totmat, 5)
            per95 = nanpercentile(totmat, 95)
            smoothper5 = ndi.filters.gaussian_filter(per5[mask], sigma)
            smoothper95 = ndi.filters.gaussian_filter(per95[mask], sigma)
            ax.fill_between(exttarray[mask], smoothper5, smoothper95, 
                            facecolor = color[i][0], edgecolor = 'none',
                            label = '90%', alpha = 0.5, zorder = 0.1)
            
            per25 = nanpercentile(totmat, 25)
            per75 = nanpercentile(totmat, 75)
            smoothper25 = ndi.filters.gaussian_filter(per25[mask], sigma)
            smoothper75 = ndi.filters.gaussian_filter(per75[mask], sigma)
            ax.fill_between(exttarray[mask], smoothper25, smoothper75, 
                            facecolor = color[i][1], edgecolor = 'none',
                            label = '50%', alpha = 0.5, zorder = 0.1)
            
            per50 = nanpercentile(totmat, 50)
            smoothper50 = ndi.filters.gaussian_filter(per50[mask], sigma)
            l2, = ax.plot(exttarray[mask], smoothper50, color[i][2], 
                           linestyle = '--', linewidth = 1.5)
            
            smootharray = ndi.filters.gaussian_filter(meanarray[mask], sigma)
            l1, = ax.plot(exttarray[mask], smootharray, color[i][2], 
                           linewidth = 1.5)
            leglist.append(l1)
            l3, = ax.plot(exttarray[mask], smoothper5, 
                color = color[i][3], linestyle = ':', linewidth = 1)
            l4, = ax.plot(exttarray[mask], smoothper95, color = color[i][3], linestyle = ':', 
                        linewidth = 1)
            l5, = ax.plot(exttarray[mask], smoothper25, color = color[i][3], linestyle = '-.', 
                        linewidth = 1)
            l6, = ax.plot(exttarray[mask], smoothper75, color = color[i][3], linestyle = '-.', 
                        linewidth = 1)
            l3.set_dashes(dashes1)
            l4.set_dashes(dashes1)
            l5.set_dashes(dashes2)
            l6.set_dashes(dashes2)
            # Get filled colors as legends
            #r1 = plt.Rectangle((0, 0), 1, 1, fc="lightgrey")
            #r2 = plt.Rectangle((0, 0), 1, 1, fc="darkgrey")
            #plt.legend([l1, l2, r2, r1], ['mean', 'median', '50%', '90%'])
        else:
            print 'ERROR'
        
        maxbin = np.max(countlist)
        print maxbin
        del totmat
        relcountlist = countlist[mask] / float(maxbin) * 100.
        returnlist.append(exttarray[mask])
        returnlist.append(relcountlist)
        if ax2:
            zero = np.zeros(relcountlist.shape[0])
            #ax2.bar(exttarray[mask], countlist[mask], linewidth = 1, 
                    #color = color[i][2], width = 5, alpha = 0.5)
            ax2.fill_between(exttarray[mask], relcountlist, zero, 
                            color = 'none', edgecolor = color[i][2], 
                            alpha = 0.5, hatch = hatch[i])
            #ax2.plot(exttarray[mask], relcountlist, c = color[i][2])
            ax.set_zorder(2)
            ax2.set_zorder(1)
            #if maxbin > 15000:
                #inc = 5000
            #else:
                #inc = 1000
            inc = 10
            ax2.set_yticks(np.arange(inc, 100 + inc, inc))
            ax2.set_ylabel('Percent of Trajectories', position = (0.1, 0.175))
            ax2.set_xlim(xlim)
            ax2.set_ylim((0, 400))
            
    if legnames != None:
        print legnames, leglist
        ax.legend(leglist, legnames, loc = legpos, fontsize = 8)
    ax.set_ylim(ylim)
    if tracer[0] == 'P':
        ax.invert_yaxis()
    if (xlim[1] - xlim[0]) < 24:
        dtick = 1
    elif (xlim[1] - xlim[0]) < 48:
        dtick = 6
    else:
        dtick = 12
    ax.xaxis.set_ticks(np.arange(-120, 120, dtick))
    ax.set_xlim(xlim)
    ax.grid(color = 'dimgrey', linestyle = '-', zorder = 0.01)
    #ax.set_frame_on(False)
    plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    plt.text(0.93, 0.93, letter, transform = plt.gca().transAxes,
            fontsize = 8)
    if tracer == 'var4':
        tracer == 'PVU'
    ax.set_ylabel(ylabel)
    ax.set_xlabel('time [h] relative to center') 
    
    
    
    if savename != None:
        print 'Save figure as', savename
        plt.subplots_adjust(left = 0.17, right = 0.95, top = 0.95, bottom = 0.15)
        plt.savefig(savename, dpi = 400)
        plt.close('all')
    #returnlist = [returnlist[i] for i in [0, 1, 3]]
    return returnlist
    
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



def draw_centered_integr(obj, loclist, idlist, tracer, carray, savename = None,
                       idtext = '', ylim = None, 
                       xlim = None, sigma = 1, middleval = 0):
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
    
    totmat, exttarray = utils._centered_mat(obj, loclist, idlist, tracer, 
                                            carray)
    
    meanarray = nanpercentile(totmat, 50)
    meanarray = meanarray * 60 * obj.dtrj
    meanarray[np.isnan(meanarray)] = 0
    middle = np.round(meanarray.shape[0]/2)
    print middle
    
    integrarray = np.cumsum(meanarray)
    
    offset = integrarray[middle] - middleval
    integrarray -= offset
    
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
    
    plt.plot(exttarray[mask], integrarray[mask], 'firebrick', linewidth = 2)
    
    
    del totmat
    
    ax.set_ylim(ylim)
    if tracer == 'P_dt':
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
    ax.set_ylabel(tracer)
    ax.set_xlabel('Time [hrs] relative to center') 
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight', dpi = 300)
        plt.close('all')

def draw_hist_2d(obj, varname1, varname2, loclist, idlist, carray, dplus,
                 savename = None):
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
    
    fig = plt.figure(figsize = (3.15, 3.15))
    plt.gca().set_aspect(1000 / 50)
    plt.hist2d(varlist1, varlist2, cmin = 1, 
               norm = clr.LogNorm(1, 1000),
               bins = [x_edges, y_edges], zorder = 10)
    cb = plt.colorbar(use_gridspec = True, shrink = 0.75)
    cb.set_label('Frequency')
    plt.tight_layout()
    plt.gcf().subplots_adjust(bottom = 0.1, left = 0.1)
    plt.xlabel('Pressure [hPa]')
    plt.ylabel('PV [pvu]')
    plt.grid()
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight', dpi = 300)
        plt.close('all')


def draw_hist_stacked(obj, datalist, namelist, ylim, idtext = '', savename = None, 
                      exp = False, xlim = None, bins = 50, realdate = False,
                      ylabel = None, legpos = 1, ylog = False):
    """
    TODO
    """
    
    fig = plt.figure(figsize = (3.15, 3))
    ax = plt.gca()
    bins = np.linspace(xlim[0], xlim[1], bins + 1)

    num1, xpos = np.histogram(datalist[0], bins)
    avg1 = np.average(datalist[0])
    med1 = np.median(datalist[0])
    
    width = np.diff(xpos)[0]
    xpos = xpos[:-1]
    
    b1 = plt.bar(xpos, num1, bottom = None, color = '#CCCCCC', 
                            width = width, linewidth = 0.3)
    l11, = plt.plot([avg1, avg1], [ylim[0], ylim[1]], linestyle = '--', 
                   color = '#CCCCCC', label = 'Mean', linewidth = 1)
    l12, = plt.plot([med1, med1], [ylim[0], ylim[1]], linestyle = ':', 
                   color = '#CCCCCC', label = 'Median', linewidth = 1)
    
    num2, xpos2 = np.histogram(datalist[1], bins)
    avg2 = np.average(datalist[1])
    med2 = np.median(datalist[1])
    
    b2 = plt.bar(xpos, num2, bottom = num1, color = '#525252', 
                            width = width, linewidth = 0.3)
    l21, = plt.plot([avg2, avg2], [ylim[0], ylim[1]], linestyle = '--', 
                   color = '#525252', label = 'Mean', linewidth = 1)
    l22, = plt.plot([med2, med2], [ylim[0], ylim[1]], linestyle = ':', 
                   color = '#525252', label = 'Median', linewidth = 1)
    
    handles = [b1, l11, l12, b2, l21, l22]
    labels = [namelist[0], 'Mean: {:.2f}h'.format(avg1), 
              'Median: {:.2f}h'.format(med1), namelist[1], 
              'Mean: {:.2f}h'.format(avg2), 'Median: {:.2f}h'.format(med2)]
    
    if not exp == False:
        xvals = np.linspace(xlim[0], xlim[1], 1000)
        x0 = exp[0]
        c = exp[1]
        yvals = x0 * np.exp(c*xvals)
        l3, = plt.plot(xvals, yvals, linestyle = '-', color = 'black', 
                       linewidth = 1) 
        handles.append(l3)
        labels.append('Exponential fit')
        
    if ylog:
        plt.gca().set_yscale('log')
        
    plt.legend(handles, labels, loc = legpos)
    ax.set_ylim(ylim)
    
    
        
    d = timedelta(1)
    if realdate:
        datelist = [(obj.date + d*i) for i in range(10)]
        plt.xticks(np.arange(0, 10 * 24, 24), datelist, rotation = 0, 
                   fontsize = 6)
    else:
        plt.xlabel('Time [h]')
    if xlim[1] == 48:
        plt.xticks(np.arange(0, 48 + 4, 4))
    ax.set_xlim(xlim)
    plt.ylabel('Frequency')
    if not ylabel == None:
        plt.xlabel(ylabel)
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    
    if savename != None:
        print 'Save figure as', savename
        plt.subplots_adjust(left = 0.175, right = 0.975, top = 0.95, bottom = 0.125)
        plt.savefig(savename, dpi = 300)
        plt.close('all')


def draw_hist_3d(datalist, namelist, idtext = '', savename = None, ylim = None):
    """
    """
    fig = plt.figure()
    ax = plt.gca()
    #ax = fig.add_subplot(111, projection = '3d')
    
    #for i in range(len(datalist)):
        #histdata, histrange = np.histogram(datalist[i], 
                                           #bins = np.arange(0, 48, 0.5))
        #ax.bar(histrange[:-1], histdata, i, zdir = 'y', color = 'b', alpha = 0.8) 
    
    #ax.set_xlabel('hrs')
    #ax.set_ylabel('criterion')
    #ax.set_zlabel('Count')
    
    #plt.yticks(range(0, len(datalist)), namelist)
    
    #ax.set_zlim(0,1000)
    newcolors = ['#A3CC52', '#CCCC00', '#FFB84D', '#FF5C33', '#8F0000']
    colors = newcolors[:len(namelist)]
    plt.hist(datalist, bins = np.linspace(0, 48, 144), 
                 histtype = 'step', linewidth = 2, label = namelist, 
                 color = colors)
    plt.legend()  
    plt.gca().set_xlim(0, 48)
    if not ylim == None:
        plt.gca().set_ylim(ylim)
    
    plt.xlabel('Ascent time [hrs]')
    plt.ylabel('Frequency')
    
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight', dpi = 300)
        plt.close('all')
    
    

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
    #fig = plt.figure(figsize = (10, 10))
    #ax = plt.gca()   
    #ax.set_aspect('equal')
    if carray == None:
        carray = 'blue'   # Set to default
        
    
    print np.nanmean(array1), np.nanmean(array2)
    #plt.scatter(array1, array2)
    # Plot scatter plot, add labels
    
    array2 = 600. / array2 / 60.
    
    #sca = ax.scatter(array1, array2, c = carray, s = 8,
                     #cmap=plt.get_cmap('Spectral'), 
                     #norm=plt.Normalize(100, 1000), linewidths = 0)
    #plt.plot([0, 5], [0, 5])
    #plt.xlim(0, xmx+0.1*xmx)
    #plt.ylim(0, xmx+0.1*xmx)
    #plt.xlabel('Maximum vertical velocity [hPa/s]')
    #plt.ylabel('Cross-tropospheric average vertical velocity [hPa/s]')
    
    fig2 = plt.figure()
    print np.max(array1), np.max(array2)
    mask = [np.isfinite(array1)]
    plt.hist(array1[mask] / array2[mask], bins = 99, range = (1, 100), 
             color = '#CCCCCC')
    plt.ylabel('Frequency')
    plt.xlabel('ratio(maximum ascent velocity / average ascent velocity)')
    fig2.gca().set_xlim(1, 100)
    fig2.gca().set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # Set limits, plot diagnonal line, set ticks
    #xmx = np.amax(array1[np.isfinite(array1)])
    #ymx = np.amax(array2[np.isfinite(array2)])
    #mx = max(xmx, ymx)
    #plt.xticks(np.arange(0, xmx+3, 3))
    #plt.yticks(np.arange(0, xmx+3, 3))
    
    # Add idtext
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    #cb = fig.colorbar(sca, shrink = 0.7)
    #cb.set_label('p')
    #cb.ax.invert_yaxis()
    #ax.set_axis_bgcolor('grey')
    
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
        


def draw_hist(array, xlim, ylim, bins = 50, idtext = '', xlabel =  None, 
              savename = None, ylog = False, exp = False, mintohrs = False,
              legend = True):
    """
    TODO
      
    """
    # Set up figure
    fig = plt.figure(figsize = (3.15, 3))
    ax = plt.gca()
    # Convert to np array
    array = np.array(array)
    
    if mintohrs:
        array = array / 60.
    
    bins = np.linspace(xlim[0], xlim[1], bins + 1)
    avg = np.nanmean(array)
    med = np.nanmedian(array)
    
    # Plot histogram, remove nans
    num, xpos = np.histogram(array, bins)
    width = np.diff(xpos)[0]
    xpos = xpos[:-1]  
    plt.bar(xpos, num, color = '#CCCCCC', width = width, linewidth = 0.3)
    #plt.hist(array[np.isfinite(array)], bins = 100)
    if ylog:
        plt.gca().set_yscale('log')
    l1, = plt.plot([avg, avg], [ylim[0], ylim[1]], linestyle = '--', 
                   color = 'black', label = 'Mean', linewidth = 1)
    l2, = plt.plot([med, med], [ylim[0], ylim[1]], linestyle = ':', 
                   color = 'black', label = 'Median', linewidth = 1)
        
    
    if not exp == False:
        xvals = np.linspace(xlim[0], xlim[1], 1000)
        x0 = exp[0]
        c = exp[1]
        yvals = x0 * np.exp(c*xvals)
        l3, = plt.plot(xvals, yvals, linestyle = '-', color = 'black', 
                       linewidth = 1)
        if legend:
            plt.legend([l1, l2, l3], ['Mean: {:.2f}h'.format(avg) , 
                                    'Median: {:.2f}h'.format(med), 
                                    'Exponential fit'])
    else:
        if legend:
            plt.legend([l1, l2], ['Mean: {:.2f}h'.format(avg) , 
                                    'Median: {:.2f}h'.format(med)])
        
    
    # Add labels and text
    # plt.title('Total Number of trajectories: ' + str( array.shape[0]))
    plt.ylabel("Frequency")
    if xlabel == None:
        plt.xlabel('Time [h]')
    else:
        plt.xlabel(xlabel)
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    if xlim[1] == 48:
        plt.xticks(np.arange(0, 48 + 4, 4))
        
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    #ax.yaxis.labelpad = 5
    
    if savename != None:
        print 'Save figure as', savename
        plt.subplots_adjust(left = 0.175, right = 0.975, top = 0.95, bottom = 0.125)
        plt.savefig(savename, dpi = 300)
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
    m = draw_contour(obj, varlist, tplot, idtext = idtext)
    
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
    norm = plt.Normalize(0, 50)
    cmap = clr.ListedColormap(['khaki'] + ['yellow'] * 2 + 
                              ['greenyellow'] * 3 + ['darksage'] * 4 + 
                              ['darkolivegreen'])
    dbin = 1   # degrees
    x_edges = np.arange(obj.xlim[0], obj.xlim[1] + dbin, dbin)
    y_edges = np.arange(obj.ylim[0], obj.ylim[1] + dbin, dbin)
    density, _, _ = np.histogram2d(latlist, lonlist, [y_edges, x_edges])

    X, Y = np.meshgrid(x_edges, y_edges)
    X, Y = utils.rcoo_2_gcoo(X, Y, obj.pollon, obj.pollat)

    X, Y = m(X, Y)
    plt.pcolormesh(X, Y, density, norm = norm, cmap = cmap, alpha = 0.5)
    plt.colorbar(shrink = 0.7, extend = 'max')
    plt.tight_layout()
    
    if savename != None:
        print 'Save figure as', savename
        plt.savefig(savename, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()
    return density

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
        

def draw_contour(obj, varlist, time, idtext, savename = None, rcoo = False, 
                 setting = None, cbar = True, title = False, diffobj = None, 
                 diffvar = None, crop = False, dot = None, line = None):
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
    fig = plt.figure(figsize = (3.15,10))
    ax = plt.gca()   
    #ax.set_aspect('equal')
    
    # Draw basemap
    m = basemap(obj.cfile, obj.xlim, obj.ylim, rcoo)
    
    # Plotting all contour fields
    for var in varlist:   
        diffind = None
        difflist = None
        # NOTE: 'CUM_PREC' not implemented right now
        #if varlist[i] == 'CUM_PREC':
            #contour(rfiles, varlist[i], cosmoind, xlim, ylim, 
                    #trjstart = trjstart)
                    
        # Get index and filelist
        print time
        if type(var) in [list, tuple]:
            if var[0] == 'THETA':
                cosmoind, filelist = obj._get_index('T', time)
            else:
                cosmoind, filelist = obj._get_index(var[0], time)
            
            if not diffobj == None and var[0] == diffvar:
                print time
                tdelta  = (obj.date - diffobj.date).total_seconds() / 60.
                print time + tdelta
                diffind, difflist = diffobj._get_index(var[0], time + tdelta)
            
            if len(var) == 2:
                leveltype = 'PS'
            else:
                leveltype = var[2]
            
            contour(obj, filelist, var[0], cosmoind, obj.xlim, obj.ylim, var[1], 
                    m = m, setting = setting, cbar = cbar, diffobj = diffobj,
                    diffind = diffind, difflist = difflist, 
                    leveltype = leveltype, crop = crop)
        else: 
            cosmoind, filelist = obj._get_index(var, time)
            if not diffobj == None and var == diffvar:
                tdelta  = (obj.date - diffobj.date).total_seconds() / 60.
                diffind, difflist = diffobj._get_index(var, time + tdelta)
                
            contour(obj, filelist, var, cosmoind, obj.xlim, obj.ylim, m = m, 
                    setting = setting, cbar = cbar, diffobj = diffobj,
                    diffind = diffind, difflist = difflist, crop = crop)
        
    ## Set plot properties
    #plt.xlim(obj.xlim)
    #plt.ylim(obj.ylim)
    
    # Set labels and title
    #plt.xlabel('longitude')
    #plt.ylabel('latitude')
    
    if not dot == None:
        dotlon = dot[0]
        dotlat = dot[1]
        dotlon, dotlat = utils.rcoo_2_gcoo(dotlon, dotlat, obj.pollon, 
                                           obj.pollat)
        dotx, doty = m(dotlon, dotlat)
        m.scatter(dotx, doty, s = 40, c = 'limegreen', linewidth = 0, 
                  zorder = 100)

    if not line == None:
        dotlon = [line[0], line[2]]
        dotlat = [line[1], line[3]]
        dotlon, dotlat = utils.rcoo_2_gcoo(dotlon, dotlat, obj.pollon, 
                                           obj.pollat)
        dotx, doty = m(dotlon, dotlat)
        m.plot(dotx, doty, linewidth = 2.5, color = 'red')
        
    
    if title:
        plt.title(obj.date + timedelta(minutes = time))
    plt.text(0.94, 1.02, idtext, transform = plt.gca().transAxes, 
             fontsize = 6)
    #plt.tight_layout()
    
        
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()
    else:
        return m


def draw_trj(obj, varlist, filelist, idlist, cfile, rfiles, pfiles, 
             savename = False,pollon = None, pollat = None, xlim = None, 
             ylim = None, onlybool = False, startarray = None, 
             stoparray = None, trjstart = None, idtext = '', linewidth = 0.7,
             carray = 'P', centdiff = False, sigma = None, thinning = 1):
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
    m = draw_contour(obj, varlist, tmpt, idtext = idtext)
    
    # Plot trajectories
    cnt = 0   # continuous counter for startarray and stoparray
    for i in range(len(filelist)):
        print 'Plotting file', i+1, 'of', len(filelist)
        rootgrp = nc.Dataset(filelist[i], 'r')
        lonmat = rootgrp.variables['longitude'][:, :]
        latmat = rootgrp.variables['latitude'][:, :]
        pmat = rootgrp.variables[carray][:, :]
    
        lonmat, latmat = utils.rcoo_2_gcoo(lonmat, latmat, obj.pollon, 
                                               obj.pollat)
        #print lonmat, latmat
        lonmat, latmat = m(lonmat, latmat)
        #if i == 0:
            #print lonmat, latmat
        
        for j in idlist[i]:
            # Filter out zero values!
            
            if onlybool:
                parray = pmat[startarray[cnt]:stoparray[cnt], j]
                lonarray = lonmat[startarray[cnt]:stoparray[cnt], j]
                latarray = latmat[startarray[cnt]:stoparray[cnt], j]
            else:
                parray = pmat[:, j][np.isfinite(pmat[:, j])]
                lonarray = lonmat[:, j][np.isfinite(pmat[:, j])]
                latarray = latmat[:, j][np.isfinite(pmat[:, j])]
            
            if centdiff:
                if sigma != None:
                    parray = ndi.filters.gaussian_filter(parray, sigma)
                parray = np.gradient(parray)
            if cnt % thinning == 0:
                single_trj(lonarray, latarray, parray, linewidth = linewidth, 
                        carray = carray)
            cnt += 1

    divider = make_axes_locatable(plt.gca())
    cax1 = divider.append_axes('bottom', size = '3%', pad = 0.30)
    cb = plt.colorbar(lc, orientation = 'horizontal', cax = cax1) 
    cb.set_label('pressure [hPa]')
    if carray == 'P':
        cb.ax.invert_yaxis()
    plt.tight_layout()
    plt.title('')
    
    
    # Save Plot
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
        plt.close('all')
        plt.clf()

def draw_trj_path(obj, filelist, idlist, stop, thinning = 10, savename = None):
    """
    TODO
    """
    tmpt = 360  # Not used!
    m = draw_contour(obj, [], tmpt, '')
    clist = ['red', 'blue']
    
    
    
    for n in range(len(filelist)):
        cnt = 0
        lonlist = []
        latlist = []
        
        for i in range(len(filelist[n])):
            print 'Plotting file', i+1, 'of', len(filelist[n])
            #print filelist[n], i
            rootgrp = nc.Dataset(filelist[n][i], 'r')
            lonmat = rootgrp.variables['longitude'][:, :]
            latmat = rootgrp.variables['latitude'][:, :]
            lonmat, latmat = utils.rcoo_2_gcoo(lonmat, latmat, obj.pollon, 
                                               obj.pollat)
            #print lonmat, latmat
            lonmat, latmat = m(lonmat, latmat)
            
            for j in idlist[n][i]:
                ind = stop[n][cnt]-(360 / obj.dtrj)
                if ind >= 0 and cnt % thinning == 0:
                    lonlist.append(lonmat[ind, j])
                    latlist.append(latmat[ind, j])
                
                #if cnt % thinning == 0:
                    #m.plot(lonarray, latarray, color = clist[n])
                cnt += 1
        m.scatter(lonlist, latlist, color = clist[n], s = 1)
                
    if savename != None:
        print "Saving figure as", savename
        plt.savefig(savename, dpi = 300, bbox_inches = 'tight')
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
            
            
            lonmat, latmat = utils.rcoo_2_gcoo(lonmat, latmat, obj.pollon, 
                                               obj.pollat)
            
            
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
                 savename = None, idtext = '', inrange = None, cafter = None,
                 thinning = False, setting = None, path = False, cbar = True,
                 ctracer = 'P', diffobj = None, diffvar = None):
    """
    tplus = time after MODEL start
    """
    
    m = draw_contour(obj, varlist, tplus, idtext = idtext, setting = setting,
                     cbar = cbar, diffobj = diffobj, diffvar = diffvar)
    cnt = 0   # continuous counter for startarray and stoparray
    if len(loclist) == 1:
        colorlist = ['gray']
    else:
        colorlist = ['red', 'mediumseagreen']
        
    didplot = False
    # Plot trajectories
    for j in range(len(loclist)):
        for i in range(len(loclist[j])):
            print 'Plotting file', i+1, 'of', len(loclist[j])
            rootgrp = nc.Dataset(loclist[j][i], 'r')
            trjstart = int(rootgrp.variables['time'][0] / 60)
            trjind = (tplus - trjstart) / obj.dtrj
            if trjind <= 0:
                break
            lonmat = rootgrp.variables['longitude'][:, :]
            lonmat = lonmat[:, idlist[j][i]]
            latmat = rootgrp.variables['latitude'][:, :]
            latmat = latmat[:, idlist[j][i]]
            lonarray = lonmat[trjind, :]
            latarray = latmat[trjind, :]
            if type(ctracer) == str:
                carray = rootgrp.variables[ctracer][:, :]
                carray = carray[trjind, idlist[j][i]]
            else:
                carray = ctracer[i]

            lonarray, latarray = utils.rcoo_2_gcoo(lonarray, latarray, obj.pollon, 
                                                obj.pollat)
            lonmat[lonmat == 0] = np.nan
            latmat[latmat == 0] = np.nan
            lonmat, latmat = utils.rcoo_2_gcoo(lonmat, latmat, obj.pollon, 
                                                obj.pollat)
            #print lonmat, latmat
            
            if not inrange == None:
                norm = plt.Normalize(inrange[0], inrange[1])
                carray = carray[(carray > inrange[0]) & (carray < inrange[1])]
                if cafter != None:
                    carray = cafter[cnt:cnt+parray.shape[0]]
                    carray = carray[(parray > inrange[0]) & (parray < inrange[1])]
                    carray = tplus - carray
                    norm = plt.Normalize(-2880, 2880)
                    
                lonarray = lonarray[(carray > inrange[0]) & (carray < inrange[1])]
                latarray = latarray[(carray > inrange[0]) & (carray < inrange[1])]
            else:
                norm = plt.Normalize(100, 1000)
            
            cmap = 'Spectral'
            cblabel = 'Pressure [hPa]'
            if 'P_dt' in ctracer:
                norm = plt.Normalize(0, 0.1)
                cmap = 'Spectral_r'
                carray = carray * (-1)
                lonarray = lonarray[carray > 0.1]
                latarray = latarray[carray > 0.1]
                carray = carray[carray > 0.1]
                
            if ctracer == 'THETA':
                cmap = 'nipy_spectral'
                norm = plt.Normalize(280, 340)
                cblabel = 'potential temperature [K]'
            if not thinning == False:
                lonmat = lonmat[:, ::thinning]
                latmat = latmat[:, ::thinning]
                lonarray = lonarray[::thinning]
                latarray = latarray[::thinning]
                carray = carray[::thinning]
            
            #cmap = clr.ListedColormap(['darkorange', 'orange', 'khaki', 'beige',
                                    #'greenyellow', 'lawngreen', 'green', 
                                    #'darkolivegreen'])
            
            if path:
                xmat, ymat = m(lonmat, latmat)
                #print xmat, ymat
                plt.plot(xmat, ymat, color = colorlist[j], linewidth = 0.5, 
                         alpha = 0.8, zorder = 0.2)
            
            x, y = m(lonarray, latarray)
            #print len(carray), x.shape
            if not type(ctracer) == str:
                carray = 1 / np.array(carray) / 60 * 300
                norm = plt.Normalize(0, 0.2)
                cmap = 'Spectral_r'
            m.scatter(x, y, c = carray, s = 5,
                        cmap = cmap, linewidth = 0.2,
                        norm = norm, edgecolor = colorlist[j])
            didplot = True
        #break
        
        
    # Set plot properties
    #plt.xlim(obj.xlim)
    #plt.ylim(obj.ylim)
    #plt.title(obj.date + timedelta(minutes = tplus))
    #plt.colorbar(shrink = 0.3)
    if cbar:
        cb = ColorbarBase(cax2, orientation = 'horizontal', cmap = cmap, 
                       norm = norm, ticks = [100, 200, 300, 400, 500, 600, 700, 
                                             800, 900])
        cb.set_label(cblabel)
        cb.ax.invert_yaxis()
    #plt.tight_layout()
    
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


def basemap(cfile, xlim, ylim, rcoo = False):
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
    #dx = 0.025   # Assumes COSMO 2.8 resolution
    #field = pwg.getfield(cfile, "FR_LAND_S", invfile = False)
    
    ## Setting up grid
    #ny, nx = field.shape
    #x = np.linspace(xlim[0], xlim[1], nx)
    #y = np.linspace(ylim[0], ylim[1], ny)

    #X, Y = np.meshgrid(x, y)
    
    ## Plot and filter field
    #field[field != 0] = 1.
    #plt.contourf(X, Y, field, levels = [0.5, 2], colors = ["0.85"], zorder=0.4)
    ##plt.contour(X, Y, field, 2, colors = "0.5")
    
    #del field
    
    
    
    cosobj = pwg.getfobj(cfile, "FR_LAND_S")
    lats = cosobj.lats
    lons = cosobj.lons
    # from cosmoutils
    jmid = len(lats[:,0])/2
    imid = len(lons[0,:])/2
    latmid = lats[jmid,imid]
    lonmid = lons[jmid,imid]
   
    lllat = lats[0,0]
    lllon = lons[0,0]
    urlat = lats[-1,-1]
    urlon = lons[-1,-1]
    
    m = Basemap(projection="stere", lon_0=lonmid, lat_0=latmid, 
                llcrnrlat=lllat, urcrnrlat=urlat, llcrnrlon=lllon, 
                urcrnrlon=urlon, resolution='i')
    m.drawcoastlines(color = 'black', linewidth = 0.5)
    m.drawcountries(color = 'grey', linewidth = 0.1)
    
    npars = 5
    nmers = 5
    if npars>0 and nmers>0:
        #parallels = np.arange(int(lllat), int(urlat)+1, (int(urlat)+1-int(lllat))/float(npars))
        parallels = np.arange(-100, 100, 10)
        m.drawparallels(parallels, labels = [0,0,0,0], linewidth = 0.3)
        #meridians = np.arange(int(lllon), int(urlon)+1, (int(urlon)+1-int(lllon))/float(nmers))
        meridians = np.arange(-100, 100, 10)
        m.drawmeridians(meridians, labels = [0,0,0,0], linewidth = 0.3)
    
    return m
    
def contour(obj, filelist, variable, cosmoind, xlim, ylim, plevel = None, 
            m = None, setting = None, cbar = True, diffobj = None, 
            diffind = None, difflist = None, leveltype = 'PS', crop = False):
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
    global cbar1, cax2
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
            field1 = pwg.getfobj(filelist[cosmoind - 1], "TOT_PREC_S")
            field2 = pwg.getfobj(filelist[cosmoind + 1], "TOT_PREC_S")
            X, Y = m(field1.lons, field1.lats)
            field = (field2.data - field1.data)/(dt*2)*3600
        except IndexError:
            # Reached end of array, use precip from previous time step
            field1 = pwg.getfobj(filelist[cosmoind - 2], "TOT_PREC_S")
            field2 = pwg.getfobj(filelist[cosmoind], "TOT_PREC_S")
            X, Y = m(field1.lons, field1.lats)
            field = (field2.data - field1.data)/(dt*2)*3600
    # Retrieve regular fields
    elif plevel != None:
        print plevel
        field, flons, flats, fobj = utils._get_level(obj, filelist[cosmoind], 
                                               variable, plevel, leveltype = 
                                               leveltype)
        X, Y = m(flons, flats)
    else:
        fobj = pwg.getfobj(filelist[cosmoind], variable)
        field = fobj.data
        X, Y = m(fobj.lons, fobj.lats)
    
    if not diffind == None:
        # for now only simple fields and plevel
        if plevel != None:
            if plevel > 600:
                print 'Error for >600hPa due to orographic data gaps'
            dfield, dflons, dflats, dfobj = utils._get_level(diffobj, difflist[diffind], 
                                                   variable, plevel, leveltype = 
                                                   leveltype)
        else:
            dfobj = pwg.getfobj(difflist[diffind], variable)
            dfield = dfobj.data

        # Interpolate to DE grid
        
        #print fobj.rlats[:, 0], fobj.rlons[0, :]
        #print dfobj.rlats[:, 0], dfobj.rlons[0, :]
        # factor = dfobj.drlat / fobj.drlat        
        rbs = RectBivariateSpline(dfobj.rlats[:, 0], dfobj.rlons[0, :], dfield)
        intfield = rbs(fobj.rlats[:, 0], fobj.rlons[0, :])
        
        if crop:
            field = intfield
        else:
            # Take difference
            field = smoothfield(field, 15)
            intfield = smoothfield(intfield, 15)
            field = field - intfield
            variable = variable + '_diff'
    
        
        
        
        
    #print X, Y
    ## Setting up grid
    #ny, nx = field.shape
    #x = np.linspace(xlim[0], xlim[1], nx)
    #y = np.linspace(ylim[0], ylim[1], ny)
    #X, Y = np.meshgrid(x, y)
    if setting == 1:
        print 'Setting 1'
        if variable == "PMSL":   # Surface pressure
            field = smoothfield(field, 3)/100
            levels = list(np.arange(900, 1100, 5))
            CS = m.contour(X, Y, field, levels = levels, colors = 'green', 
                            linewidths = 0.5, zorder = 0.5)
            plt.clabel(CS, fontsize = 4, inline = 1, fmt = "%.0f")
        elif variable == 'THETA':
            field = smoothfield(field, 5)
            levels = np.arange(200, 350, 2)
            CS = plt.contour(X, Y, field, colors = 'black', levels = levels, 
                             zorder = 0.4, alpha = 1, linewidths = 1.5, 
                             linestyles = '-')
            plt.clabel(CS, inline = 1, fontsize = 7, fmt = "%.0f")
            #cbar = plt.colorbar(shrink = 0.7, pad = 0)
        elif variable in ['var4', 'var4_diff']:   # PV
            field = smoothfield(field, 3)
            field = field * 1.e6   # PVU
            levels = [-0.2, 0.2, 1, 1.5, 2, 3, 5]
            cmap = ['#80B2CC', '#FFFFFF', '#FFE0C2', '#FF9933', '#CC3300', 
                    '#FF5C33']
            if variable == 'var4_diff':
                levels = [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
                #levels = np.array(levels) / 10.
                cmap = [ 'red', 'tomato', 'lightsalmon', 'white', 'white',  'lightblue', 'cornflowerblue', 
                        'blue']
                #levels = np.arange(-1, 1, 0.1)
                #cmap = plt.get_cmap('bwr')
                #norm = plt.Normalize(-1, 1, 0.1)
                #plt.contourf(X, Y, field, cmap = cmap, 
                                #extend = 'both', alpha = 0.5,
                                #zorder = 0.45, norm = norm)
            
            CS = plt.contourf(X, Y, field, colors = cmap, 
                            extend = 'both', levels = levels, alpha = 0.5,
                            zorder = 0.45)
            if variable == 'var4_diff':
                CS.cmap.set_under('maroon')
                CS.cmap.set_over('darkblue')
            else:
                CS.cmap.set_under('#007ACC')
                CS.cmap.set_over('#FFCC00')
            
            ax = plt.gca()
            if cbar:
                divider = make_axes_locatable(ax)
                cax1 = divider.append_axes('bottom', size = '3%', pad = 0.20)
                #if len(obj.trjfiles) != 42:
                cax2 = divider.append_axes('bottom', size = '3%', pad = 0.40)
                cbar1 = plt.colorbar(orientation = 'horizontal', cax = cax1)
                cbar1.set_label('PV [pvu]')
            else: 
                divider = make_axes_locatable(ax)
                cax1 = divider.append_axes('bottom', size = '3%', pad = 0.20)
                cbar1 = plt.colorbar(orientation = 'horizontal', cax = cax1)
                plt.gcf().delaxes(cax1)
            plt.sca(ax)
            #levels = np.arange(-5, 10, 2)
            #CS = plt.contour(X, Y, field, colors = 'black', levels = levels, 
                        #zorder = 0.6, linewidths = 2.)
            #plt.clabel(CS, inline = 1, fontsize = 7, fmt = "%.0f")

    elif setting in [2, 5]:
        print 'Setting 2'
        if variable == "PMSL":   # Surface pressure
            field = smoothfield(field, 3)/100
            levels = list(np.arange(900, 1100, 5))
            CS = m.contour(X, Y, field, levels = levels, colors = 'green', 
                            linewidths = 0.5, zorder = 0.5)
            plt.clabel(CS, fontsize = 4, inline = 1, fmt = "%.0f")
        elif variable == 'var145_S':   # CAPE
            field = smoothfield(field, 3)
            if setting == 2:
                levels = list(np.arange(100, 3000, 100))
                ticks = [100, 500, 1000, 1500, 2000, 2500, 3000]
            elif setting == 5:
                levels = list(np.arange(10, 1000, 100))
                ticks = [10, 250, 500, 750, 1000]
            plt.contourf(X, Y, field, cmap = plt.get_cmap('hot_r'), 
                        extend = 'max', levels = levels, alpha = 0.8, 
                        zorder = 0.45)
            ax = plt.gca()
            if cbar:
                divider = make_axes_locatable(ax)
                cax1 = divider.append_axes('bottom', size = '3%', pad = 0.20)
                cax2 = divider.append_axes('bottom', size = '3%', pad = 0.40)
                cbar1 = plt.colorbar(orientation = 'horizontal', cax = cax1,
                                     ticks  = ticks)
                cbar1.set_label('CAPE [J/kg]')
                cbar1.ax.set_xticklabels
            else: 
                divider = make_axes_locatable(ax)
                cax1 = divider.append_axes('bottom', size = '3%', pad = 0.20)
                cbar1 = plt.colorbar(orientation = 'horizontal', cax = cax1)
                plt.gcf().delaxes(cax1)
            plt.sca(ax)
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
            if cbar:
                cbar2 = plt.colorbar(orientation = 'horizontal', cax = cax2)
                cbar2.set_label('Precipitation [mm/h]')
        elif variable == 'THETA':
            field = smoothfield(field, 5)
            #levels = np.arange(200, 350, 2)
            #CS = plt.contourf(X, Y, field, colors = 'black', levels = levels, 
                             #zorder = 0.4, alpha = 1, linewidths = 1.5, 
                             #linestyles = '-')
            #plt.clabel(CS, inline = 1, fontsize = 7, fmt = "%.0f")
            CS = plt.contourf(X, Y, field, cmap = 'spectral', 
                             zorder = 0.4, alpha = 0.5)
            #if cbar:
            #cbar = plt.colorbar(orientation = 'horizontal', pad = 0.05, 
                                #shrink = 0.55)
            #cbar.ax.set_aspect(0.02)
    
    elif setting == 3:
        print 'Setting 3'
        if variable == "PMSL":   # Surface pressure
            field = smoothfield(field, 3)/100
            levels = list(np.arange(900, 1100, 5))
            CS = m.contour(X, Y, field, levels = levels, colors = 'green', 
                            linewidths = 0.5, zorder = 0.5)
            plt.clabel(CS, fontsize = 7, inline = 1, fmt = "%.0f")
            levels = list(np.arange(900, 1100, 1))
            CS = m.contour(X, Y, field, levels = levels, colors = 'black', 
                            linewidths = 0.2, zorder = 0.5)
    
    elif setting == 4:
        print 'Setting 4'
        if variable == "PMSL":   # Surface pressure
            field = smoothfield(field, 3)/100
            levels = list(np.arange(900, 1100, 5))
            CS = m.contour(X, Y, field, levels = levels, colors = 'green', 
                            linewidths = 1.5, zorder = 0.5)
            #plt.clabel(CS, fontsize = 7, inline = 1, fmt = "%.0f")
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
            if cbar:
                ax = plt.gca()
                divider = make_axes_locatable(ax)
                cax1 = divider.append_axes('bottom', size = '3%', pad = 0.30)
                cax2 = divider.append_axes('bottom', size = '3%', pad = 0.60)
                plt.sca(ax)
                cbar2 = plt.colorbar(orientation = 'horizontal', cax = cax1)
                
                cbar2.set_label('Precipitation [mm/h]')
    
    else:
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
        elif variable == 'THETA':
            plt.contourf(X, Y, field, cmap = 'hot_r') 
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
            CS = m.contour(X, Y, field, levels = levels, colors = 'k', 
                            linewidths = 1, zorder = 0.5, alpha = 0.5)
            #plt.clabel(CS, fontsize = 7, inline = 1, fmt = "%.0f")
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
        elif variable == 'W':
            print 'hi'
            levels = list(np.linspace(-7, 7, 100))
            plt.contourf(X, Y, field, cmap = plt.get_cmap('gist_ncar'), 
                         zorder = 0.45, levels = levels)
            cbar = plt.colorbar(shrink = 0.7)
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
    
    # Filter for zero
    x = x[p != 0]
    y = y[p != 0]
    p = p[p != 0]
    # Filter for nan
    x = x[np.isfinite(p)]
    y = y[np.isfinite(p)]
    p = p[np.isfinite(p)]

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
    lc = col.LineCollection(segments, cmap = cmap, norm = norm, alpha = 0.5, zorder = 100)
    lc.set_array(p)
    lc.set_linewidth(linewidth)
    plt.gca().add_collection(lc)
    
    
def smoothfield(field, sigma = 8):
    """
    Smoothes given field
    """
    newfield = ndi.filters.gaussian_filter(field, sigma)
    return newfield
        
    
    
    
    
    
    
    
