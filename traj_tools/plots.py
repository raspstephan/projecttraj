"""
Plotting submodule
Submodule of traj_tools
"""


from . import common
import filters
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
        plt.clf()


def draw_xy(varlist, filelist, idlist, cfile, rfiles, pfiles, savename = False, 
            realcoord = True, pollon = None, pollat = None, xlim = None, 
            ylim = None):
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
    realcoord : bool
      If True, real coordinates are used. pollon and pollat must be given.
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
    tmpt = filelist[i].split('/')[-1]
    tmpt = int(tmpt.split('_')[1].lstrip('t'))
    cosmoind = int(tmpt/outint)
    
    # Set up figure
    fig = plt.figure(figsize = (12,8))
    ax = plt.gca()   
    ax.set_aspect('equal')
    basemap(cfile)

    # Plotting all contour fields
    for i in range(len(VarList)):   # Plotting all contour fields
        try:
            Contour(pfiles, varlist[i], cosmoind)
        except:
            pass
        else:
            Contour(rfiles, varlist[i], cosmoind)
    
    # Plot trajectories
    for i in range(len(filelist)):
        rootgrp = nc.Dataset(filelist[i], 'r')
        lonmat = rootgrp.variable['longitude'][:, :]
        latmat = rootgrp.variable['latitude'][:, :]
        pmat = rootgrp.variable['P'][:, :]
    
        if realcoord:
            lonmat[:, :] += (180 - pollon)
            latmat[:, :] += (90 - pollat)

        for j in idlist[i]:
            single_traj(lonmat[:, j], latmat[:, j], pmat[:, j])
    
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
    if savebase != False:
        print "Saving figure as", savebase
        plt.savefig(savebase, dpi = 400)
        plt.close('all')
        plt.clf()
    
            
            
            
            
            
            

#######################################
### Front End Plotting Functions
#######################################


def DrawHist(Criterion = 600, Array = None, IndMatrix = None, TraceInd = 7, SaveBase = False, XAxis = "Minutes"):
    """
    Draws a histogram of ascent times
    """
    pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile, dt = LoadCaseSpec()
    if Array == None:
        Array = filters.MinXMatrix(FileList, TraceInd, Criterion, IndMatrix = WCBIndM, Flat = True)
    if XAxis == "Minutes":   # Convert x axis into minutes, for dt = 5
        Array = list(np.array(Array) * 5)
    if XAxis == "Hours":   # Convert to hours
        Array = list(np.array(Array) * 5 / 60)
        
    # Plotting parameters
    fig = plt.figure()
    plt.hist(Array, bins = 150)
    plt.title("Total number of trajectories:" + str(len(Array)))
    plt.xlabel(XAxis)
    plt.ylabel("Number of trajectories")
    plt.xlim(0,5000)
    
    
    if SaveBase != False:
        savename = SaveBase
        plt.savefig(savename, dpi = 400)

        
def DrawHistStarts(Criterion = 600, TraceInd = 7, SaveBase = False, XAxis = "Minutes"):
    """
    Draws Histograms of all starting times
    """
    pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile, dt = LoadCaseSpec()
    AscMatList = filters.MinXMatrix(FileList, TraceInd, Criterion, StartInd = True, IndMatrix = WCBIndM, Flat = True)
    StartList = filters.StartIndList()
    StartDate = datetime(2012, 10, 13, 00, 00)
    for i in range(len(AscMatList)):
        print "Plotting %i of %i" % (i+1, len(AscMatList))
        DrawHist(Array = AscMatList[i], SaveBase = SaveBase + (StartDate + timedelta(minutes = dt * (StartList[i] + TrjOffset))).isoformat("_"), XAxis = XAxis)
        print "Plotted:", SaveBase + (StartDate + timedelta(minutes = dt * (StartList[i] + TrjOffset))).isoformat("_")


def DrawXYStarts(VarList, IndMatrix = None, RealCoord = True, SaveBase = False, linewidth = 0.7):
    """
    Draws XY plots of all starting times
    """
    pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile, dt = LoadCaseSpec()
    StartIndList = filters.StartIndList()
    print "Total:", len(StartIndList)
    StartDate = datetime(2012, 10, 13, 00, 00)
    for i in range(len(StartIndList)):
        print "Plotting:", i+1
        DrawXYSingle(VarList, StartIndList[i], IndMatrix = IndMatrix, RealCoord = RealCoord, SaveBase = SaveBase + (StartDate + timedelta(minutes = dt * (StartIndList[i] + TrjOffset))).isoformat("_"), linewidth = linewidth)   


def DrawXYInterval(VarList, StartInd, Interval, IndMatrix = None, RealCoord = True, SaveBase = False):
    """
    Draws XY plots in intervals
    """
    pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile = common.LoadCaseSpec(bquiet = True)
    nplots = np.floor((len(pfiles-TrjOffset)*5/60)/Interval)
    print "Total:", nplots
    for i in range(int(nplots)):
        print "Plotting:", i+1, "t =", i*Interval
        DrawXYSingle(VarList, StartInd, IndMatrix = IndMatrix, t = (i*Interval), RealCoord = RealCoord, SaveBase = SaveBase)
    
    
def DrawXYSingle(VarList, StartInd, IndMatrix = None, t = "end", RealCoord = True, SaveBase = False, linewidth = 0.7):
    """ 
    Draws single XY plots
    rvar, pvar are lists with desired variable names
    eg rvar = ["PMSL", "TOT_PREC_S"]
    Time t in hours from start of first trajectory
    """
    # Load case specifics
    pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile, dt = common.LoadCaseSpec()
    if IndMatrix == None:
        IndMatrix = WCBIndM
    
    # Convert t in index
    if t == "end":   # Plot synoptic situation at trajectory start
        COSMOind = TrjOffset + StartInd
    else:
        trjind = t * 60 / 5
        COSMOind = trjind + TrjOffset
        
    # Plotting
    fig = plt.figure(figsize = (12,8))
    ax = plt.gca()   
    ax.set_aspect('equal')   # Keep aspect ratio
    BaseMap(cfile)
    
    for i in range(len(VarList)):   # Plotting all contour fields
        if VarList[i] in pvar:
            Contour(pfiles, VarList[i], COSMOind)
        elif VarList[i] in rvar:
            Contour(rfiles, VarList[i], COSMOind)
        else:
            print rvar, pvar, VarList[i]
            raise ValueError ('Wrong input')
    
    # Plot trajectories
    for i in range(len(IndMatrix)):
        M = filters.OpenSaveFile(FileList[i])[1]
        if RealCoord:
            M[:, 0, :] += (180 - 165)
            M[:, 1, :] += (90 - 35)
            
        for j in IndMatrix[i]:
            StartPosTmp = filters.StartPos(M[j])[0]
            if StartPosTmp == StartInd:
                # Plotting
                if t == "end":
                    XYPlot(M[j, :, StartPosTmp:], linewidth = linewidth)
                else:
                    XYPlot(M[j, :, StartPosTmp:(StartPosTmp+trjind)], linewidth = linewidth)
            elif StartPosTmp > StartInd:
                break
        else: 
            continue
        break   # Break out of both loops
    
    # Set plot properties
    if RealCoord:
        plt.xlim(-14, 44)
        plt.ylim(33, 76)
    else:
        plt.xlim(-29, 29)
        plt.ylim(-22, 21)
    cb = fig.colorbar(lc, shrink = 0.7)
    cb.set_label('p')
    cb.ax.invert_yaxis()
    plt.tight_layout()
    
    # Save Plot
    if SaveBase != False:
        print "Saving figure as", SaveBase
        plt.savefig(SaveBase, dpi = 400)
        plt.close('all')
    
    
    
    
    




#########################################
### Back End Plotting Functions
#########################################


def basemap(cfile, RealCoord = True, pollon = -165., pollat = 35.):
    """
    Draws Land-Sea Map
    """
    dx = 0.025   # Assumes COSMO 2.2 resolution
    field = pwg.getfield(cfile, "FR_LAND_S", invfile = False)
    
    # Setting up grid
    ny, nx = field.shape
    x = np.linspace(-dx * (nx -1) / 2, dx * (nx -1) / 2, nx)
    #y = np.linspace(-dx * (ny -1) / 2, dx * (ny -1) / 2, ny)   Domain is slightly off center
    y = np.linspace(-22, 21, ny)
    if RealCoord:
        roteqlon = 180 + pollon
        roteqlat = 90 - pollat
        x += roteqlon
        y += roteqlat
    X, Y = np.meshgrid(x, y)
    
    # Plot and filter field
    field[field != 0] = 1.
    plt.contourf(X, Y, field, levels = [0.5, 2], colors = ["0.85"])
    #plt.contour(X, Y, field, 2, colors = "0.5")
    
    del field
    
    
def contour(files, Variable, COSMOind, RealCoord = True, pollon = -165., pollat = 35.):
    """
    Draws contour plot of Variable
    """
    dt = 5. * 60
    dx = 0.025
    
    #Get field
    if Variable == "TOT_PREC_S":
        field1 = pwg.getfield(files[COSMOind - 1], "TOT_PREC_S", invfile = False)
        field2 = pwg.getfield(files[COSMOind + 1], "TOT_PREC_S", invfile = False)
        field = (field2 - field1)/(dt*2)*3600
    else:
        field = pwg.getfield(files[COSMOind], Variable, invfile = False)
     
    # Set up grid
    ny, nx = field.shape
    x = np.linspace(-dx * (nx -1) / 2, dx * (nx -1) / 2, nx)
    #y = np.linspace(-dx * (ny -1) / 2, dx * (ny -1) / 2, ny)   Domain is slightly off center
    y = np.linspace(-22, 21, ny)
    if RealCoord:
        roteqlon = 180 + pollon
        roteqlat = 90 - pollat
        x += roteqlon
        y += roteqlat
    X, Y = np.meshgrid(x, y)
    
    if Variable == "FI":
        field = smoothfield(field, 8)
        levels = list(np.arange(400,600,8))   # Needs smoothing?
        plt.contour(X, Y, field/100, levels = levels, colors = "k", linewidths = 2)
    elif Variable == "T":
        field = smoothfield(field, 8)
        plt.contourf(X, Y, field, alpha = 0.5)
        plt.colorbar()
        #plt.contour(X, Y, field, levels = list(np.arange(150, 350, 4)), colors = "grey", linewidths = 2)
    elif Variable == "TOT_PREC_S":
        cmPrec =( (0    , 0.627 , 1    ),
                  (0.137, 0.235 , 0.98 ),
                  (0.1  , 0.1   , 0.784),
                  (0.392, 0     , 0.627),
                  (0.784, 0     , 0.627),
                  (1    , 0.3   , 0.9  ) )   # Tobi's colormap
        levels = [0.1, 0.3, 1, 3, 10, 100]
        plt.contourf(X, Y, field, levels, colors=cmPrec, extend='max', alpha = 0.8, zorder = 10)
        plt.colorbar(shrink = 0.7)
    elif Variable == "PMSL":
        field = smoothfield(field, 8)/100
        levels = list(np.arange(900, 1100, 5))
        CS = plt.contour(X, Y, field, levels = levels, colors = 'k', linewidths = 1, zorder = 9, alpha = 0.5)
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
    x = A[0][A[7] != 0]
    y = A[1][A[7] != 0]
    p = A[7][A[7] != 0]

    points = np.array([x,y]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
    lc = col.LineCollection(segments, cmap=plt.get_cmap('Spectral'), norm=plt.Normalize(100, 1000), alpha = 0.8)
    lc.set_array(p)
    lc.set_linewidth(linewidth)
    plt.gca().add_collection(lc)
    
    
def smoothfield(field, sigma = 8):
    """
    Smoothes given field
    """
    newfield = ndi.filters.gaussian_filter(field, sigma)
    return newfield
        
    
    
    
    
    
    
    
