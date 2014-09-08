###############################################################################
#
# COSMO Planning Tool - Domain conversion into rotated coordinates
# 
###############################################################################

# Imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Callable function

def cosmo2geo(pollon, pollat, lonlim, latlim, savename = None):
    """
    Projects COSMO coordinates onto a geographical grid.
    NOTE: Projection type is hard coded at the moment!
    
    Parameters
    ----------
    pollon : float
      Longitude of rotated pole
    pollat : float
      Latitude of rotated pole
    lonlim : tuple
      x-direction limits of rotated projection, e.g. (-10, 10)
    latlim : tuple
      y-direction limits of rotated projection
    savename : str
      Path where output plot will be saved
      If None, function will return fig object
    
    """
    
    # Set up figure
    plt.figure(figsize = (10, 10)) 
    plt.text(0.5, 1.1, ('N pole lon/lat = ' + str(pollon) + '/' + str(pollat) + 
                        ' lonlim = ' + str(lonlim) + ' latlim = ' + 
                        str(latlim)), transform = plt.gca().transAxes, 
                        horizontalalignment='center')
        
    # Create boundaries 
    ns = 100   # Number of points on each side
    lon1 = np.linspace(lonlim[0], lonlim[1], ns)
    lat1 = np.linspace(latlim[1], latlim[1], ns)
    lon2 = np.linspace(lonlim[1], lonlim[1], ns)
    lat2 = np.linspace(latlim[1], latlim[0], ns)
    lon3 = np.linspace(lonlim[1], lonlim[0], ns)
    lat3 = np.linspace(latlim[0], latlim[0], ns)
    lon4 = np.linspace(lonlim[0], lonlim[0], ns)
    lat4 = np.linspace(latlim[0], latlim[1], ns)
    rlon = np.concatenate((lon1, lon2, lon3, lon4))
    rlat = np.concatenate((lat1, lat2, lat3, lat4))
    
    # Convert to radians
    rlon = np.deg2rad(rlon)   
    rlat = np.deg2rad(rlat)
    pollon = np.deg2rad(pollon + 180)   # Not sure why and if correct
    pollat = np.deg2rad(pollat)
    
    # Convert into real coordinates
    glon = np.arctan((np.cos(rlat) * np.sin(rlon)) / 
                     (np.sin(pollat) * np.cos(rlat) * np.cos(rlon) - 
                      np.sin(rlat) * np.cos(pollat))) + pollon
    glat = np.arcsin(np.sin(rlat) * np.sin(pollat) + np.cos(rlat) * 
                     np.cos(rlon) * np.cos(pollat))
    
    # Convert back to degrees
    glon = np.rad2deg(glon)   
    glat = np.rad2deg(glat)
    
    
    
    # Draw basemap
    #lonmid = (np.amax(glon) + np.amin(glon)) / 2   # Evaluate midpoints of plot
    #latmid = (np.amax(glat) + np.amin(glat)) / 2
    #dg = 5   # Extra margin around domain
    m = Basemap(projection = 'npstere', lon_0 = 0, lat_0 = 90, boundinglat = 20,
                resolution = 'l')
    m.drawcoastlines()
    #m.drawcountries()
    m.drawparallels(np.arange(0, 105, 15), 
                    labels = [0, 0, 0, 0])
    m.drawmeridians(np.arange(0, 390, 30),
                    labels = [1, 1, 1, 1])
    
    # Convert lat and lon to map projection
    mlon, mlat = m(glon, glat)
    
    # Draw converted boundaries
    plt.plot(mlon, mlat, 'r')
    
    # Save figure
    if savename != None:
        print 'Saving figure as: ', savename
        plt.savefig(savename)
    
    # Print out additional information to console
    ie = (lonlim[1] - lonlim[0]) / 0.025 + 1   # Cosmo grid points
    je = (latlim[1] - latlim[0]) / 0.025 + 1
    print '============================================'
    print 'With dx = 0.025: ie = ', ie, ' , je = ', je
    print 'Total number of grid points: ' , ie * je
    print '============================================'




# Plot test map if called directly 
if __name__ == '__main__':
    cosmo2geo(-170, 40, (-14, 14), (-12, 12), 
                   '/usr/users/stephan.rasp/tmp/latlon.png')
