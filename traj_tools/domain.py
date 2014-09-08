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

def convert_domain(pollon, pollat, lonlim, latlim, savename = None):
    """
    Projects COSMO coordinates onto a real grid
    
    Parameters
    ----------
    pollon : float
      Longitude of rotated pole
    pollat : float
      Latitude of rotated pole
    lonlim : tuple
      x-direction limits of rotated projection
    latlim : tuple
      y-direction limits of rotated projection
    savename : str
      Path where output plot will be saved
      If None, function will return fig object
      
    Returns
    -------
    fig : plt.fig object
      If savename is None
    
    """
    
    # Create boundaries 
    
    
    
    # Convert into real coordinates
    
    
    # Draw basemap
    m = Basemap(projection = 'ortho', lon_0 = pollon, lat_0 = pollat, 
                resolution = 'l')  # NOTE: Use different coordinates
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
    plt.show()
    
    # Draw converted boundaries
    


# Plot test map if called directly 
if __name__ == '__main__':
    convert_domain(0, 35, None, None, None)
