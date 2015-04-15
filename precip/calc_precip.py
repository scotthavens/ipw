'''
Created April 15, 2015

Load the distributed precipitation and do the following:
1. Calculate time since last storm
2. Albedo

@author: Scott Havens
'''

import numpy as np
import matplotlib.pyplot as plt
import isnobal as sno


#===============================================================================
# Load the precipitation data
#===============================================================================

fileName = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/wrf_data/Regrid/BRB_100m_Regrid.nc'
timeStep = range(80,81)        # timestep to load

# open the netcdf file
n = netcdf.netcdf(fileName)

for t in timeStep:
    
    # get the precip value
    p = n.netcdf.variables['RAINNC'][t,:]
    p = np.flipud(p)
    
    # get the dew point value
    d = n.netcdf.variables['dwpt'][t,:]
    d = np.flipud(d)

    # determine the precip phase


    
    plt.figure(1)
    plt.subplot(121)
    plt.imshow(p)
    plt.colorbar()
    
    plt.subplot(122)
    plt.imshow(d)
    plt.colorbar()
    
    plt.show()
    
    
n.close()
    
