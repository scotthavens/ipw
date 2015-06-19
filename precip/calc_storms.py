'''
Created May 29, 2015

Load the distributed precipitation and do the following:
1. Calculate time since last storm

@author: Scott Havens
'''

import numpy as np
import matplotlib.pyplot as plt
from isnobal import ipw, netcdf, precip
import progressbar
# import multiprocessing as mp

from datetime import datetime
startTime = datetime.now()

#===============================================================================
# NetCDF Inputs and add a new variable VP
#===============================================================================

dimensions = ('Time','south_north','west_east',)
nx = 2000
ny = 1500

dataFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/wrf_data/Regrid/BRB_100m_Regrid.nc'
psFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/precip/perc_snow.nc'
data = netcdf.netcdf(dataFile, 'r')
psdata = netcdf.netcdf(psFile, 'r')

# STORM DAYS NETCDF FILE
sdFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/precip/storm_days.nc'
sd = netcdf.netcdf(sdFile, 'w', True)

sd.netcdf.createDimension(dimensions[0],None)
sd.netcdf.createDimension(dimensions[1],ny)
sd.netcdf.createDimension(dimensions[2],nx)

varName = 'storm_days'
sd.add_variable(varName, dimensions)
setattr(sd.netcdf.variables[varName], 'units', 'time step')
setattr(sd.netcdf.variables[varName], 'description', 'Time since last storm')

# STORM PRECIP NETCDF FILE
# spFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/precip/storm_precip.nc'
# sp = netcdf.netcdf(spFile, 'w', True)
# 
# sp.netcdf.createDimension(dimensions[0],None)
# sp.netcdf.createDimension(dimensions[1],ny)
# sp.netcdf.createDimension(dimensions[2],nx)
# 
# varName = 'storm_precip'
# sp.add_variable(varName, dimensions)
# setattr(sp.netcdf.variables[varName], 'units', 'mm')
# setattr(sp.netcdf.variables[varName], 'description', 'Accumulated storm precip')
 
#===============================================================================
# Calculate time since last storm and accumulated storm precip
#===============================================================================

timeStep = range(5947)        # timestep to load

pbar = progressbar.ProgressBar(len(timeStep)).start()
j = 0

junk = data.netcdf.variables['RAINNC'][2,:];
stormDays = np.zeros(junk.shape)
stormPrecip = stormDays
for t in timeStep:
        
    # get the precip value
    p = data.netcdf.variables['RAINNC'][t,:]
    
    # get snow percent
    ps = psdata.netcdf.variables['perc_snow'][t,:]

    # determine the time since last storm
    stormDays, stormPrecip = precip.storms(p, ps, 0.5, 4, stormDays, stormPrecip)

    # output to the netcdf file
    sd.netcdf.variables['storm_days'][t,:] = stormDays
#     sp.netcdf.variables['storm_precip'][t,:] = stormPrecip

#     plt.figure(1)
#     plt.subplot(221)
#     plt.imshow(p)
#     plt.colorbar()
#        
#     plt.subplot(222)
#     plt.imshow(ps)
#     plt.colorbar()
#        
#     plt.subplot(223)
#     plt.imshow(stormDays)
#     plt.colorbar()
#        
#     plt.subplot(224)
#     plt.imshow(stormPrecip)
#     plt.colorbar()
#        
#     plt.show()
#      

    j += 1
    pbar.update(j)
 
 
pbar.finish()

data.close()
sd.close()
sp.close()
print datetime.now() - startTime
