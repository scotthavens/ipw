'''
Created April 15, 2015

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
# IPW output file
#===============================================================================
outDir = './ppt.4b/'        # direcotry to put the new file
outPrefix = 'ppt.4b_'       # prefix to use for the IPW file

u = 4770050
v = 550050
du = -100
dv = 100
units = 'm'
csys = 'UTM'

nbits = 16

#===============================================================================
# NetCDF Inputs and add a new variable VP
#===============================================================================

netcdfFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/wrf_data/Regrid/BRB_100m_Regrid.nc'
dimensions = ('Time','south_north','west_east',)

n = netcdf.netcdf(netcdfFile, 'r+')

# create the variable and add some attributes
varName = 'storm_days'
n.add_variable(varName, dimensions)
setattr(n.netcdf.variables[varName], 'units', 'time step')
setattr(n.netcdf.variables[varName], 'description', 'Time steps since last storm')

varName = 'storm_precip'
n.add_variable(varName, dimensions)
setattr(n.netcdf.variables[varName], 'units', 'mm')
setattr(n.netcdf.variables[varName], 'description', 'Accumulated storm precip')
 
#===============================================================================
# Calculate % snow and snow density
#===============================================================================

timeStep = range(100)        # timestep to load

pbar = progressbar.ProgressBar(len(timeStep)).start()
j = 0

junk = n.netcdf.variables['RAINNC'][2,:];
stormDays = np.zeros(junk.shape)
stormPrecip = stormDays
for t in timeStep:
        
    # get the precip value
    p = n.netcdf.variables['RAINNC'][t,:]
    
    # get snow percent
    ps = n.netcdf.variables['perc_snow'][t,:]

    # determine the time since last storm
    stormDays, stormPrecip = precip.storms(p, ps, 0.2, 2, stormDays, stormPrecip)

      # output to the netcdf file
    n.netcdf.variables['storm_days'][t,:] = stormDays
    n.netcdf.variables['storm_precip'][t,:] = stormPrecip

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

n.close()
print datetime.now() - startTime
