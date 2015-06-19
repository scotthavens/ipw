'''
Created April 17, 2015

Calculate the albedo

@author: Scott Havens
'''

import numpy as np
import matplotlib.pyplot as plt
from isnobal import ipw, netcdf, rad
import progressbar
# import multiprocessing as mp

from datetime import datetime
startTime = datetime.now()

#===============================================================================
# NetCDF Inputs and add a new variable VP
#===============================================================================

dataFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/wrf_data/Regrid/BRB_100m_Regrid.nc'
workFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/wrf_data/Regrid/BRB_working.nc'
dimensions = ('Time','south_north','west_east',)

data = netcdf.netcdf(dataFile, 'r+')
work = netcdf.netcdf(workFile, 'r+')

# create the variable and add some attributes
# varName = 'storm_days'
# work.add_variable(varName, dimensions)
# setattr(work.netcdf.variables[varName], 'units', 'time step')
# setattr(work.netcdf.variables[varName], 'description', 'Time steps since last storm')

# varName = 'storm_precip'
# work.add_variable(varName, dimensions)
# setattr(work.netcdf.variables[varName], 'units', 'mm')
# setattr(work.netcdf.variables[varName], 'description', 'Accumulated storm precip')
 
#===============================================================================
# Calculate % snow and snow density
#===============================================================================

timeStep = range(115,120)        # timestep to load

# pbar = progressbar.ProgressBar(len(timeStep)).start()
j = 0
zone = 7*60     # minutes west of GMT

 # lat/lon of the middle of the BRB
lat = 43 + 51/float(60) + 50/float(3600)
lon = -115 + 20/float(60) + 00/float(3600)

# constants for albeda
gsize = 150
maxgsz = 500
dirt = 2.0

for t in timeStep:
        
    # get the current time step date
    dt = data.netcdf.variables['Times'][t,:].tostring()
    date = datetime.strptime(dt,'%Y-%m-%d_%H:%M:%S')
    
   # calculate the sun angle
    cosz, azimuth = rad.sunang(date, lat, lon)
        
    # get the storm days value
    p = data.netcdf.variables['storm_days'][t,:]
    
#     telapsed = 100.0 + t

    # calculate the albedo
    alb_v, alb_ir = rad.albedo(p, cosz, gsize, maxgsz, dirt)

    
    
#     print(cosz)
#     print(alb_v, alb_ir)
    

      # output to the netcdf file
#     work.netcdf.variables['storm_days'][t,:] = stormDays

    plt.figure(1)
    plt.subplot(221)
    plt.imshow(p)
    plt.colorbar()
        
#     plt.subplot(222)
#     plt.imshow(ps)
#     plt.colorbar()
        
    plt.subplot(223)
    plt.imshow(alb_v)
    plt.colorbar()
        
    plt.subplot(224)
    plt.imshow(alb_ir)
    plt.colorbar()
        
    plt.show()
      

    j += 1
#     pbar.update(j)
 
 
# pbar.finish()

data.close()
work.close()
print datetime.now() - startTime
