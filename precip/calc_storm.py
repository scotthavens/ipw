'''
Created April 16, 2015

Load in the precip and calculate:
1. The time since last storm
2. Albedo

@author: Scott Havens
'''

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
# varName = 'storm_days'
# n.add_variable(varName, dimensions)
# setattr(n.netcdf.variables[varName], 'units', '%')
# setattr(n.netcdf.variables[varName], 'description', 'Decimal days since storm')
 
#===============================================================================
# Calculate % snow and snow density
#===============================================================================

timeStep = range(21,22)        # timestep to load

pbar = progressbar.ProgressBar(len(timeStep)).start()
j = 0

for t in timeStep:
        
    # get the precip value
    p = n.netcdf.variables['RAINNC'][t,:]
#     p = np.flipud(p)
    
    # get snow percent
    ps = n.netcdf.variables['perc_snow'][t,:]
#     ps = np.flipud(ps)

    # determine the precip phase
#     ps, sd = precip.mkprecip(p, d)


    plt.figure(1)
    plt.subplot(221)
    plt.imshow(p)
    plt.colorbar()
      
    plt.subplot(222)
    plt.imshow(ps)
    plt.colorbar()
      
    plt.subplot(223)
    plt.imshow(ps)
    plt.colorbar()
      
#     plt.subplot(224)
#     plt.imshow(sd)
#     plt.colorbar()
      
    plt.show()


    # output to the netcdf file
#     n.netcdf.variables['perc_snow'][t,:] = ps
#     n.netcdf.variables['snow_density'][t,:] = sd
    
    # write the  4-band precipitation image if there is precip
#     if np.sum(p > 0):
#         i = ipw.IPW()
#         i.new_band(p)
#         i.new_band(ps)
#         i.new_band(sd)
#         i.new_band(d)
#         i.add_geo_hdr([u, v], [du, dv], units, csys)
#         i.write('%s%s%04i' % (outDir, outPrefix, t), nbits)
    
    
    j += 1
    pbar.update(j)
 
 
pbar.finish()


# results = [pool.apply_async(calc_stuff, args=(x,)) for x in timeStep]

    

    
    
n.close()
print datetime.now() - startTime
