'''
Created May 14, 2015

Convert IPW files in a directory to a single netcdf file
This script serves a template and may get incorporated into
the netcdf module

@author: Scott Havens
'''

import numpy as np
import matplotlib.pyplot as plt
from isnobal import ipw, netcdf
import progressbar
import glob

from datetime import datetime
startTime = datetime.now()

#===============================================================================
# IPW input directory
#===============================================================================
ipwDir = '/Volumes/kormos-disk1/air_temp/wy2013/*.ipw'

ny = 3011
nx = 1602

#===============================================================================
# NetCDF Inputs and add a new variable VP
#===============================================================================

netcdfFile = '/Volumes/kormos-disk1/air_temp/air_temp.nc'
dimensions = ('time','y','x',)

n = netcdf.netcdf(netcdfFile, 'w', True)

# create the dimensions
n.netcdf.createDimension('time',None)
n.netcdf.createDimension('y',ny)
n.netcdf.createDimension('x',nx)

# create the variable and add some attributes
varName = 'air_temp'
n.add_variable(varName, dimensions)
setattr(n.netcdf.variables[varName], 'units', 'C')
setattr(n.netcdf.variables[varName], 'description', 'Air temperature')


#===============================================================================
# Get all files in the directory, open ipw file, and add to netCDF
#===============================================================================

# get all the files in the directory
d = glob.glob(ipwDir)

pbar = progressbar.ProgressBar(len(d)).start()
j = 0

for file in d:
        
    # Read the IPW file
    i = ipw.IPW(file)
        
    # output to the netcdf file
    # ensure to check that flipud is the correct thing to do
    n.netcdf.variables['air_temp'][j,:] = np.flipud(i.bands[0].data)
        
    j += 1
    pbar.update(j)
 
 
pbar.finish()
    
n.close()
print datetime.now() - startTime
