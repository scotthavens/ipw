'''
Created June 04, 2015

Convert a variable from a netcdf file into IPW files
This script serves a template and may get incorporated into
the netcdf module

@author: Scott Havens
'''

import numpy as np
import matplotlib.pyplot as plt
from isnobal import ipw, netcdf
import progressbar

from datetime import datetime
startTime = datetime.now()

#===============================================================================
# IPW output directory
#===============================================================================
ipwDir = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/precip/stormDays/'
ipwPrefix = 'storm_'

var = 'storm_days'  # variable in NetCDF file to convert

u = 4770050
v = 550050
du = -100
dv = 100
units = 'm'
csys = 'UTM'
nx = 2000
ny = 1500

nbits = 8

#===============================================================================
# NetCDF Inputs and add a new variable VP
#===============================================================================

netcdfFile = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/precip/storm_days.nc'
dimensions = ('time','y','x',)

n = netcdf.netcdf(netcdfFile, 'r', True)

timeSteps = n.netcdf.variables[var].dimSizes

#===============================================================================
# Write NetCDF variable to IPW files
#===============================================================================

n.netcdf2ipw(var, timeStep=0, outPrefix='netcdf_', outDir='./data/') 
    
n.close()
print datetime.now() - startTime
