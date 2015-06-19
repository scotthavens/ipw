'''
Created April 15, 2015

Wrapper script to add a variable to the netCDF file from a ipw image

@author: Scott Havens
'''

import numpy as np
from ipw import ipw, netcdf
# import ipw
import netCDF4 as nc
import matplotlib.pyplot as plt


#===============================================================================
# Setup the input values for the file to grab 
#===============================================================================

fileName = 'test.nc'    # file to load
# fileName = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/wrf_data/Regrid/'

var = 'test3d'              # variable to add
datatype = 'f'
dimensions = ('Time','south_north','west_east',)

timeStep = 0            # timestep to add
outDir = './data/'           # direcotry of file
outPrefix = 'ta'        # prefix to use for the IPW file

#===============================================================================
# Wrtie data to netcdf file
#===============================================================================
n = netcdf.netcdf(fileName, 'r+')

# create the variable
n.add_variable(var, dimensions)

# add some attributes
setattr(n.netcdf.variables[var], 'units', 'my units')
setattr(n.netcdf.variables[var], 'description', 'my description')

# add an ipw image to the netcdf
image = './data/ta_0000' 
n.add_image(image, var, 2, 0)



n.close()









