'''
Created on Apr 7, 2015

Playing around with loading in a netCDF file and getting infromation about it.
End goal is to take a netCDF file and convert a file, variable, time step into a IPW file

This will only deal with variables that are 3D and not the 4D variables (yet)

@author: scott
'''

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


#===============================================================================
# Setup the input values for the file to grab 
#===============================================================================

fileName = 'test.nc'    # file to load
var = 'GLW'              # variable to load
timeStep = 9            # timestep to load
outDir = './'           # direcotry to put the new file
outPrefix = 'ta'        # prefix to use for the IPW file
offset = -273.15        # offset to apply to all the data before being written

#===============================================================================
# Create the dataset object and access the variable
#===============================================================================

# list all the variables
f = nc.Dataset(fileName)

# print(f.variables.keys())

v = f.variables[var]

print(v)

#===============================================================================
# List the dimensions for the variable
#===============================================================================

# for d in f.dimensions.items(): # lists all the demsions of the file

dimNames = v.dimensions
dimSizes = v.shape

print(dimNames,dimSizes)

#===============================================================================
# Get the variable at the specified time
# flipud to ensure that (0,0) is in the upper left
#===============================================================================

tmp = np.flipud(v[timeStep,:,:] + offset)


#===============================================================================
# plot the variable
#===============================================================================

plt.imshow(tmp)
plt.colorbar()
plt.show()




# Close the file
f.close()



