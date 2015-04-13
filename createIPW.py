'''
Created April 7, 2015

Wrapper script to specify a matrix and export it into an IPW file

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
fileName = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/wrf_data/Regrid/'

var = 'T2'              # variable to load
timeStep = 9            # timestep to load
outDir = './data/'           # direcotry to put the new file
outPrefix = 'ta'        # prefix to use for the IPW file
offset = -273.15        # offset to apply to all the data before being written

#===============================================================================
# Create the dataset object and access the variable
#===============================================================================

# netCDF object
f = nc.Dataset(fileName)

# get the variable
v = f.variables[var]


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
# print(len(tmp))
# print(len(tmp[0]))
# print(tmp.shape)

# Close the file
f.close()

#===============================================================================
# plot the variable
#===============================================================================

# plt.imshow(tmp)
# plt.colorbar()
# plt.show()

#===============================================================================
# Begin to create the IPW file 
#===============================================================================

i = ipw.IPW();      # initialize a new file
i.new_band(tmp)   # add data to the band
i.new_band(tmp)   # add data to the band

i.add_geo_hdr([4886900, 570100], [-100,100], 'm', 'UTM')

# print(i)



#===============================================================================
# Wrtie netcdf file to data 
#===============================================================================
n = netcdf.netcdf(fileName)
n.netcdf2ipw('T2',0,'ta_')
n.add_geo_hdr([4886900, 570100], [-100,100], 'm', 'UTM')
n.write_ipw(16)



i = ipw.IPW('./data/ta_0000')

plt.imshow(i.bands[0].data)
plt.colorbar()

plt.show()








