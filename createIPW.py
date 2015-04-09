'''
Created April 7, 2015

Wrapper script to specify a matrix and export it into an IPW file

@author: Scott Havens
'''



import numpy as np
from ipw import ipw as i
import netCDF4 as nc
import matplotlib.pyplot as plt


#===============================================================================
# Setup the input values for the file to grab 
#===============================================================================

fileName = 'test.nc'    # file to load
var = 'GLW'              # variable to load
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


#===============================================================================
# plot the variable
#===============================================================================

# plt.imshow(tmp)
# plt.colorbar()
# plt.show()

#===============================================================================
# Begin to create the IPW file 
#===============================================================================

'''
To add a new band to the file
- new_band(data, bits)
-- get nsamps and nlines from the data, ensure that they are the same if bands
-- already exist
-- fill in some of the basic header information like bits and bytes (basic_image_

- add_geo_hdr(band number, geo data)
-- adds the geo data to the given band

 
'''

ipw = i.IPW();
ipw.new_band(100, 100)

ipw.bands = 1

print(ipw)

#             self._data_frame = None
#             self.input_file = None
#             self.file_type = None
#             self.header_dict = None
#             self.binary_data = None
#             self.bands = None
#             self.nonglobal_bands = None
#             self.geotransform = None
#             self.start_datetime = None
#             self.end_datetime = None

# Close the file
f.close()









