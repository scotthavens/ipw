'''
Created on Apr 13, 2015

@author: scott
'''

__version__ = '0.0.1'


import numpy as np
# from ipw import ipw
import ipw
import netCDF4 as nc
import os, math, operator
import progressbar
import warnings
# import sys
# import matplotlib.pyplot as plt


class netcdf:
    '''
    netcdf represents a netCDF file that can be added to or write out IPW files
    '''


    def __init__(self, fileName=None, rw='r', clobber=False):
        '''
        Args:
        fileName - file name to the existing or new netCDF file
        rw (default 'r'/read) - how to open the file 'r','w','r+'
            'r' - read only
            'w' - !!WARNING!! Will overwrite the existing file
            'r+' - append to the file
        '''
        
        self.netcdf = None
        self.dimNames = None
        self.dimSizes = None
        self.variables = None
        self.ipw = []
        
        if fileName is not None:
        
            # read in the netCDF file
            self.fileName = fileName     # file to load           
            
            # netCDF object
            self.netcdf = nc.Dataset(fileName, rw, clobber)
            
            # if it exists get some information about the netCDF file
            if os.path.isfile(fileName):
                       
                # get the variables
                self.variables = self.netcdf.variables
        
                # dimensions and sizes
#                 self.dimNames = self.variables.dimensions
#                 self.dimSizes = self.variables.shape

        else:
            # create empty class to load into
            self.fileName = None
            
            
#     def __str__(self):
#             
#         return self.netcdf.__str__()



    def netcdf2ipw(self, var, timeStep=0, outPrefix='netcdf_', outDir='./data/'):
        '''
        Take the netcdf and convert to IPW images, multiple var's will output
        all the variables to the same file. This does not write the images
                
        *** ASSUMPTION: netCDF orgin is the lower left corner ***
        *** Can only write one variable at a time ***
                
        Args:
        var - string of a 3D variable to convert to IPW
        timeStep (default 0) - which time step to convert, single value or array of values
        outPrefix (default 'netcdf') - prefix to save images
        outDir (default './data') - directory to save images 
        
        '''
        
        # check if var is in variable list
        if type(var) != str:
            raise ValueError('Variable is not specified as a string')
            
        if var not in self.variables:
            raise ValueError('Variable ( ' + var + ' ) not in netCDF file')
            
        # check the time step
        if type(timeStep) == int:
            timeStep = [timeStep]

        # check if the direcotry exists
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        
        pbar = progressbar.ProgressBar(len(timeStep)).start()
        j = 0
        for time in timeStep:
            
            j += 1
            pbar.update(j)
            
            # create the file name to save the file to
            fileName = os.path.join(outDir,outPrefix + '%04i' % time)
            
            # create an IPW instance
            i = ipw.IPW()
            i.fname = fileName
            
            # add a band for each variable
#             for v in var:
                               
            # get the data, this assumes that the orgin of the netcdf data is
            # in the lower left corner
            data = np.flipud(self.netcdf.variables[var][time,:,:])
            
            # add a new band
            i.new_band(data)
                        
            self.ipw.append(i)
            i = None
            
        pbar.finish()
        
    def add_geo_hdr(self, coordinates, d, units, csys, band='all'):
        '''
        Write a geo header to all the IPW files in netcdf
        Wrapper to call IPW.add_geo_hdr
        '''    
        
        for i in range(len(self.ipw)):
            self.ipw[i].add_geo_hdr(coordinates, d, units, csys, band)
        
        
    def write_ipw(self, nbits=8):
        '''
        Write the self.ipw's to IPW files, wrapper for IPW.write()
        Having a seperate set allows the user to make modifications to the bands
        like adding geo headers before writing           
        '''
        
        pbar = progressbar.ProgressBar(len(self.ipw)).start()
        
        for j,i in enumerate(self.ipw):
            pbar.update(j+1)
            i.write(i.fname, nbits)
            
        pbar.finish()
        
    def ipw2netcdf(self, variable, fileDir):
        '''
        Take a directory of IPW images, read in, and add to a netcdf
        '''
        
        
        
    def add_image(self, image, varName, timeStep=None, band=0):
        '''
        Add a sinle band from an IPW image specified by image to the netcdf
        
        Args:
        image - IPW image location
        varName - variable to assign data to
        timeStep (default=None) - timeStep to enter data into
            if timeStep=None, add_image expects a 2D variable
        band (default=0) - if a multi-band image, add this band
        '''
        
        # check to make sure that the image exists
        if not os.path.isfile(image):
            raise IOError('Could not find file -- %s' % image)
        
        # get the IPW file
        i = ipw.IPW(image)
        
        # add the data to the netCDF file
        if timeStep is not None:
            self.netcdf.variables[varName][timeStep,:] =  np.flipud(i.bands[band].data)
        else:
            self.netcdf.variables[varName][:] =  np.flipud(i.bands[band].data)
        
    
    def add_image_directory(self, dirName, prefix, varName, band=0):
        '''
        Read in all the ipw images in a direcotry
        
        Args:
        dirName - direcotry name
        prefix - prefix to the ipw images
        varName - variable name to load into
        band (default=0) - if a multi-band image, add this band
        '''  
    
    
    
        
        
    def add_variable(self, varName, dimensions, datatype='f'):
        '''
        Add a variable if it doesn't exist to the netCDF file with dimensions already defined
        This allows the user to repeatidly call this method with having
        netCDF4 raise an execption if the variable exists
        
        Args:
        varName - name of variable
        dimensions - tuple of dimension names
        datatype (default='float') - type of data to store
        '''
        
        # check if the variable already exists
        if varName in self.netcdf.variables.keys():
            warnings.warn('Variable -- %s -- is already in the netCDF file' % varName)
        
        else:
            
            # determine the chunking shape
#             chunk = self.chunk_shape_3D(dimensions, 8, 8096)
            
            # create the variable
            self.netcdf.createVariable(varName, datatype, dimensions)#, chunksizes=[6,100,100])
     
     
     
    def close(self):
        '''
        Close the netCDF file
        '''
        
        self.netcdf.close()
    
    #===========================================================================
    # The following came from:
    # http://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_choosing_shapes
    #===========================================================================
    def binlist(self, n, width=0):
        '''
        Return list of bits that represent a non-negative integer.
        n      -- non-negative integer
        width  -- number of bits in returned zero-filled list (default 0)
        '''
        
        return map(int, list(bin(n)[2:].zfill(width)))

    def numVals(self, shape):
        '''Return number of values in chunk of specified shape, given by a list of dimension lengths.
    
        shape -- list of variable dimension sizes'''
        if(len(shape) == 0):
            return 1
        return reduce(operator.mul, shape)
    
    def perturbShape(self, shape, onbits):
        '''Return shape perturbed by adding 1 to elements corresponding to 1 bits in onbits
    
        shape  -- list of variable dimension sizes
        onbits -- non-negative integer less than 2**len(shape)
        '''
        return map(sum, zip(shape, self.binlist(onbits, len(shape))))
    
    def chunk_shape_3D(self, varShape, valSize=4, chunkSize=4096):
        '''
        Return a 'good shape' for a 3D variable, assuming balanced 1D, 2D access
    
        varShape  -- length 3 list of variable dimension sizes
        chunkSize -- maximum chunksize desired, in bytes (default 4096)
        valSize   -- size of each data value, in bytes (default 4)
    
        Returns integer chunk lengths of a chunk shape that provides
        balanced access of 1D subsets and 2D subsets of a netCDF or HDF5
        variable var with shape (T, X, Y), where the 1D subsets are of the
        form var[:,x,y] and the 2D slices are of the form var[t,:,:],
        typically 1D time series and 2D spatial slices.  'Good shape' for
        chunks means that the number of chunks accessed to read either
        kind of 1D or 2D subset is approximately equal, and the size of
        each chunk (uncompressed) is no more than chunkSize, which is
        often a disk block size.
        '''
    
        rank = 3  # this is a special case of n-dimensional function chunk_shape
        chunkVals = chunkSize / float(valSize) # ideal number of values in a chunk
        numChunks  = varShape[0]*varShape[1]*varShape[2] / chunkVals # ideal number of chunks
        axisChunks = numChunks ** 0.25 # ideal number of chunks along each 2D axis
        cFloor = [] # will be first estimate of good chunk shape
        # cFloor  = [varShape[0] // axisChunks**2, varShape[1] // axisChunks, varShape[2] // axisChunks]
        # except that each chunk shape dimension must be at least 1
        # chunkDim = max(1.0, varShape[0] // axisChunks**2)
        if varShape[0] / axisChunks**2 < 1.0:
            chunkDim = 1.0
            axisChunks = axisChunks / math.sqrt(varShape[0]/axisChunks**2)
        else:
            chunkDim = varShape[0] // axisChunks**2
        cFloor.append(chunkDim)
        prod = 1.0  # factor to increase other dims if some must be increased to 1.0
        for i in range(1, rank):
            if varShape[i] / axisChunks < 1.0:
                prod *= axisChunks / varShape[i]
        for i in range(1, rank):
            if varShape[i] / axisChunks < 1.0:
                chunkDim = 1.0
            else:
                chunkDim = (prod*varShape[i]) // axisChunks
            cFloor.append(chunkDim)
    
        # cFloor is typically too small, (numVals(cFloor) < chunkSize)
        # Adding 1 to each shape dim results in chunks that are too large,
        # (numVals(cCeil) > chunkSize).  Want to just add 1 to some of the
        # axes to get as close as possible to chunkSize without exceeding
        # it.  Here we use brute force, compute numVals(cCand) for all
        # 2**rank candidates and return the one closest to chunkSize
        # without exceeding it.
        bestChunkSize = 0
        cBest = cFloor
        for i in range(8):
            # cCand = map(sum,zip(cFloor, binlist(i, rank)))
            cCand = self.perturbShape(cFloor, i)
            thisChunkSize = valSize * self.numVals(cCand)
            if bestChunkSize < thisChunkSize <= chunkSize:
                bestChunkSize = thisChunkSize
                cBest = list(cCand) # make a copy of best candidate so far
        return map(int, cBest)

