'''
Created on Apr 13, 2015

@author: scott
'''

__version__ = '0.0.1'


import numpy as np
# from ipw import ipw
import ipw
import netCDF4 as nc
import os
import progressbar
import warnings
# import sys
# import matplotlib.pyplot as plt


class netcdf:
    '''
    netcdf represents a netCDF file that can be added to or write out IPW files
    '''


    def __init__(self, fileName=None, rw='r'):
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
            self.netcdf = nc.Dataset(fileName, rw, clobber=False)
            
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
            # create the variable
            self.netcdf.createVariable(varName, datatype, dimensions)
     
     
     
    def close(self):
        '''
        Close the netCDF file
        '''
        
        self.netcdf.close()
