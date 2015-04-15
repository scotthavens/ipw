'''
Created on Apr 15, 2015

@author: scott
'''

__version__ = '0.0.1'


import numpy as np
# from isnobal import ipw



def mkprecip(precipitation, temperature):
    '''
    Follows the IPW command mkprecip
    
    temperature         % snow    snow density

          T < -5 C       100%        75 kg/m^3
      -5 C <= T < -3 C   100%       100 kg/m^3
      -3 C <= T < -1.5 C 100%       150 kg/m^3
    -1.5 C <= T < -0.5 C 100%       175 kg/m^3
    -0.5 C <= T < 0 C     75%       200 kg/m^3
       0 C <= T < 0.5 C   25%       250 kg/m^3
     0.5 C <= T            0%         0 kg/m^3
    
    Args:
    precipitation - array of precipitation values [mm]
    temperature - array of temperature values, use dew point temperature
        if available [degree C]
        
    Output:
    - returns the percent snow and estimated snow density
    '''
    
    # convert the inputs to numpy arrays
    precipitation = np.array(precipitation)
    temperature = np.array(temperature)
    
    # create a list from the table above
    t = []
    t.append( {'temp_min': -1e309,  'temp_max': -5,     'snow': 1,    'density':75} )
    t.append( {'temp_min': -5,      'temp_max': -3,     'snow': 1,    'density':100} )
    t.append( {'temp_min': -3,      'temp_max': -15,    'snow': 1,    'density':150} )
    t.append( {'temp_min': -1.5,    'temp_max': -0.5,   'snow': 1,    'density':175} )
    t.append( {'temp_min': -0.5,    'temp_max': 0.0,    'snow': 0.75, 'density':200} )
    t.append( {'temp_min': 0.0,     'temp_max': 0.5,    'snow': 0.25, 'density':250} )
    t.append( {'temp_min': 0.5,     'temp_max': 1e309,  'snow': 0,    'density':0} )
        
    
    # preallocate the percent snow (ps) and snow density (sd)
    ps = np.zeros(precipitation.shape)
    sd = np.zeros(ps.shape)
    
    # if no precipitation return all zeros
    if np.sum(precipitation) == 0:
        return ps, sd
    
    # determine the indicies and allocate based on the table above
    for row in t:
        
        # get the values between the temperature ranges that have precip
        ind = np.where(np.logical_and( \
            temperature >= row['temp_min'], \
            temperature <= row['temp_max']))
        
        # set the percent snow
        ps[ind] = row['snow']
        
        # set the density
        sd[ind] = row['density']
    
    # if there is no precipitation at a pixel, don't report a value
    # this may make isnobal crash, I'm not really sure
    ps[precipitation == 0] = 0
    sd[precipitation == 0] = 0
    
    return ps, sd
    
    
    
    
    
    
    
    
    
    
    
    