'''
Created on Apr 17, 2015

@author: scott
'''

__version__ = '0.0.1'


import numpy as np
import datetime as dt
import subprocess as sp
import math
import progressbar
# from isnobal import ipw

IPW = "/Users/ipw/ipw-2.1.0/bin"# IPW executables

# define some albedo parameters
MAXV = 1.0          # vis albedo when gsize = 0
MAXIR = 0.85447     # IR albedo when gsize = 0
IRFAC = -0.02123    # IR decay factor
VFAC = 500.0        # visible decay factor
VZRG = 1.375e-3     # vis zenith increase range factor
IRZRG = 2.0e-3      # ir zenith increase range factor
IRZ0 = 0.1          # ir zenith increase range, gsize=0

def growth(t):
    '''
    Calculate grain size growth
    From IPW albedo > growth
    '''
    
    a = 4.0;
    b = 3.
    c = 2.0
    d = 1.0

    factor = (a+(b*t)+(t*t))/(c+(d*t)+(t*t)) - 1.0

    return(1.0 - factor)

def albedo(telapsed, cosz, gsize, maxgsz, dirt=2):
    '''
    Calculate the abedo, adapted from IPW function albedo
    
    Inputs:
    telapsed - time since last snow storm (decimal days)
    cosz - from sunang
    gsize - gsize is effective grain radius of snow after last storm (mu m)
    maxgsz -  maxgsz is maximum grain radius expected from grain growth (mu m)
    dirt - dirt is effective contamination for adjustment to visible albedo (usually between 1.5-3.0)
    
    Output
    
    
    Created April 17, 2015
    Scott Havens
    '''
    
    telapsed = np.array(telapsed)
    
    # check inputs
    if gsize <= 0 | gsize >= 500:
        raise Exception("unrealistic input: gsize=%i", gsize)

    if (maxgsz <= gsize | maxgsz >= 2000):
        raise Exception("unrealistic input: maxgsz=%i", maxgsz)

    if dirt >= 10:
        raise Exception("unrealistic input: dirt=%i", dirt)
    if dirt <= 1:
        raise Exception("unrealistic input: dirt=%i", dirt)
    
    # set initial grain radii for vis and ir
    radius_ir = math.sqrt(gsize);
    range_ir = math.sqrt(maxgsz) - radius_ir;
    radius_v = math.sqrt(dirt * gsize);
    range_v = math.sqrt(dirt * maxgsz) -radius_v;
    
    #check to see if sun is up
    if (cosz > 0.0):
        # calc grain growth decay factor
        growth_factor = growth(telapsed + 1.0);

        # calc effective gsizes for vis & ir
        gv = radius_v + (range_v * growth_factor);
        gir = radius_ir + (range_ir * growth_factor);

        #calc albedos for cos(z)=1
        alb_v = MAXV - (gv / VFAC);
        alb_ir = MAXIR * np.exp(IRFAC * gir);

        # calculate effect of cos(z)<1

        # adjust diurnal increase range 
        dzv = gv * VZRG;
        dzir = (gir * IRZRG) + IRZ0;

        # calculate albedo
        alb_v += dzv * (1.0 - cosz);
        alb_ir += dzir * (1.0 - cosz);
    
    else:
        alb_v = np.zeros(telapsed.shape);
        alb_ir = np.zeros(telapsed.shape);
    
    return alb_v, alb_ir

def ihorizon(X, Y, Z, azm, searlen=100):
    '''
    Calculate the horizon values for an entire DEM image
    for the desired azimuth
    
    Assumes that the step size is constant
    
    Inputs:
        X - vector of x-coordinates
        Y - vector of y-coordinates
        Z - matrix of elevation data
        azm - azimuth to calculate the horizon at
        searlen - length of search (in pixels) to look for horizon
    
    Outputs:
        CZ - cosine of the local horizonal angles
        H - index along line to the point
    
    20150602 Scott Havens 
    '''
    
    # check inputs
    azm = azm*np.pi/180 # degress to radians
    m,n = Z.shape
    dx = X[2] - X[1]    # DEM step size
    
    
    # preallocate
    CZ = np.zeros(Z.shape)
    H = np.zeros(Z.shape)
    
    # iterate over the DEM
    st = 10000
    pbar = progressbar.ProgressBar(st).start()
    j = 0
    for (y0,x0), value in np.ndenumerate(Z):
        
        # create the endpoint
        x1, y1 = x0 + searlen*np.sin(azm), y0 + searlen*np.cos(azm)
        
        # create the line as indicies to the dem
#         xi, yi = np.linspace(x0, x1, searlen), np.linspace(y0, y1, searlen) 
        xi, yi = _linspace(x0, x1, searlen), _linspace(y0, y1, searlen)
        xi = xi.astype(np.int)  # indicies to the data
        yi = yi.astype(np.int)
        
        # Extract the values along the line
        ind = (xi < n) & (yi < m) # ensure that the values are on the grid
        xi = xi[ind]
        yi = yi[ind]
        xl = X[xi]
        yl = Y[yi]
        
        # distances between points
        di = np.hypot(xl[1:]-xl[:-1],yl[1:]-yl[:-1])    
        di = np.cumsum(di)
        di = np.insert(di,0,0)   # add zero to the front for the current point
        
        # remove any duplicated values
        di, ind = np.unique(di, return_index=True)
        xi = xi[ind]
        yi = yi[ind]
        
        zi = Z[yi, xi]
        
        # calculate the horizon
        CZ[y0,x0], H[y0,x0] = horizon(0,di,zi)
        
        j += 1
        pbar.update(j)
        if j == st:
            break
 
 
    pbar.finish()
    
    return CZ, H
    
    

def horizon(i,x,z):
    '''
    Calculate the horizon point for pixel i on the x,y points
    Written to mimic IPW's horizon function
    
    Inputs:
    i - pixel index on the (x,y) point
    x - horizontal distances for points
    z - elevations for the points
    
    Output:
    cz - cosine of the anlge from the zenith to the pixel's horizon
        in the forward sampling direction
    h - index to the value
    
    20150601 Scott Havens
    '''
    
    n = len(x)-1  # number of points to look at
    offset = 1      # offset from current point to start looking
    
    if i == n:
        h = i
    else:
        # origin point
        xo = x[i]
        zo = z[i]
        
        # array of other point, ignore the neighbor point
        xi = np.array(x[i+offset:])
        zi = np.array(z[i+offset:])
        
        # calculate the slope to each point
        s = (zi - zo)/(xi - float(xo))
        
        # find the maximum value
        if np.sum(s>0) == 0:
            h = i
        else:
            h = np.argmax(s) + i + offset
    
    # calculate the cosz
    cz = _cosz(x[i],z[i],x[h],z[h])
    
    
    return cz, h

def hor1f_simple(x,z):
    '''
    Calculate the horizon pixel for all x,z
    This mimics the simple algorthim from Dozier 1981
    to help understand how it's working
    
    Works backwards from the end but looks forwards for
    the horizon
    90% faster than rad.horizon
    
    Inputs:
    x - horizontal distances for points
    z - elevations for the points
    
    Output:
    h - index to the horizon point
    
    20150601 Scott Havens
    '''
    
    N = len(x)  # number of points to look at
#     offset = 1      # offset from current point to start looking
    
    # preallocate the h array
    h = np.zeros((N,1), dtype=int)
    h[N-1] = N-1
    i = N - 2
    
    # work backwarks from the end for the pixels
    while i >= 0:
        h[i] = i
        j = i + 1   # looking forward
        max_tan = 0
        
        while j < N:
            sij = _slope(i,z[i],j,z[j])
            
            if sij > max_tan:
                h[i] = j
                max_tan = sij
            
            j = j + 1
        i = i - 1
    
    return h

def hor1f(x, z, offset=2):
    '''
    BROKEN: Haven't quite figured this one out
    
    Calculate the horizon pixel for all x,z
    This mimics the algorthim from Dozier 1981 and the 
    hor1f.c from IPW
    
    Works backwards from the end but looks forwards for
    the horizon
    
    xrange stops one index before [stop]
    
    Inputs:
    x - horizontal distances for points
    z - elevations for the points
    
    Output:
    h - index to the horizon point
    
    20150601 Scott Havens
    '''
    
    N = len(x)  # number of points to look at
    
    # preallocate the h array
    h = np.zeros((N,1), dtype=int)
    h[N-1] = N-1    # the end point is it's own horizon
        
    # work backwarks from the end for the pixels
    for i in xrange(N-2,-1,-1) :
        
        zi = z[i]
        
        # Start with next-to-adjacent point in either forward or backward
        # direction, depending on which way loop is running. Note that we
        # don't consider the adjacent point; this seems to help reduce noise.
        k = i + offset
        if k >= N: k -= 1
          
        # loop until horizon is found
        # xrange will set the maximum number of iterations that can be
        # performed based on the length of the vector
        for t in xrange(k,N):
            j = k
            k = h[j][0]
            sij = _slope(i, zi, j, z[j])
            sihj = _slope(i, zi, k, z[k])
            
            # if slope(i,j) >= slope(i,h[j]), horizon has been found; otherwise
            # set j to k (=h[j]) and loop again
            # or if we are at the end of the section
            if sij < sihj: # or k == N-1:
                break
            
        
        # if slope(i,j) > slope(j,h[j]), j is i's horizon; else if slope(i,j)
        # is zero, i is its own horizon; otherwise slope(i,j) = slope(i,h[j])
        # so h[j] is i's horizon
        if sij > sihj:
            h[i] = j
        elif sij == 0:
            h[i] = i
        else:
            h[i] = k
        
    
    return h
 
def _slope(xi,zi,xj,zj):
    '''
    Slope between the two points only if greater than 0
    20150603 Scott Havens
    '''
    
    return 0 if zj <= zi else (zj - zi) / (float(xj) - xi)   

    
def _cosz(x1,z1,x2,z2):
    '''
    Cosize of the zenith between two points
    
    20150601 Scott Havens
    '''
    d = np.sqrt((x2 - x1)**2 + (z2 - z1)**2)
    
    if d == 0:
        v = 0
    else:
        v = (z2 - z1)/float(d)
    
    return v
    
    

def sunang(date, lat, lon, zone=0, slope=0, aspect=0):
    '''
    Wrapper for the IPW sunang function
    
    Inputs:
    date - date to calculate sun angle for (datetime object)
    lat - latitude in decimal degrees
    lon - longitude in decimal degrees
    zone - The  time  values  are  in the time zone which is min minutes 
        west of Greenwich (default: 0).  For example, if input times are 
        in Pacific Standard Time, then min would be 480.
    slope (default=0) - slope of surface
    aspect (default=0) - aspect of surface
    
    Output:
    cosz - cosine of the zeinith angle 
    azimuth - solar azimuth
    
    Created April 17, 2015
    Scott Havnes    
    '''
    
      
    # date string
    dstr = date.strftime('%Y,%m,%d,%H,%M,%S')
    
    # degree strings
    d, m, sd = deg_to_dms(lat)
    lat_str = str(d) + ',' + str(m) + ',' + '%02.1i' % sd
    
    d, m, sd = deg_to_dms(lon)
    lon_str = str(d) + ',' + str(m) + ',' + '%02.1i' % sd
    
    # prepare the command
    cmd_str = 'sunang -b %s -l %s -t %s -s %i -a %i -z %i' % \
        (lat_str, lon_str, dstr, slope, aspect, zone)

    p = sp.Popen(cmd_str,stdout=sp.PIPE, shell=True, env={"PATH": IPW})
    
    # get the results
    out, err = p.communicate()
    
    c = out.rstrip().split(' ')
    cosz = float(c[1])
    azimuth = c[3]
    
    return cosz, azimuth




def deg_to_dms(deg):
    '''
    Decimal degree to degree, minutes, seconds
    '''
    d = int(deg)
    md = abs(deg - d) * 60
    m = int(md)
    sd = (md - m) * 60
    return [d, m, sd]

def _linspace(start, stop, num):
    '''
    My own implementation of linspace for this work
    20150602 Scott Havens
    ''' 
    
    d = stop - start
    i = np.arange(num)
    
    return start + i*d/(num - 1)
    
    









