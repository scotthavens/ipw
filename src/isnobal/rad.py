'''
Created on Apr 17, 2015

@author: scott
'''

__version__ = '0.0.1'


import numpy as np
# import datetime as dt
import subprocess as sp
import math
# import progressbar
# from joblib import Parallel, delayed  
# import multiprocessing as mp
# import ctypes
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
STEF_BOLTZ = 5.6697e-8  # stephman boltzman constant
EMISS_TERRAIN = 0.98    # emissivity of the terrain
EMISS_VEG = 0.96    # emissivity of the vegitation

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
    cosz - cosine local solar illumination angle matrix
    gsize - gsize is effective grain radius of snow after last storm (mu m)
    maxgsz -  maxgsz is maximum grain radius expected from grain growth (mu m)
    dirt - dirt is effective contamination for adjustment to visible albedo (usually between 1.5-3.0)
    
    Output
    
    
    Created April 17, 2015
    Modified July 23, 2015 - take image of cosz and calculate albedo for one time step
    Scott Havens
    '''
    
#     telapsed = np.array(telapsed)
    
    # check inputs
    if gsize <= 0 or gsize > 500:
        raise Exception("unrealistic input: gsize=%i", gsize)

    if (maxgsz <= gsize or maxgsz > 2000):
        raise Exception("unrealistic input: maxgsz=%i", maxgsz)
    if 1 >= dirt >= 10:
        raise Exception("unrealistic input: dirt=%i", dirt)
#     if dirt <= 1:
#         raise Exception("unrealistic input: dirt=%i", dirt)
    
    # set initial grain radii for vis and ir
    radius_ir = math.sqrt(gsize);
    range_ir = math.sqrt(maxgsz) - radius_ir;
    radius_v = math.sqrt(dirt * gsize);
    range_v = math.sqrt(dirt * maxgsz) -radius_v;
    
    # calc grain growth decay factor
    growth_factor = growth(telapsed + 1.0);

    # calc effective gsizes for vis & ir
    gv = radius_v + (range_v * growth_factor);
    gir = radius_ir + (range_ir * growth_factor);

    #calc albedos for cos(z)=1
    alb_v_1 = MAXV - (gv / VFAC);
    alb_ir_1 = MAXIR * np.exp(IRFAC * gir);

    # calculate effect of cos(z)<1

    # adjust diurnal increase range 
    dzv = gv * VZRG;
    dzir = (gir * IRZRG) + IRZ0;

    # calculate albedo
    alb_v = alb_v_1;
    alb_ir = alb_ir_1;
    
    # correct if the sun is up  
    ind = cosz > 0.0;   
    alb_v[ind] += dzv[ind] * (1.0 - cosz[ind]);
    alb_ir[ind] += dzir[ind] * (1.0 - cosz[ind]);
    
    return alb_v, alb_ir

def ihorizon(x, y, Z, azm, mu=0, offset=2, ncores=0):
    '''
    Calculate the horizon values for an entire DEM image
    for the desired azimuth
    
    Assumes that the step size is constant
    
    Inputs:
        X - vector of x-coordinates
        Y - vector of y-coordinates
        Z - matrix of elevation data
        azm - azimuth to calculate the horizon at
        mu - 0 -> calculate cos(z)
             - >0 -> calculate a mask whether or not the point can see the sun
    
    Outputs:
        H   - if mask=0 cosine of the local horizonal angles
            - if mask=1 index along line to the point
    
    20150602 Scott Havens 
    '''
    
    # check inputs
    azm = azm*np.pi/180 # degress to radians
    m,n = Z.shape
     
    # transform the x,y into the azm direction xr,yr
    xi, yi = np.arange(-n/2,n/2), np.arange(-m/2,m/2)
    X, Y = np.meshgrid(xi, yi)
    xr = X*np.cos(azm) - Y*np.sin(azm)
    yr = X*np.sin(azm) + Y*np.cos(azm)
    
    # xr is the "new" column index for the profiles
    # yr is the distance along the profile
    xr = xr.round().astype(int)
    yr = (x[2] - x[1]) * yr
    
    H = np.zeros(Z.shape)
#     pbar = progressbar.ProgressBar(n).start()
#     j = 0
    
    # loop through the columns
#     if ncores == 0:
    for i in xrange(-n/2,n/2):
        find_horizon(i, H, xr, yr, Z, mu)
             
#     else:
#     shared_array_base = mp.Array(ctypes.c_double, m*n)
#     sH = np.ctypeslib.as_array(shared_array_base.get_obj())
#     sH = sH.reshape(m, n)
#     sxr = np.ctypeslib.as_array(shared_array_base.get_obj())
#     sxr = sxr.reshape(m, n)
#     syr = np.ctypeslib.as_array(shared_array_base.get_obj())
#     syr = syr.reshape(m, n)
#     sZ = np.ctypeslib.as_array(shared_array_base.get_obj())
#     sZ = sZ.reshape(m, n)
#     def wrap_horizon(i, def_param1=sH, def_parm2=sxr, def_param3=syr, def_param4=sZ):
#         find_horizon(i, sH, sxr, syr, sZ, 0.67)
 
#     pool = mp.Pool(processes=4)
#     [pool.apply(find_horizon, args=(i, shared_array, xr, yr, Z, mu, offset)) for i in range(-n/2,n/2)]
#     pool.map(wrap_horizon, range(-n/2,n/2))
    
#     print(shared_array)
                
    return H   



def find_horizon(i, H, xr, yr, Z, mu):
    # index to profile and get the elevations
    ind = xr == i
    zi = Z[ind]
    
    # distance along the profile
    di = yr[ind]
    
    # sort the y values and get the cooresponding elevation
    idx = np.argsort(di)
    di = di[idx]
    zi = zi[idx]
    
    # if there are some values in the vector
    # calculate the horizons
    if len(zi) > 0:
#             h2 = hor1f_simple(di, zi)
        h = hor1f(di, zi)
        
        cz = _cosz(di, zi, di[h], zi[h])
        
        # if we are making a mask
        if mu > 0:
#                 iz = cz == 0    # points that are their own horizon
            idx = cz > mu   # points sheltered from the sun
            cz[idx] = 0
            cz[~idx] = 1
#                 cz[iz] = 1
            
    H[ind] = cz

#         j += 1
#         pbar.update(j)
#     
#     pbar.finish()
        
    
        

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
    h = np.zeros(N, dtype=int)
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

def hor1f(x, z, offset=1):
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
    x = np.array(x)
    z = np.array(z)
    
    # preallocate the h array
    h = np.zeros(N, dtype=int)
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
            k = h[j]
                        
            sij = _slope(x[i], zi, x[j], z[j])
            sihj = _slope(x[i], zi, x[k], z[k])
                        
            # if slope(i,j) >= slope(i,h[j]), horizon has been found; otherwise
            # set j to k (=h[j]) and loop again
            # or if we are at the end of the section
            if sij > sihj: # or k == N-1:
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
    Slope between the two points only if the pixel is higher
    than the other
    20150603 Scott Havens
    '''
    
    return 0 if zj <= zi else (zj - zi) / (xj - float(xi))
    
    
def _cosz(x1,z1,x2,z2):
    '''
    Cosize of the zenith between two points
    
    20150601 Scott Havens
    '''
    d = np.sqrt((x2 - x1)**2 + (z2 - z1)**2)
    diff = z2 - z1
    
#     v = np.where(diff != 0., d/diff, 100)
    
    i = d == 0
    d[i] = 1
    v = diff/d
    v[i] = 0
    
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


def cf_cloud(beam, diffuse, cf):
    '''
    Correct beam and diffuse irradiance for cloud attenuation at a single
    time, using input clear-ky global and diffuse radiation calculations supplied by
    locally modified toporad or locally modified stoporad
    
    Inputs:
    beam - global irradiance
    diffuse - diffuse irradiance
    cf- cloud attenuation factor
    
    Outputs:
    c_grad - cloud corrected gobal irradiance
    c_drad - cloud corrected diffuse irradiance
    
    20150610 Scott Havens - adapted from cloudcalc.c
    '''
    
    # define some constants
    CRAT1 = 0.15
    CRAT2 = 0.99
    CCOEF = 1.38
    
        
    # cloud attenuation, beam ratio is reduced
    bf_c = CCOEF * (cf - CRAT1)**2
    c_grad = beam * cf
    c_brad = c_grad * bf_c
    c_drad = c_grad - c_brad
    
    # extensive cloud attenuation, no beam
    ind = cf <= CRAT1
    c_brad[ind] = 0
    c_drad[ind] = c_grad[ind]
    
    # minimal cloud attenution, no beam ratio reduction
    ind = cf > CRAT2
    c_drad[ind] = diffuse[ind] * cf[ind]
    c_brad[ind] = c_grad[ind] - c_drad[ind]
    
    return c_grad, c_drad
    

def veg_beam(data, height, cosz, k):
    '''
    Apply the vegetation correction to the beam irradiance
    using the equation from Links and Marks 1999
    
    S_b,f = S_b,o * exp[ -k h sec(theta) ] or
    S_b,f = S_b,o * exp[ -k h / cosz ]
    
    20150610 Scott Havens
    '''
    
    # ensure that the sun is visible
    cosz[cosz <= 0] = 0.01
    
    return data * np.exp(-k * height / cosz)

    
def veg_diffuse(data, tau):
    '''
    Apply the vegetation correction to the diffuse irradiance
    using the equation from Links and Marks 1999
    
    S_d,f = tau * S_d,o
    
    20150610 Scott Havens
    '''
    
    return tau * data
    

def set_min_max(data, min_val, max_val):
    '''
    Ensure that the data is in the bounds of min and max
    20150611 Scott Havens
    '''
    
    data[data <= min_val] = min_val
    data[data >= max_val] = max_val
    
    return data
    
    
def thermal_correct_terrain(th, ta, viewf):
    '''
    Correct the thermal radiation for terrain assuming that
    the terrain is at the air temperature and the pixel and 
    a sky view
    
    Inputs:
    th - thermal radiation
    ta - air temperature [C]
    viewf - sky view factor from view_f
    
    Outputs:
    th_c - correct thermal radiation
    
    20150611 Scott Havens
    '''
    
    # thermal emitted from the terrain
    terrain = STEF_BOLTZ * EMISS_TERRAIN * np.power(ta + 273.15, 4)
    
    # correct the incoming thermal
    return viewf * th + (1 - viewf) * terrain
    
    
def thermal_correct_canopy(th, ta, tau, veg_height, height_thresh=2):
    '''
    Correct thermal radiation for vegitation.  It will only correct
    for pixels where the veg height is above a threshold. This ensures
    that the open areas don't get this applied.  Vegitation temp
    is assumed to be at air temperature
    
    Inputs:
    th - thermal radiation
    ta - air temperature [C]
    tau - transmissivity of the canopy
    veg_height - vegitation height for each pixel
    height_thresh - threshold hold for height to say that there is veg in the pixel
    
    Output:
    th_c - corrected thermal radiation
    
    20150611 Scott Havens
    '''
    
    # thermal emitted from the canopy
    veg = STEF_BOLTZ * EMISS_VEG * np.power(ta + 273.15, 4)
    
    # pixels with canopy above the threshold
    ind = veg_height > height_thresh
    
    # correct incoming thermal
    th[ind] = tau[ind] * th[ind] + (1 - tau[ind]) * veg[ind]
    
    return th
    
    
    
    
    
    
    
    
    







