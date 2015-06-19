'''
20150601 Scott Havens

Calcualte the horizon angle for each point on a DEM give the sunangle

'''


import numpy as np
import matplotlib.pyplot as plt
from isnobal import rad, ipw

from datetime import datetime

# Load the DEM
dem = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/topo/dem100m.ipw'
i = ipw.IPW(dem)

# some sun angle, Oct 1, 12:00
mu = 0.670150
azm = 10.373    # angle from south
azm = 45


# process the dem file
# d = i.bands[0]
# Y = d.bline + np.arange(d.nlines)*d.dline
# # Y = Y[::-1]
# X = d.bsamp + np.arange(d.nsamps)*d.dsamp
# 
# # extract the line
# x0, y0 = 700, 390 # south fort payette
# x0, y0 = 150, 1000   # boise area
# x0, y0 = 1300, 600
# x0, y0 = 250, 1200
# x0, y0 = 1999, 10
# x1, y1 = x0 + 20*np.sin(azm*np.pi/180), y0 + 20*np.cos(azm*np.pi/180)
# 
# l = int(np.hypot(x1-x0, y1-y0))
# xi, yi = np.linspace(x0, x1, l), np.linspace(y0, y1, l) 
# xi = xi.astype(np.int)  # indicies to the data
# yi = yi.astype(np.int)
# 
# m,n = d.data.shape
# ind = (xi < n) & (yi < m) # ensure that the values are on the grid
# xi = xi[ind]
# yi = yi[ind]
# 
# # Extract the values along the line
# xl = X[xi]
# yl = Y[yi]
# di = np.hypot(xl[1:]-xl[:-1],yl[1:]-yl[:-1])    # distances between points
# di = np.cumsum(di)
# di = np.insert(di,0,0)   # add zero to the front for the current point
# 
# # remove any duplicated values
# di, ind = np.unique(di, return_index=True)
# xi = xi[ind]
# yi = yi[ind]
# 
# zi = d.data[yi, xi]
# 
# 
# 
# cz,hi = rad.horizon(0,di,zi)

# follow Dozier 1981 data for testing

x = [0, 100, 150, 200, 225, 375, 400, 450, 500, 680, 780, 850, 875, 900, 960, 1010]
y = [80, 180, 130, 50, 170, 90, 200, 95, 50, 60, 140, 160, 60, 80, 50, 75]
 
startTime = datetime.now()
h = np.empty(len(x),)
cz = np.empty(len(x),)
for j in xrange(1,1000):
    for i in range(0,len(x)):
        h[i] = rad.horizon(i, x, y)
# calculate the cosz
# cz = rad._cosz(x,y,x[h],y[h])
print datetime.now() - startTime

# hor1f_simple is 90% faster than horizon
startTime = datetime.now()
for j in xrange(1,1000):
    hi = rad.hor1f_simple(x, y)
print datetime.now() - startTime

# hor1f is 70% faster than horizon for small arrays
# It is much faster than hor1f_simple on larger arrays
startTime = datetime.now()
for j in xrange(1,1000):
    hi = rad.hor1f(x, y, 1)
print datetime.now() - startTime


print(hi)
plt.Figure
plt.plot(x,y,'ro-')
plt.show() 





# #-- Plot...
# # plt.Figure(figsize=(8,6))
# fig, axes = plt.subplots(2,1)
# # fig.set_size_inches(10,10)
# im = axes[0].imshow(d.data)#,extent=[X[0], X[-1:], Y[-1:], Y[0]])
# axes[0].plot([x0,x1],[y0,y1], 'ro-')
# axes[0].axis('image')
#  
# axes[1].plot(di,zi)
# axes[1].plot(di[hi],zi[hi],'ro')
#  
# plt.show()







