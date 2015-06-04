'''
20150602 Scott Havens

Calcualte the horizon angle for full DEM for each point on a DEM give the sunangle

'''


import numpy as np
import matplotlib.pyplot as plt
from isnobal import rad, ipw
import subprocess as sp
import os
from scipy import ndimage, misc
from datetime import datetime

IPW = "/Users/ipw/ipw-2.1.0/bin:" # IPW executables
PATH = os.environ.copy()["PATH"] + ':/usr/local/bin:' + IPW
TMPDIR = os.environ['TMPDIR'].rstrip('/')

# Load the DEM
dem = '/Volumes/Drobo1/BRB/BRB-wy09/spatial_WRF/data/topo/dem100m.ipw'
i = ipw.IPW(dem)

# some sun angle, Oct 1, 12:00
mu = 0.670150
azm = 10.373    # angle from south
# azm = 45
# azm = azm*np.pi/180


# process the dem file
d = i.bands[0]
Y = d.bline + np.arange(d.nlines)*d.dline
X = d.bsamp + np.arange(d.nsamps)*d.dsamp
Z = d.data
Z = Z[:,1:200]
# Z,f = np.meshgrid(np.arange(5),np.arange(5))
zm,zn = Z.shape    
 

# run the horizon
startTime = datetime.now()
cz = rad.ihorizon(X, Y, Z, azm, mu, 1)
print datetime.now() - startTime

# run the C code
startTime = datetime.now()
cmd = 'horizon -u %s -a %s %s > junk' % (str(mu), str(azm), dem) 
# sp.Popen(cmd, shell=True, env={"PATH": PATH, "IPW": IPW, "TMPDIR": TMPDIR}).wait()
# # horizon -u $mu -a $azm $efile > $hfile
print datetime.now() - startTime

# impor the test file
c = ipw.IPW('/Users/scott/Documents/Projects/isnobal/python/rad/junk')

#-- Plot...
# plt.Figure(figsize=(8,6))
fig, axes = plt.subplots(1,2)
# fig.set_size_inches(10,10)
im = axes[0].imshow(cz)#,extent=[X[0], X[-1:], Y[-1:], Y[0]])
plt.colorbar(im, ax=axes[0])

# im1 = axes[1].imshow(Z-Zs)
im1 = axes[1].imshow(c.bands[0].data)
# im1.set_clim(0,100)
plt.colorbar(im1, ax=axes[1])
 
plt.show()







