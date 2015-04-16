'''
Test reading an ipw file
'''


import numpy as np
from isnobal import ipw
import netCDF4 as nc
import matplotlib.pyplot as plt
import subprocess as sp


fileName = '../test_1band.ipw'


o = ipw.IPW(fileName)
# print(ipw)

o.write('out.ipw',16)


# load both and compare
# o = ipw.IPW(fileName)
n = ipw.IPW('out.ipw')
 
d = o.bands[0].data - n.bands[0].data
 
print(sum(sum(d)))
 
# print(i.header_lines)
# print(i.binary_data[:100])


plt.subplot(131)
plt.imshow(o.bands[0].data)
plt.colorbar()

plt.subplot(132)
plt.imshow(n.bands[0].data)
plt.colorbar()
 
plt.subplot(133)
plt.imshow(d)
plt.colorbar()

plt.show()

