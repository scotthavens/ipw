'''
Test reading an ipw file
'''


import numpy as np
from ipw import ipw as i
import netCDF4 as nc
import matplotlib.pyplot as plt
import subprocess as sp


fileName = '/Users/scott/Documents/Projects/isnobal/python/test_1band.ipw'


ipw = i.IPW(fileName)
# print(ipw)

ipw.write('out.ipw',16)


# load both and compare
o = i.IPW(fileName)
n = i.IPW('out.ipw')

d = o.bands[0].data - n.bands[0].data

print(sum(sum(d)))

# print(ipw.header_lines)
# print(ipw.binary_data[:100])


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

