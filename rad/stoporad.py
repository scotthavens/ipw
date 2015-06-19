'''
Test a python version of the script stoporad
Goals are to improve the ialbedo function to use the distributed
days since last storm to calculate albedo

20150603 Scott Havens
'''


import subprocess as sp
import os
from isnobal import netcdf
import progressbar
from datetime import datetime
startTime = datetime.now()

IPW = "/Users/ipw/ipw-2.1.0/bin:" # IPW executables
PATH = os.environ.copy()["PATH"] + ':/usr/local/bin:' + IPW
TMPDIR = os.environ['TMPDIR'].rstrip('/')

# Atmospheric parameters
opt_depth = 100     # elevation that the optical thickness was measured
tau = 0.2           # optical thickness of atmosphere
omega = 0.85        
gamma = 0.3

# albedo calculation parameters
grain_size = 300        # effective grain radius of snow after last storm (mu m)
max_grain = 2000        # maximum grain radius expected from grain growth (mu m)
dirt = 2.0              # effective contamination for adjustment to visible albedo (usually between 1 and 5 (integer values only))

# study site parameters
tz_min_west = 420   # minutes west of GMT
wyear = 2009        # water year
wy = datetime(wyear-1,10,1)

# stoporad DEM parameters
stoporad_in = '../../topo/stoporad_in.ipw'

# sunang file
sang = '../sun_angle'

# days since last storm
netcdfFile = '../../precip/storm_days.nc'
n = netcdf.netcdf(netcdfFile, 'r')
storm = n.netcdf.variables['min_days']

#===============================================================================
# Run stoporad
#===============================================================================

# read the sunangle file
lines = [line.rstrip('\n') for line in open(sang)]
    
# loop over each and run stoporad
pbar = progressbar.ProgressBar(len(lines)).start()
j = 0
for l in lines:
    
    # split the line and remove any whitespaces
    ln = l.split(',')
    dt = ln[0].strip()
    mu = ln[1].strip()
    azm = ln[2].strip()
    
    # calculate storm day and current day (DOY since WY start format)
    c = datetime.strptime(dt,'%Y/%m/%d+%H:%M')
    wday = c - wy
    w = wday.total_seconds()/86400  # elapsed time time decimal days
    
    s = storm[j]/24 # storm days in decimal days
    
    d,h = divmod(wday.total_seconds()/3600, 24)     # days and hours
        
    # Calculate ir radiation
    alb_f = '../albedo/ir/%03i_%02i' % (d, h)
    out = './ir/%03i_%02i' % (d, h)
    
    ir_cmd = '/bin/sh stoporad -z %i -t %s -w %s -g %s -x 0.7,2.8 -s %s'\
        ' -d %s -f %i -y %i -A %s,%s -a %i -m %i -c %i -k %s -D %s > %s' \
        % (opt_depth, str(tau), str(omega), str(gamma), str(s), str(w), tz_min_west, wyear, \
           str(mu), str(azm), grain_size, max_grain, dirt, alb_f, stoporad_in, out) 
    proc = sp.Popen(ir_cmd, shell=True, env={"PATH": PATH, "IPW": IPW, "TMPDIR": TMPDIR}).wait()
    
    # Calculate visible radiation
    alb_f = '../albedo/vis/%03i_%02i' % (d, h)
    out = './vis/%03i_%02i' % (d, h)
    
    ir_cmd = '/bin/sh stoporad -z %i -t %s -w %s -g %s -x 0.28,0.7 -s %s'\
        ' -d %s -f %i -y %i -A %s,%s -a %i -m %i -c %i -k %s -D %s > %s' \
        % (opt_depth, str(tau), str(omega), str(gamma), str(s), str(w), tz_min_west, wyear, \
           str(mu), str(azm), grain_size, max_grain, dirt, alb_f, stoporad_in, out) 
    proc = sp.Popen(ir_cmd, shell=True, env={"PATH": PATH, "IPW": IPW, "TMPDIR": TMPDIR}).wait()
    
    j += 1
    pbar.update(j)
    


pbar.finish()
print datetime.now() - startTime

