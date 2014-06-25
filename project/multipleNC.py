import fnmatch
import os
import netCDF4 as nc
import numpy as np

filename = '/array/data1/073208o/workspace_matlab/runs/2013_run'

matches = []
for root, dirnames, filenames in os.walk(filename):
  for filename in fnmatch.filter(filenames, 'dngrid_0001.nc'):
      matches.append(os.path.join(root, filename))

time = np.array([])
ua = np.array([])
for i,v in enumerate(matches):
    data = nc.Dataset(v, 'r')
    t = data.variables['time'][:]
    u = data.variables['ua'][:, 0]

    time = np.hstack((time, t))
    ua = np.hstack((ua, u))
