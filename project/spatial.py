from __future__ import division
import numpy as np
from datetime import datetime
from datetime import timedelta
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from shortest_element_path import shortest_element_path
from fvcomClass import FVCOM
import numexpr as ne
import scipy.io as sio

filename = '/EcoII/july_2012/output/dngrid_0001_03.nc'
siglay = np.array([0.98999,0.94999,0.86999,0.74999,0.58999,0.41000,0.25000,0.13000,0.05000,0.01000])

ax = [-66.35, -66.325, 44.261, 44.272]
#data = FVCOM(filename)
data = FVCOM(filename, ax=ax)
data.el_region()

u = data.u[:, :, data.region_e[:,0]]
v = data.v[:, :, data.region_e[:,0]]
#ww = data.ww[:, :, data.region_e[:,0]]

#vel = ne.evaluate('sqrt(u**2 + v**2 + ww**2)')
vel = ne.evaluate('sqrt(u**2 + v**2)')

mean_vel = np.mean(vel, axis=0)
