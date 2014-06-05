import numpy as np
from loadnc2d import loadnc2d
from ncdatasort import ncdatasort
from mjd2num import mjd2num
from closest_point import closest_point
#from cf_u_rated_turbs import *
#import matplotlib.pyplot as plt
#import scipy.io as sio
import cPickle as pickle
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv


filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'

# only variables used are: x,y, time, trinodes, lat, lon, h
(x, y, ua, va, trinodes,
 el, h, time, siglev, siglay,
 nbe, a1u, a2u, aw0, awx, awy,
 lon, lat, nele, node) = loadnc2d(filename)

(nodexy, uvnodexy, dt, deltat,
 hour, thour, TP, rho, g, period,
 nodell, uvnodell, trinodes) = ncdatasort(x, y, time*24*3600,
                                          trinodes, lon, lat)

time = mjd2num(time)

#power station location
ps_loc = np.array([[-64.4031, 45.3710], [-64.4031, 45.3710]])
psi = closest_point(ps_loc, lon, lat)

xc = np.mean(x[trinodes], axis=1)
yc = np.mean(y[trinodes], axis=1)

distance = np.sqrt((xc - x[psi[0]])**2 + (yc - y[psi[0]])**2)

x = x - x[psi[0]]
y = y - y[psi[1]]
plot_range = [-10000, 5000, -5000, 2000]

max_distance = 5000

speed = np.sqrt(ua*ua+va*va)
hc = np.mean(h[trinodes], axis=1)

#u = speed[:,0]
u_rated = np.arange(0, 8 + 0.25, 0.25)

# Skipped analysis data since its not needed besides Rayleigh
Rayleigh = np.array([0.97, 1])

ii=100000

coef = ut_solv(time, ua[:,ii], va[:,ii], uvnodell[ii,1],
                cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                nodiagn=True, linci=True, conf_int=True)

pickle.dump(coef, open( "speedCoef.p", "wb"))

coef = ut_solv(time, ua[:,ii], [], uvnodell[ii,1],
                cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                nodiagn=True, linci=True, conf_int=True)

pickle.dump(coef, open( "elevCoef.p", "wb"))
