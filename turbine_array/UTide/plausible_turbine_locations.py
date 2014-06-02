import numpy as np
from loadnc2d import loadnc2d
from ncdatasort import ncdatasort
from mjd2num import mjd2num
from closest_point import closest_point
from cf_u_rated_turbs import *
import matplotlib.pyplot as plt
import scipy.io as sio
import cPickle as pickle
import sys
sys.path.append('/home/wesley/github/UTide/')
#from UTide.ut_solv.master import *
#from UTide.ut_solv import ut_solv
#from UTide import ut_solv
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

min_depth = 30
max_depth = 80

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
turb = createTurbines(1)
turb = fillTurbines(turb, 10, 0.4, 1025, 1, 2.5, 5, 0, 0.4)


#turbines = cf_u_rated_turbs(speed, turb, u_rated)
turbines = np.load('turbines.dat')

meanP = turbines.meanP
capacity_factor = turbines.cf
capacity_factor[np.isnan(capacity_factor)] = 0

tp = turbines.Prated

# Skipped analysis data since its not needed besides Rayleigh
Rayleigh = np.array([0.97, 1])

#locate N turnbines
N = 10000

#turbine spacing
turbine_diameter = 10
spacing_along = 15 * turbine_diameter
spacing_across = 2 * turbine_diameter

turbine_power = meanP

#score=(55-capacity_factor)/10+abs(hc'-40)/10 +(max(distance',2500)-2500)/500;
#score=20*(1-meanP/1e6)+abs(hc'-40)/40 +((max(distance',3000)-3000)/250).^1;
score = 20 * (1 - meanP / 1e6)
turbine_score = score
#Find best location


loci = np.empty((N))
for ii in xrange(N):

    loci[ii] = np.argmin(turbine_score)

    # do u_tide analysis at loc
    coef = ut_solv(time, ua[:,loci[ii]], va[:,loci[ii]], uvnodell[loci[ii],1],
                  cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                  nodiagn=True, linci=True, conf_int=False)

#    coef = ut_solv(time, ua[:,loci[ii]], np.array([]), uvnodell[loci[ii],1],
#                  'auto', Rayleigh[0],'NoTrend','Rmin', 'OLS',
#                  'NoDiagn', 'LinCI')

    # for testing
    if ii == 0:
        pickle.dump(coef, open( "coef.p", "wb"))
        #coef.dump('coef.dat')

    cx = np.cos(coef['theta'][0]*np.pi/180)
    cy = np.sin(coef['theta'][0]*np.pi/180)

    #find new xy as distance from location, and in direction of M2 ellipse
    newx = (xc-xc[loci[ii]]) * cx + (yc - yc[loci[ii]]) * cy
    newy = -(xc-xc[loci[ii]]) * cy + (yc - yc[loci[ii]]) * cx

    #a=find(abs(newx)<spacing_along & abs(newy)<spacing_across);
    #a = np.argwhere(np.abs(newx)< spacing_along & np.abs(newy) < spacing_across)
    #a = np.where(np.abs(newx)< spacing_along & np.abs(newy) < spacing_across)
    a = np.argwhere(np.logical_and((np.abs(newx) < spacing_along),
                                   (np.abs(newy) < spacing_across)))[0]

    # turbine_power(a)=0;
    turbine_score[a] = 100

pickle.dump(loci, open('loci.p', 'wb'))


#loci = np.argmin(turbine_score)


#    coef= ut_solv(time, ua[:,loci], va[:,loci], uvnodell[loci,1],
#                  'auto', Rayleigh[0],'NoTrend','Rmin', 'OLS',
#                  'NoDiagn', 'LinCI')
