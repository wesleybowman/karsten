import numpy as np
from loadnc2d import loadnc2d
from ncdatasort import ncdatasort
from mjd2num import mjd2num
from closest_point import closest_point
from cf_u_rated_turbs import *
import matplotlib.pyplot as plt


filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'

(x, y, ua, va, trinodes,
 el, h, time, siglev, siglay,
 nbe, a1u, a2u, aw0, awx, awy,
 lon, lat, nele, node) = loadnc2d(filename)

#(nodexy,uvnode,dt,deltat,
# hour,thour,TP,rho,g,period,
# nodell,uvllnode)=ncdatasort(x,y,time*24*3600,trinodes,lon,lat)
(nodexy, uvnodexy, dt, deltat,
 hour, thour, TP, rho, g, period,
 nodell, uvnodell, trinodes) = ncdatasort(x,y,time*24*3600,trinodes,lon,lat)

time=mjd2num(time);

min_depth=30
max_depth=80

#power station location
#ps_loc=[-64.4031,45.3710];
#ps_loc(2,:)=[-64.4031,45.3710];
ps_loc = np.array([[-64.4031,45.3710],[-64.4031,45.3710]])
#psi=closest_point(ps_loc,[long,lat]);
psi=closest_point(ps_loc,lon,lat);
#xc=mean(x[trinodes],2);
#yc=mean(y[trinodes],2);
xc = np.mean(x[trinodes],axis=1)
yc = np.mean(y[trinodes],axis=1)
#distance=sqrt( (xc-x(psi(1))).^2+(yc-y(psi(1))).^2);
distance = np.sqrt( (xc - x[psi[0]])**2+(yc-y[psi[0]])**2)

x=x-x[psi[0]];
y=y-y[psi[1]];
plot_range=[-10000, 5000, -5000, 2000];
#xc=mean(x(trinodes),2);
#yc=mean(y(trinodes),2);

max_distance =5000;


speed = np.sqrt(ua*ua+va*va)
hc=np.mean(h[trinodes],axis=1);

#u = speed[:,0]
u_rated = np.arange(0, 8 + 0.25, 0.25)
turb = createTurbines(1)
turb = fillTurbines(turb, 10, 0.4, 1025, 1,2.5,5,0,0.4)


#turbines = cf_u_rated_turbs(speed, turb, u_rated)
turbines = np.load('../pythonTurbines')

meanP = turbines.meanP
capacity_factor=turbines.cf
#a = np.isnan(capacity_factor)
#capacity_factor[a]=0
capacity_factor[np.isnan(capacity_factor)]=0

tp = turbines.Prated

# Skipped analysis data since its not needed besides Rayleigh
Rayleigh = np.array([0.97,1])

#locate N turnbines
N=10000

#turbine spacing
turbine_diameter=10
spacing_along=15*turbine_diameter
spacing_across=2*turbine_diameter

turbine_power=meanP

#score=(55-capacity_factor)/10+abs(hc'-40)/10 +(max(distance',2500)-2500)/500;
#score=20*(1-meanP/1e6)+abs(hc'-40)/40 +((max(distance',3000)-3000)/250).^1;
score=20*(1-meanP/1e6)
turbine_score=score
#Find best location

loci = np.argmin(turbine_score)


