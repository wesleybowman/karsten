from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as net
from datetime import datetime
from datetime import timedelta
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from shortest_element_path import shortest_element_path
from fvcomClass import FVCOM

def time_index(a):
    a = a.reindex(index=a.index.to_datetime())
    return a

def magnitude(a):
    mag = np.sqrt(a['u']**2+a['v']**2+a['w']**2)
    return mag

def theta(a):
    a = np.arctan2(a['v'],a['u'])
    return a



filename = '/array2/data3/rkarsten/july_2012/output/dngrid_0001_02.nc'
siglay = np.array([0.98999,0.94999,0.86999,0.74999,0.58999,0.41000,0.25000,0.13000,0.05000,0.01000])


data = FVCOM(filename)
# North-South
ind = data.closest_point([-66.3385, -66.3385], [44.277, 44.277])
# East- West
ind = data.closest_point([-66.3412, -66.3324], [44.277, 44.277])

short_path = shortest_element_path(data.latc, data.lonc,
                                    data.lat,data.lon,
                                    data.nv, data.h)

#short_path = shortest_element_path(filename)
el, _ = short_path.getTargets([ind])

t_slice = ['2014-02-02T06:45:00','2014-02-02T07:05:00']
t_slice = np.array(t_slice,dtype='datetime64[us]')

t = data.time.shape[0]
l = []
for i in range(t):
    date = datetime.fromordinal(int(data.time[i]))+timedelta(days=data.time[i]%1)-timedelta(days=366)
    l.append(date)

time = np.array(l,dtype='datetime64[us]')
if t_slice.shape[0] != 1:
    argtime = np.argwhere((time>=t_slice[0])&(time<=t_slice[-1])).flatten()

#vel = np.sqrt(nc['u'][argtime,:,el]**2+nc['v'][argtime,:,el]**2+nc['ww'][argtime,:,el]**2)
#vel = np.sqrt(data.u[argtime,:,el]**2+data.u[argtime,:,el]**2+data.ww[argtime,:,el]**2)
vel = np.sqrt(data.u[:, :, el[0]]**2 + data.v[:, :, el[0]]**2 + data.ww[:, :, el[0]]**2)

lat = data.latc[el]
lon = data.lonc[el]
#lat = nc['latc'][el]
#lon = nc['lonc'][el]

line = lon
print vel.shape
vmax = 2.5
vmin = 0

for i in range(vel.shape[0]):
    fig,ax = plt.subplots()
    plt.rc('font',size='22')
    levels = np.linspace(0,3.3,34)
    cs = ax.contourf(line,siglay,vel[i,:],levels=levels)
    ax.contour(line,siglay,vel[i,:],cs.levels,colors='k',hold='on')
    cbar = fig.colorbar(cs,ax=ax)
    cbar.set_label(r'Velocity $(m/s)$', rotation=-90,labelpad=30)
    #plt.title(str(time[i]))
    ax.set_xlabel('Longitude')
    scale = 1
    ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
    ax.xaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_formatter(ticks)
    saveName = '/figures/figure{0}.png'.format(i)
    plt.savefig(saveName, bbox_inches=0)
