from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv
#from coefNC import coefNC2D


def ncdatasort(x, y, time, trinodes, lon=None, lat=None):

    # Note, when time is passes, it is passed as
    # time*24*3600.
    # Does this need to be done this way?

    hour = 3600
    g = 9.806
    TP = 12.42
    rho = 1026
    period = (TP*3600)/(2*np.pi)

    # time=double(time)
    time = time+678942

    dt = time[1] - time[0]
    thour = time/hour
    deltat = thour[1]-thour[0]

    size = x.shape[0]
    nodexy = np.zeros((size, 2))
    nodexy[:, 0] = x
    nodexy[:, 1] = y
    nodexy = np.array((x, y)).T

    trinodes = trinodes.T-1
    triSize = trinodes.shape[0]
    uvnodexy = np.zeros((triSize, 2))

    uvnodexy[:, 0] = (nodexy[trinodes[:, 0], 0] + nodexy[trinodes[:, 1], 0] +
                      nodexy[trinodes[:, 2], 0]) / 3

    uvnodexy[:, 1] = (nodexy[trinodes[:, 0], 1] + nodexy[trinodes[:, 1], 1] +
                      nodexy[trinodes[:, 2], 1]) / 3

    if lat != None and lon != None:
        nodell = np.array((lon, lat)).T

        uvnodell = np.zeros((triSize, 2))

        uvnodell[:, 0] = (nodell[trinodes[:, 0], 0] +
                          nodell[trinodes[:, 1], 0] +
                          nodell[trinodes[:, 2], 0]) / 3

        uvnodell[:, 1] = (nodell[trinodes[:, 0], 1] +
                          nodell[trinodes[:, 1], 1] +
                          nodell[trinodes[:, 2], 1]) / 3

    else:
        'No nodell, uvnodell set to uvnodexy'
        uvnodell = uvnodexy

    return (nodexy, uvnodexy, dt, deltat, hour, thour,
            TP, rho, g, period, nodell, uvnodell, trinodes)


def mjd2num(x):

    y = x + 678942

    return y



def closest_point(points, lon, lat):

    point_list = np.array([lon,lat]).T

    closest_dist = ((point_list[:, 0] - points[:, 0, None])**2 +
                    (point_list[:, 1] - points[:, 1, None])**2)

    closest_point_indexes = np.argmin(closest_dist, axis=1)

    return closest_point_indexes


def datetime2matlabdn(dt):
    # ordinal = dt.toordinal()
    mdn = dt + timedelta(days=366)
    frac = (dt-datetime(dt.year, dt.month, dt.day, 0, 0, 0)).seconds / \
        (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac




filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
#filename = '/array/data2/jculina/2014/GIS_project/lunar_month_sep142011/output/smallcape_force_0001.nc'

data = nc.Dataset(filename, 'r')
x = data.variables['x'][:]
y = data.variables['y'][:]
xc = data.variables['xc'][:]
yc = data.variables['yc'][:]
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
lonc = data.variables['lonc'][:]
latc = data.variables['latc'][:]
h = data.variables['h'][:]
ua = data.variables['ua']
va = data.variables['va']
time = data.variables['time'][:]
trinodes = data.variables['nv'][:]
elev = data.variables['zeta']

(nodexy, uvnodexy, dt, deltat,
 hour, thour, TP, rho, g, period,
 nodell, uvnodell, trinodes) = ncdatasort(x, y, time*24*3600,
                                          trinodes, lon, lat)

time = mjd2num(time)

Rayleigh = np.array([1.0])
lonlat = np.array([lon, lat]).T
lonclatc = np.array([lonc, latc]).T

dim1 = len(lonc)

coef = ut_solv(time, ua[:, 0], va[:, 0], uvnodell[0, 1],
                cnstit='auto', rmin=1, notrend=True, method='ols',
                nodiagn=True, linci=True, conf_int=True)

opt = coef['aux']['opt']
del coef['aux']['opt']
aux = pd.DataFrame(coef['aux'])
del coef['aux']

dim2 = len(coef['Lsmaj'])

temp = np.zeros((dim1,dim2))

Lsmaj = pd.DataFrame(temp)
Lsmaj_ci = pd.DataFrame(temp)
Lsmin = pd.DataFrame(temp)
Lsmin_ci = pd.DataFrame(temp)
theta = pd.DataFrame(temp)
theta_ci = pd.DataFrame(temp)
g = pd.DataFrame(temp)
g_ci = pd.DataFrame(temp)
name = pd.DataFrame(temp)

A = pd.DataFrame(temp)
A_ci = pd.DataFrame(temp)
gA = pd.DataFrame(temp)
gA_ci = pd.DataFrame(temp)
nameA = pd.DataFrame(temp)


for i, (lonc,latc) in enumerate(lonclatc):
    print i

    coef = ut_solv(time, ua[:, i], va[:, i], uvnodell[i, 1],
                    cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                    nodiagn=True, linci=True, conf_int=True)

    opt = pd.DataFrame(coef['aux']['opt'].items())
    del coef['aux']['opt']
    aux = pd.DataFrame(coef['aux'])
    del coef['aux']
    c = pd.DataFrame(coef)

    cat = pd.concat([c,aux], axis=1)

    Lsmaj.iloc[i] = coef['Lsmaj']
    Lsmaj_ci.iloc[i] = coef['Lsmaj_ci']
    Lsmin.iloc[i] = coef['Lsmin']
    Lsmin_ci.iloc[i] = coef['Lsmin_ci']
    theta.iloc[i] = coef['theta']
    theta_ci.iloc[i] = coef['theta_ci']
    g.iloc[i] = coef['g']
    g_ci.iloc[i] = coef['g_ci']
    name.iloc[i] = coef['name']

#    Lsmaj[i,:] = cat['Lsmaj'].values
#    Lsmin[i,:] = cat['Lsmin'].values
#    g[i,:] = cat['g'].values
#    theta[i, :] = cat['theta'].values
#    name[i, :] = cat['name'].values
#    Lsmaj_ci[i,:] = cat['Lsmaj_ci'].values
#    Lsmin_ci[i,:] = cat['Lsmin_ci'].values
#    theta_ci[i,:] = cat['theta_ci'].values
#    g_ci[i,:] = cat['g_ci'].values

    coefElev = ut_solv(time, ua[:, i], [], uvnodell[i, 1],
                    cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                    nodiagn=True, linci=True, conf_int=True)

    opt = pd.DataFrame(coefElev['aux']['opt'].items())
    del coefElev['aux']['opt']
    aux = pd.DataFrame(coefElev['aux'])
    del coefElev['aux']
    c = pd.DataFrame(coefElev)
    cat = pd.concat([c,aux], axis=1)

    A.iloc[i] = coefElev['A']
    A_ci.iloc[i] = coefElev['A_ci']
    gA.iloc[i] = coefElev['g']
    gA_ci.iloc[i] = coefElev['g_ci']
    nameA.iloc[i] = coefElev['name']


data = nc.Dataset('coef.nc', 'w', format='NETCDF4')
data.createDimension('dim', None)
data.createDimension('dimx', len(x))
data.createDimension('dimtime', len(time))
data.createDimension('dimtri', trinodes.shape[-1])

newx = data.createVariable('x', 'f8', ('dimx',))
newx[:] = x
newy = data.createVariable('y', 'f8', ('dimx',))
newy[:] = y
newxc = data.createVariable('xc', 'f8', ('dim',))
newxc[:] = xc
newyc = data.createVariable('yc', 'f8', ('dim',))
newyc[:] = yc
newlon = data.createVariable('lon', 'f8', ('dimx',))
newlon[:] = lon
newlat = data.createVariable('lat', 'f8', ('dimx',))
newlat[:] = lat
newlonc = data.createVariable('lonc', 'f8', ('dim',))
newlonc[:] = lonc
newlatc = data.createVariable('latc', 'f8', ('dim',))
newlatc[:] = latc
newh = data.createVariable('h', 'f8', ('dimx',))
newh[:] = h
newtime = data.createVariable('time', 'f8', ('dimtime',))
newtime[:] = time
newtrinodes = data.createVariable('trinodes', 'f8', ('dim','dimtri'))
newtrinodes[:] = trinodes

data.createDimension('dim2', Lsmaj.shape[-1])
#data.createDimension('optDim', len(opt))
newLsmaj = data.createVariable('Lsmaj', 'f8', ('dim','dim2'))
newLsmaj_ci = data.createVariable('Lsmaj_ci', 'f8', ('dim','dim2'))
newLsmin = data.createVariable('Lsmin', 'f8', ('dim','dim2'))
newLsmin_ci = data.createVariable('Lsmin_ci', 'f8', ('dim','dim2'))
newg = data.createVariable('g', 'f8', ('dim','dim2'))
newg_ci = data.createVariable('g_ci', 'f8', ('dim','dim2'))
newtheta = data.createVariable('theta', 'f8', ('dim','dim2'))
newtheta_ci = data.createVariable('theta_ci', 'f8', ('dim','dim2'))
newname = data.createVariable('name', 'c', ('dim','dim2'))

newA = data.createVariable('A', 'f8', ('dim','dim2'))
newA_ci = data.createVariable('A_ci', 'f8', ('dim','dim2'))
newgA = data.createVariable('gA', 'f8', ('dim','dim2'))
newgA_ci = data.createVariable('gA_ci', 'f8', ('dim','dim2'))
newnameA = data.createVariable('nameA', 'c', ('dim','dim2'))

newLsmaj[:] = Lsmaj.values
newLsmaj_ci[:] = Lsmaj_ci.values
newLsmin[:] = Lsmin.values
newLsmin_ci[:] = Lsmin_ci.values
newtheta[:] = theta.values
newtheta_ci[:] = theta_ci.values
newg[:] = g.values
newg_ci[:] = g_ci.values
newname = name.values

newA[:] = A.values
newA_ci[:] = A_ci.values
newgA[:] = gA.values
newgA_ci = gA_ci.values
newnameA = nameA.values


#
#
#for i, (lonc,latc) in enumerate(lonclatc):
#    print i
#
#    coef = ut_solv(time, ua[:, i], va[:, i], uvnodell[i, 1],
#                    cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
#                    nodiagn=True, linci=True, conf_int=True)
#
#    opt = pd.DataFrame(coef['aux']['opt'].items())
#    del coef['aux']['opt']
#    aux = pd.DataFrame(coef['aux'])
#    del coef['aux']
#    c = pd.DataFrame(coef)
#
#    cat = pd.concat([c,aux], axis=1)
#
#    Lsmaj[i,:] = cat['Lsmaj'].values
#    Lsmin[i,:] = cat['Lsmin'].values
#    g[i,:] = cat['g'].values
#    theta[i, :] = cat['theta'].values
#    name[i, :] = cat['name'].values
#    Lsmaj_ci[i,:] = cat['Lsmaj_ci'].values
#    Lsmin_ci[i,:] = cat['Lsmin_ci'].values
#    theta_ci[i,:] = cat['theta_ci'].values
#    g_ci[i,:] = cat['g_ci'].values
#
#    coefElev = ut_solv(time, ua[:, i], [], uvnodell[i, 1],
#                    cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
#                    nodiagn=True, linci=True, conf_int=True)
#
#    opt = pd.DataFrame(coefElev['aux']['opt'].items())
#    del coefElev['aux']['opt']
#    aux = pd.DataFrame(coefElev['aux'])
#    del coefElev['aux']
#    c = pd.DataFrame(coefElev)
#    cat = pd.concat([c,aux], axis=1)
#
#    A[i, :] = cat['A'].values
#    gA[i, :] = cat['g'].values
#    nameA[i, :] = cat['name'].values
#    A_ci[i, :] = cat['A_ci'].values
#    gA_ci[i, :] = cat['g_ci'].values
