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


def chunk(size, n):
    """ Yield successive n-sized chunks from l.
    """
    l = range(size)
    for i in xrange(0, size, n):
        if i+n in l:
            yield i, i+n
        else:
            yield i, l[-1]


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
#filename = '/home/abalzer/scratch/standard_run_directory/0.0015/output/dngrid_0001.nc'

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
#ua = data.variables['ua'][:]
#va = data.variables['va'][:]
time = data.variables['time'][:]
trinodes = data.variables['nv'][:]
#siglay = data.variables['siglay'][:]
#siglev = data.variables['siglev'][:]

(nodexy, uvnodexy, dt, deltat,
 hour, thour, TP, rho, g, period,
 nodell, uvnodell, trinodes) = ncdatasort(x, y, time*24*3600,
                                          trinodes, lon, lat)

time = mjd2num(time)

Rayleigh = np.array([1])

lonlat = np.array([lon, lat]).T
lonclatc = np.array([lonc, latc]).T

#size = lonclatc.shape[0]
#chunks = list(chunk(s,80000))

#coefName = 'coef{0}'.format(i)
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

for i, (lonc,latc) in enumerate(lonclatc):
#for i in range(0,10):
    print i

    coef = ut_solv(time, ua[:, i], va[:, i], uvnodell[i, 1],
                    'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
                    'NoDiagn', 'LinCI')


    opt = pd.DataFrame(coef['aux']['opt'].items())
    del coef['aux']['opt']
    aux = pd.DataFrame(coef['aux'])
    del coef['aux']
    c = pd.DataFrame(coef)

    cat = pd.concat([c,aux], axis=1)
    #coefnc = data.createVariable('coef', 'f8', ('dim',))
    #coefnc[:] = cat
    try:
        data.createDimension('dim2', cat['Lsmaj'].shape[0])
        Lsmaj = data.createVariable('Lsmaj', 'f8', ('dim','dim2'))
        Lsmin = data.createVariable('Lsmin', 'f8', ('dim','dim2'))
        g = data.createVariable('g', 'f8', ('dim','dim2'))
        theta = data.createVariable('theta', 'f8', ('dim','dim2'))
    except RuntimeError:
        pass

    Lsmaj[i,:] = c['Lsmaj'].values
    Lsmin[i,:] = c['Lsmin'].values
    g[i,:] = c['g'].values
    theta[i,:] = c['theta'].values

#coefName = 'coef{0}'.format(i)
#coefNC2D(cat, coefName)


#aux = coef['aux']
#del aux['opt']
#del coef['aux']
#
#aux = pd.DataFrame(aux)
#c = pd.DataFrame(coef)
#c = pd.concat([c, aux], axis=1)
# c['aux'] = pd.Series(c['aux'])

#runData = pd.concat([runData, nameSpacer])
#runData = pd.concat([runData, c])

# name = '{0}'.format(adcp.iloc[i,0])
# adcpData.to_hdf('adcpData.h5', name, mode='a')

#runData.to_hdf('runData.h5', 'runData', mode='a')
#runData.to_csv('runData.csv')
#adcpData.to_hdf('adcpData.h5', 'runData', mode='a')
#adcpData.to_csv('adcpData.csv')
