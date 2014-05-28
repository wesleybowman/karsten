from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv


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


#filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'
#filename = '/home/abalzer/scratch/standard_run_directory/0.0015/output/dngrid_0001.nc'
filename = '/home/abalzer/standard_run_directory/0.0015/output/dngrid_0001.nc'

data = nc.Dataset(filename, 'r')
x = data.variables['x'][:]
y = data.variables['y'][:]
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
ua = data.variables['ua']
va = data.variables['va']
time = data.variables['time'][:]
trinodes = data.variables['nv'][:]

(nodexy, uvnodexy, dt, deltat,
 hour, thour, TP, rho, g, period,
 nodell, uvnodell, trinodes) = ncdatasort(x, y, time*24*3600,
                                          trinodes, lon, lat)

time = mjd2num(time)

#Rayleigh = np.array([0.97, 1])
Rayleigh = np.array([1])

#adcpFilename = '/home/wesley/github/karsten/adcp/dngrid_adcp_2012.txt'

adcpFilename = '/home/wesleyb/github/karsten/adcp/dngrid_adcp_2012.txt'
adcp = pd.read_csv(adcpFilename)

lonlat = np.array([adcp['Longitude'], adcp['Latitude']]).T

index = closest_point(lonlat, lon, lat)

# Need to do DataFrame instead of Series
adcpData = pd.DataFrame()
runData = pd.DataFrame()
# runData = pd.DataFrame()

for i, ii in enumerate(index):

    path = adcp.iloc[i, -1]
    if path != 'None':
        ADCP = pd.read_csv(path, index_col=0)
        ADCP.index = pd.to_datetime(ADCP.index)

        adcpTime = np.empty(ADCP.index.shape)

        for j, jj in enumerate(ADCP.index):
            adcpTime[j] = datetime2matlabdn(jj)

        adcpCoef = ut_solv(time, ua[:, ii], va[:, ii], uvnodell[ii, 1],
                           'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
                           'NoDiagn', 'LinCI')

        adcpAUX = adcpCoef['aux']
        del adcpAUX['opt']
        del adcpCoef['aux']

        adcpAUX = pd.DataFrame(adcpAUX)
        a = pd.DataFrame(adcpCoef)
        a = pd.concat([a, adcpAUX], axis=1)
        # a['aux'] = pd.Series(a['aux'])

        size = a.shape[0]
        #nameSpacer = pd.DataFrame({'ADCP_Location': [adcp.iloc[i, 0]]})
        nameSpacer = pd.DataFrame({'ADCP_Location': np.repeat(adcp.iloc[i, 0],
                                                              size)})
        adcpData = pd.concat([adcpData, nameSpacer])
        adcpData = pd.concat([adcpData, a],axis=1)

        coef = ut_solv(time, ua[:, ii], va[:, ii], uvnodell[ii, 1],
                       'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
                       'NoDiagn', 'LinCI')

        aux = coef['aux']
        del aux['opt']
        del coef['aux']

        aux = pd.DataFrame(aux)
        c = pd.DataFrame(coef)
        c = pd.concat([c, aux], axis=1)
        # c['aux'] = pd.Series(c['aux'])

        runData = pd.concat([runData, nameSpacer])
        runData = pd.concat([runData, c], axis=1)

# name = '{0}'.format(adcp.iloc[i,0])
# adcpData.to_hdf('adcpData.h5', name, mode='a')

runData.to_hdf('runData.h5', 'runData', mode='a')
runData.to_csv('runData.csv')
adcpData.to_hdf('adcpData.h5', 'runData', mode='a')
adcpData.to_csv('adcpData.csv')
