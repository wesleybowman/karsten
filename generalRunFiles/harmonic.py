from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib.pyplot as plt
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


# filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'
# filename = '/home/abalzer/scratch/standard_run_directory/0.0015/output/dngrid_0001.nc'
filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
#filename = '/home/abalzer/standard_run_directory/0.0015/output/dngrid_0001.nc'

data = nc.Dataset(filename, 'r')
x = data.variables['x'][:]
y = data.variables['y'][:]
lonc = data.variables['lonc'][:]
latc = data.variables['latc'][:]
ua = data.variables['ua']
va = data.variables['va']
time = data.variables['time'][:]
trinodes = data.variables['nv'][:]

time = mjd2num(time)

Rayleigh = np.array([1])

# adcpFilename = '/home/wesley/github/karsten/adcp/dngrid_adcp_2012.txt'
adcpFilename = '/home/wesley/github/karsten/adcp/testADCP.txt'

#adcpFilename = '/home/wesleyb/github/karsten/adcp/dngrid_adcp_2012.txt'
adcp = pd.read_csv(adcpFilename)

lonclatc = np.array([adcp['Longitude'], adcp['Latitude']]).T

index = closest_point(lonclatc, lonc, latc)

adcpData = pd.DataFrame()
runData = pd.DataFrame()
bottomfriction = '{0}'.format(filename.split('/')[-3])

for i, ii in enumerate(index):

    print adcp.iloc[i, 0]
    path = adcp.iloc[i, -1]
    if path != 'None':
        ADCP = pd.read_csv(path, index_col=0)
        ADCP.index = pd.to_datetime(ADCP.index)

        adcpTime = np.empty(ADCP.index.shape)

        for j, jj in enumerate(ADCP.index):
            adcpTime[j] = datetime2matlabdn(jj)

#        adcpCoef = ut_solv(adcpTime, ADCP['u'].values, ADCP['v'].values, uvnodell[ii, 1],
#                           'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
#                           'NoDiagn', 'LinCI')

        order = ['M2','S2','N2','K2','K1','O1','P1','Q1']

        adcpCoef = ut_solv(adcpTime, ADCP['u'].values, ADCP['v'].values,
                           lonclatc[i, 1],
                           cnstit=order, rmin=Rayleigh[0], notrend=True,
                           method='ols', nodiagn=True, linci=True,
                           conf_int=False, ordercnstit='frq')

        adcpAUX = adcpCoef['aux']
        del adcpAUX['opt']
        del adcpCoef['aux']

        adcpAUX = pd.DataFrame(adcpAUX)
        a = pd.DataFrame(adcpCoef)
        size = a.shape[0]
        nameSpacer = pd.DataFrame({'ADCP_Location': np.repeat(adcp.iloc[i, 0],
                                                              size)})

        bottomName = pd.DataFrame({'bottomFriction': np.repeat(bottomfriction,
                                                              size)})

        longitude = pd.DataFrame({'lon':np.repeat(lonclatc[i, 0], 8)})

        cat = pd.concat([a, adcpAUX, longitude, nameSpacer, bottomName], axis=1)

        cat = cat.set_index('ADCP_Location')
        adcpData = pd.concat([adcpData, cat])


#        coef = ut_solv(time, ua[:, ii], va[:, ii], uvnodell[ii, 1],
#                       'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
#                       'NoDiagn', 'LinCI')


        coef = ut_solv(time, ua[:, ii], va[:, ii], lonclatc[i, 1],
                        cnstit=order, rmin=Rayleigh[0], notrend=True, method='ols',
                        nodiagn=True, linci=True, conf_int=False,
                       ordercnstit='frq')


        aux = coef['aux']
        del aux['opt']
        del coef['aux']

        aux = pd.DataFrame(aux)
        c = pd.DataFrame(coef)
        size = c.shape[0]
#        nameSpacer = pd.DataFrame({'ADCP_Location': np.repeat(adcp.iloc[i, 0],
#                                                              size)})
#
#        bottomName = pd.DataFrame({'bottomFriction': np.repeat(bottomfriction,
#                                                              size)})

        ccat = pd.concat([c, aux, longitude, nameSpacer, bottomName], axis=1)
        ccat = ccat.set_index('ADCP_Location')

        runData = pd.concat([runData, ccat])


outputName = 'runData{0}.csv'.format(filename.split('/')[-3])
runData.to_csv(outputName, index_label='ADCP_Location')
outputName = 'adcpData{0}.csv'.format(filename.split('/')[-3])
adcpData.to_csv(outputName, index_label='ADCP_Location')
