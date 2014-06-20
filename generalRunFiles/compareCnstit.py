from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import cPickle as pickle
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


filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
#filename = '/home/abalzer/standard_run_directory/0.0015/output/dngrid_0001.nc'
#filename = '/array/data1/rkarsten/dncoarse_bctest/output/dn_coarse_0001.nc'

data = nc.Dataset(filename, 'r')
x = data.variables['x'][:]
y = data.variables['y'][:]
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
lonc = data.variables['lonc'][:]
latc = data.variables['latc'][:]
ua = data.variables['ua']
va = data.variables['va']
time = data.variables['time'][:]
trinodes = data.variables['nv'][:]

#(nodexy, uvnodexy, dt, deltat,
# hour, thour, TP, rho, g, period,
# nodell, uvnodell, trinodes) = ncdatasort(x, y, time*24*3600,
#                                          trinodes, lon, lat)

time = mjd2num(time)

Rayleigh = np.array([1])

# adcpFilename = '/home/wesley/github/karsten/adcp/dngrid_adcp_2012.txt'
adcpFilename = '/home/wesley/github/karsten/adcp/testADCP.txt'

#adcpFilename = '/home/wesleyb/github/karsten/adcp/dngrid_adcp_2012.txt'
#adcpFilename = '/home/107002b/github/karsten/adcp/dngrid_adcp_2012.txt'
adcp = pd.read_csv(adcpFilename)

lonlat = np.array([adcp['Longitude'], adcp['Latitude']]).T

#index = closest_point(lonlat, lon, lat)
index = closest_point(lonlat, lonc, latc)

adcpData = pd.DataFrame()
runData = pd.DataFrame()


struct = np.array([])

for i, ii in enumerate(index):

    path = adcp.iloc[i, -1]
    if path != 'None':
        print adcp.iloc[i, 0]
        #print lonlat[i,1], uvnodell[ii,1]

        ADCP = pd.read_csv(path, index_col=0)
        ADCP.index = pd.to_datetime(ADCP.index)

        adcpTime = np.empty(ADCP.index.shape)

        for j, jj in enumerate(ADCP.index):
            adcpTime[j] = datetime2matlabdn(jj)

        adcpCoef = ut_solv(adcpTime, ADCP['u'].values, ADCP['v'].values, lonlat[i, 1],
                            cnstit='auto', rmin=Rayleigh[0], notrend=True,
                            method='ols', nodiagn=True, linci=True,
                            conf_int=False)

#        adcpAUX = adcpCoef['aux']
#        del adcpAUX['opt']
#        del adcpCoef['aux']

#        adcpAUX = pd.DataFrame(adcpAUX)
#        a = pd.DataFrame(adcpCoef)
#        size = a.shape[0]
#        nameSpacer = pd.DataFrame({'ADCP_Location': np.repeat(adcp.iloc[i, 0],
#                                                                size)})
#
#        cat = pd.concat([a, adcpAUX, nameSpacer], axis=1)
#
#        adcpData = cat.set_index('ADCP_Location')
        adcpData = adcpCoef
        #adcpData = pd.concat([adcpData, cat])

        coef = ut_solv(time, ua[:, ii], va[:, ii], lonlat[i, 1],
                        cnstit='auto', rmin=Rayleigh[0], notrend=True,
                       method='ols', nodiagn=True, linci=True, conf_int=False)

#        aux = coef['aux']
#        del aux['opt']
#        del coef['aux']
#
#        aux = pd.DataFrame(aux)
#        c = pd.DataFrame(coef)
#        size = c.shape[0]
#        nameSpacer = pd.DataFrame({'ADCP_Location': np.repeat(adcp.iloc[i, 0],
#                                                              size)})
#
#        ccat = pd.concat([c, aux, nameSpacer], axis=1)
#        runData = ccat.set_index('ADCP_Location')
        runData = coef


        mod = pd.DataFrame({'ua':ua[:, i], 'va':va[:, i]})
        obs = pd.DataFrame({'u':ADCP['u'].values, 'v':ADCP['v'].values})

        obs_loc = {'name':adcp.iloc[i,0], 'type':'ADCP', 'lat':lonlat[i,-1],
                   'lon':lonlat[0,0], 'obs_timeseries':obs,
                   'mod_timeseries':mod, 'obs_time':adcpTime,
                   'mod_time':time,'speed_obs_harmonics':adcpData,
                   'speed_mod_harmonics':runData}


        struct = np.hstack((struct, obs_loc))

pickle.dump(struct, open("struct.p", "wb"))
