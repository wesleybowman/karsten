from __future__ import division
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import cPickle as pickle
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv
import scipy.io as sio
from stationClass import station
from adcpClass import ADCP
from fvcomClass import FVCOM

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




def adcp(datafiles, debug=False):

    if debug:
        adcpFilename = '/home/wesley/github/karsten/adcp/testADCP.txt'
    else:
        adcpFilename = '/array/home/107002b/github/karsten/adcp/acadia_dngrid_adcp_2012.txt'

    #adcpFilename = '/home/wesleyb/github/karsten/adcp/dngrid_adcp_2012.txt'
    adcp = pd.read_csv(adcpFilename)

    for i,v in enumerate(adcp['Latitude']):
        path = adcp.iloc[i, -1]
        if path != 'None':
            print adcp.iloc[i, 0]
            #print lonlat[i,1], uvnodell[ii,1]

            ADCP = pd.read_csv(path, index_col=0)
            ADCP.index = pd.to_datetime(ADCP.index)

            adcpTime = np.empty(ADCP.index.shape)

            for j, jj in enumerate(ADCP.index):
                adcpTime[j] = datetime2matlabdn(jj)

            adcpCoef = ut_solv(adcpTime, ADCP['u'].values,
                               ADCP['v'].values, v,
                               cnstit='auto', rmin=0.95, notrend=True,
                               method='ols', nodiagn=True, linci=True,
                               conf_int=True)

            adcpData = adcpCoef

    obs = pd.DataFrame({'u':ADCP['u'].values, 'v':ADCP['v'].values})
    Struct = {}


    for filename in datafiles:
        print filename
        data = nc.Dataset(filename, 'r')
        #x = data.variables['x'][:]
        #y = data.variables['y'][:]
        lon = data.variables['lon'][:]
        lat = data.variables['lat'][:]
        lonc = data.variables['lonc'][:]
        latc = data.variables['latc'][:]
        ua = data.variables['ua']
        va = data.variables['va']
        time = data.variables['time'][:]
        #trinodes = data.variables['nv'][:]

        time = mjd2num(time)

        lonlat = np.array([adcp['Longitude'], adcp['Latitude']]).T

        #index = closest_point(lonlat, lon, lat)
        index = closest_point(lonlat, lonc, latc)

        adcpData = pd.DataFrame()
        runData = pd.DataFrame()

        Name = filename.split('/')[-3]
        Name = '2012_run'
        print Name
        struct = np.array([])

        for i, ii in enumerate(index):

            path = adcp.iloc[i, -1]
            if path != 'None':
                print adcp.iloc[i, 0]

                coef = ut_solv(time, ua[:, ii], va[:, ii], lonlat[i, 1],
                                cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True,
                               conf_int=True)

                runData = coef

                mod = pd.DataFrame({'ua':ua[:, ii], 'va':va[:, ii]})

                obs_loc = {'name':adcp.iloc[i,0], 'type':'ADCP', 'lat':lonlat[i,-1],
                        'lon':lonlat[0,0], 'obs_timeseries':obs,
                        'mod_timeseries':mod, 'obs_time':adcpTime,
                        'mod_time':time,'speed_obs_harmonics':adcpData,
                        'speed_mod_harmonics':runData}


                struct = np.hstack((struct, obs_loc))

        Struct[Name] = struct

    return Struct


def main(debug=False):
    if debug:
        #datafiles = ['/array/data1/rkarsten/dncoarse_bctest_old/output/dn_coarse_0001.nc',
        #            '/array/data1/rkarsten/dncoarse_bctest/output/dn_coarse_0001.nc']
        fvFiles = ['/home/wesley/ncfiles/smallcape_force_0001.nc']
    else:

        #fvFile = '/EcoII/EcoEII_server_data_tree/data/simulated/FVCOM/dngrid/june_2013_3D/output/'
        #fvFiles = ['/EcoII/EcoEII_server_data_tree/data/simulated/FVCOM/dngrid/june_2013_3D/output/']
        fvFiles = ['/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/']
        #adcpFile = '/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPa_avg5.mat'
        #adcpFile = '/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPb_avg5.mat'
        adcpFiles = ['/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPa_avg5.mat',
         '/EcoII/EcoEII_server_data_tree/data/observed/GP/ADCP/Flow_GP-130620-BPb_avg5.mat']
	fvdebug = '/EcoII/EcoEII_server_data_tree/workspace/simulated/FVCOM/dngrid/june_2013_3D/output/dngrid_0001_week2.nc'

    fvdebugData = FVCOM(fvdebug)
    saveName = 'june_2013_3D_station.p'
    Name = 'june_2013_3D_station'
    Struct = {}

    for fvFile in fvFiles:

        print fvFile
        fvData = station(fvFile)

        struct = np.array([])
        for adcpFile in adcpFiles:
            print adcpFile
            adcpData = ADCP(adcpFile)
            lonlat = np.array([adcpData.lon[0], adcpData.lat[0]]).T
            #lonlat = np.array([adcpData.x[0], adcpData.y[0]]).T
            #ind = closest_point(lonlat, fvData.lon, fvData.lat)
            newind = closest_point(lonlat, fvdebugData.lon, fvdebugData.lat)
            #ind = closest_point(lonlat, fvData.x, fvData.y)
            new = np.array([fvdebugData.xc[newind], fvdebugData.yc[newind]])
            ind = closest_point(new, fvData.x, fvData.y)

            print ind
            print adcpData.mtime.shape
            print adcpData.ua.shape
            print adcpData.va.shape
            print adcpData.surf.shape

            adcpVelCoef = ut_solv(adcpData.mtime, adcpData.ua,
                            adcpData.va, adcpData.lat[0],
                            cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, coef_int=True)

            adcpElevCoef = ut_solv(adcpData.mtime, adcpData.surf,
                            [], adcpData.lat[0],
                            cnstit='auto', rmin=0.95, notrend=True,
                            method='ols', nodiagn=True, linci=True, coef_int=True)

            adcpName = adcpFile.split('/')[-1].split('.')[0]
            #WB_COMMENT: Doesn't currently work
            obs = pd.DataFrame({'u':adcpData.ua, 'v':adcpData.va, 'elev':adcpData.surf})

            print fvData.time.shape
            print fvData.ua[:, ind].shape
            print fvData.va[:, ind].shape
            print fvData.lat[ind].shape

            fvVelCoef = ut_solv(fvData.time, fvData.ua[:, ind].flatten(),
                                fvData.va[:, ind].flatten(),
                        adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                        method='ols', nodiagn=True, linci=True, conf_int=True)

            print fvData.elev[:, ind].shape
            fvElevCoef = ut_solv(fvData.time, fvData.elev[:, ind].flatten(), [],
                        adcpData.lat[0], cnstit='auto', rmin=0.95, notrend=True,
                        method='ols', nodiagn=True, linci=True, conf_int=True)

            mod = pd.DataFrame({'ua':fvData.ua[:, ind].flatten(),
                                'va':fvData.va[:, ind].flatten(),
                                'elev':fvData.elev[:, ind].flatten()})


            obs_loc = {'name': adcpName, 'type':'ADCP', 'lat':fvdebugData.lat[newind],
                    'lon':fvdebugData.lon[newind], 'obs_timeseries':obs,
                    'mod_timeseries':mod, 'obs_time':adcpData.mtime,
                    'mod_time':fvData.time, 'vel_obs_harmonics':adcpVelCoef,
                    'elev_obs_harmonics':adcpElevCoef,
                    'vel_mod_harmonics':fvVelCoef, 'elev_mod_harmonics':fvElevCoef}

            struct = np.hstack((struct, obs_loc))

        Struct[Name] = struct

    if debug:
        pickle.dump(Struct, open("structADCP.p", "wb"))

    pickle.dump(Struct, open(saveName, "wb"))
    return Struct
    #return fvData, adcpData

if __name__ == '__main__':
    Struct = main()
