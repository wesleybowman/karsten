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



def tideGauge(datafiles, struct):
    dgFilename = '/array/rkarsten/common_tidal_files/data/observed/DG/TideGauge/DigbyWharf_015893_20140115_2221_Z.mat'

    gpFilename = '/array/rkarsten/common_tidal_files/data/observed/GP/TideGauge/Westport_015892_20140325_1212_Z.mat'

    dgtg = sio.loadmat(dgFilename, struct_as_record=False, squeeze_me=True)
    gptg = sio.loadmat(gpFilename, struct_as_record=False, squeeze_me=True)

    ut_constits = ['M2','S2','N2','K2','K1','O1','P1','Q1']

    coef_gptg = ut_solv(gptg.RBR.date_num_Z,
                        (gptg.RBR.data-np.mean(gptg.RBR.data)), [],
                        gptg.RBR.lat, cnstit=ut_constits, notrend=True,
                        rmin=0.95, method='ols', nodiagn=True, linci=True,
                        ordercnstit='frq')

    coef_dgtg = ut_solv(dgtg.RBR.date_num_Z,
                        (dgtg.RBR.data-np.mean(dgtg.RBR.data)), [],
                        dgtg.RBR.lat, cnstit=ut_constits, notrend=True,
                        rmin=0.95, method='ols', nodiagn=True, linci=True,
                        ordercnstit='frq')

    for filename in datafiles:

        data = nc.Dataset(filename, 'r')
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        time = data.variables['time'][:]

        tg_gp_id = np.argmin(np.sqrt((lon-gptg.RBR.lon)**2+(lat-gptg.RBR.lat)**2))
        tg_dg_id = np.argmin(np.sqrt((lon-dgtg.RBR.lon)**2+(lat-dgtg.RBR.lat)**2))

        elgp = data.variables['zeta'][tg_gp_id, :]
        eldg = data.variables['zeta'][tg_dg_id, :]

        time = mjd2num(time)

        obs_loc = {'mod_time':time, 'obs_time':dgtg.RBR.date_num_Z,
                  'lon':lon, 'lat':lat,
                  'dg_tidegauge_harmonics': coef_dgtg,
                  'gp_tidegauge_harmonics':coef_gptg}

        struct = np.hstack((struct, obs_loc))


    #pickle.dump(struct, open("structADCP.p", "wb"))
    return struct






def adcp(datafiles):

    for filename in datafiles:
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

        time = mjd2num(time)

        Rayleigh = np.array([1])

        #adcpFilename = '/home/wesley/github/karsten/adcp/dngrid_adcp_2012.txt'
        #adcpFilename = '/home/wesley/github/karsten/adcp/testADCP.txt'

        #adcpFilename = '/home/wesleyb/github/karsten/adcp/dngrid_adcp_2012.txt'
        adcpFilename = '/array/home/107002b/github/karsten/adcp/acadia_dngrid_adcp_2012.txt'
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

                adcpData = adcpCoef

                coef = ut_solv(time, ua[:, ii], va[:, ii], lonlat[i, 1],
                                cnstit='auto', rmin=Rayleigh[0], notrend=True,
                            method='ols', nodiagn=True, linci=True, conf_int=False)

                runData = coef

                mod = pd.DataFrame({'ua':ua[:, i], 'va':va[:, i]})
                obs = pd.DataFrame({'u':ADCP['u'].values, 'v':ADCP['v'].values})

                obs_loc = {'name':adcp.iloc[i,0], 'type':'ADCP', 'lat':lonlat[i,-1],
                        'lon':lonlat[0,0], 'obs_timeseries':obs,
                        'mod_timeseries':mod, 'obs_time':adcpTime,
                        'mod_time':time,'speed_obs_harmonics':adcpData,
                        'speed_mod_harmonics':runData}


                struct = np.hstack((struct, obs_loc))

#    pickle.dump(struct, open("structADCP.p", "wb"))
    return struct


def main():
    datafiles = ['array/data1/rkarsten/dncoarse_bctest_old/output/dn_coarse_0001.nc',
                 '/array/data1/rkarsten/dncoarse_bctest/output/dn_coarse_0001.nc',
                 '/array/data1/rkarsten/dncoarse_bctest2/output/dn_coarse_0001.nc',
                 '/array/data1/rkarsten/dncoarse_bctest_all/output/dn_coarse_0001.nc',
                 '/array/data1/rkarsten/dncoarse_bctest_EC/output/dn_coarse_0001.nc',
                 '/array/data1/rkarsten/dncoarse_bctest_timeseries/output/dn_coarse_0001.nc',
                 '/array/data1/rkarsten/dncoarse_stationtest/output/dn_coarse_0001.nc']

    struct = adcp(datafiles)
    struct = tideGauge(datafiles, struct)
