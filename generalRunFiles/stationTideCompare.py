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



def tideGauge(datafiles, Struct):
    dgFilename = '/array/home/rkarsten/common_tidal_files/data/observed/DG/TideGauge/DigbyWharf_015893_20140115_2221_Z.mat'

    gpFilename = '/array/home/rkarsten/common_tidal_files/data/observed/GP/TideGauge/Westport_015892_20140325_1212_Z.mat'

    dgtg = sio.loadmat(dgFilename, struct_as_record=False, squeeze_me=True)
    gptg = sio.loadmat(gpFilename, struct_as_record=False, squeeze_me=True)

    ut_constits = ['M2','S2','N2','K2','K1','O1','P1','Q1']

    print 'Westport TideGauge'
    coef_gptg = ut_solv(gptg['RBR'].date_num_Z,
                        (gptg['RBR'].data-np.mean(gptg['RBR'].data)), [],
                        gptg['RBR'].lat, cnstit=ut_constits, notrend=True,
                        rmin=0.95, method='ols', nodiagn=True, linci=True,
                        ordercnstit='frq')

    print 'DigbyWharf TideGauge'
    coef_dgtg = ut_solv(dgtg['RBR'].date_num_Z,
                        (dgtg['RBR'].data-np.mean(dgtg['RBR'].data)), [],
                        dgtg['RBR'].lat, cnstit=ut_constits, notrend=True,
                        rmin=0.95, method='ols', nodiagn=True, linci=True,
                        ordercnstit='frq')

    struct = np.array([])
    for filename in datafiles:

        print filename
        data = nc.Dataset(filename, 'r')
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        time_JD = data.variables['time_JD'][:]
        time_second = data.variables['time_second'][:]
        time = time_JD + 678942 + time_second / (24*3600)

        #time = mjd2num(time)

        tg_gp_id = np.argmin(np.sqrt((lon-gptg['RBR'].lon)**2+(lat-gptg['RBR'].lat)**2))
        tg_dg_id = np.argmin(np.sqrt((lon-dgtg['RBR'].lon)**2+(lat-dgtg['RBR'].lat)**2))

        #elgp = data.variables['zeta'][tg_gp_id, :]
        #eldg = data.variables['zeta'][tg_dg_id, :]
        elgp = data.variables['zeta'][:, tg_gp_id]
        eldg = data.variables['zeta'][:, tg_dg_id]

        coef_dg = ut_solv(time, eldg, [], dgtg['RBR'].lat, cnstit=ut_constits,
                          notrend=True, rmin=0.95, method='ols', nodiagn=True,
                          linci=True, ordercnstit='frq')

        coef_gp = ut_solv(time, elgp, [], gptg['RBR'].lat, cnstit=ut_constits,
                          notrend=True, rmin=0.95, method='ols', nodiagn=True,
                          linci=True, ordercnstit='frq')


        Name = filename.split('/')[-3]
        Name = '2012_station_run'

        print Name

        obs_loc = {'name':Name, 'type':'TideGauge',
                   'mod_time':time, 'dg_time':dgtg['RBR'].date_num_Z,
                   'gp_time':gptg['RBR'].date_num_Z,
                   'lon':lon, 'lat':lat,
                   'dg_tidegauge_harmonics': coef_dgtg,
                   'gp_tidegauge_harmonics':coef_gptg,
                   'dg_mod_harmonics': coef_dg,
                   'gp_mod_harmonics': coef_gp,
                   'dg_tg_data':dgtg['RBR'].data,
                   'gp_tg_data':gptg['RBR'].data,
                   'eldg':eldg, 'elgp':elgp}

        struct = np.hstack((struct, obs_loc))

        Struct[Name] = np.hstack((Struct[Name], struct))


    #pickle.dump(struct, open("structADCP.p", "wb"))
    return Struct


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
        #lat = data.variables['lat'][:]
        #lon = data.variables['lon'][:]
        time_JD = data.variables['time_JD'][:]
        time_second = data.variables['time_second'][:]
        time = time_JD + 678942 + time_second/(24*3600)

        lonc = data.variables['lon'][:]
        latc = data.variables['lat'][:]
        ua = data.variables['ua']
        va = data.variables['va']
        #trinodes = data.variables['nv'][:]

        #time = mjd2num(time)

        lonlat = np.array([adcp['Longitude'], adcp['Latitude']]).T

        #index = closest_point(lonlat, lon, lat)
        index = closest_point(lonlat, lonc, latc)

        adcpData = pd.DataFrame()
        runData = pd.DataFrame()

        Name = filename.split('/')[-3]
        Name = '2012_station_run'

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
        datafiles = ['/array/data1/rkarsten/dncoarse_bctest_old/output/dn_coarse_0001.nc',
        '/array/data1/rkarsten/dncoarse_bctest/output/dn_coarse_0001.nc']
        #datafiles = ['/home/wesley/ncfiles/smallcape_force_0001.nc']
    else:

        datafiles = ['/array2/data3/rkarsten/dncoarse_3D/output2/dn_coarse_station_timeseries.nc']
        datafiles = ['/array/home/rkarsten/common_tidal_files/data/simulated/FVCOM/dn_coarse_station_timeseries.nc']

        datafiles = ['/array/home/116822s/2012_station_run.nc']

    saveName = 'struct2012_station_run.p'
    Struct = adcp(datafiles, debug=False)

    if debug:
        pickle.dump(Struct, open("structADCP.p", "wb"))

    Struct = tideGauge(datafiles, Struct)
    pickle.dump(Struct, open(saveName, "wb"))
    return Struct

if __name__ == '__main__':
    main()
