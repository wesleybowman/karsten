from __future__ import division
import numpy as np
from rawADCPclass import rawADCP
from datetime import datetime
from datetime import timedelta
import scipy.io as sio
import scipy.interpolate as sip
import matplotlib.pyplot as plt

def date2py(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + \
        timedelta(days=matlab_datenum%1) - timedelta(days = 366)

    return python_datetime


def py2date(dt):
   mdn = dt + timedelta(days = 366)
   frac_seconds = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
   return mdn.toordinal() + frac_seconds + frac_microseconds

def calc_ensemble(x, ens, ens_dim):

    #initialize input
    ens = int(ens)
    #x = x[:, None]

    if ens_dim == 1:
        ens_size = np.floor(x.shape[0]/60)
    else:
        pass

    #x_ens = np.empty((ens_size, 1, ens))
    x_ens = np.empty((ens_size, ens))
    x_ens[:] = np.nan

    for j in xrange(ens):
        if ens_dim == 1:
            ind_ens = np.arange(j, x.shape[0] - (ens - j), ens)
            #x_ens[..., j] = x[ind_ens]
            x_ens[..., j] = x[ind_ens]

        else:
            pass

    #x_ens = np.nanmean(x_ens, axis=2)
    x_ens = np.nanmean(x_ens, axis=1)
    return x_ens


def rotate_coords(x, y, theta):
    '''
    Similar to "rotate_to_channelcoords.m" code,
    theta is now the angle
    between the old axis and the new x-axis (CCw is positive)
    '''

    xnew = x * np.cos(theta) + y * np.sin(theta)
    ynew = -x * np.sin(theta) + y * np.cos(theta)

    return xnew, ynew

def rotate_to_true(X, Y, theta=-19):
    '''
    % X,Y are the X and Y coordinates (could be speeds) relative to magnetic
    % north -- inputs can be vectors
    % x,y are the coordinates relative to true north
    % This function assumes the measured location is Nova Scotia where the
    % declination angle is -19 degrees.
    %
    % Sept 29, 2012: Changed print statement
    %
    % Sept 20, 2012: Modified the function to allow for theta to be input.
    % Default will remain at -19 degrees, but this may not be accurate for all
    % places in Nova Scotia.
    '''

    print 'Rotating velocities to be relative to true north (declination = {0})'.format(theta)

    Theta = theta * np.pi / 180

    x = X * np.cos(Theta) + Y * np.sin(Theta)
    y = -X * np.sin(Theta) + Y * np.cos(Theta)

    return x, y


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

def save_FlowFile_BPFormat(fileinfo, adcp, rbr, params, options):
    print adcp.mtime[0]
    day1 = date2py(adcp.mtime[0][0])
    print day1
    #date_time = [date2py(tval[0]) for tval in adcp.mtime[:]]
    datenum = datetime(day1.year,1,1) + timedelta(365)
    datenum = datenum.toordinal()

    yd = adcp.mtime[:].flatten() - datenum
    tind = np.where((yd > params.tmin) & (yd < params.tmax))[0]

    time = {}
    time['mtime'] = adcp.mtime[:].flatten()[tind]
    dt = np.nanmean(np.diff(time['mtime']))

    if not rbr:
        print 'Depths measured by ADCP not yet coded.'
    else:
        print 'Ensemble averaging rbr data'

        nens = round(dt/(rbr.mtime[1] - rbr.mtime[0]))



#    mini = timedelta(days=params.tmin)
#    maxi = timedelta(days=params.tmax)
#
#    nmin = datetime(day1.year,1,1) + mini
#    nmax = datetime(day1.year,1,1) + maxi
#    print yd



if __name__ == '__main__':
    filename = '140703-EcoEII_database/data/GP-120726-BPd_raw.mat'
    data = rawADCP(filename)
    #adcp = Struct(**data.adcp)
    rawADCP = data.adcp
    adcp = data.adcp
    params = Struct(**data.saveparams)
    params = data.saveparams
    rbr = Struct(**data.rbr)

#    save_FlowFile_BPFormat(data.fileinfo, data.adcp, data.rbr,
#                           data.saveparams, data.options)

    debug = False
    #day1 = date2py(adcp.mtime[0][0])
    day1 = date2py(adcp['mtime'][0][0])
    #date_time = [date2py(tval[0]) for tval in adcp.mtime[:]]
    datenum = datetime(day1.year,1,1) + timedelta(365)
    datenum = datenum.toordinal()
    #yd = adcp.mtime[:].flatten() - datenum
    yd = adcp['mtime'][:].flatten() - datenum
    #tind = np.where((yd > params.tmin) & (yd < params.tmax))[0]
    tind = np.where((yd > params['tmin']) & (yd < params['tmax']))[0]

    time = {}
    time['mtime'] = adcp['mtime'][:].flatten()[tind]
    dt = np.nanmean(np.diff(time['mtime']))
    pres = {}

    if not rbr:
        print 'Depths measured by ADCP not yet coded.'
    else:
        print 'Ensemble averaging rbr data'

        nens = round(dt/(rbr.mtime[1] - rbr.mtime[0]))
        temp = np.arange(rbr.mtime[nens/2-1], rbr.mtime[-1-nens/2], dt)
        temp2 = np.r_[rbr.mtime[nens/2-1]: rbr.mtime[-1-nens/2]: dt]

        # Load in matlab values
        filename = './140703-EcoEII_database/scripts_examples/mtime.mat'
        mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        matTimes = mat['mtimeens']
        filename = './140703-EcoEII_database/scripts_examples/dt.mat'
        mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        matdt = mat['dt']


        mtimeens = np.arange(rbr.mtime[nens/2-1], rbr.mtime[-1-nens/2], dt)
        #mtimeens = mtimeens + params.rbr_hr_offset / 24
        mtimeens = mtimeens + params['rbr_hr_offset'] / 24
        depthens = calc_ensemble(rbr.depth, nens, 1)

        filename = './140703-EcoEII_database/scripts_examples/depthens.mat'
        mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        matdepthens = mat['depthens']

        filename = './140703-EcoEII_database/scripts_examples/time.mat'
        mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        matmtime = mat['mtime']


        debug = True
        if debug:
            print matTimes.shape
            print temp - matTimes
            print temp2 - matTimes
            print dt - matdt
            print depthens - matdepthens
            print 'time'
            print time['mtime'] - matmtime

        temp = sip.interp1d(mtimeens, depthens, kind='linear')

        #pres['surf']= temp(time['mtime']) + params.dabPS
        pres['surf']= temp(time['mtime']) + params['dabPS']

        #if data.options['showRBRavg']:
        if debug:
            #plt.plot(rbr.mtime+params.rbr_hr_offset/24, rbr.depth+params.dabPS,
            #         label='RBR')
            #plt.plot(time['mtime'], pres['surf'], 'r', label='AVG')
            plt.plot(rbr.mtime+params['rbr_hr_offset']/24, rbr.depth+params['dabPS'],
                     label='RBR')
            plt.plot(time['mtime'], pres['surf'], 'r', label='AVG')
            plt.xlabel('Time')
            plt.ylabel('Elevation')
            plt.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=plt.gcf().transFigure)

            plt.show()

    ## zlevels
    data = {}
    z = adcp['config']['ranges'][:] + params['dabADCP']
    z = z.flatten()
    zind = np.where((z > params['zmin']) & (z < params['zmax']))[0]
    data['bins'] = z[zind]

    ## Currents
    #data['vert_vel'] = adcp['vert_vel'][:][zind, tind].T
    #data['error_vel'] = adcp['error_vel'][:][zind, tind].T
    data['vert_vel'] = adcp['vert_vel'][:][tind][:, zind]
    data['error_vel'] = adcp['error_vel'][:][tind][:, zind]

    # If compass wasn't calibrated
    #if isfield(params,'hdgmod'):
    if 'hdgmod' in params:
        adcp['east_vel'][:], adcp['north_vel'][:] = rotate_coords(adcp['east_vel'][:],
                                                              adcp['north_vel'][:],
                                                              params['hdgmod'])
        #Comments{end+1} = 'East and north velocity rotated by params.hdgmod';

    data['east_vel'], data['north_vel'] = \
        rotate_to_true(adcp['east_vel'][:][tind][:, zind],
                       adcp['north_vel'][:][tind][:, zind],
                       params['declination'])

    # Direction






#data.dir_vel = get_DirFromN(data.east_vel,data.north_vel);
#
#% Signed Speed
#spd_all = sqrt(data.east_vel.^2+data.north_vel.^2);
#% Determine flood and ebb based on principal direction (Polagye Routine)
#disp('Getting signed speed (Principal Direction Method) -- used all speeds')
#[s_signed_all, PA_all] = sign_speed(data.east_vel, data.north_vel, spd_all, data.dir_vel, params.flooddir);
#data.mag_signed_vel = s_signed_all;
#
#
##    save_FlowFile_BPFormat(data.fileinfo, adcp, data.rbr,
##                           params, data.options)
