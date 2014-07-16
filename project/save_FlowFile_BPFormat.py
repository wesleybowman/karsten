from __future__ import division
import numpy as np
from rawADCPclass import rawADCP
from datetime import datetime
from datetime import timedelta
import scipy.io as sio

def date2py(matlab_datenum):
    python_datetime = datetime.fromordinal(int(matlab_datenum)) + \
        timedelta(days=matlab_datenum%1) - timedelta(days = 366)

    return python_datetime


def py2date(dt):
   mdn = dt + timedelta(days = 366)
   frac_seconds = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
   return mdn.toordinal() + frac_seconds + frac_microseconds


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
    adcp = Struct(**data.adcp)
    params = Struct(**data.saveparams)
    rbr = Struct(**data.rbr)

#    save_FlowFile_BPFormat(data.fileinfo, data.adcp, data.rbr,
#                           data.saveparams, data.options)

    day1 = date2py(adcp.mtime[0][0])
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
        temp = np.arange(rbr.mtime[nens/2-1], rbr.mtime[-1-nens/2], dt)
        temp2 = np.r_[rbr.mtime[nens/2-1]: rbr.mtime[-1-nens/2]: dt]

        # Load in matlab values
        filename = './140703-EcoEII_database/scripts_examples/mtime.mat'
        mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        matTimes = mat['mtimeens']
        filename = './140703-EcoEII_database/scripts_examples/dt.mat'
        mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        matdt = mat['dt']

        print matTimes.shape
        print temp - matTimes
        print temp2 - matTimes
        print dt - matdt

#        mtimeens = rbr.mtime[nens/2]rdt:rbr.mtime[end-nens/2]
#        mtimeens = mtimeens+params.rbr_hr_offset/24
#
#        depthens = calc_ensemble(rbr.depth,nens,1);


#    save_FlowFile_BPFormat(data.fileinfo, adcp, data.rbr,
#                           params, data.options)
