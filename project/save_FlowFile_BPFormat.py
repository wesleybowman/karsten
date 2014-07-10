from __future__ import division
from rawADCPclass import rawADCP
from datetime import datetime
from datetime import timedelta

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
    yd = adcp.mtime[:] - datenum

    mini = timedelta(days=params.tmin)
    maxi = timedelta(days=params.tmax)

    nmin = datetime(day1.year,1,1) + mini
    nmax = datetime(day1.year,1,1) + maxi



    print yd



if __name__ == '__main__':
    filename = '140703-EcoEII_database/data/GP-120726-BPd_raw.mat'
    data = rawADCP(filename)
    adcp = Struct(**data.adcp)
    params = Struct(**data.saveparams)

#    save_FlowFile_BPFormat(data.fileinfo, data.adcp, data.rbr,
#                           data.saveparams, data.options)

    save_FlowFile_BPFormat(data.fileinfo, adcp, data.rbr,
                           params, data.options)
