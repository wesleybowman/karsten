from datetime import datetime, timedelta


def mjd2num(x):

    y = x + 678942

    return y


def datenum2date(matlab_datenum):

    python_datetime = datetime.fromordinal(int(matlab_datenum)) + \
        timedelta(days=matlab_datenum % 1) - timedelta(days=366)

    return python_datetime
