import netCDF4 as nc
#import pandas as pd


def coefNC2D(c, coefName):
    try:
        data = nc.Dataset('coef.nc', 'a', format='NETCDF4')
    except RuntimeError:
        data = nc.Dataset('coef.nc', 'w', format='NETCDF4')


    coefgrp = data.createGroup(coefName)
    coefgrp.createDimension('dim', None)

    Lsmaj = coefgrp.createVariable('Lsmaj', 'f8', ('dim',))
    Lsmaj[:] = c['Lsmaj'].values

    Lsmin = coefgrp.createVariable('Lsmin', 'f8', ('dim',))
    Lsmin[:] = c['Lsmin'].values

    g = coefgrp.createVariable('g', 'f8', ('dim',))
    g[:] = c['g'].values

    names = coefgrp.createVariable('names', 'c', ('dim',))
    names[:] = c['name'].values

    theta = coefgrp.createVariable('theta', 'f8', ('dim',))
    theta[:] = c['theta'].values

    frq = coefgrp.createVariable('frq', 'f8', ('dim',))
    frq[:] = c['frq'].values

    #newx = coefgrp.createVariable('x', 'f8', ('dim',))
    #newx[:] = x
    #
    #newy = coefgrp.createVariable('y', 'f8', ('dim',))
    #newy[:] = y
    #
    #newtime = coefgrp.createVariable('time', 'f8', ('dim',))
    #newtime[:] = time
    #
    #newlon = coefgrp.createVariable('x', 'f8', ('dim',))
    #newlon[:] = x
    #
    #newlat = coefgrp.createVariable('x', 'f8', ('dim',))
    #newlat[:] = x
    #
    #newlonc = coefgrp.createVariable('x', 'f8', ('dim',))
    #newlonc[:] = x
    #
    #newlatc = coefgrp.createVariable('x', 'f8', ('dim',))
    #newlatc[:] = x






    data.close()
