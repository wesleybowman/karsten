from __future__ import division
import numpy as np
import netCDF4 as nc
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr


class FVCOM:
    '''
    A class for FVCOM data.
    As of right now, only takes a filename as input. It will then load in the
    data (except for timeseries, since loading in the whole time series can be
    too large)
    '''

    def __init__(self, filename):
        self.QC = ['raw data']
        self.load(filename)

    def load(self, filename):
        self.data = nc.Dataset(filename, 'r')
        self.x = self.data.variables['x'][:]
        self.y = self.data.variables['y'][:]
        self.lon = self.data.variables['lon'][:]
        self.lat = self.data.variables['lat'][:]
        self.lonc = self.data.variables['lonc'][:]
        self.latc = self.data.variables['latc'][:]
        self.h = self.data.variables['h'][:]
        self.nbe = self.data.variables['nbe'][:]
        self.a1u = self.data.variables['a1u'][:]
        self.a2u = self.data.variables['a2u'][:]
        self.aw0 = self.data.variables['aw0'][:]
        self.awx = self.data.variables['awx'][:]
        self.awy = self.data.variables['awy'][:]
        self.time = self.data.variables['time'][:]
        self.trinodes = self.data.variables['nv'][:]
        self.siglay = self.data.variables['siglay'][:]
        self.siglev = self.data.variables['siglev'][:]

        # Need to use len to get size of dimensions
        self.nele = self.data.dimensions['nele']
        self.node = self.data.dimensions['node']

        # elev timeseries
        self.el = self.data.variables['zeta']

        try:
            self.wa = self.data.variables['w']
            self.ua = self.data.variables['u']
            self.va = self.data.variables['v']
            self.D3 = True

        except KeyError:
            self.ua = self.data.variables['ua']
            self.va = self.data.variables['va']
            self.D3 = False

    def harmonics(self, ind, twodim=True, **kwarg):

        if twodim:
            self.coef = ut_solv(self.time, self.ua[:, ind], self.va[:, ind],
                                self.lat[ind], **kwarg)

            self.QC.append('Harmonics done for velocity')

        else:
            self.coef = ut_solv(self.time, self.ua[:, ind], [],
                                self.lat[ind], **kwarg)

            self.QC.append('Harmonics done for elevation')
