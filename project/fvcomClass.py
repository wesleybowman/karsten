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

    ax can be defined as a region, i.e. a bounding box.
    An example:
        ax = [min(lon_coord), max(lon_coord), min(lat_coord), max(lat_coord)]
    '''

    def __init__(self, filename, ax=[]):
        self.QC = ['raw data']

        self.load(filename)

        if ax:
            self.ax = ax
        else:
            self.ax = [min(self.lon), max(self.lon), min(self.lat), max(self.lat)]


    def el_region(self):
        #self.lonc = np.intersect1d(np.argwhere(self.lonc>=self.ax[0]),np.argwhere(self.lonc<=self.ax[1]))
        #self.latc = np.intersect1d(np.argwhere(self.latc>=self.ax[2]),np.argwhere(self.latc<=self.ax[3]))
        #self.region_e = np.intersect1d(self.latc,self.lonc)
        #print self.region_e.shape

        self.region_e = np.argwhere((self.lonc >= self.ax[0]) &
                                    (self.lonc <= self.ax[1]) &
                                    (self.latc >= self.ax[2]) &
                                    (self.latc <= self.ax[3]))

        print self.region_e


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
        self.trinodes = self.data.variables['nv'][:]
        self.siglay = self.data.variables['siglay'][:]
        self.siglev = self.data.variables['siglev'][:]

        # Need to use len to get size of dimensions
        self.nele = self.data.dimensions['nele']
        self.node = self.data.dimensions['node']

        # elev timeseries
        self.el = self.data.variables['zeta']

        # get time and adjust it to matlab datenum
        self.time = self.data.variables['time'][:]
        self.time = self.time + 678942

        try:
            self.wa = self.data.variables['ww']
            self.u = self.data.variables['u']
            self.v = self.data.variables['v']
            self.ua = self.data.variables['ua']
            self.va = self.data.variables['va']
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

    def reconstr(self, time):
        if self.coef['aux']['opt']['twodim']:
            self.U, self.V = ut_reconstr(time, self.coef)
        else:
            self.ts_recon, _ = ut_reconstr(time, self.coef)





if __name__ == '__main__':
    filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
    test = FVCOM(filename)
    test.harmonics(0, cnstit='auto', notrend=True, nodiagn=True)
    test.reconstr(test.time)
