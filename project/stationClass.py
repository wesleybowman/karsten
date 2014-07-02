import netCDF4 as nc


class station:
    def __init__(self, filename):
        self.load(filename)

    def load(self, filename):
        self.data = nc.Dataset('filename')
        self.x = self.data.variables['x'][:]
        self.y = self.data.variables['y'][:]
        self.lon = self.data.variables['lon'][:]
        self.lat = self.data.variables['lat'][:]
        self.siglay = self.data.variables['siglay'][:]
        self.siglev = self.data.variables['siglev'][:]
        self.h = self.data.variables['h'][:]
        self.time_JD = self.data.variables['time_JD'][:]
        self.time_second = self.data.variables['time_second'][:]
        self.u = self.data.variables['u']
        self.v = self.data.variables['v']
        self.ww = self.data.variables['ww']
        self.ua = self.data.variables['ua']
        self.va = self.data.variables['va']
        self.elev = self.data.variables['elev']

if __name__ == '__main__':
    filename = '/array2/data3/rkarsten/dncoarse_3D/output/dn_coarse_station_timeseries.nc'
