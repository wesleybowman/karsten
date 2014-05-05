import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
#from mpl_toolkits.basemap import Basemap
import numpy as np


filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'
data = nc.Dataset(filename,'r')
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
nv = data.variables['nv'][:].T -1
h = data.variables['h'][:]
el = data.variables['zeta'][:]


tri = Tri.Triangulation(lon,lat,triangles=nv)

levels = np.arange(-38, -4, 1)   # depth contours to plot

fig = plt.figure(figsize=(18,10))
plt.rc('font',size='22')
ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
plt.tricontourf(tri, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
plt.tricontourf(tri, -h,shading='faceted',cmap=plt.cm.gist_earth)
#plt.tripcolor(tri, facecolors = h)
#plt.tripcolor(lon,lat,triangles=nv, facecolors=h)

plt.triplot(tri)
plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.gca().patch.set_facecolor('0.5')
cbar=plt.colorbar()
cbar.set_label('Water Depth (m)', rotation=-90,labelpad=30)

scale = 1
ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
ax.xaxis.set_major_formatter(ticks)
ax.yaxis.set_major_formatter(ticks)
plt.grid()

plt.show()



#mmap = Basemap(projection='merc', lat_0 = 44.5, lon_0 = -66.4,
#    resolution = 'h', area_thresh = 0.1,
#    llcrnrlon=-71.46, llcrnrlat=37.6,
#    urcrnrlon=-57.66, urcrnrlat=45.97)
#
#mmap.drawcoastlines()
#mmap.drawcountries()
#mmap.drawmapboundary()
#mmap.fillcontinents(color='coral')
#
#plt.show()
