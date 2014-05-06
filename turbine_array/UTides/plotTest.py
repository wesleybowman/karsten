import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as ticker
#from mpl_toolkits.basemap import Basemap
import numpy as np
import cPickle as pickle
import seaborn
import time

filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'
data = nc.Dataset(filename,'r')
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
nv = data.variables['nv'][:].T -1
h = data.variables['h'][:]
el = data.variables['zeta'][:]
x = data.variables['x'][:]
y = data.variables['y'][:]

trinodes = nv
xc = np.mean(x[trinodes], axis=1)
yc = np.mean(y[trinodes], axis=1)
hc = np.mean(h[trinodes], axis=1)

lonc = np.mean(lon[trinodes], axis=1)
latc = np.mean(lat[trinodes], axis=1)

loci = pickle.load(open('loci.p', 'rb'))
loci = loci.astype(int)

latind = np.argwhere(((45.2<lat[:]), (lat<45.4)))
lonind = np.argwhere(((-64.8<lon), (lon<-64.1)))
lat[latind]
lon[lonind]


#tri = Tri.Triangulation(lon,lat,triangles=nv)
tri = Tri.Triangulation(lon[lonind], lat[latind], triangles=nv)

levels = np.arange(-100, -8, 1)

fig = plt.figure(figsize=(18,10))
#plt.ion()
plt.rc('font',size='22')

#ax = fig.add_subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
ax = fig.add_subplot(111)

#plt.tricontourf(tri, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
plt.tricontourf(tri, -h,levels=levels,shading='faceted')

plt.triplot(tri)
plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.gca().patch.set_facecolor('0.5')
cbar = plt.colorbar()
cbar.set_label('Water Depth (m)', rotation=-90, labelpad=30)

scale = 1
ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
ax.xaxis.set_major_formatter(ticks)
ax.yaxis.set_major_formatter(ticks)
plt.grid()

#plt.plot(xc[loci], yc[loci], 'ko')
#plt.plot(lonc[loci], latc[loci], 'ko')
#plt.plot(lonc[loci[0]], latc[loci[0]], 'ko')
plt.show()


for i,v in enumerate(loci):
    print i
    plt.tricontourf(tri, -h,levels=levels,shading='faceted')
    plt.triplot(tri)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.gca().patch.set_facecolor('0.5')
    cbar = plt.colorbar()
    cbar.set_label('Water Depth (m)', rotation=-90, labelpad=30)

    scale = 1
    ticks = ticker.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
    ax.xaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_formatter(ticks)
    plt.grid()

    plt.plot(lonc[loci[0:i]], latc[loci[0:i]], 'ko')
    #time.sleep(0.5)
    #plt.draw()
    plt.show()


#plt.show()

#fig=plt.figure()
#plt.axis([0,1000,0,1])
#
#i=0
#x=list()
#y=list()
#
#while i <1000:
#    temp_y=np.random.random()
#    x.append(i)
#    y.append(temp_y)
#    plt.scatter(i,temp_y)
#    i+=1
#    plt.show()
