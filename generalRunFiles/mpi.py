from mpi4py import MPI
import numpy as np
import pandas as pd
import netCDF4 as nc
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv


def ncdatasort(x, y, time, trinodes, lon=None, lat=None):

    # Note, when time is passes, it is passed as
    # time*24*3600.
    # Does this need to be done this way?

    hour = 3600
    g = 9.806
    TP = 12.42
    rho = 1026
    period = (TP*3600)/(2*np.pi)

    # time=double(time)
    time = time+678942

    dt = time[1] - time[0]
    thour = time/hour
    deltat = thour[1]-thour[0]

    size = x.shape[0]
    nodexy = np.zeros((size, 2))
    nodexy[:, 0] = x
    nodexy[:, 1] = y
    nodexy = np.array((x, y)).T

    trinodes = trinodes.T-1
    triSize = trinodes.shape[0]
    uvnodexy = np.zeros((triSize, 2))

    uvnodexy[:, 0] = (nodexy[trinodes[:, 0], 0] + nodexy[trinodes[:, 1], 0] +
                      nodexy[trinodes[:, 2], 0]) / 3

    uvnodexy[:, 1] = (nodexy[trinodes[:, 0], 1] + nodexy[trinodes[:, 1], 1] +
                      nodexy[trinodes[:, 2], 1]) / 3

    if lat != None and lon != None:
        nodell = np.array((lon, lat)).T

        uvnodell = np.zeros((triSize, 2))

        uvnodell[:, 0] = (nodell[trinodes[:, 0], 0] +
                          nodell[trinodes[:, 1], 0] +
                          nodell[trinodes[:, 2], 0]) / 3

        uvnodell[:, 1] = (nodell[trinodes[:, 0], 1] +
                          nodell[trinodes[:, 1], 1] +
                          nodell[trinodes[:, 2], 1]) / 3

    else:
        'No nodell, uvnodell set to uvnodexy'
        uvnodell = uvnodexy

    return (nodexy, uvnodexy, dt, deltat, hour, thour,
            TP, rho, g, period, nodell, uvnodell, trinodes)

def mjd2num(x):

    y = x + 678942

    return y


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
#filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'
#filename = '/home/rkarsten/scratch/dn_coarse_2d_clean/output/dn_coarse_0001.nc'
filename = '/home/rkarsten/common_folder/dn_coarse_2d_clean/output/dn_coarse_0001.nc'

print 'Loading in Data'
data = nc.Dataset(filename, 'r')
x = data.variables['x'][:]
y = data.variables['y'][:]
xc = data.variables['xc'][:]
yc = data.variables['yc'][:]
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
lonc = data.variables['lonc'][:]
latc = data.variables['latc'][:]
h = data.variables['h'][:]
time = data.variables['time'][:]
trinodes = data.variables['nv'][:]
ua = data.variables['ua']
va = data.variables['va']
elev = data.variables['zeta']

print 'Data Loaded'

(nodexy, uvnodexy, dt, deltat,
hour, thour, TP, rho, g, period,
nodell, uvnodell, trinodes) = ncdatasort(x, y, time*24*3600,
                                        trinodes, lon, lat)

time = mjd2num(time)

print 'Done with ncdatasort'
Rayleigh = np.array([1])

lonlat = np.array([lon, lat]).T
lonclatc = np.array([lonc, latc]).T

comm.Barrier()

if rank == 0:
    print 'Creating NC layout'
    data = nc.Dataset('coef.nc', 'w', format='NETCDF4')
    data.createDimension('dim', None)
    data.createDimension('dimx', len(x))
    data.createDimension('dimtime', len(time))
    data.createDimension('dimtri', trinodes.shape[-1])

    newx = data.createVariable('x', 'f8', ('dimx',))
    newx[:] = x
    newy = data.createVariable('y', 'f8', ('dimx',))
    newy[:] = y
    newxc = data.createVariable('xc', 'f8', ('dim',))
    newxc[:] = xc
    newyc = data.createVariable('yc', 'f8', ('dim',))
    newyc[:] = yc
    newlon = data.createVariable('lon', 'f8', ('dimx',))
    newlon[:] = lon
    newlat = data.createVariable('lat', 'f8', ('dimx',))
    newlat[:] = lat
    newlonc = data.createVariable('lonc', 'f8', ('dim',))
    newlonc[:] = lonc
    newlatc = data.createVariable('latc', 'f8', ('dim',))
    newlatc[:] = latc
    newh = data.createVariable('h', 'f8', ('dimx',))
    newh[:] = h
    newtime = data.createVariable('time', 'f8', ('dimtime',))
    newtime[:] = time
    newtrinodes = data.createVariable('trinodes', 'f8', ('dim','dimtri'))
    newtrinodes[:] = trinodes

#    coef = ut_solv(time, ua[:, 0], va[:, 0], uvnodell[0, 1],
#                    'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
#                    'NoDiagn', 'LinCI')

    coef = ut_solv(time, ua[:, 0], va[:, 0], uvnodell[0, 1],
                    cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                    nodiagn=True, linci=True, conf_int=True)


    opt = pd.DataFrame(coef['aux']['opt'].items())
    del coef['aux']['opt']
    aux = pd.DataFrame(coef['aux'])
    del coef['aux']
    c = pd.DataFrame(coef)
    cat = pd.concat([c,aux], axis=1)


    data.createDimension('dim2', cat['Lsmaj'].shape[0])
    Lsmaj = data.createVariable('Lsmaj', 'f8', ('dim','dim2'))
    Lsmaj_ci = data.createVariable('Lsmaj_ci', 'f8', ('dim','dim2'))
    Lsmin = data.createVariable('Lsmin', 'f8', ('dim','dim2'))
    Lsmin_ci = data.createVariable('Lsmin_ci', 'f8', ('dim','dim2'))
    g = data.createVariable('g', 'f8', ('dim','dim2'))
    g_ci = data.createVariable('g_ci', 'f8', ('dim','dim2'))
    theta = data.createVariable('theta', 'f8', ('dim','dim2'))
    theta_ci = data.createVariable('theta_ci', 'f8', ('dim','dim2'))
    name = data.createVariable('name', 'c', ('dim','dim2'))
    A = data.createVariable('A', 'f8', ('dim','dim2'))
    A_ci = data.createVariable('A_ci', 'f8', ('dim','dim2'))
    gA = data.createVariable('gA', 'f8', ('dim','dim2'))
    gA_ci = data.createVariable('gA_ci', 'f8', ('dim','dim2'))
    nameA = data.createVariable('nameA', 'c', ('dim','dim2'))

    data.close()

else:
    data = None

comm.Barrier()

print 'Read in variables, starting calculations'
data = nc.Dataset('coef.nc', 'r', format='NETCDF4')
Lsmaj = data.variables['Lsmaj'][:]
Lsmaj_ci = data.variables['Lsmaj_ci'][:]
Lsmin = data.variables['Lsmin'][:]
Lsmin_ci = data.variables['Lsmin_ci'][:]
g = data.variables['g'][:]
g_ci = data.variables['g_ci'][:]
theta = data.variables['theta'][:]
theta_ci = data.variables['theta_ci'][:]
name = data.variables['name'][:]
A = data.variables['A'][:]
A_ci = data.variables['A_ci'][:]
gA = data.variables['gA'][:]
gA_ci = data.variables['gA_ci'][:]
nameA = data.variables['nameA'][:]

rows = [rank + size * i for i in range(int(len(lonc)/size)+1) if comm.rank + comm.size*i < len(lonc)]
#for i in xrange(rank, len(lonc), size):
for i in rows:
    print i
#    coef = ut_solv(time, ua[:, i], va[:, i], uvnodell[i, 1],
#                    'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
#                    'NoDiagn', 'LinCI')

    coef = ut_solv(time, ua[:, i], va[:, i], uvnodell[i, 1],
                    cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                    nodiagn=True, linci=True, conf_int=True)

    opt = pd.DataFrame(coef['aux']['opt'].items())
    del coef['aux']['opt']
    aux = pd.DataFrame(coef['aux'])
    del coef['aux']
    c = pd.DataFrame(coef)
    cat = pd.concat([c,aux], axis=1)

    Lsmaj[i,:] = cat['Lsmaj'].values
    Lsmin[i,:] = cat['Lsmin'].values
    g[i,:] = cat['g'].values
    theta[i, :] = cat['theta'].values
    name[i, :] = cat['name'].values
    Lsmaj_ci[i,:] = cat['Lsmaj_ci'].values
    Lsmin_ci[i,:] = cat['Lsmin_ci'].values
    theta_ci[i,:] = cat['theta_ci'].values
    g_ci[i,:] = cat['g_ci'].values

#    coefElev = ut_solv(time, elev[:, i], np.array([]), uvnodell[i, 1],
#                    'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
#                    'NoDiagn', 'LinCI')

    coefElev = ut_solv(time, ua[:, i], [], uvnodell[i, 1],
                    cnstit='auto', rmin=Rayleigh[0], notrend=True, method='ols',
                    nodiagn=True, linci=True, conf_int=True)


    opt = pd.DataFrame(coefElev['aux']['opt'].items())
    del coefElev['aux']['opt']
    aux = pd.DataFrame(coefElev['aux'])
    del coefElev['aux']
    c = pd.DataFrame(coefElev)
    cat = pd.concat([c,aux], axis=1)
    A[i, :] = cat['A'].values
    gA[i, :] = cat['g'].values
    nameA[i, :] = cat['name'].values
    A_ci[i, :] = cat['A_ci'].values
    gA_ci[i, :] = cat['g_ci'].values

comm.Barrier()
