import netCDF4 as nc


def loadnc2d(filename):

    data = nc.Dataset(filename, 'r')
    x = data.variables['x'][:]
    y = data.variables['y'][:]
    lon = data.variables['lon'][:]
    lat = data.variables['lat'][:]
    ua = data.variables['ua'][:]
    va = data.variables['va'][:]
    el = data.variables['zeta'][:]
    h = data.variables['h'][:]
    nbe = data.variables['nbe'][:]
    a1u = data.variables['a1u'][:]
    a2u = data.variables['a2u'][:]
    aw0 = data.variables['aw0'][:]
    awx = data.variables['awx'][:]
    awy = data.variables['awy'][:]
    time = data.variables['time'][:]
    trinodes = data.variables['nv'][:]
    siglay = data.variables['siglay'][:]
    siglev = data.variables['siglev'][:]
    nele = data.variables['nele'][:]
    node = data.variables['node'][:]

    return (x, y, ua, va, trinodes, el, h, time, siglev, siglay, nbe, a1u, a2u,
            nele, aw0, awx, awy, node, lon, lat)
