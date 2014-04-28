from loadnc2d import loadnc2d
from cf_u_rated_turbs import *


filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'
(x, y, ua, va, trinodes, el, h, time, siglev, siglay, nbe, a1u, a2u,
            aw0, awx, awy, lon, lat, nele, node) = loadnc2d(filename)
speed = np.sqrt(ua*ua+va*va)
u = speed[:,0]
u_rated = np.arange(0, 8 + 0.25, 0.25)
turb = createTurbines(1)
turb = createTurbines(1)
turb = fillTurbines(turb, 10, 0.4, 1025, 1,2.5,5,0,0.4)

num_t, num_els = speed.shape

turbines = createTurbines(num_els)

for j in np.arange(num_els):
    u = speed[:, j]
    max_u = np.max(u)

    turbs = createTurbines(u_rated.shape[0])

    for i, v in enumerate(u_rated):

        turbs[i] = turb
        turbs[i]['rated'] = v
        turbs[i]['Prated'] = turbs[i]['Cp'] * (1/2 * turbs[i]['rho'] *
                                                turbs[i]['area'] *
                                                turbs[i]['rated']**3)
        print turbs[i]['Prated']
        P = calculate_power(u, u, turbs, i)
