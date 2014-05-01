from loadnc2d import loadnc2d
import numpy as np
from cf_u_rated_turbs import *


filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'

(x, y, ua, va, trinodes, el, h, time, siglev, siglay, nbe, a1u, a2u,
            aw0, awx, awy, lon, lat, nele, node) = loadnc2d(filename)

speed = np.sqrt(ua*ua+va*va)
u = speed[:,0]
u_rated = np.arange(0, 8 + 0.25, 0.25)
turb = createTurbines(1)
turb = fillTurbines(turb, 10, 0.4, 1025, 1,2.5,5,0,0.4)

turbines = cf_u_rated_turbs(speed, turb, u_rated)
