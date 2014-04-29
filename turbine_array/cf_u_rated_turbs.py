from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def createTurbines(size):

    Turbine = np.zeros(size, dtype={'names':
                                    ['area', 'Cp', 'rho', 'cutin',
                                     'rated', 'cutoff', 'Prated',
                                     'meanP', 'cf'],
                                    'formats':
                                    ['f8', 'f8', 'f8', 'f8', 'f8',
                                     'f8', 'f8', 'f8', 'f8']})

    Turbine = Turbine.view(np.recarray)

    return Turbine


def fillTurbines(Turbine, radius, Cp, rho, cutin, rated, cutoff, meanP, cf):

    radius = 10
    area = np.pi * radius**2
    Cp = 0.4
    rho = 1025
    cutin = 1
    rated = 2.5
    cutoff = 5
    Prated = Cp * (1/2 * rho * area * rated**3)
    meanP = 0
    cf = 0.40

    Turbine.area = area
    Turbine.Cp = Cp
    Turbine.rho = rho
    Turbine.cutin = cutin
    Turbine.rated = rated
    Turbine.cutoff = cutoff
    Turbine.Prated = Prated
    Turbine.meanP = meanP
    Turbine.cf = cf

    return Turbine


def cf_u_rated_turbs(speed, turb, u_rated):

    num_t, num_els = speed.shape

    turbines = createTurbines(num_els)

    for j in np.arange(num_els):
        print j
        u = speed[:, j]
#        max_u = np.max(u)

        turbs = createTurbines(u_rated.shape[0])

        turbs[:] = turb
        turbs.rated = u_rated
        turbs.Prated = turbs.Cp * (1/2 * turbs.rho * turbs.area *
                                   turbs.rated**3)

        P = calculate_power(u, u, turbs)
        turbs.meanP = np.mean(P, axis=1)
        turbs.cf = turbs.meanP / turbs.Prated

        a = np.where(turbs.cf > turb.cf)[0]
        if a.shape[0] > 0:
            uri = np.argmax(turbs.rated[a])
#            ur = np.max(turbs.rated[a])
            turbines[j] = turbs[a][uri]
        else:
            turbines[j] = turbs[0]

    return turbines


def calculate_power(u, t, turbine, debug=False):

    u = np.abs(u)
    nt = u.shape[0]

    ones = np.ones((nt))

    cutinTemp = ones * turbine.cutin[:, None]
    ratedTemp = ones * turbine.rated[:, None]
    pratedTemp = ones * turbine.Prated[:, None]
    cutoffTemp = ones * turbine.cutoff[:, None]

    P = ((u - cutinTemp) / (ratedTemp - cutinTemp))**3 * pratedTemp

    P[u < cutinTemp] = 0
    P[u > ratedTemp] = pratedTemp[u > ratedTemp]
    P[u > cutoffTemp] = 0

    if debug:

        for i in xrange(len(turbine)):
            plt.plot(t - t[i], P/1000, '-')
        plt.xlabel('Time (days)')
        plt.ylabel('Turbine Power Density (kW/m^2)')

    return P
