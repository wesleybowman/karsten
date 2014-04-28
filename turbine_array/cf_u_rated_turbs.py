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

    Turbine[0] = (area, Cp, rho, cutin, rated, cutoff, Prated, meanP, cf)

    return Turbine


def cf_u_rated_turbs(speed, turb, u_rated):

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
            P = calculate_power(u, u, turbs, i)
            turbs[i]['meanP '] = np.mean(P)
            turbs[i]['cf'] = turbs[i]['meanP'] / turbs[i]['Prated']

        a = np.where(turbs['cf'] > turb['cf'])[0]
        if a.shape[0] > 0:
            ur = np.max(turbs[a]['rated'])
            uri = np.argmax(turbs[a]['rated'])
            turbines[j] = turbs[a[uri]]
        else:
            turbines[j] = turbs[0]

    return turbines


def calculate_power(u, t, turbine, index, debug=False):

    u = np.abs(u)
    P = np.zeros(u.shape)
    nt = u.shape[0]

    # Error is in here
    ones = np.ones((nt, 1))
    turbine[index]['cutin'] = ones * turbine[index]['cutin']
    turbine[index]['rated'] = ones * turbine[index]['rated']
    turbine[index]['Prated'] = ones * turbine[index]['Prated']
    turbine[index]['cutoff'] = ones * turbine[index]['cutoff']

    print turbine[index]['Prated']

    P = ((u - turbine[index]['cutin']) / (turbine[index]['rated'] - turbine[index]['cutin']))**3 * \
        turbine[index]['Prated']

    P[u < turbine[index]['cutin']] = 0

    P[u > turbine[index]['rated']] = turbine[index]['Prated'][u > turbine[index]['rated']]

    P[u > turbine[index]['cutoff']] = 0

    if debug:

        for i in xrange(len(turbine)):
            plt.plot(t - t[i], P/1000, '-')
        plt.xlabel('Time (days)')
        plt.ylabel('Turbine Power Density (kW/m^2)')

    return P
