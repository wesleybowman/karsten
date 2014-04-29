import numpy as np
import matplotlib.pyplot as plt


def createTurbines(size):

#    Turbine = np.zeros(size, dtype= [('area',float),
#                                     ('Cp',float),
#                                     ('rho',float),
#                                     ('cutin',float)
#                                     ('rated',float),
#                                     ('cutoff',float),
#                                     ('Prated',float)
#                                     ('meanP',float),
#                                     ('cf',float)])

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

    #Turbine[0] = (area, Cp, rho, cutin, rated, cutoff, Prated, meanP, cf)
    Turbine.area = np.pi*radius**2
    Turbine.Cp=0.4
    Turbine.rho=1025
    Turbine.cutin=1
    Turbine.rated=2.5
    Turbine.cutoff=5
    Turbine.Prated=Turbine.Cp*(1/2*Turbine.rho*Turbine.area*Turbine.rated**3)
    Turbine.meanP=0
    Turbine.cf=0.40


    return Turbine


def cf_u_rated_turbs(speed, turb, u_rated):

    num_t, num_els = speed.shape

    turbines = createTurbines(num_els)

    for j in np.arange(num_els):
        print j
        u = speed[:, j]
        max_u = np.max(u)

        turbs = createTurbines(u_rated.shape[0])

        for i, v in enumerate(u_rated):

            turbs[i] = turb
            turbs.rated[i] = v
            turbs.Prated[i] = turbs.Cp[i] * (1/2 * turbs.rho[i] *
                                                   turbs.area[i] *
                                                   turbs.rated[i]**3)

#            turbs.Prated[i] = turbs[i].Cp * (1/2 * turbs[i].rho *
#                                                   turbs[i].area *
#                                                   turbs[i].rated**3)
            P = calculate_power(u, u, turbs, i)
            turbs.meanP[i] = np.mean(P)
            turbs.cf[i] = turbs.meanP[i] / turbs.Prated[i]

        a = np.where(turbs.cf > turb.cf)[0]
        if a.shape[0] > 0:
            ur = np.max(turbs[a]['rated'])
            uri = np.argmax(turbs[a]['rated'])
            turbines[j] = turbs[a[uri]]
        else:
            turbines[j] = turbs[0]

    return turbines


def calculate_power(u, t, turbine, index, debug=False):

    u = np.abs(u)
    P = np.zeros(u.shape[0])
    nt = u.shape[0]

    # Error is in here
    #ones = np.ones((nt, 1))
    ones = np.ones((nt))
#    turbine[index]['cutin'] = np.tile(turbine[index]['cutin'],(nt,1))
#    turbine[index]['rated'] = np.tile(turbine[index]['rated'],(nt,1))
#    turbine[index]['Prated'] = np.tile(turbine[index]['Prated'],(nt,1))
#    turbine[index]['cutoff'] = np.tile(turbine[index]['cutoff'],(nt,1))

#    turbine.cutin = np.tile(turbine.cutin,(nt,1))
#    turbine.rated = np.tile(turbine.rated,(nt,1))
#    turbine.Prated = np.tile(turbine.Prated,(nt,1))
#    turbine.cutoff = np.tile(turbine.cutoff,(nt,1))

#    turbine[index]['cutin'] = ones * turbine[index]['cutin']
#    turbine[index]['rated'] = ones * turbine[index]['rated']
#    turbine[index]['Prated'] = ones * turbine[index]['Prated']
#    turbine[index]['cutoff'] = ones * turbine[index]['cutoff']

    cutinTemp = ones * turbine[index]['cutin']
    ratedTemp = ones * turbine[index]['rated']
    pratedTemp = ones * turbine[index]['Prated']
    cutoffTemp = ones * turbine[index]['cutoff']

    #print turbine.Prated
    P = ((u - cutinTemp) / (ratedTemp- cutinTemp))**3 * pratedTemp

#    P = ((u - turbine[index]['cutin']) / (turbine[index]['rated']- turbine[index]['cutin']))**3 * \
#        turbine[index]['Prated']

#    P = ((u - turbine.cutin) / (turbine.rated - turbine.cutin))**3 * \
#        turbine.Prated

#    P[u < turbine.cutin] = 0
#    P[u > turbine.rated] = turbine.Prated[u > turbine.rated]
#    P[u > turbine.cutoff] = 0
    #temp = pratedTemp.flatten()
    #temp2 = ratedTemp.flatten()

    P[u < cutinTemp] = 0
    P[u > ratedTemp] = pratedTemp[u>ratedTemp]
    P[u > cutoffTemp] = 0

#    P[u < turbine[index]['cutin']] = 0
#    P[u > turbine[index]['rated']] = turbine[index]['Prated'][u > turbine[index]['rated']]
#    P[u > turbine[index]['cutoff']] = 0

    if debug:

        for i in xrange(len(turbine)):
            plt.plot(t - t[i], P/1000, '-')
        plt.xlabel('Time (days)')
        plt.ylabel('Turbine Power Density (kW/m^2)')

    return P
