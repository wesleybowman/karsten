from __future__ import division
import numpy as np


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
            TP, rho, g, period, nodell, uvnodell, time, trinodes)
