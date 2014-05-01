import numpy as np

def ut_solv(tin, uin, vin, lat, cnstit, Rayleigh, *varargin):

    print varargin
    nt, t, u, v, tref, lor, elor, opt, tgd, uvgd = ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin)
    #coef = ut_solv1(tin,uin,vin,lat,cnstit,varargin)
    return nt, t, u, v, tref, lor, elor, opt, tgd, uvgd


def ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin):

    print 'ut_solv: '
    nt,t,u,v,tref,lor,elor,opt,tgd,uvgd = ut_slvinit(tin,uin,vin,cnstit,Rayleigh,varargin)

#    opt['cnstit'] = cnstit
#    [nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref, opt['rmin']/(24*lor),
#                                           opt['cnstit'], opt['infer'])
#
#    coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)
#    coef.aux.opt = opt;
#    coef.aux.lat = lat;





    return nt, t, u, v, tref, lor, elor, opt, tgd, uvgd


def ut_cnstitsel(tref,minres,incnstit,infer):

    mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)
    shallow = mat_contents['shallow']
    const = mat_contents['const']

    astro, ader = ut_astron(tref)

    ii = np.isfinite(const.ishallow)
    # might need if didn't use flatten
    const.freq[~ii] = np.dot(const.doodson[~ii, :], ader) / 24

    #k = ii.nonzero()[0]
    for k in ii.nonzero()[0]:
        ik = const.ishallow[k]+np.arange(const.nshallow[k])
        ik = ik.astype(int)-1
        const.freq[k] = np.sum(const.freq[shallow.iname[ik]]*shallow.coef[ik])

    return nNR, nR, nI, cnstit, coef


def ut_slvinit(tin,uin,vin,cnstit,Rayleigh,args):

    opt = {}
    args = list(args)
    tgd = ~np.isnan(tin)
    uin = uin[tgd]
    tin = tin[tgd]

    if vin.shape[0] == 0:
        opt['twodim'] = False
        #twodim = False
        v = np.array([])
    else:
        opt['twodim'] = True
        #twodim = True
        vin = vin[tgd]

    #if twodim:
    if opt['twodim']:
        uvgd = ~np.isnan(uin) & ~np.isnan(vin)
        v = vin[uvgd]
    else:
        uvgd = ~np.isnan(uin)

    t = tin[uvgd]
    nt = len(t)
    u = uin[uvgd]
    eps = np.finfo(np.float64).eps

    if np.var(np.unique(np.diff(tin))) < eps:
        opt['equi'] = 1 # based on times; u/v can still have nans ("gappy")
        #equi = 1 # based on times; u/v can still have nans ("gappy")
        lor = (np.max(tin)-np.min(tin))
        elor = lor*len(tin)/(len(tin)-1)
        tref = 0.5*(tin[0]+tin[-1])
    else:
        opt['equi'] = 0
        #equi = 0;
        lor = (np.max(t) - np.min(t))
        elor = lor*nt/(nt-1)
        tref = 0.5*(t[0]+t[-1])

    ## options
    opt['notrend'] = 0
    opt['prefilt'] = []
    opt['nodsatlint'] = 0
    opt['nodsatnone'] = 0
    opt['gwchlint'] = 0
    opt['gwchnone'] = 0
    opt['infer'] = []
    opt['inferaprx'] = 0
    opt['rmin'] = 1
    opt['method'] = 'cauchy'
    opt['tunrdn'] = 1
    opt['linci'] = 0
    opt['white'] = 0
    opt['nrlzn'] = 200
    opt['lsfrqosmp'] = 1
    opt['nodiagn'] = 0
    opt['diagnplots'] = 0
    opt['diagnminsnr'] = 2
    opt['ordercnstit'] = []
    opt['runtimedisp'] = 'yyy'

    #methnotset = 1
    allmethods = ['ols', 'andrews', 'bisquare', 'fair', 'huber',
                  'logistic', 'talwar', 'welsch']

    args = [string.lower() for string in args]

    if 'notrend' in args:
        #opt['notrend'] = 1
        opt['notrend'] = True

    if 'rmin' in args:
        opt['rmin'] = Rayleigh

    if 'nodiagn' in args:
        #opt['nodiagn']=1
        opt['nodiagn'] = True

    if 'linci' in args:
        #opt['linci'] = 1
        opt['linci'] = True

    if allmethods:
        methods = [i for i in allmethods if i in args]
        if len(methods) > 1:
            print 'ut_solv: Only one "method" option allowed.'
        else:
            opt['method'] = methods[0]

    if opt['method'] != 'cauchy':
        ind = np.argwhere(opt['method'] in allmethods)[0][0]
        allconst = [np.nan, 1.339, 4.685, 1.400, 1.345, 1.205, 2.795, 2.985]
        opt['tunconst'] = allconst[ind]
    else:
        opt['tunconst'] = 2.385

    opt['tunconst'] = opt['tunconst'] /opt['tunrdn']

# only needed if we sort the options
#    nf = len(opt)

    return nt, t, u, v, tref, lor, elor, opt, tgd, uvgd


def ut_astron(jd):
    '''
    UT_ASTRON()
    calculate astronomical constants
    input
    jd = time [datenum UTC] (1 x nt)
    outputs
    astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
    ader = matrix of derivatives of astro [cycles/day] (6 x nt)
    UTide v1p0 9/2011 d.codiga@gso.uri.edu
    (copy of t_astron.m from t_tide, Pawlowicz et al 2002)
    '''

    jd = np.array([jd])
    # datenum(1899,12,31,12,0,0)
    daten = 693961.500000000
    d = jd - daten
    D = d / 10000

    args = np.array([[np.ones(jd.shape[0])],[d],[D*D],[D**3]])
    sc = np.array([270.434164, 13.1763965268, -0.0000850, 0.000000039])
    hc= np.array([ 279.696678, 0.9856473354, 0.00002267,0.000000000])
    pc= np.array([ 334.329556, 0.1114040803,-0.0007739,-0.00000026])
    npc= np.array([-259.183275, 0.0529539222,-0.0001557,-0.000000050])
    ppc= np.array([281.220844, 0.0000470684, 0.0000339, 0.000000070])

    astro = np.dot(np.vstack((sc, hc, pc, npc, ppc)), args) / 360 % 1
    tau = jd%1 + astro[1,:] - astro[0,:]
    astro = np.vstack((tau,astro))
    dargs = np.array([[np.zeros(jd.shape[0])], [np.ones(jd.shape[0])],
                      [2.0e-4*D], [3.0e-4*D*D]])

    ader = np.dot(np.vstack((sc,hc,pc,npc,ppc)), dargs)/360.0
    dtau = 1.0 + ader[1,:] - ader[0, :]
    ader = np.vstack((dtau, ader))

    # might need to take out depending on shape of jd
    astro = astro.flatten()
    ader = ader.flatten()

    return astro,ader
