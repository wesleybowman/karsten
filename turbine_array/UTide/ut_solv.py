from __future__ import division
import numpy as np
import scipy.io as sio
import scipy.sparse
import scipy.signal
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

def ut_solv(tin, uin, vin, lat, cnstit, Rayleigh, *varargin):

    coef = ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin)

    return coef


def ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin):

    print 'ut_solv: '
    nt,t,u,v,tref,lor,elor,opt,tgd,uvgd = ut_slvinit(tin,uin,vin,cnstit,Rayleigh,varargin)

    opt['cnstit'] = cnstit
    [nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref, opt['rmin']/(24*lor),
                                           opt['cnstit'], opt['infer'])

    # a function we don't need
    # coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)

    coef['aux']['opt'] = opt
    coef['aux']['lat'] = lat

    print 'matrix prep ... '

    ngflgs = [opt['nodsatlint'], opt['nodsatnone'],
              opt['gwchlint'], opt['gwchnone']]

    #ngflgs = [opt.nodsatlint, opt.nodsatnone opt.gwchlint opt.gwchnone];

    E = ut_E(t,tref,cnstit['NR']['frq'],cnstit['NR']['lind'],lat,ngflgs,opt['prefilt'])

    B = np.hstack((E,E.conj()))

    # more infer stuff

    if opt['notrend']:
        B = np.hstack((B,np.ones((nt,1))))
        nm = 2 * (nNR + nR) + 1
    else:
        B = np.hstack((B, np.ones((nt,1)), (t-tref)/lor))
        nm = 2*(nNR + nR) + 2

    print 'Solution ...'

    xraw = u

    if opt['twodim']:
        #xraw = complex(u,v);
        xraw = u+v*1j

    if opt['method']=='ols':
        #m = B\xraw;
        m = np.linalg.lstsq(B, xraw)[0]
        #W = sparse(1:nt,1:nt,1);
        W = scipy.sparse.identity(nt)
#    else:
#        lastwarn('');
#        [m,solnstats] = robustfit(B,ctranspose(xraw),...
#            opt.method,opt.tunconst,'off');
#        if isequal(lastwarn,'Iteration limit reached.')
#            # nan-fill, create coef.results, reorder coef fields,
#            % do runtime display
#            coef = ut_finish(coef,nNR,nR,nI,elor,cnstit);
#            % abort remainder of calcs
#            return;
#        W = sparse(1:nt,1:nt,solnstats.w);

    xmod = B*m
    xmod = np.dot(B, m)

    if not opt['twodim']:
        xmod = np.real(xmod)

    e = W*(xraw-xmod)

    nc = nNR+nR
    ap = m[np.hstack((np.arange(nNR), 2*nNR+np.arange(nR)))]
    am = m[np.hstack((nNR+np.arange(nNR), 2*nNR+nR+np.arange(nR)))]

    Xu = np.real(ap + am)
    Yu = -np.imag(ap - am)

    if not opt['twodim']:
        #XY = np.hstack((Xu, Yu))
        coef['A'], _ , _ ,coef['g '] = ut_cs2cep(Xu, Yu)
        #coef['A'], _ , _ ,coef['g '] = ut_cs2cep(XY)

    else:
        Xv = np.imag(ap+am)
        Yv = np.real(ap-am)
        #XY = np.vstack((Xu, Yu, Xv, Yv))
        coef['Lsmaj'], coef['Lsmin'], coef['theta'], coef['g'] = ut_cs2cep(Xu, Yu, Xv, Yv)
        #coef['Lsmaj'], coef['Lsmin'], coef['theta'], coef['g '] = ut_cs2cep(XY)

    ## mean and trend
    if opt['twodim']:
        if opt['notrend']:
            coef['umean'] = np.real(m[-1])
            coef['vmean'] = np.imag(m[-1])
        else:
            coef['umean'] = np.real(m[-1-1])
            coef['vmean'] = np.imag(m[-1-1])
            coef['uslope'] = np.real(m[-1])/lor
            coef['vslope'] = np.imag(m[-1])/lor
    else:
        if opt['notrend']:
            coef['mean'] = np.real(m[-1])
        else:
            coef['mean'] = np.real(m[-1-1])
            coef['slope'] = np.real(m[-1])/lor

    # Confidence Intervals

    print 'conf. int''vls... '

    #import pdb; pdb.set_trace()

    if not opt['white']:
        # band-averaged (ba) spectral densities
        if opt['equi']:
            if np.sum(tgd) > np.sum(uvgd):
                efill = np.interp1(t,e,tin(tgd))
                if np.any(np.isnan(efill)): # fill start&/end nans w/ nearest good
                    ind = np.where(np.isnan(efill))[0]
                    #ind2 = ind(ind<find(~isnan(efill),1,'first'));
                    ind2 = ind[ind < np.where(~np.isnan(efill),1,'first')]
                    efill[ind2] = efill[np.max(ind2)+1]
                    ind2 = ind[ind>np.where(~np.isnan(efill),1,'last')]
                    efill[ind2] = efill[np.min(ind2)-1]

                ba = ut_pdgm(tin[tgd],efill,coef['aux']['frq'],1,0)
            else:
                ba = ut_pdgm(tin[tgd],e,coef['aux']['frq'],1,0)

        else:

            ba = ut_pdgm(t,e,coef.aux.frq,0,opt.lsfrqosmp);

        # power [ (e units)^2 ] from spectral density [ (e units)^2 / cph ]
        df = 1/(elor*24)
        ba['Puu']= ba['Puu']*df

        if opt['twodim']:
            ba['Pvv']= ba['Pvv']*df
            ba['Puv']= ba['Puv']*df

        # assign band-avg power values to NR & R freqs
        #Puu = np.zeros(size(coef.aux.frq));
        Puu = np.zeros(coef['aux']['frq'].shape)
        if opt['twodim']:
            Pvv = Puu
            Puv = Pvv

        #for i = 1:length(ba.Puu)
        for i in xrange(len(ba['Puu'])):
            #ind = find(coef.aux.frq>=ba.fbnd(i,1) & coef.aux.frq<=ba.fbnd(i,2))
            ind = np.where(np.logical_and(coef['aux']['frq']>=ba['fbnd'][i,0],
                           coef['aux']['frq']<=ba['fbnd'][i,0]))

            Puu[ind] = ba['Puu'][i]
            if opt['twodim']:
                Pvv[ind] = ba['Pvv'][i]
                Puv[ind] = ba['Puv'][i]



    #varMSM = real((ctranspose(xraw)*W*xraw - ctranspose(m)*ctranspose(B)*W*xraw)/(nt-nm))
#    varMSM = np.real((np.conj(xraw).T * W * xraw -
#                      np.conj(m).T[:,None] * np.conj(B).T * W * xraw)/(nt-nm))

    varMSM = np.real((np.dot(np.conj(xraw[:,None]).T * W, xraw[:,None]) -
                      np.dot(np.dot(np.conj(m[:,None]).T, np.conj(B).T) * W,
                             xraw[:,None]))/(nt-nm))

    #gamC = inv(ctranspose(B)*W*B)*varMSM
    gamC = np.linalg.inv(np.dot(np.conj(B).T * W, B)) * varMSM
    #gamP = inv(transpose(B)*W*B)*((transpose(xraw)*W*xraw - transpose(m)*transpose(B)*W*xraw)/(nt-nm))
#    gamP = np.dot(np.linalg.inv(np.dot(B.T * W, B)), ((xraw.T * W * xraw - m.T[:,None] * B.T * W *
#                                          xraw) / (nt-nm)))

    gamP = np.linalg.inv(np.dot(B.T * W, B)) * \
            (np.dot(xraw[:,None].T * W, xraw[:,None]) -
             np.dot(np.dot(m[:,None].T, B.T)*W, xraw[:,None]))/(nt-nm)


    Gall = gamC + gamP
    Hall = gamC - gamP

    coef['g_ci']= np.nan*np.ones(coef['g'].shape)
    if opt['twodim']:
        coef['Lsmaj_ci']= coef['g_ci']
        coef['Lsmin_ci']= coef['g_ci']
        coef['theta_ci']= coef['g_ci']
        varcov_mCw = np.nan*np.ones((nc,4,4))
    else:
        coef['A_ci']= coef['g_ci']
        varcov_mCw = np.nan*np.ones((nc,2,2))

    if not opt['white']:
        varcov_mCc = varcov_mCw

    #for c=1:nc
    for c in np.arange(nc):
        #G = [Gall(c,c) Gall(c,c+nc); Gall(c+nc,c) Gall(c+nc,c+nc);];
        G = np.array([[Gall[c,c], Gall[c,c+nc]],[Gall[c+nc,c], Gall[c+nc,c+nc]]])
        H = np.array([[Hall[c,c], Hall[c,c+nc]],[Hall[c+nc,c], Hall[c+nc,c+nc]]])
        #H = [Hall(c,c) Hall(c,c+nc); Hall(c+nc,c) Hall(c+nc,c+nc);];
        varXu = np.real(G[0,0]+G[1,1]+2*G[0,1])/2
        varYu = np.real(H[0,0]+H[1,1]-2*H[0,1])/2

        if opt['twodim']:
            varXv = np.real(H[0,0]+H[1,1]+2*H[0,1])/2
            varYv = np.real(G[0,0]+G[1,1]-2*G[0,1])/2
            #varXv = real(H(1,1)+H(2,2)+2*H(1,2))/2;
            #varYv = real(G(1,1)+G(2,2)-2*G(1,2))/2;

        '''Doesn't work cause Pvv1s is not working in ut_pdgm '''
        if opt['linci']: # linearized
            if not opt['twodim']:
                varcov_mCw[c,:,:] = np.diag(np.array([varXu, varYu]))
                if not opt['white']:
                    den = varXu + varYu
                    varXu = Puu[c]*varXu/den
                    varYu = Puu[c]*varYu/den
                    varcov_mCc[c,:,:] = np.diag(np.array([varXu, varYu]))
                sig1,sig2 = ut_linci(Xu[c],Yu[c],np.sqrt(varXu),np.sqrt(varYu))
                coef['A_ci'][c] = 1.96*sig1
                coef['g_ci'][c] = 1.96*sig2
            else:
                varcov_mCw[c,:,:] = np.diag(np.array([varXu, varYu, varXv, varYv]))
                if not opt['white']:
                    den = varXv + varYv
                    varXv = Pvv[c]*varXv/den
                    varYv = Pvv[c]*varYv/den
                    varcov_mCc[c,:,:] = np.diag(np.array([varXu, varYu, varXv,
                                                          varYv]))
                sig1,sig2 = ut_linci(Xu[c]+1j*Xv[c],Yu[c]+1j*Yv[c],
                                       np.sqrt(varXu)+1j*np.sqrt(varXv),
                                       np.sqrt(varYu)+1j*np.sqrt(varYv))
                coef['Lsmaj_ci'][c] = 1.96*np.real(sig1)
                coef['Lsmin_ci'][c] = 1.96*np.imag(sig1)
                coef['g_ci'][c] = 1.96*np.real(sig2)
                coef['theta_ci'][c] = 1.96*np.imag(sig2)
                #import pdb; pdb.set_trace()

        else: # monte carlo
            pass


        if opt['twodim']:
            PE = np.sum(coef['Lsmaj']**2+coef['Lsmin']**2)
            PE = 100* (coef['Lsmaj']**2+coef['Lsmin']**2)/PE
        else:
            PE = 100*coef['A']**2/np.sum(coef['A']**2)

        ind = PE.argsort()[::-1]
    if opt['twodim']:
        coef['Lsmaj'] = coef['Lsmaj'][ind]
        coef['Lsmaj_ci'] = coef['Lsmaj_ci'][ind]
        coef['Lsmin'] = coef['Lsmin'][ind]
        coef['Lsmin_ci'] = coef['Lsmin_ci'][ind]
        coef['theta'] = coef['theta'][ind]
        coef['theta_ci'] = coef['theta_ci'][ind]
    else:
        coef['A'] = coef['A'][ind]
        coef['A_ci'] = coef['A_ci'][ind]

    coef['g'] = coef['g'][ind]
    coef['name'] = coef['name'][ind]
    coef['aux']['frq'] = coef['aux']['frq'][ind]
    coef['aux']['lind'] = coef['aux']['lind'][ind]

    #import pdb; pdb.set_trace()

    return coef


def ut_pdgm(t,e,cfrq,equi,frqosmp):

    P = {}
    nt = len(e)
    hn = np.hanning(nt)

    if equi:
        # matlab pwelch
        # pwelch(x,window,noverlap,nfft)
        #[Puu1s,allfrq] = pwelch(real(e),hn,0,nt);
#        Puu1s, allfrq = scipy.signal.welch(np.real(e), window='hanning',
#                                           noverlap=0, nfft=nt, fs=2*np.pi)
        allfrq, Puu1s = scipy.signal.welch(np.real(e), window='hanning',
                                           noverlap=0, nfft=nt, fs=2*np.pi,
                                           detrend='constant',
                                           scaling='density')

#        allfrq, Puu1s = scipy.signal.periodogram(np.real(e), window='hanning',
#                                           nfft=nt, fs=2*np.pi,
#                                           detrend='constant',
#                                           scaling='density')
#        allfrq, Puu1s = scipy.signal.welch(np.real(e), window=hn, noverlap=0,
#                                           nfft=nt, fs=2*np.pi)
        #hn = mlab.window_hanning(t)
        Puu1s, allfrq = mlab.psd(np.real(e), window=hn, noverlap=0, NFFT=nt, Fs=2*np.pi)
    else:
        # ut_lmbscga
        pass

    #import pdb; pdb.set_trace()

    fac = (nt-1)/(2*np.pi*(t[-1]-t[0])*24) # conv fac: rad/sample to cph
    allfrq = allfrq*fac # to [cycle/hour] from [rad/samp]
    Puu1s = Puu1s / fac  # to [e units^2/cph] from [e units^2/(rad/samp)]

    #import pdb; pdb.set_trace()

    P['Puu'], P['fbnd'] = ut_fbndavg(Puu1s, allfrq, cfrq)

    if not np.isreal(e).all():

        if equi:
            #Pvv1s, _ = pwelch(np.imag(e),hn,0,nt)
            temp, Pvv1s = scipy.signal.welch(np.imag(e),window=hn, noverlap=0, nfft=nt, fs=2*np.pi)
            #temp, Pvv1s = scipy.signal.welch(np.imag(e),window=hn, noverlap=0, nfft=nt, fs=2*np.pi)

            # should be able to use mlab.csd
            #Puv1s, _ = cpsd(np.real(e),np.imag(e),hn,0,nt)
            #Pvv1s, temp = mlab.psd(np.imag(e), window=hn, noverlap=0, NFFT=nt, Fs=2*np.pi, sides='default')

            Puv1s, temp = mlab.csd(np.real(e),np.imag(e), noverlap=0, NFFT=nt, window=hn, Fs=2*np.pi)

        else:
            #Pvv1s, _ = ut_lmbscga(imag(e),t,hn,frqosmp);
            #Puv1s, _ = ut_lmbscgc(real(e),imag(e),t,hn,frqosmp);
            pass

        Pvv1s = Pvv1s/fac
        P['Pvv'], _ = ut_fbndavg(Pvv1s,allfrq,cfrq)
        Puv1s = np.real(Puv1s)/fac
        P['Puv'], _ = ut_fbndavg(Puv1s,allfrq,cfrq)
        P['Puv'] = np.abs(P['Puv'])


    return P

def ut_fbndavg(P,allfrq,cfrq):
    # UT_FBNDAVG()
    # line-decimate and band-average spectra
    # inputs
    # P = periodogram to treat [e units^2/cph]
    # allfrq = frequency values of (equispaced) P estimates [cph]
    # cfrq = frequencies of constituents [cph] (nc x 1)
    # outputs
    # avP = line-decimated and band-averaged spectrum [e units^2/cph] (9 x 1)
    # fbnd = frequencies [cph] at edges of averaged bands (9 x 2)
    # UTide v1p0 9/2011 d.codiga@gso.uri.edu
    # (based on residual_spectrum.m of t_tide, Pawlowicz et al 2002)

    df=allfrq[2]-allfrq[1]
    #P[np.round(cfrq/df).astype(int)+1] = np.nan
    np.round(cfrq/df).astype(int)
    P[np.round(cfrq/df).astype(int)] = np.nan

    fbnd =np.array([[.00010, .00417],
                    [.03192, .04859],
                    [.07218, .08884],
                    [.11243, .12910],
                    [.15269, .16936],
                    [.19295, .20961],
                    [.23320, .25100],
                    [.26000, .29000],
                    [.30000, .50000]])


    #nfbnd=size(fbnd,1);
    nfbnd=fbnd.shape[0]
    avP = np.zeros((nfbnd,1))

    #import pdb; pdb.set_trace()

    for k in np.arange(nfbnd-1,-1,-1):
    #for k=nfbnd:-1:1,
        b1 = np.where(allfrq>=fbnd[k,0])[0]
        b2 = np.where(allfrq<=fbnd[k,1])[0]
        b3 = np.where(np.isfinite(P))[0]
        jbnd = np.intersect1d(np.intersect1d(b1, b2), b3)

        #Issue is here, taking the mean of Pvv1s when Pvv1s is already off.
        if jbnd.any():
            #avP[k]=np.mean(P[jbnd-1])
            avP[k]=np.mean(P[jbnd])
        elif k < nfbnd:
            avP[k]=P[k+1]

    return avP, fbnd


#%---------------------------------------------------------
def ut_linci(X,Y,sigX,sigY):
    # UT_LINCI()
    # current ellipse parameter uncertainties from cosine/sine coefficient
    # uncertainties, by linearized relations w/ correlations presumed zero
    # inputs: (two-dim case complex, one-dim case real)
    # X = Xu + i*Xv
    # Y = Yu + i*Yv
    # for Xu =real(X) = u cosine coeff; Yu =real(Y) = u sine coeff
    # Xv =imag(X) = v cosine coeff; Yv =imag(Y) = v sine coeff
    # sigX = sigXu + i*sigXv
    # sigY = sigYu + i*sigYv
    # for sigXu =real(sigX) =stddev(Xu); sigYu =real(sigY) =stddev(Yu)
    # sigXv =imag(sigX) =stddev(Xv); sigYv =imag(sigY) =stddev(Yv)
    # outputs:
    # two-dim case, complex
    # sig1 = sig_Lsmaj +1i*sig_Lsmin [same units as inputs]
    # sig2 = sig_g + 1i*sig_theta [degrees]
    # one-dim case, real
    # sig1 = sig_A [same units as inputs]
    # sig2 = sig_g [degrees]
    # UTide v1p0 9/2011 d.codiga@gso.uri.edu
    # (adapted from errell.m of t_tide, Pawlowicz et al 2002)

    X = np.array([X])
    Y = np.array([Y])
    sigX = np.array([sigX])
    sigY = np.array([sigY])
    Xu = np.real(X[:])
    sigXu = np.real(sigX)
    Yu = np.real(Y[:])
    sigYu = np.real(sigY)

    Xv = np.imag(X[:])
    sigXv = np.imag(sigX[:])
    Yv = np.imag(Y[:])
    sigYv = np.imag(sigY[:])

    rp=.5*np.sqrt((Xu+Yv)**2+(Xv-Yu)**2)
    rm=.5*np.sqrt((Xu-Yv)**2+(Xv+Yu)**2)
    sigXu2=sigXu**2
    sigYu2=sigYu**2
    sigXv2=sigXv**2
    sigYv2=sigYv**2

    ex=(Xu+Yv)/rp
    fx=(Xu-Yv)/rm
    gx=(Yu-Xv)/rp
    hx=(Yu+Xv)/rm

    # major axis
    dXu2=(.25*(ex+fx))**2
    dYu2=(.25*(gx+hx))**2
    dXv2=(.25*(hx-gx))**2
    dYv2=(.25*(ex-fx))**2
    sig1 = np.sqrt(dXu2*sigXu2+dYu2*sigYu2+dXv2*sigXv2+dYv2*sigYv2)

    # phase
    rn=2*(Xu*Yu+Xv*Yv)
    rd=Xu**2-Yu**2+Xv**2-Yv**2
    den=rn**2+rd**2
    dXu2=((rd*Yu-rn*Xu)/den)**2
    dYu2=((rd*Xu+rn*Yu)/den)**2
    dXv2=((rd*Yv-rn*Xv)/den)**2
    dYv2=((rd*Xv+rn*Yv)/den)**2
    sig2 = (180/np.pi)*np.sqrt(dXu2*sigXu2+dYu2*sigYu2+dXv2*sigXv2+dYv2*sigYv2)

    #if ~isreal(X)
    if not np.isreal(X):
        # minor axis
        dXu2=(.25*(ex-fx))**2
        dYu2=(.25*(gx-hx))**2
        dXv2=(.25*(hx+gx))**2
        dYv2=(.25*(ex+fx))**2
        sig1 = sig1 + 1j*np.sqrt(dXu2*sigXu2+dYu2*sigYu2+dXv2*sigXv2+dYv2*sigYv2)

        # orientation
        rn=2.*(Xu*Xv+Yu*Yv)
        rd=Xu**2+Yu**2-(Xv**2+Yv**2)
        den=rn**2+rd**2
        dXu2=((rd*Xv-rn*Xu)/den)**2
        dYu2=((rd*Yv-rn*Yu)/den)**2
        dXv2=((rd*Xu+rn*Xv)/den)**2
        dYv2=((rd*Yu+rn*Yv)/den)**2
        sig2 = sig2 + 1j*(180/np.pi)*np.sqrt(dXu2*sigXu2+dYu2*sigYu2 + dXv2*sigXv2+dYv2*sigYv2)

    return sig1, sig2




def ut_cs2cep(Xu, Yu, Xv=np.array([False]), Yv=np.array([False])):

    #Xu = XY[:, 0]
    #Yu = XY[:, 1]

    if not Xv.all():
        Xv = np.zeros(Xu.shape)
        Yv = np.zeros(Yu.shape)

#    if XY.shape[-1] > 2:
#        Xv = XY[:, 3]
#        Yv = XY[:, 4]
#    else:
#        Xv = np.zeros(Xu.shape)
#        Yv = np.zeros(Yu.shape)

    ap = ((Xu+Yv)+1j*(Xv-Yu))/2
    am = ((Xu-Yv)+1j*(Xv+Yu))/2
    Ap = np.abs(ap)
    Am = np.abs(am)
    Lsmaj = Ap+Am
    Lsmin = Ap-Am
    epsp = np.angle(ap)*180/np.pi
    epsm = np.angle(am)*180/np.pi

    theta = ((epsp+epsm)/2) % 180
    g = (-epsp+theta) % 360

    return Lsmaj,Lsmin,theta,g


def ut_E(t,tref,frq,lind,lat,ngflgs,prefilt):

    nt = len(t)
    nc = len(lind)
    if ngflgs[1] and ngflgs[3]:
        F = np.ones((nt,nc))
        U = np.zeros((nt,nc))
        V = 24*(t-tref)*frq
    else:
        F, U, V = ut_FUV(t,tref,lind,lat,ngflgs);

    E = F * np.exp(1j*(U+V)*2*np.pi)

    #if ~isempty(prefilt)
#    if len(prefilt)!=0:
#        P=interp1(prefilt.frq,prefilt.P,frq).T
#        P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
#        E = E*P(ones(nt,1),:);

    return E


def ut_FUV(t, tref, lind, lat, ngflgs):

    nt = len(t)
    nc = len(lind)
    ## nodsat

    if ngflgs[1]:
        F = np.ones((nt,nc))
        U = np.zeros((nt,nc))
    else:
        if ngflgs[0]:
            tt = tref
        else:
            tt = t

        ntt = len(tt)

        mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)
        sat = mat_contents['sat']
        const = mat_contents['const']
        shallow = mat_contents['shallow']

        astro, ader = ut_astron(tt)

        if abs(lat) < 5:
            lat=np.sign(lat)*5;

        slat = np.sin(np.pi * lat/180)
        rr = sat.amprat
        j = np.where(sat.ilatfac==1)[0]

        rr[j] = rr[j]*0.36309*(1.0-5.0*slat*slat)/slat;

        j = np.where(sat.ilatfac==2)

        rr[j]=rr[j]*2.59808*slat;

        uu = np.dot(sat.deldood, astro[3:6, :]) + sat.phcorr[:,None]*np.ones((1,ntt)) % 1

        nfreq=len(const.isat)
        mat = rr[:,None]*np.ones((1,ntt)) * np.exp(1j*2*np.pi*uu)

        F = np.ones((nfreq, ntt)) + 0j
        ind = np.unique(sat.iconst)

        for i in xrange(len(ind)):
            F[ind[i]-1, :] = 1+np.sum(mat[sat.iconst==ind[i],:], axis=0)

        #U = imag(log(F))/(2*pi); % faster than angle(F)
        U = np.imag(np.log(F)) / (2*np.pi)
        F = np.abs(F)

        for k in np.where(np.isfinite(const.ishallow))[0]:
            ik=const.ishallow[k]+np.arange(const.nshallow[k])
            ik = ik.astype(int)
            j = shallow.iname[ik-1]
            exp1 = shallow.coef[ik-1]
            exp2 = np.abs(exp1)
            temp1 = exp1*np.ones((ntt,1))
            temp2 = exp2*np.ones((ntt,1))
            temp1 = temp1.T
            temp2 = temp2.T
            F[k,:]=np.prod(F[j-1,:]**temp2,axis=0)
            U[k,:]=np.sum(U[j-1,:]*temp1,axis=0)

        F=F[lind,:].T
        U=U[lind,:].T

        if ngflgs[1]: # nodal/satellite with linearized times
            F = F[np.ones((nt,1)),:]
            U = U[np.ones((nt,1)),:]

    ## gwch (astron arg)
    if ngflgs[3]: # none (raw phase lags not greenwich phase lags)
#        if ~exist('const','var'):
#            load('ut_constants.mat','const');
#        [~,ader] = ut_astron(tref);
#        ii=isfinite(const.ishallow);
#        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
#        for k=find(ii)'
#            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
#            const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
        V = 24*(t-tref)*const.freq(lind).T
    else:
        if ngflgs[3]: # linearized times
            tt = tref
        else:
            tt = t # exact times

        ntt = len(tt)
#        if exist('astro','var')
#            if ~isequal(size(astro,2),ntt)
#                [astro,~]=ut_astron(tt');
#            end
#        else
#            [astro,~]=ut_astron(tt');
#
#        if ~exist('const','var')
#            load('ut_constants.mat')

        mat_contents = sio.loadmat('ut_constants.mat',
                                   struct_as_record=False, squeeze_me=True)
        sat = mat_contents['sat']
        const = mat_contents['const']
        shallow = mat_contents['shallow']
        astro, ader = ut_astron(tt)

        #V = np.dot(const.doodson, astro) + const.semi[:,None]*np.ones((1,ntt)) % 1
        V = np.dot(const.doodson, astro) + const.semi[:,None]*np.ones((1,ntt))
        #V = V % 1

        for k in np.where(np.isfinite(const.ishallow))[0]:
            ik=const.ishallow[k]+np.arange(const.nshallow[k])
            ik = ik.astype(int)
            j = shallow.iname[ik-1]
            exp1 = shallow.coef[ik-1]
            temp1 = exp1[:]*np.ones((ntt,1))
            temp1 = temp1.T

            V[k,:] = np.sum(V[j-1,:]*temp1,axis=0)

        V=V[lind,:].T

#        if ngflgs(3) % linearized times
#            [~,ader] = ut_astron(tref);
#            ii=isfinite(const.ishallow);
#            const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
#            for k=find(ii)'
#                ik=const.ishallow(k)+(0:const.nshallow(k)-1);
#                const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
#                    shallow.coef(ik) );
#            end
#            V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
#        end

    return F, U, V


def ut_cnstitsel(tref,minres,incnstit,infer):

    mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)
    shallow = mat_contents['shallow']
    const = mat_contents['const']

    cnstit = {}
    coef = {}

    astro, ader = ut_astron(tref)

    ii = np.isfinite(const.ishallow)
    const.freq[~ii] = np.dot(const.doodson[~ii, :], ader) / 24

    for k in ii.nonzero()[0]:
        ik = const.ishallow[k]+np.arange(const.nshallow[k])
        ik = ik.astype(int)-1
        const.freq[k] = np.sum(const.freq[shallow.iname[ik]-1]*shallow.coef[ik])

    ## cnstit.NR
    cnstit['NR'] = {}
    if incnstit.lower()=='auto':
        cnstit['NR']['lind'] = np.where(const.df >= minres)[0]
    else:
        pass

    # skipped some stuff here cause they involve infer

    #import pdb; pdb.set_trace()
    cnstit['NR']['frq'] = const.freq[cnstit['NR']['lind']]
    cnstit['NR']['name'] = const.name[cnstit['NR']['lind']]
    nNR = len(cnstit['NR']['frq'])

    ## cnstit.R
    nR = 0
    nI = 0
    cnstit['R'] = []

    nallc = nNR+nR+nI

    coef['name'] = cnstit['NR']['name']
    coef['aux'] = {}
    coef['aux']['frq'] = cnstit['NR']['frq']
    coef['aux']['lind'] = cnstit['NR']['lind']

    # another infer if statement

    coef['aux']['reftime'] = tref


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
    d = jd[:] - daten
    D = d / 10000

    #args = np.array([[np.ones(jd.shape)],[d],[D*D],[D**3]]).flatten()[:,None]
    args = np.vstack((np.ones(jd.shape), d, D*D, D**3))
    sc = np.array([270.434164, 13.1763965268, -0.0000850, 0.000000039])
    hc= np.array([ 279.696678, 0.9856473354, 0.00002267,0.000000000])
    pc= np.array([ 334.329556, 0.1114040803,-0.0007739,-0.00000026])
    npc= np.array([-259.183275, 0.0529539222,-0.0001557,-0.000000050])
    ppc= np.array([281.220844, 0.0000470684, 0.0000339, 0.000000070])

    astro = np.dot(np.vstack((sc, hc, pc, npc, ppc)), args) / 360 % 1
    tau = jd%1 + astro[1,:] - astro[0,:]
    astro = np.vstack((tau,astro))

#    dargs = np.array([[np.zeros(jd.shape[0])], [np.ones(jd.shape[0])],
#                      [2.0e-4*D], [3.0e-4*D*D]]).flatten()[:,None]

    dargs = np.vstack((np.zeros(jd.shape), np.ones(jd.shape),
                      2.0e-4*D, 3.0e-4*D*D))

    ader = np.dot(np.vstack((sc,hc,pc,npc,ppc)), dargs)/360.0
    dtau = 1.0 + ader[1,:] - ader[0, :]
    ader = np.vstack((dtau, ader))

    # might need to take out depending on shape of jd
    #astro = astro.flatten()
    #ader = ader.flatten()

    return astro,ader


def loadMAT(filename):

    mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)

    items = []
    items = {}

    for i in mat_contents:
        name = '{0}'.format(i)
        items[name] = mat_contents[name]
        #i = mat_contents[name]
        #items.append(i)

    return items
