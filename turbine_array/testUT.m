[x,y,ua,va,trinodes,el,h,time,siglev,siglay,nbe,a1u,a2u,aw0,awx,awy,nele,node,long,lat] = loadnc2d_matlab_RK('/home/wesley/github/aidan-projects/grid','/dngrid_0001.nc');
%[x,y,ua,va,trinodes,el,h,time,siglev,siglay,nbe,a1u,a2u,aw0,awx,awy,nele,node,long,lat] = loadnc2d_matlab_RK('/home/wesley/ncfiles','/smallcape_force_0001.nc');

[nodexy,uvnode,dt,deltat,hour,thour,TP,rho,g,period,nodell,uvllnode]=ncdatasort(x,y,time*24*3600,trinodes,long,lat);
time=mjd2num(time);

%power station location
ps_loc=[-64.4031,45.3710];
ps_loc(2,:)=[-64.4031,45.3710];
psi=closest_point(ps_loc,[long,lat]);
xc=mean(x(trinodes),2);
yc=mean(y(trinodes),2);
distance=sqrt( (xc-x(psi(1))).^2+(yc-y(psi(1))).^2);

x=x-x(psi(1));
y=y-y(psi(1));
plot_range=[-10000, 5000, -5000, 2000];
xc=mean(x(trinodes),2);
yc=mean(y(trinodes),2);

max_distance =5000;

speed=sqrt(ua.^2+va.^2);

hc=mean(h(trinodes),2);

%U_tide_data
analysis.annual.run = 0;        
analysis.annual.SNR = 3;
analysis.annual.Rayleigh = [0.97,1];    %first includes annual consituents, second excludes
analysis.annual.dt = 10/(24*60);        %10 minute resolution on prediction

ii=100000
coef = ut_solv(time, double(ua(:,ii)), double(va(:,ii)), uvllnode(ii,2), 'auto','NoTrend','Rmin',analysis.annual.Rayleigh(1),'OLS','NoDiagn','LinCI');
save speedCoef.m coef
    
coef = ut_solv(time, double(ua(:,ii)), [], uvllnode(ii,2), 'auto','NoTrend','Rmin',analysis.annual.Rayleigh(1),'OLS','NoDiagn','LinCI');
save elevCoef.m coef