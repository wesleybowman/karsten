%if 1==1
%dd=cd;
    %[x,y,ua,va,trinodes,el,h,time,siglev,siglay,nbe,a1u,a2u,nele,aw0,awx,awy,node,long,lat]=loadnc2d_matlab_RK(dd,'/smallcape_0001.nc');
    [x,y,ua,va,trinodes,el,h,time,siglev,siglay,nbe,a1u,a2u,aw0,awx,awy,nele,node,long,lat] = loadnc2d_matlab_RK('/home/wesley/github/aidan-projects/grid','/dngrid_0001.nc');
[nodexy,uvnode,dt,deltat,hour,thour,TP,rho,g,period,nodell,uvllnode]=ncdatasort(x,y,time*24*3600,trinodes,long,lat);
time=mjd2num(time);
%end
%depth restrictions
min_depth=30
max_depth=80

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

%speed(:,hc<min_depth)=0;
%speed(:,hc>max_depth)=0;
%speed(:,distance>max_distance)=0;
radius=10;
turbine.area=pi*radius^2;
turbine.Cp=0.4;
turbine.rho=1025;
turbine.cutin=1;
turbine.rated=2.5;
turbine.cutoff=5;
turbine.Prated=turbine.Cp*(1/2*turbine.rho*turbine.area*turbine.rated^3);
turbine.meanP=0;
turbine.cf=0.40;
u_rated=0:0.25:8;

load('turbines.mat')
%turbines=cf_u_rated_turbs(speed,turbine,u_rated);

% if 1==1
% figure
% patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',[turbines.rated]')
% shading flat;
% axis(plot_range)
% colorbar
% caxis([0 5])
% drawnow
% end

% for jj=1:length(urated)
% %turbine characteristics
% turbine(jj).cutin=1;
% turbine(jj).rated=urated(jj);
% turbine(jj).cutoff=max(urated(jj)+10);
% turbine(jj).Prated=Cp*(1/2*rho*turbine_area*turbine(jj).rated^3);
%     P=calculate_power_RK(speed(:,jj),time,turbine(jj));
% meanP(jj)=mean(P);
% (jj)=meanP(jj)/turbine(jj).Prated*100;
% end
meanP=[turbines.meanP];
capacity_factor=[turbines.cf];
capacity_factor(isnan(capacity_factor))=0;
tp=[turbines.Prated];
% if 1==1
% figure
% patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',meanP')
% shading flat;
% colorbar
% axis(plot_range)
% 
% figure
% patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',capacity_factor')
% shading flat;
% colorbar
% axis(plot_range)
% 
% figure
% patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',tp')
% shading flat;
% colorbar
% axis(plot_range)
% end

%U_tide_data
analysis.annual.run = 0;        
analysis.annual.SNR = 3;
analysis.annual.Rayleigh = [0.97,1];    %first includes annual consituents, second excludes
analysis.annual.dt = 10/(24*60);        %10 minute resolution on prediction

%locate N turnbines
N=10000;

%turbine spacing
turbine_diameter=10;
spacing_along=15*turbine_diameter;
spacing_across=2*turbine_diameter;

turbine_power=meanP;

%score=(55-capacity_factor)/10+abs(hc'-40)/10 +(max(distance',2500)-2500)/500;
%score=20*(1-meanP/1e6)+abs(hc'-40)/40 +((max(distance',3000)-3000)/250).^1;
score=20*(1-meanP/1e6);
turbine_score=score;
%Find best location
%figure
ii=1

for ii=1:N
%    [~,loci(ii)]=max(turbine_power);
    [~,loci(ii)]=min(turbine_score);
    
    % do u_tide analysis at loc
    coef = ut_solv(time, double(ua(:,loci(ii))), double(va(:,loci(ii))), uvllnode(loci(ii),2), 'auto','NoTrend','Rmin',analysis.annual.Rayleigh(1),'OLS','NoDiagn','LinCI');
    save coef.m coef
    cx=cos(coef.theta(1)*pi/180);
    cy=sin(coef.theta(1)*pi/180);    
    %find new xy as distance from location, and in direction of M2 ellipse
    newx=(xc-xc(loci(ii)))*cx+(yc-yc(loci(ii)))*cy;
    newy=-(xc-xc(loci(ii)))*cy+(yc-yc(loci(ii)))*cx;
    a=find(abs(newx)<spacing_along & abs(newy)<spacing_across);
%    turbine_power(a)=0;    
    turbine_score(a)=100;    
if 1==1
clf
patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',turbine_score')
shading flat;
axis(plot_range)
colorbar
caxis([0 40])
drawnow
end

end


total_power=sum(meanP(loci))/1e6

if 1==1
figure
patch('Vertices',[x,y],'Faces',trinodes,'FaceVertexCdata',hc)
shading flat;
colorbar
caxis([0 100])
hold on
plot(xc(loci),yc(loci),'.m')
hold off
axis(plot_range)
end

save plausible_P.mat loci turbine_power turbines x y long lat trinodes h distance


 
