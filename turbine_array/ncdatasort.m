function [nodexy,uvnode,dt,deltat,hour,thour,TP,rho,g,period,nodell,uvnodell,time]=ncdatasort(x,y,time,trinodes,long,lat)


%takes x,y,time,trinodes in and makes some standard variables that are commonly used
%also contains constants use like hour in seconds,g, m2 tidal period,rho and period

hour=3600;
g=9.806;
TP=12.42;
rho=1026;
period=(TP*3600)/(2*pi);


time=double(time);
time=time+678942;


dt=time(2)-time(1);
thour=time/hour;
deltat=thour(2)-thour(1);


nodexy(:,1)=x;
nodexy(:,2)=y;

 uvnode(:,1)=(nodexy(trinodes(:,1),1)+...
       nodexy(trinodes(:,2),1)+...
       nodexy(trinodes(:,3),1))/3;
   uvnode(:,2)=(nodexy(trinodes(:,1),2)+...
      nodexy(trinodes(:,2),2)+...
      nodexy(trinodes(:,3),2))/3;



if nargin>4
	
nodell(:,1)=long;
nodell(:,2)=lat;

 uvnodell(:,1)=(nodell(trinodes(:,1),1)+...
       nodell(trinodes(:,2),1)+...
       nodell(trinodes(:,3),1))/3;
   uvnodell(:,2)=(nodell(trinodes(:,1),2)+...
      nodell(trinodes(:,2),2)+...
      nodell(trinodes(:,3),2))/3;
else
'No nodell, uvnodell set to uvnodexy'
uvnodell=uvnodexy;

end



