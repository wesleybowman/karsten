from oct2py import octave
from loadnc2d import loadnc2d as myfunction

filename = '/home/wesley/github/aidan-projects/grid/dngrid_0001.nc'


x,y,ua,va,trinodes,el,h,time,siglev,siglay,nbe,a1u,a2u,aw0,awx,awy,nele,node,lon,lat = octave.call('loadnc2d',filename)
x,y,ua,va,trinodes,el,h,time,siglev,siglay,nbe,a1u,a2u,aw0,awx,awy,nele,node,lon,lat = myfunction(filename)
