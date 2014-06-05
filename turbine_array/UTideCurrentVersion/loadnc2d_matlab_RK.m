function [x,y,ua,va,trinodes,el,h,time,siglev,siglay,nbe,a1u,a2u,aw0,awx,awy,nele,node,long,lat] = loadnc2d(datadir,singlename)

%Loads a 2d ncfile
%First input is a directory must include a / at the end like '/path/to/my/dir/'
%If given only a path it loads the first named file
%Second input is a way to specify which single file you want to use, uses files name.

if nargin<2
files=dir([ '' datadir '*.nc' ]);
filename=['' datadir '' files(1).name '' ];
else
 filename=['' datadir '' singlename '']
end


ncid=netcdf.open(filename,'NC_NOWRITE');

x = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x'));
y = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y'));
long = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
ua = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ua'))';
va = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'va'))';
el = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'))';
h = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'));
nbe = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nbe'));
a1u = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'a1u'));
a2u = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'a2u'));
aw0 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'aw0'));
awx = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'awx'));
awy = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'awy'));
time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
trinodes = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nv'));
siglay = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'siglay'));
siglev = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'siglev'));
nele = netcdf.getVar(ncid,netcdf.inqDimID(ncid,'nele'));
node = netcdf.getVar(ncid,netcdf.inqDimID(ncid,'node'));
netcdf.close(ncid);










