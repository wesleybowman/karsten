% Created by:  
%  Justine McMillan 
%  17-12-2012 
%  Notes: Analysis Parameters used for the SWNS report. 
%
% Modified: July 3, 2014 - Changes to paths made so that it would work for Wesley
%

stn = 'GP-120726-BPd';
%% File information
fileinfo.datadir = 	'../data/'; 	 %path to raw data files 
fileinfo.ADCP = [stn '_raw']; 	 %name of ADCP file 
fileinfo.outdir = '../data/'; 	 %path to output directory
fileinfo.flowfile =  [stn,'_Flow']; 	 %name of output file with Flow data
fileinfo.rbr = ['station4_grandPassageII_RBRSN_011857.mat'];
fileinfo.paramfile = mfilename;

%% ADCP parameters
saveparams.tmin = 209; 	 %tmin (year day)
saveparams.tmax = 240; 	 %tmax (year day)
saveparams.zmin = 0;     %minimum z to include in saves file    
saveparams.zmax = 20;
saveparams.approxdepth = 15.5; %Approximate depth
saveparams.flooddir= 0; 	 %Flood direction (relative to true north, CW is positive)
saveparams.declination = -17.25;%Declination angle
saveparams.lat =  44.2605; 	 %latitude
saveparams.lon = -66.3354; 	 %longitude
saveparams.dabADCP = 0.5; 	 %depth above bottom of ADCP
saveparams.dabPS = -0.6; 	 %depth above bottom of pressure sensor	
saveparams.rbr_hr_offset = 3;    % hour offset to convert rbr time to UTC
