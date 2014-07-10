%% Generate depth averaged below surface data
% Example script to give to Wesley
%
% JMM
% July 3, 2014

clear
P = genpath('../functions/');
addpath(P);

%%
site = 'GP-120726-BPd';

%% File info
fileinfo.flowfile = [site,'_','Flow'];
fileinfo.outdir = ['../data/'];

fid = dir([fileinfo.outdir,fileinfo.flowfile,'.mat']);
disp(['Input file created on ',datestr(fid.datenum)])

fileinfo.BSfile = [site,'_davgBS'];

%% Run function
options.showBS = 1;
depthavg_BelowSurf(fileinfo,options)
