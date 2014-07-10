%% Generate a BP format file from the raw data
% Example script to give to Wesley
%
% JMM
% July 3, 2014

clear all
P = genpath('../functions/');
addpath(P);

%% run script that sets parameters and paths to data files
run Params_Stn4_SWNSreport

%% load raw data
% ADCP
load ../data/GP-120726-BPd_raw.mat
% rbr (pressure sensor)
rbr = load_rbrdata([fileinfo.datadir fileinfo.rbr]);

%% set options
options.showPA = 1;
options.showRBRavg = 1;

%% save a flow file in BPformat
save_FlowFile_BPFormat(fileinfo,adcp,rbr,saveparams,options)
