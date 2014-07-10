function rbrout=load_rbrdata(filename,varargin)
% Loads in an rbr data file and gets into correct format (names are
% slightly different from Doug's output



load(filename)
rbrout.mtime = rbr.yd;

rbrout.temp = rbr.temperature;
rbrout.pres = rbr.pressure;
rbrout.depth = rbr.depth;