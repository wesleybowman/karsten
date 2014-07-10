%% Generate 10 minute ensemble data
% Example script to give to Wesley
%
% JMM
% July 3, 2014

clear all
P = genpath('../functions/');
addpath(P);
%%
site = 'GP-120726-BPd';

%% load flow data
load(['../data/',site,'_Flow'])

%% choose ensemble average and output file name
t_ens = 10;
outputfile = ['../data/',site,'_',num2str(t_ens),'minavg']

%% Run ensemble averaging script
options.method = 'TimeSearch';

bins = data.bins;
[time,data,pres] = EnsembleData_FlowFile(t_ens*60,time,data,pres,options);
data.bins = bins;

%% nanabove the surface
fields = fieldnames(data);
for ii=1:length(fields)
    if strcmp(fields{ii},'bins')
        continue;
    else
        field_in = data.(fields{ii});
        [data.(fields{ii}), bin_max] = nan_AboveSurf(field_in,data.bins,pres.surf,0.95);
    end
end

%% save 
Comments{end+1} = 'Data was ensemble averaged from files already in BP format';

lon = params.lon;
lat = params.lat;
disp(['Saving data to ',outputfile])
save(outputfile,'data','pres','time','lon','lat','params','Comments','-v7.3')

%% Save metadata
metadata.progname=[mfilename('fullpath')];
metadata.date = datestr(now);
metadata.paramfile = metadata.paramfile;
save(outputfile,'metadata','-append')



