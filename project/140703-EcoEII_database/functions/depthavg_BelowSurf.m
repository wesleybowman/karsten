function depthavg_BelowSurf(fileinfo,options)
% This function reads in the data (Polagye format) and then depth averages the
% variables from the bottom boundary layer to 95% below the surface signal.
%
% Justine McMillan
% Nov 14, 2013 (modified from an earlier file used for swnsRA)
%
% Updates:
% - Apr 15/14: No longer save 'lon'/'lat' variables since they are in
%              'param' structure
%- June 13/14: Added metadata
%            Added additional fields that are now in data structure
% - July 3/14: Modified the "fields" definition to be more automated


clim = 2.5;
load([fileinfo.outdir fileinfo.flowfile])
data



%% Initializing
nt = length(time.mtime);
dataavg.davglims = NaN*ones(nt,2);

%fields = {'mag_signed_vel','mag_signed_vel_std','east_vel','east_vel_std',...
   % 'north_vel','north_vel_std','vert_vel','vert_vel_std','error_vel',...
   % 'error_vel_std','dir_vel','dir_vel_std','Ualong','Ualong_std','Ucross',...
   % 'Ucross_std'}
fields = fieldnames(data);

for ii =1:length(fields)
    if strcmp(fields{ii},'bins')
        continue;
    else
        dataavg.(fields{ii}) = NaN *ones(nt,1);
    end
end

%% Depth averaging to 95 percent of surface -- calculated at every time step
dmin = min(data.bins);
for ii = 1:length(time.mtime)
    if ~isnan(pres.surf(ii))
        
        ind = find(data.bins<0.95*pres.surf(ii));
        dmax = data.bins(ind(end));
        
        dataavg.davglims(ii,:) = [dmin dmax];
        
        for jj = 1:length(fields)
            if strcmp(fields{jj},'bins')
                continue;
            else
                dataavg.(fields{jj})(ii) = calc_depthavg(data.(fields{jj})(ii,:),data.bins,dataavg.davglims(ii,:),2);
            %     dataavg.mag_signed_vel(ii) = calc_depthavg(data.mag_signed_vel(ii,:),data.bins,dataavg.davglims(ii,:),2);
            %    % dataavg.dir_vel(ii) = calc_depthavg(data.dir_vel(ii,:),data.bins,dataavg.davglims(ii,:),2);
            %     dataavg.east_vel(ii) = calc_depthavg(data.east_vel(ii,:),data.bins,dataavg.davglims(ii,:),2);
            %     dataavg.north_vel(ii) = calc_depthavg(data.north_vel(ii,:),data.bins,dataavg.davglims(ii,:),2);
            %     dataavg.vert_vel(ii) = calc_depthavg(data.vert_vel(ii,:),data.bins,dataavg.davglims(ii,:),2);
            %     dataavg.error_vel(ii) = calc_depthavg(data.error_vel(ii,:),data.bins,dataavg.davglims(ii,:),2);
            end
        end
    end
end
dataavg.bins = nanmean(dataavg.davglims(:,2))/2;
dataavg.dir_vel = get_DirFromN(dataavg.east_vel,dataavg.north_vel);

%% Plot to check surface signal
if options.showBS
    figure,clf
    ax(1) = subplot(211);
    set(gcf,'Position',[  87  132  1544  841])
    imagesc(get_yd(time.mtime),data.bins,data.mag_signed_vel',[-1 1]*clim);
    set(gca,'YDir','Normal')
    ylabel('range (m)')
    xlabel('time')
    c= colorbar;
    ylabel(c,'speed (m/s)')
    hold on
    plot(get_yd(time.mtime),pres.surf,'k','Linewidth',2)
    if time.mtime(end)-time.mtime(1)>8
        xlim([mean(get_yd(time.mtime)),mean(get_yd(time.mtime))+2])
    end
    plot(get_yd(time.mtime),dataavg.davglims(:,1),'w','Linewidth',2)
    plot(get_yd(time.mtime),dataavg.davglims(:,2),'w','Linewidth',2)
    ax(2) = subplot(212);
    plot(get_yd(time.mtime),dataavg.mag_signed_vel)
    ylabel('signed speed (m/s)')
    colorbar
    linkaxes(ax,'x')
end


%% Save
data = dataavg;
Comments{end+1} = ['Data is depth averaged between data.davglims (about 95 percent of surface signal)'];
Comments{end+1} = ['Data.bins = nanmean(dataavg.davglims(:,2))/2 (Kind of a middepth)'];

%if exist([fileinfo.outdir, fileinfo.BSfile,'.mat'],'file')
%    SAVEFLG = input(['Do you want to overwrite ',...
%        fileinfo.outdir, fileinfo.BSfile,'? (1 = yes) \n']);
%else
SAVEFLG = 1;
%end

if SAVEFLG
    disp(['Saving data to ',fileinfo.outdir, fileinfo.BSfile])
    
    clear metadata
    metadata.fullflowfile = [fileinfo.outdir fileinfo.flowfile];
    metadata.progname = mfilename('fullpath');
    metadata.date = datestr(now);
    save([fileinfo.outdir, fileinfo.BSfile],'data','pres','time','params','Comments','metadata');
else
    disp('DID NOT SAVE OUTPUT FILE')
end