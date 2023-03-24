% Extract lat and lon from coastline shape file

% coast = shaperead('/Volumes/Storage2/Ziegler_Lisa/Lisa/HousatonicRiver.shp','UseGeoCoords',true);
coast_file = shaperead('/Users/lisaziegler/Desktop/North_Atlantic_clip.shp');
coast2 = coast_file;

lat_big=[];
lon_big=[];
for i = 1:length(coast2)
    lon = coast2(i).X';
    lat = coast2(i).Y';
    lon_big = [lon_big; lon];
    lat_big = [lat_big;lat];
end 

lat_lon = [41.17 -72.96];
st_name = 'H2';
station_mod = dsearchn([lld_n(:,2) lld_n(:,3)],([lat_lon(:,2) ; lat_lon(:,1)])'); %find the model nodes that are close to the wqm station points

 tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n);
 hold on
 
 for i = 1 : length(station_mod);
    plot(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),'ro')
    text(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),st_name(i),'fontsize',24)
    
 end
 
 saveas(gcf,'AllStations_map','png');
 
%% loop start
load('/Users/lisaziegler/Desktop/Housatonic_fvcom/BLUEGREEN_CMAP.mat');
load('/Users/lisaziegler/Desktop/Housatonic_fvcom/blue_orange_div_cmap.mat');
%load('/Users/lisaziegler/Desktop/Housatonic_fvcom/fvcom43_inputs/Inputs_v5/fvcom_grid_v5.mat');
load('/Volumes/Storage2/Ziegler_Lisa/run365/Tonic_26Jul.mat')
%load('/Users/lisaziegler/Desktop/Housatonic_fvcom/fvcom_grid/tonic_grid_nanfree.mat')
%mkdir('tempsalt_plots');

lld_nan = find(lld_n(:,4) >=-140 & lld_n(:,4)<=-1);
%%
idir = dir('*.nc'); %dir('*.nc');
close all
for ii=1:length(idir)
    if ii<10    
     file_index=['000',num2str(ii)];
    elseif ii>=10 && ii<=99
     file_index=['00',num2str(ii)];
    else
     file_index=['0',num2str(ii)];
    end
 
    fname = 'tonic_0181.nc';   
    %fname=['tonic_',file_index,'.nc'];
    salt=nc_varget(fname,'salinity');
    temp=nc_varget(fname,'temp');
    time=nc_varget(fname,'time');
    MJD_epoch = 86400; % use with output from fvcom2.7
    %MJD_epoch='Nov 17, 1858,00:00'; % use with output from fvcom4.3

    
    for j = 1:24:length(time)
    date = datestr(time(j)./MJD_epoch); %use with output from fvcom2.7
    %date = datestr(time(j)+datenum(MJD_epoch)); % use with output from
                                                    %fvcom4.3
    stemp = squeeze(temp(:,1,:));
    ssalt = squeeze(salt(:,1,:));
    btemp = squeeze(temp(:,9,:));
    bsalt = squeeze(salt(:,9,:));
    %   stemp(:,lld_nan) = NaN;


    figure(200)
    %set(gcf,'Renderer', 'painters', 'Position', [0.2 0.2 1500 800])
    %set(gcf,'Renderer', 'painters', 'Position', [0.2 0.2 900 1200])
    %set(gca,'CLim',[-0.5 28])
    
    subplot(2,2,1)
    set(gca,'CLim',[0 30])
    set(gca,'color','w','fontsize',16)
    patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',stemp(j,:),'facecolor','interp','edgecolor','interp');
    hold on
    plot(lon_big,lat_big,'-k','linewidth',1)
    xlim([-73.18 -72.9])
    ylim([41.06 41.4])
    %colormap(matlab_colormap);
    colormap(jet(15))
    %caxis([0 15])
    hbar = colorbar('Ticks',[0,2,4:2:26:2:30],'YColor',[0 0 0]);
    % Now create a label for each tick mark (you can modify these however you want)
    ylabel(hbar,'Temperature (degC)','fontsize',13)
    title(['Surface Temperature on day ',date]);
    % 
    hold on;
    % 
     subplot(2,2,2)
    set(gca,'CLim', [0 32]);
    set(gca,'color','w','fontsize',16)
    patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',ssalt(j,:),'facecolor','interp','edgecolor','interp');
    hold on
    plot(lon_big,lat_big,'-k','linewidth',1)
    xlim([-73.18 -72.9])
    ylim([41.06 41.4])
    %colormap(BlueSpec_colormap);
     colormap(jet(15))
    %caxis([0 30])
    hbar = colorbar('Ticks',[0:2:32],'YColor',[0 0 0]);
    ylabel(hbar,'salinity (PSU)','fontsize',13)
    title(['Surface Salt on day ',date]);
    % 
hold on;
    
    subplot(2,2,3)
    set(gca,'CLim',[0 30])
    set(gca,'color','w','fontsize',16)
    patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',btemp(j,:),'facecolor','interp','edgecolor','interp');
    hold on
    plot(lon_big,lat_big,'-k','linewidth',1)
    xlim([-73.18 -72.9])
    ylim([41.06 41.4])
    %colormap(matlab_colormap);

    %caxis([0 15])
    hbar = colorbar('Ticks',[0,2,4:2:26:2:30],'YColor',[0 0 0]);
    % Now create a label for each tick mark (you can modify these however you want)
    ylabel(hbar,'Temperature (degC)','fontsize',13)
    title(['Near Bottom Temperature on day ',date]);
    % 
    hold on;
    % 
     subplot(2,2,4)
    set(gca,'CLim', [0 32]);
    set(gca,'color','w','fontsize',16)
    patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',ssalt(j,:),'facecolor','interp','edgecolor','interp');
    hold on
    plot(lon_big,lat_big,'-k','linewidth',1)
    xlim([-73.18 -72.9])
    ylim([41.06 41.4])
    %colormap(BlueSpec_colormap);
     colormap(jet(15))
    %caxis([0 30])
    hbar = colorbar('Ticks',[0:2:32],'YColor',[0 0 0]);
    ylabel(hbar,'salinity (PSU)','fontsize',13)
    title(['Near Bottom Salt on day ',date]);
     drawnow
    % % 
% saveas(gcf,['tempsalt_plots/plot_ts' num2str(j,'%02d'), '.png']);



end 

end
%%
% salinity

for ii=1:length(idir)
%     if ii<10    
%      file_index=['000',num2str(ii)];
%     elseif ii>=10 && ii<=99
%      file_index=['00',num2str(ii)];
%     else
%      file_index=['0',num2str(ii)];
%     end
%     
% fname=['tonic_',file_index,'.nc'];
% salt=nc_varget(fname,'salinity');
% temp=nc_varget(fname,'temp');
% time=nc_varget(fname,'time');
% %MJD_epoch = 'Jan 01, 2017,00:00';%MJD_epoch='Nov 17, 1858,00:00';
% depth = nc_varget(fname,'zeta');

for j = 1:length(time)

%date = datestr(time(j)+datenum(MJD_epoch))
% stemp = squeeze(temp(:,1,:));
% ssalt = squeeze(salt(:,1,:));

date = datestr(time(j)+datenum(MJD_epoch));
clf
%set(gcf,'Renderer', 'painters', 'Position', [0.2 0.2 1500 800])
set(gcf,'Renderer', 'painters', 'Position', [0.2 0.2 900 1200])
%set(gca,'CLim',[-0.5 28])
% set(gca,'CLim',[0 30])
% %subplot(1,2,1)
% set(gca,'color','w','fontsize',16)
% patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',stemp(j,:),'facecolor','interp','edgecolor','interp');
% hold on
% plot(lon_big,lat_big,'-k','linewidth',1)
% xlim([-73.18 -72.9])
% ylim([41.06 41.4])
% colormap(matlab_colormap);
% %%colormap(jet(15))
% %caxis([0 15])
% hbar = colorbar('Ticks',[-0.5,2,4:2:26],'YColor',[0 0 0]);
% % Now create a label for each tick mark (you can modify these however you want)
% ylabel(hbar,'Temperature (degC)','fontsize',13)
% title(['Temperature on day ',date]);
% 
% hold on;
% % 
%  subplot(1,2,2)
set(gca,'color','w','fontsize',16)
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',ssalt(j,:),'facecolor','interp','edgecolor','interp');
hold on
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%%colormap(BlueSpec_colormap);
 colormap(jet(15))
%caxis([0 30])
hbar = colorbar('Ticks',[0:2:30],'YColor',[0 0 0]);
ylabel(hbar,'salinity (PSU)','fontsize',13)
title(['Salt on day ',date]);
% % 
 drawnow
% % 
% saveas(gcf,['tempsalt_plots/plot_ts' num2str(j,'%02d'), '.png']);



end 

end
%% BGC
idir = dir('*.nc');
close all
for ii=1:length(idir)
    if ii<10    
     file_index=['000',num2str(ii)];
    elseif ii>=10 && ii<=99
     file_index=['00',num2str(ii)];
    else
     file_index=['0',num2str(ii)];
    end
    
fname=['tonic_',file_index,'.nc'];
NO3=nc_varget(fname,'NO3');
poc = nc_varget(fname,'LPOC')+nc_varget(fname,'RPOC'); %+(nc_varget(fname,'B1')+nc_varget(fname,'B2'))/5.68;
KD = nc_varget(fname,'KD');
alg1 = nc_varget(fname,'B1');
alg2 = nc_varget(fname,'B2');
time=nc_varget(fname,'time');

chla = (alg1./50 + alg2./50).*1000; % Convert mg C/l to ug chl/l based on chl:C ratio from model
MJD_epoch = 86400;
%MJD_epoch='Nov 17, 1858,00:00';
date = time./MJD_epoch;
date_round = round(date);

for j = 1:24:length(time)
date = round(time(j)./MJD_epoch);
schla = squeeze(chla(:,1,:));
skd = squeeze(KD(:,1,:));
sno3 = squeeze(NO3(:,1,:));
spoc = squeeze(poc(:,1,:));
%date = datestr(time(j)+datenum(MJD_epoch));
%stemp = squeeze(temp(:,1,:));
%ssalt = squeeze(salt(:,1,:));

clf
%set(gcf,'Renderer', 'painters', 'Position', [0.2 0.2 1500 800])
set(gcf,'Renderer', 'painters', 'Position', [0.2 0.2 900 1200])
%set(gca,'CLim',[-0.5 28])
set(gca,'CLim',[0 20])
%subplot(1,2,1)
set(gca,'color','w','fontsize',16)
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',schla(j,:),'facecolor','interp','edgecolor','interp');
hold on
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
colormap(matlab_colormap);
%%colormap(jet(15))
%caxis([0 15])
hbar = colorbar;%('Ticks',[-0.5,2,4:2:26],'YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'POC mg Car','fontsize',13)
title(['poc on day ',num2str(date)]);
% 
% hold on;
% % 
%  subplot(1,2,2)
% set(gca,'color','w','fontsize',16)
% patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',ssalt(j,:),'facecolor','interp','edgecolor','interp');
% hold on
% plot(lon_big,lat_big,'-k','linewidth',1)
% xlim([-73.18 -72.9])
% ylim([41.06 41.4])
% colormap(BlueSpec_colormap);
%  colormap(jet(15))
% %caxis([0 30])
% hbar = colorbar('Ticks',[0:2:30],'YColor',[0 0 0]);
% ylabel(hbar,'salinity (PSU)','fontsize',13)
% title(['Salt on day ',date]);
% % 
 drawnow
% % 
% saveas(gcf,['tempsalt_plots/plot_ts' num2str(j,'%02d'), '.png']);



end 

end

%%
for ii=1:length(idir)
    if ii<10    
     file_index=['000',num2str(ii)];
    elseif ii>=10 && ii<=99
     file_index=['00',num2str(ii)];
    else
     file_index=['0',num2str(ii)];
    end
    
fname=['tonic_',file_index,'.nc'];
s_hybrid = nc_varget(fname,'siglay');
salt=nc_varget(fname,'salinity');
temp=nc_varget(fname,'temp');
time=nc_varget(fname,'time');
MJD_epoch='Nov 17, 1858,00:00';
depth = nc_varget(fname,'zeta');
JD = datestr((time)+datenum(MJD_epoch));
ntime = length(JD);

depth_out=[];
        for i =1:ntime;
            depth_out(i,:,:) = s_hybrid.*depth(i,:);
        end


 contourf(repmat(mod_time,1,nlayers),squeeze(depth_out(:,:,342),TEMP_mod(:,:),'LineStyle','none')
 hold on
 [cmin,cmax]=caxis;

end 