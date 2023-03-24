% Load BGC NetCDF files and extract only the variables interested in 
% L Ziegler. HPL Sep 2020

idir = dir('*.nc');

kb=11;
for k = 1 : kb;
    s_hybrid(k,1:length(xyd_n))=-((k-1)/(kb-1))^1.5;
end

for i = 2:11;
   layer_thickness(i-1,:)=mean(s_hybrid(i-1:i,:),1)*-1;
end

alg1_mod=[];
alg2_mod=[];
DO_mod=[];
NH4_mod=[];
NO3_mod=[];
depth1=[]; 
DON_mod=[];
DOC_mod=[];
POC_mod=[];
PON_mod=[];
SALT_mod=[];
TEMP_mod=[];
KD_mod=[];
mod_time=[];
JWNCDOC1_mod=[];
JWNCDOC2_mod=[];
JWNCDOC3_mod=[];
JWCDOC1_mod=[];
JWCDOC2_mod=[];
JWCDOC3_mod=[]; 

numfile = length(idir);
for idir = 1 : numfile;   
    numfile_count = [1:numfile];   
    kstr          = num2str(numfile_count(idir),'%04d');    
    fname         = [wqm_hisdir '/' nameseg,'_',kstr,'.nc'];
    nc            = netcdf([fname]);
    
    idir
    
    alg1_in  = nc{'B1'}(:);
    alg2_in  = nc{'B2'}(:);
    DO_in    = nc{'DOXG'}(:);
    NH4_in   = nc{'NH4'}(:);
    NO3_in   = nc{'NO3'}(:);
    depth_in = nc{'depth'}(:);
    DON_in   = nc{'WC_CDON1'}(:) + ...
                                nc{'WC_NCDON1'}(:) + ...
                                nc{'WC_CDON2'}(:) + ...
                                nc{'WC_NCDON2'}(:) + ...
                                nc{'WC_CDON3'}(:) + ...
                                nc{'WC_NCDON3'}(:);
    DOC_in   = nc{'WC_CDOC1'}(:) + ...
                                nc{'WC_NCDOC1'}(:) + ...
                                nc{'WC_CDOC2'}(:) + ...
                                nc{'WC_NCDOC2'}(:) + ...
                                nc{'WC_CDOC3'}(:) + ...
                                nc{'WC_NCDOC3'}(:);
                            
    POC_in  = nc{'LPOC'}(:) + nc{'RPOC'}(:);
    PON_in  = nc{'LPON'}(:) + nc{'RPON'}(:);
    
    SALT_in = nc{'salinity'}(:);
    TEMP_in = nc{'temp'}(:);
    KD_in   = nc{'KD'}(:);
    
    
    depth1   =[depth1 ; depth_in];
    mod_time =[mod_time;nc{'time'}(:)./86400];
    
    alg1_mod =[alg1_mod ; alg1_in];
    alg2_mod =[alg2_mod ; alg2_in];
   
    DO_mod   =[DO_mod ; DO_in];
    NH4_mod  =[NH4_mod ; NH4_in];
    NO3_mod  =[NO3_mod ; NO3_in];
    
    
    DON_mod  =[DON_mod ; DON_in];
    DOC_mod  =[DOC_mod ; DOC_in];
    
    chla_mod = (alg1_mod./50+alg2_mod./50).*1000;
    PON_mod  =[PON_mod ; PON_in]+(alg1_mod+alg2_mod)./5.68;
    POC_mod  =[POC_mod ; POC_in]+alg1_mod+alg2_mod;
    
    SALT_mod =[SALT_mod ; SALT_in];
    TEMP_mod =[TEMP_mod ; TEMP_in];
    
    KD_mod   =[KD_mod ; KD_in];
        
    clear nc
    
end

[ntimes,nlayers,nstations]=size(big_chla);

    hr_table = struct2table(hr_exo);
    hr_name = hr_table.station_description;
    hr_lon = hr_table.lon;
    hr_lat = hr_table.lat;
    unique_stations = unique(hr_name);
    hr_time = hr_table.datetime;
    hr_time = datestr(hr_time, 'yyyy-mm-dd');

    temp_exo = zeros(38,1);
    salt_exo = zeros(38,1);
    do_exo = zeros(38,1);
    chl_exo = zeros(38,1);
    turb_exo = zeros(38,1);

    for i = 1:38;
    temp_mean = mean(hr_table.exo(i,1).temp_C);
    salt_mean = mean(hr_table.exo(i,1).sal_psu);
    do_mean = mean(hr_table.exo(i,1).do_mgL);
    chl_mean = mean(hr_table.exo(i,1).chl_ugL);
    turb_mean = mean(hr_table.exo(i,1).turb_FNU);
    temp_exo(i,1) = temp_mean;
    salt_exo(i,1) = salt_mean;
    do_exo(i,1) = do_mean;
    chl_exo(i,1) = chl_mean;
    turb_exo(i,1) = turb_mean;
    end

idir = dir('*.nc');
numfile = 365;
nameseg = 'tonic';
temp = [];
salt = [];
mod_time = [];
for idir = 1 : numfile;
numfile_count = [1:numfile];
kstr = num2str(numfile_count(idir),'%04d');
fname = [nameseg,'_',kstr,'.nc'];
nc = netcdf([fname]);
temp_mod= squeeze(nc{'temp'}(:,1,:));
salt_mod= squeeze(nc{'salinity'}(:,1,:));
time= nc{'time'}(:)./24/86400;
clear nc
temp = [temp; temp_mod];
salt = [salt; salt_mod];
mod_time = [mod_time; time];
end
% daily average
temp_daily = zeros(365,7430);
salt_daily = zeros(365,7430);
time_daily = zeros(365,1);
for ii = 1:365;
temp_daily(ii,:) = mean(temp((ii-1)*24+1:end,:));
salt_daily(ii,:) = mean(salt((ii-1)*24+1:end,:));
time_daily(ii,1) = mean(mod_time((ii-1)*24+1:end,:));
end

figure(2)
close
set(gcf,'position',[604 629 1469 716],'color','w');
% 2018
for i = 1 : 22;
set(gca,'color','w','fontsize',14);
subplot(1,3,1)
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on ;
plot(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),'r.','markersize',25)
text(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),hr_name(i),'fontsize',12)
title('2018', 'fontweight','normal')
end
hold on;
% 2019
for i = 23:32;
set(gca,'color','w','fontsize',14);
subplot(1,3,2)
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on ;
plot(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),'r.','markersize',25)
text(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),hr_name(i),'fontsize',12)
title('2019','fontweight','normal')
end
hold on;
% 2020
for i = 33:39;
set(gca,'color','w','fontsize',14);
subplot(1,3,3)
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on ;
plot(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),'r.','markersize',25)
text(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),hr_name(i),'fontsize',12)
title('2020','fontweight','normal')
end
shg
saveas('HRsample_stations_2018_2020.png');
saveas(gcf,'sample_stations_2018_2020.png');

figure(3)
close
%set(gcf,'Renderer', 'painters', 'Position', [0.2 0.2 900 1200])
%set(gca,'CLim',[-0.5 28])
%set(gca,'CLim',[0 12])
subplot(2,3,1)
set(gca,'color','w','fontsize',16)
% patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',salt_daily(178,:),'facecolor','interp','edgecolor','interp');
% hold on
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%colormap(matlab_colormap);
colormap(jet(15))
caxis([0 15])
hbar = colorbar('YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'Chlorophyll (ug/l)','fontsize',13)
hold on
scatter(lld_n(stations_mod(1:7),2),lld_n(stations_mod(1:7),3),50, chl_exo(1:7,1),'o','filled', 'markeredgecolor','k');
hold on
%set(gca,'CLim',[0 30])
subplot(2,3,2)
set(gca,'color','w','fontsi ze',16)
%patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',salt_daily(277,:),'facecolor','interp','edgecolor','interp');
hold on
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%colormap(matlab_colormap);
colormap(jet(15))
caxis([0 15])
hbar = colorbar('YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'Chlorophyll (ug/l)','fontsize',13)
hold on
scatter(lld_n(stations_mod(8:12),2),lld_n(stations_mod(8:12),3),50, chl_exo(8:12,1),'o','filled', 'markeredgecolor','k');
hold on
%set(gca,'CLim',[0 30])
subplot(2,3,3)
set(gca,'color','w','fontsize',16)
%patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',salt_daily(303,:),'facecolor','interp','edgecolor','interp');
hold on
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%colormap(matlab_colormap);
colormap(jet(15))
caxis([0 15])
hbar = colorbar('YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'Chlorophyll (ug/l)','fontsize',13)
hold on
scatter(lld_n(stations_mod(13:22),2),lld_n(stations_mod(13:22),3),50, chl_exo(13:22,1),'o','filled', 'markeredgecolor','k');
hold on
%set(gca,'CLim',[0 30])
subplot(2,3,4)
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%colormap(matlab_colormap);
colormap(jet(15))
caxis([0 15])
hbar = colorbar('YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'Chlorophyll (ug/l)','fontsize',13)
hold on
scatter(lld_n(stations_mod(23),2),lld_n(stations_mod(23),3),50, chl_exo(23,1),'o','filled', 'markeredgecolor','k');
hold on
%set(gca,'CLim',[0 30])
subplot(2,3,5)
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%colormap(matlab_colormap);
colormap(jet(15))
caxis([0 15])
hbar = colorbar('YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'Chlorophyll (ug/l)','fontsize',13)
hold on
scatter(lld_n(stations_mod(24:32),2),lld_n(stations_mod(24:32),3),50, chl_exo(24:32,1),'o','filled', 'markeredgecolor','k');
hold on
%set(gca,'CLim',[0 30])
subplot(2,3,6)
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%colormap(matlab_colormap);
colormap(jet(15))
caxis([0 15])
hbar = colorbar('YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'Chlorophyll (ug/l)','fontsize',13)
hold on
scatter(lld_n(stations_mod(33:38),2),lld_n(stations_mod(33:38),3),50, chl_exo(33:38,1),'o','filled', 'markeredgecolor','k');
shg
