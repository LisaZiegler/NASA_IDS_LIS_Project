%%% Plot daily and monthly wqm variables 
% load the following
%load('/Users/lisaziegler/Desktop/Tonic_BGC/wqm_run2_2020/bgc_variables.mat')
load('/Volumes/LaCie/Housatonic_BGC/Data/LIS_Nutrients_matlab/ct_doc_chl2018.mat');
load('/Volumes/LaCie/Housatonic_BGC/Data/LIS_Nutrients_matlab/river_doc_chl2018.mat');
load('/Volumes/LaCie/Housatonic_BGC/Data/LIS_Nutrients_matlab/ct_riv_notes.mat');
load('/Volumes/LaCie/Housatonic_BGC/Data/LIS_Nutrients_matlab/s1_data.mat');
load('/Volumes/LaCie/Housatonic_BGC/Data/LIS_ALLnutrients_monthlysurface2018.mat');
load('/Volumes/LaCie/Housatonic_BGC/Data/LIS_ALLphysics_monthlysurface2018.mat');
load('/Volumes/LaCie/Housatonic_BGC/Housatonicwqm_grid,mat');

coast_file = shaperead('/Volumes/LaCie/Housatonic_BGC/Data/North_Atlantic_clip.shp');
coast2 = coast_file;

lat_big=[];
lon_big=[];
for i = 1:length(coast2)
    lon = coast2(i).X';
    lat = coast2(i).Y';
    lon_big = [lon_big; lon];
    lat_big = [lat_big;lat];
end 
%%
%%%%Extract all node points corresponding to sample sites lat/lon coordinates

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
lat_lis = [41.12,41.05,41.16,41.08,41.14,41.15,41.08,41.17,41.10];
lon_lis = [-73.09,-73.08,-73.01,-73.02,-72.94,-72.84,-73.16,-72.96,-72.93];
lisnames = {'18','19','21','22','23','27','F2','H2','H4'};
stations_mod_lis = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_lis ; lat_lis])');

lat_nodes=[];
lon_nodes=[];

for i = 1 : length(stations_mod_lis);
    lat_node = mean(lld_n(stations_mod_lis(i),2));
    lon_node = mean(lld_n(stations_mod_lis(i),3));
    lat_nodes(i)=lat_node;
    lon_nodes(i)=lon_node;
end 

chl_st = [];

for i=2:8;
chl_st{i-1} = riv_chl_6{i};
end
lat_r = riv_chl_6{9,1}';
lon_r = riv_chl_6{10,1}';
stations_mod_riv = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_r ; lat_r]'));

for i = 1 : length(stations_mod_riv);
    lat_node = mean(lld_n(stations_mod_riv(i),2));
    lon_node = mean(lld_n(stations_mod_riv(i),3));
    lat_nodes1(i)=lat_node;
    lon_nodes1(i)=lon_node;
end 

chl_ct_st = [];

for i=2:4;
chl_ct_st{i-1} = ct_chl_6{i};
end

lat_ct = ct_chl_6{5,1}';
lon_ct = ct_chl_6{6,1}';
stations_mod_ct = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_ct ; lat_ct]'));

for i = 1 : length(stations_mod_ct);
    lat_node = mean(lld_n(stations_mod_ct(i),2));
    lon_node = mean(lld_n(stations_mod_ct(i),3));
    lat_nodes2(i)=lat_node;
    lon_nodes2(i)=lon_node;
end 

chl_ct_st8 = [];

for i=2:4;
chl_ct_st8{i-1} = ct_chl_8{i};
end

lat_ct = ct_chl_8{5,1}';
lon_ct = ct_chl_8{6,1}';
stations_mod_ct8 = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_ct ; lat_ct]'));

for i = 1 : length(stations_mod_ct8);
    lat_node = mean(lld_n(stations_mod_ct8(i),2));
    lon_node = mean(lld_n(stations_mod_ct8(i),3));
    lat_nodes3(i)=lat_node;
    lon_nodes3(i)=lon_node;
end 

lat_ct = 41.1555;
lon_ct = -73.0894;
stations_mod_S1 = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_ct ; lat_ct]'));

for i = 1 : length(stations_mod_S1);
    lat_node = mean(lld_n(stations_mod_S1(i),2));
    lon_node = mean(lld_n(stations_mod_S1(i),3));
    lat_nodesS1(i)=lat_node;
    lon_nodesS1(i)=lon_node;
end 

lat_ct = riv_chl_12{8,1}';
lon_ct = riv_chl_12{9,1}';
stations_mod_riv12 = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_ct ; lat_ct]'));

for i = 1 : length(stations_mod_riv12);
    lat_node = mean(lld_n(stations_mod_riv12(i),2));
    lon_node = mean(lld_n(stations_mod_riv12(i),3));
    lat_nodes4(i)=lat_node;
    lon_nodes4(i)=lon_node;
end 

nut_stat = LISALLnutrientsmonthlysurface2018;
phy_stat = LISALLphysicsmonthlysurface2018;

%% MONTHLY CHLOROPHYLL PLOTS
%%%% Monthly Averaged Surface Chlorophyll

chl_tbl = nut_stat(nut_stat{:,4} == 'Chlorophyll a',:);
chl_stat = table2cell(chl_tbl);

close; 
set(gcf,'color','w');
set(gca, 'Fontsize',13)
chl_p = phy_stat(1:3,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);
chl1 = chl_stat(1,5);
chl2 = chl_stat(15,5);
chl3 = chl_stat(29,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;



subplot(4,3,1);

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), squeeze(temperature(:,1,15));
%cmin = min(mon_chla_pt(1,:));cmax = max(mon_chla_pt(1,:));
caxis([0 16]);
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on ;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
title('Jan')

hold on

chl_p = phy_stat(4:6,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);
chl1 = chl_stat(2,5);
chl2 = chl_stat(16,5);
chl3 = chl_stat(30,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;

subplot(4,3,2);

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(2,:));
%cmin = min(mon_chla_pt(2,:));cmax = max(mon_chla_pt(2,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
title('Feb')

hold on;


chl_p = phy_stat(7:9,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);
chl1 = chl_stat(3,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(31,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;

subplot(4,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(3,:));
%cmin = min(mon_chla_pt(3,:));cmax = max(mon_chla_pt(3,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
title('Mar');

hold on;

chl_p = phy_stat(10:12,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(4,5);
chl2 = chl_stat(18,5);
chl3 = chl_stat(32,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;
subplot(4,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(4,:));
%cmin = min(mon_chla_pt(4,:));cmax = max(mon_chla_pt(4,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
title('Apr')

hold on

chl_p = phy_stat(13:15,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(5,5);
chl2 = chl_stat(19,5);
chl3 = chl_stat(33,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;

subplot(4,3,5);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(5,:));
%cmin = min(mon_chla_pt(5,:));cmax = max(mon_chla_pt(5,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
title('May')

hold on

chl_p = phy_stat(16,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(6,5);

chl_val = chl1;
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;

chl_r = riv_chl_6{1,1};
chl_r_c = chl_r(~isnan(chl_r));
chl_ct = ct_chl_6{1,1};

subplot(4,3,6);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(6,:));
%cmin = min(mon_chla_pt(6,:));cmax = max(mon_chla_pt(6,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes1(:,[1:2,4:5,7]), lon_nodes1(:,[1:2,4:5,7]),45, chl_r_c ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes2, lon_nodes2,45, chl_ct' ,'o','filled', 'markeredgecolor','k');
hold on; 
scatter(lld_n(stations_mod(1:7),2),lld_n(stations_mod(1:7),3),50, chl_exo(1:7,1),'o','filled','markeredgecolor','k');
title('Jun')

hold on;

subplot(4,3,7);
chl_p = phy_stat(17:24,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(7,5);
chl2 = chl_stat(21,5);
chl3 = chl_stat(35,5);

chl_val = [0 0 0 0 0 chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


chl_ct = ct_chl_7{1,1};
chl_ct18 = chl_ct(29,1);
chl_ct19 = chl_ct(39,1);
chl_ct21 = chl_ct([13,28],1);
chl_ct22 = chl_ct(40,1);
chl_ct23 = chl_ct(14,1);
chl_ctf2 = chl_ct(31,1);
chl_cth2 = chl_ct([15,27],1);
chl_cth4 = chl_ct(43,1);

chl_ct_val = [chl_ct18;mean(chl_ct21);chl_ct22;chl_ct23;0;chl_ctf2;mean(chl_cth2);chl_cth4];

chl_avg = (chl_val + chl_p + chl_ct_val)./2;

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(7,:));
%cmin = min(mon_chla_pt(7,:));cmax = max(mon_chla_pt(7,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes(:,2), lon_nodes(:,2),45, chl_ct19 ,'o','filled', 'markeredgecolor','k');
title('Jul')

hold on;


chl_p = phy_stat(25:32,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(8,5);
chl2 = chl_stat(24,5);
chl3 = chl_stat(38,5);
chl_val = [0 0 0 0 0 chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

chl_avg = (chl_val + chl_p)./2;

chl_ct = ct_chl_8{1,1};

subplot(4,3,8)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(8,:));
%cmin = min(mon_chla_pt(8,:));cmax = max(mon_chla_pt(8,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes3, lon_nodes3,45, chl_ct ,'o','filled', 'markeredgecolor','k');
title('Aug')

hold on;
subplot(4,3,9)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(9,:));
%cmin = min(mon_chla_pt(9,:));cmax = max(mon_chla_pt(9,:));
caxis([0 16]); 
hold on;
colormap('winter');
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

hbar = colorbar;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Sep')

hold on;


chl_p = phy_stat(33:35,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(12,5);
chl2 = chl_stat(26,5);
chl3 = chl_stat(40,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



chl_ct = ct_chl_10{1,1};
chl_ct10 = chl_ct(1,1);
chl_ct10s = [chl_ct10; 0; 0];

chl_avg = (chl_val + chl_p + chl_ct10s)./2;

chl_s1 = s1_data_oct302018(1,1);
subplot(4,3,10)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(10,:));
%cmin = min(mon_chla_pt(10,:));cmax = max(mon_chla_pt(10,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodesS1, lon_nodesS1,45, chl_s1 ,'o','filled', 'markeredgecolor','k');
hold on; 
scatter(lld_n(stations_mod(8:22),2),lld_n(stations_mod(8:22),3),50, chl_exo(8:22,1),'o','filled','markeredgecolor','k');
title('Oct')

hold on;


chl_p = phy_stat(36:38,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(13,5);
chl2 = chl_stat(27,5);
chl3 = chl_stat(41,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

chl_ct = ct_chl_11{1,1};
chl_ct11 = chl_ct(8,1);

chl_ct11s = [chl_ct11; 0; 0];

chl_avg = (chl_val + chl_p + chl_ct11s)./2;

subplot(4,3,11)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(11,:));
%cmin = min(mon_chla_pt(11,:));cmax = max(mon_chla_pt(11,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
title('Nov');

hold on;
subplot(4,3,12)
chl_p = phy_stat(39,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(14,5);
chl_val = cell2mat(chl1);


chl_ct = ct_chl_12{1,1};
chl_ct12 = chl_ct(4,1);

chl_avg = (chl_val + chl_p + chl_ct12)./2;

riv_chl12 = riv_chl_12{1,1};

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(12,:));
%cmin = min(mon_chla_pt(12,:));cmax = max(mon_chla_pt(12,:));
caxis([0 16]); 
hold on;
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes4, lon_nodes4,45, riv_chl12 ,'o','filled', 'markeredgecolor','k');
title('Dec')

%% SPRING VERSUS SUMMER

close 
set(gcf,'color','w');
set(gca, 'Fontsize',13)
chl_p = phy_stat(7:9,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);
chl1 = chl_stat(3,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(31,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;

subplot(2,3,1);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla_pt(3,:));
%cmin = min(mon_chla_pt(3,1,:));cmax = max(mon_chla_pt(3,1,:));
caxis([0 16]); 
hold on;
%colormap(jet(15))
%Cmap = cmocean('deep');
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Mar');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');

hold on;

chl_p = phy_stat(10:12,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(4,5);
chl2 = chl_stat(18,5);
chl3 = chl_stat(32,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;
subplot(2,3,2);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla_pt(4,:));
%cmin = min(mon_chla_pt(3,1,:));cmax = max(mon_chla_pt(3,1,:));
caxis([0 16]);
hold on;
%Cmap = cmocean('deep');
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Apr');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');


hold on

chl_p = phy_stat(13:15,5);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(5,5);
chl2 = chl_stat(19,5);
chl3 = chl_stat(33,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;

subplot(2,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla_pt(5,:));
%cmin = min(mon_chla_pt(3,1,:));cmax = max(mon_chla_pt(3,1,:));
caxis([0 16]);
hold on;
%Cmap = cmocean('deep');
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('May');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');

hold on

chl_p = phy_stat(16,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(6,5);

chl_val = chl1;
chl_val = cell2mat(chl_val)';
chl_avg = (chl_val + chl_p)./2;

chl_r = riv_chl_6{1,1};
chl_r_c = chl_r(~isnan(chl_r));
chl_ct = ct_chl_6{1,1};

subplot(2,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla_pt(6,:));
%cmin = min(mon_chla_pt(3,1,:));cmax = max(mon_chla_pt(3,1,:));
caxis([0 16]);
hold on;
%Cmap = cmocean('deep');
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Jun');
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes1(:,[1:2,4:5,7]), lon_nodes1(:,[1:2,4:5,7]),45, chl_r_c ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes2, lon_nodes2,45, chl_ct' ,'o','filled', 'markeredgecolor','k');

hold on;

subplot(2,3,5);
chl_p = phy_stat(17:24,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(7,5);
chl2 = chl_stat(21,5);
chl3 = chl_stat(35,5);

chl_val = [0 0 0 0 0 chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


chl_ct = ct_chl_7{1,1};
chl_ct18 = chl_ct(29,1);
chl_ct19 = chl_ct(39,1);
chl_ct21 = chl_ct([13,28],1);
chl_ct22 = chl_ct(40,1);
chl_ct23 = chl_ct(14,1);
chl_ctf2 = chl_ct(31,1);
chl_cth2 = chl_ct([15,27],1);
chl_cth4 = chl_ct(43,1);

chl_ct_val = [chl_ct18;mean(chl_ct21);chl_ct22;chl_ct23;0;chl_ctf2;mean(chl_cth2);chl_cth4];

chl_avg = (chl_val + chl_p + chl_ct_val)./2;

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla_pt(7,:));
%cmin = min(mon_chla_pt(3,1,:));cmax = max(mon_chla_pt(3,1,:));
caxis([0 16]); 
hold on;
%Cmap = cmocean('deep');
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Jul');
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes(:,2), lon_nodes(:,2),45, chl_ct19 ,'o','filled', 'markeredgecolor','k');

hold on;


chl_p = phy_stat(25:32,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

chl1 = chl_stat(8,5);
chl2 = chl_stat(24,5);
chl3 = chl_stat(38,5);
chl_val = [0 0 0 0 0 chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

chl_avg = (chl_val + chl_p)./2;

chl_ct = ct_chl_8{1,1};

subplot(2,3,6)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla_pt(8,:));
%cmin = min(mon_chla_pt(3,1,:));cmax = max(mon_chla_pt(3,1,:));
caxis([0 16]);
hold on;
%Cmap = cmocean('deep');
colormap('winter');
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(ug chl/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Aug');
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes3, lon_nodes3,45, chl_ct ,'o','filled', 'markeredgecolor','k');
