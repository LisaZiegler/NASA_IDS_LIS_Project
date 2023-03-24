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
lat_nodes=[];lon_nodes=[];
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

%cpb_surf_h2 = [1; 4.1; 1.3; 3.1; 4.1; 1.5; 3.4; 3.1; 1.6; 5.3; 3.5; 6.2; 3.6; 3.2];

%%
%================== CHLOROPHYLL ========================================%

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

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_chla1(1,:));
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


%% DOC
%================== DOC ========================================%

%% MONTHLY DOC PLOTS
%%%% Monthly Averaged Surface DOC
chl_tbl = nut_stat(nut_stat{:,4} == 'Dissolved Organic Carbon',:);
chl_stat = table2cell(chl_tbl);


%
close;
set(gcf,'color','w');
set(gca, 'Fontsize',13)
chl1 = chl_stat(1,5);
chl2 = chl_stat(13,5);
chl3 = chl_stat(25,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


 
subplot(4,3,1);


tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(1,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hbar = colorbar;
Cmap = cmocean('-haline');
colormap(Cmap);
hold on;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Jan')

hold on

chl1 = chl_stat(2,5);
chl2 = chl_stat(14,5);
chl3 = chl_stat(26,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,2);

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(2,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Feb')

hold on;


chl1 = chl_stat(3,5);
chl2 = chl_stat(15,5);
chl3 = chl_stat(27,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(3,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Mar')

hold on;


chl1 = chl_stat(4,5);
chl2 = chl_stat(16,5);
chl3 = chl_stat(28,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

subplot(4,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(4,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Apr')

hold on


chl1 = chl_stat(5,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(29,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,5);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(5,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('May')

hold on

chl1 = chl_stat(6,5);
chl_val = cell2mat(chl1);

chl_r = riv_doc_6{1,1};
chl_ct = ct_doc_6{1,1};

subplot(4,3,6);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(6,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes1, lon_nodes1,45, chl_r ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes2, lon_nodes2,45, chl_ct' ,'o','filled', 'markeredgecolor','k');
hold on;

title('Jun')

hold on;

subplot(4,3,7);


chl1 = chl_stat(7,5);
chl2 = chl_stat(19,5);
chl3 = chl_stat(31,5);

chl_val = [0 0 0 0 0 chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


chl_ct = ct_doc_7{1,1};
chl_ct18 = chl_ct(29,1);
chl_ct19 = chl_ct(39,1);
chl_ct21 = chl_ct([13,28],1);
chl_ct22 = chl_ct(40,1);
chl_ct23 = chl_ct(14,1);
chl_ctf2 = chl_ct(31,1);
chl_cth2 = chl_ct([15,27],1);
chl_cth4 = chl_ct(43,1);

chl_ct_val = [chl_ct18;mean(chl_ct21);chl_ct22;chl_ct23;0;chl_ctf2;mean(chl_cth2);chl_cth4];

chl_avg = (chl_val  + chl_ct_val)./2;

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(7,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes(:,2), lon_nodes(:,2),45, chl_ct19 ,'o','filled', 'markeredgecolor','k');
title('Jul')

hold on;



chl1 = chl_stat(9,5);
chl2 = chl_stat(21,5);
chl3 = chl_stat(33,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



chl_ct = ct_doc_8{1,1};

subplot(4,3,8)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(8,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]);  
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes3, lon_nodes3,45, chl_ct ,'o','filled', 'markeredgecolor','k');
title('Aug')

hold on;
subplot(4,3,9)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(9,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Sep')

hold on;



chl1 = chl_stat(10,5);
chl2 = chl_stat(22,5);
chl3 = chl_stat(34,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



chl_ct = ct_doc_10{1,1};
chl_ct10 = chl_ct(1,1);
chl_ct10s = [chl_ct10; 0; 0];

chl_avg = (chl_val + chl_ct10s)./2;


subplot(4,3,10)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(10,:));
%cmin = min(mon_doc_pt(3,1,:));cmax = max(mon_doc_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;

ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
title('Oct')
hold on;




chl1 = chl_stat(11,5);
chl2 = chl_stat(23,5);
chl3 = chl_stat(35,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

chl_ct = ct_doc_11{1,1};
chl_ct11 = chl_ct(8,1);

chl_ct11s = [chl_ct11; 0; 0];

chl_avg = (chl_val + chl_ct11s)./2;
subplot(4,3,11)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(11,:));
%cmin = min(mon_doc_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_avg ,'o','filled', 'markeredgecolor','k');
title('Nov')

hold on;
subplot(4,3,12)


chl1 = chl_stat(12,5);
chl_val = cell2mat(chl1);


tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(12,:));
%cmin = min(mon_doc_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Dec')



%% SPRING VERSUS SUMMER
chl_tbl = nut_stat(nut_stat{:,4} == 'Dissolved Organic Carbon',:);
chl_stat = table2cell(chl_tbl);

chl1 = chl_stat(3,5);
chl2 = chl_stat(15,5);
chl3 = chl_stat(27,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(2,3,1);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(3,:));
%cmin = min(mon_doc_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Mar');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');

hold on;


chl1 = chl_stat(4,5);
chl2 = chl_stat(16,5);
chl3 = chl_stat(28,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

subplot(2,3,2);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(4,:));
%cmin = min(mon_doc_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Apr');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');


hold on


chl1 = chl_stat(5,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(29,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(2,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(5,:));
%cmin = min(mon_doc_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('May');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');

hold on

chl1 = chl_stat(6,5);
chl_val = cell2mat(chl1);

chl_r = riv_doc_6{1,1};
chl_ct = ct_doc_6{1,1};

subplot(2,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(6,:));
%cmin = min(mon_doc_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Jun');
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes1, lon_nodes1,45, chl_r ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes2, lon_nodes2,45, chl_ct' ,'o','filled', 'markeredgecolor','k');

hold on;

subplot(2,3,5);


chl1 = chl_stat(7,5);
chl2 = chl_stat(19,5);
chl3 = chl_stat(31,5);

chl_val = [0 0 0 0 0 chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


chl_ct = ct_doc_7{1,1};
chl_ct18 = chl_ct(29,1);
chl_ct19 = chl_ct(39,1);
chl_ct21 = chl_ct([13,28],1);
chl_ct22 = chl_ct(40,1);
chl_ct23 = chl_ct(14,1);
chl_ctf2 = chl_ct(31,1);
chl_cth2 = chl_ct([15,27],1);
chl_cth4 = chl_ct(43,1);

chl_ct_val = [chl_ct18;mean(chl_ct21);chl_ct22;chl_ct23;0;chl_ctf2;mean(chl_cth2);chl_cth4];

chl_avg = (chl_val  + chl_ct_val)./2;

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_doc1(7,:));
%cmin = min(mon_doc_pt_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt_pt(3,1,:));
caxis([0 3]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Jul');
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_avg ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes(:,2), lon_nodes(:,2),45, chl_ct19 ,'o','filled', 'markeredgecolor','k');

hold on;



chl1 = chl_stat(9,5);
chl2 = chl_stat(21,5);
chl3 = chl_stat(33,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



chl_ct = ct_doc_8{1,1};

subplot(2,3,6)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),mon_doc1(8,:));
%cmin = min(mon_doc_pt_pt_pt_pt(3,1,:));cmax = max(mon_doc_pt_pt_pt_pt(3,1,:));
caxis([0 3]);  
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Aug');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodes3, lon_nodes3,45, chl_ct ,'o','filled', 'markeredgecolor','k');


%% NO3

%================== NO3 ========================================%

%% MONTHLY NO3 PLOTS
%%%% Monthly Averaged Surface NO3
chl_tbl = nut_stat(nut_stat{:,4} == 'Nitrate + Nitrite',:);
chl_stat = table2cell(chl_tbl);

set(gcf,'color','w');
set(gca, 'Fontsize',13)
chl1 = chl_stat(1,5);
chl2 = chl_stat(13,5);
chl3 = chl_stat(25,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



subplot(4,3,1);


tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(1,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]);
hbar = colorbar;
cmocean('-haline');
hold on;

ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Jan')

hold on

chl1 = chl_stat(2,5);
chl2 = chl_stat(14,5);
chl3 = chl_stat(26,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,2);

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(2,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Feb')
hold on;


chl1 = chl_stat(3,5);
chl2 = chl_stat(15,5);
chl3 = chl_stat(27,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(3,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Mar')

hold on;


chl1 = chl_stat(4,5);
chl2 = chl_stat(16,5);
chl3 = chl_stat(28,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

subplot(4,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(4,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Apr')

hold on


chl1 = chl_stat(5,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(29,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,5);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(5,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('May')

hold on

chl1 = chl_stat(6,5);
chl_val = cell2mat(chl1);

subplot(4,3,6);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(6,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Jun')

hold on;

subplot(4,3,7);


chl1 = chl_stat(7,5);
chl2 = chl_stat(19,5);
chl3 = chl_stat(31,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(7,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Jul')

hold on;



chl1 = chl_stat(9,5);
chl2 = chl_stat(21,5);
chl3 = chl_stat(33,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,8)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(8,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Aug')

hold on;
subplot(4,3,9)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(9,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Sep')
hold on;



chl1 = chl_stat(10,5);
chl2 = chl_stat(22,5);
chl3 = chl_stat(34,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



subplot(4,3,10)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(10,:));
%cmin = min(mon_no3_pt(3,1,:));cmax = max(mon_no3_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Oct')
hold on;




chl1 = chl_stat(11,5);
chl2 = chl_stat(23,5);
chl3 = chl_stat(35,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,11)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(11,:));
%cmin = min(mon_no3_pt_pt(3,1,:));cmax = max(mon_no3_pt_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Nov')

hold on;
subplot(4,3,12)


chl1 = chl_stat(12,5);
chl_val = cell2mat(chl1);


tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(12,:));
%cmin = min(mon_no3_pt_pt(3,1,:));cmax = max(mon_no3_pt_pt(3,1,:));
caxis([0 0.6]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Dec')

%% SPRING VERSUS SUMMER

chl_tbl = nut_stat(nut_stat{:,4} == 'Nitrate + Nitrite',:);
chl_stat = table2cell(chl_tbl);

chl1 = chl_stat(3,5);
chl2 = chl_stat(15,5);
chl3 = chl_stat(27,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

set(gcf,'color','w');
set(gca, 'Fontsize',13)

subplot(2,3,1);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(3,:));
%cmin = min(mon_no3(3,1,:));cmax = max(mon_no3(3,1,:));
caxis([0 0.5]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Mar');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');

hold on;


chl1 = chl_stat(4,5);
chl2 = chl_stat(16,5);
chl3 = chl_stat(28,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

subplot(2,3,2);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(4,:));
%cmin = min(mon_no3(3,1,:));cmax = max(mon_no3(3,1,:));
caxis([0 0.5]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Apr');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');


hold on


chl1 = chl_stat(5,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(29,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(2,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(5,:));
%cmin = min(mon_no3(3,1,:));cmax = max(mon_no3(3,1,:));
caxis([0 0.5]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('May');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');

hold on

chl1 = chl_stat(6,5);
chl_val = cell2mat(chl1);

subplot(2,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(6,:));
%cmin = min(mon_no3(3,1,:));cmax = max(mon_no3(3,1,:));
caxis([0 0.5]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Jun');
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');

hold on;

subplot(2,3,5);


chl1 = chl_stat(7,5);
chl2 = chl_stat(19,5);
chl3 = chl_stat(31,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(7,:));
%cmin = min(mon_no3(3,1,:));cmax = max(mon_no3(3,1,:));
caxis([0 0.5]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Jul');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');


hold on;



chl1 = chl_stat(9,5);
chl2 = chl_stat(21,5);
chl3 = chl_stat(33,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(2,3,6)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_no31(8,:));
%cmin = min(mon_no3(3,1,:));cmax = max(mon_no3(3,1,:));
caxis([0 0.5]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Aug');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');

%% NH4 

%================== NH4 ========================================%

%% MONTHLY NH4 PLOTS
%%%% Monthly Averaged Surface NH4
chl_tbl = nut_stat(nut_stat{:,4} == 'Ammonia',:);
chl_stat = table2cell(chl_tbl);

close; 
set(gcf,'color','w');
set(gca, 'Fontsize',13)
chl1 = chl_stat(1,5);
chl2 = chl_stat(13,5);
chl3 = chl_stat(25,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



subplot(4,3,1);


tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(3,:));
%cmin = min(mon_nh4_pt(3,1,:)); 
%cmax = max(mon_nh4_pt(3,1,:));
caxis([0 0.1]);
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Mar');
hold on;
scatter(lat_nodes(:,[7,9]), lon_nodes(:,[7,9]),45, chl_val([1,3],1) ,'o','filled', 'markeredgecolor','k');

hold on

chl1 = chl_stat(2,5);
chl2 = chl_stat(14,5);
chl3 = chl_stat(26,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,2);

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(4,:));
%cmin = min(mon_nh4_pt(4,1,:));
%cmax = max(mon_nh4_pt(4,1,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Apr');
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');

hold on;

subplot(4,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(5,:));
%cmin = min(mon_nh4_pt(5,1,:));
%cmax = max(mon_nh4_pt(5,1,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('May');

hold on;


chl1 = chl_stat(4,5);
chl2 = chl_stat(16,5);
chl3 = chl_stat(28,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

subplot(4,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(6,:));
%cmin = min(mon_nh4_pt(6,1,:));
%cmax = max(mon_nh4_pt(6,1,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Jun');
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val(1,1) ,'o','filled', 'markeredgecolor','k');


hold on


chl1 = chl_stat(5,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(29,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

subplot(4,3,5);


chl1 = chl_stat(7,5);
chl2 = chl_stat(19,5);
chl3 = chl_stat(31,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(7,:));
%cmin = min(big_nh4_mon(12,:));cmax = max(big_nh4_mon(12,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val(1,1) ,'o','filled', 'markeredgecolor','k');
title('Jul')

hold on;


chl1 = chl_stat(9,5);
chl2 = chl_stat(21,5);
chl3 = chl_stat(33,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,6)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(8,:));
%cmin = min(big_nh4_mon(12,:));cmax = max(big_nh4_mon(12,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val(1,1) ,'o','filled', 'markeredgecolor','k');
title('Aug')

hold on;
subplot(4,3,7)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(9,:));
%cmin = min(big_nh4_mon(12,:));cmax = max(big_nh4_mon(12,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Sep')
hold on;



chl1 = chl_stat(10,5);
chl2 = chl_stat(22,5);
chl3 = chl_stat(34,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';



subplot(4,3,8)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(10,:));
%cmin = min(big_nh4_mon(12,:));cmax = max(big_nh4_mon(12,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Oct')
hold on;




chl1 = chl_stat(11,5);
chl2 = chl_stat(23,5);
chl3 = chl_stat(35,5);
chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(4,3,9)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(11,:));
%cmin = min(big_nh4_mon(12,:));cmax = max(big_nh4_mon(12,:));
caxis([0 0.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,8), lon_nodes(:,8),45, chl_val(2,1) ,'o','filled', 'markeredgecolor','k');
title('Nov')

hold on;

subplot(4,3,10)


chl1 = chl_stat(12,5);
chl_val = cell2mat(chl1);


tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh41(12,:));
%cmin = min(big_nh4_mon(12,:));cmax = max(big_nh4_mon(12,:));
caxis([0 0.1]);  
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');
title('Dec')


%% SPRING VERSUS SUMMER

chl_tbl = nut_stat(nut_stat{:,4} == 'Ammonia',:);
chl_stat = table2cell(chl_tbl);

chl1 = chl_stat(1,5);
chl2 = chl_stat(13,5);
chl3 = chl_stat(25,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


%close;
subplot(2,3,1);
set(gcf,'color','w');
set(gca, 'Fontsize',13)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh4_pt(3,:));
%cmin = min(mon_nh4_pt(3,1,:)); cmax = max(mon_nh4_pt(3,1,:));
%caxis([0 1.1]);
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,[7,9]), lon_nodes(:,[7,9]),45, chl_val([1,3],1) ,'o','filled', 'markeredgecolor','k');

hold on

chl1 = chl_stat(2,5);
chl2 = chl_stat(14,5);
chl3 = chl_stat(26,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(2,3,2);

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh4_pt(4,:));
%cmin = min(mon_nh4_pt(3,1,:)); cmax = max(mon_nh4_pt(3,1,:));
caxis([0 1.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_val ,'o','filled', 'markeredgecolor','k');

hold on;

subplot(2,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh4_pt(5,:));
%cmin = min(mon_nh4_pt(3,1,:)); cmax = max(mon_nh4_pt(3,1,:));
caxis([0 1.1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);

hold on;


chl1 = chl_stat(4,5);
chl2 = chl_stat(16,5);
chl3 = chl_stat(28,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';

subplot(2,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh4_pt(6,:));
%cmin = min(mon_nh4_pt(3,1,:)); cmax = max(mon_nh4_pt(3,1,:));
caxis([0 1.1]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val(1,1) ,'o','filled', 'markeredgecolor','k');


hold on


chl1 = chl_stat(5,5);
chl2 = chl_stat(17,5);
chl3 = chl_stat(29,5);

chl_val = [chl1 chl2 chl3];
chl_val = cell2mat(chl_val)';


subplot(2,3,5);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh4_pt(7,:));
%cmin = min(mon_nh4_pt(3,1,:)); cmax = max(mon_nh4_pt(3,1,:));
caxis([0 1.1]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,[8,9]), lon_nodes(:,[8,9]),45, chl_val([2,3],1) ,'o','filled', 'markeredgecolor','k');

hold on

chl1 = chl_stat(6,5);
chl_val = cell2mat(chl1);

subplot(2,3,6);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_nh4_pt(8,:));
%cmin = min(mon_nh4_pt(3,1,:)); cmax = max(mon_nh4_pt(3,1,:));
caxis([0 1.1]);
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_val ,'o','filled', 'markeredgecolor','k');

%% DO

%================== DO ========================================%

%% MONTHLY DO PLOTS
%%%% Monthly Averaged Surface DO

close; 
set(gcf,'color','w');
set(gca, 'Fontsize',13)
chl_p = phy_stat(1:3,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);



subplot(4,3,1);


tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(1,:));
%cmin = min(mon_do(1,:)); cmax = max(mon_do(1,:));
caxis([0 41]);
hbar = colorbar;
Cmap = cmocean('-haline');
colormap(Cmap);
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Jan')
hold on

chl_p = phy_stat(4:6,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

subplot(4,3,2);

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(2,:));
%cmin = min(mon_do(2,:));cmax = max(mon_do(2,:));
caxis([0 41]); 
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Feb')
hold on;


chl_p = phy_stat(7:9,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

subplot(4,3,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(3,:));
%cmin = min(mon_do(3,:));cmax = max(mon_do(3,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Mar')
hold on;


chl_p = phy_stat(10:12,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);


subplot(4,3,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(4,:));
%cmin = min(mon_do(4,:));cmax = max(mon_do(4,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Apr')

hold on


chl_p = phy_stat(13:15,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

subplot(4,3,5);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(5,:));
%cmin = min(mon_do(5,:));cmax = max(mon_do(5,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('May')
hold on

chl_p = phy_stat(15,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);


subplot(4,3,6);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(6,:));
%cmin = min(mon_do(6,:));cmax = max(mon_do(6,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_p ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lld_n(stations_mod(1:7),2),lld_n(stations_mod(1:7),3),50, do_exo(1:7,1),'o','filled','markeredgecolor','k');
title('Jun')
hold on;


subplot(4,3,7);


chl_p = phy_stat(17:24,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);



tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(7,:));
%cmin = min(mon_do(7,:));cmax = max(mon_do(7,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Jul')

hold on;



chl_p = phy_stat(25:32,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);

subplot(4,3,8)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(8,:));
%cmin = min(mon_do(8,:));cmax = max(mon_do(8,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,[1,3:9]), lon_nodes(:,[1,3:9]),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Aug')

hold on;
subplot(4,3,9)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(9,:));
%cmin = min(mon_do(9,:));cmax = max(mon_do(9,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title('Sep')
hold on;


chl_p = phy_stat(33:35,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);
chl_exo = s1_data_oct302018(1,2);

subplot(4,3,10)
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(10,:));
%cmin = min(mon_do(10,:));cmax = max(mon_do(10,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_p ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lat_nodesS1, lon_nodesS1,45, chl_exo ,'o','filled', 'markeredgecolor','k');
hold on;
scatter(lld_n(stations_mod(8:22),2),lld_n(stations_mod(8:22),3),50, do_exo(8:22,1),'o','filled','markeredgecolor','k');
title('Oct')
hold on;

chl_p = phy_stat(36:38,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);


subplot(4,3,11)

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(11,:));
%cmin = min(mon_do(11,:));cmax = max(mon_do(11,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7:9), lon_nodes(:,7:9),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Nov')

hold on;
subplot(4,3,12)


chl_p = phy_stat(39,6);
chl_p = table2cell(chl_p); chl_p = cell2mat(chl_p);



tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_do1(12,:));
%cmin = min(mon_do(12,:));cmax = max(mon_do(12,:));
caxis([0 41]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg O2/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
hold on;
scatter(lat_nodes(:,7), lon_nodes(:,7),45, chl_p ,'o','filled', 'markeredgecolor','k');
title('Dec')

%% POC
set(gcf,'color','w');
set(gca, 'Fontsize',13)
for ifile = 1:12
    
    subplot(4,3,ifile)
    
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_poc1(ifile,:));
%cmin = min(big_poc_mon(ifile,:));cmax = max(big_poc_mon(ifile,:));
caxis([0 1]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg C/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title(months(ifile))
end

% PON
set(gcf,'color','w');
set(gca, 'Fontsize',13)
for ifile = 1:12
    
    subplot(4,3,ifile)
    
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_pon1(ifile,:));
%cmin = min(big_pon_mon(ifile,:));cmax = max(big_pon_mon(ifile,:));
caxis([0 0.2]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title(months(ifile))
end 

% DON
set(gcf,'color','w');
set(gca, 'Fontsize',13)
for ifile = 1:12
    
    subplot(4,3,ifile)
    
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_don1(ifile,:));
%cmin = min(big_don_mon(ifile,:));cmax = max(big_don_mon(ifile,:));
caxis([0 0.8]); 
hold on;
Cmap = cmocean('-haline');
colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(mg N/l)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);
title(months(ifile))
end  

% KD
for ifile = 1:12
    
    subplot(4,3,ifile)
    
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), mon_kd(ifile,:));
%cmin = min(big_don_mon(ifile,:));cmax = max(big_don_mon(ifile,:));
%caxis([0]); 
hold on;
%Cmap = %cmocean('-haline');
%colormap(Cmap);
hbar = colorbar;
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on;
ylabel(hbar, '(m-1)','Fontsize',12);
xlabel('Longitude', 'Fontsize', 13); ylabel('Latitude', 'Fontsize', 13);

end

%% LINE PLOTS

close;

figure(3)
chl_h2 = LISALLnutrientssurface2018(LISALLnutrientssurface2018{:,7} == 'Chlorophyll a',:);
chl_2 = chl_h2(chl_h2{:,5} == 2,:);

chl_h2_t = datenum(chl_2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,1)
%plot(big_chla(:,342),'linewidth',1.5)
hold on;
plot(big_chla1(:,342), 'linewidth',1.5)
hold on;
plot(chl_h2_t,chl_2{:,8}, 'r.','markersize',20)
hold on;
%plot(chl_phy_t, chl_phy, 'r.', 'markersize',20)
ylabel('chlorophyll (ug chl/l)', 'fontsize',14)
%xlabel('days in 2018', 'fontsize', 14)
%legend('old run','new run','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on;

doc_h2 = LISALLnutrientssurface2018(LISALLnutrientssurface2018{:,7} == 'Dissolved Organic Carbon',:);
doc_2 = doc_h2(doc_h2{:,5} == 2,:);

doc_h2_t = datenum(doc_2{:,1})-datenum([2018,01,01,00,00,00]);
subplot(5,1,2)
%plot(big_doc(:,342), 'linewidth',1.5)
hold on;
plot(big_doc1(:,342), 'linewidth',1.5)
hold on;
plot(doc_h2_t,doc_2{:,8}, 'r.','markersize',20)
ylabel('DOC (mg C/l)', 'fontsize',14)
%xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on

no3_h2 = LISALLnutrientssurface2018(LISALLnutrientssurface2018{:,7} == 'Nitrate + Nitrite',:);
no3_2 = no3_h2(no3_h2{:,5} == 2,:);

no3_h2_t = datenum(no3_2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,3)
%plot(big_no3(:,342), 'linewidth',1.5)
hold on;
plot(big_no31(:,342), 'linewidth',1.5)
hold on;
plot(no3_h2_t,no3_2{:,8}, 'r.','markersize',20)
ylabel('nitrate (mg N/l)', 'fontsize',14)
%xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])

hold on;

nh4_h2 = LISALLnutrientssurface2018(LISALLnutrientssurface2018{:,7} == 'Ammonia',:);
nh4_2 = nh4_h2(nh4_h2{:,5} == 2,:);

nh4_h2_t = datenum(no3_2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,4)
%plot(big_nh4(:,342), 'linewidth',1.5)
hold on;
plot(big_nh41(:,342), 'linewidth',1.5)
hold on;
plot(nh4_h2_t,nh4_2{:,8}, 'r.','markersize',20)
ylabel('ammonium (mg N/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on


%do_phy = LISphyh2mean{:,2};
 %do_phy_t =  datenum(LISphyh2mean{:,1})-datenum([2018,01,01,00,00,00]);
 do_phy = surface1mphysics{:,2};
 do_phy_t =  datenum(surface1mphysics{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,5)
%plot(big_do(:,342), 'linewidth',1.5)
hold on;
plot(big_do1(:,342), 'linewidth',1.5)
hold on;
plot(do_phy_t,do_phy, 'r.','markersize',20)
ylabel('DO (mg O2/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
legend('new run','H2 CT DEEP', 'Orientation','horizontal', 'location','southoutside')
xlim([1 365])


close;
figure(4)

subplot(3,1,1)
%plot(big_pon1(:,342), 'linewidth',1.5)
%hold on;
plot(big_pon1(:,342), 'linewidth',1.5)
hold on;
ylabel('PON (mg N/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')

hold on


subplot(3,1,2)
%plot(big_poc1(:,342), 'linewidth',1.5)
%hold on;
plot(big_poc1(:,342), 'linewidth',1.5)
hold on;
ylabel('POC (mg C/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')

hold on


subplot(3,1,3)
%plot(big_don1(:,342), 'linewidth',1.5)
%hold on;
plot(big_don1(:,342), 'linewidth',1.5)
hold on;
ylabel('DON (mg N/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')

hold on

subplot(4,1,4)
%plot(big_kd1(:,342), 'linewidth',1.5)
%hold on;
plot(big_kd(:,342), 'linewidth',1.5)
hold on;
ylabel('KD (m-1)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')

%%
subplot(9,1,1)
plot((11).date, EE61_variables(11).value, '-k') %chla
hold on
plot(mean_EE61_variables(11).dates, mean_EE61_variables(11).values, 'or', 'linewidth',1)
ylabel('discharge (m3/s)')
hold on
subplot(9,1,2)
plot(EE61_variables(1).date, EE61_variables(1).value, '-k') %chla
hold on
plot(mean_EE61_variables(1).dates, mean_EE61_variables(1).values, 'or', 'linewidth',1)
ylabel('chla (ug chl/L)')
hold on
subplot(9,1,3)
plot(EE61_variables(2).date, EE61_variables(2).value, '-k') % DO
hold on
plot(mean_EE61_variables(2).dates, mean_EE61_variables(2).values, 'or','linewidth',1)
ylabel('DO2 (mg O2/l)')
hold on
subplot(9,1,4)
plot(EE61_variables(5).date, EE61_variables(5).value, '-k') % NH4
hold on
plot(mean_EE61_variables(5).dates, mean_EE61_variables(5).values, 'or','linewidth',1)
ylabel('NH4 (mg N/l)')
hold on
subplot(9,1,5)
plot(EE61_variables(6).date, EE61_variables(6).value, '-k') % NO3
hold on
plot(mean_EE61_variables(6).dates, mean_EE61_variables(6).values, 'or','linewidth',1)
ylabel('NO3 (mg N/l)')
hold on
subplot(9,1,6)
plot(EE61_variables(7).date, EE61_variables(7).value, '-k') % PO4
hold on
plot(mean_EE61_variables(7).dates, mean_EE61_variables(7).values, 'or','linewidth',1)
hold on
ylabel('PO4 (mg P/l)')
hold on
subplot(9,1,7)
plot(EE61_variables(4).date, EE61_variables(4).value, '-k') % DOP
hold on
plot(mean_EE61_variables(4).dates, mean_EE61_variables(4).values, 'or','linewidth',1)
ylim([0 0.15])
ylabel('DOP (mg P/l)')
hold on
subplot(9,1,8)
plot(EE61_variables(3).date, EE61_variables(3).value, '-k') % DON
hold on
plot(mean_EE61_variables(3).dates, mean_EE61_variables(3).values, 'or','linewidth',1)
ylabel('DON (mg N/l)')
hold on
subplot(9,1,9)
plot(EE61_variables(8).date, EE61_variables(8).value, '-k') % PON
hold on
plot(mean_EE61_variables(8).dates, mean_EE61_variables(8).values, 'or','linewidth',1)
ylabel('PON (mg N/l)')

legend('obs','mean','Orientation','horizontal','location','southoutside')
