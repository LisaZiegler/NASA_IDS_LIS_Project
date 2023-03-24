% Create timeseries plot of bgc model outputs and validate with observation
%

readtable('/Volumes/LaCie/Housatonic_BGC/Data/LIS_Nutrients_matlab/LIS_ALLnutrients_surface2018.csv');
readtable('/Volumes/LaCie/Housatonic_BGC/Data/LIS_Nutrients_matlab/LIS_ALLphysics_surface2018.csv');
readtable('/Volumes/LaCie/Housatonic_BGC/Data/surface_1m_physics.csv');
load('/Users/lisaziegler/Downloads/CT_DEEP/housatonic_river_exo2_data.mat');
load('/Volumes/LaCie/Housatonic_BGC/tonic_grid16Sep.mat');
load('/Volumes/LaCie/Housatonic_BGC/Model_Runs/run_Jul821/daily_bgc_Jul821.mat');
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

lat_lis = [41.12,41.05,41.16,41.08,41.14,41.15,41.08,41.17,41.10];
lon_lis = [-73.09,-73.08,-73.01,-73.02,-72.94,-72.84,-73.16,-72.96,-72.93];
lisnames = {'18','19','21','22','23','27','F2','H2','H4'};
stations_mod_lis = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_lis ; lat_lis])');

hr_table = struct2table(hr_exo);
hr_name = hr_table.station_description;
hr_lon = hr_table.lon';
hr_lat = hr_table.lat';
unique_stations = unique(hr_name);
hr_time = hr_table.datetime;
hr_time = datestr(hr_time, 'yyyy-mm-dd');
stations_mod = dsearchn([lld_n(:,2) lld_n(:,3)],([hr_lon ; hr_lat])');

% Plot only stations that have data collected in 2018
for i = 1
close;

set(gcf,'color','w'); 

figure(i)
set(gca,'color','w','fontsize',14);
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4]) 

hold on ;

for ii = 1 : 22;
plot(lld_n(stations_mod(ii),2),lld_n(stations_mod(ii),3),'r.','markersize',25)
text(lld_n(stations_mod(ii),2),lld_n(stations_mod(ii),3),hr_name(ii),'fontsize',12)
title('2018', 'fontweight','normal')
end

hold on

for j = 1:length(stations_mod_lis)
plot(lld_n(stations_mod_lis(j),2),lld_n(stations_mod_lis(j),3),'r.','markersize',25)
text(lld_n(stations_mod_lis(j),2),lld_n(stations_mod_lis(j),3),lisnames(j),'fontsize',12)
end

end 

% Extract EXO data from structure
temp_exo = zeros(22,1);
salt_exo = zeros(22,1);
do_exo = zeros(22,1);
chl_exo = zeros(22,1);

for i = 1:22;

do_mean = mean(hr_table.exo{i,1}.do_mgL);
chl_mean = mean(hr_table.exo{i,1}.chl_ugL);
temp_mean = mean(hr_table.exo{i,1}.temp_C);
salt_mean = mean(hr_table.exo{i,1}.sal_psu);

do_exo(i,1) = do_mean;
chl_exo(i,1) = chl_mean;
temp_exo(i,1) = temp_mean;
salt_exo(i,1) = salt_mean;

end

%%
close;
set(gca, 'fontsize',13);
figure(4)
chl_h2 = LISALLnutrientssurface2018(LISALLnutrientssurface2018{:,7} == 'Chlorophyll a',:);
chl_2 = chl_h2(chl_h2{:,5} == 2,:);

chl_h2_t = datenum(chl_2{:,1})-datenum([2018,01,01,00,00,00]);

%chl_h2 = nutsurfF2(nutsurfF2{:,2} == 'Chlorophyll a',:);

%chl_h2_t = datenum(chl_2{:,1})-datenum([2018,01,01,00,00,00]);

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
%xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on


 do_phy = surface1mphysics{:,2};
 do_phy_t =  datenum(surface1mphysics{:,1})-datenum([2018,01,01,00,00,00]);
 %do_phy = surface1mphysicsH4{:,4};
 %do_phy_t =  datenum(surface1mphysicsH4{:,1})-datenum([2018,01,01,00,00,00]);

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


%%

close;

figure

chl_h2 = nutsurfF2(nutsurfF2{:,2} == 'Chlorophyll a',:);

chl_h2_t = datenum(chl_2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,1)
%plot(big_chla(:,342),'linewidth',1.5)
hold on;
plot(big_chla1(:,7), 'linewidth',1.5)
hold on;
plot(chl_h2_t,chl_2{:,8}, 'r.','markersize',20)
hold on;
%plot(chl_phy_t, chl_phy, 'r.', 'markersize',20)
ylabel('chlorophyll (ug chl/l)', 'fontsize',14)
%xlabel('days in 2018', 'fontsize', 14)
%legend('old run','new run','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on;

doc_h2 = nutsurfF2(nutsurfF2{:,2} == 'Dissolved Organic Carbon',:);


doc_h2_t = datenum(doc_2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,2)
%plot(big_doc(:,342), 'linewidth',1.5)
hold on;
plot(big_doc1(:,7), 'linewidth',1.5)
hold on;
plot(doc_h2_t,doc_2{:,8}, 'r.','markersize',20)
ylabel('DOC (mg C/l)', 'fontsize',14)
%xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on

no3_h2 = nutsurfF2(nutsurfF2{:,2} == 'Nitrate + Nitrite',:);

no3_h2_t = datenum(no3_2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,3)
%plot(big_no3(:,342), 'linewidth',1.5)
hold on;
plot(big_no31(:,7), 'linewidth',1.5)
hold on;
plot(no3_h2_t,no3_2{:,8}, 'r.','markersize',20)
ylabel('nitrate (mg N/l)', 'fontsize',14)
%xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])

hold on;

nh4_h2 = nutsurfF2(nutsurfF2{:,2} == 'Ammonia',:);

nh4_h2_t = datenum(no3_2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,4)
%plot(big_nh4(:,342), 'linewidth',1.5)
hold on;
plot(big_nh41(:,7), 'linewidth',1.5)
hold on;
plot(nh4_h2_t,nh4_2{:,8}, 'r.','markersize',20)
ylabel('ammonium (mg N/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
%legend('MODEL w/rivON','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on


 %do_phy = LISphyh2mean{:,2};
 %do_phy_t =  datenum(LISphyh2mean{:,1})-datenum([2018,01,01,00,00,00]);
 do_phy = surface1mphysicsF2{:,4};
 do_phy_t =  datenum(surface1mphysicsF2{:,1})-datenum([2018,01,01,00,00,00]);

subplot(5,1,5)
%plot(big_do(:,342), 'linewidth',1.5)
hold on;
plot(big_do1(:,7), 'linewidth',1.5)
hold on;
plot(do_phy_t,do_phy, 'r.','markersize',20)
ylabel('DO (mg O2/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
legend('new run','F2 CT DEEP', 'Orientation','horizontal', 'location','southoutside')
xlim([1 365])

%%
figure

chl_h2 = LISphy18mean{:,3};
chl_h2_t =  datenum(LISphy18mean{:,1})-datenum([2018,01,01,00,00,00]);

hr = datenum(hr_time)-datenum([2018,01,01,00,00,00]);

subplot(2,1,1)
%plot(big_chla(:,342),'linewidth',1.5)
hold on;
plot(big_chla1(:,456), 'k','linewidth',1.5)
hold on;
plot(chl_h2_t,chl_h2, 'r.','markersize',20)
 hold on;
plot(hr(12,1), chl_exo(12,1), 'b.', 'markersize',20)
hold on;
plot(hr(17,1), chl_exo(17,1), 'b.', 'markersize',20)
ylabel('chlorophyll (ug chl/l)', 'fontsize',14)
%xlabel('days in 2018', 'fontsize', 14)
%legend('old run','new run','H2 CT DEEP', 'Orientation','horizontal')
xlim([1 365])
hold on;

do_phy = LISphy18mean{:,2};
do_phy_t =  datenum(LISphy18mean{:,1})-datenum([2018,01,01,00,00,00]);

subplot(2,1,2)
%plot(big_do(:,342), 'linewidth',1.5)
hold on;
plot(big_do1(:,456), 'k','linewidth',1.5)
hold on;;
plot(do_phy_t,do_phy, 'r.','markersize',20)
hold on
plot(hr(12,1),do_exo(12,1), 'b.','markersize',20)
hold on;
plot(hr(17,1),do_exo(17,1), 'b.','markersize',20)
ylabel('DO (mg O2/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
legend('new run','18 CT DEEP', 'Sample5','Orientation','horizontal', 'location','southoutside')
xlim([1 365])


%%
subplot(3,1,1)
plot(big_chla1(:,5865), 'k','linewidth',1.5)
hold on;
plot(hr(8,1), chl_exo(8,1), 'b.', 'markersize',20)
hold on;
plot(hr(13,1), chl_exo(13,1), 'b.', 'markersize',20)
ylabel('chlorophyll (ug chl/l)', 'fontsize',14)
xlabel('days in 2018', 'fontsize', 14)
legend('new run','Sample1', 'Orientation','horizontal')
xlim([1 365])


%%
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

