load('/Users/lisaziegler/Downloads/CT_DEEP/housatonic_river_exo2_data.mat');
load('/Volumes/LaCie/Housatonic_BGC/tonic_grid16Sep.mat');
load('/Volumes/LaCie/fvcom_runs/LIS_allstations_2018.mat');
readtable('/Volumes/LaCie/fvcom_runs/LIS_phy_F2.csv');

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

marsh_file = shaperead('/Volumes/LaCie/CT_DEM/meshgrid/wetlands/WheelerMarsh.shp');
marsh2 = marsh_file;
lat_big1=[];
lon_big1=[];
for i = 1:length(marsh2)
    lon = marsh2(i).X';
    lat = marsh2(i).Y';
    lon_big1 = [lon_big1; lon];
    lat_big1 = [lat_big1;lat];
end


utm = repmat('18 T', [2387 1]);

[Y,X] = utm2deg(lon_big1,lat_big1,utm);

stations_mod = dsearchn([lld_n(:,2) lld_n(:,3)],([stations.lon ; stations.lat])'); %find the model nodes that are close to the wqm station points

 tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n);
 hold on
 
 for i = 1 : length(stations_mod);
    plot(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),'ro')
    text(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),stations(i).name,'fontsize',24)
    
 end
 
 saveas(gcf,'AllStations_map','png');
 
%% 

hr_table = struct2table(hr_exo);
hr_name = hr_table.station_description;
hr_lon = hr_table.lon';
hr_lat = hr_table.lat';
unique_stations = unique(hr_name);
hr_time = hr_table.datetime;
hr_time = datestr(hr_time, 'yyyy-mm-dd');
stations_mod = dsearchn([lld_n(:,2) lld_n(:,3)],([hr_lon ; hr_lat])');


% Extract EXO data from structure
temp_exo = zeros(22,1);
salt_exo = zeros(22,1);


for i = 1:22;

temp_mean = mean(hr_table.exo{i,1}.temp_C);
salt_mean = mean(hr_table.exo{i,1}.sal_psu);

temp_exo(i,1) = temp_mean;
salt_exo(i,1) = salt_mean;

end

set(gca,'color','w','fontsize',16)
% patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',salt_daily(178,:),'facecolor','interp','edgecolor','interp');
% hold on
plot(lon_big,lat_big,'-k','linewidth',1)
hold on 
plot(X, Y,'-k','linewidth',1);
hold on
xlim([-73.18 -72.9])
ylim([41.06 41.4])
%colormap(matlab_colormap);
colormap(jet(15))
caxis([0 15])
hbar = colorbar('YColor',[0 0 0]);
% Now create a label for each tick mark (you can modify these however you want)
ylabel(hbar,'Chlorophyll (ug/l)','fontsize',13)
hold on
scatter(lld_n(stations_mod(2:4),2),lld_n(stations_mod(2:4),3),50, squeeze(big_temp(4142,1,2:4)),'o','filled', 'markeredgecolor','k');
hold on
scatter(lld_n(stations_mod(2:4),2),lld_n(stations_mod(2:4),3),60, temp_exo(2:4,1),'o','filled', 'markeredgecolor','k');





