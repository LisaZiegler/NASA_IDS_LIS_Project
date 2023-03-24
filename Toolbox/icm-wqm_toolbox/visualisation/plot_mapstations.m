%%%%%%%%%% Match model points with observational points
%% Search for latitude and longitude points located near the closest node 

% This is if you just want to make a quick plot of where your stations are
% and you do not feel like putting your data into a structure

stations_node = {'st18','st19','st21','st22','st23','st27','F2','H2','H4'};
lat_station = [41.12,41.05,41.16,41.08,41.14,41.15,41.08,41.17,41.10];
lon_station = [-73.09,-73.08,-73.01,-73.02,-72.94,-72.84,-73.16,-72.96,-72.93];

stations_mod_hr = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_station ; lat_station])');


% Plot grid 
close all;

set(gcf,'color','w');
set(gca, 'Fontsize',14);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n);
cmocean('thermal');
colorbar off
 hold on
% Plot stations 
 for i = 1 : length(stations_mod_hr);
    plot(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),'r.','Markersize',20)
    text(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),stations_node{i},'fontsize',18)   
 end
 
 export_fig('Allstations.png','zbuffer','-r600');
%% 
% Chla
lat_f2 = [41.08];
lon_f2 = [-73.16];
st_hr = {'F2'};

chla_hr = [3.9950];

stations_mod_hr = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_f2 ; lat_f2])');
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n); colorbar off
 hold on
 
 for i = 1 : length(stations_mod_hr);
    plot(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),'ro')
    text(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),st_hr{i},'fontsize',24)
    
 end
 
lat_nodes=[];lon_nodes=[];
for i = 1 : length(stations_mod_hr);
    lat_node = mean(lld_n(stations_mod_hr(i),2));
    lon_node = mean(lld_n(stations_mod_hr(i),3));
    lat_nodes(i)=lat_node';
    lon_nodes(i)=lon_node';
end 
DO_surf = squeeze(chla_mod(:,1,:));
DO1 = mean(DO_surf(1:31,:));
DO2 = mean(DO_surf(32:59,:));
DO3 = mean(DO_surf(60:89,:));
DO4 = mean(DO_surf(90:120,:));
DO5 = mean(DO_surf(121:151,:));
DO6 = mean(DO_surf(152:181,:));
DO7 = mean(DO_surf(182:212,:));
DO8 = mean(DO_surf(213:243,:));
DO9 = mean(DO_surf(244:273,:));
DO10 = mean(DO_surf(274:304,:));
DO11 = mean(DO_surf(305:334,:));
DO12 = mean(DO_surf(335:end,:));

DO_monthly = {DO1, DO2, DO3, DO4, DO5, DO6,DO7, DO8, DO9, DO10, DO11, DO12};
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for ifile = 11;

    monthly = months{ifile};

figure(1001)
set(gcf,'color','w');
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',DO_monthly{ifile},'facecolor','interp','edgecolor','interp');
caxis([0 3])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '(ug chl/l)','Fontsize',13);
title(['monthly averaged Surface Chlorophyll --' monthly], 'Fontsize',14);
hold on

scatter(lat_nodes, lon_nodes,30, chla_hr,'o','filled', 'markeredgecolor','k')
drawnow

export_fig(sprintf('chlanovcomparison_000%d', ifile),'-r600', '-png');

end

% DOC
lat_hr = [41.2272,41.192817,41.18495,41.179817,41.171617];
lon_hr = [-73.110517,-73.113467,-73.112,-73.117217,-73.1113];

st_hr = {'HR1','HR2','HR3','HR4','HR5'};
doc_hr = [2.12,2.22,2.03,1.96,2.21];
 
lat_samples = [41.16751,41.14,41.13,41.111,41.0959,41.172,41.1342,41.1591];
lon_samples = [-73.106,-73.092,-73.0953,-73.0957,-73.09081,-73.11355,-73.063,-73.09288];
st_hr = {'S1','S2','S3','S4','S5','S6','S7','S8','S9'};

doc_hr = [2.65,2.02,1.61,1.49,1.51,2.51,1.489,1.479];

stations_mod_hr = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_samples ; lat_samples])');
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n); colorbar off
 hold on
 
 for i = 1 : length(stations_mod_hr);
    plot(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),'ro')
    text(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),st_hr{i},'fontsize',24)
    
 end
 
lat_nodes=[];lon_nodes=[];
for i = 1 : length(stations_mod_hr);
    lat_node = mean(lld_n(stations_mod_hr(i),2));
    lon_node = mean(lld_n(stations_mod_hr(i),3));
    lat_nodes(i)=lat_node';
    lon_nodes(i)=lon_node';
end 
DO_surf = squeeze(DOC_mod(:,1,:));
DO1 = mean(DO_surf(1:31,:));
DO2 = mean(DO_surf(32:59,:));
DO3 = mean(DO_surf(60:89,:));
DO4 = mean(DO_surf(90:120,:));
DO5 = mean(DO_surf(121:151,:));
DO6 = mean(DO_surf(152:181,:));
DO7 = mean(DO_surf(182:212,:));
DO8 = mean(DO_surf(213:243,:));
DO9 = mean(DO_surf(244:273,:));
DO10 = mean(DO_surf(274:304,:));
DO11 = mean(DO_surf(305:334,:));
DO12 = mean(DO_surf(335:end,:));

DO_monthly = {DO1, DO2, DO3, DO4, DO5, DO6,DO7, DO8, DO9, DO10, DO11, DO12};
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for ifile = 10;

    monthly = months{ifile};

figure(1001)
set(gcf,'color','w');
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',DO_monthly{ifile},'facecolor','interp','edgecolor','interp');
caxis([0 2])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '(mg C l^-^1)','Fontsize',13);
title(['monthly averaged Surface DOC --' monthly], 'Fontsize',14);
hold on

scatter(lat_nodes, lon_nodes,30, doc_hr,'o','filled', 'markeredgecolor','k')
drawnow

export_fig(sprintf('DOCoctobercomparison_000%d', ifile),'-r600', '-png');

end



 