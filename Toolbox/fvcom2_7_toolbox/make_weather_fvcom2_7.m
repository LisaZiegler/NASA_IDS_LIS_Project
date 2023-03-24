% write _mc_air.dat and _mc.dat forcings for FVCOM2.7
%first run the shell script to pull the files from the FTP server

%get_NARR.sh

%pull out the lat and long from the first file, all files will be the same

% find the closest NARR point for each of the grid cells in the model
% domain, and use the indices for the subsequet files

% open first netcdf file to get the grid dimensions
%below is the *.mat grid file created from sms 2dm file using sms2dm2mat
%load('/Volumes/Storage2/Ziegler_Lisa/run365/Tonic_26Jul.mat');
load('/Volumes/LaCie/tonicfvcom_casestudy/fvcom2.7fine_grid/meshgrid/tonic_grid16Sep.mat');
%get coordinates
%ncflag=input('Do you want to write out the FVCOM netcdf files for weather, 0 for no 1 for yes? ---> ');
%myyears_in = input('What year is this data for? --->');%;2016; 2017;2018];
myyears_in = 2018;
%find the closest NARR Grid to each element in the model
% if(exist('MyNARR_Zones.mat'));
%     load('MyNARR_Zones.mat');
% else
    disp('Finding the closest NARR zones');
    tic
    nc = netcdf([num2str(myyears_in(1)) '/air.2m.' num2str(myyears_in(1)) '.nc']);
    lat = nc{'lat'}(:);
    lon=nc{'lon'}(:);
    [ydim,xdim] = size(lon);
    %get one long vector for long and lat for the NARR grid;
    
    lat_long_out_bigrid(:,1)=reshape(lon,1,xdim*ydim);
    lat_long_out_bigrid(:,2)=reshape(lat,1,xdim*ydim);
    clear nc
    
    myzone_ele = dsearchn(lat_long_out_bigrid,[ll_e(:,1) ll_e(:,2)]);
    myzone_node= dsearchn(lat_long_out_bigrid,[lld_n(:,2) lld_n(:,3)]);
    
    myzone_ele = myzone_ele';
    myzone_node = myzone_node';
 
    toc
% end

disp('Found the closest NARR zones');

%%

% plot locations

% node_1 = lat_long_out_bigrid(71875,:);
% node_2 = lat_long_out_bigrid(71876,:);
% 
% values = [1:1:96673];
% v = strsplit(num2str(values))';
% %struc = gshhs('gshhs_f.b',[40.5 41.5],[-74 -72]);
% close;
% figure;
% 
% %geoshow(struc,'FaceColor',[1 1 1]*0.82);
% load coast;
% hold on;
% plot(long,lat);
% % plot(long+360,lat);
% hold on;
% grid on
% 
% for i = 1:96673;
% plot(lat_long_out_bigrid(i,1),lat_long_out_bigrid(i,2),'.k','markersize',10);
% text(lat_long_out_bigrid(i,1),lat_long_out_bigrid(i,2),v{i},'fontsize',12)
% end 
% hold on;
% ylim([36 41.5]); xlim([-74 -68]);
% 
% plot(lld_n(:,2),lld_n(:,3),'.r','markersize',5);
% hold on;
% plot(node_1(:,1),node_1(:,2),'*k','markersize',10)
% hold on;
% plot(node_2(:,1),node_2(:,2),'*k','markersize',10)

%%
% U WIND, V WIND, TEMPERATURE, REL. HUM, PRESSURE, DOWNWARD LWRAD, UPWARD
% LWRAD, DW SWRAD, UW SWRAD, EVAPORATION, PRECIPITATION,LATENT HEAT FLUX,
% SENSIBLE HEAT FLUX
% Open all of the files

myyear=num2str(myyears_in);

%  %see if the data exists, if not then can download it
% if(~exist('air.2m.2016.nc'))
%     myftp=ftp('ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/');
%     mget(myftp,['air.2m.' myyear '.nc']);
%
% end
%
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/dlwrf.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/dswrf.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/evap.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/lhtfl.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/prate.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/pres.sfc.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/shtfl.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/ulwrf.sfc.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/uswrf.sfc.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/uwnd.10m.$year.nc"
% wget "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Dailies/monolevel/vwnd.10m.$year.nc"


disp(['Running year ' myyear])

disp('Opening and extracting all of the data');
nc = netcdf([myyear '/uwnd.10m.' myyear '.nc']); % uwind (m/s)
uwind = nc{'uwnd'}(:);
clear nc
nc = netcdf([myyear '/vwnd.10m.' myyear '.nc']); % vwind (m/s)
vwind = nc{'vwnd'}(:);
clear nc
nc = netcdf([myyear '/air.2m.' myyear '.nc']); % air temperature  K
air = nc{'air'}(:);
clear nc
nc = netcdf([myyear '/rhum.2m.' myyear '.nc']); % relative humidity %
rhum = nc{'rhum'}(:);
clear nc
nc = netcdf([myyear '/pres.sfc.' myyear '.nc']); % surface pressure Pa
pres = nc{'pres'}(:);
clear nc
nc = netcdf([myyear '/dlwrf.' myyear '.nc']); % downward long wave radiation flux (W m^2)
dlwrf = nc{'dlwrf'}(:);
clear nc
nc = netcdf([myyear '/ulwrf.sfc.' myyear '.nc']); % upward long wave radiation flux  (W m^2)
ulwrf = nc{'ulwrf'}(:);
clear nc
nc = netcdf([myyear '/dswrf.' myyear '.nc']); % downward short wave flux (W m^2)
dswrf = nc{'dswrf'}(:);
clear nc
nc = netcdf([myyear '/uswrf.sfc.' myyear '.nc']); % upward short wave flux (W m^2)
uswrf = nc{'uswrf'}(:);
clear nc
nc = netcdf([myyear '/evap.' myyear '.nc']); % evaporation mm
evap = nc{'evap'}(:);
clear nc
nc = netcdf([myyear '/prate.' myyear '.nc']); % precipitation mm
prate = nc{'prate'}(:);
clear nc
nc = netcdf([myyear '/lhtfl.' myyear '.nc']); % latent heat flux
lhtfl = nc{'lhtfl'}(:);
clear nc
nc = netcdf([myyear '/shtfl.' myyear '.nc']); % sensible heat flux
shtfl = nc{'shtfl'}(:);
nc = netcdf([myyear '/tcdc.' myyear '.nc']); % total cloud cover
tcdc = nc{'tcdc'}(:);
nc = netcdf([myyear '/shum.2m.' myyear '.nc']); % specific humidity
shum = nc{'shum'}(:);
clear nc
disp('Extracted all the data, now putting into some big vectors');
[ntimes,~,~]=size(rhum);
%%

% need to find the coordinates in the grid that match the zone for the Chesapeake Bay
% dont want to loop through the entire grid for each time point, it takes
% too long

% now reshape into a 3 dimensional vector with the 9 met variables at each
% location in a column vector, similar to the met zones found previously

uwind_vector = reshape(uwind,ntimes,xdim*ydim)';
vwind_vector = reshape(vwind,ntimes,xdim*ydim)';
air_vector   = reshape(air,ntimes,xdim*ydim)';
rhum_vector   = reshape(rhum,ntimes,xdim*ydim)';
pres_vector   = reshape(pres,ntimes,xdim*ydim)';
dlwrf_vector   = reshape(dlwrf,ntimes,xdim*ydim)';
ulwrf_vector   = reshape(ulwrf,ntimes,xdim*ydim)';
dswrf_vector   = reshape(dswrf,ntimes,xdim*ydim)';
uswrf_vector   = reshape(uswrf,ntimes,xdim*ydim)';
evap_vector   = reshape(evap,ntimes,xdim*ydim)';
prate_vector   = reshape(prate,ntimes,xdim*ydim)';
lhtfl_vector   = reshape(lhtfl,ntimes,xdim*ydim)';
shtfl_vector   = reshape(shtfl,ntimes,xdim*ydim)';
%tcdc_vector   = reshape(tcdc,ntimes,xdim*ydim)';
shum_vector   = reshape(shum,ntimes,xdim*ydim)';

%correct bad values
uwind_vector(uwind_vector < -9E30)=0;
vwind_vector(vwind_vector < -9E30)=0;
air_vector(air_vector < -9E30)=0;
rhum_vector(rhum_vector < -9E30)=0;
pres_vector(pres_vector < -9E30)=0;
dlwrf_vector(dlwrf_vector < -9E30)=0;
ulwrf_vector(ulwrf_vector < -9E30)=0;
dswrf_vector(dswrf_vector < -9E30)=0;
uswrf_vector(uswrf_vector < -9E30)=0;
evap_vector(evap_vector < -9E30)=0;
prate_vector(prate_vector < -9E30)=0;
lhtfl_vector(lhtfl_vector < -9E30)=0;
shtfl_vector(shtfl_vector < -9E30)=0;
%tcdc_vector(tcdc_vector < -9E30)=0;
shum_vector(shum_vector < -9E30)=0;


uwind_mygrid=uwind_vector(myzone_node,:);
vwind_mygrid=vwind_vector(myzone_node,:);
%downward long wave radiation (W m^-2)
dlwrf_mygrid=dlwrf_vector(myzone_node,:);
%upward long wave radiation (W m^-2)
ulwrf_mygrid=ulwrf_vector(myzone_node,:);
%downward short wave radiation (W m^-2)
dswrf_mygrid=dswrf_vector(myzone_node,:);
%upward short wave radiation (W m^-2)
uswrf_mygrid=uswrf_vector(myzone_node,:);
%accumulated evaporation (convert kg m^-2 3hr^-1 --> m s^-1), negative is a
%loss of water
evap_mygrid=evap_vector(myzone_node,:)./3./1000./3600*-1;
%precipation rate (convert from kg m^-2 s^-1 to m s^-1);
prate_mygrid=prate_vector(myzone_node,:)./1000;
%latent heat flux (W m^-2
lhtfl_mygrid=lhtfl_vector(myzone_node,:);
%sensible heat flux (W m^-2)
shtfl_mygrid=shtfl_vector(myzone_node,:);
%air temperature, Celsius
air_mygrid=air_vector(myzone_node,:)-273.15;
%Pressure
pres_mygrid=pres_vector(myzone_node,:);
%Relative Humidity
rhum_mygrid=rhum_vector(myzone_node,:);
%specific Humidity
shum_mygrid=shum_vector(myzone_node,:);
%total cloud cover
% tcdc_mygrid=tcdc_vector(myzone_node,:);

% calculate net heat flux
% SWDOWN,GLW ---- lw and sw downward radiaton W/m^2
%SHTFL,LHTFL ----- sensible and latent heat flux W/m^2
%LWUP,SWUP -----short and long wave upward rad W/m^2
%PRATE ---- precipitation rate (+/-) kg/m^2/s
% ================================

net_Shortwave=dswrf_mygrid-uswrf_mygrid;
net_Longwave=dlwrf_mygrid-ulwrf_mygrid;
net_Heatflux=net_Shortwave+net_Longwave+lhtfl_mygrid+shtfl_mygrid;
disp('Assigned the data to the grid, now need to load into a data structure');

%%
%loop to calculate the wind velocity for points
%calculates the wind velocity and direction for all points at each time step

%following is used to convert from atan2 to Dir met
%
% % % Dirmet = atan2(-Umet,-Vmet) * DperR = 270 - ( atan2(Vmet,Umet) * DperR )
% from http://www.eol.ucar.edu/content/wind-direction-quick-reference
my_unique_zones = unique(myzone_node,'stable');%get the uinique
%my_unique_zones = my_unique_zones(1,1); 

for i= 1: length(my_unique_zones); % fill a big array
  % 1st dimensino is variable, second dimension is time, 3rd dimension is
  % zone
     met_array_out(:,1,i)=uwind_mygrid(my_unique_zones(i),:);
     met_array_out(:,2,i)=vwind_mygrid(my_unique_zones(i),:);
     met_array_out(:,3,i)=air_mygrid(my_unique_zones(i),:);
     met_array_out(:,4,i)=rhum_mygrid(my_unique_zones(i),:);
     met_array_out(:,5,i)=pres_mygrid(my_unique_zones(i),:);
     met_array_out(:,6,i)=dlwrf_mygrid(my_unique_zones(i),:);
     met_array_out(:,7,i)=ulwrf_mygrid(my_unique_zones(i),:);
     met_array_out(:,8,i)=dswrf_mygrid(my_unique_zones(i),:);
     met_array_out(:,9,i)=uswrf_mygrid(my_unique_zones(i),:);
     met_array_out(:,10,i)=evap_mygrid(my_unique_zones(i),:);
     met_array_out(:,11,i)=prate_mygrid(my_unique_zones(i),:);
     met_array_out(:,12,i)=lhtfl_mygrid(my_unique_zones(i),:);
     met_array_out(:,13,i)=shtfl_mygrid(my_unique_zones(i),:);              
end

wind_vel=[];
wind_direction_met=[];
for t = 1 : ntimes
                               
for igrid = 1:9441;

         wind_vel(igrid,t)=sqrt(abs(uwind_mygrid(igrid,t)^2+vwind_mygrid(igrid,t)^2)) ;
           %find the wind direction in radians 
           
         wind_direction_met(igrid,t)= 270-(atan2(vwind_mygrid(igrid,t),uwind_mygrid(igrid,t)).*180/pi) ; %-180 to 180 range
%          wind_direction(igrid,t)= atan2(v(igrid,t),u(igrid,t)).*180/pi ; %-180 to 180 range, meteorlogical;

        if wind_direction_met(igrid,t) > 360
            
            wind_direction_met(igrid,t) = wind_direction_met(igrid,t) - 360;
            
        end
        
    end
       
end


% ========================
% THE order of our met array is

% U WIND, V WIND, TEMPERATURE, REL. HUM, PRESSURE, DOWNWARD LWRAD, UPWARD
% LWRAD, DW SWRAD, UW SWRAD, EVAPORATION, PRECIPITATION,LATENT HEAT FLUX,
% SENSIBLE HEAT FLUX

% now write the zone files with the met data
 if(~isdir('Zone_files'));
    
    mkdir('Zone_files')
    
end
   for i = 1 :length(my_unique_zones);
     
     if(my_unique_zones(i)<1000)
         
       filename = ['Zone_files/hst_wind_non_uniform_00', num2str(my_unique_zones(i)),'.dat'];
     else
         filename = ['Zone_files/hst_wind_non_uniform_', num2str(my_unique_zones(i)),'.dat'];
     end
           fid = fopen(filename,'w');
           
           fprintf(fid,'%s\n','Time Variable Meteorlogical Data');
           

           
    
    
       for t = 1 : time;
         
            fprintf(fid,'%f\n',(t-1)*3);
         
            fprintf(fid,'%f\t %f\t %f\t %f\t %f\t %f\n', prate_mygrid(i,t)', evap_mygrid(i,t)',...
              wind_vel(i,t),wind_direction_met(igrid,t), Q_net(i,t)',...
              shortwave(i,t));
        
        end
       
       fclose(fid);
     
%  

%      if (typeflag == 0 )
   if(my_unique_zones(i)<1000)
       fid2=fopen(['Zone_files/hst_hfx_non_uniform_00', num2str(my_unique_zones(i)),'.dat'],'w');
    else
      fid2=fopen(['Zone_files/hst_hfx_non_uniform_', num2str(my_unique_zones(i)),'.dat'],'w');
    end

%    fid2 = fopen('tonic_mc_air.dat','w')
       fprintf(fid2,'%f\n',time);
    
       for t =1:time ;
        
         fprintf(fid2,'%f\n',(t-1)*3);
        
          fprintf(fid2,'%f\t %f\t %f\t %f\t %f\t %f\t %f\n',t10(i,t), RH(i,t),...
              Pressure(i,t), dwlw(i,t), dwsw(i,t),...
              net_sw(i,t),net_lw(i,t)) ;
            
       end
    
     fclose(fid2)
    
     disp(' met file is created');

%      end
  
   end
       
    disp('congratulations, now run prep_meteorology.f90 to get the non-uniform forcing file')



