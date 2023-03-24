y%%%%%%% ICM post processing
%%%%% L. Zielger, 2020, HPL

% This script is used to read in ICM outputs and validate it against
% observational data. 

% load grid file
load ('/Users/lisaziegler/Desktop/BGCoutputs/Housatonicwqm_grid.mat');

kb=11;
for k = 1 : kb;
    s_hybrid(k,1:length(xyd_n))=-((k-1)/(kb-1))^1.5;
end

for i = 2:11;
   layer_thickness(i-1,:)=mean(s_hybrid(i-1:i,:),1)*-1;
end
% 
% % Pull out the data and depth integrate to get mass per unit are
% % for each layer (allows for a better average)
% fname = 'tonic_0001.nc';
idir = dir('*.nc');
numfile = length(idir);
nameseg = 'tonic';

DO_mod1 = [];
alg1_mod1 =[];
alg2_mod1 = [];
NH4_mod1  = [];
NO3_mod1  = [];
lpoc_mod1 = [];
rpoc_mod1 = [];
lpon_mod1 = [];
rpon_mod1 = [];
wcdoc1_1 = [];
wcdoc2_1 = [];
wcdoc3_1 = [];
wcncdoc1_1 = [];
wcncdoc2_1 = [];
wcncdoc3_1 = [];
wc_cdon1_1 = [];
wc_cdon2_1 = [];
wc_cdon3_1 = [];
wc_ncdon1_1 = [];
wc_ncdon2_1 = [];
wc_ncdon3_1 = [];
kd1 = [];
mod_time1= [];
depth1_1 = [];
temp = [];
salt = [];
chla_mod1  = [];
PON_mod1 = [];
POC_mod1 = [];
DOC_mod1 = [];
DON_mod1 = [];
alg1_spring = [];
alg2_summer = [];
for idir = 1 : numfile;
    
    numfile_count = [1:numfile];
    
    kstr = num2str(numfile_count(idir),'%04d');
    
    fname = [nameseg,'_',kstr,'.nc'];
    
    nc = netcdf([fname]);
    
DO_mod= squeeze(nc{'DOXG'}(:,1,:));
alg1_mod = squeeze(nc{'B1'}(:,1,:));
alg2_mod = squeeze(nc{'B2'}(:,1,:));
NH4_mod  = squeeze(nc{'NH4'}(:,1,:));
NO3_mod  = squeeze(nc{'NO3'}(:,1,:));
lpoc_mod = squeeze(nc{'LPOC'}(:,1,:));
rpoc_mod = squeeze(nc{'RPOC'}(:,1,:));
lpon_mod = squeeze(nc{'LPON'}(:,1,:));
rpon_mod = squeeze(nc{'RPON'}(:,1,:));
wcdoc1 = squeeze(nc{'WC_CDOC1'}(:,1,:));
wcdoc2 = squeeze(nc{'WC_CDOC2'}(:,1,:));
wcdoc3 = squeeze(nc{'WC_CDOC3'}(:,1,:));
wcncdoc1 = squeeze(nc{'WC_NCDOC1'}(:,1,:));
wcncdoc2 = squeeze(nc{'WC_NCDOC2'}(:,1,:));
wcncdoc3 = squeeze(nc{'WC_NCDOC3'}(:,1,:));
wc_cdon1 = squeeze(nc{'WC_CDON1'}(:,1,:));
wc_cdon2 = squeeze(nc{'WC_CDON2'}(:,1,:));
wc_cdon3 = squeeze(nc{'WC_CDON3'}(:,1,:));
wc_ncdon1 = squeeze(nc{'WC_NCDON1'}(:,1,:));
wc_ncdon2 = squeeze(nc{'WC_NCDON2'}(:,1,:));
wc_ncdon3 = squeeze(nc{'WC_NCDON3'}(:,1,:));
kd = squeeze(nc{'KD'}(:,1,:));
%mod_time= nc{'time'}(:)./24/86400;
%depth1 = nc{'depth'}(:);
temp_in = squeeze(nc{'temp'}(:,1,:));
salt_in = squeeze(nc{'salinity'}(:,1,:));

    clear nc
    
PON_in = lpon_mod + rpon_mod;
POC_in = lpoc_mod + rpoc_mod;

alg1_mod_spr = (alg1_mod./50).*1000;
alg2_mod_sum = (alg2_mod./50).*1000;

alg1_spring = [alg1_spring; alg1_mod_spr];
alg2_summer = [alg2_summer; alg2_mod_sum];

chla_mod = (alg1_mod./50+alg2_mod./50).*1000;   % convert mg C/l to ug chl/l based on chl:C ratio from model              
PON_mod =PON_in + (alg1_mod+alg2_mod)./5.68; % Use Redfield to convert algae C to N
POC_mod = POC_in + alg1_mod+alg2_mod;
DOC_mod = (wcdoc1+wcdoc2+wcdoc3) + (wcncdoc1+wcncdoc2+wcncdoc3);
DON_mod = (wc_cdon1 + wc_cdon2 + wc_cdon3) + (wc_ncdon1 + wc_ncdon2 + wc_ncdon3);

    DO_mod1 = [DO_mod1; DO_mod];
    %alg1_mod1 =[alg1_mod1; alg1_mod];
    %alg2_mod1 = [alg2_mod1; alg2_mod];
    chla_mod1 = [chla_mod1; chla_mod];
    NH4_mod1  = [NH4_mod1; NH4_mod];
    NO3_mod1  = [NO3_mod1; NO3_mod];
    %lpoc_mod1 = [lpoc_mod1; lpoc_mod];
    %rpoc_mod1 = [rpoc_mod1; rpoc_mod];
    %lpon_mod1 = [lpon_mod1; lpon_mod];
    %rpon_mod1 = [rpon_mod1; rpon_mod];
    %wcdoc1_1 = [wcdoc1_1; wcdoc1];
    %wcdoc2_1 = [wcdoc2_1; wcdoc2];
    %wcdoc3_1 = [wcdoc3_1; wcdoc3];
    %wcncdoc1_1 = [wcncdoc1_1; wcncdoc1];
    %wcncdoc2_1 = [wcncdoc2_1; wcncdoc2];
    %wcncdoc3_1 = [wcncdoc3_1; wcncdoc3];
    %wc_cdon1_1 = [wc_cdon1_1; wc_cdon1];
    %wc_cdon2_1 = [wc_cdon2_1; wc_cdon2];
    %wc_cdon3_1 = [wc_cdon3_1; wc_cdon3];
    %wc_ncdon1_1 = [wc_ncdon1_1; wc_ncdon1];
    %wc_ncdon2_1 = [wc_ncdon2_1; wc_ncdon2];
    %wc_ncdon3_1 = [wc_ncdon3_1; wc_ncdon3];
    PON_mod1 = [PON_mod1; PON_mod];
    POC_mod1 = [POC_mod1; POC_mod];
temp = [temp; temp_in];
salt = [salt; salt_in];

%kd1 = [kd1; kd];
%mod_time1= [mod_time1; mod_time];
%depth1_1 = [depth1_1; depth1];



DOC_mod1 = [DOC_mod1; DOC_mod];
DON_mod1 = [DON_mod1; DON_mod];
    
end

% DO_mod = squeeze(nc{'DOXG'}(:,1,:));
% alg1_mod = squeeze(nc{'B1'}(:,1,:));
% alg2_mod = squeeze(nc{'B2'}(:,1,:));
% NH4_mod  = squeeze(nc{'NH4'}(:,1,:));
% NO3_mod  = squeeze(nc{'NO3'}(:,1,:));
% lpoc_mod = squeeze(nc{'LPOC'}(:,1,:));
% rpoc_mod = squeeze(nc{'RPOC'}(:,1,:));
% lpon_mod = squeeze(nc{'LPON'}(:,1,:));
% rpon_mod = squeeze(nc{'RPON'}(:,1,:));
% wcdoc1 = squeeze(nc{'WC_CDOC1'}(:,1,:));
% wcdoc2 = squeeze(nc{'WC_CDOC2'}(:,1,:));
% wcdoc3 = squeeze(nc{'WC_CDOC3'}(:,1,:));
% wcncdoc1 = squeeze(nc{'WC_NCDOC1'}(:,1,:));
% wcncdoc2 = squeeze(nc{'WC_NCDOC2'}(:,1,:));
% wcncdoc3 = squeeze(nc{'WC_NCDOC3'}(:,1,:));
% wc_cdon1 = squeeze(nc{'WC_CDON1'}(:,1,:));
% wc_cdon2 = squeeze(nc{'WC_CDON2'}(:,1,:));
% wc_cdon3 = squeeze(nc{'WC_CDON3'}(:,1,:));
% wc_ncdon1 = squeeze(nc{'WC_NCDON1'}(:,1,:));
% wc_ncdon2 = squeeze(nc{'WC_NCDON2'}(:,1,:));
% wc_ncdon3 = squeeze(nc{'WC_NCDON3'}(:,1,:));
% kd = squeeze(nc{'KD'}(:,1,:));
% mod_time= nc{'time'}(:)./86400;
% depth1 = nc{'depth'}(:);


% DO_mod = nc_varget('tonic_0001.nc', 'DOXG');
% alg1_mod = nc_varget('tonic_0001.nc', 'B1');
% alg2_mod = nc_varget('tonic_0001.nc', 'B2');
% NH4_mod  = nc_varget('tonic_0001.nc', 'NH4');
% NO3_mod  = nc_varget('tonic_0001.nc', 'NO3');
% lpoc_mod = nc_varget('tonic_0001.nc', 'LPOC');
% rpoc_mod = nc_varget('tonic_0001.nc', 'RPOC');
% lpon_mod = nc_varget('tonic_0001.nc', 'LPON');
% rpon_mod = nc_varget('tonic_0001.nc', 'RPON');
% wcdoc1 = nc_varget('tonic_0001.nc', 'WC_CDOC1');
% wcdoc2 = nc_varget('tonic_0001.nc', 'WC_CDOC2');
% wcdoc3 = nc_varget('tonic_0001.nc', 'WC_CDOC3');
% wcncdoc1 = nc_varget('tonic_0001.nc', 'WC_NCDOC1');
% wcncdoc2 = nc_varget('tonic_0001.nc', 'WC_NCDOC2');
% wcncdoc3 = nc_varget('tonic_0001.nc', 'WC_NCDOC3');
% wc_cdon1 = nc_varget('tonic_0001.nc', 'WC_CDON1');
% wc_cdon2 = nc_varget('tonic_0001.nc', 'WC_CDON2');
% wc_cdon3 = nc_varget('tonic_0001.nc', 'WC_CDON3');
% wc_ncdon1 = nc_varget('tonic_0001.nc', 'WC_NCDON1');
% wc_ncdon2 = nc_varget('tonic_0001.nc', 'WC_NCDON2');
% wc_ncdon3 = nc_varget('tonic_0001.nc', 'WC_NCDON3');
% kd = nc_varget('tonic_0001.nc','KD');
% 
% mod_time= nc_varget('tonic_0001.nc','time')./86400;
% depth1 = nc_varget('tonic_0001.nc','depth');
% PON_mod = [];
% POC_mod = [];
% PON_in = lpon_mod + rpon_mod;
% POC_in = lpoc_mod + rpoc_mod;
% 
% chla_mod = (alg1_mod./50+alg2_mod./50).*1000;   % convert mg C/l to ug chl/l based on chl:C ratio from model              
% PON_mod =PON_in + (alg1_mod+alg2_mod)./5.68; % Use Redfield to convert algae C to N
% POC_mod = POC_in + alg1_mod+alg2_mod;
% DOC_mod = (wcdoc1+wcdoc2+wcdoc3) + (wcncdoc1+wcncdoc2+wcncdoc3);
% DON_mod = (wc_cdon1 + wc_cdon2 + wc_cdon3) + (wc_ncdon1 + wc_ncdon2 + wc_ncdon3);

big_chla1 = zeros(365,7430);
big_nh41 = zeros(365,7430);
big_no31 = zeros(365,7430);
big_pon1 = zeros(365,7430);
big_poc1 = zeros(365,7430);
big_doc1 = zeros(365,7430);
big_don1 = zeros(365,7430);
big_do1 = zeros(365,7430);
big_temp = zeros(365,7430);
big_salt = zeros(365,7430);
% big_chla1 = zeros(203,7430);
% big_nh41 = zeros(203,7430);
% big_no31 = zeros(203,7430);
% big_pon1 = zeros(203,7430);
% big_poc1 = zeros(203,7430);
% big_doc1 = zeros(203,7430);
% big_don1 = zeros(203,7430);
% big_do1 = zeros(203,7430);
% big_temp = zeros(203,7430);
% big_salt = zeros(203,7430);
%big_depth1 = zeros(365,7430);
%big_time1 = zeros(365,1);
%big_kd1 = zeros(365,7430);
big_spring_alg = zeros(365,7430);
big_summer_alg = zeros(365,7430);
    for ii = 1:365;
    if ii*24<8758
    big_chla1(ii,:) = mean(chla_mod1((ii-1)*24+1:ii*24,:));
    big_nh41(ii,:) = mean(NH4_mod1((ii-1)*24+1:ii*24,:));
    big_no31(ii,:) = mean(NO3_mod1((ii-1)*24+1:ii*24,:));
    big_pon1(ii,:) = mean(PON_mod1((ii-1)*24+1:ii*24,:));
    big_poc1(ii,:) = mean(POC_mod1((ii-1)*24+1:ii*24,:));
    big_doc1(ii,:) = mean(DOC_mod1((ii-1)*24+1:ii*24,:));
    big_don1(ii,:) = mean(DON_mod1((ii-1)*24+1:ii*24,:));
    big_do1(ii,:) = mean(DO_mod1((ii-1)*24+1:ii*24,:));
    big_temp(ii,:) = mean(temp((ii-1)*24+1:ii*24,:));
    big_salt(ii,:) = mean(salt((ii-1)*24+1:ii*24,:));
    big_spring_alg(ii,:) = mean(alg1_spring((ii-1)*24+1:ii*24,:));
    big_summer_alg(ii,:) = mean(alg2_summer((ii-1)*24+1:ii*24,:));
    %big_depth1(ii,:) = mean(depth1_1((ii-1)*24+1:ii*24,:));
    %big_time1(ii,1) = mean(mod_time1((ii-1)*24+1:ii*24,1));
    %big_kd1(ii,:) = mean(kd1((ii-1)*24+1:ii*24,:));
    else
    big_chla1(ii,:) = mean(chla_mod1((ii-1)*24+1:end,:));
    big_nh41(ii,:) = mean(NH4_mod1((ii-1)*24+1:end,:));
    big_no31(ii,:) = mean(NO3_mod1((ii-1)*24+1:end,:));
    big_pon1(ii,:) = mean(PON_mod1((ii-1)*24+1:end,:));
    big_poc1(ii,:) = mean(POC_mod1((ii-1)*24+1:end,:));
    big_doc1(ii,:) = mean(DOC_mod1((ii-1)*24+1:end,:));
    big_don1(ii,:) = mean(DON_mod1((ii-1)*24+1:end,:));
    big_do1(ii,:) = mean(DO_mod1((ii-1)*24+1:end,:));
    big_temp(ii,:) = mean(temp((ii-1)*24+1:end,:));
    big_salt(ii,:) = mean(salt((ii-1)*24+1:end,:));
    big_spring_alg(ii,:) = mean(alg1_spring((ii-1)*24+1:end,:));
    big_summer_alg(ii,:) = mean(alg2_summer((ii-1)*24+1:end,:));
    %big_depth1(ii,:) = mean(depth1_1((ii-1)*24+1:end,:));
    %big_time1(ii,1) = mean(mod_time1((ii-1)*24+1:end,1));
    %big_kd1(ii,:) = mean(kd1((ii-1)*24+1:end,:));
    end
    end
    
%save('daily_bgc_rivOn.mat', 'big_do', 'big_chla', 'big_nh4', 'big_no3', 'big_doc', 'big_poc','big_pon','mod_time', 'big_time', 'big_depth', 'depth1','big_kd');    
save('daily_bgc_tempfunc.mat', 'big_do1', 'big_chla1', 'big_nh41', 'big_no31', 'big_doc1', 'big_poc1','big_pon1', 'big_don1', 'big_temp', 'big_salt');%'mod_time1', 'big_time1', 'big_depth1', 'depth1','big_kd1', 'big_poc1','big_pon1', 'big_don1'); 

d1 = datetime(2018,01,01);
d2 = datetime(2018,31,12);
date = (d1:d2)';
date1 = datevec(date);

% mon_doc1 = zeros(12,7430);
% mon_do1 = zeros(12,7430);
% mon_chla1 = zeros(12,7430);
% mon_nh41 = zeros(12,7430);
% mon_no31 = zeros(12,7430);
% mon_kd1 = zeros(12,7430);
% mon_poc1 = zeros(12,7430);
% mon_pon1 = zeros(12,7430);
% mon_don1 = zeros(12,7430);
% mon_temp = zeros(12, 7430);
% mon_salt = zeros(12, 7430);

mon_doc1 = zeros(12,7430);
mon_do1 = zeros(12,7430);
mon_chla1 = zeros(12,7430);
mon_nh41 = zeros(12,7430);
mon_no31 = zeros(12,7430);
mon_kd1 = zeros(12,7430);
mon_poc1 = zeros(12,7430);
mon_pon1 = zeros(12,7430);
mon_don1 = zeros(12,7430);
mon_temp = zeros(12, 7430);
mon_salt = zeros(12, 7430);
 for i=2018                                    %Years Looping
for j=1:12                                      %Months Looping
idx=find((date1(:,1)==i) & (date1(:,2)==j));  %Find all indeces for the specific month
meandata1=mean(big_doc1(idx,:));
meandata2=mean(big_do1(idx,:));
meandata3=mean(big_chla1(idx,:));
meandata4=mean(big_nh41(idx,:));
meandata5=mean(big_no31(idx,:));
%meandata7 = mean(big_kd1(idx,:));
meandata8 = mean(big_pon1(idx,:));
meandata9 = mean(big_poc1(idx,:));
meandata10 = mean(big_don1(idx,:));
meandata11 = mean(big_temp(idx,:));
meandata12 = mean(big_salt(idx,:));
mon_doc1(j,:)=meandata1;          %Store in a new matrix the mean values per month
mon_do1(j,:)=meandata2;
mon_chla1(j,:)=meandata3;
mon_nh41(j,:)=meandata4;
mon_no31(j,:)=meandata5;
%mon_kd1(j,:)=meandata7;
mon_pon1(j,:)=meandata8;
mon_poc1(j,:)=meandata9;
mon_don1(j,:)=meandata10;
mon_temp(j,:)=meandata11;
mon_salt(j,:)=meandata12;
end
 end
%save('monthly_bgc_nonptOn.mat', 'mon_do', 'mon_chla', 'mon_nh4', 'mon_no3', 'mon_doc','mon_kd'); 
save('monthly_bgc_tempfunc.mat', 'mon_do1', 'mon_chla1', 'mon_nh41', 'mon_no31', 'mon_doc1', 'mon_poc1','mon_pon1','mon_don1','mon_temp','mon_salt');    

% big_chla = zeros(365,1,7478);
% big_nh4 = zeros(365,1,7478);
% big_no3 = zeros(365,1,7478);
% big_pon = zeros(365,1,7478);
% big_poc = zeros(365,1,7478);
% big_doc = zeros(365,1,7478);
% big_don = zeros(365,1,7478);
% big_do = zeros(365,1,7478);
% big_depth = zeros(365,1,7478);
% big_time = zeros(365,1);
% big_kd = zeros(365,1,7478);
% 
%     for ii = 1:365;
%     if ii*24<8758
%     big_chla(ii,:,:) = mean(chla_mod((ii-1)*24+1:ii*24,:,:));
%     big_nh4(ii,:,:) = mean(NH4_mod((ii-1)*24+1:ii*24,:,:));
%     big_no3(ii,:,:) = mean(NO3_mod((ii-1)*24+1:ii*24,:,:));
%     big_pon(ii,:,:) = mean(PON_mod((ii-1)*24+1:ii*24,:,:));
%     big_poc(ii,:,:) = mean(POC_mod((ii-1)*24+1:ii*24,:,:));
%     big_doc(ii,:,:) = mean(DOC_mod((ii-1)*24+1:ii*24,:,:));
%     big_don(ii,:,:) = mean(DON_mod((ii-1)*24+1:ii*24,:,:));
%     big_do(ii,:,:) = mean(DO_mod((ii-1)*24+1:ii*24,:,:));
%     big_depth(ii,:,:) = mean(depth1((ii-1)*24+1:ii*24,:,:));
%     big_time(ii,1) = mean(mod_time((ii-1)*24+1:ii*24,1));
%     big_kd(ii,:,:) = mean(kd((ii-1)*24+1:ii*24,1));
%     else
%     big_chla(ii,:,:) = mean(chla_mod((ii-1)*24+1:end,:,:));
%     big_nh4(ii,:,:) = mean(NH4_mod((ii-1)*24+1:end,:,:));
%     big_no3(ii,:,:) = mean(NO3_mod((ii-1)*24+1:end,:,:));
%     big_pon(ii,:,:) = mean(PON_mod((ii-1)*24+1:end,:,:));
%     big_poc(ii,:,:) = mean(POC_mod((ii-1)*24+1:end,:,:));
%     big_doc(ii,:,:) = mean(DOC_mod((ii-1)*24+1:end,:,:));
%     big_don(ii,:,:) = mean(DON_mod((ii-1)*24+1:end,:,:));
%     big_do(ii,:,:) = mean(DO_mod((ii-1)*24+1:end,:,:));
%     big_depth(ii,:,:) = mean(depth1((ii-1)*24+1:end,:,:));
%     big_time(ii,1) = mean(mod_time((ii-1)*24+1:end,1));
%     big_kd(ii,:,:) = mean(kd((ii-1)*24+1:end,1));
%     end
%     end
%     
% save('daily_bgc.mat', 'big_do', 'big_chla', 'big_nh4', 'big_no3', 'big_doc', 'big_poc','big_pon''mod_time', 'big_time', 'big_depth', 'depth1','big_kd');    
% 
% 
% d1 = datetime(2018,01,01);
% d2 = datetime(2018,31,12);
% date = (d1:d2)';
% date1 = datevec(date);
% 
% mon_doc = zeros(12,1,7478);
% mon_do = zeros(12,1,7478);
% mon_chla = zeros(12,1,7478);
% mon_nh4 = zeros(12,1,7478);
% mon_no3 = zeros(12,1,7478);
% mon_kd = zeros(12,1,7478);
% mon_poc = zeros(12,1,7478);
% mon_pon = zeros(12,1,7478);
% mon_cpb = zeros(12,1);
% 
%  for i=2018                                    %Years Looping
% for j=1:12                                      %Months Looping
% idx=find((date1(:,1)==i) & (date1(:,2)==j));  %Find all indeces for the specific month
% meandata1=mean(big_doc(idx,:,:));
% meandata2=mean(big_do(idx,:,:));
% meandata3=mean(big_chla(idx,:,:));
% meandata4=mean(big_nh4(idx,:,:));
% meandata5=mean(big_no3(idx,:,:));
% meandata7 = mean(big_kd(idx,:,:));
% meandata8 = mean(big_pon(idx,:,:));
% meandata9 = mean(big_poc(idx,:,:));
% mon_doc(j,:,:)=meandata1;          %Store in a new matrix the mean values per month
% mon_do(j,:,:)=meandata2;
% mon_chla(j,:,:)=meandata3;
% mon_nh4(j,:,:)=meandata4;
% mon_no3(j,:,:)=meandata5;
% mon_kd(j,:,:)=meandata7;
% mon_pon(j,:,:)=meandata8;
% mon_poc(j,:,:)=meandata9;
% end
%  end
% 
% save('monthly_bgc.mat', 'mon_do', 'mon_chla', 'mon_nh4', 'mon_no3', 'mon_doc', 'mon_poc','mon_pon','mon_kd');    


first_day = 1;
starting_date=first_day; 
data_int=ICM_int; %1/24 -- hourly output

mod_time=mod_time+1;

[ntimes,nlayers,nstations]=size(chla_mod);

for iz=1:nlayers;
    
    depth(:,iz,:)=depth1.*repmat(layer_thickness(iz),ntimes,1);

end

%%
day = 1:365;
day = day';
rounded_mod_times=day;
mod_times=round(mod_time);
%%
stations_mod = dsearchn([lld_n(:,2) lld_n(:,3)],([stations.lon ; stations.lat])'); %find the model nodes that are close to the wqm station points

day = 1:351;
day = day';
rounded_mod_times=day;
% search window for calculating statistics
data_int = 1;
Swindow=3;
close all;

% Chla
for ifile = 1:351;

days = day(ifile);

figure(1000)
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',squeeze(chla_mod(ifile,1,:)),'facecolor','interp','edgecolor','interp');
caxis([0 3])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '[ug chl/l]','Fontsize',13);
title(['Surface Chlorophyll on day' num2str(days)], 'Fontsize',14);

drawnow

export_fig(sprintf('chl_000%d', ifile),'-a4','-q100', '-jpeg');

end

% monthly values
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

for ifile = 10:12;

    monthly = months{ifile};

figure(1001)
set(gcf,'color','w');
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',DO_monthly{ifile},'facecolor','interp','edgecolor','interp');
caxis([0 2])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '(ug chl/l)','Fontsize',13);
title(['monthly averaged Surface Chlorophyll --' monthly], 'Fontsize',14);

drawnow

export_fig(sprintf('fallmonthly_000%d', ifile),'-r600', '-png');

end

% DO
for ifile = 1:351;

days = day(ifile);

figure(1000)
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',squeeze(DO_mod(ifile,1,:)),'facecolor','interp','edgecolor','interp');
caxis([0 15])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, 'O_2 (mg O_2 l^-^1)','Fontsize',13);
title(['Surface DO on day' num2str(days)], 'Fontsize',14);

drawnow

export_fig(sprintf('do_000%d', ifile),'-a4','-q100', '-jpeg');

end
% monthly values
DO_surf = squeeze(DO_mod(:,1,:));
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

for ifile = 10:12;

    monthly = months{ifile};

figure(1001)
set(gcf,'color','w');
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',DO_monthly{ifile},'facecolor','interp','edgecolor','interp');
caxis([0 15])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, 'O_2 (mg O_2 l^-^1)','Fontsize',13);
title(['monthly averaged Surface DO --' monthly], 'Fontsize',14);

drawnow

export_fig(sprintf('fallmonthly_000%d', ifile),'-r600', '-png');

end
%NO3
% monthly values
DO_surf = squeeze(NO3_mod(:,1,:));
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

for ifile = 10:12;

    monthly = months{ifile};

figure(1001)
set(gcf,'color','w');
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',DO_monthly{ifile},'facecolor','interp','edgecolor','interp');
caxis([0 0.2])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '(mg N l^-^1)','Fontsize',13);
title(['monthly averaged Surface NO_3^2^- --' monthly], 'Fontsize',14);

drawnow

export_fig(sprintf('fallmonthly_000%d', ifile),'-r600', '-png');


end

%NH4
% monthly values
DO_surf = squeeze(NH4_mod(:,1,:));
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

for ifile = 10:12;

    monthly = months{ifile};

figure(1001)
set(gcf,'color','w');
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',DO_monthly{ifile},'facecolor','interp','edgecolor','interp');
caxis([0 0.1])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '(mg N l^-^1)','Fontsize',13);
title(['monthly averaged Surface NH_4^+ --' monthly], 'Fontsize',14);

drawnow

export_fig(sprintf('fallmonthly_000%d', ifile),'-r600', '-png');

end

%DOC
% monthly values
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

for ifile = 10:12;

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

%scatter(lat_nodes, lon_nodes,30, doc_hr,'o','filled', 'markeredgecolor','k')
drawnow

export_fig(sprintf('fallmonthly_000%d', ifile),'-r600', '-png');

end
lat_hr = [41.2272,41.192817,41.18495,41.179817,41.171617];
lon_hr = [-73.110517,-73.113467,-73.112,-73.117217,-73.1113];
st_hr = ['HR1','HR2','HR3','HR4','HR5'];
st_hr = {'HR1','HR2','HR3','HR4','HR5'};
doc_hr = [2.12,2.22,2.03,1.96,2.21];
open FVCOM_plotsalt1
date = 178;
stations_mod_hr = dsearchn([lld_n(:,2) lld_n(:,3)],([lon_hr ; lat_hr])');
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), -h_n);
 hold on
 
 for i = 1 : length(stations_mod_hr);
    plot(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),'ro')
    text(lld_n(stations_mod_hr(i),2),lld_n(stations_mod_hr(i),3),st_hr{i},'fontsize',24)
    
 end
 lat_nodes=[];lon_nodes=[];for i = 1 : length(stations_mod_hr);lat_node = mean(lld_n(stations_mod_hr(i),2));lon_node = mean(lld_n(stations_mod_hr(i),3));lat_nodes(i)=lat_node';lon_nodes(i)=lon_node';end
close all;
 for ifile = 27;

days = day(ifile);

figure(1003)
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',squeeze(DOC_mod(ifile,1,:)),'facecolor','interp','edgecolor','interp');
caxis([0 2])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '(mg C l^-^1)','Fontsize',13);
title(['Surface DO on day' num2str(days)], 'Fontsize',14);

hold on

scatter(lat_nodes, lon_nodes,30, doc_hr,'o','fill')

drawnow

%export_fig(sprintf('doc_000%d', ifile),'-a4','-q100', '-jpeg');

end

% DON: DON (mg N l^-^1)'

% monthly values
DO_surf = squeeze(DON_mod(:,1,:));
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

for ifile = 10:12;

    monthly = months{ifile};

figure(1001)
set(gcf,'color','w');
patch('Vertices',lld_n(:,2:3),'Faces',e2n(:,2:4),'Cdata',DO_monthly{ifile},'facecolor','interp','edgecolor','interp');
caxis([0 0.2])
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hbar = colorbar;
ylabel(hbar, '(mg N l^-^1)','Fontsize',13);
title(['monthly averaged Surface DON --' monthly], 'Fontsize',14);

drawnow

export_fig(sprintf('fallmonthly_000%d', ifile),'-r600', '-png');

end


for istat=1:length(stations_mod);
    
    myvar_id=1; % change this ID based on where the variable is located within the CBP data structure 
                % generated by read_cbp.m
        
    statstring = ['station' num2str(istat)];
    mkdir(statstring)
   
    figure
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),-h_n);
    hold on
    plot(lld_n(stations_mod(istat),2),lld_n(stations_mod(istat),3),'*','color','k','markersize',24);

    saveas(gcf,[statstring '/station_map'],'png');
    
  figure  

  % get the depth averaged model output at our station
  for i = 1:351;
  depth_averaged = squeeze(mean(chla_mod(i,:,stations_mod(istat))));
  mean_chla(i,:) = depth_averaged;
  end
  p1 = plot(day,mean_chla,'d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
         
  hold on
  %depth_averaged=squeeze(sum(chla_mod(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
    %p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
            %hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date); 
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      
    % now for each unique date, find corresponding values and take the
    
      % average
      mean_values=[];
      
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate));
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

      
    p2 = plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end
h = [p1;p2];
      title([statstring ' Chl a'],'Fontsize',20);
      %legend(h,'Modeled','Observed');
      xlabel('Days in 2018','Fontsize',20);ylabel('Chl a (ug l^-^1)','Fontsize',20)  
      saveas(gcf,[statstring '/chla'],'png');
       saveas(gcf,[statstring '/chla'],'eps');     
       saveas(gcf,[statstring '/chla'],'fig'); 
end
       figure;
%  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
       end
       
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    
axis([axis_min axis_max axis_min axis_max ]); 

plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed chla (ug chla l^-^1)');ylabel('Modeled chla (ug chla l^-^1)');
title([ 'Depth Averaged chla ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/chla_stats'],'png');
       saveas(gcf,[statstring '/chla_stats'],'eps');     
       saveas(gcf,[statstring '/chla_stats'],'fig');  
end        
%
% Dissolved Oxygen
for istat=1:length(stations_mod);
myvar_id=2;
  figure  
  
    for i = 1:351;
  depth_averaged = squeeze(mean(DO_mod(i,:,stations_mod(istat))));
  mean_chla(i,:) = depth_averaged;
  end
  plot(day,mean_chla,'d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
         
  hold on
  
  % get the depth averaged model output at our station
  %depth_averaged=squeeze(sum(DO_mod(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
   % p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
            hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates]=unique(stations(istat).variables(myvar_id).date); %ic refernces where the dates came from
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      

      % now for each unique date, find corresponding values and take the
      % average
      mean_values=[];
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate))
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

            % plot the real world measured data
     plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end

      %legend([p1,p2],'Modeled','Observed');
      xlabel('Days in 2018','Fontsize',20);ylabel('O_2 (mg O_2 l^-^1)','Fontsize',20)  
      saveas(gcf,[statstring '/DO'],'png');
      saveas(gcf,[statstring '/DO'],'eps');
       saveas(gcf,[statstring '/DO'],'fig'); 
end
%  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
       end
       
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;

axis([axis_min axis_max axis_min axis_max ]); 

p4=plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed O_2 (mg O_2 l^-^1)');ylabel('Modeled O_2 (mg O_2 l^-^1)');
title([ 'Depth Averaged O_2^- ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/DO_stats'],'png');
       saveas(gcf,[statstring '/DO_stats'],'eps');     
       saveas(gcf,[statstring '/DO_stats'],'fig');  
% NITRATE
for istat=1:length(stations_mod); figure  
  
 myvar_id=8;
 
  for i = 1:351;
  depth_averaged = squeeze(mean(NO3_mod(i,:,stations_mod(istat))));
  mean_chla(i,:) = depth_averaged;
  end
  plot(day,mean_chla,'d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
         
  hold on
  % get the depth averaged model output at our station
  %depth_averaged=squeeze(sum(NO3_mod(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
    %p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
           % hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates]=unique(stations(istat).variables(myvar_id).date); %ic refernces where the dates came from
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      

        
      % now for each unique date, find corresponding values and take the
      % average
      mean_values=[];
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate))
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

    plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end

      %legend([p1,p2],'Modeled','Observed');
      xlabel('Days in 2018','Fontsize',20);ylabel('NO_3^2^- (mg N l^-^1)','Fontsize',20)  
      saveas(gcf,[statstring '/NO3'],'png');
      saveas(gcf,[statstring '/NO3'],'eps');
       saveas(gcf,[statstring '/NO3'],'fig');
       end
       %  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
       end
       
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    
axis([axis_min axis_max axis_min axis_max ]); 

p4=plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed NO_3^- (mg N l^-^1)');ylabel('Modeled NO_3^- (mg N l^-^1)');
title([ 'Depth Averaged NO_3^- ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/NO3_stats'],'png');
       saveas(gcf,[statstring '/NO3_stats'],'eps');     
       saveas(gcf,[statstring '/NO3_stats'],'fig');  
% Ammonium
for istat=1:length(stations_mod)
    close all;
 figure  
  
 myvar_id=6;
  for i = 1:351;
  depth_averaged = squeeze(mean(NH4_mod(i,:,stations_mod(istat))));
  mean_chla(i,:) = depth_averaged;
  end
  plot(day,mean_chla,'d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
         
  hold on
  % get the depth averaged model output at our station
  %depth_averaged=squeeze(sum(NH4_mod(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
   % p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
   %         hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates]=unique(stations(istat).variables(myvar_id).date); %ic refernces where the dates came from
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      
     
      % now for each unique date, find corresponding values and take the
      % average
      mean_values=[];
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate))
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

     plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end

     % legend([p1,p2],'Modeled','Observed');
      xlabel('Days in 2018','Fontsize',20);ylabel('NH_4^+ (mg N l^-^1)','Fontsize',20)  
      saveas(gcf,[statstring '/NH4'],'png');
      saveas(gcf,[statstring '/NH4'],'eps');
       saveas(gcf,[statstring '/NH4'],'fig'); 
end      
       %  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
       end
       
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    
axis([axis_min axis_max axis_min axis_max ]); 

p4=plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed NH_4^+ (mg N l^-^1)');ylabel('Modeled NH_4^+ (mg N l^-^1)');
title([ 'Depth Averaged NH_4^+ ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/NH4_stats'],'png');
       saveas(gcf,[statstring '/NH4_stats'],'eps');     
       saveas(gcf,[statstring '/NH4_stats'],'fig');         
% Phosphate

 figure  
  
 myvar_id=9;
  % get the depth averaged model output at our station
  depth_averaged=squeeze(sum(PO4_mod(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
    p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
            hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates]=unique(stations(istat).variables(myvar_id).date); %ic refernces where the dates came from
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      

        
      % now for each unique date, find corresponding values and take the
      % average
      mean_values=[];
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate))
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

     p2=  plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end

      legend([p1,p2],'Modeled','Observed');
      xlabel('Days in 2005','Fontsize',20);ylabel('PO_4^3^- (mg P l^-^1)','Fontsize',20)  
      saveas(gcf,[statstring '/PO4'],'png');
      saveas(gcf,[statstring '/PO4'],'eps');
       saveas(gcf,[statstring '/PO4'],'fig');
       %  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
       end
       
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;

axis([axis_min axis_max axis_min axis_max ]); 

p4=plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed PO_4^3+ (mg P l^-^1)');ylabel('Modeled PO_4^3+ (mg P l^-^1)');
title([ 'Depth Averaged PO_4^3+ ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/PO4_stats'],'png');
       saveas(gcf,[statstring '/PO4_stats'],'eps');     
       saveas(gcf,[statstring '/PO4_stats'],'fig');             
       
% DON
for istat=1:length(stations_mod)
 figure  
  
 myvar_id=4;
 for i = 1:351;
  depth_averaged = squeeze(mean(DON_mod(i,:,stations_mod(istat))));
  mean_chla(i,:) = depth_averaged;
  end
  plot(day,mean_chla,'d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
         
  hold on
  
  % get the depth averaged model output at our station
  %depth_averaged=squeeze(sum(DON(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
   % p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
           % hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates]=unique(stations(istat).variables(myvar_id).date); %ic refernces where the dates came from
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      

        
      % now for each unique date, find corresponding values and take the
      % average
      mean_values=[];
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate))
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

     plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end

      %legend([p1,p2],'Modeled','Observed');
      xlabel('Days in 2018','Fontsize',20);ylabel('DON (mg N l^-^1)','Fontsize',20)  
      saveas(gcf,[statstring '/DON'],'png');
      saveas(gcf,[statstring '/DON'],'eps');
       saveas(gcf,[statstring '/DON'],'fig'); 
end 
       %  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
       end
       
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    
axis([axis_min axis_max axis_min axis_max ]); 

p4=plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed DON (mg N l^-^1)');ylabel('Modeled DON (mg N l^-^1)');
title([ 'Depth Averaged DON ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/DON_stats'],'png');
       saveas(gcf,[statstring '/DON_stats'],'eps');     
       saveas(gcf,[statstring '/DON_stats'],'fig');           

 % DOC
 for istat=1:length(stations_mod)
  figure  
  
 myvar_id=3;
 for i = 1:351;
  depth_averaged = squeeze(mean(DOC_mod(i,:,stations_mod(istat))));
  mean_chla(i,:) = depth_averaged;
  end
  plot(day,mean_chla,'d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
         
  hold on
  % get the depth averaged model output at our station
 % depth_averaged=squeeze(sum(DOC(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
   % p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
   %         hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates]=unique(stations(istat).variables(myvar_id).date); %ic refernces where the dates came from
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      

        
      % now for each unique date, find corresponding values and take the
      % average
      mean_values=[];
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate))
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

      plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end

      %legend([p1,p2],'Modeled','Observed');
      xlabel('Days in 2018','Fontsize',20);ylabel('DOC (mg C l^-^1)','Fontsize',20)  
      saveas(gcf,[statstring '/DOC'],'png');
      saveas(gcf,[statstring '/DOC'],'eps');
       saveas(gcf,[statstring '/DOC'],'fig'); 
 end       
       
        %  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
           
       end
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;

axis([axis_min axis_max axis_min axis_max ]); 

p4=plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed DOC (mg C l^-^1)');ylabel('Modeled DOC (mg C l^-^1)');
title([ 'Depth Averaged DOC ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/DOC_stats'],'png');
       saveas(gcf,[statstring '/DOC_stats'],'eps');     
       saveas(gcf,[statstring '/DOC_stats'],'fig');    
       
       % Salinity
   for istat=1:length(stations_mod)
  figure  
  
 myvar_id=11;
  for i = 1:351;
  depth_averaged = squeeze(mean(salt(i,:,stations_mod(istat))));
  mean_chla(i,:) = depth_averaged;
  end
  plot(day,mean_chla,'d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
         
  hold on
  % get the depth averaged model output at our station
  %depth_averaged=squeeze(sum(SALT_mod(:,:,stations_mod(istat)),1))./depth(:,stations_mod(istat))';
  %  p1=  plot(mod_time,depth_averaged,'-d','color','k')%,'markersize')%,mod_station_data(ilay,days,6,istat)^1.2+10);  % make a plot where the size of the dot is relative to the depth
  %          hold on
% get the dates for the real world sampling and average them if there were
% duplicates
            [unique_dates]=unique(stations(istat).variables(myvar_id).date); %ic refernces where the dates came from
            
      cbp_station_date= stations(istat).variables(myvar_id).date ; 
      cbp_station_value= stations(istat).variables(myvar_id).value;      

        
      % now for each unique date, find corresponding values and take the
      % average
      mean_values=[];
            for idate=1:length(unique_dates);                
                mydate=find(cbp_station_date==unique_dates(idate))
                mean_values(idate)= mean(cbp_station_value(mydate));      
            end

     plot(unique_dates,mean_values,'*','color','r','markersize',10);
     
%       end

      %legend([p1,p2],'Modeled','Observed');
      xlabel('Days in 2018','Fontsize',20);ylabel('Salinity (ppt)','Fontsize',20)  
      saveas(gcf,[statstring '/SALT'],'png');
      saveas(gcf,[statstring '/SALT'],'eps');
       saveas(gcf,[statstring '/SALT'],'fig'); 
        end       
       
        %  Statistics 
       % run some stats
          % get the matching model times for statistics;   
            [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
             X_out=mean_values(ib)'; % observational data;
       
             Y_out=[]; % take the average model data over a 3 day window instead of a direct comparison
       
       for i = 1 : length(ia)
           Y_out=[Y_out; mean(depth_averaged(ia(i)-Swindow:ia(i)+Swindow))];
       end
       
       
[MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;

axis([axis_min axis_max axis_min axis_max ]); 

p4=plot(X_out,Y_out,'kd','markersize',12)
xlabel('Observed Salinity (ppt)');ylabel('Modeled Salinity (ppt)');
title([ 'Depth Averaged Salinity ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);


      saveas(gcf,[statstring '/SALT_stats'],'png');
       saveas(gcf,[statstring '/SALT_stats'],'eps');     
       saveas(gcf,[statstring '/SALT_stats'],'fig');        
             
end





%         %algae1 (g C m^-2)
%         alg1_mod(ilay,:,:) = squeeze(hisdata(:,:,find(strcmp(varnames,'ALG1')))).*dzz;
%          %algae2 (g C m^-2)
%         alg2_mod(ilay,:,:)  = squeeze(hisdata(:,:,find(strcmp(varnames,'ALG2')))).*dzz;
%         %Dissolved O2 (g O2 m^-2)
%         DO_mod(ilay,:,:)  = squeeze(hisdata(:,:,find(strcmp(varnames,'DO')))).*dzz;
%         %NH4 (g N m^-2)
%         NH4_mod(ilay,:,:)  = squeeze(hisdata(:,:,find(strcmp(varnames,'NH4')))).*dzz;
%         %NO3 (g N m^-2)
%         NO3_mod(ilay,:,:)  = squeeze(hisdata(:,:,find(strcmp(varnames,'NO3')))).*dzz;
%         %PO4 (g P m^-2)        
%         PO4_mod(ilay,:,:)  = squeeze(hisdata(:,:,find(strcmp(varnames,'PO4')))).*dzz;
%         %DOC (g C m^-2)
%          DOC(ilay,:,:) = squeeze(sum(hisdata(:,:,17:22),3)).*dzz;
%          %POC (g C m^-2)
%          POC(ilay,:,:) = squeeze(sum(hisdata(:,:,15:16),3)).*dzz;
%         %DON (g N m^-2)         
%         DON(ilay,:,:) = squeeze(sum(hisdata(:,:,23:28),3)).*dzz;
%         %DOP (g P m^-2)
%         DOP(ilay,:,:) = squeeze(sum(hisdata(:,:,29:34),3)).*dzz;      
%         
%         SALT_mod(ilay,:,:)  = squeeze(hisdata(:,:,find(strcmp(varnames,'S')))).*dzz; 