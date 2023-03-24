%%%%% ICM Post Processing
%%% L Ziegler, 2021, HPL

% This script reads in the .nc model outputs and extract variables of
% interest and puts them into a daily and monthly timeseries data frame

% idir = dir('/Lacie/Housatonic_BGC/tonic_grid16Sep.mat/','*.nc');
gridfilein = input('What is the full path and file name to the *.mat grid file? --->'); 
load(gridfilein);
numfile = input('What is the number of files you want to be read? --->');


DO_mod    = [];
alg1_mod  = [];
alg2_mod  = [];
NH4_mod   = [];
NO3_mod   = [];
lpoc_mod  = [];
rpoc_mod  = [];
lpon_mod  = [];
rpon_mod  = [];
wcdoc1    = [];
wcdoc2    = [];
wcdoc3    = [];
wcncdoc1  = [];
wcncdoc2  = [];
wcncdoc3  = [];
wc_cdon1  = [];
wc_cdon2  = [];
wc_cdon3  = [];
wc_ncdon1 = [];
wc_ncdon2 = [];
wc_ncdon3 = [];
chla_mod  = [];
PON_mod   = [];
POC_mod   = [];
DON_mod   = [];
DOC_mod   = [];


kd       = [];
mod_time = [];
depth   = [];

for idir = 1:numfile;
    
    numfile_count = [1:numfile];
    
    kstr = num2str(numfile_count(idir),'%04d');
    
    fname = [nameseg, '_', kstr, '.nc'];
    
    nc = netcdf([fname]);
    
    
    DO_mod    = squeeze(nc{'DOXG'}(:,1,:));
    alg1_mod  = squeeze(nc{'B1'}(:,1,:));
    alg2_mod  = squeeze(nc{'B2'}(:,1,:));
    NH4_mod   = squeeze(nc{'NH4'}(:,1,:));
    NO3_mod   = squeeze(nc{'NO3'}(:,1,:));
    lpoc_mod  = squeeze(nc{'LPOC'}(:,1,:));
    rpoc_mod  = squeeze(nc{'RPOC'}(:,1,:));
    lpon_mod  = squeeze(nc{'LPON'}(:,1,:));
    rpon_mod  = squeeze(nc{'RPON'}(:,1,:));
    wcdoc1    = squeeze(nc{'WC_CDOC1'}(:,1,:));
    wcdoc2    = squeeze(nc{'WC_CDOC2'}(:,1,:));
    wcdoc3    = squeeze(nc{'WC_CDOC3'}(:,1,:));
    wcncdoc1  = squeeze(nc{'WC_NCDOC1'}(:,1,:));
    wcncdoc2  = squeeze(nc{'WC_NCDOC2'}(:,1,:));
    wcncdoc3  = squeeze(nc{'WC_NCDOC3'}(:,1,:));
    wc_cdon1  = squeeze(nc{'WC_CDON1'}(:,1,:));
    wc_cdon2  = squeeze(nc{'WC_CDON2'}(:,1,:));
    wc_cdon3  = squeeze(nc{'WC_CDON3'}(:,1,:));
    wc_ncdon1 = squeeze(nc{'WC_NCDON1'}(:,1,:));
    wc_ncdon2 = squeeze(nc{'WC_NCDON2'}(:,1,:));
    wc_ncdon3 = squeeze(nc{'WC_NCDON3'}(:,1,:));
    
    kd        = squeeze(nc{'KD'}(:,1,:));
    mod_time  = nc{'time'}(:)./24/86400;
    depth1    = nc{'depth'}(:);

    
    
    PON_in    = lpon_mod + rpon_mod;
    POC_in    = lpoc_mod + rpoc_mod;

    chla_mod  = (alg1_mod./50+alg2_mod./50).*1000;   % convert mg C/l to ug chl/l based on chl:C ratio from model              
    
    PON_mod   = PON_in + (alg1_mod+alg2_mod)./5.68; % Use Redfield to convert algae C to N
    POC_mod   = POC_in + alg1_mod+alg2_mod;
    
    DOC_mod   = (wcdoc1+wcdoc2+wcdoc3) + (wcncdoc1+wcncdoc2+wcncdoc3);
    DON_mod   = (wc_cdon1 + wc_cdon2 + wc_cdon3) + (wc_ncdon1 + wc_ncdon2 + wc_ncdon3);

    DO_mod    = [DO_mod; DO_mod];

    chla_mod  = [chla_mod; chla_mod];
    NH4_mod   = [NH4_mod; NH4_mod];
    NO3_mod   = [NO3_mod; NO3_mod];

    PON_mod   = [PON_mod; PON_mod];
    POC_mod   = [POC_mod; POC_mod];

    DOC_mod   = [DOC_mod; DOC_mod];
    DON_mod   = [DON_mod; DON_mod];
    
    kd        = [kd; kd];
    mod_time  = [mod_time; mod_time];
    depth    = [depth; depth1];
    
end

% daily timeseries
big_chla = zeros(numfile, length);
big_nh4  = zeros(numfile, 1);
big_no3  = zeros(numfile, 1);
big_pon  = zeros(numfile, 1);
big_poc  = zeros(numfile, 1);
big_doc  = zeros(numfile, 1);
big_don  = zeros(numfile, 1);
%big_depth1 = zeros(numfile, );
%big_time1 = zeros(365,1);
%big_kd1 = zeros(365,7430);


    for ii = 1:length(numfile);
    if ii*24<mod_time
    big_chla1(ii,:) = mean(chla_mod1((ii-1)*24+1:ii*24,:));
    big_nh41(ii,:) = mean(NH4_mod1((ii-1)*24+1:ii*24,:));
    big_no31(ii,:) = mean(NO3_mod1((ii-1)*24+1:ii*24,:));
    big_pon1(ii,:) = mean(PON_mod1((ii-1)*24+1:ii*24,:));
    big_poc1(ii,:) = mean(POC_mod1((ii-1)*24+1:ii*24,:));
    big_doc1(ii,:) = mean(DOC_mod1((ii-1)*24+1:ii*24,:));
    big_don1(ii,:) = mean(DON_mod1((ii-1)*24+1:ii*24,:));
    big_do1(ii,:) = mean(DO_mod1((ii-1)*24+1:ii*24,:));
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
    %big_depth1(ii,:) = mean(depth1_1((ii-1)*24+1:end,:));
    %big_time1(ii,1) = mean(mod_time1((ii-1)*24+1:end,1));
    %big_kd1(ii,:) = mean(kd1((ii-1)*24+1:end,:));
    end
    end



