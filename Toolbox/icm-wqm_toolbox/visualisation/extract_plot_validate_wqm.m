% script to plot the WQM output at CB stations to get a comparison and do
% some basic statistics

% B Clark sep 2015
% UMCES, HPL
% ========================================================================
% Updated by L Ziegler 2020 for Long Island Sound

% Load Grid file
matfilein  = input('What is the full path and file name to the *.mat grid file? ---> ','s');
% '/Volumes/LaCie/Toolkit/Mirror/Housatonic_BGC/Housatonicwqm_grid.mat';
load(matfilein);
% Load ICM nc files
wqm_hisdir = input('What is full path to the WQM history directory? --->  ','s');
% /Volumes/LaCie/Toolkit/Mirror/Housatonic_BGC/run1_Jun2020
nlayers    = input('How many layers are in the model? --> ');
outdir     = 'nutrient_plots';
ICM_int    = input('What is the interval that the data was output at, in days?  ---> ');

% Load observational data
load('/data/users/bclark/CBAY_DATA/Sediment_stations.mat');
CBP_struct = input('What is the full path and file name to the .mat file with the CBP data structure ---> ','s');
load(CBP_struct);
first_day = input('What is the first day of the model? --->  ');

%domflag = input('Is this the WQM with the new DOM or old DOM formulas? enter 0 for no and 1 for yes ---->');

numfile=input('How many files do you want to compare? ---> '); 
nameseg=input('What is the prefix for the ICM files? ---> ','s');

kb=11;
for k = 1 : kb;
    s_hybrid(k,1:length(xyd_n))=-((k-1)/(kb-1))^1.5;
end

for i = 2:11;
   layer_thickness(i-1,:)=mean(s_hybrid(i-1:i,:),1)*-1;
end

clear nc
%
if(~exist(outdir));
    mkdir(outdir);
end

% get the station indices by matching the CBP data coordinates
% with the FVCOM grid coordinates
stations_mod = dsearchn([lld_n(:,2) lld_n(:,3)],([stations.lon ; stations.lat])'); %find the model nodes that are close to the wqm station points

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n);
hold on
hold on;
plot(lon_big,lat_big,'-k','linewidth',1)
xlim([-73.18 -72.9])
ylim([41.06 41.4])
hold on ;

for i = 1 : length(stations_mod);
   text(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),stations(i).name,'fontsize',24)
    
end

saveas(gcf,'AllStations_map','png');

%% Extract data from model output and save as a variable
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
for idir=1 : numfile;   
    numfile_count = [1:numfile];   
    kstr = num2str(numfile_count(idir),'%04d');    
    fname = [wqm_hisdir '/' nameseg,'_',kstr,'.nc'];
    nc = netcdf([fname]);
    idir
    
    alg1_in= nc{'B1'}(:);
    alg2_in  = nc{'B2'}(:);
    DO_in  = nc{'DOXG'}(:);
    NH4_in  = nc{'NH4'}(:);
    NO3_in = nc{'NO3'}(:);
    depth_in= nc{'depth'}(:);
    DON_in = nc{'WC_CDON1'}(:) + ...
                                nc{'WC_NCDON1'}(:) + ...
                                nc{'WC_CDON2'}(:) + ...
                                nc{'WC_NCDON2'}(:) + ...
                                nc{'WC_CDON3'}(:) + ...
                                nc{'WC_NCDON3'}(:);
    DOC_in = nc{'WC_CDOC1'}(:) + ...
                                nc{'WC_NCDOC1'}(:) + ...
                                nc{'WC_CDOC2'}(:) + ...
                                nc{'WC_NCDOC2'}(:) + ...
                                nc{'WC_CDOC3'}(:) + ...
                                nc{'WC_NCDOC3'}(:);
                            
    POC_in= nc{'LPOC'}(:) + nc{'RPOC'}(:);
    PON_in= nc{'LPON'}(:) + nc{'RPON'}(:);
    
    SALT_in = nc{'salinity'}(:);
    TEMP_in  = nc{'temp'}(:);
    KD_in = nc{'KD'}(:);
    
    
    
    mod_time=[mod_time;nc{'time'}(:)./86400];
    
    alg1_mod=[alg1_mod ; alg1_in];
    alg2_mod=[alg2_mod ; alg2_in];
    DO_mod=[DO_mod ; DO_in];
    NH4_mod=[NH4_mod ; NH4_in];
    NO3_mod=[NO3_mod ; NO3_in];
    
     depth1=[depth1 ; depth_in];
    DON_mod=[DON_mod ; DON_in];
    DOC_mod=[DOC_mod ; DOC_in];
    
    PON_mod=[PON_mod ; PON_in]+(alg1_mod+alg2_mod)./5.68;
    POC_mod=[POC_mod ; POC_in]+alg1_mod+alg2_mod;
    
    SALT_mod=[SALT_mod ; SALT_in];
    TEMP_mod=[TEMP_mod ; TEMP_in];
    
    KD_mod=[KD_mod ; KD_in];
    

%     alg1_mod=[alg1_mod ; alg1_in(:,:,stations_mod)];
%     alg2_mod=[alg2_mod ; alg2_in(:,:,stations_mod)];
%     DO_mod=[DO_mod ; DO_in(:,:,stations_mod)];
%     NH4_mod=[NH4_mod ; NH4_in(:,:,stations_mod)];
%     NO3_mod=[NO3_mod ; NO3_in(:,:,stations_mod)];
%     
%      depth1=[depth1 ; depth_in(:,stations_mod)];
%     DON_mod=[DON_mod ; DON_in(:,:,stations_mod)];
%     DOC_mod=[DOC_mod ; DOC_in(:,:,stations_mod)];
%     
%     PON_mod=[PON_mod ; PON_in(:,:,stations_mod)]+(alg1_mod+alg2_mod)./5.68;
%     POC_mod=[POC_mod ; POC_in(:,:,stations_mod)]+alg1_mod+alg2_mod;
%     
%     SALT_mod=[SALT_mod ; SALT_in(:,:,stations_mod)];
%     TEMP_mod=[TEMP_mod ; TEMP_in(:,:,stations_mod)];
%     
%     KD_mod=[KD_mod ; KD_in(:,:,stations_mod)];
    
    clear nc
    
end


starting_date=first_day; 
data_int=ICM_int; %1/24 -- hourly output

% convert mg C/l to ug chl/l based on chl:C ratio from model
chla_mod = (alg1_mod./50+alg2_mod./50).*1000;

mod_time=mod_time+starting_date;

[ntimes,nlayers,nstations]=size(big_chla);

for iz=1:nlayers;
    
    depth(:,iz,:)=depth1.*repmat(layer_thickness(iz,:,1));

end

%%
rounded_mod_times=round(mod_time);
%%

for istat=1:length(stations_mod);
    
    statstring = ['station' num2str(istat)];
    mkdir(statstring)
    
    figure
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),-h_n);
    hold on
    plot(lld_n(stations_mod(istat),2),lld_n(stations_mod(istat),3),'*','color','k','markersize',24);
    
    saveas(gcf,[statstring '/station_map'],'png');
    
    %     for ivar=1: length(stations(istat).variables)
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'CHLA'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
   
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,chla_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    xlabel('Days in 2005');ylabel('Depth (m)')
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        for iday=1:length(ib)
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,stations_mod(2))',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,stations_mod(2));
        modelValues=big_chla(ia(iday),model_depthIds,stations_mod(2));
        end  
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=mean(myUniqueValues);
        Y_out_mean=mean(modelValues);
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    h.Label.String='chla (ug chla l^-^1)'   ;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    title([ 'Chl a ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    subplot(2,2,3)
    
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    
    xlabel('Observed Chl a (ug l^-^1)');ylabel('Modeled Chl a (ug l^-^1)')  ;
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(chla_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('Chl a (ug l^-^1)');%ylabel('Probability')
    saveas(gcf,[statstring '/chla_depth'],'png');
    saveas(gcf,[statstring '/chla_depth'],'eps');
    saveas(gcf,[statstring '/chla_depth'],'fig');
    
    
    
    figure;     % depth averaged
    subplot(2,2,[1,2]);
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    
    plot(mod_time,mean(chla_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged Chl a ' ,'  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Chl a (ug l^-^1)');
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed Chl a (ug l^-^1)');ylabel('Modeled Chl a (ug l^-^1)')  ;
    
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(chla_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('Chl a (ug l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/chla_depth_stats'],'png');
    saveas(gcf,[statstring '/chla_depth_stats'],'eps');
    saveas(gcf,[statstring '/chla_depth_stats'],'fig');
    %
    %========================= Dissolved Oxygen ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DO'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,DO_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=DO_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='DO (mg O_2 l^-^1)'       ;
    
  
    title([ 'DO ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DO (mg O_2 l^-^1)');ylabel('Modeled DO (mg O_2 l^-^1)')  ;
    
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(DO_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    ylabel('DO (mg O_2 l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/DO_depth'],'png');
    saveas(gcf,[statstring '/DO_depth'],'eps');
    saveas(gcf,[statstring '/DO_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(DO_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged DO ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Observed DO (mg O_2 l^-^1)');
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DO (mg O_2 l^-^1)');ylabel('Modeled DO (mg O_2 l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(DO_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('DO (mg O_2 l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/DO_depth_stats'],'png');
    saveas(gcf,[statstring '/DO_depth_stats'],'eps');
    saveas(gcf,[statstring '/DO_depth_stats'],'fig');
    
    %
    %========================= Nitrate ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'NO23F'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2])
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,NO3_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=NO3_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='NO_3^- (mg N l^-^1)'       ;
    
    title([ 'NO_3^- ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('NO_3^- (mg N l^-^1)')
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NO_3^- (mg N l^-^1)');ylabel('Modeled NO_3^- (mg N l^-^1)');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(NO3_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('NO_3(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/NO3_depth'],'png');
    saveas(gcf,[statstring '/NO3_depth'],'eps');
    saveas(gcf,[statstring '/NO3_depth'],'fig');
    
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(NO3_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged NO_3^- ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('NO_3^- (mg N l^-^1)')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NO_3^- (mg N l^-^1)');ylabel('Modeled NO_3^- (mg N l^-^1)');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(NO3_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('NO_3(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/NO3_depth_stats'],'png');
    saveas(gcf,[statstring '/NO3_depth_stats'],'eps');
    saveas(gcf,[statstring '/NO3_depth_stats'],'fig');
    %
    %========================= Ammonium ============================
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'NH4F'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,NH4_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=NH4_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='NH_4^+ (mg N l^-^1)'       ;
    
  
    
    title([ 'NH_4^+ ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NH_4^+ (mg N l^-^1)');ylabel('Modeled NH_4^+ (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(NH4_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('NH_4(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/NH4_depth'],'png');
    saveas(gcf,[statstring '/NH4_depth'],'eps');
    saveas(gcf,[statstring '/NH4_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(NH4_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged NH_4^+' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('NH_4^+ (mg N l^-^1)')  ;
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NH_4^+ (mg N l^-^1)');ylabel('Modeled NH_4^+ (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(NH4_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('NH_4(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/NH4_depth_stats'],'png');
    saveas(gcf,[statstring '/NH4_depth_stats'],'eps');
    saveas(gcf,[statstring '/NH4_depth_stats'],'fig');
    
    
    %
    %  %========================= PHOSPHATE ============================
    %
    %        myvar_id=9; % change this ID based on where the variable is located within the CBP data structure
    %                 % generated by read_cbp.m
    %
    %     % pull out the data for this station
    %         [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    %       cbp_station_date= stations(istat).variables(myvar_id).date ;
    %       cbp_station_value= stations(istat).variables(myvar_id).value;
    %       cbp_station_depths=stations(istat).variables(myvar_id).depth;
    %     % now for each unique date, find corresponding values and take the
    %
    %       % get the matching model times for statistics and plotting;
    %        [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    %
    %
    %            % now plot date by date
    %         figure;
    %
    %        X_out=[];
    %        Y_out=[];
    %        X_out_mean=[];
    %        Y_out_mean=[];
    %        unique_dates_out=[];
    %        subplot(2,1,1);
    %
    %           contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,PO4_mod(:,:,istat),'LineStyle','none')
    %           hold on
    %           [cmin,cmax]=caxis;
    %
    %             for iday=1:length(ib)
    %
    %                 % get the good day from the data back out so we can get
    %                 % each value and depth
    %                 myIds=find(cbp_station_date==unique_dates(ib(iday)));
    %
    %                 myDates=cbp_station_date(myIds);
    %                 myDepths=cbp_station_depths(myIds);
    %                 myvalues=cbp_station_value(myIds);
    %
    %                 % get rid of duplicates;
    %                 [myUniqueDepths,Zids]=unique(myDepths);
    %                 myUniqueValues=myvalues(Zids)
    %                 myUniqueDates=unique(myDates);
    %
    %                 % now match up the model
    %                 model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
    %                 modelDepths=depth(ia(iday),model_depthIds,istat);
    %                 modelValues=PO4_mod(ia(iday),model_depthIds,istat);
    %
    %                         for iplot=1:length(modelDepths);
    %                                 scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
    %                         end
    %                 X_out=[X_out;myUniqueValues];
    %                 Y_out=[Y_out;modelValues'];
    %
    %                 % take the depth averaged quantity
    %                 X_out_mean=[X_out_mean;mean(myUniqueValues)];
    %                 Y_out_mean=[Y_out_mean;mean(modelValues)];
    %
    %                 unique_dates_out=[unique_dates_out;myUniqueDates];
    %
    %             end
    %          h= colorbar;
    %         caxis([cmin cmax])
    %         h.Label.String='PO_4^3^- (mg P l^-^1)'       ;
    %
    %  subplot(2,1,2);
    %      [MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    %      title([ 'PO_4^3^+ ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    % %          legend([p1(1),p2,3],'Modeled','Modeled','Observed');
    %         xlabel('Days in 2005');ylabel('PO_4^3^+(mg P l^-^1)')
    %
    %         plot(X_out,Y_out,'kd','markersize',8);
    %
    %        axis([axis_min axis_max axis_min axis_max]);
    %        title([ 'PO_4^3^+' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    %                xlabel('Observed PO_4^3^+ (mg P l^-^1)');ylabel('Modeled PO_4^3^+ (mg P l^-^1)')
    %         saveas(gcf,[statstring '/PO4_depth'],'png');
    %         saveas(gcf,[statstring '/PO4_depth'],'eps');
    %         saveas(gcf,[statstring '/PO4_depth'],'fig');
    %   figure;     % depth averaged
    %
    %   [MEF,r,WMS,RMSE,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    %      subplot(2,1,1);
    %      plot(mod_time,mean(PO4_mod(:,:,istat),2),'k-o');
    %      hold on
    %      plot(unique_dates_out,X_out_mean,'r-*');
    %      title([ 'Depth Averaged PO_4^3^+ ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    %         xlabel('Days in 2005');ylabel('PO_4^3^+(mg P l^-^1)');
    %      legend('Modeled','Observed')
    %      subplot(2,1,2);
    %
    %            plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    %            axis([axis_min axis_max axis_min axis_max]);
    %       xlabel('Observed PO_4^3^+ (mg P l^-^1)');ylabel('Modeled PO_4^3^+ (mg P l^-^1)')
    %         saveas(gcf,[statstring '/PO4_depth_stats'],'png');
    %         saveas(gcf,[statstring '/PO4_depth_stats'],'eps');
    %         saveas(gcf,[statstring '/PO4_depth_stats'],'fig');
    %
    %========================= DON ============================
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DON'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,DON_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=DON_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='DON (mg N l^-^1)'       ;
    
    title([ 'DON ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')
    
    subplot(2,2,3);
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);;
    %     title([ 'DON ' , 'MEF =' num2str(MEF) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed DON (mg N l^-^1)');ylabel('Modeled DON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(DON_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('DON (mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/DON_depth'],'png');
    saveas(gcf,[statstring '/DON_depth'],'eps');
    saveas(gcf,[statstring '/DON_depth'],'fig');
    
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(DON_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged DON ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('DON(mg N l^-^1)')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DON (mg N l^-^1)');ylabel('Modeled DON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(DON_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('DON (mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/DON_depth_stats'],'png');
    saveas(gcf,[statstring '/DON_depth_stats'],'eps');
    saveas(gcf,[statstring '/DON_depth_stats'],'fig');
    %
    %========================= DOC ============================
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DOC'))
            
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.
            break
        end
    end
    CtoN=7.2
    % pull out the data for this station
    
    if(isempty(stations(istat).variables(myvar_id).value));  % pull out the DON data and use a C to N ratio if there is no data
        myvar_id=[];
        ivar=0;
        while isempty(myvar_id)
            ivar=ivar+1;
            if(strcmp(stations(istat).variables(ivar).name,'DON'))
                myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
                % generated by read_cbp.m
            end
        end
        [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
        cbp_station_date= stations(istat).variables(myvar_id).date ;
        cbp_station_value= stations(istat).variables(myvar_id).value*CtoN;
        cbp_station_depths=stations(istat).variables(myvar_id).depth;
        % now for each unique date, find corresponding values and take the
    else
        
        [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
        cbp_station_date= stations(istat).variables(myvar_id).date ;
        cbp_station_value= stations(istat).variables(myvar_id).value;
        cbp_station_depths=stations(istat).variables(myvar_id).depth;
        % now for each unique date, find corresponding values and take the
    end
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,DOC_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=DOC_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='DOC (mg C l^-^1)'       ;
   
    title([ 'DOC ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    %     title([ 'DOC ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed DOC (mg C l^-^1)');ylabel('Modeled DOC (mg C l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(DOC_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('DOC (mg C l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/DOC_depth'],'png');
    saveas(gcf,[statstring '/DOC_depth'],'eps');
    saveas(gcf,[statstring '/DOC_depth'],'fig');
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(DOC_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged DOC ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('DOC(mg C l^-^1)')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DOC (mg C l^-^1)');ylabel('Modeled DOC (mg C l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(DOC_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('DOC (mg C l^-^1)');%ylabel('Probability')
    
    saveas(gcf,[statstring '/DOC_depth_stats'],'png');
    saveas(gcf,[statstring '/DOC_depth_stats'],'eps');
    saveas(gcf,[statstring '/DOC_depth_stats'],'fig');
    %
    %========================= Salinity ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'SALINITY'))
            
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value ;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,SALT_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);;
        modelValues=SALT_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='Salinity'       ;
  
    title([ 'Salinity' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    %     title([ 'Salinity ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed Salinity');ylabel('Modeled Salinity ');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(SALT_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('Salinity');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Salt_depth'],'png');
    saveas(gcf,[statstring '/Salt_depth'],'eps');
    saveas(gcf,[statstring '/Salt_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(SALT_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged Salinity ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Salinity ')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed Salinity');ylabel('Modeled Salinity')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(SALT_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('Salinity');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Salt_depth_stats'],'png');
    saveas(gcf,[statstring '/Salt_depth_stats'],'eps');
    saveas(gcf,[statstring '/Salt_depth_stats'],'fig');
    %========================= Temperature ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'WTEMP'))
            
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value ;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,TEMP_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);;
        modelValues=TEMP_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='Temperature'       ;

    title([ 'Temperature' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    %     title([ 'Salinity ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed Temperature');ylabel('Modeled Temperature');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(TEMP_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('Temperature');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Temp_depth'],'png');
    saveas(gcf,[statstring '/Temp_depth'],'eps');
    saveas(gcf,[statstring '/Temp_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(TEMP_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged Temperature ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Temperature')
    
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed Temperature');ylabel('Modeled Temperature')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(TEMP_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('Temperature');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Temp_depth_stats'],'png');
    saveas(gcf,[statstring '/Temp_depth_stats'],'eps');
    saveas(gcf,[statstring '/Temp_depth_stats'],'fig');
    %
    %========================= Kd ============================
    %
    
%     myvar_id=[];
%     ivar=0;
%     while isempty(myvar_id)
%         ivar=ivar+1;
%         if(strcmp(stations(istat).variables(ivar).name,'KD'))
%             
%             myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
%             % generated by read_cbp.m
%         end
%     end
%     
%     if(~isempty(stations(istat).variables(myvar_id).value));  % pull out the KD data and use a C to N ratio if there is no data
%         
%         % pull out the data for this station
%         [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
%         cbp_station_date= stations(istat).variables(myvar_id).date ;
%         cbp_station_value= stations(istat).variables(myvar_id).value ;
%         cbp_station_depths=stations(istat).variables(myvar_id).depth;
%         % now for each unique date, find corresponding values and take the
%         
%         % get the matching model times for statistics and plotting;
%         [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
%         
%         
%         % now plot date by date
%         figure;
%         
%         X_out=[];
%         Y_out=[];
%         X_out_mean=[];
%         Y_out_mean=[];
%         unique_dates_out=[];
%         subplot(2,1,1);
%         
%         subplot(2,2,[1,2])
%         contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,KD_mod(:,:,istat),'LineStyle','none')
%         hold on
%         [cmin,cmax]=caxis;
%         
%         for iday=1:length(ib)
%             
%             % get the good day from the data back out so we can get
%             % each value and depth
%             myIds=find(cbp_station_date==unique_dates(ib(iday)));
%             
%             myDates=cbp_station_date(myIds);
%             myDepths=cbp_station_depths(myIds);
%             myvalues=cbp_station_value(myIds);
%             
%             % get rid of duplicates;
%             [myUniqueDepths,Zids]=unique(myDepths);
%             myUniqueValues=myvalues(Zids)
%             myUniqueDates=unique(myDates);
%             
%             % now match up the model
%             model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
%             modelDepths=depth(ia(iday),model_depthIds,istat);;
%             modelValues=KD_mod(ia(iday),model_depthIds,istat);
%             
%             for iplot=1:length(modelDepths);
%                 scatter(myUniqueDates,myUniqueDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
%             end
%             X_out=[X_out;myUniqueValues];
%             Y_out=[Y_out;modelValues'];
%             
%             % take the depth averaged quantity
%             X_out_mean=[X_out_mean;mean(myUniqueValues)];
%             Y_out_mean=[Y_out_mean;mean(modelValues)];
%             
%             unique_dates_out=[unique_dates_out;myUniqueDates];
%             
%         end
%         h= colorbar;
%         caxis([cmin cmax])
%         h.Label.String='Kd (m^-^1)'       ;
%         
%         
%         [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
%         title([ 'KD' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
%         %          legend([p1(1),p2,3],'Modeled','Modeled','Observed');
%         xlabel('Days in 2005');ylabel('kd (m^-^1)')  ;
%         
%         subplot(2,2,3)
%         plot(X_out,Y_out,'kd','markersize',8);
%         
%         axis([axis_min axis_max axis_min axis_max]);
%         title([ 'KD ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
%         xlabel('Observed KD (m^-^1)');ylabel('Modeled KD (m^-^1)');
%         
%         
%         
%         subplot(2,2,4);
%         %     histogram(Y_out,10,'normalization','probability');
%         %     hold on
%         %     histogram(X_out,10,'normalization','probability');
%         g1={};g2={};g3={};
%         VAR_cat=reshape(KD_mod(:,:,istat),1,nlayers*ntimes);
%         g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
%         g2(1:length(Y_out))={'Modeled'};
%         g3(1:length(VAR_cat))={'Modeled-ALL'};
%         grouping=[g1 g2 g3];
%         boxplot([X_out ; Y_out; VAR_cat'],grouping);
%         ylabel('Kd (m^-^1)');%ylabel('Probability')
%         
%         saveas(gcf,[statstring '/KD_depth'],'png');
%         saveas(gcf,[statstring '/KD_depth'],'eps');
%         saveas(gcf,[statstring '/KD_depth'],'fig');
%         
%         
%     end% if loop to make sure there are values for kd, some stations dont have it
    
    %========================= POC ============================
    myvar_id=[];
    myvar_id2=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DON'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    ivar=0;
    while isempty(myvar_id2)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'TON'))
            myvar_id2=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    [unique_dates2,unique_ids2]=unique(stations(istat).variables(myvar_id2).date);
    
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_date2= stations(istat).variables(myvar_id2).date ;
    
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_value2= stations(istat).variables(myvar_id2).value;
    
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    cbp_station_depths2=stations(istat).variables(myvar_id2).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    if(length(unique_dates)~=length(unique_dates2));
        disp('Different days for PON and DON, exiting');
        return
    end
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,POC_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        myIds2=find(cbp_station_date2==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        
        myvalues=cbp_station_value(myIds);
        myvalues2=cbp_station_value2(myIds2);
        
        
        POC_in=(myvalues2-myvalues)*5.67; % TON-PON multiplied by redfield (g C g N^-1)
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=POC_in(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=POC_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='POC(mg C l^-^1)'       ;
    
    title([ 'POC ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')
    
    subplot(2,2,3);
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);;
    %     title([ 'DON ' , 'MEF =' num2str(MEF) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed POC (mg C l^-^1)');ylabel('Modeled POC (mg c l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(POC_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/POC_depth'],'png');
    saveas(gcf,[statstring '/POC_depth'],'eps');
    saveas(gcf,[statstring '/POC_depth'],'fig');
    
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(POC_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged POC ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('POC (mg C l^-^1)')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed POC (mg C l^-^1)');ylabel('Modeled POC (mg C l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(POC_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %      set(gca,'FontWeight','bold','LineWidth',2);
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/POC_depth_stats'],'png');
    saveas(gcf,[statstring '/POC_depth_stats'],'eps');
    saveas(gcf,[statstring '/POC_depth_stats'],'fig');
    
    close all
    
    
    
       %========================= PON ============================
    myvar_id=[];
    myvar_id2=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DON'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    ivar=0;
    while isempty(myvar_id2)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'TON'))
            myvar_id2=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    [unique_dates2,unique_ids2]=unique(stations(istat).variables(myvar_id2).date);
    
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_date2= stations(istat).variables(myvar_id2).date ;
    
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_value2= stations(istat).variables(myvar_id2).value;
    
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    cbp_station_depths2=stations(istat).variables(myvar_id2).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    if(length(unique_dates)~=length(unique_dates2));
        disp('Different days for PON and DON, exiting');
        return
    end
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,PON_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        myIds2=find(cbp_station_date2==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        
        myvalues=cbp_station_value(myIds);
        myvalues2=cbp_station_value2(myIds2);
        
        
        PON_in=(myvalues2-myvalues); % TON-PON 
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=PON_in(Zids)
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=PON_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='PON(mg C l^-^1)'       ;
    
    title([ 'PON ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')
    
    subplot(2,2,3);
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);;
    %     title([ 'DON ' , 'MEF =' num2str(MEF) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed PON (mg N l^-^1)');ylabel('Modeled PON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(PON_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/PON_depth'],'png');
    saveas(gcf,[statstring '/PON_depth'],'eps');
    saveas(gcf,[statstring '/PON_depth'],'fig');
    
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(PON_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged PON ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('PON (mg N l^-^1)')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed PON (mg N l^-^1)');ylabel('Modeled PON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(PON_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %      set(gca,'FontWeight','bold','LineWidth',2);
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/PON_depth_stats'],'png');
    saveas(gcf,[statstring '/PON_depth_stats'],'eps');
    saveas(gcf,[statstring '/PON_depth_stats'],'fig');
    
    close all
    
    
    
end
