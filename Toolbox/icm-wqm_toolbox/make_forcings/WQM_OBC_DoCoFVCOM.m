% Script to make the OBC file for the FVCOM-ICM model

% B Clark, UMCES, HPL, March 2018



         % NOW HAVE A STRUCT WITH ALL THE INFORMATION WE NEED
         % need to extract each date and interpolate to the model depth
         %first need to know the model depth
         
% structin=input('What is the input *.mat file with the station structure ---> ');
structin='/Users/bclark/Desktop/Water_Quality_forcing/CBP_mainstem_Dorchester_2016.mat';
% structin='/home/bclark/Desktop/Water_Quality_forcing/CB52_only.mat';
load(structin);
% gridfilein=input('What is the full path and file name to the *.mat grid file? ---> ');
gridfilein='/Users/bclark/Desktop/Water_Quality_forcing/Blackwater_grid_v4_5_11a_fineMESH.mat';

load(gridfilein);

temp_obc_file=importdata('/Users/bclark/Desktop/Water_Quality_forcing/Blackwater_obc.dat')
boundary_ids = temp_obc_file(:,1);

% get the coordinates of the boundary nodes
model_boundary_coords=lld_n(boundary_ids,2:4);

% now find the closest CBP stations to each model node   
station_coords=([stations.lon; stations.lat])';     
obc_stations=dsearchn(station_coords,model_boundary_coords(:,1:2));


% plot the 
%  figure;
%  tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n);
%  hold on
%  plot(station_coords(:,1),station_coords(:,2),'kd','markersize',15);
 
 
% ncfilein=('/Users/bclark/Desktop/test1/Blackwater_0001.nc') ;
% nc=netcdf(ncfilein);

model_depth=model_boundary_coords(:,3);%nc{'h'}(boundary_ids);
model_slevels=[-0, -0.03162277, -0.08944272, -0.1643168, -0.2529822, -0.3535534, ...
    -0.464758, -0.585662, -0.7155418, -0.853815, -1] ; %nc{'s_hybrid'}(:,boundary_ids);

% clear nc

 
write_flag = input('Do you want to write the obc input file? 0 for no 1 for yes ---> ');
plot_flag=input('Do you want to plot contours of the CBP data?  0 for no 1 for yes ---> ');

 
%%


% loop over each node, each variable and each time and interpolate to the
% depth of the model

for inode = 1: length(model_boundary_coords);
    
    mystation=obc_stations(inode);% get the station for this node
    
    nvars=length(stations(mystation).variables);% get the number of variables
    
    % set up an array with each depth layer specified, evenly in the
    % vertical taken from the hydrodynamic model output
    depth_array=model_depth(inode).*model_slevels(2:end).*-1;
    
    for ivar= 1 : nvars;
        
    unique_dates=unique(stations(mystation).variables(ivar).date);% get the dates for this variable
   
    if(length(unique_dates)==15) % get a big array for interpolating to when data is lacking
        big_unique_dates=unique_dates;
    end
    
%             pause
    good_ids=[];
     depth_interped_out=[];
     date_counter=[];
         for jtime=1:length(unique_dates); % now go into each date and interpolate to the model node depth
                         
             good_ids=find(stations(mystation).variables(ivar).date==unique_dates(jtime));
             
             good_dates=stations(mystation).variables(ivar).date(good_ids);
             
              myvalues = stations(mystation).variables(ivar).value(good_ids);
              
              % sometimes the depths are out of order, need to sort
              mydepth = stations(mystation).variables(ivar).depth(good_ids);
              [sorted_depth,sorted_ids]=sort(mydepth);
              sorted_values=myvalues(sorted_ids);
              
              % also need to eliminate duplicates, lets just take the mean
              % value for all duplicates 
              [unique_sorted_depths,unique_sorted_ids]=unique(sorted_depth);     
              
                  good_sorted_values=[];

                for iz=1:length(unique_sorted_ids);
                    
                    good_sorted_ids= find(sorted_depth==sorted_depth(unique_sorted_ids(iz)));
                    good_sorted_values(iz)=nanmean(sorted_values(good_sorted_ids));
                                          
                end
             if(length(unique_sorted_depths)>1);  % interpolate in depth range      
                 
                depth_interped=interp1(unique_sorted_depths,good_sorted_values,depth_array,'linear',median(good_sorted_values));
                
                date_counter=[date_counter jtime];

                depth_interped_out=[depth_interped_out depth_interped'];
             end  
             

             
         end

              time_array_out=1:14:365 ; % biweekly time series

          for iz=1:length(depth_array);
              
              time_interped_out(iz,:,ivar,inode)=interp1(unique_dates(date_counter),depth_interped_out(iz,:),time_array_out,'linear',median(good_sorted_values));
              
          end

         var_names_out(ivar)=stations(mystation).variables(ivar).name;
         

    end
      
end
%%
% plot some time points to double check

% % set up some matrices for building contours near the boundarys;
% date_mat=repmat(unique_dates',length(depth_array),1);
% depth_mat=repmat(depth_array,1,length(unique_dates));
% closest_model_nodes=dsearchn(model_boundary_coords(:,1:2),station_coords);

if(plot_flag);
    mkdir('Contours');
    closest_model_nodes=dsearchn(model_boundary_coords(:,1:2),station_coords);
    
    station_names=stations(:).name;
    
    for j = 1 : length(closest_model_nodes);
        
        station_names=stations(j).name;
        
        
        mynode= closest_model_nodes(j);
        depth_array=model_depth(mynode).*model_slevels(2:end).*-1;
        
        % set up some matrices for building contours near the boundarys;
        date_mat=repmat(time_array_out,length(depth_array),1);
        depth_mat=repmat(depth_array',1,length(time_array_out));
        
        for i = 1 : nvars;
            
            plotting_vars=squeeze(squeeze(time_interped_out(:,:,i,closest_model_nodes(j))));
            figure;
            contourf(date_mat,depth_mat*-1,plotting_vars);
            title([var_names_out(i) ' at station ' station_names])
            colorbar;
            saveas(gcf,['./Contours/' , char(var_names_out(i)) '_node_' num2str(j)] ,'png')
            
        end
        
    end
end

%%
write_flag=1;

if(write_flag)
    
    
    
    Redfield_CtoNmass_ratio = (106*12.011)/(16*14.007); %(mols C*12.011g/molC)/(mols N*14g/mol N) = g C / g N = 5.6786
    Redfield_NtoPmass_ratio = (16*14.007)/(1*30.974); %(mols C*12.011g/molC)/(mols P*30.974 g/mol P)= g C/ g P = 41.1043
    
    carbon2chl1=75; %carbon to chlorophyll a ratio micro g C/ micro chl
    carbon2chl2=60;
    %what fraction of the organics is labile and refractory??
    % sum of L and R must equal 1
    RPOM_frac = 0.5;
    LPOM_frac = 0.5;
    
    % RDOM_FRAC = 0.5;
    % LDOM_FRAC = 0.5;
    DOM1_FRAC = 0.4; % Keller and Hood 2011
    DOM2_FRAC = 0.59;%
    DOM3_FRAC = 0.01;
    Colored_frac = 0.46; %following data collected on NASA GEOCAPE cruise when regressing DOC vs a355
    nonColored_frac = 0.54;%
    
    %fill out another structure for ICM;
    % wqm_vars_openboundary_ICM!=wqm_vars_openboundary;
    
    alg1_frac=0.6; %fraction of algal biomass that is diatoms;
    alg2_frac=1-alg1_frac; %fraction of algal biomass that is dinoflagellates
    
    % loop through both boundaries once the formulas are figured out;
    
    chla_id=find(strcmp(var_names_out,'CHLA'));
    TON_id=find(strcmp(var_names_out,'TON'));
    DON_id=find(strcmp(var_names_out,'DON'));
    DOP_id=find(strcmp(var_names_out,'DOP'));
    DO_id=find(strcmp(var_names_out,'DO'));
    NH4_id=find(strcmp(var_names_out,'NH4'));
    NO3_id=find(strcmp(var_names_out,'NO3'));
    TSS_id=find(strcmp(var_names_out,'TSS'));
    PO4_id=find(strcmp(var_names_out,'PO4'));
    TP_id=find(strcmp(var_names_out,'TP'));
    PP_id=find(strcmp(var_names_out,'PP'));
    TEMP_id=find(strcmp(var_names_out,'Temperature'));
    SALT_id=find(strcmp(var_names_out,'Salinity'));
    
    %algae
    ALG1=squeeze(time_interped_out(:,:,chla_id,:))*alg1_frac.*carbon2chl1/1000; %convert from micrograms to mg C /l
    ALG2=squeeze(time_interped_out(:,:,chla_id,:))*alg2_frac.*carbon2chl2/1000; %convert from micrograms to mg C /l
    
    % ========== Organic Nitrogen ==========
    PON=squeeze(time_interped_out(:,:,TON_id,:))-squeeze(time_interped_out(:,:,DON_id,:)); % mg N l^-1
    DON=squeeze(time_interped_out(:,:,DON_id,:));  % mg N l^-1
    %fractionate the PON and DON
    LPON=PON*LPOM_frac;
    RPON=PON*RPOM_frac;
    CDON1=DON*Colored_frac*DOM1_FRAC;
    CDON2=DON*Colored_frac*DOM2_FRAC;
    CDON3=DON*Colored_frac*DOM3_FRAC;
    NCDON1=DON*nonColored_frac*DOM1_FRAC;
    NCDON2=DON*nonColored_frac*DOM2_FRAC;
    NCDON3=DON*nonColored_frac*DOM3_FRAC;
    
    % ========= organic Carbon ========
    
    %DOC (Redfield)
    % DOC C:N ratio from historical data at CB52
    DOM_CtoN=14; % mg C mg N^-^1 CLIS(14±3) Vlahos & Whitney 2020
    DOC=DON*14;
    
    CDOC1=DOC*Colored_frac*DOM1_FRAC;
    CDOC2=DOC*Colored_frac*DOM2_FRAC;
    CDOC3=DOC*Colored_frac*DOM3_FRAC;
    NCDOC1=DOC*nonColored_frac*DOM1_FRAC;
    NCDOC2=DOC*nonColored_frac*DOM2_FRAC;
    NCDOC3=DOC*nonColored_frac*DOM3_FRAC;
    %POC (Redfield)
    POC=PON.*Redfield_CtoNmass_ratio;
    LPOC=POC*LPOM_frac;
    RPOC=POC*RPOM_frac;
    
    
    % ========== organic Phosphorus =============
    %POP (Redfield from PON)
    POP=PON./Redfield_NtoPmass_ratio;
    LPOP=POP*LPOM_frac;
    RPOP=POP*RPOM_frac;
    % DOP
    DOP=squeeze(time_interped_out(:,:,DOP_id,:));
    CDOP1=DOP*Colored_frac*DOM1_FRAC;
    CDOP2=DOP*Colored_frac*DOM2_FRAC;
    CDOP3=DOP*Colored_frac*DOM3_FRAC;
    NCDOP1=DOP*nonColored_frac*DOM1_FRAC;
    NCDOP2=DOP*nonColored_frac*DOM2_FRAC;
    NCDOP3=DOP*nonColored_frac*DOM3_FRAC;
    
    
    % inorganic phosphorous
    
    %particulate inorganic Phosphorus
    PIP = time_interped_out(:,:,PP_id,:); % PP from bay program is particualte Phosphate
    %dissolved PO4
    PO4=time_interped_out(:,:,PO4_id,:);
    %total PO4
    totalPO4=PIP+PO4;
    
    
    %inorganic Nitrogen
    
    NH4=time_interped_out(:,:,NH4_id,:);
    NO3=time_interped_out(:,:,NO3_id,:);
    
    %make sure units are correct for ICM; all bay program data comes in mg/l
    %accept chla
    
    % temperatuer and salinity, not used but mine as well add it in there
    Temp=time_interped_out(:,:,TEMP_id,:);
    Salt = time_interped_out(:,:,SALT_id,:);
    TP_id;
    % tss
    
    TSS=time_interped_out(:,:,TSS_id,:);
    
    % oxygen
    
    DO = time_interped_out(:,:,DO_id,:);
    
end
    %%
    
    [nlayers,ntimes,nvars,nnodes]=size(time_interped_out);
    
    nvars=44; % we use 44 variables for the ICM model
    
    
    %now there are all of the variables, print out in the format the WQM
    %accepts ;
    %
    % %the list of variables to call into the ICM printing loop
    % ICM_varlist=char('WTEMP','SAL','TSS',...
    %                 'ALG1','ALG2','ALG3','ZOO1','ZOO2',...  %plankton
    %                 'LDOC','RDOC','LPOC','RPOC',...                %carbon
    %                 'NH4F','NO23F','UREA','LDON','RDON','LPON','RPON',... %nitrogen
    %                 'totalPO4','LDOP','RDOP','LPOP','RPOP','PIP',... %phosphorus
    %                 'COD','DO',...                                  %oxygen
    %                 'PSI','SIF',...                                         %silica
    %                 'tz1','tz2','tz3');             % round out the list with matrices with zeros;
    %
    ICM_varlist=char('WTEMP','SAL','TSS',...
        'ALG1','ALG2','ALG3','ZOO1','ZOO2',...  %plankton
        'CDOC1','NCDOC1','LPOC','RPOC',...                %carbon
        'NH4F','NO23F','UREA','CDON1','NCDON1','LPON','RPON',... %nitrogen
        'totalPO4','CDOP1','NCDOP1','LPOP','RPOP','PIP',... %phosphorus
        'COD','DO',...                                  %oxygen
        'PSI','SIF',...                                         %silica
        'tz1','tz2','tz3',...
        'CDOC2','NCDOC2','CDOC3','NCDOC3','CDON2','NCDON2','CDON3','NCDON3',...
        'CDOP2','NCDOP2','CDOP3','NCDOP3' );             % round out the list with matrices with zeros;
    %
    % now fill the array
    big_output_array=zeros(nlayers,ntimes,nvars,nnodes);
    big_output_array(:,:,1,:)=Temp;
    big_output_array(:,:,2,:)=Salt;
    big_output_array(:,:,3,:)=TSS;
    big_output_array(:,:,4,:)=ALG1;
    big_output_array(:,:,5,:)=ALG2;
    6
    7
    8
    big_output_array(:,:,9,:)=CDOC1;
    big_output_array(:,:,10,:)=NCDOC1;
    big_output_array(:,:,11,:)=LPOC;
    big_output_array(:,:,12,:)=RPOC;
    big_output_array(:,:,13,:)=NH4;
    big_output_array(:,:,14,:)=NO3;
    15
    big_output_array(:,:,16,:)=CDON1;
    big_output_array(:,:,17,:)=NCDON1;
    big_output_array(:,:,18,:)=LPON;
    big_output_array(:,:,19,:)=RPON;
    
    big_output_array(:,:,20,:)=totalPO4;
    big_output_array(:,:,21,:)=CDOP1;
    big_output_array(:,:,22,:)=NCDOP1;
    big_output_array(:,:,23,:)=LPOP;
    big_output_array(:,:,24,:)=RPOP;
    big_output_array(:,:,25,:)=PIP;
    26
    big_output_array(:,:,27,:)=DO;
    28
    29
    30
    31
    32
    big_output_array(:,:,33,:)=CDOC2;
    big_output_array(:,:,34,:)=NCDOC2;
    big_output_array(:,:,35,:)=CDOC3;
    big_output_array(:,:,36,:)=NCDOC3;
    
    big_output_array(:,:,37,:)=CDON2;
    big_output_array(:,:,38,:)=NCDON2;
    big_output_array(:,:,39,:)=CDON3;
    big_output_array(:,:,40,:)=NCDON3;
    
    big_output_array(:,:,41,:)=CDOP2;
    big_output_array(:,:,42,:)=NCDOP2;
    big_output_array(:,:,43,:)=CDOP3;
    big_output_array(:,:,44,:)=NCDOP3;
    
    
    
    %%
    obc_file_out='tonic_obc_wq_CtoNDOCinc.dat'; %printing file
    fid = fopen(obc_file_out,'w');
    %
    fprintf(fid,'%d\n',length(boundary_ids));
    
    for ifp=1:length(boundary_ids); % print the layer information at the top
        
        fprintf(fid,'%16.4f %16.4f\n',temp_obc_file(ifp,1:2));
        
    end
    %
    for tt=1 : length(time_array_out)
        
        t_hours=time_array_out(tt).*24-24; % time in hours from the first day
        
        fprintf(fid,'%10.4f\n',t_hours); %output time
        
        for ivar = 1 : nvars;
            
            for inode= 1 : nnodes;
                
                fprintf(fid,'%d' ,boundary_ids(inode));
                fprintf(fid,'%16.8f',big_output_array(:,tt,ivar,inode));%,varname_print); %print the values at each node
                fprintf(fid,'\n');
                
            end
        end
        
    end
    fclose(fid);
    
    











