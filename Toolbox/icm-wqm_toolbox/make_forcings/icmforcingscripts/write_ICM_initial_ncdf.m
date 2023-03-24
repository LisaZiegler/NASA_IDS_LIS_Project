% SCript to write the initial condition file for
% Restart of ICM
%using netcdf outputs
%B Clark Nov 2018

wqm_hisdir = input('What is full path to the WQM history directory? --->  ','s');
nlayers = input('How many layers are in the model? --> ');
% dom_switch = input('Is this the old DOM formula or new DOM formula?  enter 0 for old and 1 for new ---> ')
matfilein = input('What is the full path and file name to the *.mat grid file? ---> ','s');
load(matfilein);


% if dom_switch ==0
%     nvars = 32;
% else
     nvars = 44;
% end

myday=input('What day do you want to write the ICM input for? -->  ');
first_day=input('What is the first day of the model? --->   ');
fname= [wqm_hisdir '/tonic_0001.nc']
%fname= [wqm_hisdir '/rhd_0001.nc']
myinfo=ncinfo(fname);
ngrids=myinfo.Dimensions(2).Length;
nlayers=myinfo.Dimensions(4).Length;
ntimes=myinfo.Dimensions(end).Length;
nc=netcdf(fname);
  
time_in=nc{'time'}(:);
time_int=(time_in(2)-time_in(1))/86400;
mytime=round((myday-first_day)/time_int);

Big_array_out=zeros(nlayers,ngrids,nvars);
  
        Big_array_out(:,:,1)=squeeze(nc{'temp'}(mytime,:,:));%temperature
        Big_array_out(:,:,2)=squeeze(nc{'salinity'}(mytime,:,:));%Salinity
        Big_array_out(:,:,3)=10;%SSI
        Big_array_out(:,:,4)=0.2;%squeeze(nc{'B1'}(mytime,:,:)); %ALG1
        Big_array_out(:,:,5)=0.2;%squeeze(nc{'B2'}(mytime,:,:)); %ALG2
        Big_array_out(:,:,6)=0.0 ;%AL3
        Big_array_out(:,:,7)=0.0 ; % SZ
        Big_array_out(:,:,8)=0.0;  % LZ
        Big_array_out(:,:,9)=0.1663; %squeeze(nc{'WC_CDOC1'}(mytime,:,:)); %CDOC1
        Big_array_out(:,:,10)=0.1952; %squeeze(nc{'WC_NCDOC1'}(mytime,:,:)); %NCDOC1

        Big_array_out(:,:,11)= 0.1140;%squeeze(nc{'LPOC'}(mytime,:,:)); %LPOC
        Big_array_out(:,:,12)= 0.1140;%squeeze(nc{'RPOC'}(mytime,:,:)); %RPOC
        Big_array_out(:,:,13)= 0.0102;%squeeze(nc{'NH4'}(mytime,:,:)); %NH4
        Big_array_out(:,:,14)= 0.02;%squeeze(nc{'NO3'}(mytime,:,:)); %N03
        Big_array_out(:,:,15)=0.0; %UREA
        
        Big_array_out(:,:,16)= 0.0221;%squeeze(nc{'WC_CDON1'}(mytime,:,:)); %CDON1
        Big_array_out(:,:,17)= 0.0259;%squeeze(nc{'WC_NCDON1'}(mytime,:,:)); %NCDON1
        
        Big_array_out(:,:,18)= 0.0201;%squeeze(nc{'LPON'}(mytime,:,:)); %LPOC; %LPON
        Big_array_out(:,:,19)= 0.0201;%squeeze(nc{'RPON'}(mytime,:,:)); %RPOC
        
        Big_array_out(:,:,20)=0.0;%hisdata(mytime,:,11); %PO4
        
        Big_array_out(:,:,21)=0.0;%hisdata(mytime,:,29); %CDOP1
        Big_array_out(:,:,22)=0.0;%hisdata(mytime,:,30); %NCDOP1
        
        Big_array_out(:,:,23)=0.0;%0.01; %LPOP
        Big_array_out(:,:,24)=0.0;%0.01; %RPOP
        
        Big_array_out(:,:,25)=0.0; %PIP
        Big_array_out(:,:,26)=0.0; %SOD
        Big_array_out(:,:,27)= 9.2271;%squeeze(nc{'DOXG'}(mytime,:,:)); %DO2
        
        Big_array_out(:,:,28:32)=0.0; %rest are zeros
        
        Big_array_out(:,:,33)= 0.2453;%squeeze(nc{'WC_CDOC2'}(mytime,:,:)); %CDOC2
        Big_array_out(:,:,34)= 0.2880;%squeeze(nc{'WC_NCDOC2'}(mytime,:,:)); %NCDOC2
        Big_array_out(:,:,35)= 0.0042;%squeeze(nc{'WC_CDOC3'}(mytime,:,:)); %CDOC3
        Big_array_out(:,:,36)= 0.0049;%squeeze(nc{'WC_NCDOC3'}(mytime,:,:)); %NCDOC3
        
        Big_array_out(:,:,37)= 0.0325;%squeeze(nc{'WC_CDON2'}(mytime,:,:)); %CDON2
        Big_array_out(:,:,38)= 0.0382;%squeeze(nc{'WC_NCDON2'}(mytime,:,:)); %NCDON2
        Big_array_out(:,:,39)= 0.0;%squeeze(nc{'WC_CDON3'}(mytime,:,:)); %CDON3
        Big_array_out(:,:,40)= 0.0;%squeeze(nc{'WC_NCDON3'}(mytime,:,:)); %NCDON3
        
                
        Big_array_out(:,:,41)=0.0;%hisdata(mytime,:,31); %CDOC2
        Big_array_out(:,:,42)=0.0;%hisdata(mytime,:,32); %NCDOC2
        Big_array_out(:,:,43)=0.0;%hisdata(mytime,:,33); %CDOC3
        Big_array_out(:,:,44)=0.0;%hisdata(mytime,:,34); %NCDOC3
        
         Big_array_out(nlayers+1,:,:)=Big_array_out(nlayers,:,:);

       % Big_array_out(:,nlayers+1,:)=Big_array_out(:,nlayers,:);

%%
fid = fopen('tonic_initial_wq_vert.dat','w');

for igrid = 1:ngrids
    for jlay=1:nlayers+1
        for kvars = 1:nvars
        fprintf(fid,'%8.4f',Big_array_out(jlay,igrid,kvars));
        end
        fprintf(fid,'\n');
    end
end


