% cut up and make the forcing for the Nanticoke River

% for use in the DoCoFVCOM model

% B Clark UMCES, HPL, March 2018


% get the dates for this location
unique_dates = stations.variables(3).date;
EE61_variables = stations.variables;
% get mean value for each date

mean_EE61_variables =struct;

% get the mean values for each date
for ivar = 1 : length(EE61_variables);
    
    mean_EE61_variables(ivar).name = EE61_variables(ivar).name;
    
    for i = 1 : length(unique_dates);
        
        my_ids= find(EE61_variables(ivar).date==unique_dates(i));
        
        mean_EE61_variables(ivar).dates(i)=unique_dates(i);
        
        mean_EE61_variables(ivar).values(i)=nanmean(EE61_variables(ivar).value(my_ids));
        
        
    end
    
    mean_EE61_variables(ivar).values(isnan(mean_EE61_variables(ivar).values))=0.0;
end



% now process the data and build a big cumbersome array that has all the
% data for each date


% first get a big array that has time as the 3rd dimension
forestPercent= 0.4; 

CtoN =(4.785*forestPercent+11.273)*0.857; % Mass C to N Ratio from regression Lu et al., 2014
CDOMFrac = (59.85*forestPercent-2.74)/100; % CDOM fraction, regression from Lu et al., 2013
DOM1frac=0.163 ; %labile DOC fraction from Lu et al., 2013 for all watersheds
DOM2frac = 0.3; %fraction of remaining DOM that is semilabile after removing the labile fraction
DOM3frac = 0.70;

%will use the Redfield for POM P:N ratio

PtoN=30.974./(16*14.0067);

% get the DON and fractionate it out, can do the same with DOP as well
% will use the DON data to calculate DOC
DON = mean_EE61_variables(3).values;
CDON=DON*CDOMFrac;
NCDON=DON*(1-CDOMFrac);
CDON1=CDON*DOM1frac;
NCDON1=NCDON*DOM1frac;
CDON2 = (CDON-CDON1)*0.3;
NCDON2=(NCDON-NCDON1)*0.3;
CDON3 = (CDON-CDON1)*0.7;
NCDON3=(NCDON-NCDON1)*0.7;

% Get the DOP and give same fractionation as DON
DOP=mean_EE61_variables(4).values;
CDOP=DOP*CDOMFrac;
NCDOP=DOP*(1-CDOMFrac);
CDOP1=CDOP*DOM1frac;
NCDOP1=NCDOP*DOM1frac;
CDOP2 = (CDOP-CDOP1)*0.3;
NCDOP2=(NCDOP-NCDOP1)*0.3;
CDOP3 = (CDOP-CDOP1)*0.7;
NCDOP3=(NCDOP-NCDOP1)*0.7;

% get the particulate organic nitrogen by difference between totoal organic
% nitrogen and DON
PON = mean_EE61_variables(8).values-DON;
%fractionate by 20 % labile and 80% refractory following Lu et al., 2014
LPON=PON*0.2;
RPON=PON*0.8;

% if(CDON1+CDON2+CDON3~=CDON)
%     write('DON is bad');
% end
% 

% get total CDOC and NCDOC

myRiver_array = zeros(44,1,20); % 44 -> river source

my_dates = unique_dates*24;
%%
%%%%dischargre
myRiver_array(1,:,:)=repmat(mean_EE61_variables(11).values,1,1);
%%%%temperature
myRiver_array(2,:,:)=repmat(mean_EE61_variables(10).values,1,1);
%Salinity
3
%silica
4
%%%%Algae 1 50C to chl, convert to mg l^-^1 and take half for algae 1
myRiver_array(5,:,:)=repmat(mean_EE61_variables(1).values,1,1).*50./1000*0.5; % 
%%%%Algae 2 50C to chl, convert to mg l^-^1 and take half for algae 2
myRiver_array(6,:,:)=repmat(mean_EE61_variables(1).values,1,1).*50./1000*0.5 ;% 
%Algae 3
7
%Small Zooplankton
8
%large Zooplankton
9
%%%%All DOC values take the calculated above DON fractions and the
%%%%regressed CtoN value
myRiver_array(10,:,:)= repmat(CDON1*CtoN,1,1);
%%%%NCDOC1
myRiver_array(11,:,:)= repmat(NCDON1*CtoN,1,1);
%%%%LPOC
myRiver_array(12,:,:)=repmat(LPON*CtoN,1,1);%uses same regression as the DOC and DON
%%%%RPOC
myRiver_array(13,:,:)=repmat(RPON*CtoN,1,1);
%%%%NH4
myRiver_array(14,:,:)=repmat(mean_EE61_variables(5).values,1,1);
%%%%NO3
myRiver_array(15,:,:)=repmat(mean_EE61_variables(6).values,1,1);
% Urea
16
%%%%CDON1
myRiver_array(17,:,:)=repmat(CDON1,1,1);
%%%%NCDON1
myRiver_array(18,:,:)=repmat(NCDON1,1,1);
%%%%LPON
myRiver_array(19,:,:)=repmat(LPON,1,1);
%%%%RPON 
myRiver_array(20,:,:)=repmat(RPON,1,1);
%%%%PO4
myRiver_array(21,:,:)=repmat(mean_EE61_variables(7).values,1,1);
%%%%CDOP1
myRiver_array(22,:,:)=repmat(CDOP1,1,1);
%%%%NCDOP1
myRiver_array(23,:,:)=repmat(NCDOP1,1,1);
%%%%LPOP
myRiver_array(24,:,:)=repmat(LPON*PtoN,1,1); % convert measured PON to POP using Redfield
%%%%RPOP
myRiver_array(25,:,:)=repmat(RPON*PtoN,1,1); % convert measured PON to POP using Redfield
%PIP
26
%COD
27
%%%%DOXG
myRiver_array(28,:,:)=repmat(mean_EE61_variables(2).values,1,1)
%SIUPB
29
%SIAT
30
%PIB1
31
%PIB2
32
%PIB3
33
%%%%CDOC2
myRiver_array(34,:,:)=repmat(CDON2*CtoN,1,1);
%%%%NCDOC2
myRiver_array(35,:,:)=repmat(NCDON2*CtoN,1,1);
%%%%CDOC3
myRiver_array(36,:,:)=repmat(CDON3*CtoN,1,1);
%%%%NCDOC3
myRiver_array(37,:,:)=repmat(NCDON3*CtoN,1,1);
%%%%CDON2
myRiver_array(38,:,:)=repmat(CDON2,1,1);
%%%%NCDON2
myRiver_array(39,:,:)=repmat(NCDON2,1,1);
%%%%CDON3
myRiver_array(40,:,:)=repmat(CDON3,1,1);
%%%%NCDON3
myRiver_array(41,:,:)=repmat(NCDON3,1,1);
%%%%CDOP2
myRiver_array(42,:,:)=repmat(CDOP2,1,1);
%%%%NCDOP2
myRiver_array(43,:,:)=repmat(NCDOP2,1,1);
%%%%CDOP3
myRiver_array(44,:,:)=repmat(CDOP3,1,1);
%%%%NCDOP3
myRiver_array(45,:,:)=repmat(NCDOP3,1,1);




%NOW CAN WRITE THE RIVER FILE OUTPUT
[r,c,t]=size(myRiver_array);


fid=fopen('Blackwater_river_pnt.dat','w');

for i = 1: t;
    fprintf(fid,'%8.4f',my_dates(i));
    fprintf(fid,'\n');
    for j= 1 : r ;
        fprintf(fid,'%16.8f',squeeze(myRiver_array(j,:,i)));
        fprintf(fid,'\n');
    end


end
%%
%%%%temperature
myRiver_array(1,:,:)=repmat(mean_EE61_variables(10).values,1,1);
%Salinity
2
%silica
3
%%%%Algae 1 50C to chl, convert to mg l^-^1 and take half for algae 1
myRiver_array(4,:,:)=repmat(mean_EE61_variables(1).values,1,1).*50./1000*0.5; % 
%%%%Algae 2 50C to chl, convert to mg l^-^1 and take half for algae 2
myRiver_array(5,:,:)=repmat(mean_EE61_variables(1).values,1,1).*50./1000*0.5 ;% 
%Algae 3
6
%Small Zooplankton
7
%large Zooplankton
8
%%%%All DOC values take the calculated above DON fractions and the
%%%%regressed CtoN value
myRiver_array(9,:,:)= repmat(CDON1*CtoN,1,1);
%%%%NCDOC1
myRiver_array(10,:,:)= repmat(NCDON1*CtoN,1,1);
%%%%LPOC
myRiver_array(11,:,:)=repmat(LPON*CtoN,1,1);%uses same regression as the DOC and DON
%%%%RPOC
myRiver_array(12,:,:)=repmat(RPON*CtoN,1,1);
%%%%NH4
myRiver_array(13,:,:)=repmat(mean_EE61_variables(5).values,1,1);
%%%%NO3
myRiver_array(14,:,:)=repmat(mean_EE61_variables(6).values,1,1);
% Urea
15
%%%%CDON1
myRiver_array(16,:,:)=repmat(CDON1,1,1);
%%%%NCDON1
myRiver_array(17,:,:)=repmat(NCDON1,1,1);
%%%%LPON
myRiver_array(18,:,:)=repmat(LPON,1,1);
%%%%RPON 
myRiver_array(19,:,:)=repmat(RPON,1,1);
%%%%PO4
myRiver_array(20,:,:)=repmat(mean_EE61_variables(7).values,1,1);
%%%%CDOP1
myRiver_array(21,:,:)=repmat(CDOP1,1,1);
%%%%NCDOP1
myRiver_array(22,:,:)=repmat(NCDOP1,1,1);
%%%%LPOP
myRiver_array(23,:,:)=repmat(LPON*PtoN,1,1); % convert measured PON to POP using Redfield
%%%%RPOP
myRiver_array(24,:,:)=repmat(RPON*PtoN,1,1); % convert measured PON to POP using Redfield
%PIP
25
%COD
26
%%%%DOXG
myRiver_array(27,:,:)=repmat(mean_EE61_variables(2).values,1,1)
%SIUPB
28
%SIAT
29
%PIB1
30
%PIB2
31
%PIB3
32
%%%%CDOC2
myRiver_array(33,:,:)=repmat(CDON2*CtoN,1,1);
%%%%NCDOC2
myRiver_array(34,:,:)=repmat(NCDON2*CtoN,1,1);
%%%%CDOC3
myRiver_array(35,:,:)=repmat(CDON3*CtoN,1,1);
%%%%NCDOC3
myRiver_array(36,:,:)=repmat(NCDON3*CtoN,1,1);
%%%%CDON2
myRiver_array(37,:,:)=repmat(CDON2,1,1);
%%%%NCDON2
myRiver_array(38,:,:)=repmat(NCDON2,1,1);
%%%%CDON3
myRiver_array(39,:,:)=repmat(CDON3,1,1);
%%%%NCDON3
myRiver_array(40,:,:)=repmat(NCDON3,1,1);
%%%%CDOP2
myRiver_array(41,:,:)=repmat(CDOP2,1,1);
%%%%NCDOP2
myRiver_array(42,:,:)=repmat(NCDOP2,1,1);
%%%%CDOP3
myRiver_array(43,:,:)=repmat(CDOP3,1,1);
%%%%NCDOP3
myRiver_array(44,:,:)=repmat(NCDOP3,1,1);




%NOW CAN WRITE THE RIVER FILE OUTPUT
[r,c,t]=size(myRiver_array);


fid=fopen('tonic0000_riv_chldec.dat','w');
for i = 1: t;
    fprintf(fid,'%8.4f',my_dates(i));
    fprintf(fid,'\n');
    for j= 1 : r ;
        fprintf(fid,'%16.8f',squeeze(myRiver_array(j,:,i)));
        fprintf(fid,'\n');
    end

end




