% This file is used to create the river boundary forcing for FVCOM2.7
% The model requires daily discharge, temp and salinity. 
% In this case salinity is zero because the river is fresh

% load data file
%T= readtable('sumQ_Wtemp_Dailymean2018.dat');
%save('meandaily_discharge_temp.mat', 'T', '-v7.3');

% PATH = '/Volumes/LaCie/tonicfvcom_casestudy/fvcom2.7fine_grid/River'
load('meandaily_discharge_temp.mat');

% test to see if increasing temperature will improve temperature in river

t_inc = 4;

temp = [];
for i = 1:365
    
    d = T.wtemp(i) + t_inc;
    
    temp(i,1) = d;
end

%start_time = datenum([2018 01 01 00 00 00]); % data and model start time
%t = datenum(T.Date)-start_time;

 t = (0:24:8736)'; % Julian days but in hours


% making forcing for smaller domain with a single river input
 fid = fopen('hr_riv.dat','w');
 for i = 1 : length(temp);
 fprintf(fid,'%8.4f\n',t(i));
 fprintf(fid, '%8.4f\n', T.Q_sum_m3(i));
 fprintf(fid,'%8.4f\n', temp(i));
 fprintf(fid, '%8.4f\n', T.salt(i));
 end

% making forcing for Large model domain with multiple river inputs
 fid = fopen('lis_riv.dat','w');
 for i = 1 : length(t);
 fprintf(fid,'%8.4f\n',t(i));
 fprintf(fid, '%8.4f %8.4f %8.4f\n', T.Q_sum_m3(i), T.Q_sum_m3(i), T.Q_sum_m3(i));
 fprintf(fid,'%8.4f %8.4f %8.4f\n', temp(i), temp(i), temp(i));
 fprintf(fid, '%8.4f %8.4f %8.4f\n', T.salt(i), T.salt(i), T.salt(i));
 end