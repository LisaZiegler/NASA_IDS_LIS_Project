%% Plot velocity vector
clear all;
clc;
% Load files

load('../Documents/CT_DEM/Housatonic_Forcings/Meshgrid/Old_mesh/Tonic_26Jul.mat');
nc = 'tonic_0005.nc';

%% Extract variables


lat = ll_e(:,2); % degrees coordinates from mesh info file
lon = ll_e(:,1); % degree coordinates from mesh info file

u = ncread(nc,'u');
v = ncread(nc,'v');
%nv = ncread(nc,'nv');
%x = ncread(nc,'x');
%y = ncread(nc,'y');

% Another way to get lat and lon of elements
%xc = mean(x(nv),2); %decimal degree
%yc = mean(y(nv),2); %decimal degree

latlim = [40.5 41.3];
lonlim = [-73 -72.5];

struc = gshhs('../Documents/MATLAB/Apps/gshhg-bin-2.3.7/gshhs_f.b',latlim,lonlim);
%% Plot 
close all;

scale = 3;
interval = 10;
figure(1)

geoshow(struc,'FaceColor','w');
hold on;
quiver(lon(1:interval:end),lat(1:interval:end),...
    squeeze(u(1:interval:end,1,1)),squeeze(v(1:interval:end,1,1)), scale,...
    'MaxHeadSize',1,'AutoScale','off');






