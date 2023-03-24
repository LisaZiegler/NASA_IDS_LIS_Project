clear all; clc
load('node.mat');%load the nodes data to draw the x axis
load('mesh.mat');%load nodes coordinate and depth data
a=mesh.nodexy(node(:),1);%extract the x coordinate from mesh.mat file
depth=mesh.depth(node(:),1);%extract depth data from mesh.met file
depth=0:-10/56:-10; %interpolate the depth ranging from 0 to -10 m, into 57 depth values, which can match the extracted nodes number.
depth=depth';
nc='tst_0.50001.nc';
tp=ncread(nc,'temp');%extract the output data from nc output
tp1=tp(node(:),:,:);%get the temperature data on the selected nodes
%%
tp3=squeeze(tp1(:,:,121));%select the time you want to plot
tp3=tp3';
%%
for i=1:57; %create the salinity data for contour plot, match the salinity dimensions with x, depth coordinate
    x=tp3(1,i);
    y=tp3(10,i);
    if x==y;
        z(1:57,i)=x; % this loop needs to be changed if your cases is more complicated, for example, like if you have some warmer or cooler temperature in the middle layer while the surface and bottom temperature are the same.
    else
    vq=x:((y-x)/55):y;
    k=1*a;
    ss=interp1(a,k,vq);%this equation needs changing if the variable doesn't drop monotonically.
    z(1,i)=x;
    z(57,i)=y;
    for j=1:55
    z(j+1,i)=vq(j);
    end
    end
end
%%
% pcolor(a,depth,z);shading flat
% hold on 
% [C,hfigc] = contour(a,depth,z,[10:2:20]);
% set(hfigc, ...
%     'LineWidth',1.0, ...
%     'Color', [1 1 1]);
subplot(1,2,2)%change the plot order
contourf(a,depth,z,7,'ShowText','on')
caxis([10 16]);
colormap(winter);

%you can change the caxis to change the range of the colorbar
%%
title('Temperature Vertical Distribution/Mtide:0.5m/after 5 days');
%%
%%
c = colorbar;
c.Label.String = 'Temperature/Celcius';
c.Label.FontSize = 16;