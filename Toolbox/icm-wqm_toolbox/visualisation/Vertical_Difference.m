clear all; clc
load('node.mat');
load('mesh.mat');
a=mesh.nodexy(node(:),1);
depth=mesh.depth(node(:),1);
depth=0:-10/56:-10;
depth=depth';
%%
nc='tst_0.50001.nc';
tp=ncread(nc,'temp');
tp1=tp(node(:),:,:);
tp3=squeeze(tp1(:,:,121));
tp3=tp3';%read the first file data in
%%
nc2='tst_0001.nc';
tp0=ncread(nc2,'temp');
tp01=tp0(node(:),:,:);
tp03=squeeze(tp01(:,:,121));
tp03=tp03';%read the second file data in
tp4=tp03-tp3;%do the subtraction
%%
for i=1:57; %create the salinity data for contour plot, match the dimensions to x, y coordinate
    x=tp4(1,i);
    y=tp4(10,i);
    if x==y;
        z(1:57,i)=x;
    else
    vq=x:((y-x)/55):y;
    k=1*a;
    ss=interp1(a,k,vq);%%change the loop and equation if the variable doesn't drop or increase motonically
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
% subplot(1,2,2)%change the plot order
contourf(a,depth,z,7,'ShowText','on')
caxis([-1 0.3]);
colormap(autumn);

%%
title('Temperature Vertical Difference/TMtide(1-0.5m)/after 5 days');
%%
%%
c = colorbar;
c.Label.String = 'Temperature/Celcius';
c.Label.FontSize = 16;