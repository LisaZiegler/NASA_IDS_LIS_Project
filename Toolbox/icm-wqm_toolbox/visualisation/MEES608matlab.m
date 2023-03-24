%%read the variable data from the output file, you can use command ncdump -h
%%filename to check the variable names in your nc output file
nc=('chn_0001.nc');
xc=ncread(nc,'xc');
yc=ncread(nc,'yc');
zeta=ncread(nc,'zeta');
temp=ncread(nc,'temp');
sa=ncread(nc,'salinity');
depth=ncread(nc,'h');
xc=ncread(nc,'xc');
yc=ncread(nc,'yc');

%%
load('mesh.mat');
patch('Vertices',mesh.nodexy,'Faces',mesh.trinodes,...
           'facecolor','non','edgecolor','blue');%% plot the triangle grids
hold on 
plot(mesh.nodexy(:,1),mesh.nodexy(:,2),'r.');%%plot the nodes
plot(xc,yc,'g*');%% plot the elements
%%
patch('Vertices',mesh.nodexy,'Faces',mesh.trinodes,'Cdata',depth,...
           'facecolor','interp','edgecolor','interp')
axis equal
colorbar
%%
hold on
sa1=squeeze(sa(:,1,48)); %% salinity data in nc file is in three dimensions, we should change it into the data lower than two dimensions


axis equal
colorbar
%%
temp1=squeeze(temp(:,1,48));
patch('Vertices',mesh.nodexy,'Faces',mesh.trinodes,'Cdata',temp1,...
           'facecolor','interp','edgecolor','interp')
axis equal
colorbar
%%
patch('Vertices',mesh.nodexy,'Faces',mesh.trinodes,'Cdata',zeta(:,48),...
           'facecolor','interp','edgecolor','interp')
axis equal
colorbar
%%
%%we can also plot the time sries data, for example:
plot(zeta(635,:));
sa2=squeeze(sa(635,2,:));
plot(sa2);