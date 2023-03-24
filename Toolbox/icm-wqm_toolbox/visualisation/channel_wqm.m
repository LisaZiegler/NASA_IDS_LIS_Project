% nc=('chn_0001.nc');
% nodexy(:,1)=ncread(nc,'x');
% nodexy(:,2)=ncread(nc,'y');
% depth=ncread(nc,'h');
% trinodes=ncread(nc,'nv');
% temp=ncread(nc,'temp');
% save channel.mat;
%%
load('mesh.mat');
figure(1)
patch('Vertices',mesh.nodexy,'Faces',mesh.trinodes,...
           'facecolor','non','edgecolor','blue');%% plot the triangle grids
hold on 
plot(mesh.nodexy(:,1),mesh.nodexy(:,2),'r.');%%plot the nodes
plot(mesh.nodexy(1569,1),mesh.nodexy(1569,2),'b.')

M=505;%nodes number
KBM1=10;%sigma layers
a= 6;
b=10.02;
n=num2str(i,'%04d');
filename1='DO0005.opt';
simvar0=textread(filename1);
for j=1:KBM1
for k=1:M
simvar(k,j)=simvar0(M*(j-1)+k,1);
end
end
figure
patch('Vertices',nodexy,'Faces',trinodes,'Cdata',simvar(1:M,10),...
'facecolor','interp','edgecolor','interp')
axis equal
colorbar
caxis([a,b])

