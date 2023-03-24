% Converting sms to mat

meshfile = input('What is the mesh.mat file you want FVCOM to use? ----->  ');
load(meshfile);

% Make tonnic_cor forcing file
fid = fopen('lis_cor.dat','w');
for i = 1 : length(lld_n)
fprintf(fid,'%16.4f',xyd_n(i,2:3),lld_n(i,3))
fprintf(fid,'\n');
end

% Plot grid depth 
figure
tricolor(e2n(:,2:4),xyd_n(:,2),xyd_n(:,3), -h_n)
shading flat 

% Make obc forcing file
% in file tonic_obc.dat

%global_node_number  node_number(1:end)  type (use 1)  ini_temp ini_sal
 
type = 1;
ini_temp = 2.3;
ini_sal = 28;
salt = [];
temp = [];
type1 =[];
for i=1:47;
    temp(i,1)=ini_temp;
    salt(i,1) = ini_sal;
    type1(i,1)=type;
end

fid = fopen('lis_obc.dat','w');
for i = 1 : 47;
fprintf(fid,'%8.f %8.f %8.4f %8.4f %8.4f', [i i], type1(i), temp(i), salt(i));
fprintf(fid, '\n');
end




