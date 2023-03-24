 fid = fopen('tonic_met_good.dat','w');
for i = 1: 366;
fprintf(fid,'%16.8f',tonicmetgrid(i,:));
fprintf(fid,'\n');
end