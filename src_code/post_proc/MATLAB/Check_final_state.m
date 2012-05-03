clear all
close all
fid = fopen('/home/es205/codes/issued_codes/svn_lucian/MD_dCSE/src_code/results/final_state','r','n');

a = fread(fid,'double');
b = a(1:2048*3*2);
c=reshape(b,[3,2,2048]);
scatter3(c(1,1,:),c(2,1,:),c(3,1,:));
hold on
quiver3(c(1,1,:),c(2,1,:),c(3,1,:),c(1,2,:),c(2,2,:),c(3,2,:));
T = mean(c(1,2,:).^2+c(2,2,:).^2+c(3,2,:).^2)