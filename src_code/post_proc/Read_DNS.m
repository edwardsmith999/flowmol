
%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx
%==========================================================================

clear all
close all
fclose('all')

pdir = pwd;

%--- grid size ----
ngx = 32+1;
ngy = 32+1;
ngz = 8+1;

%--- domain size ----
Lx = 1.;
Ly = 52.1;
Lz = 1.;

%Axis for plots
x = linspace(0, Lx, ngx);
y = linspace(0, Ly, ngy);
z = linspace(0, Lz, ngz);

xa = linspace(0, Lx+1/ngx, ngx+1);
ya = linspace(0, Ly+1/ngy, ngy+1);
za = linspace(0, Lz+1/ngz, ngz+1);

%Reynols number
Re = 19.5;

cd '/home/es205/codes/coupled/CFD_dCSE/src_code/results'
%cd '/home/es205/codes/coupled/coupler_dCSE/src_code/couette_data'

%--- Read grid ----
%read_grid
fid = fopen('grid.data','r','ieee-le.l64');
xpg = fread(fid,[ngx ngy],'double');
ypg = fread(fid,[ngx ngy],'double');
fclose(fid);

files = dir('Sub*');

dt = 0;
filenames = filenames_const_dt(files,dt);

skipk = 1;
skipi = 1;
skipj = 1;

pz = 1;
px = 1;
py = 1:ngy;

% read subdoms one at a time
figure
hold on

m = 1
for n = 1:2^2:100 %length(filenames)
    
    %Analytical solution
	t = (n-1)*2.
    analy = couette_analytical_fn(t,Re,[1.0,0],Ly,ngy-1,'top');
    plot(y,analy,'r');
    drawnow

    %Read from DNS files
    V = read_sub(filenames(n).name,ngz,ngx,ngy,pz,px,py,skipk,skipi,skipj,3);
    u(:,m) = V{1};
    v(:,m) = V{2};
    w(:,m) = V{3};
    scatter(y,u(:,m),'s')
    drawnow
    m = m + 1;
end

% plot(squeeze(mean(mean(uc,1),2)))