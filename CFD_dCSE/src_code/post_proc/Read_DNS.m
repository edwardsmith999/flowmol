%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx
%==========================================================================

function[u,v,w] = Read_DNS(filename,resultfile_dir)

pdir = pwd;

%--- grid size ----
ngx = 7+1;
ngy = 7;
ngz = 8+1;

%--- domain size ----
Lx = 35.90949;
Ly = 35.90949;
Lz = 20.0;

%Axis for plots
x = linspace(0, Lx, ngx);
y = linspace(0, Ly, ngy);
z = linspace(0, Lz, ngz);

xa = linspace(0, Lx+1/ngx, ngx+1);
ya = linspace(0, Ly+1/ngy, ngy+1);
za = linspace(0, Lz+1/ngz, ngz+1);

%Reynolds number
Re = 5.0;

if (exist('filename') == 0)
    filename = 'grid.data';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    %resultfile_dir = '/home/es205/codes/coupled/CFD_dCSE/src_code/results';
    resultfile_dir = '/home/es205/codes/coupled/coupler_dCSE/src_code/couette_data/';
    display('setting results file to default');
end

cd(resultfile_dir)

%--- Read grid ----
%read_grid

fid = fopen(strcat(resultfile_dir,filename),'r','ieee-le.l64');
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
%figure

m = 1
for n = 1:length(filenames)
    
    %Analytical solution
    %t = (n-1)*15;
    %analy = couette_analytical_fn(t,Re,[1.0,0],Ly,ngy-1,'top');
    %plot(y,analy,'r');
    %hold on
    
    %axis([0 60 -0.1 1.1])
    
    %Read from DNS files
    V = read_sub(filenames(n).name,ngz,ngx,ngy,pz,px,py,skipk,skipi,skipj,3);
    u(:,m) = V{1};
    v(:,m) = V{2};
    w(:,m) = V{3};
    %scatter(y,u(:,m),'s')
    %drawnow
    m = m + 1;
    %hold off
end

end

% plot(squeeze(mean(mean(uc,1),2)))