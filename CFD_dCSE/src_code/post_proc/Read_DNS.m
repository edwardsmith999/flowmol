%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx
%==========================================================================

function[u,v,w] = Read_DNS(filename,resultfile_dir,ngx,ngy,ngz,Lx,Ly,Lz)

pdir = pwd;

%--- grid size ----
if (exist('ngz','var') == 0)
    ngx = 128+1;
    ngy = 128+1;
    ngz = 8+1;
end

%--- domain size ----
if (exist('Lz','var') == 0)
    Lx = 181.1;
    Ly = 181.1;
    Lz = 12.0;
end

%Axis for plots
x = linspace(0, Lx, ngx);
y = linspace(0, Ly, ngy);
z = linspace(0, Lz, ngz);

xa = linspace(0, Lx+1/ngx, ngx+1);
ya = linspace(0, Ly+1/ngy, ngy+1);
za = linspace(0, Lz+1/ngz, ngz+1);

%Reynolds number
Re = 0.375;

if (exist('filename') == 0)
    filename = 'grid.data';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = '/home/es205/results/MD_continuum_results/code/coupled_couette/varying_processor_study/CFD_dCSE/src_code/results/';
    %resultfile_dir = '/home/es205/codes/coupled/coupler_dCSE/src_code/couette_data/';
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

pz = 1:ngz;
px = 1:ngx;
py = 1:ngy;

% read subdoms one at a time
%figure

m = 1
for n = 1:length(filenames)
    n
    %Analytical solution
    %t = (n-1)*5;
    %analy = couette_analytical_fn(t,Re,[1.0,0],Ly,ngy-1,'top');
    %plot(y,analy,'k');
    %hold on
    
    %axis([-1 58 -0.1 1.1])
    
    %Read from DNS files
    V = read_sub(filenames(n).name,ngz,ngx,ngy,pz,px,py,skipk,skipi,skipj,3);
    u(:,m) = squeeze(mean(mean(V{1},1),2));
    v(:,m) = squeeze(mean(mean(V{2},1),2));
    w(:,m) = squeeze(mean(mean(V{3},1),2));
    %scatter(ypg(4,:)+5.12/2,u(:,m),'s')
    %drawnow
    m = m + 1;
    %hold off
    %savefig(strcat('couette',num2str(m)),'eps')
end

end

% plot(squeeze(mean(mean(uc,1),2)))