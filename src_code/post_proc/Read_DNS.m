%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx
%
% filename - grid.data file to read DNS grid from
% resultfile_dir - Location of result files 
% ngx,ngy,ngz - number of grid points (cell vertices)
% Lx,Ly,Lz - Domain size
% nvar - number of variable in Subdomain snapshot files
% cc - cell centred data logical flag (true or flase)
%==========================================================================

function[u,v,w,p,stress] = Read_DNS(filename,resultfile_dir,ngx,ngy,ngz,Lx,Ly,Lz,nvar,cc)

pdir = pwd;

%Default filename
if (exist('filename') == 0)
    filename = 'grid.data';
end

%Only return 3 velocity components unless specified
if (exist('nvar','var') == 0)
    nvar = 3;
end

%--- grid size ----
if (exist('ngz','var') == 0)
    error('Read_DNS error - Grid data ngx,ngy or ngz missing')
end

%--- domain size ----
if (exist('Lz','var') == 0)
    error('Read_DNS error - Domain size data Lx, Ly or Lz missing')
end

%Axis for plots
x = linspace(0, Lx, ngx);
y = linspace(0, Ly, ngy);
z = linspace(0, Lz, ngz);

xa = linspace(0, Lx+1/ngx, ngx+1);
ya = linspace(0, Ly+1/ngy, ngy+1);
za = linspace(0, Lz+1/ngz, ngz+1);

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results/';
    display(strcat('setting results file to default',resultfile_dir));
end

%Read grid data
cd(resultfile_dir)
fid = fopen(strcat(resultfile_dir,filename),'r','ieee-le.l64');
xpg = fread(fid,[ngx ngy],'double');
ypg = fread(fid,[ngx ngy],'double');
fclose(fid);

%Get list of Subdomain files
files = dir('Sub*'); dt = 0;
filenames = filenames_const_dt(files,dt);

skipk = 1;
skipi = 1;
skipj = 1;

pz = 1:ngz;
px = 1:ngx;
py = 1:ngy;

%  USE 1 LESS THAN ngx, ngy, ngz if cell centred 
if (exist('cc','var'))
   if (cc == true)
        nn = 1;
    else
        nn = 0;
    end
else
    nn = 0;
end

for n = 1:length(filenames)
   
    %Read from DNS files
    V = read_sub(filenames(n).name,ngz,ngx,ngy,pz,px,py,skipk,skipi,skipj,nvar);

    u(:,n) = squeeze(mean(mean(V{1}(1:end-nn,1:end-nn,:),1),2));
%    G3DtoDX(xpg(1:end,1),ypg(1,1:end),z(1:end), ... 
%            permute(V{1}(:,:,2:end-1),[2 3 1]),strcat('./vmd_volumes/DNS',num2str(n),'.dx'),-Lx/2,2.8499,-Lz/2)
    v(:,n) = squeeze(mean(mean(V{2}(1:end-nn,1:end-nn,:),1),2));
    w(:,n) = squeeze(mean(mean(V{3}(1:end-nn,1:end-nn,:),1),2));
    if (nvar >= 4)
        p = V{4}(1:end-nn,1:end-nn,:);
    end
    %Pressure and Stress if required
    if (nvar > 4)
        if (nvar == 13)
            for jxyz=1:3
            for ixyz=1:3
                stress(ixyz,jxyz,:,n) = squeeze(mean(mean(V{4+ixyz+(jxyz-1)*3}(1:end-nn,1:end-nn,:),1),2));
            end
            end
        else
            disp('Error in Read DNS - Stress is 9 components so nvar should be <4 or 13')
        end
    end
end

end