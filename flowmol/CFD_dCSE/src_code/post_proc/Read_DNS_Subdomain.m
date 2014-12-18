%==========================================================================
% Read VELOCITY and FLUX field from DNS subdomain files 
%
% n - timestep
% filename - grid.data file to read DNS grid from
% resultfile_dir - Location of result files 
% ngx,ngy,ngz - number of grid points (cell vertices)
% Lx,Ly,Lz - Domain size
% nvar - number of variable in Subdomain snapshot files
% cc - cell centred data logical flag (true or flase)
%==========================================================================

function[u,v,w,p,stress] = Read_DNS_Subdomain(n,filename,resultfile_dir,ngx,ngy,ngz,skipk,skipi,skipj,nvar,cc)

spwd = pwd;

%Only return 3 velocity components unless specified
if (exist('nvar','var') == 0)
    nvar = 3;
end

%--- grid size ----
if (exist('ngz','var') == 0)
    error('Read_DNS error - Grid data ngx,ngy or ngz missing')
end

if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results/';
    display(strcat('setting results file to default',resultfile_dir));
end

%Go to results directory
cd(resultfile_dir)

%Get list of Subdomain files
files = dir('Sub*'); dt = 0;
filenames = filenames_const_dt(files,dt);

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

 
%Read from DNS files
V = read_sub(filenames(n).name,ngz,ngx,ngy,pz,px,py,skipk,skipi,skipj,nvar);
u = V{1}(1:end-nn,1:end-nn,:);
v = V{2}(1:end-nn,1:end-nn,:);
w = V{3}(1:end-nn,1:end-nn,:);
if (nvar >= 4)
    p = V{4}(1:end-nn,1:end-nn,:);
end
%Pressure and Stress if required
if (nvar > 4)
    if (nvar == 13)
        for jxyz=1:3
        for ixyz=1:3
            stress(ixyz,jxyz,:,n) = V{4+ixyz+(jxyz-1)*3}(1:end-nn,1:end-nn,:);
        end
        end
    else
        disp('Error in Read DNS - Stress is 9 components so nvar should be =< 4 or 13')
    end
end

cd(spwd);

end