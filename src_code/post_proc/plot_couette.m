%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx
%==========================================================================

clear all
close all
fclose('all')

%--- grid size ----
ngx = 8+1;
ngy = 7+1;
ngz = 8+1;

%--- domain size ----
Lx = 1.;
Ly = 1.;
Lz = 1.;

cd '/home/es205/codes/coupled/CFD_dCSE/src_code/results'

%--- Read grid ----
fid = fopen('grid.data','r','ieee-le.l64');
xpg = fread(fid,[ngx ngy],'double');
ypg = fread(fid,[ngx ngy],'double');
fclose(fid);

figure
hold all

files = dir('Sub*');

dt = 40;
filenames = filenames_const_dt(files,dt);

Nskip = 20; Nsteps=10000;
for ntime=0:Nskip:Nsteps
    if(ntime<=9                          )
        file_prefix= 'SubDom_dble.00000';
    elseif(ntime>=10     && ntime<=99    )
         file_prefix='SubDom_dble.0000';       
    elseif(ntime>=100    && ntime<=999   )
        file_prefix= 'SubDom_dble.000';
    elseif(ntime>=1000   && ntime<=9999  )
        file_prefix= 'SubDom_dble.00';
    elseif(ntime>=10000  && ntime<=99999 )
        file_prefix= 'SubDom_dble.0';
    elseif(ntime>=100000 && ntime<=999999)
        file_prefix= 'SubDom_dble.';
    end

    %---- Read in velocity field ----
    filename=strcat(file_prefix, num2str(ntime));
    fid = fopen(filename,'r','ieee-le.l64');
    %Try either side of current value
    if (fid == -1) 
        filename=strcat(file_prefix, num2str(ntime+1));
        fid = fopen(filename,'r','ieee-le.l64');
    end
    if (fid == -1) 
        filename=strcat(file_prefix, num2str(ntime-1));
        fid = fopen(filename,'r','ieee-le.l64');
    end
    if (fid == -1) 
        errorstr=strcat('Error, filename --', ...
                        file_prefix, num2str(ntime-1),'/', ...
                        file_prefix, num2str(ntime  ),'/', ...
                        file_prefix, num2str(ntime+1), ...
                        '  Not found in code output'                 );
        continue
        %error(errorstr)
    end
    uc = fread(fid, (ngz+1)*(ngx+2)*(ngy+1),'double');
    %vc = fread(fid, (ngz+1)*(ngx+1)*(ngy+2),'double');
    %wc = fread(fid, (ngz+2)*(ngx+1)*(ngy+1),'double');
    fclose(fid);

    %--- Reshape velocities and fluxes into matrices ---
    uc = reshape(uc, ngz+1, ngx+2, ngy+1);
    %vc = reshape(vc, ngz+1, ngx+1, ngy+2);
    %wc = reshape(wc, ngz+2, ngx+1, ngy+1);

	plot(squeeze(mean(mean(uc(:,:,1:end-2),1),2)))
    drawnow
end

