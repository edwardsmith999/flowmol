%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx 
%==========================================================================

clear all
close all
fclose('all')

	%--- grid size ----
        ngx = 32+1;
        ngy = 32+1;
        ngz = 8+1;

	%--- domain size ----
	Lx = 1.;
	Ly = 1.;
	Lz = 1.;

        cd './../results'

        %--- Read grid ----
        fid = fopen('grid.data','r','ieee-le.l64');
        xpg = fread(fid,[ngx ngy],'double');
        ypg = fread(fid,[ngx ngy],'double');
        fclose(fid);


        %---- Read in velocity field ----
        fid = fopen('SubDom_dble.000055','r','ieee-le.l64');
        uc = fread(fid, (ngz+1)*(ngx+2)*(ngy+1),'double');
        vc = fread(fid, (ngz+1)*(ngx+1)*(ngy+2),'double');
        wc = fread(fid, (ngz+2)*(ngx+1)*(ngy+1),'double');
        fclose(fid);


        %--- Reshape velocities and fluxes into matrices ---
        uc = reshape(uc, ngz+1, ngx+2, ngy+1);
        vc = reshape(vc, ngz+1, ngx+1, ngy+2);
        wc = reshape(wc, ngz+2, ngx+1, ngy+1);

	x = linspace(0, Lx, ngx);
	y = linspace(0, Ly, ngy);
	z = linspace(0, Lz, ngz);


contour_plots()

