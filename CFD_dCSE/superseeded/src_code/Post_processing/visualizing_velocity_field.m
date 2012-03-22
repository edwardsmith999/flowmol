%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx 
%==========================================================================

clear all
close all
fclose('all')

	%--- grid size ----
        ngx = 256+1;
        ngy = 128+1;
        ngz = 4+1;

	%--- domain size ----
	Lx = 100.;
	Ly =  20.;
	Lz =   1.;

        %--- Read grid ----
        fid = fopen('grid.data','r','ieee-le.l64');
        xpg = fread(fid,[ngx ngy],'double');
        ypg = fread(fid,[ngx ngy],'double');
        fclose(fid);


        %---- Read in velocity field ----
        fid = fopen('ucvcwc.dble.002000','r','ieee-le.l64');
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

        %=======================================
        %    Contour plots of uc, vc and wc
        %=======================================

        %---------------------------------------
	%	pick k-plane to visualize
        %---------------------------------------
	kpick = (ngz-1)/2;

	%--- uc ---
        vel = zeros(ngx,ngy);
        vel(:,:)  = uc(kpick,1:ngx,1:ngy);
        figure; 
        	map  = colormap('gray');
        	contourf(xpg,ypg,vel);
        	colorbar
        	xlabel('x'); ylabel('y')
        	title('Mean Streamwise Velocity U - xy plane')

	%--- vc ---
         vel = zeros(ngx,ngy);
         vel(:,:)  = vc(kpick,1:ngx,1:ngy);
         figure; 
         	map  = colormap('gray');
         	n    = length(map);
         	pam(1:n , :) = map(n:-1:1 , :);
         	colormap(pam);
         	contourf(xpg,ypg,vel);
         	colorbar
         	xlabel('x'); ylabel('y')
         	title('Mean Wall-Normal Velocity V - xy plane')
 
         
	%--- wc ---
	 vel = zeros(ngx,ngy);
         vel(:,:)  = wc(kpick,1:ngx,1:ngy);
         figure;
         	map  = colormap('gray');
         	n    = length(map);
         	pam(1:n , :) = map(n:-1:1 , :);
         	colormap(pam);
         	contourf(xpg,ypg,vel);
         	colorbar
         	xlabel('x'); ylabel('y')
         	title('Mean Spanwise Velocity W - xy plane')

        %---------------------------------------
	%	pick i-plane to visualize
        %---------------------------------------
        ipick = (ngx-1)/2;

	!--- uc ---
        vel = zeros(ngz,ngy);
        vel(:,:)  = uc(1:ngz,ipick,1:ngy);
        figure;
        	map  = colormap('gray');
        	contourf(z,y,vel');
        	colorbar
        	xlabel('z'); ylabel('y')
        	title('Mean Streamwise Velocity U - zy plane')

	!--- vc ---
         vel = zeros(ngz,ngy);
         vel(:,:)  = vc(1:ngz,ipick,1:ngy);
         figure;
         	map  = colormap('gray');
         	n    = length(map);
         	pam(1:n , :) = map(n:-1:1 , :);
         	colormap(pam);
         	contourf(z,y,vel');
         	colorbar
         	xlabel('z'); ylabel('y')
         	title('Mean Wall-Normal Velocity V - zy plane')

	!--- wc ---
         vel = zeros(ngz,ngy);
         vel(:,:)  = wc(1:ngz,ipick,1:ngy);
         figure;
         	map  = colormap('gray');
         	n    = length(map);
         	pam(1:n , :) = map(n:-1:1 , :);
         	colormap(pam);
         	contourf(z,y,vel');
         	colorbar
         	xlabel('z'); ylabel('y')
         	title('Mean Spanwise Velocity W - zy plane')


        %---------------------------------------
	%	pick j-plane to visualize
        %---------------------------------------
        jpick = (ngy-1)/2;

	%--- uc ---
        vel = zeros(ngz,ngx);
        vel(:,:)  = uc(1:ngz, 1:ngx, jpick);
        figure;
        	map  = colormap('gray');
        	contourf(z,x,vel');
        	colorbar
        	xlabel('z'); ylabel('x')
        	title('Mean Streamwise Velocity U - zx plane')

	%--- vc ---
         vel = zeros(ngz,ngx);
         vel(:,:)  = vc(1:ngz, 1:ngx, jpick);
         figure;
         	map  = colormap('gray');
         	n    = length(map);
         	pam(1:n , :) = map(n:-1:1 , :);
         	colormap(pam);
         	contourf(z,x,vel');
         	colorbar
         	xlabel('z'); ylabel('x')
         	title('Mean Wall-Normal Velocity V - zx plane')


	%--- wc ---
         vel = zeros(ngz,ngx);
         vel(:,:)  = wc(1:ngz, 1:ngx, jpick);
         figure;
         	map  = colormap('gray');
         	n    = length(map);
         	pam(1:n , :) = map(n:-1:1 , :);
         	colormap(pam);
         	contourf(z,x,vel');
         	colorbar
         	xlabel('z'); ylabel('x')
         	title('Mean Spanwise Velocity W - zx plane')

