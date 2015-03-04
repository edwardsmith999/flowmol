%=======================================
%    Contour plots of uc, vc and wc
%=======================================
function [] = contour_plots()
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
end