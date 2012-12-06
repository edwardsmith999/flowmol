%=============================================================================
% Load a DNS grid.data file and output xpg, ypg, xsg, ysg, ngx, ngy, nsx, nsy
%=============================================================================
% Usage:
%
% read_grid(filename,skip,option)
%
% read_grid(filename,skip,'plot') reads the grid file from the specified
% directory and plots the grid
%
% If no grid file name is specified the default name grid.data is used
% Either of the following are valid:
% read_grid([],[2 2 2])
% read_grid('',[2 2 2])
%
% similarly if no skip is specified the default of [2 2 2] is used
% read_grid; with no inputs will read the file grid.data with skip [2 2 2]

function read_grid(filename,skip,varargin)

% WD = cd;
show_plot = 0;

if nargin == 0
    filename = [];
    skip = [2 2 2];
end

if isempty(filename)
filename = 'grid.data';
end

disp(['Reading ',filename ,' in the current directory using kskip = ',num2str(skip(1)),...
                                                          ' iskip = ',num2str(skip(2)),...
                                                          ' jskip = ',num2str(skip(3))])

% plot grid on exit
if ~isempty(varargin)
    if strcmp(varargin{1}, 'plot')
        show_plot = 1;
    else
        error('Too many input arguements. See ''help read_grid''')
    end
end

% Close any open files
fclose('all');

% Stored snapshot grid start point
istart = 1;
jstart = 1;

% Stored snapshot grid skip
% iskip = 2;
% jskip = 2;
% kskip = 2;
iskip = skip(2);
jskip = skip(3);
kskip = skip(1);

% Read grid file
fid = fopen(filename,'r','ieee-le.l64');

switch fid
    case -1
        error(['Specified grid file ''',filename,''' not found in current directory'])
end
xy_grid = fread(fid,inf,'double');
fclose(fid);

% Split into x and y components
x_grid = xy_grid(1:end/2);
y_grid = xy_grid(end/2+1:end);

ngy = sum(x_grid == max(x_grid));
ngx = length(x_grid)./ngy;
ngz = ngy;

% Stored grid size
nsx = (ngx-1)/iskip + 1;
nsy = (ngy-1)/jskip + 1;

xpg = reshape(x_grid,ngx,ngy);
ypg = reshape(y_grid,ngx,ngy);

xsg = xpg(istart:iskip:end,jstart:jskip:end);
ysg = ypg(istart:iskip:end,jstart:jskip:end);


%   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% 	! calculate the location of p points in physical domain
% 	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

alz = 10; % Z-width of the domain

    for k=1:ngz
        zeta=double(k-1);
        delta=2.0;
        x1=delta/2.*(zeta/double(ngz-1)-1.);
        x2=delta/2.;
        tanh1=(exp(x1)-exp(-x1))/(exp(x1)+exp(-x1));
        tanh2=(exp(x2)-exp(-x2))/(exp(x2)+exp(-x2));
%     !	coef3=1.+tanh1/tanh2
        coef3=zeta/double(ngz-1);
        zpw(k)=(1.0-coef3)*0.0+coef3*alz;
    end

	for k=1:ngz
		zpg(k)=zpw(k);
    end

	for k=1:ngz-1 
		zpp(k)=0.5*(zpw(k)+zpw(k+1));
		zpu(k)=0.5*(zpw(k)+zpw(k+1));
		zpv(k)=0.5*(zpw(k)+zpw(k+1));
    end

% 	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% 	! calculate the location of u points in physical domain
% 	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	for j=1:ngy-1
	for i=1:ngx 
		xpu(i,j)=0.5*(xpg(i,j)+xpg(i,j+1));
		ypu(i,j)=0.5*(ypg(i,j)+ypg(i,j+1));
    end
    end

% 	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% 	! calculate the location of v points in physical domain
% 	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	for j=1:ngy
	for i=1:ngx-1
		xpv(i,j)=0.5*(xpg(i,j)+xpg(i+1,j));
		ypv(i,j)=0.5*(ypg(i,j)+ypg(i+1,j));
    end
    end

	for j=1:ngy-1
	for i=1:ngx-1
		xpp(i,j)=0.5*(xpu(i,j)+xpu(i+1,j));
		ypp(i,j)=0.5*(ypv(i,j)+ypv(i,j+1));
    end
    end

	for j=1:ngy-1
	for i=1:ngx-1
		xpw(i,j)=xpp(i,j);
		ypw(i,j)=ypp(i,j);
    end
    end
    
    % control volume size in z
	pzk=alz/double(ngz-1);
    
% !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% ! calculate the area vector for the p control volume
% ! faces (e,w,n,s,f,b): spxi(e,w),speta(n,s), spz(f,b)
% !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

% 	!cccccccccccccccccccccccccc
% 	! for point e,w
% 	!cccccccccccccccccccccccccc
	for j=1:ngy-1
	for i=1:ngx
		petax=xpg(i,j+1)-xpg(i,j);
		petay=ypg(i,j+1)-ypg(i,j);
		spxix(i,j) = petay*pzk;
		spxiy(i,j) =-petax*pzk;
    end
    end
% 	!cccccccccccccccccccccccccc
% 	! for point n,s  
% 	!cccccccccccccccccccccccccc
	for j=1:ngy
	for i=1:ngx-1
		pxix=xpg(i+1,j)-xpg(i,j);
		pxiy=ypg(i+1,j)-ypg(i,j);
		spetax(i,j)=-pzk*pxiy;
		spetay(i,j)= pzk*pxix;
    end
    end
% 	!cccccccccccccccccccccccccc
% 	! for point f,b
% 	!cccccccccccccccccccccccccc
	for j=1:ngy-1   
	for i=1:ngx-1
		pxix =(xpg(i+1,j)+xpg(i+1,j+1))/2.-(xpg(i,j)+xpg(i,j+1))/2.;
		pxiy =(ypg(i+1,j)+ypg(i+1,j+1))/2.-(ypg(i,j)+ypg(i,j+1))/2.;
		petax=(xpg(i,j+1)+xpg(i+1,j+1))/2.-(xpg(i,j)+xpg(i+1,j))/2.;
		petay=(ypg(i,j+1)+ypg(i+1,j+1))/2.-(ypg(i,j)+ypg(i+1,j))/2.;
		spz(i,j)=pxix*petay-pxiy*petax;
    end
    end

% 	for j=1,ngy-1;
% 		spz(0,j)  =spz(1    ,j);
% 		spz(ngx,j)=spz(ngx-1,j);
%     end
% 	for i=0,ngx;
% 		spz(i,0)  =spz(i,1);
% 		spz(i,ngy)=spz(i,ngy-1);
%     end
     
    

%   !ccccccccccccccccccccccccccccccccccccccccccccccccc
% 	! calcuate the volume size for p control volume
% 	!ccccccccccccccccccccccccccccccccccccccccccccccccc
	for j=1:ngy-1
	for i=1:ngx-1
% 		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% 		! first find the positions of fne and bsw for p control volume
% 		! calculate rfne-rbsw
% 		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		rdifx=xpg(i+1,j+1)-xpg(i,j);
		rdify=ypg(i+1,j+1)-ypg(i,j);
		rdifz=pzk;
% 		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% 		! calculate spxiw+spseta+spbz
% 		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		sumx=spxix(i,j)+spetax(i,j);
		sumy=spxiy(i,j)+spetay(i,j);
		sumz=spz(i,j);
		sumx=sumx/3.;
		sumy=sumy/3.;
		sumz=sumz/3.;
% 		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% 		! perform inner product to obtain Vp
% 		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		vp(i,j)=sumx*rdifx + sumy*rdify + sumz*rdifz;
    end
    end
    
assignin('caller','xpg',xpg)
assignin('caller','ypg',ypg)
assignin('caller','xsg',xsg)
assignin('caller','ysg',ysg)

assignin('caller','ngx',ngx)
assignin('caller','ngy',ngy)
assignin('caller','nsx',nsx)
assignin('caller','nsy',nsy)
assignin('caller','vpg' ,vp )

% Plot grid
if show_plot
%figure
plot(xpp,ypp,'ro');hold all
plot(xpg,ypg,'k.-')
plot(xpg',ypg','k.-')

plot(xpu,ypu,'>b')
plot(xpv,ypv,'^b')
end

%Check cell ratio in y
delta_y = diff(ypg(floor(ngx/2.),:));
for i = 1:ngy-2
    ratio(i)=(delta_y(i)/delta_y(i+1));
end
maxratio = max(ratio);
if (maxratio < 1) 
    maxratio = min(ratio);
    maxratio = 100*(1/maxratio-1);
else
    maxratio = (maxratio-1)*100;
end
if (maxratio > 3) 
   disp(strcat('**WARNING** Ratio of cell sizes in y direction is greater than 3%. Largest ratio = ',num2str(maxratio),'%'))
end

end
