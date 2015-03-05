%close all

xmin =   -2.56496  ;
xmax =    2.57  ;
deltax =  2*0.85499  ;
ymin =   xmin  ;
ymax =   xmax  ;
deltay = deltax ;
zmin =   xmin  ;
zmax =   xmax  ;
deltaz = deltax ; 

[gridx,gridy,gridz]=meshgrid(xmin:deltax:xmax,ymin:deltay:ymax,zmin:deltaz:zmax);

%Setup Grid
for i = 1:max(size(gridx))
    scatter3 (gridz(i,:), gridx(i,:), gridy(i,:), 'r');
    hold on
end

% scatter3(a(1),a(2),a(3))
% scatter3(b(1),b(2),b(3))
% scatter3(c(1),c(2),c(3))
% scatter3(d(1),d(2),d(3))

% quiver3 (a(1),a(2),a(3), -a(1)+b(1),-a(2)+b(2),-a(3)+b(3)); figure(gcf)
