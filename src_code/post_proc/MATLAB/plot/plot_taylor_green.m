 Vo = 1; L = pi; res = pi/8; H =0;
 x = -L-H:res:L+H; y = -L-H:res:L+H; z = -L-H:res:L+H;
 [X,Y,Z] = meshgrid(x,y,z);
 u =  Vo.*sin(X).*cos(Y).*cos(Z);
 v = -Vo.*cos(X).*sin(Y).*cos(Z);
 w =  Vo.*cos(X).*cos(Y).*sin(Z);
% 
 curlanalx = -Vo*2*cos(X).*sin(Y).*sin(Z);
 [curlx,curly,curlz] = curl(u,v,w);

%Make object files
[f,v]=isosurface(curlx, 0.2);
vertface2obj(v,f,'TG_xvort.obj')
[f,v]=isosurface(curly, 0.01);
vertface2obj(v,f,'TG_yvort.obj')
[f,v]=isosurface(curlz, 0.2);
vertface2obj(v,f,'TG_zvort.obj')

 
%sliceomatic(curlx); sliceomatic(curly); sliceomatic(curlz); 

% uhat = fft(u);
% E = 0.5*uhat.*conj(uhat)
% loglog(E(1,:,1))
% figure
% slice(u(:,:,:),[15 0 0],[],[1 0 0]);
% figure
%quiver3(X,Y,Z,u,v,w)

clear all
close all
resultfile_dir = '/home/es205/results/md_results/fortran/3D_code/parallel/results/taylor_green/taylor_green_Re100_2/';
read_header
Nrecords = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
%Start
%    mbins = read_mbins('./mbins',resultfile_dir,1);
%     velbins = read_vbins('./vbins',resultfile_dir,1);
%     vbins = squeeze(velbins(:,:,:,1)./(mbins(:,:,:)+1));
%    sliceomatic(vbins);

% %Curl
% mbins = read_mbins(1,'./mbins',resultfile_dir);
% velbins = read_vbins(1,'./vbins',resultfile_dir);
% U = squeeze(velbins(:,:,:,1)./(mbins(:,:,:)+1));
% V = squeeze(velbins(:,:,:,2)./(mbins(:,:,:)+1));
% W = squeeze(velbins(:,:,:,3)./(mbins(:,:,:)+1));
% [curlx,curly,curlz,cav] = curl(U,V,W);
% %Finish
%   mbins = read_mbins('./mbins',resultfile_dir,Nrecords-1);
%   velbins = read_vbins('./vbins',resultfile_dir,Nrecords-1);
%   vbins = squeeze(velbins(:,:,:,1)./(mbins(:,:,:)+1));
%   sliceomatic(vbins);

mbins = read_mbins('./mbins',resultfile_dir,1);
velbins = read_vbins('./vbins',resultfile_dir,1);
ubins = squeeze(velbins(:,:,:,2)./(mbins(:,:,:)+1));
vbins = squeeze(velbins(:,:,:,1)./(mbins(:,:,:)+1));
wbins = squeeze(velbins(:,:,:,3)./(mbins(:,:,:)+1));
skip =16;
[curlxi,curlyi,curlzi] = curl(ubins(1:skip:end,1:skip:end,1:skip:end), ...
                              vbins(1:skip:end,1:skip:end,1:skip:end), ...
                              wbins(1:skip:end,1:skip:end,1:skip:end));

%sliceomatic(curlxi); sliceomatic(curlyi) ; sliceomatic(curlzi)

meansqrt_vorti = mean(mean(mean(curlxi.^2 + curlyi.^2+ curlzi.^2)));
%meansqrt_vorti = mean(mean(mean(abs(curlzi))));

h = figure;
camangle_theta = 0; camangel_phi = 0;

for i=1:1:Nrecords-1
    mbins = read_mbins('./mbins',resultfile_dir,i);
    velbins = read_vbins('./vbins',resultfile_dir,i);
    ubins = squeeze(velbins(:,:,:,2)./(mbins(:,:,:)+1));
    vbins = squeeze(velbins(:,:,:,1)./(mbins(:,:,:)+1));
    wbins = squeeze(velbins(:,:,:,3)./(mbins(:,:,:)+1));

    skip = 16;
    [curlx,curly,curlz] = curl(ubins(1:skip:end,1:skip:end,1:skip:end), ...
                               vbins(1:skip:end,1:skip:end,1:skip:end), ...
                               wbins(1:skip:end,1:skip:end,1:skip:end));

    %sliceomatic(curlx); sliceomatic(curly); sliceomatic(curlz); 

	meansqrt_vort(i) =  mean(mean(mean(curlx.^2 + curly.^2+ curlz.^2)))/meansqrt_vorti;
    %meansqrt_vort(i) = mean(mean(mean(abs(curlz))))/meansqrt_vorti;


 	%figure(h);
 	%isosurface(curly, 0.8);% isosurface(curly, 0.5); isosurface(curlz, 0.5);
 	%hold all
 	%isosurface(curly,-0.8); %isosurface(curly,-0.5); isosurface(curlz,-0.5);
    %view(3); axis vis3d; daspect([1 1 1]); box on; camorbit(camangle_theta,camangel_phi);
    %colormap jet; camlight; lighting gouraud
 	%box on; 

%     uslice = ubins(:,15,:);
%     vslice = vbins(:,15,:);
%     uslice = reshape(uslice,size(ubins,1),size(ubins,3));
%     vslice = reshape(vslice,size(vbins,1),size(vbins,3));
%     skip = 4;
%     quiver(uslice(1:skip:end,1:skip:end),vslice(1:skip:end,1:skip:end));
    %hold on
    %contour(sqrt(uslice.^2 + vslice.^2));
%     pause(0.1)
%     close all
%    drawnow
%    clf

    %hold all
end

figure
plot(meansqrt_vort)


% Get k and eps
for i = 1:Nrecords-1
    mbins = read_mbins('./mbins',resultfile_dir,i);
    velbins = read_vbins('./vbins',resultfile_dir,i);
    vbins(:,:,:,1) = velbins(:,:,:,1)./(mbins(:,:,:)+1);
    vbins(:,:,:,2) = velbins(:,:,:,2)./(mbins(:,:,:)+1);
    vbins(:,:,:,3) = velbins(:,:,:,3)./(mbins(:,:,:)+1);
    %eps = diff(vbins(:,:,:,1,i),1).^2 + diff(vbins(:,:,:,2,i),2).^2 + diff(vbins(:,:,:,3,i),3).^2;
    %slice(velbins(:,:,:,2,i),[15 0 0],[],[1 0 0]);
    %pause();

    k(i) = squeeze( mean(mean(mean(vbins(:,:,:,1).^2 + ...
                                   vbins(:,:,:,2).^2+ ...
                                   vbins(:,:,:,3).^2,1),2),3));
end


delta_t = 0.005;
eps = diff(k)/delta_t;

plot(k)
figure
plot(eps)
figure
cases=[1 floor(Nrecords/16) floor(Nrecords/8) floor(Nrecords/4) floor(Nrecords/3) floor(Nrecords/2) Nrecords-1];
%cases = 1:Nrecords-1;
for i=1:size(cases,2)
	mbins = read_mbins('./mbins',resultfile_dir,cases(i));
	velbins = read_vbins('./vbins',resultfile_dir,cases(i));
	vbins = squeeze(velbins(:,:,:,1)./(mbins(:,:,:)+1));
    vbinshat = fft(vbins(:,:,:,1));
    E = 0.5*vbinshat.*conj(vbinshat);
    loglog(mean(mean(E(1:size(E,1)/2,:,:),2),3))
    hold all
end
alpha = 100;
k53 = alpha*(1:size(E,1)/2).^(-5/3);
%loglog(k53,'k--')

legend(num2str(cases(1)),num2str(cases(2)),num2str(cases(3)),num2str(cases(4)),num2str(cases(5)),num2str(cases(6)),num2str(cases(7)))%,'k^{-5/3}')