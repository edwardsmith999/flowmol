clear all
close all
fig1 = figure


% =================================
% READ COUPLED velocity streamlines
% =================================
resultfile_dir = '/home/es205/results/CPL_runs/1NCER_bump_allzthermostat/results/no_bottom_thermostat/'
resultfile_dir_MD = strcat(resultfile_dir,'md_data/results/');
resultfile_dir_CFD = strcat(resultfile_dir,'couette_data/');
resultfile_dir_CPL = strcat(resultfile_dir,'results/');

%Read MD Header
resultfile_dir = resultfile_dir_MD;
read_header
Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));

%Wall Normal Direction
wall_normal = 2;
tstart = 1;
tend = Nvel_records-2;

% = = = Read MD header = = =
%Check output flags and read data accordingly
if (velocity_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
elseif ((velocity_outflag > 0) && (velocity_outflag < 4) )
    read_vslice
elseif (mass_outflag == 4)
    mass_bins = read_mbins('mbins',resultfile_dir_MD);
elseif ((mass_outflag > 0) && (mass_outflag < 4) )
    read_mslice
end
%%MD domain set-up
Domain_setup
MD_domain = globaldomain;
MD_liquiddomain = liquiddomain;
MD_gnbins = gnbins;
MD_binsize = binsize;
clear liquiddomain gnbins binsize

%---Get CFD grid size ----
[ngx, ngy, ngz, Lx, Ly, Lz, dx, dy, dz] = read_report(strcat(resultfile_dir_CFD,'report'));

%Continuum_Domain_setup
read_grid(strcat(resultfile_dir_CFD,'grid.data'),[1 1 1])
CFD_domain = [Lx,Ly,Lz];
CFD_binsize= [dx,dy,dz];
CFD_gnbins = [ngx,ngy,ngz]-1;
clear Lx Ly Lz dx dy dz ngx ngy ngz

%Read Coupled data
resultfile_dir = resultfile_dir_CPL;
CPL_read_header
CPL_olap_nbins = [icmax_olap-icmin_olap+1, jcmax_olap-jcmin_olap+1,kcmax_olap-kcmin_olap+1];
CPL_olap_domain= CFD_binsize.*CPL_olap_nbins;

%Get cell ratio
MD_cells_per_CFD = CFD_binsize./MD_binsize;
coupleddomain = CFD_domain + MD_domain - CPL_olap_domain;
%Check domain size agrees
domainratio  = CFD_domain./MD_domain;
CFD_nbinswallbot= round(nbinswallbot./MD_cells_per_CFD);

%Get orthogonal direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    %[null,ixyz]=max(domainratio);
    ixyz = 2;
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
end


m = 0;
for i = 1:Nvel_records-2
    i;
    m = m + 1

    %Read velocity
    filename = strcat(resultfile_dir_MD,'/mbins');
    mass_bins = read_mbins(filename,resultfile_dir_MD,m);
    filename = strcat(resultfile_dir_MD,'/vbins');
    vel_bins = read_vbins(filename,resultfile_dir_MD,m);
    

        % ====== COUPLED VELOCITY ======
        %Load CFD data
        [xi,yi,zi] = meshgrid(1:size(vel_bins,2), ...
                              1:size(vel_bins,1), ...
                              1:size(vel_bins,3));
        [xrange,yrange,zrange] = meshgrid(1:MD_cells_per_CFD(2):size(vel_bins,2), ...
                                          1:MD_cells_per_CFD(1):size(vel_bins,1), ... 
                                          1:MD_cells_per_CFD(3):size(vel_bins,3));
        [ngx, ngy, ngz, Lx, Ly, Lz, dx, dy, dz] = read_report(strcat(resultfile_dir_CFD,'report'));
        [u,v,w,P] = Read_DNS_Subdomain(m,'grid.data',resultfile_dir_CFD,ngx-2,ngy-1,ngz-2,1,1,1,4,true);
        
        %Coarse grain MD data and combine with CFD grid -- note CFD data is assumed to include halo
        u_MD = interp3(xi,yi,zi,vel_bins(:,:,:,1)./mass_bins(:,:,:),xrange,yrange,zrange); 
        u_CPL = cat(2,u_MD(:,1:end-CPL_olap_nbins(2)/2,:),permute(u(:,:,(CPL_olap_nbins(2)/2+2):end-1),[2,3,1]));
        v_MD = interp3(xi,yi,zi,vel_bins(:,:,:,2)./mass_bins(:,:,:),xrange,yrange,zrange); 
        v_CPL = cat(2,v_MD(:,1:end-CPL_olap_nbins(2)/2,:),permute(v(:,:,(CPL_olap_nbins(2)/2+2):end-1),[2,3,1]));
        w_MD = interp3(xi,yi,zi,vel_bins(:,:,:,3)./mass_bins(:,:,:),xrange,yrange,zrange); 
        w_CPL = cat(2,w_MD(:,1:end-CPL_olap_nbins(2)/2,:),permute(w(:,:,(CPL_olap_nbins(2)/2+2):end-1),[2,3,1]));
               
        cumulative_vel_CPL(:,:,:,1,m) = u_CPL;
        cumulative_vel_CPL(:,:,:,2,m) = v_CPL;
        cumulative_vel_CPL(:,:,:,3,m) = w_CPL;
end

% =================================
% READ Full MD velocity streamlines
% =================================
clearvars -except cumulative_vel_CPL MD_cells_per_CFD wall_normal fig1 tstart tend u_CPL v_CPL x_CPL y_CPL coupleddomain
resultfile_dir = '/home/es205/results/CPL_runs/0NCER_bump_MD_only/51p2_MD_unit_domain__Matched_to_CPL_Case/z_thermostat/results/results/'
cd(resultfile_dir)

%Read MD Header file
read_header

% = = = Read MD header = = =
%Check output flags and read data accordingly
if (velocity_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
elseif ((velocity_outflag > 0) && (velocity_outflag < 4) )
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
elseif (mass_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nmass_ave));
    mass_bins = read_mbins('mbins',resultfile_dir_md);
elseif ((mass_outflag > 0) && (mass_outflag < 4) )
    read_mslice
end

if (tend ~= Nvel_records)
    disp(['Warning coupled records number ', num2str(tend), ...
          ' while uncoupled records only ',  num2str(Nvel_records)])
   tend =  Nvel_records-2;
end


%%MD domain set-up
Domain_setup

%Get orthogonal direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    ixyz = wall_normal;
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
end



m = 0;
for i = 1:Nvel_records-2
    i;
    m = m + 1

    %Read velocity
    filename = strcat(resultfile_dir,'/mbins');
    mass_bins = read_mbins(filename,resultfile_dir,m);
    filename = strcat(resultfile_dir,'/vbins');
    vel_bins = read_vbins(filename,resultfile_dir,m);
    

        % ====== COUPLED VELOCITY ======
        %Load CFD data
        [xi,yi,zi] = meshgrid(1:size(vel_bins,2), ...
                              1:size(vel_bins,1), ...
                              1:size(vel_bins,3));
        [xrange,yrange,zrange] = meshgrid(1:MD_cells_per_CFD(2):size(vel_bins,2), ...
                                          1:MD_cells_per_CFD(1):size(vel_bins,1), ... 
                                          1:MD_cells_per_CFD(3):size(vel_bins,3));

        %Coarse grain MD data and combine with CFD grid
        u_MD = interp3(xi,yi,zi,vel_bins(:,:,:,1)./mass_bins(:,:,:),xrange,yrange,zrange); 
        v_MD = interp3(xi,yi,zi,vel_bins(:,:,:,2)./mass_bins(:,:,:),xrange,yrange,zrange); 
        w_MD = interp3(xi,yi,zi,vel_bins(:,:,:,3)./mass_bins(:,:,:),xrange,yrange,zrange); 
               
        cumulative_vel_MD(:,:,:,1,m) = u_MD;
        cumulative_vel_MD(:,:,:,2,m) = v_MD;
        cumulative_vel_MD(:,:,:,3,m) = w_MD;
end




% =================================
% PLOT ALL velocity streamlines
% =================================

%COUPLED 
u_CPL = mean(cumulative_vel_CPL(:,:,:,1,tstart:tend),5);
v_CPL = mean(cumulative_vel_CPL(:,:,:,2,tstart:tend),5);
w_CPL = mean(cumulative_vel_CPL(:,:,:,3,tstart:tend),5);

xaxis   = linspace(0,coupleddomain(1),size(u_CPL,1));
yaxis   = linspace(0,coupleddomain(2),size(u_CPL,2));
[x_CPL,y_CPL]   = meshgrid(xaxis,yaxis);

%contourf(x_CPL,y_CPL,mean(u_CPL,3)',cres,'LineStyle','None')
%vmagnitude = sqrt(u_CPL.^2+v_CPL.^2+w_CPL.^2);
%contourf(x_CPL,y_CPL,squeeze(mean(vmagnitude(:,:,:),3))',cres,'LineStyle','None')
U = squeeze(mean(u_CPL,3))';
V = squeeze(mean(v_CPL,3))';

hold on
%h = quiver(x_CPL,y_CPL,U,V);
%set(h,'Color','w'); 

figure(fig1);
yaxis   = linspace(0,coupleddomain(2),2*size(u_MD,2));
[sx,sy] = meshgrid(0,yaxis(8:end));
h = streamline(x_CPL,y_CPL,U,V,sx,sy);
set(h,'Color','k','LineWidth',0.01);  hold off


%MD DOMAIN
u_MD = mean(cumulative_vel_MD(:,:,:,1,tstart:tend),5);
v_MD = mean(cumulative_vel_MD(:,:,:,2,tstart:tend),5);
w_MD = mean(cumulative_vel_MD(:,:,:,3,tstart:tend),5);

xaxis   = linspace(0,globaldomain(1),size(u_MD,1));
yaxis   = linspace(0,globaldomain(2),size(u_MD,2));
[x_MD,y_MD]   = meshgrid(xaxis,yaxis);

%contourf(x_MD,y_MD,mean(u_MD,3)',cres,'LineStyle','None')
%vmagnitude = sqrt(u_MD.^2+v_MD.^2+w_MD.^2);
%contourf(x_MD,y_MD,squeeze(mean(vmagnitude(:,:,:),3))',cres,'LineStyle','None')
U = squeeze(mean(u_MD,3))';
V = squeeze(mean(v_MD,3))';

figure(fig1); hold on
yaxis   = linspace(0,globaldomain(2),2*size(u_MD,2));
[sx,sy] = meshgrid(0,yaxis(8:end));
h = streamline(x_MD,y_MD,U,V,sx,sy);
set(h,'Color','r','LineStyle','--','LineWidth',0.01);  hold off


% = = = Get plot of v profiles per cell = = =

%COUPLED CASE
yaxis = linspace(0,1,size(v_CPL,2));
figure; hold all
%Plot coupled results
colist=hsv(size(v_CPL,1));
for i=1:size(v_CPL,1)
    lgnd{i} = num2str(0.06 + (i-1)*0.11);
    plot(yaxis,(mean( v_CPL(i,:,:),3)),'col',colist(i,:))

end


%MD CASE
yaxis = linspace(0,1,size(v_MD,2));
hold all
%Plot coupled results
colist=hsv(size(v_MD,1));
for i=1:size(v_MD,1)
    lgnd{i} = num2str(0.06 + (i-1)*0.11);
    plot(yaxis,(mean( v_MD(i,:,:),3)),'--','col',colist(i,:))

end


axis([0 1 -0.03 0.03]); 
xlabel('y/H','FontSize',16); ylabel('v/U','FontSize',16)
set(gca,'FontSize',16);
set(gca,'FontName','Times'); box on;
legend(lgnd,'orientation','verticle','location','EastOutside')
savefig('./v_velocity','png','eps')




