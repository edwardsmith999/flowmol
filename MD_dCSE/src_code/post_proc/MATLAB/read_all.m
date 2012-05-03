%Read all outputs present in results
clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read simulation properties from header file
read_header;
%Read mass in each bin
mass_bins_s = read_mbins('./mbins_s');
mass_bins_p = read_mbins('./mbins_p');
a = mass_bins_s-mass_bins_p;
for i=1:size(size(mass_bins_s),2)
    a = squeeze(max(a));
end
disp(strcat('error in mass serial vs parallel = ',num2str(a,10)))
%Read velocity in each bin
vel_bins_s = read_vbins('./vbins_s');
vel_bins_p = read_vbins('./vbins_p');
a = vel_bins_s-vel_bins_p;
for i=1:size(size(vel_bins_s),2)
    a = squeeze(max(a));
end
disp(strcat('error in momentum serial vs parallel = ',num2str(a,10)))
%Read VA stress
pressure_VA_s = read_pVA('./pVA_s');
pressure_VA_p = read_pVA('./pVA_p');
a = pressure_VA_s-pressure_VA_p;
for i=1:size(size(pressure_VA_s),2)
    a = squeeze(max(a));
end
disp(strcat('error in VA Pressure serial vs parallel = ',num2str(a,10)))
ixyz = 3; kxyz = mod(ixyz,3)+1;	jxyz = mod(ixyz+1,3)+1
scatter((1:10),mean(mean(mean(pressure_VA_s(:,:,:,ixyz,ixyz),jxyz),kxyz),5),'s')
set(gca,'XTick',0:1:11)
grid on
hold all
scatter((1:10),mean(mean(mean(pressure_VA_p(:,:,:,ixyz,ixyz),jxyz),kxyz),5),'x')

%Read mass flux
[mass_flux_s,mass_snapshot_s]=read_mflux('mflux_s','msnap_s');
[mass_flux_p,mass_snapshot_p]=read_mflux('mflux_p','msnap_p');
%Serial
ixyz = 3; kxyz = mod(ixyz,3)+1;	jxyz = mod(ixyz+1,3)+1
scatter((1:10)-0.5,mean(mean(mean(mass_flux_s(:,:,:,ixyz,:),jxyz),kxyz),5),'bs')
set(gca,'XTick',0:1:11)
grid on
hold all
scatter((1:10)+0.5,-mean(mean(mean(mass_flux_s(:,:,:,ixyz+3,:),jxyz),kxyz),5),'rx')
%Parallel
plot((1:10)-0.5,squeeze(mean(mean(mean(mass_flux_p(:,:,:,ixyz,:),jxyz),kxyz),5)),'--k')
hold all
plot((1:10)+0.5,-squeeze(mean(mean(mean(mass_flux_p(:,:,:,ixyz+3,:),jxyz),kxyz),5)),'c--')

a = mass_flux_s-mass_flux_p;
for i=1:size(size(mass_flux_s),2)
    a = squeeze(max(a));
end
disp(strcat('error in mass flux serial vs parallel = ',num2str(a,10)))
a = mass_snapshot_s-mass_snapshot_p;
for i=1:size(size(mass_snapshot_s),2)
    a = squeeze(max(a));
end
disp(strcat('error in mass snapshot serial vs parallel = ',num2str(a,10)))


%Read momentum flux
b = [0 0 0];
ixyz = 1; kxyz = mod(ixyz,3)+1;	jxyz = mod(ixyz+1,3)+1;
nxyz = 1;
for t=2:(Nsteps-initialstep)
	[velocity_snapshot_s,velocity_flux_s,pressure_surface_s] = ...
        read_vflux(t-1,resultfile_dir,gnbins,nd,'vflux_s','psurface_s','vsnap_s');
	[velocity_snapshot_p,velocity_flux_p,pressure_surface_p] = ...
        read_vflux(t-1,resultfile_dir,gnbins,nd,'vflux_p','psurface_p','vsnap_p');
    a = velocity_snapshot_s-velocity_snapshot_p;
%     ixyz = 2; kxyz = mod(ixyz,3)+1;	jxyz = mod(ixyz+1,3)+1;
%     scatter((1:10),mean(mean(velocity_snapshot_s(:,:,:,ixyz),jxyz),kxyz),'bs')
%     set(gca,'XTick',0:1:11)
%     grid on
%     hold all
%     scatter((1:10),mean(mean(velocity_snapshot_p(:,:,:,ixyz),jxyz),kxyz),'rx')
%     pause()
    for i=1:size(size(velocity_snapshot_s),2)
        a = squeeze(max(a));
    end
    b(1) = max(a,b(1));
    c(1,:,t) = mean(mean(velocity_snapshot_s(:,:,:,ixyz),jxyz),kxyz);
    c(2,:,t) = mean(mean(velocity_snapshot_p(:,:,:,ixyz),jxyz),kxyz);
    a = velocity_flux_s-velocity_flux_p;
    for i=1:size(size(velocity_flux_s),2)
        a = squeeze(max(a));
    end
    b(2) = max(a,b(2));
    c(3,:,t) = mean(mean(mean(velocity_flux_s(:,:,:,nxyz,ixyz  ),jxyz),kxyz),5);
    c(4,:,t) = mean(mean(mean(velocity_flux_s(:,:,:,nxyz,ixyz+3),jxyz),kxyz),5);
    c(5,:,t) = mean(mean(mean(velocity_flux_p(:,:,:,nxyz,ixyz  ),jxyz),kxyz),5);
    c(6,:,t) = mean(mean(mean(velocity_flux_p(:,:,:,nxyz,ixyz+3),jxyz),kxyz),5);

    a = pressure_surface_s-pressure_surface_p;
    for i=1:size(size(pressure_surface_s),2)
        a = squeeze(max(a));
    end
    b(3) = max(a,b(3));
    c(7,:,t) = mean(mean(mean(pressure_surface_s(:,:,:,nxyz,ixyz  ),jxyz),kxyz),5);
    c(8,:,t) = mean(mean(mean(pressure_surface_s(:,:,:,nxyz,ixyz+3),jxyz),kxyz),5);
    c(9,:,t) = mean(mean(mean(pressure_surface_p(:,:,:,nxyz,ixyz  ),jxyz),kxyz),5);
    c(10,:,t)= mean(mean(mean(pressure_surface_p(:,:,:,nxyz,ixyz+3),jxyz),kxyz),5);

% 
%     ixyz = 2; kxyz = mod(ixyz,3)+1;	jxyz = mod(ixyz+1,3)+1;
%     scatter((1:10)-0.5,mean(mean(mean(velocity_flux_s(:,:,:,1,ixyz),jxyz),kxyz),5),'bs')
%     set(gca,'XTick',0:1:11)
%     grid on
%     hold all
%     scatter((1:10)+0.5,mean(mean(mean(velocity_flux_s(:,:,:,1,ixyz+3),jxyz),kxyz),5),'rx')
%     %Parallel
%     plot((1:10)-0.5,squeeze(mean(mean(mean(velocity_flux_p(:,:,:,1,ixyz),jxyz),kxyz),5)),'--k')
%     hold all
%     plot((1:10)+0.5,squeeze(mean(mean(mean(velocity_flux_p(:,:,:,1,ixyz+3),jxyz),kxyz),5)),'c--')
%     pause()
end
disp(strcat('error in velocity_snapshot serial vs parallel = ',num2str(b(1),10)))
disp(strcat('error in velocity_flux_s serial vs parallel = ',  num2str(b(2),10)))
disp(strcat('error in pressure_surface serial vs parallel = ', num2str(b(3),10)))

figure
set(gca,'XTick',0:1:11)
grid on
hold all
plot(1:10, squeeze(mean(c(1,:,:),3)),'rh','MarkerSize',15)
plot(1:10, squeeze(mean(c(2,:,:),3)),'r-')
plot((1:10)-0.5, squeeze(mean(c(3,:,:),3)),'bx')
plot((1:10)+0.5,-squeeze(mean(c(4,:,:),3)),'bs','MarkerSize',10)
plot((1:10)-0.5, squeeze(mean(c(5,:,:),3)),'b--')
plot((1:10)+0.5,-squeeze(mean(c(6,:,:),3)),'b-')
plot((1:10)-0.5, squeeze(mean(c(7,:,:),3)),'k^','MarkerSize',10)
plot((1:10)+0.5, squeeze(mean(c(8,:,:),3)),'k*')
plot((1:10)-0.5, squeeze(mean(c(9,:,:),3)),'k:')
plot((1:10)+0.5, squeeze(mean(c(10,:,:),3)),'k-.')

%Read Energy flux
%read_eflux


