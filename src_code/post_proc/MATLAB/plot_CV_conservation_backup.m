%Plot all output files from simulation
close all
clear all


ibin = 4; jbin = 4; kbin=4;
resultfile_dir = './../../results/120201_CV_conservation';
%resultfile_dir = './../../results/'

%========CV Mass Conservation=======
read_mflux
%Calculate total CV flux and change in mass
totalmflux = squeeze(sum(mass_flux,4));
for i =1:Nmflux_records
    dmassdt(:,:,:,i) = mass_snapshot(:,:,:,i+1) - mass_snapshot(:,:,:,i);
end
%Verify that CV mass is exactly conservative
display(strcat('Maximum error in CV mass conservation =', ... 
   num2str(max(max(max(max(squeeze(totalmflux(:,:,:,:)) - squeeze(dmassdt(:,:,:,:))))))*...
           min(min(min(min(squeeze(totalmflux(:,:,:,:)) - squeeze(dmassdt(:,:,:,:))))))*100),'%'));
       
%========CV Momementum Conservation=======         
%Find results files
resultfile_dir = './../../results/120201_CV_conservation';

%Check CV momentum has been recorded
cd(resultfile_dir);
fid = fopen('./vflux','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    return % stop rest if not...
end

%Read Header file
read_header
%Establish and Plot domain set-up
Domain_setup       
%setup stress direction and face
ixyz = 1;
jxyz = 2;
Nvflux_records = Nsteps / (Nvflux_ave);
%Check CV are satisfied
for m =1:1:Nvflux_records-3
    m
    %Load momentum flux values for current timestep
    [velocity_snapshot(:,:,:,:), ...
     velocity_flux(:,:,:,:,:),   ... 
     pressure_surface(:,:,:,:,:)    ] = read_vflux(m,resultfile_dir,gnbins,nd);

     [velocity_snapshot_tplus1(:,:,:,:)] = read_vflux(m+1,resultfile_dir,gnbins,nd);
     [velocity_snapshot_tplus2(:,:,:,:)] = read_vflux(m+2,resultfile_dir,gnbins,nd);
     
     % %Calculate total CV flux and change in mass
     totalflux =((velocity_flux(:,:,:,:,1)+velocity_flux(:,:,:,:,4)))/(binsize(1)) ...
               +((velocity_flux(:,:,:,:,2)+velocity_flux(:,:,:,:,5)))/(binsize(2)) ...
               +((velocity_flux(:,:,:,:,3)+velocity_flux(:,:,:,:,6)))/(binsize(3));
     totalpressure =((pressure_surface(:,:,:,:,1)-pressure_surface(:,:,:,:,4)))/(binsize(1)) ...
                   +((pressure_surface(:,:,:,:,2)-pressure_surface(:,:,:,:,5)))/(binsize(2)) ...
                   +((pressure_surface(:,:,:,:,3)-pressure_surface(:,:,:,:,6)))/(binsize(3));
                  
     %totalpressure = totalpressure*delta_t
     dvelocitydt(:,:,:,:) =  (velocity_snapshot_tplus2(:,:,:,:) - velocity_snapshot_tplus1(:,:,:,:))/(delta_t*Nvflux_ave);
     
     
     domain_ave_pressure_tensor(:,1,m) = mean(mean(mean((pressure_surface(:,:,:,:,1)+pressure_surface(:,:,:,:,4)),1),2),3)/2;
     domain_ave_pressure_tensor(:,2,m) = mean(mean(mean((pressure_surface(:,:,:,:,2)+pressure_surface(:,:,:,:,5)),1),2),3)/2; 
     domain_ave_pressure_tensor(:,3,m) = mean(mean(mean((pressure_surface(:,:,:,:,3)+pressure_surface(:,:,:,:,6)),1),2),3)/2;

     a(m,1,:) =totalpressure(ibin,jbin,kbin,:);
     a(m,2,:) =dvelocitydt(ibin,jbin,kbin,:);
     a(m,3,:) =totalflux(ibin,jbin,kbin,:);
     
     Error(m) =  max(max(max(((squeeze(totalpressure(:,:,:,ixyz)) ...
                             -squeeze(totalflux(:,:,:,ixyz)) )    ...
                             -squeeze(dvelocitydt(:,:,:,ixyz)))*10000))) ...
               + min(min(min(((squeeze(totalpressure(:,:,:,ixyz)) ...
                             -squeeze(totalflux(:,:,:,ixyz)))     ...
                             -squeeze(dvelocitydt(:,:,:,ixyz)))*10000)));
     
     %sliceomatic(( squeeze(totalpressure(:,:,:,ixyz)) ...
     %             -squeeze(totalflux(:,:,:,ixyz))     ...
     %             -squeeze(dvelocitydt(:,:,:,ixyz)))*10000)
    % pause()
    
    
end
%Adjust as multiplied by delta_t before write out
a = a*delta_t;

%========CV Energy Conservation=======  
for m =1:1:Nvflux_records-2
    m
    %Load energy flux values for current timestep
    [energy_snapshot,energy_flux,energy_surfaces]= read_eflux(m,resultfile_dir,gnbins);
       
    [energy_snapshot_tplus1(:,:,:,:)] = read_eflux(m,resultfile_dir,gnbins);
    [energy_snapshot_tplus2(:,:,:,:)] = read_eflux(m+1,resultfile_dir,gnbins);
     
     % %Calculate total CV flux and change in mass
	totalflux =((energy_flux(:,:,:,1)+energy_flux(:,:,:,4)))/(binsize(1)) ...
               +((energy_flux(:,:,:,2)+energy_flux(:,:,:,5)))/(binsize(2)) ...
               +((energy_flux(:,:,:,3)+energy_flux(:,:,:,6)))/(binsize(3));
    totalpower    =((energy_surfaces(:,:,:,1)-energy_surfaces(:,:,:,4)))/(binsize(1)) ...
                   +((energy_surfaces(:,:,:,2)-energy_surfaces(:,:,:,5)))/(binsize(2)) ...
                   +((energy_surfaces(:,:,:,3)-energy_surfaces(:,:,:,6)))/(binsize(3));
    
     denergydt(:,:,:,:) =  (energy_snapshot_tplus2(:,:,:,:) - energy_snapshot_tplus1(:,:,:,:))/(delta_t*Nvflux_ave);
     
     Error(m) =  max(max(max(((squeeze(totalpower(:,:,:,ixyz)) ...
                             -squeeze(totalflux(:,:,:,ixyz)) )    ...
                             -squeeze(denergydt(:,:,:,ixyz)))*10000))) ...
               + min(min(min(((squeeze(totalpower(:,:,:,ixyz)) ...
                             -squeeze(totalflux(:,:,:,ixyz)))     ...
                             -squeeze(denergydt(:,:,:,ixyz)))*10000)));
                         
     b(m,1) =totalpower(ibin,jbin,kbin);
     b(m,2) =totalflux(ibin,jbin,kbin);
     b(m,3) =denergydt(ibin,jbin,kbin);             
end
%Adjust as multiplied by delta_t before write out
b = b*delta_t;

%Plot time evolution graphs
%======Mass=======
scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
xaxis = [0:0.005:0.005*(size(dmassdt(ibin,jbin,kbin,:),4)-1)]; 
plot(xaxis,squeeze(totalmflux(ibin,jbin,kbin,:)),'k--', 'markersize', 5,'lineWidth', 3);
hold on
plot(xaxis,-squeeze(dmassdt(ibin,jbin,kbin,:)),'--','Color',[.6 .6 .6],'lineWidth', 5);
Residual = totalmflux(ibin,jbin,kbin,:) - dmassdt(ibin,jbin,kbin,:);
plot(xaxis,squeeze(Residual),'k-','lineWidth', 2);
set(gca,'FontSize',16)
%xlabel('Time'); 
ylabel('|Mass|')
set(gca,'xtick',[])
axis([1.3 2.1 -1.1 1.1])
%legend ('MdS','\Delta M','Residual','location','BestOutside')
%legend('boxoff')
hold off

%======Momentum=======
subplot(1, 1, 1)
%main plot
scrsz = get(0,'ScreenSize');
fig2 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
xaxis = [0:0.005:0.005*(size(dmassdt(ibin,jbin,kbin,:),4)-1)]; 
plot(xaxis(3:size(a,1)+2),a(:,2,1),'-','Color',[.6 .6 .6],'lineWidth', 4);
hold on
plot(xaxis(3:size(a,1)+2),a(:,3,1),'k--', 'markersize', 5,'lineWidth', 4);
plot(xaxis(3:4:size(a,1)+2),a(1:4:end,1,1), 'kx','lineWidth', 2,'MarkerSize',10);
Residual = squeeze(a(:,1,:))-squeeze(a(:,2,:))-squeeze(a(:,3,:));
plot(xaxis(1:size(a,1)),squeeze(Residual),'k-','lineWidth', 4);
set(gca,'FontSize',16)
xlabel('Time');
ylabel('Momentum')
%set(gca,'xtick',[])
axis([1.3 2.01 -0.31 0.31])
%legend ('Accumulation','Advection','Forcing','Residual','location','SouthWest')
%legend('boxoff')
hold off
%Plot insert
 hax = axes('Position', [.2, .15, .23, .23]);
 plot(xaxis(3:size(a,1)+2),a(:,2,1),'-','Color',[.6 .6 .6],'lineWidth', 4);
 hold on
 plot(xaxis(3:size(a,1)+2),a(:,3,1),'k--', 'markersize', 5,'lineWidth', 4);
 plot(xaxis(3:4:size(a,1)+2),a(1:4:end,1,1), 'kx','lineWidth', 2,'MarkerSize',2);
 Residual = squeeze(a(:,1,:))-squeeze(a(:,2,:))-squeeze(a(:,3,:));
 plot(xaxis(1:size(a,1)),squeeze(Residual),'k-','lineWidth', 4);
 set(hax,'XTick',[]); 
 axis([1.3 2.01 -1.3 1.3])

%======Energy=======
scrsz = get(0,'ScreenSize');
fig3 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
plot(xaxis(2:size(b,1)+1),b(:,3)','-','Color',[.6 .6 .6],'lineWidth', 4)
hold on
%plot(xaxis(2:size(b,1)+1),0.5*(squeeze(b(:,1))+circshift(squeeze(b(:,1)),-1)), 'b-','lineWidth', 7)
plot(xaxis(2:size(b,1)+1),b(:,2)','k--', 'markersize', 5,'lineWidth', 4)
plot(xaxis(2:4:size(b,1)+1)+0.5*delta_t,b(1:4:end,1), 'kx','lineWidth', 2,'MarkerSize',10)
Residual = squeeze(b(:,3)) - 0.5*(squeeze(b(:,1))+circshift(squeeze(b(:,1)),1)) + squeeze(b(:,2));
plot(xaxis(2:size(b,1)+1),squeeze(Residual),'-k','lineWidth', 4);
xaxis = [0:0.005:0.005*(size(dmassdt(ibin,jbin,kbin,:),4)-1)]; 
set(gca,'FontSize',16)
xlabel('Time'); ylabel('Energy')
%legend ('Accumulation','Advection','Forcing','Residual','location','SouthWest')
%legend('boxoff')
hold off
axis([1.3 2.01 -0.31 0.31])

hax = axes('Position', [.2, .15, .23, .23]);
plot(xaxis(2:size(b,1)+1),b(:,3)','-','Color',[.6 .6 .6],'lineWidth', 4)
hold on
%plot(xaxis(2:size(b,1)+1),0.5*(squeeze(b(:,1))+circshift(squeeze(b(:,1)),-1)), 'b-','lineWidth', 7)
plot(xaxis(2:size(b,1)+1),b(:,2)','k--', 'markersize', 5,'lineWidth', 4)
plot(xaxis(2:4:size(b,1)+1)+0.5*delta_t,b(1:4:end,1), 'kx','lineWidth', 2,'MarkerSize',2)
Residual = squeeze(b(:,3)) - 0.5*(squeeze(b(:,1))+circshift(squeeze(b(:,1)),1)) + squeeze(b(:,2));
plot(xaxis(2:size(b,1)+1),squeeze(Residual),'-k','lineWidth', 4);
axis([1.3 2.01 -2.5 2.5])
 set(hax,'XTick',[]); 
%Adjust overall suplot spacing
% p = get(h(2), 'pos');
% p(2) = p(2) + 0.054;
% set(h(2), 'pos', p);
% 
% p = get(h(3), 'pos');
% p(2) = p(2) + 0.108;
% set(h(3), 'pos', p);

minx = 0;
maxx = 200;
%subplot(3,1,1), axis([minx maxx -1 1])

%subplot(3,1,2), axis([minx maxx -20 20])
%subplot(3,1,3), axis([minx maxx -50 50])

%semilogy(abs(b(:,3)),'b')
%hold on
%semilogy(1.5:size(b(:,1),1)+0.5,abs(b(:,1)+b(:,2)),'--r')
%legend('totalpower','totalflux','denergydt','location','NorthWest')
%legend('boxoff')