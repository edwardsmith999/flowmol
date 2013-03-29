%Plot all output files from simulation
close all
clear all

scrsz = get(0,'ScreenSize');
set(0,'DefaultFigureRenderer','OpenGL')
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/4 scrsz(4)/2]);

%pwdir = '/home/es205/codes/coupled/MD_dCSE/src_code/post_proc/MATLAB';
%resultfile_dir = './../../results/';
%pwdir='/home/es205/results/md_results/fortran/3D_code/parallel/results/converge_diverge';
%resultfile_dir = '/home/es205/results/md_results/fortran/3D_code/parallel/results/converge_diverge/';
pwdir='/home/es205/codes/coupled/MD_dCSE/src_code/';
resultfile_dir = '/home/es205/codes/coupled/MD_dCSE/src_code/results/';
%Read Header
read_header

%ibin = floor(gnbins(1)/2.0); jbin = floor(gnbins(2)/2.0); kbin=floor(gnbins(3)/2.0);
ibin = floor(gnbins(2)/2.0); jbin = 2 ; kbin=floor(gnbins(3)/2.0);

%Check if cv conservation is employed - exit if not
if (exist('cv_conserve'))
    cv_conserve = cv_conserve;
else
    cv_conserve = 0;
end
if (cv_conserve == 1)
    Nmflux_records = (Nsteps-initialstep) / (Nmflux_ave);
else
    error('CV used for averages only so not conserved')
end

%========CV Mass Conservation=======
n = 1;
skip = 10;
for m =1:skip:Nmflux_records-2
    m
    %Spacial evolution of domain at time half way from start to finish
    [mass_flux,mass_snapshot] = read_mflux('./mflux','./msnap',resultfile_dir,m);
    
    %Calculate total CV flux and change in mass
    totalmflux = squeeze(sum(mass_flux,4));
    
    %Spacial evolution of domain between time m and m+1
    [mass_flux_tp1,mass_snapshot_tp1] = read_mflux('./mflux','./msnap',resultfile_dir,m+1);
    dmassdt(:,:,:) = mass_snapshot_tp1 - mass_snapshot;
    
    
    
    %Verify that CV mass is exactly conservative
    conserved = (squeeze(totalmflux(:,:,:)) - squeeze(dmassdt(:,:,:)));
    if (max(conserved(:)) > 0.0 || min(conserved(:)) < 0.0)
        %Display message
        display(strcat('Error in CV mass conservation =', ...
            num2str(max(conserved(:))+abs(min(conserved(:))))*100),'% - beginning debug');
        
        %Plot slice of regions where not conserved
        h=slice(conserved,[],[],[10]);
        view([10,20,5]); axis 'equal';
        set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
        colorbar; drawnow; pause()
        
        %Log temporal evolution over 100 timesteps
        skip = 1;
        %Save conservation time plot for a cell
        dt(n) = squeeze(dmassdt(ibin,jbin,kbin));
        dm(n) = squeeze(totalmflux(ibin,jbin,kbin));
        n = n + 1;
        if (n == 100)
            %Plot 100 steps of a single cell
            plot(dt)
            hold all
            plot(dm,'r--')
            break
        end
    else
        display(strcat('CV mass conservation OK'));
    end
    
end

%========CV Momementum Conservation=======
clear mass_flux mass_snapshot mass_flux_tp1 mass_snapshot_tp1 conserved

%Read Header file
read_header
%Establish and Plot domain set-up
Domain_setup
%setup stress direction and face
ixyz = 1;
jxyz = 2;
Nvflux_records = Nsteps / (Nvflux_ave);
skip =1; n = 1; tol = 10*eps;
external_force_flag = 1;
%Check CV are satisfied
for m =1:skip:Nvflux_records-3
    m
    %Load momentum flux values for current timestep -- include external
    %forces if relevant
    if (external_force_flag)
        [velocity_snapshot(:,:,:,:), ...
            velocity_flux(:,:,:,:,:),   ...
            pressure_surface(:,:,:,:,:) ...
            F_ext(:,:,:,:)                  ] = read_vflux(m,resultfile_dir,gnbins,nd);
    else
        [velocity_snapshot(:,:,:,:), ...
            velocity_flux(:,:,:,:,:),   ...
            pressure_surface(:,:,:,:,:) ] = read_vflux(m,resultfile_dir,gnbins,nd);
        F_ext = zeros(size(velocity_snapshot));
    end
    
    [velocity_snapshot_tplus1(:,:,:,:)] = read_vflux(m+1,resultfile_dir,gnbins,nd);
   
    % Calculate total CV flux and change in mass
    totalflux =((velocity_flux(:,:,:,:,1)+velocity_flux(:,:,:,:,4)))/(binsize(1)) ...
              +((velocity_flux(:,:,:,:,2)+velocity_flux(:,:,:,:,5)))/(binsize(2)) ...
              +((velocity_flux(:,:,:,:,3)+velocity_flux(:,:,:,:,6)))/(binsize(3));
    totalpressure =((pressure_surface(:,:,:,:,1)-pressure_surface(:,:,:,:,4)))/(binsize(1)) ...
                  +((pressure_surface(:,:,:,:,2)-pressure_surface(:,:,:,:,5)))/(binsize(2)) ...
                  +((pressure_surface(:,:,:,:,3)-pressure_surface(:,:,:,:,6)))/(binsize(3));
    
    %totalpressure = totalpressure*delta_t
    dvelocitydt(:,:,:,:) =  (velocity_snapshot_tplus1(:,:,:,:)  ...
        - velocity_snapshot(:,:,:,:))/(delta_t*Nvflux_ave);
    
    
    %domain_ave_pressure_tensor(:,1,m) = mean(mean(mean((pressure_surface(:,:,:,:,1)+pressure_surface(:,:,:,:,4)),1),2),3)/2;
    %domain_ave_pressure_tensor(:,2,m) = mean(mean(mean((pressure_surface(:,:,:,:,2)+pressure_surface(:,:,:,:,5)),1),2),3)/2;
    %domain_ave_pressure_tensor(:,3,m) = mean(mean(mean((pressure_surface(:,:,:,:,3)+pressure_surface(:,:,:,:,6)),1),2),3)/2;
    
    %Verify that CV momentum is exactly conservative
    conserved = ( squeeze(sum(totalpressure(:,:,:,:),4)) ...
                 -squeeze(sum(totalflux(:,:,:,:),4))     ...
                 -squeeze(sum(dvelocitydt(:,:,:,:),4))   ...
                 +squeeze(sum(F_ext(:,:,:,:),4)/prod(binsize)));
    
    if (max(conserved(:)) > tol || min(conserved(:)) < -tol)
        Error(n) =  max(conserved(:))+abs(min(conserved(:)));
        
        display(strcat('Error in CV mass conservation =', ...
            num2str(Error(n)*100),'% - beginning debug'));
        
        h=slice(conserved(:,:,:),[],[],[5]);
        view([2]); axis 'tight';
        set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
        title(['d\rhou/dt - \nabla \cdot \Pi',num2str(ibin),',',num2str(jbin),',',num2str(kbin)])
        colorbar;  drawnow; pause(0.1)
        

        %Log temporal evolution over 100 timesteps
        skip = 1;
        %Save conservation time plot for a cell
        a(n,1,:) = totalpressure(ibin,jbin,kbin,:);
        a(n,2,:) = dvelocitydt(ibin,jbin,kbin,:);
        a(n,3,:) = totalflux(ibin,jbin,kbin,:);
        a(n,4,:) = F_ext(ibin,jbin,kbin,:);
        n = n + 1;
        if (n == 100 || m == Nvflux_records-3)

            %Adjust as multiplied by delta_t before write out
            a = a*delta_t;

            %Plot 100 steps of a single cell
            plot(sum(a(1:end,2,:),3))   %du/dt =
            hold all
            plot(-sum(a(1:end,3,:),3) + sum(a(1:end,1,:),3)+sum(a(1:end,4,:),3)/prod(binsize),'r--') % - d\rho uu/dr - dP/dr + dsigma/dr
            title(['d\rhou/dt and \nabla \cdot \Pi vs time in CV ',num2str(ibin),',',num2str(jbin),',',num2str(kbin)])
            %Plot error
            figure
            plot(Error)
            hold all
            plot(0:size(Error,2),ones(1,size(Error,2)+1)*eps('double'))
            title(['Maximum Error in whole domain vs time'])

            break
        end
        
    else
        display(strcat('CV momentum conservation OK'));
    end
    %
    %             figure(fig2)
    %             h=slice(-squeeze(sum(F_ext(:,:,:,:),4)),[],[],[5]);
    %             view([2]); axis 'tight';
    %             set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
    %             colorbar; caxis([-1 1]); drawnow; pause(0.1)
    
    %sliceomatic(( squeeze(sum(totalpressure(:,:,:,:),4)) ...
    %             -squeeze(sum(totalflux(:,:,:,:),4))     ...
    %             -squeeze(sum(dvelocitydt(:,:,:,:),4)))*10000)
    
end


%Plot all components
% plot(sum(a(1:end,2,:),3),'g')   %du/dt
% hold on
% plot(sum(a(2:end-1,4,:),3)/prod(binsize),'k') %F_ext
% plot(-sum(a(1:end,3,:),3),'r') % - d\rho uu/dr - dP/dr
% plot(-sum(a(1:end,1,:),3),'b') %  dsigma/dr



%pause()
%
% %========CV Energy Conservation=======
% for m =1:1:Nvflux_records-2
%     m
%     %Load energy flux values for current timestep
%     [energy_snapshot,energy_flux,energy_surfaces]= read_eflux(m,resultfile_dir,gnbins);
%
%     [energy_snapshot_tplus1(:,:,:,:)] = read_eflux(m,resultfile_dir,gnbins);
%     [energy_snapshot_tplus2(:,:,:,:)] = read_eflux(m+1,resultfile_dir,gnbins);
%
%      % %Calculate total CV flux and change in mass
% 	totalflux =((energy_flux(:,:,:,1)+energy_flux(:,:,:,4)))/(binsize(1)) ...
%                +((energy_flux(:,:,:,2)+energy_flux(:,:,:,5)))/(binsize(2)) ...
%                +((energy_flux(:,:,:,3)+energy_flux(:,:,:,6)))/(binsize(3));
%     totalpower    =((energy_surfaces(:,:,:,1)-energy_surfaces(:,:,:,4)))/(binsize(1)) ...
%                    +((energy_surfaces(:,:,:,2)-energy_surfaces(:,:,:,5)))/(binsize(2)) ...
%                    +((energy_surfaces(:,:,:,3)-energy_surfaces(:,:,:,6)))/(binsize(3));
%
%      denergydt(:,:,:,:) =  (energy_snapshot_tplus2(:,:,:,:) - energy_snapshot_tplus1(:,:,:,:))/(delta_t*Nvflux_ave);
%
%      Error(m) =  max(max(max(((squeeze(totalpower(:,:,:,ixyz)) ...
%                              -squeeze(totalflux(:,:,:,ixyz)) )    ...
%                              -squeeze(denergydt(:,:,:,ixyz)))*10000))) ...
%                + min(min(min(((squeeze(totalpower(:,:,:,ixyz)) ...
%                              -squeeze(totalflux(:,:,:,ixyz)))     ...
%                              -squeeze(denergydt(:,:,:,ixyz)))*10000)));
%
%      b(m,1) =totalpower(ibin,jbin,kbin);
%      b(m,2) =totalflux(ibin,jbin,kbin);
%      b(m,3) =denergydt(ibin,jbin,kbin);
% end
% %Adjust as multiplied by delta_t before write out
% b = b*delta_t;

%Plot time evolution graphs
% h = tight_subplot(3,1);
% %======Mass=======
% xaxis = [0:0.005:0.005*(size(dmassdt(ibin,jbin,kbin,:),4)-1)];
% subplot(h(1)),
% plot(xaxis,squeeze(totalmflux(ibin,jbin,kbin,:)),'k--', 'markersize', 5,'lineWidth', 3);
% hold on
% plot(xaxis,-squeeze(dmassdt(ibin,jbin,kbin,:)),'--','Color',[.6 .6 .6],'lineWidth', 5);
% Residual = totalmflux(ibin,jbin,kbin,:) - dmassdt(ibin,jbin,kbin,:);
% plot(xaxis,squeeze(Residual),'k-','lineWidth', 2);
% set(gca,'FontSize',16)
% %xlabel('Time');
% ylabel('|Mass|')
% set(gca,'xtick',[])
% axis([1.2 1.6 -1.1 1.1])
% %legend ('MdS','\Delta M','Residual','location','BestOutside')
% %legend('boxoff')
% hold off
%
% %======Momentum=======
% subplot(h(2)),
% xaxis = [0:0.005:0.005*(size(dmassdt(ibin,jbin,kbin,:),4)-1)];
% plot(xaxis(3:size(a,1)+2),a(:,1,1), 'k-','lineWidth', 7);
% hold on
% plot(xaxis(3:size(a,1)+2),a(:,3,1),'k--', 'markersize', 5,'lineWidth', 3);
% plot(xaxis(3:size(a,1)+2),a(:,2,1),'--','Color',[.6 .6 .6],'lineWidth', 5);
% Residual = squeeze(a(:,1,:))-squeeze(a(:,2,:))-squeeze(a(:,3,:));
% plot(xaxis(1:size(a,1)),squeeze(Residual),'k-','lineWidth', 2);
% set(gca,'FontSize',16)
% %xlabel('Time');
% ylabel('|Momentum|')
% set(gca,'xtick',[])
% axis([1.2 2.0 -0.3 0.3])
% %legend ('\Delta F','MdS','\Delta M','Residual','location','BestOutside')
% %legend('boxoff')
% hold off
%
% %======Energy=======
% subplot(h(3)),
% plot(xaxis(2:size(b,1)+1)+0.5*delta_t,b(:,1), 'k-','lineWidth', 7)
% hold on
% %plot(xaxis(2:size(b,1)+1),0.5*(squeeze(b(:,1))+circshift(squeeze(b(:,1)),-1)), 'b-','lineWidth', 7)
% plot(xaxis(2:size(b,1)+1),b(:,2)','k--', 'markersize', 5,'lineWidth', 3)
% plot(xaxis(2:size(b,1)+1),b(:,3)','--','Color',[.6 .6 .6],'lineWidth', 5)
% Residual = squeeze(b(:,3)) - 0.5*(squeeze(b(:,1))+circshift(squeeze(b(:,1)),1)) + squeeze(b(:,2));
% plot(xaxis(2:size(b,1)+1),squeeze(Residual),'-k','lineWidth', 2);
% xaxis = [0:0.005:0.005*(size(dmassdt(ibin,jbin,kbin,:),4)-1)];
% set(gca,'FontSize',16)
% xlabel('Time'); ylabel('|Energy|')
% %legend('\Delta Fv','EdS','\Delta E','Residual','location','BestOutside')
% %legend('boxoff')
% hold off
% axis([1.2 2.0 -0.2 0.4])

%Adjust overall suplot spacing
% p = get(h(2), 'pos');
% p(2) = p(2) + 0.054;
% set(h(2), 'pos', p);
%
% p = get(h(3), 'pos');
% p(2) = p(2) + 0.108;
% set(h(3), 'pos', p);

%minx = 0;
%maxx = 200;
%subplot(3,1,1), axis([minx maxx -1 1])

%subplot(3,1,2), axis([minx maxx -20 20])
%subplot(3,1,3), axis([minx maxx -50 50])

%semilogy(abs(b(:,3)),'b')
%hold on
%semilogy(1.5:size(b(:,1),1)+0.5,abs(b(:,1)+b(:,2)),'--r')
%legend('totalpower','totalflux','denergydt','location','NorthWest')
%legend('boxoff')