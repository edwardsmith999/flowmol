%Plot all output files from simulation
close all
clear all

%Find results files
resultfile_dir = './../results';

%Read Header file
read_header
read_continuum_header

%Establish and Plot domain set-up
Domain_setup

%Read continuum output
if (continuum_vflag == 3)
    read_continuum_vbins
else
    read_continuum_vslice
end

%Convert to stress
continuum_stress=continuum_vslice_to_stress(continuum_velslice,meu,ly,ny);
Nvflux_records = Nsteps / (tplot * Nvflux_ave);
t_ave = 100;

%setup stress direction and face
ixyz = 1;
jxyz = 2;

%Plot Stress
overlap = 3;
coupleddomain = ly + globaldomain(2) - 2*overlap*binsize(2);
coupledliquiddomain = coupleddomain - wallbot(2);
couplednbins = (globalnbins(ixyz)-1)/2+ny-overlap;
xaxis = 0.05:0.1:1.05;
xaxis2 = 1/nbins(2):1/nbins(2):1;
spectral_res = 6; %Ratio of points for spectral analytical solution to domain bins
analy_points = spectral_res*(couplednbins); % Number of spectral points
xaxis_analy =  0:1/(analy_points):1;
u_wall = 1;
%xloc_MDCFDhalocell = 0.5*(xaxis(couplednbins-ny) + xaxis(couplednbins-ny+1));
timeratio = delta_t/continuum_delta_t;

scrsz = get(0,'ScreenSize');
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3) scrsz(4)]);

% set(gca,'nextplot','replacechildren','visible','off')
% f = getframe;
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,20) = 0;
%    hold on
   
% %Check CV are satisfied
% for m =1:1:Nvflux_records-1
%     [velocity_snapshot(:,:,:,:), ...
%      velocity_flux(:,:,:,:,:),   ... 
%      pressure_surface(:,:,:,:,:)    ] = read_vflux(m,resultfile_dir,globalnbins,nd);
% 
%      [velocity_snapshot_tplus1(:,:,:,:)] = read_vflux(m+1,resultfile_dir,globalnbins,nd);
%  
%      % %Calculate total CV flux and change in mass
%      totalflux =((velocity_flux(:,:,:,:,1)+velocity_flux(:,:,:,:,4)))/(binsize(1)) ...
%                +((velocity_flux(:,:,:,:,2)+velocity_flux(:,:,:,:,5)))/(binsize(2)) ...
%                +((velocity_flux(:,:,:,:,3)+velocity_flux(:,:,:,:,6)))/(binsize(3));
%      totalpressure =((pressure_surface(:,:,:,:,1)-pressure_surface(:,:,:,:,4)))/(binsize(1)) ...
%                    +((pressure_surface(:,:,:,:,2)-pressure_surface(:,:,:,:,5)))/(binsize(2)) ...
%                    +((pressure_surface(:,:,:,:,3)-pressure_surface(:,:,:,:,6)))/(binsize(3));
%                   
%      %totalpressure = totalpressure*delta_t
%      dvelocitydt(:,:,:,:) =  (velocity_snapshot_tplus1(:,:,:,:) - velocity_snapshot(:,:,:,:))/(delta_t*Nvflux_ave);
%      
%      sliceomatic(( squeeze(totalpressure(:,:,:,ixyz)) ...
%                   -squeeze(totalflux(:,:,:,ixyz))     ...
%                   -squeeze(dvelocitydt(:,:,:,ixyz)))*10000)
%      pause()
% end

%Get steady state stress for normalisation
[temp, velocity_flux, pressure_surface ] = read_vflux(Nvflux_records,resultfile_dir,globalnbins,nd);
final_pressure =-0.5*(squeeze(mean(mean(mean(pressure_surface(:,:,:,ixyz,jxyz),1),2),3)) ...
                     +squeeze(mean(mean(mean(velocity_flux(:,:,:,ixyz,jxyz),1),2),3)) ...
                     +squeeze(mean(mean(mean(pressure_surface(:,:,:,ixyz,jxyz+3),1),2),3)) ...
                 	 +squeeze(mean(mean(mean(velocity_flux(:,:,:,ixyz,jxyz+3),1),2),3)));
                 
final_pressure = 1; %0.041;

%vidObj = VideoWriter('couette_stress.avi');
%vidObj.FrameRate = 10;
%open(vidObj);

for m = 1:10:Nvflux_records-1
    
    t = (m-0.5)*delta_t*Nvflux_ave*tplot;
    analy = couette_analytical_stress_fn(t,Re,1,couplednbins*2*binsize(2),spectral_res*couplednbins);
    plot(xaxis_analy,-analy*meu/final_pressure,'r','LineWidth',5);
    set(gca,'FontSize',20)
    hold on
    axis([-0.1 1.1 -0.1 0.6])
        
    %Coarse grain slices as 2 bins in MD = 1 bin in CFD and in time
    ave_pressure_surface = zeros(nbins(1),floor(nbins(2)/2),nbins(3),nd);
    for j=1:t_ave
        [velocity_snapshot(:,:,:,:), ...
           velocity_flux(:,:,:,:,:), ... 
         pressure_surface(:,:,:,:,:)     ] = read_vflux(m+(j-1),resultfile_dir,globalnbins,nd);
    
        n = 1;
        for i=2:2:globalnbins(jxyz)
            ave_pressure_surface(:,n,:,:) = ave_pressure_surface(:,n,:,:) + pressure_surface(:,i  ,:,:,jxyz  ) ...
                                                                          + velocity_flux   (:,i  ,:,:,jxyz  ) ...
                                                                          + pressure_surface(:,i  ,:,:,jxyz+3) ...
                                                                          + velocity_flux   (:,i  ,:,:,jxyz+3) ...
                                                                          + pressure_surface(:,i+1,:,:,jxyz+3) ...
                                                                          + velocity_flux   (:,i+1,:,:,jxyz+3);
            n = n+1;
        end
    end
    ave_pressure_surface = ave_pressure_surface/(3*t_ave);
    
    %plot(0.025:0.05:0.575,squeeze(mean(mean(pressure_surface(:,2:end,:,ixyz,jxyz),1),3)),'.')   
    
    plot(xaxis(1:(globalnbins(jxyz)-1)/2), -squeeze(mean(mean(ave_pressure_surface(:,:,:,ixyz),1),3)), ...
         'x','LineWidth',5,'Color',[.5 .5 .5],'MarkerSize',20);
          
    %Plot continuum stress profile
    plot(xaxis(((globalnbins(jxyz)-1)/2)+1-overlap:end-1),continuum_stress(:,m),'s','Color',[.5 .5 .5],'MarkerSize',20,'LineWidth',5);
    
    %Plot lines to show continuous stress of CV
    %plot(xaxis,squeeze(mean(mean(pressure_surface(:,:,:,ixyz,jxyz,m),1),3)),'LineWidth',5,'Color',[.5 .5 .5]);
    %plot(xaxis2,squeeze(mean(mean(pressure_surface(:,:,:,ixyz,jxyz+3,m),1),3)),'LineWidth',5,'Color',[.5 .5 .5]);

    legend ('Analytical','MD','CFD','location','NorthWest'); legend('boxoff')
    xlabel('y/H'); ylabel('P_{xy}')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    hold off
    pause(0.5)
    %f = getframe(gcf);
    %im(:,:,3,m) = rgb2ind(f.cdata,map,'nodither');

    %Store pictures and videos
    %savefig(strcat('Velocity_',num2str(i)),'png')
    %currFrame = getframe(gcf);
    %writeVideo(vidObj,currFrame);
   
    if (m > 2000)
        break
    end
end

% Close the file.
%close(vidObj)
%movefile ./couette_stress.avi ./../results
% 
% timestep = 5;
% [velocity_snapshot(:,:,:,:,timestep)] = read_vflux(timestep,resultfile_dir,globalnbins,nd);
% 
%   
% %Verify adjacent surface values are equal
% %Random Bin
% for i=1:nd
%     randbin(i) = round(rand*(globalnbins(i)-1))+1;
% end
% scatter(0:globalnbins(1)-1,velocity_flux(:,randbin(2),randbin(3),1,1)+...
%                         pressure_surface(:,randbin(2),randbin(3),1,1),'rx')
% hold on
% scatter(1:globalnbins(1),-velocity_flux(:,randbin(2),randbin(3),1,4)+...
%                        pressure_surface(:,randbin(2),randbin(3),1,4),'bo')
%                    
% % %Calculate total CV flux and change in mass
%  totalflux     =((velocity_flux(:,:,:,:,1,:)+velocity_flux(:,:,:,:,4,:)))/(binsize(1)) ...
%                +((velocity_flux(:,:,:,:,2,:)+velocity_flux(:,:,:,:,5,:)))/(binsize(2)) ...
%                +((velocity_flux(:,:,:,:,3,:)+velocity_flux(:,:,:,:,6,:)))/(binsize(3));
%  totalpressure =((pressure_surface(:,:,:,:,1,:)-pressure_surface(:,:,:,:,4,:)))/(binsize(1)) ...
%                +((pressure_surface(:,:,:,:,2,:)-pressure_surface(:,:,:,:,5,:)))/(binsize(2)) ...
%                +((pressure_surface(:,:,:,:,3,:)-pressure_surface(:,:,:,:,6,:)))/(binsize(3));
% %   
% %totalpressure = totalpressure*delta_t           
%   for i =1:Nvflux_records-1
%       dvelocitydt(:,:,:,:,i) =  velocity_snapshot(:,:,:,:,i+1) - velocity_snapshot(:,:,:,:,i);
%   end
%   
% %Verify that CV is exactly conservative
% figure
% %Random Bin
%   for i=1:nd
%      randbin(i) = round(rand*(globalnbins(i)-1))+1;
%  end
%  plot(squeeze(dvelocitydt(randbin(1),randbin(2),randbin(3),1,:)));
%  hold on
%  plot(squeeze(dvelocitydt(randbin(1),randbin(2),randbin(3),2,:)));
% %plot(squeeze(dvelocitydt(randbin(1),randbin(2),randbin(3),3,:)));
% 
% plot(squeeze(totalpressure(randbin(1),randbin(2),randbin(3),1,:)) ...
%     -squeeze(totalflux(randbin(1),randbin(2),randbin(3),1,:)),'-r');
% 
% plot(-squeeze(totalflux(randbin(1),randbin(2),randbin(3),1,:)),'--r');
% %plot(-squeeze(totalflux(randbin(1),randbin(2),randbin(3),2,:)),'--r');
% %plot(-squeeze(totalflux(randbin(1),randbin(2),randbin(3),3,:)),'--r');
% 
% plot(squeeze(totalpressure(randbin(1),randbin(2),randbin(3),1,:)),'-.c');
% %plot(squeeze(totalpressure(randbin(1),randbin(2),randbin(3),2,:)),'-.c');
% %plot(squeeze(totalpressure(randbin(1),randbin(2),randbin(3),3,:)),'-.c');
% 
% display(strcat('Maximum error in CV conservation =', ... 
%    num2str(max(max(max(max(max(squeeze(totalpressure)-squeeze(totalflux)-squeeze(dvelocitydt))))))*...
%            min(min(min(min(min(squeeze(totalpressure)-squeeze(totalflux)-squeeze(dvelocitydt))))))*100),'%'));

% % 
% % scrsz = get(0,'ScreenSize');
% % %fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
% % fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3) scrsz(4)]);
% % 
% % %Calculate spanwise direction
% % if (velocity_outflag < 4)
% %     ixyz = velocity_outflag;
% % else
% %     %Find location of smallest dimension and take slice through this
% %     ixyz = find(nbinsliquid==min(nbinsliquid));
% %     jxyz = mod(ixyz+1,3)+1;
% %     kxyz = mod(ixyz,3)+1;
% %     mass_slice = squeeze(sum(sum(mass_bins(:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
% %     vel_slice =  squeeze(sum(sum(vel_bins(:,:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
% % end
% % 
% % 
% % %Coarse grain slices as 2 bins in MD = 1 bin in CFD
% % n = 1;
% % for i=2:2:globalnbins(ixyz)
% %     ave_vel_slice(n,:,:) = vel_slice(i,:,:) + vel_slice(i+1,:,:);
% %     ave_mass_slice(n,:)  = mass_slice(i,:)  + mass_slice(i+1,:);
% %     n = n+1;
% % end
% % 
% % %Average velocity per molecule
% % for n =1:3
% %     vel_slice(:,n,:) = squeeze(vel_slice(:,n,:))./mass_slice;
% %     ave_vel_slice(:,n,:) = squeeze(ave_vel_slice(:,n,:))./ave_mass_slice;
% % end
% % 
% % %Plot mass and velocity
% % overlap = 3;
% % coupleddomain = ly + globaldomain(2) - 2*overlap*binsize(2);
% % coupledliquiddomain = coupleddomain - wallbot(2);
% % couplednbins = (globalnbins(ixyz)-1)/2+ny-overlap;
% % xaxis = 0.05:0.1:1.05; %binsize(2)/coupledliquiddomain:1/(couplednbins):1-binsize(2)/coupledliquiddomain;
% % xaxis2 = -0.025:0.05:0.975; %-1/(4*(couplednbins-1)):1/(2*(couplednbins-1)):1-1/(4*(couplednbins-1));
% % spectral_res = 6; %Ratio of points for spectral analytical solution to domain bins
% % analy_points = spectral_res*(couplednbins); % Number of spectral points
% % xaxis_analy =  0:1/(analy_points):1;
% % u_wall = 0.5*(continuum_velslice(ny+1,1) + continuum_velslice(ny+2,1));
% % xloc_MDCFDhalocell = 0.5*(xaxis(couplednbins-ny) + xaxis(couplednbins-ny+1));
% % timeratio = delta_t/continuum_delta_t;
% % 
% % %DELETE THESE LINES
% % %DELETE THESE LINES
% % continuum_velslice = continuum_velslice/u_wall;
% % u_wall = 1;
% % %DELETE THESE LINES
% % %DELETE THESE LINES
% % 
% % %Write a sequence of frames to a compressed AVI file, couette.avi:
% % % Prepare the new file.
% % %%vidObj = VideoWriter('couette.avi');
% % %vidObj.FrameRate = 30;
% % %open(vidObj);
% % 
% % m = 100; t_ave = 100; %Average over timesteps
% % set(0,'currentfigure',fig2)
% % for i = 1:Nvel_records
% %     
% %     %Plot density fluctuations in domain
% %     %set(0,'currentfigure',fig1)
% %     %plot(xaxis,ave_density_slice(:,i),'--');
% % 
% %     %Time -0.5 to record value at half time interval
% %     t = (m-0.5)*delta_t*Nmass_ave*tplot;
% %     analy = couette_analytical_fn(t,Re,[u_wall,0], ...
% %         couplednbins*2*binsize(2), ...
% %         spectral_res*couplednbins,'top');
% %     plot(xaxis_analy,analy,'r','LineWidth',5);
% %     set(gca,'FontSize',20)
% %     hold on
% %     %plot molecular velocity profile
% %     plot(xaxis(1:(globalnbins(ixyz)-1)/2),ave_vel_slice(:,1,m)/u_wall,'x','LineWidth',5,'Color',[.5 .5 .5],'MarkerSize',20);
% %     %plot(xaxis(1:(globalnbins(ixyz)-1)/2),mean(ave_vel_slice(:,1,m-t_ave/2:m+t_ave/2),3)/u_wall,'x','LineWidth',5,'Color',[.0 .5 .0],'MarkerSize',20);
% % 
% %     %plot(xaxis2(1:(globalnbins(ixyz))),vel_slice(:,1,m)/u_wall,'.','Color',[.5 .5 .5]);
% %     axis([-0.1 1.1 -0.1 1.1]);
% %     
% %     %Plot continuum velocity profile
% %     plot(xaxis(((globalnbins(ixyz)-1)/2)+1-overlap:end-1),continuum_velslice(2:ny+1,m)/u_wall,'s','Color',[.5 .5 .5],'MarkerSize',20,'LineWidth',5);
% %     %Plot Halo values
% %     plot(xaxis(end),continuum_velslice(end,m+1)/u_wall,'o','Color',[.5 .5 .5],'MarkerSize',20,'LineWidth',3,'HandleVisibility','off');
% %     plot(xaxis(((globalnbins(ixyz)-1)/2)-overlap),continuum_velslice(1,m+1)/u_wall,'o','Color',[.5 .5 .5],'MarkerSize',20,'LineWidth',3);
% % 
% %     legend ('Analytical','MD','CFD','CFD_{halo}','location','NorthWest'); legend('boxoff')
% %     xlabel('y/H'); ylabel('U_x/U')
% %     plottitle=num2str(t,'%10.6f');
% %     title(strcat('Plot after  ',plottitle,' time units'));
% %     hold off
% %     pause(0.5);
% %     
% %     %Store pictures and videos
% %    %savefig(strcat('Velocity_',num2str(i)),'png')
% %    %currFrame = getframe(gcf);
% %    %writeVideo(vidObj,currFrame);
% %     
% %     %Increase plot counter and check for end of file
% %     if (m == Nvel_records-3)
% %         break
% %     end
% %     m = m+20
% %     if (m > Nvel_records-3)
% %         m = Nvel_records-3
% %     end
% % end
% % 
% % % Close the file.
% % %close(vidObj)
% % 
% % %movefile ./couette.avi ./../results
% % % 
% % % for m = 1:Ncontinuum_records
% % %      z = [continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, continuum_velslice(2:ny+1,m)/u_wall, ...
% % %           continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, ...
% % %           continuum_velslice(2:ny+1,m)/u_wall];
% % %     imagesc(z);
% % %     colormap(gray)
% % %   
% % % 
% % %     [x,y] = meshgrid(xaxis(((globalnbins(ixyz)-1)/2)+1-overlap:end-1),0:1/6:1);
% % %     z = [continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, continuum_velslice(2:ny+1,m)/u_wall, ...
% % %          continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, ...
% % %          continuum_velslice(2:ny+1,m)/u_wall];
% % %     mesh(x,y,z')
% % %     gray
% % %     axis([0 1 0 1 0 1])
% % %      pause(0.01)  
% % %  end


