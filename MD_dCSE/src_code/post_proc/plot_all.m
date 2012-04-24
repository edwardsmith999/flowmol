%Plot all output files from simulation
close all
clear all

%Read Header file
read_header

%Find results files
resultfile_dir = './../results';

%Check output flags and read data accordingly
if (velocity_outflag == 4)
    read_vbins
elseif ((velocity_outflag > 0) & (velocity_outflag < 4) )
    read_vslice
elseif (mass_outflag == 4)
    read_mbins
elseif ((mass_outflag > 0) & (mass_outflag < 4) )
    read_mslice
end

%Establish and Plot domain set-up
Domain_setup

%Plot walls, thermostatted region and sliding vectors
%plot_domain

%Plot Virial Stress
%read_pvirial

%Calculate Reynolds number
visc = 1.6;
L = 1; %liquiddomain(velocity_outflag)/globaldomain(velocity_outflag);
Re = 0.8*L*max(wallslidev)/visc;

scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);

%Calculate spanwise direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    %Find location of smallest dimension and take slice through this
    ixyz = find(nbinsliquid==min(nbinsliquid));
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
    density_slice = squeeze(sum(sum(mass_bins(:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
    vel_slice =  squeeze(sum(sum(vel_bins(:,:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
end

%Plot mass and velocity     
m = 1;
xaxis = 0:1/(gnbins(ixyz)-1):1;
spectral_res = 10; %Ratio of points for spectral analytical solution to domain bins
analy_points = spectral_res*nbinsliquid(ixyz); % Number of spectral points
xaxis_analy =  wallbot(ixyz) ...
              :liquiddomain(ixyz)/analy_points ...
              :globaldomain(ixyz)-walltop(ixyz); %Liquid part of domain
xaxis_analy = xaxis_analy / globaldomain(ixyz); %Normalised by whole domain width

for i = 1:Nvel_records
    
     %Plot density fluctuations in domain    
     set(0,'currentfigure',fig1)
     plot(xaxis,density_slice(:,i),'--');
 
     %plot molecular velocity profile
     set(0,'currentfigure',fig2)
     plot(xaxis,vel_slice(:,1,i),'--x','Color',[.5 .5 .5]);
     %axis([-0.1 1.1 -0.1 1.1]);
     hold on
     %Time -0.5 to record value at half time interval
     t = (i-0.5)*delta_t*Nmass_ave*tplot;
     %analy = couette_analytical_fn(t,Re,[topwallslidev,botwallslidev], ...
     %                               liquiddomain(ixyz), ...
    %                                spectral_res*nbinsliquid(ixyz),sliding_wall);
     %plot(xaxis_analy,analy,'r');
     legend ('MD simulation','Analytical','location','NorthWest'); legend('boxoff')
     xlabel('y/H'); ylabel('U_x/U')
     plottitle=num2str(t,'%10.6f');
     title(strcat('Plot after ',plottitle,' seconds'));
%     savefig(strcat('Velocity_',num2str(i)),'png')
     hold off
     pause(0.2);
%     %M(n) = getframe(fig1); %Store Frames for video
%       if (m == Nresults-3)
%           break
%       end
%       m = 2*m
%       if (m > Nresults-3)
%           m = Nresults-3
%       end
% 
end
 

