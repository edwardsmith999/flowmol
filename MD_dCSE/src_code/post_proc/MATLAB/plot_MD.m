%Plot all output files from simulation
close all
clear all

%Find results files
resultfile_dir = './../../results';

%Read Header file
read_header

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
%continuum_Domain_setup


%Plot walls, thermostatted region and sliding vectors
%plot_domain

%Plot Virial Stress
%read_pvirial

%Calculate Reynolds number
%visc = 1.6;
%L = 1; %liquiddomain(velocity_outflag)/globaldomain(velocity_outflag);
%Re = 0.8*L*max(wallslidev)/visc;

scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3) scrsz(4)]);

%Calculate spanwise direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    %Find location of smallest dimension and take slice through this
    ixyz = find(nbinsliquid==min(nbinsliquid));
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
    mass_slice = squeeze(sum(sum(mass_bins(:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
    vel_slice =  squeeze(sum(sum(vel_bins(:,:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
end

%Average same as slice
t_ave = 200;
%ave_vel_slice(:,:,1) = mean(vel_slice(:,:,    0 + 1:  t_ave),3);
ave_vel_slice(:,:,1) = mean(vel_slice(:,:,  t_ave+1:2*t_ave),3);
ave_vel_slice(:,:,2) = mean(vel_slice(:,:,2*t_ave+1:3*t_ave),3);
ave_vel_slice(:,:,3) = mean(vel_slice(:,:,3*t_ave+1:4*t_ave),3);
ave_vel_slice(:,:,4) = mean(vel_slice(:,:,4*t_ave+1:5*t_ave),3);
%ave_mass_slice(:,1) = mean(mass_slice(:,    0 + 1:  t_ave),2);
ave_mass_slice(:,1) = mean(mass_slice(:,  t_ave+1:2*t_ave),2);
ave_mass_slice(:,2) = mean(mass_slice(:,2*t_ave+1:3*t_ave),2);
ave_mass_slice(:,3) = mean(mass_slice(:,3*t_ave+1:4*t_ave),2);
ave_mass_slice(:,4) = mean(mass_slice(:,4*t_ave+1:5*t_ave),2);

%Average velocity per molecule
for n =1:3
    vel_slice(:,n,:) = squeeze(vel_slice(:,n,:))./mass_slice;
    ave_vel_slice(:,n,:) = squeeze(ave_vel_slice(:,n,:))./ave_mass_slice;
end

%Plot mass and velocity
overlap = 3;
xaxis = -1/8:1/8:1+1/8; %binsize(2)/coupledliquiddomain:1/(couplednbins):1-binsize(2)/coupledliquiddomain;
xaxis2 = -0.025:0.05:0.975; %-1/(4*(couplednbins-1)):1/(2*(couplednbins-1)):1-1/(4*(couplednbins-1));
u_0 = 2; t_0 = 160;

%Analytical Solution
spectral_res = 6;
viscosity = 2.14;
Re = density*u_0*globaldomain(ixyz)/viscosity;
analy_points = 20;%spectral_res*(nbins(ixyz)); % Number of spectral points
xaxis_analy = 0:1/21:0.97; %As the liquid solid interaction is still 1, the domain is moved in slightly due to molecular sticking at the wall

%Write a sequence of frames to a compressed AVI file, couette.avi:
% Prepare the new file.
%vidObj = VideoWriter('couette.avi');
%vidObj.FrameRate = 10;
%open(vidObj);

%Plot density fluctuations in domain
set(0,'currentfigure',fig1)
plot(xaxis,mean(density_slice(:,:),2),'--');

m = 1; %Average over timesteps
set(0,'currentfigure',fig2)
load 'digitised_oconnell'
for i = 1:Nvel_records
    
    %Time -0.5 to record value at half time interval
    set(0,'currentfigure',fig2)
    t = (m-0.5)*delta_t*Nmass_ave*tplot*t_ave;
    analy = startup_plate_analytical_fn(t,t_0,viscosity,density,u_0,liquiddomain(ixyz),analy_points);
    plot(xaxis_analy,analy/u_0,'r','LineWidth',5);
    hold on
    %plot(digitised_soln(:,1),digitised_soln(:,m+1),'--r','LineWidth',5);
    set(gca, 'YMinorTick', 'on')
    %plot molecular velocity profile
    plot(xaxis(1:(gnbins(ixyz))),ave_vel_slice(:,1,m)/u_0,'x','LineWidth',5,'Color',[.5 .5 .5],'MarkerSize',20);
    %plot(xaxis(1:(gnbins(ixyz)-1)/2),mean(ave_vel_slice(:,1,m-t_ave/2:m+t_ave/2),3)/u_wall,'x','LineWidth',5,'Color',[.0 .5 .0],'MarkerSize',20);
    set(gca,'FontSize',20)
    %plot(xaxis2(1:(gnbins(ixyz))),vel_slice(:,1,m)/u_wall,'.','Color',[.5 .5 .5]);
    axis([-0.2 1.2 -0.1 1.01]);
    
    %Plot continuum velocity profile
    xlabel('y/H'); ylabel('U_x/U')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    %hold off

    pause(0.1)
    
    %Store pictures and videos
   %savefig(strcat('Velocity_',num2str(i)),'png')
   %currFrame = getframe(gcf);
   %writeVideo(vidObj,currFrame);
    
    %Increase plot counter and check for end of file
    step = 1;
    if (m == Nvel_records/t_ave-1)
        break
    end
    m = m+step
    if (m > Nvel_records/t_ave-1)
        m = Nvel_records/t_ave-1
    end

end

% Close the file.
%close(vidObj)

%movefile ./couette.avi ./../../results
