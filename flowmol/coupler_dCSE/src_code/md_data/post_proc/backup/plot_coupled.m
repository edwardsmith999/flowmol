%Plot all output files from simulation
close all
clear all

%Find results files
resultfile_dir = './../results';

%Read Header file
read_header
%read_continuum_header

Read_DNS

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

if (continuum_vflag == 3)
    read_continuum_vbins
else
    read_continuum_vslice
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
%fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
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


%Coarse grain slices as 2 bins in MD = 1 bin in CFD
n = 1;
for i=2:2:globalnbins(ixyz)
    ave_vel_slice(n,:,:) = vel_slice(i,:,:) + vel_slice(i+1,:,:);
    ave_mass_slice(n,:)  = mass_slice(i,:)  + mass_slice(i+1,:);
    n = n+1;
end

%Average velocity per molecule
for n =1:3
    vel_slice(:,n,:) = squeeze(vel_slice(:,n,:))./mass_slice;
    ave_vel_slice(:,n,:) = squeeze(ave_vel_slice(:,n,:))./ave_mass_slice;
end

%Plot mass and velocity
overlap = 3;
coupleddomain = ly + globaldomain(2) - 2*overlap*binsize(2);
coupledliquiddomain = coupleddomain - wallbot(2);
couplednbins = (globalnbins(ixyz)-1)/2+ny-overlap;
xaxis = 0.05:0.1:1.05; %binsize(2)/coupledliquiddomain:1/(couplednbins):1-binsize(2)/coupledliquiddomain;
xaxis2 = -0.025:0.05:0.975; %-1/(4*(couplednbins-1)):1/(2*(couplednbins-1)):1-1/(4*(couplednbins-1));
spectral_res = 6; %Ratio of points for spectral analytical solution to domain bins
analy_points = spectral_res*(couplednbins); % Number of spectral points
xaxis_analy =  0:1/(analy_points):1;
u_wall = 0.5*(continuum_velslice(ny+1,1) + continuum_velslice(ny+2,1));
xloc_MDCFDhalocell = 0.5*(xaxis(couplednbins-ny) + xaxis(couplednbins-ny+1));
timeratio = delta_t/continuum_delta_t;

%Write a sequence of frames to a compressed AVI file, couette.avi:
% Prepare the new file.
%vidObj = VideoWriter('couette.avi');
%vidObj.FrameRate = 10;
%open(vidObj);

m = 1; t_ave = 1; %Average over timesteps
set(0,'currentfigure',fig2)
for i = 1:Nvel_records
    
    %Plot density fluctuations in domain
    %set(0,'currentfigure',fig1)
    %plot(xaxis,ave_density_slice(:,i),'--');

    %Time -0.5 to record value at half time interval
    t = (m-0.5)*delta_t*Nmass_ave*tplot;
    analy = couette_analytical_fn(t,Re,[u_wall,0], ...
        couplednbins*2*binsize(2), ...
        spectral_res*couplednbins,'top');
    plot(xaxis_analy,analy,'r','LineWidth',5);
    set(gca,'FontSize',20)
    hold on
    %plot molecular velocity profile
    plot(xaxis(1:(globalnbins(ixyz)-1)/2),ave_vel_slice(:,1,m)/u_wall,'x','LineWidth',5,'Color',[.5 .5 .5],'MarkerSize',20);
    %plot(xaxis(1:(globalnbins(ixyz)-1)/2),mean(ave_vel_slice(:,1,m-t_ave/2:m+t_ave/2),3)/u_wall,'x','LineWidth',5,'Color',[.0 .5 .0],'MarkerSize',20);

    %plot(xaxis2(1:(globalnbins(ixyz))),vel_slice(:,1,m)/u_wall,'.','Color',[.5 .5 .5]);
    axis([-0.1 1.1 -0.1 1.1]);
    
    %Plot continuum velocity profile
    plot(xaxis(((globalnbins(ixyz)-1)/2)+1-overlap:end-1),continuum_velslice(2:ny+1,m)/u_wall,'s','Color',[.5 .5 .5],'MarkerSize',20,'LineWidth',5);
    %Plot Halo values
    plot(xaxis(end),continuum_velslice(end,m+1)/u_wall,'o','Color',[.5 .5 .5],'MarkerSize',20,'LineWidth',3,'HandleVisibility','off');
    plot(xaxis(((globalnbins(ixyz)-1)/2)-overlap),continuum_velslice(1,m+1)/u_wall,'o','Color',[.5 .5 .5],'MarkerSize',20,'LineWidth',3);

    legend ('Analytical','MD','CFD','CFD_{halo}','location','NorthWest'); legend('boxoff')
    xlabel('y/H'); ylabel('U_x/U')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    hold off

    pause(0.5)
    
    %Store pictures and videos
   %savefig(strcat('Velocity_',num2str(i)),'png')
   %currFrame = getframe(gcf);
   %writeVideo(vidObj,currFrame);
    
    %Increase plot counter and check for end of file
    if (m > 2200)
        break
    end
    m = m+10
    if (m > Nvel_records-3)
        m = Nvel_records-3
    end
end

% Close the file.
%close(vidObj)

%movefile ./couette.avi ./../results
% 
% for m = 1:Ncontinuum_records
%      z = [continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, continuum_velslice(2:ny+1,m)/u_wall, ...
%           continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, ...
%           continuum_velslice(2:ny+1,m)/u_wall];
%     imagesc(z);
%     colormap(gray)
%   
% 
%     [x,y] = meshgrid(xaxis(((globalnbins(ixyz)-1)/2)+1-overlap:end-1),0:1/6:1);
%     z = [continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, continuum_velslice(2:ny+1,m)/u_wall, ...
%          continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall,continuum_velslice(2:ny+1,m)/u_wall, ...
%          continuum_velslice(2:ny+1,m)/u_wall];
%     mesh(x,y,z')
%     gray
%     axis([0 1 0 1 0 1])
%      pause(0.01)  
%  end


