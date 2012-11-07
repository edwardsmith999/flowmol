%Plot all output files from simulation
close all
clear all

%Select diffusive solver (0) or full DNS (1)
CFD = 1;

%Find results files
%resultfile_dir = './../results/';
resultfile_dir = '/home/es205/results/CX1_data/CPL_testing/multiple_cells/';
%resultfile_dir = '/home/es205/results/MD_continuum_results/code/coupled_couette/varying_processor_study/coupler_dCSE/src_code/';
resultfile_dir_md = strcat(resultfile_dir,'md_data/results/');
resultfile_dir_cfd = strcat(resultfile_dir,'couette_data/'); 

%Read MD Header file
resultfile_dir = resultfile_dir_md;
read_header

if (CFD == 0) 
    read_continuum_header
elseif(CFD == 1)
    %---Get CFD grid size ----
    [ngx, ngy, ngz, Lx, Ly, Lz, dx, dy, dz] = read_report(strcat(resultfile_dir_cfd,'report'));
end

% = = = Read MD header = = = 
%Check output flags and read data accordingly
if (velocity_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
elseif ((velocity_outflag > 0) & (velocity_outflag < 4) )
    read_vslice
elseif (mass_outflag == 4)
    mass_bins = read_mbins('mbins',resultfile_dir_md);
elseif ((mass_outflag > 0) & (mass_outflag < 4) )
    read_mslice
end

% = = = Read CFD data = = =
if (CFD == 0) 
	if (continuum_vflag == 3)
	    read_continuum_vbins
	else
	    read_continuum_vslice
	end
elseif(CFD == 1)
    u = Read_DNS('grid.data',resultfile_dir_cfd,ngx-2,ngy-1,ngz-2,Lx,Ly,Lz);
    continuum_velslice = u;
end



%Continuum_Domain_setup
%set(0,'currentfigure',fig1)
%read_grid(strcat(resultfile_dir_cfd,'grid.data'),[1 1 1],'plot')
%%MD domain set-up
Domain_setup 
%Plot walls, thermostatted region and sliding vectors
%set(0,'currentfigure',fig1)
%plot_domain

%Plot Virial Stress
%read_pvirial

% %Average multiple MD cells into one
% MD_cells_per_CFD = 1;
% 
% switch MD_cells_per_CFD
% case 1
%     ave_vel_slice = vel_slice;
%     ave_mass_slice = mass_slice;
% case 2
%     ave_vel_slice = zeros(size(vel_slice,1)/MD_cells_per_CFD,size(vel_slice,2),size(vel_slice,3));
%     ave_mass_slice = zeros(size(mass_slice,1)/MD_cells_per_CFD,size(mass_slice,2),size(mass_slice,3));
%     n=1;
%     for i=1:MD_cells_per_CFD:size(vel_slice,1)
%         ave_vel_slice(n,:,:) = 0.5*(vel_slice(i,:,:) + vel_slice(i+1,:,:));
%         ave_mass_slice(n,:,:) = 0.5*(mass_slice(i,:,:)+mass_slice(i,:,:));
%         n = n + 1;
%     end
% end

MD_cells_per_CFD = 1

%Plot mass and velocity
if (CFD == 0) 
    overlap = 3;
    coupleddomain(2) = ly + globaldomain(2) - 2*overlap*binsize(2);
    coupledliquiddomain = coupleddomain - wallbot(2);
    couplednbins = (gnbins(ixyz)-1)/2+ny-overlap;
    xaxis = -0.05:1/17:1.05; %binsize(2)/coupledliquiddomain:1/(couplednbins):1-binsize(2)/coupledliquiddomain;
    xaxis2 = -0.025:0.05:0.975; %-1/(4*(couplednbins-1)):1/(2*(couplednbins-1)):1-1/(4*(couplednbins-1));

    %xloc_MDCFDhalocell = 0.5*(xaxis(couplednbins-ny) + xaxis(couplednbins-ny+1));
    timeratio = delta_t/continuum_delta_t;
elseif(CFD == 1)
    cfd_binsize = binsize*MD_cells_per_CFD;
    if (abs(cfd_binsize(1)-dx)>0.001 || abs(cfd_binsize(2)-dy)>0.001)
       %error('CFD and MD cells do not agree')
       disp('CFD and MD cells do not agree')
    end
    wallsize = zeros(1,3);
    wallsize(2) = MD_cells_per_CFD*binsize(2);
    MD_domain = globaldomain - wallsize;
    CFD_domain = [Lx,Ly,Lz];
    overlap = 40*MD_cells_per_CFD*binsize;
    coupleddomain = CFD_domain + MD_domain - overlap;
    timeratio = 4; %delta_t/1.2;
end

%Get orthogonal direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    ixyz = 2;
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
end

%Analytical Solution
u_0 = 1; t_0 = 160;
spectral_res = 6;
viscosity = 1.6;
Re = density*u_0*1/viscosity;
analy_points = 20; % Number of spectral points
xaxis_analy = 0:1/20:1.00; %As the liquid solid interaction is still 1, the domain is moved in slightly due to molecular sticking at the wall

%setup figures
scrsz = get(0,'ScreenSize');
%fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/5 scrsz(4)/2]);
gax = axes('Position', [0.1, 0.1, 0.8, 0.8]);
hax = axes('Position', [.18, .48, .23, .35]);

%Write a sequence of frames to a compressed AVI file, couette.avi:
% Prepare the new file.
% vidObj = VideoWriter('couette.avi');
% vidObj.FrameRate = 10;
% open(vidObj);

intial_offset = 53.33;
t_ave = 1;  %Average over timesteps
m = 0+t_ave; %Initial Timestep
for i = 1:Nvel_records
    i

    %Read velocity
    filename = strcat(resultfile_dir_md,'/mbins');
    mass_bins = read_mbins(filename,resultfile_dir_md,m);
    filename = strcat(resultfile_dir_md,'/vbins');
    vel_bins = read_vbins(filename,resultfile_dir_md,m);

    %Calculate spanwise direction
    if (velocity_outflag < 4)
        ixyz = velocity_outflag;
    else
        %Find location of smallest dimension and take slice through this
        %ixyz = find(nbinsliquid==min(nbinsliquid));
        ixyz = 2;
        jxyz = mod(ixyz+1,3)+1;
        kxyz = mod(ixyz,3)+1;
        mass_slice = squeeze(sum(sum(mass_bins(:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
        vel_slice =  squeeze(sum(sum(vel_bins(:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
    end

    %Average velocity per molecule
    for n =1:3
        ave_vel_slice(:,n) = vel_slice(:,n)./mass_slice';
    end

    %Plot density fluctuations in domain
%     set(0,'currentfigure',fig1)
%     plot(xaxis(1:(gnbins(ixyz))),mean(density_slice(:,:),2),'--');
%     hold on
%     plot(xaxis(1:(gnbins(ixyz))),density_slice(:,m),'r-');
%     hold off

    %Time -0.5 to record value at half time interval
    t =(intial_offset+m-0.5)*2*delta_t*Nmass_ave*tplot;
    
    % =================================
    %    Plot MD region on full axis
	% =================================
    set(0,'currentfigure',fig2)
    set(fig2,'CurrentAxes',gax)

    
    % =================================
    %    Plot MD region
	% =================================
    set(0,'currentfigure',fig2)
    set(fig2,'CurrentAxes',gax)

    %Plot CFD velocity profile
    plot(linspace(MD_domain(ixyz)-overlap(ixyz)+0.5*cfd_binsize(ixyz),coupleddomain(ixyz)-0.5*cfd_binsize(ixyz),ngy-1)/coupleddomain(ixyz),continuum_velslice(:,m)/u_0,'s','Color',[.5 .5 .5],'MarkerSize',10,'LineWidth',5);
    hold on
    %plot molecular velocity profile
    plot(linspace(-0.5*cfd_binsize(ixyz),MD_domain(ixyz)-0.5*cfd_binsize(ixyz),gnbins(ixyz)/MD_cells_per_CFD)/coupleddomain(ixyz),ave_vel_slice(:,1),'x','LineWidth',5,'Color',[.2 .2 .2],'MarkerSize',10);
      
    %Plot anayltical solution
    analy = couette_analytical_fn(t,Re,[1,0],coupleddomain(ixyz), ...
                                   analy_points,'top');
    plot(xaxis_analy,analy/u_0,'k','LineWidth',4);

    %Make plot look nice
    %legend ('MD','CFD','CFD_{halo}','location','NorthEast'); legend('boxoff')
    set(gca,'FontSize',20)
    axis([-0.1 1.1 -0.1 1.1]);
    set(gca,'FontSize',20)
    axis([-0.1 0.65 -0.05 0.25]);
    xlabel('y/H'); ylabel('U_x/U')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    hold off

    % =================================
    %    Plot Insert 
	% =================================
    set(fig2,'CurrentAxes',hax)

    plot(linspace(MD_domain(ixyz)-overlap(ixyz)+0.5*cfd_binsize(ixyz),coupleddomain(ixyz)-0.5*cfd_binsize(ixyz),ngy-1)/coupleddomain(ixyz),continuum_velslice(:,m)/u_0,'s','Color',[.5 .5 .5],'MarkerSize',5,'LineWidth',2);
    hold on
	plot(linspace(-0.5*cfd_binsize(ixyz),MD_domain(ixyz)-0.5*cfd_binsize(ixyz),gnbins(ixyz)/MD_cells_per_CFD)/coupleddomain(ixyz),ave_vel_slice(:,1),'x','LineWidth',5,'Color',[.2 .2 .2],'MarkerSize',5);


    %Fix axis position to averaged continuum BC
    tave_vbc = 10;
    if (m < tave_vbc+1)
        ave_vbc = mean(continuum_velslice(1,m:m+tave_vbc));
    elseif(m > Nvel_records-tave_vbc-1)
        ave_vbc = mean(continuum_velslice(1,m-tave_vbc:m));
    else
        ave_vbc = mean(continuum_velslice(1,m-tave_vbc:m+tave_vbc));
    end
	axis([0.35 0.6 ave_vbc-0.025 ave_vbc+0.1 ])
    analy = couette_analytical_fn(t,Re,[1,0], ...
        coupleddomain(ixyz), ...
        analy_points,'top');
    %analy = startup_plate_analytical_fn(t,t_0,viscosity,density,u_0,liquiddomain(ixyz),analy_points);
    plot(xaxis_analy,analy/u_0,'k','LineWidth',4);
 	set(hax,'XTick',[]); 
    hold off

    drawnow
    pause(0.1)

    %Store pictures and videos
    %if (mod(m,Nvel_records/10) == 0)
    %    savefig(strcat('Velocity_',num2str(m)),'png')
    %end
% 	currFrame = getframe(gcf);
% 	writeVideo(vidObj,currFrame);
%     
    %Increase plot counter and check for end of file
    if (m > Nvel_records)
        break
    end
    m = m+1
    if (m > Nvel_records)
        m = Nvel_records-30
    end
end

% Close the file.
% close(vidObj)

%movefile ./couette.avi ./../results

