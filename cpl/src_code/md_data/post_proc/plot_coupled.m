%Plot all output files from simulation
close all
clear all

%Select diffusive solver (0) or full DNS (1)
CFD = 1;

%Turn on/off video and picture output
savevid = 0;
savepic = 1;

%Find results files
%resultfile_dir = './../results/';
%resultfile_dir = '/home/es205/results/CX1_data/CPL_testing/slice/';
%resultfile_dir = '/home/es205/codes/coupled/coupler_dCSE/src_code/';
%resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/flekkoy/50CFDMDratio/';
%resultfile_dir = '/home/es205/results/MD_continuum_results/code/coupled_couette/varying_processor_study/coupler_dCSE/src_code/';
resultfile_dir = '/home/djt06/Documents/Academia/PhD/Code/Development/branch/coupler_dCSE/src_code/';
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
    [u,v,w] = Read_DNS('grid.data',resultfile_dir_cfd,ngx-2,ngy-1,ngz-2,Lx,Ly,Lz,3,true);
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


MD_cells_per_CFD = 2;

%Calculate properties
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
    wallsize(2) = tethdistbot(2); %2*MD_cells_per_CFD*binsize(2);
    MD_domain = globaldomain - wallsize;
    CFD_domain = [Lx,Ly,Lz];
    overlap = 4*MD_cells_per_CFD*binsize;
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
viscosity = 2.14;
Re = density*u_0*1/viscosity;
analy_points = 20; % Number of spectral points

xaxis_md  = linspace(-1.5*cfd_binsize(ixyz),MD_domain(ixyz)-0.5*cfd_binsize(ixyz),gnbins(ixyz)/MD_cells_per_CFD)/coupleddomain(ixyz);
%xaxis_md  = linspace(0,MD_domain(ixyz),gnbins(ixyz)/MD_cells_per_CFD)/coupleddomain(ixyz);
xaxis_cfd = linspace(MD_domain(ixyz)-overlap(ixyz)-0.5*cfd_binsize(ixyz),coupleddomain(ixyz)+0.5*cfd_binsize(ixyz),ngy-1)/coupleddomain(ixyz);
xaxis_analy = 0:1/20:1.00; %As the liquid solid interaction is still 1, the domain is moved in slightly due to molecular sticking at the wall

xaxis_md2  = linspace(-1.75*cfd_binsize(ixyz),MD_domain(ixyz)-0.25*cfd_binsize(ixyz),gnbins(ixyz))/coupleddomain(ixyz);

%setup figures
scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/5 scrsz(4)/2]);
gax = axes('Position', [0.1, 0.1, 0.8, 0.8]);
%hax = axes('Position', [.18, .48, .23, .35]);

%Check if stress is avilable to plot
files = dir(resultfile_dir_md);
plot_pVA = 0;
for i=1:size(files,1)
    if (strmatch('pVA',files(i).name) == 1)
        plot_pVA=1;
    end
end

%Write a sequence of frames to a compressed AVI file, couette.avi:
% Prepare the new file.
if (savevid == 1)
    vidObj = VideoWriter('velocity.avi');
    vidObj.FrameRate = 10;
    open(vidObj);
    if (plot_pVA == 1)
        vidObj2 = VideoWriter('stress.avi');
        vidObj2.FrameRate = 10;
        open(vidObj2);
    end
end

last_record = 0;
intial_offset = -1; %53.33;
t_ave = 10;  %Average over timesteps
m = 0+ceil(t_ave/2); %Initial Timestep
snaps = [2,4,8, 16, 700];
for i = 1:Nvel_records
    i
    
    % =================================
    %    Plot COUPLED velocity
    % =================================
    %Read velocity
    filename = strcat(resultfile_dir_md,'/mbins');
    mass_bins = read_mbins(filename,resultfile_dir_md,m-2); %NOTE WE NEED minus 2 here not sure why yet!
    filename = strcat(resultfile_dir_md,'/vbins');
    vel_bins = read_vbins(filename,resultfile_dir_md,m-2); %NOTE WE NEED minus 2 here not sure why yet!
    
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
    
    %Average multiple MD cells into one
    switch MD_cells_per_CFD
        case 1
            ave_vel_slice = vel_slice;
            ave_mass_slice = mass_slice;
        case 2
            ave_vel_slice = zeros(size(vel_slice,1)/MD_cells_per_CFD,size(vel_slice,2));
            ave_mass_slice = zeros(size(mass_slice,2)/MD_cells_per_CFD,size(mass_slice,1));
            n=1;
            for icell=1:MD_cells_per_CFD:size(vel_slice,1)
                ave_vel_slice(n,:)  = 0.5*(vel_slice(icell,:) + vel_slice(icell+1,:));
                ave_mass_slice(n)   = 0.5*(mass_slice(icell)  + mass_slice(icell+1));
                n = n + 1;
            end
        case 4
            ave_vel_slice = zeros(size(vel_slice,1)/MD_cells_per_CFD,size(vel_slice,2));
            ave_mass_slice = zeros(size(mass_slice,2)/MD_cells_per_CFD,size(mass_slice,1));
            n=1;
            for icell=1:MD_cells_per_CFD:size(vel_slice,1)
                ave_vel_slice(n,:) = 0.25*( vel_slice(icell  ,:) ...
                    +vel_slice(icell+1,:) ...
                    +vel_slice(icell+2,:) ...
                    +vel_slice(icell+3,:));
                ave_mass_slice(n) = 0.25*(mass_slice(icell  ) ...
                    +mass_slice(icell+1) ...
                    +mass_slice(icell+2) ...
                    +mass_slice(icell+3));
                n = n + 1;
            end
    end
    
    
    %Average velocity per molecule
    for n =1:3
        ave_vel_slice(:,n) = ave_vel_slice(:,n)./ave_mass_slice;
        vel_slice(:,n) = vel_slice(:,n)./mass_slice';
    end

    %Plot density fluctuations in domain
    %set(0,'currentfigure',fig1)
    %plot(xaxis_md,2*ave_mass_slice/(Nmass_ave*dx*dy*dz),'x','Color',[.5 .5 .5],'MarkerSize',8,'LineWidth',5);
    %axis([-0.1 1.1 0.7 1.0 ])

    %Time -0.5 to record value at half time interval
    %if (i  == 1)
    %    t = 0;
    %else
        t =(intial_offset+m-0.5)*delta_t*Nmass_ave*tplot;
    %end
    
    % =================================
    %    Plot MD region on full axis
    % =================================
    set(0,'currentfigure',fig2)
    set(fig2,'CurrentAxes',gax)
    
    %Plot CFD velocity profile
    plot(xaxis_cfd,continuum_velslice(:,m)/u_0,'s','Color',[.5 .5 .5],'MarkerSize',8,'LineWidth',5);
    hold on

    %plot molecular velocity profile
    plot(xaxis_md(1:end),ave_vel_slice(1:end,1),'x','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10);
    %plot(xaxis_md2,vel_slice(:,1),'^-')
    
    %Plot anayltical solution
    analy = couette_analytical_fn(t,Re,[1,0],coupleddomain(ixyz),analy_points,'top');
    plot(xaxis_analy,analy/u_0,'k','LineWidth',5);
    
    %Make plot look nice
    legend ('CFD','MD','location','NorthWest'); legend('boxoff')
    set(gca,'FontSize',20)
    axis([-0.1 1.1 -0.1 1.1]);
    xlabel('y/H'); ylabel('U_x/U')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    hold off
    
    % =================================
    %    Plot Insert
    % =================================
    %set(fig2,'CurrentAxes',hax)
    %set(fig2,'CurrentAxes',hax)
    %plot(abs(ave_vel_slice(:,2)));
    %axis([0 12 -0.01 0.2])
    %
    %     plot(linspace(MD_domain(ixyz)-overlap(ixyz)+0.5*cfd_binsize(ixyz),coupleddomain(ixyz)-0.5*cfd_binsize(ixyz),ngy-1)/coupleddomain(ixyz),continuum_velslice(:,m)/u_0,'s','Color',[.5 .5 .5],'MarkerSize',5,'LineWidth',2);
    %     hold on
    % 	plot(linspace(-0.5*cfd_binsize(ixyz),MD_domain(ixyz)-0.5*cfd_binsize(ixyz),gnbins(ixyz)/MD_cells_per_CFD)/coupleddomain(ixyz),ave_vel_slice(:,1),'x','LineWidth',5,'Color',[.2 .2 .2],'MarkerSize',5);
    %
    %
    %     %Fix axis position to averaged continuum BC
    %     tave_vbc = 10;
    %     if (m < tave_vbc+1)
    %         ave_vbc = mean(continuum_velslice(1,m:m+tave_vbc));
    %     elseif(m > Nvel_records-tave_vbc-1)
    %         ave_vbc = mean(continuum_velslice(1,m-tave_vbc:m));
    %     else
    %         ave_vbc = mean(continuum_velslice(1,m-tave_vbc:m+tave_vbc));
    %     end
    % 	axis([0.35 0.6 ave_vbc-0.025 ave_vbc+0.1 ])
    %  	set(hax,'XTick',[]);
    %hold off

    %Store pictures and videos
    if (savepic == 1)
        if (mod(m,Nvel_records/10) == 0)
            savefig(strcat('Velocity_',num2str(m)),'png')
        end
    elseif (savevid == 1)
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end
    
    drawnow
    %pause(0.5)

    % =================================
    %    Plot COUPLED stress
    % =================================

    if ( plot_pVA == 1)
        set(0,'currentfigure',fig1)

        cumulative_PVA = 0;
        filename = strcat(resultfile_dir_md,'/pVA');
        for ave = ceil(-floor(t_ave/2)):floor(t_ave/2)
            [pressure_VA]=read_pVA(m+ave-2,resultfile_dir_md,filename);
            cumulative_PVA = cumulative_PVA + pressure_VA;
        end
        pressure_VA = cumulative_PVA/t_ave;

        %Calculate spanwise direction
        if (velocity_outflag < 4)
            ixyz = velocity_outflag;
        else
            %Find location of smallest dimension and take slice through this
            %ixyz = find(nbinsliquid==min(nbinsliquid));
            ixyz = 2;
            jxyz = mod(ixyz+1,3)+1;
            kxyz = mod(ixyz,3)+1;
            P_slice =  squeeze(sum(sum(pressure_VA(:,:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
        end

        %Average multiple MD cells into one
        switch MD_cells_per_CFD
            case 1
                ave_P_slice = P_slice;
            case 2
                ave_P_slice = zeros(size(P_slice,1)/MD_cells_per_CFD,size(P_slice,2),size(P_slice,3));
                n=1;
                for icell=1:MD_cells_per_CFD:size(P_slice,1)
                    ave_P_slice(n,:,:)  = 0.5*(P_slice(icell,:,:) + P_slice(icell+1,:,:));
                    n = n + 1;
                end
            case 4
                ave_P_slice = zeros(size(P_slice,1)/MD_cells_per_CFD,size(P_slice,2),size(P_slice,3));
                n=1;
                for icell=1:MD_cells_per_CFD:size(P_slice,1)
                    ave_P_slice(n,:,:) = 0.25*( P_slice(icell  ,:,:) ...
                        +P_slice(icell+1,:,:) ...
                        +P_slice(icell+2,:,:) ...
                        +P_slice(icell+3,:,:));
                    n = n + 1;
                end
        end


        %Average velocity per molecule
        for n =1:3
            ave_P_slice(:,n) = ave_P_slice(:,n)./ave_mass_slice;
            P_slice(:,n) = P_slice(:,n)./mass_slice';
        end


        %Plot CFD shear stress
        plot(xaxis_cfd(2:end),squeeze(stress(1,2,2:end,m)),'s','Color',[.5 .5 .5],'MarkerSize',8,'LineWidth',5)
        hold on

        %Plot MD shear stress
        plot(xaxis_md,-squeeze(ave_P_slice(:,1,2)),'x','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10)
        plot(xaxis_md2,-squeeze(P_slice(:,1,2)),'^-')

        %Plot anayltical solution
        clear analy
        analy = couette_analytical_stress_fn(t,Re,1,coupleddomain(ixyz),analy_points);
        plot(xaxis_analy,analy/Re,'k','LineWidth',5);

        %Make plot look nice
        legend ('CFD','MD','location','NorthWest'); legend('boxoff')
        set(gca,'FontSize',20)
        axis([-0.1 1.1 -0.01 0.1]);
        xlabel('y/H'); ylabel('\Pi')
        plottitle=num2str(t,'%10.6f');
        title(strcat('Plot after  ',plottitle,' time units'));
        hold off

        drawnow

        %Store pictures and videos
        if (savepic == 1)
            if (mod(m,Nvel_records/10) == 0)
                savefig(strcat('Stress_',num2str(m)),'png')
            end
        elseif (savevid == 1)
            currFrame = getframe(gcf);
            writeVideo(vidObj2,currFrame);
        end

    end


    % =================================
    %    Plot COUPLED results
    % =================================
    %Increase plot counter and check for end of file
    if (last_record == 1)
        break
    end
    m = m + 1 %snaps(i)
    %m = t_array(i)
    if (m > Nvel_records-ceil(t_ave/2))
        m = Nvel_records-ceil(t_ave/2)
        last_record = 1;
    end



end

% Close the Video file.
if (savevid == 1)
    close(vidObj)
    if (plot_pVA == 1)
        close(vidObj2)
    end
 end

%movefile ./couette.avi ./../results

