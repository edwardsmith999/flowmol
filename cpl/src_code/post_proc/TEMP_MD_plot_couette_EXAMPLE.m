%Plot all output files from simulation
close all
clear all

%Turn on/off video and picture output
savevid = 1;
savepic = 0;

%Wall Normal Direction
wall_normal = 2;

%Choose level of detail to plot
% 1 - density and velocity
% 2 - density, velocity and VA stresses/fluxes
% 3 - density, velocity and VA/CV stresses/fluxes
plot_level = 3;

%Find results files
resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/flekkoy/MD_only_testing';

fid = fopen(strcat(resultfile_dir,'/couette_stress_analy'),'r','n');
applied_shear = fread(fid,'double');
fclose(fid);


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


%%MD domain set-up
Domain_setup

%Calculate properties
MD_cells_per_ave = 1;
wallsize = zeros(1,3);
wallsize(2) = binsize(2);
MD_domain = globaldomain - wallsize;

%Get orthogonal direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    ixyz = wall_normal;
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
end

%Analytical Solution
u_0 = 1; t_0 = 160;
spectral_res = 6;
viscosity = 1.6;
Re = density*u_0*1/viscosity;
analy_points = 24; % Number of spectral points

xaxis_md  = linspace(0,0.5*globaldomain(ixyz),gnbins(ixyz)/MD_cells_per_ave)/globaldomain(ixyz);
xaxis_analy = linspace(0,1,analy_points); %As the liquid solid interaction is still 1, the domain is moved in slightly due to molecular sticking at the wall

%setup figures
scrsz = get(0,'ScreenSize');
fig3 = figure('Position',[     1     scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
gax = axes('Position', [0.1, 0.1, 0.8, 0.8]);
hax = axes('Position', [.18, .48, .23, .35]);

if (plot_level > 1)
	fig1 = figure('Position',[scrsz(3)*2/6 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
end


%Check if stress is avilable to plot
if (plot_level > 1)
    files = dir(resultfile_dir);
    plot_pVA = 0;
    for i=1:size(files,1)
        if (strmatch('pVA',files(i).name) == 1)
            plot_pVA=1;
        end
    end
    if (savevid == 1)
        if (plot_pVA == 1)
            vidObj2 = VideoWriter('stress.avi');
            vidObj2.FrameRate = 10;
            open(vidObj2);
        end
    end

	if (plot_level > 2)
        %Check if CV is avilable to plot
        files = dir(resultfile_dir);
        plot_vflux = 0;
        for i=1:size(files,1)
            if (strmatch('vflux',files(i).name) == 1)
                plot_vflux=1;
            end
        end
	end
end

%Write a sequence of frames to a compressed AVI file, couette.avi:
% Prepare the new file.
if (savevid == 1)
    vidObj = VideoWriter('velocity.avi');
    vidObj.FrameRate = 10;
    open(vidObj);
end

last_record = 0;
intial_offset = 50; %53.33;
t_ave = 100;  %Average over timesteps
m = 0+ceil(t_ave/2); %Initial Timestep
snaps = [2,4,8, 16, 700];
for i = 1:Nvel_records
    i
    
    % =================================
    %    Plot velocity
    % =================================
    %Find location of smallest dimension and take slice through this
    ixyz = wall_normal;
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
    %Read MD velocity
    if (velocity_outflag == 4)
        vel_slice = read_vslice('./vslice',resultfile_dir);
        filename = strcat(resultfile_dir,'/mbins');
        mass_bins = read_mbins(filename,resultfile_dir,m-2); %NOTE WE NEED minus 2 here not sure why yet!
        filename = strcat(resultfile_dir,'/vbins');
        vel_bins = read_vbins(filename,resultfile_dir,m-2); %NOTE WE NEED minus 2 here not sure why yet!
        mass_slice = squeeze(sum(sum(mass_bins(:,:,:)  ,jxyz),kxyz)/(gnbins(jxyz)*gnbins(kxyz)));
        vel_slice  = squeeze(sum(sum(vel_bins(:,:,:,:) ,jxyz),kxyz)/(gnbins(jxyz)*gnbins(kxyz)));
    elseif ((velocity_outflag > 0) && (velocity_outflag < 4) )
        ixyz = velocity_outflag;
        mass_slice = read_mslice('mslice',resultfile_dir,m-2);
        vel_slice  = read_vslice('vslice',resultfile_dir,m-2);
    end

    %Average multiple MD cells into one
    switch MD_cells_per_ave
        case 1
            ave_vel_slice = vel_slice;
            ave_mass_slice = mass_slice;
        case 2
            ave_vel_slice  = zeros(size(vel_slice,1)/MD_cells_per_ave,size(vel_slice,2));
            ave_mass_slice = zeros(size(mass_slice,2)/MD_cells_per_ave,size(mass_slice,1));
            n=1;
            for icell=1:MD_cells_per_ave:size(vel_slice,1)
                ave_vel_slice(n,:)  = 0.5*(vel_slice(icell,:) + vel_slice(icell+1,:));
                ave_mass_slice(n)   = 0.5*(mass_slice(icell)  + mass_slice(icell+1));
                n = n + 1;
            end
        case 4
            ave_vel_slice = zeros(size(vel_slice,1)/MD_cells_per_ave,size(vel_slice,2));
            ave_mass_slice = zeros(size(mass_slice,2)/MD_cells_per_ave,size(mass_slice,1));
            n=1;
            for icell=1:MD_cells_per_ave:size(vel_slice,1)
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
        vel_slice(:,n) = vel_slice(:,n)./mass_slice;
    end

    %Plot density fluctuations in domain
    set(0,'currentfigure',fig3)
    plot(xaxis_md,mass_slice/(gnbins(jxyz)*gnbins(kxyz)*Nmass_ave*prod(binsize)),'-x','Color',[.5 .5 .5],'MarkerSize',8,'LineWidth',5);
    axis([-0.1 1.1 0.5 1.1 ])

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
    
    %plot molecular velocity profile
    plot(xaxis_md(1:end),ave_vel_slice(1:end,1),'x','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10);
    hold on
    
    %Plot anayltical solution
    analy = couette_analytical_fn(t,Re,[u_0,0],2*globaldomain(ixyz),analy_points,'top');
    plot(xaxis_analy,analy,'k','LineWidth',5);
    
    %Make plot look nice
    legend ('MD','analy','location','NorthEast'); legend('boxoff')
    set(gca,'FontSize',20)
    axis([-0.1 1.1 -0.1 1.1]);
    xlabel('y/H'); ylabel('U_x/U')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    hold off
    
        % =================================
        %  Plot normal velocity on insert
        % =================================
        set(fig2,'CurrentAxes',hax)
        plot(xaxis_md,abs(ave_vel_slice(:,2)));
        axis([-0.1 1.1 -0.01 0.1])
        set(hax,'XTick',[]);

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
    %    Plot stress
    % =================================
    if (plot_level > 1)
        if ( plot_pVA == 1)
            set(0,'currentfigure',fig1)

            %Get cumulative pressure
            cumulative_PVA = 0;
            cumulative_convection = 0;
            filename = strcat(resultfile_dir,'/pVA');
            for ave = ceil(-floor(t_ave/2)):floor(t_ave/2)
                %Add the pressures
                [pressure_VA]=read_pVA(m+ave-2,resultfile_dir,filename);
                cumulative_PVA = cumulative_PVA + pressure_VA;  
                %Cumulative Convection
                mass_slice = read_mslice('mslice',resultfile_dir,m+ave-2);
                vel_slice  = read_vslice('vslice',resultfile_dir,m+ave-2);
                for n=1:size(vel_slice,1)
                    rhovv(n,:,:)= density*(vel_slice(n,:)./mass_slice(n))'*(vel_slice(n,:)./mass_slice(n));
                end
                cumulative_convection =  cumulative_convection + rhovv;
            end
            pressure_VA = cumulative_PVA/t_ave;
            cumulative_convection = cumulative_convection/t_ave;

            %Find location of smallest dimension and take slice through this
            %ixyz = find(nbinsliquid==min(nbinsliquid));
            ixyz = wall_normal;
            jxyz = mod(ixyz+1,3)+1;
            kxyz = mod(ixyz,3)+1;
            P_slice =  squeeze(sum(sum(pressure_VA(:,:,:,:,:),jxyz),kxyz)/(gnbins(jxyz)*gnbins(kxyz)));


            %Average multiple MD cells into one
            switch MD_cells_per_ave
                case 1
                    ave_P_slice = P_slice;
                    ave_cumulative_convection = cumulative_convection;
                case 2
                    ave_P_slice = zeros(size(P_slice,1)/MD_cells_per_ave,size(P_slice,2),size(P_slice,3));
                    n=1;
                    for icell=1:MD_cells_per_ave:size(P_slice,1)
                        ave_P_slice(n,:,:)  = 0.5*(P_slice(icell,:,:) + P_slice(icell+1,:,:));
                        n = n + 1;
                    end
                case 4
                    ave_P_slice = zeros(size(P_slice,1)/MD_cells_per_ave,size(P_slice,2),size(P_slice,3));
                    n=1;
                    for icell=1:MD_cells_per_ave:size(P_slice,1)
                        ave_P_slice(n,:,:) = 0.25*( P_slice(icell  ,:,:) ...
                            +P_slice(icell+1,:,:) ...
                            +P_slice(icell+2,:,:) ...
                            +P_slice(icell+3,:,:));
                        n = n + 1;
                    end
            end


            %Average stress per molecule
            %for n1 =1:3
            %for n2 =1:3
            %    ave_P_slice(:,n1,n2) = ave_P_slice(:,n1,n2)./ave_mass_slice;
            %    ave_cumulative_convection(:,n1,n2) = ave_cumulative_convection(:,n1,n2)./ave_mass_slice;
            %    P_slice(:,n1,n2) = P_slice(:,n1,n2)./mass_slice;
            %end
            %end


            %Plot MD shear stress
            %plot(xaxis_md,-squeeze(density*ave_vel_slice(:,1).*ave_vel_slice(:,2)),'o','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10)
            plot(xaxis_md,-squeeze(ave_P_slice(:,1,2)),'s','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10)
            hold on
            plot(xaxis_md,-ave_cumulative_convection(:,1,2),'o','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10)

            %plot(xaxis_md,,'x','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10)
            %plot(xaxis_md,-squeeze(ave_P_slice(:,1,2))+ave_cumulative_convection(:,1,2),'s','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10)


            %Plot anayltical solution
            clear analy
            analy = couette_analytical_stress_fn(t,Re,u_0,2*globaldomain(ixyz),analy_points);
            plot(xaxis_analy,analy,'k','LineWidth',5);

            %Plot applied stress
            plot(xaxis_md(end),-2*applied_shear((intial_offset+m-0.5)*Nmass_ave*tplot),'x','MarkerSize',20)

            %Make plot look nice
            legend ('flux','stress','analy','location','NorthWest'); legend('boxoff')
            set(gca,'FontSize',20)
            axis([-0.3 1.1 -0.1 0.8]);
            xlabel('y/H'); ylabel('\Pi')
            plottitle=num2str(t,'%10.6f');
            title(strcat('Plot after  ',plottitle,' time units'));
            hold off


            if (plot_vflux == 0)
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

        end
    end


    % =================================
    %    Plot CV stress
    % =================================
    if (plot_level > 2)
        if ( plot_vflux == 1)
            set(0,'currentfigure',fig1)

            cumulative_vflux = 0;
            filename1 = strcat(resultfile_dir,'/vflux');
            filename2 = strcat(resultfile_dir,'/psurface');
            filename3 = strcat(resultfile_dir,'/vsnap');
            cumulative_vflux    = 0;
            cumulative_psurface = 0;
            for ave = ceil(-floor(t_ave/2)):floor(t_ave/2)
                [velocity_snapshot,velocity_flux,pressure_surface] = ...
                    read_vflux(m+ave-2,resultfile_dir,gnbins,nd,filename1,filename2,filename3);
                cumulative_vflux    = cumulative_vflux + velocity_flux; 
                cumulative_psurface = cumulative_psurface + pressure_surface;  
            end
            cumulative_vflux    = cumulative_vflux/(t_ave);
            cumulative_psurface = cumulative_psurface/(t_ave);
            [velocity_snapshot,velocity_flux,pressure_surface] = ...
                    read_vflux(m-2,resultfile_dir,gnbins,nd,filename1,filename2,filename3);

            %Find location of smallest dimension and take slice through this
            ixyz = wall_normal;
            jxyz = mod(ixyz+1,3)+1;
            kxyz = mod(ixyz,3)+1;
            vflux_slice    =  squeeze(sum(sum(cumulative_vflux(:,:,:,:,:),jxyz),kxyz)/(gnbins(jxyz)*gnbins(kxyz)));
            psurface_slice =  squeeze(sum(sum(cumulative_psurface(:,:,:,:,:),jxyz),kxyz)/(gnbins(jxyz)*gnbins(kxyz)));

            %Average multiple MD cells into one
            switch MD_cells_per_ave
                case 1
                    ave_vflux_slice = vflux_slice;
                    ave_psurface_slice = psurface_slice;
                case 2
                    ave_vflux_slice = zeros(size(vflux_slice,1)/MD_cells_per_ave,size(vflux_slice,2),size(vflux_slice,3));
                    ave_psurface_slice = zeros(size(psurface_slice,1)/MD_cells_per_ave,size(psurface_slice,2),size(psurface_slice,3));
                    n=1;
                    for icell=1:MD_cells_per_ave:size(vflux_slice,1)
                        ave_vflux_slice(n,:,:)  = 0.5*(vflux_slice(icell,:,:) + vflux_slice(icell+1,:,:));
                        ave_psurface_slice(n,:,:)  = 0.5*(psurface_slice(icell,:,:) + psurface_slice(icell+1,:,:));
                        n = n + 1;
                    end
                case 4
                    ave_vflux_slice = zeros(size(vflux_slice,1)/MD_cells_per_ave,size(vflux_slice,2),size(vflux_slice,3));
                    n=1;
                    for icell=1:MD_cells_per_ave:size(vflux_slice,1)
                        ave_vflux_slice(n,:,:) = 0.25*( vflux_slice(icell  ,:,:) ...
                            +vflux_slice(icell+1,:,:) ...
                            +vflux_slice(icell+2,:,:) ...
                            +vflux_slice(icell+3,:,:));
                        ave_psurface_slice(n,:,:) = 0.25*( psurface_slice(icell  ,:,:) ...
                            +psurface_slice(icell+1,:,:) ...
                            +psurface_slice(icell+2,:,:) ...
                            +psurface_slice(icell+3,:,:));
                        n = n + 1;
                    end
            end


            %Average stress per molecule
            for n =1:3
                ave_vflux_slice(:,n) = ave_vflux_slice(:,n)./(binsize(jxyz)*binsize(kxyz));
                ave_psurface_slice(:,n) = ave_psurface_slice(:,n)./(binsize(jxyz)*binsize(kxyz));
            end


            %Plot MD shear stress
            hold on
            plot(xaxis_md,-squeeze(ave_vflux_slice(:,1,2)),   '--','LineWidth',3,'Color',[.5 .5 .5],'MarkerSize',10)
            plot(xaxis_md,-squeeze(ave_psurface_slice(:,1,2)),'-','LineWidth',3,'Color',[.5 .5 .5],'MarkerSize',10)
            %plot(xaxis_md,-squeeze(ave_psurface_slice(:,1,2))+squeeze(ave_vflux_slice(:,1,2)),'-','LineWidth',3,'Color',[.5 .5 .5],'MarkerSize',10)

            %Plot anayltical solution
            clear analy
            analy = couette_analytical_stress_fn(t,Re,u_0,2*globaldomain(ixyz),analy_points);
            plot(xaxis_analy,analy,'k','LineWidth',5);

            %Make plot look nice
            legend ('flux','stress','analy','location','NorthWest'); legend('boxoff')
            set(gca,'FontSize',20)
            axis([-0.3 1.1 -0.01 0.075]);
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
    end


    % =================================
    %    Increase current output
    % =================================
    %Increase plot counter and check for end of file
    if (last_record == 1)
        break
    end
    m = m + 100 %snaps(i)
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

