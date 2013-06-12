%Plot all output files from simulation
close all
clear all

set(0,'DefaultAxesFontName', 'Times')

%Select diffusive solver (0) or full DNS (1)
CFD = 1;

%Turn on/off video and picture output
savevid = 0;
savepic = 0;

%Include line plots from MD data
% 0 = off
% 1 = on
plot_profiles = 1;


%Plot contours


%Include plot using CV data
% 0 = off
% 1 = 2D plots
% 2 = 3D plots
plot_CV = 0;
cres = 20;

%Analytical Solution parameters
u_0 = 1; t_0 = 160;
spectral_res = 6;
viscosity = 2.14; density = 0.8;
Re = density*u_0*1/viscosity;
analy_points = 20; % Number of spectral points

%Find results files
%resultfile_dir = './../results/';
%resultfile_dir = '/home/es205/results/MD_continuum_results/results/CPL_testing/130430_NCER_mdws/';
resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/NCER_wall_bump/basecase_no_bump/';
%resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/flekkoy/Flekkoy_meu10/';
%resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/flekkoy/Inc_specular_walls_large/';
%resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/flekkoy/50CFDMDratio/';
%resultfile_dir = '/home/djt06/Documents/Academia/PhD/Code/Development/branch/coupler_dCSE/src_code/';
resultfile_dir_MD = strcat(resultfile_dir,'md_data/results/');
resultfile_dir_CFD = strcat(resultfile_dir,'couette_data/');
resultfile_dir_CPL = strcat(resultfile_dir,'results/');

%Read MD Header
resultfile_dir = resultfile_dir_MD;
read_header

% = = = Read MD header = = =
%Check output flags and read data accordingly
if (velocity_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
elseif ((velocity_outflag > 0) && (velocity_outflag < 4) )
    read_vslice
elseif (mass_outflag == 4)
    mass_bins = read_mbins('mbins',resultfile_dir_MD);
elseif ((mass_outflag > 0) && (mass_outflag < 4) )
    read_mslice
end
%%MD domain set-up
Domain_setup
MD_domain = globaldomain;
MD_liquiddomain = liquiddomain;
MD_gnbins = gnbins;
MD_binsize = binsize;
clear liquiddomain gnbins binsize


if (CFD == 0)
    read_continuum_header
elseif(CFD == 1)
    %---Get CFD grid size ----
    [ngx, ngy, ngz, Lx, Ly, Lz, dx, dy, dz] = read_report(strcat(resultfile_dir_CFD,'report'));
end
% = = = Read CFD data = = =
if (CFD == 0)
    if (continuum_vflag == 3)
        read_continuum_vbins
    else
        read_continuum_vslice
    end
elseif(CFD == 1)
    [u,v,w] = Read_DNS_slice('grid.data',resultfile_dir_CFD,ngx-2,ngy-1,ngz-2,Lx,Ly,Lz,3,true);
    continuum_velslice = u;
end
%Continuum_Domain_setup
read_grid(strcat(resultfile_dir_CFD,'grid.data'),[1 1 1])
CFD_domain = [Lx,Ly,Lz];
CFD_binsize= [dx,dy,dz];
CFD_gnbins = [ngx,ngy,ngz]-1;
clear Lx Ly Lz dx dy dz ngx ngy ngz

%Read Coupled data
resultfile_dir = resultfile_dir_CPL;
CPL_read_header
CPL_olap_nbins = [icmax_olap-icmin_olap+1, jcmax_olap-jcmin_olap+1,kcmax_olap-kcmin_olap+1];
CPL_olap_domain= CFD_binsize.*CPL_olap_nbins;

%Calculate coupler properties
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
    %Get cell ratio
    MD_cells_per_CFD = CFD_binsize./MD_binsize;
    coupleddomain = CFD_domain + MD_liquiddomain - CPL_olap_domain;
    %Check domain size agrees
    domainratio  = CFD_domain./MD_domain;
    CFD_nbinswallbot= round(nbinswallbot./MD_cells_per_CFD);
end

%Get orthogonal direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    %[null,ixyz]=max(domainratio);
    ixyz = 2;
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
end

xaxis_MD  = linspace(0.5*CFD_binsize(ixyz)-tethdistbot(2),MD_liquiddomain(ixyz)-0.5*CFD_binsize(ixyz), ...
                    (ceil(gnbinsliquid(ixyz)/MD_cells_per_CFD(ixyz)))) ...
                    /coupleddomain(ixyz);
xaxis_MD2  = linspace(0,MD_liquiddomain(ixyz),gnbinsliquid(ixyz))/coupleddomain(ixyz);
xaxis_CFD = linspace(MD_liquiddomain(ixyz)-CPL_olap_domain(ixyz)-0.5*CFD_binsize(ixyz), ...
                    coupleddomain(ixyz)+0.5*CFD_binsize(ixyz),CFD_gnbins(ixyz)+2)/coupleddomain(ixyz);
xaxis_analy = linspace(0.0,1.0,analy_points); %As the liquid solid interaction is still 1, the domain is moved in slightly due to molecular sticking at the wall

%setup figures
scrsz = get(0,'ScreenSize');

if (plot_CV == 1)
    %Old multiple figure 
    %fig3 = figure('Position',[0 0 scrsz(3)/6 scrsz(4)/3]);
    %fig4 = figure('Position',[scrsz(3)/6 0 scrsz(3)/6 scrsz(4)/3]);
    %fig5 = figure('Position',[scrsz(3)/3 0 scrsz(3)/6 scrsz(4)/3]);

    %New subplot array
    close all
    fig1 = figure('Position',[1 1 scrsz(3)/2 scrsz(4)]);

else
    fig1 = figure('Position',[1 scrsz(4) scrsz(3)/4 scrsz(4)/2]);
    fig2 = figure('Position',[scrsz(3)/4 scrsz(4) scrsz(3)/4 scrsz(4)/2]);
    gax = axes('Position', [0.1, 0.1, 0.8, 0.8]);
    fig3 = figure('Position',[0 0 scrsz(3)/6 scrsz(4)/3]);
    fig4 = figure('Position',[scrsz(3)/6 0 scrsz(3)/6 scrsz(4)/3]);
    fig5 = figure('Position',[scrsz(3)/3 0 scrsz(3)/6 scrsz(4)/3]);
end

%Plot overlap domain schematic
%plot_domain
%read_grid(strcat(resultfile_dir_CFD,'grid.data'),[1 1 1],'plot')

%hax = axes('Position', [.18, .48, .23, .35]);

%Check if stress is avilable to plot
files = dir(resultfile_dir_MD);
plot_pVA = 0;
for i=1:size(files,1)
    if (strmatch('pVA',files(i).name) == 1)
        plot_pVA=1;
    end
end

%If MD stress is plotted, ensure CFD stress is available
if (plot_pVA == 1 && exist('stress','var') == 0) 
    for m =1:size(u,2)
        stress(1,2,:,m) = (1/Re) * cdiff(u(:,m),CFD_binsize(ixyz));
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
    if (plot_CV == 1)
        vidObj3 = VideoWriter('CV.avi');
        vidObj3.FrameRate = 10;
        open(vidObj3);
    end
end

%Contour plot axis
xaxis = linspace(0,globaldomain(1),globalncells(1));
yaxis = linspace(0,globaldomain(2),globalncells(2));
[X,Y] = meshgrid(xaxis,yaxis);

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
    filename = strcat(resultfile_dir_MD,'/mbins');
    mass_bins = read_mbins(filename,resultfile_dir_MD,m); 
    filename = strcat(resultfile_dir_MD,'/vbins');
    vel_bins = read_vbins(filename,resultfile_dir_MD,m); 
%     filename = strcat(resultfile_dir_MD,'/Tbins');
%     T_bins = read_Tbins(filename,resultfile_dir_MD,m); 

    if (plot_profiles == 1)

        %Calculate spanwise direction
        if (velocity_outflag == 4)
            %Find location of smallest dimension and take slice through this
            mass_slice = squeeze(sum(sum(mass_bins(:,:,:) ,jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
            vel_slice =  squeeze(sum(sum(vel_bins(:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
        end

        %Average multiple MD cells into one
        ave_mass_slice = interp1(mass_slice(:),1:MD_cells_per_CFD(ixyz):size(mass_slice(:)));

        %Average velocity per molecule
        for n =1:3
            ave_vel_slice(:,n) = interp1(vel_slice(:,n),1:MD_cells_per_CFD(ixyz):size(vel_slice(:,n)));
            ave_vel_slice(:,n) = ave_vel_slice(:,n)./ave_mass_slice';
            vel_slice(:,n) = vel_slice(:,n)./mass_slice';
        end

        %Plot density fluctuations in domain
        %set(0,'currentfigure',fig1)
        %plot(xaxis_MD,ave_mass_slice/(Nmass_ave*prod(CFD_binsize)),'x','Color',[.5 .5 .5],'MarkerSize',8,'LineWidth',5);
        %axis([-0.1 1.1 0.7 1.0 ])

        %Time -0.5 to record value at half time interval
        %if (i  == 1)
        %    t = 0;
        %else
            t =(intial_offset+m-0.5)*delta_t*Nmass_ave*tplot;
        %end
        disp(['step number ', num2str(m*Nmass_ave*tplot)])
        % =================================
        %    Plot MD region on full axis
        % =================================
        set(0,'currentfigure',fig2)
        set(fig2,'CurrentAxes',gax)

        %Plot anayltical solution
        analy = couette_analytical_fn(t,Re,[1,0],coupleddomain(ixyz),analy_points,'top');
        plot(xaxis_analy,analy,'k','LineWidth',5);
        hold on

        %Plot CFD velocity profile
        plot(xaxis_CFD,continuum_velslice(:,m),'s','Color',[.5 .5 .5],'MarkerSize',12,'LineWidth',3);

        %plot molecular velocity profile
        plot(xaxis_MD(:),ave_vel_slice(2:end,1),'x','LineWidth',3,'Color',[.5 .5 .5],'MarkerSize',14);
        %plot(xaxis_MD2,vel_slice(5:end,1),'^-')


        %Make plot look nice
        legend ('CFD','MD','location','NorthWest'); legend('boxoff')
        set(gca,'FontSize',16)
        axis([-0.1 1.1 -0.1 1.1]);
        xlabel('y/H'); ylabel('U_x/U')
        plottitle=num2str(t,'%10.6f');
        %title(strcat('Plot after  ',plottitle,' time units'));
        hold off

        % =================================
        %    Plot Insert
        % =================================
        %set(fig2,'CurrentAxes',hax)
        %set(fig2,'CurrentAxes',hax)
        %plot(abs(ave_vel_slice(:,2)));
        %axis([0 12 -0.01 0.2])
        %
        %     plot(linspace(MD_domain(ixyz)-overlap(ixyz)+0.5*CFD_binsize(ixyz),coupleddomain(ixyz)-0.5*CFD_binsize(ixyz),ngy-1)/coupleddomain(ixyz),continuum_velslice(:,m)/u_0,'s','Color',[.5 .5 .5],'MarkerSize',5,'LineWidth',2);
        %     hold on
        % 	plot(linspace(-0.5*CFD_binsize(ixyz),MD_domain(ixyz)-0.5*CFD_binsize(ixyz),gnbins(ixyz)/MD_cells_per_CFD)/coupleddomain(ixyz),ave_vel_slice(:,1),'x','LineWidth',5,'Color',[.2 .2 .2],'MarkerSize',5);
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
            %if (mod(m,Nvel_records/10) == 0)
                savefig(strcat('Velocity_',num2str(m)),'png')
            %end
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

             %Check stress and velocity evolve at same pace
             if (Nvel_ave ~= Nstress_ave)
                disp(['Stress is averaged at ', num2str(Nstress_ave), ...
                       ' step intervals but velocity at ',num2str(Nvel_ave), ' step intervals.' ...
                        ,'Attempting to compensate' ])
                ratio_stress2vel = Nvel_ave/Nstress_ave;
                rec = m * ratio_stress2vel;
             else
                rec = m;
             end

             cumulative_PVA = 0;
             filename = strcat(resultfile_dir_MD,'/pVA');
             for ave = ceil(-floor(t_ave/2)):floor(t_ave/2)
                 [pressure_VA]=read_pVA(rec+ave,resultfile_dir_MD,filename);
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
                 P_slice =  squeeze(sum(sum(pressure_VA(:,:,:,:,:),jxyz),kxyz)/(MD_gnbins(jxyz)*MD_gnbins(kxyz)));
             end

            %Average velocity per molecule
             for nxyz =1:3
             for mxyz =1:3
                 ave_P_slice(:,nxyz,mxyz) = interp1(P_slice(:,nxyz,mxyz),1:MD_cells_per_CFD(ixyz):size(P_slice,1));
                 %ave_P_slice(:,nxyz,mxyz) = ave_P_slice(:,nxyz,mxyz)./ave_mass_slice';
                 %P_slice(:,nxyz,mxyz) = P_slice(:,nxyz,mxyz)./mass_slice';
             end
             end


             %Stress figure
             figure(fig1)

             %Plot anayltical solution
             clear analy
             analy = couette_analytical_stress_fn(t,Re,1,coupleddomain(ixyz),analy_points);
             plot(xaxis_analy,analy,'k','LineWidth',5);
             hold on

             %Plot CFD shear stress
             plot(xaxis_CFD(2:end-1),squeeze(stress(1,2,:,m)),'s','Color',[.5 .5 .5],'MarkerSize',12,'LineWidth',3)

             %Plot MD shear stress
             %plot(xaxis_MD,-squeeze(density*ave_vel_slice(2:end,1).*ave_vel_slice(2:end,2)),'o','LineWidth',3,'Color',[.2 .2 .2],'MarkerSize',10)
             plot(xaxis_MD,-squeeze(ave_P_slice(2:end,1,2)),'x','LineWidth',3,'Color',[.5 .5 .5],'MarkerSize',10)
             %plot(xaxis_MD2,-squeeze(P_slice(5:end,1,2)),'^-')

             %Make plot look nice
             legend ('CFD','MD','location','NorthWest'); legend('boxoff')
             set(gca,'FontSize',16)
             axis([-0.1 1.1 -0.01 0.1]);
             xlabel('y/H'); ylabel('P_{xy}')
             plottitle=num2str(t,'%10.6f');
             %title(strcat('Plot after  ',plottitle,' time units'));
             hold off

             drawnow

             %Store pictures and videos
             if (savepic == 1)
                 %if (mod(m,Nvel_records/10) == 0)
                     savefig(strcat('Stress_',num2str(m)),'png')
                 %end
             elseif (savevid == 1)
                 currFrame = getframe(gcf);
                 writeVideo(vidObj2,currFrame);
             end

            figure(fig5)
            %Direct Stress
%             subplot(1,2,1); 
%             pVA = -squeeze(pressure_VA(:,:,:,1,1) + pressure_VA(:,:,:,2,2)+pressure_VA(:,:,:,3,3))/3;
%             contourf(X,Y,squeeze(mean(pVA(:,:,:),3))',cres,'LineStyle','None');
%             title('Direct Pressure','FontSize',16);  %colorbar('location','southoutside')
%             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
%             set(gca,'FontName','Times'); box on; %caxis([0.0 2.0])
%             set(gca,'FontSize',16);  colormap('hot')

            %Shear Stress
            %subplot(1,2,2); 
            contourf(X,Y,squeeze(mean(pressure_VA(:,:,:,1,2),3))',cres,'LineStyle','None');
            title('Shear Pressure','FontSize',16); colorbar('location','southoutside')
            xlabel('x','FontSize',16); ylabel('y','FontSize',16)
            set(gca,'FontName','Times'); box on; %caxis([0.0 2.0])
            set(gca,'FontSize',16);  colormap('hot')

        end

        %Mass plots
        figure(fig3)
        contourf(X,Y,squeeze(mean(mass_bins(:,:,:),3)/(prod(MD_binsize)*Nmass_ave))',cres,'LineStyle','None')
        axis equal; axis 'tight';  colorbar; title('velocity','FontSize',16)
        xlabel('x','FontSize',16); ylabel('y','FontSize',16)
        set(gca,'FontName','Times'); box on; caxis([0.7 0.85])
        set(gca,'FontSize',16);  colormap('hot')

        %Velocity plots
        figure(fig4)
        vmagnitude = sqrt(vel_bins(:,:,:,1).^2+vel_bins(:,:,:,2).^2+vel_bins(:,:,:,3).^2);
        contourf(X,Y,squeeze(mean(vmagnitude(:,:,:),3)./mean(mass_bins(:,:,:),3))',cres,'LineStyle','None')
        axis equal; axis 'tight';  colorbar; title('velocity','FontSize',16)
        xlabel('x','FontSize',16); ylabel('y','FontSize',16)
        set(gca,'FontName','Times'); box on; 
        set(gca,'FontSize',16);  colormap('hot')

        U = squeeze(mean(vel_bins(:,:,:,1),3))';
        V = squeeze(mean(vel_bins(:,:,:,2),3))';

        figure(fig4); hold on
        h = quiver(X,Y,U,V); 
        set(h,'Color','w'); hold off

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
    if (m > Nvel_records-ceil(t_ave/2)-2)
        m = Nvel_records-ceil(t_ave/2)-1
        last_record = 1;
    end


    %CV fluxes

    if (plot_CV ~= 0)

        filename = strcat(resultfile_dir_MD,'/pVA');
        [pressure_VA]=read_pVA(m,resultfile_dir_MD,filename);

        [mass_flux,mass_snapshot] = read_mflux('./mflux','./msnap',resultfile_dir_MD,m);
        [velocity_snapshot,velocity_flux,   ...
            pressure_surface,F_ext] = read_vflux(m,resultfile_dir_MD,MD_gnbins,nd);
        %Average just top surfaces
        Pk_CV = (velocity_flux(:,:,:,1,1) + velocity_flux(:,:,:,2,2)+velocity_flux(:,:,:,3,3))/3;
        Pc_CV = -squeeze(pressure_surface(:,:,:,1,1) + pressure_surface(:,:,:,2,2)+pressure_surface(:,:,:,3,3))/3;
        %Average all 6 surfaces
%         Pk_CV = (  velocity_flux(:,:,:,1,1) + velocity_flux(:,:,:,2,2)+velocity_flux(:,:,:,3,3) ...
%                  + velocity_flux(:,:,:,1,4) + velocity_flux(:,:,:,2,5)+velocity_flux(:,:,:,3,6))/6;
%         Pc_CV = -squeeze(  pressure_surface(:,:,:,1,1) + pressure_surface(:,:,:,2,2)+pressure_surface(:,:,:,3,3) ...
%                          + pressure_surface(:,:,:,1,4) + pressure_surface(:,:,:,2,5)+pressure_surface(:,:,:,3,6))/6;
% 

        pVA = -squeeze(pressure_VA(:,:,:,1,1) + pressure_VA(:,:,:,2,2)+pressure_VA(:,:,:,3,3))/3;
        totalflux =(velocity_flux(:,:,:,:,1)+velocity_flux(:,:,:,:,4))/MD_binsize(1) ...
                  +(velocity_flux(:,:,:,:,2)+velocity_flux(:,:,:,:,5))/MD_binsize(2) ...
                  +(velocity_flux(:,:,:,:,3)+velocity_flux(:,:,:,:,6))/MD_binsize(3);

        cum_Pk_CV = cum_Pk_CV + Pk_CV;
        cum_count = cum_count + 1;

        % ---------- 2 D Plots --------------------

        %Contour Plots
        if (plot_CV == 1)
            figure(fig1)

            rhouu = zeros(size(velocity_flux));
            rhouu(:,:,:,1,1) = (velocity_flux(:,:,:,1,1)-cum_Pk_CV/cum_count);
            rhouu(:,:,:,2,2) = (velocity_flux(:,:,:,2,2)-cum_Pk_CV/cum_count);
            rhouu(:,:,:,3,3) = (velocity_flux(:,:,:,3,3)-cum_Pk_CV/cum_count);

             xaxis   = linspace(0,globaldomain(1),globalncells(1));
             yaxis   = linspace(0,globaldomain(2),globalncells(2));
             [X,Y]   = meshgrid(xaxis,yaxis);

             %figure(fig1)
            subplot(2,8,1:2)
            contourf(X,Y,mean(mass_bins,3)'/(prod(MD_binsize)*Nmass_ave),cres,'LineStyle','None')
             axis equal; axis 'tight';  colorbar; title('density','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; caxis([0.4 0.85])
             set(gca,'FontSize',16);  colormap('hot')


             %figure(fig2)
            subplot(2,8,3:4)
            vmagnitude = sqrt(vel_bins(:,:,:,1).^2+vel_bins(:,:,:,2).^2+vel_bins(:,:,:,3).^2);
            contourf(X,Y,squeeze(mean(vmagnitude(:,:,:),3)./mean(mass_bins(:,:,:),3))',cres,'LineStyle','None')
             axis equal; axis 'tight';  colorbar; title('velocity','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; 
             set(gca,'FontSize',16);  colormap('hot')
            U = squeeze(mean(vel_bins(:,:,:,1),3)./mean(mass_bins(:,:,:),3))';
            V = squeeze(mean(vel_bins(:,:,:,2),3)./mean(mass_bins(:,:,:),3))';
            hold on
            h = quiver(X,Y,U,V); 
            set(h,'Color','w'); hold off
%             hold on
%             [sx,sy] = meshgrid(0,linspace(tethdistbot(2)+0.5,globaldomain(2),20));
%             h = streamline(X,Y,U,V,sx,sy);
%             set(h,'Color','w','LineWidth',0.01);  hold off
%              [sx,sy] = meshgrid(26,linspace(tethdistbot(2)+0.5,tethdistbot(2)+10,10));
%              h = streamline(X,Y,U,V,sx,sy);
%              set(h,'Color','w','LineWidth',0.01);  hold off
             %vf_slice_hist(i,:,:,:,:) = velocity_flux(:,6,:,:,:);

          
            %figure(fig3)
            subplot(2,8,5:6)
            contourf(X,Y,(mean(T_bins,3)./(3*mean(mass_bins,3)))',cres,'LineStyle','None')
             axis equal; axis 'tight';  colorbar; title('Temperature','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; caxis([0.0 1.1])
             set(gca,'FontSize',16);  colormap('hot')

            %figure(fig4)
%             subplot(2,8,5:6)
%              contourf(X,Y,squeeze(mean(pVA(:,:,:),3))',cres,'LineStyle','None');
%              axis equal; axis 'tight';  colorbar; title('Pressure','FontSize',16)
%              xlabel('x','FontSize',16); ylabel('y','FontSize',16)
%              set(gca,'FontName','Times'); box on; caxis([0.0 2.0])
%              set(gca,'FontSize',16);  colormap('hot')

            subplot(2,8,7:8)
             %contourf(X,Y,-squeeze(mean(cum_Pk_CV(:,:,:)/cum_count,3))',cres,'LineStyle','None');
             %contourf(X,Y,squeeze(mean(Pk_CV(:,:,:),3)+mean(Pc_CV(:,:,:),3))',cres,'LineStyle','None');  
             contourf(X,Y,squeeze(mean(pVA(:,:,:),3))'-squeeze(mean(Pk_CV(:,:,:),3)+mean(Pc_CV(:,:,:),3))',cres,'LineStyle','None');             axis equal; axis 'tight';  colorbar; title('Pressure','FontSize',16)

           axis equal; axis 'tight';  colorbar; title('Pressure','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; %caxis([0.0 2.0])
             set(gca,'FontSize',16);  colormap('hot')


            %figure(fig5)
            subplot(2,8,9:10)
             contourf(X,Y,-squeeze(mean(rhouu(:,:,:,1,1),3))',cres,'LineStyle','None');
             axis equal; axis 'tight';  colorbar; title('\int_S_x ( \rho u_x u_x )  dS_x','FontSize',16)

             %contourf(X,Y,-squeeze(mean(velocity_flux(:,:,:,1,1),3))',cres,'LineStyle','None');
             %axis equal; axis 'tight';  colorbar; title('\int_S_x ( \rho u_x u_x + P ) dS_x','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; caxis([-0.2 0.25])
             set(gca,'FontSize',16);  colormap('hot')

            %figure(fig6)
            subplot(2,8,11:12)
             contourf(X,Y,-squeeze(mean(rhouu(:,:,:,2,2),3))',cres,'LineStyle','None');
             axis equal; axis 'tight';  colorbar; title('\int_S_y ( \rho u_y u_y )  dS_y','FontSize',16)

             %contourf(X,Y,-squeeze(mean(velocity_flux(:,:,:,2,2),3))',cres,'LineStyle','None');
             %axis equal; axis 'tight';  colorbar; title('\int_S_y ( \rho u_y u_y + P )  dS_y','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; caxis([-0.2 0.25])
             set(gca,'FontSize',16);  colormap('hot')


            %figure(fig7)
            subplot(2,8,13:14)
            contourf(X,Y,mean(sum(totalflux(:,:,:,:),4),3)',cres,'LineStyle','None')
             axis equal; axis 'tight';  colorbar; title('\int_S ( \rho u u + P ) \cdot dS','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; %caxis([0.0 1.1])
             set(gca,'FontSize',16);  colormap('hot')

            %figure(fig8);
            subplot(2,8,15:16)
            %plot(0.333*mean(mean(T_bins,1),3)./mean(mean(mass_bins,1),3),1:size(T_bins,2))
            contourf(X,Y,mean(mass_flux(:,:,:,1,1),3)',cres,'LineStyle','None')
             axis equal; axis 'tight';  colorbar; title('\int_S \rho u  \cdot dS','FontSize',16)
             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
             set(gca,'FontName','Times'); box on; %caxis([0.0 1.1])
             set(gca,'FontSize',16);  colormap('hot')


            if (savevid == 1)
                currFrame = getframe(gcf);
                writeVideo(vidObj3,currFrame);
            end

%             CV_contourf_plot(globaldomain,globalncells,-rhouu,1,1,fig1);
%             CV_contourf_plot(globaldomain,globalncells,rhouu,2,2,fig2);
%             CV_contourf_plot(globaldomain,globalncells,-rhouu,3,3,fig3);
% 
%             %Streamline plots
%             %Apparently streamlines require seeding everywhere!
%             xaxis = linspace(0,globaldomain(1),globalncells(1));
%             yaxis = linspace(0,globaldomain(2),globalncells(2));
%             [X,Y] = meshgrid(xaxis,yaxis);
%             [sx,sy] = meshgrid(xaxis,yaxis);
% 
%             U = -squeeze(mean(rhouu(:,:,:,1,1),3))';
%             V = squeeze(mean(rhouu(:,:,:,2,2),3))';
% 
%             figure(fig1); hold on
%             h = quiver(X,Y,U,V); 
%             %h = streamline(X,Y,U,V,sx(1:end),sy(1:end));
%             set(h,'Color','k'); hold off
% 
%             figure(fig4)
%             xaxis = linspace(0,globaldomain(1),globalncells(1));
%             yaxis = linspace(0,globaldomain(2),globalncells(2));
%             [X,Y] = meshgrid(xaxis,yaxis);
%             contourf(X,Y,squeeze(mean(totalflux(:,:,:,1),3))',cres,'LineStyle','None');
%             axis equal; axis 'tight';  colorbar; title('velocity','FontSize',16)
%             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
%             set(gca,'FontName','Times'); box on; %caxis([-0.1 1.0])
%             set(gca,'FontSize',16);  colormap('hot')
% 
%             %Streamline plots
%             %Apparently streamlines require seeding everywhere!
%             xaxis = linspace(0,globaldomain(1),globalncells(1));
%             yaxis = linspace(0,globaldomain(2),globalncells(2));
%             [X,Y] = meshgrid(xaxis,yaxis);
%             [sx,sy] = meshgrid(xaxis,yaxis);
% 
%             U = squeeze(mean(vel_bins(:,:,:,1),3))';
%             V = squeeze(mean(vel_bins(:,:,:,2),3))';
% 
%             figure(fig4); hold on
%             h = quiver(X,Y,U,V); 
%             set(h,'Color','w'); hold off
% 
%             figure(fig5)
%             xaxis = linspace(0,globaldomain(1),globalncells(1));
%             yaxis = linspace(0,globaldomain(2),globalncells(2));
%             [X,Y] = meshgrid(xaxis,yaxis);
%             contourf(X,Y,squeeze(mean(Pk_CV(:,:,:),3))',cres,'LineStyle','None');
%             axis equal; axis 'tight';  colorbar; title('kinetic pressure','FontSize',16)
%             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
%             set(gca,'FontName','Times'); box on; %caxis([-0.1 1.0])
%             set(gca,'FontSize',16);  colormap('hot')

%             figure(fig1); hold on
%             h = streamline(X,Y,U,V,sx(1:2:end),sy(1:2:end));
%             set(h,'Color','w'); hold off

    elseif (plot_CV == 2)
        % ---------- 3 D Plots --------------------

        %Isosurface plots
        xaxis = linspace(0,globaldomain(1),globalncells(1));
        yaxis = linspace(0,globaldomain(2),globalncells(2));
        zaxis = linspace(0,globaldomain(3),globalncells(3));
        [X,Y,Z] = meshgrid(yaxis,xaxis,zaxis);

        %Fixes a 3D plotting bug in matlab running on 2 screens
        figure(fig3);  clf;
        set(0,'DefaultFigureRenderer','OpenGL');
        %Red kinetic part
        [faces,vertices] = isosurface(X,Y,Z,Pk_CV, 0.3);
        p=patch('Faces',faces,'Vertices',vertices);
        isonormals(X,Y,Z,Pk_CV,p);
        set(p,'FaceColor','red','EdgeColor','none');
        %Blue configurational
        [faces,vertices] = isosurface(X,Y,Z,(squeeze(vel_bins(:,:,:,1))./mass_bins(:,:,:)), 0.5);
        p=patch('Faces',faces,'Vertices',vertices);
        isonormals(X,Y,Z,Pc_CV,p);
        set(p,'FaceColor','blue','EdgeColor','none');
        daspect([1 1 1]); box on;
        view([72 -68]); axis tight
        camlight 
        lighting gouraud

    end
 
%          U = -squeeze(vel_bins(:,:,:,1));
%          V = -squeeze(vel_bins(:,:,:,2));
%          W = -squeeze(vel_bins(:,:,:,3));
% 
%         [sx,sy,sz] = meshgrid(16,xaxis,zaxis);

% 
%         figure(fig3); hold on
%         h = streamline(X,Y,Z,U,V,W,sx(1:end),sy(1:end),sz(1:end));
%         set(h,'Color','k'); hold off
        

%          CV_contourf_plot(globaldomain,globalncells,velocity_flux,1,1,fig3)
%          if (savevid == 1)
%             currFrame = getframe(gcf);
%             writeVideo(vidObj3,currFrame);
%          end
%          CV_contourf_plot(globaldomain,globalncells,velocity_flux,2,2,fig4)
%          %CV_contourf_plot(globaldomain,globalncells,velocity_flux,2,1)
%          CV_contourf_plot(globaldomain,globalncells,velocity_flux,3,3,fig5)
% 
%         CV_contourf_plot(globaldomain,globalncells,pressure_surface,1,1)
%         CV_contourf_plot(globaldomain,globalncells,pressure_surface,1,2)
%         CV_contourf_plot(globaldomain,globalncells,pressure_surface,2,2)
%         CV_contourf_plot(globaldomain,globalncells,pressure_surface,2,1)
%         CV_contourf_plot(globaldomain,globalncells,pressure_surface,3,3)



%             contourf(X,Y,-squeeze(mean(temp,3))',cres,'LineStyle','None');

%             axis equal; axis 'tight';  colorbar; title([inputname(3),' Direction =',num2str(ixyz), ' Face =',num2str(face)],'FontSize',16)
%             xlabel('x','FontSize',16); ylabel('y','FontSize',16)
%             set(gca,'FontName','Times'); box on;
%             set(gca,'FontSize',16);  colormap('hot')
%             ierr = 'success';

%         figure
%         xaxis = linspace(0,globaldomain(1),globalncells(1));
%         yaxis = linspace(0,globaldomain(2),globalncells(2));
%         [X,Y] = meshgrid(xaxis,yaxis);
%         contourf(X,Y,-squeeze(mean(velocity_flux(:,:,:,1,1),3))',cres,'LineStyle','None');
%         axis equal; axis 'tight';  colorbar; title('[\rho u u ]_{xx} Momentum flux')
%         xlabel('x','FontSize',16); ylabel('y','FontSize',16)
%         set(gca,'FontSize',16);  colormap('hot')
%         set(gca,'FontName','Times'); box on;


%         figure
%         contourf(X,Y,-squeeze(mean(pressure_surface(:,:,:,1,2),3))',cres,'LineStyle','None');
%         axis equal; axis 'tight';  colorbar
 
    end







end

% Close the Video file.
if (savevid == 1)
    close(vidObj)
    if (plot_pVA == 1)
        close(vidObj2)
    end
    if (plot_CV == 1)
        close(vidObj3);
    end
 end

%movefile ./couette.avi ./../results

