clear all
close all

%Turn on/off video and picture output
savevid = 0;
savepic = 0;

%Plot type
% 1 = mass/velocity/temperature in bins
% 2 = CV analysis
plot_type = 2;

%Plot options
mass_splines = 0; %Add splines to mass results

%Slicetype
% 0 = off
% 1 = slice
% 2 = contourf
slicetype = 2;
cres = 100; %Contour resolution

%Shitfting cooefficients - may need some tuning to give correct colorbar
%placements
figup = 0.04;
clrbardown = 0.06;

%Nozzle mask
% Average only fluid regions inside nozzle
nozzle_mask = 20;


%Starting time for steady state averages
start = 1;

resultfile_dir_MD = '/home/es205/results/md_results/fortran/3D_code/parallel/results/converge_diverge/high_temporal_resolution/'
%resultfile_dir_MD = '/home/es205/codes/coupled/MD_dCSE/src_code/results/'
cd(resultfile_dir_MD)

%Read Header file
resultfile_dir = resultfile_dir_MD;
read_header
N_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));

%Range of timesteps
tstart = 1;
tskip = 1;
tend = N_records-3;


%PLOTS STUFF
if (plot_type == 1)

    %Setup figures
    scrsz = get(0,'ScreenSize');
    set(0,'DefaultFigureRenderer','OpenGL')
    if (slicetype ~= 0)
        fig1 = figure('Position',[1 scrsz(4) scrsz(3)/4 scrsz(4)/3]);
        fig2 = figure('Position',[scrsz(3)/4 scrsz(4) scrsz(3)/4 scrsz(4)/3]);
        fig3 = figure('Position',[1 1 scrsz(3)/4 scrsz(4)/3]);
    end

    if (savevid == 1)
        vidObjm = VideoWriter('mass.avi');
        vidObjv = VideoWriter('velocity.avi');
        vidObjT = VideoWriter('temperature.avi');
        vidObjm.FrameRate = 10;
        vidObjv.FrameRate = 10;
        vidObjT.FrameRate = 10;
        open(vidObjm); open(vidObjv); open(vidObjT);
    end

    %Mass
    minlt_t = zeros(gnbins(2),tend-start+1);
    mqrtr_t = zeros(gnbins(2),tend-start+1);
    mhalf_t = zeros(gnbins(2),tend-start+1);

    %Velocity
    uinlt_t = zeros(gnbins(2),tend-start+1);
    uqrtr_t = zeros(gnbins(2),tend-start+1);
    uhalf_t = zeros(gnbins(2),tend-start+1);

    %Temperature
    Tinlt_t = zeros(gnbins(2),tend-start+1);
    Tqrtr_t = zeros(gnbins(2),tend-start+1);
    Thalf_t = zeros(gnbins(2),tend-start+1);
    

    for i=tstart:tskip:tend
        i

        filename = strcat(resultfile_dir_MD,'/mbins');
        mass_bins = read_mbins(filename,resultfile_dir_MD,i);
        filename = strcat(resultfile_dir_MD,'/vbins');
        vel_bins = read_vbins(filename,resultfile_dir_MD,i);
        filename = strcat(resultfile_dir,'/Tbins');
        T_bins = read_Tbins(filename,resultfile_dir,i); 

        %Mass
        if (slicetype == 0)
            %Do nothing
        elseif (slicetype == 1)
            figure(fig1); clf
            h=slice(mass_bins(:,:,:)/(prod(binsize)*Nmass_ave),[],[],[10]);
            axis 'equal';  view([90,90]);
            set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
            caxis([0.6 1.2])
            colorbar('location','southoutside'); 
        elseif (slicetype == 2)   
            figure(fig1); clf
            xaxis = linspace(0,globaldomain(1),gnbins(1));
            yaxis = linspace(0,globaldomain(2),gnbins(2));
            [X,Y] = meshgrid(xaxis,yaxis);
            contourf(X,Y,(squeeze(mean(mass_bins(:,:,:,1),3))./(prod(binsize)*Nmass_ave))',cres,'LineStyle','None')
            axis equal; axis 'tight';  caxis([0.4 0.8]); 

            %Include colourbar and shift down slightly
            colorbar('delete'); c = colorbar('location','southoutside');
            posc=get(c,'Position'); posf=get(gca,'Position');
            set(gca,'position',[posf(1),posf(2)+figup,posf(3),posf(4)]);
            set(c,'Position',[posc(1),posc(2)-clrbardown,posc(3),posc(4)]);

            %Plot nozzzle bounding lines
            hold on
            a = 0.25*globaldomain(2)/2;
            b = tethdistbot(2) + 0.5*binsize(2);
            plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'k','LineWidth', 2)
            plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'k','LineWidth', 2)
            plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'w:','LineWidth', 2)
            plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'w:','LineWidth', 2)

            %labels etc
            xlabel('x','FontSize',16); ylabel('y','FontSize',16)
            set(gca,'FontSize',16);  colormap('hot')
            set(gca,'FontName','Times'); box on;

            if (savepic == 1)
                savefig('./nozzle_mcontour','png','eps')
            end
            if (savevid == 1)
                currFrame = getframe(gcf);
                writeVideo(vidObjm,currFrame);
            end

            %pause(0.1)
            drawnow

        end


        
        %Velocity
        if (slicetype == 0)
            %Do nothing
        elseif (slicetype == 1)
            figure(fig2);  clf
            h=slice(vel_bins(:,:,:,1)./mass_bins(:,:,:),[],[],[10]);
            axis 'equal'; view([90,90]);
            set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
            caxis([-0.2 0.51])
            colorbar('location','southoutside'); 
        elseif (slicetype == 2)  
            figure(fig2);  clf
            xaxis = linspace(0,globaldomain(1),gnbins(1));
            yaxis = linspace(0,globaldomain(2),gnbins(2));
            [X,Y] = meshgrid(xaxis,yaxis);
            U = (squeeze(mean(vel_bins(:,:,:,1),3))./squeeze(mean(mass_bins(:,:,:),3)))';
            V = (squeeze(mean(vel_bins(:,:,:,2),3))./squeeze(mean(mass_bins(:,:,:),3)))';
            contourf(X,Y,U,cres,'LineStyle','None')
            %hold on
            %qs = 4
            %quiver(X(1:qs:end,1:qs:end,1:qs:end), ...
            %       Y(1:qs:end,1:qs:end,1:qs:end), ...
            %       U(1:qs:end,1:qs:end,1:qs:end), ...
            %       V(1:qs:end,1:qs:end,1:qs:end))
            axis equal; axis 'tight';  caxis([0.0 1.8]); %ylim([8 42]);
            %Include colourbar and shift down slightly
            colorbar('delete'); c = colorbar('location','southoutside');
            posc=get(c,'Position'); posf=get(gca,'Position');
            set(gca,'position',[posf(1),posf(2)+figup,posf(3),posf(4)]);
            set(c,'Position',[posc(1),posc(2)-clrbardown,posc(3),posc(4)]);

            %Plot nozzzle bounding lines
            hold on
            a = 0.25*globaldomain(2)/2;
            b = tethdistbot(2) + 0.5*binsize(2);
            plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'k','LineWidth', 2)
            plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'k','LineWidth', 2)
            plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'w:','LineWidth', 2)
            plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'w:','LineWidth', 2)

            %labels etc
            xlabel('x','FontSize',16); ylabel('y','FontSize',16)
            set(gca,'FontSize',16);  colormap('hot')
            set(gca,'FontName','Times'); box on;

            if (savepic == 1)
                savefig('./nozzle_vcontour','png','eps')
            end

            if (savevid == 1)
                currFrame = getframe(gcf);
                writeVideo(vidObjv,currFrame);
            end

            %pause(0.1)
            drawnow

        end

        %Temperature
        if (slicetype == 0)
            %Do nothing
        elseif (slicetype == 1)
            figure(fig3);  clf
            h=slice(T_bins(:,:,:)./(3*mass_bins(:,:,:)),[],[],[10]);
            axis 'equal'; view([90,90]);
            set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
            caxis([-0.2 2.01])
            colorbar('location','southoutside'); 
        elseif (slicetype == 2)
            figure(fig3);  clf
            xaxis = linspace(0,globaldomain(1),gnbins(1));
            yaxis = linspace(0,globaldomain(2),gnbins(2));
            [X,Y] = meshgrid(xaxis,yaxis);
            contourf(X,Y,(squeeze(mean(T_bins(:,:,:,1),3))./squeeze(3*mean(mass_bins(:,:,:),3)))',cres,'LineStyle','None')
            axis equal; axis 'tight';  caxis([0.8 2.0]); 

            %Include colourbar and shift down slightly
            colorbar('delete'); c = colorbar('location','southoutside');
            posc=get(c,'Position'); posf=get(gca,'Position');
            set(gca,'position',[posf(1),posf(2)+figup,posf(3),posf(4)]);
            set(c,'Position',[posc(1),posc(2)-clrbardown,posc(3),posc(4)]);

            %Plot nozzzle bounding lines
            hold all
            a = 0.25*globaldomain(2)/2;
            b = tethdistbot(2) + 0.5*binsize(2);
            plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'k','LineWidth', 2)
            plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'k','LineWidth', 2)
            plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'w:','LineWidth', 2)
            plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'w:','LineWidth', 2)

            %labels etc
            xlabel('x','FontSize',16); ylabel('y','FontSize',16)
            set(gca,'FontSize',16);  colormap('hot')
            set(gca,'FontName','Times'); box on;

            if (savepic == 1)
                savefig('./nozzle_Tcontour','png','eps')
            end

            if (savevid == 1)
                currFrame = getframe(gcf);
                writeVideo(vidObjT,currFrame);
            end


        elseif (slicetype == 3)
            %Line plots
            T = mean(T_bins(:,2,:)./(3*mass_bins(:,2,:)),3);
            T(isnan(T)) = 0;
            plot(T,'rs')
            hold all
            plot(mean(mass_bins(:,2,10)/(prod(binsize)*Nmass_ave),3),'bx')
            hold off

            %pause(0.1)
            drawnow
        end

        %Save time averaged values
        if (i > start)

            %Mass
            minlt_t(:,i) = mean(mean(squeeze(mass_bins(  1  ,:,:)/(prod(binsize)*Nmass_ave)),2),3);
            mqrtr_t(:,i) = mean(mean(squeeze(mass_bins(end/4,:,:)/(prod(binsize)*Nmass_ave)),2),3);
            mhalf_t(:,i) = mean(mean(squeeze(mass_bins(end/2,:,:)/(prod(binsize)*Nmass_ave)),2),3);

            %Velocity
            uinlt_t(:,i) = mean(mean(squeeze(vel_bins(  1  ,:,:,1))./squeeze(mass_bins(  1  ,:,:)),2),3);
            uqrtr_t(:,i) = mean(mean(squeeze(vel_bins(end/4,:,:,1))./squeeze(mass_bins(end/4,:,:)),2),3);
            uhalf_t(:,i) = mean(mean(squeeze(vel_bins(end/2,:,:,1))./squeeze(mass_bins(end/2,:,:)),2),3);

            %Temperature
            Tinlt_t(:,i) = mean(mean(squeeze(T_bins(  1  ,:,:))./squeeze(3*mass_bins(  1  ,:,:)),2),3);
            Tqrtr_t(:,i) = mean(mean(squeeze(T_bins(end/4,:,:))./squeeze(3*mass_bins(end/4,:,:)),2),3);
            Thalf_t(:,i) = mean(mean(squeeze(T_bins(end/2,:,:))./squeeze(3*mass_bins(end/2,:,:)),2),3);


            %Centreline values
            centreline = round(size(vel_bins,2)/2);
            vcl(i) =   squeeze(mean(mean(squeeze(vel_bins(:,centreline,:,1,:)) ...
                     ./squeeze(mass_bins(:,centreline,:,:)),1),2));

            Tcl(i) =   squeeze(mean(mean(squeeze(T_bins(:,centreline,:,:)) ...
                     ./squeeze(3*mass_bins(:,centreline,:,:)),1),2));
        end

    end
   
    %Plot density variation along channel
    figure
    %plot(mean(mean(squeeze(mass_bins(:,:,:,i)/(prod(binsize)*Nmass_ave)),2),3),'x')
    
    % Plot velcity profile at 3 locations
    y = linspace(binsize(2)/2,globaldomain(2)-binsize(2)/2,gnbins(2));
    ycf = binsize(2)/2:binsize(2)/100:globaldomain(2)-binsize(2)/2;

    %Mass bins
    %filename = strcat(resultfile_dir_MD,'/mbins');
    %mass_bins = read_mbins(filename,resultfile_dir_MD);
    minlt = mean(minlt_t(:,start:end),2); %mean(mean(squeeze(mass_bins(  1  ,:,:,start:end)/(prod(binsize)*Nmass_ave)),2),3);
    mqrtr = mean(mqrtr_t(:,start:end),2); %mean(mean(squeeze(mass_bins(end/4,:,:,start:end)/(prod(binsize)*Nmass_ave)),2),3);
    mhalf = mean(mhalf_t(:,start:end),2); %mean(mean(squeeze(mass_bins(end/2,:,:,start:end)/(prod(binsize)*Nmass_ave)),2),3);
	hold all

    %Fit interplolated splines to density data
    if (mass_splines == 1)
        ft = fittype('splineinterp'); 
        ps = 1;

        %Halfway
        cf = fit(y(:),mhalf(:),ft);    %Get spline data
        mhalf_cf = feval(cf,ycf);       %Evaluate and convert to array
        plot(mhalf_cf(ps:end-ps),ycf(ps:end-ps),'color', [0.7,0.7,0.7],'LineWidth',2) %Plot spline

        %Quarter
        cf = fit(y(:),mqrtr(:),ft);    %Get spline data
        mqrtr_cf = feval(cf,ycf);       %Evaluate and convert to array
        plot(mqrtr_cf(ps:end-ps),ycf(ps:end-ps),'--','color', [0.3,0.3,0.3],'LineWidth',2) %Plot spline

        %Initial
        cf = fit(y(:),minlt(:),ft);    %Get spline data
        minlt_cf = feval(cf,ycf);       %Evaluate and convert to array
        plot(minlt_cf(ps:end-ps),ycf(ps:end-ps),'-','color', [0.0,0.0,0.0],'LineWidth',2) %Plot spline

    end

	plot(mhalf,y,'-x','color', [0.7,0.7,0.7],'MarkerSize',8) %Plot data
    plot(mqrtr,y,'-s','color', [0.3,0.3,0.3],'MarkerSize',10) %Plot data
    plot(minlt,y,'-o','color', [0.0,0.0,0.0],'MarkerSize',6) %Plot data

    %labels etc
    axis([0.4 1.55  -2.0 58.0 ])
    xlabel('Density','FontSize',16); ylabel('y','FontSize',16)
    set(gca,'FontSize',16); 
    set(gca,'FontName','Times'); box on;
    if (savepic == 1)
       savefig('./conv_divg_mprofiles','png','eps')
    end

    %Vbins
    figure
    %filename = strcat(resultfile_dir_MD,'/vbins');
    %vel_bins = read_vbins(filename,resultfile_dir_MD);
    uinlt = mean(uinlt_t(:,start:end),2); %mean(mean(squeeze(vel_bins(  1  ,:,:,1,start:end))./squeeze(mass_bins(  1  ,:,:,start:end)),2),3);
    uqrtr = mean(uqrtr_t(:,start:end),2); %mean(mean(squeeze(vel_bins(end/4,:,:,1,start:end))./squeeze(mass_bins(end/4,:,:,start:end)),2),3);
    uhalf = mean(uhalf_t(:,start:end),2); %mean(mean(squeeze(vel_bins(end/2,:,:,1,start:end))./squeeze(mass_bins(end/2,:,:,start:end)),2),3);

    plot(uhalf,y,'color', [0.7,0.7,0.7],'LineWidth',2)
    hold all
    plot(uqrtr,y,'--','color', [0.3,0.3,0.3],'LineWidth',2)
    plot(uinlt,y,'-','color', [0.0,0.0,0.0],'LineWidth',2)

    %labels etc
    axis([-0.1 1.8  -2.0 58.0 ])
    xlabel('Velocity','FontSize',16); ylabel('y','FontSize',16)
    set(gca,'FontSize',16); 
    set(gca,'FontName','Times'); box on;
    if (savepic == 1)
        savefig('./conv_divg_vprofiles','png','eps')
    end

    %Obtain time evolution of centreline velocity
    figure
    %vcl = squeeze(mean(mean(squeeze(vel_bins(:,15,:,1,:))./squeeze(mass_bins(:,15,:,:)),1),2));
    taxis = (1:size(vcl,2))*delta_t*Nvel_ave*tplot;
    plot(taxis,vcl,'color',[0.6,0.6,0.6])
    hold on
    %svcl = smooth(vcl,0.5);
    %plot(taxis(10:end-10),svcl(10:end-10),'color',[0.0,0.0,0.0],'LineWidth',2)
    %axis([0 5500 0.02 0.05])
    xlabel('Time','FontSize',16); ylabel('Velocity','FontSize',16)
    set(gca,'FontSize',16); 
    set(gca,'FontName','Times'); box on;
    if (savepic == 1)
        savefig('./vcl_tevo','png','eps')
    end

    %axis([0 250 0 0.05])

    clear vel_bins

    %Tbins
    %filename = strcat(resultfile_dir,'/Tbins');
    %T_bins = read_Tbins(filename,resultfile_dir); 
    Tinlt = mean(Tinlt_t(:,start:end),2); %mean(mean(squeeze(T_bins(  1  ,:,:,start:end))./squeeze(3*mass_bins(  1  ,:,:,start:end)),2),3);
    Tqrtr = mean(Tqrtr_t(:,start:end),2); %mean(mean(squeeze(T_bins(end/4,:,:,start:end))./squeeze(3*mass_bins(end/4,:,:,start:end)),2),3);
    Thalf = mean(Thalf_t(:,start:end),2); %mean(mean(squeeze(T_bins(end/2,:,:,start:end))./squeeze(3*mass_bins(end/2,:,:,start:end)),2),3);
    
    figure
    plot(Thalf,y,'color', [0.7,0.7,0.7],'LineWidth',2)
    hold all
    plot(Tqrtr,y,'--','color', [0.3,0.3,0.3],'LineWidth',2)
    plot(Tinlt,y,'-','color', [0.0,0.0,0.0],'LineWidth',2)

    %labels etc
    axis([0.85 1.75  -2.0 58.0 ])
    xlabel('Temperature','FontSize',16); ylabel('y','FontSize',16)
    set(gca,'FontSize',16); 
    set(gca,'FontName','Times'); box on;
    if (savepic == 1)
        savefig('./conv_divg_Tprofiles','png','eps')
    end

    %Obtain time evolution of centreline Temperature
    figure
    %Tcl = squeeze(mean(mean(squeeze(T_bins(:,15,:,:))./squeeze(3*mass_bins(:,15,:,:)),1),2));
    taxis = (1:size(Tcl,2))*delta_t*Nvel_ave*tplot;
    plot(taxis,Tcl,'color',[0.6,0.6,0.6])
    hold on
    sTcl = smooth(Tcl,0.5);
    plot(taxis(10:end-10),sTcl(10:end-10),'color',[0.0,0.0,0.0],'LineWidth',2)
    %axis([0 5500 0.985 1.01])
    xlabel('Time','FontSize',16); ylabel('Temperature','FontSize',16)
    set(gca,'FontSize',16); 
    set(gca,'FontName','Times'); box on;
    if (savepic == 1)
        savefig('./Tcl_tevo','png','eps')
    end

% 	figure
% 	plot(y,mqrtr.*uqrtr)
% 	hold all
% 	plot(y,mhalf.*uhalf)
% 	plot(y,minlt.*uinlt)
%     
%     trapz(minlt.*uinlt)
%     trapz(mqrtr.*uqrtr)
%     trapz(mhalf.*uhalf)
    
elseif(plot_type == 2)

    %Establish and Plot domain set-up
    Domain_setup
    %setup stress direction and face
    ixyz = 1;
    jxyz = 2;
    Nvflux_records = Nsteps / (Nvflux_ave);
    n = 1; tol = 10000*eps;
    external_force_flag = 1
    ave_count = 0;
    ave_dv = zeros(gnbins(1),1);
    ave_mf = zeros(gnbins(1),1);
    ave_vf = zeros(gnbins(1),1);
    ave_vf1 = zeros(gnbins(1),1);

    if (nozzle_mask == 0)
        %Plot entire domain values
        nozzleCV = ones(nbins(1),nbins(2),nbins(3));
    else
        %Plot only values inside nozzle
        [nozzleCV,ytop,ybot,xstepup,xstepdown] = ...
                   nozzle_CV(   -globaldomain(1)/2,globaldomain(1)/2, ...
                                -globaldomain(2)/2,globaldomain(2)/2, ...
                                -globaldomain(3)/2,globaldomain(3)/2, ...
                                -globaldomain(1)/2,globaldomain(1)/2, ...
                                -globaldomain(2)/2+cellsidelength(2)+tethdistbot(2), ...
                                 globaldomain(2)/2-cellsidelength(2)+tethdisttop(2), ...
                                -globaldomain(3)/2,globaldomain(3)/2, ...
                                 gnbins(1),gnbins(2),gnbins(3),0.125,0);

        %Plot and check
        xaxis = linspace(0,globaldomain(1),gnbins(1));
        yaxis = linspace(0,globaldomain(2),gnbins(2));
        [X,Y] = meshgrid(xaxis,yaxis);
        contourf(X,Y,squeeze(mean(nozzleCV(:,:,:),3))',cres,'LineStyle','None')
        axis equal; axis 'tight'; colormap(gray)
        hold on

        %Plot steps
        scatter(xstepup*binsize(1),ytop(xstepup)*binsize(2),'rx')
        scatter(xstepdown*binsize(1),ytop(xstepdown)*binsize(2),'rx')
        scatter(xstepup*binsize(1),ybot(xstepup)*binsize(2),'rx')
        scatter(xstepdown*binsize(1),ybot(xstepdown)*binsize(2),'rx')

        %Plot lines
        a = 0.25*globaldomain(2)/2;
        b = tethdistbot(2) + 0.5*binsize(2);
        plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'r','LineWidth', 2)
        plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'r','LineWidth', 2)
        plot(xaxis, a*(1-cos(2*(pi*xaxis/globaldomain(1))))+ b      ,'b:','LineWidth', 2)
        plot(xaxis,-a*(1-cos(2*(pi*xaxis/globaldomain(1))))+globaldomain(2)-b,'b:','LineWidth', 2)
    end

    %Plot time evolution of velocity and fluxes
    for m =tstart:tskip:tend
        m

        [mass_flux,mass_snapshot]=read_mflux('./mflux','./msnap',resultfile_dir,m);
        [velocity_snapshot, ...
        velocity_flux,   ...
        pressure_surface ...
        F_ext                  ] = read_vflux(m,resultfile_dir,gnbins,nd);
        [velocity_snapshot_pt  ] = read_vflux(m+1,resultfile_dir,gnbins,nd);
        dvdt = velocity_snapshot_pt(:,:,:,1) - velocity_snapshot(:,:,:,1);


        xaxis = linspace(0,globaldomain(1),gnbins(1));
        dv = trapz(mean(nozzleCV.*dvdt(:,:,:,1),3),2);
        %mf = trapz(mean(nozzleCV.*mass_flux(:,:,:,1),3),2);
        mf = sum(sum(mass_flux(:,:,:,1),3),2)/(globaldomain(2)*globaldomain(3));
        vf = trapz(mean(nozzleCV.*velocity_flux(:,:,:,1,1),3),2);
        vf1 = trapz(mean(velocity_flux(:,:,:,1,1),3),2);

        %plot instantanous values
        plot(xaxis,dv,'LineWidth',1,'Color',[0.6 0.6 0.6])
        hold all
        plot(xaxis,mf,'--','LineWidth',1,'Color',[0.6 0.6 0.6])
        plot(xaxis,-vf,'-.','LineWidth',1,'Color',[0.0 0.0 0.0])

        %Collect time average
        if (m > start)
            ave_count = ave_count + 1;
            ave_dv = ave_dv + dv;
            ave_mf = ave_mf + mf;
            ave_vf = ave_vf + vf;
            ave_vf1 = ave_vf1 + vf1;

            %Plot averaged values
            plot(xaxis,ave_dv/ave_count,'LineWidth',2.5,'Color',[0.7 0.7 0.7])
            hold all
            plot(xaxis,ave_mf/ave_count,'-','LineWidth',2.5,'Color',[0.5 0.5 0.5])
            %plot(xaxis,smooth(ave_mf/ave_count,0.5),'-','LineWidth',2,'Color',[0.5 0.5 0.5])
            %plot(xaxis,ybot/5+4.7,'r','LineWidth',3)

            xlabel('x','FontSize',16); ylabel('Mass Flux','FontSize',16)
            set(gca,'FontSize',16);
            set(gca,'FontName','Times'); box on;

            hold all
            %plot(xaxis,-ave_vf/ave_count,'x','LineWidth',1,'Color',[0.4 0.4 0.4])
            plot(xaxis,-ave_vf1/ave_count,'x','LineWidth',1,'Color',[0.4 0.4 0.4])

            plot(xaxis,-smooth(ave_vf1/ave_count,0.08),'-','LineWidth',2,'Color',[0.0 0.0 0.0])
            %plot(xaxis,1.6*(ybot-13),'r','LineWidth',3)

            %Attmpt to fit curve
%             xstart = 1; xend = 200;
%             a = 0.25*globaldomain(2)/2;
%             b = 42.4 %40  %tethdistbot(2) + 0.5*binsize(2);
%             %b = 42; 
%             c = 3.886 %0.25;

            %Linear
%             %dissipation = (1-c*xaxis(xstart:xend)/globaldomain(1));
%             %Exponential
%             dissipation = exp(-c*xaxis(xstart:xend)/globaldomain(1));
%             %Cosine wall function
%             wall_profile = (-a*(1-cos(2*pi*xaxis(xstart:xend)/globaldomain(1)))) ;
%             %Sine function (derivative of cosine)
%             wall_profile = ((2*pi/globaldomain(1))*a*(1-sin(2*pi*xaxis(xstart:xend)/globaldomain(1)))) ;
% 
%             plot(xaxis(xstart:xend),dissipation.*wall_profile +b,'k','LineWidth',2)

            xlabel('x','FontSize',16); ylabel('Momentum Flux','FontSize',16)
            set(gca,'FontSize',16);
            set(gca,'FontName','Times'); box on;

        end

        xlabel('x','FontSize',16); ylabel('mass Flux','FontSize',16)
        set(gca,'FontSize',16);
        set(gca,'FontName','Times'); box on;
        %axis([-5.0 305 -2 25])
        %legend('dvdt', 'mass flux', 'mom flux','Orientation','horizontal','location','South','boxoff')
        hold off; drawnow; pause(0.01);

        %Store Temporal evolution of mass/momentum flux at inlet and throat
        %of nozzle
        dvdt_CV(m)     = mean(trapz(mean(dvdt(8:100 ,:,:,1),3),2),1);
        mf_inlet(m)    = trapz(mean(mass_flux(8 ,:,:,1),3),2);
        mf_throat(m)   = trapz(mean(mass_flux(100,:,:,1),3),2);
        vf_inlet(m)    = trapz(mean(velocity_flux(8 ,:,:,1,1),3),2);
        vf_throat(m)   = trapz(mean(velocity_flux(100,:,:,1,1),3),2);
        ps_inlet(m)    = trapz(mean(pressure_surface(8 ,:,:,1,1),3),2);
        ps_throat(m)   = trapz(mean(pressure_surface(100,:,:,1,1),3),2);
        Fext_CV(m)     = mean(trapz(mean(F_ext(8:100 ,:,:,1),3),2),1);
    end


    figure
    %Plot momentum Fluxes
    hold all
    plot(xaxis,-ave_vf1/ave_count,'x','LineWidth',1,'Color',[0.4 0.4 0.4])

    %Attmpt to fit cubic polynomial
    p3 = -3.455e-6;
    p2 = 0.002244;
    p1 =-0.4253;
    p0 =50.02;
    f = p3*xaxis.^3 + p2*xaxis.^2 + p1*xaxis.^1 +p0;

    plot(xaxis(18:end-10),f(18:end-10),'k','LineWidth',2)
    %smth_ave_vf1 = -smooth(ave_vf1/ave_count,0.1)
    %plot(xaxis(15:end-10),smth_ave_vf1(15:end-10),'-','LineWidth',2,'Color',[0.0 0.0 0.0])

    xlabel('x','FontSize',16); ylabel('Momentum Flux','FontSize',16)
    set(gca,'FontSize',16);
    set(gca,'FontName','Times'); box on;
	savefig('./nozzle_steamwise_flux','png','eps')



    figure
    plot(dvdt_CV,'-','Color',[.6 .6 .6],'lineWidth', 2)
    hold all
    %plot((mf_inlet(1:tskip:end)-mf_throat(1:tskip:end)))
    plot((vf_inlet(1:tskip:end)-vf_throat(1:tskip:end)),'k--', 'markersize', 5,'lineWidth', 2)
    plot((ps_inlet(1:tskip:end)-ps_throat(1:tskip:end))/100, 'kx','lineWidth', 0.5,'MarkerSize',5)
    plot(Fext_CV, 'k-','lineWidth', 0.5)
    legend('dvdt', ...            %'\int \rho u dS^+ - \int \rho u dS^- ',  ...
           '\int \rho u u dS^+ - \int \rho u u dS^- ',    ....
            '\int\sigma dS^+ - \int\sigma dS^- ', ...
            ' F_{ext}', 'Orientation','Vertical','Location','NorthEast')
    legend('boxoff')


    clear all
    close all

    %Read Header file
	resultfile_dir     = '/home/es205/results/md_results/fortran/3D_code/parallel/results/converge_diverge/'
	resultfile_dir_htr = '/home/es205/results/md_results/fortran/3D_code/parallel/results/converge_diverge/high_temporal_resolution/'

    read_header
    N_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));

    %Range of timesteps
    tstart = 1;
    tskip = 1;
    tend = N_records-3;

    read_header

    for m =tstart:tskip:tend
        m

        [velocity_snapshot,velocity_flux] = read_vflux(m,resultfile_dir,gnbins,nd);
        vf_inlet(m)    = trapz(mean(velocity_flux(8 ,:,:,1,1),3),2);
        vf_throat(m)   = trapz(mean(velocity_flux(100,:,:,1,1),3),2);

        [velocity_snapshot,velocity_flux] = read_vflux(m,resultfile_dir_htr,gnbins,nd);
        vf_inlet_htr(m)    = trapz(mean(velocity_flux(8 ,:,:,1,1),3),2);
        vf_throat_htr(m)   = trapz(mean(velocity_flux(100,:,:,1,1),3),2);

        taxis(m) = m * Nvflux_ave * delta_t;
    end

    figure 

    hold all
    plot(taxis/10,(vf_inlet_htr(1:tskip:end)-vf_throat_htr(1:tskip:end)),'-', 'markersize', 5,'lineWidth', 2,'Color',[0.7 0.7 0.7])
    plot(taxis,(vf_inlet(1:tskip:end)-vf_throat(1:tskip:end)),'-', 'markersize', 5,'lineWidth', 2,'Color',[0.0 0.0 0.0])

 
end













%       %Loops through and collect data from many CV
%         for n = 1:size(xstepup,2)-1
% 
%             LV_botx = xstepup(n);           LV_topx = xstepup(n+1);
%             LV_boty = ytop(xstepup(n+1));   LV_topy = ybot(xstepup(n)); 
%             LV_botz = 1;                    LV_topz = gnbins(3);
% 
%             [VS(m,n,:),VF(m,n,:,:),PS(m,n,:,:),facesize,Fext(m,n,:)]  ...
%                 = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
%                                            LV_boty, LV_topy, ...
%                                            LV_botz, LV_topz,m,resultfile_dir);
% 
%             [VSp1(m,n,:)]...
%                 = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
%                                            LV_boty, LV_topy, ...
%                                            LV_botz, LV_topz,m+1,resultfile_dir);
% 
%             xaxis = linspace(0,globaldomain(1),size(xstepup,2)-1);
%             DV(m,n,:) =  (VSp1(m,n,:) - VS(m,n,:))/(delta_t*Nvflux_ave);
%             
%         end
% 
% 
%         %plot instantanous values
%         plot(xaxis,DV(m,:,1),'LineWidth',2,'Color',[0.6 0.6 0.6])
%         hold all
%         plot(xaxis,VF(m,:,1,1),'-.','LineWidth',2,'Color',[0.6 0.6 0.6])
%         plot(xaxis,PS(m,:,1,1),'--','LineWidth',2,'Color',[0.6 0.6 0.6])
% 
% 
% 
%     
%     %Establish and Plot domain set-up
%     Domain_setup
%     %setup stress direction and face
%     ixyz = 1;
%     jxyz = 2;
%     Nvflux_records = Nsteps / (Nvflux_ave);
%     skip =2; n = 1; tol = 10000*eps;
%     external_force_flag = 1
%     %Check CV are satisfied
%     for m =1:skip:Nvflux_records-3
%         m
%         
%         [mass_flux,mass_snapshot]=read_mflux('./mflux','./msnap',resultfile_dir,m);
% 
%         %plot(ytop)
%         %hold all
%         %plot(ybot)
%         %scatter(xstepup,ytop(xstepup))
%         %scatter(xstepdown,ytop(xstepdown))
%         %scatter(xstepup,ybot(xstepup))
%         %scatter(xstepdown,ybot(xstepdown))
%         %hold off
% 
%         %Loops through and collect data from many CV
%         for n = 1:size(xstepup,2)-1
%             %scatter(xstepup(n),ytop(xstepup(n)))
%             %hold all
%             %scatter(xstepup(n+1),ytop(xstepup(n+1)))
% 
%             LV_botx = xstepup(n);           LV_topx = xstepup(n+1);
%             LV_boty = ytop(xstepup(n+1));   LV_topy = ybot(xstepup(n)); 
%             LV_botz = 1;                    LV_topz = gnbins(3);
% 
%             [VS(m,n,:),VF(m,n,:,:),PS(m,n,:,:),facesize,Fext(m,n,:)]  ...
%                 = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
%                                            LV_boty, LV_topy, ...
%                                            LV_botz, LV_topz,m,resultfile_dir);
% 
%             [VSp1(m,n,:)]...
%                 = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
%                                            LV_boty, LV_topy, ...
%                                            LV_botz, LV_topz,m+1,resultfile_dir);
% 
%              totalflux  =  ((VF(m,n,:,1)+VF(m,n,:,4)))/(binsize(1)) ...
%                           +((VF(m,n,:,2)+VF(m,n,:,5)))/(binsize(2)) ...
%                           +((VF(m,n,:,3)+VF(m,n,:,6)))/(binsize(3));
%              totalpressure  =  ((PS(m,n,:,1)-PS(m,n,:,4)))/(binsize(1)) ...
%                               +((PS(m,n,:,2)-PS(m,n,:,5)))/(binsize(2)) ...
%                               +((PS(m,n,:,3)-PS(m,n,:,6)))/(binsize(3));
%          
%              dvelocitydt =  (VSp1(m,n,:) - VS(m,n,:))/(delta_t*Nvflux_ave);
%              %Verify that CV momentum is exactly conservative
%              conserved(m,n) = ( squeeze(sum(totalpressure(:))) ...
%                                -squeeze(sum(totalflux(:)))     ...
%                                -squeeze(sum(dvelocitydt(:)))   ...
%                                +squeeze(sum(Fext(:))/prod(binsize)));
% 
% 
%             a(m,n,1) = squeeze(sum(totalpressure(:))) ;
%             a(m,n,2) = squeeze(sum(dvelocitydt(:)));
%             a(m,n,3) = squeeze(sum(totalflux(:)));
%             a(m,n,4) = squeeze(sum(Fext(:)/prod(binsize)));
% 
%         end      
% 
% %         plot(VSp1(m,:,1)-VS(m,:,1),'ko')
% %         hold all
% %         plot(VF(m,:,1,1),'bx'); plot(-VF(m,:,1,4),'bs')
% %         plot(PS(m,:,1,1)/10,'rx'); plot(PS(m,:,1,4)/10,'rs')
% %         plot(conserved(m,:),'k-')
% %         drawnow; pause(0.2)
% %         hold off
% 
% %         %Plot all components
% %         plot(sum(a(:,:),3),'g')   %du/dt
% %         hold on
% %         plot(sum(a(:,:),3)/prod(binsize),'k') %F_{ext}
% %         plot(-sum(a(:,:),3),'r') % - d\rho uu/dr - dP/dr
% %         plot(-sum(a(:,:),3),'b') %  dsigma/dr
% %         legend('du/dt','F_{ext}','- d\rhouu/dr - dP/dr','d\sigma/dr')
%     end
% 
%     %Adjust as multiplied by delta_t before write out
%     a = a*delta_t;
% 
%     figure
% 	for n = 1:size(xstepup,2)-1
%         plot(a(1:end,n,2)/delta_t)   %du/dt =
%         hold all
%         plot(-a(1:end,n,3) + a(1:end,n,1)+a(1:end,n,4),'r--') % - d\rho uu/dr - dP/dr + dsigma/dr + F_ext
%         title(['d\rhou/dt and \nabla \cdot \Pi vs time in CV '])
%         hold off; pause(2.0)
%     end
        
% Check CV conservation on the large volume ensconcing the nozzle
% x = 1 to end/2
% y = 1 to end
% z = 1 to end
read_header
%Establish and Plot domain set-up
Domain_setup
%setup stress direction and face
ixyz = 1;
jxyz = 2;
Nvflux_records = Nsteps / (Nvflux_ave);
skip =2; n = 1; tol = 100*eps;
external_force_flag = 1
clear a

%Choose nozzle area excluding non-conserved thermostatted regions
LV_botx = 8; LV_topx = 100;%floor(gnbins(1)/2);
LV_boty = 4; LV_topy = gnbins(2)-3;%gnbins(2);
LV_botz = 1; LV_topz = gnbins(3);%gnbins(3);

for m=1:skip:Nvflux_records-3
    [LV_velocity_snapshot,   LV_velocity_flux, ...
    LV_pressure_surface,LV_facesize,LV_F_ext]  ...
        = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
                                   LV_boty, LV_topy, ...
                                   LV_botz, LV_topz,m,resultfile_dir);

    [LV_velocity_snapshot_tplus1] ...
        = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
                                   LV_boty, LV_topy, ...
                                   LV_botz, LV_topz,m+1,resultfile_dir);

    % Calculate total CV flux and change in mass - NOTE binsize used
    % as during sum all dx/dy/dz are identical and cancel leaving
    % a single dx,dy or dz
    totalflux =((LV_velocity_flux(:,1)+LV_velocity_flux(:,4)))/(binsize(1)) ...
              +((LV_velocity_flux(:,2)+LV_velocity_flux(:,5)))/(binsize(2)) ...
              +((LV_velocity_flux(:,3)+LV_velocity_flux(:,6)))/(binsize(3));
    totalpressure =((LV_pressure_surface(:,1)-LV_pressure_surface(:,4)))/(binsize(1)) ...
                  +((LV_pressure_surface(:,2)-LV_pressure_surface(:,5)))/(binsize(2)) ...
                  +((LV_pressure_surface(:,3)-LV_pressure_surface(:,6)))/(binsize(3));

    %totalpressure = totalpressure*delta_t
    dvelocitydt =  (LV_velocity_snapshot_tplus1(:)  ...
        - LV_velocity_snapshot(:))/(delta_t*Nvflux_ave);

    %Verify that CV momentum is exactly conservative
    conserved = ( squeeze(sum(totalpressure(:))) ...
                 -squeeze(sum(totalflux(:)))     ...
                 -squeeze(sum(dvelocitydt(:)))   ...
                 +squeeze(sum(LV_F_ext(:))/prod(binsize)));

    if (conserved > tol || conserved < -tol)
        Error(n) =  conserved;

        display(strcat('Error in CV momentum conservation =', ...
            num2str(Error(n)*100),'% - beginning debug'));

        %Log temporal evolution over 100 timesteps
        skip = 1;
        %Save conservation time plot for a cell
        a(n,1,:) = totalpressure(:);
        a(n,2,:) = dvelocitydt(:);
        a(n,3,:) = totalflux(:);
        a(n,4,:) = LV_F_ext(:);
        n = n + 1;
        if (n == 100 || m == Nvflux_records-3)

            %Adjust as multiplied by delta_t before write out
            a = a*delta_t;

            %Plot all components
            plot(sum(a(1:end,2,:),3),'g')   %du/dt
            hold on
            plot(sum(a(1:end,4,:),3)/prod(binsize),'k') %F_{ext}
            plot(-sum(a(1:end,3,:),3),'r') % - d\rho uu/dr - dP/dr
            plot(-sum(a(1:end,1,:),3),'b') %  dsigma/dr
            legend('du/dt','F_{ext}','- d\rhouu/dr - dP/dr','d\sigma/dr')

            %Plot 100 steps of a single cell
            plot(sum(a(1:end,2,:),3))   %du/dt =
            hold all
            plot(-sum(a(1:end,3,:),3) + sum(a(1:end,1,:),3) ...
                 +sum(a(1:end,4,:),3)/prod(binsize),'r--') % - d\rho uu/dr - dP/dr + dsigma/dr
            title(['d\rhou/dt and \nabla \cdot \Pi vs time in CV ', ...
                num2str(LV_botx),'to',num2str(LV_topx),',', ...
                num2str(LV_boty),'to',num2str(LV_topy),',', ...
                num2str(LV_botz),'to',num2str(LV_topz)])
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
end









%Plot only values inside nozzle
[nozzleCV,ytop,ybot,xstepup,xstepdown] = ...
           nozzle_CV(   -globaldomain(1)/2,globaldomain(1)/2, ...
                        -globaldomain(2)/2,globaldomain(2)/2, ...
                        -globaldomain(3)/2,globaldomain(3)/2, ...
                        -globaldomain(1)/2,globaldomain(1)/2, ...
                        -globaldomain(2)/2+cellsidelength(2)+tethdistbot(2), ...
                         globaldomain(2)/2-cellsidelength(2)+tethdisttop(2), ...
                        -globaldomain(3)/2,globaldomain(3)/2, ...
                         gnbins(1),gnbins(2),gnbins(3),0.125,0);



clear a b
%Choose nozzle area excluding wall
for m=1:skip:Nvflux_records-3
    m
    for n=2:size(xstepup,2)-1
        LV_botx = xstepup(n);           LV_topx = xstepup(n+1);
        LV_boty = ytop(xstepup(n+1))+1; LV_topy = ybot(xstepup(n))-1; 
        LV_botz = 1;                    LV_topz = gnbins(3);

        %for m=1:skip:Nvflux_records-3
            [LV_velocity_snapshot,   LV_velocity_flux, ...
            LV_pressure_surface,LV_facesize,LV_F_ext]  ...
                = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
                                           LV_boty, LV_topy, ...
                                           LV_botz, LV_topz,m,resultfile_dir);

            [LV_velocity_snapshot_tplus1] ...
                = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
                                           LV_boty, LV_topy, ...
                                           LV_botz, LV_topz,m+1,resultfile_dir);

            % Calculate total CV flux and change in mass - NOTE binsize used
            % as during sum all dx/dy/dz are identical and cancel leaving
            % a single dx,dy or dz
            totalflux =((LV_velocity_flux(:,1)+LV_velocity_flux(:,4)))/(binsize(1)) ...
                      +((LV_velocity_flux(:,2)+LV_velocity_flux(:,5)))/(binsize(2)) ...
                      +((LV_velocity_flux(:,3)+LV_velocity_flux(:,6)))/(binsize(3));
            totalpressure =((LV_pressure_surface(:,1)-LV_pressure_surface(:,4)))/(binsize(1)) ...
                          +((LV_pressure_surface(:,2)-LV_pressure_surface(:,5)))/(binsize(2)) ...
                          +((LV_pressure_surface(:,3)-LV_pressure_surface(:,6)))/(binsize(3));

            b(m,n,1,:,:) = LV_velocity_flux(:,:);
            b(m,n,2,:,:) = LV_pressure_surface(:,:);

            %totalpressure = totalpressure*delta_t
            dvelocitydt =  (LV_velocity_snapshot_tplus1(:)  ...
                - LV_velocity_snapshot(:))/(delta_t*Nvflux_ave);

            %Verify that CV momentum is exactly conservative
            conserved = ( squeeze(sum(totalpressure(:))) ...
                         -squeeze(sum(totalflux(:)))     ...
                         -squeeze(sum(dvelocitydt(:)))   ...
                         +squeeze(sum(LV_F_ext(:))/prod(binsize)));

            %if (conserved > tol || conserved < -tol)
            %    Error(n) =  conserved;

            %    display(strcat('Error in CV momentum conservation =', ...
            %        num2str(Error(n)*100),'% - beginning debug'));

                %Log temporal evolution over 100 timesteps
            %    skip = 1;
                %Save conservation time plot for a cell
                a(m,n,1,:) = totalpressure(:);
                a(m,n,2,:) = dvelocitydt(:);
                a(m,n,3,:) = totalflux(:);
                a(m,n,4,:) = LV_F_ext(:);


%                 if (n == 100 || m == Nvflux_records-3)
% 
%                     %Adjust as multiplied by delta_t before write out
%                     a = a*delta_t;
% 

% 
%                     %Plot 100 steps of a single cell
%                     plot(sum(a(1:end,2,:),3))   %du/dt =
%                     hold all
%                     plot(-sum(a(1:end,3,:),3) + sum(a(1:end,1,:),3) ...
%                          +sum(a(1:end,4,:),3)/prod(binsize),'r--') % - d\rho uu/dr - dP/dr + dsigma/dr
%                     title(['d\rhou/dt and \nabla \cdot \Pi vs time in CV ', ...
%                         num2str(LV_botx),'to',num2str(LV_topx),',', ...
%                         num2str(LV_boty),'to',num2str(LV_topy),',', ...
%                         num2str(LV_botz),'to',num2str(LV_topz)])
%                     %Plot error
%                     figure
%                     plot(Error)
%                     hold all
%                     plot(0:size(Error,2),ones(1,size(Error,2)+1)*eps('double'))
%                     title(['Maximum Error in whole domain vs time'])
% 
%                     break
%                 end

            %else
            %    display(strcat('CV momentum conservation OK'));
            %end
        %end

    end


    %Plot all components each timestep
    xCV_loc = xstepup(2:end-1) + diff(xstepup(2:end))/2;

    CV_area = -(ytop(xstepup(n+1))+1-ybot(xstepup(n))-1)*binsize(2)*binsize(3); 
    plot(xCV_loc,sum(a(m,2:end,2,:),4),'g')   %du/dt
    hold on
    plot(xCV_loc,sum(a(m,2:end,4,:),4)/prod(binsize),'k') %F_{ext}
    plot(xCV_loc,-sum(a(m,2:end,3,:),4),'r') % - d\rho uu/dr - dP/dr
    plot(xCV_loc,-sum(a(m,2:end,1,:),4),'b') %  dsigma/dr
    legend('du/dt','F_{ext}','- d\rhouu/dr - dP/dr','d\sigma/dr')
    hold off; drawnow; pause(0.01)
end

%Plot average change along domain
plot(-mean(b(:,2:end,1,1,1),1)/CV_area,'-','color',[0.0 0.0 0.0],'LineWidth',2)
hold all
plot( mean(b(:,2:end,1,1,4),1)/CV_area,'-','color',[0.6 0.6 0.6],'LineWidth',2)
plot( mean(b(:,2:end,2,1,1),1)/CV_area,'--','color',[0.0 0.0 0.0],'LineWidth',2)
plot( mean(b(:,2:end,2,1,4),1)/CV_area,'--','color',[0.6 0.6 0.6],'LineWidth',2)
legend('- \int ( \rho u_x u_x- P ) \cdot dS_x^+', ...
       '- \int ( \rho u_x u_x- P ) \cdot dS_x^-', ...
       '- \int \sigma \cdot dS_x^+', ...
       '- \int \sigma \cdot dS_x^-','orientation','horizontal')
xlabel('x location','FontSize',16); ylabel('Stress','FontSize',16)
set(gca,'FontName','Times'); box on; %caxis([0.0 1.1])
set(gca,'FontSize',16);  


plot(mean(sum(a(:,:,2,1),4),1))
hold all
plot(mean(sum(a(:,:,4,1),4),1))
plot(mean(sum(a(:,:,3,1),4),1))
plot(mean(sum(a(:,:,1,1),4),1))
legend('du/dt','F_{ext}','- d\rhouu/dr - dP/dr','d\sigma/dr')

%Plot time history of average 
CV_no_in = 9; CV_no_out = 9;
plot(-mean(b(:,CV_no_in,1,1,1),4)/CV_area,'-','color',[0.0 0.0 0.0],'LineWidth',2)
hold all
plot( mean(b(:,CV_no_out,1,1,4),4)/CV_area,'-','color',[0.6 0.6 0.6],'LineWidth',2)
plot( mean(b(:,CV_no_in,2,1,1),4)/CV_area,'--','color',[0.0 0.0 0.0],'LineWidth',2)
plot( mean(b(:,CV_no_out,2,1,4),4)/CV_area,'--','color',[0.6 0.6 0.6],'LineWidth',2)
legend('- \int ( \rho u_x u_x- P ) \cdot dS_x^+', ...
       '- \int ( \rho u_x u_x- P ) \cdot dS_x^-', ...
       '- \int \sigma \cdot dS_x^+', ...
       '- \int \sigma \cdot dS_x^-','orientation','horizontal')
xlabel('Time','FontSize',16); ylabel('Stress','FontSize',16)
set(gca,'FontName','Times'); box on; %caxis([0.0 1.1])
set(gca,'FontSize',16);  

CV_no = 9
%dudt
hLine=plot(sum(a(:,CV_no,2,1),4),'-.','color',[0.8 0.8 0.8])
hold all
plot(smooth(sum(a(:,CV_no,2,1),4),0.7,'loess'),'-.','color',[0.4 0.4 0.4],'LineWidth',3)
set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');

%F_ext
plot(sum(a(:,CV_no,4,1),4),'color',[0.2 0.2 0.2])

%rho u u dS
hLine=plot(sum(a(:,CV_no,3,1),4),'-','color',[0.8 0.8 0.8])
hold all
plot(smooth(sum(a(:,CV_no,3,1),4),0.7,'loess'),'-','color',[0.4 0.4 0.4],'LineWidth',3)
set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');

%sigma dS
hLine=plot(sum(a(:,CV_no,1,1),4),'--','color',[0.6 0.6 0.6])
hold all
plot(smooth(sum(a(:,CV_no,1,1),4),0.7,'loess'),'--','color',[0.0 0.0 0.0],'LineWidth',3)
set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');

legend('du/dt','F_{ext}','- d\rhouu/dr - dP/dr','d\sigma/dr','orientation','horizontal')  

xlabel('Time','FontSize',16); ylabel('Flux','FontSize',16)
set(gca,'FontName','Times'); box on; %caxis([0.0 1.1])
set(gca,'FontSize',16);  

axis([0 170 -5 30])
savefig('./nozzle_CV_time_evolution','png','eps')











       
        
%%
        %Load momentum flux values for current timestep -- include external
        %forces if relevant
%         if (external_force_flag)
%             [velocity_snapshot(:,:,:,:), ...
%                 velocity_flux(:,:,:,:,:),   ...
%                 pressure_surface(:,:,:,:,:) ...
%                 F_ext(:,:,:,:)                  ] = read_vflux(m,resultfile_dir,gnbins,nd);
%         else
%             [velocity_snapshot(:,:,:,:), ...
%                 velocity_flux(:,:,:,:,:),   ...
%                 pressure_surface(:,:,:,:,:) ] = read_vflux(m,resultfile_dir,gnbins,nd);
%             F_ext = zeros(size(velocity_snapshot));
%         end
        
       % [velocity_snapshot_tplus1(:,:,:,:)] = read_vflux(m+1,resultfile_dir,gnbins,nd);

        %Plot x dependant properties
%         plot(squeeze(sum(sum(pressure_surface(:,:,:,1,1)/10,2),3)))
%         hold all
%         plot(squeeze(sum(sum(velocity_flux(:,:,:,1,1),2),3)))
%         plot(squeeze(sum(sum(sum(F_ext(:,:,:,:)/prod(binsize),2),3),4)))
%         plot(squeeze(sum(sum(sum(velocity_snapshot_tplus1(:,:,:,:)-velocity_snapshot,2),3),4)))
%         plot(sum(sum(mass_flux(:,:,:,1),2),3))
%         legend('pressure surface/10', 'velocity flux','F ext','dvelocitydt','mass flux')
%         savefig(['nozzle_fluxes', num2str(m), '.png'],'png')
%         hold off




%%
        
        % Calculate total CV flux and change in mass - NOTE binsize used
        % as during sum all dx/dy/dz are identical and cancel leaving
        % a single dx,dy or dz
%        totalflux =((LV_velocity_flux(:,1)+LV_velocity_flux(:,4)))/(binsize(1)) ...
%                   +((LV_velocity_flux(:,2)+LV_velocity_flux(:,5)))/(binsize(2)) ...
%                   +((LV_velocity_flux(:,3)+LV_velocity_flux(:,6)))/(binsize(3));
%         totalpressure =((LV_pressure_surface(:,1)-LV_pressure_surface(:,4)))/(binsize(1)) ...
%                       +((LV_pressure_surface(:,2)-LV_pressure_surface(:,5)))/(binsize(2)) ...
%                       +((LV_pressure_surface(:,3)-LV_pressure_surface(:,6)))/(binsize(3));
%         
%         %totalpressure = totalpressure*delta_t
%         dvelocitydt =  (LV_velocity_snapshot_tplus1(:)  ...
%             - LV_velocity_snapshot(:))/(delta_t*Nvflux_ave);
%         
        %Plot flux in and out
        %         figure(fig2);
        %          h=slice(squeeze(velocity_flux(:,:,:,1,1)),[],[],[5]);
        %          view([90 90]); axis 'tight';
        %          set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
        %          colorbar;  drawnow; pause(0.1)
        %          figure(fig1);
        %         h=slice(squeeze(pressure_surface(:,:,:,1,1)),[],[],[5]);
        %         view([90 90]); axis 'tight';
        %         set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
        %         colorbar;  drawnow; pause(0.1)
               
%         figure(fig1);
%         nozzleCV(nozzleCV == 0) = NaN;
%         h=slice(squeeze(nozzleCV.*(-mass_flux(:,:,:,1))),[],[],[5]);
%         view([90 90]); axis 'tight';
%         set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
%         colorbar;  drawnow; pause(0.1)
%         figure(fig2);
        %F_ext(abs(F_ext) < 0.001) = NaN;
        %h=slice(squeeze(F_ext(:,:,:,1)),[],[],[5]);
        %view([90 90]); axis 'tight';
        %set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)
        %colorbar;  drawnow; pause(0.1)
        %plot(squeeze(sum(sum(mass_flux(:,:,:,1),2),3)/(gnbins(2)*binsize(2)*gnbins(3)*binsize(3))))

        
