% General purpose routine to plot various statistics from output of
% Transflow DNS code. Optional out flags can be tuned to output
% slices, (velcoity/vorticity) isosurface, energy spectra, velocity profiles
% autocorrelations and wall region plots. Also allows output as .dx files



clear all
close all

%Optional outputs
contourplots = 4; sliceplane = 1;
channel_profile = 1;
spectrum = 1;

%Parameters
Re_b = 4.0;
visc = 1/Re_b;
spectrum_yloc = 16;

%Directory names
%resultfile_dir = '/home/es205/results/CFD_results/results/Transflow/minimal_channel_laminar/';
%resultfile_dir = '/home/es205/results/CFD_results/results/Transflow/minimal_channel_turbulent_scaled/'
resultfile_dir = '/home/es205/results/CFD_results/results/Transflow/turbulent_couette/couette_scaled/CFD_dCSE/src_code/results/'

%---Get CFD grid size ----
read_grid(strcat(resultfile_dir,'grid.data'),[1 1 1])
[ngx, ngy, ngz, Lx, Ly, Lz, dx, dy, dz] = read_report(strcat(resultfile_dir,'report'));
z = linspace(0,Lz,ngz-2);
[Yz Zz] = meshgrid(ypg(1,:),z);
[Xx Zx] = meshgrid(xpg(:,1),z);
[ugx,ugy,ugz]=meshgrid(linspace(0,Lx,ngx-3),linspace(0,Lz,ngz-3),linspace(0,Ly,ngy-1));




%Setup figures
scrsz = get(0,'ScreenSize');
if (contourplots ~= 0)
    fig1=figure('Position',[    1         scrsz(4)     scrsz(3)/6 scrsz(4)/2]);
    fig5=figure('Position',[    1            1         scrsz(3)/6 scrsz(4)/3]);
end
if (channel_profile == 1)
    fig2=figure('Position',[ scrsz(3)/6   scrsz(4)     scrsz(3)/6 scrsz(4)/2]);
elseif (channel_profile > 1)
    fig2=figure('Position',[ scrsz(3)/6   scrsz(4)     scrsz(3)/6 scrsz(4)/2.6]);
    fig4=figure('Position',[ scrsz(3)/6      1         scrsz(3)/6 scrsz(4)/2.6]);
end
if (spectrum == 1)
    fig3=figure('Position',[ scrsz(3)*2/6 scrsz(4)     scrsz(3)/6 scrsz(4)/2]);
elseif (spectrum == 2)
    fig3=figure('Position',[ scrsz(3)*2/6 scrsz(4)     scrsz(3)/6 scrsz(4)/2.6]);
    fig6=figure('Position',[ scrsz(3)*2/6     1        scrsz(3)/6 scrsz(4)/2.6]);
end


%Main loop
u_sum = 0; v_sum = 0; w_sum = 0; uv_fluct_sum=0; dudy_sum=0;
u_fluct_rmsum=0; v_fluct_rmsum=0; Rxu_fft = zeros(2*(ngx-3),1);
Rxu = zeros(ngx-3,1); Rxv = zeros(ngx-3,1); Rxw = zeros(ngx-3,1);
Rzu = zeros(ngz-3,1); Rzv = zeros(ngz-3,1); Rzw = zeros(ngz-3,1);
cd(resultfile_dir); files = dir('Sub*');
for n=50:size(files,1)

    disp(strcat('Iter_', int2str(n), '_of_', int2str(size(files,1))))
    
    %Read domain
    [u,v,w,p] =  Read_DNS_Subdomain(n,'grid.data',resultfile_dir,ngx-2,ngy-1,ngz-2,1,1,1,4,true);

    %Write velocity field to .dx format which can be read by VMD
    %G3DtoDX(xpg(1:end-1,1),ypg(1,1:end-1),z(1:end-1), permute(u(:,:,2:end-1),[2 3 1]),strcat('./vmd_volumes/DNS',num2str(n),'.dx'),-Lx/2,-Ly/2,-Lz/2)

    %sliceomatic(cav,linspace(0,Lx,ngx-3),linspace(0,Lz,ngz-3),linspace(0,Ly,ngy-1))

    % --- Plot contours ---
        % Plot of u,v,w, and pressure on xy slice 
    if (contourplots == 1)
        figure(fig1);
        subplot(2,2,1); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(u(sliceplane,:,1:end-1)),10); colorbar
        subplot(2,2,2); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(v(sliceplane,:,1:end-1)),10); colorbar
        subplot(2,2,3); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(w(sliceplane,:,1:end-1)),10); colorbar
        subplot(2,2,4); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(p(sliceplane,:,1:end-1)),10); colorbar
        drawnow; pause(0.0001)
    elseif (contourplots == 2)
        % Plot of u,v,w, and pressure on yz slice 
        figure(fig1);
        subplot(2,2,1); contourf(Zz(1:end-1,:),Yz(1:end-1,:),squeeze(u(:,sliceplane,1:end-1)),10); colorbar
        subplot(2,2,2); contourf(Zz(1:end-1,:),Yz(1:end-1,:),squeeze(v(:,sliceplane,1:end-1)),10); colorbar
        subplot(2,2,3); contourf(Zz(1:end-1,:),Yz(1:end-1,:),squeeze(w(:,sliceplane,1:end-1)),10); colorbar
        subplot(2,2,4); contourf(Zz(1:end-1,:),Yz(1:end-1,:),squeeze(p(:,sliceplane,1:end-1)),10); colorbar
        drawnow; pause(0.0001)
    elseif (contourplots == 3)
        % Plot of u,v,w, and pressure on xz slice 
        figure(fig1);
        subplot(2,2,1); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(u(:,:,sliceplane)),10); colorbar
        subplot(2,2,2); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(v(:,:,sliceplane)),10); colorbar
        %hold on; plot(Xx, Zx,'k-'); plot(Xx',Zx', 'k-'); hold off
        subplot(2,2,3); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(w(:,:,sliceplane)),10); colorbar
        subplot(2,2,4); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(p(:,:,sliceplane)),10); colorbar
        drawnow; pause(0.0001)
        figure(fig5); 
        subplot(2,2,1); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(u(:,:,sliceplane+1)),10); colorbar
        subplot(2,2,2); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(v(:,:,sliceplane+1)),10); colorbar
        %hold on; plot(Xx, Zx,'k-'); plot(Xx',Zx', 'k-'); hold off
        subplot(2,2,3); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(w(:,:,sliceplane+1)),10); colorbar
        subplot(2,2,4); contourf(Xx(1:end-1,1:end-1),Zx(1:end-1,1:end-1),squeeze(p(:,:,sliceplane+1)),10); colorbar
        drawnow; pause(0.0001)
        %feather(squeeze(mean(u(:,:,sliceplane),1)),squeeze(mean(v(:,:,sliceplane),1)))
        
    elseif (contourplots == 4)
        % Plot of angular momentum, on xy slice 
        [curlx,curly,curlz,cav] = curl(u,v,w);
        figure(fig1);
        subplot(2,2,1); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(curlx(sliceplane,:,1:end-1)),10); colorbar
        subplot(2,2,2); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(curly(sliceplane,:,1:end-1)),10); colorbar
        subplot(2,2,3); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(curlz(sliceplane,:,1:end-1)),10); colorbar
        subplot(2,2,4); contourf(xpg(1:end-1,:),ypg(1:end-1,:),squeeze(cav(sliceplane,:,1:end-1)),10); colorbar
        drawnow; pause(0.0001)
        figure(fig5); clf;
        %Fixes a 3D plotting bug in matlab running on 2 screens
        set(0,'DefaultFigureRenderer','OpenGL');
        [faces,vertices] = isosurface(ugx,ugy,ugz,curlz, 0.1);
        p=patch('Faces',faces,'Vertices',vertices);
        isonormals(ugx,ugy,ugz,curlz,p);
        set(p,'FaceColor','red','EdgeColor','none');
        [faces,vertices] = isosurface(ugx,ugy,ugz,curlz,-0.1);
        p=patch('Faces',faces,'Vertices',vertices);
        isonormals(ugx,ugy,ugz,curlz,p);
        set(p,'FaceColor','blue','EdgeColor','none');
        daspect([1 1 1]); box on;
        view(3); axis tight
        camlight 
        lighting gouraud
        hold off
    end


    %Energy spectrum at centreline
    if (n>50)

        % Get u fluctuations
        u_sum = u_sum + mean(mean(u,1),2);
        u_mean = u_sum / (n-50) ;
        dudy_sum = dudy_sum + diff(squeeze(mean(mean(u(:,:,1:end-1),1),2)))./diff(ypg(1,1:end))';
        for i=1:size(u,1)
        for j=1:size(u,2)
            u_fluct(i,j,:) = squeeze(u(i,j,:)) - squeeze(u_mean(:));
        end
        end

        % Get v fluctuations
        v_sum = v_sum + mean(mean(v,1),2);
        v_mean = v_sum / (n-50) ;
        for i=1:size(v,1)
        for j=1:size(v,2)
            v_fluct(i,j,:) = squeeze(v(i,j,:)) - squeeze(v_mean(:));
        end
        end

        %Get w fluctuations
        w_sum = w_sum + mean(mean(w,1),2);
        w_mean = w_sum / (n-50) ;
        for i=1:size(w,1)
        for j=1:size(w,2)
            w_fluct(i,j,:) = squeeze(w(i,j,:)) - squeeze(w_mean(:));
        end
        end

        %Collect turbulent statistics
        u_fluct_rmsum = u_fluct_rmsum + abs(u_fluct);
        v_fluct_rmsum = v_fluct_rmsum + abs(v_fluct);
        uv_fluct_sum  = uv_fluct_sum + u_fluct.*v_fluct;

        %Calculate scaling parameter
        dudy = dudy_sum / (n-50);
        dudy_top =  max(dudy); dudy_bot = -min(dudy);
        dudy_ave = 0.5*(dudy_top + dudy_bot);
        u_tau = (dudy_ave*visc)^0.5;
        delta_tau = (visc/dudy_ave)^0.5;

        if (channel_profile > 0)
            figure(fig2);
            plot(ypg(1,:),squeeze(u_mean(1:end-1)),'Color',[.5 .5 .5],'LineWidth',3)
            hold all
            plot(ypg(1,:),squeeze(u(4,4,1:end-1)),'k-.','LineWidth',3)
            plot(ypg(1,:),squeeze(u_fluct(4,4,1:end-1)),'k')
            plot(ypg(1,:),squeeze(v_fluct(4,4,1:end-1)),'k--')
            plot(ypg(1,:),squeeze(w_fluct(4,4,1:end-1)),'k:')
            hold off
            %axis([-0.1 2.1 -0.4 1.5])
            drawnow
            %savefig('channel','png')
        	if (channel_profile == 2)
                figure(fig4);
                uv_fluct = uv_fluct_sum / (n-50);
                plot(ypg(1,:),-squeeze(uv_fluct(4,4,1:end-1))./u_tau^2,'-')
                hold all
                plot(ypg(1,1:end-1),visc*dudy(1:end)./u_tau^2,'o')
                plot(ypg(1,1:end-1),(-squeeze(uv_fluct(4,4,2:end-1))+visc*dudy(1:end))./u_tau^2,'--')
                axis([-0.1 2.1 -1.1 1.1])
                hold off
                drawnow
            elseif (channel_profile == 3)
                u_fluct_rms = u_fluct_rmsum / (n-50);
                v_fluct_rms = v_fluct_rmsum / (n-50);
                plot(ypg(1,:),squeeze(u_fluct_rms(4,4,1:end-1)),'x')
                hold all
                plot(ypg(1,:),squeeze(v_fluct_rms(4,4,1:end-1)),'s')
                hold off
                drawnow
            end
        end
    
        %Collect spectrum information
        uhatx = fft(squeeze(mean(u_fluct(:,:,spectrum_yloc),1)));
        vhatx = fft(squeeze(mean(v_fluct(:,:,spectrum_yloc),1)));
        whatx = fft(squeeze(mean(w_fluct(:,:,spectrum_yloc),1)));
        Exu(:,n-50) =  0.5*uhatx.*conj(uhatx);
        Exv(:,n-50) =  0.5*vhatx.*conj(vhatx);
        Exw(:,n-50) =  0.5*whatx.*conj(whatx);

        %Collect spatial autocorrelation functions at yplus ~7.2 (7.335294532208900) 
        % x autocorrelation
        for i=1:size(u_fluct,2)
            Rxu(i) = Rxu(i) + mean(u_fluct(:,1,21).*u_fluct(:,i,21),1);
            Rxv(i) = Rxv(i) + mean(v_fluct(:,1,21).*v_fluct(:,i,21),1);
            Rxw(i) = Rxw(i) + mean(w_fluct(:,1,21).*w_fluct(:,i,21),1);
        end

        % Calculate autocorrelation using FFT
        len = size(u_fluct,2);
        nfft = 2^nextpow2(2*len-1);
        Rxu_fft(:) = Rxu_fft(:) + ifft(     fft(squeeze(mean(u_fluct(:,:,spectrum_yloc),1)),nfft) .* ...
                                       conj(fft(squeeze(mean(u_fluct(:,:,spectrum_yloc),1)),nfft)) )';

        % z autocorrelation
        for i=1:size(u_fluct,1)
            Rzu(i) = Rzu(i) + mean(u_fluct(1,:,21).*u_fluct(i,:,21),2);
            Rzv(i) = Rzv(i) + mean(v_fluct(1,:,21).*v_fluct(i,:,21),2);
            Rzw(i) = Rzw(i) + mean(w_fluct(1,:,21).*w_fluct(i,:,21),2);
        end

        if (spectrum == 1)
            figure(fig3);
            loglog(xpg(1:end/2,1),Exu(1:end/2,n-50),'k')    
            hold on
            loglog(xpg(1:end/2,1),Exv(1:end/2,n-50),'k--')    
            loglog(xpg(1:end/2,1),Exw(1:end/2,n-50),'k:')  
            hold off
            %axis([100 10000 10^-7 10^-2  ])
            %uhatz = fft(u(:,:,spectrum_yloc),1);
            %Ez = Ez + 0.5*uhatz.*conj(uhatz);
        elseif (spectrum == 2)
            %Plot spatial autocorrelation functions ay yplus ~7.2 (7.335294532208900) 
            figure(fig3);
            plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxu(1:end/2)/(n-50))/(Rxu(1)/(n-50)),'k')
            hold on
            plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxu_fft(1:end/4)/(n-50))/(Rxu_fft(1)/(n-50)),'r')
            plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxv(1:end/2)/(n-50))/(Rxv(1)/(n-50)),'k--')
            plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxw(1:end/2)/(n-50))/(Rxw(1)/(n-50)),'k:')
            plot(xpg(1:(end-1)/2)/delta_tau,zeros(size(xpg(1:(end-1)/2),2),1),'k-.')
            axis([0 300 -0.5 1.1])
            drawnow; pause(0.0001)
            hold off
            figure(fig6);
            plot(z(1:(end-1)/2)/delta_tau,(Rzu(1:end/2)/(n-50))/(Rzu(1)/(n-50)),'k')
            hold on
            plot(z(1:(end-1)/2)/delta_tau,(Rzv(1:end/2)/(n-50))/(Rzv(1)/(n-50)),'k--')
            plot(z(1:(end-1)/2)/delta_tau,(Rzw(1:end/2)/(n-50))/(Rzw(1)/(n-50)),'k:')
            plot(z(1:(end-1)/2)/delta_tau,zeros(size(z(1:(end-1)/2),2),1),'k-.')
            axis([0 85 -1 1.1])
            drawnow; pause(0.0001)
            hold off
        end

    end
end

close all
scrsz = get(0,'ScreenSize');

%Scaling parameters
dudy = dudy_sum / (n-50) ;
dudy_top =  max(dudy); dudy_bot = -min(dudy);
dudy_ave = 0.5*(dudy_top + dudy_bot);
u_tau = (dudy_ave*visc)^0.5;
delta_tau = (visc/dudy_ave)^0.5;
%Get plus units
y_plus = ypg/delta_tau;
u_plus = u_mean/u_tau;

%Plot time averaged Reynolds shear stress
u_fluct_rms = u_fluct_rmsum / (n-50);
v_fluct_rms = v_fluct_rmsum / (n-50);
uv_fluct    = uv_fluct_sum  / (n-50);
figure('Position',[     1     scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
plot(ypg(1,:),-squeeze(uv_fluct(4,4,1:end-1))/u_tau^2,'ks')
hold all
plot(ypg(1,1:end-1),(-squeeze(uv_fluct(4,4,2:end-1))+visc*dudy(1:end))./u_tau^2,'k--','LineWidth',5)
set(gca,'FontSize',20)
savefig('stresses','eps')
figure('Position',[     1     scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
plot(ypg(1,:),squeeze(u_fluct_rms(4,4,1:end-1))/u_tau,'kx')
hold all
plot(ypg(1,:),squeeze(v_fluct_rms(4,4,1:end-1))/u_tau,'ko')
set(gca,'FontSize',20)
savefig('RMS_vel','eps')

%Plot time average energy spectra
figure('Position',[     1     scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);;
C = 100;
loglog(xpg(1:end/2,1)/delta_tau,mean(Exu(1:end/2,1:end).^0.5,2),'k')    
hold on
loglog(xpg(1:end/2,1)/delta_tau,mean(Exv(1:end/2,1:end).^0.5,2),'k--')   
loglog(xpg(1:end/2,1)/delta_tau,mean(Exw(1:end/2,1:end).^0.5,2),'k:') 
loglog(xpg(1:end/2,1)/delta_tau,C*(xpg(1:end/2,1)/delta_tau).^(-5/3),'k.-')
%axis([7 400 1e-6 5 ]); set(gca,'FontSize',20)
savefig('x_energy_spectra','eps')

% Plot correlations
% x correlation
figure('Position',[     1     scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxu(1:end/2)/(n-50))/(Rxu(1)/(n-50)),'k')
hold on
%plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxu_fft(1:end/4)/(n-50))/(Rxu_fft(1)/(n-50)),'r')
plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxv(1:end/2)/(n-50))/(Rxv(1)/(n-50)),'k--')
plot(xpg(1:(end-1)/2,1)/delta_tau,(Rxw(1:end/2)/(n-50))/(Rxw(1)/(n-50)),'k:')
plot(xpg(1:(end-1)/2)/delta_tau,zeros(size(xpg(1:(end-1)/2),2),1),'k-.')
%axis([0 300 -0.5 1.1]); set(gca,'FontSize',20)
drawnow; pause(0.0001)
hold off
savefig('x_correlation','eps')

%z correlation
figure('Position',[     1     scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);;
plot(z(1:(end-1)/2)/delta_tau,(Rzu(1:end/2)/(n-50))/(Rzu(1)/(n-50)),'k')
hold on
plot(z(1:(end-1)/2)/delta_tau,(Rzv(1:end/2)/(n-50))/(Rzv(1)/(n-50)),'k--')
plot(z(1:(end-1)/2)/delta_tau,(Rzw(1:end/2)/(n-50))/(Rzw(1)/(n-50)),'k:')
plot(z(1:(end-1)/2)/delta_tau,zeros(size(z(1:(end-1)/2),2),1),'k-.')
%axis([0 85 -1 1.1]); set(gca,'FontSize',20)
drawnow; pause(0.0001)
hold off
savefig('z_correlation','eps')

%Plot wall layer
figure('Position',[     1     scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);;
visc_sublayer = u_tau*ypg(1,1:end-1)/visc;
loglaw = 2.41*log(y_plus(1,1:end-1)) + 5.4;
semilogx(y_plus(1,1:end-1),squeeze(mean(mean(u_plus(:,:,2:end-1),1),2)),'ko')
hold all
semilogx(y_plus(1,1:end-1),visc_sublayer,'k-')
semilogx(y_plus(1,1:end-1),loglaw,'k--')
axis([ 0 200 0 25 ]); set(gca,'FontSize',20)
savefig('wall_layer','eps')
