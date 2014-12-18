%Plot continuum velocity slice vs analytical solution
clear all
close all

%Find results files
resultfile_dir = './../results';

%Read data for continuum plot
[continuum_velslice] = read_continuum_vslice(resultfile_dir);
read_continuum_header
continuum_tau_xx=read_continuum_tauslice(resultfile_dir,'continuum_tauslice_xx');
continuum_tau_xy=read_continuum_tauslice(resultfile_dir,'continuum_tauslice_xy');
continuum_tau_yx=read_continuum_tauslice(resultfile_dir,'continuum_tauslice_yx');
continuum_tau_yy=read_continuum_tauslice(resultfile_dir,'continuum_tauslice_yy');
Ncontinuum_records = floor(continuum_Nsteps / continuum_tplot);

scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);

u_wall = 1.0;%0.5*(continuum_velslice(ny+1,1) + continuum_velslice(ny+2,1));

spectral_res = 1.0;
dy = 1/(ny+2);
xaxis = linspace(0,1,ny+2); %1/(ny+1):1/(ny+1):1-1/(ny+1);
%xaxis = linspace(0,1,ny+2); %1/(ny+1):1/(ny+1):1-1/(ny+1);
%xaxis = -1/(ny+2):1/(ny):1+1/(ny+2);
xaxis_analy = linspace(dy/2,1-dy/2,spectral_res*ny);
m = 1; %range = [4, 16, 32, 64,  512];
range = [2^0, 2.^[2:6] ]/(continuum_delta_t*continuum_tplot);
Re = (u_wall * rho * 1 ) /meu ;
m = range(1);
for i = 1:Ncontinuum_records

    %Time -0.5 to record value at half time interval
    t = (m)*continuum_delta_t*continuum_tplot;
    
    %Plot Velocity Contour
    set(0,'currentfigure',fig1)
    analy = couette_analytical_fn(t,Re,[1,0],ly, ...
                                    spectral_res*ny,'top');
    h(1) = plot(xaxis_analy,analy,'k','LineWidth',3);

    hold on 
    
    %Plot analytical velocity contour
    h(2) = plot(xaxis,continuum_velslice(:,m),'s','Color',[.5 .5 .5],'MarkerSize',7,'LineWidth',1);
    
    %Labels etc
    c = {'CFD','Analytical'}
    order=[2 1];
    legend (h(:),c{order},'location','NorthWest'); legend('boxoff')
    xlabel('y/H','FontSize',16); ylabel('U_x/U','FontSize',16)
    set(gca,'FontSize',16);
    set(gca,'FontName','Times'); box on;
    %plottitle=num2str(t,'%10.6f');
    %title(strcat('Plot after ',plottitle,' seconds'));
    axis([-0.1 1.1 -0.1 1.1])
    %savefig(strcat('Velocity_',num2str(i)),'png','eps')
    %hold off
  
    
    %Plot Stress Contour
    set(0,'currentfigure',fig2)
    continuum_stress=continuum_vslice_to_stress(continuum_velslice,meu,ly,ny);
    plot(xaxis(2:ny+1),continuum_stress(1:ny,m),'s','Color',[.5 .5 .5]);
    axis([-0.1 1.1 -0.05 1.5]);
    hold on
    plot(xaxis(2:ny+1),-meu*continuum_tau_xy(1:ny,:,m),'o','Color',[.3 .3 .3]);

    analy=couette_analytical_stress_fn(t,Re,1,ly,spectral_res*ny,'top');
    plot(xaxis_analy,analy*rho,'r');
    legend ('CFD','Analytical','location','NorthWest'); legend('boxoff')
    xlabel('y/H'); ylabel('stress')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after ',plottitle,' seconds'));
    hold off
    %savefig(strcat('Stress_',num2str(i)),'png')
    pause(0.5);
    %     %M(n) = getframe(fig1); %Store Frames for video
    if (m == Ncontinuum_records-3)
        set(0,'currentfigure',fig1)
        savefig('CFD_couette_velocity','png','eps')
        break
    end
    m = range(i);
    if (m > Ncontinuum_records-3 || i == size(range,2))
        m = Ncontinuum_records-3
    end
    
end
