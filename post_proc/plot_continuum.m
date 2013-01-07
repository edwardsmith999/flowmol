%Plot continuum velocity slice vs analytical solution
clear all
close all

%Find results files
resultfile_dir = './../../results';

%Read data for continuum plot
read_continuum_vslice

scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);

u_wall = 1;%0.5*(continuum_velslice(ny+1,1) + continuum_velslice(ny+2,1));

spectral_res = 10;
xaxis = 1/(ny+1):1/(ny+1):1-1/(ny+1);
xaxis = -1/(ny+2):1/(ny):1+1/(ny+2);
xaxis_analy = 0:1/(spectral_res*ny):1;
m = 1;
for i = 1:Ncontinuum_records
    
    %Plot Velocity Contour
    set(0,'currentfigure',fig1)
    hold off
    plot(xaxis,continuum_velslice(1:ny+2,m)/u_wall,'s','Color',[.5 .5 .5]);
    %axis([-0.1 1.1 -0.1 1.1]);
    hold on
    
    %Time -0.5 to record value at half time interval
    t = (m)*continuum_delta_t*continuum_tplot;

    analy = couette_analytical_fn(t,Re,[1,0],ly, ...
                                    spectral_res*ny,'top');
    plot(xaxis_analy,analy,'r');
    legend ('CFD','Analytical','location','NorthWest'); legend('boxoff')
    xlabel('y/H'); ylabel('U_x/U')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after ',plottitle,' seconds'));
    %savefig(strcat('Velocity_',num2str(i)),'png')
    hold off
    
    
    %Plot Stress Contour
    set(0,'currentfigure',fig2)
    continuum_stress=continuum_vslice_to_stress(continuum_velslice,meu,ly,ny);
    plot(xaxis(2:ny+1),continuum_stress(1:ny,m)/meu,'s','Color',[.5 .5 .5]);
    axis([-0.1 1.1 -0.5 0.5]);
    hold on

    analy=couette_analytical_stress_fn(t,Re,1,ly,spectral_res*ny);
    plot(xaxis_analy,-analy,'r');
    legend ('CFD','Analytical','location','NorthWest'); legend('boxoff')
    xlabel('y/H'); ylabel('stress')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after ',plottitle,' seconds'));
    hold off
    %savefig(strcat('Stress_',num2str(i)),'png')
    pause(0.5);
    %     %M(n) = getframe(fig1); %Store Frames for video
    if (m == Ncontinuum_records-3)
        break
    end
    m = m*2
    if (m > Ncontinuum_records-3)
        m = Ncontinuum_records-3
    end
    
end
