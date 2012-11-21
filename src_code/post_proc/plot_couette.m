%==========================================================================
% Visualization routine for VELOCITY and FLUX field from DNS files
%  ucvcwc.dble.xxxxxx and uuvvww.dble.xxxxxx
%==========================================================================

clear all
close all
fclose('all')

resultfile_dir = '/home/es205/codes/coupled/CFD_dCSE/src_code/results/';

%---Get CFD grid size ----
[ngx, ngy, ngz, Lx, Ly, Lz, dx, dy, dz] = read_report(strcat(resultfile_dir,'report'));

%--- Read grid ----
fid = fopen(strcat(resultfile_dir,'grid.data'),'r','ieee-le.l64');
xpg = fread(fid,[ngx ngy],'double');
ypg = fread(fid,[ngx ngy],'double');
fclose(fid);


%Analytical Solution
u_0 = 1;
spectral_res = 6; 
Re = 0.5; dt = 15;
analy_points = 33; % Number of spectral points
xaxis_analy = 0:1/analy_points:1.00;

%setup figures
scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
fig2 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/5 scrsz(4)/2]);

files = dir('Sub*');

%Load DNS data
[u,v,w,p,stress] = Read_DNS('grid.data',resultfile_dir,ngx-2,ngy-1,ngz-2,Lx,Ly,Lz,13,true);
xaxis = linspace(-dy/2,Ly+dy/2,ngy-1)/Ly;
dlmwrite(strcat(resultfile_dir,'../analystress'),zeros(1,analy_points), 'delimiter', ' ')
for ntime=1:size(files,1)

    % =================================
    %    Plot CFD velocity
    % =================================
    set(0,'currentfigure',fig2)

    %Time -0.5 to record value at half time interval
    if (ntime  == 1)
        t = 0;
    else
        t =(ntime-1.5)*dt;
    end

    %Plot CFD velocity profile
    plot(xaxis,u(:,ntime)/u_0,'s','Color',[.5 .5 .5],'MarkerSize',8,'LineWidth',5);
    hold on
    
    %Plot anayltical solution
    analy = couette_analytical_fn(t,Re,[1,0],Ly,analy_points,'top');
    plot(xaxis_analy,analy/u_0,'k','LineWidth',5);
  
    %Make plot look nice
    %legend ('MD','CFD','CFD_{halo}','location','NorthEast'); legend('boxoff')
    set(gca,'FontSize',20)
    axis([-0.1 1.1 -0.1 1.1]);
    xlabel('y/H'); ylabel('U_x/U')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    hold off

	drawnow
    pause(0.1)

    % =================================
    %    Plot CFD stress
    % =================================
    set(0,'currentfigure',fig1)
    %disp(strcat('Maximum Pressure in domain = ',num2str(max(abs(p(:))))))

    %Plot CFD shear stress
    plot(xaxis,squeeze(stress(1,2,:,ntime)),'s','Color',[.5 .5 .5],'MarkerSize',8,'LineWidth',5)
    hold on

    %Plot anayltical solution
    clear analy
    analy = couette_analytical_stress_fn(t,Re,1,Ly,analy_points);
    %dlmwrite(strcat(resultfile_dir,'../analystress'),analy', 'delimiter', ' ','-append')
    plot(xaxis_analy,analy/Re,'k','LineWidth',5);

    %Make plot look nice
    %legend ('MD','CFD','CFD_{halo}','location','NorthEast'); legend('boxoff')
    set(gca,'FontSize',20)
    axis([-0.1 1.1 -0.01 0.1]);
    xlabel('y/H'); ylabel('\Pi')
    plottitle=num2str(t,'%10.6f');
    title(strcat('Plot after  ',plottitle,' time units'));
    hold off

	drawnow
    pause(0.1)

end
