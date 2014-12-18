clear all
close all

dirname = '/home/es205/results/MD_continuum_results/code/coupled_couette/111223_CV_FEM/';
%Open results file
filename = strcat(dirname,'fort.99');
fid = fopen(filename);
A = fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g',inf);

eqib_steps = 9;

%Check simulation progress so far
progressfile = strcat(dirname,'results/simulation_progress');
fid = fopen(progressfile);
nsteps = fscanf(fid,'%g');
nsteps = floor(nsteps/100)-eqib_steps;

%Create output arrays
nbins = size(A,1)/(13*nsteps);
b=reshape(A,[13,nbins,nsteps]);

xplots = ceil(sqrt(nbins)); 
yplots = ceil(nbins/xplots);

for i =1:nbins
    %subplot(nbins,i,ceil(i/xplots))
    figure
	col = get(gca,'ColorOrder');
	c(:,i) = smooth(b(13,i,:),200);
	%c(:,i) = squeeze(cumsum(b(6,i,:)))./[1:nsteps]';
	plot(squeeze(abs(c(:,i))),'color',[0,0,0],'LineWidth',5);
    axis([1,nsteps,-0.1,0.5]);
    hold on
    %c(:,i) = smooth(b(5,i,:),200);
    c(:,i) =  squeeze(cumsum(b(5,i,:)))./[1:nsteps]';
    plot(squeeze(abs(c(:,i))),'color',col(mod(i-1,7)+1,:),'LineWidth',3);
    c(:,i) = squeeze(cumsum(b(13,i,:)))./[1:nsteps]';
    plot(squeeze(abs(c(:,i))),'--','color',[0,0,0],'LineWidth',5);
    %c(:,i) = smooth(b(12,i,:),200);
    c(:,i) =  squeeze(cumsum(b(12,i,:)))./[1:nsteps]';
    plot(squeeze(abs(c(:,i))),'color',col(mod(i,7)+1,:),'LineWidth',3);
   legend ('Bottom Setpoint','Bottom Traction','Top Setpoint','Top Traction','location','NorthEast'); legend('boxoff')
   xlabel('time'); ylabel('P_{xy}')
   plottitle=strcat('Cell Number ',num2str(i));
   title(plottitle)
%    savefig(strcat('Control_',num2str(i)),'png')
end
pause()
t_ave = 100; t_ave = ceil(t_ave/2);
figure
for i=0+t_ave:nsteps-t_ave+1
	plot(1:6,mean(b(5,:,i-t_ave+1:i+t_ave-1),3),'-x')
    hold on
    plot(1:6,mean(b(6,:,i-t_ave+1:i+t_ave-1),3),'-s')
    axis([1,7,-0.5,0.5])
    plot(2:7,mean(b(12,:,i-t_ave+1:i+t_ave-1),3),'-x')
    plot(2:7,mean(b(13,:,i-t_ave+1:i+t_ave-1),3),'-s')
    hold off
    pause(0.5)
end
