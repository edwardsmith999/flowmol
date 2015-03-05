%Read mass average output from simulation

%Store Present Working directory
pwdir = pwd;
if (exist(resultfile_dir) == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_header
Nvel_records = Nsteps / (tplot * Nvel_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./pvirial','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('pvirial file does not exist in results')
end
pvirial = fread(fid,'double');
pressure_virial = reshape(pvirial,nd,nd,Nvel_records);
fclose(fid);

scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);

set(0,'CurrentFigure',fig1)
ColOrd = get(gca,'ColorOrder');

for i=1:nd
    for j=1:i
        % Get the color
        Col = ColOrd(i+j,:);
        plot(squeeze(pressure_virial(i,j,:)),'Color',Col)
        hold on
    end
end

% xlabel('x'); ylabel('y'); zlabel('z'); 
% s=3.0;
% x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*s;
% y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*s;
% z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*s;

% for i = 1:Nvel_records
%     i
%     % Plot CUBE
%     for i=1:6
%         plot3(x(:,i),y(:,i),z(:,i),'k');
%         hold on
%         %axis([-0.1 0.1 -0.1 0.1 -0.1 0.1])
%     end
%     
%     %Plot surface stresses
%     quiver3 (s/2,0,0,pressure_virial(1,1,i),0,0, 'DisplayName', 'stress');
%     quiver3 (s/2,0,0,0,pressure_virial(1,2,i),0, 'DisplayName', 'stress');
%     quiver3 (s/2,0,0,0,0,pressure_virial(1,3,i), 'DisplayName', 'stress');
%     quiver3 (0,s/2,0,pressure_virial(2,1,i),0,0, 'DisplayName', 'stress');
%     quiver3 (0,s/2,0,0,pressure_virial(2,2,i),0, 'DisplayName', 'stress');
%     quiver3 (0,s/2,0,0,0,pressure_virial(2,3,i), 'DisplayName', 'stress');
%     quiver3 (0,0,s/2,pressure_virial(3,1,i),0,0, 'DisplayName', 'stress');
%     quiver3 (0,0,s/2,0,pressure_virial(3,2,i),0, 'DisplayName', 'stress');
%     quiver3 (0,0,s/2,0,0,pressure_virial(3,3,i), 'DisplayName', 'stress');
%     hold off
%     pause(0.5)
% end
