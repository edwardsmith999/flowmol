%Read Volume Averaged (VA) pressure output from simulation
function[pressure_VA,pressure_VA_k,pressure_VA_c]=read_pVA(filename1,filename2)

if (exist('filename1') == 0)
    filename1 = './pVA';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_header
Nvel_records = (Nsteps-initialstep) / (tplot * Nstress_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen(filename1,'r','n');
%Check file exists and if not then check for 2 k/c files
if (fid == -1)
    filename = filename1;
    if (exist('filename1') == 0)
        filename1 = './pVA_k'
        filename2 = './pVA_c'
    end
    if (exist('filename2') == 0)
        filename2 = './pVA_c'
    end
    %Read volume averaged kinetic file
    fid = fopen(filename1,'r','n');
    if (fid == -1)
        error(strcat('neither ',filename,' or ', filename1,' file exists in results'))
    end
    pVA_k = fread(fid,'double');
    pressure_VA_k = reshape(pVA_k,gnbins(1),gnbins(2),gnbins(3),nd,nd,Nvel_records);
    %Read volume averaged configurational file
    fid = fopen(filename2,'r','n');
    if (fid == -1)
        error(strcat('neither ',filename,' or ', filename2,' file exists in results'))
    end
    pVA_c = fread(fid,'double');
    pressure_VA_c = reshape(pVA_c,gnbins(1),gnbins(2),gnbins(3),nd,nd,Nvel_records);
    fclose(fid);
    cd (pwdir);
    %Calculate pressure_VA from kinetic and confiurational
    pressure_VA = pressure_VA_k + pressure_VA_c;
else
    %Read pressure_VA from file
    pVA = fread(fid,'double');
    pressure_VA = reshape(pVA,gnbins(1),gnbins(2),gnbins(3),nd,nd,Nvel_records);
    fclose(fid);
    cd (pwdir);
end


% scrsz = get(0,'ScreenSize');
% fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
% 
% set(0,'CurrentFigure',fig1)
% ColOrd = get(gca,'ColorOrder');
% 
% for i=1:nd
%     for j=1:i
%         % Get the color
%         Col = ColOrd(i+j,:);
%         plot(squeeze(pressure_virial(i,j,:)),'Color',Col)
%         hold on
%     end
% end

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
