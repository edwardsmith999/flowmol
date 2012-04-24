%Read mass average output from simulation
clear all

%Store Present Working directory
pwdir = pwd;
if (exist(resultfile_dir) == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
read_header;
Nvflux_records = Nsteps / (tplot * Nvflux_ave);
Ncubeface = 6;

%%
%Load mass flux CV data
cd(resultfile_dir);
fid = fopen('./vflux','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vflux file does not exist in results')
end
vflux = fread(fid,'double');
velocity_flux = reshape(vflux,globalnbins(1),globalnbins(2),globalnbins(3),nd,Ncubeface,Nvflux_records);
fclose(fid);
%%
%Load surface pressure CV data
cd(resultfile_dir);
fid = fopen('./psurface','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('psurface file does not exist in results')
end
psurface = fread(fid,'double');
pressure_surface = reshape(psurface,globalnbins(1),globalnbins(2),globalnbins(3),nd,Ncubeface,Nvflux_records);
fclose(fid);
%%
%Load mass snapshot CV data
cd(resultfile_dir);
fid = fopen('./vsnap','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vsnap file does not exist in results')
end
vsnap = fread(fid,'double');
velocity_snapshot = reshape(vsnap,globalnbins(1),globalnbins(2),globalnbins(3),nd,Nvflux_records+1);
fclose(fid);
%%

%Verify adjacent surface values are equal
%Random Bin
for i=1:nd
    randbin(i) = round(rand*(globalnbins(i)-1))+1;
end
scatter(0:globalnbins(1)-1,velocity_flux(:,randbin(2),randbin(3),1,1,4)+...
                        pressure_surface(:,randbin(2),randbin(3),1,1,4),'rx')
hold on
scatter(globalnbins(1),velocity_flux(1,randbin(2),randbin(3),1,1,4)+...
                       pressure_surface(1,randbin(2),randbin(3),1,1,4),'rx')
scatter(1:globalnbins(1),-velocity_flux(:,randbin(2),randbin(3),1,4,4)+...
                       pressure_surface(:,randbin(2),randbin(3),1,4,4),'bo')
scatter(0,-velocity_flux(end,randbin(2),randbin(3),1,4,4)+...
                       pressure_surface(end,randbin(2),randbin(3),1,4,4),'bo')
                   
% %Calculate total CV flux and change in mass
 totalflux     =((velocity_flux(:,:,:,:,1,:)+velocity_flux(:,:,:,:,4,:)))/(binsize(1)) ...
               +((velocity_flux(:,:,:,:,2,:)+velocity_flux(:,:,:,:,5,:)))/(binsize(2)) ...
               +((velocity_flux(:,:,:,:,3,:)+velocity_flux(:,:,:,:,6,:)))/(binsize(3));
 totalpressure =((pressure_surface(:,:,:,:,1,:)-pressure_surface(:,:,:,:,4,:)))/(binsize(1)) ...
               +((pressure_surface(:,:,:,:,2,:)-pressure_surface(:,:,:,:,5,:)))/(binsize(2)) ...
               +((pressure_surface(:,:,:,:,3,:)-pressure_surface(:,:,:,:,6,:)))/(binsize(3));
  
%totalpressure = totalpressure*delta_t           
 for i =1:Nvflux_records
     dvelocitydt(:,:,:,:,i) =  velocity_snapshot(:,:,:,:,i+1) - velocity_snapshot(:,:,:,:,i);
 end
 
 %Verify that CV is exactly conservative
 figure
 %Random Bin
 for i=1:nd
     randbin(i) = round(rand*(globalnbins(i)-1))+1;
 end
plot(squeeze(dvelocitydt(randbin(1),randbin(2),randbin(3),1,:)));
hold on
%plot(squeeze(dvelocitydt(randbin(1),randbin(2),randbin(3),2,:)));
%plot(squeeze(dvelocitydt(randbin(1),randbin(2),randbin(3),3,:)));

plot(-squeeze(totalflux(randbin(1),randbin(2),randbin(3),1,:)),'--r');
%plot(-squeeze(totalflux(randbin(1),randbin(2),randbin(3),2,:)),'--r');
%plot(-squeeze(totalflux(randbin(1),randbin(2),randbin(3),3,:)),'--r');

plot(squeeze(totalpressure(randbin(1),randbin(2),randbin(3),1,:)),'-.c');
%plot(squeeze(totalpressure(randbin(1),randbin(2),randbin(3),2,:)),'-.c');
%plot(squeeze(totalpressure(randbin(1),randbin(2),randbin(3),3,:)),'-.c');

display(strcat('Maximum error in CV conservation =', ... 
   num2str(max(max(max(max(max(squeeze(totalpressure)-squeeze(totalflux)-squeeze(dvelocitydt))))))*...
           min(min(min(min(min(squeeze(totalpressure)-squeeze(totalflux)-squeeze(dvelocitydt))))))*100),'%'));
% %sliceomatic(mass_flux(:,:,:,1,8))
