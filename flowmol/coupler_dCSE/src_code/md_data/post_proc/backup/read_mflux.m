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
Nmflux_records = Nsteps / (tplot * Nmflux_ave);
Ncubeface = 6;

%Load mass flux CV data
cd(resultfile_dir);
fid = fopen('./mflux','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('mflux file does not exist in results')
end
mflux = fread(fid,'int32');
mass_flux = reshape(mflux,globalnbins(1),globalnbins(2),globalnbins(3),Ncubeface,Nmflux_records);
fclose(fid);

%Load mass snapshot CV data
cd(resultfile_dir);
fid = fopen('./msnap','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('msnap file does not exist in results')
end
msnap = fread(fid,'int32');
mass_snapshot = reshape(msnap,globalnbins(1),globalnbins(2),globalnbins(3),Nmflux_records+1);
fclose(fid);

%Calculate total CV flux and change in mass
totalflux = squeeze(sum(mass_flux,4));
for i =1:Nmflux_records
    dmassdt(:,:,:,i) = mass_snapshot(:,:,:,i+1) - mass_snapshot(:,:,:,i);
end

%Verify that CV is exactly conservative

display(strcat('Maximum error in CV conservation =', ... 
   num2str(max(max(max(max(squeeze(totalflux(:,:,:,:)) - squeeze(dmassdt(:,:,:,:))))))*...
           min(min(min(min(squeeze(totalflux(:,:,:,:)) - squeeze(dmassdt(:,:,:,:))))))*100),'%'));
%sliceomatic(mass_flux(:,:,:,1,8))
