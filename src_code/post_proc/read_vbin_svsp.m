%Read velocity average output from simulation
clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
%Read simulation properties from header file and calculate simulation
%properties
read_header;
Nmass_records = (Nsteps-initialstep) / (tplot * Nmass_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./mbins_s','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('mbins file does not exist in results')
end
mbins = fread(fid,'int32');
mass_bins_s = reshape(mbins,gnbins(1),gnbins(2),gnbins(3),Nmass_records);
fclose(fid);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./mbins_p','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('mbins file does not exist in results')
end
mbins = fread(fid,'int32');
mass_bins_p = reshape(mbins,gnbins(1),gnbins(2),gnbins(3),Nmass_records);
fclose(fid);

max(max(max(mass_bins_p(:,:,:,1) - mass_bins_s(:,:,:,1))))


Nvel_records = (Nsteps-initialstep) / (tplot * Nvel_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./vbins_s','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vbins file does not exist in results')
end
velbins = fread(fid,'double');
vel_bins_s = reshape(velbins,gnbins(1),gnbins(2),gnbins(3),nd,Nvel_records);
fclose(fid);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./vbins_p','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vbins file does not exist in results')
end
velbins = fread(fid,'double');
vel_bins_p = reshape(velbins,gnbins(1),gnbins(2),gnbins(3),nd,Nvel_records);
fclose(fid);

max(max(max(vel_bins_p(:,:,:,1,1) - vel_bins_s(:,:,:,1,1))))

