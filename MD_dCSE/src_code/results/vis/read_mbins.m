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
Nmass_records = Nsteps / (tplot * Nmass_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./mbins','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('mbins file does not exist in results')
end
mbins = fread(fid,'int32');
mass_bins = reshape(mbins,globalnbins(1),globalnbins(2),globalnbins(3),Nmass_records);
fclose(fid);

%sliceomatic(mass_bins(:,:,:,2))
