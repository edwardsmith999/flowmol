%Read mass average output from simulation
function[mass_bins]=read_mbins(filename)

if (exist('filename') == 0)
    filename = './mbins';
end
%clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
read_header;
Nmass_records = (Nsteps-initialstep) / (tplot * Nmass_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen(filename,'r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error(strcat(filename,'file does not exist in results'))
end
mbins = fread(fid,'int32');
mass_bins = reshape(mbins,gnbins(1),gnbins(2),gnbins(3),Nmass_records);
fclose(fid);

%sliceomatic(mass_bins(:,:,:,2))
