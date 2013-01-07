%Read mass average output from simulation
function[mass_bins,density_bins]=read_mbins(filename,resultfile_dir,read_time)

if (exist('filename') == 0)
    filename = './mbins';
end
%clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%Check mass CV file exists an open
cd(resultfile_dir);
fid = fopen(filename,'r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error(strcat(filename,'file does not exist in results'))
end

% Check if snapshot or all simulation should be read then
% read simulation properties from header file, calculate 
% datasize to read and read required data
read_header;
if (exist('read_time') == 0)
    Nmass_records = floor((Nsteps-initialstep) / (tplot * Nmass_ave));
    mbins = fread(fid,'int32');
    mass_bins = reshape(mbins,gnbins(1),gnbins(2),gnbins(3),Nmass_records);
else
    datasize = gnbins(1)*gnbins(2)*gnbins(3);
    bytes = 4;
    fseek(fid, bytes*datasize*read_time, 'bof');
    mbins = fread(fid,datasize,'int32');
    mass_bins = reshape(mbins,gnbins(1),gnbins(2),gnbins(3));
end
fclose(fid);
density_bins = mass_bins*((tplot)/(Nmass_ave*binsize(1)*binsize(2)*binsize(3)));

end

%sliceomatic(mass_bins(:,:,:,2))
