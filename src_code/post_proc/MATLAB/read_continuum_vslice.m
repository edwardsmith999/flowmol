%Read velocity average output from simulation
%clear all
%close all

%Store Present Working directory
pwdir = pwd;
if (exist(resultfile_dir) == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_continuum_header
Ncontinuum_records = floor(continuum_Nsteps / continuum_tplot);

%Load velocity in x direction CV data
cd(resultfile_dir);
fid = fopen('./continuum_vslice','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('continuum_vslice file does not exist in results')
end

continuum_vslice = fread(fid,'double');
continuum_velslice = reshape(continuum_vslice,(ny+2),Ncontinuum_records);
fclose(fid);

