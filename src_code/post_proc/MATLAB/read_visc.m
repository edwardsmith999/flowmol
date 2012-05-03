%Read velocity average output from simulation
clear all

%Store Present Working directory
pwdir = pwd;
if (exist(resultfile_dir) == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_header
Nvisc_records = Nsteps / (tplot * Nstress_ave * Nvisc_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./visc','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('visc file does not exist in results')
end
visc = fread(fid,'double');
fclose(fid);
