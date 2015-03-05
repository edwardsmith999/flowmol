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
read_mbins
Nvel_records = Nsteps / (tplot * Nvel_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./vbins_s','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vbins file does not exist in results')
end
velbins = fread(fid,'double');
vel_bins = reshape(velbins,globalnbins(1),globalnbins(2),globalnbins(3),nd,Nvel_records);
fclose(fid);
