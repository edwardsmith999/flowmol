%Read velocity average output from simulation
function[vel_bins]=read_vbins(filename)

if (exist('filename') == 0)
    filename = './vbins';
end
%clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_header
%mass_bins = read_mbins;
Nvel_records = (Nsteps-initialstep) / (tplot * Nvel_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen(filename,'r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error(strcat(filename,' file does not exist in results'))
end
velbins = fread(fid,'double');
vel_bins = reshape(velbins,gnbins(1),gnbins(2),gnbins(3),nd,Nvel_records);
fclose(fid);

