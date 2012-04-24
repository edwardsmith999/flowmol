%Read Temperature average output from simulation
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
Ntemp_records = Nsteps / (tplot * NTemp_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./Tbins','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vbins file does not exist in results')
end
Tempbins = fread(fid,'double');
Temp_bins = reshape(Tempbins,gnbins(1),gnbins(2),gnbins(3),Ntemp_records);
fclose(fid);
