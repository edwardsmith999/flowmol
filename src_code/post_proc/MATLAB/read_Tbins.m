%Read Temperature average output from simulation
function[Temp_bins]=read_Tbins(filename,resultfile_dir)

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_header
mass_bins = read_mbins(strcat(resultfile_dir,'/mbins'));
NTemp_records = (Nsteps-initialstep) / (tplot * NTemp_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./Tbins','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('Tbins file does not exist in results')
end
Tempbins = fread(fid,'double');
Temp_bins = reshape(Tempbins,gnbins(1),gnbins(2),gnbins(3),NTemp_records-1);
fclose(fid);

end
