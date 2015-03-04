%Read velocity average output from simulation
function [continuum_velbins] = read_continuum_vbins(resultfile_dir)


%Store Present Working directory
pwdir = pwd;
if (exist(resultfile_dir) == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_continuum_header
Ncontinuum_records = continuum_Nsteps / continuum_tplot;

%Load velocity in x direction CV data
cd(resultfile_dir);
fid = fopen('./continuum_vxbins','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('continuum_vxbins file does not exist in results')
end
continuum_vxbins = fread(fid,'double');
continuum_velbins = reshape(continuum_vxbins,(nx+2),(ny+2),Ncontinuum_records);
fclose(fid);

%Load velocity in y direction CV data
% cd(resultfile_dir);
% fid = fopen('./continuum_vybins','r','n');
% cd (pwdir);
% %Check file exists
% if (fid == -1)
%     error('continuum_vybins file does not exist in results')
% end
% continuum_vybins = fread(fid,'double');
% fclose(fid);


