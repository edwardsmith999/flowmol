%Read mass flux output from simulation
function[mass_flux,mass_snapshot,Nmflux_records]=read_mflux(filename1,filename2,resultfile_dir)

if (exist('filename1') == 0)
    filename1 = './mflux';
end
if (exist('filename2') == 0)
    filename2 = './msnap';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
read_header;
Nmflux_records = (Nsteps-initialstep) / (tplot * Nmflux_ave);
Ncubeface = 6;

%Load mass flux CV data
cd(resultfile_dir);
fid = fopen(filename1,'r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error(strcat(filename1,'file does not exist in results'))
end
mflux = fread(fid,'int32');
if (exist('globalnbins','var') == 1)
    gnbins = globalnbins;
end
mass_flux = reshape(mflux,gnbins(1),gnbins(2),gnbins(3),Ncubeface,Nmflux_records);
fclose(fid);

%Load mass snapshot CV data
cd(resultfile_dir);
fid = fopen(filename2,'r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error(strcat(filename2 ,'file does not exist in results'))
end
msnap = fread(fid,'int32');
mass_snapshot = reshape(msnap,gnbins(1),gnbins(2),gnbins(3),Nmflux_records+1);
fclose(fid);
