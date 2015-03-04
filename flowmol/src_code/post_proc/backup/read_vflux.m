%Read stress measurements at time "read_time"

function[velocity_snapshot_t,velocity_flux_t,pressure_surface_t] = read_vflux(read_time,resultfile_dir,globalnbins,nd)

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
Ncubeface = 6;
datasize = globalnbins(1)*globalnbins(2)*globalnbins(3)*nd*Ncubeface;
bytes = 8;

%Load momentum flux CV data
cd(resultfile_dir);
fid = fopen('./vflux','r','ieee-le');
%Check file exists
if (fid == -1)
    error('vflux file does not exist in results')
end
fseek(fid, bytes*datasize*read_time, 'bof');
vflux = fread(fid,datasize,'double');
velocity_flux_t = reshape(vflux,globalnbins(1),globalnbins(2),globalnbins(3),nd,Ncubeface);
fclose(fid);


%Load surface pressure CV data
fid = fopen('./psurface','r','ieee-le');
%Check file exists
if (fid == -1)
    error('psurface file does not exist in results')
end
fseek(fid, bytes*datasize*read_time, 'bof');
psurface = fread(fid,datasize,'double');
pressure_surface_t = reshape(psurface,globalnbins(1),globalnbins(2),globalnbins(3),nd,Ncubeface);
fclose(fid);


%Load momentum snapshot CV data
fid = fopen('./vsnap','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vsnap file does not exist in results')
end
fseek(fid, bytes*datasize*read_time/Ncubeface, 'bof');
vsnap = fread(fid,datasize/Ncubeface,'double');
velocity_snapshot_t = reshape(vsnap,globalnbins(1),globalnbins(2),globalnbins(3),nd);
fclose(fid);