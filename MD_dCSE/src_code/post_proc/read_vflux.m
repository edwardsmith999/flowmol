%Read stress measurements at time "read_time"
function[velocity_snapshot_t,velocity_flux_t,pressure_surface_t] = read_vflux(read_time,resultfile_dir,gnbins,nd,filename1,filename2,filename3)

if (exist('filename1') == 0)
    filename1 = './vflux';
end
if (exist('filename2') == 0)
    filename2 = './psurface';
end

if (exist('filename3') == 0)
    filename3 = './vsnap';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
Ncubeface = 6;
datasize = gnbins(1)*gnbins(2)*gnbins(3)*nd*Ncubeface;
bytes = 8;

%Load momentum flux CV data
cd(resultfile_dir);
fid = fopen(filename1,'r','ieee-le');
%Check file exists
if (fid == -1)
    error(strcat(filename1,'vflux file does not exist in results'))
end
fseek(fid, bytes*datasize*read_time, 'bof');
vflux = fread(fid,datasize,'double');
velocity_flux_t = reshape(vflux,gnbins(1),gnbins(2),gnbins(3),nd,Ncubeface);
fclose(fid);


%Load surface pressure CV data
fid = fopen(filename2,'r','ieee-le');
%Check file exists
if (fid == -1)
    error(strcat(filename2,'file does not exist in results'))
end
fseek(fid, bytes*datasize*read_time, 'bof');
psurface = fread(fid,datasize,'double');
pressure_surface_t = reshape(psurface,gnbins(1),gnbins(2),gnbins(3),nd,Ncubeface);
fclose(fid);


%Load momentum snapshot CV data
fid = fopen(filename3,'r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error(strcat(filename3,'file does not exist in results'))
end
fseek(fid, bytes*datasize*read_time/Ncubeface, 'bof');
vsnap = fread(fid,datasize/Ncubeface,'double');
velocity_snapshot_t = reshape(vsnap,gnbins(1),gnbins(2),gnbins(3),nd);
fclose(fid);