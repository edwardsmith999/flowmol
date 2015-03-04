%Read stress measurements at time "read_time"
%Input of form
%read_vflux(read_time,resultfile_dir,gnbins,nd,filename1,filename2,filename3)
% read_time - value of record to read. note that snapshot 1 is taken before
% evolution so all snapshots are one ahead of fluxes (read_time+1 is read
% by this routine). Forces are taken at time before current flux reading so
% (read_time-1 is read by this routine)
% resultfile_dir - location of files
% gnbins - number of global bins
% nd - number of dimensions
% filename1,filename2,filename3 - vflux, psurface and vsnap files
%Output variables
%[velocity_snapshot_t,velocity_flux_t,pressure_surface_t]


function[velocity_snapshot_t,velocity_flux_t,pressure_surface_t,F_ext_t] = read_vflux(read_time,resultfile_dir,gnbins,nd,filename1,filename2,filename3,filename4)

if (exist('filename1') == 0)
    filename1 = './vflux';
end
if (exist('filename2') == 0)
    filename2 = './psurface';
end

if (exist('filename3') == 0)
    filename3 = './vsnap';
end

if (exist('filename4') == 0)
    filename4 = './Fext';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
nd = 3;
Ncubeface = 6;
datasize = gnbins(1)*gnbins(2)*gnbins(3)*nd*Ncubeface;
bytes = 8;

%Load momentum flux CV data
flux_read_time = read_time;
cd(resultfile_dir);
fid = fopen(filename1,'r','ieee-le');
%Check file exists
if (fid == -1)
    error(strcat(filename1,' does not exist in results'))
end
fseek(fid, bytes*datasize*flux_read_time, 'bof');
vflux = fread(fid,datasize,'double');
velocity_flux_t = reshape(vflux,gnbins(1),gnbins(2),gnbins(3),nd,Ncubeface);
fclose(fid);


%Load surface pressure CV data
stress_read_time = read_time - 1;
fid = fopen(filename2,'r','ieee-le');
%Check file exists
if (fid == -1)
    error(strcat(filename2,' does not exist in results'))
end
fseek(fid, bytes*datasize*stress_read_time, 'bof');
psurface = fread(fid,datasize,'double');
pressure_surface_t = reshape(psurface,gnbins(1),gnbins(2),gnbins(3),nd,Ncubeface);
fclose(fid);


%Load momentum snapshot CV data
snap_read_time = read_time ;
fid = fopen(filename3,'r','ieee-le');
%Check file exists
if (fid == -1)
    error(strcat(filename3,' does not exist in results'))
end
fseek(fid, bytes*datasize*snap_read_time/Ncubeface, 'bof');
vsnap = fread(fid,datasize/Ncubeface,'double');
velocity_snapshot_t = reshape(vsnap,gnbins(1),gnbins(2),gnbins(3),nd);
fclose(fid);

%Load External force CV data (if it is requested)
if (nargout == 4)
    force_read_time = read_time ;
    fid = fopen(filename4,'r','ieee-le');
    %Check file exists
    if (fid == -1)
        error(strcat(filename4,' does not exist in results'))
    end
    fseek(fid, bytes*datasize*force_read_time/Ncubeface, 'bof');
    F_ext = fread(fid,datasize/Ncubeface,'double');
    F_ext_t = reshape(F_ext,gnbins(1),gnbins(2),gnbins(3),nd);
    fclose(fid);
end