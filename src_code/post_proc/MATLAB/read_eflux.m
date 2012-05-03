%Read energy measurements at time "read_time"

function[energy_snapshot_t,energy_flux_t,energy_surface_t] = read_eflux(read_time,resultfile_dir,gnbins)

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
Ncubeface = 6;
datasize = gnbins(1)*gnbins(2)*gnbins(3)*Ncubeface;
bytes = 8;

%Load mass flux CV data
cd(resultfile_dir);
fid = fopen('./eflux','r','ieee-le');
%Check file exists
if (fid == -1)
    error('eflux file does not exist in results')
end
fseek(fid, bytes*datasize*read_time, 'bof');
eflux = fread(fid,datasize,'double');
energy_flux_t = reshape(eflux,gnbins(1),gnbins(2),gnbins(3),Ncubeface);
fclose(fid);


%Load surface power CV data
fid = fopen('./esurface','r','ieee-le');
%Check file exists
if (fid == -1)
    error('esurface file does not exist in results')
end
fseek(fid, bytes*datasize*read_time, 'bof');
esurface = fread(fid,datasize,'double');
energy_surface_t = reshape(esurface,gnbins(1),gnbins(2),gnbins(3),Ncubeface);
fclose(fid);


%Load energy snapshot CV data
fid = fopen('./esnap','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('esnap file does not exist in results')
end
fseek(fid, bytes*datasize*read_time/Ncubeface, 'bof');
esnap = fread(fid,datasize/Ncubeface,'double');
energy_snapshot_t = reshape(esnap,gnbins(1),gnbins(2),gnbins(3));
fclose(fid);