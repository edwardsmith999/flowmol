%Read energy measurements at time "read_time"

function[energy_snapshot_t,energy_flux_t,energy_surface_t,Fvext_t] = read_eflux(read_time,resultfile_dir,gnbins,filename1,filename2,filename3,filename4)

if (exist('filename1') == 0)
    filename1 = './eflux';
end
if (exist('filename2') == 0)
    filename2 = './esurface';
end

if (exist('filename3') == 0)
    filename3 = './esnap';
end

if (exist('filename4') == 0)
    filename4 = './Fvext';
end


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
fid = fopen(filename1,'r','ieee-le');
%Check file exists
if (fid == -1)
    error('eflux file does not exist in results')
end
fseek(fid, bytes*datasize*read_time, 'bof');
eflux = fread(fid,datasize,'double');
energy_flux_t = reshape(eflux,gnbins(1),gnbins(2),gnbins(3),Ncubeface);
fclose(fid);


%Load surface power CV data
fid = fopen(filename2,'r','ieee-le');
%Check file exists
if (fid == -1)
    error('esurface file does not exist in results')
end
fseek(fid, bytes*datasize*read_time, 'bof');
esurface = fread(fid,datasize,'double');
energy_surface_t = reshape(esurface,gnbins(1),gnbins(2),gnbins(3),Ncubeface);
fclose(fid);


%Load energy snapshot CV data
fid = fopen(filename3,'r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('esnap file does not exist in results')
end
fseek(fid, bytes*datasize*read_time/Ncubeface, 'bof');
esnap = fread(fid,datasize/Ncubeface,'double');
energy_snapshot_t = reshape(esnap,gnbins(1),gnbins(2),gnbins(3));
fclose(fid);


%Load external power CV data
fid = fopen(filename4,'r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('Fvext file does not exist in results')
end
fseek(fid, bytes*datasize*read_time/Ncubeface, 'bof');
Fvext = fread(fid,datasize/Ncubeface,'double');
Fvext_t = reshape(Fvext,gnbins(1),gnbins(2),gnbins(3));
fclose(fid);

