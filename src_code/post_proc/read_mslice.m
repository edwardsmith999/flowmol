%Read mass average output from simulation
%clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
read_header;
Nmass_records = floor(Nsteps / (tplot * Nmass_ave));

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./mslice','r','ieee-le');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('mslice file does not exist in results')
end
mbins = fread(fid,'int32');
mass_slice = reshape(mbins,gnbins(mass_outflag),Nmass_records);
fclose(fid);

%Get two directions orthogonal to slice direction
jxyz = mod(mass_outflag+1,3)+1;
kxyz = mod(mass_outflag,3)+1;

density_slice = mass_slice/(Nvel_ave*nbins(jxyz)*nbins(kxyz)*binsize(1)*binsize(2)*binsize(3));

%Plot evolution
% for i = 1:Nmass_records
% %   plot(mass_slice(:,i));
%    plot(density_slice(:,i));
%    pause(0.5);
% end
