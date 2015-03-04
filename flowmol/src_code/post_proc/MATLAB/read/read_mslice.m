%Read mass slice average output from simulation
%Input of form
%read_mslice(filename,resultfile_dir,read_time)
% read_time - value of record to read.
% resultfile_dir - location of files
% filename - mslice file
%Output variables
%[m_slice]

function[m_slice]=read_mslice(filename,resultfile_dir,read_time,normalise)

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
read_header;

%Load mass CV data
cd(resultfile_dir);
fid = fopen(filename,'r','ieee-le');
cd (pwdir);

%Check file exists
if (fid == -1)
    error('mslice file does not exist in results')
end


if (exist('read_time','var') == 0)
    % Read everything
    Nmass_records = floor((Nsteps-initialstep) / (tplot * Nmass_ave));
    mbins = fread(fid,'int32');
    m_slice = reshape(mbins,gnbins(mass_outflag),Nmass_records);
else
    %Read data at specified time
    datasize = gnbins(mass_outflag);
    bytes = 4;
    bytes*datasize*(read_time-1);
    fseek(fid, bytes*datasize*(read_time-1), 'bof');
    m_slice = fread(fid,datasize,'int32');
end
fclose(fid);

% - - - Divide mass by area - - -
if (exist('normalise','var') == 1)
    %Get two directions orthogonal to slice direction
    jxyz = mod(mass_outflag+1,3)+1;
    kxyz = mod(mass_outflag,3)+1;
    m_slice = m_slice/(Nvel_ave*gnbins(jxyz)*gnbins(kxyz)*binsize(1)*binsize(2)*binsize(3));
end
