%Read mass flux output from simulation at optional time "read_time"
%Input of form
% read_mflux(filename1,filename2,resultfile_dir,read_time)
%   - filename1,filename2 - mflux and msnap files
%   - resultfile_dir - location of files
%   - read_time - value of record to read. note that snapshot 1 is taken 
%                 before evolution so all snapshots are one ahead of fluxes 

%Output variables
%[mass_flux,mass_snapshot]

function[mass_flux,mass_snapshot]=read_mflux(filename1,filename2,resultfile_dir,read_time)

if (exist('filename1') == 0)
    filename1 = './mflux';
end
if (exist('filename2') == 0)
    filename2 = './msnap';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%Read simulation properties from header file and calculate simulation
%properties
cd(resultfile_dir);
read_header;

%Open files for mass flux and mass snapshot
%Mass flux CV data
fid_mflux = fopen(filename1,'r','ieee-le');
%Check file exists
if (fid_mflux == -1)
    error([filename1,' file does not exist in results'])
end
fid_msnap = fopen(filename2,'r','ieee-le');
%Check file exists
if (fid_msnap == -1)
    error([filename2 ,' file does not exist in results'])
end

%Check if all data to be read or just single time
Ncubeface = 6;
bytes = 4; %Integers

if (exist('read_time'))     % - - Read single time- -

    %Load read_time mass flux CV data
    datasize = gnbins(1)*gnbins(2)*gnbins(3)*Ncubeface;
    fseek(fid_mflux, bytes*datasize*read_time, 'bof');
    mflux = fread(fid_mflux,datasize,'int32');
    if (exist('globalnbins','var') == 1)
        gnbins = globalnbins; %Backward compatibility
    end
    mass_flux = reshape(mflux,gnbins(1),gnbins(2),gnbins(3),Ncubeface);
    fclose(fid_mflux);

    %Load read_time mass snapshot CV data
    datasize = gnbins(1)*gnbins(2)*gnbins(3);
    fseek(fid_msnap, bytes*datasize*read_time, 'bof');
    msnap = fread(fid_msnap,datasize,'int32');
    mass_snapshot = reshape(msnap,gnbins(1),gnbins(2),gnbins(3));
    fclose(fid_msnap);

else    % - - Read all data - -

    %Check if cv conservation is employed
    if (exist('cv_conserve'))
        cv_conserve = cv_conserve;
    else
        cv_conserve = 0;
    end
    if (cv_conserve == 1)
        Nmflux_records = (Nsteps-initialstep) / (Nmflux_ave);
    else
        Nmflux_records = (Nsteps-initialstep) / (tplot*Nmflux_ave);
    end


    %Check file size is not too large - stop
    maxGB = 2;
    file = bytes*gnbins(1)*gnbins(2)*gnbins(3)*Ncubeface*Nmflux_records;
    if (file > maxGB*1024^3)
        disp(['Filesize in read_mflux is non-crompulently large -- ' ,num2str(file/1024^3)  , 'GB!'])
        reply = input('Do you want to continue y/n? ', 's');
        if (reply == 'y')
            %continue
        else
            error('Stopping')
        end
    end

    %Load all mass flux CV data
    mflux = fread(fid_mflux,'int32');
    if (exist('globalnbins','var') == 1)
        gnbins = globalnbins; %Backward compatibility
    end
    mass_flux = reshape(mflux,gnbins(1),gnbins(2),gnbins(3),Ncubeface,Nmflux_records);
    fclose(fid_mflux);

    %Load all mass snapshot CV data
    msnap = fread(fid_msnap,'int32');
    mass_snapshot = reshape(msnap,gnbins(1),gnbins(2),gnbins(3),Nmflux_records+1);
    fclose(fid_msnap);

end


