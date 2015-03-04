%Read velocity slice average output from simulation
%Input of form
%read_vslice(filename,resultfile_dir,read_time)
% read_time - value of record to read.
% resultfile_dir - location of files
% filename - vslice file
%Output variables
%[vel_slice]

function[vel_slice]=read_vslice(filename,resultfile_dir,read_time)

%Store Present Working directory and check for results directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_header

%Load mass CV data
cd(resultfile_dir);
fid = fopen(filename,'r','n');
cd (pwdir);

%Check file exists
if (fid == -1)
    error('vslice file does not exist in results')
end

if (exist('read_time','var') == 0)
    % Read everything
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
    velslice = fread(fid,'double');
    vel_slice = reshape(velslice,gnbins(velocity_outflag),nd,Nvel_records);
else
    %Read data at specified time
    datasize = gnbins(velocity_outflag)*nd;
    bytes = 8;
    fseek(fid, bytes*datasize*(read_time-1), 'bof');
    velslice = fread(fid,datasize,'double');
    vel_slice = reshape(velslice,gnbins(velocity_outflag),nd);
end
fclose(fid);

%Set Nan values to zero
vel_slice(isnan(vel_slice)) = 0;

% - - - Divide velocity by mass and area - - -
if (exist('normalise','var') == 1)
    if (normalise == 1)
        m_slice = mslice('./mslice',resultfile_dir,read_time);
        for ixyz = 1:nd
            vel_slice(:,ixyz,:)  = squeeze(vel_slice(:,ixyz,:))/m_slice(:,:);
        end
    elseif (normalise == 2)
        %Get two directions orthogonal to slice direction
        jxyz = mod(velocity_outflag+1,3)+1;
        kxyz = mod(velocity_outflag,3)+1;
        for ixyz = 1:nd
            vel_slice(:,ixyz,:) = squeeze(vel_slice(:,ixyz,:))...
                              ./(nbins(jxyz)*gnbins(kxyz)*Nvel_ave*density*binsize(1)*binsize(2)*binsize(3));
        end
    else
        disp(strcat('normalise value',num2str(normalise),' not recognised'))
    end
end


