%Read velocity average output from simulation
function[vel_bins]=read_vbins(filename,resultfile_dir,read_time)

    if (exist('filename') == 0)
        filename = './vbins';
    end

    %Store Present Working directory
    pwdir = pwd;
    if (exist('resultfile_dir') == 0)
        resultfile_dir = './../../results';
        display('setting results file to default "./../../results"');
    end

    %Load velocity data
    cd(resultfile_dir);
    fid = fopen(filename,'r','n');
    cd (pwdir);
    %Check file exists
    if (fid == -1)
        error(strcat(filename,' file does not exist in results'))
    end


    % Check if snapshot or all simulation should be read then
    % read simulation properties from header file, calculate 
    % datasize to read and read required data
    read_header
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
    if (exist('read_time') == 0)
        velbins = fread(fid,'double');
        % Check that data is expected size 
        if (size(velbins,1)/(nd*prod(gnbins)) ~= Nvel_records)
            disp(['Expected Records = ', num2str(Nvel_records), ...
              '. Actual Records = ', num2str(size(velbins,1)/(nd*prod(gnbins))), ...
              '. Using Actual Records for Velocity bins'])
            Nvel_records = size(velbins,1)/(nd*prod(gnbins));
        end
        vel_bins = reshape(velbins,gnbins(1),gnbins(2),gnbins(3),nd,Nvel_records);
    else
        datasize = gnbins(1)*gnbins(2)*gnbins(3)*nd;
        bytes = 8;
        status = fseek(fid, bytes*datasize*read_time, 'bof');
        if (status ~= -1)
            velbins = fread(fid,datasize,'double');
            vel_bins = reshape(velbins,gnbins(1),gnbins(2),gnbins(3),nd);
        else
            error(['Requested record = ' num2str(read_time) ' is greater than maximum = ' num2str(Nvel_records)])
        end
    end
    fclose(fid);

end
