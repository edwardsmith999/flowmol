%Read Temperature average output from simulation
function[Temp_bins]=read_Tbins(filename,resultfile_dir,read_time)

    if (exist('filename') == 0)
        filename = './Tbins';
    end

    %Store Present Working directory
    pwdir = pwd;
    if (exist('resultfile_dir') == 0)
        resultfile_dir = './../../results';
        display('setting results file to default "./../../results"');
    end

    %Load temperature data
    cd(resultfile_dir);
    fid = fopen('./Tbins','r','n');
    cd (pwdir);

    %Check file exists
    if (fid == -1)
        error('Tbins file does not exist in results')
    end

    % Check if snapshot or all simulation should be read then
    % read simulation properties from header file, calculate 
    % datasize to read and read required data
    read_header
    if (exist('read_time','var') == 0)
        NTemp_records = floor((Nsteps-initialstep) / (tplot * NTemp_ave));
        Tempbins = fread(fid,'double');
        Temp_bins = reshape(Tempbins,gnbins(1),gnbins(2),gnbins(3),NTemp_records);
    else
        datasize = gnbins(1)*gnbins(2)*gnbins(3);
        bytes = 8;
        fseek(fid, bytes*datasize*read_time, 'bof');
        Tempbins = fread(fid,datasize,'double');
        Temp_bins = reshape(Tempbins,gnbins(1),gnbins(2),gnbins(3));
    end

    fclose(fid);

end
