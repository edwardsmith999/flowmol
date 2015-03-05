%Read Volume Averaged (VA) pressure output from simulation
function[pressure_VA,pressure_VA_k,pressure_VA_c]=read_pVA(read_time,resultfile_dir,filename1,filename2)

if (exist('filename1') == 0)
    filename1 = './pVA';
end

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_header
datasize = gnbins(1)*gnbins(2)*gnbins(3)*nd*nd;
bytes = 8;
%Nvel_records = (Nsteps-initialstep) / (tplot * Nstress_ave);

%Load mass CV data
cd(resultfile_dir);
fid = fopen(filename1,'r','n');
fclose(fid);
%Check file exists and if not then check for 2 k/c files
if (fid == -1 || strcmp(filename1,strcat(resultfile_dir,'/pVA_k')))
    if (exist('filename1') == 0)
        filename1 = './pVA_k'
        filename2 = './pVA_c'
    end
    if (exist('filename2') == 0)
        filename2 = './pVA_c'
    end
    %Read volume averaged kinetic file
    fid = fopen(filename1,'r','n');
    if (fid == -1)
        error([filename1,' doesnt file exists in results'])
    end
    fseek(fid, bytes*datasize*read_time, 'bof');
    pVA_k = fread(fid,datasize,'double');
    pressure_VA_k = reshape(pVA_k,gnbins(1),gnbins(2),gnbins(3),nd,nd);
    fclose(fid);
    %Read volume averaged configurational file
    fid = fopen(filename2,'r','n');
    if (fid == -1)
        error(strcat('neither ',filename,' or ', filename2,' file exists in results'))
    end
    fseek(fid, bytes*datasize*read_time, 'bof');
    pVA_c = fread(fid,datasize,'double');
    pressure_VA_c = reshape(pVA_c,gnbins(1),gnbins(2),gnbins(3),nd,nd);
    fclose(fid);
    cd (pwdir);
    %Calculate pressure_VA from kinetic and confiurational
    pressure_VA = pressure_VA_k + pressure_VA_c;
else
    %Read pressure_VA from file
    fid = fopen(filename1,'r','n');
    fseek(fid, bytes*datasize*read_time, 'bof');
    pVA = fread(fid,datasize,'double');
    pressure_VA = reshape(pVA,gnbins(1),gnbins(2),gnbins(3),nd,nd);
    fclose(fid);
    cd (pwdir);
end
