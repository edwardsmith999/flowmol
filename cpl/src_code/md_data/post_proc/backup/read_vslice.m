%Read mass average output from simulation
%clear all

%Store Present Working directory and check for results directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_mslice
Nvel_records = floor(Nsteps / (tplot * Nvel_ave));

%Load mass CV data
cd(resultfile_dir);
fid = fopen('./vslice','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('vslice file does not exist in results')
end
vslice = fread(fid,'double');
vel_slice = reshape(vslice,globalnbins(velocity_outflag),nd,Nvel_records);
fclose(fid);

%Get two directions orthogonal to slice direction
jxyz = mod(velocity_outflag+1,3)+1;
kxyz = mod(velocity_outflag,3)+1;

%Divide velocity by mass
for ixyz = 1:nd
    vel_sliceV(:,ixyz,:) = squeeze(vel_slice(:,ixyz,:))...
                      ./(nbins(jxyz)*nbins(kxyz)*Nvel_ave*density*binsize(1)*binsize(2)*binsize(3));
    vel_slice(:,ixyz,:)  = squeeze(vel_slice(:,ixyz,:));%./(mass_slice);
end

%Set Nan values to zero
vel_slice(isnan(vel_slice)) = 0;

