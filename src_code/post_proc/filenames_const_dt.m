%====================================================
% Return file numbers of DNS snapshots at regular dt
% removing filenames samller than dt while ignoring
% any gaps greater than dt
%====================================================
%
% Usage:
% files = dir('SubDom*'); dt = 40;
% filenames = filenames_const_dt(files,dt);
%

function filenames = filenames_const_dt(files,dt)
% Extract file numbers from filenames
len_name   = regexp(files(1).name,'\D');
len_number = regexp(files(1).name,'\d');

files = struct2cell(files);
files = files(1,:)';

files = cell2mat(files);
fname = files(:,len_name);
fnumber = files(:,len_number);

files = str2num(fnumber);

% Find files at inconsistent dt due to simulation restart
restart_files = files(diff(diff(files) < dt)==-1);

% Remove the restart files from the file list
[files ind] = setdiff(files,restart_files);
fname = fname(ind,:);
fnumber = fnumber(ind,:);

% Construct corrected filenames structure for output
filenames = [fname,fnumber];
filenames = mat2cell(filenames,ones(size(files,1),1),size(filenames,2));
filenames = cell2struct(filenames,'name',2);