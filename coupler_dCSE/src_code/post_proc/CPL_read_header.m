%Open simulation header file and load all simulation variables to workspace
%clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../../results"');
end

%Import date from header file
cd(resultfile_dir);
coupler_header = importdata('./coupler_header',';');
cd (pwdir);

%Loop through all elements and extract variables
for i = 1:size(coupler_header.data,1)
    varname =  coupler_header.textdata{i,2};
    evalc([varname '= coupler_header.data(i)']);
end
