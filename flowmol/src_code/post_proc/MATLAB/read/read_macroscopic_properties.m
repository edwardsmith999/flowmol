%Reads macroscopic properties from macroscopic_properties
%file and stores in an array

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../../results';
    display('setting results file to default "./../../results"');
end

%Import date from header file
cd(resultfile_dir);
macroscopic_properties = importdata('./macroscopic_properties',';');
cd (pwdir);