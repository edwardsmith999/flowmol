%Open simulation header file and load all simulation variables to workspace
%clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

cd('resultfile_dir');
continuum_header = importdata('./continuum_header',';');
cd (pwdir);

%Loop through all elements and extract variables
for i = 1:size(continuum_header.data,1)
    varname =  continuum_header.textdata{i,2};
    eval([varname '= continuum_header.data(i)']);
end

%If simulation is not finished, set number of steps to last
%output step of the simulation
cd(resultfile_dir);
if (exist('./simulation_progress') == 1)
    continuum_Nsteps = ceil((importdata('./simulation_progress'))/50);
end
cd (pwdir);
