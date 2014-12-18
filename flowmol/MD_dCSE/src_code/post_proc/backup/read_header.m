%Open simulation header file and load all simulation variables to workspace
%clear all

%Store Present Working directory
pwdir = pwd;
if (exist('resultfile_dir') == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%Import date from header file
cd(resultfile_dir);
simulation_header = importdata('./simulation_header',';');
cd (pwdir);

%Loop through all elements and extract variables
for i = 1:size(simulation_header.data,1)
    varname =  simulation_header.textdata{i,2};
    evalc([varname '= simulation_header.data(i)']);
end

%If simulation is not finished, set number of steps to last
%output step of the simulation
cd(resultfile_dir); 
if (exist('./simulation_progress') > 0)
    Nsteps = importdata('./simulation_progress');
end
cd (pwdir);

%Set mass flag to same as velocity flag if zero
if velocity_outflag ~= 0
    mass_outflag = velocity_outflag;
end
