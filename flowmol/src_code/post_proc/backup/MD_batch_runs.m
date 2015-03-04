%Batch Processing Script to run fortran code and process results
clear all
close all

%Store Present Working directory
codefile = '/home/es205/codes/issued_codes/svn_lucian/MD_dCSE/src_code/';
cd (codefile) %Go to code file

%Write starting date to csv
dlmwrite('./post_proc/CV_virial_out.csv',clock,'-append');

%=========Run Md executable================
%for xcells = 5:10
%for ycells = 5:20
%for zcells = 3:10
for i = 1:5

    density = 0.5 + (i-1)*0.05
    %zcells = ycells;

    cd (codefile) %Go to code file

    %Change size of domain in input file
    %[status] = change_line(codefile,'INITIALNUNITS',[xcells,ycells,zcells]);
    [status] = change_line(codefile,'DENSITY',density,'double');

    %Setup restart file in liquid state
    [status] = change_line(codefile,'NSTEPS',[10000],'integer');
    [status] = system('./md.exe');
    
    %Run main case
    %[status] = change_line(codefile,'NSTEPS',[5000],'integer');
    %[status] = system('./md.exe -r ./results/final_state');
    
    % Output pressure values to csv file
    cd('./post_proc')
    plot_CV_virial_pressure
    
    clearvars -except i xcells ycells zcells codefile

%end
end

