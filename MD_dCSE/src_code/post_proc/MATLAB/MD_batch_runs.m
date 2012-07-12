%Batch Processing Script to run fortran code and process results
clear all
close all

%Store Present Working directory
codefile = '/home/es205/codes/coupled/MD_dCSE/src_code/';
cd (codefile) %Go to code file

%Write starting date to csv
dlmwrite('./post_proc/MATLAB/CV_virial_out.csv',clock,'-append');

%=========Run Md executable================
%for xcells = 5:10
%for ycells = 5:20
%for zcells = 3:10
for i = 1:13

	density = 0.828 + (i-1)*0.002
    %zcells = ycells;

    cd (codefile) %Go to code file

    %Change size of domain in input file
    %[status] = change_line(codefile,'INITIALNUNITS',[xcells,ycells,zcells]);
    [status] = change_line(codefile,'DENSITY',density,'double');

    %Setup restart file in liquid state
    %[status] = change_line(codefile,'NSTEPS',[10000],'integer');
    [status] = system('./md.exe');
    
    %Run main case
    %[status] = change_line(codefile,'NSTEPS',[5000],'integer');
    %[status] = system('./md.exe -r ./results/final_state');
    
    % Output pressure values to csv file
    cd('./post_proc/MATLAB')
    %plot_CV_virial_pressure
    read_macroscopic_properties

    %Copy to results array
    range = 6:size(macroscopic_properties.data,1);
    result(1)  = density;
    result(2)  = mean(macroscopic_properties.data(range,2));
    result(3)  = mean(macroscopic_properties.data(range,3));
    result(4)  = mean(macroscopic_properties.data(range,4));
    result(5)  = mean(macroscopic_properties.data(range,5));
    result(6)  = mean(macroscopic_properties.data(range,6));
    result(7)  = mean(macroscopic_properties.data(range,7));
    result(8)  = mean(macroscopic_properties.data(range,8));

    %Append CSV file
    dlmwrite('CV_virial_out.csv',result,'-append');
    
    clearvars -except i xcells ycells zcells codefile

%end
end

