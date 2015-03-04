%Read velocity average output from simulation
%clear all
%close all

%Store Present Working directory
pwdir = pwd;
if (exist(resultfile_dir) == 0)
    resultfile_dir = './../results';
    display('setting results file to default "./../results"');
end

%%Read mass ouput incluing simulation properties from header file and calculate simulation
%properties
read_continuum_header
Ncontinuum_records = continuum_Nsteps / continuum_tplot;

%Load velocity in x direction CV data
cd(resultfile_dir);
fid = fopen('./continuum_vxbins','r','n');
cd (pwdir);
%Check file exists
if (fid == -1)
    error('continuum_vxbins file does not exist in results')
end
continuum_vxbins = fread(fid,'double');

% %Load velocity in y direction CV data
% cd(resultfile_dir);
% fid = fopen('./continuum_vybins','r','n');
% cd (pwdir);
% %Check file exists
% if (fid == -1)
%     error('continuum_vybins file does not exist in results')
% end

continuum_vybins = fread(fid,'double');
continuum_velbins = reshape(continuum_vxbins,(nx+2),(ny+2),Ncontinuum_records);
fclose(fid);

continuum_velslice = squeeze(mean(continuum_velbins,1));

% %Plot evolving y profile
% xaxis = 0:1/(ny-1):1;
% n = 1;
% for i=1:Ncontinuum_records
%     n = n*2;
%     if (n > Ncontinuum_records) 
%         break 
%     end
%     plot(xaxis,continuum_velslice(2:ny+1,n),'s');
%     axis([-0.1 1.1 -0.1 1.1]);
%     hold on
%     %Time -0.5 to record value at half time interval
%     t = (n)*continuum_delta_t*continuum_tplot;
%     analy = couette_analytical_fn(t,Re,[1,0],ly,ny-1,'top');
%     plot(xaxis,analy,'r');
% 	legend ('CFD simulation','Analytical','location','NorthWest'); legend('boxoff')
% 	xlabel('y/H'); ylabel('U_x/U')
% end

% %Plot contour plot as it evolves in time
% figure
% for i=1:Ncontinuum_records
%     imagesc(continuum_velbins(:,:,i)')
%     pause(0.01)
% end
% 
% %Plot slicomatic of velocity array
% sliceomatic(continuum_velbins)

