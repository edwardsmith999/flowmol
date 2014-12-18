%Plot all output files from simulation
close all
clear all

%Find results files
resultfile_dir = './../results';

%Read Header file
read_header
read_continuum_header

%Establish and Plot domain set-up
Domain_setup

%Read continuum output
if (continuum_vflag == 3)
    read_continuum_vbins
else
    read_continuum_vslice
end

%Convert to stress
continuum_stress=continuum_vslice_to_stress(continuum_velslice,meu,ly,ny);
Nvflux_records = Nsteps / (tplot * Nvflux_ave);
t_ave = 100;

%setup stress direction and face
ixyz = 1;
jxyz = 2;

%Plot Stress
overlap = 3;
coupleddomain = ly + globaldomain(2) - 2*overlap*binsize(2);
coupledliquiddomain = coupleddomain - wallbot(2);
couplednbins = (globalnbins(ixyz)-1)/2+ny-overlap;
xaxis = 0.05:0.1:1.05;
xaxis2 = 1/nbins(2):1/nbins(2):1;
spectral_res = 6; %Ratio of points for spectral analytical solution to domain bins
analy_points = spectral_res*(couplednbins); % Number of spectral points
xaxis_analy =  0:1/(analy_points):1;
u_wall = 1;
%xloc_MDCFDhalocell = 0.5*(xaxis(couplednbins-ny) + xaxis(couplednbins-ny+1));
timeratio = delta_t/continuum_delta_t;
   
%Check CV are satisfied
for m =1:1:Nvflux_records-1
    [velocity_snapshot(:,:,:,:), ...
     velocity_flux(:,:,:,:,:),   ... 
     pressure_surface(:,:,:,:,:)    ] = read_vflux(m,resultfile_dir,globalnbins,nd);

     [velocity_snapshot_tplus1(:,:,:,:)] = read_vflux(m+1,resultfile_dir,globalnbins,nd);
 
     % %Calculate total CV flux and change in mass
     totalflux =((velocity_flux(:,:,:,:,1)+velocity_flux(:,:,:,:,4)))/(binsize(1)) ...
               +((velocity_flux(:,:,:,:,2)+velocity_flux(:,:,:,:,5)))/(binsize(2)) ...
               +((velocity_flux(:,:,:,:,3)+velocity_flux(:,:,:,:,6)))/(binsize(3));
     totalpressure =((pressure_surface(:,:,:,:,1)-pressure_surface(:,:,:,:,4)))/(binsize(1)) ...
                   +((pressure_surface(:,:,:,:,2)-pressure_surface(:,:,:,:,5)))/(binsize(2)) ...
                   +((pressure_surface(:,:,:,:,3)-pressure_surface(:,:,:,:,6)))/(binsize(3));
                  
     %totalpressure = totalpressure*delta_t
     dvelocitydt(:,:,:,:) =  (velocity_snapshot_tplus1(:,:,:,:) - velocity_snapshot(:,:,:,:))/(delta_t*Nvflux_ave);
     
     sliceomatic(( squeeze(totalpressure(:,:,:,ixyz)) ...
                  -squeeze(totalflux(:,:,:,ixyz))     ...
                  -squeeze(dvelocitydt(:,:,:,ixyz)))*10000)
     pause()
end



