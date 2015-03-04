%Plot continuum velocity slice vs analytical solution
clear all
close all

%Find results files
resultfile_dir = './../../results';

%Read data for continuum plot
read_continuum_header

scrsz = get(0,'ScreenSize');
fig1 = figure('Position',[1 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);

u_wall = 2;%0.5*(continuum_velslice(ny+1,1) + continuum_velslice(ny+2,1));

spectral_res = 1;
xaxis = 1/(ny+1):1/(ny+1):1-1/(ny+1);
xaxis = -1/(ny+2):1/(ny):1+1/(ny+2);
xaxis_analy = 0:1/(spectral_res*ny):1;
m = 1;

fid = fopen('./coute_strs_anayl','w','n');
for i = 1:1000  %Ncontinuum_records
    
    i    
    %Time -0.5 to record value at half time interval
    t = (i)*continuum_delta_t*continuum_tplot;
    
    analy=couette_analytical_stress_fn(t,Re,1,ly,spectral_res*ny-1);
    fwrite(fid,analy,'double');
  
    %plot(analy)
    %axis([0,ny+2,-0.4,0.0])
    %pause(0.5)
    
end
