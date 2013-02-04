%Plot all output files from simulation
close all
clear all

%Select diffusive solver (0) or full DNS (1)
CFD = 1;

%Wall Normal Direction
wall_normal = 2;

%Find results files
resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/flekkoy/Inc_specular_walls_large/';
resultfile_dir_MD = strcat(resultfile_dir,'md_data/results/');
resultfile_dir_CFD = strcat(resultfile_dir,'couette_data/');
resultfile_dir_CPL = strcat(resultfile_dir,'results/');

%Create directory for volume data files
system(['mkdir -p ',resultfile_dir_CPL, 'vol_data'])

% = = = Read MD header = = =
resultfile_dir = resultfile_dir_MD;
read_header

%Get interval which correspond to VMD data
Max_interval = 20;
for i=1:Max_interval
    if (i < 10)
        var_start = ['vmd_start_0', int2str(i)];
        var_end   = ['vmd_end_0', int2str(i)];
    else
        var_start = ['vmd_start_', int2str(i)];
        var_end   = ['vmd_end_', int2str(i)];
    end
    if (exist(var_start,'var'))
        interval(1,i) = eval(var_start);
        interval(2,i) = eval(var_end);
    end
end

%Write interval to file for vmd to read
fid = fopen(strcat(resultfile_dir_CPL,'vol_data/vmd_intervals'),'w');
fprintf(fid, '%u \n', interval)
fclose(fid);

%Intervals divided by velocity output frequncy
interval = floor(interval/ (tplot * Nvel_ave));

% = = = Read MD data = = =
%Check output flags and read data accordingly
if (velocity_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
elseif ((velocity_outflag > 0) && (velocity_outflag < 4) )
    read_vslice
elseif (mass_outflag == 4)
    mass_bins = read_mbins('mbins',resultfile_dir_MD);
elseif ((mass_outflag > 0) && (mass_outflag < 4) )
    read_mslice
end
%%MD domain set-up
Domain_setup
MD_domain = globaldomain;
MD_liquiddomain = liquiddomain;
MD_gnbins = gnbins;
MD_binsize = binsize;
clear liquiddomain gnbins binsize
%plot_domain

% = = = Read CFD data = = =
if (CFD == 0)
    read_continuum_header
    if (continuum_vflag == 3)
        read_continuum_vbins
    else
        read_continuum_vslice
    end
elseif(CFD == 1)
    %---Get CFD grid size ----
    [ngx, ngy, ngz, Lx, Ly, Lz, dx, dy, dz] = read_report(strcat(resultfile_dir_CFD,'report'));
    if (velocity_outflag == 4)
        %Get data step by step
    elseif ((velocity_outflag > 0) && (velocity_outflag < 4) )
        [u,v,w] = Read_DNS_slice('grid.data',resultfile_dir_CFD,ngx-2,ngy-1,ngz-2,Lx,Ly,Lz,3,true);
        continuum_velslice = u;
    end
end

%Continuum_Domain_setup
%input_CFD_gnbins = [ngx,ngy,ngz];
%read_grid(strcat(resultfile_dir_CFD,'grid.data'),[1 1 1])
CFD_domain = [Lx,Ly,Lz];
CFD_binsize= [dx,dy,dz];
CFD_gnbins = [ngx,ngy,ngz]-1;
clear Lx Ly Lz dx dy dz ngx ngy ngz

%Read Coupled data
resultfile_dir = resultfile_dir_CPL;
CPL_read_header
CPL_olap_nbins = [icmax_olap-icmin_olap+1, jcmax_olap-jcmin_olap+1,kcmax_olap-kcmin_olap+1];
CPL_olap_domain= CFD_binsize.*CPL_olap_nbins;

%Calculate coupler properties
if (CFD == 0)
    overlap = 3;
    coupleddomain(2) = ly + globaldomain(2) - 2*overlap*binsize(2);
    coupledliquiddomain = coupleddomain - wallbot(2);
    couplednbins = (gnbins(ixyz)-1)/2+ny-overlap;
    xaxis = -0.05:1/17:1.05; %binsize(2)/coupledliquiddomain:1/(couplednbins):1-binsize(2)/coupledliquiddomain;
    xaxis2 = -0.025:0.05:0.975; %-1/(4*(couplednbins-1)):1/(2*(couplednbins-1)):1-1/(4*(couplednbins-1));
    
    %xloc_MDCFDhalocell = 0.5*(xaxis(couplednbins-ny) + xaxis(couplednbins-ny+1));
    timeratio = delta_t/continuum_delta_t;
elseif(CFD == 1)
    %Get cell ratio
    MD_cells_per_CFD = CFD_binsize./MD_binsize;
    coupleddomain = CFD_domain + MD_liquiddomain - CPL_olap_domain;
    %Check domain size agrees
    domainratio  = CFD_domain./MD_domain;
    CFD_nbinswallbot= round(nbinswallbot./MD_cells_per_CFD);
end

%Get orthogonal direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    [null,ixyz]=min(domainratio);
    jxyz = mod(ixyz+1,3)+1;
    kxyz = mod(ixyz,3)+1;
end

t_ave = 10;
m = 1; %Initial Timestep
out_no = 1; %Number of output
interval_num = 1; inside_interval = 0; % Interval number
for i = 1:Nvel_records
    i

    if (m >= interval(1,interval_num) && m <= interval(2,interval_num))
        just_inside_interval = 1;
        % =================================
        %    Plot velocity
        % =================================
        %Find location of smallest dimension and take slice through this
        ixyz = wall_normal;
        jxyz = mod(ixyz+1,3)+1;
        kxyz = mod(ixyz,3)+1;
        %Read MD velocity
        if (velocity_outflag == 4)
            mass_bins = read_mbins('/mbins',resultfile_dir_MD,m-2); %NOTE WE NEED minus 2 here not sure why yet!
            vel_bins = read_vbins('/vbins',resultfile_dir_MD,m-2); %NOTE WE NEED minus 2 here not sure why yet!
            [xi,yi,zi] = meshgrid(1:size(mass_bins,2), ...
                                  1:size(mass_bins,1), ...
                                  1:size(mass_bins,3));
            [xrange,yrange,zrange] = meshgrid(1:MD_cells_per_CFD(2):size(mass_bins,2), ...
                                              1:MD_cells_per_CFD(1):size(mass_bins,1), ... 
                                              1:MD_cells_per_CFD(3):size(mass_bins,3));
            ave_mass_bins = interp3(xi,yi,zi,mass_bins,xrange,yrange,zrange);
            ave_vel_bins  = interp3(vel_bins(:,:,:,1), xrange,yrange,zrange);
            ave_vel_bins(:,:,:) = ave_vel_bins(:,:,:)./ave_mass_bins;

            %mass_slice = squeeze(sum(sum(mass_bins(:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
            %vel_slice =  squeeze(sum(sum(vel_bins(:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
            %Read CFD velocity
            [u,v,w] = Read_DNS_Subdomain(m,'grid.data',resultfile_dir_CFD, ...
                                    CFD_gnbins(1)-1,CFD_gnbins(2),CFD_gnbins(3)-1, ...
                                    1,1,1,3,true);
            u = permute(u,[2,3,1]);
            %Combine both coupled arrays
            halo =1;
            CPL_bins(:,:,:) = cat(2, ave_vel_bins(:,:,:) ...
                                   , u(:,1+CPL_olap_nbins(2)+halo:size(u,2),:) );


        elseif ((velocity_outflag > 0) & (velocity_outflag < 4) )
            ixyz = velocity_outflag;
            cumulative_m = 0; cumulative_v = 0;
            for ave = ceil(-floor(t_ave/2)):floor(t_ave/2)
                %Cumulative
                if (m+ave-2 < 1) 
                    read = 1;
                else
                    read = m;
                end
                mass_slice = read_mslice('mslice',resultfile_dir,read);
                vel_slice  = read_vslice('vslice',resultfile_dir,read);
                cumulative_m = cumulative_m + mass_slice;  
                cumulative_v = cumulative_v + vel_slice;  
            end
            mass_slice = cumulative_m/t_ave;
            vel_slice  = cumulative_v/t_ave;
            %Average multiple MD cells into one
            ave_mass_slice = interp1(mass_slice(:),1:MD_cells_per_CFD(ixyz):size(mass_slice(:)));
            %Average velocity per molecule
            for n =1:3
                ave_vel_slice(:,n) = ave_vel_slice(:,n)./ave_mass_slice;
                vel_slice(:,n) = vel_slice(:,n)./mass_slice;
            end

        end
 
        %Pad slice array to bins and get domain size vectors
        if ((velocity_outflag > 0) & (velocity_outflag < 4) )
            xdx = globaldomain(1); %linspace(-globaldomain(1)/2,globaldomain(1)/2,1) ;
            ydx = linspace(-globaldomain(2)/2,globaldomain(2)/2,size(ave_vel_slice,1)+1);
            zdx = globaldomain(3); %linspace(-globaldomain(3)/2,globaldomain(3)/2,1) ;
            vel_bins = reshape(squeeze(ave_vel_slice(:,1)),[1,size(ave_vel_slice,1),1]);
            vel_bins(:,end+1,:) = maxv;
        else
            xdx = linspace(0,coupleddomain(1)+wallbot(1),size(ave_vel_bins,1)+size(u,1)-CPL_olap_nbins(1)) ;
            ydx = linspace(0,coupleddomain(2)+wallbot(2),size(ave_vel_bins,2)+size(u,2)-CPL_olap_nbins(2)-halo);
            zdx = linspace(0,coupleddomain(3)+wallbot(3),size(ave_vel_bins,3)+size(u,3)-CPL_olap_nbins(3)) ;
        end
        G3DtoDX(xdx,ydx,zdx,CPL_bins(:,:,:,1),strcat(resultfile_dir_CPL,'./vol_data/CPL',num2str(out_no),'.dx'),-globaldomain(1)/2,-globaldomain(2)/2,-globaldomain(3)/2)
        out_no = out_no +1;

    elseif (just_inside_interval == 1)
        interval_num = interval_num + 1;
        just_inside_interval = 0;
    end

    m = m +1;

end
