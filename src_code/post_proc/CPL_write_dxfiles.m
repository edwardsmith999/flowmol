%Plot all output files from simulation
close all
clear all

system('mkdir -p ./vol_data')

%Wall Normal Direction
wall_normal = 2;

%Find results files
resultfile_dir = '/home/es205/results/MD_continuum_results/results/coupled_couette/flekkoy/MD_only_testing';
cd(resultfile_dir)

% = = = Read MD header = = =
read_header
%Check output flags and read data accordingly
if (velocity_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
elseif ((velocity_outflag > 0) & (velocity_outflag < 4) )
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nvel_ave));
    %Get maximum value to normalise everything by
    %Read all slice values
    v=read_vslice('vslice',resultfile_dir);
    m=read_mslice('mslice',resultfile_dir);
    maxv = max(max(squeeze(v(:,1,:))./m));
elseif (mass_outflag == 4)
    Nvel_records = floor((Nsteps-initialstep) / (tplot * Nmass_ave));
    mass_bins = read_mbins('mbins',resultfile_dir_md);
elseif ((mass_outflag > 0) & (mass_outflag < 4) )
    read_mslice
end

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
%Intervals divided by velocity output frequncy
interval = floor(interval/ (tplot * Nvel_ave));

%%MD domain set-up
Domain_setup

%Calculate properties
MD_cells_per_ave = 1;
wallsize = zeros(1,3);
wallsize(2) = binsize(2);
MD_domain = globaldomain - wallsize;

%Get orthogonal direction
if (velocity_outflag < 4)
    ixyz = velocity_outflag;
else
    ixyz = wall_normal;
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
            vel_slice = read_vslice('./vslice',resultfile_dir);
            filename = strcat(resultfile_dir,'/mbins');
            mass_bins = read_mbins(filename,resultfile_dir,m-2); %NOTE WE NEED minus 2 here not sure why yet!
            filename = strcat(resultfile_dir,'/vbins');
            vel_bins = read_vbins(filename,resultfile_dir,m-2); %NOTE WE NEED minus 2 here not sure why yet!
            mass_slice = squeeze(sum(sum(mass_bins(:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
            vel_slice =  squeeze(sum(sum(vel_bins(:,:,:,:),jxyz),kxyz)/(nbins(jxyz)*nbins(kxyz)));
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
        end

        %Average multiple MD cells into one
        switch MD_cells_per_ave
            case 1
                ave_vel_slice = vel_slice;
                ave_mass_slice = mass_slice;
            case 2
                ave_vel_slice = zeros(size(vel_slice,1)/MD_cells_per_ave,size(vel_slice,2));
                ave_mass_slice = zeros(size(mass_slice,2)/MD_cells_per_ave,size(mass_slice,1));
                n=1;
                for icell=1:MD_cells_per_ave:size(vel_slice,1)
                    ave_vel_slice(n,:)  = 0.5*(vel_slice(icell,:) + vel_slice(icell+1,:));
                    ave_mass_slice(n)   = 0.5*(mass_slice(icell)  + mass_slice(icell+1));
                    n = n + 1;
                end
            case 4
                ave_vel_slice = zeros(size(vel_slice,1)/MD_cells_per_ave,size(vel_slice,2));
                ave_mass_slice = zeros(size(mass_slice,2)/MD_cells_per_ave,size(mass_slice,1));
                n=1;
                for icell=1:MD_cells_per_ave:size(vel_slice,1)
                    ave_vel_slice(n,:) = 0.25*( vel_slice(icell  ,:) ...
                        +vel_slice(icell+1,:) ...
                        +vel_slice(icell+2,:) ...
                        +vel_slice(icell+3,:));
                    ave_mass_slice(n) = 0.25*(mass_slice(icell  ) ...
                        +mass_slice(icell+1) ...
                        +mass_slice(icell+2) ...
                        +mass_slice(icell+3));
                    n = n + 1;
                end
        end

        %Average velocity per molecule
        for n =1:3
            ave_vel_slice(:,n) = ave_vel_slice(:,n)./ave_mass_slice;
            vel_slice(:,n) = vel_slice(:,n)./mass_slice;
        end

        if ((velocity_outflag > 0) & (velocity_outflag < 4) )
            xdx = globaldomain(1); %linspace(-globaldomain(1)/2,globaldomain(1)/2,1) ;
            ydx = linspace(-globaldomain(2)/2,globaldomain(2)/2,size(ave_vel_slice,1)+1);
            zdx = globaldomain(3); %linspace(-globaldomain(3)/2,globaldomain(3)/2,1) ;
            vel_bins = reshape(squeeze(ave_vel_slice(:,1)),[1,size(ave_vel_slice,1),1]);
            vel_bins(:,end+1,:) = maxv;
        end
        G3DtoDX(xdx,ydx,zdx,vel_bins,strcat('./vol_data/MD',num2str(out_no),'.dx'),-globaldomain(1)/2,-globaldomain(2)/2,-globaldomain(3)/2)
        out_no = out_no +1;

    elseif (just_inside_interval == 1)
        interval_num = interval_num + 1;
        just_inside_interval = 0;
    end

    m = m +1;

end
