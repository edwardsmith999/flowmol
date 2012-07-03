
clear all
close all
time = 10;

runs = 3;
names{1} = '/home/es205/codes/coupled/MD_dCSE/src_code/results/1bin';
names{2} =  '/home/es205/codes/coupled/MD_dCSE/src_code/results/4bins';
names{3} = '/home/es205/codes/coupled/MD_dCSE/src_code/results/4bin_parallel';

for o = 1:runs
    resultfile_dir = names{o};
    read_header

    %Test types integer flag
    % mass_bins     = 1
    % mass_flux     = 2
    % momentum_flux = 3

    test_type = 4;
    switch(test_type)
    case(1) 
        %First iteration set up arrays
        if (o==1) 
            mass_bins_sum = zeros(globalncells(1),globalncells(2),globalncells(3),(Nsteps-initialstep)/(Nmass_ave*tplot),runs); 
        end
        bpc = gnbins./globalncells;
        filename = './mbins';
        [mass_bins,density_bins]=read_mbins(filename,resultfile_dir);

        %Avergae multiple bins to give same as number of cells
        l = 1; m = 1; n = 1; 
        for i=1:bpc(1):size(mass_bins,1)
        for j=1:bpc(2):size(mass_bins,2)
        for k=1:bpc(3):size(mass_bins,3)
            for p = 0:bpc(1)-1
            for q = 0:bpc(2)-1
            for r = 0:bpc(3)-1
                mass_bins_sum(l,m,n,:,o) = mass_bins_sum(l,m,n,:,o) + mass_bins(i+p,j+q,k+r,:);
            end
            end
            end
            n = n + 1;
        end
            n = 1;
            m = m + 1;
        end
            m = 1;
            l = l + 1;
        end
        %Check error between current and previous runs
        if (o == 2)
            error = mass_bins_sum(:,:,:,:,o) - mass_bins_sum(:,:,:,:,o-1);
        elseif (o > 2)
            error = error + mass_bins_sum(:,:,:,:,o) - mass_bins_sum(:,:,:,:,o-1);
        end

    case(2)

        resultfile_dir = names{o};
        read_header
        bpc = gnbins./globalncells;
        nhb = bpc;

        filename1 = './mflux';
        filename2 = './msnap';
        [mass_flux,mass_snapshot,Nmflux_records]=read_mflux(filename1,filename2,resultfile_dir);
        read_header

        if (o==1) 
            mass_flux_sum = zeros(size(mass_flux,1)/bpc(1),size(mass_flux,2)/bpc(2),size(mass_flux,3)/bpc(3),6,(Nsteps-initialstep)/(Nmflux_ave*tplot),runs);        
            mass_snap_sum = zeros(globalncells(1),globalncells(2),globalncells(3),(Nsteps-initialstep)/(Nmflux_ave*tplot)+1,runs); 
        end
        l = 1; m = 1; n = 1;
        for i=1:bpc(1):size(mass_flux,1)
        for j=1:bpc(2):size(mass_flux,2)
        for k=1:bpc(3):size(mass_flux,3)
            for p = 0:bpc(1)-1
            for q = 0:bpc(2)-1
            for r = 0:bpc(3)-1
                mass_flux_sum(l,m,n,:,:,o) = mass_flux_sum(l,m,n,:,:,o) + mass_flux(i+p,j+q,k+r,:,:);
                mass_snap_sum(l,m,n,:,o)   = mass_snap_sum(l,m,n,:,o)   + mass_snapshot(i+p,j+q,k+r,:);
            end
            end
            end
            n = n + 1;
        end
            n = 1;
            m = m + 1;
        end
            m = 1;
            l = l + 1;
        end
        %sliceomatic(sum(mass_flux_sum(:,:,:,:,time,o),4))
        %sliceomatic(sum(mass_snap_sum(:,:,:,time,o),4))
        %Check error between current and previous runs
        if (o == 2)
            error = mass_flux_sum(:,:,:,:,:,o) - mass_flux_sum(:,:,:,:,:,o-1);
            error2= mass_snap_sum(:,:,:,:,o) - mass_snap_sum(:,:,:,:,o-1);
        elseif (o > 2)
            error = error + mass_flux_sum(:,:,:,:,:,o) - mass_flux_sum(:,:,:,:,:,o-1);
            error2= error2 + mass_snap_sum(:,:,:,:,o) - mass_snap_sum(:,:,:,:,o-1);
        end

        figure
        plot(0:size(mass_flux,1)-1, sum(sum(sum(mass_flux(:,:,:,1,:),2),3),5),'s')
        hold all
        plot(1:size(mass_flux,1)  ,-sum(sum(sum(mass_flux(:,:,:,4,:),2),3),5),'x')
        hold off

        %sliceomatic(sum(mass_flux(2:end,:,:,1,:),5)+sum(mass_flux(1:end-1,:,:,4,:),5))

    case(3)

        resultfile_dir = names{o};
        read_header
        bpc = gnbins./globalncells;
        nhb = bpc;
        surface_plot            = zeros(gnbins(1),runs);

        filename1 = './vflux';
        filename2 = './psurface';
        filename3 = './vsnap';

        for t = 1:(Nsteps-initialstep)/(Nvflux_ave*tplot)-1
            t
            [velocity_snapshot,velocity_flux,pressure_surface] = read_vflux(t,resultfile_dir,gnbins,nd,filename1,filename2,filename3);
            read_header

            if (o==1) 
                velocity_flux_sum       = zeros(globalncells(1),globalncells(2),globalncells(3),3,6,(Nsteps-initialstep)/(Nvflux_ave*tplot),runs);        
                velocity_snap_sum       = zeros(globalncells(1),globalncells(2),globalncells(3),3,(Nsteps-initialstep)/(Nvflux_ave*tplot)+1,runs); 
                pressure_surface_sum    = zeros(globalncells(1),globalncells(2),globalncells(3),3,6,(Nsteps-initialstep)/(Nvflux_ave*tplot),runs);
            end

            l = 1; m = 1; n = 1;
            for i=1:bpc(1):size(velocity_flux,1)
            for j=1:bpc(2):size(velocity_flux,2)
            for k=1:bpc(3):size(velocity_flux,3)
                for p = 0:bpc(1)-1
                for q = 0:bpc(2)-1
                for r = 0:bpc(3)-1
                    velocity_flux_sum(l,m,n,:,:,t,o)      = velocity_flux_sum(l,m,n,:,:,t,o) + velocity_flux(i+p,j+q,k+r,:,:);
                    velocity_snap_sum(l,m,n,:,t,o)        = velocity_snap_sum(l,m,n,:,t,o)   + velocity_snapshot(i+p,j+q,k+r,:);
                    pressure_surface_sum(l,m,n,:,:,t,o)   = pressure_surface_sum(l,m,n,:,:,t,o) + pressure_surface(i+p,j+q,k+r,:,:);
         
                end
                end
                end
                n = n + 1;
            end
                n = 1;
                m = m + 1;
            end
                m = 1;
                l = l + 1;
            end

            surface_plot(:,1) = surface_plot(:,1)  ...
                                + squeeze(sum(sum(   velocity_flux(:,:,:,1,1),2),3)) ...
                                + squeeze(sum(sum(pressure_surface(:,:,:,1,1),2),3));
            surface_plot(:,2) = surface_plot(:,2)  ...
                                - squeeze(sum(sum(   velocity_flux(:,:,:,1,4),2),3)) ...
                                + squeeze(sum(sum(pressure_surface(:,:,:,1,4),2),3));
        end

        figure
        plot(0:size(surface_plot,1)-1, surface_plot(:,1),'s')
        hold all
        plot(1:size(surface_plot,1)  , surface_plot(:,2),'x')
        hold off

        %sliceomatic(sum(velocity_flux_sum(:,:,:,:,time,o),4))
        %sliceomatic(sum(velocity_snap_sum(:,:,:,time,o),4))
        %Check error between current and previous runs
        if (o == 2)
            error =    velocity_flux_sum(:,:,:,:,:,:,o) -    velocity_flux_sum(:,:,:,:,:,:,o-1);
            error2=    velocity_snap_sum(:,:,:,:,:,o)   -    velocity_snap_sum(:,:,:,:,:,o-1);
            error3= pressure_surface_sum(:,:,:,:,:,:,o) - pressure_surface_sum(:,:,:,:,:,:,o-1);
        elseif (o > 2)
            error = error  +    velocity_flux_sum(:,:,:,:,:,:,o) -    velocity_flux_sum(:,:,:,:,:,:,o-1);
            error2= error2 +    velocity_snap_sum(:,:,:,:,:,  o) -    velocity_snap_sum(:,:,:,:,:,  o-1);
            error3= error3 + pressure_surface_sum(:,:,:,:,:,:,o) - pressure_surface_sum(:,:,:,:,:,:,o-1);
        end

    case(4)

        resultfile_dir = names{o};
        read_header
        bpc = gnbins./globalncells;
        nhb = bpc;
        surface_plot            = zeros(gnbins(1),runs);

        filename1 = './eflux';
        filename2 = './esurface';
        filename3 = './esnap';

        for t = 3:(Nsteps-initialstep)/(Neflux_ave*tplot)-1
            t
            [energy_snapshot,energy_flux,energy_surface] = read_eflux(t,resultfile_dir,gnbins,filename1,filename2,filename3);
            read_header

            if (o==1) 
                energy_flux_sum       = zeros(globalncells(1),globalncells(2),globalncells(3),6,(Nsteps-initialstep)/(Nvflux_ave*tplot),runs);        
                energy_snap_sum       = zeros(globalncells(1),globalncells(2),globalncells(3),(Nsteps-initialstep)/(Nvflux_ave*tplot)+1,runs); 
                energy_surface_sum    = zeros(globalncells(1),globalncells(2),globalncells(3),6,(Nsteps-initialstep)/(Nvflux_ave*tplot),runs);
            end

            l = 1; m = 1; n = 1;
            for i=1:bpc(1):size(energy_flux,1)
            for j=1:bpc(2):size(energy_flux,2)
            for k=1:bpc(3):size(energy_flux,3)
                for p = 0:bpc(1)-1
                for q = 0:bpc(2)-1
                for r = 0:bpc(3)-1
                    energy_flux_sum(l,m,n,:,t,o)      = energy_flux_sum(l,m,n,:,t,o) + energy_flux(i+p,j+q,k+r,:);
                    energy_snap_sum(l,m,n,t,o)        = energy_snap_sum(l,m,n,t,o)   + energy_snapshot(i+p,j+q,k+r);
                    energy_surface_sum(l,m,n,:,t,o)   = energy_surface_sum(l,m,n,:,t,o) + energy_surface(i+p,j+q,k+r,:);
         
                end
                end
                end
                n = n + 1;
            end
                n = 1;
                m = m + 1;
            end
                m = 1;
                l = l + 1;
            end

            surface_plot(:,1) = surface_plot(:,1)  ...
                                + squeeze(sum(sum(   energy_flux(:,4,4,1),2),3)) ...
                                + squeeze(sum(sum(energy_surface(:,4,4,1),2),3));
            surface_plot(:,2) = surface_plot(:,2)  ...
                                - squeeze(sum(sum(   energy_flux(:,4,4,4),2),3)) ...
                                + squeeze(sum(sum(energy_surface(:,4,4,4),2),3));
        end

        figure
        plot(0:size(surface_plot,1)-1, surface_plot(:,1),'s')
        hold all
        plot(1:size(surface_plot,1)  , surface_plot(:,2),'x')
        hold off

        %sliceomatic(sum(energy_flux_sum(:,:,:,:,time,o),4))
        %sliceomatic(sum(energy_snap_sum(:,:,:,time,o),4))
        %Check error between current and previous runs
        if (o == 2)
            error =    energy_flux_sum(:,:,:,:,:,o) -    energy_flux_sum(:,:,:,:,:,o-1);
            error2=    energy_snap_sum(:,:,:,:,o)   -    energy_snap_sum(:,:,:,:,o-1);
            error3= energy_surface_sum(:,:,:,:,:,o) - energy_surface_sum(:,:,:,:,:,o-1);
        elseif (o > 2)
            error = error  +    energy_flux_sum(:,:,:,:,:,o) -    energy_flux_sum(:,:,:,:,:,o-1);
            error2= error2 +    energy_snap_sum(:,:,:,:,  o) -    energy_snap_sum(:,:,:,:,  o-1);
            error3= error3 + energy_surface_sum(:,:,:,:,:,o) - energy_surface_sum(:,:,:,:,:,o-1);
        end

    end

end

disp(strcat('Error is cell refinements over ', num2str(runs), ' runs = ', ...
    num2str(sum(error(:))+mean(error(:))+max(error(:))+min(error(:)))));

disp(strcat('Error is cell refinements over ', num2str(runs), ' runs = ', ...
    num2str(sum(error2(:))+mean(error2(:))+max(error2(:))+min(error2(:)))));

disp(strcat('Error is cell refinements over ', num2str(runs), ' runs = ', ...
    num2str(sum(error3(:))+mean(error3(:))+max(error3(:))+min(error3(:)))));

if (sum(error(:))+mean(error(:))+max(error(:))+min(error(:)) > 0)
    sliceomatic(sum(error,4))
end