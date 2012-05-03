%Plots domain evolution in times as contour, isosurfaces or 2d plots

clear all
close all

%Output = 1 for contour, 2 for isosurfaces or 3 for 2d plots
output = 3;

%Contour plots
switch(output)
case(1) 
    read_header
    [vel_bins]=read_vbins('./vbins');

    for i =1:size(vel_bins,5)-1
        subplot(2,1,1)
            h=slice(vel_bins(:,:,:,1,i),[],[],5,'cubic');
            set(h,'FaceColor','interp',...
              'EdgeColor','none',...
              'DiffuseStrength',.8)
            axis fill 
            view(90,90)
            colorbar; caxis([min(vel_bins(:)),500])
        subplot(2,1,2)
            [velocity_snapshot_t,velocity_flux_t,pressure_surface_t] = ...
                read_vflux(i,'./../../results',gnbins,nd);
            h=slice(pressure_surface_t(:,:,:,1,1),[],[],5,'cubic');
            set(h,'FaceColor','interp',...
              'EdgeColor','none',...
              'DiffuseStrength',.8)
            view(90,90)
            colorbar;caxis([-20,60])
            pause(0.1)
            drawnow

    end

%3D isosurfaces
case(2)
    read_header
    [vel_bins]=read_vbins('./vbins');
    [x,y,z] = meshgrid(1:9,1:120,1:9);

    for i =1:size(vel_bins,5)-1
        subplot(2,1,1)
             [v] = vel_bins(:,:,:,1,i);
             [faces,verts,colors] = isosurface(x,y,z,v,200,x); 
             patch('Vertices', verts, 'Faces', faces, ... 
                 'FaceVertexCData', colors, ... 
                 'FaceColor','interp', ... 
                 'edgecolor', 'interp');
            view(45,45)
            colorbar; caxis([min(vel_bins(:)),500])
        subplot(2,1,2)
            [velocity_snapshot_t,velocity_flux_t,pressure_surface_t] = ...
                read_vflux(i,'./../../results',gnbins,nd);
            [v] = pressure_surface_t(:,:,:,1,1);
            [faces,verts,colors] = isosurface(x,y,z,v,15,x); 
            patch('Vertices', verts, 'Faces', faces, ... 
                 'FaceVertexCData', colors, ... 
                 'FaceColor','interp', ... 
                 'edgecolor', 'interp');
            view(45,45)
            colorbar;caxis([-20,60])
            pause(0.1)
            drawnow

    end

%Averaged Plots
case(3)
    read_header
    [mass_bins,density_bins] = read_mbins('./mbins');
    vel_bins=read_vbins('./vbins');
    T_bins = read_Tbins('./Tbins');

    for i =1:size(vel_bins,5)-1
            subplot(4,1,1)
                plot(mean(mean(density_bins(:,:,:,i),2),3))

            subplot(4,1,2)
                plot(mean(mean(vel_bins(:,:,:,1,i),2),3))

            subplot(4,1,3)
                [velocity_snapshot_t,velocity_flux_t,pressure_surface_t] = ...
                    read_vflux(i,'./../../results',gnbins,nd);
                plot(mean(mean(pressure_surface_t(:,:,:,1,1),2),3))

            subplot(4,1,4)
                plot(mean(mean(T_bins(:,:,:,i),2),3))
            pause(0.1)
            drawnow

            %[energy_snapshot_t,energy_flux_t,energy_surface_t] = ...
            %    read_eflux(i,'./../../results',gnbins)

           % subplot(2,2,2)
           %     plot(mean(mean(energy_surface_t(:,:,:,1),2),3))
    end

end
