%Plots domain evolution in times as contour, isosurfaces or 2d plots

clear all
close all

%Output = 1 for contour, 2 for isosurfaces or 3 for 2d plots
output = 3;

%Contour plots
switch(output)
case(1) 

    %Write a sequence of frames to a compressed AVI file, couette.avi:
    % Prepare the new file.
    vidObj = VideoWriter('shockcontour.avi');
    vidObj.FrameRate = 2;
    open(vidObj);

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
                read_vflux(i,'./../results',gnbins,nd);
            h=slice(pressure_surface_t(:,:,:,1,1),[],[],5,'cubic');
            set(h,'FaceColor','interp',...
              'EdgeColor','none',...
              'DiffuseStrength',.8)
            view(90,90)
            colorbar;caxis([-20,60])
            pause(0.1)
            drawnow

            %Store pictures and videos
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);

    end

    % Close the file.
    close(vidObj)

%3D isosurfaces
case(2)

    read_header
    [vel_bins]=read_vbins('./vbins');
    [x,y,z] = meshgrid(1:9,1:gnbins(1),1:9);

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
                read_vflux(i,'./../results',gnbins,nd);
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

        scrsz = get(0,'ScreenSize');
        h=figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)]);
		daspect([1 1 1]) %Fix aspect ratio

        %Write a sequence of frames to a compressed AVI file, couette.avi:
        % Prepare the new file.
        vidObj = VideoWriter('./shockplot.avi');
        vidObj.FrameRate = 2;
        open(vidObj);

        %resultfile_dir = './../results';
        resultfile_dir = '/home/es205/results/md_results/fortran/3D_code/parallel/code/shock_wave/MD_dCSE/src_code/results/ensemble';

        read_header
        [mass_bins,density_bins] = read_mbins('./mbins',resultfile_dir);
        [mass_flux,mass_snapshot,Nmflux_records] = read_mflux;
        vel_bins=read_vbins('./vbins',resultfile_dir);
        v_bins(:,:,:,1,:) = squeeze(vel_bins(:,:,:,1,:));%./(mass_bins(:,:,:,:)+1);
        %T_bins = read_Tbins('./Tbins');


        %Scaling and tings
        binvol = prod(binsize);
        sigmatoA = 3.4; MDPressuretoPa = 4.214e6;
        Pa2kbar = 1/100000000;
		
        xaxis_scale = binsize(1)*sigmatoA;
        rho_0 = 0.8;

        %Averaging window sizes
        t_win_sze = 10;
        win_sze = 10;
        win_srt = 10;
        Wave_average = zeros(win_sze+1,5);
        for i =180:1:Nsteps-t_win_sze-1

            subplot(5,1,1)
                xaxis = (1:gnbins)*xaxis_scale;
                plot(xaxis,(mean(mean(mass_snapshot(:,:,:,i+1),2),3)/binvol)/rho_0)
                hold all
                %plot(mean(mean(mass_bins(:,:,:,i),2),3)/(binvol))   
                axis([0 max(xaxis) 0.8 1.5])
                xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                xlabel('x'); ylabel('\rho / \rho_0'); 
                figure = h;
                title(strcat('Iteration number ',num2str(i)))   
                hold off
     
            subplot(5,1,2)
                xaxis = (1:gnbins)*xaxis_scale;
                plot(xaxis, sum(sum(sum(mass_flux(:,:,:,:,i),2),3),4),'r')
                hold all
                dmassdt(:,:,:,i) = mass_snapshot(:,:,:,i+1) - mass_snapshot(:,:,:,i);
                plot(xaxis,sum(sum(dmassdt(:,:,:,i),2),3),'b--')
                plot(xaxis,sum(sum(mass_flux(:,:,:,1,i),2),3) ... 
                     +sum(sum(mass_flux(:,:,:,4,i),2),3),'b--')
                xlabel('x'); ylabel('d\rho /dt ','Interpreter','tex'); 
                axis([0 max(xaxis) -5 5])
                xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                hold off

            subplot(5,1,3)
                [velocity_snapshot_t] = ...
                    read_vflux(i,resultfile_dir,gnbins,nd);
                xaxis = (1:gnbins)*xaxis_scale;
                plot(xaxis,mean(mean(velocity_snapshot_t(:,:,:,1),2),3))
                hold all
                xaxis = (0.5:gnbins(1)-0.5)*xaxis_scale;
                plot(xaxis,mean(sum(sum(mass_flux(:,:,:,1,i:i+t_win_sze),2),3),5))
                %plot(mean(mean(v_bins(:,:,:,1,i),2),3))
                axis([0 max(xaxis) -0.1 max(velocity_snapshot_t(:))])
                xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                xlabel('x'); ylabel('u_x'); 
                hold off

            subplot(5,1,4)
            	[velocity_snapshot_t,velocity_flux_t,pressure_surface_t] = ...
                    read_vflux(i,resultfile_dir,gnbins,nd);
%                 pressure_surface_t(1*gnbins(1)/npx  ,:,:,:,4) = pressure_surface_t(1*gnbins(1)/npx  ,:,:,:,4)/2;
%                 pressure_surface_t(1*gnbins(1)/npx+1,:,:,:,1) = pressure_surface_t(1*gnbins(1)/npx+1,:,:,:,1)/2;
%                 pressure_surface_t(2*gnbins(1)/npx  ,:,:,:,4) = pressure_surface_t(1*gnbins(1)/npx  ,:,:,:,4)/2;
%                 pressure_surface_t(2*gnbins(1)/npx+1,:,:,:,1) = pressure_surface_t(2*gnbins(1)/npx+1,:,:,:,1)/2;
%                 pressure_surface_t(3*gnbins(1)/npx  ,:,:,:,4) = pressure_surface_t(1*gnbins(1)/npx  ,:,:,:,4)/2;
%                 pressure_surface_t(3*gnbins(1)/npx+1,:,:,:,1) = pressure_surface_t(3*gnbins(1)/npx+1,:,:,:,1)/2;
                xaxis = (0:gnbins-1)*xaxis_scale;
                plot(xaxis,mean(mean(pressure_surface_t(:,:,:,1,1),2),3)*MDPressuretoPa*Pa2kbar,'b')
                hold all
                plot(xaxis,-mean(mean(  velocity_flux_t(:,:,:,1,1),2),3)*MDPressuretoPa*Pa2kbar,'r')
                xaxis = (1:gnbins)*xaxis_scale;
                plot(xaxis,mean(mean(pressure_surface_t(:,:,:,1,4),2),3)*MDPressuretoPa*Pa2kbar,'bx')
                plot(xaxis, mean(mean(  velocity_flux_t(:,:,:,1,4),2),3)*MDPressuretoPa*Pa2kbar,'rx')
                legend('configurational','kinetic')
                axis([0 max(xaxis) min(mean(mean(pressure_surface_t(:,:,:,1,4),2),3))*MDPressuretoPa*Pa2kbar ...
                                   max(mean(mean(pressure_surface_t(:,:,:,1,4),2),3))*MDPressuretoPa*Pa2kbar    ])
                xaxis = (1.0:gnbins(1)+1)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                xlabel('x'); ylabel('\Pi'); 
                hold off

             subplot(5,1,5)
                [velocity_snapshot(:,:,:,:)] = read_vflux(i,resultfile_dir,gnbins,nd);
                velocity_flux = velocity_flux_t;
                pressure_surface = pressure_surface_t;

                [velocity_snapshot_tplus1(:,:,:,:)] = read_vflux(i+1,resultfile_dir,gnbins,nd);
                [velocity_snapshot_tplus2(:,:,:,:)] = read_vflux(i+2,resultfile_dir,gnbins,nd);
                 % %Calculate total CV flux and change in mass
                 totalflux =((velocity_flux(:,:,:,:,1)+velocity_flux(:,:,:,:,4)))/(binsize(1)) ...
                           +((velocity_flux(:,:,:,:,2)+velocity_flux(:,:,:,:,5)))/(binsize(2)) ...
                           +((velocity_flux(:,:,:,:,3)+velocity_flux(:,:,:,:,6)))/(binsize(3));
                 totalpressure =((pressure_surface(:,:,:,:,1)-pressure_surface(:,:,:,:,4)))/(binsize(1)) ...
                               +((pressure_surface(:,:,:,:,2)-pressure_surface(:,:,:,:,5)))/(binsize(2)) ...
                               +((pressure_surface(:,:,:,:,3)-pressure_surface(:,:,:,:,6)))/(binsize(3));
                 dvelocitydt(:,:,:,:) =  ( velocity_snapshot_tplus2(:,:,:,:) ...
                                          -velocity_snapshot_tplus1(:,:,:,:))/(delta_t*Nvflux_ave);
                 xaxis = (1:gnbins)*xaxis_scale;
                 plot(xaxis,delta_t*squeeze(sum(sum((totalpressure(:,:,:,1) ...
                                       -  totalflux(:,:,:,1)),2),3)),'r')
                 hold all
                 plot(xaxis,delta_t*squeeze(sum(sum(   dvelocitydt(:,:,:,1),2),3)),'b--')

                %disp(strcat(' iter= ',num2str(i), ...
                %            ' pressure= ',num2str(totalpressure(50,4,4,1)), ...
                %            ' flux= ', num2str(totalflux(50,4,4,1)), ...
                %            ' dt= ', num2str(dvelocitydt(50,4,4,1))))

                 %plot(squeeze(sum(sum((totalpressure(:,:,:,1) ...
                 %                       -  totalflux(:,:,:,1) ...
                 %                       -dvelocitydt(:,:,:,1)),2),3)),'r')
                 axis([0 max(xaxis) -7 7])
                 xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                 gridxy(xaxis,'LineStyle',':')
                 xlabel('x'); ylabel('d (\rho u_x) / dt');
                 hold off


            %Calculate lagrangian averages!
            c(:,i) = delta_t*squeeze(sum(sum((totalpressure(:,:,:,1) ...
                 -  totalflux(:,:,:,1)),2),3))./sum(sum(sum(mass_flux(:,:,:,:,i),2),3),4);

            Wave_average(:,1) = Wave_average(:,1) + mean(mean(mass_snapshot(win_srt:win_srt+win_sze,:,:,i+1),2),3);
            Wave_average(:,2) = Wave_average(:,2) + sum(sum(mass_flux(win_srt:win_srt+win_sze,:,:,1,i),2),3) + ...
                                                    sum(sum(mass_flux(win_srt:win_srt+win_sze,:,:,4,i),2),3);
            Wave_average(:,3) = Wave_average(:,3) + 4*mean(mean(velocity_snapshot_t(win_srt:win_srt+win_sze,:,:,1),2),3);
            Wave_average(:,4) = Wave_average(:,4) + mean(mean(pressure_surface_t(win_srt:win_srt+win_sze,:,:,1,4),2),3) ... 
                                                   -mean(mean(   velocity_flux_t(win_srt:win_srt+win_sze,:,:,1,4),2),3);
            %Wave_average(5,:) = Wave_average(5,:) + delta_t*squeeze(sum(sum((totalpressure(win_srt:win_srt+win_sze,:,:,1) ...
            %                                                                  -  totalflux(win_srt:win_srt+win_sze,:,:,1)),2),3));


            %hold all
%             subplot(5,1,5)
%                 plot(mean(mean((T_bins(:,:,:,i)./(nd*mass_bins(:,:,:,i))),2),3))
%                 grid on
%                 %hold all
%                 axis([0 gnbins(1) 0 0.5])
%                 xlabel('x'); ylabel('Temperature'); 

            %pause(0.8)
            whitebg(h,'w')
            drawnow
             

            %savefig(strcat('shock',num2str(i)),'png')

            %Store pictures and videos
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
            
%             [energy_snapshot_t,energy_flux_t,energy_surface_t] = ...
%                 read_eflux(i,'./../results',gnbins);
%             
%              subplot(5,1,5)
%              plot(mean(mean(energy_surface_t(:,:,:,1),2),3))
%              hold all
%              plot(mean(mean(energy_flux_t(:,:,:,1),2),3))
%              hold off
%              grid on
%              %hold all
%              %axis([0 gnbins(1) 0 1])
%              xlabel('x'); ylabel('Energy'); 
        end

    % Close the file.
    close(vidObj)

        
end
