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
        resultfile_dir = '/home/es205/results/md_results/fortran/3D_code/parallel/code/shock_wave/MD_dCSE/src_code/results/';
        %resultfile_dir = '/home/es205/results/md_results/fortran/3D_code/parallel/code/shock_wave/MD_dCSE/src_code/results/ensemble/run_of_50';

        read_header
        [mass_bins,density_bins] = read_mbins('./mbins',resultfile_dir);
        [mass_flux,mass_snapshot,Nmflux_records] = read_mflux('./mflux','./msnap',resultfile_dir);
        vel_bins=read_vbins('./vbins',resultfile_dir);
        v_bins(:,:,:,1,:) = squeeze(vel_bins(:,:,:,1,:));%./(mass_bins(:,:,:,:)+1);
        %T_bins = read_Tbins('./Tbins');


        %Scaling and tings
        binvol = prod(binsize);
        sigmatoA = 3.4; MDPressuretoPa = 4.214e6;
        Pa2kbar = 1/100000000;
		
        xaxis_scale = binsize(1)*sigmatoA;
        rho_0 = 0.707;

        %Averaging window sizes
        t_win_sze = 10;
        win_sze = 10;
        win_srt = 10;
        Wave_average = zeros(win_sze+1,5);
        for i =1:5:(Nsteps)/tplot-5
            i

            g(1)=subplot(6,1,1);
                figure = h;
                title(strcat('Iteration number ',num2str(i)))  
                if (i<= 9) 
                    filename = strcat('/home/es205/results/md_results/fortran/3D_code/parallel/code/shock_wave/MD_dCSE/src_code/results/ensemble/run_of_50/untitled.000',num2str(i),'.ppm');
                end
                if(i>=10 && i<=99) 
                   filename = strcat('/home/es205/results/md_results/fortran/3D_code/parallel/code/shock_wave/MD_dCSE/src_code/results/ensemble/run_of_50/untitled.00',num2str(i),'.ppm');
                end
                if(i>=100 && i<=999)
                   filename = strcat('/home/es205/results/md_results/fortran/3D_code/parallel/code/shock_wave/MD_dCSE/src_code/results/ensemble/run_of_50/untitled.0',num2str(i),'.ppm');
                end
                A = imread(filename);
                B = A(30:220,40:1400);
                %BW = im2bw(A(30:220,40:1400),0.01);
                %imshow(BW,'Parent', g(1))
                image(B); colormap(gray); 
                set(gca,'xtick',[])

            g(2)=subplot(6,1,2);
                xaxis = (1:gnbins)*xaxis_scale;
                plot(xaxis,(mean(mean(mass_snapshot(:,:,:,i+1),2),3)/binvol)/rho_0)
                hold all
                plot(xaxis,mean(mean(mass_bins(:,:,:,i),2),3)/(binvol)/rho_0,'--')   
                axis([0 max(xaxis) 0.8 1.5])
                xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                xlabel('x'); ylabel('\rho / \rho_0'); 
                hold off
                set(gca,'xtick',[])
            g(3)=subplot(6,1,3);
                xaxis = (1:gnbins)*xaxis_scale;
                %plot(xaxis,sum(sum(mass_flux(:,:,:,1,i),2),3) ... 
                %          +sum(sum(mass_flux(:,:,:,4,i),2),3),'r-')
                plot(xaxis, sum(sum(sum(mass_flux(:,:,:,:,i),2),3),4),'r-')
                hold all
                dmassdt(:,:,:,i) = mass_snapshot(:,:,:,i+1) - mass_snapshot(:,:,:,i);
                plot(xaxis,sum(sum(dmassdt(:,:,:,i),2),3),'b--')
                xlabel('x'); ylabel('d\rho /dt ','Interpreter','tex'); 
                axis([0 max(xaxis) -5 5])
                xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                hold off

                set(gca,'xtick',[])
            g(4)=subplot(6,1,4);
                [velocity_snapshot_t] = ...
                    read_vflux(i,resultfile_dir,gnbins,nd);
                xaxis = (1:gnbins)*xaxis_scale;
                plot(xaxis,mean(mean(velocity_snapshot_t(:,:,:,1),2),3),'s')
                hold all
                xaxis = (0.5:gnbins(1)-0.5)*xaxis_scale;
                plot(xaxis,mean(sum(sum(mass_flux(:,:,:,1,i),2),3)/(0.1*Nvflux_ave),5))
                %plot(mean(mean(v_bins(:,:,:,1,i),2),3))
                uin = 1.8; % max(mean(mean(velocity_snapshot_t(:,:,:,1),2),3));
                uout =0.0;
                [u_analy,dudt_analy] = burger_analy(uout,uin,2.7,xaxis,(i)*delta_t*tplot*sigmatoA,-7,'2nd');
                %vu = burgers_solution ( 2.5, size(xaxis), xaxis, 10, 0:8);
                plot(xaxis,u_analy,'k')%*(max(mean(mean(velocity_snapshot_t(:,:,:,1),2),3))/uin))
                axis([0 max(xaxis) -0.1 2.1])
                legend('MD velocity','Burgers eqn')
                xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                xlabel('x'); ylabel('u_x'); 
                hold off

                set(gca,'xtick',[])
            g(5)=subplot(6,1,5);
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
                %axis([0 max(xaxis) min(mean(mean(pressure_surface_t(:,:,:,1,4),2),3))*MDPressuretoPa*Pa2kbar ...
                %                   max(mean(mean(pressure_surface_t(:,:,:,1,4),2),3))*MDPressuretoPa*Pa2kbar    ])
                axis([0 max(xaxis) 0 0.9 ])
                xaxis = (1.0:gnbins(1)+1)*xaxis_scale;
                gridxy(xaxis,'LineStyle',':')
                xlabel('x'); ylabel('\Pi'); 
                hold off
                set(gca,'xtick',[])
             g(6)=subplot(6,1,6);
            	[velocity_snapshot_t,velocity_flux,pressure_surface] = ...
                    read_vflux(i,resultfile_dir,gnbins,nd);
                %velocity_flux = velocity_flux_t;
                %pressure_surface = pressure_surface_t;

                           % %Calculate total CV flux and change in mass
                 totalflux =((velocity_flux(:,:,:,:,1)+velocity_flux(:,:,:,:,4)))/(binsize(1)) ...
                           +((velocity_flux(:,:,:,:,2)+velocity_flux(:,:,:,:,5)))/(binsize(2)) ...
                           +((velocity_flux(:,:,:,:,3)+velocity_flux(:,:,:,:,6)))/(binsize(3));
                 totalpressure =((pressure_surface(:,:,:,:,1)-pressure_surface(:,:,:,:,4)))/(binsize(1)) ...
                               +((pressure_surface(:,:,:,:,2)-pressure_surface(:,:,:,:,5)))/(binsize(2)) ...
                               +((pressure_surface(:,:,:,:,3)-pressure_surface(:,:,:,:,6)))/(binsize(3));
                [velocity_snapshot_tplus1(:,:,:,:)] =  read_vflux(i-1,resultfile_dir,gnbins,nd);
                [velocity_snapshot_tplus2(:,:,:,:)] = velocity_snapshot_t %read_vflux(i+1,resultfile_dir,gnbins,nd);
                 dvelocitydt(:,:,:,:) =  ( velocity_snapshot_tplus2(:,:,:,:) ...
                                          -velocity_snapshot_tplus1(:,:,:,:))/(delta_t*Nvflux_ave);
                 xaxis = (1:gnbins)*xaxis_scale;
                 Error(:,1)=((squeeze(totalpressure(:,:,:,1)) ...
                           -squeeze(totalflux(:,:,:,1)) )    ...
                           -squeeze(dvelocitydt(:,:,:,1)) ...
                        + ((squeeze(totalpressure(:,:,:,1)) ...
                           -squeeze(totalflux(:,:,:,1)))     ...
                           -squeeze(dvelocitydt(:,:,:,1))));

                 sum(Error(10:90))


                 plot(xaxis,squeeze(sum(sum((totalpressure(:,:,:,1)-totalflux(:,:,:,1)),2),3)),'r')
                 hold all
                 plot(xaxis,squeeze(sum(sum(   dvelocitydt(:,:,:,1),2),3)),'b--')
                 plot(xaxis,-dudt_analy*400,'k')%*(max(mean(mean(velocity_snapshot_t(:,:,:,1),2),3))/uin))

                    
                %disp(strcat(' iter= ',num2str(i), ...
                %            ' pressure= ',num2str(totalpressure(50,4,4,1)), ...
                %            ' flux= ', num2str(totalflux(50,4,4,1)), ...
                %            ' dt= ', num2str(dvelocitydt(50,4,4,1))))

                 %plot(squeeze(sum(sum((totalpressure(:,:,:,1) ...
                 %                       -  totalflux(:,:,:,1) ...
                 %                       -dvelocitydt(:,:,:,1)),2),3)),'r')
                 axis([0 max(xaxis) -0.1 20.0])
                 xaxis = (0.5:gnbins+0.5)*xaxis_scale;
                 gridxy(xaxis,'LineStyle',':')
                 xlabel('x'); ylabel('d (\rho u_x) / dt');
                 hold off


            %Calculate lagrangian averages!
            c(:,i) = delta_t*squeeze(sum(sum((totalpressure(:,:,:,1) ...
                 -  totalflux(:,:,:,1)),2),3))./sum(sum(sum(mass_flux(:,:,:,:,i),2),3),4);

            %Wave_average(:,1) = Wave_average(:,1) + mean(mean(mass_snapshot(win_srt:win_srt+win_sze,:,:,i+1),2),3);
            %Wave_average(:,2) = Wave_average(:,2) + sum(sum(mass_flux(win_srt:win_srt+win_sze,:,:,1,i),2),3) + ...
            %                                        sum(sum(mass_flux(win_srt:win_srt+win_sze,:,:,4,i),2),3);
            %Wave_average(:,3) = Wave_average(:,3) + 4*mean(mean(velocity_snapshot_t(win_srt:win_srt+win_sze,:,:,1),2),3);
            %Wave_average(:,4) = Wave_average(:,4) + mean(mean(pressure_surface_t(win_srt:win_srt+win_sze,:,:,1,4),2),3) ... 
            %                                       -mean(mean(   velocity_flux_t(win_srt:win_srt+win_sze,:,:,1,4),2),3);
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
            %Set background to white
            whitebg(h,'w')
            %Move plots closer together
            for n=1:6
                p = get(g(n), 'pos');
                p(2) = p(2) -0.0 + 0.00*(i-1);
                p(4) = p(4) + 0.03;
                set(g(n), 'pos', p);
            end


            drawnow
             

            if (mod(i,50))
                savefig(strcat('shock',num2str(i)),'png')
            end

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
