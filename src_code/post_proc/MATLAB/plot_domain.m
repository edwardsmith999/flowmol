%Plot walls, thermostatted region and sliding walls

%Setup Domain
resolution = (2*nbins);
[x,y,z] = meshgrid(-0.6*globaldomain:globaldomain/resolution:0.6*globaldomain);

scrsz = get(0,'ScreenSize');
%fig_cube = figure('Position',[2000 scrsz(4)/4 scrsz(3)/6 scrsz(4)/2]);
%Plot grid of bins size
%plot_cube(x,y,z,0.5*binsize,-0.5*binsize,'red',fig_cube)
grid on;
set(gca,'xtick',[-0.5*globaldomain(1):binsize(1):0.5*globaldomain(1)]);
set(gca,'ytick',[-0.5*globaldomain(2):binsize(2):0.5*globaldomain(2)]);
set(gca,'ztick',[-0.5*globaldomain(3):binsize(3):0.5*globaldomain(3)]);

%Plot wall regions
for ixyz = 1:3
    direction = zeros(1,3);
    direction(ixyz) = 1;
    % Bottom wall
    if (wallbot(ixyz) ~= 0)
        bot_wall_top =  not(direction) .* globaldomain *0.5 - direction .* (globaldomain* 0.5 - wallbot);
        bot_wall_bot = -not(direction) .* globaldomain *0.5 - direction .* (globaldomain ) * 0.5;
        plot_cube(x,y,z,bot_wall_top,bot_wall_bot,[1, 0, 0],fig1)
    end
    % Top wall
    if (walltop(ixyz) ~= 0)
        top_wall_top =  not(direction) .* globaldomain *0.5 + direction .* (globaldomain* 0.5 );
        top_wall_bot = -not(direction) .* globaldomain *0.5 + direction .* (globaldomain* 0.5 - walltop);
        plot_cube(x,y,z,top_wall_top,top_wall_bot,'blue',fig1)
    end
    %Bottom Thermostatted region
    if (thermstatbot(ixyz) ~= 0)
        bot_therm_top =  not(direction) .* globaldomain *0.5 - direction .* (globaldomain* 0.5 - thermstatbot);
        bot_therm_bot = -not(direction) .* globaldomain *0.5 - direction .* (globaldomain ) * 0.5;
        plot_cube(x,y,z,bot_therm_top,bot_therm_bot,'yellow',fig1) 
    end
    %Top Thermostatted region
    if (thermstattop(ixyz) ~= 0)
        top_therm_top =  not(direction) .* globaldomain *0.5 + direction .* (globaldomain* 0.5 );
        top_therm_bot = -not(direction) .* globaldomain *0.5 + direction .* (globaldomain* 0.5 - thermstattop);
        plot_cube(x,y,z,top_therm_top,top_therm_bot,'yellow',fig1) 
    end
    %Bottom Sliding regions
    if (slidedistbot(ixyz) ~= 0)
        bot_slide_top =  not(direction) .* globaldomain *0.5 - direction .* (globaldomain* 0.5 - slidedistbot);
        bot_slide_bot = -not(direction) .* globaldomain *0.5 - direction .* (globaldomain ) * 0.5;
        plot_slide(x,y,z,bot_slide_top,bot_slide_bot,wallslidev,fig1) 
    end
    %Top Sliding regions
    if (slidedisttop(ixyz) ~= 0)
        top_slide_top =  not(direction) .* globaldomain *0.5 + direction .* (globaldomain* 0.5 );
        top_slide_bot = -not(direction) .* globaldomain *0.5 + direction .* (globaldomain* 0.5 - slidedisttop);
        plot_slide(x,y,z,top_slide_top,top_slide_bot,wallslidev,fig1) 
    end
    
    %legend('Bottom wall','Top wall','Bottom Thermostatted region','Top Thermostatted region','Sliding','location', 'NorthOutside')

end