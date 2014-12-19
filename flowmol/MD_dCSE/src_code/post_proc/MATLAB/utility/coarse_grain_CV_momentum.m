% Return surface values from a large control volume

function[LV_velocity_snap, ...
         LV_velocity_flux, ...
         LV_pressure_surface, ...
         LV_facesize,         ...
         LV_Fext ] = coarse_grain_CV_momentum(LV_botx, LV_topx, ...
                                              LV_boty, LV_topy, ...
                                              LV_botz, LV_topz,m,resultfile_dir)

    %Read Header file
    read_header

    %Read CV field
    if (nargout == 5) %exist('LV_Fext','var'))
        [velocity_snapshot, ...
        velocity_flux,   ...
        pressure_surface, ...
        F_ext                  ] = read_vflux(m,resultfile_dir,gnbins,nd);
    else
        [velocity_snapshot, ...
        velocity_flux,   ...
        pressure_surface ]     = read_vflux(m,resultfile_dir,gnbins,nd);
    end

    %==============Get size of CV face==============
    LV_facesize(1) = (LV_topx-LV_botx+1)*binsize(1);
    LV_facesize(2) = (LV_topy-LV_boty+1)*binsize(2);
    LV_facesize(3) = (LV_topz-LV_botz+1)*binsize(3);

    %==============velocity snapshot=================
    LV_velocity_snap        = squeeze(sum(sum(sum( ...
                                velocity_snapshot(LV_botx:LV_topx, ...
                                                  LV_boty:LV_topy, ...
                                                  LV_botz:LV_topz,:) ...
                             ,1),2),3));

    %==============surface flux==============
    % X Surfaces
    LV_velocity_flux(:,1) = sum(sum(sum(velocity_flux(    LV_botx, ...
                                                       LV_boty:LV_topy, ...
                                                       LV_botz:LV_topz,:,1) ...
                             ,1),2),3);
    LV_velocity_flux(:,4) = sum(sum(sum(velocity_flux(    LV_topx, ...
                                                       LV_boty:LV_topy, ...
                                                       LV_botz:LV_topz,:,4) ...
                             ,1),2),3);
    % Y Surfaces
    LV_velocity_flux(:,2) = sum(sum(sum(velocity_flux(LV_botx:LV_topx, ...
                                                            LV_boty, ...
                                                       LV_botz:LV_topz,:,2) ...
                             ,1),2),3);
    LV_velocity_flux(:,5) = sum(sum(sum(velocity_flux( LV_botx:LV_topx, ...
                                                            LV_topy, ...
                                                       LV_botz:LV_topz,:,5) ...
                             ,1),2),3);
    % Z Surfaces
    LV_velocity_flux(:,3) = sum(sum(sum(velocity_flux( LV_botx:LV_topx, ...
                                                       LV_boty:LV_topy, ...
                                                            LV_botz,:,3) ...
                             ,1),2),3);
    LV_velocity_flux(:,6) = sum(sum(sum(velocity_flux( LV_botx:LV_topx, ...
                                                       LV_boty:LV_topy, ...
                                                            LV_topz,:,6) ...
                             ,1),2),3);


    % ==============surface forces==============
    % X Surfaces
    LV_pressure_surface(:,1) = sum(sum(sum(pressure_surface(    LV_botx, ...
                                                            LV_boty:LV_topy, ...
                                                            LV_botz:LV_topz,:,1) ...
                             ,1),2),3);
    LV_pressure_surface(:,4) = sum(sum(sum(pressure_surface(    LV_topx, ...
                                                            LV_boty:LV_topy, ...
                                                            LV_botz:LV_topz,:,4) ...
                             ,1),2),3);
    % Y Surfaces
    LV_pressure_surface(:,2) = sum(sum(sum(pressure_surface(LV_botx:LV_topx, ...
                                                                LV_boty, ...
                                                            LV_botz:LV_topz,:,2) ...
                            ,1),2),3);
    LV_pressure_surface(:,5) = sum(sum(sum(pressure_surface( LV_botx:LV_topx, ...
                                                                LV_topy, ...
                                                             LV_botz:LV_topz,:,5) ...
                             ,1),2),3);
    % Z Surfaces
    LV_pressure_surface(:,3) = sum(sum(sum(pressure_surface( LV_botx:LV_topx, ...
                                                             LV_boty:LV_topy, ...
                                                                LV_botz,:,3) ...
                             ,1),2),3);
    LV_pressure_surface(:,6) = sum(sum(sum(pressure_surface( LV_botx:LV_topx, ...
                                                             LV_boty:LV_topy, ...
                                                                LV_topz,:,6) ...
                             ,1),2),3);

    %==============external forces=================
    if (nargout == 5) %exist('LV_Fext','var'))
        LV_Fext = squeeze(sum(sum(sum( ...
                                    F_ext(LV_botx:LV_topx, ...
                                          LV_boty:LV_topy, ...
                                          LV_botz:LV_topz,:) ...
                       ,1),2),3));
    end




end