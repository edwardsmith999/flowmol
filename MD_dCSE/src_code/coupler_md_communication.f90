! Dummy module needed for stand alone runs
! of MD code
module coupler_md_communication

contains

        subroutine boundary_box_average(send_data)
               implicit none
               logical,intent(in) :: send_data      

        end subroutine boundary_box_average

        subroutine simulation_apply_continuum_forces(iter)
               implicit none
               integer, intent(in) :: iter 
        end subroutine simulation_apply_continuum_forces

        subroutine coupler_constrain_forces(np,pressure,r,a)
                implicit none
                integer, intent(in)                  :: np
                real(kind=kind(0.d0)), intent(inout) :: a(:,:)
                real(kind=kind(0.d0)), intent(in)    :: pressure, r(:,:)
        end subroutine coupler_constrain_forces

end module coupler_md_communication
