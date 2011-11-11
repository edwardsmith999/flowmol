! Dummy module needed for stand alone runs
! of MD code
! It sets several parameter to nedeed for MD run
module coupler_md_global_data
        implicit none
        save

	! basic coupling data
        logical, parameter :: use_coupling = .false.

        integer MD_COMM                  ! MD communicator
        character(len=*), parameter :: code_name = "MD"
        integer, parameter :: CFD = 1, MD = 2

        integer :: nsteps = 1, average_period = 1

        real xL_md, yL_md, zL_md ! macroscopic sizes of MD domain. 

        integer ibmin_md(1), ibmax_md(1),jbmin_md(1), jbmax_md(1), kbmin_md(1), kbmax_md(1)

contains 

        subroutine create_communicators
                implicit none

        end subroutine create_communicators

end module coupler_md_global_data
