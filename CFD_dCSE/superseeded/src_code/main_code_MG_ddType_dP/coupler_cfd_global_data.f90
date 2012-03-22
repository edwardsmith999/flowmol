module coupler_cfd_global_data
        implicit none
        save

! basic coupling data
        logical, parameter :: use_coupling = .false.

        integer CFD_COMM                  ! CF communicator
        character(len=*), parameter :: code_name = "CFD"
        integer, parameter :: CFD = 1, MD = 2

contains 

        subroutine create_communicators
                implicit none

        end subroutine create_communicators

end module coupler_cfd_global_data
