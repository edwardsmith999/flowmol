!=============================================================================
!
!                                  Coupler parameters
!
! Data accessible from application ( molecular or continuum ).
! No need to have a separate use stamentment for this module as it is used in coupler module
!
! Lucian Anton, November 2011
!
!=============================================================================

module coupler_parameters
	implicit none
	save
     
	! realms
	integer, parameter :: COUPLER_CFD = 1, COUPLER_MD = 2
	character(len=*),parameter :: code_name(2) = (/ "CFD", "MD " /)
 	integer COUPLER_REALM 

	! error codes
	integer, parameter :: COUPLER_ERROR_REALM  = 1, &  ! wrong realm value
                          COUPLER_ERROR_ONE_REALM  = 2, &  ! one realm missing
                          COUPLER_ERROR_INIT       = 3, &  ! initialisation error
                          COUPLER_ERROR_INPUT_FILE = 4, &  ! wrong value in input file
                          COUPLER_ERROR_READ_INPUT = 5     ! error in processing input file or data transfers

end module coupler_parameters
