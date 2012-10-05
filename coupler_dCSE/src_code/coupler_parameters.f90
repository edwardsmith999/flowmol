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


  ! tags values
	integer, parameter :: & 
		VOID=-666,         &
		CPL=777,		  &
		CFD=888 			! tag values marking the input file for domain data

	! realms
	integer, parameter :: &
		cfd_realm = 1,    &
		md_realm  = 2
	character(len=*),parameter :: & 
		realm_name(2) = (/ "CFD", "MD " /)

	! error codes
	integer, parameter :: & 
		COUPLER_ERROR_REALM  = 1,     &  ! wrong realm value
		COUPLER_ERROR_ONE_REALM  = 2, &  ! one realm missing
		COUPLER_ERROR_INIT       = 3, &  ! initialisation error
		COUPLER_ERROR_INPUT_FILE = 4, &  ! wrong value in input file
		COUPLER_ERROR_READ_INPUT = 5, &  ! error in processing input file or data transfers
		COUPLER_ERROR_CONTINUUM_FORCE = 6, & !the region in which the continuum constrain force is apply spans over two MD domains
		COUPLER_ABORT_ON_REQUEST = 7, & ! used in request_abort 
		COUPLER_ABORT_SEND_CFD   = 8, &	! error in coupler_cfd_send
		COUPLER_ERROR_CART_COMM   = 9 	! Wrong comm value in CPL_Cart_coords

    ! CFD code ids useful for MD
    integer, parameter ::   & 
		couette_serial=101, & 
		couette_parallel=102

    !derived type that describes the local CFD grid for MD subroutine 
    type cfd_grid_info
        !indices for the cfd grid box covering the local domain, outside and inside
        integer icmino,icmaxo,jcmino,jcmaxo,kcmino,kcmaxo
        integer icmin, icmax,  jcmin, jcmax, kcmin, kcmax
        !global limits of CFD grid
        integer gicmin, gicmax, gjcmin, gjcmax, gkcmin,gkcmax 
        
        real(kind(0.d0)) xmin,xmax,dx,zmin,zmax,dz,ymin,ymax
        real(kind(0.d0)), allocatable :: y(:) 
        
    end type cfd_grid_info


end module coupler_parameters
