!-----------------------------------------------------------------------------
!
!                                Move Particles
! Move particles as a result of forces using the verlet time integration
!
!-----------------------------------------------------------------------------

module module_move_particles_lfv

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD

end module module_move_particles_lfv

!=========================================================================================!
! simulation_move_particles_lfv
! authors: Ed Smith & David Trevelyan
!
! description:
!     simulation_move_particles_lfv performs time integration on the equations of motion
!     for each particle. For a full description of the algorithm, please see "Computer
!     simulation of liquids" by Allen & Tildesley. 
!
!     basic algorithm:
!        - v(t+dt/2) = v(t-dt/2) + a(t)dt
!        - r(t+dt)   = r(t)      + v(t+dt/2)dt
!
!     computations that are only required by extended system ensembles (for example,
!     evaluating the time derivative of the damping parameter "zeta" in the Nosé-Hoover
!     algorithm) are, in the interests of clarity, computed in subroutines contained
!     within simulation_move_particles_lfv. 
!
!-----------------------------------------------------------------------------------------!
subroutine simulation_move_particles_lfv
use module_move_particles_lfv
implicit none
	
	integer :: n
	double precision :: ascale, bscale
	double precision, save :: zeta=0.d0
	double precision, dimension(np,nd) :: U

	select case(ensemble)

	case(nve)
		do n=1,np
			v(n,:) = v(n,:) + delta_t*a(n,:)
			r(n,:) = r(n,:) + delta_t*v(n,:)
		end do

	case(nvt_NH)
		call evaluate_NH_params
		do n=1,np
	        v(n,:) = v(n,:)*ascale + a(n,:)*delta_t*bscale
			r(n,:) = r(n,:)        + v(n,:)*delta_t			
		end do

	case(nvt_GIK)
		stop 'GIK thermostat not available with leap-frog Verlet'

	case(nvt_PUT_NH)
		call evaluate_U_PUT
		call evaluate_NH_params_PUT
		do n=1,np
	        v(n,:) = v(n,:)*ascale + a(n,:)*delta_t*bscale + zeta*U(n,:)*delta_t*bscale
			r(n,:) = r(n,:)        + v(n,:)*delta_t			
		end do
	
	case(nvt_pwa_NH)
		stop 'pwa_NH not available for leap-frog Verlet.'
				
	case default
		stop 'Unrecognised move option in simulation_move_particles_lfv.f90'

	end select

contains
	
	!-----------------------------------------------------------------------
	!Evaluate Nosé-Hoover parameters
	subroutine evaluate_NH_params
		implicit none
	
		double precision, dimension (nd) :: vel
		double precision :: v2sum
		double precision :: Q
		double precision :: dzeta_dt
		
		v2sum = 0.d0
		do n=1,np
			vel(:) = v(n,:) - 0.5d0*a(n,:)*delta_t	
			v2sum = v2sum + dot_product(vel,vel)
		end do
		
		Q        = np*delta_t
		dzeta_dt = (v2sum - (nd*np + 1)*inputtemperature)/Q
		zeta     = zeta + delta_t*dzeta_dt
		bscale   = 1.0/(1.0+0.5*delta_t*zeta)
		ascale   = (1-0.5*delta_t*zeta)*bscale

	end subroutine evaluate_NH_params

	!-----------------------------------------------------------------------
	!Evaluate N-H parameters for profile unbiased thermostat
	subroutine evaluate_NH_params_PUT
		implicit none
		
		double precision :: pec_v2sum
		double precision :: Q
		double precision :: dzeta_dt
		double precision, dimension(nd) :: pec_v

		pec_v2sum = 0.d0
		do n=1,np
			pec_v(:)  = v(n,:) - U(n,:) - 0.5d0*a(n,:)*delta_t      ! PUT: Find peculiar velocity
			pec_v2sum = pec_v2sum + dot_product(pec_v,pec_v)        ! PUT: Sum peculiar velocities squared
		end do
		Q        = np*delta_t                                       ! PUT: Thermal inertia
		dzeta_dt = (pec_v2sum - (np*nd+1)*inputtemperature)/Q       ! PUT: dzeta_dt(t-dt)
		zeta     = zeta + delta_t*dzeta_dt                          ! PUT: zeta(t)
		bscale   = 1.0/(1.0+0.5*zeta*delta_t)                       
		ascale   = (1.0-0.5*zeta*delta_t)*bscale

	end subroutine evaluate_NH_params_PUT

	!-----------------------------------------------------------------------
	!Evaluate streaming velocity for each particle
	subroutine evaluate_U_PUT
		use calculated_properties_MD,   only: nbins, get_mass_slices, get_velo_slices
		use shear_info_MD,              only: shear_plane
		implicit none
		
		integer :: slicebin
		integer, dimension(:), allocatable :: m_slice
		double precision, dimension(nd) :: slicebinsize
		double precision, dimension(:,:), allocatable :: v_slice
		double precision, dimension(:,:), allocatable :: v_avg
	
		allocate(m_slice(nbins(shear_plane)))                                   ! PUT: Allocate instantaneous mass slices
		allocate(v_slice(nbins(shear_plane),nd))                                ! PUT: Allocate instantaneous velocity slices
		allocate(v_avg(nbins(shear_plane),nd))                                  ! PUT: Allocate instantaneous velocity averages
		slicebinsize(:) = domain(:)/nbins(:)                                    ! PUT: Get bin size for PUT
		m_slice = get_mass_slices(shear_plane)                                  ! PUT: Get total mass in all slices
		v_slice = get_velo_slices(shear_plane)                                  ! PUT: Get total velocity in all slices
		do slicebin=1,nbins(shear_plane)                                        ! PUT: Loop through all slices
			v_avg(slicebin,:) = v_slice(slicebin,:)/m_slice(slicebin)           ! PUT: average velocity
		end do
		
		do n=1,np
			slicebin = ceiling((r(n,shear_plane)+halfdomain(shear_plane))/&
			                    slicebinsize(shear_plane))
			if (slicebin > nbins(shear_plane)) slicebin = nbins(shear_plane)    ! PUT: Prevent out-of-range values
			if (slicebin < 1) slicebin = 1                                      ! PUT: Prevent out-of-range values
			U(n,:) = v_avg(slicebin,:)
		end do

		deallocate(m_slice)
		deallocate(v_slice)
		deallocate(v_avg)

	end subroutine evaluate_U_PUT

end subroutine simulation_move_particles_lfv

!======================================================================================
!======================================================================================
!--------------------------------------------------------------------------------------

