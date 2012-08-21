!-----------------------------------------------------------------------------
!
!                                Move Particles
! Move particles as a result of forces using the verlet time integration
!
!-----------------------------------------------------------------------------

module module_move_particles_lfv

	use interfaces
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

			v(1,n) = v(1,n) + delta_t*a(1,n)
			v(2,n) = v(2,n) + delta_t*a(2,n)
			v(3,n) = v(3,n) + delta_t*a(3,n)

			r(1,n) = r(1,n) + delta_t*v(1,n)
			r(2,n) = r(2,n) + delta_t*v(2,n)
			r(3,n) = r(3,n) + delta_t*v(3,n)

		end do

	case(nvt_NH)
		call evaluate_NH_params
		do n=1,np
	        v(:,n)     = v(:,n)*ascale + a(:,n)*delta_t*bscale
			r(:,n)     = r(:,n)        + v(:,n)*delta_t			
		end do

	case(nvt_GIK)
		call error_abort('GIK thermostat not available with leap-frog Verlet')

	case(nvt_PUT_NH)
		call evaluate_U_PUT
		call evaluate_NH_params_PUT
		do n=1,np
	        v(:,n)     = v(:,n)*ascale + a(:,n)*delta_t*bscale + zeta*U(n,:)*delta_t*bscale
			r(:,n)     = r(:,n)        + v(:,n)*delta_t			
		end do
	
	case(nvt_pwa_NH)
		call error_abort('pwa_NH not available for leap-frog Verlet.')

	case(tag_move)
		call simulation_move_particles_lfv_tag
						
	case default
		call error_abort('Unrecognised move option in simulation_move_particles_lfv.f90')

	end select
	
	if (rtrue_flag .eq. 1) then
		call simulation_move_particles_true_lfv
	endif

contains

	!-----------------------------------------------------------------------
	!Evaluate "true" positions and velocities without periodic wrapping
	subroutine simulation_move_particles_true_lfv
		use shear_info_MD, only: le_sp, le_sv, le_sd
		implicit none

		vtrue = v

		do n=1,np
			vtrue(n,le_sd) = v(n,le_sd) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
			rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n)
		end do

	end subroutine simulation_move_particles_true_lfv
	
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
			vel(:) = v(:,n) - 0.5d0*a(:,n)*delta_t	
			v2sum = v2sum + dot_product(vel,vel)
		end do
		call globalSum(v2sum)		
		Q        = globalnp*delta_t
		dzeta_dt = (v2sum - (nd*globalnp + 1)*inputtemperature)/Q
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
			pec_v(:)  = v(:,n) - U(n,:) - 0.5d0*a(:,n)*delta_t      ! PUT: Find peculiar velocity
			pec_v2sum = pec_v2sum + dot_product(pec_v,pec_v)        ! PUT: Sum peculiar velocities squared
		end do
		call globalSum(pec_v2sum)
		Q        = globalnp*delta_t                                 ! PUT: Thermal inertia
		dzeta_dt = (pec_v2sum - (globalnp*nd+1)*inputtemperature)/Q ! PUT: dzeta_dt(t-dt)
		zeta     = zeta + delta_t*dzeta_dt                          ! PUT: zeta(t)
		bscale   = 1.0/(1.0+0.5*zeta*delta_t)                       
		ascale   = (1.0-0.5*zeta*delta_t)*bscale

	end subroutine evaluate_NH_params_PUT

	!-----------------------------------------------------------------------
	!Evaluate streaming velocity for each particle
	subroutine evaluate_U_PUT
		use calculated_properties_MD,   only: nbins, get_mass_slices, get_velo_slices
		use shear_info_MD,              only: le_sp
		implicit none
		
		integer :: slicebin
		integer, dimension(:), allocatable :: m_slice
		double precision, dimension(nd) :: slicebinsize
		double precision, dimension(:,:), allocatable :: v_slice
		double precision, dimension(:,:), allocatable :: v_avg
	
		allocate(m_slice(nbins(le_sp)))                                   ! PUT: Allocate instantaneous mass slices
		allocate(v_slice(nbins(le_sp),nd))                                ! PUT: Allocate instantaneous velocity slices
		allocate(v_avg(nbins(le_sp),nd))                                  ! PUT: Allocate instantaneous velocity averages
		slicebinsize(:) = domain(:)/nbins(:)                                    ! PUT: Get bin size for PUT
		m_slice = get_mass_slices(le_sp)                                  ! PUT: Get total mass in all slices
		v_slice = get_velo_slices(le_sp)                                  ! PUT: Get total velocity in all slices
		do slicebin=1,nbins(le_sp)                                        ! PUT: Loop through all slices
			v_avg(slicebin,:) = v_slice(slicebin,:)/m_slice(slicebin)           ! PUT: average velocity
		end do
		
		do n=1,np
			slicebin = ceiling((r(le_sp,n)+halfdomain(le_sp))/&
			                    slicebinsize(le_sp))
			if (slicebin > nbins(le_sp)) slicebin = nbins(le_sp)    ! PUT: Prevent out-of-range values
			if (slicebin < 1) slicebin = 1                                      ! PUT: Prevent out-of-range values
			U(n,:) = v_avg(slicebin,:)
		end do

		deallocate(m_slice)
		deallocate(v_slice)
		deallocate(v_avg)

	end subroutine evaluate_U_PUT

	!---------------------------------------------------------------------
	!Tag move routine
	subroutine simulation_move_particles_lfv_tag
		use interfaces
		implicit none
				
		integer :: maxtag, n, thermostatnp 
		double precision :: dzeta_dt, pec_v2sum, v2sum,Q
		double precision :: ascale, bscale
		double precision, dimension(nd)	:: vel
		logical :: PUT

		if (all(tag.eq.8)) then		
			PUT = .true.
			call evaluate_U_PUT										! Only works for serial code.
			pec_v2sum = 0.d0
		else 
			PUT = .false.
		end if

		!Check if any molecules are thermostatted and calculate appropriate coefficients
		maxtag = maxval(tag)
		call globalMaxInt(maxtag)
		if (maxtag .ge. 4) then

			v2sum = 0.d0    									  	! Reset total v2sum
			thermostatnp = 0										! Reset count of thermostatted mols
			do n = 1, np   											! Loop all molecules
				if (tag(n) .lt. 4) cycle							! Only include thermostatted molecules - DO YOU WANT THIS LINE UNCOMMENTED?
				if (PUT) then										! PUT: If using PUT find peculiar v2sum
					vel(:) = v(:,n) - U(n,:) - 0.5d0*a(:,n)*delta_t			! PUT: Find peculiar velocity
					pec_v2sum = pec_v2sum + dot_product(vel,vel)						! PUT: Sum peculiar velocities squared
				else
					vel(:) = v(:,n) - 0.5d0*a(:,n)*delta_t	
					v2sum = v2sum + dot_product(vel,vel)
				end if
				thermostatnp = thermostatnp + 1
			enddo

			!Obtain global sums for all parameters
			call globalSumInt(thermostatnp)
			call globalSum(v2sum)	

			Q        = thermostatnp * delta_t
			dzeta_dt = (v2sum - (nd*thermostatnp + 1)*inputtemperature) / Q
			if (PUT) dzeta_dt = (pec_v2sum - (nd*thermostatnp + 1)*inputtemperature) / Q
			zeta 	 = zeta + delta_t*dzeta_dt
			bscale	 = 1.0/(1.0+0.5*delta_t*zeta)
			ascale	 = (1-0.5*delta_t*zeta)*bscale

		endif

		!Step through each particle n
		do n = 1,np        
			select case (tag(n))
			case (0)
				!Leapfrog mean velocity calculated here at v(t+0.5delta_t) = v(t-0.5delta_t) + a*delta_t 
				!Leapfrog mean position calculated here at r(t+delta_t) = r(t) + v(t+0.5delta_t)*delta_t
				v(:,n) = v(:,n) + delta_t * a(:,n) 	!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n)	!Position calculated from velocity
			case (1)
				!Fixed Molecules - no movement r(n+1) = r(n)
			case (2)
				!Fixed with constant sliding speed
				r(:,n) = r(:,n) + delta_t*slidev(:,n)	!Position calculated from velocity

				!Moving piston for shock wave with 1000 equilibrate, 100 piston moving and 
				!no moving wall for next 2000*0.005 time units taken for wave to cover whole domain 
				!at which point the simulation blows up!
				!if (iter .lt. 200) then !Initialisation
				!	!Fixed Molecules - no movement r(n+1) = r(n)
				!elseif (iter .ge. 200 .and. iter .lt. 300) then
				!	r(:,n) = r(:,n) + delta_t*slidev(:,n)	!Position calculated from velocity
				!else
				!	!Fixed Molecules - no movement r(n+1) = r(n)
				!endif
			case (3)
				!Tethered molecules
				call tether_force(n)
				v(:,n) = v(:,n) + delta_t * a(:,n) 	!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n)	!Position calculated from velocity
			case (4)
				!Nose Hoover Thermostatted Molecule
				v(1,n) = v(1,n)*ascale + a(1,n)*delta_t*bscale
				r(1,n) = r(1,n)    +     v(1,n)*delta_t			
				v(2,n) = v(2,n)*ascale + a(2,n)*delta_t*bscale
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t				
				v(3,n) = v(3,n)*ascale + a(3,n)*delta_t*bscale
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	
			case (5)
				!Thermostatted Tethered molecules unfixed with no sliding velocity
				call tether_force(n)
				v(1,n) = v(1,n)*ascale + a(1,n)*delta_t*bscale
				r(1,n) = r(1,n)    +     v(1,n)*delta_t			
				v(2,n) = v(2,n)*ascale + a(2,n)*delta_t*bscale
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t				
				v(3,n) = v(3,n)*ascale + a(3,n)*delta_t*bscale
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	
			case (6)
				!Tethered molecules with sliding velocity
				call tether_force(n)
				v(:,n) = v(:,n) + delta_t * a(:,n) 							!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n) + delta_t*slidev(:,n)	!Position calculated from velocity+slidevel
			case (7)
				!Thermostatted Tethered molecules unfixed with sliding velocity
				call tether_force(n)
				v(1,n) = v(1,n)*ascale + a(1,n)*delta_t*bscale
				r(1,n) = r(1,n)    +     v(1,n)*delta_t	+ slidev(1,n)*delta_t		
				v(2,n) = v(2,n)*ascale + a(2,n)*delta_t*bscale
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t	+ slidev(2,n)*delta_t			
				v(3,n) = v(3,n)*ascale + a(3,n)*delta_t*bscale
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	+ slidev(3,n)*delta_t
			case (8)
				!Profile unbiased thermostat (Nose-Hoover)
	        	v(:,n) = v(:,n)*ascale + a(:,n)*delta_t*bscale + zeta*U(n,:)*delta_t*bscale
				r(:,n) = r(:,n)        + v(:,n)*delta_t			
			case default
				call error_abort("Invalid molecular Tag")
			end select
		enddo
		
	end subroutine simulation_move_particles_lfv_tag

end subroutine simulation_move_particles_lfv

!======================================================================================

subroutine simulation_move_particles_lfv_opt
	use arrays_MD, only	: r,v,a
	use physical_constants_MD, only : np
	use computational_constants_MD, only : delta_t
	implicit none

	integer		:: n

	do n=1,np

		v(1,n) = v(1,n) + delta_t*a(1,n)
		v(2,n) = v(2,n) + delta_t*a(2,n)
		v(3,n) = v(3,n) + delta_t*a(3,n)

		r(1,n) = r(1,n) + delta_t*v(1,n)
		r(2,n) = r(2,n) + delta_t*v(2,n)
		r(3,n) = r(3,n) + delta_t*v(3,n)

	enddo

end subroutine simulation_move_particles_lfv_opt
!======================================================================================
!--------------------------------------------------------------------------------------

