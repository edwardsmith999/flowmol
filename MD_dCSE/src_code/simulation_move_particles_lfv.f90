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
	double precision, dimension(:,:),allocatable :: U

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
		allocate(U(nd,np))
		call evaluate_U_PUT
		call evaluate_NH_params_PUT
		do n=1,np
	        v(:,n)     = v(:,n)*ascale + (a(:,n)+zeta*U(:,n))*delta_t*bscale
			r(:,n)     = r(:,n)        + v(:,n)*delta_t			
		end do
		deallocate(U)
	
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

	simtime = simtime + delta_t

contains

	!-----------------------------------------------------------------------
	!Evaluate "true" positions and velocities without periodic wrapping
	subroutine simulation_move_particles_true_lfv
		use shear_info_MD, only: le_sp, le_sv, le_sd
		implicit none

		vtrue = v

		do n=1,np
			vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
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

		print*, 'Warning: PUT evaluation only applicable to Lees-Edwards systems'

		pec_v2sum = 0.d0
		do n=1,np
			pec_v(:)  = v(:,n) - U(:,n) - 0.5d0*a(:,n)*delta_t      ! PUT: Find peculiar velocity
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

		print*, 'Warning: PUT evaluation only applicable to Lees-Edwards systems'
	
		allocate(m_slice(nbins(le_sp)))                                   ! PUT: Allocate instantaneous mass slices
		allocate(v_slice(nbins(le_sp),nd))                                ! PUT: Allocate instantaneous velocity slices
		allocate(v_avg(nbins(le_sp),nd))                                  ! PUT: Allocate instantaneous velocity averages
		slicebinsize(:) = domain(:)/nbins(:)                              ! PUT: Get bin size for PUT
		m_slice = get_mass_slices(le_sp)                                  ! PUT: Get total mass in all slices
		v_slice = get_velo_slices(le_sp)                                  ! PUT: Get total velocity in all slices
		do slicebin=1,nbins(le_sp)                                        ! PUT: Loop through all slices
			v_avg(slicebin,:) = v_slice(slicebin,:)/m_slice(slicebin)     ! PUT: average velocity
		end do
		
		do n=1,np
			slicebin = ceiling((r(le_sp,n)+halfdomain(le_sp))/&
			                    slicebinsize(le_sp))
			if (slicebin > nbins(le_sp)) slicebin = nbins(le_sp)          ! PUT: Prevent out-of-range values
			if (slicebin < 1) slicebin = 1                                ! PUT: Prevent out-of-range values
			U(:,n) = v_avg(slicebin,:)
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
				
		integer	:: n, thermostatnp
		double precision :: freq, dzeta_dt, v2sum, Q
		double precision :: ascale, bscale
		double precision, dimension(nd)	:: vel

		if (tag_thermostat_active) then
		
			! -------------   W  A  R  N  I  N  G  -------------------------	
			! --------------------------------------------------------------
			! 
			!  This is only applicable for systems under a homogeneous
			!  temperature field. It is NOT intended to be used as a tool
			!  for applying temperature gradients or "regional" thermostats. 
			!  If this functionality is desired, the routine will need to 
			!  be redeveloped with multiple "v2sum"s for each region that 
			!  is thermostatted separately.
			! 
			! --------------------------------------------------------------

			! PUT only works for serial code.
			if ( any(tag(:).eq.PUT_thermo) ) call evaluate_U_PUT  

			! Warn user of lack of development
			if ( any(tag(:).eq.z_thermo) ) then
				if (mod(iter,tplot).eq.0 .and. irank .eq. iroot) then
				print*, 'Warning: thermostatting only in z direction. This &
				       & has not been fully developed and requires checking.'
				end if
			end if

			v2sum = 0.d0
			thermostatnp = 0
			do n = 1, np

				select case ( tag(n) )
				case ( thermo, teth_thermo, teth_thermo_slide )
					vel(:) = v(:,n) - 0.5d0*a(:,n)*delta_t

				case ( PUT_thermo )
					vel(:) = v(:,n) - U(:,n) - 0.5d0*a(:,n)*delta_t

				case ( z_thermo )
					vel(:) = v(:,n) - 0.5d0*a(:,n)*delta_t

				case default
					cycle ! Don't include non-thermostatted molecules

				end select

				v2sum = v2sum + dot_product(vel,vel)
				thermostatnp = thermostatnp + 1

			enddo


			!Obtain global sums for all parameters
			call globalSumInt(thermostatnp)
			call globalSum(v2sum)	

			Q        = thermostatnp * delta_t
			dzeta_dt = (v2sum - (nd*thermostatnp + 1)*inputtemperature) / Q
			zeta 	 = zeta + delta_t*dzeta_dt
			bscale	 = 1.0/(1.0+0.5*delta_t*zeta)
			ascale	 = (1-0.5*delta_t*zeta)*bscale

		else
			!Reduces to the un-thermostatted equations
			ascale = 1
			bscale = 1
		endif

		!call pointsphere((/ 0.0, 0.0, -6.34 /),2.d0)

		!Step through each particle n
		do n = 1,np        
			select case (tag(n))
			case (free)
				!Leapfrog mean velocity calculated here at v(t+0.5delta_t) = v(t-0.5delta_t) + a*delta_t 
				!Leapfrog mean position calculated here at r(t+delta_t) = r(t) + v(t+0.5delta_t)*delta_t
				v(:,n) = v(:,n) + delta_t * a(:,n) 	!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n)	!Position calculated from velocity
			case (fixed)
				!Fixed Molecules - no movement r(n+1) = r(n)
			case (fixed_slide)
				!Fixed with constant sliding speed
				r(:,n) = r(:,n) + delta_t*slidev(:,n)	!Position calculated from velocity

				!Moving piston for shock wave with 1000 equilibrate, 100 piston moving and 
				!no moving wall for next 2000*0.005 time units taken for wave to cover whole domain 
				!at which point the simulation blows up!
				!if (iter .lt. 6000) then !Initialisation
				!	!Fixed Molecules - no movement r(n+1) = r(n)
				!elseif (iter .ge. 6000 .and. iter .lt. 8000) then
				!	if (irank .eq. iroot) then
				!		delta_t = 0.0005 !Reduce timestep
				!		call globalbroadcast(delta_t,1,irank)
				!	endif
				!	freq = 10
					!t = freq*iter*delta_t
				!	r(:,n) = r(:,n) + delta_t*slidev(:,n)!*sin(t)	!Position calculated from velocity
				!elseif (iter .ge. 8000) then
					!Fixed Molecules - no movement r(n+1) = r(n)
				!else
					!Fixed Molecules - no movement r(n+1) = r(n)
				!endif

			case (teth)
				!Tethered molecules
				call tether_force(n)
				v(:,n) = v(:,n) + delta_t * a(:,n) 	!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n)	!Position calculated from velocity
			case (thermo)
				!Nose Hoover Thermostatted Molecule
				v(1,n) = v(1,n)*ascale + a(1,n)*delta_t*bscale
				r(1,n) = r(1,n)    +     v(1,n)*delta_t			
				v(2,n) = v(2,n)*ascale + a(2,n)*delta_t*bscale
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t				
				v(3,n) = v(3,n)*ascale + a(3,n)*delta_t*bscale
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	
			case (teth_thermo)
				!Thermostatted Tethered molecules unfixed with no sliding velocity
				call tether_force(n)
				v(1,n) = v(1,n)*ascale + a(1,n)*delta_t*bscale
				r(1,n) = r(1,n)    +     v(1,n)*delta_t			
				v(2,n) = v(2,n)*ascale + a(2,n)*delta_t*bscale
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t				
				v(3,n) = v(3,n)*ascale + a(3,n)*delta_t*bscale
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	
			case (teth_slide)
				!Tethered molecules with sliding velocity
				call tether_force(n)
				v(:,n) = v(:,n) + delta_t * a(:,n) 							!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n) + delta_t*slidev(:,n)	!Position calculated from velocity+slidevel
			case (teth_thermo_slide)
				!Thermostatted Tethered molecules unfixed with sliding velocity
				call tether_force(n)
				v(1,n) = v(1,n)*ascale + a(1,n)*delta_t*bscale
				r(1,n) = r(1,n)    +     v(1,n)*delta_t	+ slidev(1,n)*delta_t		
				v(2,n) = v(2,n)*ascale + a(2,n)*delta_t*bscale
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t	+ slidev(2,n)*delta_t			
				v(3,n) = v(3,n)*ascale + a(3,n)*delta_t*bscale
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	+ slidev(3,n)*delta_t
			case (PUT_thermo)
				!Profile unbiased thermostat (Nose-Hoover)
	        	v(:,n) = v(:,n)*ascale + (a(:,n)+zeta*U(:,n))*delta_t*bscale
				r(:,n) = r(:,n)        + v(:,n)*delta_t			
			case (z_thermo)
				!Thermostat in the z direction only (Nose-Hoover)
				v(1,n) = v(1,n) + delta_t * a(1,n) 	
				r(1,n) = r(1,n) + delta_t * v(1,n)	
				v(2,n) = v(2,n) + delta_t * a(2,n) 	
				r(2,n) = r(2,n) + delta_t * v(2,n)	
				v(3,n) = v(3,n)*ascale + a(3,n)*delta_t*bscale
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	
			case (10)

			case default
				call error_abort("Invalid molecular Tag")
			end select

			!Specular walls
			if (specular_wall(1) .ne. 0.0) call specular_flat_wall(n, ascale, bscale, 1, globaldomain(1)/2.d0-specular_wall(1))
			if (specular_wall(2) .ne. 0.0) call specular_flat_wall(n, ascale, bscale, 2, globaldomain(2)/2.d0-specular_wall(2))
			if (specular_wall(3) .ne. 0.0) call specular_flat_wall(n, ascale, bscale, 3, globaldomain(3)/2.d0-specular_wall(3))
		enddo



	end subroutine simulation_move_particles_lfv_tag

end subroutine simulation_move_particles_lfv

!======================================================================================
! Minimal form of the move particles subroutine
!

subroutine simulation_move_particles_lfv_basic
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

end subroutine simulation_move_particles_lfv_basic
!======================================================================================



!--------------------------------------------------------------------------------------

subroutine specular_flat_wall(molno, ascale, bscale, dir, spec_pos)
	use module_molecule_properties
	use arrays_MD
	use interfaces
	implicit none

	integer,intent(in)			   :: molno, dir
	double precision,intent(in)	   :: ascale, bscale, spec_pos

	double precision, dimension(3) :: r_glob
	integer                        :: normal
	double precision               :: newxd

	!Get position in global co-ordinates
	r_glob(1) = r(1,molno) - (halfdomain(1)*(npx-1)) + domain(1)*(iblock-1)
	r_glob(2) = r(2,molno) - (halfdomain(2)*(npy-1)) + domain(2)*(jblock-1)
	r_glob(3) = r(3,molno) - (halfdomain(3)*(npz-1)) + domain(3)*(kblock-1)

	if (abs(r_glob(dir)) .gt. spec_pos) then

		!print'(a,i8,4f10.5)', 'Greater than spec_pos', molno,r_glob(dir),abs(r_glob(dir)),spec_pos,r(dir,molno)

		!Get normal direction of molecule by checking if top or bottom
		if ( r_glob(dir) .lt. 0) then
			normal = 1   !normal should point outwards
	 	else
			normal = -1  !normal should point inwards
		endif

		! Reflect normal velocity.
		v(dir,molno) = normal*abs(v(dir,molno))
		!Calculate distance over the spectral barrier
		newxd = abs(r_glob(dir)) - spec_pos
		!Move molecule same distance back into the domain on other side of spectral barrier
		r(dir, molno) = r(dir, molno) + normal*2.d0*newxd

		!write(7777,'(4i8,4f10.5)'), irank,dir,normal, molno,newxd,r_glob(dir), r(dir,molno),spec_pos

	endif


end subroutine specular_flat_wall 
