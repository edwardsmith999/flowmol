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
    use messenger, only: globalise
    use boundary_MD, only: specular_flag, specular_flat, specular_wall, &
                           specular_radial, specular_wall_flag
	implicit none
	
	integer :: n
	real(kind(0.d0)) :: ascale, bscale, ry_old, vy_old, shear
	real(kind(0.d0)), dimension(:), allocatable, save  :: zeta    	!Parameter used in Nose Hoover thermostat

	!Apply external force field to regions of spaces
	select case(external_force_flag)
	case(0)
		!Do nothing - no force applied
	case(1)
		call simulation_apply_global_force(F_ext_ixyz,F_ext)
	case(2)
		call simulation_apply_local_force(F_ext_ixyz,F_ext, & 
										  F_ext_limits(1),F_ext_limits(2), & 
										  F_ext_limits(3),F_ext_limits(4), & 
										  F_ext_limits(5),F_ext_limits(6))
	case default
		call error_abort("Error - incorrectly specified external_force_flag")
	end select

	!Select case and evolve system in time
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


	case(SLLOD)

		!Designed to apply a linear x velocity profile in the y direction
		shear  = 2.d0/globaldomain(2)

		do n=1,np

			vy_old = v(2,n); ry_old = r(2,n)
			v(2,n) = v(2,n) + delta_t*a(2,n)
			r(2,n) = r(2,n) + delta_t*v(2,n)

			v(3,n) = v(3,n) + delta_t*a(3,n)
			r(3,n) = r(3,n) + delta_t*v(3,n)

			!Apply to bottom half of domain only
			if (globalise(r(2,n),2) .lt. 0.d0) then
				v(1,n) = v(1,n) + delta_t*a(1,n) - 0.5d0*shear*delta_t*(v(2,n)+vy_old)
				r(1,n) = r(1,n) + delta_t*v(1,n) + 0.5d0*shear*delta_t*globalise(r(2,n)+ry_old,ixyz=2)
			else
				v(1,n) = v(1,n) + delta_t*a(1,n)
				r(1,n) = r(1,n) + delta_t*v(1,n)
			endif
		end do

!		shear = (/ 2.d0/globaldomain(2), 0.d0, 0.d0 /)
!		ascale3(:) = 1.d0+0.5*delta_t*shear(:)
!		bscale3(:) = 1.d0-0.5*delta_t*shear(:)
!		do n=1,np
!		    v(:,n) = (1.d0/ascale3(:))*(a(:,n)*delta_t + v(:,n)*bscale3(:))
!			r(:,n) = (1.d0/bscale3(:))*(v(:,n)*delta_t + r(:,n)*ascale3(:))		
!		enddo

	case default
		call error_abort('Unrecognised move option in simulation_move_particles_lfv.f90')

	end select

	if (specular_flag .eq. specular_flat) then
		
		!if (specular_wall(1) .ne. 0.0) call specular_flat_wall(1, globaldomain(1)/2.d0-specular_wall(1),specular_wall_flag)
		!if (specular_wall(2) .ne. 0.0) call specular_flat_wall(2, globaldomain(2)/2.d0-specular_wall(2),specular_wall_flag)
		!if (specular_wall(3) .ne. 0.0) call specular_flat_wall(3, globaldomain(3)/2.d0-specular_wall(3),specular_wall_flag)
		if (specular_wall(1) .ne. 0.0) call specular_flat_walls(1)
		if (specular_wall(2) .ne. 0.0) call specular_flat_walls(2)
		if (specular_wall(3) .ne. 0.0) call specular_flat_walls(3)

	else if (specular_flag .eq. specular_radial) then
	
		call specular_walls_cylinders

	end if
		
	if (rtrue_flag .eq. 1) then
		call simulation_move_particles_true_lfv()
	endif

	simtime = simtime + delta_t

contains

	!-----------------------------------------------------------------------
	!Evaluate Nosé-Hoover parameters
	subroutine evaluate_NH_params
    	use messenger_data_exchange, only : globalSum
        use module_set_parameters, only : mass
		implicit none
	
		real(kind(0.d0)), dimension (nd) :: vel
		real(kind(0.d0)) :: mv2sum, thermostatnp
		real(kind(0.d0)) :: Q
		real(kind(0.d0)) :: dzeta_dt
		
		mv2sum = 0.d0
        thermostatnp = 0.d0
		do n=1,np
			vel(:) = v(:,n) + 0.5d0*a(:,n)*delta_t	
			mv2sum = mv2sum + mass(n) * dot_product(vel,vel)
        	thermostatnp = thermostatnp + mass(n)          
		end do
		call globalSum(mv2sum)		
		call globalSum(thermostatnp)		
		Q        = thermostatnp*delta_t
		dzeta_dt = (mv2sum - (real(nd*thermostatnp + 1,kind(0.d0)))*thermostattemperature(1))/Q
		zeta     = zeta + delta_t*dzeta_dt
		bscale   = 1.0/(1.0+0.5*delta_t*zeta(1))
		ascale   = (1-0.5*delta_t*zeta(1))*bscale

	end subroutine evaluate_NH_params

	!-----------------------------------------------------------------------
	!Evaluate N-H parameters for profile unbiased thermostat
	subroutine evaluate_NH_params_PUT
	    use messenger_data_exchange, only : globalSum
        use module_set_parameters, only : mass
		implicit none
		
		real(kind(0.d0)) :: pec_mv2sum, thermostatnp
		real(kind(0.d0)) :: Q
		real(kind(0.d0)) :: dzeta_dt
		real(kind(0.d0)), dimension(nd) :: pec_v

		print*, 'Warning: PUT evaluation only applicable to Lees-Edwards systems'

        call error_abort("Error -- This routine cannot work as U is not longer defined!")

		pec_mv2sum = 0.d0
		do n=1,np
			!pec_v(:)  = v(:,n) + 0.5d0*a(:,n)*delta_t - U(:,n)      ! PUT: Find peculiar velocity
			pec_mv2sum = pec_mv2sum + mass(n) * dot_product(pec_v,pec_v)        ! PUT: Sum peculiar velocities squared
        	thermostatnp = thermostatnp + mass(n)          
		end do
		call globalSum(pec_mv2sum)
		call globalSum(thermostatnp)		
		Q        = thermostatnp*delta_t                                 ! PUT: Thermal inertia
		dzeta_dt = (pec_mv2sum - (real(thermostatnp*nd+1,kind(0.d0)))*thermostattemperature(1))/Q ! PUT: dzeta_dt(t-dt)
		zeta     = zeta + delta_t*dzeta_dt                          ! PUT: zeta(t)
		bscale   = 1.0/(1.0+0.5*zeta(1)*delta_t)                       
		ascale   = (1.0-0.5*zeta(1)*delta_t)*bscale

	end subroutine evaluate_NH_params_PUT

	!-----------------------------------------------------------------------
	!Evaluate streaming velocity for each particle
	subroutine evaluate_U_PUT
		use calculated_properties_MD,   only: nbins, get_mass_slices, get_velo_slices
		use shear_info_MD,              only: le_sp
		implicit none
		
		integer :: slicebin
		integer, dimension(:), allocatable :: m_slice
		real(kind(0.d0)), dimension(nd) :: slicebinsize
		real(kind(0.d0)), dimension(:,:), allocatable :: v_slice
		real(kind(0.d0)), dimension(:,:), allocatable :: v_avg

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
			!U(:,n) = v_avg(slicebin,:)
		end do

		deallocate(m_slice)
		deallocate(v_slice)
		deallocate(v_avg)

	end subroutine evaluate_U_PUT


    function get_therm_region(n)
		use messenger, only: globalise
        use computational_constants_MD, only : nthermo
    	implicit none

        integer, intent(in) :: n
        integer :: get_therm_region
        real(kind(0.d0)), dimension(3) :: rglob

        if (nthermo .eq. 1) then
            get_therm_region = 1 
        else
            rglob =globalise(r(:,n))
            if (rglob(2) .gt. 0.d0) then
                get_therm_region = 1
            else
                get_therm_region = 2
            endif
        endif

    end function get_therm_region

	!---------------------------------------------------------------------
	!Tag move routine
	subroutine simulation_move_particles_lfv_tag
		use interfaces
		use librarymod, only: cartesianiser, cartesianisev, &
		                      cpolariser, cpolarisev
		use messenger, only: localise, globalise
	    use messenger_data_exchange, only : globalSum
		use concentric_cylinders, only: omega
        use module_set_parameters, only : mass, nthermo
		implicit none
				
		integer	:: n, tr
		real(kind(0.d0)), dimension(:), allocatable :: thermostatnp, dzeta_dt, mv2sum, Q
		real(kind(0.d0)), dimension(:), allocatable :: ascale, bscale
		real(kind(0.d0)), dimension(nd)	:: vel
		real(kind(0.d0)), dimension(nd)	:: rpol, rglob

		!Dynamically reassign tags based on spatial location
		if (dynamically_update_tags) then
			call reset_location_tags
		endif
	
        allocate(ascale(nthermo), bscale(nthermo))
		if (tag_thermostat_active) then
		
			! --------------------------------------------------------------!
			! -------------   W  A  R  N  I  N  G  -------------------------!	
			! --------------------------------------------------------------!
			! 
			!  This is only applicable for systems under a homogeneous
			!  temperature field. It is NOT intended to be used as a tool
			!  for applying temperature gradients or "regional" thermostats. 
			!  If this functionality is desired, the routine will need to 
			!  be redeveloped with multiple "mv2sum"s for each region that 
			!  is thermostatted separately.
			! 
			! --------------------------------------------------------------!

			! PUT only works for serial code.
			if ( any(tag(:).eq.PUT_thermo) ) call evaluate_U_PUT  

			! Warn user of lack of development
			if ( any(tag(:).eq.z_thermo) ) then
				if (mod(iter,tplot).eq.0 .and. irank .eq. iroot) then
				print*, 'Warning: thermostatting only in z direction. This &
				       & has not been fully developed and requires checking.'
				end if
			end if

            allocate(mv2sum(nthermo))
            allocate(thermostatnp(nthermo))
			mv2sum = 0.d0
			thermostatnp = 0.d0
			do n = 1, np
				if ( tag(n) .eq. PUT_thermo ) then
					vel(:) = v(:,n) - U(:,n) + 0.5d0*a(:,n)*delta_t
				else if ( any( thermo_tags .eq. tag(n) ) ) then
					vel(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
				else
					! Don't include non-thermostatted molecules in calculation
					cycle
				end if

                !multiple thermostat regions
                tr = get_therm_region(n)
				mv2sum(tr) = mv2sum(tr) + mass(n) * dot_product(vel,vel)
				thermostatnp(tr) = thermostatnp(tr) + mass(n)                

			enddo

			!Obtain global sums for all parameters
			call globalSum(thermostatnp(:), size(thermostatnp))
			call globalSum(mv2sum(:), size(mv2sum))

			!Nose Hoover thermostat coefficients
            allocate(Q(nthermo))
            allocate(dzeta_dt(nthermo))
            if (.not.allocated(zeta)) then
                allocate(zeta(nthermo))
                zeta = 0.d0
            endif

            do n=1,nthermo
			    Q(n)        = 0.1*thermostatnp(n) * delta_t
                if (thermostatnp(n) .gt. 0.d0) then
    			    dzeta_dt(n) = (mv2sum(n) - (nd*thermostatnp(n) + 1) & 
                                   *thermostattemperature(n)) / Q(n)
                else
                    dzeta_dt(n) = 0.d0
                endif
			    zeta(n) 	 = zeta(n) + delta_t*dzeta_dt(n)
			    bscale(n)	 = 1.d0/(1.d0+0.5d0*delta_t*zeta(n))
			    ascale(n)	 = (1.d0-0.5d0*delta_t*zeta(n))*bscale(n)
            enddo

		else

			!Reduces to the un-thermostatted equations
			ascale(:) = 1.d0
			bscale(:) = 1.d0
		endif

		!Step through each particle n
		do n = 1,np

			select case (tag(n))
			case (free)
				!Leapfrog mean velocity calculated here at v(t+0.5delta_t) = v(t-0.5delta_t) + a*delta_t 
				v(:,n) = v(:,n) + delta_t * a(:,n) 	!Velocity calculated from acceleration
				!Leapfrog mean position calculated here at r(t+delta_t) = r(t) + v(t+0.5delta_t)*delta_t
				r(:,n) = r(:,n) + delta_t * v(:,n)	!Position calculated from velocity
			case (fixed)
				!Fixed Molecules - no movement r(t+dt) = r(t)
			case (fixed_slide)
				!Fixed with constant sliding speed
				r(:,n) = r(:,n) + delta_t*slidev(:,n)	!Position calculated from velocity
			case (teth)
				!Tethered molecules
				call tether_force(n)
				v(:,n) = v(:,n) + delta_t * a(:,n) 	!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n)	!Position calculated from velocity
			case (thermo)
				!Nose Hoover Thermostatted Molecule
                                tr = get_therm_region(n)
				v(1,n) = v(1,n)*ascale(tr) + a(1,n)*delta_t*bscale(tr)
				r(1,n) = r(1,n)    +     v(1,n)*delta_t			
				v(2,n) = v(2,n)*ascale(tr) + a(2,n)*delta_t*bscale(tr)
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t				
				v(3,n) = v(3,n)*ascale(tr) + a(3,n)*delta_t*bscale(tr)
				r(3,n) = r(3,n)    +     v(3,n)*delta_t
			case (teth_thermo)
                                tr = get_therm_region(n)
				!Thermostatted Tethered molecules unfixed with no sliding velocity
				call tether_force(n)
				v(1,n) = v(1,n)*ascale(tr) + a(1,n)*delta_t*bscale(tr)
				r(1,n) = r(1,n)    +     v(1,n)*delta_t			
				v(2,n) = v(2,n)*ascale(tr) + a(2,n)*delta_t*bscale(tr)
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t				
				v(3,n) = v(3,n)*ascale(tr) + a(3,n)*delta_t*bscale(tr)
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	
				!write(1000000+iter*100+irank,'(8f10.5)'), globalise(r(:,n)),v(:,n),thermostattemperature,mv2sum/(nd*thermostatnp)
			case (teth_slide)
				!Tethered molecules with sliding velocity
				call tether_force(n)
				v(:,n) = v(:,n) + delta_t * a(:,n) 							!Velocity calculated from acceleration
				r(:,n) = r(:,n) + delta_t * v(:,n) + delta_t*slidev(:,n)	!Position calculated from velocity+slidevel
			case (teth_thermo_slide)
                                tr = get_therm_region(n)
				!Thermostatted Tethered molecules unfixed with sliding velocity
				call tether_force(n)
				v(1,n) = v(1,n)*ascale(tr) + a(1,n)*delta_t*bscale(tr)
				r(1,n) = r(1,n)    +     v(1,n)*delta_t	+ slidev(1,n)*delta_t		
				v(2,n) = v(2,n)*ascale(tr) + a(2,n)*delta_t*bscale(tr)
				r(2,n) = r(2,n)    + 	 v(2,n)*delta_t	+ slidev(2,n)*delta_t			
				v(3,n) = v(3,n)*ascale(tr) + a(3,n)*delta_t*bscale(tr)
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	+ slidev(3,n)*delta_t
			case (PUT_thermo)
                                tr = get_therm_region(n)
				!Profile unbiased thermostat (Nose-Hoover)
	        	v(:,n) = v(:,n)*ascale(tr) + (a(:,n)+zeta*U(:,n))*delta_t*bscale(tr)
				r(:,n) = r(:,n)        + v(:,n)*delta_t			
			case (z_thermo)
                tr = get_therm_region(n)
				!Thermostat in the z direction only (Nose-Hoover)
				v(1,n) = v(1,n) + delta_t * a(1,n) 	
				r(1,n) = r(1,n) + delta_t * v(1,n)	
				v(2,n) = v(2,n) + delta_t * a(2,n) 	
				r(2,n) = r(2,n) + delta_t * v(2,n)	
				v(3,n) = v(3,n)*ascale(tr) + a(3,n)*delta_t*bscale(tr)
				r(3,n) = r(3,n)    +     v(3,n)*delta_t	
			case (cyl_teth_thermo_rotate)

				! Tether force
				call tether_force(n)

				! Thermostatted move of peculiar vel/pos
				v(:,n) = v(:,n)*ascale + a(:,n)*delta_t*bscale
				r(:,n) = r(:,n)        + v(:,n)*delta_t			
	
				! Ad-hoc rotate cylinder molecules
				rglob(:) = globalise(r(:,n))
				rpol(:)  = cpolariser(rglob(:))
				rpol(2)  = rpol(2) + omega*delta_t
				rglob(:) = cartesianiser(rpol(:))
				r(:,n)   = localise(rglob(:))

				! Ad-hoc rotate cylinder tether sites
				rglob(:) = globalise(rtether(:,n))
				rpol(:)  = cpolariser(rglob(:))
				rpol(2)  = rpol(2) + omega*delta_t
				rglob(:) = cartesianiser(rpol(:))
				rtether(:,n) = localise(rglob(:))

			case default
                print'(3(a,i8))', ' Error in simulation_move_lfv at Iter =', iter, ' Mol no = ', n, ' Tag number = ', tag(n)
				call error_abort("Invalid molecular Tag")
			end select

		enddo

		if (config_special_case .eq. 'rotate_cylinders') call update_omega

	end subroutine simulation_move_particles_lfv_tag

	subroutine update_omega
		use concentric_cylinders, only: omega_f, omega_i, omega, &
		                                omega_rampiters
		use computational_constants_MD, only: iter, initialstep, delta_t
		implicit none

		if ( iter - initialstep .lt. omega_rampiters ) then
			
			omega = omega_i + (omega_f-omega_i)*0.5* &
			                  (1 - cos(pi*(iter-initialstep)/omega_rampiters))
		
		else 
			
			omega = omega_f
		
		end if

	end subroutine update_omega


	!-----------------------------------------------------------------------
	!Evaluate "true" positions and velocities without periodic wrapping
	subroutine simulation_move_particles_true_lfv
		use shear_info_MD, only: le_sp, le_sv, le_sd
		implicit none

        if (any(periodic .gt. 1)) then
		vtrue = v
		do n=1,np
			vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
			rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n)
		end do
        else
		    do n=1,np
			    rtrue(:,n) = rtrue(:,n) + delta_t*v(:,n)
		    end do
        endif

	end subroutine simulation_move_particles_true_lfv

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
! Specular walls for cylinder setup
! Author: David Trevelyan 2013
subroutine specular_walls_cylinders
    use concentric_cylinders
    use computational_constants_MD, only: delta_rneighbr 
    use physical_constants_MD, only: np
    use arrays_MD, only: r, v
    use messenger, only: localise, globalise
    use librarymod, only: cpolariser, cartesianiser, cpolarisev, cartesianisev
    implicit none

    real(kind(0.d0)) :: rcpol(3), vcpol(3), rglob(3)
    real(kind(0.d0)) :: tol, rr
    integer :: n,cyl
        
    tol = delta_rneighbr

    do n=1,np

        rglob(:) = globalise(r(:,n))
        rcpol(:) = cpolariser(rglob(:))
        vcpol(:) = cpolarisev(v(:,n),rcpol(2))
        rr = rcpol(1)

        if ( rr .lt. r_oi + tol) cyl = cyl_inner
        if ( rr .gt. r_io - tol) cyl = cyl_outer

        if ( cyl .eq. cyl_outer ) then
    
            call speculate(rcpol,vcpol,r_oo,r_io)
    
        else if ( cyl .eq. cyl_inner ) then

            call speculate(rcpol,vcpol,r_oi,r_ii)
        
        end if 
    
        rglob(:) = cartesianiser(rcpol(:))
        r(:,n) = localise(rglob(:))
        v(:,n) = cartesianisev(vcpol(:),rcpol(2))

    end do

contains

    subroutine speculate(rcpol,vcpol,r_o,r_i)

        real(kind(0.d0)), intent(in)    :: r_o, r_i
        real(kind(0.d0)), intent(inout) :: rcpol(3),vcpol(3)

        real(kind(0.d0)) :: dr

        if ( rcpol(1) .gt. r_o ) then

            dr = rcpol(1) - r_o	
            rcpol(1) = rcpol(1) - 2.0*dr
            vcpol(1) = -vcpol(1)

        else if ( rcpol(1) .lt. r_i ) then

            dr = r_i - rcpol(1)	
            rcpol(1) = rcpol(1) + 2.0*dr
            vcpol(1) = -vcpol(1)

        end if			

    end subroutine speculate

end subroutine specular_walls_cylinders

!subroutine specular_radial_wall(molno, r_prev, r_in, r_out) 
!Author: Musab Khawaja 2012
!	use module_molecule_properties
!	use arrays_MD
!    use interfaces
!	implicit none
!
!	integer                        :: molno, normal_dirn
!	real(kind(0.d0))               :: r_in, r_out, wall, perp_dist, vr_mag, tol
!    real(kind(0.d0)), dimension(3) :: r_prev, r_prev_pol, r_prev_glob, r_pol, &
!                                      r_glob, z_hat, theta_hat, normal, r2_X, &
!                                      rX, vr, vtheta, perp_vec
!
!    z_hat = (/ 0.d0, 0.d0, 1.d0 /)
!    tol = 0.0000001d0
!
!    r_prev_glob = globalize(r_prev)
!    r_prev_pol = cpolarize(r_prev_glob)
!    r_glob = globalize(r(:,molno))
!    r_pol = cpolarize(r_glob)
!	
!    !check if molecule will leave cell: 
!    if (r_pol(1) < r_in .or. r_pol(1) > r_out) then
!      !print*, 'specular collision'
!      if (r_prev_pol(1) < r_in .or. r_prev_pol(1) > r_out) print*, 'Already out'
!      
!      if (r_pol(1) < r_in) then 
!		  wall = r_in         !inner wall:
!          normal_dirn = 1   !  normal should point outwards
!	  else
!		  wall = r_out      !outer wall:
!          normal_dirn = -1  !  normal should point inwards
!	  endif
!
!      rX = find_cylinder_intersection(r_prev_glob,r_glob,wall,tol)
!      normal(:) = normal_dirn*normalise(rX(:))
!
!      ! Resolve velocity relative to normal and tangent
!!      theta_hat(:) = crossprod(z_hat(:), normal(:))
!      vr_mag = dot_product(v(:,molno), normal(:))
!      vtheta(:) = dot_product(v(:,molno), theta_hat(:))*theta_hat(:)
!      ! Force vr to point in normal direction
!      vr(1) = abs(vr_mag)*normal(1)
!      vr(2) = abs(vr_mag)*normal(2)
!
!      ! Reflect normal velocity.
!      v(1,molno) = vr(1) + vtheta(1)
!      v(2,molno) = vr(2) + vtheta(2)
!
!      ! Relative displacement vector of new point r2 from intersection
!      r2_X(1) = r_glob(1) - rX(1)
!      r2_X(2) = r_glob(2) - rX(2)
!      r2_X(3) = 0.d0
!      
!      ! Shortest perpendicular distance of r2_X from tangent
!      perp_vec(:) = crossprod(r2_X(:), theta_hat(:))
!      perp_dist = magnitude3(perp_vec(:))
!
!      !reflected new position 
!      r_glob(1) = r_glob(1) + 2.d0*perp_dist*normal(1)
!      r_glob(2) = r_glob(2) + 2.d0*perp_dist*normal(2)
!
!      !convert back to local processor co-ordinates.
!      r(:,molno) = localize(r_glob)	
!
!	endif
!
!end subroutine specular_radial_wall

! New routine written by DT during debugging, now slightly cleaner syntactically 
! and safer with globalise/localise routines
subroutine specular_flat_walls(ixyz)
    use physical_constants_MD, only: np
    use computational_constants_MD, only: globaldomain, tether_tags
    use arrays_MD, only: r, v, tag
    use messenger, only: localise, globalise
    use boundary_MD, only: specular_wall
    implicit none

    integer, intent(in) :: ixyz 

    real(kind(0.d0)) :: spec_top, spec_bot
    real(kind(0.d0)) :: rglob(3), molvel(3)
    integer :: n

    spec_top =  globaldomain(ixyz)/2.d0 - specular_wall(ixyz)
    spec_bot = -globaldomain(ixyz)/2.d0 + specular_wall(ixyz)

    do n = 1, np

        if (any(tag(n) .eq. tether_tags)) cycle

        rglob(:) = globalise(r(:,n))
        molvel(:) = v(:,n) 

        call speculate(rglob,molvel)
    
        r(:, n) = localise(rglob)
        v(:, n) = molvel(:) 
    
    end do

contains

    subroutine speculate(rin, vin)

        real(kind(0.d0)), intent(inout) :: rin(3), vin(3)
        
        if (rin(ixyz) .gt. spec_top) then
            
            rin(ixyz) = rin(ixyz) - 2.d0*(rin(ixyz) - spec_top)
            vin(ixyz) = vin(ixyz) * (-1.d0)

        else if (rin(ixyz) .lt. spec_bot) then

            rin(ixyz) = rin(ixyz) + 2.d0*(spec_top - rin(ixyz))
            vin(ixyz) = vin(ixyz) * (-1.d0)

        end if

    end subroutine speculate

end subroutine specular_flat_walls


!--------------------------------------------------------------------------------------
!subroutine specular_flat_wall(dir, spec_pos, flag)
!    use module_molecule_properties
!    use arrays_MD
!    use interfaces
!    use librarymod, only :  Maxwell_Boltzmann_vel3
!    implicit none
!
!    integer,intent(in) :: dir
!    integer,intent(in) :: flag
!    real(kind(0.d0)),intent(in)	:: spec_pos
!
!    real(kind(0.d0)), dimension(3) :: r_glob
!    integer                        :: n, ixyz
!    integer                        :: normal
!    real(kind(0.d0))               :: newxd
!
!    do n = 1,np
!
!        !Get position in global co-ordinates
!        newxd = 0.d0
!        r_glob(1) = r(1,n) - (halfdomain(1)*(npx-1)) + domain(1)*(iblock-1)
!        r_glob(2) = r(2,n) - (halfdomain(2)*(npy-1)) + domain(2)*(jblock-1)
!        r_glob(3) = r(3,n) - (halfdomain(3)*(npz-1)) + domain(3)*(kblock-1)
!
!        ! Skip tethered mols
!        if (any(tag(n) .eq. tether_tags)) cycle
!
!        if (abs(r_glob(dir)) .gt. spec_pos) then
!
!            !Get normal direction of molecule by checking if top or bottom
!            if ( r_glob(dir) .lt. 0) then
!                normal = 1   !normal should point outwards
!            else
!                normal = -1  !normal should point inwards
!            endif
!
!            ! Reflect normal velocity.
!            v(dir,n) = normal*abs(v(dir,n))
!            !Calculate distance over the spectral barrier
!            newxd = abs(r_glob(dir)) - spec_pos
!            !Move molecule same distance back into the domain on other side of spectral barrier
!            r(dir,n) = r(dir,n) + normal*2.d0*newxd
!
!            select case(flag)
!            case(0)
!                !Do nothing - control case
!            case(1)
!                !Pick specified temperature and velocity from Maxwell Boltzmann style distribution
!                v(:,n) = Maxwell_Boltzmann_vel3(thermostattemperature,wallslidev(:)*sign(1.d0,r_glob(dir)))
!!                    do ixyz = 1,3
!!                        v(ixyz,n) = Maxwell_Boltzmann_vel(thermostattemperature,wallslidev(ixyz)*sign(1.d0,r_glob(dir)))
!!                    enddo
!                !if (newxd .gt. 0.0000d0) print'(a,2i7,7f10.5)', 'specular wall', iter, n, v(:,n), thermostattemperature,wallslidev(:)
!            end select
!
!            !write(7777,'(4i8,4f10.5)'), irank,dir,normal, n,newxd,r_glob(dir), r(dir,n),spec_pos
!
!        endif
!
!        !if (newxd .gt. 0.0000d0) print'(a,i8,6f10.5,i8)', 'Greater than spec_pos', n,r_glob(dir),abs(r_glob(dir)),spec_pos,r(dir,n),newxd, periodic(dir)
!
!    end do
!
!end subroutine specular_flat_wall 
