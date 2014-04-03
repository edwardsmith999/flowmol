!=========================================================================================!
!-----------------------------------------------------------------------------------------!
!
!                             M O V E    P A R T I C L E S
!
!       Move particles as a result of forces using the velocity-Verlet algorithm
!
!-----------------------------------------------------------------------------------------!

module module_move_particles_vv

	use interfaces
	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD

end module module_move_particles_vv

!=========================================================================================!
! simulation_move_particles_vv(pass_num)
! authors: David Trevelyan & Ed Smith
!
! description:
!     simulation_move_particles_vv contains both parts of the velocity-Verlet algorithm
!     in which particle positions or velocities are updated. For a full description of
!     the basic algorithm, please see "Computer Simulation of Liquids" by Allen &
!     Tildesley. The subroutine must first be called with input pass_num equal to 1 and,
!     after applying boundary conditions and a force computation, should be subsequently
!     called with pass_num equal to 2. 
!
!     first pass:
!        - r(t+dt)   = r(t) + v(t)*dt + 0.5*a(t)*dt^2
!        - v(t+dt/2) = v(t) + 0.5*a(t)*dt
!     (apply BCs, compute forces)
!     second pass:
!        - v(t+dt)   = v(t+dt/2) + 0.5*a(t+dt)*dt
!
!     computations that are only required by extended system ensembles (for example,
!     evaluating the time derivative of the damping parameter "zeta" in the Nosé-Hoover
!     algorithm) are, in the interests of clarity, computed in subroutines contained
!     within simulation_move_particles_vv. 
!
!-----------------------------------------------------------------------------------------!
subroutine simulation_move_particles_vv(pass_num)
	use module_move_particles_vv
	implicit none

	integer                :: n,i
	integer, intent(in)    :: pass_num
	double precision       :: alpha
	double precision       :: zeta_old
	double precision, save :: dzeta_dt
	double precision, save :: zeta=0.d0
	double precision, dimension(nd,np) :: v_old
	double precision, dimension(nd,np) :: vrelsum
	double precision, dimension(nd,np) :: Utrue


	!--------First half of velocity-Verlet algorithm. Finds r(t+dt) and v(t+dt/2).--------!
	if (pass_num.eq.1) then


		select case(ensemble)
		case(nve)
			if (rtrue_flag.eq.1) call simulation_move_particles_true_vv
			do n=1,np
				r(:,n)     = r(:,n)     + delta_t*v(:,n) + 0.5d0*(delta_t**2.d0)*a(:,n)
				v(:,n)     = v(:,n)     + 0.5d0*delta_t*a(:,n)
			end do

		case(nvt_NH)
			call evaluate_dzeta_dt
			if (rtrue_flag.eq.1) call simulation_move_particles_true_vv
			do n=1,np
				r(:,n)     = r(:,n)     + delta_t*v(:,n) + 0.5d0*(delta_t**2.d0)*(a(:,n)-zeta*v(:,n))
				v(:,n)     = v(:,n)     + 0.5d0*delta_t*(a(:,n)-zeta*v(:,n))
			end do

		case(nvt_GIK)
			if (rtrue_flag.eq.1) call simulation_move_particles_true_vv
			do n=1,np
				r(:,n)     = r(:,n)     + delta_t*v(:,n) + 0.5d0*(delta_t**2.d0)*(a(:,n)-zeta*v(:,n))
				v(:,n)     = v(:,n)     + 0.5d0*delta_t*(a(:,n)-zeta*v(:,n))
			end do			

		case(nvt_PUT_NH)	
			call evaluate_U_PUT
			call evaluate_dzeta_dt_PUT
			if (rtrue_flag.eq.1) call simulation_move_particles_true_vv
			do n=1,np
				r(:,n)     = r(:,n)     + delta_t*v(:,n) + 0.5d0*(delta_t**2.d0)*(a(:,n)-zeta*(v(:,n)-U(:,n)))
				v(:,n)     = v(:,n)     + 0.5d0*delta_t*(a(:,n)-zeta*(v(:,n)-U(:,n)))
			end do
		
		case(nvt_pwa_NH)
			call evaluate_pwa_terms_pwaNH
			if (rtrue_flag.eq.1) call simulation_move_particles_true_vv
			do n=1,np
				v(:,n)     = v(:,n) + 0.5d0*delta_t*(a(:,n) - zeta*vrelsum(:,n))
				r(:,n)     = r(:,n) + delta_t*v(:,n)
			end do

		case(nvt_DPD)
			if (iter .eq. initialstep) then
				call evaluate_DPD(1)
			else
				call messenger_updateborders(0)    
				call evaluate_DPD(0)
			endif
			if (rtrue_flag.eq.1) call simulation_move_particles_true_vv
			do n=1,np
				v(:,n)     = v(:,n)     + 0.5d0*delta_t*(a(:,n)+aD(:,n)+aR(:,n))
				r(:,n)     = r(:,n)     + delta_t*v(:,n) 
			end do
			!Regenerate random number for next iteration
			call random_number(theta)
			theta = sqrt(3.d0)*(2.d0*theta-1.d0)

		case(tag_move)
			call error_abort('Tag mode for velocity-Verlet not yet implemented.')
	
		case default
			call error_abort('Unrecognised ensemble, stopping.')

		end select	
	
	!-------Second half of velocity-Verlet algorithm.-------------------------------------!
	else if (pass_num.eq.2) then

		select case(ensemble)

			case(nve)
				do n=1,np
					v(:,n) = v(:,n) + 0.5d0*delta_t*a(:,n)
				end do
	
			case(nvt_NH)
				zeta  = zeta + delta_t*dzeta_dt
				alpha = 1.d0 + 0.5d0*delta_t*zeta
				do n=1,np
					v(:,n) = (v(:,n) + 0.5d0*delta_t*a(:,n))/alpha
				end do
			
			case(nvt_GIK)
				v_old = v
				do i=1,3
					call evaluate_zeta_GIK
					alpha = 1.d0 + 0.5d0*delta_t*zeta
					do n=1,np
						v(:,n) = (v_old(:,n) + 0.5d0*delta_t*a(:,n))/alpha
					end do
				end do

			case(nvt_PUT_NH)
				call evaluate_U_PUT
				zeta  = zeta + delta_t*dzeta_dt
				alpha = 1.d0 + 0.5d0*delta_t*zeta
				do n=1,np
					v(:,n) = (v(:,n) + 0.5d0*delta_t*(a(:,n)+zeta*U(:,n)))/alpha
				end do

			case(nvt_pwa_NH)
				call evaluate_pwa_terms_pwaNH
				v_old = v
				zeta_old = zeta + 0.5d0*delta_t*dzeta_dt
				do i = 1,5
					call messenger_updateborders(0)
					call evaluate_pwa_terms_pwaNH
					do n=1,np
						v(:,n) = v_old(:,n) + 0.5d0*delta_t*(a(:,n)- zeta*vrelsum(:,n))
					end do
					zeta = zeta_old + 0.5d0*delta_t*dzeta_dt
				end do

			case(nvt_DPD)	
				call evaluate_DPD(1)
				do n=1,np
					v(:,n) = v(:,n) + 0.5d0*delta_t*(a(:,n)+aD(:,n)+aR(:,n))
				end do

			case(tag_move)
				call error_abort('Tag mode for velocity-Verlet not yet implemented.')
			
			case default
				call error_abort('Unrecognised ensemble, stopping.')

		end select
			
		if (rtrue_flag.eq.1) call simulation_move_particles_true_vv	
		
		simtime = simtime + delta_t
		
	endif

contains

	!-----------------------------------------------------------------------
	!Evaluate "true" positions and velocities without periodic wrapping
	subroutine simulation_move_particles_true_vv
		use shear_info_MD, only: le_sv, le_sp, le_sd
		implicit none

		vtrue = v

		if (pass_num .eq. 1) then

			select case (ensemble)
			case (nve)

				do n=1,np
					vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
					rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n) &
					                            + 0.5d0*(delta_t**2.d0)*a(:,n)
					
				end do

			case (nvt_NH)

				do n=1,np
					vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
					rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n) &
					                            + 0.5d0*(delta_t**2.d0)*(a(:,n)-zeta*vtrue(:,n))
				end do

			case (nvt_GIK)
				
				do n=1,np
					vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
					rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n) &
					                            + 0.5d0*(delta_t**2.d0)*(a(:,n)-zeta*vtrue(:,n))
				end do
	
			case(nvt_PUT_NH)	
				do n=1,np
					vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
					Utrue(:,n)     = U(:,n)
					Utrue(le_sd,n) = U(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
					rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n) + &
					                 0.5d0*(delta_t**2.d0)*(a(:,n)-zeta*(v(:,n)-U(:,n)))
				end do
			
			case(nvt_pwa_NH)
				do n=1,np
					vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
					rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n)
				end do

			case(nvt_DPD)
				do n=1,np
					vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
					rtrue(:,n)     = rtrue(:,n) + delta_t*vtrue(:,n) 
				end do

			case default

			end select
		
		else if (pass_num .eq.2) then
	
			do n=1,np
				vtrue(le_sd,n) = v(le_sd,n) + anint(rtrue(le_sp,n)/domain(le_sp))*le_sv
			end do
	
		end if

	end subroutine simulation_move_particles_true_vv

	!---------------------------------------------------------------------------
	!Evaluate dzeta_dt for global Nose-Hoover equations of motion
	subroutine evaluate_dzeta_dt	
		implicit none

		double precision :: v2sum
		double precision :: Q
	
		v2sum=0.d0
		do n=1,np
			v2sum = v2sum + dot_product(v(:,n),v(:,n))
		end do
		call globalSum(v2sum)

		Q        = globalnp*delta_t
		dzeta_dt = (v2sum - (globalnp*nd+1)*inputtemperature)/Q

	end subroutine evaluate_dzeta_dt

	!---------------------------------------------------------------------------
	!Evaluate zeta for global Gaussian isokinetic thermostat
	subroutine evaluate_zeta_GIK
		implicit none

		double precision :: avsum
		double precision :: v2sum

		avsum = 0.d0
		v2sum = 0.d0
		do n=1,np
			avsum = avsum + dot_product(a(:,n),v(:,n))
			v2sum = v2sum + dot_product(v(:,n),v(:,n))
		end do
		call globalSum(avsum)
		call globalSum(v2sum)

		zeta = avsum/v2sum
		
	end subroutine evaluate_zeta_GIK
	
	!---------------------------------------------------------------------------
	!Evaluate dzeta_dt for profile unbiased Nose-Hoover equations of motion
	subroutine evaluate_dzeta_dt_PUT
		use calculated_properties_MD,  only: nbins, get_mass_slices, get_velo_slices
		use shear_info_MD,             only: le_sp	
		implicit none
		
		double precision :: pec_v2sum
		double precision :: Q
		double precision, dimension(nd) :: pec_v

		pec_v2sum = 0.d0
		do n=1,np
			pec_v(:)  = v(:,n) - U(:,n)                                         ! PUT: Find peculiar velocity
			pec_v2sum = pec_v2sum + dot_product(pec_v,pec_v)                    ! PUT: Sum peculiar velocities squared
		end do
		call globalSum(pec_v2sum)

		Q = globalnp*delta_t
		dzeta_dt = (pec_v2sum - (globalnp*nd+1)*inputtemperature)/Q

	end subroutine evaluate_dzeta_dt_PUT

	!---------------------------------------------------------------------------
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
		v_slice = get_velo_slices(le_sp)                                  ! PUT: Get total velocity in all slices (note that on the second pass this is half a timestep behind.)
		do slicebin=1,nbins(le_sp)                                        ! PUT: Loop through all slices
			v_avg(slicebin,:) = v_slice(slicebin,:)/m_slice(slicebin)           ! PUT: average velocity
		end do
		
		do n=1,np
			slicebin = ceiling((r(le_sp,n)+halfdomain(le_sp))/&
			                    slicebinsize(le_sp))
			if (slicebin > nbins(le_sp)) slicebin = nbins(le_sp)    ! PUT: Prevent out-of-range values
			if (slicebin < 1) slicebin = 1                                      ! PUT: Prevent out-of-range values
			U(:,n) = v_avg(slicebin,:)
		end do

		deallocate(m_slice)
		deallocate(v_slice)
		deallocate(v_avg)

	end subroutine evaluate_U_PUT

	!----------------------------------------------------------------------------
	!Evaluate pairwise terms for Allen & Schmid thermostat
	subroutine evaluate_pwa_terms_pwaNH_AP
		use linked_list
		implicit none
	
		integer :: i,j
		double precision :: rij2,wsq,vr,Q
		double precision, dimension(nd) :: ri,rj,rij,rijhat
		double precision, dimension(nd) :: vi,vj,vij

		vrelsum = 0.d0
		dzeta_dt = 0.d0
		Q = np*0.02*delta_t

		do i=1,np
			
			ri(:) = r(:,i)
			vi(:) = v(:,i)
	
			do j=i+1,np
	
				rj(:)     = r(:,j)
				rij(:)    = ri(:) - rj(:)
				rij2      = dot_product(rij,rij)
				if (rij2.ge.rcutoff2) cycle
				vj(:)     = v(:,j)
				vij(:)    = vi(:) - vj(:)
				rijhat(:) = rij(:)/sqrt(rij2)
				wsq       = (1.d0-(sqrt(rij2)/rcutoff))*(1.d0-(sqrt(rij2)/rcutoff))
				vr        = dot_product(vij,rijhat)

				vrelsum(:,i) = vrelsum(:,i) + wsq*vr*rijhat(:)
				vrelsum(:,j) = vrelsum(:,j) - wsq*vr*rijhat(:)
				
				dzeta_dt = dzeta_dt + wsq*(vr**2.d0 - inputtemperature*2.d0)/Q 

			end do

		end do
	
	end subroutine evaluate_pwa_terms_pwaNH_AP

	!----------------------------------------------------------------------------
	!Evaluate pairwise terms for Allen & Schmid thermostat
	subroutine evaluate_pwa_terms_pwaNH
		use linked_list
		implicit none
	
		integer :: noneighbrs
		integer :: j,molnoi,molnoj
		double precision :: rij2,wsq,vr,Q
		double precision, dimension(nd) :: ri,rj,rij,rijhat
		double precision, dimension(nd) :: vi,vj,vij
		type(neighbrnode), pointer :: old, current

		vrelsum = 0.d0
		dzeta_dt = 0.d0
		Q = np*0.02*delta_t

		do molnoi=1,np
 
	    	noneighbrs = neighbour%noneighbrs(molnoi)   !Determine number of elements in neighbourlist
			old        => neighbour%head(molnoi)%point  !Set old to head of neighbour list
			ri(:)      = r(:,molnoi)
			vi(:)      = v(:,molnoi)
	
			do j=1,noneighbrs
	
				molnoj    = old%molnoj
				if (molnoj.eq.molnoi) call error_abort("Self interaction in pwaNH thermostat")
				rj(:)     = r(:,molnoj)
				vj(:)     = v(:,molnoj)
				rij(:)    = ri(:) - rj(:)
				vij(:)    = vi(:) - vj(:)
				rij2      = dot_product(rij,rij)
				rijhat(:) = rij(:)/sqrt(rij2)
				wsq       = (1.d0-(sqrt(rij2)/rcutoff))*(1.d0-(sqrt(rij2)/rcutoff))
				if (rij2.ge.rcutoff2) wsq = 0.d0
				vr        = dot_product(vij,rijhat)

				vrelsum(:,molnoi) = vrelsum(:,molnoi) + wsq*vr*rijhat(:)
				vrelsum(:,molnoj) = vrelsum(:,molnoj) - wsq*vr*rijhat(:)
				
				dzeta_dt = dzeta_dt + wsq*(vr**2.d0 - inputtemperature*2.d0)/Q 

				current => old	
				old => current%next

			end do

		end do
		
		nullify(current)
		nullify(old)
	
	end subroutine evaluate_pwa_terms_pwaNH

	!----------------------------------------------------------------------------
	!Evaluate pairwise terms for DPD thermostat by 
	!Soddemann, Dunweg an Kremer Phys Rev E 68, 046702 (2003)
	!From this paper : typical 0.5 < zeta < 1.5 
	!Random numbers do not need to be Gaussian

	subroutine evaluate_DPD(flag)
		use interfaces
		use linked_list
		implicit none
	
		integer,intent(in)				:: flag
		integer 						:: noneighbrs
		integer 						:: j,molnoi,molnoj
		double precision 				:: rij2,vr,wR,wD,sigma
		double precision, dimension(nd) :: ri,rj,rij,rijhat
		double precision, dimension(nd)	:: randi,randj,theta_ij
		double precision, dimension(nd)	:: vi,vj,vij
		type(neighbrnode), pointer 		:: old, current

		zeta = 10.d0; sigma = sqrt(2.d0*inputtemperature*zeta)

		!Either dissipative terms only or both dissipative and random terms
		select case(flag)
		case(0)
			!Evaluate dissipative terms only
			aD = 0.d0
			do molnoi=1,np
	 	    	noneighbrs = neighbour%noneighbrs(molnoi)   !Determine number of elements in neighbourlist
				old        => neighbour%head(molnoi)%point  !Set old to head of neighbour list
				ri(:)      = r(:,molnoi)
				vi(:)      = v(:,molnoi)
	
				do j=1,noneighbrs
					molnoj    = old%molnoj
					if (molnoj.eq.molnoi) call error_abort("Self interaction in DPD vv")	!Self interactions are unacceptable!
					rj(:)     = r(:,molnoj)
					rij(:)    = ri(:) - rj(:)
					rij2      = dot_product(rij,rij)

					!Thermostat force only local for molecules in cutoff range
					if (rij2 .lt. rcutoff2) then
						vj(:)     = v(:,molnoj)
						vij(:)    = vi(:)-vj(:)
						rijhat(:) = rij(:)/sqrt(rij2)
						wD        = -(1.d0-sqrt(rij2)/rcutoff)**2.d0
						vr        = dot_product(rijhat,vij)

						aD(:,molnoi) = aD(:,molnoi) + zeta*wD*vr*rijhat(:)
						aD(:,molnoj) = aD(:,molnoj) - zeta*wD*vr*rijhat(:)

					endif
					
					current => old	
					old => current%next

				enddo
			enddo
			
			nullify(current)
			nullify(old)
		case(1)
			!Evaluate random and dissipative terms
			aD = 0.d0; aR=0.d0
			do molnoi=1,np
	 	    	noneighbrs = neighbour%noneighbrs(molnoi)   !Determine number of elements in neighbourlist
				old        => neighbour%head(molnoi)%point  !Set old to head of neighbour list
				ri(:)      = r(:,molnoi)
				vi(:)      = v(:,molnoi)
				randi(:)   = theta(:,molnoi)
		
				do j=1,noneighbrs
					molnoj    = old%molnoj
					if (molnoj.eq.molnoi) stop "Self interaction in DPD vv"	!self interactions are unacceptable!
					rj(:)     = r(:,molnoj)
					rij(:)    = ri(:)-rj(:)
					rij2      = dot_product(rij,rij)

					!Thermostat force only local for molecules in cutoff range
					if (rij2 .lt. rcutoff2) then
						vj(:)     = v(:,molnoj)
						vij(:)    = vi(:)-vj(:)
						rijhat(:) = rij(:)/sqrt(rij2)
						wR        = 1.d0-sqrt(rij2)/rcutoff
						wD        = -wR**2.d0
						vr        = dot_product(rijhat,vij)
						randj(:)  = theta(:,molnoj)			

						!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
						!We have N random numbers theta_i and need to generate the theta_ij interaction
						!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -		
						!Random noise sum of theta_i+theta_j normalised by std - This gives the correct 
						!distribution but results in correlation as theta_i is used in multiple interactions
						!theta_ij(:)= (((randi(:)+randj(:))) - meantheta )/sqrt(vartheta)

						!Use product of random numbers - this no longer gives the same distribution
						!but satisfies the Weiner criterion, <W>=0 and <W W> = kron_delta*kron_delta*kron_delta
						theta_ij(:)= randi(:)*randj(:)

						!Divide by sqrt of dt so delta t can be used for aD and aR 
						!aR has sqrt of dt as required by ito calculus
						aD(:,molnoi) = aD(:,molnoi) + zeta*wD*vr*rijhat(:)
						aR(:,molnoi) = aR(:,molnoi) + sigma*wR*theta_ij(:)*rijhat(:)/sqrt(delta_t)

						aD(:,molnoj) = aD(:,molnoj) - zeta*wD*vr*rijhat(:)
						aR(:,molnoj) = aR(:,molnoj) - sigma*wR*theta_ij(:)*rijhat(:)/sqrt(delta_t)

					endif
					
					current => old	
					old => current%next

				enddo
			enddo
			
			nullify(current)
			nullify(old)

		case default
			call error_abort("Flag input to evaluate DPD incorrect")
		end select

	end subroutine evaluate_DPD

	!----------------------------------------------------------------------------
	!ALL PAIRS Evaluate pairwise terms for DPD thermostat by 
	!Soddemann, Dunweg an Kremer Phys Rev E 68, 046702 (2003)
	!From this paper : typical 0.5 < zeta < 1.5 
	!Random numbers do not need to be Gaussian

	subroutine evaluate_DPD_ap(flag)
		use linked_list
		implicit none
	
		integer,intent(in)				:: flag
		integer 						:: ixyz, tempi
		integer 						:: molnoi,molnoj
		double precision 				:: rij2,vr,wR,wD,sigma, temp, temp2
		double precision, dimension(nd)	:: ri,rj,rij,rijhat
		double precision, dimension(nd)	:: rand,randi,randj,theta_ij,meantheta,vartheta
		double precision, dimension(nd)	:: vi,vj,vij

		zeta = 10.d0; sigma = sqrt(2.d0*inputtemperature*zeta)
		select case(flag)
		case(0)
			!Evaluate dissipative terms only
			aD = 0.d0
			do molnoi = 1,np					!Step through each particle in list 
				ri = r(:,molnoi)         	!Retrieve ri
				vi = v(:,molnoi)

				do molnoj = molnoi+1,np				!Step through all j for each i
					rj = r(:,molnoj)				!Retrieve rj

					!Calculate rij using nearest image convention, i.e. if more than half
					! a domain betwen molecules, must be closer over periodic boundaries  
					rij2 = 0.d0
					do ixyz=1,nd
						rij(ixyz) = r(ixyz,molnoi) - r(ixyz,molnoj)          !Evaluate distance between particle i and j
		    				if (abs(rij(ixyz)) > halfdomain(ixyz)) then	
								rij(ixyz) = rij(ixyz) - sign(domain(ixyz),rij(ixyz)) 
							endif
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					if (rij2 < rcutoff2) then

						vj        = v(:,molnoj)
						vij(:)    = vi(:)-vj(:)
						rijhat(:) = rij(:)/sqrt(rij2)
						wD        = -(1.d0-sqrt(rij2)/rcutoff)**2.d0
						vr        = dot_product(rijhat,vij)

						aD(:,molnoi) = aD(:,molnoi) + zeta*wD*vr*rijhat(:)
						aD(:,molnoj) = aD(:,molnoj) - zeta*wD*vr*rijhat(:)
					endif

				enddo
			enddo

		case(1)
			!Evaluate random and dissipative terms
			aD = 0.d0; aR=0.d0
			!Calculate mean and variance of random number array
			meantheta(:) = sum(theta(:,1:np),1)/real(np,kind(0.d0))
			do ixyz = 1,nd
				vartheta(ixyz)  = sum((theta(ixyz,1:np)-meantheta(ixyz))**2.d0)/real(np,kind(0.d0))
			enddo
			!Get mean and variance for theta_ij = (theta_i + theta_j)
			!using mean_ij = mean_i + mean_j and var_ij = (var_i^2 + var_j^2)^0.5
			meantheta   = 2.d0*meantheta
			vartheta(:) = 2.d0*vartheta(:)
			temp = 0.d0; tempi = 0; temp2 = 0.d0

			do molnoi = 1,np					!Step through each particle in list 
				ri = r(:,molnoi)         	!Retrieve ri
				vi 		  = v(:,molnoi)
				randi(:)  = theta(:,molnoi)	

				do molnoj = molnoi+1,np				!Step through all j for each i
					rj = r(:,molnoj)				!Retrieve rj
					!Calculate rij using nearest image convention, i.e. if more than half
					! a domain betwen molecules, must be closer over periodic boundaries  
					rij2 = 0.d0
					do ixyz=1,nd
						rij(ixyz) = r(ixyz,molnoi) - r(ixyz,molnoj)          !Evaluate distance between particle i and j
		    				if (abs(rij(ixyz)) > halfdomain(ixyz)) then	
								rij(ixyz) = rij(ixyz) - sign(domain(ixyz),rij(ixyz)) 
							endif
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					if (rij2 < rcutoff2) then
						vi 		  = v(:,molnoi)
						vj        = v(:,molnoj)
						vij(:)    = vi(:)-vj(:)
						rijhat(:) = rij(:)/sqrt(rij2)
						wR        = 1.d0-sqrt(rij2)/rcutoff
						wD        = -wR**2.d0
						vr        = dot_product(rijhat,vij)
						randj(:)  = theta(:,molnoj)


						!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
						!We have N random numbers theta_i and need to generate the theta_ij interaction
						!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
						!Entirely regenerate the random numbers to prevent any correlation
						!Use random number to get mean zero, std 1
						theta_ij = sqrt(3.d0)*(2.d0*rand-1.d0)

						!Entirely regenerate the random numbers to prevent any correlation (else use theta_i array)
						!call random_number(randi); randi=sqrt(3.d0)*(2.d0*randi-1.d0)
						!call random_number(randj); randj=sqrt(3.d0)*(2.d0*randj-1.d0)

						!Random noise sum of theta_i+theta_j normalised by std - This gives the correct 
						!distribution but results in correlation as theta_i is used in multiple interactions
						!theta_ij(:)= (((randi(:)+randj(:))) - meantheta )/sqrt(vartheta)

						!Use product of random numbers - this no longer gives the same distribution
						!but satisfies the Weiner criterion, <W>=0 and <W W> = kron_delta*kron_delta*kron_delta	
						!theta_ij = randi(:)*randj(:)
						
						if (molnoi .eq. 150) write(200,'(3f10.5)') theta_ij(:)

						temp  = temp + theta_ij(1)
						temp2 = temp2+ theta_ij(1)**2
						tempi = tempi + 1

						!Divide by sqrt of dt so delta t can be used for aD and aR 
						!aR has sqrt of dt as required by ito calculus
						aD(:,molnoi) = aD(:,molnoi) + zeta*wD*vr*rijhat(:)
						aR(:,molnoi) = aR(:,molnoi) + sigma*wR*theta_ij(:)*rijhat(:)/sqrt(delta_t)

						aD(:,molnoj) = aD(:,molnoj) - zeta*wD*vr*rijhat(:)
						aR(:,molnoj) = aR(:,molnoj) - sigma*wR*theta_ij(:)*rijhat(:)/sqrt(delta_t)

					endif

				enddo
			enddo

		case default
			stop "Flag input to evaluate DPD incorrect"
		end select

		if(mod(iter,100) .eq. 0) then
			write(1000,'(i8,2f10.5)') iter, temp/tempi, sqrt(temp2/tempi - (temp/tempi)**2)
		!	print'(a,2f10.5)', 'Fluctuation dissipation required 2 zeros here:', & 
		!				sigma**2.d0-2.d0*inputtemperature*zeta, wD + wR**2.d0
		!	print'(a,2i5,2f18.5)', 'Sum of D and R Forces',iter,flag,sum(aD(:,1:np)), sum(aR(:,1:np))
		endif

	end subroutine evaluate_DPD_ap

end subroutine simulation_move_particles_vv
