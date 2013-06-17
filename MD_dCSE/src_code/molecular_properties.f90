!----------------------------------------------------------------------------------
!                                Molecular Properties
! Routines for setting and adjusting the properties of the molecules in the domain
!
!---------------------------------------------------------------------------------

module module_molecule_properties

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use librarymod

end module module_molecule_properties

!----------------------------------------------------------------------------------
!Create an array of molecular tags

	!For initialunitsize "a"
	!                [  o     o ]
	!a (1 cell size) |     o    ]  a/2 (distance between molcules)	
	!                |  o     o
	!                 __________]  a/4 (distance from bottom of domain)
	!
	!So use (0.20+0.5d0*mol_layers)*initialunitsize(ixyz)
subroutine setup_cylinder_tags
	use module_molecule_properties
	use concentric_cylinders
	use interfaces, only: error_abort
	use messenger, only: globalise
	implicit none

	integer :: n
	real(kind(0.d0)) :: rglob(3), rpol(3)

	do n = 1,np
		
		rglob = globalise(r(:,n))
		rpol  = cpolariser(rglob)
		
		if (rpol(1) .lt. r_oi) tag(n) = cyl_teth_thermo_rotate
		if (rpol(1) .gt. r_io) tag(n) = teth_thermo

	end do	

end subroutine setup_cylinder_tags

subroutine setup_location_tags
	use module_molecule_properties
	use interfaces, only: error_abort
	use messenger, only: globalise
	use computational_constants_MD, only : texture_type
	implicit none

	integer :: n
	real(kind(0.d0)) :: rglob(3)
	logical :: l_thermo
	logical :: l_teth
	logical :: l_fixed
	logical :: l_slide

	!Initialise all molecules free of thermostats, etc (i.e. tag free=0)
	tag(:) = free 

	if (ensemble .eq. tag_move) then

		!Setup fixed wall and location dependent thermostat tags
		do n = 1,np

			rglob(:) = globalise(r(:,n))

			l_thermo = get_tag_status(rglob,'thermo') 
			l_teth   = get_tag_status(rglob,'teth') 
			l_fixed  = get_tag_status(rglob,'fixed')
			l_slide  = get_tag_status(rglob,'slide')

			! Checks
			if ( l_teth .and. l_fixed ) then
				call error_abort("User attempted to fix AND tether a molecule. Aborting.")
			end if
			if ( l_thermo .and. l_fixed ) then
				call error_abort("User attempted to fix AND thermostat a molecule. Aborting.")
			end if

			! Teth
			if ( l_teth .and. .not. l_thermo .and. .not. l_slide  ) tag(n) = teth
			if ( l_teth .and. .not. l_thermo .and.       l_slide  ) tag(n) = teth_slide
			if ( l_teth .and.       l_thermo .and. .not. l_slide  ) tag(n) = teth_thermo
			if ( l_teth .and.       l_thermo .and.       l_slide  ) tag(n) = teth_thermo_slide

			! Fixed 
			if ( l_fixed .and. .not. l_slide ) tag(n) = fixed
			if ( l_fixed .and.       l_slide ) tag(n) = fixed_slide

			! Thermo only
			if ( l_thermo .and. .not. l_teth .and. .not. l_slide ) tag(n) = thermo

		enddo

		!Check if any molecules are thermostatted and switch on flag
		call get_tag_thermostat_activity(tag_thermostat_active)

	end if

contains

	function get_tag_status(rg,status_type) result(tag_status)
		implicit none

		real(kind(0.d0)) :: rg(3)   ! Global position
		character(*) :: status_type ! Specifies thermo, tethered, fixed, etc

		integer :: ixyz
		real(kind(0.d0)) :: bottom(3) ! Global position of bottom of domain
		real(kind(0.d0)) :: top(3)
		real(kind(0.d0)) :: tagdistbottom(3) ! Distance specified by user
		real(kind(0.d0)) :: tagdisttop(3)
		logical :: tag_status 

		bottom = (/ -globaldomain(1)/2.d0, -globaldomain(2)/2.d0, -globaldomain(3)/2.d0 /)
		top    = (/  globaldomain(1)/2.d0,  globaldomain(2)/2.d0,  globaldomain(3)/2.d0 /)

		select case (status_type)
		case ('thermo')
			tagdistbottom(:) = thermstatbottom(:)
			tagdisttop(:)    = thermstattop(:)
			!Thermostat complicated wall texture if specified, otherwise thermostat
			!is based on thermstatbottom/thermstattop
			if (texture_type .ne. 0 .and. texture_therm .eq. 1) then
				call wall_textures(texture_type,rg,tagdistbottom,tagdisttop)
			endif
		case ('teth')
			tagdistbottom(:) = tethereddistbottom(:)
			tagdisttop(:)    = tethereddisttop(:)
			!Apply a complicated wall texture if specified
			if (texture_type .ne. 0) call wall_textures(texture_type,rg,tagdistbottom,tagdisttop)
		case ('fixed')
			tagdistbottom(:) = fixdistbottom(:)
			tagdisttop(:)    = fixdisttop(:)
		case ('slide')
			tagdistbottom(:) = slidedistbottom(:)
			tagdisttop(:)    = slidedisttop(:)
		case default
			call error_abort("Unrecognised tag status type")
		end select

		tag_status = .false.
		! If within a tagged region then mark tag status as true and return
		do ixyz = 1,3
			if (rg(ixyz) .le. bottom(ixyz) + tagdistbottom(ixyz) .or. &
				rg(ixyz) .ge. top(ixyz)    - tagdisttop(ixyz)) then
				tag_status = .true.
				return

			end if	
		end do

	end function get_tag_status

end subroutine setup_location_tags

!----------------------------------------------------------------------------------
! Check if thermostat is active anywhere in world
subroutine get_tag_thermostat_activity(active)
	use physical_constants_MD, only: np
	use computational_constants_MD, only: thermo_tags, ensemble, tag_move
	use arrays_MD, only: tag
	implicit none

	logical, intent(out) :: active
	integer :: i, n

	if (ensemble.ne.tag_move) return 

	i = 0
	do n = 1, np
		if (any(tag(n).eq.thermo_tags)) then
			i = 1	
			exit
		end if
	end do
	call globalMaxInt(i)

	if (i.eq.1) then
		active = .true.
	else
		active = .false.
	end if
	
end subroutine get_tag_thermostat_activity

!----------------------------------------------------------------------------------

subroutine wall_textures(texture_type,rg,tagdistbottom,tagdisttop)
	use physical_constants_MD, only : pi,tethereddistbottom,tethereddisttop
	use computational_constants_MD, only : posts,roughness,converge_diverge,texture_intensity, &
										   globaldomain, cellsidelength,texture_therm,nh,halfdomain,ncells
	implicit none

	integer,intent(in)		:: texture_type
	double precision,dimension(3),intent(in) :: rg
	double precision,dimension(3),intent(out):: tagdistbottom,tagdisttop

	integer 				:: n,icell,jcell,kcell
	double precision		:: xlocation,ylocation,zlocation,rand,fraction_domain,postheight

	select case (texture_type)
	case(posts)

		postheight = 5.12d0
		tagdistbottom = tethereddistbottom
		ylocation = rg(2) + 0.5*globaldomain(2)

		! Post is above wall
		if ((ylocation .gt. tethereddistbottom(2) ) .and. & 
			(ylocation .lt. tethereddistbottom(2) + postheight)) then

			!Single strip in the middle of the domain
			if (rg(1) .gt. -postheight/2.d0 .and. &
				rg(1) .lt.  postheight/2.d0) then
				tagdistbottom(2) = tethereddistbottom(2) + postheight 
			else
				tagdistbottom(2) = tethereddistbottom(2)
			endif
		endif

	case(roughness)

		!rough wall
		call random_number(rand)
		tagdistbottom = 0.d0; tagdisttop=0.d0
		xlocation = rg(1)/globaldomain(1) + 0.5
		zlocation = rg(3)/globaldomain(3) + 0.5
		!X roughness
		tagdistbottom(2) =   0.5d0*texture_intensity *globaldomain(2) &  
  						   + 0.25d0*texture_intensity*globaldomain(2)*sin(2*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*sin(5*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*sin(20*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*2.d0*(rand-1)
		tagdisttop(2)    =   0.5d0*texture_intensity*globaldomain(2) &  !Minimum height
  						   + 0.25d0*texture_intensity*globaldomain(2)*sin(2*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*sin(5*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*sin(20*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*2.d0*(rand-1)
		!Z roughness
		!tagdistbottom(2) = tagdistbottom(2)  &
  		!				   + 0.1d0*globaldomain(2)*sin(2*pi*zlocation) + &
		!				   + 0.1d0*globaldomain(2)*sin(5*pi*zlocation) + &
		!				   + 0.1d0*globaldomain(2)*sin(20*pi*zlocation) + &
		!				   + 0.1d0*globaldomain(2)*2.d0*(rand-1)
		!tagdisttop(2)    =   tagdisttop(2)   &
  		!				   + 0.1d0*globaldomain(2)*sin(2*pi*zlocation) + &
		!				   + 0.1d0*globaldomain(2)*sin(5*pi*zlocation) + &
		!				   + 0.1d0*globaldomain(2)*sin(20*pi*zlocation) + &
		!				   + 0.1d0*globaldomain(2)*2.d0*(rand-1)

		!Ensure Minimum height
		if (tagdistbottom(2) .lt. 0.05d0*globaldomain(2)) then
			tagdistbottom(2) = 0.05d0*globaldomain(2)
		endif
	
		if (tagdisttop(2) .lt. 0.05d0*globaldomain(2)) then
			tagdisttop(2) = 0.05d0*globaldomain(2)
		endif

	case(converge_diverge)
		!A converging diverging channel
		tagdistbottom = 0.d0; tagdisttop=0.d0
		xlocation = rg(1)/globaldomain(1) + 0.5
		fraction_domain =  0.5d0*texture_intensity
		!Ensure Minimum height
		if (tethereddistbottom(2) .lt. 0.05d0*globaldomain(2)) then
			tethereddistbottom(2) = 0.05d0*globaldomain(2)
		endif
		if (tethereddisttop(2) .lt. 0.05d0*globaldomain(2)) then
			tethereddisttop(2) = 0.05d0*globaldomain(2)
		endif
		!N.B. One minus cosine used as periodic and has a periodic derivative (sine does not)
		tagdistbottom(2) =  0.5d0*fraction_domain*globaldomain(2) & 
								 *(1-cos(2*pi*xlocation)) & 
						  	+ tethereddistbottom(2) + 0.5d0*cellsidelength(2)
		tagdisttop(2)    =  0.5d0*fraction_domain*globaldomain(2)   & 
								 *(1-cos(2*pi*xlocation)) & 
							+ tethereddisttop(2)    + 0.5d0*cellsidelength(2)
	end select

end subroutine wall_textures

!----------------------------------------------------------------------------------
!Read molecular tags using properties defined by different place values

subroutine decode_tag(molno)
	use module_molecule_properties
	implicit none

	integer			:: l,n,m
	integer 		:: molno, no_2_decode
	integer,dimension(5) :: check, place_value

	no_2_decode = tag(molno)
	m = 1
	!Maximum 32 bit signed integer represented number is ~4x10^9
	do n = 8,0,-2
		!Move through with 100 values for each criterion
		place_value(m) = 10**n
		check(m) = no_2_decode/place_value(m)
		do l = 1,m
			check(m) = check(m) - check(l-1) * place_value(l-1) / place_value(m)
		enddo
		m = m + 1
	enddo

end subroutine decode_tag

!----------------------------------------------------------------------------------
!Read molecular tags and assign values accordingly

subroutine read_tag(molno)
	use interfaces
	use module_molecule_properties
	implicit none

	integer :: molno

	if (ensemble .eq. tag_move) then

		!Check tag and assign properties accordingly
		select case (tag(molno))
		case (free)
			!Set molecules to unfixed with no sliding velocity
			fix(:,molno) = 1
			slidev(:,molno) = 0.d0
		case (fixed)
			!Fixed Molecules
			fix(:,molno) = 0
			slidev(:,molno) = 0.d0
		case (fixed_slide)
			!Fixed with constant sliding speed
			fix(:,molno) = 0
			slidev(:,molno) = wallslidev
		case (teth)
			!Tethered molecules unfixed with no sliding velocity
			fix(:,molno) = 1
			slidev(:,molno) = 0.d0
		case (thermo)
			!Thermostatted molecules 
			fix(:,molno) = 1
			slidev(:,molno) = 0.d0
		case (teth_thermo)
			!Thermostatted Tethered molecules unfixed with no sliding velocity
			fix(:,molno) = 1
			slidev(:,molno) = 0.d0
		case (teth_slide)
			!Tethered molecules with sliding velocity
			fix(:,molno) = 1
			slidev(:,molno) = wallslidev
		case (teth_thermo_slide)
			!Thermostatted Tethered molecules unfixed with sliding velocity
			fix(:,molno) = 1
			slidev(:,molno) = wallslidev
		case (PUT_thermo)
			!Profile unbiased thermostat (Nose-Hoover)
			fix(:,molno) = 1
			slidev(:,molno) = 0.d0
		case (z_thermo)
			!Thermostat in the z direction only (Nose-Hoover) 
			fix(:,molno) = 1
			slidev(:,molno) = 0.d0
		case (cyl_teth_thermo_rotate)
			fix(:,molno) = 1
			slidev(:,molno) = 0.d0
		case default
			call error_abort("Invalid molecular Tag")
		end select

	end if

end subroutine read_tag

!------------------------------------------------------------------------------
subroutine tether_force(molno)
	use module_molecule_properties
	use arrays_MD
	implicit none

	integer                        :: molno
	double precision               :: acctmag
	double precision, dimension(3) :: at, rio

	!COEFFICIENTS MOVED TO INPUT FILE
	!Define strength of tethering potential ~ phi= k2*rio^2 
	!Force constants (k2 = 0.5*57.15) from Liem Brown and Clarke via B. D. Todd, Peter J. Daivis, and Denis J. Evans (1995) PRE. 52, 5
	!double precision, parameter    :: teth_k2=28.575
	!Define strength of tethering potential ~ phi= -k4*rio^4 - k6*rio^6
	!Force constants (k4 = 5,000, k6 = 5,000,000) from Petravich and Harrowell (2006) J. Chem. Phys.124, 014103.
	!double precision, parameter    :: teth_k4=5000.d0    
	!double precision, parameter    :: teth_k6=5000000.d0

	!Obtain displacement from initial position
	rio(:) = r(:,molno) - rtether(:,molno)

	!Apply tethering forces
	acctmag = -2.d0*teth_k2*magnitude(rio) 		 & 
			  -4.d0*teth_k4*magnitude(rio)**2.d0 & 
			  -6.d0*teth_k6*magnitude(rio)**4.d0
	at(:) = acctmag * rio(:)

	!Adjust molecular acceleration accordingly
	a(:,molno) = a(:,molno) + at(:)

	!Adjust initial postion if molecule is sliding
	rtether(:,molno) = rtether(:,molno) + slidev(:,molno)*delta_t

	!Add tethered force to stress calculation
	if (vflux_outflag .eq. 4) then
		if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
			call record_external_forces(at(:),r(:,molno))
			!call control_volume_stresses(at(:),r(:,molno),rtether(:,molno),molno)
		endif
	endif

end subroutine tether_force




! TEST ADDING SOME MOLECULES


module particle_insertion

contains

	!Inserts molecule with specified positions, velocity and optionally tag
	! Code Author: Edward Smith 2013
	subroutine insert_molecule(rin,vin,tagin)
		use arrays_MD, only : r, v, tag
		use physical_constants_MD, only : np
		use computational_constants_MD, only :halfdomain, cellsidelength, nh,ensemble,tag_move
		implicit none

		integer, intent(in),optional			:: tagin
		double precision,dimension(3),intent(in):: rin, vin

		integer									:: icell,jcell,kcell
		double precision,parameter				:: tol = 0.1d0

		!Add new molecule position and velocity to top of array
		r(:,np+1)  = rin(:)
		v(:,np+1)  = vin(:)
		if (present(tagin)) then
			tag(np+1)  = tagin
			!Read molecular tag and assign correct properties to reordered molecules
			call read_tag(np+1)
		else
			if (ensemble.eq.tag_move) tag(np+1) = 0
		endif

		!Update processor number of molecules
		np = np + 1

		!Get new molecule's cell
		icell = ceiling((r(1,np)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add nh due to halo(s)
		jcell = ceiling((r(2,np)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add nh due to halo(s)
		kcell = ceiling((r(3,np)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add nh due to halo(s)

		!Add new molecule to current cell list
		call linklist_checkpush(icell, jcell, kcell, np)

	end subroutine insert_molecule


	!Remove molecules from the array and reorder to fill the gap
	! Code Author: Edward Smith 2013
	subroutine remove_molecule(molno)
		use arrays_MD, only : r, v, tag, rtrue, rtether
		use physical_constants_MD, only : np
		use computational_constants_MD, only : ensemble,tag_move,tether_tags,rtrue_flag, &
											   halfdomain, cellsidelength, nh
		implicit none

		integer, intent(in)		:: molno

		integer					:: icell,jcell,kcell
		integer					:: icell_top,jcell_top,kcell_top

		!Get removed molecule's cell
		icell = ceiling((r(1,molno)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add nh due to halo(s)
		jcell = ceiling((r(2,molno)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add nh due to halo(s)
		kcell = ceiling((r(3,molno)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add nh due to halo(s)

		!Get top molecule's cell
		icell_top = ceiling((r(1,np)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add nh due to halo(s)
		jcell_top = ceiling((r(2,np)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add nh due to halo(s)
		kcell_top = ceiling((r(3,np)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add nh due to halo(s)

		!Pop removed molecule and top molecule from cell lists
		call linklist_pop(icell, jcell, kcell, molno)
		call linklist_pop(icell_top, jcell_top, kcell_top, np) 

		!Replace specified molecule with top one from array
		if (ensemble.eq.tag_move) then
			tag(molno)  = tag(np)
		endif
		r(:,molno)       = r(:,np)
		v(:,molno)       = v(:,np)
		if (rtrue_flag.eq.1) then
			rtrue(:,molno) = rtrue(:,np)
		endif
		if (ensemble.eq.tag_move) then
			if (any(tag(molno).eq.tether_tags)) then
				rtether(:,molno) = rtether(:,np)
			endif
		endif

		!Add top molecule to current cell list
		call linklist_checkpush(icell_top, jcell_top, kcell_top, molno) 
		
		!Update processor number of molecules
		np = np - 1

	end subroutine remove_molecule

	!Create molecular position randomly in specified cell
	subroutine create_position(icell,jcell,kcell,rout,dir,u,check)
		use librarymod, only : Maxwell_Boltzmann_vel
#if USE_COUPLER
		use CPL, only : CPL_get
#endif
		use computational_constants_MD, only : halfdomain,globaldomain,domain,delta_t,cellsidelength,nhb, &
											   iblock,jblock,kblock,npx,npy,npz
		use physical_constants_MD, only : density,np
		use calculated_properties_MD, only : temperature
		use linked_list
		use arrays_MD, only : r
		implicit none


		logical,intent(in),optional				 	:: check
		integer,intent(in) 						 	:: icell,jcell,kcell,dir
		double precision,dimension(3),intent(in) 	:: u
		double precision,dimension(3),intent(out)	:: rout

		logical							:: insert_ok
		integer							:: ixyz,i,j,k,cellnp,molnoi
		integer,dimension(3)			:: cells,block,procs
		double precision,dimension(3)	:: rand,dxyz,ri,rij
		double precision				:: rtest,dx,dy,dz
		type(node), pointer 	        :: oldi, currenti

#if USE_COUPLER
		call CPL_get(dx=dx,dy=dy,dz=dz)
#else
		dx = cellsidelength(1)
		dy = cellsidelength(2)
		dz = cellsidelength(3)
#endif

		! Define 3 arrays 
		block = (/ iblock,jblock,kblock /)
		procs = (/  npx  ,  npy ,  npz  /)
		cells = (/ icell ,jcell ,kcell  /)
		dxyz  = (/  dx  ,   dy ,   dz   /)

		
		select case(dir)
		case(1:3)
			! Surface insertion -- chooses position in 2D space a small distance from cell surface

		case(4)
			! Volume insertion -- chooses position in 3D space inside cell
			!Get 3D position
			insert_ok = .false.
			do while (insert_ok .eq. .false.) 
				!Get three random numbers
				call random_number(rand)
				do ixyz = 1,3
					rout(ixyz) = (cells(ixyz)-1)*dxyz(ixyz) + rand(ixyz)*dxyz(ixyz) & 
								  - 0.5d0*globaldomain(ixyz) &
								  - domain(ixyz)*(block(ixyz)-1)+0.5d0*domain(ixyz)*(procs(ixyz)-1)
					!print'(3i8,6f10.5)', cells,rout,(cells(ixyz)+nhb(ixyz))*dxyz(ixyz) + rand(ixyz)*dxyz(ixyz) & 
					!			  , 0.5d0*globaldomain(ixyz) &
					!			  , domain(ixyz)*(block(ixyz)-1)+0.5d0*domain(ixyz)*(procs(ixyz)-1)
				enddo

				if (present(check)) then
					insert_ok = .true.
				   !Get new molecule's cell
				    i = ceiling((rout(1)+halfdomain(1)) &
				    /cellsidelength(1))+nhb(1) !Add nh due to halo(s)
			    	j = ceiling((rout(2)+halfdomain(2)) &
				    /cellsidelength(2))+nhb(2) !Add nh due to halo(s)
				    k = ceiling((rout(3)+halfdomain(3)) &
				    /cellsidelength(3))+nhb(3) !Add nh due to halo(s)

				    !Check it against cell's current occupants
				    cellnp = cell%cellnp(i,j,k)
				    oldi => cell%head(i,j,k)%point !First molecule in list
			    	do i = 1,cellnp                    !Step through all particle in list
			        	molnoi = oldi%molno          !Get molecule number
				        ri(:) = r(:,molnoi)
				        rij = rout - ri
		
						if (dot_product(rij,rij) .lt. 1.0d0) then
							!print'(4f10.5)', rij, dot_product(rij,rij)
							insert_ok = .false.
						endif 

			    	    currenti => oldi
			        	oldi => currenti%next !obtain next item in list
				    enddo

				else
					insert_ok = .true.
				endif
			enddo
		case default 
			stop "Error in create_position - Insertion flag not specified"
		end select


		!if (present(u)) then
			!Distance using maxwell boltzmann distribution

		!	rout(2) = globaldomain(2) - delta_t * Maxwell_Boltzmann_vel(temperature,u(2))
		!else
			!Random distance from top of domain
		!	rout(2) = globaldomain(2) - 0.01d0 * rand(2) * dxyz(2)
		!endif
		!rout(3) = kcell * dz + rand(3) * dz 


		!Average position molecule would move into domain from
		!integral of Maxwell Boltzmann distribution
		!N.B. erf() is an intrinsic in fortran 2008 - for earlier version this should be defined
		!if (present(u)) then
		!	rtest = 0.5d0*delta_t*density*(1 + erf(sqrt(1/(2*temperature))*(u(2)-(globaldomain(2)-rout(2))/delta_t)))
			!if (rtest .gt. 0.99) then
			!	print*, rtest
			!endif
		!endif

	end subroutine create_position

	! Return a velocity vector for a new molecule
	subroutine create_velocity(u,vout)
		use librarymod, only : Maxwell_Boltzmann_vel3
		use calculated_properties_MD, only : temperature
		implicit none

		double precision,dimension(3),intent(out) :: vout
		double precision,dimension(3),intent(in)  :: u

		vout = Maxwell_Boltzmann_vel3(temperature,u)

	end subroutine create_velocity


	!=============================================================================
	! Particle insertion and removal algorithms
	!=============================================================================
	! USHER algorithm: Delgado-Buscalioni and Coveney JCP 119 2 July 2003
	! Code Author: David Trevelyan 2013
	subroutine usher_get_insertion_pos(startpos, Utarget, finishpos, succeed, &
	                                   Ufinish, maxiter_op, dSmax_op, tol_op)
		use physical_constants_MD, only: density
		implicit none

		real(kind(0.d0)) , intent(in) :: startpos(3), Utarget
		real(kind(0.d0)), intent(out) :: finishpos(3), Ufinish
		logical, intent(out)          :: succeed

		integer, optional, intent(in) :: maxiter_op	
		real(kind(0.d0)),optional, intent(in) :: dSmax_op 
		real(kind(0.d0)),optional, intent(in) :: tol_op

		integer :: n, maxiter                     
		real(kind(0.d0)) :: dS, dSmax            ! Step lengths
		real(kind(0.d0)) :: dir                  ! Desired direction up/down in U
		real(kind(0.d0)) :: tol                  ! U toler for stopping algorithm
		real(kind(0.d0)) :: steppos(3) 
		real(kind(0.d0)) :: U, Uold, Uovlp       
		real(kind(0.d0)) :: f(3), fhat(3), fmag

		! Init	
		succeed = .false.	

		! Set default parameters: see pg. 986 of JCP reference for figures of
		! optimal values 
		maxiter = 200
		dSmax   = 0.1*(density**-1.5)
		dS      = dSmax
		Uovlp   = 10.0**4.0
		tol     = 0.01 
		if (present(maxiter_op)) maxiter = maxiter_op
		if (present(dSmax_op)) dSmax = dSmax_op 
		if (present(tol_op)) tol = tol_op 

		! Start steepest descent
		steppos = startpos
		Uold = Uovlp
		do n = 1, maxiter	

			call compute_force_and_potential_at(steppos, U, f)

			!fmag = norm2(f) 
			fmag = sqrt(dot_product(f,f)) 
			dir  = sign(1.d0,(U-Utarget))

			! Evaluate conditions for stopping the algorithm
			if ( abs((U-Utarget)/Utarget) .lt. tol) then
				succeed = .true.
				finishpos = steppos
				Ufinish = U
				return	

			else if ( fmag .eq. 0.0 ) then
				! Hit a region of 0 force
				succeed = .false.
				finishpos = 0.d0 !sqrt(-1.d0)
				Ufinish = U
				return	
		
			else if ( dir*(U-Uold) > 0.0 ) then
				! Else if step has gone the wrong way
				succeed = .false.
				finishpos = 0.d0 !sqrt(-1.d0)
				Ufinish = U
				return	
			
			end if 

			! Determine step size
			if ( U .gt. Uovlp ) then
				dS = 0.9 - (4.0/U)**(1.0/12.0)
			else
				dS = min(dSmax,abs(U-Utarget)/fmag)	
			end if

			! Perform steepest descent step
			fhat = f / fmag 
			steppos(:) = steppos(:) + dir*fhat(:)*dS
			Uold = U

		end do

	end subroutine usher_get_insertion_pos

end module particle_insertion

!Insert Molecules and call USHER algorithm
! Code Author: Edward Smith 2013
subroutine usher_insert(nparticles)
	use particle_insertion
	use computational_constants_MD, only: iter, halfdomain,ncells
	use physical_constants_MD, only: Potential_sLRC, np
	implicit none

	integer, intent(in) :: nparticles

	integer :: n,i,maxattempts=10000,icell,jcell,kcell
	logical :: pos_found 
	real(kind(0.d0)) :: Utarget, Ufinish
	real(kind(0.d0)) :: startpos(3), insertpos(3)
	real(kind(0.d0)),dimension(3) :: avevel,vnew,rand

	avevel = 0.d0

	Utarget  = get_Utarget() 

	do n = 1,nparticles	

		call random_number(rand)
		icell = ceiling(rand(1) * ncells(1))
		jcell = ceiling(rand(2) * ncells(2))
		kcell = ceiling(rand(3) * ncells(3))

		call create_velocity(avevel,vnew)

		do i = 1, maxattempts

			call create_position(icell,jcell,kcell,startpos,4,avevel)	

			!startpos = get_randompos()	
			call usher_get_insertion_pos(startpos,Utarget,insertpos,pos_found,Ufinish)

			if (pos_found) then

				!print('(a,2f9.3,a,3f9.3,a,i4,a,i6)'), ' USHER SUCCEEDED: Utarget, Ufinish: ', &
				!                                    Utarget, Ufinish, ', r:',  &
				!                                    insertpos, ', after ', i,  &
				!                                    ' attempts , new np', np
				!write(40000+iter,*) insertpos
				call insert_molecule(insertpos,vnew)	
				call messenger_updateborders(1) !Update borders ready for next insertion - halo rebuilt

				!Write local xy field @ z
				!call simulation_write_potential_field(max(-halfdomain(1),insertpos(1)-3.0), &
				!                                      min( halfdomain(1),insertpos(1)+3.0), &
				!                                      max(-halfdomain(2),insertpos(2)-3.0), &
				 !                                     min( halfdomain(2),insertpos(2)+3.0), &
				!                                      insertpos(3), 300, 30000+iter)
				!Write entire xy field @ z
				!call simulation_write_potential_field(-halfdomain(1),halfdomain(1), &
				!                                      -halfdomain(2),halfdomain(2), &
				!                                       insertpos(3), 400, 30000+iter)

				exit

			end if

		end do

	end do

	if (.not. pos_found) then
		print*, ' usher failed  Utarget =', Utarget
	end if

contains

	function get_Utarget() result(U)
		use calculated_properties_MD, only: potenergy
		implicit none
	
		real(kind(0.d0)) :: U 
		U = potenergy + potential_sLRC	
	
	end function get_Utarget
	
	function get_randompos() result(pos)
		use computational_constants_MD, only: domain
		implicit none
	
		integer :: ixyz
		real(kind(0.d0)) :: pos(3), rand
		
		do ixyz = 1,3
			call random_number(rand)
			rand = (rand-0.5)*domain(ixyz)
			pos(ixyz) = rand
		end do
	
	end function get_randompos

end subroutine

subroutine remove_mols(nmols)
	use particle_insertion, only : remove_molecule
	use physical_constants_MD, only: np
	implicit none

	integer,intent(in)	:: nmols

	integer			 :: n,molno
	double precision :: rand

	do n=1,nmols
		call random_number(rand)
		molno = ceiling(rand * np)
		call remove_molecule(molno)
	enddo

end subroutine remove_mols

