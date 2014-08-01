!----------------------------------------------------------------------------------
!								Molecular Properties
! Routines for setting and adjusting the properties of the molecules in the domain
!
!---------------------------------------------------------------------------------

module module_molecule_properties

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD

contains

    function get_tag_status(rg,status_type) result(tag_status)
        use interfaces, only: error_abort
        use computational_constants_MD, only : texture_type, domain
        use calculated_properties_MD, only : nbins
        implicit none

        double precision :: rg(3)   ! Global position
        character(*) :: status_type ! Specifies thermo, tethered, fixed, etc

        integer :: ixyz
        double precision :: bottom(3) ! Global position of bottom of domain
        double precision :: top(3)
        double precision :: tagdistbottom(3) ! Distance specified by user
        double precision :: tagdisttop(3)
        double precision :: Mbinsize(3)
        logical :: tag_status 

        bottom = (/ -globaldomain(1)/2.d0, -globaldomain(2)/2.d0, -globaldomain(3)/2.d0 /)
        top	= (/  globaldomain(1)/2.d0,  globaldomain(2)/2.d0,  globaldomain(3)/2.d0 /)

        select case (status_type)
        case ('thermo')
            tagdistbottom(:) = thermstatbottom(:)
            tagdisttop(:)	 = thermstattop(:)
            !Thermostat complicated wall texture if specified, otherwise thermostat
            !is based on thermstatbottom/thermstattop
            if (texture_type .ne. 0 .and. texture_therm .eq. 1) then
                call wall_textures(texture_type,rg,tagdistbottom,tagdisttop)
            endif
        case ('teth')
            tagdistbottom(:) = tethereddistbottom(:)
            tagdisttop(:)	= tethereddisttop(:)
            !Apply a complicated wall texture if specified
            if (texture_type .ne. 0) call wall_textures(texture_type,rg,tagdistbottom,tagdisttop)
        case ('fixed')
            tagdistbottom(:) = fixdistbottom(:)
            tagdisttop(:)	= fixdisttop(:)
        case ('slide')
            tagdistbottom(:) = slidedistbottom(:)
            tagdisttop(:)	= slidedisttop(:)
        case default
            call error_abort("Unrecognised tag status type")
        end select

        tag_status = .false.
        ! If within a tagged region then mark tag status as true and return
        Mbinsize = domain/nbins
        do ixyz = 1,3
            if (rg(ixyz) .le. bottom(ixyz) + tagdistbottom(ixyz) .or. &
               (rg(ixyz) .ge. top(ixyz)	   - tagdisttop(ixyz)) ) then
                tag_status = .true.
                return

            end if  
        end do

    end function get_tag_status

end module module_molecule_properties

!----------------------------------------------------------------------------------
!Create an array of molecular tags

	!For initialunitsize "a"
	!				[  o	 o ]
	!a (1 cell size) |	 o	]  a/2 (distance between molcules)  
	!				|  o	 o
	!				 __________]  a/4 (distance from bottom of domain)
	!
	!So use (0.20+0.5d0*mol_layers)*initialunitsize(ixyz)
subroutine setup_cylinder_tags
	use module_molecule_properties
	use concentric_cylinders
	use interfaces, only: error_abort
	use messenger, only: globalise
	use librarymod, only : cpolariser 
	implicit none

	integer :: n
	double precision :: rglob(3), rpol(3)

	do n = 1,np
		
		rglob = globalise(r(:,n))
		rpol  = cpolariser(rglob)
		
		if (rpol(1) .lt. r_oi) then
			tag(n) = cyl_teth_thermo_rotate
		else if (rpol(1) .gt. r_io) then
			tag(n) = teth_thermo
		else
			tag(n) = free 
		end if

	end do  

end subroutine setup_cylinder_tags

subroutine setup_cylinder_tags_equilibrate
	use module_molecule_properties
	use concentric_cylinders
	use interfaces, only: error_abort
	use messenger, only: globalise
	use librarymod, only : cpolariser 
	implicit none

	integer :: n
	double precision :: rglob(3), rpol(3)

	do n = 1,np
		
		rglob = globalise(r(:,n))
		rpol  = cpolariser(rglob)
		
		if (rpol(1) .lt. r_oi) then
			tag(n) = cyl_teth_thermo_rotate
		else if (rpol(1) .gt. r_io) then
			tag(n) = teth_thermo
		else
			tag(n) = thermo
		end if

	end do  

	tag_thermostat_active = .true.

end subroutine setup_cylinder_tags_equilibrate


subroutine setup_location_tags
	use module_molecule_properties
	use interfaces, only: error_abort
	use messenger, only: globalise
	use computational_constants_MD, only : texture_type, domain
	use calculated_properties_MD, only : nbins
	implicit none

	integer :: n
	double precision :: rglob(3)
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
			if ( l_teth .and. .not. l_thermo .and.	   l_slide  ) tag(n) = teth_slide
			if ( l_teth .and.	   l_thermo .and. .not. l_slide  ) tag(n) = teth_thermo
			if ( l_teth .and.	   l_thermo .and.	   l_slide  ) tag(n) = teth_thermo_slide

			! Fixed 
			if ( l_fixed .and. .not. l_slide ) tag(n) = fixed
			if ( l_fixed .and.	   l_slide ) tag(n) = fixed_slide

			! Thermo only
			if ( l_thermo .and. .not. l_teth .and. .not. l_slide ) tag(n) = thermo

		enddo

		!Check if any molecules are thermostatted and switch on flag
		call get_tag_thermostat_activity(tag_thermostat_active)

	end if

end subroutine setup_location_tags




subroutine reset_location_tags
	use module_molecule_properties
	use interfaces, only: error_abort
	use messenger, only: globalise
	use computational_constants_MD, only : texture_type, domain
	use calculated_properties_MD, only : nbins
	implicit none

	integer :: n
	double precision :: rglob(3), Mbinsize(3)
	logical :: l_thermo

	if (ensemble .eq. tag_move) then

		!Setup fixed wall and location dependent thermostat tags
		do n = 1,np

			rglob(:) = globalise(r(:,n))

            if (any(tag(n) .eq. tether_tags)) then
                cycle
            else if (tag(n) .eq. fixed .or. tag(n) .eq. fixed_slide) then
                cycle 
            end if

			l_thermo = get_tag_status(rglob,'thermo') 
			if (l_thermo) tag(n) = thermo

		enddo

		!Check if any molecules are thermostatted and switch on flag
		call get_tag_thermostat_activity(tag_thermostat_active)

	end if

end subroutine reset_location_tags

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

	integer,intent(in)	  :: texture_type
	double precision,dimension(3),intent(in) :: rg
	double precision,dimension(3),intent(out):: tagdistbottom,tagdisttop

	double precision		:: xlocation,ylocation,zlocation,rand,fraction_domain,postheight

	select case (texture_type)
	case(0)
		!Case of simply tethered walls
		tagdistbottom = tethereddistbottom
		tagdisttop	  = tethereddisttop
	case(posts)

        tagdisttop = tethereddisttop

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
		tagdisttop(2)	=   0.5d0*texture_intensity*globaldomain(2) &  !Minimum height
						   + 0.25d0*texture_intensity*globaldomain(2)*sin(2*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*sin(5*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*sin(20*pi*xlocation) + &
						   + 0.25d0*texture_intensity*globaldomain(2)*2.d0*(rand-1)
		!Z roughness
		!tagdistbottom(2) = tagdistbottom(2)  &
		!				  + 0.1d0*globaldomain(2)*sin(2*pi*zlocation) + &
		!				  + 0.1d0*globaldomain(2)*sin(5*pi*zlocation) + &
		!				  + 0.1d0*globaldomain(2)*sin(20*pi*zlocation) + &
		!				  + 0.1d0*globaldomain(2)*2.d0*(rand-1)
		!tagdisttop(2)	=   tagdisttop(2)   &
		!				  + 0.1d0*globaldomain(2)*sin(2*pi*zlocation) + &
		!				  + 0.1d0*globaldomain(2)*sin(5*pi*zlocation) + &
		!				  + 0.1d0*globaldomain(2)*sin(20*pi*zlocation) + &
		!				  + 0.1d0*globaldomain(2)*2.d0*(rand-1)

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
		tagdisttop(2)	=  0.5d0*fraction_domain*globaldomain(2)   & 
								 *(1-cos(2*pi*xlocation)) & 
							+ tethereddisttop(2)	+ 0.5d0*cellsidelength(2)
	end select

end subroutine wall_textures

!----------------------------------------------------------------------------------
!Read molecular tags using properties defined by different place values

subroutine decode_tag(molno)
	use module_molecule_properties
	implicit none

	integer		 :: l,n,m
	integer		 :: molno, no_2_decode
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
	use messenger, only : globalise
	implicit none

	integer :: molno
	double precision   :: r_global(3)

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
			!r_global = globalise(r(:,molno))
			slidev(:,molno) = wallslidev!*sign(1.d0,r_global(2))
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
			r_global = globalise(r(:,molno))
			slidev(:,molno) = wallslidev*sign(1.d0,r_global(2))
		case (teth_thermo_slide)
			!Thermostatted Tethered molecules unfixed with sliding velocity
			fix(:,molno) = 1
			r_global = globalise(r(:,molno))
			slidev(:,molno) = wallslidev*sign(1.d0,r_global(2))
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
	use module_record_external_forces, only : record_external_forces
	use librarymod, only : magnitude 
	implicit none

	integer						:: molno
	double precision			   :: acctmag
	double precision, dimension(3) :: at,mat, rio

	!COEFFICIENTS MOVED TO INPUT FILE
	!Define strength of tethering potential ~ phi= k2*rio^2 
	!Force constants (k2 = 0.5*57.15) from Liem Brown and Clarke via B. D. Todd, Peter J. Daivis, and Denis J. Evans (1995) PRE. 52, 5
	!double precision, parameter	:: teth_k2=28.575
	!Define strength of tethering potential ~ phi= -k4*rio^4 - k6*rio^6
	!Force constants (k4 = 5,000, k6 = 5,000,000) from Petravich and Harrowell (2006) J. Chem. Phys.124, 014103.
	!double precision, parameter	:: teth_k4=5000.d0	
	!double precision, parameter	:: teth_k6=5000000.d0

	!Obtain displacement from initial position
	rio(:) = r(:,molno) - rtether(:,molno)

	!Apply tethering forces
	acctmag = -2.d0*teth_k2*magnitude(rio)	   & 
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
			mat = -at
			call record_external_forces(mat(:),r(:,molno))
			! There was a time when I though that the CV conservation would require interactions
			! which cross the surface. However, it appears this is not the case...
			!call control_volume_stresses(at(:),r(:,molno),rtether(:,molno),molno)
		endif
	endif

end subroutine tether_force




! TEST ADDING SOME MOLECULES


module particle_insertion

contains

	subroutine insert_molecules(rin,vin,tagin)
		implicit none

		integer,dimension(:), intent(in),optional	:: tagin
		double precision,dimension(:,:),intent(in)	:: rin, vin

		integer	:: n		

		if (present(tagin)) then
			do n =1,size(tagin)
				call insert_molecule(rin(:,n),vin(:,n),tagin(n))
			enddo
		else
			do n =1,size(rin,2)
				call insert_molecule(rin(:,n),vin(:,n))
			enddo
		endif

	end subroutine insert_molecules

	!Inserts molecule with specified positions, velocity and optionally tag
	! Code Author: Edward Smith & DT 2014
	subroutine insert_molecule(rin,vin,tagin)
		use arrays_MD, only : r, v, tag
		use physical_constants_MD, only : np, halo_np, globalnp
		use computational_constants_MD, only :halfdomain,cellsidelength, nh,&
											  ensemble,tag_move, ncells
		implicit none

		integer, intent(in),optional			:: tagin
		double precision,dimension(3),intent(in):: rin, vin

		integer :: icell, jcell, kcell
		integer :: icell_halo, jcell_halo, kcell_halo
		double precision,parameter :: tol = 0.1d0

		!Copy np + 1 halo molecule to (np + halo_np) + 1 position
		r(:,np+halo_np+1)  = r(:,np+1)
		v(:,np+halo_np+1)  = v(:,np+1)

		!Linklist update -- pop np+1 from its cell
		call linklist_findandpop(np+1, icell_halo, jcell_halo, kcell_halo)

		!print*, 'halo after pop and before push:'
		!call linklist_print(icell_halo, jcell_halo, kcell_halo)

		call linklist_checkpush(icell_halo, jcell_halo, kcell_halo, np+1+halo_np)

		!print*, 'halo after push:'
		!call linklist_print(icell_halo, jcell_halo, kcell_halo)
   
		!Add new molecule position and velocity to top of array
		r(:,np+1)  = rin(:)
		v(:,np+1)  = vin(:)
		if (present(tagin)) then
			tag(np+1)  = tagin
			!Read molecular tag and assign correct properties to reordered molecules
			call read_tag(np+1)
			!todo rtether etc
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

		!print*, 'cell before push:'
		!call linklist_print(icell, jcell, kcell)

		!Add new molecule to current cell list
		call linklist_checkpush(icell, jcell, kcell, np)

		!print*, 'cell after push:'
		!call linklist_print(icell, jcell, kcell)

	end subroutine insert_molecule


	!Remove molecules from the array and reorder to fill the gap
	! Code Author: Edward Smith & DT 2014
	subroutine remove_molecule(molno)
		use arrays_MD, only : r, v, tag, rtrue, rtether
		use physical_constants_MD, only : np, globalnp, halo_np
		use computational_constants_MD, only : ensemble,tag_move,tether_tags,rtrue_flag, &
											   halfdomain, cellsidelength, nh, ncells
		implicit none

		integer, intent(in) :: molno

		integer :: icell,jcell,kcell
		integer :: icell_top,jcell_top,kcell_top
		integer :: icell_halo,jcell_halo,kcell_halo

		!Get removed molecule's cell
		icell = ceiling((r(1,molno)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add nh due to halo(s)
		jcell = ceiling((r(2,molno)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add nh due to halo(s)
		kcell = ceiling((r(3,molno)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add nh due to halo(s)


		!print*, 'molno cell linklist before pop:'
		!call linklist_print(icell,jcell,kcell)
		!print*, 'np cell linklist before pop:'
		!call linklist_print(icell_top,jcell_top,kcell_top)

		!Pop removed molecule and top molecule from cell lists
		call linklist_findandpop(molno,icell,jcell,kcell)
		call linklist_findandpop(np,icell,jcell,kcell) !TODO new cell return,
													   !to be developed

		!print*, 'molno cell linklist after pop:'
		!call linklist_print(icell,jcell,kcell)
		!print*, 'np cell linklist after pop:'
		!call linklist_print(icell_top,jcell_top,kcell_top)

		!Replace specified molecule with top one from array
		if (ensemble.eq.tag_move) then
			tag(molno)  = tag(np)
		endif
		r(:,molno)	   = r(:,np)
		v(:,molno)	   = v(:,np)
		if (rtrue_flag.eq.1) then
			rtrue(:,molno) = rtrue(:,np)
		endif
		if (ensemble.eq.tag_move) then
			if (any(tag(np).eq.tether_tags)) then
				rtether(:,molno) = rtether(:,np)
			endif
		endif

		!Add top molecule to current cell list
		!Get top molecule's cell
		icell_top = ceiling((r(1,np)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add nh due to halo(s)
		jcell_top = ceiling((r(2,np)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add nh due to halo(s)
		kcell_top = ceiling((r(3,np)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add nh due to halo(s)
		call linklist_checkpush(icell_top, jcell_top, kcell_top, molno) 

		!print*, 'molno cell linklist after push again:'
		!call linklist_print(icell_top,jcell_top,kcell_top)

		!Replace void at np with top halo molecule TODO
		!Get removed molecule's cell
		icell_halo = ceiling((r(1,np+halo_np)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add nh due to halo(s)
		jcell_halo = ceiling((r(2,np+halo_np)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add nh due to halo(s)
		kcell_halo = ceiling((r(3,np+halo_np)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add nh due to halo(s)

		!Replace specified molecule with top one from array
		if (ensemble.eq.tag_move) then
			tag(np)  = tag(np+halo_np)
		endif
		r(:,np)	   = r(:,np+halo_np)
		v(:,np)	   = v(:,np+halo_np)
		if (rtrue_flag.eq.1) then
			rtrue(:,np) = rtrue(:,np+halo_np)
		endif
		if (ensemble.eq.tag_move) then
			if (any(tag(np+halo_np).eq.tether_tags)) then
				rtether(:,np) = rtether(:,np+halo_np)
			endif
		endif

		!print*, 'np+halo_np cell linklist before pop:'
		!call linklist_print(icell_halo,jcell_halo,kcell_halo)

		call linklist_findandpop(np+halo_np)

		!print*, 'np+halo_np cell linklist after pop and before push:'
		!call linklist_print(icell_halo,jcell_halo,kcell_halo)

		! If halo cell molecule is out of the domain 
		if (icell_halo .gt. ncells(1)+2) stop "Error -- mol outside halo!"  ! icell_halo = ncells(1) + 2 
		if (jcell_halo .gt. ncells(2)+2) stop "Error -- mol outside halo!"  ! jcell_halo = ncells(2) + 2 
		if (kcell_halo .gt. ncells(3)+2) stop "Error -- mol outside halo!"  ! kcell_halo = ncells(3) + 2 
		if (icell_halo .lt. 1) stop "Error -- mol outside halo!"  ! icell_halo = 1
		if (jcell_halo .lt. 1) stop "Error -- mol outside halo!"  ! jcell_halo = 1
		if (kcell_halo .lt. 1) stop "Error -- mol outside halo!"  ! kcell_halo = 1
		call linklist_checkpush(icell_halo, jcell_halo, kcell_halo, np) 

		!print*, 'np+halo_np cell linklist after push:'
		!call linklist_print(icell_halo,jcell_halo,kcell_halo)

		! Update borders to be safe
		!call messenger_updateborders(1)	   

		!Update processor number of molecules
		np = np - 1
		!globalnp = globalnp - 1

	end subroutine remove_molecule

	! Remove a specified molecule and stick a new one at the top
	! of the stack
	subroutine replace_molecule(molno,rin,vin,tagin)
		implicit none

		integer, intent(in) :: molno
		integer, intent(in),optional			:: tagin
		double precision,dimension(3),intent(in):: rin, vin

		call remove_molecule(molno)
		call insert_molecule(rin,vin,tagin)

	end subroutine replace_molecule

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


		logical,intent(in),optional				 :: check
		integer,intent(in)						  :: icell,jcell,kcell,dir
		double precision,dimension(3),intent(in)	:: u
		double precision,dimension(3),intent(out)   :: rout

		logical						 :: insert_ok
		integer						 :: ixyz,i,j,k,cellnp,molnoi
		integer,dimension(3)			:: cells,block,procs
		double precision,dimension(3)   :: rand,dxyz,ri,rij
		double precision				:: dx,dy,dz
		type(node), pointer			 :: oldi, currenti

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
					!			 , 0.5d0*globaldomain(ixyz) &
					!			 , domain(ixyz)*(block(ixyz)-1)+0.5d0*domain(ixyz)*(procs(ixyz)-1)
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
					do i = 1,cellnp					!Step through all particle in list
						molnoi = oldi%molno		  !Get molecule number
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

		!   rout(2) = globaldomain(2) - delta_t * Maxwell_Boltzmann_vel(temperature,u(2))
		!else
			!Random distance from top of domain
		!   rout(2) = globaldomain(2) - 0.01d0 * rand(2) * dxyz(2)
		!endif
		!rout(3) = kcell * dz + rand(3) * dz 


		!Average position molecule would move into domain from
		!integral of Maxwell Boltzmann distribution
		!N.B. erf() is an intrinsic in fortran 2008 - for earlier version this should be defined
		!if (present(u)) then
		!   rtest = 0.5d0*delta_t*density*(1 + erf(sqrt(1/(2*temperature))*(u(2)-(globaldomain(2)-rout(2))/delta_t)))
			!if (rtest .gt. 0.99) then
			!   print*, rtest
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
									   Ufinish, maxiter_op, dSmax_op, tol_op, &
									   extra_pos)
		use physical_constants_MD, only: density
		use module_compute_forces, only : compute_force_and_potential_at
		implicit none

		double precision , intent(in) :: startpos(3), Utarget
		double precision, intent(out) :: finishpos(3), Ufinish
		logical, intent(out)		  :: succeed
		double precision,dimension(:,:),allocatable,optional,intent(in)  :: extra_pos

		integer, optional, intent(in) :: maxiter_op 
		double precision,optional, intent(in) :: dSmax_op 
		double precision,optional, intent(in) :: tol_op

		integer :: n, maxiter					 
		double precision :: dS, dSmax			! Step lengths
		double precision :: dir				  ! Desired direction up/down in U
		double precision :: ratio,tol			! U toler for stopping algorithm
		double precision :: steppos(3) 
		double precision :: U, Uold, Uovlp	   
		double precision :: f(3), fhat(3), fmag

		!print*, 'usher_get_insertion_pos', present(extra_pos)

		! Init  
		succeed = .false.   

		! Set default parameters: see pg. 986 of JCP reference for figures of
		! optimal values 
		maxiter = 200
		dSmax   = 0.1*(density**-1.5)
		dS	  = dSmax
		Uovlp   = 10.0**4.0
		tol	 = 0.001 
		if (present(maxiter_op)) maxiter = maxiter_op
		if (present(dSmax_op)) dSmax = dSmax_op 
		if (present(tol_op)) tol = tol_op 

		! Start steepest descent
		steppos = startpos
		Uold = Uovlp
		do n = 1, maxiter   

			if (present(extra_pos)) then
				call compute_force_and_potential_at(input_pos=steppos, Usum=U, f=f,extra_pos=extra_pos)
			else
				call compute_force_and_potential_at(input_pos=steppos, Usum=U, f=f)
			endif
			!print('(a,4x,i4,8x,3f10.4,8x,f10.4,8x,3f10.4)'), 'n,steppos,U,f', n, steppos, U, f

			!fmag = norm2(f) 
			fmag = sqrt(dot_product(f,f)) 
			dir  = sign(1.d0,(U-Utarget))

			! Evaluate conditions for stopping the algorithm
			if (abs(Utarget) .ne. 0.d0) then
				ratio = abs((U-Utarget)/Utarget)
			else
				ratio = abs(U-Utarget)
			end if

			if ( ratio .lt. tol) then
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

subroutine reinsert_molecules
	use particle_insertion
	use computational_constants_MD, only: halfdomain, nproc, VOID, & 
										  irank,iroot,iter,iblock,jblock,kblock,npx,npy,npz
	use physical_constants_MD, only: np, tethereddistbottom, insertnp, rcutoff2
    use boundary_MD, only: open_boundary
	use arrays_MD, only: r, v
	use  messenger_data_exchange, only : GlobalGatherv, planeGatherv
	use messenger, only : globalise,  plane_nproc, planerank
	implicit none

	integer :: n, i,j,ixyz
	integer :: reinsertnp,reinsert_attemp,procmolcount,insert_flag, block(3),npxyz(3)
	integer :: subcomminsertnp, insert_temp(1),proc, nproc_insert, commrank
	integer,dimension(:),allocatable ::rdisps, subcommmol_to_proc
	
	integer,dimension(:,:),allocatable :: insertnp_array,reinsert_molno
	double precision :: x,y,z,dx,dy,dz,vx,vy,vz
	double precision :: Ufinish, fdummy(3)
	double precision :: startpos(3), insertpos(3),rij(3),rij2
	double precision,dimension(:),allocatable :: sendbuf, recvbuf
	double precision,dimension(:,:),allocatable :: insert_locs,insert_vels,reinsert_locs,subcomminsert_locs

	insert_flag = 4

	! If no open boundaries then return
	if (all(open_boundary.eq.0)) then
		return
	end if

	!Insert np is obtained from all molecules which leave during timestep
	insert_temp(1) = insertnp

	select case (insert_flag)
	case(0)

		!Get global number of molecules to insert from leaving molecules
		nproc_insert = nproc; commrank = irank
		allocate(insertnp_array(1,nproc_insert))
		call globalGather(insert_temp,insertnp_array,1)
		subcomminsertnp = sum(insertnp_array)

		call calculate_mol_proc_redistribution_random(subcomminsertnp,nproc_insert,commrank, & 
														insertnp, insertnp_array,subcommmol_to_proc)

	case(1,3,5) !Bottom x,y,z
		ixyz = nint(0.5*(insert_flag+1))
		nproc_insert = plane_nproc(ixyz); commrank = planerank(ixyz)
		block = (/ iblock, jblock, kblock /)
		if (block(ixyz) .eq. 1 .and. open_boundary(insert_flag) .eq. 1) then

			!Share number of molecules to insert across subcomm
			allocate(insertnp_array(1,nproc_insert))
			call planeGatherInt(insert_temp,insertnp_array,1,ixyz)
			subcomminsertnp = sum(insertnp_array)

			call calculate_mol_proc_redistribution_random(subcomminsertnp,nproc_insert,commrank, & 
															insertnp, insertnp_array,subcommmol_to_proc)

		elseif (open_boundary(insert_flag) .eq. 0 ) then
			stop "Error - open boundary off but insertion requested"
		else
			!Outside of reinsertion region
			insertnp = 0; subcomminsertnp = 0
			allocate(insertnp_array(1,nproc_insert)); insertnp_array = 0
		endif

	case(2,4,6) !Top x,y,z
		ixyz = nint(0.5*insert_flag)
		nproc_insert = plane_nproc(ixyz); commrank = planerank(ixyz)
		block = (/ iblock, jblock, kblock /)
		npxyz = (/   npx,   npy,    npz   /)
		if (block(ixyz) .eq. npxyz(ixyz) .and. open_boundary(insert_flag) .eq. 1) then
			!Share number of molecules to insert across subcomm
			allocate(insertnp_array(1,nproc_insert))
			call planeGatherInt(insert_temp,insertnp_array,1,ixyz)
			subcomminsertnp = sum(insertnp_array)

			call calculate_mol_proc_redistribution_random(subcomminsertnp,nproc_insert, commrank, & 
															insertnp, insertnp_array,subcommmol_to_proc)

			!print*, 'values',irank,iblock,jblock,kblock,insertnp_array,insertnp,subcomminsertnp

		elseif (open_boundary(insert_flag) .eq. 0 ) then
			stop "Error - open boundary off but insertion requested"
		else
			!Outside of reinsertion region
			insertnp = 0; subcomminsertnp = 0
			allocate(insertnp_array(1,nproc_insert)); insertnp_array = 0
		endif

	case default
		stop "Error in insert flag type"
	end select

	if (subcomminsertnp .eq. 0) return

    !print'(a,7i8)', 'Insert in parallel ',irank,iblock,jblock,kblock,nproc, insertnp,subcomminsertnp 
	!print*, 'Insert array ',irank,iblock,jblock,kblock,insertnp_array

	!Insert molecules using the USHER algorithm
	call usher_get_positions(flag=insert_flag,insertnp=insertnp,insert_locs=insert_locs)

	!write(100+irank,*), 'insert locs',irank,iblock,jblock,kblock,iter, insert_locs

	!Gatherv takes array of different size from each processes and exchanges
	allocate(subcomminsert_locs(3,subcomminsertnp))
	allocate(rdisps(nproc_insert))
	rdisps = 0
	do i=1,nproc_insert-1
		rdisps(i+1) = rdisps(i) + 3*insertnp_array(1,i)
		!print'(9i5)', irank,iblock,jblock,kblock,i,rdisps(i+1),rdisps(i),3*insertnp_array(1,i),nproc_insert
	enddo

	!Exchange all inserted positions and syncronise
	allocate(recvbuf(3*subcomminsertnp))		
	allocate(sendbuf(3*insertnp))
	!print'(a,13i8)', 'BUFFER SHAPES', irank,insertnp,insertnp_array(1,:),subcomminsertnp, shape(sendbuf),shape(recvbuf),shape(insert_locs),shape(subcomminsert_locs)
	sendbuf = reshape(globalise(insert_locs),(/3*insertnp/))
	recvbuf = reshape(subcomminsert_locs,(/3*subcomminsertnp/))

	select case (insert_flag)
	case(0)
		call globalGatherv(sendbuf,3*insertnp, & 
						   recvbuf,3*insertnp_array(1,:),rdisps)
	case(1:6) !Bottom/Top x,y,z
		call planeGatherv( sendbuf,3*insertnp, & 
						   recvbuf,3*insertnp_array(1,:),rdisps,ixyz)
	case default
		stop "Error in insert flag type"
	end select

	subcomminsert_locs = reshape(recvbuf,(/3,nint(size(recvbuf)/3.d0)/))

!	!write(1000+irank,*), 'subcomminsert locs',irank,iter, subcomminsert_locs
!	if (irank .eq. 1) then
!		proc = 1
!		do i=1,subcomminsertnp
!			proc = subcommmol_to_proc(i)
!			!print*, 'aaa', i, proc, nproc_insert, rdisps(proc), insertnp_array(1,proc), sum(insertnp_array(1,1:proc))
!		enddo
!	endif


	call check_insert_overlap(gnip=subcomminsertnp,inp=insertnp,gi_locs=subcomminsert_locs, & 
							  gm_to_p = subcommmol_to_proc, rank_=commrank, ri_molno = reinsert_molno,rinp=reinsertnp)

	!Exit if nothing left to reinsert
	if (reinsertnp .eq. 0) then
		!do n=1,size(insert_locs,2)
		!	write(5000+irank,'(4i6,6f14.5)'), iter, iblock,jblock,kblock,globalise(insert_locs(:,n)), insert_vels(:,n)
		!enddo
		call get_random_velocity(insertnp,insert_vels)
		call insert_molecules(insert_locs,insert_vels)
	else
		do reinsert_attemp = 1,100
			print'(2(a,i8),i6,a,2i8)', 'At iter = ',  iter, ' attempt = ',  reinsert_attemp,reinsertnp, ' molecule(s) overlap -- attempt to reinsert ',  irank, insertnp
			call usher_get_positions(flag=insert_flag,insertnp=reinsertnp,insert_locs=reinsert_locs)
			!print'(a,4i6,3f10.5,2i7)', 'details', iter, irank, reinsert_attemp, reinsertnp, reinsert_locs, reinsert_molno
			do n =1,reinsertnp
				insert_locs(:,reinsert_molno(1,n)) = reinsert_locs(:,n)
				subcomminsert_locs(:,reinsert_molno(2,n)) = globalise(reinsert_locs(:,n))
			enddo
			call check_insert_overlap(gnip=subcomminsertnp,inp=insertnp,gi_locs=subcomminsert_locs, & 
							  		  gm_to_p = subcommmol_to_proc, rank_=commrank, ri_molno = reinsert_molno,rinp=reinsertnp)
			if (reinsertnp .eq. 0) exit
		enddo
		do n=1,size(insert_locs,2)
			write(5000+irank,'(4i6,6f14.5)'), iter, iblock,jblock,kblock,globalise(insert_locs(:,n)), insert_vels(:,n)
		enddo
		call get_random_velocity(insertnp,insert_vels)
		call insert_molecules(insert_locs,insert_vels)
	endif

	deallocate(insertnp_array)
	deallocate(insert_locs)
	deallocate(subcomminsert_locs)
	deallocate(rdisps)
	deallocate(recvbuf)
	deallocate(sendbuf)

contains 


	subroutine usher_get_positions(flag,insertnp,insert_locs,Utarget_in,extra_pos)
		use interfaces, only : error_abort
		implicit none

		integer,intent(in)	:: flag,insertnp
		double precision,intent(in),optional :: Utarget_in
		double precision,dimension(:,:),allocatable,intent(out) :: insert_locs
		!An array of extra molecular positions to include in potential calculation
		double precision,dimension(:,:),allocatable,optional,intent(in)  :: extra_pos

		integer				:: maxattempts=10000
		double precision 	:: Utarget
		logical 			:: pos_found

		!print*, 'usher_get_positions', present(extra_pos)

		!Allocate array to record all inserted locations
		if (allocated(insert_locs)) deallocate(insert_locs)
		allocate(insert_locs(3,insertnp))

		!Loop through and get position of all molecules to insert
		do n = 1, insertnp

			do i = 1, maxattempts

				!Choose starting guess for input position
				select case(flag)
				case(0)
					startpos = get_randompos_procdomain() 
				case(1,3,5)
					startpos = get_randompos_procbottom(nint(0.5*(flag+1)))
				case(2,4,6)
					startpos = get_randompos_proctop(nint(0.5*flag))
				case default
					call error_abort("Error - Usher flag not recognised ")
				end select
				
				!Get target potential energy
				if (present(Utarget_in)) then
					Utarget = Utarget_in
				else
					Utarget = get_Utarget() 
				endif

				!Call USHER to get molecule position
				if (present(extra_pos)) then
					call usher_get_insertion_pos(startpos,Utarget,insertpos,pos_found,Ufinish,extra_pos=extra_pos)
				else
					call usher_get_insertion_pos(startpos,Utarget,insertpos,pos_found,Ufinish)
				endif

				!Insert molecule if suitable position is found
				if (pos_found) then
!					print('(a,2f9.3,a,3f9.3,a,i4,a,i6)'), ' USHER SUCCEEDED: Utarget, Ufinish: ', &
!														Utarget, Ufinish, ', r:',  &
!														insertpos, ', after ', i,  &
!														' attempts , new np', np
					exit
				endif

			enddo

			if (.not. pos_found) then 
				print*, ' usher failed  Utarget =', Utarget
				STOP "ERROR -- For USHER, failure is not an option"
			endif

			!Save inserted location and return
			insert_locs(:,n) = insertpos

		enddo

	end subroutine usher_get_positions

	!Return an array of randomly chosen molecular velocities 
	subroutine get_random_velocity(insertnp,insert_vels)
		implicit none

		integer,intent(in)	:: insertnp
		double precision,dimension(:,:),allocatable,intent(out) :: insert_vels

		integer :: n

		if (allocated(insert_vels)) deallocate(insert_vels)
		allocate(insert_vels(3,insertnp))

		do n = 1,insertnp
			insert_vels(:,n) = sample_MB_vel3()
		enddo

	end subroutine get_random_velocity

	! For re-insert insert, divide evenly between all specified processors
	! and add remainder processorwise with shift used to prevent preferential insertion

	subroutine calculate_mol_proc_redistribution_random(subcomminsertnp,nproc_,rank_,insertnp,  & 
														insertnp_array, subcommmol_to_proc)
		implicit none

		integer, intent(in)									:: subcomminsertnp, nproc_,rank_
		integer, intent(out)								:: insertnp
		integer,allocatable,dimension(:), intent(out)		:: subcommmol_to_proc
		integer,allocatable,dimension(:,:),intent(inout)	:: insertnp_array

		integer 		:: n, molno, proc, extra_insert, insertproc
		integer,save 	:: currentproc = 0

		! For subcomm insert, divide evenly between all processors
		! and add remainder processorwise
		insertnp = floor(subcomminsertnp/dble(nproc_))
		insertnp_array(1,:) = insertnp
		extra_insert = mod(subcomminsertnp,nproc_)
		!Keep shifting the processor used for insert to prevent preferential insertion
		currentproc = currentproc + 1
		do n=1,extra_insert
			insertproc = mod(currentproc+n-1,nproc_)+1
			!print*,rank_,n, 'shiftproc', insertproc,extra_insert
			if (rank_ .eq. insertproc) then
				insertnp = insertnp + 1
			endif
			insertnp_array(1,insertproc) = insertnp_array(1,insertproc) + 1
			!print*, 'EXTRA', rank_, n, insertnp, insertnp_array(1,:)
		enddo

		!Define a mapping between subcomm inserted molecule number and the processor its added to
		allocate(subcommmol_to_proc(subcomminsertnp))
		do proc=1,size(insertnp_array,2)
			do n=1,insertnp_array(1,proc)
				molno = sum(insertnp_array(1,1:proc-1)) + n
				subcommmol_to_proc(molno) = proc
				!print*, n, molno, proc, subcommmol_to_proc(molno),size(subcommmol_to_proc)
			enddo
		enddo


	end subroutine calculate_mol_proc_redistribution_random

	!Check molecular positions to see if any should be reinserted
	subroutine check_insert_overlap(gnip,inp,gi_locs,gm_to_p,rank_,ri_molno,rinp)
		use computational_constants_MD, only : irank, globaldomain, VOID
		use physical_constants_MD, only : rcutoff2
		implicit none

		integer ,intent(in)										:: gnip, inp,rank_
		integer ,intent(out)									:: rinp
		integer,dimension(:),allocatable,intent(in) 			:: gm_to_p
		integer,dimension(:,:),allocatable,intent(out) 			:: ri_molno 
		double precision,dimension(:,:),allocatable,intent(in)	:: gi_locs

		integer	:: i,j,proc,procmolcount
		double precision	:: rij(3), rij2

		!Check on every processor if any inserted molecules is too close to any other
		allocate(ri_molno(2,inp)); ri_molno = VOID
		rinp = 0; proc = 1; procmolcount = 0
		do i=1,gnip
			!Keep track of which processor the current molecule is located
			!if (i .gt. sum(insertnp_array(1,1:proc)) .and. proc .ne. nproc_) proc = proc + 1
			proc = gm_to_p(i)
			!if (3*i .gt. rdisps(proc) .and. proc .ne. nproc_) proc = proc + 1
			if (rank_ .eq. proc) then
				procmolcount = procmolcount + 1
				!print'(a,9i7)', 'PROCESSOR COUNT', rank_, iter, i,inp,gnip, proc, shape(ri_molno), procmolcount
			endif
			do j=i+1,gnip
				rij(:) = gi_locs(:,i)-gi_locs(:,j)
				!If more than half a domain betwen molecules, must be closer over periodic boundaries  
				do ixyz = 1,3    
					if (abs(rij(ixyz)) .gt. 0.5d0*globaldomain(ixyz)) then  !N.B array operator-corresponding element compared 
						rij(ixyz) = rij(ixyz) - sign(globaldomain(ixyz),rij(ixyz))    !Sign of rij applied to domain
					endif
				enddo
				rij2 = dot_product(rij,rij)
				write(300+irank,'(5i5,7f18.9)'), iter, irank,  i,j, gm_to_p(j),  gi_locs(:,i), gi_locs(:,j),rij2
				if (rij2 .lt. rcutoff2) then
					if (rank_ .eq. proc) then
						print'(a,i8,3i5,7f12.7,i6)', 'Inserted Molecular overlap',iter,irank, i,j,sqrt(rij2),gi_locs(:,i),gi_locs(:,j),proc
						rinp = rinp + 1
						ri_molno(1,rinp) = procmolcount	!Local number
						ri_molno(2,rinp) = i			!subcomm number
					endif
				endif
			enddo
		enddo


	end subroutine check_insert_overlap

	

	!Create a random positions on a given processors local domain
	function get_randompos_procdomain() result(pos)
		use computational_constants_MD, only: domain
		implicit none
	
		integer :: ixyz
		double precision :: pos(3), rand
		
		do ixyz = 1,3
			call random_number(rand)
			rand = (rand-0.5)*domain(ixyz)
			pos(ixyz) = rand
		enddo
	
	end function get_randompos_procdomain

	!Create a random positions from the top of a processors local domain
	function get_randompos_proctop(dir) result(pos)
		use computational_constants_MD, only: domain
		use calculated_properties_MD, only : nbins
		implicit none
	
		integer,intent(in) :: dir
		integer			   :: ixyz
		double precision   :: pos(3), rand, binsize
	
		do ixyz = 1,3
			call random_number(rand)
			if (ixyz == dir) then
				binsize = domain(ixyz)/nbins(ixyz)
				rand = 0.5d0*domain(ixyz) - rand*binsize
				pos(ixyz) = rand
			else
				rand = (rand-0.5)*domain(ixyz)
				pos(ixyz) = rand
			endif
		enddo
	
	end function get_randompos_proctop

	!Create a random positions from the bottom of a processors local domain
	function get_randompos_procbottom(dir) result(pos)
		use computational_constants_MD, only: domain
		use calculated_properties_MD, only : nbins
		implicit none
	
		integer,intent(in) :: dir
		integer			   :: ixyz
		double precision   :: pos(3), rand, binsize
	
		do ixyz = 1,3
			call random_number(rand)
			if (ixyz == dir) then
				binsize = domain(ixyz)/nbins(ixyz)
				rand = -0.5d0*domain(ixyz) + rand*binsize
				pos(ixyz) = rand
			else
				rand = (rand-0.5)*domain(ixyz)
				pos(ixyz) = rand
			endif
		enddo
	
	end function get_randompos_procbottom

	! Get maxwell Boltzmann speed
	function sample_MB_vel3() result(vsample)
		use librarymod, only : Maxwell_Boltzmann_vel3
		use calculated_properties_MD, only : temperature
		use physical_constants_MD, only: inputtemperature
		implicit none

		double precision :: vsample(3)
		double precision :: zeromean(3)

		zeromean(:) = 0.d0
		!Insert molecules with the same value as the current temperature
		vsample = Maxwell_Boltzmann_vel3(temperature,zeromean)
		!Insert molecules with the setpoint temperature
		!vsample = Maxwell_Boltzmann_vel3(inputtemperature,zeromean)

	end function sample_MB_vel3 

	! Get potential energy of system
	function get_Utarget() result(U)
		use calculated_properties_MD, only: potenergy
		use physical_constants_MD, only: potential_sLRC
		implicit none
	
		double precision :: U 

		!Insert molecules with the same value as the current potential
		U = potenergy + potential_sLRC  
		!Insert molecules with the setpoint potential
		!U = inputpotential

	end function get_Utarget


!	subroutine usher_insert(flag,insertnp,insert_locs,Utarget_in,v_ins,skipinsert)
!		use interfaces, only : error_abort
!		implicit none

!		integer,intent(in)	:: flag,insertnp
!		double precision,intent(in),optional :: Utarget_in
!		double precision,dimension(:,:),intent(in),optional :: v_ins
!		double precision,dimension(:,:),allocatable,optional,intent(out) :: insert_locs
!		logical,optional,intent(in)	:: skipinsert

!		integer				:: maxattempts=1000
!		double precision 	:: Utarget, insertvel(3)
!		logical 			:: pos_found

!		!Allocate array to record all inserted locations
!		if (allocated(insert_locs)) deallocate(insert_locs)
!		allocate(insert_locs(3,insertnp))

!		!Loop through and insert all molecules
!		do n = 1, insertnp

!			do i = 1, maxattempts

!				!Choose starting guess for input position
!				select case(flag)
!				case(0)
!					startpos = get_randompos_procdomain() 
!				case(2,4,6)
!					startpos = get_randompos_proctop(nint(0.5*flag))
!				case(1,3,5)
!					startpos = get_randompos_procbottom(nint(0.5*(flag+1)))
!				case default
!					call error_abort("Error - Usher flag not recognised ")
!				end select

!				!Get velocity to insert molecule with
!				if (present(v_ins)) then
!					insertvel = v_ins(:,n)
!				else
!					insertvel = sample_MB_vel3()
!				endif
!				
!				!Get target potential energy
!				if (present(Utarget_in)) then
!					Utarget = Utarget_in
!				else
!					Utarget = get_Utarget() 
!				endif

!				!Call USHER to get molecule position
!				call usher_get_insertion_pos(startpos,Utarget,insertpos,pos_found,Ufinish)

!				!Insert molecule if suitable position is found
!				if (pos_found) then
!					call insert_molecule(insertpos,insertvel)	
!					print('(a,2f9.3,a,3f9.3,a,i4,a,i6)'), ' USHER SUCCEEDED: Utarget, Ufinish: ', &
!														Utarget, Ufinish, ', r:',  &
!														insertpos, ', after ', i,  &
!														' attempts , new np', np
!					exit
!				endif

!			enddo

!			if (.not. pos_found) then 
!				print*, ' usher failed  Utarget =', Utarget
!				STOP "ERROR -- For USHER, failure is not an option"
!			endif
!			!Save inserted location and return
!			if (present(insert_locs)) insert_locs(:,n) = insertpos
!		enddo

!	end subroutine usher_insert




end subroutine reinsert_molecules




!! Remove Molecules and then re-insert with same potential energy
!! using the USHER algorithm
!! Code Author: Edward Smith 2013
!subroutine usher_teleport(nparticles)
!	use particle_insertion
!	use computational_constants_MD, only: iter, halfdomain,ncells
!	use physical_constants_MD, only: Potential_sLRC, np
!	use arrays_MD, only: potenergymol_LJ, r, a, v
!	implicit none

!	integer, intent(in) :: nparticles

!	integer :: n,i,maxattempts=1000,icell,jcell,kcell
!	integer :: molno
!	double precision :: Ufinish, rand1, fdummy(3)
!	double precision :: startpos(3), insertpos(3)
!	double precision,dimension(3) :: avevel,vnew,vold,rand
!	double precision,dimension(:),allocatable :: Utarget
!	logical :: pos_found

!	avevel = 0.d0
!	allocate(Utarget(nparticles))
!	do n = 1,nparticles 

!		!Choose molecule at random and store its potential energy as target
!		call random_number(rand1)
!		molno = ceiling(rand1 * np)
!		startpos = r(:,molno)
!		vold = v(:,molno)
!		call compute_force_and_potential_at(startpos,Utarget(n),fdummy)

!		! Remove the molecule 
!		call remove_molecule(molno)

!	!enddo
!	!do n = 1,nparticles 

!		do i = 1, maxattempts

!			startpos = get_randompos() 
!			call usher_get_insertion_pos(startpos,Utarget(n),insertpos,pos_found,Ufinish)

!			if (pos_found) then

!				call insert_molecule(insertpos,vold)	

!				print('(a,2f9.3,a,3f9.3,a,i4,a,i6)'), ' USHER SUCCEEDED: Utarget, Ufinish: ', &
!													Utarget(n), Ufinish, ', r:',  &
!													insertpos, ', after ', i,  &
!													' attempts , new np', np

!				!Write local xy field @ z
!				!call simulation_write_potential_field(max(-halfdomain(1),insertpos(1)-3.0), &
!				!									  min( halfdomain(1),insertpos(1)+3.0), &
!				!									  max(-halfdomain(2),insertpos(2)-3.0), &
!				 !									 min( halfdomain(2),insertpos(2)+3.0), &
!				!									  insertpos(3), 300, 30000+iter)
!				!Write entire xy field @ z
!				!call simulation_write_potential_field(-halfdomain(1),halfdomain(1), &
!				!									  -halfdomain(2),halfdomain(2), &
!				!									   insertpos(3), 400, 30000+iter)

!				exit

!			end if

!		end do

!		if (.not. pos_found) then 
!			print*, ' usher failed  Utarget =', Utarget(n)
!			! Put molecule back where it was if failed
!			call insert_molecule(startpos,vold)
!		end if

!	end do
!	deallocate(Utarget)

!contains

!	function get_Utarget() result(U)
!		use calculated_properties_MD, only: potenergy
!		implicit none
!	
!		double precision :: U 
!		U = potenergy + potential_sLRC  
!	
!	end function get_Utarget
!	
!	function get_randompos() result(pos)
!		use computational_constants_MD, only: domain
!		implicit none
!	
!		integer :: ixyz
!		double precision :: pos(3), rand
!		
!		do ixyz = 1,3
!			call random_number(rand)
!			rand = (rand-0.5)*domain(ixyz)
!			pos(ixyz) = rand
!		end do
!	
!	end function get_randompos

!end subroutine usher_teleport






!! Remove Molecules and then re-insert with same potential energy
!! using the USHER algorithm
!! Code Author: Edward Smith 2013
!subroutine insert_mols(nparticles,Utarget,vold)
!	use particle_insertion
!	use computational_constants_MD, only: iter, halfdomain,ncells
!	use physical_constants_MD, only: Potential_sLRC, np
!	use arrays_MD, only: potenergymol_LJ, r, a, v
!	implicit none

!	integer, intent(in) :: nparticles
!	double precision,dimension(:),allocatable,intent(in) :: Utarget,vold

!	logical  :: pos_found
!	integer :: n,i,maxattempts=10000
!	double precision :: Ufinish 
!	double precision :: startpos(3), insertpos(3)

!	do n = 1,nparticles 
!		do i = 1, maxattempts

!			startpos = get_randompos() 
!			call usher_get_insertion_pos(startpos,Utarget(n),insertpos,pos_found,Ufinish)

!			if (pos_found) then

!				call insert_molecule(insertpos,vold)	

!				print('(a,2f9.3,a,3f9.3,a,i4,a,i6)'), ' USHER SUCCEEDED: Utarget, Ufinish: ', &
!													Utarget(n), Ufinish, ', r:',  &
!													insertpos, ', after ', i,  &
!													' attempts , new np', np
!				exit

!			end if
!		end do

!		if (.not. pos_found) then 
!			print*, ' usher failed  Utarget =', Utarget(n)
!			! Put molecule back where it was if failed
!			call insert_molecule(startpos,vold)
!		end if

!	end do

!contains

!	function get_Utarget() result(U)
!		use calculated_properties_MD, only: potenergy
!		implicit none
!	
!		double precision :: U 
!		U = potenergy + potential_sLRC  
!	
!	end function get_Utarget
!	
!	function get_randompos() result(pos)
!		use computational_constants_MD, only: domain
!		implicit none
!	
!		integer :: ixyz
!		double precision :: pos(3), rand
!		
!		do ixyz = 1,3
!			call random_number(rand)
!			rand = (rand-0.5)*domain(ixyz)
!			pos(ixyz) = rand
!		end do
!	
!	end function get_randompos

!end subroutine insert_mols


!subroutine remove_mols(nparticles,vold,Utarget)
!	use particle_insertion, only : remove_molecule
!	use physical_constants_MD, only: np
!	use arrays_MD, only: r, v
!	implicit none

!	integer, intent(in) :: nparticles
!	double precision,dimension(:),allocatable,intent(out) :: Utarget
!	double precision,dimension(:,:),allocatable,intent(out) :: vold

!	integer :: n
!	integer :: molno
!	double precision :: fdummy(3)
!	double precision :: startpos(3)
!	double precision,dimension(3) :: vnew,rand

!	if (allocated(Utarget)) deallocate(Utarget)

!	allocate(Utarget(nparticles),vold(3,nparticles))
!	do n = 1,nparticles 

!		!Choose molecule at random and store its potential energy as target
!		call random_number(rand)
!		molno = ceiling(rand(1) * np)
!		startpos = r(:,molno)
!		vold(:,n) = v(:,molno)
!		call compute_force_and_potential_at(startpos,Utarget(n),fdummy)

!		! Remove the molecule 
!		call remove_molecule(molno)

!	enddo

!end subroutine remove_mols
