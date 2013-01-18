!----------------------------------------------------------------------------------
!                                Molecular Properties
! Routines for setting and adjusting the properties of the molecules in the domain
!
! Notes : Where statements used for conciseness (slower than do/if combo)
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

subroutine setup_tag
	use module_molecule_properties
	use interfaces, only: error_abort
	use messenger, only: globalise
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
		case ('teth')
			tagdistbottom(:) = tethereddistbottom(:)
			tagdisttop(:)    = tethereddisttop(:)
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

end subroutine setup_tag

!----------------------------------------------------------------------------------
! Check if thermostat is active anywhere in world
subroutine get_tag_thermostat_activity(active)
	use physical_constants_MD, only: np
	use computational_constants_MD, only: thermo_tags
	use arrays_MD, only: tag
	implicit none

	logical, intent(out) :: active
	integer :: i, n

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

subroutine wall_textures(texture_type)
	use module_molecule_properties
	implicit none

	integer 				:: n, icell,jcell
	character				:: texture_type
	double precision		:: xlocation,zlocation

	select case (texture_type)
	case('square')
	!Square ridges
		jcell = 1
		do icell=1,ncells(1),3
			!print*, 'lb=',icell*cellsidelength(1)-halfdomain(1) ,'ub=',(icell+1)*cellsidelength(1)-halfdomain(1)
			do n = 1, np
				if(r(2,n) .gt. (jcell*cellsidelength(2)-halfdomain(2))) then
				if(r(2,n) .lt. ((jcell+1)*cellsidelength(2)-halfdomain(2))) then
					if(r(1,n) .gt. (icell*cellsidelength(1)-halfdomain(1))) then
					if(r(1,n) .lt. ((icell+1)*cellsidelength(1)-halfdomain(1))) then
					!print*, n
						tag(n) = free
					endif
					endif
				endif
				endif
			enddo
		enddo

	case('spikes')
		!Spikes
		do n = 1, np
			xlocation = pi*(r(1,n)+halfdomain(1))/domain(1)
			if (r(2,n) .gt. cellsidelength(2)-halfdomain(2)) then
				!--x direction--
				if (r(2,n) .gt. (0.3d0+sin(5*xlocation)*heaviside(sin(5*xlocation))) & 
						 *5.d0*cellsidelength(2)-halfdomain(2)) tag(n) = free
				if (r(2,n) .gt. (3.d0+sin(10*xlocation)) & 
						 *2.5d0*cellsidelength(2)-halfdomain(2)) tag(n) = free 
				!Add Noise
				if (r(2,n) .gt. (0.3d0+sin(5*xlocation)*heaviside(sin(5*xlocation))+0.3*sin(100*xlocation)) & 
						 *4.d0*cellsidelength(2)-halfdomain(2)) tag(n) = free

				!--z direction--
				if (r(2,n) .gt. (0.3d0+sin(5*zlocation)*heaviside(sin(5*zlocation))) & 
						 *5.d0*cellsidelength(2)-halfdomain(2)) tag(n) = free
				if (r(2,n) .gt. (3.d0+sin(10*zlocation)) & 
						 *2.5d0*cellsidelength(2)-halfdomain(2)) tag(n) = free
				!Add Noise
				if (r(2,n) .gt. (0.3d0+sin(5*zlocation)*heaviside(sin(5*zlocation))+0.3*sin(100*zlocation)) & 
						 *4.d0*cellsidelength(2)-halfdomain(2)) tag(n) = free
			endif
		enddo
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
	case default
		call error_abort("Invalid molecular Tag")
	end select

end subroutine read_tag

!------------------------------------------------------------------------------
subroutine tether_force(molno)
	use module_molecule_properties
	use arrays_MD
	implicit none

	integer                        :: molno
	double precision               :: acctmag
	double precision, dimension(3) :: at, rio
	double precision, parameter    :: k4=5000.d0    !Force constants from...
	double precision, parameter    :: k6=5000000.d0 !...Petravich and Harrowell

	!Obtain displacement from initial position
	rio(:) = r(:,molno) - rtether(:,molno)

	!Calculate applied tethering force
	acctmag = -4.d0*k4*magnitude(rio)**2 - 6.d0*k6*magnitude(rio)**4.d0
	at(:) = acctmag * rio(:)

	!Adjust molecular acceleration accordingly
	a(:,molno) = a(:,molno) + at(:)

	!Adjust initial postion if molecule is sliding
	rtether(:,molno) = rtether(:,molno) + slidev(:,molno)*delta_t

end subroutine tether_force

