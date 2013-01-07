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
	implicit none

	integer 				:: ixyz, n, mol_layers
	integer,dimension(3)	:: block, npxyz

	block = (/ iblock, jblock, kblock /)
	npxyz = (/  npx  ,  npy  ,  npz   /)
	mol_layers = 2

	!Initialise
	tag(:) = free 

	if (ensemble .eq. tag_move) then
		!Setup fixed wall and location dependent thermostat tags
		do n = 1,np
			do ixyz = 1,3
				!Bottom
				if (block(ixyz) .eq. 1) then
					if(r(ixyz,n).lt.-halfdomain(ixyz)+thermstatbottom(ixyz)) 	tag(n) = thermo 
					if(r(ixyz,n).lt.-halfdomain(ixyz)+tethereddistbottom(ixyz))	tag(n) = teth
					if(r(ixyz,n).lt.-halfdomain(ixyz)+tethereddistbottom(ixyz) 		& 
				 .and. r(ixyz,n).lt.-halfdomain(ixyz)+thermstatbottom(ixyz)) 	tag(n) = teth_thermo 
					if(r(ixyz,n).lt.-halfdomain(ixyz)+fixdistbottom(ixyz))		tag(n) = fixed 
					if(r(ixyz,n).lt.-halfdomain(ixyz)+slidedistbottom(ixyz)) 	tag(n) = fixed_slide 
					if(r(ixyz,n).lt.-halfdomain(ixyz)+tethereddistbottom(ixyz) 		& 
				 .and. r(ixyz,n).lt.-halfdomain(ixyz)+slidedistbottom(ixyz))	tag(n) = teth_slide 
					if(r(ixyz,n).lt.-halfdomain(ixyz)+tethereddistbottom(ixyz) 		& 
				 .and. r(ixyz,n).lt.-halfdomain(ixyz)+thermstatbottom(ixyz)	&
				 .and. r(ixyz,n).lt.-halfdomain(ixyz)+slidedistbottom(ixyz))	tag(n) = teth_thermo_slide
				endif

				!Top	
				if (block(ixyz) .eq. npxyz(ixyz)) then
					if(r(ixyz,n).ge. halfdomain(ixyz)-thermstattop(ixyz)) 		tag(n) = thermo 
					if(r(ixyz,n).ge. halfdomain(ixyz)-tethereddisttop(ixyz)) 	tag(n) = teth 
					if(r(ixyz,n).gt. halfdomain(ixyz)-tethereddisttop(ixyz) 	& 
					 .and. r(ixyz,n).gt. halfdomain(ixyz)-thermstattop(ixyz)) 	tag(n) = teth_thermo 
					if(r(ixyz,n).ge. halfdomain(ixyz)-fixdisttop(ixyz)) 		tag(n) = fixed 
					if(r(ixyz,n).ge. halfdomain(ixyz)-slidedisttop(ixyz)) 		tag(n) = fixed_slide 
					if(r(ixyz,n).gt. halfdomain(ixyz)-tethereddisttop(ixyz)		& 
					 .and. r(ixyz,n).gt. halfdomain(ixyz)-slidedisttop(ixyz))	tag(n) = teth_slide 
					if(r(ixyz,n).gt. halfdomain(ixyz)-tethereddisttop(ixyz)		& 
					 .and. r(ixyz,n).gt. halfdomain(ixyz)-thermstattop(ixyz)   	&
					 .and. r(ixyz,n).gt. halfdomain(ixyz)-slidedisttop(ixyz))	tag(n) = teth_thermo_slide
				endif
			enddo
		enddo
	end if

!	do n = 1,np
		!x bottom
!		if (iblock .eq. 1) then
!			if(r(1,n).lt.-halfdomain(1)+   thermstatbottom(1)) 	tag(n) = 4
!			if(r(1,n).lt.-halfdomain(1)+tethereddistbottom(1)) 	tag(n) = 3
!			if(r(1,n).lt.-halfdomain(1)+tethereddistbottom(1) & 
!		     .and. r(1,n).lt.-halfdomain(1)+   thermstatbottom(1)) 	tag(n) = 5
!			if(r(1,n).lt.-halfdomain(1)+     fixdistbottom(1))	tag(n) = 1
!			if(r(1,n).lt.-halfdomain(1)+   slidedistbottom(1)) 	tag(n) = 2
!			if(r(1,n).lt.-halfdomain(1)+tethereddistbottom(1) & 
!		     .and. r(1,n).lt.-halfdomain(1)+   slidedistbottom(1))	tag(n) = 6
!			if(r(1,n).lt.-halfdomain(1)+tethereddistbottom(1) & 
!		     .and. r(1,n).lt.-halfdomain(1)+   thermstatbottom(1)    &
!		     .and. r(1,n).lt.-halfdomain(1)+   slidedistbottom(1))	tag(n) = 7
!		endif

!		!x top	
!		if (iblock .eq. npx) then
!			if(r(1,n).ge. halfdomain(1)-thermstattop(1)) 		tag(n) = 4
!			if(r(1,n).ge. halfdomain(1)-tethereddisttop(1)) 	tag(n) = 3
!			if(r(1,n).gt. halfdomain(1)-tethereddisttop(1) 	& 
!		     .and. r(1,n).gt. halfdomain(1)-thermstattop(1)) 		tag(n) = 5
!			if(r(1,n).ge. halfdomain(1)-fixdisttop(1)) 		tag(n) = 1
!			if(r(1,n).ge. halfdomain(1)-slidedisttop(1)) 		tag(n) = 2
!			if(r(1,n).gt. halfdomain(1)-tethereddisttop(1)	& 
!		     .and. r(1,n).gt. halfdomain(1)-slidedisttop(1))		tag(n) = 6
!			if(r(1,n).gt. halfdomain(1)-tethereddisttop(1)	& 
!		     .and. r(1,n).gt. halfdomain(1)-thermstattop(1)   	&
!		     .and. r(1,n).gt. halfdomain(1)-slidedisttop(1))		tag(n) = 7
!		endif

		!y bottom
!		if (jblock .eq. 1) then
!			if(r(2,n).lt.-halfdomain(2)+thermstatbottom(2)) 	tag(n) = 4
!			if(r(2,n).lt.-halfdomain(2)+tethereddistbottom(2)) 	tag(n) = 3
!			if(r(2,n).lt.-halfdomain(2)+tethereddistbottom(2) & 
!		     .and. r(2,n).lt.-halfdomain(2)+thermstatbottom(2) ) 	tag(n) = 5
!			if(r(2,n).lt.-halfdomain(2)+fixdistbottom(2))		tag(n) = 1
!			if(r(2,n).lt.-halfdomain(2)+slidedistbottom(2)) 	tag(n) = 2
!			if(r(2,n).lt.-halfdomain(2)+tethereddistbottom(2) & 
!		     .and. r(2,n).lt.-halfdomain(2)+slidedistbottom(2))		tag(n) = 6
!			if(r(2,n).lt.-halfdomain(2)+tethereddistbottom(2) & 
!		     .and. r(2,n).lt.-halfdomain(2)+thermstatbottom(2)    &
!		     .and. r(2,n).lt.-halfdomain(2)+slidedistbottom(2))		tag(n) = 7
	!	endif
	
		!y top
	!	if (jblock .eq. npy) then
	!		if(r(2,n).ge. halfdomain(2)-thermstattop(2)) 		tag(n) = 4
	!		if(r(2,n).ge. halfdomain(2)-tethereddisttop(2)) 	tag(n) = 3
	!		if(r(2,n).gt. halfdomain(2)-tethereddisttop(2) 	& 
	!	     .and. r(2,n).gt. halfdomain(2)-thermstattop(2) ) 	tag(n) = 5
	!		if(r(2,n).ge. halfdomain(2)-fixdisttop(2)) 			tag(n) = 1
	!		if(r(2,n).ge. halfdomain(2)-slidedisttop(2)) 		tag(n) = 2
!!			if(r(2,n).gt. halfdomain(2)-tethereddisttop(2) 	& 
!		     .and. r(2,n).gt. halfdomain(2)-slidedisttop(2))	tag(n) = 6 
!			if(r(2,n).gt. halfdomain(2)-tethereddisttop(2) 	& 
!		     .and. r(2,n).gt. halfdomain(2)-thermstattop(2) &
!		     .and. r(2,n).gt. halfdomain(2)-slidedisttop(2))	tag(n) = 7 
!		endif

		!z bottom
!		if (kblock .eq. 1) then
!			if(r(3,n).lt.-halfdomain(3)+thermstatbottom(3)) 	tag(n) = 4
!			if(r(3,n).lt.-halfdomain(3)+tethereddistbottom(3)) 	tag(n) = 3
!			if(r(3,n).lt.-halfdomain(3)+tethereddistbottom(3) & 
!		     .and. r(3,n).lt.-halfdomain(3)+thermstatbottom(3) ) 	tag(n) = 5
!			if(r(3,n).lt.-halfdomain(3)+fixdistbottom(3))		tag(n) = 1
!			if(r(3,n).lt.-halfdomain(3)+slidedistbottom(3)) 	tag(n) = 2
!			if(r(3,n).lt.-halfdomain(3)+tethereddistbottom(3) & 
!		     .and. r(3,n).lt.-halfdomain(3)+slidedistbottom(3))		tag(n) = 6
!			if(r(3,n).lt.-halfdomain(3)+tethereddistbottom(3) & 
!		     .and. r(3,n).lt.-halfdomain(3)+thermstatbottom(3)    &
!		     .and. r(3,n).lt.-halfdomain(3)+slidedistbottom(3))		tag(n) = 7
!		endif

		!z top
!		if (kblock .eq. npz) then
!			if(r(3,n).ge. halfdomain(3)-thermstattop(3)) 		tag(n) = 4
!			if(r(3,n).ge. halfdomain(3)-tethereddisttop(3)) 	tag(n) = 3
!			if(r(3,n).gt. halfdomain(3)-tethereddisttop(3) & 
!		     .and. r(3,n).gt. halfdomain(3)-thermstattop(3) ) 		tag(n) = 5
!			if(r(3,n).ge. halfdomain(3)-fixdisttop(3)) 		tag(n) = 1
!			if(r(3,n).ge. halfdomain(3)-slidedisttop(3)) 		tag(n) = 2
!			if(r(3,n).gt. halfdomain(3)-tethereddisttop(3) & 
!		     .and. r(3,n).gt. halfdomain(3)-slidedisttop(3))		tag(n) = 6
!			if(r(3,n).gt. halfdomain(3)-tethereddisttop(3) & 
!		     .and. r(3,n).gt. halfdomain(3)-thermstattop(3)    &
!		     .and. r(3,n).gt. halfdomain(3)-slidedisttop(3))		tag(n) = 7
!		endif
!
!	enddo

end subroutine setup_tag

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

