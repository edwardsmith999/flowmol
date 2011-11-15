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

end module module_molecule_properties

!----------------------------------------------------------------------------------
!Create an array of molecular tags

subroutine setup_tag
	use module_molecule_properties
	implicit none

	integer 			:: n, mol_layers
	!double precision		:: xlocation,zlocation, heaviside


	!For initialunitsize "a"
	!		 [  o     o ]
	!a (1 cell size) [     o    ]  a/2 (distance between molcules)	
	!		 [  o     o
	!		  __________]  a/4 (distance from bottom of domain)
	!
	!So use (0.20+0.5d0*mol_layers)*initialunitsize(ixyz)

	mol_layers = 2

	wallslidev(1) = 0.d0
	wallslidev(2) = 0.d0
	wallslidev(3) = 0.d0

	!Setup fixed molecules
	!+0.2 to include gap between wall and 1st molecule
	fixdistbottom(1) = 0.d0*cellsidelength(1) !initialunitsize(1)
	fixdistbottom(2) = 1.d0*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	fixdistbottom(3) = 0.d0*cellsidelength(3) !initialunitsize(3)
	fixdisttop(1) = 0.d0*cellsidelength(1) !initialunitsize(1)
	fixdisttop(2) = 0.d0*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	fixdisttop(3) = 0.d0*cellsidelength(3) !initialunitsize(3)

	!Setup sliding molecules
	slidedistbottom(1) = 0.d0*cellsidelength(1) !initialunitsize(1)
	slidedistbottom(2) = 0.d0*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	slidedistbottom(3) = 0.d0*cellsidelength(3) !initialunitsize(3)
	slidedisttop(1) = 0.d0 !initialunitsize(1)
	slidedisttop(2) = 0.d0*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	slidedisttop(3) = 0.d0 !initialunitsize(3)

	!Setup molecules with tethered potentials
	tethereddistbottom(1) = 0.d0 !initialunitsize(1)
	tethereddistbottom(2) = 0.d0*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	tethereddistbottom(3) = 0.d0 !initialunitsize(3)
	tethereddisttop(1) = 0.d0 !initialunitsize(1)
	tethereddisttop(2) = 0.d0*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	tethereddisttop(3) = 0.d0 !initialunitsize(3)

	!Setup thermostatted molecules
	thermstatbottom(1) = 0.d0 !initialunitsize(1)
	thermstatbottom(2) = 0.d0*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	thermstatbottom(3) = 0.d0 !initialunitsize(3)
	thermstattop(1) = 0.d0 !initialunitsize(1)
	thermstattop(2) = 0.d0*(ncells(2)-1)*cellsidelength(2) !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	thermstattop(3) = 0.d0 !initialunitsize(3)

	!Set all molecules tag equal to zero (normal)
	tag = 0

	do n = 1,np
		!x bottom
		if (iblock .eq. 1) then
			if(r(n,1).lt.-halfdomain(1)+thermstatbottom(1)) 	tag(n) = 4
			if(r(n,1).lt.-halfdomain(1)+tethereddistbottom(1)) 	tag(n) = 3
			if(r(n,1).lt.-halfdomain(1)+tethereddistbottom(1) & 
		     .and. r(n,1).lt.-halfdomain(1)+thermstatbottom(1)) 	tag(n) = 5
			if(r(n,1).lt.-halfdomain(1)+fixdistbottom(1))		tag(n) = 1
			if(r(n,1).lt.-halfdomain(1)+slidedistbottom(1)) 	tag(n) = 2
			if(r(n,1).lt.-halfdomain(1)+tethereddistbottom(1) & 
		     .and. r(n,1).lt.-halfdomain(1)+slidedistbottom(1))		tag(n) = 6
			if(r(n,1).lt.-halfdomain(1)+tethereddistbottom(1) & 
		     .and. r(n,1).lt.-halfdomain(1)+thermstatbottom(1)    &
		     .and. r(n,1).lt.-halfdomain(1)+slidedistbottom(1))		tag(n) = 7
		endif

		!x top	
		if (iblock .eq. npx) then
			if(r(n,1).ge. halfdomain(1)-thermstattop(1)) 		tag(n) = 4
			if(r(n,1).ge. halfdomain(1)-tethereddisttop(1)) 	tag(n) = 3
			if(r(n,1).gt. halfdomain(1)-tethereddisttop(1) 	& 
		     .and. r(n,1).gt. halfdomain(1)-thermstattop(1)) 		tag(n) = 5
			if(r(n,1).ge. halfdomain(1)-fixdisttop(1)) 		tag(n) = 1
			if(r(n,1).ge. halfdomain(1)-slidedisttop(1)) 		tag(n) = 2
			if(r(n,1).gt. halfdomain(1)-tethereddisttop(1)	& 
		     .and. r(n,1).gt. halfdomain(1)-slidedisttop(1))		tag(n) = 6
			if(r(n,1).gt. halfdomain(1)-tethereddisttop(1)	& 
		     .and. r(n,1).gt. halfdomain(1)-thermstattop(1)   	&
		     .and. r(n,1).gt. halfdomain(1)-slidedisttop(1))		tag(n) = 7
		endif

		!y bottom
		if (jblock .eq. 1) then
			if(r(n,2).lt.-halfdomain(2)+thermstatbottom(2)) 	tag(n) = 4
			if(r(n,2).lt.-halfdomain(2)+tethereddistbottom(2)) 	tag(n) = 3
			if(r(n,2).lt.-halfdomain(2)+tethereddistbottom(2) & 
		     .and. r(n,2).lt.-halfdomain(2)+thermstatbottom(2) ) 	tag(n) = 5
			if(r(n,2).lt.-halfdomain(2)+fixdistbottom(2))		tag(n) = 1
			if(r(n,2).lt.-halfdomain(2)+slidedistbottom(2)) 	tag(n) = 2
			if(r(n,2).lt.-halfdomain(2)+tethereddistbottom(2) & 
		     .and. r(n,2).lt.-halfdomain(2)+slidedistbottom(2))		tag(n) = 6
			if(r(n,2).lt.-halfdomain(2)+tethereddistbottom(2) & 
		     .and. r(n,2).lt.-halfdomain(2)+thermstatbottom(2)    &
		     .and. r(n,2).lt.-halfdomain(2)+slidedistbottom(2))		tag(n) = 7
		endif
	
		!y top
		if (jblock .eq. npy) then
			if(r(n,2).ge. halfdomain(2)-thermstattop(2)) 		tag(n) = 4
			if(r(n,2).ge. halfdomain(2)-tethereddisttop(2)) 	tag(n) = 3
			if(r(n,2).gt. halfdomain(2)-tethereddisttop(2) & 
		     .and. r(n,2).gt. halfdomain(2)-thermstattop(2) ) 		tag(n) = 5
			if(r(n,2).ge. halfdomain(2)-fixdisttop(2)) 		tag(n) = 1
			if(r(n,2).ge. halfdomain(2)-slidedisttop(2)) 		tag(n) = 2
			if(r(n,2).gt. halfdomain(2)-tethereddisttop(2) & 
		     .and. r(n,2).gt. halfdomain(2)-slidedisttop(2))		tag(n) = 6 
			if(r(n,2).gt. halfdomain(2)-tethereddisttop(2) & 
		     .and. r(n,2).gt. halfdomain(2)-thermstattop(2)    &
		     .and. r(n,2).gt. halfdomain(2)-slidedisttop(2))		tag(n) = 7 
		endif

		!z bottom
		if (kblock .eq. 1) then
			if(r(n,3).lt.-halfdomain(3)+thermstatbottom(3)) 	tag(n) = 4
			if(r(n,3).lt.-halfdomain(3)+tethereddistbottom(3)) 	tag(n) = 3
			if(r(n,3).lt.-halfdomain(3)+tethereddistbottom(3) & 
		     .and. r(n,3).lt.-halfdomain(3)+thermstatbottom(3) ) 	tag(n) = 5
			if(r(n,3).lt.-halfdomain(3)+fixdistbottom(3))		tag(n) = 1
			if(r(n,3).lt.-halfdomain(3)+slidedistbottom(3)) 	tag(n) = 2
			if(r(n,3).lt.-halfdomain(3)+tethereddistbottom(3) & 
		     .and. r(n,3).lt.-halfdomain(3)+slidedistbottom(3))		tag(n) = 6
			if(r(n,3).lt.-halfdomain(3)+tethereddistbottom(3) & 
		     .and. r(n,3).lt.-halfdomain(3)+thermstatbottom(3)    &
		     .and. r(n,3).lt.-halfdomain(3)+slidedistbottom(3))		tag(n) = 7
		endif

		!z top
		if (kblock .eq. npz) then
			if(r(n,3).ge. halfdomain(3)-thermstattop(3)) 		tag(n) = 4
			if(r(n,3).ge. halfdomain(3)-tethereddisttop(3)) 	tag(n) = 3
			if(r(n,3).gt. halfdomain(3)-tethereddisttop(3) & 
		     .and. r(n,3).gt. halfdomain(3)-thermstattop(3) ) 		tag(n) = 5
			if(r(n,3).ge. halfdomain(3)-fixdisttop(3)) 		tag(n) = 1
			if(r(n,3).ge. halfdomain(3)-slidedisttop(3)) 		tag(n) = 2
			if(r(n,3).gt. halfdomain(3)-tethereddisttop(3) & 
		     .and. r(n,3).gt. halfdomain(3)-slidedisttop(3))		tag(n) = 6
			if(r(n,3).gt. halfdomain(3)-tethereddisttop(3) & 
		     .and. r(n,3).gt. halfdomain(3)-thermstattop(3)    &
		     .and. r(n,3).gt. halfdomain(3)-slidedisttop(3))		tag(n) = 7
		endif

	enddo

	!print*, tag

	!Square ridges
	!jcell = 1
	!do icell=1,ncells(1),3
	!	print*, 'lb=',icell*cellsidelength(1)-halfdomain(1) ,'ub=',(icell+1)*cellsidelength(1)-halfdomain(1)
	!	do n = 1, np
	!		if(r(n,2) .gt. (jcell*cellsidelength(2)-halfdomain(2))) then
	!		if(r(n,2) .lt. ((jcell+1)*cellsidelength(2)-halfdomain(2))) then
	!			if(r(n,1) .gt. (icell*cellsidelength(1)-halfdomain(1))) then
	!			if(r(n,1) .lt. ((icell+1)*cellsidelength(1)-halfdomain(1))) then
	!				!print*, n
	!				tag(n) = 0
	!			endif
	!			endif
	!		endif
	!		endif
	!	enddo
	!enddo
	!VMD colours y>-5 and y <-3.5 and x <X

	!Spikes
	!do n = 1, np
	!	xlocation = pi*(r(n,1)+halfdomain(1))/domain(1)
	!	if (r(n,2) .gt. cellsidelength(2)-halfdomain(2)) then
			!No negative
			!if (r(n,2) .gt. (0.3d0+sin(5*xlocation)*heaviside(sin(5*xlocation))) & 
			!		 *5.d0*cellsidelength(2)-halfdomain(2)) tag(n) = 0
			!if (r(n,2) .gt. (3.d0+sin(10*xlocation)) & 
			!		 *2.5d0*cellsidelength(2)-halfdomain(2)) tag(n) = 0
			!Noise
			!if (r(n,2) .gt. (0.3d0+sin(5*xlocation)*heaviside(sin(5*xlocation))+0.3*sin(100*xlocation)) & 
			!		 *4.d0*cellsidelength(2)-halfdomain(2)) tag(n) = 0
	!	endif
	!enddo

	!do n = 1, np
	!	zlocation = pi*(r(n,3)+halfdomain(3))/domain(3)
	!	if (r(n,2) .gt. cellsidelength(2)-halfdomain(2)) then
			!No negative
			!if (r(n,2) .gt. (0.3d0+sin(5*zlocation)*heaviside(sin(5*zlocation))) & 
			!		 *5.d0*cellsidelength(2)-halfdomain(2)) tag(n) = 0
			!if (r(n,2) .gt. (3.d0+sin(10*zlocation)) & 
			!		 *2.5d0*cellsidelength(2)-halfdomain(2)) tag(n) = 0
			!Noise
			!if (r(n,2) .gt. (0.3d0+sin(5*zlocation)*heaviside(sin(5*zlocation))+0.3*sin(100*zlocation)) & 
			!		 *4.d0*cellsidelength(2)-halfdomain(2)) tag(n) = 0
	!	endif
	!enddo
	!VMD colours y >4.38*sin(3*3.141*(x+10.22)/20.435)-10.22 and y<-4.0 and y>-9

end subroutine setup_tag

!----------------------------------------------------------------------------------
!Read molecular tags using properties defined by different place values

subroutine decode_tag(molno)
	use module_molecule_properties
	implicit none

	integer			:: l,n,m
	integer 		:: molno, no_2_decode
	integer			:: fixed, forced, thermostated
	integer,dimension(5)	:: check, place_value

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
	use module_molecule_properties
	implicit none

	integer :: molno

	!Check tag and assign properties accordingly
	select case (tag(molno))
	case (0)
		!Set molecules to unfixed with no sliding velocity
		fix(molno,:) = 1
		slidev(molno,:) = 0.d0
	case (1)
		!Fixed Molecules
		fix(molno,:) = 0
		slidev(molno,:) = 0.d0
	case (2)
		!Fixed with constant sliding speed
		fix(molno,:) = 0
		slidev(molno,:) = wallslidev
	case (3)
		!Tethered molecules unfixed with no sliding velocity
		fix(molno,:) = 1
		slidev(molno,:) = 0.d0
	case (4)
		!Thermostatted molecules 
		fix(molno,:) = 1
		slidev(molno,:) = 0.d0
		thermostat(molno,:) = 1
	case (5)
		!Thermostatted Tethered molecules unfixed with no sliding velocity
		fix(molno,:) = 1
		slidev(molno,:) = 0.d0
		thermostat(molno,:) = 1
	case (6)
		!Tethered molecules with sliding velocity
		fix(molno,:) = 1
		slidev(molno,:) = wallslidev
	case (7)
		!Thermostatted Tethered molecules unfixed with sliding velocity
		fix(molno,:) = 1
		slidev(molno,:) = wallslidev
		thermostat(molno,:) = 1
	case default
		stop "Invalid molecular Tag"
	end select

end subroutine read_tag

!-----------------------------------------------------------------------

subroutine tether_force(molno)
	use module_molecule_properties
	use arrays_MD
	implicit none

	integer 			:: molno, ixyz
	double precision		:: k4, k6
	double precision		:: magnitude,acctmag
	double precision, dimension(3)	:: at, rio

	!Define strength of tethering potential ~ phi= -k4*rio^4 - k6*rio^6
	!Constants from Petravic and Harrowell (k4 = 5,000, k6 = 5,000,000)
	k4 = 5000.d0
	k6 = 5000000.d0

	!Obtain displacement from initial position
	rio(:) = r(molno,:) - rinitial(molno,:)
	!For serial code, the nearest neighbour convention could be used. This becomes more complicated in parallel
	!requiring rinitial to be passed to the adjacent processor.
	do ixyz = 1,3
		if(abs(rio(ixyz)) .gt. halfdomain(ixyz)) rio(ixyz)=rio(ixyz)-sign(domain(ixyz),rio(ixyz))
	enddo

	!Calculate applied tethering force
	acctmag = -4.d0*k4*magnitude(rio)**2 - 6.d0*k6*magnitude(rio)**4.d0
	at(:) = acctmag * rio(:)

	!Adjust molecular acceleration accordingly
	a(molno,:) = a(molno,:) + at(:)

	!Adjust initial postion if molecule is sliding
	rinitial(molno,:) = rinitial(molno,:) + slidev(molno,:)*delta_t
	do ixyz = 1,3
		if(abs(rinitial(molno,ixyz)) .gt. halfdomain(ixyz)) then
			rinitial(molno,ixyz)=rinitial(molno,ixyz)-sign(domain(ixyz),rinitial(molno,ixyz))
		endif
	enddo

	!print'(i8,9f10.5)', molno, rinitial(molno,:), r(molno,:), at(:)
	!if (molno .eq. 1311) print'(i8,9f10.5)', molno, rinitial(molno,:), r(molno,:),  rio(:)
	!if (maxval(a(molno,:)) .gt. 1000.d0) print'(i8,9f10.5)', molno, rinitial(molno,:), r(molno,:),  rio(:)

end subroutine tether_force

