!-----------------------------------------------------------------------------
!
!                                  Apply BC     
! Apply periodic boundary conditions 
!                             
!-----------------------------------------------------------------------------

module module_exchange_particles

	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use linked_list

end module module_exchange_particles
!----------------------------------------------------------------------------------

subroutine simulation_exchange_particles
use module_exchange_particles
implicit none

	integer,dimension(3,3,3) :: passarray, recvarray

	!Reset passlist
	call linklist_deallocatepasslist
	passarray = 0 !Zero number of sent molecules
	recvarray = 0 !Zero number of recieved molecules

	call buildpasslist(passarray)
	call buildrecvlist(passarray,recvarray)
	call molsendrecv(recvarray)
	call reorderdata(recvarray)
	call messenger_syncall()

end subroutine simulation_exchange_particles

!------------------------------------------------------------------------------------
!Builds a list of molecules to pass from one processor to adjacent processor

subroutine buildpasslist(passarray)
use module_exchange_particles
implicit none

	integer :: n
	integer :: ipass, jpass, kpass
	integer,dimension(3,3,3) :: passarray

	!Build list of molecules to pass
	do n=1,np 
		!Define an array pass with central element for stationary
		!molecules and surrounding elements detailing direction of passing
		ipass = 0 !Set x passing index to zero - no passing
		jpass = 0 !Set y passign index to zero - no passing
		kpass = 0 !Set z passign index to zero - no passing

		!Check to see if molecules have left the domain and should be passed
		ipass = floor((r(n,1)+halfdomain(1))/domain(1))
		jpass = floor((r(n,2)+halfdomain(2))/domain(2))
		kpass = floor((r(n,3)+halfdomain(3))/domain(3))

		!Add current molecule number to linklist of sent molecules
		if (ipass .ne. 0 .or. jpass .ne. 0 .or. kpass .ne. 0) then
			!print*,irank, ipass, jpass, kpass
			call linklist_checkpushmol(n,ipass,jpass,kpass)
			passarray(2+ipass,2+jpass,2+kpass)= &
			passarray(2+ipass,2+jpass,2+kpass) + 1
		endif

	enddo

end subroutine buildpasslist


!------------------------------------------------------------------------------------
!Builds a list of molecules to pass from one processor to adjacent processor

subroutine buildpasslistcells(passarray)
use module_exchange_particles
implicit none

	integer :: i, n
	integer :: molno,cellnp,nopassed
	integer :: icell, jcell, kcell, ipass, jpass, kpass
	integer,dimension(3,3,3) :: passarray
	type(node), pointer :: old, current

	do n =1,nsurfacecells
		icell = surfacecell(n,1)
		jcell = surfacecell(n,2) 
		kcell = surfacecell(n,3)

		old => cell%head(icell,jcell,kcell)%point 
		cellnp = cell%cellnp(icell,jcell,kcell)
		current => old     !make current point to head of list
	
		!Build list of molecules to pass
		do i=1,cellnp 
			molno = old%molno !Number of molecule

			!Define an array pass with central element for stationary
			!molecules and surrounding elements detailing direction of passing
			ipass = 0 !Set x passing index to zero - no passing
			jpass = 0 !Set y passign index to zero - no passing
			kpass = 0 !Set z passign index to zero - no passing

			!Check to see if molecules have left the domain and should be passed
			ipass = floor((r(molno,1)+halfdomain(1))/domain(1))
			jpass = floor((r(molno,2)+halfdomain(2))/domain(2))
			kpass = floor((r(molno,3)+halfdomain(3))/domain(3))

			!Add current molecule number to linklist of sent molecules
			if (ipass .ne. 0 .or. jpass .ne. 0 .or. kpass .ne. 0) then
				!print*,irank, ipass, jpass, kpass
				call linklist_checkpushmol(molno,ipass,jpass,kpass)
				passarray(2+ipass,2+jpass,2+kpass)= &
				passarray(2+ipass,2+jpass,2+kpass) + 1
			endif

			old => current%next  !make old point to next node of current
			current  => old      !Set current item to old ready for next loop
		enddo
	enddo
	nullify (current)   !Set current pointer to null
	nullify (old)       !Set old pointer to null

end subroutine buildpasslistcells


!------------------------------------------------------------------------------------
!Builds a list of molecules to pass from one processor to adjacent processor

subroutine buildrecvlist(passarray,recvarray)
use module_exchange_particles
implicit none

	integer :: ipass, jpass, kpass
	integer,dimension(3,3,3) :: passarray, recvarray

	!Receive molecules
	do ipass = -1,1
	do jpass = -1,1
	do kpass = -1,1
		call arraysendrecv(ipass,jpass,kpass,passarray,recvarray)
	enddo
	enddo
	enddo

end subroutine buildrecvlist

!------------------------------------------------------------------------------------
!Move through molecule list of passed molecule and fill with new arrivals starting 
!from the end of the new molecules added to r this timestep

subroutine reorderdata(recvarray)
use module_exchange_particles
implicit none

	integer			   :: i
	integer			   :: molno,nopassed,recvcount
	integer,dimension(3,3,3)   :: recvarray
	type(passnode), pointer    :: oldp, currentp

	recvcount = sum(recvarray)

	!Work through passed molecules to fill gaps
	oldp => pass%head
	nopassed = pass%nopassed
	currentp => oldp     !make current point to head of list

	!Loop through all empty locations
	do i=1,nopassed

		molno = oldp%molno	!Position of empty location 

		r(molno,:) = r(np+recvcount,:)
		v(molno,:) = v(np+recvcount,:)

		!print*, irank, 'mol no',np+recvcount, 'into gap', molno

		!Update number of new molecules
		recvcount = recvcount - 1

		oldp => currentp%next  !make old point to next node of current
		currentp  => oldp      !Set current item to old ready for next loop

	enddo

	!Adjust total number of molecules to reflect difference between
	!recieved molecules and passed molecules
	np = np + recvcount

	nullify (currentp)   !Set current pointer to null
	nullify (oldp)       !Set old pointer to null

	!print*, irank, 'number of particles', np

end subroutine reorderdata

!------------------------------------------------------------------------------------
!Routine for passing molecules from one processor to adjacent processor

!subroutine passmolecules()
!use module_exchange_particles
!implicit none

!	integer :: i, nopassed
!	integer :: ipass, jpass, kpass

	!Allocate an array pass with central element for stationary
	!molecules and surrounding elements detailing direction of passing
!	do i=1,np 

!		ipass = 0 !Set x passing index to zero - no passing
!		jpass = 0 !Set y passign index to zero - no passing
!		kpass = 0 !Set z passign index to zero - no passing

		!Check to see if molecules have left the domain and should be passed
!		ipass = floor((r(i,1)+halfdomain(1)))
!		jpass = floor((r(i,2)+halfdomain(2)))
!		kpass = floor((r(i,3)+halfdomain(3)))

		!Record passing
		!pass%nopassed(2+ipass,2+jpass,2+kpass)= pass%nopassed(2+ipass,2+jpass,2+kpass)+1

!		if (ipass .ne. 0 .or. jpass .ne. 0 .or. kpass .ne. 0) then
!		!Add current molecule number to linklist of sent molecules
!			call linklist_checkpushmol(molno, nopassed)
!		endif

!	enddo

!	recvcount = 0 !Zero number of recieved molecules

!	do i=1,nopassed

		!Update number of new molecules
!		recvcount = recvcount + 1

		!Send molecule to correct processors
!		call molsendrecv(molno,ipass,jpass,kpass,r(np+recvcount,:),v(np+recvcount,:))

!	enddo
!		endif

!	enddo

!end subroutine passmolecules

