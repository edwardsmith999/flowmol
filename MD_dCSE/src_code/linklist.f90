!-----------------------------------------------------------------------
! 				Linklist
!                                
! Linklist subroutines used to build, manipulate and output linked lists
!
!======================================================================
!		Cell and neighbour building subroutines               =
!======================================================================
! Subroutines used to build cell and neighbour lists
! assigntoneighbourlist()
! assigntocell()
! assigntohalocell()
!!======================================================================
!		Linklist manipulation Subroutines                     =
!======================================================================
! Subroutines used throughout code for all linklist operations
!
! linklist_build(start,finish) !Build mol no. start to mol no. finish
! linklist_circbuild(start,finish) !As build but a circular link list 
! linklist_pop(icell, jcell, kcell)
! linklist_push(icell, jcell, kcell, molno)
! linklist_checkpush(icell, jcell, kcell, molno)
! linklist_checkpushneighbr(molnoi, molnoj)
! linklist_print(icell, jcell, kcell)
! linklist_printallcells()
! linklist_printalldomaincells()
! linklist_deallocate(icell, jcell, kcell)
! linklist_deallocateall()
! linklist_gototop(icell, jcell, kcell)
! linklist_gotobottom(icell, jcell, kcell)
! linklist_gotomolecule(icell, jcell, kcell, n) !N.B. molecule must be in cell
! linklist_movethrough(icell, jcell, kcell, n,dir) !dir='goup' OR 'down'	 
! linklist_printcurrent(icell, jcell, kcell)
!
!-------------------------------------------------------------------------------

module module_linklist

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use linked_list

	integer, dimension(6,6,6)	:: intcount

end module module_linklist

!======================================================================
!		Cell and neighbour building subroutines               =
!======================================================================

!--------------------------------------------------------
!Routine used to rebuild linked list for each cell

subroutine assign_to_cell
use module_linklist
implicit none

	integer		:: n
	integer		:: icell, jcell, kcell

	do n=1,np
		icell = ceiling((r(n,1)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add 1 due to halo
		jcell = ceiling((r(n,2)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add 1 due to halo
		kcell = ceiling((r(n,3)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add 1 due to halo
		!if (irank .eq. iroot) print'(5i5,3f10.5)', irank, &
		!n, icell, jcell, kcell, r(n,1), r(n,2), r(n,3)
		call linklist_checkpush(icell, jcell, kcell, n)
	enddo

end subroutine assign_to_cell

!--------------------------------------------------------
!Routine used to rebuild linked list for each halo cell

subroutine assign_to_halocell
use module_linklist
implicit none

	integer		:: n, start, finish
	integer		:: icell, jcell, kcell

	start = np+1
	finish = np+halo_np

	do n=start,finish
		icell = ceiling((r(n,1)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add 1 due to halo
		jcell = ceiling((r(n,2)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add 1 due to halo
		kcell = ceiling((r(n,3)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add 1 due to halo
		!if (irank .eq. iroot) print'(4i5,3f10.5)', &
		!n, icell, jcell, kcell, r(n,1), r(n,2), r(n,3)
		call linklist_checkpush(icell, jcell, kcell, n)
	enddo

end subroutine assign_to_halocell

!----------------------------------------------------------------------------------
! Assign to Neighbourlist                               
! Build a list connecting molecules which are within interacting distance of
! each other ~ re-ordered so that each all i and all j interactions are calculated
! for each cell pair before moving to next cell pair

subroutine assign_to_neighbourlist
use module_linklist
implicit none

	integer                         :: i, j, ixyz   !Define dummy index
	integer				:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp
	integer				:: molnoi, molnoj
	double precision		:: rij2   !magnitude^2 between i and j
	double precision,dimension(3)   :: ri, rj !Position of molecule i and j
	double precision,dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	!Create Neighbourlist array based on current np
	allocate(neighbour%noneighbrs(np))
	allocate(neighbour%head(np))
	do i = 1,np
		neighbour%noneighbrs(i) = 0	!Zero number of molecules in neighbour list
		nullify(neighbour%head(i)%point)!Nullify neighbour list head pointer 
	enddo

	!if (maxval(cell%cellnp(:,:,:)) .gt. 9) stop "ERROR - greater than 10 per cell"

	do icell=2, ncells(1)+1
	do jcell=2, ncells(2)+1
	do kcell=2, ncells(3)+1

		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp                  !Step through each particle in list 
			molnoi = oldi%molno 	 !Number of molecule
			ri = r(molnoi,:)         !Retrieve ri

			do icellshift = -1,1
			do jcellshift = -1,1
			do kcellshift = -1,1
				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				!print*, icell+icellshift,jcell+jcellshift,kcell+kcellshift

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	 !Number of molecule
					rj = r(molnoj,:)         !Retrieve rj
					currentj => oldj
					oldj => currentj%next    !Use pointer in datatype to obtain next item in list

					if(molnoi==molnoj) cycle !Check to prevent interaction with self

					rij2=0                   !Set rij^2 to zero
					rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

					do ixyz=1,nd
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					if (rij2 < rneighbr2) call linklist_checkpushneighbr(molnoi, molnoj)
				enddo
			enddo
			enddo
			enddo
			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine assign_to_neighbourlist

!----------------------------------------------------------------------------------
! Assign to Neighbourlist count each interaction only once
! Build a list connecting molecules which are within interacting distance of
! each other with each interaction counted only once

subroutine assign_to_neighbourlist_halfint
use module_linklist
implicit none

	integer                         :: i, j,k, ixyz   !Define dummy index
	integer				:: icell, jcell, kcell
	integer                         :: cellnp, adjacentcellnp
	integer				:: molnoi, molnoj
	integer, dimension(13)          :: icellshift, jcellshift, kcellshift
	double precision		:: rij2   !magnitude^2 between i and j
	double precision,dimension(3)   :: ri, rj !Position of molecule i and j
	double precision,dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldihead, oldi, currenti, oldjhead, oldj, currentj

	!Create Neighbourlist array based on current np
	allocate(neighbour%noneighbrs(np))
	allocate(neighbour%head(np))
	do i = 1,np
		neighbour%noneighbrs(i) = 0	!Zero number of molecules in neighbour list
		nullify(neighbour%head(i)%point)!Nullify neighbour list head pointer 
	enddo

	!Assign cell offsets
	icellshift = (/ 1, 1, 0,-1, 0, 1, 1, 0,-1,-1,-1, 0, 1/)
	jcellshift = (/ 0, 1, 1, 1, 0, 0, 1, 1, 1, 0,-1,-1,-1/)
	kcellshift = (/ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1/)

	!Calculate interactions between all cells in domain and halos
	do icell=2, ncells(1)+1
	do jcell=2, ncells(2)+1
	do kcell=2, ncells(3)+1

		!Retrieve cell np and set old to first molecule in list
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point 

		!Check interaction within own cell once
		do i = 1,cellnp                 !Step through each particle in list 
			molnoi = oldi%molno 	!Number of molecule
			ri = r(molnoi,:)        !Retrieve ri
			oldj => oldi%next	!Point j molecule to next molecule to i

			do j = i+1,cellnp          !Step through all j for each i

				molnoj = oldj%molno 	!Number of molecule
				rj = r(molnoj,:)        !Retrieve rj
				currentj => oldj
				oldj => currentj%next   !Use pointer in datatype to obtain next item in list

				if(molnoi==molnoj) cycle!Check to prevent interaction with self
	
				rij2=0                  !Set rij^2 to zero
				rij(:) = ri(:) - rj(:)  !Evaluate distance between particle i and j

				do ixyz=1,nd
					rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
				enddo

				if (rij2 < rneighbr2) call linklist_checkpushneighbr(molnoi, molnoj)
				!if (rij2 < rneighbr2) print*,'own_cell', molnoi, molnoj
			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo

		!Reset old to first molecule in list
		oldihead => cell%head(icell,jcell,kcell)%point

		do k = 1,13
			oldi => oldihead
			oldjhead => cell%head(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))%point
			adjacentcellnp = cell%cellnp(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))

			do i = 1,cellnp				!Step through each particle in list 
				molnoi = oldi%molno		!Number of molecule
				ri = r(molnoi,:)		!Retrieve ri

				oldj => oldjhead		!Reset j to head of linked list

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	!Number of molecule
					rj = r(molnoj,:)	!Retrieve rj
					currentj => oldj
					oldj => currentj%next	!Use pointer in datatype to obtain next item in list

					rij2=0			!Set rij^2 to zero
					rij(:) = ri(:) - rj(:)	!Evaluate distance between particle i and j

					do ixyz=1,nd
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					!if (rij2 < rneighbr2) print*,'neighbr_cells',  molnoi, molnoj
					if (rij2 < rneighbr2) call linklist_checkpushneighbr(molnoi, molnoj)

				enddo

				currenti => oldi
				oldi => currenti%next !Use pointer in datatype to obtain next item in list

			enddo


		enddo
	enddo
	enddo
	enddo

	!Build up interactions using unchecked halo cells

	!-----------------------------
	! Bottom xy plane domain face-
		   kcell=1
	!       9 Interactions	     -
	!-----------------------------

	!Perpendicular in z direction
	k = 5
	do icell=2, ncells(1)+1
	do jcell=2, ncells(2)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 6
	do icell=1, ncells(1)
	do jcell=2, ncells(2)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 10
	do icell=3,ncells(1)+2
	do jcell=2, ncells(2)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in yz plane
	k = 8
	do icell=2, ncells(1)+1
	do jcell=1, ncells(2)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 12
	do icell=2, ncells(1)+1
	do jcell=3, ncells(2)+2
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k =7
	do icell=1, ncells(1)
	do jcell=1, ncells(2)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 9
	do icell=3, ncells(1)+2
	do jcell=1, ncells(2)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 11
	do icell=3, ncells(1)+2
	do jcell=3, ncells(2)+2
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 13
	do icell=1, ncells(1)
	do jcell=3, ncells(2)+2
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!------------------------------
	! Bottom yz plane domain face -
		  icell = 1
	!	5 Interactions
	!------------------------------

	!Perpendicular in x direction
	k = 1
	do jcell = 2, ncells(2)+1
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 6
	do jcell = 2, ncells(2)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xy plane
	k = 2
	do jcell = 1, ncells(2)
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 7
	do jcell = 1, ncells(2)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 13
	do jcell = 3, ncells(2)+2
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!------------------------------
	! Top yz plane domain face    -
	    icell = ncells(1) + 2
	!	4 Interactions	      -
	!------------------------------

	!Edge cells diagonals in xy plane
	k = 4
	do jcell = 1, ncells(2)
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 10
	do jcell = 2, ncells(2)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 9
	do jcell = 1, ncells(2)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 11
	do jcell = 3, ncells(2)+2
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!------------------------------
	! Bottom xy plane domain face -
		  jcell = 1
	!	6 Interactions
	!------------------------------

	!Perpendicular in y direction
	k = 3
	do icell = 2, ncells(1)+1
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xy plane
	k = 2
	do icell = 2, ncells(1)
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 4
	do icell = 3,ncells(1)+1
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 8
	do icell = 2, ncells(1)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 7
	do icell = 2, ncells(1)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 9
	do icell = 3, ncells(1)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!------------------------------
	! Top xy plane domain face    -
	     jcell = ncells(2)+2
	!	3 Interactions
	!------------------------------

	!Edge cells diagonals in xz plane
	k = 12
	do icell = 2, ncells(1)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 11
	do icell = 3, ncells(1)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	k = 13
	do icell = 2, ncells(1)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(icell, jcell, kcell, k)
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine assign_to_neighbourlist_halfint


!------------------------------------------------------------------------------
!Routine to calculate all molecular interactions between two specified cells
!Used only for halo cell so works out interaction of j with i but adds to
!molecule i's neighbourlist

subroutine calculate_cell_interactions(icell, jcell, kcell, k)
use module_linklist
implicit none

	integer				:: i, j, ixyz
	integer				:: icell, jcell, kcell, k
	integer                         :: cellnp, adjacentcellnp
	integer				:: molnoi, molnoj
	integer, dimension(13)		:: icellshift, jcellshift, kcellshift
	double precision		:: rij2   !magnitude^2 between i and j
	double precision,dimension(3)   :: ri, rj !Position of molecule i and j
	double precision,dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldi, currenti, oldjhead, oldj, currentj

	!Assign cell offsets
	icellshift = (/ 1, 1, 0,-1, 0, 1, 1, 0,-1,-1,-1, 0, 1/)
	jcellshift = (/ 0, 1, 1, 1, 0, 0, 1, 1, 1, 0,-1,-1,-1/)
	kcellshift = (/ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1/)

		!Retrieve cell np and set old to first molecule in list
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list
		adjacentcellnp = cell%cellnp(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))
		oldjhead => cell%head(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))%point

		!Count number of interactions per cell
		!intcount(icell, jcell, kcell) = intcount(icell, jcell, kcell) + 1
		!intcount(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k)) = & 
		!intcount(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k)) + 1

		!write(100,'(6(i8,a))'), icell,';',jcell,';',kcell,';', &
		!icell+icellshift(k),';',jcell+jcellshift(k),';',kcell+kcellshift(k),';'

		!Check interaction with neighbouring cells
		do i = 1,cellnp                 !Step through each particle in list 
			molnoi = oldi%molno 	!Number of molecule
			ri = r(molnoi,:)        !Retrieve ri

			oldj => oldjhead	!Reset j to head of linked list

			do j = 1,adjacentcellnp          !Step through all j for each i

				molnoj = oldj%molno !Number of molecule
				rj = r(molnoj,:)         !Retrieve rj
				currentj => oldj
				oldj => currentj%next    !Use pointer in datatype to obtain next item in list

				rij2=0                   !Set rij^2 to zero
				rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

				do ixyz=1,nd
					rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
				enddo

				!if (rij2 < rneighbr2) print*,'neighbr_cells',  molnoi, molnoj
				!Used for halo cell so molnoj and molnoi swapped over!
				if (rij2 < rneighbr2) call linklist_checkpushneighbr(molnoj, molnoi)

			enddo

		currenti => oldi
		oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo

end subroutine calculate_cell_interactions

!======================================================================
!		Linklist manipulation Subroutines                     =
!======================================================================

!-----------------------------------------------------------------------
!Build linklist from array of positions of molecules

subroutine linklist_build(start, finish)
use module_linklist
implicit none

	integer            :: j
	integer            :: cellnp
	integer            :: icell, jcell, kcell	
	integer            :: start, finish
	type(node), pointer:: old, current

	print*, 'building linklist'
	allocate (old)            !Allocate storage space for pointer old to point at
	old%rp => r(start,:)          !Point to first molecule's r
	old%vp => v(start,:)          !Set up pointer to current molecule's v
	old%ap => a(start,:)          !Set up pointer to current molecule's a
	old%molno  = start       !Assign first molecule number 1
	nullify (old%next)        !Points to null as first in linked list
	
	do j=start+1,finish
		allocate (current)          !Allocate storage space for current to point at
		current%rp  => r(j,:)       !Set up pointer to current molecule's r
		current%vp  => v(j,:)       !Set up pointer to current molecule's v
		current%ap  => a(j,:)       !Set up pointer to current molecule's a
		current%molno  = j     !Assign each molecule a number
		old%previous => current     !Set up forward pointer to allow movement forward
		current%next => old         !Set current's link pointer to point to previous item
		old     => current          !Set current item to old ready for next loop
	enddo

	cellnp = finish - start+1
	cell%cellnp(icell,jcell,kcell) = cellnp
	nullify(current)      !Nullify current as no longer required
	nullify(old%previous) !Points top to null as last in linked list


end subroutine linklist_build

!===================================================================================
!Build circular linklist from array of positions of molecules

subroutine linklist_circbuild(start, finish)
use module_linklist
implicit none

	integer              :: j
	integer              :: cellnp
	integer              :: start, finish
	type(node), pointer  :: old, current
	type(node), pointer  :: bottom
        
	print*, 'building circular linklist'
	allocate (bottom)          !Allocate storage space for pointer bottom to point at
	bottom%rp => r(start,:)    !Point to first molecule's r
	bottom%vp => v(start,:)    !Set up pointer to current molecule's v
	bottom%ap => a(start,:)    !Set up pointer to current molecule's a
	bottom%molno  = start !Assign first molecule number

	allocate (old)            !Allocate storage space for pointer old to point at
	old%rp => r(start+1,:)    !Point to second molecule's r
	old%vp => v(start+1,:)    !Point to second molecule's v
	old%ap => a(start+1,:)    !Point to second molecule's a
	old%molno  = start+1 !Assign second molecule number
	bottom%previous => old    !Set up forward pointer to allow movement forward 
	old%next => bottom 	  !Set up first elements pointer to bottom of list

	do j=start+2,finish
		allocate (current)          !Allocate storage space for current to point at
		current%rp  => r(j,:)       !Set up pointer to current molecule's r
		current%vp  => v(j,:)       !Set up pointer to current molecule's v
		current%ap  => a(j,:)       !Set up pointer to current molecule's a
		current%molno   = j    !Assign each molecule a number
		old%previous => current     !Set up forward pointer to allow movement forward
		current%next => old         !Set current's link pointer to point to previous item
		old   => current            !Set current item to old ready for next loop
	enddo

	bottom%next => old  !Point next pointer of bottom to top of list - circular
	old%previous => bottom 	!Point previous pointer of top of list to bottom - circular
	cellnp = finish - start+1

	nullify(current) !Nullify current as no longer required
	
end subroutine linklist_circbuild

!===================================================================================
!Remove ('pop') a molecule from the stack an return its number, position and velocity

subroutine linklist_pop(icell, jcell, kcell, molnopop)
use module_linklist
implicit none

	integer            	       :: cellnp, flag
	integer, intent(in)    	       :: icell, jcell, kcell
	integer, intent(in)    	       :: molnopop
	type(node), pointer 	       :: old, current
	type(node), pointer 	       :: pop

	call linklist_gotomolecule(icell, jcell, kcell, molnopop,flag)
	old => cell%head(icell,jcell, kcell)%point

	if (associated(old) .eqv. .true.) then    !Exit if null
	if (molnopop == old%molno) then !Exit if not correct molno
		print*, 'pop', old%molno, 'from cell', icell, jcell, kcell

		cellnp = cell%cellnp(icell,jcell, kcell)

		pop => old
		old     => old%next !Set old to next item in list
		current => old%previous !Set current to previous item in list

		deallocate(pop)

		old%previous => current     !Previous pointer connects to old (pop missed out)
		current%next => old         !Next pointer connects to current (pop missed out)
	
		cellnp = cellnp - 1	          !Reduce number of molecules by one
		cell%cellnp(icell,jcell, kcell) = cellnp !Update cell molecular number
		print*, 'cell np', cellnp

		!Check molecules in cell, if none then cell pointer is nullified, otherwise
		!cell pointer is set to old to ensure it points to a molecule in the list
		if(cellnp==0) then
			nullify(cell%head(icell,jcell, kcell)%point) 
			print*, 'cell pointer set to null'
		else
			cell%head(icell,jcell, kcell)%point => old 
		endif
	endif
	endif

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required
	nullify(pop)  			    !Nullify pop as no longer required

end subroutine linklist_pop

!===================================================================================
!Adds molecule specified by passed variables moleculenpush to linked list

subroutine linklist_push(icell, jcell, kcell, molnopush)
use module_linklist
implicit none

	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	integer, target     :: molnopush
	type(node), pointer :: old, current
	type(node), pointer :: push

	!Obtain molecular number and top item of link list from cell head
	old => cell%head(icell,jcell, kcell)%point
	cellnp = cell%cellnp(icell,jcell, kcell)

	print*, 'push', molnopush

	allocate(push) !Allocate type to add to stack
	push%molno = molnopush !Build type from inputs
	push%rp => r(molnopush,:)   !Build type from inputs
	push%vp => v(molnopush,:)   !Build type from inputs
	push%ap => a(molnopush,:)   !Build type from inputs

	current => old           !Set current to first item in list
	old     => current%next  !Set old to next item in list
	
	old%previous => push     !Old previous pointer connects to new push
	current%next => push     !Current next pointer connects to new push
	push%next    => old      !Push next pointer connects to old
	push%previous=> current  !Push previous pointer connects to current
	old => push

	cell%head(icell,jcell,kcell)%point => old !Set cell pointer to top of cell list
	cellnp = cellnp + 1                 !Increase number of molecules by one
	cell%cellnp(icell,jcell,kcell) = cellnp   !Update cell molecular number

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required
	nullify(push)                       !Nullify old as no longer required

end subroutine linklist_push


!===================================================================================
!Adds molecule specified by passed variables to linked list with check to see if
!linklist is empty so new linklist can be established

subroutine linklist_checkpush(icell, jcell, kcell, molnopush)
use module_linklist
implicit none

	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	integer	            :: molnopush
	type(node), pointer :: old, current

	cellnp = cell%cellnp(icell,jcell, kcell)
	allocate(current) 		!Allocate type to add to stack
	current%molno = molnopush 	!Build type from inputs
	current%rp => r(molnopush,:)   	!Build type from inputs
	current%vp => v(molnopush,:)   	!Build type from inputs
	current%ap => a(molnopush,:)  	!Build type from inputs

	select case (cellnp)
	case(0)
		!print*, 'empty cell - cell contains', cellnp, 'molecules'
		nullify(current%previous) !Nullify pointer at top of list
		nullify(current%next)     !Nullify pointer at bottom of list
	case(1:)
		old => cell%head(icell,jcell,kcell)%point !Set old to top of list
		current%next => old                	!Set old to next item in list
		old%previous => current            	!Old previous pointer connects to new current
		nullify(current%previous)          	!Nullify pointer at top of list
	end select

	cell%head(icell,jcell,kcell)%point => current 	!Set cell pointer to top of cell list
	cellnp = cellnp + 1                     	!Increase number of particles by one
	cell%cellnp(icell,jcell,kcell) = cellnp       	!Update cell molecular number

	!Nullify pointers current and old so cell head is left pointing to top of the
	!required linked list
	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_checkpush

!===================================================================================
!Adds molecule specified by passed variables to linked list with check to see if
!linklist is empty so new linklist can be established

subroutine linklist_checkpush_bin(ibin, jbin, kbin, molnopush)
use module_linklist
implicit none

	integer             :: binnp		
	integer             :: ibin, jbin, kbin
	integer		    :: molnopush
	type(node), pointer :: old, current

	binnp = bin%cellnp(ibin,jbin,kbin)
	allocate(current) 		!Allocate type to add to stack
	current%molno = molnopush 	!Build type from inputs
	nullify(current%rp)    		!Build type from inputs
	nullify(current%vp)    		!Build type from inputs
	nullify(current%ap) 	  	!Build type from inputs

	select case (binnp)
	case(0)
		!print*, 'empty bin - bin contains', binnp, 'molecules'
		nullify(current%previous) !Nullify pointer at top of list
		nullify(current%next)     !Nullify pointer at bottom of list
	case(1:)
		old => bin%head(ibin,jbin,kbin)%point !Set old to top of list
		current%next => old                	!Set old to next item in list
		old%previous => current            	!Old previous pointer connects to new current
		nullify(current%previous)          	!Nullify pointer at top of list
	end select

	bin%head(ibin,jbin,kbin)%point => current 	!Set bin pointer to top of bin list
	binnp = binnp + 1                     		!Increase number of particles by one
	bin%cellnp(ibin,jbin,kbin) = binnp       	!Update bin molecular number

	!Nullify pointers current and old so bin head is left pointing to top of the
	!required linked list
	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_checkpush_bin

!===================================================================================
!Adds molecule specified by passed variables to neighbour linked list with check if
!linklist is empty so new linklist can be established

subroutine linklist_checkpushneighbr(molnoi, molnoj)
use module_linklist
implicit none

	integer			  	:: noneighbrs
	integer                    	:: molnoi, molnoj
	type(neighbrnode), pointer 	:: old, current

	noneighbrs = neighbour%noneighbrs(molnoi)
	allocate(current)                         !Allocate type to add to stack
	current%molnoj = molnoj                   !Build type from inputs

	select case (noneighbrs)
	case(0)
		nullify(current%previous)         !Nullify pointer at top of list
		nullify(current%next)             !Nullify pointer at bottom of list
	case(1:)
		old => neighbour%head(molnoi)%point!Set old to top of list
		current%next => old			!Set old to next item in list
		old%previous => current			!Old previous pointer connects to new current
		nullify(current%previous)		!Nullify pointer at top of list
	end select

	neighbour%head(molnoi)%point => current	!Set cell pointer to top of cell list
	noneighbrs = noneighbrs + 1               	!Increase number of particles by one
	neighbour%noneighbrs(molnoi) = noneighbrs	!Update neighbour list molecular number

	!Nullify pointers current and old so neighbour head is left pointing to top
	nullify(current)                          !Nullify current as no longer required
	nullify(old)                              !Nullify old as no longer required

end subroutine linklist_checkpushneighbr

!===================================================================================
!Adds molecule specified by passed variables to passed molecule linked list with check if
!linklist is empty so new linklist can be established

subroutine linklist_checkpushmol(molno,ipass,jpass,kpass)
use module_linklist
implicit none

	integer,intent(in)         :: molno
	integer			   :: sendnp
	integer			   :: ipass, jpass, kpass
	type(passnode), pointer    :: old, current

	sendnp = pass%sendnp
	allocate(current)                         !Allocate type to add to stack
	current%molno = molno                     !Build type from inputs
	current%ipass = ipass			  !Build type from inputs
	current%jpass = jpass			  !Build type from inputs
	current%kpass = kpass			  !Build type from inputs

	select case (sendnp)
	case(0)
		nullify(current%previous)         !Nullify pointer at top of list
		nullify(current%next)             !Nullify pointer at bottom of list
	case(1:)
		old => pass%head                  !Set old to top of list
		current%next => old               !Set old to next item in list
		old%previous => current           !Old previous pointer connects to new current
		nullify(current%previous)         !Nullify pointer at top of list
	end select

	pass%head => current                      !Set cell pointer to top of cell list
	sendnp = sendnp + 1                   	  !Increase number of particles by one
	pass%sendnp = sendnp                  	  !Update neighbour list molecular number

	!Nullify pointers current and old so neighbour head is left pointing to top
	nullify(current)                          !Nullify current as no longer required
	nullify(old)                              !Nullify old as no longer required

end subroutine linklist_checkpushmol

!===================================================================================
!Move backwards through linked list and print out results
	
subroutine linklist_print(icell, jcell, kcell)
use module_linklist
implicit none

	integer             :: j
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current

	!Obtain molecular number and top item of link list from cell head
	old => cell%head(icell,jcell, kcell)%point
	cellnp = cell%cellnp(icell,jcell, kcell)

	print*, 'linklist print called for cell', icell, jcell, kcell
	if(cellnp == 0) print*, 'linklist empty'

	current => old ! make current point to head of list
	do j=1,cellnp
		!print*, 'more items in linked list?: ', associated(old%next)
		print'(i8,3f10.5)', current%molno, current%rp
		if (associated(old%next) .eqv. .true. ) then !Exit if null
			old => current%next ! Use pointer in datatype to obtain next item in list
			current => old      ! make current point to old - move alone one
		endif
	enddo

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_print

!===================================================================================
!Move through all cells linked list and print out results reseting back to the top
!upon completion 
	
subroutine linklist_printallcells
use module_linklist
implicit none

	integer             :: j
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current

	print*, 'linklist print called for all cells'

	do icell=1,ncells(1)+2
	do jcell=1,ncells(2)+2
	do kcell=1,ncells(3)+2

		!Obtain molecular number and top item of link list from cell head
		old => cell%head(icell,jcell, kcell)%point
		cellnp = cell%cellnp(icell,jcell, kcell)

		print*, 'Cell', icell, jcell, kcell
		if(cellnp == 0) print*, 'linklist empty'

		current => old ! make current point to head of list
		do j=1,cellnp
			!print*, 'more items in linked list?: ', associated(old%next)
			print*, current%molno, current%rp
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next ! Use pointer in datatype to obtain next item in list
				current => old      ! make current point to old - move alone one
			endif
		enddo

	enddo
	enddo
	enddo

	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_printallcells

!===================================================================================
!Move through all cells linked list not including the halo cells and print out 
!results reseting back to the top upon completion 
	
subroutine linklist_printalldomaincells
use module_linklist
implicit none

	integer             :: j
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current

	print*, 'linklist print called for all domain cells'

	do icell=2,ncells(1)+1
	do jcell=2,ncells(2)+1
	do kcell=2,ncells(3)+1

		!Obtain molecular number and top item of link list from cell
		old => cell%head(icell,jcell, kcell)%point
		cellnp = cell%cellnp(icell,jcell, kcell)

		print*, 'Cell', icell, jcell, kcell
		if(cellnp == 0) print*, 'linklist empty'

		current => old ! make current point to head of list
		do j=1,cellnp
			!print*, 'more items in linked list?: ', associated(old%next)
			print*, current%molno, current%rp
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next ! Use pointer in datatype to obtain next item in list
				current => old      ! make current point to old - move alone one
			endif
		enddo

	enddo
	enddo
	enddo
	
	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_printalldomaincells


!===================================================================================
!Move backwards through linked list and print out results
	
subroutine linklist_printneighbourlist
use module_linklist
implicit none

	integer                    :: i, j
	integer			   :: noneighbrs
	type(neighbrnode), pointer :: old, current

	do i = 1,np

		!Obtain molecular number and top item of link list from cell head
 		noneighbrs = neighbour%noneighbrs(i)  	!Determine number of elements in neighbourlist
		old => neighbour%head(i)%point		!Set old to head of neighbour list

		if(noneighbrs == 0) print*, i, 'has 0 neighbours'

		current => old ! make current point to head of list
		do j=1,noneighbrs
			!print*, 'more items in linked list?: ', associated(old%next)
			print'(2(a,i8),6f10.5)', 'i = ', i,' j = ', current%molnoj, r(i,:), r(current%molnoj,:)
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next ! Use pointer in datatype to obtain next item in list
				current => old      ! make current point to old - move alone one
			endif
		enddo

	enddo

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_printneighbourlist

!===================================================================================
!Move backwards through linked list and print out results
	
subroutine linklist_printpassedlist
use module_linklist
implicit none

	integer                    :: j
	integer			   :: sendnp
	type(passnode), pointer    :: old, current

	!Obtain molecular number and top item of link list from cell head
        sendnp = pass%sendnp	!Determine number of elements in neighbourlist
	old => pass%head	!Set old to head of neighbour list

	print*, 'linklist print called for pass list with', sendnp, 'elements'
	if(sendnp == 0) print*, 'linklist empty'

	current => old ! make current point to head of list
	do j=1,sendnp
		!print*, 'more items in linked list?: ', associated(old%next)
		print*, current%molno, current%ipass, current%jpass, current%kpass
		if (associated(old%next) .eqv. .true. ) then !Exit if null
			old => current%next ! Use pointer in datatype to obtain next item in list
			current => old      ! make current point to old - move alone one
		endif
	enddo

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_printpassedlist


!===================================================================================
!Print array and linklist to compare value 
	
subroutine linklist_compareprint(icell, jcell, kcell)
use module_linklist
implicit none

	integer            :: j         !Dummy counter
	integer		   :: n         !Molecule number
	integer            :: cellnp		
	integer            :: icell, jcell, kcell
	type(node), pointer:: old, current

	!Obtain molecular number and top item of link list from cell
	old => cell%head(icell,jcell, kcell)%point
	cellnp = cell%cellnp(icell,jcell, kcell)

	!call linklist_gototop(icell, jcell, kcell) 
	!print*, 'Print array and linklist for comparison has been called - Should be identical'
	current => old
	do j=1,cellnp
		n = old%molno
		print'(i0,a,2f20.15,a,2f10.5,a,2f10.5)',old%molno, ' ap = ',old%ap, &
			 ' vp = ', old%vp, ' rp = ', old%rp
		print'(i0,a,2f20.15,a,2f10.5,a,2f10.5)',n, ' a  = ',a(n,:), ' v  = ', &
			 v(n,:), ' r  = ', r(n,:)
		old => current%next !Use pointer in datatype to obtain next item in list
		current => old          !make current point to old - move along one
	enddo

	!cell%head(icell,jcell, kcell)%point => old

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_compareprint

!===================================================================================
!Deallocate cell linklist

subroutine linklist_deallocate(icell, jcell, kcell)
use module_linklist
implicit none

	integer            :: j
	integer            :: cellnp
	integer            :: icell, jcell, kcell
	type(node), pointer:: old, current		

	if (associated(cell%head(icell,jcell, kcell)%point) .eqv. .true. ) then !Exit if null

		!Obtain molecular number and top item of link list from cell head
		old => cell%head(icell,jcell, kcell)%point
		cellnp = cell%cellnp(icell,jcell, kcell)

		current => old ! make current point to head of list
		do j=1,cellnp-1
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next       !make list point to next node of old
				nullify(current%next)     !Remove pointer to next
				nullify(current%previous) !Remove pointer to previous
				deallocate(current)       !Deallocate current entry
				current => old            !Make current point to new previous
			endif
		enddo
		nullify(old%next)     !Remove pointer to next
		nullify(old%previous) !Remove pointer to previous
		deallocate(old)       !Deallocate final entry
		nullify(cell%head(icell,jcell, kcell)%point) !Set cell head pointer to null
		cell%cellnp(icell,jcell, kcell) = 0          !Zero cell molecule number
	endif

end subroutine linklist_deallocate

!===================================================================================
!Deallocate bin linklist

subroutine linklist_deallocate_bins
use module_linklist
use calculated_properties_MD
implicit none

	integer            	:: j
	integer            	:: binnp
	integer           	:: ibin, jbin, kbin
	type(node), pointer	:: old, current

	do ibin = 1,1
	do jbin = 1,nbins(1)
	do kbin = 1,1	

		if (associated(bin%head(ibin,jbin,kbin)%point) .eqv. .true. ) then !Exit if null

			!Obtain molecular number and top item of link list from bin head
			old => bin%head(ibin,jbin,kbin)%point
			binnp = bin%cellnp(ibin,jbin,kbin)

			current => old ! make current point to head of list
			do j=1,binnp-1
				if (associated(old%next) .eqv. .true. ) then !Exit if null
					old => current%next       !make list point to next node of old
					nullify(current%next)     !Remove pointer to next
					nullify(current%previous) !Remove pointer to previous
					deallocate(current)       !Deallocate current entry
					current => old            !Make current point to new previous
				endif
			enddo
			nullify(old%next)     !Remove pointer to next
			nullify(old%previous) !Remove pointer to previous
			deallocate(old)       !Deallocate final entry
			nullify(bin%head(ibin,jbin,kbin)%point)  !Set bin head pointer to null
			bin%cellnp(ibin,jbin, kbin) = 0          !Zero bin molecule number
		endif

	enddo
	enddo
	enddo

end subroutine linklist_deallocate_bins


!===================================================================================
!Deallocate all linklists

subroutine linklist_deallocatepasslist
use module_linklist
implicit none

	integer            :: j
	integer            :: sendnp
	type(passnode), pointer :: old, current

	if (associated(pass%head) .eqv. .true. ) then !Exit if null
        	sendnp = pass%sendnp  		!Determine number of elements in passlist
		old => pass%head		!Set old to head of pass list
		current => old 			!make current point to head of list
		do j=1,sendnp-1
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next       !Make list point to next node of old
				nullify(current%next)     !Remove pointer to next
				nullify(current%previous) !Remove pointer to previous
				deallocate(current)       !Deallocate current entry
				current => old            !Make current point to new previous
			endif
		enddo

		nullify(old%next)	!Remove pointer to next
		nullify(old%previous)	!Remove pointer to previous
		deallocate(old)         !Deallocate final entry
		nullify(pass%head)	!Set pass head pointer to null
		pass%sendnp = 0		!Zero cell molecule number
	endif

end subroutine linklist_deallocatepasslist

!===================================================================================
!Deallocate all linklists

subroutine linklist_deallocateall
use module_linklist
implicit none

	integer            :: i, j
	integer            :: cellnp, noneighbrs
	integer            :: icell, jcell, kcell
	type(node), pointer:: old, current
	type(neighbrnode), pointer :: oldn, currentn

	do icell=1,ncells(1)+2
	do jcell=1,ncells(2)+2
	do kcell=1,ncells(3)+2

		if (associated(cell%head(icell,jcell,kcell)%point) .eqv. .true. ) then !Exit if null
			old => cell%head(icell,jcell,kcell)%point
			cellnp = cell%cellnp(icell,jcell, kcell)
			current => old ! make current point to head of list
			do j=1,cellnp-1
				if (associated(old%next) .eqv. .true. ) then !Exit if null
					old => current%next       !Make list point to next node of old
					nullify(current%next)     !Remove pointer to next
					nullify(current%previous) !Remove pointer to previous
					deallocate(current)       !Deallocate current entry
					current => old            !Make current point to new previous
				endif
			enddo

			nullify(old%next)     !Remove pointer to next
			nullify(old%previous) !Remove pointer to previous
			deallocate(old)       !Deallocate final entry
			nullify(cell%head(icell,jcell,kcell)%point) !Set cell head pointer to null
			cell%cellnp(icell,jcell, kcell) = 0         !Zero cell molecule number
		endif
	enddo
	enddo
	enddo

	!do icell=1,ncells(1)+2
	!do jcell=1,ncells(2)+2
	!do kcell=1,ncells(3)+2
	!	nullify(cell%head(icell,jcell,kcell)%point)  !Set cell head pointer to null
	!	cell%cellnp(icell,jcell,kcell) = 0           !Zero cell molecule number
	!enddo
	!enddo
	!enddo

	do i = 1, np
		if (associated(neighbour%head(i)%point) .eqv. .true. ) then !Exit if null
        		noneighbrs = neighbour%noneighbrs(i)  !Determine number of elements in neighbourlist
			oldn => neighbour%head(i)%point		   !Set old to head of neighbour list
			currentn => oldn ! make current point to head of list
			do j=1,noneighbrs-1
				if (associated(oldn%next) .eqv. .true. ) then !Exit if null
					oldn => currentn%next       !Make list point to next node of old
					nullify(currentn%next)      !Remove pointer to next
					nullify(currentn%previous)  !Remove pointer to previous
					deallocate(currentn)        !Deallocate current entry
					currentn => oldn            !Make current point to new previous
				endif
			enddo

			nullify(oldn%next)         !Remove pointer to next
			nullify(oldn%previous)     !Remove pointer to previous
			deallocate(oldn)           !Deallocate final entry
			nullify(neighbour%head(i)%point)   !Set neighbour head pointer to null
			neighbour%noneighbrs(i) = 0  !Zero cell molecule number
		endif
	enddo

	!Deallocate array of molecules neighbourlist pointers
	deallocate(neighbour%noneighbrs)
	deallocate(neighbour%head)

end subroutine linklist_deallocateall

!=============================================================================
!Moves to the bottom of the linklist

subroutine linklist_gotobottom(icell, jcell, kcell)
use module_linklist
implicit none

	integer            :: j
	integer            :: cellnp		
	integer            :: icell, jcell, kcell
	type(node), pointer:: old, current

	old => cell%head(icell,jcell, kcell)%point
	cellnp = cell%cellnp(icell,jcell, kcell)

	print*, 'linklist go to bottom called'

	current => old ! make current point to head of list
	do j=1,cellnp
		print*, 'more items in linked list?: ', associated(old%next)
		print*, current%molno, current%rp
		if (associated(old%next) .eqv. .true. ) then !Exit if null
			old => current%next ! make list point to next node of old
			current => old      ! make current point to new previous
		endif
	enddo

	cell%head(icell,jcell, kcell)%point => old

	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_gotobottom

!============================================================================
!Moves to the top of the linklist

subroutine linklist_gototop(icell, jcell, kcell)
use module_linklist
implicit none

	integer            :: j
	integer            :: cellnp		
	integer            :: icell, jcell, kcell
	type(node), pointer:: old, current

	old => cell%head(icell,jcell, kcell)%point
	cellnp = cell%cellnp(icell,jcell, kcell)

	print*, 'linklist go to top called'

	current => old ! make current point to head of list
	do j=1,cellnp
		!print*, 'more items in linked list?: ', associated(old%previous)
		!print*, current%molno, current%rp
		if (associated(old%previous) .eqv. .true. ) then !Exit if null
			old => current%previous     !Set up forward pointer to allow movement forward
			current  => old             !Set current item to old ready for next loop
		endif
	enddo

	cell%head(icell,jcell, kcell)%point => old

	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_gototop	

!============================================================================
!Finds a specified molecule and sets the cellpointer to that molecule
!If molecule is not found, flag is set to zero

subroutine linklist_gotomolecule(icell, jcell, kcell, n, flag)
use module_linklist
implicit none

	integer             :: j
	integer, intent(in) :: n
	integer, intent(out):: flag
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current

	old => cell%head(icell,jcell, kcell)%point
	cellnp = cell%cellnp(icell,jcell, kcell)

	!print*, 'linklist go to molecule', n, '  called'

	current => old ! make current point to head of list
	do j=1,cellnp
		old => current%previous  !Set up forward pointer to allow movement forward
		current  => old          !Set current item to old ready for next loop
		if (old%molno .eq. n ) then !Exit if null
			print*, 'molecule', n, '  found'
			cell%head(icell,jcell, kcell)%point => old
			flag = 1 !Set flag to 1 to indicate molecule found
			return
		endif
	enddo

	flag = 0 !Set flag to zero to indicate molecule not in list
	!print*, '**ERROR - molecule', n, 'is not in this cell**'
	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_gotomolecule


!==========================================================================
!Move a specified number of places upwards or downwards through the linklist
!determine by input n and direction dir which should be either 'goup' or
!down

subroutine linklist_movethrough(icell, jcell, kcell, n, dir)
use module_linklist
implicit none

	integer            :: n, j
	character(len=4)   :: dir
	integer            :: icell, jcell, kcell
	type(node), pointer:: old, current

	old => cell%head(icell,jcell, kcell)%point

	print*, 'linklist move through called'
	current => old ! make current point to head of list

	if (dir == 'goup') then
		do j=1,n 
			!print*, 'more items in linked list?: ', associated(old%previous)
			!print*, current%molno, current%rp
			if (associated(old%previous) .eqv. .true. ) then !Exit if null
				old => current%previous  !make old point to previous node of current
				current  => old          !Set current item to old ready for next loop
			endif
		enddo
	elseif (dir == 'down') then
		do j=1,n
			!print*, 'more items in linked list?: ', associated(old%next)
			!print*, current%molno, current%rp
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next      !make old point to next node of current
				current => old           !Set current item to old ready for next loop
			endif
		enddo
	else
	print*, 'error - direction specified should be of type character and either *goup* or *down*'
	endif

	cell%head(icell,jcell, kcell)%point => old !Assign cell pointer to top of list

	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_movethrough	


!===================================================================================
!Prints the details of the current element of the linked list as well as those either
!side in the linked list

subroutine linklist_printcurrent(icell, jcell, kcell)
use module_linklist
implicit none

	integer            :: icell, jcell, kcell
	type(node), pointer :: old

	old => cell%head(icell,jcell, kcell)%point
	if (associated(old) .eqv. .true.) then !exit if empty

		if (associated(old%next) .eqv. .true. ) then
			print'(a,i8,3f10.5)', 'element below =', old%next%molno, old%next%rp
		else
			print*, 'no element below'
		endif

		print'(a,i8,3f10.5)', 'current element =', old%molno, old%rp

		if (associated(old%previous) .eqv. .true. ) then
			print'(a,i8,3f10.5)', 'element above = ', old%previous%molno, old%previous%rp
		else
			print*, 'no element above'
		endif
	endif

end subroutine linklist_printcurrent
