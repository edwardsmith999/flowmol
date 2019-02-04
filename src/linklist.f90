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

module linked_list

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD

	integer, dimension(6,6,6)	:: intcount

	type listinfo
        !Dummy list type
    end type listinfo	

	!3D Information for each cell and pointer to cell/bin linked list
	type cellinfo
		integer  , dimension(:,:,:), allocatable :: cellnp
		type(ptr), dimension(:,:,:), pointer     :: head
	end type cellinfo
	
	!Cannot create array of pointers in Fortran so create an array of
	!type ptr which contains a pointer
	type ptr
		type(node), pointer :: point
	end type ptr

	!A Node of the linked list which is built for each cell/bin
	type node
		integer             :: molno
		type(node), pointer :: next, previous
	end type node

	! Neighbourlist with cell head 
	!Information for pointer to neighbour list
	type neighbrinfo
		integer,dimension(:),allocatable    :: Nlist
		type(ptr), dimension(:), pointer	:: head
	end type neighbrinfo
	
	!Information for pass list
	type passinfo
		integer					:: sendnp
		type(passnode), pointer :: head
	end type passinfo

	!A Node of the passed molecule list
	type passnode
		integer                 :: molno
		integer					:: ipass
		integer					:: jpass
		integer					:: kpass
		type(passnode), pointer :: next, previous
	end type passnode


	! Neighbourlist with cell head 
	!Information for pointer to neighbour list
	type, extends(neighbrinfo) :: clusterinfo
        integer :: Nclust
        integer :: maxclusts=10000
        integer, dimension(:),allocatable :: inclust, clusterngbrs
	end type clusterinfo
	

	!Global type cell used to record location of all cells linked lists
	!N.B. Type cell is not a pointer & nothing points to it
	type(cellinfo)    	:: cell
	type(cellinfo)    	:: bin	!Used to bin molecules for statistics
	type(neighbrinfo) 	:: neighbour
	type(clusterinfo)   :: cluster
	type(passinfo)    	:: pass


contains


	!============================================================================
	!Finds a specified molecule and sets the cellpointer to that molecule
	!If molecule is not found, flag is set to zero

	subroutine linklist_gotomolecule(self, icell, jcell, kcell, n, mol, success)
		implicit none

        type(cellinfo),intent(inout)    :: self
		integer, intent(in) 			:: n
		type(node),pointer,intent(out)  :: mol
        logical, intent(out)            :: success

		integer             			:: j
		integer             			:: cellnp		
		integer             			:: icell, jcell, kcell
		type(node), pointer 			:: old, current

        !Initialise success to be false
        success = .false.

        !If cells are out of range, return null and fail
        if (icell .gt. ncells(1)+2 .or. &
            icell .lt. 1           .or. &
            jcell .gt. ncells(2)+2 .or. &
            jcell .lt. 1           .or. &
            kcell .gt. ncells(3)+2 .or. &
            kcell .lt. 1                 ) then
            
            nullify(mol)
            return

        end if 

		!If celllist is empty, return null and fail
		cellnp = self%cellnp(icell,jcell, kcell)
		if (cellnp .eq. 0) then
			nullify(mol)	 !Nullify to indicate molecule not in list
			return
		endif

		old => self%head(icell,jcell, kcell)%point

		!print'(a,i8,a,3i6)', 'linklist go to molecule', n, '  called in cell ', icell,jcell,kcell

		current => old ! make current point to head of list
		do j=1,cellnp
			!print*, 'mol ',j, ' of ', cellnp, 'number', old%molno
			if (old%molno .eq. n ) then 
				!print*, 'molecule', n, '  found'
				mol => old
				!self%head(icell,jcell, kcell)%point => old
                success = .true.
				return
			endif
			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo
				
		nullify(mol)	 !Nullify to indicate molecule not in list
		nullify(current) !Nullify current as no longer required
		nullify(old)     !Nullify old as no longer required

	end subroutine linklist_gotomolecule



!----------------------------------------------------------------------------------
! Assign to Neighbourlist each molecules only once so each interaction is
! counted twice.

subroutine assign_to_neighbourlist_allint(self, cellin)
	implicit none

    type(neighbrinfo),intent(out) :: self
    type(cellinfo),intent(in)       :: cellin

    integer                         :: i, j !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	real(kind(0.d0))				:: rij2   !magnitude^2 between i and j
	real(kind(0.d0)),dimension(3)   :: ri, rj !Position of molecule i and j
	real(kind(0.d0)),dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	!Create Neighbourlist array based on current np
	allocate(self%Nlist(np))
	allocate(self%head(np))
	do i = 1,np
		self%Nlist(i) = 0	!Zero number of molecules in neighbour list
		nullify(self%head(i)%point)!Nullify neighbour list head pointer 
	enddo

	!if (maxval(cellin%cellnp(:,:,:)) .gt. 9) call error_abort("ERROR - greater than 10 per cell")

	do icell=2, ncells(1)+1
	do jcell=2, ncells(2)+1
	do kcell=2, ncells(3)+1

		cellnp = cellin%cellnp(icell,jcell,kcell)
		oldi => cellin%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp                  !Step through each particle in list 
			molnoi = oldi%molno 	 !Number of molecule
			ri = r(:,molnoi)         !Retrieve ri

			do icellshift = -1,1
			do jcellshift = -1,1
			do kcellshift = -1,1
				oldj => cellin%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cellin%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				!print*, icell+icellshift,jcell+jcellshift,kcell+kcellshift

				do j = 1,adjacentcellnp         !Step through all j for each i

					molnoj = oldj%molno 	 	!Number of molecule
					rj = r(:,molnoj)         	!Retrieve rj

					currentj => oldj
					oldj => currentj%next    	!Use pointer in datatype to obtain next item in list
					
					if(molnoi .eq. molnoj) cycle 	!Check to prevent interaction with self

					rij(:) = ri(:) - rj(:)   	!Evaluate distance between particle i and j
					rij2 = dot_product(rij,rij) !Square of rij
					
					if (potential_flag.eq.1) call check_update_adjacentbeadinfo_allint(molnoi,molnoj)	

					if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoi, molnoj)


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

	nullify(oldi)      		!Nullify as no longer required
	nullify(oldj)      		!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine assign_to_neighbourlist_allint

!----------------------------------------------------------------------------------
! Assign to Neighbourlist count each interaction only once
! Build a list connecting molecules which are within interacting distance of
! each other with each interaction counted only once

subroutine assign_to_neighbourlist_halfint(self, cellin)
	implicit none

    type(neighbrinfo),intent(inout) :: self
    type(cellinfo),intent(in)       :: cellin

	integer                         :: i, j,k !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	integer, dimension(13)          :: icellshift, jcellshift, kcellshift
	real(kind(0.d0))				:: rij2   !magnitude^2 between i and j
	real(kind(0.d0)),dimension(3)   :: ri, rj !Position of molecule i and j
	real(kind(0.d0)),dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldihead, oldi, currenti, oldjhead, oldj, currentj

	!Create Neighbourlist array based on current np
	allocate(self%Nlist(np))
	allocate(self%head(np))
	do i = 1,np
		self%Nlist(i) = 0	!Zero number of molecules in neighbour list
		nullify(self%head(i)%point)!Nullify neighbour list head pointer 
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
		cellnp = cellin%cellnp(icell,jcell,kcell)
		oldi => cellin%head(icell,jcell,kcell)%point 

		!Check interaction within own cell once
		do i = 1,cellnp                 !Step through each particle in list 
			molnoi = oldi%molno 	!Number of molecule
			ri = r(:,molnoi)        !Retrieve ri
			oldj => oldi%next	!Point j molecule to next molecule to i

			do j = i+1,cellnp          !Step through all j for each i

				molnoj = oldj%molno 	!Number of molecule
				rj = r(:,molnoj)        !Retrieve rj
				currentj => oldj
				oldj => currentj%next   !Use pointer in datatype to obtain next item in list

				if(molnoi==molnoj) cycle!Check to prevent interaction with self
	
				rij(:) = ri(:) - rj(:)  !Evaluate distance between particle i and j

				rij2 = dot_product(rij,rij)	!Square of vector calculated

				if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
				if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoi, molnoj)
				!if (rij2 < rneighbr2) print*,'own_cell', molnoi, molnoj
			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo

		!Reset old to first molecule in list
		oldihead => cellin%head(icell,jcell,kcell)%point

		do k = 1,13
			oldi => oldihead
			oldjhead => cellin%head(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))%point
			adjacentcellnp = cellin%cellnp(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))

			do i = 1,cellnp				!Step through each particle in list 
				molnoi = oldi%molno		!Number of molecule
				ri = r(:,molnoi)		!Retrieve ri

				oldj => oldjhead		!Reset j to head of linked list

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	!Number of molecule
					rj = r(:,molnoj)	!Retrieve rj
					currentj => oldj
					oldj => currentj%next	!Use pointer in datatype to obtain next item in list

					rij(:) = ri(:) - rj(:)	!Evaluate distance between particle i and j
					rij2 = dot_product(rij,rij)	!Square of vector calculated

					if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
					!if (rij2 < rneighbr2) print*,'neighbr_cells',  molnoi, molnoj
					if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoi, molnoj)

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
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 6
	do icell=1, ncells(1)
	do jcell=2, ncells(2)+1
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 10
	do icell=3,ncells(1)+2
	do jcell=2, ncells(2)+1
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in yz plane
	k = 8
	do icell=2, ncells(1)+1
	do jcell=1, ncells(2)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 12
	do icell=2, ncells(1)+1
	do jcell=3, ncells(2)+2
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k =7
	do icell=1, ncells(1)
	do jcell=1, ncells(2)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 9
	do icell=3, ncells(1)+2
	do jcell=1, ncells(2)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 11
	do icell=3, ncells(1)+2
	do jcell=3, ncells(2)+2
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 13
	do icell=1, ncells(1)
	do jcell=3, ncells(2)+2
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
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
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 6
	do jcell = 2, ncells(2)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xy plane
	k = 2
	do jcell = 1, ncells(2)
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 7
	do jcell = 1, ncells(2)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 13
	do jcell = 3, ncells(2)+2
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
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
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 10
	do jcell = 2, ncells(2)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 9
	do jcell = 1, ncells(2)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 11
	do jcell = 3, ncells(2)+2
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
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
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xy plane
	k = 2
	do icell = 2, ncells(1)
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 4
	do icell = 3,ncells(1)+1
	do kcell = 2, ncells(3)+1
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Edge cells diagonals in xz plane
	k = 8
	do icell = 2, ncells(1)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 7
	do icell = 2, ncells(1)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 9
	do icell = 3, ncells(1)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
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
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	!Corner cell xyz diagonals
	k = 11
	do icell = 3, ncells(1)+1
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	k = 13
	do icell = 2, ncells(1)
	do kcell = 2, ncells(3)
		call calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

contains

    subroutine all_cell_interactions

        integer               :: icell,jcell,kcell,x,y,z,x1,y1,z1,x2,y2,z2
        integer,dimension(27) :: imin,imax,jmin,jmax,kmin,kmax,shift

        x = ncells(1); x1 = x+1; x2 = x+2
        y = ncells(2); y1 = y+1; y2 = y+2
        z = ncells(3); z1 = z+1; z2 = z+2
        imin = [2,1,3,2,2,1,3,3,1,1,1,1,1,1,x2,x2,x2,x2,2,2,3,2,2,3,2,3,2]
        imax = [x1,x,x2,x1,x1,x,x2,x2,x,1,1,1,1,1,x2,x2,x2,x2,x1,x,x1,x1,x,x1,x1,x1,x]
        jmin = [2,2,2,1,3,1,1,3,3,2,2,1,1,3,1,2,1,3,1,1,1,1,1,1,y2,y2,y2]
        jmax = [y1,y1,y1,y,y2,y,y,y2,y2,y1,y1,y,y,y2,y,y1,y,y2,1,1,1,1,1,1,y2,y2,y2]
        kmin = [1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
        kmax = [1,1,1,1,1,1,1,1,1,z1,z,z1,z,z,z1,z,z,z,z1,z1,z1,z,z,z,z,z,z]
        shift= [5,6,10,8,12,7,9,11,13,1,6,2,7,13,4,10,9,11,3,2,4,8,7,9,12,11,13]

        do i = 1,size(shift)
            do icell = imin(i),imax(i)
            do jcell = jmin(i),jmax(i)
            do kcell = kmin(i),kmax(i)
                call calculate_cell_interactions(cellin, icell, jcell, kcell, shift(i))
            enddo
            enddo
            enddo
        enddo

    end subroutine all_cell_interactions

end subroutine assign_to_neighbourlist_halfint




!----------------------------------------------------------------------------------
! Assign to Neighbourlist count each interaction only once
! Build a list connecting molecules which are within interacting distance of
! each other with each interaction counted only once

subroutine assign_to_neighbourlist_halfint_opt(self, cellin)
	implicit none

    type(neighbrinfo),intent(inout) :: self
    type(cellinfo),intent(in)       :: cellin

	integer                         :: i, j,k !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	integer, dimension(13)          :: icellshift, jcellshift, kcellshift
	real(kind(0.d0))				:: rij2   !magnitude^2 between i and j
	real(kind(0.d0)),dimension(3)   :: ri, rj !Position of molecule i and j
	real(kind(0.d0)),dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldihead, oldi, currenti, oldjhead, oldj, currentj


	!Create Neighbourlist array based on current np
	allocate(self%Nlist(np))
	allocate(self%head(np))
	do i = 1,np
		self%Nlist(i) = 0	!Zero number of molecules in neighbour list
		nullify(self%head(i)%point)!Nullify neighbour list head pointer 
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
		cellnp = cellin%cellnp(icell,jcell,kcell)
		oldi => cellin%head(icell,jcell,kcell)%point 
		
		!print('(a,3i6,3i6,i6)'), 'ncells, icell, jcell, kcell, cellnp', ncells, icell, jcell, kcell, cellnp

		!Check interaction within own cell once
		do i = 1,cellnp                 !Step through each particle in list 
			molnoi = oldi%molno 	!Number of molecule
			ri = r(:,molnoi)        !Retrieve ri
			oldj => oldi%next	!Point j molecule to next molecule to i

			do j = i+1,cellnp          !Step through all j for each i

				molnoj = oldj%molno 	!Number of molecule
				rj = r(:,molnoj)        !Retrieve rj
				currentj => oldj
				oldj => currentj%next   !Use pointer in datatype to obtain next item in list

				if(molnoi==molnoj) cycle!Check to prevent interaction with self
	
				rij(:) = ri(:) - rj(:)  !Evaluate distance between particle i and j
				rij2 = dot_product(rij,rij)	!Square of vector calculated
                
				if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
				if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoi, molnoj)
				!if (rij2 < rneighbr2) print*,'own_cell', molnoi, molnoj
			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo

		!Reset old to first molecule in list
		oldihead => cellin%head(icell,jcell,kcell)%point

		do k = 1,13
			oldi => oldihead
			oldjhead => cellin%head(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))%point
			adjacentcellnp = cellin%cellnp(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))

			do i = 1,cellnp				!Step through each particle in list 
				molnoi = oldi%molno		!Number of molecule
				ri = r(:,molnoi)		!Retrieve ri

				oldj => oldjhead		!Reset j to head of linked list

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	!Number of molecule
					rj = r(:,molnoj)	!Retrieve rj
					currentj => oldj
					oldj => currentj%next	!Use pointer in datatype to obtain next item in list

					rij(:) = ri(:) - rj(:)	!Evaluate distance between particle i and j
					rij2 = dot_product(rij,rij)	!Square of vector calculated

					if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
					!if (rij2 < rneighbr2) print*,'neighbr_cells',  molnoi, molnoj
					if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoi, molnoj)

				enddo

				currenti => oldi
				oldi => currenti%next !Use pointer in datatype to obtain next item in list

			enddo

		enddo
	enddo
	enddo
	enddo

	!Build up interactions using unchecked halo cells
    call all_cell_interactions()

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

contains

    subroutine all_cell_interactions()
	    implicit none

        integer               :: x,y,z,x1,y1,z1,x2,y2,z2
        integer,dimension(27) :: imin,imax,jmin,jmax,kmin,kmax,shift

        x = ncells(1); x1 = x+1; x2 = x+2
        y = ncells(2); y1 = y+1; y2 = y+2
        z = ncells(3); z1 = z+1; z2 = z+2
        imin = [2,1,3,2,2,1,3,3,1,1,1,1,1,1,x2,x2,x2,x2,2,2,3,2,2,3,2,3,2]
        imax = [x1,x,x2,x1,x1,x,x2,x2,x,1,1,1,1,1,x2,x2,x2,x2,x1,x,x1,x1,x,x1,x1,x1,x]
        jmin = [2,2,2,1,3,1,1,3,3,2,2,1,1,3,1,2,1,3,1,1,1,1,1,1,y2,y2,y2]
        jmax = [y1,y1,y1,y,y2,y,y,y2,y2,y1,y1,y,y,y2,y,y1,y,y2,1,1,1,1,1,1,y2,y2,y2]
        kmin = [1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
        kmax = [1,1,1,1,1,1,1,1,1,z1,z,z1,z,z,z1,z,z,z,z1,z1,z1,z,z,z,z,z,z]
        shift= [5,6,10,8,12,7,9,11,13,1,6,2,7,13,4,10,9,11,3,2,4,8,7,9,12,11,13]

        do i = 1,size(shift)
            call calculate_cell_interactions_opt(self, cellin, &
                                                 imin(i),imax(i),icellshift(shift(i)), &
                                                 jmin(i),jmax(i),jcellshift(shift(i)), &
                                                 kmin(i),kmax(i),kcellshift(shift(i))   )
        enddo

    end subroutine all_cell_interactions

end subroutine assign_to_neighbourlist_halfint_opt




!----------------------------------------------------------------------------------
! Assign to Neighbourlist each molecules only once so each interaction is
! counted twice. Halo molecules are also included which is essential for a
! number of calculated properties

subroutine assign_to_neighbourlist_allint_halo(self, cellin)
	implicit none

    type(neighbrinfo),intent(inout) :: self
    type(cellinfo),intent(in)       :: cellin

	integer                         :: i, j !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	real(kind(0.d0))				:: rij2   !magnitude^2 between i and j
	real(kind(0.d0)),dimension(3)   :: ri, rj !Position of molecule i and j
	real(kind(0.d0)),dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	!Create Neighbourlist array based on current np
	allocate(self%Nlist(np+halo_np))
	allocate(self%head(np+halo_np))
	do i = 1,np+halo_np
		self%Nlist(i) = 0			!Zero number of molecules in neighbour list
		nullify(self%head(i)%point)	!Nullify neighbour list head pointer 
	enddo

	!if (maxval(cellin%cellnp(:,:,:)) .gt. 9) call error_abort("ERROR - greater than 10 per cell")

	do icell=1, ncells(1)+2
	do jcell=1, ncells(2)+2
	do kcell=1, ncells(3)+2

		cellnp = cellin%cellnp(icell,jcell,kcell)
		oldi => cellin%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp				!Step through each particle in list 
			molnoi = oldi%molno		!Number of molecule
			ri = r(:,molnoi)		!Retrieve ri

			do icellshift = -1,1
			do jcellshift = -1,1
			do kcellshift = -1,1

				!Prevents out of range values in i
				if (icell+icellshift .lt. 	  1		 ) cycle
				if (icell+icellshift .gt. ncells(1)+2) cycle
				!Prevents out of range values in j
				if (jcell+jcellshift .lt. 	  1		 ) cycle
				if (jcell+jcellshift .gt. ncells(2)+2) cycle
				!Prevents out of range values in k
				if (kcell+kcellshift .lt. 	  1		 ) cycle
				if (kcell+kcellshift .gt. ncells(3)+2) cycle

				oldj => cellin%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cellin%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				!print*, icell+icellshift,jcell+jcellshift,kcell+kcellshift

				do j = 1,adjacentcellnp         !Step through all j for each i

					molnoj = oldj%molno 	 	!Number of molecule
					rj = r(:,molnoj)         	!Retrieve rj

					currentj => oldj
					oldj => currentj%next    	!Use pointer in datatype to obtain next item in list
					
					if(molnoi==molnoj) cycle 	!Check to prevent interaction with self

					rij(:) = ri(:) - rj(:)   	!Evaluate distance between particle i and j
					rij2 = dot_product(rij,rij) !Square of rij
					
					if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)
					if (rij2 < rneighbr2) 	 call linklist_checkpushneighbr(self, molnoi, molnoj)

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

	nullify(oldi)      		!Nullify as no longer required
	nullify(oldj)      		!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine assign_to_neighbourlist_allint_halo


!----------------------------------------------------------------------------------
! Assign to Neighbourlist count each interaction only once
! Build a list connecting molecules which are within interacting distance of
! each other with each interaction counted only once.
! Halo molecules are also included which is essential for a
! number of calculated properties

subroutine assign_to_neighbourlist_halfint_halo(self, cellin)
	implicit none

    type(neighbrinfo),intent(inout) :: self
    type(cellinfo),intent(in)       :: cellin

	integer                         :: i,j,k !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	integer, dimension(13)          :: icellshift, jcellshift, kcellshift
	real(kind(0.d0))				:: rij2   !magnitude^2 between i and j
	real(kind(0.d0)),dimension(3)   :: ri, rj !Position of molecule i and j
	real(kind(0.d0)),dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldihead, oldi, currenti, oldjhead, oldj, currentj

	!Create Neighbourlist array based on current np
	allocate(self%Nlist(np+halo_np))
	allocate(self%head(np+halo_np))
	do i = 1,np+halo_np
		self%Nlist(i) = 0	!Zero number of molecules in neighbour list
		nullify(self%head(i)%point)!Nullify neighbour list head pointer 
	enddo

	!Assign cell offsets
	icellshift = (/ 1, 1, 0,-1, 0, 1, 1, 0,-1,-1,-1, 0, 1/)
	jcellshift = (/ 0, 1, 1, 1, 0, 0, 1, 1, 1, 0,-1,-1,-1/)
	kcellshift = (/ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1/)

	!Calculate interactions between all cells in domain and halos
	do icell=1, ncells(1)+2
	do jcell=1, ncells(2)+2
	do kcell=1, ncells(3)+2

		!Retrieve cell np and set old to first molecule in list
		cellnp = cellin%cellnp(icell,jcell,kcell)
		oldi => cellin%head(icell,jcell,kcell)%point 

		!Check interaction within own cell once
		do i = 1,cellnp                 !Step through each particle in list 
			molnoi = oldi%molno 		!Number of molecule
			ri = r(:,molnoi)        	!Retrieve ri
			oldj => oldi%next			!Point j molecule to next molecule to i

			do j = i+1,cellnp			!Step through all j for each i

				molnoj = oldj%molno 	!Number of molecule
				rj = r(:,molnoj)        !Retrieve rj
				currentj => oldj
				oldj => currentj%next   !Use pointer in datatype to obtain next item in list

				if(molnoi==molnoj) cycle!Check to prevent interaction with self
	
				rij2=0                  !Set rij^2 to zero
				rij(:) = ri(:) - rj(:)  !Evaluate distance between particle i and j

				rij2 = dot_product(rij,rij)	!Square of vector calculated

				if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
				if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoi, molnoj)
				!if (rij2 < rneighbr2) print*,'own_cell', molnoi, molnoj
			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo

		!Reset old to first molecule in list
		oldihead => cellin%head(icell,jcell,kcell)%point

		do k = 1,13

			!Prevents out of range values in i
			if (icell+icellshift(k) .lt. 	  1		) cycle
			if (icell+icellshift(k) .gt. ncells(1)+2) cycle
			!Prevents out of range values in j
			if (jcell+jcellshift(k) .lt. 	  1		) cycle
			if (jcell+jcellshift(k) .gt. ncells(2)+2) cycle
			!Prevents out of range values in k
			if (kcell+kcellshift(k) .lt. 	  1		) cycle
			if (kcell+kcellshift(k) .gt. ncells(3)+2) cycle

			oldi => oldihead
			oldjhead => cellin%head(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))%point
			adjacentcellnp = cellin%cellnp(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))

			do i = 1,cellnp				!Step through each particle in list 
				molnoi = oldi%molno		!Number of molecule
				ri = r(:,molnoi)		!Retrieve ri

				oldj => oldjhead		!Reset j to head of linked list

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	!Number of molecule
					rj = r(:,molnoj)	!Retrieve rj
					currentj => oldj
					oldj => currentj%next	!Use pointer in datatype to obtain next item in list

					rij2=0			!Set rij^2 to zero
					rij(:) = ri(:) - rj(:)	!Evaluate distance between particle i and j

					rij2 = dot_product(rij,rij)	!Square of vector calculated

					if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
					!if (rij2 < rneighbr2) print*,'neighbr_cells',  molnoi, molnoj
					if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoi, molnoj)

				enddo

				currenti => oldi
				oldi => currenti%next !Use pointer in datatype to obtain next item in list

			enddo


		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      		!Nullify as no longer required
	nullify(oldj)      		!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine assign_to_neighbourlist_halfint_halo

!------------------------------------------------------------------------------
!Routine to calculate all molecular interactions between two specified cells
!Used only for halo cell so works out interaction of j with i but adds to
!molecule i's neighbourlist

subroutine calculate_cell_interactions(cellin, icell, jcell, kcell, k)
	implicit none

    type(cellinfo),intent(in)       :: cellin
	integer,intent(in)				:: icell, jcell, kcell, k

	integer							:: i,j
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	integer, dimension(13)			:: icellshift, jcellshift, kcellshift
	real(kind(0.d0))				:: rij2   !magnitude^2 between i and j
	real(kind(0.d0)),dimension(3)   :: ri, rj !Position of molecule i and j
	real(kind(0.d0)),dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldi, currenti, oldjhead, oldj, currentj

	!Assign cell offsets
	icellshift = (/ 1, 1, 0,-1, 0, 1, 1, 0,-1,-1,-1, 0, 1/)
	jcellshift = (/ 0, 1, 1, 1, 0, 0, 1, 1, 1, 0,-1,-1,-1/)
	kcellshift = (/ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1/)

		!Retrieve cell np and set old to first molecule in list
		cellnp = cellin%cellnp(icell,jcell,kcell)
		oldi => cellin%head(icell,jcell,kcell)%point !Set old to first molecule in list
		adjacentcellnp = cellin%cellnp(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))
		oldjhead => cellin%head(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k))%point

		!Count number of interactions per cell
		!intcount(icell, jcell, kcell) = intcount(icell, jcell, kcell) + 1
		!intcount(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k)) = & 
		!intcount(icell+icellshift(k),jcell+jcellshift(k),kcell+kcellshift(k)) + 1

		!write(100,'(6(i8,a))'), icell,';',jcell,';',kcell,';', &
		!icell+icellshift(k),';',jcell+jcellshift(k),';',kcell+kcellshift(k),';'

		!Check interaction with neighbouring cells
		do i = 1,cellnp                 	!Step through each particle in list 
			molnoi = oldi%molno 			!Number of molecule
			ri = r(:,molnoi)        		!Retrieve ri

			oldj => oldjhead				!Reset j to head of linked list

			do j = 1,adjacentcellnp         !Step through all j for each i

				molnoj = oldj%molno 		!Number of molecule
				rj = r(:,molnoj)         	!Retrieve rj
				currentj => oldj
				oldj => currentj%next    	!Use pointer in datatype to obtain next item in list

				rij(:) = ri(:) - rj(:)   	!Evaluate distance between particle i and j
				rij2 = dot_product(rij,rij)	!Square of vector calculated

				if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
				!Used for halo cell so molnoj and molnoi swapped over!
				if (rij2 < rneighbr2) call linklist_checkpushneighbr(neighbour, molnoj, molnoi)

			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo

    end subroutine calculate_cell_interactions

!------------------------------------------------------------------------------
! * * * Optimised by performing loop inside routine * * *
!Routine to calculate all molecular interactions between two specified cells
!Used only for halo cell so works out interaction of j with i but adds to
!molecule i's neighbourlist


subroutine calculate_cell_interactions_opt(self, cellin, &
                                           istart,iend,ishift, & 
                                           jstart,jend,jshift, & 
                                           kstart,kend,kshift)
	implicit none

    type(neighbrinfo),intent(inout) :: self
    type(cellinfo),intent(in)       :: cellin

	integer,intent(in)				:: istart,iend,ishift,jstart,jend,jshift,kstart,kend,kshift

	integer							:: i,j
	integer							:: icell, jcell, kcell
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	real(kind(0.d0))				:: rij2   !magnitude^2 between i and j
	real(kind(0.d0)),dimension(3)   :: ri, rj !Position of molecule i and j
	real(kind(0.d0)),dimension(3)   :: rij    !vector between particles i and j
	type(node), pointer 	        :: oldi, currenti, oldjhead, oldj, currentj


	do kcell=kstart,kend
	do jcell=jstart,jend
	do icell=istart,iend

		cellnp = cellin%cellnp(icell,jcell,kcell)
		oldi => cellin%head(icell,jcell,kcell)%point !Set old to first molecule in list

		!Retrieve cell np and set old to first molecule in list
		adjacentcellnp = cellin%cellnp(icell+ishift,jcell+jshift,kcell+kshift)
		oldjhead => cellin%head(icell+ishift,jcell+jshift,kcell+kshift)%point

		!Check interaction with neighbouring cells
		do i = 1,cellnp                 	!Step through each particle in list 
			molnoi = oldi%molno 			!Number of molecule
			ri = r(:,molnoi)        		!Retrieve ri

			oldj => oldjhead				!Reset j to head of linked list

			do j = 1,adjacentcellnp         !Step through all j for each i

				molnoj = oldj%molno 		!Number of molecule
				rj = r(:,molnoj)         	!Retrieve rj
				currentj => oldj
				oldj => currentj%next    	!Use pointer in datatype to obtain next item in list

				rij(:) = ri(:) - rj(:)   	!Evaluate distance between particle i and j
				rij2 = dot_product(rij,rij)	!Square of vector calculated

				if (potential_flag.eq.1) call check_update_adjacentbeadinfo(molnoi,molnoj)	
				if (rij2 < rneighbr2) call linklist_checkpushneighbr(self, molnoj, molnoi)

			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo

	enddo
	enddo
	enddo

end subroutine calculate_cell_interactions_opt


!------------------------------------------------------------------------
! Test routine which allows a local cell and neighbour list to be 
! built using a user specified (i.e. subset or test) array of molecular
! positions

subroutine build_cell_and_neighbourlist_using_debug_positions(rdebug, celldebug, neighbourdebug)
    use physical_constants_MD, only : np
    use computational_constants_MD, only : halfdomain, cellsidelength, nh, ncells, rneighbr2
    use arrays_MD, only : r
    implicit none

    type(neighbrinfo),intent(out) :: neighbourdebug
    type(cellinfo),intent(out)   :: celldebug

    real(kind(0.d0)),dimension(:,:),allocatable,intent(in)    :: rdebug

    integer		:: i, j, n, m
    integer		:: cellnp, adjacentcellnp,  testmols, molcount, molnoi, molnoj
    integer		:: icell, jcell, kcell, icellshift, jcellshift, kcellshift
    real(kind(0.d0))                :: rd2, rij2
    real(kind(0.d0)),dimension(3)   :: ri, rj, rij

	type(node), pointer 	        :: oldi, currenti, oldj, currentj

    testmols = size(rdebug,2)

    !Allocate and zero new cell lists
    allocate(celldebug%head(ncells(1)+2,ncells(2)+2,ncells(3)+2))
    allocate(celldebug%cellnp(ncells(1)+2,ncells(2)+2,ncells(3)+2))
    do icell =1,ncells(1)+2
    do jcell =1,ncells(2)+2
    do kcell =1,ncells(3)+2
	    celldebug%cellnp(icell,jcell,kcell)=0 !Zero number of molecules per cell
	    nullify(celldebug%head(icell,jcell,kcell)%point) !Nullify cell head pointers
    enddo
    enddo
    enddo

    !Build cell list
    do n=1,testmols
	    icell = ceiling((rdebug(1,n)+halfdomain(1)) &
	    /cellsidelength(1))+nh !Add 1 due to halo
	    jcell = ceiling((rdebug(2,n)+halfdomain(2)) &
	    /cellsidelength(2))+nh !Add 1 due to halo
	    kcell = ceiling((rdebug(3,n)+halfdomain(3)) &
	    /cellsidelength(3))+nh !Add 1 due to halo

	    call linklist_checkpush(celldebug, icell, jcell, kcell, n)
    enddo

    !Build neighbour list
    !Create Neighbourlist array based on current np
    allocate(neighbourdebug%Nlist(np))
    allocate(neighbourdebug%head(np))
    do i = 1,np
	    neighbourdebug%Nlist(i) = 0	!Zero number of molecules in neighbour list
	    nullify(neighbourdebug%head(i)%point)!Nullify neighbour list head pointer 
    enddo

    do icell=2, ncells(1)+1
    do jcell=2, ncells(2)+1
    do kcell=2, ncells(3)+1

	    cellnp = celldebug%cellnp(icell,jcell,kcell)
	    oldi => celldebug%head(icell,jcell,kcell)%point !Set old to first molecule in list

	    do i = 1,cellnp                  !Step through each particle in list 
		    molnoi = oldi%molno 	 !Number of molecule
		    ri = rdebug(:,molnoi)         !Retrieve ri

		    do icellshift = -1,1
		    do jcellshift = -1,1
		    do kcellshift = -1,1
			    oldj => celldebug%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
			    adjacentcellnp = celldebug%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

			    do j = 1,adjacentcellnp         !Step through all j for each i

				    molnoj = oldj%molno 	 	!Number of molecule
				    rj = rdebug(:,molnoj)         	!Retrieve rj

				    currentj => oldj
				    oldj => currentj%next    	!Use pointer in datatype to obtain next item in list
				
				    if(molnoi .eq. molnoj) cycle 	!Check to prevent interaction with self

				    rij(:) = ri(:) - rj(:)   	!Evaluate distance between particle i and j
				    rij2 = dot_product(rij,rij) !Square of rij

				    if (rij2 < rneighbr2) call linklist_checkpushneighbr(neighbourdebug, molnoi, molnoj)

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

    nullify(oldi)      		!Nullify as no longer required
    nullify(oldj)      		!Nullify as no longer required
    nullify(currenti)      	!Nullify as no longer required
    nullify(currentj)      	!Nullify as no longer required

end subroutine  build_cell_and_neighbourlist_using_debug_positions

!======================================================================
!		Linklist manipulation Subroutines                     =
!======================================================================

!-----------------------------------------------------------------------
!Build linklist from array of positions of molecules
! NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED 
subroutine linklist_build(cellin, start, finish)
	implicit none

    type(cellinfo),intent(out)       :: cellin

	integer            :: j
	integer            :: cellnp
	integer            :: icell, jcell, kcell	
	integer            :: start, finish
	type(node), pointer:: old, current

	print*, 'building linklist'
	allocate (old)            !Allocate storage space for pointer old to point at
	!old%rp => r(start,:)          !Point to first molecule's r
	!old%vp => v(start,:)          !Set up pointer to current molecule's v
	!old%ap => a(start,:)          !Set up pointer to current molecule's a
	old%molno  = start       !Assign first molecule number 1
	nullify (old%next)        !Points to null as first in linked list
	
	do j=start+1,finish
		allocate (current)          !Allocate storage space for current to point at
		!current%rp  => r(:,j)       !Set up pointer to current molecule's r
		!current%vp  => v(:,j)       !Set up pointer to current molecule's v
		!current%ap  => a(:,j)       !Set up pointer to current molecule's a
		current%molno  = j     !Assign each molecule a number
		old%previous => current     !Set up forward pointer to allow movement forward
		current%next => old         !Set current's link pointer to point to previous item
		old     => current          !Set current item to old ready for next loop
	enddo

	cellnp = finish - start+1
	cellin%cellnp(icell,jcell,kcell) = cellnp
	nullify(current)      !Nullify current as no longer required
	nullify(old%previous) !Points top to null as last in linked list


end subroutine linklist_build
! NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED 

!===================================================================================
!Build circular linklist from array of positions of molecules
! NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED 
subroutine linklist_circbuild(start, finish)
	implicit none

	integer              :: j
	integer              :: cellnp
	integer              :: start, finish
	type(node), pointer  :: old, current, bottom
        
	print*, 'building circular linklist'
	allocate (bottom)         !Allocate storage space for pointer bottom to point at
	bottom%molno  = start     !Assign first molecule number

	allocate (old)            !Allocate storage space for pointer old to point at
	old%molno  = start+1      !Assign second molecule number
	bottom%previous => old    !Set up forward pointer to allow movement forward 
	old%next => bottom 	      !Set up first elements pointer to bottom of list

	do j=start+2,finish
		allocate (current)          !Allocate storage space for current to point at
		current%molno   = j         !Assign each molecule a number
		old%previous => current     !Set up forward pointer to allow movement forward
		current%next => old         !Set current's link pointer to point to previous item
		old   => current            !Set current item to old ready for next loop
	enddo

	bottom%next => old  !Point next pointer of bottom to top of list - circular
	old%previous => bottom 	!Point previous pointer of top of list to bottom - circular
	cellnp = finish - start+1

	nullify(current) !Nullify current as no longer required
	
end subroutine linklist_circbuild
! NOT USED NOT USED NOT USED NOT USED NOT USED NOT USED 

!===================================================================================
! Remove ('pop') a molecule from the stack and returns its cell

! Find which cell list molecule "molno" is associated with (it might not necessarily
! be calculatable with integer division if it has strayed between rebuilds)
subroutine linklist_findandpop(self, molno, icell_return, jcell_return, kcell_return)
    implicit none

    type(cellinfo),intent(inout)    :: self

    integer, intent(in) :: molno
    integer, intent(out) :: icell_return, jcell_return, kcell_return

    integer :: icell, jcell, kcell
    integer :: icell_guess, jcell_guess, kcell_guess
    integer :: ishift, jshift, kshift
	type(node), pointer :: pop
    logical :: found=.false.

    ! Guess cell by integer division first
    icell_guess = ceiling((r(1,molno)+halfdomain(1)) &
    /cellsidelength(1))+nh !Add nh due to halo(s)
    jcell_guess = ceiling((r(2,molno)+halfdomain(2)) &
    /cellsidelength(2))+nh !Add nh due to halo(s)
    kcell_guess = ceiling((r(3,molno)+halfdomain(3)) &
    /cellsidelength(3))+nh !Add nh due to halo(s)

    icell = icell_guess
    jcell = jcell_guess
    kcell = kcell_guess

    call linklist_gotomolecule(self, icell, jcell, kcell, molno, pop, found) !TODO Make this optional

    ! If molecule in guessed cell (should be most of the time)
    if (found) then

        call linklist_pop(self, icell, jcell, kcell, molno)
        icell_return = icell
        jcell_return = jcell
        kcell_return = kcell   
        return

    ! If molecule has strayed out of the cell region it is associated with
    ! (we don't know this yet) then search neighbouring cells to the guess.
    else

        do ishift = -1,1
        do jshift = -1,1
        do kshift = -1,1

            icell = icell_guess + ishift
            jcell = jcell_guess + jshift
            kcell = kcell_guess + kshift

            ! If neighbour cell is out of range ignore
            if (icell .gt. ncells(1)+2) cycle
            if (jcell .gt. ncells(2)+2) cycle
            if (kcell .gt. ncells(3)+2) cycle
            if (icell .lt. 1) cycle
            if (jcell .lt. 1) cycle
            if (kcell .lt. 1) cycle

            call linklist_gotomolecule(self, icell, jcell, kcell, molno, pop, found)

            if (found) then

                call linklist_pop(self, icell, jcell, kcell, molno)
                icell_return = icell
                jcell_return = jcell
                kcell_return = kcell   
                return 

            end if

        end do
        end do
        end do

        if (.not. found) then
            print'(a,i8,6(a,i6))', 'Molecule', molno, ' not found in x cells from ', icell_guess-1, ' to ', icell_guess+1, & 
                                                                   ' y cells from ', jcell_guess-1, ' to ', jcell_guess+1, & 
                                                                   ' z cells from ', kcell_guess-1, ' to ', kcell_guess+1
             stop "ERROR in linklist_findandpop - molecule not found "
        endif

    end if

 

end subroutine linklist_findandpop

subroutine linklist_pop(self, icell, jcell, kcell, molnopop)
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer, intent(in)    	       :: icell, jcell, kcell
	integer, intent(in)    	       :: molnopop

	integer            	           :: cellnp
	type(node), pointer 	       :: old, current
	type(node), pointer 	       :: pop
    logical :: foundmol 

	call linklist_gotomolecule(self, icell, jcell, kcell, molnopop, pop, foundmol)
	if (foundmol .eqv. .true.) then    !Exit if null

        if (molnopop .eq. pop%molno) then !Exit if not correct molno

            !print*, 'pop', pop%molno, 'from cell', icell, jcell, kcell

            cellnp = self%cellnp(icell,jcell, kcell)
            
            !Check if popped molecule is the one head pointer is pointing to
            if (associated(self%head(icell,jcell, kcell)%point,pop)) then

                !Check there are other molecules in cell
                if (associated(pop%next)) then
                    !Set head pointer to next in list and remove top
                    self%head(icell,jcell, kcell)%point => pop%next 
                    old     => pop%next         !Set old to next item in list 
                    nullify(old%previous)
                else
                    !If none then cell pointer is nullified and the list is empty
                    nullify(self%head(icell,jcell,kcell)%point)
                endif
                deallocate(pop)				!Destroy pop

            else

                !Check if popped molecule is last in list
                if (associated(pop%next)) then
                    !If just an element in list, remove it
                    old     => pop%next         !Set old to next item in list
                    current => pop%previous     !Set current to previous item in list
                    deallocate(pop)				!Destroy pop

                    old%previous => current     !Previous pointer connects to old (pop missed out)
                    current%next => old     	!Next pointer connects to current (pop missed out)

                else
                    !If last in list
                    current => pop%previous     !Set current to previous item in list
                    deallocate(pop)				!Destroy pop

                    nullify(current%next)		!Next pointer connects to null (now last in list)
                endif

            endif
        
            cellnp = cellnp - 1	          !Reduce number of molecules by one
            self%cellnp(icell,jcell, kcell) = cellnp !Update cell molecular number
            !print*, 'new cell np', cellnp

            nullify(current) !Nullify current as no longer required
            nullify(old)     !Nullify old as no longer required

        endif

        nullify(pop)  	 !Nullify pop as no longer required
        return

    else
       
        print*, 'linklist_pop failed to pop molecule ', molnopop, 'from cell',&
                                                        icell, jcell, kcell 
        stop 

	endif

end subroutine linklist_pop

!===================================================================================
!Adds molecule specified by passed variables moleculenpush to linked list

subroutine linklist_push(self, icell, jcell, kcell, molnopush)
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	integer, target     :: molnopush
	type(node), pointer :: old, current
	type(node), pointer :: push

	!Obtain molecular number and top item of link list from cell head
	old => self%head(icell,jcell, kcell)%point
	cellnp = self%cellnp(icell,jcell, kcell)

	print*, 'push', molnopush

	allocate(push) !Allocate type to add to stack
	push%molno = molnopush !Build type from inputs

	current => old           !Set current to first item in list
	old     => current%next  !Set old to next item in list
	
	old%previous => push     !Old previous pointer connects to new push
	current%next => push     !Current next pointer connects to new push
	push%next    => old      !Push next pointer connects to old
	push%previous=> current  !Push previous pointer connects to current
	old => push

	self%head(icell,jcell,kcell)%point => old !Set cell pointer to top of cell list
	cellnp = cellnp + 1                 !Increase number of molecules by one
	self%cellnp(icell,jcell,kcell) = cellnp   !Update cell molecular number

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required
	nullify(push)                       !Nullify old as no longer required

end subroutine linklist_push

!===================================================================================
!Adds molecule specified by passed variables to linked list with check to see if
!linklist is empty so new linklist can be established

subroutine linklist_checkpush(self, icell, jcell, kcell, molnopush)
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	integer	            :: molnopush
	type(node), pointer :: old, current

	cellnp = self%cellnp(icell,jcell,kcell)
	allocate(current) 		        !Allocate type to add to stack
	current%molno = molnopush 	    !Build type from inputs

	select case (cellnp)
	case(0)
		nullify(current%previous) !Nullify pointer at top of list
		nullify(current%next)     !Nullify pointer at bottom of list
	case(1:)
		old => self%head(icell,jcell,kcell)%point !Set old to top of list
		current%next => old                	!Set old to next item in list
		old%previous => current            	!Old previous pointer connects to new current
		nullify(current%previous)          	!Nullify pointer at top of list
	end select

	self%head(icell,jcell,kcell)%point => current 	!Set cell pointer to top of cell list
	cellnp = cellnp + 1                     		!Increase number of particles by one
	self%cellnp(icell,jcell,kcell) = cellnp       	!Update cell molecular number

	!Nullify pointers current and old so cell head is left pointing to top of the
	!required linked list
	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_checkpush

!===================================================================================
!Adds molecule specified by passed variables to neighbour linked list with check if
!linklist is empty so new linklist can be established

subroutine linklist_checkpushneighbr(self, molnoi, molnoj)
	implicit none

    class(neighbrinfo),intent(inout) :: self

	integer			  			:: noneighbrs
	integer                    	:: molnoi, molnoj
	type(node), pointer 	    :: old, current

	noneighbrs = self%Nlist(molnoi)
	allocate(current)                         	!Allocate type to add to stack
	current%molno = molnoj                   	!Build type from inputs

	select case (noneighbrs)
	case(0)
		nullify(current%previous)         		!Nullify pointer at top of list
		nullify(current%next)             		!Nullify pointer at bottom of list
	case(1:)
		old => self%head(molnoi)%point		    !Set old to top of list
		current%next => old						!Set old to next item in list
		old%previous => current					!Old previous pointer connects to new current
		nullify(current%previous)				!Nullify pointer at top of list
	end select

	self%head(molnoi)%point => current		    !Set cell pointer to top of cell list
	noneighbrs = noneighbrs + 1               	!Increase number of particles by one
	self%Nlist(molnoi) = noneighbrs	            !Update neighbour list molecular number

	!Nullify pointers current and old so neighbour head is left pointing to top
	nullify(current)                          !Nullify current as no longer required
	nullify(old)                              !Nullify old as no longer required

end subroutine linklist_checkpushneighbr

!===================================================================================
!Adds molecule specified by passed variables to passed molecule linked list with check if
!linklist is empty so new linklist can be established

subroutine linklist_checkpushmol(self, molno,ipass,jpass,kpass)
	implicit none

    type(passinfo),intent(inout)    :: self

	integer,intent(in)			:: molno
	integer			   			:: sendnp
	integer						:: ipass, jpass, kpass
	type(passnode), pointer		:: old, current

	sendnp = self%sendnp
	allocate(current)                 	!Allocate type to add to stack
	current%molno = molno				!Build type from inputs
	current%ipass = ipass			  	!Build type from inputs
	current%jpass = jpass			  	!Build type from inputs
	current%kpass = kpass			  	!Build type from inputs

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

	self%head => current                      !Set cell pointer to top of cell list
	sendnp = sendnp + 1                   	  !Increase number of particles by one
	self%sendnp = sendnp                  	  !Update neighbour list molecular number

	!Nullify pointers current and old so neighbour head is left pointing to top
	nullify(current)                          !Nullify current as no longer required
	nullify(old)                              !Nullify old as no longer required

end subroutine linklist_checkpushmol

!===================================================================================
! Merge linklist delete into linklist keep 

subroutine linklist_merge(self, keep, delete)
	implicit none

    integer,intent(in)               :: keep, delete
    class(neighbrinfo),intent(inout) :: self

    integer             :: m, b4
	type(node), pointer :: old, current

    !Save number of molecules before
    b4 = self%Nlist(delete)+self%Nlist(keep)

    current => self%head(delete)%point
    do m = 1,self%Nlist(delete)

        !Add current smaller cluster molecule to the bigger cluster
        call linklist_checkpushneighbr(self, keep, current%molno)

        !Pop molecule from list
        if (associated(current%next) .eqv. .true. ) then 
            old => current%next       !Make list point to next node
        else
            nullify(old)
        endif
        !Deallocate current
        nullify(current%next)     !Remove pointer to next
        nullify(current%previous) !Remove pointer to previous
        deallocate(current)       !Deallocate current entry
        current => old
    enddo

    !Set delete to zero
    nullify(self%head(delete)%point)
    self%Nlist(delete) = 0

    !Check total number of molecules after
    if (self%Nlist(delete)+self%Nlist(keep) .ne. b4) then
        stop "Error in linklist_merge -- Number of elements is not conserved"
    endif

end subroutine linklist_merge

!===================================================================================
!Move backwards through linked list and print out results
	
subroutine linklist_print(self, icell, jcell, kcell)
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer             :: j
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current

	!Obtain molecular number and top item of link list from cell head
	old => self%head(icell,jcell, kcell)%point
	cellnp = self%cellnp(icell,jcell, kcell)

	print*, 'linklist print called for cell', icell, jcell, kcell
	if(cellnp == 0) print*, 'linklist empty'

	current => old ! make current point to head of list
	do j=1,cellnp
		!print*, 'more items in linked list?: ', associated(old%next)
		print'(i8,3f10.5)', current%molno, r(:,current%molno)
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
	
subroutine linklist_printallcells(self, print_empty)
	implicit none

    type(cellinfo),intent(inout)    :: self
    logical,intent(in),optional     :: print_empty

    logical             :: print_empty_
	integer             :: j
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current

    if (present(print_empty)) then
    	print*, 'linklist print called for only cells with molecules'
        print_empty_ = print_empty
    else
    	print*, 'linklist print called for all cells'
        print_empty_ = .true.
    endif

	do icell=1,ncells(1)+2
	do jcell=1,ncells(2)+2
	do kcell=1,ncells(3)+2

		!Obtain molecular number and top item of link list from cell head
		old => self%head(icell,jcell, kcell)%point
		cellnp = self%cellnp(icell,jcell, kcell)

        if(cellnp == 0) then
            if (print_empty_) then
        		print*, 'Cell', icell, jcell, kcell
	        	print*, 'linklist empty'
            endif
        else
      		print*, 'Cell', icell, jcell, kcell
        endif

		current => old ! make current point to head of list
		do j=1,cellnp
			!print*, 'more items in linked list?: ', associated(old%next)
			print'(i5,3f10.5)', current%molno, r(:,current%molno)
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


subroutine linklist_writeallcells(self)
    use librarymod
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer             :: j, f
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current
    character(len=8) :: rankstr

    f = get_new_fileunit()
    write(rankstr,'(i8.8)') irank 
    open(unit=f,file=trim(prefix_dir)//"all_celllists_"//rankstr, &
    status='replace')
    
	print('(a,3i5,a,3f10.4)'), 'linklist_writeallcells, ncells = ', ncells, ', halfdomain = ', halfdomain

	do icell=1,ncells(1)+2
	do jcell=1,ncells(2)+2
	do kcell=1,ncells(3)+2

		!Obtain molecular number and top item of link list from cell head
		old => self%head(icell,jcell, kcell)%point
		cellnp = self%cellnp(icell,jcell, kcell)

		write(f, '(a,3i10)') 'Cell', icell, jcell, kcell
		if(cellnp == 0) write(f,'(a)') 'linklist empty'

		current => old ! make current point to head of list
		do j=1,cellnp
			!print*, 'more items in linked list?: ', associated(old%next)
			write(f, '(i10, 3f15.4)') current%molno, r(:,current%molno)
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

    close(f,status='keep')

end subroutine linklist_writeallcells

!===================================================================================
!Move through all cells linked list not including the halo cells and print out 
!results reseting back to the top upon completion 
	
subroutine linklist_printalldomaincells(self)
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer             :: j
	integer             :: cellnp		
	integer             :: icell, jcell, kcell
	type(node), pointer :: old, current

	print*, 'linklist print called for all domain cells'

	do icell=2,ncells(1)+1
	do jcell=2,ncells(2)+1
	do kcell=2,ncells(3)+1

		!Obtain molecular number and top item of link list from cell
		old => self%head(icell,jcell, kcell)%point
		cellnp = self%cellnp(icell,jcell, kcell)

		print*, 'Cell', icell, jcell, kcell
		if(cellnp == 0) print*, 'linklist empty'

		current => old ! make current point to head of list
		do j=1,cellnp
			!print*, 'more items in linked list?: ', associated(old%next)
			print*, current%molno, r(:,current%molno)
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
	
subroutine linklist_printneighbourlist(self, molno, unitno)
	implicit none

    class(neighbrinfo),intent(inout) :: self
	integer,intent(in)			    :: molno
	integer,intent(in),optional		:: unitno

	integer						:: noneighbrs, j
	type(node), pointer	:: old, current

	!Obtain molecular number and top item of link list from cell head
	noneighbrs = self%Nlist(molno)  	!Determine number of elements in neighbourlist
	old => self%head(molno)%point		!Set old to head of neighbour list

	if(noneighbrs == 0) print*, molno, 'has 0 neighbours'

	current => old ! make current point to head of list
	do j=1,noneighbrs
        
		!print*, 'more items in linked list?: ', associated(old%next)
        if (present(unitno)) then
            write(unitno,'(i6, 2(a,i8),6f10.5)') j, &
            ' Linklist print called for i = ', molno,' j = ', &
            current%molno, r(:,molno), r(:,current%molno)
        else
    		print'(i6, 2(a,i8),6f10.5)', j, &
            ' Linklist print called for i = ', molno,' j = ', &
            current%molno, r(:,molno), r(:,current%molno)
        endif
		if (associated(old%next) .eqv. .true. ) then !Exit if null
			old => current%next ! Use pointer in datatype to obtain next item in list
			current => old      ! make current point to old - move alone one
		endif
	enddo

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_printneighbourlist

!===================================================================================
!Move backwards through linked list and print out results
	
subroutine linklist_printpassedlist(self)
	implicit none

    type(passinfo),intent(inout)    :: self

	integer                    :: j
	integer			   :: sendnp
	type(passnode), pointer    :: old, current

	!Obtain molecular number and top item of link list from cell head
    sendnp = self%sendnp	!Determine number of elements in neighbourlist
	old => self%head	    !Set old to head of neighbour list

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
	
subroutine linklist_compareprint(self, icell, jcell, kcell)
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer            	:: j         !Dummy counter
	integer		   		:: n         !Molecule number
	integer            	:: cellnp		
	integer            	:: icell, jcell, kcell
	type(node), pointer	:: old, current

	!Obtain molecular number and top item of link list from cell
	old => self%head(icell,jcell, kcell)%point
	cellnp = self%cellnp(icell,jcell, kcell)

	!call linklist_gototop(icell, jcell, kcell) 
	!print*, 'Print array and linklist for comparison has been called - Should be identical'
	current => old
	do j=1,cellnp
		n = old%molno
		print'(i0,a,2f20.15,a,2f10.5,a,2f10.5)',old%molno, ' ap = ', a(:,old%molno), &
			 ' vp = ', v(:,old%molno), ' rp = ', r(:,old%molno)
		print'(i0,a,2f20.15,a,2f10.5,a,2f10.5)',n, ' a  = ',a(:,n), ' v  = ', &
			 v(:,n), ' r  = ', r(:,n)
		old => current%next !Use pointer in datatype to obtain next item in list
		current => old          !make current point to old - move along one
	enddo

	!cellin%head(icell,jcell, kcell)%point => old

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine linklist_compareprint

!===================================================================================
!Deallocate cell linklist

subroutine linklist_deallocate(self, icell, jcell, kcell)
	implicit none


    type(cellinfo),intent(inout)    :: self
	integer, intent(in)             :: icell, jcell, kcell

	integer            :: j
	integer            :: cellnp
	type(node), pointer:: old, current		

	if (associated(self%head(icell,jcell, kcell)%point) .eqv. .true. ) then !Exit if null

		!Obtain molecular number and top item of link list from cell head
		old => self%head(icell,jcell, kcell)%point
		cellnp = self%cellnp(icell,jcell, kcell)

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
		nullify(self%head(icell,jcell, kcell)%point) !Set cell head pointer to null
		self%cellnp(icell,jcell, kcell) = 0          !Zero cell molecule number
	endif

end subroutine linklist_deallocate

!===================================================================================
!Deallocate bin linklist

subroutine linklist_deallocate_bins(self)
	use calculated_properties_MD
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer            	:: j
	integer            	:: binnp
	integer           	:: ibin, jbin, kbin
	type(node), pointer	:: old, current

	do ibin = 1,1
	do jbin = 1,nbins(1)
	do kbin = 1,1	

        call linklist_deallocate(self, ibin, jbin, kbin)

!		if (associated(bin%head(ibin,jbin,kbin)%point) .eqv. .true. ) then !Exit if null

!			!Obtain molecular number and top item of link list from bin head
!			old => bin%head(ibin,jbin,kbin)%point
!			binnp = bin%cellnp(ibin,jbin,kbin)

!			current => old ! make current point to head of list
!			do j=1,binnp-1
!				if (associated(old%next) .eqv. .true. ) then !Exit if null
!					old => current%next       !make list point to next node of old
!					nullify(current%next)     !Remove pointer to next
!					nullify(current%previous) !Remove pointer to previous
!					deallocate(current)       !Deallocate current entry
!					current => old            !Make current point to new previous
!				endif
!			enddo
!			nullify(old%next)     !Remove pointer to next
!			nullify(old%previous) !Remove pointer to previous
!			deallocate(old)       !Deallocate final entry
!			nullify(bin%head(ibin,jbin,kbin)%point)  !Set bin head pointer to null
!			bin%cellnp(ibin,jbin, kbin) = 0          !Zero bin molecule number
!		endif

	enddo
	enddo
	enddo

end subroutine linklist_deallocate_bins


!===================================================================================
!Deallocate all linklists

subroutine linklist_deallocatepasslist
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
!Deallocate all cell linklists

subroutine linklist_deallocate_cells(self)
	use polymer_info_MD, only: bond, bondcount
	implicit none

    type(cellinfo),intent(inout)    :: self

	integer            			:: i, j
	integer           			:: cellnp, noneighbrs, np_neigbrs
	integer            			:: icell, jcell, kcell
	type(node), pointer			:: old, current

	do icell=1,ncells(1)+2
	do jcell=1,ncells(2)+2
	do kcell=1,ncells(3)+2


        call linklist_deallocate(self, icell, jcell, kcell)

!		if (associated(self%head(icell,jcell,kcell)%point) .eqv. .true. ) then !Exit if null
!			old => self%head(icell,jcell,kcell)%point
!			cellnp = self%cellnp(icell,jcell, kcell)
!			current => old ! make current point to head of list
!			do j=1,cellnp-1
!				if (associated(old%next) .eqv. .true. ) then !Exit if null
!					old => current%next       !Make list point to next node of old
!					nullify(current%next)     !Remove pointer to next
!					nullify(current%previous) !Remove pointer to previous
!					deallocate(current)       !Deallocate current entry
!					current => old            !Make current point to new previous
!				endif
!			enddo

!			nullify(old%next)     !Remove pointer to next
!			nullify(old%previous) !Remove pointer to previous
!			deallocate(old)       !Deallocate final entry
!			nullify(self%head(icell,jcell,kcell)%point) !Set cell head pointer to null
!			self%cellnp(icell,jcell, kcell) = 0         !Zero cell molecule number
!		endif
	enddo
	enddo
	enddo

    !We don't deallocate array of cell pointers as this will never change!!
    !deallocate(self%cellnp)
    !deallocate(self%head)

end subroutine linklist_deallocate_cells

!===================================================================================
!Deallocate all neighbour linklists

subroutine linklist_deallocate_neighbour(self)
	use polymer_info_MD, only: bond, bondcount
	implicit none

    type(neighbrinfo),intent(inout)    :: self

	integer            			:: i, j
	integer           			:: noneighbrs, np_neigbrs
	type(node), pointer 	    :: oldn, currentn

	if (force_list .gt. 1) then
		if (size(self%Nlist).eq.np) then
		!if (vflux_outflag.ne.4 ) then
			np_neigbrs = np
		else
			np_neigbrs = np + halo_np
		endif
		do i = 1, np_neigbrs
			if (associated(self%head(i)%point) .eqv. .true. ) then !Exit if null
	       		noneighbrs = self%Nlist(i)  !Determine number of elements in neighbourlist
				oldn => self%head(i)%point		   !Set old to head of neighbour list
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
				nullify(self%head(i)%point)   !Set neighbour head pointer to null
				self%Nlist(i) = 0  !Zero cell molecule number
			endif
		enddo

	    !Deallocate array of molecules neighbourlist pointers
	    deallocate(self%Nlist)
	    deallocate(self%head)
	endif

end subroutine linklist_deallocate_neighbour


subroutine linklist_deallocate_cluster(self)
	implicit none

    type(clusterinfo),intent(inout)    :: self

    integer             :: i, j
	integer             :: noclust
	type(node), pointer :: old, current

	do i = 1, self%Nclust
		if (associated(self%head(i)%point) .eqv. .true. ) then !Exit if null
       		noclust = self%Nlist(i)         !Determine number of elements in clusterlist
			old => self%head(i)%point		!Set old to head of cluster list
		    current => old                  ! make current point to head of list
		    do j=1,noclust-1
			    if (associated(old%next) .eqv. .true. ) then !Exit if null
				    old => current%next        !Make list point to next node of old
				    nullify(current%next)      !Remove pointer to next
				    nullify(current%previous)  !Remove pointer to previous
				    deallocate(current)        !Deallocate current entry
				    current => old             !Make current point to new previous
			    endif
		    enddo
		    nullify(old%next)         !Remove pointer to next
		    nullify(old%previous)     !Remove pointer to previous
		    deallocate(old)           !Deallocate final entry
			nullify(self%head(i)%point)   !Set cluster head pointer to null
			self%Nlist(i) = 0             !Zero cell molecule number
		endif
	enddo

    !Check all other possible molecules
    !E.s. March 2018. this should not be needed, Nclusts is 
    !the number of allocated clusters so arrays should not be np+extralloc
    !large, generally only ~300 clusters
!    do i = 1, np+extralloc
!        if (associated(self%head(i)%point) .eqv. .true.) then
!            old => self%head(i)%point
!            if (associated(old%next) .eqv. .false.) then
!                nullify(self%head(i)%point)   !Set cluster head pointer to null
!            else
!                print*, iter, i, self%Nlist(i), self%inclust(i), self%Nclust, associated(old%next)
!                !stop "ERROR in linklist_deallocate_cluster"
!            endif
!        endif
!    enddo

    !Deallocate array of molecules neighbourlist pointers
    self%Nclust = 0
    deallocate(self%Nlist)
    deallocate(self%head)
    deallocate(self%inclust)

end subroutine linklist_deallocate_cluster

!=============================================================================
!Moves to the bottom of the linklist

subroutine linklist_gotobottom(self, icell, jcell, kcell)
	implicit none

    type(cellinfo),intent(inout)  :: self
	integer,intent(in)         :: icell, jcell, kcell

	integer            :: j
	integer            :: cellnp		

	type(node), pointer:: old, current

	old => self%head(icell,jcell, kcell)%point
	cellnp = self%cellnp(icell,jcell, kcell)

	print*, 'linklist go to bottom called'

	current => old ! make current point to head of list
	do j=1,cellnp
		print*, 'more items in linked list?: ', associated(old%next)
		print*, current%molno, r(:,current%molno)
		if (associated(old%next) .eqv. .true. ) then !Exit if null
			old => current%next ! make list point to next node of old
			current => old      ! make current point to new previous
		endif
	enddo

	self%head(icell,jcell, kcell)%point => old

	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_gotobottom

!============================================================================
!Moves to the top of the linklist

subroutine linklist_gototop(self, icell, jcell, kcell)
	implicit none

    type(cellinfo),intent(inout)  :: self
	integer,intent(in)         :: icell, jcell, kcell

	integer            :: j
	integer            :: cellnp		
	type(node), pointer:: old, current

	old => self%head(icell,jcell, kcell)%point
	cellnp = self%cellnp(icell,jcell, kcell)

	print*, 'linklist go to top called'

	current => old ! make current point to head of list
	do j=1,cellnp
		!print*, 'more items in linked list?: ', associated(old%previous)
		!print*, current%molno, r(:,current%molno)
		if (associated(old%previous) .eqv. .true. ) then !Exit if null
			old => current%previous     !Set up forward pointer to allow movement forward
			current  => old             !Set current item to old ready for next loop
		endif
	enddo

	self%head(icell,jcell, kcell)%point => old

	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_gototop	

!============================================================================
!Finds a specified molecule and sets the cellpointer to that molecule
!If molecule is not found, flag is set to zero

!subroutine linklist_gotomolecule(self, icell, jcell, kcell, n, flag)
!!	implicit none

!	integer             :: j
!	integer, intent(in) :: n
!	integer, intent(out):: flag
!	integer             :: cellnp		
!	integer             :: icell, jcell, kcell
!	type(node), pointer :: old, current

!	old => self%head(icell,jcell, kcell)%point
!	cellnp = self%cellnp(icell,jcell, kcell)

	!print*, 'linklist go to molecule', n, '  called'

!	current => old ! make current point to head of list
!	do j=1,cellnp
!		old => current%previous  !Set up forward pointer to allow movement forward
!		current  => old          !Set current item to old ready for next loop
!		if (old%molno .eq. n ) then !Exit if null
!			print*, 'molecule', n, '  found'
!			self%head(icell,jcell, kcell)%point => old
!			flag = 1 !Set flag to 1 to indicate molecule found
!			return
!		endif
!	enddo

!	flag = 0 !Set flag to zero to indicate molecule not in list
	!print*, '**ERROR - molecule', n, 'is not in this cell**'
!	nullify(current) !Nullify current as no longer required
!	nullify(old)     !Nullify old as no longer required

!end subroutine linklist_gotomolecule


!==========================================================================
!Move a specified number of places upwards or downwards through the linklist
!determine by input n and direction dir which should be either 'goup' or
!down

subroutine linklist_movethrough(self, icell, jcell, kcell, n, dir)
	implicit none

    type(cellinfo),intent(inout)  :: self

	integer            :: n, j
	character(len=4)   :: dir
	integer            :: icell, jcell, kcell
	type(node), pointer:: old, current

	old => self%head(icell,jcell, kcell)%point

	print*, 'linklist move through called'
	current => old ! make current point to head of list

	if (dir == 'goup') then
		do j=1,n 
			!print*, 'more items in linked list?: ', associated(old%previous)
			!print*, current%molno, r(:,current%molno)
			if (associated(old%previous) .eqv. .true. ) then !Exit if null
				old => current%previous  !make old point to previous node of current
				current  => old          !Set current item to old ready for next loop
			endif
		enddo
	elseif (dir == 'down') then
		do j=1,n
			!print*, 'more items in linked list?: ', associated(old%next)
			!print*, current%molno, r(:,current%molno)
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next      !make old point to next node of current
				current => old           !Set current item to old ready for next loop
			endif
		enddo
	else
	print*, 'error - direction specified should be of type character and either *goup* or *down*'
	endif

	self%head(icell,jcell, kcell)%point => old !Assign cell pointer to top of list

	nullify(current) !Nullify current as no longer required
	nullify(old)     !Nullify old as no longer required

end subroutine linklist_movethrough	


!===================================================================================
!Prints the details of the current element of the linked list as well as those either
!side in the linked list

subroutine linklist_printcurrent(self, icell, jcell, kcell)
	implicit none

    type(cellinfo),intent(inout)  :: self

	integer            :: icell, jcell, kcell
	type(node), pointer :: old

	old => self%head(icell,jcell, kcell)%point
	if (associated(old) .eqv. .true.) then !exit if empty

		if (associated(old%next) .eqv. .true. ) then
			print'(a,i8,3f10.5)', 'element below =', old%next%molno, r(:,old%next%molno)
		else
			print*, 'no element below'
		endif

		print'(a,i8,3f10.5)', 'current element =', old%molno, r(:,old%molno)

		if (associated(old%previous) .eqv. .true. ) then
			print'(a,i8,3f10.5)', 'element above = ', old%previous%molno, r(:,old%previous%molno)
		else
			print*, 'no element above'
		endif
	endif

end subroutine linklist_printcurrent

!----------------------------------------------------------------------------------
!- check_update_adjacentbeadinfo
!  updates left and right molnos during rebuild
!----------------------------------------------------------------------------------
subroutine check_update_adjacentbeadinfo(molnoi,molnoj)
	use polymer_info_MD
	implicit none

	integer              :: chaindiff
	integer              :: jscID
	integer, intent(in)  :: molnoi, molnoj
	integer              :: bflag

	chaindiff = monomer(molnoi)%chainID - monomer(molnoj)%chainID      !Check for same chainID

	if (chaindiff.eq.0) then                                           !If same chain
		jscID = monomer(molnoj)%subchainID                             !Get subchainID of molnoj
		bflag = get_bondflag(molnoi,jscID)                             !Check if molnoi/j are connected

		select case (bflag)
		case(1)                                                        !If connected, add to bond list
			bondcount(molnoi) = bondcount(molnoi) + 1
			bondcount(molnoj) = bondcount(molnoj) + 1

			if (bondcount(molnoi).gt.monomer(molnoi)%funcy) call too_many_bonds_error(molnoi)
			if (bondcount(molnoj).gt.monomer(molnoj)%funcy) call too_many_bonds_error(molnoj)

			bond(bondcount(molnoi),molnoi) = molnoj
			bond(bondcount(molnoj),molnoj) = molnoi

		case default                                                   !If not connected, pass
		end select
	end if

contains 
	
	subroutine too_many_bonds_error(molno)
	use interfaces
	implicit none
		
		integer, intent(in) :: molno
        r(0,-1) = 1
		call error_abort('Error: too many bonds for molno ', molno)
	
	end subroutine too_many_bonds_error 

end subroutine check_update_adjacentbeadinfo

subroutine check_update_adjacentbeadinfo_allint(molnoi,molnoj)
	use polymer_info_MD
	implicit none

	integer              :: chaindiff
	integer              :: jscID
	integer, intent(in)  :: molnoi, molnoj
	integer              :: bflag

	chaindiff = monomer(molnoi)%chainID - monomer(molnoj)%chainID      !Check for same chainID

!    if (irank .eq. 1) then
!        if (molnoi .eq. 19118 .or. molnoj .eq. 19118) then
!            write(8001,*) molnoi, molnoj 
!        end if
!    end if

	if (chaindiff.eq.0) then                                           !If same chain
		jscID = monomer(molnoj)%subchainID                             !Get subchainID of molnoj
		bflag = get_bondflag(molnoi,jscID)                             !Check if molnoi/j are connected
		select case (bflag)
		case(1)                                                        !If connected, add to bond list
			bondcount(molnoi) = bondcount(molnoi) + 1
			if (bondcount(molnoi).gt.monomer(molnoi)%funcy) call too_many_bonds_error(molnoi)
			bond(bondcount(molnoi),molnoi) = molnoj
		case default                                                   !If not connected, pass
		end select
	end if

contains 
	
	subroutine too_many_bonds_error(molno)
	use interfaces
	implicit none
		
		integer, intent(in) :: molno
        r(0,-1) = 1
		call error_abort('Error: too many bonds for molno ', molno)
	
	end subroutine too_many_bonds_error 

end subroutine check_update_adjacentbeadinfo_allint


end module linked_list





!======================================================================
!		Cell and neighbour building subroutines               =
!======================================================================

!--------------------------------------------------------
!Routine used to rebuild linked list for each cell

subroutine assign_to_cell
	use physical_constants_MD, only : np
	use computational_constants_MD, only : halfdomain, cellsidelength, &
                                           ncells, nh
	use arrays_MD, only : r
	use interfaces, only: error_abort
    use linked_list, only : linklist_checkpush, cell
	implicit none

	integer		:: n
	integer		:: icell, jcell, kcell

	do n=1,np
		icell = ceiling((r(1,n)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add 1 due to halo
		jcell = ceiling((r(2,n)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add 1 due to halo
		kcell = ceiling((r(3,n)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add 1 due to halo
!#if DEBUG
		if (any((/icell,jcell,kcell/).lt.1) .or. &
		    icell.gt.ncells(1)+2            .or. &
			jcell.gt.ncells(2)+2            .or. &
			kcell.gt.ncells(3)+2          ) then
			call print_mol_escape_error(n)
			call error_abort("Aborting due to escaped mol.")
		end if
!#endif
		call linklist_checkpush(cell, icell, jcell, kcell, n)
	enddo

end subroutine assign_to_cell

!--------------------------------------------------------
!Routine used to rebuild linked list for each halo cell

subroutine assign_to_halocell(start,finish)
	use computational_constants_MD, only : halfdomain, cellsidelength, &
                                           ncells, nh
	use arrays_MD, only : r
    use linked_list, only : linklist_checkpush, cell
	implicit none

	integer		:: n, start, finish
	integer		:: icell, jcell, kcell

	do n=start,finish
		icell = ceiling((r(1,n)+halfdomain(1)) &
		/cellsidelength(1))+nh !Add 1 due to halo
		jcell = ceiling((r(2,n)+halfdomain(2)) &
		/cellsidelength(2))+nh !Add 1 due to halo
		kcell = ceiling((r(3,n)+halfdomain(3)) &
		/cellsidelength(3))+nh !Add 1 due to halo
		call linklist_checkpush(cell, icell, jcell, kcell, n)
	
	enddo
	
end subroutine assign_to_halocell

!----------------------------------------------------------------------------------
! Assign to Neighbourlist                               
! Build a list connecting molecules which are within interacting distance of
! each other ~ re-ordered so that each all i and all j interactions are calculated
! for each cell pair before moving to next cell pair

subroutine assign_to_neighbourlist!(self)
	use computational_constants_MD, only : force_list
    use interfaces
    use linked_list, only : assign_to_neighbourlist_allint, &
                            assign_to_neighbourlist_halfint_opt, &
                            neighbour, cell
    implicit none

    !type(neighbrinfo),intent(inout) :: self

	!Choose between all pairs, cell list or neighbour list
	select case(force_list)
	case(-1:1)
		!Forces calculated using all pairs or cell list so nothing to build
		! ********* Do nothing *******
	case(2)
		!Forces calculated using neighbour lists with all interactions
		call assign_to_neighbourlist_allint(neighbour, cell)
	case(3)
		!Forces calculated using neighbour lists optimised using 
		!Newton's 3rd law to count only half of the interactions
		call assign_to_neighbourlist_halfint_opt(neighbour, cell)
	case default
		call error_abort("Error in force_list flag")
	end select	


end subroutine assign_to_neighbourlist



!Sort molecular locations to improve cache efficiecy using
!blocks of multiple cells

subroutine sort_mols
	use arrays_MD, only : r, v, tag, rtether
	use physical_constants_MD
	use computational_constants_MD
	use interfaces, only : error_abort
	implicit none

	integer											:: iblk,jblk,kblk,start,finish
	integer											:: n,i,ave_molperblock,blocks
	integer,dimension(3)							:: nblocks
	integer,allocatable,dimension(:)				:: molperblock
	integer,allocatable,dimension(:,:)				:: tagtemp
	real(kind(0.d0)),dimension(3)					:: blocksidelength
	real(kind(0.d0)),allocatable,dimension(:,:,:)	:: rtemp,vtemp
	real(kind(0.d0)),allocatable,dimension(:,:,:)	:: rtethertemp
	!Safety factors in allocations
	integer,save									:: sf2=5
	real(kind(0.d0)),save							:: sf1=1.5d0
	!Number of rebuilds since last sort
	integer,save									::  nrebuilds=1

	! Determine if rebuild is required yet
	if (nrebuilds .ne. sort_freq) then
		nrebuilds = nrebuilds + 1
		return
	else
		nrebuilds = 1
	endif

	!Choose between the various sorting methodolgies
	select case(sort_flag)
	case(0)
		return 	!No sort - return
	case(1)
		!Use ordered array generated during setup
		if (potential_flag .eq. 1) call error_abort("Sort should be turned off - Not developed for polymers")
        if (Mie_potential .eq. 1) call error_abort("Sort should be turned off - Not developed for Mie potential")
	case(2)
		!Use Hilbert curve generated during setup
		if (potential_flag .eq. 1) call error_abort("Sort should be turned off - Not developed for polymers")
        if (Mie_potential .eq. 1) call error_abort("Sort should be turned off - Not developed for Mie potential")
	case default
		call error_abort('Incorrect value of sort_flag')
	end select

	!Work out number of cell-blocks required
	blocksidelength = sortblocksize*cellsidelength
	nblocks = ceiling(domain/blocksidelength)
	blocks = product(nblocks)	

	!Average number of molecules per block of sortblocksize cell
	ave_molperblock = np/blocks
	ave_molperblock = ceiling(ave_molperblock*sf1+sf2) 	!Add safety factor

	!All variables associated with molecules
	if (ensemble.eq.tag_move) then
		allocate(tagtemp( ave_molperblock,blocks))
		allocate(rtethertemp(nd,ave_molperblock,blocks))
		tagtemp=0; rtethertemp=0.d0
	endif
	allocate(rtemp(nd,ave_molperblock,blocks))
	allocate(vtemp(nd,ave_molperblock,blocks))
	rtemp=0.d0; vtemp=0.d0
	!if (rtrue_flag.eq.1) then
	!	allocate(rtruetemp(nd,ave_molperblock,blocks))
	!	rtruetemp=0.d0
	!endif
	allocate(molperblock(blocks)); molperblock = 0

	!Copy all molecules to temp arrays in order of blocks
	do n = 1,np
		iblk = ceiling((r(1,n)+halfdomain(1))/blocksidelength(1)) 
		jblk = ceiling((r(2,n)+halfdomain(2))/blocksidelength(2)) 
		kblk = ceiling((r(3,n)+halfdomain(3))/blocksidelength(3)) 

		!i = iblk + nblocks(1)*(jblk-1) + nblocks(1)*nblocks(2)*(kblk-1)
		i = Hcurve(iblk,jblk,kblk)
		molperblock(i) = molperblock(i) + 1
		! Error handeling
		if (molperblock(i) .gt. ave_molperblock) then
			sf1 = sf1 + 0.1d0; sf2 = sf2 + 5;
			!print'(a,i6,a,3i4,a,i4,2(a,i6))', 'Sort_mols Warning -- There are ',molperblock(i), ' Mols in block  ', & 
			!				 iblk,jblk,kblk, ' on proc ', irank, & 
			!				' which is greater than the expected', ave_molperblock,' increased to ',nint(ave_molperblock*sf1+sf2)
			return	!Miss this round of sorting 
		endif
		!print'(7i8)', n,iblk,jblk,kblk,i, molperblock(i), ave_molperblock
		if (ensemble.eq.tag_move) then
			tagtemp(molperblock(i),i) = tag(n)
			rtethertemp(:,molperblock(i),i) = rtether(:,n)
		endif
		rtemp(:,molperblock(i),i) = r(:,n)
		vtemp(:,molperblock(i),i) = v(:,n)
		!if (rtrue_flag.eq.1) then
		!	rtruetemp(:,n) = rtrue(:,n)
		!endif

		!print*, 'b4',n, r(:,n)
	enddo

	!Copy back sorted molecles
	start = 1
	do i=1,blocks
		finish = start + molperblock(i)-1
		!print*, i, molperblock(i), start, finish!, rtemp(:,i,start:finish)
		if (ensemble.eq.tag_move) then
			tag(start:finish) = tagtemp(1:molperblock(i),i)
			rtether(:,start:finish) = rtethertemp(:,1:molperblock(i),i)
		endif
		r(:,start:finish) = rtemp(:,1:molperblock(i),i)
		v(:,start:finish) = vtemp(:,1:molperblock(i),i)

		!Read molecular tag and assign correct properties to reordered molecules
		if (ensemble.eq.tag_move) then
			do n=start,finish
				call read_tag(n)
			enddo
		endif
		start = start + molperblock(i)
	enddo

end subroutine sort_mols



!===================================================================================
!Deallocate all linklists

subroutine linklist_deallocateall
    use computational_constants_MD, only : potential_flag
	use polymer_info_MD, only: bond, bondcount
    use linked_list, only : linklist_deallocate_bins, & 
                            linklist_deallocate_cells, &
                            linklist_deallocate_neighbour, &
                            bin, cell, neighbour
	implicit none

	!integer            			:: i, j
	!integer           			:: cellnp, noneighbrs, np_neigbrs
	!integer            			:: icell, jcell, kcell
	!type(node), pointer 	    :: old, current, oldn, currentn

    !Reset all polymer information
	select case(potential_flag)
	case(0)
	case(1)
		bond = 0
		bondcount = 0
	end select

!	do icell=1,ncells(1)+2
!	do jcell=1,ncells(2)+2
!	do kcell=1,ncells(3)+2

!		if (associated(cell%head(icell,jcell,kcell)%point) .eqv. .true. ) then !Exit if null
!			old => cell%head(icell,jcell,kcell)%point
!			cellnp = cell%cellnp(icell,jcell, kcell)
!			current => old ! make current point to head of list
!			do j=1,cellnp-1
!				if (associated(old%next) .eqv. .true. ) then !Exit if null
!					old => current%next       !Make list point to next node of old
!					nullify(current%next)     !Remove pointer to next
!					nullify(current%previous) !Remove pointer to previous
!					deallocate(current)       !Deallocate current entry
!					current => old            !Make current point to new previous
!				endif
!			enddo

!			nullify(old%next)     !Remove pointer to next
!			nullify(old%previous) !Remove pointer to previous
!			deallocate(old)       !Deallocate final entry
!			nullify(cell%head(icell,jcell,kcell)%point) !Set cell head pointer to null
!			cell%cellnp(icell,jcell, kcell) = 0         !Zero cell molecule number
!		endif
!	enddo
!	enddo
!	enddo
    call linklist_deallocate_bins(bin)
    call linklist_deallocate_cells(cell)
    call linklist_deallocate_neighbour(neighbour)


!	if (force_list .gt. 1) then
!		if (size(neighbour%Nlist).eq.np) then
!		!if (vflux_outflag.ne.4 ) then
!			np_neigbrs = np
!		else
!			np_neigbrs = np + halo_np
!		endif
!		do i = 1, np_neigbrs
!			if (associated(neighbour%head(i)%point) .eqv. .true. ) then !Exit if null
!	       		noneighbrs = neighbour%Nlist(i)  !Determine number of elements in neighbourlist
!				oldn => neighbour%head(i)%point		   !Set old to head of neighbour list
!				currentn => oldn ! make current point to head of list
!				do j=1,noneighbrs-1
!					if (associated(oldn%next) .eqv. .true. ) then !Exit if null
!						oldn => currentn%next       !Make list point to next node of old
!						nullify(currentn%next)      !Remove pointer to next
!						nullify(currentn%previous)  !Remove pointer to previous
!						deallocate(currentn)        !Deallocate current entry
!						currentn => oldn            !Make current point to new previous
!					endif
!				enddo

!				nullify(oldn%next)         !Remove pointer to next
!				nullify(oldn%previous)     !Remove pointer to previous
!				deallocate(oldn)           !Deallocate final entry
!				nullify(neighbour%head(i)%point)   !Set neighbour head pointer to null
!				neighbour%Nlist(i) = 0  !Zero cell molecule number
!			endif
!		enddo

!	    !Deallocate array of molecules neighbourlist pointers
!	    deallocate(neighbour%Nlist)
!	    deallocate(neighbour%head)
!	endif



end subroutine linklist_deallocateall

