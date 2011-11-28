!=======================================================================
! Message-passing routines for parallel processing using MPI
!
!===============================
!  Key Messenger Subroutines   =
!===============================
! messenger_invoke()
! messenger_init()
! messenger_syncall()
! messenger_free()
!===============================
!  Border Update Subroutines   =
!===============================
! messenger_updateborders()
! updateface(ixyz)
! updateedge(face1,face2)
! updatecorners()
!===============================
! Molecule Sending Subroutines =
!===============================
! sendmols()
! sendface(ixyz)
! sendedge(face1,face2)
! sendcorners()
!!===============================
! Send-Probe-Receive Subroutines=
!================================
! sendproberecv(recvsize,sendsize,sendbuffer,pos,isource,idest)
! pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,isource,idest,pairsplit)
! NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
!===============================
! Data gathering subroutines   =
!===============================
! globalsyncreduce(A, na, meanA, maxA, minA)
! globalSum(A)
! globalSumInt(A)
! globalMax(A)
! globalMin(A)
! globalSumVect(A,na)
! globalSumIntVect(A,na)
! globalAverage(A, na)
!

module messenger
	use mpi
	use computational_constants_MD
        use coupler

        integer MD_COMM                      ! global communicator
	integer myid                         ! my process rank
	integer idroot                       ! rank of root process

	! Grid topology
	integer icomm_grid                   ! comm for grid topology
	integer, allocatable :: icoord(:,:)  ! proc grid coordinates
	integer	icomm_xyz(3) 		     ! Directional subcomms
	integer, dimension(8,2) 	:: proc_topology_corner
	integer, dimension(4,3,2) 	:: proc_topology_edge

	logical :: Lperiodic(3)

	double precision wallTime

end module

!======================================================================
!			Key Messenger Subroutines                     =
!======================================================================
subroutine messenger_invoke()
	use messenger
 

        call MPI_init(ierr)
        
        if (coupler_is_active) then
                call coupler_create_comm(COUPLER_MD,MD_COMM,ierr)
                prefix_dir = "./md_data/"
        else
                MD_COMM = MPI_COMM_WORLD 
        endif
         
end


subroutine messenger_init()
	use messenger
	use physical_constants_MD
	use librarymod, only : locate
	implicit none
	!include "mpif.h"

	integer idims(nd)
        integer  ndims, ip, ixyz, iper(3)
	logical  Lremain_dims(nd)

	! Initialize MPI
	call MPI_comm_size (MD_COMM, nproc, ierr)
	call MPI_comm_rank (MD_COMM, myid, ierr)

        ! Get the periodic constrains form MD.in
        ! and the processor topology description
        open(1,file=trim(prefix_dir)//'MD.in')

        call locate(1,'PERIODIC',.true.)
	read(1,*) iper(1)
	read(1,*) iper(2)
	if (nd.eq.3) read(1,*) iper(3)

        call locate(1,'PROCESSORS',.true.)
        read(1,*) npx
        read(1,*) npy
        read(1,*) npz

	close(1,status='keep')      !Close input file

        ! set Lperiodic
        Lperiodic(1:nd) = .true.
        where(iper(1:nd) == 0) Lperiodic(1:nd) =.false.  
        write(0,*) 'Lperiodic ', Lperiodic

	! Grid topology

        ! if npz == 0 in MD.in it means that npz = nproc/(npx*npy)
        if (npz == 0) then 
                npz = nproc/(npx*npy)
        endif

        !check if npx*npy*npz=nproc
        if (npx * npy * npz /= nproc ) then
                write(*,*) ' Wrong  specification for processor topology, nproc /= npx*npy*npz'
                call MPI_Abort(MPI_COMM_WORLD,1,ierr)
        endif

        ! allocate arrays that depend on topology parameters

        allocate(ibmin(npx), ibmax(npx), ibmino(npx), ibmaxo(npx), &
		jbmin(npy), jbmax(npy), jbmino(npy), jbmaxo(npy), &
		kbmin(npz), kbmax(npz), kbmino(npz), kbmaxo(npz), stat=ierr)
        if (ierr /= 0) then 
                write(*,*) "Error allocating topology arrays in messenger_init"
                call MPI_Abort(MD_COMM,2,ierr)
        endif

        allocate(icoord(3,nproc),stat=ierr)
        if (ierr /= 0) then 
                write(*,*) "Error allocating icoord in messenger_init"
                call MPI_Abort(MD_COMM,2,ierr)
        endif

	ndims = nd
	idims(1) = npx
	idims(2) = npy
	idims(3) = npz

	call MPI_Cart_create(MD_COMM, ndims, idims, Lperiodic, .true., &
	                     icomm_grid, ierr)

	do ip=1,nproc
		call MPI_Cart_coords(icomm_grid, ip-1, ndims, icoord(1,ip), ierr)
	end do
	icoord = icoord + 1
	call MPI_comm_rank (icomm_grid, irank, ierr)
	irank  = irank + 1
	iblock = icoord(1, irank)
	jblock = icoord(2, irank)
	kblock = icoord(3, irank)

	! Directional subcomms
	do ixyz=1,3
		Lremain_dims(:) = .false.
		Lremain_dims(ixyz) = .true.
		call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm_xyz(ixyz), ierr)
	end do

	call MPI_comm_rank (icomm_xyz(1), irankx, ierr)
	call MPI_comm_rank (icomm_xyz(2), iranky, ierr)
	call MPI_comm_rank (icomm_xyz(3), irankz, ierr)

	!print *, " Old comm=",irank-1,"new x,y,z:",irankx,iranky,irankz

	! Root process at coordinates (0,0,0)
	idims = 0
	call MPI_Cart_rank(icomm_grid, idims, idroot, ierr)
	iroot = idroot + 1

	! Molecules per processor and i/o writing offset
	allocate(procnp(nproc))	

	! Save current time
	wallTime = MPI_wtime()

	return
end

subroutine messenger_proc_topology()
	use messenger
	implicit none
	!include "mpif.h"

	integer				:: i, ixyz,n, idest, isource
	integer, dimension(3)  		:: pcoords, pshiftcoords
	integer, dimension(8)   	:: icornercell, jcornercell, kcornercell
	integer, dimension(3,4)  	:: edge1, edge2

	!---Setup edge topology--

	!Set up all 12 edges
	edge1(1,:) = (/2, 2, ncells(2)+1, ncells(2)+1/)
	edge2(1,:) = (/2, ncells(3)+1, 2, ncells(3)+1/)
	edge1(2,:) = (/2, 2, ncells(1)+1, ncells(1)+1/)
	edge2(2,:) = (/2, ncells(3)+1, 2, ncells(3)+1/)
	edge1(3,:) = (/2, 2, ncells(1)+1, ncells(1)+1/)
	edge2(3,:) = (/2, ncells(2)+1, 2, ncells(2)+1/)

	!Obtain processor coordinates
	call MPI_Cart_coords (icomm_grid,myid,3,pcoords,ierr)

	do i = 1,4 !Counter for each of the 4 edges

		!Along the 4 x edges
		pshiftcoords(1)=pcoords(1)
		pshiftcoords(2)=pcoords(2)-sign(1,ncells(2)-edge1(1,i))
		pshiftcoords(3)=pcoords(3)-sign(1,ncells(3)-edge2(1,i))

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)

		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
                if (Lperiodic(2)) then
		        pshiftcoords(2)=modulo(pshiftcoords(2),npy)
                endif
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

		!Obtain destination processor id
                if (pshiftcoords(2) < 0 .or. pshiftcoords(2) >= npy )then
                        idest = MPI_PROC_NULL
                else
		        call MPI_Cart_rank(icomm_grid,pshiftcoords,idest,ierr)
                endif
		!Add to array
		proc_topology_edge(i,1,1) = idest

		!isource = - idest
		pshiftcoords(1)=pcoords(1)
		pshiftcoords(2)=pcoords(2)+sign(1,ncells(2)-edge1(1,i))
		pshiftcoords(3)=pcoords(3)+sign(1,ncells(3)-edge2(1,i))

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)
		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
                if (Lperiodic(2)) then
		        pshiftcoords(2)=modulo(pshiftcoords(2),npy)
                endif
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

		!Obtain source processor id
                if (pshiftcoords(2) < 0 .or. pshiftcoords(2) >= npy )then
                        isource = MPI_PROC_NULL
                else
		        call MPI_Cart_rank(icomm_grid,pshiftcoords,isource,ierr)
                endif
		!Add to array
		proc_topology_edge(i,1,2) = isource

		!Along the 4 y edges
		pshiftcoords(1)=pcoords(1)-sign(1,ncells(1)-edge1(2,i))
		pshiftcoords(2)=pcoords(2)
		pshiftcoords(3)=pcoords(3)-sign(1,ncells(3)-edge2(2,i))

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)
		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
                if (Lperiodic(2)) then
		        pshiftcoords(2)=modulo(pshiftcoords(2),npy)
                endif
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

		!Obtain destination processor id
                if (pshiftcoords(2) < 0 .or. pshiftcoords(2) >= npy )then
                        idest = MPI_PROC_NULL
                else
		        call MPI_Cart_rank(icomm_grid,pshiftcoords,idest,ierr)
                endif
		!Add to array
		proc_topology_edge(i,2,1) = idest

		!isource = - idest
		pshiftcoords(1)=pcoords(1)+sign(1,ncells(1)-edge1(2,i))
		pshiftcoords(2)=pcoords(2)
		pshiftcoords(3)=pcoords(3)+sign(1,ncells(3)-edge2(2,i))

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)
		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
		pshiftcoords(2)=modulo(pshiftcoords(2),npy)
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

		!Obtain source processor id
		call MPI_Cart_rank(icomm_grid,pshiftcoords,isource,ierr)
		proc_topology_edge(i,2,2) = isource

		!Along the 4 z edges
		pshiftcoords(1)=pcoords(1)-sign(1,ncells(1)-edge1(3,i))
		pshiftcoords(2)=pcoords(2)-sign(1,ncells(2)-edge2(3,i))
		pshiftcoords(3)=pcoords(3)

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)
		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
                if (Lperiodic(2)) then
		        pshiftcoords(2)=modulo(pshiftcoords(2),npy)
                endif
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

		!Obtain destination processor id
                if (pshiftcoords(2) < 0 .or. pshiftcoords(2) >= npy )then
                        idest = MPI_PROC_NULL
                else
		        call MPI_Cart_rank(icomm_grid,pshiftcoords,idest,ierr)
                endif
		!Add to array
		proc_topology_edge(i,3,1) = idest

		!isource = - idest
		pshiftcoords(1)=pcoords(1)+sign(1,ncells(1)-edge1(3,i))
		pshiftcoords(2)=pcoords(2)+sign(1,ncells(2)-edge2(3,i))
		pshiftcoords(3)=pcoords(3)

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)
		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
                if (Lperiodic(2)) then
		        pshiftcoords(2)=modulo(pshiftcoords(2),npy)
                endif
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

		!Obtain source processor id
                if (pshiftcoords(2) < 0 .or. pshiftcoords(2) >= npy )then
                        isource = MPI_PROC_NULL
                else
		        call MPI_Cart_rank(icomm_grid,pshiftcoords,isource,ierr)
                endif
		!Add to array
		proc_topology_edge(i,3,2) = isource

	enddo

	!---Setup corner topology--

 	icornercell = (/ 2, 2, 2, 2, ncells(1)+1, ncells(1)+1, ncells(1)+1, ncells(1)+1/)
 	jcornercell = (/ 2, 2, ncells(2)+1, ncells(2)+1, 2, 2, ncells(2)+1, ncells(2)+1/) 
 	kcornercell = (/ 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1/)

	!Obtain processor coordinates
	call MPI_Cart_coords (icomm_grid,myid,3,pcoords,ierr)

	!Halo Corner Cells
	do i = 1,8

		!Obtain rank of diagonal processor
		pshiftcoords(1)=pcoords(1)-sign(1,ncells(1)-icornercell(i))
		pshiftcoords(2)=pcoords(2)-sign(1,ncells(2)-jcornercell(i))
		pshiftcoords(3)=pcoords(3)-sign(1,ncells(3)-kcornercell(i))

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)
		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
                if (Lperiodic(2)) then 
             	        pshiftcoords(2)=modulo(pshiftcoords(2),npy)
                endif
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

                if (pshiftcoords(2) < 0 .or. pshiftcoords(2) >= npy) then
                        idest = MPI_PROC_NULL
                else
		        call MPI_Cart_rank(icomm_grid,pshiftcoords,idest,ierr)
                endif
		!Add to array
		proc_topology_corner(i,1) = idest

		!isource = - idest
		pshiftcoords(1)=pcoords(1)+sign(1,ncells(1)-icornercell(i))
		pshiftcoords(2)=pcoords(2)+sign(1,ncells(2)-jcornercell(i))
		pshiftcoords(3)=pcoords(3)+sign(1,ncells(3)-kcornercell(i))

		!Adjust for periodic boundaries (works using -ve numbers and periodic 
		!topology but causes problems with some mpi flags and compilers)
		pshiftcoords(1)=modulo(pshiftcoords(1),npx)
                if (Lperiodic(2)) then
          	        pshiftcoords(2)=modulo(pshiftcoords(2),npy)
                endif
		pshiftcoords(3)=modulo(pshiftcoords(3),npz)

                if (pshiftcoords(2) < 0 .or. pshiftcoords(2) >= npy) then
                        isource = MPI_PROC_NULL
                else
		        call MPI_Cart_rank(icomm_grid,pshiftcoords,isource,ierr)
                endif
		!Add to array
		proc_topology_corner(i,2) = isource

	enddo

	return
end

subroutine messenger_syncall()
	use messenger
	!include "mpif.h"

	call MPI_Barrier(MD_COMM,ierr)

	return
end

subroutine messenger_lasterrorcheck
	use messenger
	!include "mpif.h"

	integer resultlen
	character*12 err_buffer

	call MPI_Error_string(ierr,err_buffer,resultlen,ierr)
	print*, err_buffer

end subroutine messenger_lasterrorcheck


subroutine messenger_free()
	use messenger
	!include "mpif.h"

	! Report time used
	print "(a,f8.2)", "time: ", MPI_wtime() - wallTime

	! Finalize MPI
         call MPI_finalize (ierr)

	return
end

!======================================================================
!			Border Update Subroutines                     =
!======================================================================

subroutine messenger_updateborders(rebuild)
use messenger
implicit none

	integer				 	:: rebuild
	
	if (all(periodic.lt.2)) then
	 	call messenger_updateborders_quiescent(rebuild)
	else
		stop "CANNOT USE LEES EDWARDS IN PARALLEL (YET!!)"
	end if

end subroutine messenger_updateborders


subroutine messenger_updateborders_quiescent(rebuild)
	use messenger
	use physical_constants_MD
	use arrays_MD
	implicit none

	integer				 	:: rebuild
	
	halo_np = 0

	!Update faces of domain
	call updatefacedown(1)
        call updatefacedown(2)
	call updatefacedown(3)

	call updatefaceup(1)
        call updatefaceup(2)
	call updatefaceup(3)

	!Update edges of domain
	call updateedge(1,2)
	call updateedge(2,3)
	call updateedge(1,3)

	!Update Corners of domain
	call updatecorners

	if (rebuild.eq.1) call assign_to_halocell(np+1,np+halo_np)

	return
end

!-----------------------------------------------------------------
! 		      Send to lower neighbor 	        	 -
!-----------------------------------------------------------------

!Update face halo cells by passing to neighbours
subroutine updatefacedown(ixyz)
	use physical_constants_MD
	use messenger
	use arrays_MD
	use linked_list
	implicit none
	!include "mpif.h"

	integer :: i, n, ixyz
	integer :: icell,jcell,kcell
	integer :: molno,cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,datasize,buffsize
	integer :: isource,idest
	double precision, dimension(nd) :: rpack	!Temporary array used to pack position
	double precision, dimension(:), allocatable :: sendbuffer
	type(node), pointer 	        :: old, current

	!Obtain processor ID of lower neighbour
	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)

	!print*, irank, 'facedown', ixyz, 'values', isource+1, idest+1

	!Obtain amount of data to send
	sendnp = 0
	select case (ixyz)
	case (1)
		icell = 2
		do jcell=2,ncells(2)+1
		do kcell=2,ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			sendnp = sendnp + cellnp
		enddo
		enddo
	case (2)
		jcell = 2
		do icell=2, ncells(1)+1
		do kcell=2, ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			sendnp = sendnp + cellnp
		enddo
		enddo
	case (3)
		kcell = 2
		do icell=2, ncells(1)+1
		do jcell=2, ncells(2)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			sendnp = sendnp + cellnp
		enddo
		enddo
	case default
		stop "updateBorder: invalid value for ixyz"
	end select

	!Three dimensions of data for each molecules
	sendsize = nd*sendnp
	allocate(sendbuffer(sendsize))
	call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
	icomm_grid,buffsize,ierr)

	!Package data ready to send
	pos = 0
	select case (ixyz)
        case (1)
		icell = 2
		do jcell=2, ncells(2)+1
		do kcell=2, ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point

			do i = 1,cellnp    !Step through each molecule in list 
				molno = old%molno !Number of molecule
				rpack(:) = r(molno,:)	!Load into temp array
				call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
				current => old
				old => current%next
			enddo
		enddo
		enddo
       	case (2)
		jcell = 2
		do icell=2, ncells(1)+1
		do kcell=2, ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point

			do i = 1,cellnp    !Step through each molecule in list 
				molno = old%molno !Number of molecule
				rpack(:) = r(molno,:)	!Load into temp array
				call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
				current => old
				old => current%next
			enddo
		enddo
		enddo
        case (3)
		kcell = 2
		do icell=2, ncells(1)+1
		do jcell=2, ncells(2)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point

			do i = 1,cellnp    !Step through each molecule in list 
				molno = old%molno !Number of molecule
				rpack(:) = r(molno,:)	!Load into temp array
				call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
				current => old
				old => current%next
			enddo
		enddo
		enddo
        case default
                stop "updateBorder: invalid value for ixyz"
        end select

	!If processor is its own neighbour - no passing required
!	if (idest+1 .eq. irank) then
!		recvsize = sendsize
!		call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)
!		length = recvsize*datasize
!		allocate(recvbuffer(recvsize))
!		recvbuffer = sendbuffer
!	else
		!Send, probe for size and then receive data
		call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
!	endif

	!call pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,ixyz)
	!call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)

	!Three dimensions per molecule
	recvnp = recvsize/real(nd,kind(0.d0))

	pos = 0
	do n=halo_np+1,halo_np+recvnp
		call MPI_Unpack(recvbuffer,length,pos,rpack, &
		nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
		r(np+n,:) = rpack
	enddo

	!Correct positions in new processor to halo cells
	do n=halo_np+1,halo_np+recvnp
		r(np+n,ixyz) = r(np+n,ixyz) + domain(ixyz)
		!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(np+n,:)
	enddo

	!Update number of molecules in halo to include number recieved
	halo_np = halo_np + recvnp

	!call MPI_Barrier(icomm_grid,ierr)

	deallocate(recvbuffer)
	deallocate(sendbuffer)
	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end

!-----------------------------------------------------------------
! 		      Send to upper neighbor 	        	 -
!-----------------------------------------------------------------

!Update face halo cells by passing to neighbours
subroutine updatefaceup(ixyz)
	use physical_constants_MD
	use messenger
	use arrays_MD
	use linked_list
	implicit none
	!include "mpif.h"

	integer :: i, n, ixyz
	integer :: icell,jcell,kcell
	integer :: molno,cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,datasize,buffsize
	integer :: isource,idest
	double precision, dimension(nd) :: rpack	!Temporary array used to pack/unpack position
	double precision, dimension(:), allocatable :: sendbuffer
	type(node), pointer 	        :: old, current

	!Obtain processor ID of upper neighbour
        call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)

	!print*, irank, 'faceup', ixyz,'values', isource+1, idest+1

	!Obtain amount of data to send
	sendnp = 0
	select case (ixyz)
        case (1)
		icell = ncells(1)+1
		do jcell=2, ncells(2)+1
		do kcell=2, ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			sendnp = sendnp + cellnp
		enddo
		enddo
       	case (2)
		jcell = ncells(2)+1
		do icell=2, ncells(1)+1
		do kcell=2, ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			sendnp = sendnp + cellnp
		enddo
		enddo
        case (3)
		kcell = ncells(3)+1
		do icell=2, ncells(1)+1
		do jcell=2, ncells(2)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			sendnp = sendnp + cellnp
		enddo
		enddo
        case default
                stop "updateBorder: invalid value for ixyz"
        end select

	!Three dimensions of data for each molecules
	sendsize = nd*sendnp
	allocate(sendbuffer(sendsize))
	call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
	icomm_grid,buffsize,ierr)

	!Package data ready to send
	pos = 0
	select case (ixyz)
        case (1)
		icell = ncells(1)+1
		do jcell=2, ncells(2)+1
		do kcell=2, ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point

			do i = 1,cellnp    !Step through each molecule in list 
				molno = old%molno !Number of molecule
				rpack(:) = r(molno,:)	!Load into temp array
				call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
				current => old
				old => current%next
			enddo
		enddo
		enddo
       	case (2)
		jcell = ncells(2)+1
		do icell=2, ncells(1)+1
		do kcell=2, ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point

			do i = 1,cellnp    !Step through each molecule in list 
				molno = old%molno !Number of molecule
				rpack(:) = r(molno,:)	!Load into temp array
				call MPI_Pack(rpack(:),nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
				current => old
				old => current%next 
			enddo
		enddo
		enddo
        case (3)
		kcell = ncells(3)+1
		do icell=2, ncells(1)+1
		do jcell=2, ncells(2)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point

			do i = 1,cellnp    !Step through each molecule in list 
				molno = old%molno !Number of molecule
				rpack(:) = r(molno,:)	!Load into temp array
				call MPI_Pack(rpack(:),nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
				current => old
				old => current%next
			enddo
		enddo
		enddo
        case default
                stop "updateBorder: invalid value for ixyz"
        end select

	!If processor is its own neighbour - no passing required
	if (idest+1 .eq. irank) then
		recvsize = sendsize
		call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)
		length = recvsize*datasize
		allocate(recvbuffer(recvsize))
		recvbuffer = sendbuffer
	else
		!Send, probe for size and then receive data
		call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
	endif

	!Send, probe for size and then receive data
	!call pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,ixyz)
	!call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)

	!Three dimensions per molecule
	recvnp = recvsize/real(nd,kind(0.d0))

	pos = 0
	do n=halo_np+1,halo_np+recvnp
		call MPI_Unpack(recvbuffer,length,pos,rpack, &
		nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
		r(np+n,:) = rpack
	enddo

	!Correct positions in new processor to halo cells
	do n=halo_np+1,halo_np+recvnp
		r(np+n,ixyz) = r(np+n,ixyz) - domain(ixyz)
		!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(np+n,:)
	enddo

	!Update number of molecules in halo to include number recieved
	halo_np = halo_np + recvnp

	!call MPI_Barrier(icomm_grid,ierr)

	deallocate(recvbuffer)
	deallocate(sendbuffer)
	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return
end


!Update edge halo cells by passing to neighbours

subroutine updateedge(face1,face2)
	use messenger
	use physical_constants_MD
	use arrays_MD
	use linked_list
	implicit none
	!include "mpif.h"

	integer :: i, n, ixyz,face1,face2
	integer :: icell,jcell,kcell
	integer :: molno,cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,datasize,buffsize
	integer :: isource,idest
	integer, dimension(3,4)   :: edge1, edge2
	double precision, dimension(nd) :: rpack	!Temporary array used to pack position
	double precision, dimension(:), allocatable :: sendbuffer
	type(node), pointer 	        :: old, current

	!Set up all 12 edges
	edge1(1,:) = (/2, 2, ncells(2)+1, ncells(2)+1/)
	edge2(1,:) = (/2, ncells(3)+1, 2, ncells(3)+1/)
	edge1(2,:) = (/2, 2, ncells(1)+1, ncells(1)+1/)
	edge2(2,:) = (/2, ncells(3)+1, 2, ncells(3)+1/)
	edge1(3,:) = (/2, 2, ncells(1)+1, ncells(1)+1/)
	edge2(3,:) = (/2, ncells(2)+1, 2, ncells(2)+1/)

	ixyz = 6 - face1 - face2 !Determine coordinate along edge

	do i = 1,4 !Counter for each of the 4 edges

		!Obtain rank of diagonal processor
		idest = proc_topology_edge(i,ixyz,1)
		isource = proc_topology_edge(i,ixyz,2)

		!print*, irank, 'edge', i, 'values:',isource+1, idest+1

		!Obtain amount of data to send
		sendnp = 0
		select case (ixyz)
        	case (1)
			do icell = 2, ncells(1)+1 !Move along x-axis
				cellnp = cell%cellnp(icell,edge1(1,i),edge2(1,i))
				sendnp = sendnp + cellnp
			enddo
       		case (2)
			do jcell = 2, ncells(2)+1 !Move along y-axis
				cellnp = cell%cellnp(edge1(2,i),jcell,edge2(2,i))
				sendnp = sendnp + cellnp
			enddo
 		case (3)
			do kcell = 2, ncells(3)+1 !Move along z-axis
				cellnp = cell%cellnp(edge1(3,i),edge2(3,i),kcell)
				sendnp = sendnp + cellnp
			enddo
		case default
                	stop "updateBorder: invalid value for ixyz"
		end select

		!Three dimensions of data for each molecules
		sendsize = nd*sendnp
		allocate(sendbuffer(sendsize))
		call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
		icomm_grid,buffsize,ierr)

		!Package data ready to send
		pos = 0
		select case (ixyz)
        	case (1)
			do icell = 2, ncells(1)+1 !Move along x-axis
				cellnp = cell%cellnp(icell,edge1(1,i),edge2(1,i))
				old => cell%head(icell,edge1(1,i),edge2(1,i))%point !Set old to top of link list
				do n=1,cellnp
					molno = old%molno		    !Obtain molecule number
					rpack(:) = r(molno,:)	!Load into temp array
					call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
					sendbuffer,buffsize,pos,icomm_grid,ierr)
					current => old			    !Use current to move to next
					old => current%next 		    !Use pointer to obtain next item in list
				enddo
			enddo
       		case (2)
			do jcell = 2, ncells(2)+1 !Move along y-axis
				cellnp = cell%cellnp(edge1(2,i),jcell,edge2(2,i))
				old => cell%head(edge1(2,i),jcell,edge2(2,i))%point !Set old to top of link list
				do n=1,cellnp
					molno = old%molno		    !Obtain molecule number
					rpack(:) = r(molno,:)	!Load into temp array
					call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
					sendbuffer,buffsize,pos,icomm_grid,ierr)
					current => old			    !Use current to move to next
					old => current%next 		    !Use pointer to obtain next item in list
				enddo
			enddo

 		case (3)
			do kcell = 2, ncells(3)+1 !Move along z-axis
				cellnp = cell%cellnp(edge1(3,i),edge2(3,i),kcell)
				old => cell%head(edge1(3,i),edge2(3,i),kcell)%point !Set old to top of link list
				do n=1,cellnp
					molno = old%molno		    !Obtain molecule number
					rpack(:) = r(molno,:)	!Load into temp array
					call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
					sendbuffer,buffsize,pos,icomm_grid,ierr)
					current => old			    !Use current to move to next
					old => current%next 		    !Use pointer to obtain next item in list
				enddo
			enddo

		case default
                	stop "updateBorder: invalid value for ixyz"
		end select

		!If processor is its own neighbour - no passing required
		if (idest+1 .eq. irank) then
			recvsize = sendsize
			call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)
			length = recvsize*datasize
			allocate(recvbuffer(recvsize))
			recvbuffer = sendbuffer
		else
			!Send, probe for size and then receive data
			call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
		endif

		!call pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,face1)
		!call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)

		!Three dimensions per molecule
		recvnp = recvsize/real(nd,kind(0.d0))

		pos = 0
		do n=halo_np+1,halo_np+recvnp
			call MPI_Unpack(recvbuffer,length,pos,rpack, &
			nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
			r(np+n,:) = rpack
		enddo

		!Correct positions in new processor to halo cells
		select case (ixyz)
        	case (1)
			do n=halo_np+1,halo_np+recvnp 
				r(np+n,2) = r(np+n,2) &  !Move to other side of domain
				+ sign(1,ncells(2)-edge1(1,i))*domain(2)
				r(np+n,3) = r(np+n,3) &  !Move to other side of domain
				+ sign(1,ncells(3)-edge2(1,i))*domain(3)
				!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(np+n,:)
			enddo
       		case (2)
			do n=halo_np+1,halo_np+recvnp
				r(np+n,1) = r(np+n,1) &  !Move to other side of domain
				+ sign(1,ncells(1)-edge1(2,i))*domain(1)
				r(np+n,3) = r(np+n,3) &  !Move to other side of domain
				+ sign(1,ncells(3)-edge2(2,i))*domain(3)
				!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(np+n,:)
			enddo

 		case (3)
			do n=halo_np+1,halo_np+recvnp
				r(np+n,1) = r(np+n,1) &  !Move to other side of domain
				+ sign(1,ncells(1)-edge1(3,i))*domain(1)
				r(np+n,2) = r(np+n,2) &  !Move to other side of domain
				+ sign(1,ncells(2)-edge2(3,i))*domain(2)
				!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(np+n,:)
			enddo
		case default
                	stop "updateBorder: invalid value for ixyz"
		end select

		!Update number of molecules in halo to include number recieved
		halo_np = halo_np + recvnp

		!call MPI_Barrier(MD_COMM,ierr)

		deallocate(recvbuffer)
		deallocate(sendbuffer)
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return
end


subroutine updatecorners()
	use messenger
	use physical_constants_MD
	use arrays_MD
	use linked_list
	implicit none
	!include "mpif.h"

	integer :: i, n
	integer :: molno,cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,buffsize
	integer :: isource,idest
	integer, dimension(8)   :: icornercell, jcornercell, kcornercell
	double precision, dimension(nd) :: rpack	!Temporary array used to pack position
	double precision, dimension(:), allocatable :: sendbuffer
	type(node), pointer 	        :: old, current

 	icornercell = (/ 2, 2, 2, 2, ncells(1)+1, ncells(1)+1, ncells(1)+1, ncells(1)+1/)
 	jcornercell = (/ 2, 2, ncells(2)+1, ncells(2)+1, 2, 2, ncells(2)+1, ncells(2)+1/) 
 	kcornercell = (/ 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1/)

	!Halo Corner Cells
	do i = 1,8 !Counter for each of the 4 edges

		!Obtain rank of diagonal processor

		idest = proc_topology_corner(i,1)
		isource = proc_topology_corner(i,2)

		!print*, irank, 'corner',i,'values', isource+1, idest+1

		!Obtain amount of data to send
		sendnp = 0
		cellnp = cell%cellnp(icornercell(i),jcornercell(i),kcornercell(i))
		sendnp = sendnp + cellnp

		!Three dimensions of data for each molecules
		sendsize = nd*sendnp
		allocate(sendbuffer(sendsize))
		call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
		icomm_grid,buffsize,ierr)

		!Package data ready to send
		pos = 0
		cellnp = cell%cellnp(icornercell(i),jcornercell(i),kcornercell(i))
		old => cell%head(icornercell(i),jcornercell(i),kcornercell(i))%point !Set old to top of link list
		do n=1,cellnp
			molno = old%molno		    !Obtain molecule number
			rpack(:) = r(molno,:)	!Load into temp array
			call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
			sendbuffer,buffsize,pos,icomm_grid,ierr)
			current => old			    !Use current to move to next
			old => current%next 		    !Use pointer to obtain next item in list
		enddo

		!Send, probe for size and then receive data
		!call pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,1)
		call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)

		!Three dimensions per molecule
		recvnp = recvsize/real(nd,kind(0.d0))

		pos = 0
		do n=halo_np+1,halo_np+recvnp
			call MPI_Unpack(recvbuffer,length,pos,rpack, &
			nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
			r(np+n,:) = rpack
		enddo

		!Correct positions in new processor to halo cells
		do n=halo_np+1,halo_np+recvnp
			r(np+n,1) = r(np+n,1) &  !Move to other side of domain
			+ sign(1,ncells(1)-icornercell(i))*domain(1)
			r(np+n,2) = r(np+n,2) &  !Move to other side of domain
			+ sign(1,ncells(2)-jcornercell(i))*domain(2)
			r(np+n,3) = r(np+n,3) &  !Move to other side of domain
			+ sign(1,ncells(3)-kcornercell(i))*domain(3)
			!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(np+n,:)
		enddo

		!Update number of molecules in halo to include number recieved
		halo_np = halo_np + recvnp

		!call MPI_Barrier(MD_COMM,ierr)

		deallocate(sendbuffer)
		deallocate(recvbuffer)
	enddo

	!do n=np+1,np+halo_np
	!	if(irank .eq. iroot) print'(i5,3f10.5)', iter, r(n,:)
	!enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return
end

!======================================================================
!			Molecule Transfer Subroutines                 =
!======================================================================

subroutine sendmols()
	use messenger
	use physical_constants_MD
	use linked_list
	implicit none

	integer		:: i,ixyz,dir,maxnew_np,sendnp,new_np

	pass%sendnp = 0

	!Three loops to cover possibility of more than one transfer
	do i = 1,nd
		do ixyz = 1,nd
		do dir= -1,1,2

			new_np = 0
			!Check if molecules to send
			call checksendbuild(ixyz,sendnp,dir)
			!Send molecules
			call sendrecvface(ixyz,sendnp,new_np,dir)
			!Reorder recieved molecules
			call reorderdata(new_np)
			!Deallocated passed list
			call linklist_deallocatepasslist

		enddo
		enddo

		!If no new molecules have been passed or received
		!then no need to check for further transfers or 
		!reorder new molecules
		maxnew_np = new_np
		call globalMaxInt(maxnew_np)
		!if(irank .eq. iroot) print*,'transfer number',i,'maximum number passed on any processor',maxnew_np
		!if (maxnew_np .eq. 0 ) exit
	enddo

	return
end

!-----------------------------------------------------------------------
! Check number of molecules to send to neighbour 	        	 

subroutine checksendbuild(ixyz,sendnp,dir)
	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use linked_list
	implicit none

	integer :: n, ixyz,dir
	integer :: sendnp

	!Obtain amount of data to send
	sendnp = 0
	select case(dir)
	case(-1)
		do n = 1,np
			if(r(n,ixyz) < -halfdomain(ixyz)) then
				call linklist_checkpushmol(n,0,0,0)
				!print*, 'proc id', irank, 'passes mol', n
				sendnp = sendnp + 1
			endif
		enddo
	case(1)
		do n = 1,np
			if(r(n,ixyz) >= halfdomain(ixyz)) then
				call linklist_checkpushmol(n,0,0,0)
				!print*, 'proc id', irank, 'passes mol', n
				sendnp = sendnp + 1
			endif
		enddo
	end select

end

!-----------------------------------------------------------------------
!Update face halo cells by passing to neighbours
subroutine sendrecvface(ixyz,sendnp,new_np,dir)
	use physical_constants_MD
	use messenger
	use arrays_MD
	use linked_list
	implicit none
	!include "mpif.h"

	integer :: i, n, ixyz, dir
	integer :: molno,sendnp,sendsize,recvnp,recvsize
	integer :: new_np,pos,length,datasize,buffsize
	integer :: isource,idest
	double precision			    :: tagpack	!Temporay packing buffer
	double precision, dimension(nd)		    :: Xpack 	!Temporay packing buffer
	double precision, dimension(:), allocatable :: sendbuffer
	type(passnode), pointer :: old, current

	!Obtain processor ID of lower neighbour
	call MPI_Cart_shift(icomm_grid, ixyz-1, dir, isource, idest, ierr)

	!print'(a,i5,a,i5,a,i5,a,i5,a,i5)', ' proc ', irank, ' coordinate ' &
	!,ixyz, ' dir ',dir, ' passes to ', idest+1, ' receives from ', isource+1

	!One Tag, three position and three velocity components for each molecules
	sendsize = (6 + 1)*sendnp
	allocate(sendbuffer(sendsize))
	call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
				icomm_grid,buffsize,ierr)

	!Package data ready to send
	old => pass%head 
	current => old     !make current point to head of list
	pos = 0
	do i=1,sendnp

		molno = old%molno	!Number of molecule
		!Positions
		Xpack(:) = r(molno,:)
		call MPI_Pack(Xpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
		!Velocity
		Xpack(:) = v(molno,:)
		call MPI_Pack(Xpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
		!Molecular Tag - Convert to double precision so all passed variables same type
		tagpack = real(tag(molno),kind(0.d0)) 
		call MPI_Pack(tagpack,1,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)

		!if(irank .eq. 6) print*, 'proc', irank, 'dest', idest+1, 'sending', sendbuffer

		old => current%next  !make old point to next node of current
		current  => old      !Set current item to old ready for next loop

	enddo

	!If processor is its own neighbour - no passing required
	if (idest+1 .eq. irank) then
		recvsize = sendsize
		call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)
		length = recvsize*datasize
		allocate(recvbuffer(recvsize))
		recvbuffer = sendbuffer
	else
		!Send, probe for size and then receive data
		call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
	endif

	!Send, probe for size and then receive data
	!call pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,ixyz)
	!call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)

	!One tag, three spatial and Three velocity components per molecule
	recvnp = recvsize/(2.d0*nd + 1)

	!Unpack data into correct arrays
	pos = 0
	do n=new_np+1,new_np+recvnp
		call MPI_Unpack(recvbuffer,length,pos,Xpack, &
				nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
		r(np+n,:) = Xpack
		call MPI_Unpack(recvbuffer,length,pos,Xpack, &
				nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
		v(np+n,:) = Xpack
		call MPI_Unpack(recvbuffer,length,pos,tagpack, &
				1,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
		tag(np+n) = nint(tagpack)	!Convert back to integer

		!if(irank .eq. 6) print*, 'proc', irank, 'from', isource+1, 'receiving', recvbuffer
	enddo

	!Correct positions in new processor
	do n=new_np+1,new_np+recvnp
		r(np+n,ixyz) = r(np+n,ixyz) - dir * domain(ixyz)
	enddo

	!Update number of molecules in halo to include number recieved
	new_np = new_np + recvnp

	!print*, 'proc', irank, 'new molecules', new_np

	deallocate(recvbuffer)
	deallocate(sendbuffer)
	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end

!------------------------------------------------------------------------------------
!Move through molecule list of passed molecule and fill with new arrivals starting 
!from the end of the new molecules added to r this timestep

subroutine reorderdata(new_np)
	use computational_constants_MD
 	use physical_constants_MD
	use arrays_MD
	use linked_list
	implicit none

	integer			   :: i
	integer			   :: molno,sendnp,new_np
	type(passnode), pointer    :: old, current

	!Work through passed molecules to fill gaps
	old => pass%head
	sendnp = pass%sendnp
	current => old     !make current point to head of list

	!Loop through all empty locations
	do i=1,sendnp

		molno = old%molno	!Position of empty location 
		r(molno,:) = r(np+new_np,:)
		v(molno,:) = v(np+new_np,:)
		tag(molno) = tag(np+new_np)

		!Read molecular tag and assign correct properties to reordered molecules
		call read_tag(molno)

		!if (molno .eq. np+new_np)
		!print*, irank, 'mol no',np+new_np, 'into gap', molno, 'of', np

		!Update number of new molecules
		new_np = new_np - 1

		old => current%next  !make old point to next node of current
		current  => old      !Set current item to old ready for next loop
	enddo

	!Read molecular tag and assign correct properties to new molecules
	do i=np+1, np+new_np
		call read_tag(i)
	enddo
	!Adjust total number of molecules to reflect difference between
	!recieved molecules and passed molecules
	np = np + new_np

	nullify (current)   !Set current pointer to null
	nullify (old)       !Set old pointer to null

	!print*, irank, 'number of particles', np

end subroutine reorderdata

!======================================================================
!		ALTERNATE Molecule Transfer Subroutines               =
!======================================================================
!FASTER BUT MORE COMPLICATED - ONLY SENDS WHEN IT NEEDS TO!
!----------------------------------------------------------------------
!Send an array detailing molecules which will be sent and recieved by
!each of the other processors

subroutine arraysendrecv(ipass,jpass,kpass,passarray,recvarray)
	use messenger
	use physical_constants_MD
	use linked_list
	use arrays_MD
	implicit none
	!include "mpif.h"

	integer :: ipass, jpass, kpass
	integer :: status(MPI_STATUS_SIZE)
	integer :: isource,idest,icount
	integer :: buf1, buf2
	integer, dimension(nd,nd,nd) :: passarray, recvarray

	!Obtain processor coordinates
	call procneighbours(ipass,jpass,kpass,isource,idest)

	!Obtain number of molecules to receive on all processes
	icount = 1
	buf1 = passarray(2+ipass,2+jpass,2+kpass)

	call MPI_sendrecv(buf1, icount, MPI_integer, & 
		idest, 0, buf2, icount, MPI_integer, & 
		isource, 0, icomm_grid, status, ierr)

	recvarray(2+ipass,2+jpass,2+kpass) = buf2

	return
end

!------------------------------------------------------------------------------------
!Send molecule from current processor and receive on destination processor

subroutine molsendrecv(recvarray)
	use messenger
	use physical_constants_MD
	use linked_list
	use arrays_MD
	implicit none
	!include "mpif.h"

	integer :: i , n
	integer :: molno,norecvd,sendnp,ipass,jpass,kpass
	integer :: status(MPI_STATUS_SIZE),req
	integer :: isource,idest,icount
	double precision, dimension(6):: buf1, buf2
	integer,dimension(nd,nd,nd)   :: recvarray
	type(passnode), pointer    :: oldp, currentp

	icount = 6

	! --------------- SEND --------------------
	!Work through list and send molecules
	oldp => pass%head 
	sendnp = pass%sendnp
	currentp => oldp     !make current point to head of list

	do i=1,sendnp

		molno = oldp%molno	!Number of molecule
		ipass = oldp%ipass	!Direction to pass
		jpass = oldp%jpass	!Direction to pass
		kpass = oldp%kpass	!Direction to pass

		!Build Buffer to send
		buf1(1) = r(molno,1) - ipass*domain(1)
		buf1(2) = r(molno,2) - jpass*domain(2)
		buf1(3) = r(molno,3) - kpass*domain(3)
		buf1(4) = v(molno,1)
		buf1(5) = v(molno,2)
		buf1(6) = v(molno,3)

		!Obtain processors neighbours
		call procneighbours(ipass,jpass,kpass,isource,idest)

		!print*, irank, 'sending', r(molno,2), 'to', idest+1

		!Non blocking send
		call MPI_isend(buf1, icount, MPI_DOUBLE_PRECISION, &
		idest,0,icomm_grid,req,ierr)

		oldp => currentp%next  !make old point to next node of current
		currentp  => oldp      !Set current item to old ready for next loop

	enddo

	nullify (currentp)   !Set current pointer to null
	nullify (oldp)       !Set old pointer to null

	n = 0 !Set new molecule counter to zero

	! --------------- RECEIVE --------------------
	!Work through recvarray and receive molecules
	do ipass = -1,1
	do jpass = -1,1
	do kpass = -1,1

		if(recvarray(2+ipass,2+jpass,2+kpass) .ne. 0) then

			!Obtain processors neighbours
			call procneighbours(ipass,jpass,kpass,isource,idest)
	
			norecvd = recvarray(2+ipass,2+jpass,2+kpass)
			do i = 1,norecvd

				!Recieve molecules
				call MPI_Recv(buf2,icount,MPI_DOUBLE_PRECISION, &
				isource,0,icomm_grid,status,ierr) 
	
				!Increment new molecule counter by one
				n = n + 1

				!Store new molecule position and velocity
				r(np+n,1) = buf2(1)
				r(np+n,2) = buf2(2)
				r(np+n,3) = buf2(3)
				v(np+n,1) = buf2(4)
				v(np+n,2) = buf2(5)
				v(np+n,3) = buf2(6)

				!print*, irank, 'receiving', r(np+n,2), 'from', isource+1

				!Correct position on receiving processor
				!r(np+n,1) = r(np+n,1) + ipass*domain(1)
				!r(np+n,2) = r(np+n,2) + jpass*domain(2)
				!r(np+n,3) = r(np+n,3) + kpass*domain(3)

			enddo
		endif
	enddo
	enddo
	enddo

	return
end

!--------------------------------------------------------------------------
!Obtain processor's neighbours from ipass, jpass and kpass
subroutine procneighbours(ipass,jpass,kpass,isource,idest)
	use messenger
	implicit none
	!include "mpif.h"

	integer,intent(in) :: ipass, jpass, kpass
	integer,intent(out):: isource,idest
	integer, dimension(3)   :: pcoords, pshiftcoords

	!Obtain processor coordinates
	call MPI_Cart_coords (icomm_grid,myid,3,pcoords,ierr)

	pshiftcoords(1)=pcoords(1)-ipass
	pshiftcoords(2)=pcoords(2)-jpass
	pshiftcoords(3)=pcoords(3)-kpass
	call MPI_Cart_rank(icomm_grid,pshiftcoords,isource,ierr)
	pshiftcoords(1)=pcoords(1)+ipass
	pshiftcoords(2)=pcoords(2)+jpass
	pshiftcoords(3)=pcoords(3)+kpass
	call MPI_Cart_rank(icomm_grid,pshiftcoords,idest,ierr)

	return
end

!======================================================================
!	        	Send-Probe-Receive Subroutines                =
!======================================================================

!Send packaged data; test for size of package; receive data.  
!THIS ROUTINE RELIES ON BUFFERING SO MAY FAIL FOR LARGE DATA VOLUMES

subroutine sendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
	use messenger
	use physical_constants_MD
	use arrays_MD
	implicit none
	!include "mpif.h"

	integer :: status(MPI_STATUS_SIZE,2), datasize
	integer,intent(in) :: pos,isource,idest
	integer,intent(in) :: sendsize
	integer,intent(out):: recvsize,length
	double precision,intent(in) :: sendbuffer(sendsize)

	!Send data to neighbour
	call MPI_Send(sendbuffer, pos, MPI_PACKED, &
	idest, 0, icomm_grid, ierr)	

	!Find out how many particles are being recieved
	call MPI_probe(isource,0,icomm_grid,status,ierr)
	call MPI_Get_count(status,MPI_PACKED,length,ierr)
	call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)

	recvsize = length/datasize
	allocate(recvbuffer(recvsize))

	!receive particles
	call MPI_Recv(recvbuffer,length,MPI_PACKED, &
	isource,0,icomm_grid,status,ierr)

	return
end

!Send packaged data using paired send and receive to prevent dependance on buffering;
!test for size of package; receive data

subroutine pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,pairsplit)
	use messenger
	use physical_constants_MD
	use arrays_MD
	implicit none
	!include "mpif.h"

	integer :: pairsplit
	integer :: status(MPI_STATUS_SIZE,2), datasize
	integer,intent(in) :: pos,isource,idest
	integer,intent(in) :: sendsize
	integer,intent(out):: recvsize,length
	double precision,intent(in) :: sendbuffer(sendsize)
	integer,dimension(nd) :: coords

	call MPI_cart_coords(icomm_grid,myid,3,coords,ierr)

	if (mod(coords(pairsplit)+1,2).eq.0) then

		!Send data to neighbour
		call MPI_sSend(sendbuffer, pos, MPI_PACKED, &
		idest, 0, icomm_grid, ierr)	

		!Find out how many particles are being recieved
		call MPI_probe(isource,0,icomm_grid,status,ierr)
		call MPI_Get_count(status,MPI_PACKED,length,ierr)
		call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)

		recvsize = length/datasize
		allocate(recvbuffer(recvsize))

		!receive particles
		call MPI_Recv(recvbuffer,length,MPI_PACKED, &
		isource,0,icomm_grid,status,ierr)

	else

		!Find out how many particles are being recieved
		call MPI_probe(isource,0,icomm_grid,status,ierr)
		call MPI_Get_count(status,MPI_PACKED,length,ierr)
		call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)

		recvsize = length/datasize
		allocate(recvbuffer(recvsize))

		!receive particles
		call MPI_Recv(recvbuffer,length,MPI_PACKED, &
		isource,0,icomm_grid,status,ierr)

		!Send data to neighbour
		call MPI_sSend(sendbuffer, pos, MPI_PACKED, &
		idest, 0, icomm_grid, ierr)

	endif

	return
end

!Non-blocking Send of packaged data; test for size of package; receive data 

subroutine NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
	use messenger
	use physical_constants_MD
	use arrays_MD
	implicit none
	!include "mpif.h"

	integer	:: datasize, request
	integer,intent(in) :: pos,isource,idest
	integer,intent(in) :: sendsize
	integer,intent(out):: recvsize, length
	integer, dimension(:), allocatable :: status
	double precision,intent(in) :: sendbuffer(sendsize)

	allocate(status(MPI_STATUS_SIZE))

	!Send data to neighbour
	call MPI_isend(sendbuffer, pos, MPI_PACKED, &
	idest,0,icomm_grid,request,ierr)

	!Find out how many particles are being recieved
	call MPI_probe(isource,0,icomm_grid,status(:),ierr)
	call MPI_Get_count(status(:),MPI_PACKED,length,ierr)
	call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)

	recvsize = length/datasize
	allocate(recvbuffer(recvsize))

	!It appears blocking probe works before send is complete
	!so wait for send to complete before recv called


	!Receive particles
	call MPI_Recv(recvbuffer,length,MPI_PACKED, &
	isource,0,icomm_grid,status(:),ierr)

        call MPI_wait(request, status(:), ierr)

	return
end

!======================================================================
!			Data gathering subroutines                    =
!======================================================================

subroutine globalbroadcast(A,na,broadprocid)
	use messenger
	implicit none

	integer			:: na, broadprocid
	double precision	:: A

	call MPI_BCAST(A,na,MPI_DOUBLE_PRECISION,broadprocid-1,MD_COMM,ierr)

	return
end

subroutine globalsyncreduce(A, na, meanA, maxA, minA)
	use messenger
	implicit none
	!include "mpif.h"

        integer na, nprocs
	double precision, intent(in)  :: A(na)
	double precision, intent(out) :: meanA(na), maxA(na), minA(na)
	double precision buf(na)

	call MPI_Reduce(A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_MAX, iroot-1, MD_COMM, ierr)

	maxA = buf

	call MPI_Reduce(A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_MIN, iroot-1, MD_COMM, ierr)

	minA = buf

	call MPI_Reduce(A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, iroot-1, MD_COMM, ierr)
	call MPI_comm_size (MD_COMM, nprocs, ierr)

	meanA = buf/nprocs

	return
end

subroutine globalSum(A)
	use messenger


	double precision :: A, buf

	call MPI_AllReduce (A, buf, 1, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalSumInt(A)
	use messenger

	integer :: A, buf

	call MPI_AllReduce (A, buf, 1, MPI_INTEGER, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalMaxInt(A)
	use messenger

	integer :: A, buf

	call MPI_AllReduce (A, buf, 1, MPI_INTEGER, &
	                    MPI_MAX, MD_COMM, ierr)
	A = buf

	return
end


subroutine globalSumVect(A, na)
	use messenger
	!include "mpif.h"

        integer, intent(in) :: na
	double precision A(na)
	double precision buf(na)

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalSumIntVect(A, na)
	use messenger

        integer, intent(in) :: na
	integer A(na)
	integer buf(na)

	call MPI_AllReduce (A, buf, na, MPI_INTEGER, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalMaxVect(A, na)
	use messenger

        integer, intent(in) :: na
	double precision A(na)
	double precision buf(na)

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_MAX, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalMaxIntVect(A, na)
	use messenger

        integer, intent(in) :: na
	integer A(na)
	integer buf(na)

	call MPI_AllReduce (A, buf, na, MPI_INTEGER, &
	                    MPI_MAX, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalMinVect(A, na)
	use messenger

        integer, intent(in) :: na
	double precision A(na)
	double precision buf(na)

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_MIN, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalAverage(A, na)
	use messenger

        integer, intent(in) :: na
	double precision A(na)
	double precision buf(na)

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, MD_COMM, ierr)
	call MPI_comm_size (MD_COMM, np, ierr)

	buf = buf / np
	
	A = buf

	return
end

subroutine globalGathernp()
	use physical_constants_MD
	use messenger
	!include "mpif.h"

	call MPI_Allgather (np, 1, MPI_INTEGER, procnp, 1, &
			    MPI_INTEGER,icomm_grid, ierr)

	return
end


!----Sum routines over global sub communitcators

subroutine SubcommSumInt(A, ixyz)
	use messenger

        integer, intent(in) :: ixyz !Direction of sub-comm
	integer	A
	integer buf

	call MPI_AllReduce (A, buf, 1, MPI_INTEGER, &
	                    MPI_SUM, icomm_xyz(ixyz), ierr)
	A = buf

	return
end

subroutine SubcommSumVect(A, na, ixyz)
	use messenger

        integer, intent(in) :: na, ixyz !Direction of sub-comm
	double precision A(na)
	double precision buf(na)

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, icomm_xyz(ixyz), ierr)
	A = buf

	return
end

subroutine SubcommSumIntVect(A, na, ixyz)
	use messenger

        integer, intent(in) :: na, ixyz !Direction of sub-comm
	integer	A(na)
	integer buf(na)

	call MPI_AllReduce (A, buf, na, MPI_INTEGER, &
	                    MPI_SUM, icomm_xyz(ixyz), ierr)
	A = buf

	return
end

subroutine error_abort(msg)
        use mpi
        implicit none
        
        character(len=*), intent(in), optional :: msg

        integer errcode,ierr

        if (present(msg)) then 
                write(*,*) msg
        endif

        call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)

end subroutine error_abort
