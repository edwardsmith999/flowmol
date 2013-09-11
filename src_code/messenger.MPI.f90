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

    integer	::  MD_COMM                      ! global communicator
	integer	::  myid                         ! my process rank
	integer	::  idroot                       ! rank of root process

	! Grid topology
	integer 					:: icomm_grid	! comm for grid topology
	integer, allocatable 		:: icoord(:,:)  ! proc grid coordinates
	integer						:: icomm_xyz(3)	! Directional "line" subcomms
	integer						:: plane_comm(3)! Directional "plane" subcomms
	integer, dimension(8,2) 	:: proc_topology_corners
	integer, dimension(4,3,2) 	:: proc_topology_edge

	integer :: planerankx, planeranky, planerankz


	logical :: Lperiodic(3)

	double precision wallTime

contains

	!=============================================================================
	! Get molecule's global position from position local to processor.
	!-----------------------------------------------------------------------------
	function globalise(rloc) result(rglob)
		implicit none
		
		real(kind(0.d0)), intent(in)  :: rloc(3)
		real(kind(0.d0))              :: rglob(3)

		rglob(1) = rloc(1)-halfdomain(1)*(npx-1)+domain(1)*(iblock-1)
		rglob(2) = rloc(2)-halfdomain(2)*(npy-1)+domain(2)*(jblock-1)
		rglob(3) = rloc(3)-halfdomain(3)*(npz-1)+domain(3)*(kblock-1)

	end function globalise

	function localise(rglob) result(rloc)
		implicit none
		
		real(kind(0.d0)), intent(in)  :: rglob(3)
		real(kind(0.d0))              :: rloc(3)

		rloc(1) = rglob(1)+(halfdomain(1)*(npx-1))-domain(1)*(iblock-1)
		rloc(2) = rglob(2)+(halfdomain(2)*(npy-1))-domain(2)*(jblock-1)
		rloc(3) = rglob(3)+(halfdomain(3)*(npz-1))-domain(3)*(kblock-1)

	end function localise 

end module messenger

!======================================================================
!					Key Messenger Subroutines                         =
!======================================================================
subroutine messenger_invoke()
	use messenger

	!Initialise MPI
	call MPI_init(ierr)

#if (USE_COUPLER == 0)
    MD_COMM = MPI_COMM_WORLD 
#endif
        
end subroutine messenger_invoke


subroutine messenger_init()
	use messenger
	use physical_constants_MD
	use librarymod, only : locate, find3factors
	use interfaces, only : error_abort
	implicit none
	!include "mpif.h"

	logical					:: found_in_input
	integer 				:: ndims, ip, ixyz
	integer,dimension(3)	:: idims
	logical,dimension(3)    :: Lremain_dims

	! Initialize MPI
	call MPI_comm_size (MD_COMM, nproc, ierr)
	call MPI_comm_rank (MD_COMM, myid, ierr)

    ! Get the periodic constrains form MD.in
    ! and the processor topology description
    open(1,file=trim(input_file))

	call locate(1,'PERIODIC',.true.)
	read(1,*) periodic(1)
	read(1,*) periodic(2)
	read(1,*) periodic(3)

    call locate(1,'PROCESSORS',.false.,found_in_input)
	if (found_in_input) then
	    read(1,*) npx
   		read(1,*) npy
	    read(1,*) npz
	    !check if npx*npy*npz=nproc
	    if (npx * npy * npz .ne. nproc ) then
			print*, npx, npy, npz , nproc
			call error_abort(' Wrong specification for processor topology, nproc not equal to npx*npy*npz')
	    endif
	else
		!Assign (using prime factors) to each dimension if not specified
		call find3factors(nproc,npx,npy,npz)
		print*, 'WARNING - Number of processors not specified - Arbitrarily assigned as follows:'
		print*, 'npx = ', npx, 'npy = ', npy, 'npz = ', npz

	endif

	close(1,status='keep')      !Close input file

    ! set Lperiodic
    Lperiodic(1:nd) = .true.
    where(periodic(1:nd) .eq. 0) Lperiodic(1:nd) = .false.
 
	! ==== Grid topology ====
    allocate(icoord(3,nproc),stat=ierr)
    if (ierr .ne. 0) then 
		call error_abort('Error allocating icoord in messenger_init')
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

	! Directional line subcomms
	do ixyz=1,3
		Lremain_dims(:) = .false.
		Lremain_dims(ixyz) = .true.
		call MPI_Cart_sub(icomm_grid, Lremain_dims, icomm_xyz(ixyz), ierr)
	end do
	if (npx .lt. 2) icomm_xyz(1) = MPI_COMM_SELF
	if (npy .lt. 2) icomm_xyz(2) = MPI_COMM_SELF
	if (npz .lt. 2) icomm_xyz(3) = MPI_COMM_SELF
	!if (any(icomm_xyz .eq. MPI_COMM_NULL)) icomm_xyz(:)=MPI_COMM_SELF

	call MPI_comm_rank (icomm_xyz(1), irankx, ierr)
	call MPI_comm_rank (icomm_xyz(2), iranky, ierr)
	call MPI_comm_rank (icomm_xyz(3), irankz, ierr)

	! Directional plane subcomms
	Lremain_dims = (/.false.,.true.,.true./)
	call MPI_Cart_sub(icomm_grid,Lremain_dims,plane_comm(1),ierr)
	Lremain_dims = (/.true.,.false.,.true./)
	call MPI_Cart_sub(icomm_grid,Lremain_dims,plane_comm(2),ierr)
	Lremain_dims = (/.true.,.true.,.false./)
	call MPI_Cart_sub(icomm_grid,Lremain_dims,plane_comm(3),ierr)
	
	call MPI_comm_rank (plane_comm(1), planerankx, ierr)
	call MPI_comm_rank (plane_comm(2), planeranky, ierr)
	call MPI_comm_rank (plane_comm(3), planerankz, ierr)

!	call MPI_comm_size (plane_comm(1), planenprocx, ierr)
!	call MPI_comm_size (plane_comm(2), planenprocy, ierr)
!	call MPI_comm_size (plane_comm(3), planenprocz, ierr)

	! Root process at coordinates (0,0,0)
	idims = 0
	call MPI_Cart_rank(icomm_grid, idims, idroot, ierr)
	iroot = idroot + 1

	! Molecules per processor for i/o writing offset
	allocate(procnp(nproc))	
	! Molecules per processor for i/o writing offset
	allocate(proctethernp(nproc))	

	! Save current time
	wallTime = MPI_wtime()

end subroutine messenger_init

subroutine messenger_proc_topology()
	use messenger
	implicit none
	!include "mpif.h"

	integer						:: i, idest, isource
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
		if (Lperiodic(2)) pshiftcoords(2)=modulo(pshiftcoords(2),npy)
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
		proc_topology_corners(i,1) = idest

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
		proc_topology_corners(i,2) = isource

	enddo

end subroutine messenger_proc_topology

subroutine messenger_syncall()
	use messenger
	!include "mpif.h"

	call MPI_Barrier(MD_COMM,ierr)

	return
end subroutine messenger_syncall

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

	call MPI_Barrier(MD_COMM,ierr)

	! Report time used
	if (irank.eq.iroot) then
		print "(1x,a,f8.2,a)", "MPI run-time: ", MPI_wtime() - wallTime, ' seconds.'
	end if

	! Finalize MPI
    call MPI_finalize (ierr)

end subroutine messenger_free

!======================================================================
!			Border Update Subroutines                     =
!======================================================================

subroutine messenger_updateborders(rebuild)
	use interfaces
	use messenger
	use arrays_MD
	use physical_constants_MD, only: halo_np, np
	implicit none

	integer				 	:: rebuild
	
	if (all(periodic.ge.2)) then
		call error_abort( "CANNOT USE LEES EDWARDS IN PARALLEL")
	end if

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

end subroutine messenger_updateborders


!=======================================================================
! 	Routine to pack/unpack halo variables per cells using MPI commands =
!=======================================================================

module pack_unpack_cell
	use physical_constants_MD, only : nd
	use computational_constants_MD, only :	potential_flag
	use arrays_MD, only : r, v
	use polymer_info_MD, only: nsdmi
	use linked_list
	use messenger
	implicit none

contains

!Wrapper for the routine to Pack cells using MPI_pack
subroutine pack_cell(icell,jcell,kcell,sendbuffer,buffsize,pos)

	integer, intent(in)								:: icell,jcell,kcell,buffsize
	integer, intent(inout)							:: pos
	double precision, dimension(:), intent(out) 	:: sendbuffer

	integer 										:: i, molno,cellnp
	double precision, dimension(nd) 				:: rpack, vpack	!Temporary arrays used to pack
	double precision, dimension(nsdmi) 				:: FENEpack

	type(node), pointer 	        				:: old, current

	cellnp = cell%cellnp(icell,jcell,kcell)
	old => cell%head(icell,jcell,kcell)%point

	!print*, associated(cell%head(icell,jcell, kcell)%point),cellnp
	!call linklist_print(icell, jcell, kcell)

	do i = 1,cellnp    !Step through each molecule in list 
		molno = old%molno !Number of molecule
		select case (potential_flag)
		case(0)
			rpack(:) = r(:,molno)	!Load into temp array
			call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
			sendbuffer,buffsize,pos,icomm_grid,ierr)
			if (pass_vhalo .ne. 0) then
				vpack(:) = v(:,molno)	!Load into temp array
				call MPI_Pack(vpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
			endif
		case(1)
			rpack(:) = r(:,molno)	!Load into temp array
			call MPI_Pack(rpack,nd,MPI_DOUBLE_PRECISION,& 
			sendbuffer,buffsize,pos,icomm_grid,ierr)	
			if (pass_vhalo .ne. 0) then
				vpack(:) = v(:,molno)	!Load into temp array
				call MPI_Pack(vpack,nd,MPI_DOUBLE_PRECISION,& 
				sendbuffer,buffsize,pos,icomm_grid,ierr)
			endif
			call prepare_FENEbuffer(molno,FENEpack)
			call MPI_Pack(FENEpack,nsdmi,MPI_DOUBLE_PRECISION,&
			sendbuffer,buffsize,pos,icomm_grid,ierr)
		end select
		current => old
		old => current%next
	enddo

end subroutine pack_cell


!Wrapper for the routine to unpack cells using MPI_unpack
subroutine unpack_recvbuffer(halo_np,recvnp,length,recvbuffer)
	use physical_constants_MD, only : nd, np
	use computational_constants_MD, only :	potential_flag
	use arrays_MD, only : r, v
	use messenger

	integer, intent(in)								:: halo_np,recvnp,length
	double precision, dimension(:), intent(in) 		:: recvbuffer

	integer 										:: n, pos
	double precision, dimension(nd) 				:: rpack, vpack	!Temporary arrays used to pack
	double precision, dimension(nsdmi)  				:: FENEpack

	pos = 0
	do n=halo_np+1,halo_np+recvnp
		select case(potential_flag)
		case(0)
			call MPI_Unpack(recvbuffer,length,pos,rpack, &
			nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
			r(:,np+n) = rpack
			if (pass_vhalo .ne. 0) then
				call MPI_Unpack(recvbuffer,length,pos,vpack, &
				nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
				v(:,np+n) = vpack
			endif
		case(1)
			call MPI_Unpack(recvbuffer,length,pos,rpack, &
			nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
			r(:,np+n) = rpack
			if (pass_vhalo .ne. 0) then
				call MPI_Unpack(recvbuffer,length,pos,vpack, &
				nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
				v(:,np+n) = vpack
			endif
			call MPI_Unpack(recvbuffer,length,pos,FENEpack, &
			nsdmi,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
			call assign_FENEbuffer(np+n,FENEpack)
		end select
	enddo

end subroutine unpack_recvbuffer


!Get size of array to send from number of molecules
subroutine get_sendsize(sendnp,sendsize)
	use physical_constants_MD, only : nd
	use computational_constants_MD, only :	potential_flag

	integer, intent(in)		:: sendnp
	integer, intent(out)	:: sendsize

	select case (potential_flag)
	case(0) !LJ only
		sendsize = nd*sendnp
		if (pass_vhalo .ne. 0) sendsize = sendsize + nd*sendnp
	case(1) !+FENE info 
		sendsize = nd*sendnp + (nsdmi)*sendnp
		if (pass_vhalo .ne. 0) sendsize = sendsize + nd*sendnp
	end select

end subroutine get_sendsize

!Get number of molecules from size of recieved array

subroutine get_recvnp(recvsize,recvnp)
	use physical_constants_MD, only : nd
	use computational_constants_MD, only :	potential_flag

	integer, intent(in)		:: recvsize
	integer, intent(out)	:: recvnp

	select case(potential_flag)
	case(0)
		if (pass_vhalo .eq. 0) then
			recvnp = recvsize/real(nd,kind(0.d0))
		else
			recvnp = recvsize/real(2*nd,kind(0.d0))
		endif
	case(1)
		if (pass_vhalo .eq. 0) then
			recvnp = recvsize/real(nd+nsdmi,kind(0.d0))
		else
			recvnp = recvsize/real(2*nd+nsdmi,kind(0.d0))
		endif
	end select

end subroutine get_recvnp

end module pack_unpack_cell


!======================================================================
! Buffer preparation and unpacking routines for polymer simulation    =
!======================================================================

subroutine prepare_FENEbuffer(molno,FENEpack)
	use messenger
	use polymer_info_MD, only: monomer,nmonomers
	implicit none
	
	integer, intent(in) :: molno
	double precision, dimension(*), intent(out) :: FENEpack
		
	FENEpack(1)   = real(monomer(molno)%chainID,        kind(0.d0))
	FENEpack(2)   = real(monomer(molno)%subchainID,     kind(0.d0))
	FENEpack(3)   = real(monomer(molno)%funcy,          kind(0.d0))
	FENEpack(4)   = real(monomer(molno)%glob_no,        kind(0.d0))
	FENEpack(5:8) = real(monomer(molno)%bin_bflag(1:4), kind(0.d0))

end subroutine prepare_FENEbuffer

subroutine assign_FENEbuffer(molno,FENEpack)
	use messenger
	use polymer_info_MD, only: monomer,nmonomers
	implicit none

	integer, intent(in) :: molno
	double precision, dimension(*), intent(in) :: FENEpack
	
	monomer(molno)%chainID         = nint(FENEpack(1))	
	monomer(molno)%subchainID      = nint(FENEpack(2))	
	monomer(molno)%funcy           = nint(FENEpack(3))	
	monomer(molno)%glob_no         = nint(FENEpack(4))	
	monomer(molno)%bin_bflag(1:4)  = nint(FENEpack(5:8))

end subroutine assign_FENEbuffer

!-----------------------------------------------------------------
! 		      Send to lower neighbor 	        	 			 -
!-----------------------------------------------------------------

!Update face halo cells by passing to neighbours
subroutine updatefacedown(ixyz)
 	use interfaces
	use physical_constants_MD
	use polymer_info_MD
	use messenger
	use arrays_MD
	use linked_list
	use pack_unpack_cell
	implicit none
	!include "mpif.h"

	integer :: n, ixyz
	integer :: icell,jcell,kcell
	integer :: cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,datasize,buffsize
	integer :: isource,idest
	double precision, dimension(:), allocatable :: sendbuffer
	type(node), pointer :: old, current

	!Obtain processor ID of lower neighbour
	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)

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
		call error_abort( "updateBorder: invalid value for ixyz")
	end select


	!Determine size of send buffer
	call get_sendsize(sendnp,sendsize)

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
				call pack_cell(icell,jcell,kcell,sendbuffer,buffsize,pos)
			enddo
			enddo
       	case (2)
			jcell = 2
			do icell=2, ncells(1)+1
			do kcell=2, ncells(3)+1
				call pack_cell(icell,jcell,kcell,sendbuffer,buffsize,pos)
			enddo
			enddo
        case (3)
			kcell = 2
			do icell=2, ncells(1)+1
			do jcell=2, ncells(2)+1
				call pack_cell(icell,jcell,kcell,sendbuffer,buffsize,pos)
			enddo
			enddo
        case default
		call error_abort("updateBorder: invalid value for ixyz")
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

	!call pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,ixyz)
	!call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)

	!Get number of molecules from recieved data size
	call get_recvnp(recvsize,recvnp)

	!Unpack recieved halo data 
	call unpack_recvbuffer(halo_np,recvnp,length,recvbuffer)

	!Correct positions in new processor to halo cells
	do n=halo_np+1,halo_np+recvnp
		r(ixyz,np+n) = r(ixyz,np+n) + domain(ixyz)
		!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(:,np+n)
	enddo

	!Update number of molecules in halo to include number recieved
	halo_np = halo_np + recvnp

	!call MPI_Barrier(icomm_grid,ierr)

	deallocate(recvbuffer)
	deallocate(sendbuffer)
	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end subroutine updatefacedown

!-----------------------------------------------------------------
! 		      Send to upper neighbor 	        	 -
!-----------------------------------------------------------------

!Update face halo cells by passing to neighbours
subroutine updatefaceup(ixyz)
	use interfaces
	use physical_constants_MD
	use polymer_info_MD
	use messenger
	use arrays_MD
	use linked_list
	use pack_unpack_cell
	implicit none
	!include "mpif.h"

	integer :: n, ixyz
	integer :: icell,jcell,kcell
	integer :: cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,datasize,buffsize
	integer :: isource,idest
	double precision, dimension(:), allocatable :: sendbuffer
	type(node), pointer 	        :: old, current

	!Obtain processor ID of upper neighbour
	call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
	
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
		call error_abort("updateBorder: invalid value for ixyz")
        end select

	!Determine size of send buffer
	call get_sendsize(sendnp,sendsize)

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
				call pack_cell(icell,jcell,kcell,sendbuffer,buffsize,pos)
			enddo
			enddo
       	case (2)
			jcell = ncells(2)+1
			do icell=2, ncells(1)+1
			do kcell=2, ncells(3)+1
				call pack_cell(icell,jcell,kcell,sendbuffer,buffsize,pos)
			enddo
			enddo
        case (3)
			kcell = ncells(3)+1
			do icell=2, ncells(1)+1
			do jcell=2, ncells(2)+1
				call pack_cell(icell,jcell,kcell,sendbuffer,buffsize,pos)
			enddo
			enddo
        case default
		call error_abort("updateBorder: invalid value for ixyz")
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

	!Get number of molecules from recieved data size
	call get_recvnp(recvsize,recvnp)

	!Unpack recieved halo data 
	call unpack_recvbuffer(halo_np,recvnp,length,recvbuffer)

	!Correct positions in new processor to halo cells
	do n=halo_np+1,halo_np+recvnp
		r(ixyz,np+n) = r(ixyz,np+n) - domain(ixyz)
		!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(:,np+n)
	enddo

	!Update number of molecules in halo to include number recieved
	halo_np = halo_np + recvnp

	!call MPI_Barrier(icomm_grid,ierr)
	
	deallocate(recvbuffer)
	deallocate(sendbuffer)
	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return
end subroutine updatefaceup


!Update edge halo cells by passing to neighbours

subroutine updateedge(face1,face2)
	use interfaces
	use messenger
	use physical_constants_MD
	use polymer_info_MD
	use arrays_MD
	use linked_list
	use pack_unpack_cell
	implicit none
	!include "mpif.h"

	integer :: i, n, ixyz,face1,face2
	integer :: icell,jcell,kcell
	integer :: cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,datasize,buffsize
	integer :: isource,idest
	integer, dimension(3,4)   :: edge1, edge2
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
			call error_abort("updateBorder: invalid value for ixyz")
		end select

		!Determine size of send buffer
		call get_sendsize(sendnp,sendsize)

		allocate(sendbuffer(sendsize))
		call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
		icomm_grid,buffsize,ierr)

		!Package data ready to send
		pos = 0
		select case (ixyz)
    	case (1)
			do icell = 2, ncells(1)+1 !Move along x-axis
				call pack_cell(icell,edge1(1,i),edge2(1,i),sendbuffer,buffsize,pos)
			enddo
		case (2)
			do jcell = 2, ncells(2)+1 !Move along y-axis
				call pack_cell(edge1(2,i),jcell,edge2(2,i),sendbuffer,buffsize,pos)
			enddo
		case (3)
			do kcell = 2, ncells(3)+1 !Move along z-axis
				call pack_cell(edge1(3,i),edge2(3,i),kcell,sendbuffer,buffsize,pos)
			enddo
		case default
			call error_abort("updateBorder: invalid value for ixyz")
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
	
		!Get number of molecules from recieved data size
		call get_recvnp(recvsize,recvnp)

		!Unpack recieved halo data 
		call unpack_recvbuffer(halo_np,recvnp,length,recvbuffer)

		!Correct positions in new processor to halo cells
		select case (ixyz)
        	case (1)
				do n=halo_np+1,halo_np+recvnp 
					r(2,np+n) = r(2,np+n) &  !Move to other side of domain
					+ sign(1,ncells(2)-edge1(1,i))*domain(2)
					r(3,np+n) = r(3,np+n) &  !Move to other side of domain
					+ sign(1,ncells(3)-edge2(1,i))*domain(3)
				enddo
       		case (2)
				do n=halo_np+1,halo_np+recvnp
					r(1,np+n) = r(1,np+n) &  !Move to other side of domain
					+ sign(1,ncells(1)-edge1(2,i))*domain(1)
					r(3,np+n) = r(3,np+n) &  !Move to other side of domain
					+ sign(1,ncells(3)-edge2(2,i))*domain(3)
				enddo
 			case (3)
				do n=halo_np+1,halo_np+recvnp
					r(1,np+n) = r(1,np+n) &  !Move to other side of domain
					+ sign(1,ncells(1)-edge1(3,i))*domain(1)
					r(2,np+n) = r(2,np+n) &  !Move to other side of domain
					+ sign(1,ncells(2)-edge2(3,i))*domain(2)
					!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(:,np+n)
				enddo
		case default
			call error_abort("updateBorder: invalid value for ixyz")
		end select

		!Update number of molecules in halo to include number recieved
		halo_np = halo_np + recvnp

		deallocate(recvbuffer)
		deallocate(sendbuffer)

		
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end subroutine updateedge


subroutine updatecorners()
	use messenger
	use physical_constants_MD
	use polymer_info_MD
	use arrays_MD
	use linked_list
	use pack_unpack_cell
	implicit none
	!include "mpif.h"

	integer :: i, n
	integer :: cellnp,sendnp,sendsize,recvnp,recvsize,pos,length,buffsize
	integer :: isource,idest
	integer, dimension(8)   :: icornercell, jcornercell, kcornercell
	double precision, dimension(:), allocatable :: sendbuffer
	type(node), pointer 	        :: old, current


 	icornercell = (/ 2, 2, 2, 2, ncells(1)+1, ncells(1)+1, ncells(1)+1, ncells(1)+1/)
 	jcornercell = (/ 2, 2, ncells(2)+1, ncells(2)+1, 2, 2, ncells(2)+1, ncells(2)+1/) 
 	kcornercell = (/ 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1/)

	!Halo Corner Cells
	do i = 1,8 !Counter for each of the 4 edges

		!Obtain rank of diagonal processor

		idest = proc_topology_corners(i,1)
		isource = proc_topology_corners(i,2)

		!print*, irank, 'corner',i,'values', isource+1, idest+1

		!Obtain amount of data to send
		sendnp = 0
		cellnp = cell%cellnp(icornercell(i),jcornercell(i),kcornercell(i))
		sendnp = sendnp + cellnp

		!Determine size of send buffer
		call get_sendsize(sendnp,sendsize)

		allocate(sendbuffer(sendsize))
		call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
		icomm_grid,buffsize,ierr)

		!Package data ready to send
		pos = 0
		call pack_cell(icornercell(i),jcornercell(i),kcornercell(i),sendbuffer,buffsize,pos)

		!Send, probe for size and then receive data
		!call pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest,1)
		call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)

		!Get number of molecules from recieved data size
		call get_recvnp(recvsize,recvnp)

		!Unpack recieved halo data 
		call unpack_recvbuffer(halo_np,recvnp,length,recvbuffer)

		!Correct positions in new processor to halo cells
		do n=halo_np+1,halo_np+recvnp
			r(1,np+n) = r(1,np+n) &  !Move to other side of domain
			+ sign(1,ncells(1)-icornercell(i))*domain(1)
			r(2,np+n) = r(2,np+n) &  !Move to other side of domain
			+ sign(1,ncells(2)-jcornercell(i))*domain(2)
			r(3,np+n) = r(3,np+n) &  !Move to other side of domain
			+ sign(1,ncells(3)-kcornercell(i))*domain(3)
			!if(irank .eq. 18 .and. iter .gt. 33) print'(i5,3f10.5)', iter, r(:,np+n)
		enddo

		!Update number of molecules in halo to include number recieved
		halo_np = halo_np + recvnp

		!call MPI_Barrier(MD_COMM,ierr)

		deallocate(sendbuffer)
		deallocate(recvbuffer)
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine updatecorners

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

	integer 		:: n
	integer,intent(in)	:: ixyz,dir
	integer,intent(out) 	:: sendnp

	!Obtain amount of data to send
	sendnp = 0
	select case(dir)
	case(-1)
		do n = 1,np
			if(r(ixyz,n) < -halfdomain(ixyz)) then
				call linklist_checkpushmol(n,0,0,0)
				sendnp = sendnp + 1
			endif
		enddo
	case(1)
		do n = 1,np
			if(r(ixyz,n) >= halfdomain(ixyz)) then
				call linklist_checkpushmol(n,0,0,0)
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
	use polymer_info_MD
	use arrays_MD
	use linked_list
	use interfaces, only: error_abort
	implicit none

	integer, intent(inout) :: new_np
	integer, intent(in)    :: ixyz
	integer, intent(in)    :: sendnp
	integer, intent(in)    :: dir
	
	integer :: i, n
	integer :: molno,sendsize,recvnp,recvsize
	integer :: pos,length,datasize,buffsize
	integer :: isource,idest
	double precision :: dppack !Temporay packing buffer
	double precision, dimension(nd) :: Xpack !Temporay packing buffer
	double precision, dimension(nsdmi)  :: FENEpack
	double precision, dimension(:), allocatable :: sendbuffer
	type(passnode), pointer :: old, current
	character(len=1024) :: string

	!Obtain processor ID of lower neighbour
	call MPI_Cart_shift(icomm_grid, ixyz-1, dir, isource, idest, ierr)

	!Calculate required buffer size and allocate enough memory
	call get_sendmols_buffersize
	allocate(sendbuffer(sendsize))
	
	!Provide MPI with upper bound of pack size
	call MPI_Pack_size(sendsize,MPI_DOUBLE_PRECISION, &
	                   icomm_grid,buffsize,ierr)
	!Initialise position in packed buffer
	pos = 0

	!Pack buffer header (sendnp) first
	dppack = real(sendnp,kind(0.d0))
	call MPI_Pack(dppack,1,MPI_DOUBLE_PRECISION,sendbuffer,&
	              buffsize,pos,icomm_grid,ierr)

	!Pack rest of data -----------------------------------------!
	old => pass%head 
	current => old     !make current point to head of list
	do i=1,sendnp

		molno = old%molno!Number of molecule

		if (ensemble.eq.tag_move) then
			dppack = real(tag(molno),kind(0.d0))
			call MPI_Pack(dppack,1,MPI_DOUBLE_PRECISION, &
		              sendbuffer,buffsize,pos,icomm_grid,ierr)
		endif

		Xpack(:) = r(:,molno)
		call MPI_Pack(Xpack,nd,MPI_DOUBLE_PRECISION,& 
		              sendbuffer,buffsize,pos,icomm_grid,ierr)

		Xpack(:) = v(:,molno)
		call MPI_Pack(Xpack,nd,MPI_DOUBLE_PRECISION,& 
		              sendbuffer,buffsize,pos,icomm_grid,ierr)

		if (rtrue_flag .eq. 1) then
			Xpack(:) = rtrue(:,molno)
			call MPI_Pack(Xpack,nd,MPI_DOUBLE_PRECISION,& 
			              sendbuffer,buffsize,pos,icomm_grid,ierr)
		end if

		if (ensemble.eq.tag_move) then
			if (any(tag(molno).eq.tether_tags)) then
				Xpack(:) = rtether(:,molno)
				call MPI_Pack(Xpack,nd,MPI_DOUBLE_PRECISION,& 
			              sendbuffer,buffsize,pos,icomm_grid,ierr)
			end if
		endif

		if (potential_flag .eq. 1) then
			call prepare_FENEbuffer(molno,FENEpack)
			call MPI_Pack(FENEpack,nsdmi,MPI_DOUBLE_PRECISION,& 
			              sendbuffer,buffsize,pos,icomm_grid,ierr)
		end if	
		
		old => current%next  !make old point to next node of current
		current => old      !Set current item to old ready for next loop

	enddo
	!------------------------------------------------------------!

	!If processor is its own neighbour - no passing required
	if (idest+1 .eq. irank) then
		recvsize = sendsize
		call MPI_type_size(MPI_DOUBLE_PRECISION,datasize,ierr)
		length = recvsize*datasize
		allocate(recvbuffer(recvsize))
		recvbuffer = sendbuffer
	elseif (idest .eq. MPI_PROC_NULL .and. sendnp .ne. 0) then
		!Loop through escaping molecules and print
		old => pass%head 
		current => old     !make current point to head of list
		do i=1,sendnp
			call print_mol_escape_error(old%molno)
			old => current%next	!make old point to next node of current
			current => old		!Set current item to old ready for next loop
		enddo
		write(string,'(a,i6,a,i6,a)') "sendrecvface Error: Processor rank ", irank, &
		" is attempting to send ", sendnp, " molecules to MPI_PROC_NULL."
		call error_abort(string)
	else
		!Send, probe for size and then receive data
		call NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
	endif

	!Unpack header data (recvnp)
	pos = 0
	if (length .eq. 0) then
		recvnp = 0
	else
		call MPI_Unpack(recvbuffer,length,pos,dppack,1,MPI_DOUBLE_PRECISION, &
						icomm_grid,ierr)
		recvnp = nint(dppack)
		!Unpack rest of data into correct arrays --------------------!
		do n=new_np+1,new_np+recvnp

			if (ensemble.eq.tag_move) then
				call MPI_Unpack(recvbuffer,length,pos,dppack, &
							1,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
				tag(np+n)     = nint(dppack)
			endif

			call MPI_Unpack(recvbuffer,length,pos,Xpack, &
							nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
			r(:,np+n)     = Xpack

			call MPI_Unpack(recvbuffer,length,pos,Xpack, &
							nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
			v(:,np+n)     = Xpack

			if (rtrue_flag.eq.1) then
				call MPI_Unpack(recvbuffer,length,pos,Xpack, &
								nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
				rtrue(:,np+n) = Xpack
			end if
			
			if (ensemble.eq.tag_move) then
				if (any(tag(np+n).eq.tether_tags)) then
					call MPI_Unpack(recvbuffer,length,pos,Xpack, &
									nd,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
					rtether(:,np+n) = Xpack
				end if
			endif

			if (potential_flag .eq. 1 ) then
				call MPI_Unpack(recvbuffer,length,pos,FENEpack, &
								nsdmi,MPI_DOUBLE_PRECISION,icomm_grid,ierr)
				call assign_FENEbuffer(np+n,FENEpack)
			end if

		enddo
		!-------------------------------------------------------------!

		!Correct positions in new processor
		do n=new_np+1,new_np+recvnp
			r(ixyz,np+n) = r(ixyz,np+n) - dir * domain(ixyz)
			if (ensemble.eq.tag_move) then
				if (any(tag(np+n).eq.tether_tags)) then
					rtether(ixyz,np+n) = rtether(ixyz,np+n) - dir * domain(ixyz)
				endif
			endif
		enddo

	end if



	!Update number of molecules in halo to include number recieved
	new_np = new_np + recvnp
	!new_tethernp = new_tethernp + recvtethernp

	deallocate(recvbuffer)
	deallocate(sendbuffer)
	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

contains

	subroutine get_sendmols_buffersize
		implicit none

		!Include sendnp as double at head of buffer
		sendsize = 1

		!Add 2*nd for r, v and 1 for tag
		if (ensemble.eq.tag_move) then
			sendsize = sendsize + 2*nd*sendnp + 1*sendnp
		else
			sendsize = sendsize + 2*nd*sendnp
		endif
	
		!Add rtrue if necessary	
		if (rtrue_flag.eq.1) then
			sendsize = sendsize + nd*sendnp
		end if

		!Check if any molecules are tethered
		if (ensemble.eq.tag_move) then
			old => pass%head
			current => old
			do i=1,sendnp 
				molno = old%molno
				if (any(tag(molno).eq.tether_tags)) then
					sendsize = sendsize + nd
				end if
				old => current%next
				current => old	
			end do
		endif

		!Add space for polymer info if required
		if (potential_flag.eq.1) then
			sendsize = sendsize + nsdmi*sendnp
		end if

	end subroutine get_sendmols_buffersize

end

!------------------------------------------------------------------------------------
!Move through molecule list of passed molecule and fill with new arrivals starting 
!from the end of the new molecules added to r this timestep

subroutine reorderdata(new_np)
	use computational_constants_MD
 	use physical_constants_MD
	use polymer_info_MD
	use arrays_MD
	use linked_list
	implicit none

	integer			        :: i
	integer			        :: molno,sendnp
	integer, intent(inout)  :: new_np
	type(passnode), pointer :: old, current

	!Work through passed molecules to fill gaps
	old => pass%head
	sendnp = pass%sendnp
	current => old     !make current point to head of list

	!Loop through all empty locations
	do i=1,sendnp

		molno = old%molno	!Position of empty location 

		if (ensemble.eq.tag_move) then
			tag(molno)       = tag(np+new_np)
		endif
		r(:,molno)       = r(:,np+new_np)
		v(:,molno)       = v(:,np+new_np)
		if (rtrue_flag.eq.1) then
			rtrue(:,molno) = rtrue(:,np+new_np)
		endif
		if (ensemble.eq.tag_move) then
			if (any(tag(molno).eq.tether_tags)) then
				rtether(:,molno) = rtether(:,np+new_np)
			endif
		endif
		
		if (potential_flag .eq. 1) then	
			monomer(molno) = monomer(np+new_np)
		end if

		!Read molecular tag and assign correct properties to reordered molecules
		if (ensemble.eq.tag_move) call read_tag(molno)

		!Update number of new molecules
		new_np = new_np - 1

		old => current%next  !make old point to next node of current
		current  => old      !Set current item to old ready for next loop

	enddo

	!Read molecular tag and assign correct properties to new molecules
	if (ensemble.eq.tag_move) then
		do i=np+1, np+new_np
			call read_tag(i)
		enddo
	endif

	!Adjust total number of molecules to reflect difference between
	!recieved molecules and passed molecules
	np = np + new_np

	nullify (current)   !Set current pointer to null
	nullify (old)       !Set old pointer to null

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
		buf1(1) = r(1,molno) - ipass*domain(1)
		buf1(2) = r(2,molno) - jpass*domain(2)
		buf1(3) = r(3,molno) - kpass*domain(3)
		buf1(4) = v(1,molno)
		buf1(5) = v(2,molno)
		buf1(6) = v(3,molno)

		!Obtain processors neighbours
		call procneighbours(ipass,jpass,kpass,isource,idest)

		!print*, irank, 'sending', r(2,molno), 'to', idest+1

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
				r(1,np+n) = buf2(1)
				r(2,np+n) = buf2(2)
				r(3,np+n) = buf2(3)
				v(1,np+n) = buf2(4)
				v(2,np+n) = buf2(5)
				v(3,np+n) = buf2(6)

				!print*, irank, 'receiving', r(2,np+n), 'from', isource+1

				!Correct position on receiving processor
				!r(1,np+n) = r(1,np+n) + ipass*domain(1)
				!r(2,np+n) = r(2,np+n) + jpass*domain(2)
				!r(3,np+n) = r(3,np+n) + kpass*domain(3)

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
end subroutine sendproberecv

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
end subroutine pairedsendproberecv

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

end subroutine NBsendproberecv

!======================================================================
!			Bin averaged handling subroutines        	              =
!======================================================================
! Contains routines:
! swaphalos				--   TOP LEVEL Swap halos with adjacent processors
!							 calls updatefaces, pack and unpack routines
! pack_bins_into_cells  --   For multiple bin per computational cell pack
! 					         into cell ready to exchange with adjacent procs
! unpack_cells_into_bins --  Unpack bin data from cells received from 
!							 adjacent processors
! updatefaces			 --  Facilitate the MPI based exchange of data

!Swap halos of bins between processors
module messenger_bin_handler
	implicit none

	!Generic interface so pack/unpack can be used with both integers and reals
	interface swaphalos
		module procedure iswaphalos, rswaphalos
	end interface swaphalos
	private iswaphalos, rswaphalos

	interface updatefaces
		module procedure iupdatefaces, rupdatefaces
	end interface updatefaces
	private  iupdatefaces, rupdatefaces

	interface pack_bins_into_cells
		module procedure ipack_bins_into_cells, rpack_bins_into_cells
	end interface pack_bins_into_cells
	private  ipack_bins_into_cells, rpack_bins_into_cells

	interface unpack_cells_into_bins
		module procedure iunpack_cells_into_bins, runpack_cells_into_bins
	end interface unpack_cells_into_bins
	private  iunpack_cells_into_bins, runpack_cells_into_bins

contains


!Update face halo cells by passing integers to neighbours 
subroutine iswaphalos(A,n1,n2,n3,nresults)
	use computational_constants_MD, only : ncells,nhb,nhalocells,halocells
	use librarymod, only : heaviside
	implicit none

	integer,intent(in)			:: n1,n2,n3,nresults
	integer,intent(inout)		:: A(:,:,:,:)

	integer									:: n,i,j,k,ic,jc,kc,nresultscell
	integer,dimension(:,:,:,:),allocatable	:: buf

	!Pack bins into array of cells
	nresultscell = nresults * nhb(1) * nhb(2) * nhb(3) 
	allocate(buf(ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell)); buf = 0
	call pack_bins_into_cells(buf,A,nresults)

	!Exchange faces with adjacent processors
	call iupdatefaces(buf,ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell,1)
	call iupdatefaces(buf,ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell,2)
	call iupdatefaces(buf,ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell,3)

	!halo values to correct cells in array
	do n = 1, nhalocells
		i = halocells(n,1); j = halocells(n,2); k = halocells(n,3)

		!Change in number of Molecules in halo cells
		ic = i + heaviside(ncells(1)+1-i)-heaviside(i-2)
		jc = j + heaviside(ncells(2)+1-j)-heaviside(j-2)
		kc = k + heaviside(ncells(3)+1-k)-heaviside(k-2)

		buf(ic,jc,kc,:) = buf(ic,jc,kc,:) + buf(i,j,k,:)

	enddo

	!Unpack array of cells into bins
	call unpack_cells_into_bins(A,buf,nresults)
	deallocate(buf)

end subroutine iswaphalos

!Update face halo cells by passing to neighbours (double precision version)
subroutine rswaphalos(A,n1,n2,n3,nresults)
	use computational_constants_MD, only : ncells,nhb,nhalocells,halocells
	use librarymod, only : heaviside
	implicit none

	integer,intent(in)								:: n1,n2,n3,nresults
	double precision,intent(inout)					:: A(:,:,:,:)

	integer											:: n,i,j,k,ic,jc,kc,nresultscell
	double precision,dimension(:,:,:,:),allocatable	:: buf

	!Pack bins into array of cells
	nresultscell = nresults * nhb(1) * nhb(2) * nhb(3) 
	allocate(buf(ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell)); buf = 0.d0
	call pack_bins_into_cells(buf,A,nresults)

	call rupdatefaces(buf,ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell,1)
	call rupdatefaces(buf,ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell,2)
	call rupdatefaces(buf,ncells(1)+2,ncells(2)+2,ncells(3)+2,nresultscell,3)

	!halo values to correct cells in array
	do n = 1, nhalocells
		i = halocells(n,1); j = halocells(n,2); k = halocells(n,3)

		!Change in number of Molecules in halo cells
		ic = i + heaviside(ncells(1)+1-i)-heaviside(i-2)
		jc = j + heaviside(ncells(2)+1-j)-heaviside(j-2)
		kc = k + heaviside(ncells(3)+1-k)-heaviside(k-2)

		buf(ic,jc,kc,:) = buf(ic,jc,kc,:) + buf(i,j,k,:)

	enddo

	!Unpack array of cells into bins
	call unpack_cells_into_bins(A,buf,nresults)
	deallocate(buf)

end subroutine rswaphalos

!Pack bin data into array of sizes ncells to pass efficiently
subroutine ipack_bins_into_cells(cells,bins,nresults)
	use computational_constants_MD, only : ncells, nhb
	implicit none

	integer,intent(in)						:: nresults
	integer,dimension(:,:,:,:),intent(in)	:: bins
	integer,dimension(:,:,:,:),intent(out)	:: cells

	integer									:: result,icell,jcell,kcell,ibin,jbin,kbin,n

	do icell = 1,ncells(1)+2
	do jcell = 1,ncells(2)+2
	do kcell = 1,ncells(3)+2
		do ibin = 1,nhb(1)
		do jbin = 1,nhb(2)
		do kbin = 1,nhb(3)
			do result = 1,nresults
				n = result + (ibin-1)*nresults + (jbin-1)*nhb(1)*nresults + (kbin-1)*nhb(1)*nhb(2)*nresults
				cells(icell,jcell,kcell,n) = bins((icell-1)*nhb(1)+ibin,(jcell-1)*nhb(2)+jbin,(kcell-1)*nhb(3)+kbin,result)
			enddo

		enddo
		enddo
		enddo
	enddo
	enddo
	enddo

end subroutine ipack_bins_into_cells

!Unpack data from array of sizes ncells into bins 
subroutine iunpack_cells_into_bins(bins,cells,nresults)
	use computational_constants_MD, only : ncells, nhb
	implicit none

	integer,intent(in)						:: nresults
	integer,dimension(:,:,:,:),intent(out)	:: bins
	integer,dimension(:,:,:,:),intent(in)	:: cells

	integer									:: result,icell,jcell,kcell,ibin,jbin,kbin,n

	do icell = 1,ncells(1)+2
	do jcell = 1,ncells(2)+2
	do kcell = 1,ncells(3)+2
		do ibin = 1,nhb(1)
		do jbin = 1,nhb(2)
		do kbin = 1,nhb(3)
			do result = 1,nresults
				n = result + (ibin-1)*nresults + (jbin-1)*nhb(1)*nresults + (kbin-1)*nhb(1)*nhb(2)*nresults
				bins((icell-1)*nhb(1)+ibin,(jcell-1)*nhb(2)+jbin,(kcell-1)*nhb(3)+kbin,result) = cells(icell,jcell,kcell,n) 
			enddo
		enddo
		enddo
		enddo
	enddo
	enddo
	enddo

end subroutine iunpack_cells_into_bins


!Pack bin data into array of sizes ncells to pass efficiently
subroutine rpack_bins_into_cells(cells,bins,nresults)
	use computational_constants_MD, only : ncells, nhb
	implicit none

	integer,intent(in)								:: nresults
	double precision,dimension(:,:,:,:),intent(in)	:: bins
	double precision,dimension(:,:,:,:),intent(out)	:: cells

	integer									:: result,icell,jcell,kcell,ibin,jbin,kbin,n

	do icell = 1,ncells(1)+2
	do jcell = 1,ncells(2)+2
	do kcell = 1,ncells(3)+2
		do ibin = 1,nhb(1)
		do jbin = 1,nhb(2)
		do kbin = 1,nhb(3)
			do result = 1,nresults
				n = result + (ibin-1)*nresults + (jbin-1)*nhb(1)*nresults + (kbin-1)*nhb(1)*nhb(2)*nresults
				cells(icell,jcell,kcell,n) = bins((icell-1)*nhb(1)+ibin,(jcell-1)*nhb(2)+jbin,(kcell-1)*nhb(3)+kbin,result)
			enddo

		enddo
		enddo
		enddo
	enddo
	enddo
	enddo

end subroutine rpack_bins_into_cells

!Unpack data from array of sizes ncells into bins 
subroutine runpack_cells_into_bins(bins,cells,nresults)
	use computational_constants_MD, only : ncells, nhb
	implicit none

	integer,intent(in)								:: nresults
	double precision,dimension(:,:,:,:),intent(out)	:: bins
	double precision,dimension(:,:,:,:),intent(in)	:: cells

	integer											:: result,icell,jcell,kcell,ibin,jbin,kbin,n

	do icell = 1,ncells(1)+2
	do jcell = 1,ncells(2)+2
	do kcell = 1,ncells(3)+2
		do ibin = 1,nhb(1)
		do jbin = 1,nhb(2)
		do kbin = 1,nhb(3)
			do result = 1,nresults
				n = result + (ibin-1)*nresults + (jbin-1)*nhb(1)*nresults + (kbin-1)*nhb(1)*nhb(2)*nresults
				bins((icell-1)*nhb(1)+ibin,(jcell-1)*nhb(2)+jbin,(kcell-1)*nhb(3)+kbin,result) = cells(icell,jcell,kcell,n) 
			enddo
		enddo
		enddo
		enddo
	enddo
	enddo
	enddo

end subroutine runpack_cells_into_bins

!Update face halo cells by passing to neighbours (integer version)
subroutine iupdatefaces(A,n1,n2,n3,nresults,ixyz)
	use messenger, only : icomm_grid, ierr
	use mpi
	implicit none

	integer,intent(in)						:: n1,n2,n3,nresults
	integer,intent(inout)					:: A(:,:,:,:)

	integer 								:: ixyz
	integer 								:: icount,isource,idest
	integer,dimension(:,:,:,:),allocatable	:: buf1, buf2

	!Determine size of send buffer and copy to buffer data to pass to lower neighbour
	select case (ixyz)
	case (1)
		allocate(buf1(1,n2,n3,nresults), buf2(1,n2,n3,nresults))
		icount = 1*n2*n3*nresults
		buf1(1,:,:,:) = A(1,:,:,:)
		buf2 = 0
	case (2)
		allocate(buf1(n1,1,n3,nresults), buf2(n1,1,n3,nresults))
		icount = n1*1*n3*nresults
		buf2 = 0
		buf1(:,1,:,:) = A(:,1,:,:)
	case (3)
		allocate(buf1(n1,n2,1,nresults), buf2(n1,n2,1,nresults))
		icount = n1*n2*1*nresults
		buf2 = 0
		buf1(:,:,1,:) = A(:,:,1,:)
	case default
		stop "updateBorder: invalid value for ixyz"
	end select

	! Send to lower neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_integer, idest, 0, &
	                  buf2, icount, MPI_integer, isource, 0, &
	                  icomm_grid, MPI_STATUS_IGNORE, ierr)

	!Save recieved data from upper neighbour and copy to buffer data to pass to upper neighbour
	select case (ixyz)
	case (1)
		buf1(1,:,:,:) = A(n1,:,:,:)
		A(n1,:,:,:) = buf2(1,:,:,:)
	case (2)
		buf1(:,1,:,:) = A(:,n2,:,:)
		A(:,n2,:,:) = buf2(:,1,:,:)
	case (3)
		buf1(:,:,1,:) = A(:,:,n3,:)
		A(:,:,n3,:) = buf2(:,:,1,:)
	end select

	! Send to upper neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_integer, idest, 0, &
	                  buf2, icount, MPI_integer, isource, 0, &
	                  icomm_grid, MPI_STATUS_IGNORE, ierr)


	!Save recieved data from lower neighbour
	select case (ixyz)
	case (1)
		A(1,:,:,:) = buf2(1,:,:,:)
	case (2)
		A(:,1,:,:) = buf2(:,1,:,:) 
	case (3)
		A(:,:,1,:) = buf2(:,:,1,:)
	end select

	deallocate(buf1, buf2)

end subroutine iupdatefaces


!Update face halo cells by passing to neighbours
subroutine rupdatefaces(A,n1,n2,n3,nresults,ixyz)
	use messenger, only : icomm_grid, ierr
	use mpi
	implicit none

	integer,intent(in)								:: n1,n2,n3,nresults
	double precision,intent(inout)					:: A(:,:,:,:)

	integer 										:: ixyz
	integer 										:: icount,isource,idest
	double precision,dimension(:,:,:,:),allocatable	:: buf1, buf2

	!Determine size of send buffer and copy to buffer data to pass to lower neighbour
	select case (ixyz)
	case (1)
		allocate(buf1(1,n2,n3,nresults), buf2(1,n2,n3,nresults))
		icount = 1*n2*n3*nresults
		buf2 = 0.d0
		buf1(1,:,:,:) = A(1,:,:,:)
	case (2)
		allocate(buf1(n1,1,n3,nresults), buf2(n1,1,n3,nresults))
		icount = n1*1*n3*nresults
		buf2 = 0.d0
		buf1(:,1,:,:) = A(:,1,:,:)
	case (3)
		allocate(buf1(n1,n2,1,nresults), buf2(n1,n2,1,nresults))
		icount = n1*n2*1*nresults
		buf2 = 0.d0
		buf1(:,:,1,:) = A(:,:,1,:)
	case default
		stop "updateBorder: invalid value for ixyz"
	end select

	! Send to lower neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_double_precision, idest, 0, &
	                  buf2, icount, MPI_double_precision, isource, 0, &
	                  icomm_grid, MPI_STATUS_IGNORE, ierr)

	!Save recieved data from upper neighbour and copy to buffer data to pass to upper neighbour
	select case (ixyz)
	case (1)
		buf1(1,:,:,:)= A(n1,:,:,:)
		A(n1,:,:,:)  = buf2(1,:,:,:)
	case (2)
		buf1(:,1,:,:)= A(:,n2,:,:)
		A(:,n2,:,:)  = buf2(:,1,:,:)

	case (3)
		buf1(:,:,1,:)= A(:,:,n3,:)
		A(:,:,n3,:)  = buf2(:,:,1,:)
	end select

	! Send to upper neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_double_precision, idest, 0, &
	                  buf2, icount, MPI_double_precision, isource, 0, &
	                  icomm_grid, MPI_STATUS_IGNORE, ierr)


	!Save recieved data from lower neighbour
	select case (ixyz)
	case (1)
		A(1,:,:,:) = buf2(1,:,:,:)
	case (2)
		A(:,1,:,:) = buf2(:,1,:,:) 
	case (3)
		A(:,:,1,:) = buf2(:,:,1,:)
	end select

	deallocate(buf1, buf2)

end subroutine rupdatefaces

end module messenger_bin_handler


!======================================================================
!					Data gathering subroutines	                      =
!======================================================================

subroutine globalbroadcast(A,na,broadprocid)
	use messenger
	implicit none

	integer				:: na, broadprocid
	double precision	:: A

	call MPI_BCAST(A,na,MPI_DOUBLE_PRECISION,broadprocid-1,MD_COMM,ierr)

	return
end

subroutine globalsyncreduce(A, na, meanA, maxA, minA)
	use messenger
	implicit none
	!include "mpif.h"

 	integer						  :: na, nprocs
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
	implicit none

	double precision :: A, buf

	call MPI_AllReduce (A, buf, 1, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalMax(A)
	use messenger
	implicit none

	double precision :: A, buf

	call MPI_AllReduce (A, buf, 1, MPI_DOUBLE_PRECISION, &
	                    MPI_MAX, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalSumInt(A)
	use messenger
	implicit none

	integer :: A, buf

	call MPI_AllReduce (A, buf, 1, MPI_INTEGER, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalMaxInt(A)
	use messenger
	implicit none

	integer :: A, buf

	call MPI_AllReduce (A, buf, 1, MPI_INTEGER, &
	                    MPI_MAX, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalSumVectReal(A, na)
	use messenger
	implicit none

	integer, intent(in) :: na
	real A(na)
	real buf(na)

	call MPI_AllReduce (A, buf, na, MPI_REAL, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf
	return
end

subroutine globalSumVect(A, na)
	use messenger
	!include "mpif.h"
	implicit none

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
	implicit none

    integer, intent(in) :: na
	integer A(na)
	integer buf(na)

	call MPI_AllReduce (A, buf, na, MPI_INTEGER, &
	                    MPI_SUM, MD_COMM, ierr)
	A = buf

	return
end

subroutine PlaneSumIntVect(PLANE_COMM_IN, A, na)
	use messenger
	implicit none

    integer, intent(in) :: na
	integer, intent(in) :: PLANE_COMM_IN
	integer, intent(inout) :: A(na)
	integer, allocatable :: buf(:)

	allocate(buf(na))

	call MPI_AllReduce (A, buf, na, MPI_INTEGER, &
	                    MPI_SUM, PLANE_COMM_IN, ierr)
	A = buf

	deallocate(buf)

	return

end

subroutine PlaneSumVect(PLANE_COMM_IN, A, na)
	use messenger
	implicit none

    integer, intent(in) :: na
	integer, intent(in) :: PLANE_COMM_IN
	real(kind(0.d0)), intent(inout) :: A(na)
	real(kind(0.d0)), allocatable :: buf(:)

	allocate(buf(na))

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, PLANE_COMM_IN, ierr)
	A = buf

	deallocate(buf)

	return
end

subroutine globalMaxVect(A, na)
	use messenger
	implicit none

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
	implicit none

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
	implicit none

	integer, intent(in) :: na
	double precision A(na)
	double precision buf(na)

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_MIN, MD_COMM, ierr)
	A = buf

	return
end

subroutine globalSumTwoDim(A,na1,na2)
	use messenger
	implicit none

	integer, intent(in) :: na1,na2
	double precision 	:: A(na1,na2)
	double precision 	:: buf(na1,na2)
	
	call MPI_AllReduce (A, buf, na1*na2, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, MD_COMM, ierr)
	
	A = buf

	return
end

subroutine globalSumIntTwoDim(A,na1,na2)
	use messenger
	implicit none

	integer, intent(in) :: na1,na2
	integer A(na1,na2)
	integer buf(na1,na2)
	
	call MPI_AllReduce (A, buf, na1*na2, MPI_INTEGER, &
	                    MPI_SUM, MD_COMM, ierr)
	
	A = buf

	return
end

subroutine globalAverage(A, na)
	use messenger
	implicit none

	integer, intent(in) :: na

	integer	:: nprocs
	double precision A(na)
	double precision buf(na)

	call MPI_AllReduce (A, buf, na, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, MD_COMM, ierr)
	call MPI_comm_size (MD_COMM, nprocs, ierr)

	buf = buf / nprocs
	
	A = buf

	return
end

subroutine globalGather(A,B,na)
	use messenger
	implicit none

	integer, intent(in) :: na

	integer, intent(in) :: A(na)
	integer, intent(out):: B(na,nproc)

	call MPI_Allgather (A, na, MPI_INTEGER, B, na, &
			    MPI_INTEGER,icomm_grid, ierr)

end subroutine

subroutine globalGathernp()
	use physical_constants_MD
	use messenger
	implicit none
	!include "mpif.h"

	call MPI_Allgather (np, 1, MPI_INTEGER, procnp, 1, &
			    MPI_INTEGER,icomm_grid, ierr)

	return
end


!----Sum routines over global sub communitcators
subroutine SubcommGather(A,B,na,ixyz,npixyz)
	use messenger
	implicit none

	integer, intent(in) :: na, ixyz, npixyz
	integer, intent(in) :: A(na)
	integer, intent(out):: B(na,npixyz)

	call MPI_Allgather (A, na, MPI_INTEGER, B, na, &
			    MPI_INTEGER,icomm_xyz(ixyz), ierr)
	
end subroutine

subroutine SubcommSum(A, ixyz)
	use messenger
	implicit none

	integer, intent(in) :: ixyz !Direction of sub-comm
	double precision	A
	double precision buf

	call MPI_AllReduce (A, buf, 1, MPI_DOUBLE_PRECISION, &
	                    MPI_SUM, icomm_xyz(ixyz), ierr)
	A = buf

	return
end

subroutine SubcommSumInt(A, ixyz)
	use messenger
	implicit none

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
	implicit none

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
	implicit none

    integer, intent(in) :: na, ixyz !Direction of sub-comm

	integer	A(na)
	integer buf(na)

	call MPI_AllReduce (A, buf, na, MPI_INTEGER, &
         MPI_SUM, icomm_xyz(ixyz), ierr)
	A = buf

	return
end subroutine SubcommSumIntVect

subroutine error_abort_s(msg)
    use mpi
    implicit none

    character(len=*), intent(in), optional :: msg
   
	integer errcode,ierr

    if (present(msg)) then 
        write(*,*) msg
    endif

    call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)

end subroutine error_abort_s


subroutine error_abort_si(msg,i)
    use mpi
    implicit none

    character(len=*), intent(in) :: msg
    integer, intent(in) :: i

    integer errcode,ierr

    write(*,*) msg,i

    call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)

end subroutine error_abort_si

