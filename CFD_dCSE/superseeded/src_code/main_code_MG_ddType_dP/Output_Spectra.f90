!====================================================================
!     ******************************************************************
!     *   Read in COMPLEX uhat(0:ngx+1,0:ngy)			       *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!=======================================================================


module Output_uhat_mod
	use messenger
	use statistics
	implicit none
	include "mpif.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(4)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(4)      :: global_indices, local_indices
	INTEGER, DIMENSION(3)      :: gsizes_nd, lsizes_nd, memsizes_nd
	INTEGER, DIMENSION(3)      :: global_indices_nd, local_indices_nd
	INTEGER, DIMENSION(3)      :: gsizes_ez, lsizes_ez, memsizes_ez
	INTEGER, DIMENSION(3)      :: global_indices_ez, local_indices_ez
	INTEGER                    :: filetype, memtype

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0
	REAL aan
end module


subroutine Write_uhat()
use Output_uhat_mod
use data, only : file_dir
	! real	:: uhatR(0:nlx+1, 0:nly+1)
	! real	:: uhatC(0:nlx+1, 0:nly+1)
	real	:: uhatR(ngzm, nfreq_, inx, iny)
	real	:: uhatC(ngzm, nfreq_, inx, iny)
	real, dimension(2,inx,iny) :: ndx_loc
        character(len=100) local_fname

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

     !==============================================================
     !			Output of Uhat(:,:,:,:)	
     !==============================================================
        if (irank.eq.iroot) write(*,*) 'writting  ', FNAME
        local_fname=trim(file_dir)//FNAME
        call MPI_FILE_OPEN(icomm_grid, local_fname, &
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

        !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
        global_indices = (/ 0 , 0 , inx*(iblock-1) , iny*(jblock-1) /)
        local_indices  = (/ 0 , 0 , 0 , 0 /)
        memsizes       = (/ ngzm , nfreq_ , inx, iny /)

     !==============================================================
     !				U_HAT(:,:,:,:)
     !==============================================================
     !			real( uhat ), aimag( uhat )
     !==============================================================
        global_cnt = 0
        offset     = 0
        
	uhatR = real (uhat)
	uhatC = aimag(uhat)

        gsizes = (/ ngzm , nfreq_ , inx*npx, iny*npy /)
        lsizes = (/ ngzm , nfreq_ ,   inx  ,   iny   /)

        CALL Array_File_Type_uhat() 

        !------------------------ Real Part  ---------------------------
        CALL MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, filetype, &
                               'native', MPI_INFO_NULL, ierr)
	CALL MPI_FILE_WRITE_ALL(fh, uhatR, 1, memtype, status, ierr)

	global_cnt = global_cnt + product(gsizes)
	offset = global_cnt * FloatSize

        !------------------------ Imaginary Part ---------------------------
        CALL MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, filetype, &
                               'native', MPI_INFO_NULL, ierr)
	CALL MPI_FILE_WRITE_ALL(fh, uhatC, 1, memtype, status, ierr)

	global_cnt = global_cnt + product(gsizes)
	offset = global_cnt * FloatSize

        !------------------------ -------------- ---------------------------
	!		Write the indices of (inx,iny) probes
        !------------------------ -------------- ---------------------------
	! write(*,'(a,20(i5))')'irank,iblock,jblock,ibmap,jbmap = ',irank,iblock,jblock,&
	! 	(ibmap_1(indx(i)),i=1,inx), &
	! 	(jbmap_1(indy(j)),j=1,iny)

	!---- Store them in (ndx_loc) ------
	do j = 1,iny
		ndx_loc(2,:,j) = real( jbmap_1(indy(j)) )
		do i = 1,inx
			ndx_loc(1,i,j) = real( ibmap_1(indx(i)) )
		end do
	end do

        global_indices_nd = (/ 0 , inx*(iblock-1) , iny*(jblock-1) /)
        local_indices_nd  = (/ 0 , 0 , 0 /)
        memsizes_nd     = (/ 2 , inx, iny /)

        gsizes_nd = (/ 2 , inx*npx, iny*npy /)
        lsizes_nd = (/ 2 ,   inx  ,   iny   /)

        CALL Array_File_Type_nd() 

        !------------------------ Append to EOF  ---------------------------
        CALL MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, filetype, &
                               'native', MPI_INFO_NULL, ierr)
	CALL MPI_FILE_WRITE_ALL(fh, ndx_loc, 1, memtype, status, ierr)

	global_cnt = global_cnt + product(gsizes)
	offset = global_cnt * FloatSize

	CALL MPI_FILE_CLOSE(fh,ierr)
	call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
	call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0


	return
end


subroutine Array_File_Type_uhat()
use Output_uhat_mod
implicit none

	if (FILE_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(filetype,ierr)
		FILE_FLAG = 0
	end if
        !----- Create File-Type --------
        CALL MPI_TYPE_CREATE_SUBARRAY(4, gsizes, lsizes, global_indices, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  filetype, ierr)
        CALL MPI_TYPE_COMMIT(filetype, ierr)
	FILE_FLAG = 1

	if (MEM_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(memtype,ierr)
		MEM_FLAG = 0
	end if
        !----- Create Memory-Type -----
        CALL MPI_TYPE_CREATE_SUBARRAY(4, memsizes, lsizes, local_indices, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  memtype, ierr)
        CALL MPI_TYPE_COMMIT(memtype, ierr)
	MEM_FLAG = 1

	return
end


subroutine Array_File_Type_nd()
use Output_uhat_mod
implicit none

	if (FILE_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(filetype,ierr)
		FILE_FLAG = 0
	end if
        !----- Create File-Type --------
        CALL MPI_TYPE_CREATE_SUBARRAY(3, gsizes_nd, lsizes_nd, global_indices_nd, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  filetype, ierr)
        CALL MPI_TYPE_COMMIT(filetype, ierr)
	FILE_FLAG = 1

	if (MEM_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(memtype,ierr)
		MEM_FLAG = 0
	end if
        !----- Create Memory-Type -----
        CALL MPI_TYPE_CREATE_SUBARRAY(3, memsizes_nd, lsizes_nd, local_indices_nd, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  memtype, ierr)
        CALL MPI_TYPE_COMMIT(memtype, ierr)
	MEM_FLAG = 1

	return
end


!======================== FINISHED (uhat) and STARTING (Ez) ======================


subroutine Write_Ez()
use Output_uhat_mod
use data, only : file_dir

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!---------------------------------------------------
	double precision, allocatable :: OutBuffer(:,:,:)
        character(len=100) :: local_fname

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

     !==============================================================
     !			Output of Ez(:,:,:)
     !==============================================================
        if (irank.eq.iroot) write(*,*) 'writting  ', FNAME
        local_fname=trim(file_dir)//FNAME
        call MPI_FILE_OPEN(icomm_grid, local_fname, &
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

        !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
        global_indices_ez = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
        if (iblock.eq. 1 ) global_indices_ez(2) = 0
        if (jblock.eq. 1 ) global_indices_ez(3) = 0

        !ALC_BUG	local_indices_ez = (/ 0 , nox1 , noy1 /)
        !ALC_BUG	if (iblock.eq. 1 ) local_indices_ez(2) = 0
        !ALC_BUG	if (jblock.eq. 1 ) local_indices_ez(3) = 0
			local_indices_ez = (/ 0 , 0, 0 /)				!ALC_FIX
			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
			if (iblock.eq. 1 ) is = 0					!ALC_FIX
			if (jblock.eq. 1 ) js = 0					!ALC_FIX
        !ALC_BUG	memsizes_ez = (/ ngzm/2+1 , nlx+2 , nly+2 /)

     !==============================================================
     ! (iuvw = 1)	Ez(0:ngzm/2, 0:ngx  , 0:ngy  )
     ! (iuvw = 2)	Ez(0:ngzm/2, 0:ngx  , 0:ngy  )
     ! (iuvw = 3)	Ez(0:ngzm/2, 0:ngx  , 0:ngy  )
     !==============================================================
     !			Ez(0:ngzm/2, 0:nlx+1, 0:nly+1)
     !==============================================================
        global_cnt = 0
        offset     = 0

       	gsizes_ez = (/ ngzm/2+1 , ngx+1 , ngy+1 /)
       	lsizes_ez = (/ ngzm/2+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
       	if (iblock.eq. 1 )   lsizes_ez(2)=(nlxb-nox1)-  0  +1
       	if (iblock.eq.npx)   lsizes_ez(2)=(nlxb+  0 )-nox1 +1
	if ( npx  .eq. 1 )   lsizes_ez(2)=gsizes_ez(2)
       	if (jblock.eq. 1 )   lsizes_ez(3)=(nlyb-noy1)-  0  +1
       	if (jblock.eq.npy)   lsizes_ez(3)=(nlyb+  0 )-noy1 +1
	if ( npy  .eq. 1 )   lsizes_ez(3)=gsizes_ez(3)

		allocate(OutBuffer(lsizes_ez(1),lsizes_ez(2),lsizes_ez(3)))		!ALC_FIX
		memsizes_ez = lsizes_ez							!ALC_FIX
		ke = ks+lsizes_ez(1)-1;	ie = is+lsizes_ez(2)-1;	je = js+lsizes_ez(3)-1	!ALC_FIX
		OutBuffer = Ez(ks:ke, is:ie, js:je)					!ALC_FIX

        CALL Array_File_Type_Ez() 

        !------------------------ Write Ez  ---------------------------
        CALL MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, filetype, &
                               'native', MPI_INFO_NULL, ierr)
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		deallocate(OutBuffer)							!ALC_FIX

	global_cnt = global_cnt + product(gsizes_ez)
	offset = global_cnt * FloatSize


	CALL MPI_FILE_CLOSE(fh,ierr)
	call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
	call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	return
end


subroutine Array_File_Type_Ez()
use Output_uhat_mod
implicit none

	if (FILE_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(filetype,ierr)
		FILE_FLAG = 0
	end if
        !----- Create File-Type --------
        CALL MPI_TYPE_CREATE_SUBARRAY(3, gsizes_ez, lsizes_ez, global_indices_ez, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  filetype, ierr)
        CALL MPI_TYPE_COMMIT(filetype, ierr)
	FILE_FLAG = 1

	if (MEM_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(memtype,ierr)
		MEM_FLAG = 0
	end if
        !----- Create Memory-Type -----
        CALL MPI_TYPE_CREATE_SUBARRAY(3, memsizes_ez, lsizes_ez, local_indices_ez, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  memtype, ierr)
        CALL MPI_TYPE_COMMIT(memtype, ierr)
	MEM_FLAG = 1

	return
end


