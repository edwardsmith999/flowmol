
!====================================================================
!			STATISTICS at Center of C.V.
!====================================================================
!     ******************************************************************
!     *   Read in Statistics					       *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!=======================================================================

module Input_Stats_mod
	use messenger
	use statistics
	implicit none
	include "mpif.h"
	! include "mpiof.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)      :: global_indices, local_indices
	INTEGER                    :: filetype, memtype

	CHARACTER(LEN=10)          :: StatNAME
	CHARACTER(LEN=16)          :: NAME

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0
end module


subroutine Read_Stats()
use Input_Stats_mod
	!----------- RECALL THE SIZES ----------
        ! real Ui (0:nlx+1, 0:nly+1, 6)		! velocity statistics
        ! real Rij(0:nlx+1, 0:nly+1, 6)		! Reynolds stresses
        ! real UiP(0:nlx+1, 0:nly+1, 3)		! Pressure Corrolations
        ! real UiUj2(0:nlx+1, 0:nly+1, 3)	! TKE transport (u*TKE, v*TKE, w*TKE)
        ! real eps(0:nlx+1, 0:nly+1, 1)		! Dissipation (\epsilon)
	!---------------------------------------

	!----------------------------------------------------
	! (InBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!---------------------------------------------------
	double precision, allocatable :: InBuffer(:,:,:)

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

     !=======================================================
	NAME='Stats_TKE.'
	StatNAME = NAME(1:10)
	call In_Stats_FileName()
        if (irank.eq.iroot) write(*,*) 'reading  ', NAME

        call MPI_FILE_OPEN(icomm_grid, NAME, &
                           MPI_MODE_RDONLY, &
                           MPI_INFO_NULL, fh, ierr)

        !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	!  Limit don't change from one statistics array to another
	!  Only the number of quantities calculated (not the mesh)
	!-------------------------------------------------------
        global_indices = (/ iTmin_1(iblock) , jTmin_1(jblock) , 0/)
        if (iblock.eq. 1 ) global_indices(1) = 0
        if (jblock.eq. 1 ) global_indices(2) = 0

        !ALC_BUG	local_indices = (/ nox1 , noy1 , 0 /)
        !ALC_BUG	if (iblock.eq. 1 ) local_indices(1) = 0
        !ALC_BUG	if (jblock.eq. 1 ) local_indices(2) = 0
			local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
			is = nox1;	js = noy1;		ks = 1			!ALC_FIX
			if (iblock.eq. 1 ) is = 0					!ALC_FIX
			if (jblock.eq. 1 ) js = 0					!ALC_FIX
     
        gsizes = (/ ngx+1 , ngy+1 , 1 /)
	!ALC_BUG	memsizes = (/ nlx+2 , nly+2 , 1 /)

        lsizes = (/ nlxb-2*nox1+1 , nlyb-2*noy1+1 , 1/)
        if (iblock.eq. 1 )   lsizes(1)=nlxb-nox1-0+1
        if (iblock.eq.npx)   lsizes(1)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(1)=ngx+1
        if (jblock.eq. 1 )   lsizes(2)=nlyb-noy1-0+1
        if (jblock.eq.npy)   lsizes(2)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(2)=ngy+1

        		memsizes = lsizes						!ALC_FIX

        !------------------------ Ui ---------------------------
        global_cnt = 0
        offset     = 0

	gsizes(3) = 6	;	lsizes(3) = 6	;	memsizes(3) = 6
        
		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
        CALL Create_commit_Stats_inputview()
        CALL Create_commit_Stats_in_subarray()
        CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
		ie = is+lsizes(1)-1;	je = js+lsizes(2)-1;	ke = ks+lsizes(3)-1	!ALC_FIX
		Ui(is:ie, js:je, ks:ke) = InBuffer					!ALC_FIX
		deallocate(InBuffer)							!ALC_FIX
        
        !------------------------ Rij ---------------------------
        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        offset = global_cnt * FloatSize
        
	gsizes(3) = 6	;	lsizes(3) = 6	;	memsizes(3) = 6

		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
        CALL Create_commit_Stats_inputview()
        CALL Create_commit_Stats_in_subarray()
        CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
		ie = is+lsizes(1)-1;	je = js+lsizes(2)-1;	ke = ks+lsizes(3)-1	!ALC_FIX
		Rij(is:ie, js:je, ks:ke) = InBuffer					!ALC_FIX
		deallocate(InBuffer)							!ALC_FIX

        !------------------------ UiP ---------------------------
        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        offset = global_cnt * FloatSize
        
	gsizes(3) = 3	;	lsizes(3) = 3	;	memsizes(3) = 3

		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
        CALL Create_commit_Stats_inputview()
        CALL Create_commit_Stats_in_subarray()
        CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
		ie = is+lsizes(1)-1;	je = js+lsizes(2)-1;	ke = ks+lsizes(3)-1	!ALC_FIX
		UiP(is:ie, js:je, ks:ke) = InBuffer					!ALC_FIX
		deallocate(InBuffer)							!ALC_FIX

        !------------------------ UiUj2 ---------------------------
        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        offset = global_cnt * FloatSize
        
	gsizes(3) = 3	;	lsizes(3) = 3	;	memsizes(3) = 3

		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
        CALL Create_commit_Stats_inputview()
        CALL Create_commit_Stats_in_subarray()
        CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
		ie = is+lsizes(1)-1;	je = js+lsizes(2)-1;	ke = ks+lsizes(3)-1	!ALC_FIX
		UiUj2(is:ie, js:je, ks:ke) = InBuffer					!ALC_FIX
		deallocate(InBuffer)							!ALC_FIX

        !------------------------ eps ---------------------------
        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        offset = global_cnt * FloatSize
        
	gsizes(3) = 1	;	lsizes(3) = 1	;	memsizes(3) = 1

		allocate(InBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
        CALL Create_commit_Stats_inputview()
        CALL Create_commit_Stats_in_subarray()
        CALL MPI_FILE_READ_ALL(fh, InBuffer, 1, memtype, status, ierr)
		ie = is+lsizes(1)-1;	je = js+lsizes(2)-1;	ke = ks+lsizes(3)-1	!ALC_FIX
		eps(is:ie, js:je, ks:ke) = InBuffer					!ALC_FIX
		deallocate(InBuffer)							!ALC_FIX

        CALL MPI_FILE_CLOSE(fh, ierr)
	CALL MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
	CALL MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

return
end



subroutine  Create_commit_Stats_inputview()
        use Input_Stats_mod

	if (FILE_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(filetype,ierr)
		FILE_FLAG = 0
	end if
        CALL MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, global_indices, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  filetype, ierr)
        CALL MPI_TYPE_COMMIT(filetype, ierr)
        CALL MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, filetype, &
                               'native', MPI_INFO_NULL, ierr)
	FILE_FLAG = 1
        return
end

subroutine Create_commit_Stats_in_subarray()
        use Input_Stats_mod

	if (MEM_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(memtype,ierr)
		MEM_FLAG = 0
	end if
        CALL MPI_TYPE_CREATE_SUBARRAY(3, memsizes, lsizes, local_indices, &
                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  memtype, ierr)
        CALL MPI_TYPE_COMMIT(memtype, ierr)
	MEM_FLAG = 1
        return
end

subroutine In_Stats_FileName()
	use Input_Stats_mod

        IF(ntime.LE.9                        ) &
        WRITE(NAME,'(A10,A5,I1)') StatNAME,'00000',ntime
        IF(ntime.GE.10     .AND. ntime.LE.99    ) &
        WRITE(NAME,'(A10,A4,I2)') StatNAME,'0000' ,ntime
        IF(ntime.GE.100    .AND. ntime.LE.999   ) &
        WRITE(NAME,'(A10,A3,I3)') StatNAME,'000'  ,ntime
        IF(ntime.GE.1000   .AND. ntime.LE.9999  ) &
        WRITE(NAME,'(A10,A2,I4)') StatNAME,'00'   ,ntime
        IF(ntime.GE.10000  .AND. ntime.LE.99999 ) &
        WRITE(NAME,'(A10,A1,I5)') StatNAME,'0'    ,ntime
        IF(ntime.GE.100000 .AND. ntime.LE.999999) &
        WRITE(NAME,'(A10,   I6)') StatNAME,        ntime

	return
end
