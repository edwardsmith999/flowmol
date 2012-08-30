!====================================================================
!     ******************************************************************
!     *                                                                *
!     *   Output Z-planes of the instantaneous field                   *
!     *   Each plane is (ngx,ngy):                                     * 
!     *                   No halos in order to ease post-processing.   *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!=======================================================================


module Output_OneZplane_mod
      use messenger
      use data_export
      implicit none
	include "mpif.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)	   :: global_indices, local_indices
	INTEGER			   :: filetype, memtype
	CHARACTER(LEN=18)          :: FN2
	CHARACTER(LEN=12)          :: NAME

	INTEGER, PARAMETER            :: n_Kplanes = 4
	INTEGER, DIMENSION(n_Kplanes), PARAMETER :: Kplane = (/ 1 , ngzm/4 , 2*ngzm/4 , 3*ngzm/4 /)
	!INTEGER, DIMENSION(n_Kplanes), PARAMETER :: Kplane = (/ 38, 56, 57, 58 /)

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0

	REAL aan
end module



subroutine WriteZplane()
use Output_OneZplane_mod

	REAL	, DIMENSION( 1 , 1:nlx, 1:nly) :: u_Zplane
	REAL	, DIMENSION( 1 , 1:nlx, 1:nly) :: v_Zplane
	REAL	, DIMENSION( 1 , 1:nlx, 1:nly) :: w_Zplane

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!---------------------------------------------------
	integer	:: is, js, ks
	integer	:: ie, je, ke
	double precision, allocatable	:: OutBuffer(:,:,:)

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

     !=======================================================
     !         CARTESIAN VELOCITIES
     !=======================================================
        NAME='Zplane.dble.'
        call OneZplane_FileName()
        if (irank.eq.iroot) write(*,*) FN2

        call MPI_FILE_OPEN(icomm_grid, FN2, &
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/  0  , iTmin_1(iblock)-1 , jTmin_1(jblock)-1 /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/  0  , nox1-1 , noy1-1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices = (/  0  , 0  , 0 /)				!ALC_FIX
			ks = 1;		is = nox1;	js = noy1;			!ALC_FIX
			if (iblock.eq. 1 ) is = 1					!ALC_FIX
			if (jblock.eq. 1 ) js = 1					!ALC_FIX


	!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
	gsizes = (/  1  , ngx , ngy /)
	lsizes = (/  1  , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=(nlxb-nox1)-(1)   +1
	if (iblock.eq.npx)   lsizes(2)=(nlxb+0)   -(nox1)+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx
	if (jblock.eq. 1 )   lsizes(3)=(nlyb-noy1)-(1)   +1
	if (jblock.eq.npy)   lsizes(3)=(nlyb+0)   -(noy1)+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy

	!ALC_BUG	!--------- MEMORY SIZE OF (u_Zplane, v_Zplane, w_Zplane) -------
	!ALC_BUG	memsizes = (/  1  , nlx , nly /)

        global_cnt = 0
        offset     = 0

                allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
                ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
                memsizes = lsizes							!ALC_FIX

	!------------------------ uc ---------------------------
	do k = 1,n_Kplanes
		kndx = Kplane(k)
		do j=1,nly
			u_Zplane(1,:,j) = 0.25*( uc( kndx-1, 1:nlx, j-1)+uc( kndx  , 1:nlx, j-1) &
						+uc( kndx  , 1:nlx, j  )+uc( kndx-1, 1:nlx, j  )  )
		end do

		OutBuffer = u_Zplane(ks:ke, is:ie, js:je)				!ALC_FIX

		CALL Create_commit_OneZplane_fileview()
		CALL Create_commit_OneZplane_subarray() 
		CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)

        	global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        	offset = global_cnt * FloatSize
	end do
	!------------------------ vc ---------------------------
	do k = 1,n_Kplanes
		kndx = Kplane(k)
		do i=1,nlx
			v_Zplane(1,i,:) = 0.25*( vc( kndx-1, i-1, 1:nly)+vc( kndx  , i-1, 1:nly) &
						+vc( kndx  , i  , 1:nly)+vc( kndx-1, i  , 1:nly)  )
		end do

		OutBuffer = v_Zplane(ks:ke, is:ie, js:je)				!ALC_FIX

		CALL Create_commit_OneZplane_fileview()
		CALL Create_commit_OneZplane_subarray() 
		CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)

        	global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        	offset = global_cnt * FloatSize
	end do

       !--------------- wc ------------------
	do k = 1,n_Kplanes
		kndx = Kplane(k)
		do j=1,nly
		do i=1,nlx
			w_Zplane(1,i,j)= 0.25*(  wc(kndx, i-1, j-1)+wc(kndx, i-1, j  ) &
						+wc(kndx, i  , j  )+wc(kndx, i  , j-1)  )
		end do
		end do

		OutBuffer = w_Zplane(ks:ke, is:ie, js:je)				!ALC_FIX

		CALL Create_commit_OneZplane_fileview()
		CALL Create_commit_OneZplane_subarray() 
        	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)

        	global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        	offset = global_cnt * FloatSize
	end do

        CALL MPI_FILE_CLOSE(fh, ierr)
                deallocate(OutBuffer)							!ALC_FIX

       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
        CALL MPI_BARRIER(icomm_grid,IERR)
        if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from processor=',irank

	call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
	call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0
return
end


subroutine Create_commit_OneZplane_fileview()
	use Output_OneZplane_mod

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

subroutine Create_commit_OneZplane_subarray()
	use Output_OneZplane_mod

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

subroutine OneZplane_FileName()
        use Output_OneZplane_mod

        IF(ntime.LE.9                        ) &
        WRITE(FN2,'(A12,A5,I1)') NAME,'00000',ntime
        IF(ntime.GE.10     .AND. ntime.LE.99    ) &
        WRITE(FN2,'(A12,A4,I2)') NAME,'0000' ,ntime
        IF(ntime.GE.100    .AND. ntime.LE.999   ) &
        WRITE(FN2,'(A12,A3,I3)') NAME,'000'  ,ntime
        IF(ntime.GE.1000   .AND. ntime.LE.9999  ) &
        WRITE(FN2,'(A12,A2,I4)') NAME,'00'   ,ntime
        IF(ntime.GE.10000  .AND. ntime.LE.99999 ) &
        WRITE(FN2,'(A12,A1,I5)') NAME,'0'    ,ntime
        IF(ntime.GE.100000 .AND. ntime.LE.999999) &
        WRITE(FN2,'(A12,   I6)') NAME,        ntime

        return
end
