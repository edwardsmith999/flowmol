!====================================================================
!     ******************************************************************
!     *                                                                *
!     *   SAVE THE RESTART IMAGE TO ONE SINGLE FILE USING              *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!=======================================================================


module Output_mod
      use messenger
      use data_export
      implicit none
	include "mpif.h"
	! include "mpiof.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)	   :: global_indices, local_indices
	INTEGER			   :: filetype, memtype
	CHARACTER(LEN=18)          :: FN2
	CHARACTER(LEN=12)          :: NAME

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0

	REAL aan
end module



subroutine WriteRstrt()
use Output_mod

	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	double precision, allocatable :: OutBuffer(:,:,:)

	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!************************************************
	! Note (ntime, ntime_, stime, stime_) are global variables
	! ntime  :  time step at end of this simulation
	! ntime_ :  time step at begining of this simulation

	!=========================================================
	!                 DATA
	!=========================================================
	IF (irank.EQ.iroot) THEN
		open(28,file='data',form='formatted')
		write(28,*) ntime
		close(28)
	END IF

	aan=1.0*ntime
	!========================================================




	!===================================
	!   Output (ntime, stime, aan)
	!===================================
	IF (irank.EQ.iroot) THEN 
		NAME='ucvcwc.data.'
		call OutFileName()
		write(*,*) FN2
		open(29,file=FN2,form='formatted')
		write(29,'(i10,2e16.8)') ntime,stime,aan
		close(29)
	END IF

     !=======================================================
     !         CARTESIAN VELOCITIES
     !=======================================================
        NAME='ucvcwc.dble.'
        call OutFileName()
        if (irank.eq.iroot) write(*,*) FN2

        call MPI_FILE_OPEN(icomm_grid, FN2, &
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
			if (iblock.eq. 1 ) is = 0					!ALC_FIX
			if (jblock.eq. 1 ) js = 0					!ALC_FIX

	!------------------------ uc ---------------------------
        global_cnt = 0
        offset     = 0

	gsizes = (/ ngz+1 , ngx+2 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+1-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx+2
	if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy+1

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = uc(ks:ke, is:ie, js:je)					!ALC_FIX
		memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+2 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		deallocate(OutBuffer)							!ALC_FIX

	!------------------------ vc ---------------------------
        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        offset = global_cnt * FloatSize

	gsizes = (/ ngz+1 , ngx+1 , ngy+2 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx+1
	if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+1-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy+2

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = vc(ks:ke, is:ie, js:je)					!ALC_FIX
		memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+2 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		deallocate(OutBuffer)							!ALC_FIX

       !--------------- wc ------------------
        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
	offset = global_cnt*FloatSize

	gsizes = (/ ngz+2 , ngx+1 , ngy+1 /)
        lsizes = (/ ngz+2 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
        if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
        if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx+1
        if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
        if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy+1
     
		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = wc(ks:ke, is:ie, js:je)					!ALC_FIX
        	memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+2 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray() 
        CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
        CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)							!ALC_FIX

       !=========================================================
       !                 FLUXES
       !=========================================================
        NAME='uuvvww.dble.'
        call OutFileName()
        if (irank.eq.iroot) write(*,*) FN2

        CALL MPI_FILE_OPEN(icomm_grid, FN2, & 
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
			if (iblock.eq. 1 ) is = 0					!ALC_FIX
			if (jblock.eq. 1 ) js = 0					!ALC_FIX

	!------------------------ u  ---------------------------
        global_cnt = 0
        offset     = 0

	gsizes = (/ ngz+1 , ngx+2 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+1-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx+2
	if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy+1

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = u (ks:ke, is:ie, js:je)					!ALC_FIX
		memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+2 , nly+1 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		deallocate(OutBuffer)							!ALC_FIX


       !---------------- v -----------------
        global_cnt = global_cnt + (ngx+2)*(ngy+1)*(ngz+1) 
        offset = global_cnt * FloatSize

	gsizes = (/ ngz+1 , ngx+1 , ngy+2 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
	if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx+1
	if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
	if (jblock.eq.npy)   lsizes(3)=nlyb+1-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy+2

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = v (ks:ke, is:ie, js:je)					!ALC_FIX
		memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+2 /)
	CALL Create_commit_subarray() 
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		deallocate(OutBuffer)							!ALC_FIX

       !--------------- w  ------------------
	global_cnt = global_cnt + (ngx+1)*(ngy+2)*(ngz+1)
	offset = global_cnt*FloatSize

	gsizes = (/ ngz+2 , ngx+1 , ngy+1 /)
        lsizes = (/ ngz+2 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
        if (iblock.eq. 1 )   lsizes(2)=nlxb-nox1-0+1
        if (iblock.eq.npx)   lsizes(2)=nlxb+0-nox1+1
	if ( npx  .eq. 1 )   lsizes(2)=ngx+1
        if (jblock.eq. 1 )   lsizes(3)=nlyb-noy1-0+1
        if (jblock.eq.npy)   lsizes(3)=nlyb+0-noy1+1
	if ( npy  .eq. 1 )   lsizes(3)=ngy+1

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = w (ks:ke, is:ie, js:je)					!ALC_FIX
        	memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+2 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray() 
        CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
        CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)							!ALC_FIX

       !=========================================================
       !                  CONVECTIVE TERMS 
       !=========================================================
        NAME='conold.dble.'
        call OutFileName()
        if (irank.eq.iroot) write(*,*) FN2

        call MPI_FILE_OPEN(icomm_grid, FN2, & 
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	!  Note:  MPI assumes here that numbering starts from zero
	!  Since my numbering starts from (one), I have to subtract
	!  (one) from every indeex
	!-------------------------------------------------------
	global_indices = (/ 0 , iTmin_1(iblock)-1 , jTmin_1(jblock)-1 /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1-1 , noy1-1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
			ks = 1;		is = nox1;		js = noy1		!ALC_FIX
			if (iblock.eq. 1 ) is = 1					!ALC_FIX
			if (jblock.eq. 1 ) js = 1					!ALC_FIX

       !--------------- conx ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz-1 , ngx-1 , ngy-1 /)
	lsizes = (/ ngz-1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=(nlxb-nox1)-(1)   +1
	if (iblock.eq.npx)   lsizes(2)=(nlxb-1)   -(nox1)+1
	if ( npx  .eq. 1 )   lsizes(2)=(ngx-1)
	if (jblock.eq. 1 )   lsizes(3)=(nlyb-noy1)-(1)   +1
	if (jblock.eq.npy)   lsizes(3)=(nlyb-1)   -(noy1)+1
	if ( npy  .eq. 1 )   lsizes(3)=(ngy-1)

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = conx(ks:ke, is:ie, js:je)					!ALC_FIX
		memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz-1 , nlx-1 , nly-1 /)
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)

       !---------------- cony -----------------
	global_cnt = global_cnt + (ngx-1)*(ngy-1)*(ngz-1)
	offset = global_cnt * FloatSize

		OutBuffer = cony(ks:ke, is:ie, js:je)					!ALC_FIX
	CALL Create_commit_fileview()
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)

       !--------------- conz ------------------
	global_cnt = global_cnt + (ngx-1)*(ngy-1)*(ngz-1)
	offset = global_cnt*FloatSize

	!OR  global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
	!OR  global_cnt = global_cnt + product(gsizes)

		OutBuffer = conz(ks:ke, is:ie, js:je)					!ALC_FIX
	CALL Create_commit_fileview()
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		deallocate(OutBuffer)							!ALC_FIX

	!----------- codetime ------------
	!---- Should all be in Archive ---
	!  global_cnt = global_cnt + (ngx-1)*(ngy-1)*(ngz-1)
	!  offset = global_cnt*FloatSize
	!  IF (irank.EQ.iroot) THEN
	!     call MPI_FILE_SET_VIEW(fh, offset, MPI_DOUBLE_PRECISION, &
	!                            MPI_DOUBLE_PRECISION, 'native',   &
	!                            MPI_INFO_NULL, ierr)
	!     print *, 'stime = ', stime
	!     call MPI_FILE_WRITE(fh, stime, 1, MPI_DOUBLE_PRECISION, &
	!                             status, ierr)
	!     print *, 'aan      = ', aan
	!     call MPI_FILE_WRITE(fh, aan, 1, MPI_DOUBLE_PRECISION, &
	!                             status, ierr)
	!  END IF



	CALL MPI_FILE_CLOSE(fh, ierr)

       !=========================================================
       !                 PRESSURE
       !=========================================================
        NAME='pres_p.dble.'
        call OutFileName()
        if (irank.eq.iroot) write(*,*) FN2

        call MPI_FILE_OPEN(icomm_grid, FN2, & 
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices = (/ 0 , nox1 , noy1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices = (/ 0 , 0, 0 /)					!ALC_FIX
			ks = 0;		is = nox1;		js = noy1		!ALC_FIX
			if (iblock.eq. 1 ) is = 0					!ALC_FIX
			if (jblock.eq. 1 ) js = 0					!ALC_FIX

       !--------------- p ------------------
	global_cnt = 0
	offset     = 0

	gsizes = (/ ngz+1 , ngx+1 , ngy+1 /)
	lsizes = (/ ngz+1 , nlxb-2*nox1+1 , nlyb-2*noy1+1 /)
	if (iblock.eq. 1 )   lsizes(2)=(nlxb-nox1)-(0)   +1
	if (iblock.eq.npx)   lsizes(2)=(nlxb+0)   -(nox1)+1
	if ( npx  .eq. 1 )   lsizes(2)=(ngx+1)
	if (jblock.eq. 1 )   lsizes(3)=(nlyb-noy1)-(0)   +1
	if (jblock.eq.npy)   lsizes(3)=(nlyb+0)   -(noy1)+1
	if ( npy  .eq. 1 )   lsizes(3)=(ngy+1)

		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
		ke = ks+lsizes(1)-1;    ie = is+lsizes(2)-1;    je = js+lsizes(3)-1	!ALC_FIX
		OutBuffer = p(ks:ke, is:ie, js:je)					!ALC_FIX
		memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngz+1 , nlx+1 , nly+1 /)
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
	CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)							!ALC_FIX

       !==========================================================
       !     SYNC ALL THE PROCESSORS
       !----------------------------------------------------------
        CALL MPI_BARRIER(icomm_grid,IERR)
        if (irank.eq.iroot) PRINT *,'restart file(1p) dumped from processor=',irank

       !==========================================================
       !     Dump Out Pressure (phatr) to check I/O working properly
       !----------------------------------------------------------
       !
       !  NAME='pressure_ph.'
       !  WRITE(FN2,'(A12,A5,I1)') NAME,'00000',irank
       !  write(*,*) FN2
       !  open(29,file=FN2,form='unformatted')
       !  write(29) phatr
       !  close(29)
        !=======================================================
        !                  Phatr
        !=======================================================
        NAME='pressure_ph.'
        call OutFileName()
        if (irank.eq.iroot) write(*,*) FN2

        call MPI_FILE_OPEN(icomm_grid, FN2, &
                           MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, ierr)

        !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	global_indices = (/ 0 , iTmin_1(iblock) , jTmin_1(jblock) /)
	if (iblock.eq. 1 ) global_indices(2) = 0
	if (jblock.eq. 1 ) global_indices(3) = 0

	!ALC_BUG	local_indices  = (/ 0 , 1 ,  1 /)
	!ALC_BUG	if (iblock.eq. 1 ) local_indices(2) = 0
	!ALC_BUG	if (jblock.eq. 1 ) local_indices(3) = 0
			local_indices  = (/ 0 , 0 ,  0 /)					!ALC_FIX
			ks = 1;			is = 1;			js = 1			!ALC_FIX
			if (iblock.eq. 1 ) is = 0						!ALC_FIX
			if (jblock.eq. 1 ) js = 0						!ALC_FIX

	!--------------- phatr ------------------
        global_cnt = 0
        offset     = 0

	gsizes = (/ ngzm , ngx+1 , ngy+1 /)
	lsizes = (/ ngzm , nixp  , niyp  /)
	if (iblock.eq. 1 )   lsizes(2) = nixp + 1
	if (iblock.eq.npx)   lsizes(2) = nixp + 1
	if ( npx  .eq. 1 )   lsizes(2) = (ngx+1)
	if (jblock.eq. 1 )   lsizes(3) = niyp + 1
	if (jblock.eq.npy)   lsizes(3) = niyp + 1
	if ( npy  .eq. 1 )   lsizes(3) = (ngy+1)

			allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))			!ALC_FIX
			ke = ks+lsizes(1)-1;	ie = is+lsizes(2)-1;	je = js+lsizes(3)-1	!ALC_FIX
			OutBuffer = phatr(ks:ke, is:ie, js:je)					!ALC_FIX
			memsizes = lsizes							!ALC_FIX
	CALL Create_commit_fileview()
	!ALC_BUG	memsizes = (/ ngzm , nixp+2 , niyp+2 /)
	CALL Create_commit_subarray()
	CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		deallocate(OutBuffer)							!ALC_FIX

        CALL MPI_FILE_CLOSE(fh, ierr)

	call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
	call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0
return
end



subroutine Create_commit_fileview()
	use Output_mod

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

subroutine Create_commit_subarray()
	use Output_mod

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

subroutine OutFileName()
        use Output_mod

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
