!====================================================================
!     ******************************************************************
!     *                                                                *
!     *   Output Y-planes of the instantaneous field                   *
!     *   Each plane is (ngx,ngz):                                     * 
!     *                   No halos in order to ease post-processing.   *
!     *   PARALLEL I/O OF MPI-2                                        *
!     *                                                                *
!     ******************************************************************
!     *   This code only works when each processor makes the same      *
!     *   number of write, or file access, calls -- otherwise          *
!     *   the code doesn't work                                        *
!     ******************************************************************
!=======================================================================


module Output_SubDomain_mod
      use messenger
      use data_export
      use computation_parameters
      implicit none
	include "mpif.h"

	INTEGER , DIMENSION(nproc) :: DataPnts
	INTEGER                    :: global_cnt, FloatSize, fh
	INTEGER                    :: status(MPI_STATUS_SIZE)
	INTEGER (kind=MPI_OFFSET_KIND) :: offset
	INTEGER, DIMENSION(3)      :: gsizes, lsizes, memsizes
	INTEGER, DIMENSION(3)	   :: global_indices, local_indices
	INTEGER					   :: filetype, memtype
	CHARACTER(LEN=18)          :: FN2
	CHARACTER(LEN=12)          :: NAME

	INTEGER, PARAMETER            :: n_Jplanes = 10
	! INTEGER, DIMENSION(n_Jplanes), PARAMETER :: Jplane = (/ 5 , 10 , 15 , ngy-4 , ngy-9 , ngy-14 /)
	INTEGER, DIMENSION(n_Jplanes), PARAMETER :: Jplane = (/ 15 , 25 , 35 , 45 , 55 , 60 , 65 , 70 , 75 , 80 /)

	INTEGER	:: WRITER_COMM

	INTEGER :: MEM_FLAG = 0
	INTEGER :: FILE_FLAG = 0

	REAL aan

end module

!============================================================================================
! Write subdomain with co-located velocities

subroutine WriteSubDomain(FNAME, istart,iend,iskip, jstart,jend,jskip, kstart,kend,kskip)
use Output_SubDomain_mod

	implicit none
	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	CHARACTER(LEN=12)          :: FNAME
	integer	:: istart, iend, iskip
	integer	:: jstart, jend, jskip
	integer	:: kstart, kend, kskip
	double precision, allocatable	:: OutBuffer(:,:,:)

	INTEGER	:: loc_istart, loc_iend
	INTEGER	:: loc_jstart, loc_jend
	INTEGER	:: my_istart,  my_iend,  global_indx,  loc_indx
	INTEGER	:: my_jstart,  my_jend,  global_jndx,  loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	::  i, j, k
	INTEGER	:: ii,jj,kk

	NAME = FNAME
	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output box on processes) -------
	!--------------------------------------------------------------------------------------------
	!my_istart = ibmin	;	if (iblock ==  1 ) my_istart =  1
	!my_iend   = ibmax	;	if (iblock == npx) my_iend   = ngx
	!my_jstart = jbmin	;	if (jblock ==  1 ) my_jstart =  1
	!my_jend   = jbmax	;	if (jblock == npy) my_jend   = ngy

	my_istart = iTmin_1(iblock)	;	if (iblock ==  1 ) my_istart =  1
	my_iend   = iTmax_1(iblock)	;	if (iblock == npx) my_iend   = ngx
	my_jstart = jTmin_1(jblock)	;	if (jblock ==  1 ) my_jstart =  1
	my_jend   = jTmax_1(jblock)	;	if (jblock == npy) my_jend   = ngy

	!------ Create WRITER_COMM -------
	FLAG_RPOC_WRITES = 1
	if (iend.lt.my_istart  .or.  istart.gt.my_iend  .or.  jend.lt.my_jstart  .or.  jstart.gt.my_jend) then
		FLAG_RPOC_WRITES = 0
	else
		FLAG_RPOC_WRITES = 1
	end if

	CALL MPI_COMM_SPLIT(CFD_COMM, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
		call SubDomain_FileName()

		call MPI_FILE_OPEN(WRITER_COMM, trim(prefix_dir)//FN2, &
                           	MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , 0 , 0 /)
		local_indices  = (/  0  , 0 , 0 /)

		!------- x-direction -------
		if (my_istart .le. istart) then
			global_indices(2) = 0
			loc_istart	  = imap_1(istart)
		else
			global_indices(2) = ( my_istart - istart -1 )/iskip + 1
			!loc_istart	  = imap_1( my_istart + iskip - mod(my_istart-istart,iskip) )
			loc_istart	  = imap_1( istart + ((my_istart-istart-1)/iskip)*iskip + iskip)
		end if

		if (my_iend .ge. iend) then
			loc_iend	  = imap_1(iend)
		else
			loc_iend	  = imap_1( my_iend - mod(my_iend-istart,iskip) )
		end if

		!------- y-direction -------
		if (my_jstart .le. jstart) then
			global_indices(3) = 0
			loc_jstart	  = jmap_1(jstart)
		else
			global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
			!loc_jstart	  = jmap_1( my_jstart + jskip - mod(my_jstart-jstart,jskip) )
			loc_jstart	  = jmap_1( jstart + ((my_jstart-jstart-1)/jskip)*jskip + jskip)
		end if

		if (my_jend .ge. jend) then
			loc_jend	  = jmap_1(jend)
		else
			loc_jend	  = jmap_1( my_jend - mod(my_jend-jstart,jskip) )
		end if

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)
		memsizes = lsizes
		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))

		global_cnt = 0
		offset     = 0

		!------------------------ uc ---------------------------
		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

			OutBuffer (kk,ii,jj)= 0.25*( uc( k-1, i, j-1)+uc( k  , i, j-1) &
						    			+uc( k  , i, j  )+uc( k-1, i, j  )  )

		end do
		end do
		end do

		CALL Create_commit_SubDomain_fileview()
		CALL Create_commit_SubDomain_subarray() 
		CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
		
       		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
       		offset = global_cnt * FloatSize

		!------------------------ vc ---------------------------
		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

			OutBuffer(kk,ii,jj) = 0.25*( vc( k-1, i-1, j)+vc( k  , i-1, j) &
						    +vc( k  , i  , j)+vc( k-1, i  , j)  )

		end do
		end do
		end do

		CALL Create_commit_SubDomain_fileview()
		CALL Create_commit_SubDomain_subarray() 
		CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
		
       		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
       		offset = global_cnt * FloatSize

		!------------------------ wc ---------------------------
		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

			OutBuffer(kk,ii,jj) = 0.25*( wc( k, i-1, j-1)+wc( k, i-1, j  ) &
						    +wc( k, i  , j  )+wc( k, i  , j-1)  )

		end do
		end do
		end do

		CALL Create_commit_SubDomain_fileview()
		CALL Create_commit_SubDomain_subarray() 
		CALL MPI_FILE_WRITE(fh, OutBuffer, 1, memtype, status, ierr)
		
       		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
       		offset = global_cnt * FloatSize

		!------ close file and clean up -----
	
		CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if
       	!==========================================================
       	!     SYNC ALL THE PROCESSORS
       	!----------------------------------------------------------
       	CALL MPI_BARRIER(icomm_grid,IERR)
       	if (irank.eq.iroot) PRINT *,'Snapshot (Y-plane) dumped by processes'

        CALL MPI_COMM_FREE(WRITER_COMM, ierr)
return
end subroutine WriteSubDomain

!============================================================================================
! Write subdomain with cell centered velocities and pressures

subroutine WriteSubDomain_cc(FNAME, istart,iend,iskip, jstart,jend,jskip, kstart,kend,kskip)
	use Output_SubDomain_mod
	use Stress_strain

	implicit none
	!----------------------------------------------------
	! (OutBuffer) is allocated of exact size as input data.
	! It circumvents a bug in ALC MPI-IO
	!----------------------------------------------------
	CHARACTER(LEN=12),intent(in)	:: FNAME
	integer,intent(in)				:: istart, iend, iskip
	integer,intent(in)				:: jstart, jend, jskip
	integer,intent(in)				:: kstart, kend, kskip
	double precision, allocatable	:: OutBuffer(:,:,:)
	real(kind(0.d0)),dimension(:,:,:,:,:),allocatable	:: stress

	INTEGER	:: loc_istart, loc_iend
	INTEGER	:: loc_jstart, loc_jend
	INTEGER	:: my_istart,  my_iend,  global_indx,  loc_indx
	INTEGER	:: my_jstart,  my_jend,  global_jndx,  loc_jndx
	INTEGER	:: FLAG_RPOC_WRITES = 0

	INTEGER	::  i, j, k
	INTEGER	::  ii,jj,kk
	INTEGER	::  ixyz, jxyz

	NAME = FNAME
	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, FloatSize, ierr)

	!=======================================================
	!		CREATE (WRITER_COMM):
	!  The communicator which writes the output plane
	!=======================================================
	!--------- DEFINE PROCESS LIMITS (in order to locate desired output box on processes) -------
	!--------------------------------------------------------------------------------------------
	my_istart = iTmin_1(iblock)	;	if (iblock ==  1 ) my_istart =  1
	my_iend   = iTmax_1(iblock)	;	if (iblock == npx) my_iend   = ngx
	my_jstart = jTmin_1(jblock)	;	if (jblock ==  1 ) my_jstart =  0
	my_jend   = jTmax_1(jblock)	;	if (jblock == npy) my_jend   = ngy

	!------ Create WRITER_COMM -------
	FLAG_RPOC_WRITES = 1
	if (iend.lt.my_istart  .or.  istart.gt.my_iend  .or.  jend.lt.my_jstart  .or.  jstart.gt.my_jend) then
		FLAG_RPOC_WRITES = 0
	else
		FLAG_RPOC_WRITES = 1
	end if

	CALL MPI_COMM_SPLIT(CFD_COMM, FLAG_RPOC_WRITES, 0, WRITER_COMM, ierr)

	if (FLAG_RPOC_WRITES .eq. 1) then
		!=======================================================
		!         CARTESIAN VELOCITIES
		!=======================================================
		!---------- Open File for writing -------
		call SubDomain_FileName()

		call MPI_FILE_OPEN(WRITER_COMM, trim(prefix_dir)//FN2, &
                           	MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                           	MPI_INFO_NULL, fh, ierr)

		!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
		global_indices = (/  0  , 0 , 0 /)
		local_indices  = (/  0  , 0 , 0 /)

		!------- x-direction -------
		if (my_istart .le. istart) then
			global_indices(2) = 0
			loc_istart	  = imap_1(istart)
		else
			global_indices(2) = ( my_istart - istart -1 )/iskip + 1
			loc_istart	  = imap_1( istart + ((my_istart-istart-1)/iskip)*iskip + iskip)
		end if

		if (my_iend .ge. iend) then
			loc_iend	  = imap_1(iend)
		else
			loc_iend	  = imap_1( my_iend - mod(my_iend-istart,iskip) )
		end if

		!------- y-direction -------
		if (my_jstart .le. jstart) then
			global_indices(3) = 0
			loc_jstart	  = jmap_1(jstart)
		else
			global_indices(3) = ( my_jstart - jstart -1 )/jskip + 1
			loc_jstart	  = jmap_1( jstart + ((my_jstart-jstart-1)/jskip)*jskip + jskip)
		end if

		if (my_jend .ge. jend) then
			loc_jend	  = jmap_1(jend)
		else
			loc_jend	  = jmap_1( my_jend - mod(my_jend-jstart,jskip) )
		end if

		!--------- DEFINE SIZES (GLOBAL & LOCAL PORTION OF ARRAY) -------
		gsizes = (/  (kend-kstart)/kskip +1 ,         (iend-istart)/iskip +1 ,         (jend-jstart)/jskip +1  /)
		lsizes = (/  (kend-kstart)/kskip +1 , (loc_iend-loc_istart)/iskip +1 , (loc_jend-loc_jstart)/jskip +1  /)
		memsizes = lsizes
		allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))

		global_cnt = 0
		offset     = 0

		!------------------------ uc ---------------------------
		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

			OutBuffer (kk,ii,jj)= 0.5*( uc(k,i,j)+uc(k,i+1,j) )

		end do
		end do
		end do

		CALL Create_commit_SubDomain_fileview()
		CALL Create_commit_SubDomain_subarray() 
		CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		
		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
 		offset = global_cnt * FloatSize

		!------------------------ vc ---------------------------
		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

			OutBuffer(kk,ii,jj) = 0.5*( vc(k,i,j  )+vc(k,i,j+1) )

		end do
		end do
		end do

		CALL Create_commit_SubDomain_fileview()
		CALL Create_commit_SubDomain_subarray() 
		CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		
		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
		offset = global_cnt * FloatSize

		!------------------------ wc ---------------------------
		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

			OutBuffer(kk,ii,jj) = 0.5*( wc(k  ,i,j)+wc(k+1,i,j) )

		end do
		end do
		end do

		CALL Create_commit_SubDomain_fileview()
		CALL Create_commit_SubDomain_subarray() 
		CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		
 		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
 		offset = global_cnt * FloatSize

		!------------------------  P ---------------------------
		do j = loc_jstart, loc_jend, jskip
		do i = loc_istart, loc_iend, iskip
		do k =     kstart,     kend, kskip
			jj = (j-loc_jstart)/jskip + 1
			ii = (i-loc_istart)/iskip + 1
			kk = (k-    kstart)/kskip + 1

			OutBuffer (kk,ii,jj)= p(k,i,j)

		end do
		end do
		end do

		CALL Create_commit_SubDomain_fileview()
		CALL Create_commit_SubDomain_subarray() 
		CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
		
		global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
		offset = global_cnt * FloatSize

		!---------------------- Stress -------------------------
		call Evaluate_stress(uc,vc,wc,P,stress)
		! Write all nine components of stress
		do jxyz = 1,3
		do ixyz = 1,3

			do j = loc_jstart, loc_jend, jskip
			do i = loc_istart, loc_iend, iskip
			do k =     kstart,     kend, kskip
				jj = (j-loc_jstart)/jskip + 1
				ii = (i-loc_istart)/iskip + 1
				kk = (k-    kstart)/kskip + 1

				OutBuffer (kk,ii,jj)= stress(ixyz,jxyz,k,i,j)

			end do
			end do
			end do

			CALL Create_commit_SubDomain_fileview()
			CALL Create_commit_SubDomain_subarray() 
			CALL MPI_FILE_WRITE_ALL(fh, OutBuffer, 1, memtype, status, ierr)
			
			global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
			offset = global_cnt * FloatSize

		enddo
		enddo

		!------ close file and clean up -----
	
		CALL MPI_FILE_CLOSE(fh, ierr)
		deallocate(OutBuffer)

		call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
		call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

	end if

	!==========================================================
	!     SYNC ALL THE PROCESSORS
	!----------------------------------------------------------
	CALL MPI_BARRIER(icomm_grid,IERR)
	if (irank.eq.iroot) PRINT *,'Snapshot (Y-plane) dumped by processes'

    CALL MPI_COMM_FREE(WRITER_COMM, ierr)

return
end subroutine WriteSubDomain_cc




subroutine Create_commit_SubDomain_fileview()
	use Output_SubDomain_mod

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

subroutine Create_commit_SubDomain_subarray()
	use Output_SubDomain_mod

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

subroutine SubDomain_FileName()
        use Output_SubDomain_mod

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
