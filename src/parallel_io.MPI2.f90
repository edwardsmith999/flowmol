!======================================================================

! Parallel i/o

! --- INPUT ROUTINES ---
! setup_inputs              Read input file
! setup_restart_inputs      Read restart files and input file
! setup_restart_microstate  Read Initial configuration from restart file

! --- OUTPUT ROUTINES ---
! parallel_io_final_state   Write final configuration and simulation details for restart
! parallel_io_vmd           Output for VMD
! parallel_io_vmd_sl        Output for VMD with separate solid/lquid regions
! parallel_io_vmd_halo      Output Halos
! CV AVERAGING
! mass_averaging            slice/CV
! cumulative_mass           slice/CV
! mass_snapshot             CV
! momentum_averaging        slice/CV
! cumulative_momentum       slice/CV
! momentum_snapshot         CV
! pressure_averaging        domain/CV
! cumulative_pressure       domain/CV
! FLUX AVERAGING
! mass_flux_averaging       CV_surface
! cumulative_mass_flux      CV_surface
! momentum_flux_averaging   plane(MOP)/CV_surface
! cumulative_momentum_flux  plane(MOP)/CV_surface
! surface_pressure          plane(MOP)/CV_surface
! cumulative_pressure       plane(MOP)/CV_surface
! 
!======================================================================

module module_parallel_io
    use mpi
    use computational_constants_MD
    use physical_constants_MD
    use arrays_MD
    use polymer_info_MD
    use concentric_cylinders, only: cyl_file
    use shear_info_MD
    use calculated_properties_MD
    use messenger, only : MD_COMM
    use interfaces
    implicit none

    integer     :: restartfileid, fileid, fileidtrue !File name used for parallel i/o

    !Generic interface so write arrays can be used with both integers and reals
    interface write_arrays
        module procedure iwrite_arrays_1, rwrite_arrays_1, iwrite_arrays, rwrite_arrays 
    end interface write_arrays
    private  iwrite_arrays_1, rwrite_arrays_1, iwrite_arrays, rwrite_arrays
    
    interface write_zplane 
        module procedure iwrite_zplane_1, iwrite_zplane, rwrite_zplane_1, rwrite_zplane
    end interface write_zplane
    private iwrite_zplane_1, iwrite_zplane, rwrite_zplane_1, rwrite_zplane

contains


!====================================================================
!           Array writing subroutines
!--------------------------------------------------------------------


! --- Write integer arrays ---

!   1D array wrapper
subroutine iwrite_arrays_1(temp,nresults,outfile,outstep)

    integer, intent(in)                     :: nresults,outstep
    integer, dimension(:,:,:),intent(in)    :: temp
    character(*),intent(in)                 :: outfile

    integer, dimension(:,:,:,:),allocatable :: some_array

    allocate(some_array(size(temp,1),size(temp,2),size(temp,3),1))
    some_array(:,:,:,1) = temp(:,:,:)
    call iwrite_arrays(some_array,nresults,outfile,outstep)
    deallocate(some_array)

end subroutine iwrite_arrays_1

!   Main iwrite Routine
subroutine iwrite_arrays(some_array,nresults,outfile,outstep)
    use messenger, only : icomm_grid
    use computational_constants_MD, only : separate_outfiles
    use librarymod, only : get_Timestep_FileName

    integer, intent(in)                     :: nresults,outstep
    integer, dimension(:,:,:,:),intent(in)  :: some_array!(nx,ny,nz,nresults)
    character(*),intent(in)                 :: outfile

    integer                                 :: n, fh, outstep_
    integer                                 :: MEM_FLAG = 0
    integer                                 :: FILE_FLAG = 0
    integer                                 :: int_size,datatype
    integer                                 :: status(mpi_status_size)
    integer                                 :: filetype, memtype
    integer (kind=MPI_offset_kind)          :: offset, global_cnt
    integer, dimension(3)                   :: gsizes, lsizes, memsizes
    integer, dimension(3)                   :: global_indices, local_indices
    integer, allocatable,dimension(:,:)     :: proc_lsizes 
    integer, allocatable,dimension(:,:,:)   :: OutBuffer
    character(200)                          :: outfile_t

    datatype = MPI_INTEGER
    call MPI_TYPE_SIZE(datatype, int_size, ierr)

    !Error if outstep less than one
    if (outstep .le. 0) then
        !print*, "Warning -- outstep ", outstep, " for filename ", & 
        !        trim(outfile) , " results in ", trim(outfile_t)
        outstep_ = 1
        !call error_abort("Error in write_arrays -- requires outstep > 1 ")
    else
        outstep_ = outstep
    endif

    !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
    !  Note:  MPI assumes here that numbering starts from zero
    !  Since numbering starts from (one), subtract (one) from every index
    !-------------------------------------------------------
    if (separate_outfiles) then
        !Either write a separate output files for each timestep
        call get_Timestep_FileName(outstep_-1,outfile,outfile_t)
        global_cnt = 0; offset = 0
    else
        !or a single file for whole run
        outfile_t = outfile
        global_cnt  = (outstep_-1)*gnbins(1)*gnbins(2)*gnbins(3)*nresults
        offset      = global_cnt * int_size
    endif

    gsizes      = gnbins
    lsizes      = nbins
    local_indices(:) = (/  0  , 0 , 0 /)
    !Calculate global_indices
    global_indices(:)= (/  0  , 0 , 0 /)
    allocate(proc_lsizes(3,nproc))
    call globalGather(lsizes,proc_lsizes,3)
    global_indices(1) = sum(proc_lsizes(1,1:iblock-1))
    global_indices(2) = sum(proc_lsizes(2,1:jblock-1))
    global_indices(3) = sum(proc_lsizes(3,1:kblock-1))
    deallocate(proc_lsizes)

    !Allocate ouput buffer
    allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))
    memsizes = lsizes

    !Open file to write
    call MPI_file_open(icomm_grid, outfile_t, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                                MPI_INFO_NULL, fh, ierr)

    do n =1,nresults
        !Copy to outbuffer
        !print*, n,1+nhb(1),nbins(1)+nhb(1),1+nhb(2),nbins(2)+nhb(2),1+nhb(3),nbins(3)+nhb(3)
        OutBuffer =  some_array(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),n)
        !Update array datatype and reset fileview to correct section of array
        CALL Create_commit_fileview(gsizes,lsizes,global_indices,offset,datatype,FILE_FLAG,filetype,fh)
        !Update local array datatype to ignore halo cells
        CALL Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
        !Write to file
        CALL MPI_file_write_all(fh, OutBuffer, 1, memtype, status, ierr)

        !Calculate global cnt offset
        global_cnt = global_cnt + gnbins(1)*gnbins(2)*gnbins(3)
        offset = global_cnt * int_size
    
    enddo

    deallocate(OutBuffer)
    CALL MPI_FILE_CLOSE(fh, ierr)

    !==========================================================
    !     FREE DATA TYPES
    !----------------------------------------------------------
    call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
    call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

end subroutine iwrite_arrays


! --- Double precision arrays ---

!   1D array wrapper
subroutine rwrite_arrays_1(temp,nresults,outfile,outstep)

    integer, intent(in)                             :: nresults,outstep
    real(kind(0.d0)), dimension(:,:,:),intent(in)   :: temp
    character(*),intent(in)                         :: outfile

    real(kind(0.d0)), dimension(:,:,:,:),allocatable    :: some_array

    allocate(some_array(size(temp,1),size(temp,2),size(temp,3),1))
    some_array(:,:,:,1) = temp(:,:,:)
    call rwrite_arrays(some_array,nresults,outfile,outstep)
    deallocate(some_array)

end subroutine rwrite_arrays_1

!   Main rwrite Routine
subroutine rwrite_arrays(some_array,nresults,outfile,outstep)
    use messenger, only : icomm_grid
    use computational_constants_MD, only : separate_outfiles
    use librarymod, only : get_Timestep_FileName

    integer, intent(in)                             :: nresults,outstep
    real(kind(0.d0)), dimension(:,:,:,:),intent(in) :: some_array
    character(*),intent(in)                         :: outfile

    integer                             :: n, fh, outstep_
    integer                             :: MEM_FLAG = 0
    integer                             :: FILE_FLAG = 0
    integer                             :: datatype
    integer                             :: status(mpi_status_size)
    integer                             :: filetype, memtype, dp_size_temp
    integer (kind=MPI_offset_kind)      :: offset, global_cnt, dp_size
    integer, dimension(3)               :: gsizes, lsizes, memsizes
    integer, dimension(3)               :: global_indices, local_indices
    integer, dimension(:,:),allocatable :: proc_lsizes 
    real(kind(0.d0)), allocatable,dimension(:,:,:)      :: OutBuffer
    character(200)                      :: outfile_t

    datatype = MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(datatype, dp_size_temp, ierr)
    dp_size = dp_size_temp

    !Error if outstep less than one
    if (outstep .le. 0) then
        !print*, "Warning -- outstep ", outstep, " for filename ", & 
        !        trim(outfile) , " results in ", trim(outfile_t)
        outstep_ = 1
        !call error_abort("Error in write_arrays -- requires outstep > 1 ")
    else
        outstep_ = outstep
    endif

    !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
    !  Note:  MPI assumes here that numbering starts from zero
    !  Since numbering starts from (one), subtract (one) from every index
    !-------------------------------------------------------
    if (separate_outfiles) then
        !Either write a separate output files for each timestep
        call get_Timestep_FileName(outstep_-1,outfile,outfile_t)
        global_cnt = 0; offset = 0
    else
        !or a single file for the whole run
        outfile_t = outfile
        global_cnt  = (outstep_-1)*gnbins(1)*gnbins(2)*gnbins(3)*nresults !Global number of items written so far
        offset      = global_cnt * dp_size                               !Bytes written already 
    endif

    gsizes      = gnbins                                             !Global "sizes", i.e. bins (need a better name for this)
    lsizes      = nbins                                              !Local "sizes", i.e. bins (need a better name for this)
    local_indices(:) = (/  0  , 0 , 0 /)                             !Not sure, goes into MPI_TYPE_CREATE_SUBARRAY !todo
    !Calculate global_indices
    global_indices(:)= (/  0  , 0 , 0 /)                             !Same, goes into MPI_TYPE_CREATE_SUBARRAY
    allocate(proc_lsizes(3,nproc))                                   !Number of bins on each processor 
    call globalGather(lsizes,proc_lsizes,3)                          !(populated)
    global_indices(1) = sum(proc_lsizes(1,1:iblock-1))               !Seems to be a set of offsets in x,y,z
    global_indices(2) = sum(proc_lsizes(2,1:jblock-1))
    global_indices(3) = sum(proc_lsizes(3,1:kblock-1))
    deallocate(proc_lsizes)

    !Allocate ouput buffer
    allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))               !Buffer that's written
    memsizes = lsizes                                                !Not sure why redefining?

    !Open file to write
    call MPI_file_open(icomm_grid, outfile_t, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                                MPI_INFO_NULL, fh, ierr)

    do n =1,nresults
        !Copy to outbuffer
        OutBuffer =  some_array(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),n)
        !Update array datatype and reset fileview to correct section of array
        CALL Create_commit_fileview(gsizes,lsizes,global_indices,offset,datatype,FILE_FLAG,filetype,fh)
        !Update local array datatype to ignore halo cells
        CALL Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
        !print'(a,9i6,2i18)',"rwrite",irank, n, shape(OutBuffer), shape(some_array),global_cnt,offset
        !Write to file
        CALL MPI_file_write_all(fh, OutBuffer, 1, memtype, status, ierr)

        !Calculate global cnt offset
        global_cnt = global_cnt + gnbins(1)*gnbins(2)*gnbins(3)
        offset = global_cnt * dp_size
    
    enddo

    deallocate(OutBuffer)
    CALL MPI_FILE_CLOSE(fh, ierr)

    !==========================================================
    !     FREE DATA TYPES
    !----------------------------------------------------------
    CALL MPI_BARRIER(icomm_grid,IERR)
    call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
    call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

end subroutine rwrite_arrays

subroutine rwrite_zplane_1(temp,nresults,outfile,outstep)

    integer, intent(in) :: nresults,outstep
    real(kind(0.d0)), dimension(:,:,:), intent(in) :: temp 
    character(*),intent(in) :: outfile

    real(kind(0.d0)), dimension(:,:,:,:), allocatable :: cpol_array

    allocate(cpol_array(size(temp,1),size(temp,2),size(temp,3),1))
    cpol_array(:,:,:,1) = temp(:,:,:)
    call rwrite_zplane(cpol_array,nresults,outfile,outstep)
    deallocate(cpol_array)

end subroutine rwrite_zplane_1

! Write array to file
subroutine rwrite_zplane(cpol_array,nresults,outfile,outstep)
    use mpi
    use messenger, only : icomm_xyz, kblock
    use concentric_cylinders, only: gcpol_bins, cpol_bins, cpol_nhbz
    use computational_constants_MD, only: npz, separate_outfiles
    use librarymod, only : get_Timestep_FileName
    implicit none

    integer, intent(in) :: nresults,outstep
    real(kind(0.d0)), dimension(:,:,:,:),intent(in) :: cpol_array
    character(*),intent(in) :: outfile

    integer :: n, fh
    integer :: i, j, k, l
    integer :: dp_size,datatype
    integer :: status(mpi_status_size)
    integer :: filetype, memtype
    integer (kind=MPI_offset_kind) :: offset, global_cnt
    integer, dimension(3) :: gsizes, lsizes, memsizes
    integer, dimension(3) :: global_indices, local_indices
    integer, dimension(:,:),allocatable :: kblock_lsizes 
    real(kind(0.d0)), allocatable,dimension(:,:,:) :: OutBuffer
    character(200) :: outfile_t

    datatype = MPI_DOUBLE_PRECISION
    call MPI_Type_size(datatype, dp_size, ierr)

    if (separate_outfiles) then
        !Either write a separate output files for each timestep
        call get_Timestep_FileName(outstep-1, outfile, outfile_t)
        global_cnt = 0; offset = 0
    else
        !or a single file for the whole run
        outfile_t = outfile
        ! Define offsets, etc. Remember MPI is written in C so indices are
        ! counted from 0, rather than 1 as in Fortran.
        global_cnt  = int8(int8((outstep-1))*&
                      int8(product(int8(gcpol_bins(:))))*int8(nresults))
        offset      = int8(global_cnt) * int8(dp_size)
    endif

    gsizes      = gcpol_bins
    lsizes      = cpol_bins
    local_indices(:) = (/ 0 , 0 , 0 /)
    global_indices(:)= (/ 0 , 0 , 0 /)

    ! Number of bins on each kblock 
    allocate(kblock_lsizes(3,npz))
    call SubcommGather(lsizes,kblock_lsizes,3,3,npz)
    global_indices(1) = 0
    global_indices(2) = 0 
    global_indices(3) = sum(kblock_lsizes(3,1:kblock-1))
    deallocate(kblock_lsizes)

    !Allocate ouput buffer
    allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))
    memsizes = lsizes

    ! Open file, create and commit MPI subarray and datatype 
    call MPI_File_open(icomm_xyz(3), outfile_t, &
                       MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                       MPI_INFO_NULL, fh, ierr)
    call MPI_Type_create_subarray(3, gsizes, lsizes, global_indices, &
                                  MPI_ORDER_FORTRAN, datatype, &
                                  filetype, ierr)
    call MPI_Type_commit(filetype, ierr)
    call MPI_Type_create_subarray(3, memsizes, lsizes, local_indices,&
                                  MPI_ORDER_FORTRAN, datatype, memtype,&
                                  ierr)
    call MPI_Type_commit(memtype, ierr)

    do n =1,nresults

        OutBuffer = cpol_array(:,:,cpol_nhbz+1:cpol_bins(3)+cpol_nhbz,n)
        
        call MPI_File_set_view(fh, offset, datatype, filetype, &
                               'native', MPI_INFO_NULL, ierr)
        call MPI_File_write_all(fh, OutBuffer, 1, memtype, status, ierr)

        !Calculate global count offset
        global_cnt = int8(global_cnt) + int8(product(int8(gcpol_bins(:))))
        offset = int8(global_cnt) * int8(dp_size)
    
    enddo

    deallocate(OutBuffer)
    call MPI_File_close(fh, ierr)

    ! Free data types
    call MPI_Barrier(icomm_xyz(3),ierr)
    call MPI_Type_free(memtype,ierr)
    call MPI_Type_free(filetype,ierr)

end subroutine rwrite_zplane

subroutine iwrite_zplane_1(temp,nresults,outfile,outstep)

    integer, intent(in) :: nresults,outstep
    integer, dimension(:,:,:), intent(in) :: temp 
    character(*),intent(in) :: outfile

    integer, dimension(:,:,:,:), allocatable :: cpol_array

    allocate(cpol_array(size(temp,1),size(temp,2),size(temp,3),1))
    cpol_array(:,:,:,1) = temp(:,:,:)
    call iwrite_zplane(cpol_array,nresults,outfile,outstep)
    deallocate(cpol_array)

end subroutine iwrite_zplane_1

subroutine iwrite_zplane(cpol_array,nresults,outfile,outstep)
    use mpi
    use messenger, only : icomm_xyz, kblock
    use concentric_cylinders, only: gcpol_bins, cpol_bins, cpol_nhbz
    use computational_constants_MD, only: npz, separate_outfiles
    use librarymod, only : get_Timestep_FileName
    implicit none

    integer, intent(in) :: outstep
    integer, intent(in) :: nresults 
    integer, dimension(:,:,:,:),intent(in) :: cpol_array
    character(*),intent(in) :: outfile

    integer :: n, fh, ierr
    integer :: MEM_FLAG = 0
    integer :: FILE_FLAG = 0
    integer :: int_size,datatype
    integer :: status(mpi_status_size)
    integer :: filetype, memtype
    integer (kind=MPI_offset_kind) :: offset,global_cnt
    integer, dimension(3) :: gsizes, lsizes, memsizes
    integer, dimension(3) :: global_indices, local_indices
    integer, dimension(:,:),allocatable :: kblock_lsizes 
    integer, allocatable,dimension(:,:,:) :: OutBuffer
    character(200) :: outfile_t

    datatype = MPI_INTEGER
    call MPI_TYPE_SIZE(datatype, int_size, ierr)

    if (separate_outfiles) then
        !Either write a separate output files for each timestep
        call get_Timestep_FileName(outstep-1, outfile, outfile_t)
        global_cnt = 0; offset = 0
    else
        !or a single file for the whole run
        outfile_t = outfile
        ! Define offsets, etc. Remember MPI is written in C so indices are
        ! counted from 0, rather than 1 as in Fortran.
        global_cnt  = int8(int8((outstep-1))*&
                      int8(product(int8(gcpol_bins(:))))*int8(nresults))
        offset      = int8(global_cnt) * int8(int_size)
    endif

    gsizes      = gcpol_bins
    lsizes      = cpol_bins
    local_indices(:) = (/ 0 , 0 , 0 /)
    global_indices(:)= (/ 0 , 0 , 0 /)

    !Number of bins on each kblock processor group 
    allocate(kblock_lsizes(3,npz))
    call SubcommGather(lsizes,kblock_lsizes,3,3,npz)
    global_indices(1) = 0
    global_indices(2) = 0 
    global_indices(3) = sum(kblock_lsizes(3,1:kblock-1))
    deallocate(kblock_lsizes)

    allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))
    memsizes = lsizes

    ! Open file, create and commit MPI subarray and datatype 
    call MPI_file_open(icomm_xyz(3), outfile_t, &
                       MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                       MPI_INFO_NULL, fh, ierr)

    call MPI_Type_create_subarray(3, gsizes, lsizes, global_indices, &
                                  MPI_ORDER_FORTRAN, datatype, &
                                  filetype, ierr)
    call MPI_Type_commit(filetype, ierr)
    call MPI_Type_create_subarray(3, memsizes, lsizes, local_indices,&
                                  MPI_ORDER_FORTRAN, datatype, memtype,&
                                  ierr)
    call MPI_Type_commit(memtype, ierr)

    do n =1,nresults

        OutBuffer = cpol_array(:,:,cpol_nhbz+1:cpol_bins(3)+cpol_nhbz,n)

        call MPI_File_set_view(fh, offset, datatype, filetype, &
                               'native', MPI_INFO_NULL, ierr)
        call MPI_File_write_all(fh, OutBuffer, 1, memtype, status, ierr)

        !Calculate global count offset
        global_cnt = int8(global_cnt) + int8(product(int8(gnbins(:))))
        offset = int8(global_cnt) * int8(int_size)
    
    enddo

    deallocate(OutBuffer)
    CALL MPI_FILE_CLOSE(fh, ierr)

    ! Free data types
    CALL MPI_BARRIER(icomm_xyz(3),IERR)
    call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
    call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

end subroutine iwrite_zplane

subroutine Create_commit_fileview(gsizes,lsizes,global_indices,offset,datatype,FILE_FLAG,filetype,fh)
    implicit none

    integer,intent(inout)                       :: filetype
    integer,intent(in)                          :: datatype,fh
    integer (kind=MPI_offset_kind),intent(in)   :: offset
    integer, dimension(3),intent(in)            :: gsizes, lsizes, global_indices
    integer,intent(inout)                       :: FILE_FLAG

    if (FILE_FLAG .eq. 1) then
        CALL MPI_TYPE_FREE(filetype,ierr)
        FILE_FLAG = 0
    end if
    CALL MPI_TYPE_CREATE_SUBARRAY(nd, gsizes, lsizes, global_indices, &
                                    MPI_ORDER_FORTRAN, datatype,  filetype, ierr)
    CALL MPI_TYPE_COMMIT(filetype, ierr)
    CALL MPI_FILE_SET_VIEW(fh, offset, datatype, filetype, &
                            'native', MPI_INFO_NULL, ierr)
    FILE_FLAG = 1
    
end subroutine Create_commit_fileview


subroutine Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
    implicit none

    integer,intent(inout)           :: memtype
    integer,intent(in)              :: datatype
    integer, dimension(3),intent(in):: memsizes, lsizes, local_indices
    integer,intent(inout)           :: MEM_FLAG


    if (MEM_FLAG .eq. 1) then
        CALL MPI_TYPE_FREE(memtype,ierr)
        MEM_FLAG = 0
    end if
    CALL MPI_TYPE_CREATE_SUBARRAY(nd, memsizes, lsizes, local_indices, &
            MPI_ORDER_FORTRAN, datatype,  memtype, ierr)
    CALL MPI_TYPE_COMMIT(memtype, ierr)
    MEM_FLAG = 1
    
end subroutine Create_commit_subarray


! Get the current number for output file 

function get_iter()
    use computational_constants_MD, only : restart_numbering

    integer :: get_iter

    if (restart_numbering) then
        get_iter = iter-initialstep+1
    else
        get_iter = iter
    endif

end function get_iter

end module module_parallel_io

!======================================================================
!                   INPUTS                        =
!======================================================================


!=============================================================================
! setup_command_arguments
! Checks for command-line arguments passed to the program and assigns the 
! relevant values.
!-----------------------------------------------------------------------------
subroutine setup_command_arguments()
use module_parallel_io
implicit none
    
    integer :: i,argcount
    logical :: restart_file_exists, input_file_exists
    character(len=200) :: arg,nextarg
    character(len=270) :: errorstr

    !Set default values in case they aren't specified by user
    restart = .false.           
    input_file_exists = .false.
    restart_file_exists = .false.
    input_file = trim(prefix_dir)//'MD.in'
    initial_microstate_file = trim(prefix_dir)//'final_state'
    cyl_file = trim(prefix_dir)//'cylinders'

    argcount = command_argument_count()

    if (argcount.gt.0) then         !If more than 0 arguments   
                
        do i=1,argcount         !Loop through all arguments...
                
            call get_command_argument(i,arg)    !Reading two at once
            if ( i < argcount) then 
                call get_command_argument(i+1,nextarg)
            else
                nextarg=''
            endif
                
            if (trim(arg).eq.'-r' .and. nextarg(1:1).ne.'-') then
                initial_microstate_file = trim(nextarg)
                inquire(file=initial_microstate_file, exist=restart) !Check file exists
                if (restart .eqv. .false.) then
                    write(errorstr,'(3a)') "Restart file not found, please check file ", & 
                                        trim(initial_microstate_file)," exists or remove -r flag"
                    call error_abort(errorstr)
                endif
            end if

            if (trim(arg).eq.'-i' .and. nextarg(1:1).ne.'-') then
                input_file = trim(nextarg)
            end if

            if (trim(arg).eq.'-c' .and. nextarg(1:1).ne.'-') then
                cyl_file = trim(nextarg)
            end if

        end do

    end if

    inquire(file=input_file, exist=input_file_exists)   !Check file exists
    if (input_file.eq.'./MD.in'.and..not. input_file_exists) input_file = './default.in'
    inquire(file=input_file, exist=input_file_exists)
    
    if(.not. input_file_exists) then
        print*, 'Input file ', trim(input_file), ' not found. Stopping simulation.'
        call error_abort
    end if

end subroutine setup_command_arguments

!-----------------------------------------------------------------------------
! Subroutine:   setup_inputs
! Author(s):    David Trevelyan & Edward Smith
! Description:
!       The input file MD.in contains capitalised keywords followed by
!       numerical values. The "locate" subroutine rewinds to the beginning
!       of MD.in and scans each line until the keyword is matched.
!
!       Consequently, the file position is set for the next statement to read
!       the line underneath the previously "located" keyword. 
!-----------------------------------------------------------------------------

subroutine setup_inputs()
    use module_parallel_io
    use librarymod, only : locate
    implicit none
    
    integer             :: n, tvalue(8)

    call random_seed(size=n)
    allocate(seed(n))

    !Read input file
    call setup_read_input

    rcutoff2= rcutoff**2        !Useful definition to save computational time
    initialstep = 0             !Set initial step to one to start
    
    if (seed(1)==seed(2)) then
        ! Randomisations 
        call random_seed(get=seed(1:n))
        call date_and_time(values=tvalue)
        seed=IEOR(tvalue(8)+irank,seed)
    else 
        !Assign different random number seed to each processor
        seed =  irank
    endif

    !Assign seed to random number generator
    call random_seed(put=seed(1:n))

    elapsedtime = delta_t*Nsteps !Set elapsed time to end of simualtion

end subroutine setup_inputs

!-------------------------------------------------------------------------------------
!                    Restart Simulation inputs
! Set up inputs on every processor, based on the final state of a previous simulation
!-------------------------------------------------------------------------------------

subroutine setup_restart_inputs()
    use module_parallel_io
    use librarymod, only : locate
    implicit none

    integer                         :: n, ixyz
    integer                         :: prev_nproc
    integer                         :: extrasteps
    integer                         :: checkint
    integer(kind=MPI_OFFSET_KIND)   :: ofs, header_ofs
    integer(selected_int_kind(18))  :: header_pos
    real(kind(0.d0))                :: checkdp
    character(400)                  :: error_message

    !Allocate random number seed
    call random_seed(size=n)
    allocate(seed(n))

    !Read input file
    call setup_read_input
    extrasteps = Nsteps

    !=====================================================================================================!
    !========================   R E A D    R E S T A R T    H E A D E R   ================================!
    !Check if values from input file are different and alert user - all processors have
    !read the same file so only need to check on one processor

    !Open on a single process and broadcast if different
    if (irank .eq. iroot) then
            
        call MPI_File_open(MPI_COMM_SELF, initial_microstate_file, & 
        MPI_MODE_RDONLY , MPI_INFO_NULL, restartfileid, ierr)

        ! read the size of offset, just in case we hit a system that still uses 32 bits addresses
        ! read in a 8 bit integer
        ofs = -8
        call mpi_file_seek(restartfileid,ofs,mpi_seek_end,ierr)
        call mpi_file_read(restartfileid,header_pos    ,1,mpi_integer8,MPI_STATUS_IGNORE,ierr)
        header_ofs = header_pos
        call MPI_FILE_SEEK(restartfileid,header_ofs,MPI_SEEK_SET,ierr)

        call MPI_File_read(restartfileid,globalnp        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        if (checkint .ne. initialnunits(1)) then
            print*, 'Discrepancy between x initialnunits in input and restart file:', &
                    'restart file will be used. Resetting initialnunits(1) from', initialnunits(1), &
                    'to', checkint 
            initialnunits(1) = checkint
        endif
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        if (checkint .ne. initialnunits(2)) then
            print*, 'Discrepancy between y initialnunits in input and restart file:', &
                    'restart file will be used. Resetting initialnunits(2) from', initialnunits(2), &
                    'to', checkint 
            initialnunits(2) = checkint
        endif
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        if (checkint .ne. initialnunits(3)) then
            print*, 'Discrepancy between z initialnunits in input and restart file:', &
                    'restart file will be used. Resetting initialnunits(3) from', initialnunits(3), &
                    'to', checkint 
            initialnunits(3) = checkint
        endif
        call MPI_File_read(restartfileid,Nsteps          ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !tplot
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !seed
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !seed
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !periodic
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !periodic
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !periodic
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !potential_flag
        if (checkint .ne. potential_flag) then
            print*, 'Discrepancy between potential_flag', &
                      'in input & restart file - restart file will be used'
            potential_flag = checkint
        endif
        call MPI_File_read(restartfileid,prev_rtrue_flag ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !rtrue_flag
        if (prev_rtrue_flag .ne. rtrue_flag) then
            print*, 'Discrepancy between rtrue_flag', &
                    'in input & restart file - current file will be used,', &
                    ' but rtrue will still be read from restart file.'  
            print*, 'prev_rtrue_flag:', prev_rtrue_flag
            print*, 'rtrue_flag:', rtrue_flag
        endif
        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !solvent_flag
        if (checkint .ne. solvent_flag) then
            print*, 'Discrepancy between solvent_flag', &
                      'in input & restart file - restart file will be used'
            solvent_flag = checkint
        endif

        call MPI_File_read(restartfileid,checkint       ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !nmonomers
        if (checkint.ne.nmonomers) then
            print*, 'Discrepancy between nmonomers', &
                  'in input & restart file - restart file will be used'
            nmonomers = checkint
        endif

        procnp = 0; proc_reorder = 0;   prev_nproc = 1
        !Small debugging run (nproc<27) - if proc mismatch use serial reordering (all read everything and discard)
        if(npx .le. 3 .and. npy .le. 3 .and. npz .le. 3) then
            error_message = 'Small debug run (less than 3 x 3 x 3 processors). &
                            &Molecules will be assigned to correct processors - all read everything and discard'
            call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            if (checkint .ne. npx)  proc_reorder = 1; prev_nproc = prev_nproc*checkint
            call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            if (checkint .ne. npy)  proc_reorder = 1; prev_nproc = prev_nproc*checkint
            call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            if (checkint .ne. npz)  proc_reorder = 1; prev_nproc = prev_nproc*checkint
            if (proc_reorder.eq.0) then
                do n=1,prev_nproc           !Loop through all processors and store for restart
                    call MPI_File_read(restartfileid,procnp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
                enddo
                do n=1,prev_nproc           !Loop through all processors and store for restart
                    call MPI_File_read(restartfileid,proctethernp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
                enddo
            else
                do n=1,prev_nproc           !Loop through all processors and discard
                    call MPI_File_read(restartfileid,checkint,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
                enddo
                do n=1,prev_nproc           !Loop through all processors and discard
                    call MPI_File_read(restartfileid,checkint,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
                enddo
            endif
        !Large run - if proc mismatch, stop code and suggest reordering in serial
        else
            error_message = 'Number of processors in input file does not match the restart file.            &
                            &Options:                                                                        & 
                            &1) Set processors in input file to zeros to use restart proc topology.         & 
                            &2) Reorder restart file in serial for current proc topology                    & 
                            &3) Run in serial or with less than 3x3x3 processors                        '
            call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            if (checkint .ne. npx) call error_abort(error_message)
            call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            if (checkint .ne. npy) call error_abort(error_message)
            call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            if (checkint .ne. npz) call error_abort(error_message)
            prev_nproc = npx*npy*npz    !If parallel run, number of molecules per processor written
            do n=1,prev_nproc           !Loop through all processors and store for restart
                call MPI_File_read(restartfileid,procnp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            enddo
            do n=1,prev_nproc           !Loop through all processors and store for restart
                call MPI_File_read(restartfileid,proctethernp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            enddo
        endif


        ! - - - Check here if record is globaldomain or density- - -
        !Included for backwards compatibility -- before revision 673 finalstate did not include globaldomain
        !The assumption here is that globaldomain is always greater than 3.0 while density is always less than 3
        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        if (checkdp .lt. 3.d0) then
            print*,  'Warning -- Since revision 673, globaldomain included in finalstate ', &
                     ' this appears to be an older restart file.', & 
                      '(or globaldomain in x = ',checkdp,'). ', &
                      'It is assumed that global domain is not specified and density = ', checkdp

            !Read density
            if (checkdp .ne. density) then
                print*, 'Discrepancy between system density', &
                        'in input & restart file - restart file will be used'
                density = checkdp
            endif
            globaldomain = (initialnunits(:)/((density/4.d0)**(1.d0/nd)))
        else
        
            !Read 2 other global domain values followed by density and check they match
            do ixyz = 1,nd
                if (checkdp .ne. globaldomain(ixyz)) then
                    print*, 'Discrepancy between globaldomain(', ixyz, ')', globaldomain(ixyz), checkdp, &
                            'in input & restart file - restart file will be used'
                    globaldomain(ixyz) = checkdp
                endif
                call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            end do
            if (checkdp .ne. density) then
                print*, 'Discrepancy between system density', &
                        'in input & restart file - restart file will be used'
                density = checkdp
            endif

            ! Check units match globaldomain
            do ixyz=1,nd
                if ( globaldomain(ixyz) .ne.  (initialnunits(ixyz)/((density/4.d0)**(1.d0/nd))) ) then
                    print('(a,i1,a,f15.5,a,f15.5)'),'Warning: globaldomain(',ixyz,') = ',globaldomain(ixyz),&
                    ' does not match initialunits and density calculation value: ', &
                    initialnunits(ixyz)/((density/4.d0)**(1.d0/nd))
                end if
            enddo

        endif

        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        if (checkdp .ne. rcutoff) then
            print*, 'Discrepancy between cut off radius', &
                    'in input & restart file - restart file will be used'
            rcutoff = checkdp
        endif
        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)    !delta_t 
        call MPI_File_read(restartfileid,elapsedtime     ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_read(restartfileid,simtime         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)    !k_c 
        if (checkdp.ne.k_c) then
            print*, 'Discrepancy between k_c', &
                'in input & restart file - restart file will be used'
            k_c = checkdp
        endif
        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)    !R_0
        if (checkdp.ne.R_0) then
            print*, 'Discrepancy between R_0', &
                      'in input & restart file - restart file will be used'
            R_0 = checkdp
        endif
        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)    !eps_pp
        if (checkdp.ne.eps_pp) then
            print*, 'Discrepancy between eps_pp', &
                      'in input & restart file - restart file will be used'
            eps_pp = checkdp
        endif
        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)    !eps_ps
        if (checkdp.ne.eps_ps) then
            print*, 'Discrepancy between eps_ps', &
                      'in input & restart file - restart file will be used'
            eps_ps = checkdp
        endif
        call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)    !eps_ss
        if (checkdp.ne.eps_ss) then
            print*, 'Discrepancy between eps_ss', &
                      'in input & restart file - restart file will be used'
            eps_ss = checkdp
        endif
        call MPI_File_read(restartfileid,delta_rneighbr  ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)    !delta_rneighbr

        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !Mie_potential flag
        if (checkint .eq. 0 .and. Mie_potential .eq. 1) then
            print*, 'Discrepancy between potential_flag', &
                    'in input & restart file - Mie flags ', & 
                    'will attempt to be assigned to restart config'
            Mie_potential = 2
        endif

        call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !global_numbering flag
        if (checkint .eq. 1 .and. global_numbering .eq. 0) then
            print*, 'global_numbering used in restart file but not requested in', &
                    'input file - global_numbering will be used '
            global_numbering = checkint
        elseif (checkint .eq. 0 .and. global_numbering .eq. 1) then
            print*, 'global_numbering not in restart file but requested in', &
                    'input file - global_numbering info will be created '
            global_numbering = 2
        endif
        call MPI_File_close(restartfileid,ierr)
    endif

    !---------------Broadcast data read by root to all other processors-------------------------!
    ! temporary np, exact value will be fixed after reading r and v arrays
    ! np is needed in set_parameters_outputs that is called before microstate initialisation
    ! when restart is true
    call MPI_BCAST(globalnp,          1,MPI_integer,iroot-1,MD_COMM,ierr)
    np = globalnp/nproc
    call MPI_BCAST(rcutoff,           1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
    rcutoff2 = rcutoff**2             !Useful definition to save computational time
    call MPI_BCAST(initialnunits,     3,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(Nsteps,            1,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(potential_flag,    1,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(prev_rtrue_flag,   1,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(solvent_flag,      1,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(proc_reorder,      1,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(procnp, size(procnp),MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(proctethernp, size(proctethernp),MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(nmonomers,         1,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(globaldomain,      3,MPI_double_precision,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(density,           1,MPI_double_precision,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(inputtemperature,  1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
    call MPI_BCAST(elapsedtime,       1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
    call MPI_BCAST(simtime,           1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
    call MPI_BCAST(k_c,               1,MPI_double_precision,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(R_0,               1,MPI_double_precision,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(eps_pp,            1,MPI_double_precision,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(eps_ps,            1,MPI_double_precision,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(eps_ss,            1,MPI_double_precision,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(Mie_potential,     1,MPI_integer,iroot-1,MD_COMM,ierr)
    call MPI_BCAST(global_numbering,  1,MPI_integer,iroot-1,MD_COMM,ierr)

    elapsedtime = elapsedtime + delta_t*extrasteps  !Set elapsed time to end of simualtion
    initialstep = Nsteps                            !Set plot count to final plot of last
    iter = initialstep                              !Set iter to initialstep
    Nsteps = Nsteps + extrasteps                    !Establish final iteration step based on previous

end subroutine setup_restart_inputs

!------------------------------------------------------------------------------
subroutine setup_restart_microstate()
    use module_parallel_io
    use messenger_data_exchange, only : globalSum
    implicit none

    logical                                      :: tag_off=.false.
    integer                                      :: i,n,nl,procassign
    integer                                      :: pos
    integer                                      :: dp_datasize
    integer(kind=MPI_OFFSET_KIND)                :: disp, procdisp
    integer, dimension(:), allocatable           :: bufsize
    real(kind(0.d0))                             :: tagtemp, moltypetemp, globnotemp
    real(kind(0.d0)), dimension (nd)             :: rtemp,vtemp,rtruetemp,rtethertemp
    real(kind(0.d0)), dimension (nsdmi)          :: monomertemp
    real(kind(0.d0)), dimension (:), allocatable :: buf !Temporary variable

    allocate(bufsize(nproc))
    bufsize = 0

    !Allocate temporary tag array
    if (allocated(tag) .neqv. .true.) then
        allocate(tag(np+extralloc)); tag = free
        tag_off = .true.
    endif
    
    !Determine size of datatypes
    call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

    !Open restart file on all processor
    call MPI_FILE_OPEN(MD_COMM, initial_microstate_file, & 
                       MPI_MODE_RDONLY , MPI_INFO_NULL, restartfileid, ierr)

    select case(proc_reorder)
    case(0)

        bufsize(irank) = 2*nd*procnp(irank) + 1*procnp(irank)
        if (prev_rtrue_flag.eq.1) then
            bufsize(irank) = bufsize(irank) + nd*procnp(irank)
        end if
        bufsize(irank) = bufsize(irank) + nd*proctethernp(irank)
        if (Mie_potential .eq. 1) then
            bufsize(irank) = bufsize(irank) + procnp(irank)
        end if
        !Add global number if required
        if (global_numbering .eq. 1) then
            bufsize(irank) = bufsize(irank) + procnp(irank)
        endif
        if (potential_flag .eq. 1) then
            bufsize(irank) = bufsize(irank) + nsdmi*procnp(irank)
        end if
        call globalSum(bufsize,nproc)

        !Obtain displacement of each processor using procs' np from restart file
        !with 6 position (3 wrapped, 3 unwrapped) and 3 velocity components
        !for each molecule
        procdisp = 0
        do i=1,irank-1
            procdisp = procdisp + bufsize(i)*dp_datasize
        end do
        !Obtain location to write in file
        disp =  procdisp

        allocate(buf(bufsize(irank)))
        !Set each processor to that location and read particlewise
        call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
                               MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_READ_ALL(restartfileid, buf, bufsize(irank), MPI_double_precision, & 
                               MPI_STATUS_IGNORE, ierr) !Read position from file
    
        pos = 1
        do nl=1,procnp(irank)

            tag(nl) = nint(buf(pos))
            pos = pos + 1
            !Correct to local coordinates
            r(1,nl) = buf(pos)  -domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
            r(2,nl) = buf(pos+1)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
            r(3,nl) = buf(pos+2)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)
            pos = pos + 3
            !Read velocities
            v(1,nl) = buf(pos)
            v(2,nl) = buf(pos+1)
            v(3,nl) = buf(pos+2)
            pos = pos + 3
            if (prev_rtrue_flag.eq.1) then

                if (rtrue_flag .eq. 1) then
                    !Store true positions
                    rtrue(1,nl) = buf(pos)
                    rtrue(2,nl) = buf(pos+1)
                    rtrue(3,nl) = buf(pos+2)
                end if
                pos = pos + 3

            end if
            if (any(tag(nl).eq.tether_tags)) then
                !Read tether position, corrected to local coords
                rtether(1,nl) = buf(pos)  -domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
                rtether(2,nl) = buf(pos+1)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
                rtether(3,nl) = buf(pos+2)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)
                pos = pos + 3
            end if
            if (Mie_potential .eq. 1) then
                moltype(nl) = nint(buf(pos))
                pos = pos + 1
            endif
            !Load global number if required
            if (global_numbering .eq. 1) then
                glob_no(nl) = nint(buf(pos))
                pos = pos + 1
            !Create global number if not present in restart file
            elseif (global_numbering .eq. 2) then
                glob_no(nl) = nl + sum(procnp(1:irank-1))
                write(2530,*) nl, glob_no(nl)
            endif
            if (potential_flag.eq.1) then
                !Read monomer data
                monomer(nl)%chainID        = nint(buf(pos))
                monomer(nl)%subchainID     = nint(buf(pos+1))
                monomer(nl)%funcy          = nint(buf(pos+2))
                monomer(nl)%glob_no        = nint(buf(pos+3))
                monomer(nl)%bin_bflag(1:4) = nint(buf(pos+4:pos+7))
                pos = pos + 8
            end if

        enddo

        np = procnp(irank)
    
        deallocate(buf)

    case(1) !Reorder flag triggered

        nl = 0      !Reset local molecules count nl

        !---------- For all molecule positions ------------
        !Move through location of position co-ordinates
        do n=1,globalnp

            !---------------  READ ONE MOL -------------------------!
            !Read tag
            call MPI_FILE_READ_ALL(restartfileid, tagtemp, 1, MPI_DOUBLE_PRECISION, &
                                   MPI_STATUS_IGNORE, ierr)
            !Read position  
            call MPI_FILE_READ_ALL(restartfileid, rtemp, 3, MPI_DOUBLE_PRECISION, &
                                   MPI_STATUS_IGNORE, ierr)
            !Read velocity
            call MPI_FILE_READ_ALL(restartfileid, vtemp, 3, MPI_DOUBLE_PRECISION, &
                                   MPI_STATUS_IGNORE, ierr)
            if (prev_rtrue_flag .eq. 1) then
                call MPI_FILE_READ_ALL(restartfileid, rtruetemp, 3, MPI_DOUBLE_PRECISION, &
                                       MPI_STATUS_IGNORE, ierr)
            end if
            if (any(nint(tagtemp).eq.tether_tags)) then
                call MPI_FILE_READ_ALL(restartfileid, rtethertemp, 3, MPI_DOUBLE_PRECISION, &
                                       MPI_STATUS_IGNORE, ierr)
            end if
            if (Mie_potential .eq. 1) then
                call MPI_FILE_READ_ALL(restartfileid, moltypetemp, 1, MPI_DOUBLE_PRECISION, &
                                       MPI_STATUS_IGNORE, ierr)
            endif
            if (global_numbering .eq. 1) then
                call MPI_FILE_READ_ALL(restartfileid, globnotemp, 1, MPI_DOUBLE_PRECISION, &
                                       MPI_STATUS_IGNORE, ierr)
            endif
            if (potential_flag.eq.1) then
                call MPI_FILE_READ_ALL(restartfileid, monomertemp, nsdmi, MPI_DOUBLE_PRECISION, &
                                       MPI_STATUS_IGNORE, ierr)
            end if
            !------------ END READ ONE MOL -------------------------!

            !Use integer division to determine which processor to assign molecule to
            procassign = ceiling((rtemp(1)+globaldomain(1)/2.d0)/domain(1))
            if (procassign .ne. iblock) cycle
            procassign = ceiling((rtemp(2)+globaldomain(2)/2.d0)/domain(2))
            if (procassign .ne. jblock) cycle
            procassign = ceiling((rtemp(3)+globaldomain(3)/2.d0)/domain(3))
            if (procassign .ne. kblock) cycle

            !If molecules is in the domain then add to processor's total
            nl = nl + 1 !Local molecule count

            tag(nl) = nint(tagtemp)
            !Correct to local coordinates
            r(1,nl) = rtemp(1)-domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
            r(2,nl) = rtemp(2)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
            r(3,nl) = rtemp(3)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)

            v(:,nl) = vtemp(:)

            if (prev_rtrue_flag .eq. 1) then
                !Read true unwrapped positions
                rtrue(1,nl) = rtruetemp(1)
                rtrue(2,nl) = rtruetemp(2)
                rtrue(3,nl) = rtruetemp(3)
            end if  

            if (any(tag(nl).eq.tether_tags)) then
                !Read tethered positions
                rtether(1,nl) = rtethertemp(1)-domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
                rtether(2,nl) = rtethertemp(2)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
                rtether(3,nl) = rtethertemp(3)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)
            end if
            if (Mie_potential .eq. 1) then
                moltype(nl) = moltypetemp
            endif
            !Add global number if required
            if (global_numbering .eq. 1) then
                glob_no(nl) = globnotemp
            !Create global number if not present in restart file
            elseif (global_numbering .eq. 2) then
                glob_no(nl) = nl + sum(procnp(1:irank-1))
            endif
            if (potential_flag.eq.1) then
                monomer(nl)%chainID        = nint(monomertemp(1))
                monomer(nl)%subchainID     = nint(monomertemp(2))
                monomer(nl)%funcy          = nint(monomertemp(3))
                monomer(nl)%glob_no        = nint(monomertemp(4))
                monomer(nl)%bin_bflag(1:4) = nint(monomertemp(5:8))
            end if
            if (mod(n,1000) .eq. 0) print'(a,f10.2)', & 
                ' Redistributing molecules to input processor topology - % complete =', (100.d0*n/globalnp)
        enddo
    
        np = nl !Correct local number of particles on processor

    case default

        call error_abort('processor re-ordering flag incorrect in restart microstate')

    end select

    if (global_numbering .eq. 2) then
        global_numbering = 1
    endif

    ! Mie moltype should be from restart file!! If mie_potential was zero in restart
    ! but now input requests one, setup as if new run (based on location, etc).
    if (Mie_potential .eq. 2) then
        call setup_moltypes_wall                  !Setup type of molecules
        mie_potential = 1
    endif

    ! Determine number of chains by global maximum of chainID
    if (potential_flag .eq. 1) then
        nchains = maxval(monomer(:)%chainID)
        call globalMaxInt(nchains)
    end if

    !Close file used to load initial state and remove if called "final_state" 
    !to prevent confusion with final state of current run
    call MPI_FILE_CLOSE(restartfileid, ierr)
    if (initial_microstate_file .eq. './results/final_state') then
        call MPI_FILE_DELETE(initial_microstate_file, MPI_INFO_NULL, ierr)
    endif

    if (reset_tags_on_restart .eq. 1) then
        if (irank.eq.iroot) print*, 'Molecular tags reset based on input file.'
        call setup_location_tags(0)              !Setup locn of fixed mols
    elseif (reset_tags_on_restart .eq. 2) then
        if (irank.eq.iroot) print*, 'Molecular tags reset based on input file.'
        call setup_location_tags(1)              !Setup locn of thermostats only
    else
        if (irank.eq.iroot) print*, 'Molecular tags have been read from restart file.'
    endif
    call get_tag_thermostat_activity(tag_thermostat_active)
    if (irank.eq.iroot) print*, 'Thermostat thermostat activity = ', tag_thermostat_active
    do n = 1,np
        call read_tag(n)        !Read tag and assign properties
    enddo
    if (irank.eq.iroot) print*, 'Molecular tags have been imported and processed.'

    deallocate(bufsize)
    if (tag_off) deallocate(tag)    !Tags off so tag info not necessary
    
    !Choose initial molecular velocities using velocity flag
    select case(initial_velocity_flag)
    case(0)
        !Do nothing and use restart file velocities
    case(1)
        select case (trim(velocity_special_case))
        case('debug')
            call setup_initialise_velocities_test
        case('taylor_green')
            call setup_initialise_velocities_TG_parallel
        case('dns')
            call set_velocity_field_from_DNS_restart(trim(DNS_filename),DNS_ngx,DNS_ngy,DNS_ngz)
        case default
            call error_abort('Unidentified initial velocities_special_case')    
        end select
    case default
        call error_abort('Unidentified initial velocity flag')  
    end select

end subroutine setup_restart_microstate


#if __INTEL_COMPILER > 1200
    subroutine load_bforce_pdf()
        use boundary_MD
        use librarymod, only: get_new_fileunit
        implicit none

        integer :: nperbin, f
        real(kind(0.d0)) :: histbinsize

        f = get_new_fileunit()
        open(unit=f, file='bforce.input',action='read',form='unformatted',&
             access='stream')
        read(f) bforce_pdf_nsubcells
        read(f) bforce_pdf_nbins
        read(f) bforce_pdf_min
        read(f) bforce_pdf_max
        read(f) nperbin 
        allocate(bforce_pdf_input_data(bforce_pdf_nsubcells, bforce_pdf_nbins, nperbin))
        read(f) bforce_pdf_input_data 
        close(f,status='keep')

        bforce_pdf_binsize = (bforce_pdf_max - bforce_pdf_min)/real(bforce_pdf_nbins)

    end subroutine load_bforce_pdf
#endif

!=============================================================================
! Import cylinders to be filled 
subroutine parallel_io_cyl_footer(infile)
    use mpi
    use computational_constants_MD, only: globaldomain
    use physical_constants_MD, only: globalnp
    use concentric_cylinders, only: cyl_file, r_oo, r_io, r_oi, r_ii, cyl_np
    implicit none

    character(*), intent(in) :: infile

    integer :: ierr, fileid
    integer(kind=MPI_OFFSET_KIND) :: ofs, header_ofs
    integer(selected_int_kind(18)) :: header_pos

    call MPI_FILE_OPEN(MPI_COMM_SELF, infile, MPI_MODE_RDONLY, &
                       MPI_INFO_NULL, fileid, ierr)

    if (ierr .ne. 0) stop 'Error opening cylinder file'

    ! Read size of offset, in case we hit a system that uses 32 bit
    ! addresses
    ofs = -8 ! Read an 8-bit integer
    call MPI_FILE_SEEK(fileid, ofs, MPI_SEEK_END, ierr)
    call MPI_FILE_READ(fileid, header_pos, 1, MPI_INTEGER8, &
                       MPI_STATUS_IGNORE,ierr)
    header_ofs = header_pos
    call MPI_FILE_SEEK(fileid,header_ofs,MPI_SEEK_SET,ierr)

    call MPI_FILE_READ(fileid,cyl_np,1,MPI_INTEGER, &
                       MPI_STATUS_IGNORE,ierr)

    call MPI_FILE_READ(fileid,globaldomain,3,MPI_DOUBLE_PRECISION, &
                       MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_READ(fileid,r_oo,1,MPI_DOUBLE_PRECISION, &
                       MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_READ(fileid,r_io,1,MPI_DOUBLE_PRECISION, &
                       MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_READ(fileid,r_oi,1,MPI_DOUBLE_PRECISION, &
                       MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_READ(fileid,r_ii,1,MPI_DOUBLE_PRECISION, &
                       MPI_STATUS_IGNORE,ierr)

    call MPI_FILE_CLOSE(fileid,ierr)
    !print*, 'Read cylinder footer information:'
    !print*, 'globalnp:', globalnp
    !print*, 'globaldomain:', globaldomain
    !print*, 'r_oo:', r_oo
    !print*, 'r_io:', r_io 
    !print*, 'r_oi:', r_oi 
    !print*, 'r_ii:', r_ii 

end subroutine parallel_io_cyl_footer

subroutine parallel_io_import_cylinders()
    use mpi
    use concentric_cylinders, only: cyl_file, cyl_np
    use physical_constants_MD, only: globalnp, np
    use computational_constants_MD, only: iblock, jblock, kblock, irank,&
                                          globaldomain, domain, teth
    use arrays_MD, only: r, tag
    use messenger, only: localise, MD_COMM
    implicit none

    integer :: fileid, ierr
    integer :: nl, n, procassign
    real    :: rtemp(3)
   

    !Open restart file on all processor
    call MPI_FILE_OPEN(MPI_COMM_SELF, cyl_file, & 
                       MPI_MODE_RDONLY , MPI_INFO_NULL, fileid, ierr)

    nl = 0      !Reset local molecules count nl
    !---------- For all molecule positions ------------
    !Move through location of position co-ordinates
    do n=1,cyl_np

        !Read global position   
        call MPI_FILE_READ(fileid, rtemp, 3, MPI_REAL, &
                           MPI_STATUS_IGNORE, ierr)

        !Use integer division to determine which processor to assign molecule to
        procassign = ceiling((rtemp(1)+globaldomain(1)/2.d0)/domain(1))
        if (procassign .ne. iblock) cycle
        procassign = ceiling((rtemp(2)+globaldomain(2)/2.d0)/domain(2))
        if (procassign .ne. jblock) cycle
        procassign = ceiling((rtemp(3)+globaldomain(3)/2.d0)/domain(3))
        if (procassign .ne. kblock) cycle

        !If molecules is in the domain then add to processor's total
        nl = nl + 1 !Local molecule count

        !Correct to local coordinates
        r(:,nl) = localise(real(rtemp(:),kind(0.d0)))

        if (mod(n,1000) .eq. 0) print'(a,f10.2)', & 
            ' Redistributing cylinder molecules to input processor topology - % complete =', (100.d0*dble(n)/dble(cyl_np))

    enddo

    np = nl

    call MPI_FILE_CLOSE(fileid, ierr)

end subroutine parallel_io_import_cylinders

!======================================================================
!=                                      OUTPUTS                       =
!======================================================================

!------------------------------------------------------------------------
!Write positions and velocities of molecules to a file for restart

subroutine parallel_io_final_state(interim_output)
    use module_parallel_io
    use polymer_info_MD
    use messenger_data_exchange, only : globalSum
    use librarymod, only : get_Timestep_FileName
    implicit none

    logical,intent(in) :: interim_output

    integer                                         :: n, i
    integer                                         :: pos
    integer                                         :: dp_datasize
    integer, dimension(:), allocatable              :: bufsize
    integer(kind=MPI_OFFSET_KIND)                   :: disp, procdisp, filesize
    integer(kind=selected_int_kind(18))             :: header_pos
    real(kind(0.d0)), dimension(:,:), allocatable   :: rglobal,rtetherglobal
    real(kind(0.d0)), dimension(:)  , allocatable   :: buf
    character(32)                                   :: filename_incno

    allocate(bufsize(nproc))
    bufsize = 0
    
    !Rebuild simulation before recording final state
    call linklist_deallocateall         !Deallocate all linklist components
    call sendmols                       !Exchange particles between processors
    call sort_mols                      !Improved spatial locality helps restart on different proc topology
    call assign_to_cell                 !Re-build linklist every timestep
    call messenger_updateborders(1)     !Update borders between processors
    call assign_to_neighbourlist        !Setup neighbourlist

    !Build array of number of particles on neighbouring
    !process' subdomains on current proccess
    call globalGathernp

    !Determine size of datatypes
    call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

    !Adjust r according to actual location for storage according
    !to processor topology with r = 0 at centre
    allocate(rglobal(3,np))
    rglobal(1,:) = r(1,1:np)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
    rglobal(2,:) = r(2,1:np)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
    rglobal(3,:) = r(3,1:np)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

    !Allocate tag array to write for restart so as to maximise compatibility
    if (ensemble.ne.tag_move) then
        !allocate(tag(np));     tag(:) = free 
    else
        !Convert any tethered molecules to global coordinates ready to write out
        allocate(rtetherglobal(3,np))
        rtetherglobal(1,:) = rtether(1,1:np)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
        rtetherglobal(2,:) = rtether(2,1:np)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
        rtetherglobal(3,:) = rtether(3,1:np)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
    endif

    if (.not. interim_output) then
        !Remove previous final state file
        call MPI_FILE_DELETE(trim(prefix_dir)//'results/final_state', MPI_INFO_NULL, ierr)

        !Open file on all processors
        call MPI_FILE_OPEN( MD_COMM,trim(prefix_dir)//'results/final_state', & 
                            MPI_MODE_RDWR + MPI_MODE_CREATE, & 
                            MPI_INFO_NULL, restartfileid, ierr)
    else
        call get_Timestep_FileName(iter, trim(prefix_dir)//'results/interim_state', filename_incno)

        !Open file on all processors
        call MPI_FILE_OPEN( MD_COMM, filename_incno, & 
                            MPI_MODE_RDWR + MPI_MODE_CREATE, & 
                            MPI_INFO_NULL, restartfileid, ierr)
    endif

    ! Calculate buffer size ---------------------------------------!
    ! Attention: np is changed inside reorder_data%sendmols call
    ! Add 2*nd for r and v, 1 for tag
    bufsize(irank) = 2*procnp(irank)*nd + 1*procnp(irank)
    ! If rtrue on, add nd
    if (rtrue_flag .eq. 1) then
        bufsize(irank) = bufsize(irank) + procnp(irank)*nd
    end if
    ! For any tethered, add nd for rtether and count proctethernp
    proctethernp = 0
    do n=1,np
        if (any(tag(n).eq.tether_tags)) then
            proctethernp(irank) = proctethernp(irank) + 1
            bufsize(irank) = bufsize(irank) + nd
        end if
    end do
    ! If mie potential, add space for moltypes
    if (Mie_potential .eq. 1) then
        bufsize(irank) = bufsize(irank) + procnp(irank)
    end if
    !If global molecular numbers, add space
    if (global_numbering .eq. 1) then
        bufsize(irank) = bufsize(irank) + procnp(irank)
    endif
    ! If polymer sim, add space for polymer info
    if (potential_flag .eq. 1) then
        bufsize(irank) = bufsize(irank) + nsdmi*procnp(irank)
    end if
    !--------------------------------------------------------------!

    !Collect all bufsizes and proctethernps
    call globalSum(bufsize,nproc)
    call globalSum(proctethernp,nproc)

    ! Allocate buffer to be written to final_state file
    allocate(buf(bufsize(irank)))
    !Populate buffer -----------------------------------------------------!
    pos = 1
    do n = 1,np
        buf(pos) = real(tag(n),kind(0.d0)); pos = pos + 1
        buf(pos:pos+2) = rglobal(:,n);      pos = pos + 3
        buf(pos:pos+2) = v(:,n);            pos = pos + 3
        if (rtrue_flag .eq. 1) then
            buf(pos:pos+2) = rtrue(:,n);    pos = pos + 3
        end if
        if (any(tag(n) .eq. tether_tags)) then
            buf(pos:pos+2) = rtetherglobal(:,n);  pos = pos + 3
        end if
        if (mie_potential .eq. 1) then
            buf(pos) = moltype(n);  pos = pos + 1
        end if
        !If global molecular numbers, add space
        if (global_numbering .eq. 1) then
            buf(pos) = glob_no(n);  pos = pos + 1
        endif
        if (potential_flag .eq. 1) then
            buf(pos)     = real(monomer(n)%chainID,kind(0.d0))
            buf(pos+1)   = real(monomer(n)%subchainID,kind(0.d0))
            buf(pos+2)   = real(monomer(n)%funcy,kind(0.d0))
            buf(pos+3)   = real(monomer(n)%glob_no,kind(0.d0))
            buf(pos+4:pos+7) = real(monomer(n)%bin_bflag(1:4),kind(0.d0))
            pos = pos + 8
        end if
    end do
    !---------------------------------------------------------------------!

    ! Obtain file displacement for each process     
    procdisp = 0
    do i=1,irank-1
        procdisp = procdisp + bufsize(i)*dp_datasize
    enddo
    disp = procdisp
    !Set each processor to that location
    call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
                           MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
    ! Write buffer to own space in file
    call MPI_FILE_WRITE(restartfileid, buf,bufsize(irank), & 
                            MPI_double_precision, MPI_STATUS_IGNORE, ierr)
    !Close file on all processors
    call MPI_FILE_CLOSE(restartfileid, ierr)
    !This barrier is needed in order to get the correct file size in the next write
    call MPI_BARRIER(MD_COMM, ierr) 
 

    !----------------Write header-----------------------
    !Written at the end for performance and simplicity reasons 
    !(See Gropp, lusk & Thakur Using MPI-2)

    !Write the header with one processor only
    if (irank .eq. iroot) then

        if (.not. interim_output) then
            call MPI_file_open(MPI_COMM_SELF,trim(prefix_dir)//'results/final_state', & 
                 MPI_MODE_WRONLY, MPI_INFO_NULL, restartfileid, ierr)
            print*, "Output final"
        else
            call get_Timestep_FileName(iter, trim(prefix_dir)//'results/interim_state', filename_incno)
            call MPI_file_open(MPI_COMM_SELF,filename_incno, & 
                 MPI_MODE_WRONLY, MPI_INFO_NULL, restartfileid, ierr)

        endif

                
        if (ierr .ne. 0) then 
            write(0,*) "MD parallel_io: error in MPI_File open"
        endif

        call MPI_File_get_size(restartfileid,filesize,ierr)
        disp = filesize

        call MPI_file_set_view(restartfileid,disp,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierr)

        call MPI_File_write(restartfileid,sum(procnp)   ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,initialnunits ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        if (iter .lt. Nsteps ) then
            call MPI_File_write(restartfileid,iter      ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        else
            call MPI_File_write(restartfileid,Nsteps    ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        endif
        call MPI_File_write(restartfileid,tplot         ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,seed          ,2,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,periodic      ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,potential_flag,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,rtrue_flag    ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,solvent_flag  ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,nmonomers     ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,npx           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,npy           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,npz           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,procnp,size(procnp),MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,proctethernp,size(proctethernp),MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,globaldomain  ,3,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,density       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,rcutoff       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,delta_t       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,elapsedtime   ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,simtime       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,k_c           ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,R_0           ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,eps_pp        ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,eps_ps        ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,eps_ss        ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,delta_rneighbr,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,mie_potential,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,global_numbering,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

        header_pos = filesize ! just in case offset kind is 32 bit, rather improbable these days  !!!
        call MPI_File_write(restartfileid,header_pos,1,MPI_INTEGER8,MPI_STATUS_IGNORE,ierr)
        call MPI_file_close(restartfileid, ierr)

    endif

    if (allocated(rglobal)) deallocate(rglobal) 
    if (allocated(rtetherglobal)) deallocate(rtetherglobal) 
    if (allocated(buf)) deallocate(buf)
    if (allocated(bufsize)) deallocate(bufsize)

end subroutine parallel_io_final_state

!------------------------------------------------------------------------
! Call appropriate vmd writing routine
subroutine parallel_io_write_vmd(intervalno,reccount)
    use computational_constants_MD, only : vmd_outflag, vmd_intervals, & 
                                           iter, initialstep, & 
                                           vmd_skip
	use interfaces, only : error_abort
    implicit none

    integer,intent(out)  :: intervalno, reccount

	integer			:: vmd_iter

    vmd_iter = iter-initialstep+1

    if (vmd_iter.ge.vmd_intervals(1,intervalno) .and. & 
        vmd_iter.lt.vmd_intervals(2,intervalno)) then
        select case(vmd_outflag)
        case(1)
            call parallel_io_vmd(reccount)
        case(2)
            call parallel_io_vmd_sl(reccount)
        case(3)
            call parallel_io_vmd(reccount)
            call parallel_io_vmd_halo(reccount)
        case(4)
            call parallel_io_vmd_true(reccount)
        case default
            call error_abort('Unrecognised vmd_outflag in simulation_record')
        end select
        reccount = reccount + 1
    else if (vmd_iter.ge.vmd_intervals(2,intervalno)) then
        intervalno = intervalno + 1			
    endif

end subroutine parallel_io_write_vmd

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd(recno)
    use module_parallel_io
    use messenger_data_exchange, only : globalSum
    implicit none

    integer, intent(in)             :: recno

    integer                         :: i, datasize
    integer                         :: n,globmolno,ordered_write=0
    integer(kind=MPI_OFFSET_KIND)   :: disp, procdisp
    integer,dimension(:),pointer    :: globalno
    real,dimension(:),allocatable   :: Xbuf, Ybuf, Zbuf
    real,dimension(:),allocatable   :: Xbufglob,Ybufglob,Zbufglob

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp()

    !Determine size of real datatype
    call MPI_type_size(MPI_real,datasize,ierr)

    if (global_numbering .eq. 1) then
        ordered_write = 1
        globalno => glob_no
    endif

    if (potential_flag .eq. 1) then
        ordered_write = 1
        globalno => monomer(:)%glob_no
    endif


    !Load buffers with single precision r and adjust according
    !to processor topology with r = 0 at centre
    select case(ordered_write)
    case(0)

        !Allocate buffers
        allocate(Xbuf(np),Ybuf(np),Zbuf(np))

        Xbuf(:) = r(1,1:np)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
        Ybuf(:) = r(2,1:np)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
        Zbuf(:) = r(3,1:np)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

        procdisp = 0
        !Obtain displacement of each processor using all other procs' np
        do i = 1, irank -1
            procdisp = procdisp + procnp(i)*datasize
        enddo

        !Open file on all processors
        call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp.dcd',      & 
                           MPI_MODE_RDWR + MPI_MODE_CREATE, & 
                           MPI_INFO_NULL, fileid,     ierr)

        !-------------Write X coordinates--------------------
        !Obtain location to write in file

        !If intervals set to zero then full simulation recorded
        !if (Nvmd_intervals.eq.0) then
        !    disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
        !        + procdisp  
        !else
        !Otherwise, calculate number of previous intervals
        disp = recno * nd * globalnp * datasize & !Current iteration
            + procdisp  
        !endif

        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
                               MPI_REAL, 'native', MPI_INFO_NULL, ierr)
        
        !Write information to file
        call MPI_FILE_WRITE_ALL(fileid,     Xbuf,     np, MPI_REAL, & 
                                MPI_STATUS_IGNORE, ierr) 

        !-------------Write Y coordinates--------------------
        disp = disp + globalnp * datasize
        !Obtain location to write in file
        !If intervals set to zero then full simulation recorded
        !if (Nvmd_intervals.eq.0) then
        !    disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
        !        + procdisp  &
        !        + globalnp * datasize               !Y Coordinate location
        !else
            !Otherwise, calculate number of previous intervals
        !    disp = recno * nd * globalnp * datasize & !Current iteration
        !        + procdisp      &
        !        + globalnp * datasize               !Y Coordinate location
        !endif

        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
                               MPI_REAL, 'native', MPI_INFO_NULL, ierr)
        
        !Write information to file
        call MPI_FILE_WRITE_ALL(fileid,     Ybuf,     np, MPI_REAL, & 
                                MPI_STATUS_IGNORE, ierr) 

        !-------------Write Z coordinates--------------------

        !If intervals set to zero then full simulation recorded
        disp = disp + globalnp * datasize
!        if (Nvmd_intervals.eq.0) then
!            disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
!                + procdisp  &
!                + 2 * globalnp * datasize                !Z Coordinate location
!        else
!            !Otherwise, calculate number of previous intervals
!            disp = recno * nd * globalnp * datasize & !Current iteration
!                + procdisp      &
!                + 2 * globalnp * datasize                !Z Coordinate location
!        endif

        !Obtain location to write in file
        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
                               MPI_REAL, 'native', MPI_INFO_NULL, ierr)
        
        !Write information to file
        call MPI_FILE_WRITE_ALL(fileid,     Zbuf,     np, MPI_REAL, & 
                                MPI_STATUS_IGNORE, ierr) 

        !-------------- CLOSE -------------------------------   
        call MPI_FILE_CLOSE(fileid, ierr)
    
    case(1)

        allocate(Xbufglob(globalnp),Ybufglob(globalnp),Zbufglob(globalnp))

        Xbufglob = 0.0           !Initialise to zero so that global array...
        Ybufglob = 0.0           !...may be found by summation in parallel
        Zbufglob = 0.0           !------------------------------------------

        !Build sparse individual "global" buffers according to global molecular ID of each monomer
        do n=1,np
            globmolno           = globalno(n)
            Xbufglob(globmolno) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
            Ybufglob(globmolno) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
            Zbufglob(globmolno) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
        end do

        call globalSum(Xbufglob,globalnp)  !Global summation to complete global buffer
        call globalSum(Ybufglob,globalnp)
        call globalSum(Zbufglob,globalnp)

        if (irank .eq. iroot) then  

            call MPI_FILE_OPEN(MPI_COMM_SELF,trim(prefix_dir)//'results/vmd_temp.dcd', & 
                               MPI_MODE_RDWR + MPI_MODE_CREATE,       & 
                               MPI_INFO_NULL, fileid, ierr)

            !disp =((iter-initialstep+1)/(tplot-1)) * nd * globalnp * datasize   !Current iteration
            disp = recno * nd * globalnp * datasize !& !Current iteration

            !Write X positions---------------------------------------------
            call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
                                   'native', MPI_INFO_NULL, ierr)
            call MPI_FILE_WRITE(fileid,     Xbufglob,     globalnp, MPI_REAL, &      !Write buffer
                                    MPI_STATUS_IGNORE, ierr)

            disp = disp + globalnp*datasize                                          !Update file disp

            !Write Y positions---------------------------------------------
            call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
                                   'native', MPI_INFO_NULL, ierr)
            call MPI_FILE_WRITE(fileid,     Ybufglob,     globalnp, MPI_REAL, &      !Write buffer
                                    MPI_STATUS_IGNORE, ierr)

            disp = disp + globalnp*datasize

            !Write Z positions---------------------------------------------
            call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
                                   'native', MPI_INFO_NULL, ierr)
            call MPI_FILE_WRITE(fileid,     Zbufglob,     globalnp, MPI_REAL, &      !Write buffer
                                    MPI_STATUS_IGNORE, ierr)
        
            call MPI_FILE_CLOSE(fileid, ierr) 

        end if
    
    case default
    end select

    if (allocated(Xbuf)) deallocate(Xbuf)
    if (allocated(Ybuf)) deallocate(Ybuf)
    if (allocated(Zbuf)) deallocate(Zbuf)
    if (allocated(Xbufglob)) deallocate(Xbufglob)
    if (allocated(Ybufglob)) deallocate(Ybufglob)
    if (allocated(Zbufglob)) deallocate(Zbufglob)

end subroutine parallel_io_vmd


!------------------------------------------------------------------------
!Write true (no periodic wrapping) positions of molecules to a file

subroutine parallel_io_vmd_true(recno)
    use module_parallel_io
    use messenger_data_exchange, only : globalSum
    implicit none

    integer, intent(in)             :: recno

    integer                         :: procdisp
    integer                         :: i, datasize
    real,dimension(:),allocatable   :: Xbuf, Ybuf, Zbuf
    real,dimension(:),allocatable   :: Xbufglob,Ybufglob,Zbufglob
    integer                         :: n,globmolno
    integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

    !Determine size of real datatype
    call MPI_type_size(MPI_real,datasize,ierr)

    !Load buffers with single precision r and adjust according
    !to processor topology with r = 0 at centre
    select case(potential_flag)
    case(0)

        !Allocate buffers
        allocate(Xbuf(np),Ybuf(np),Zbuf(np))

        Xbuf(:) = rtrue(1,1:np)
        Ybuf(:) = rtrue(2,1:np)
        Zbuf(:) = rtrue(3,1:np)

        procdisp = 0
        !Obtain displacement of each processor using all other procs' np
        do i = 1, irank -1
            procdisp = procdisp + procnp(i)*datasize
        enddo

        !Open file on all processors
        call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp_true.dcd', & 
                           MPI_MODE_RDWR + MPI_MODE_CREATE, & 
                           MPI_INFO_NULL, fileid, ierr)

        !-------------Write X coordinates--------------------
        !Obtain location to write in file

        !If intervals set to zero then full simulation recorded
        !if (Nvmd_intervals.eq.0) then
        !    disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
        !        + procdisp  
        !else
        !Otherwise, calculate number of previous intervals
        disp = recno * nd * globalnp * datasize & !Current iteration
            + procdisp  
        !endif

        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
                               MPI_REAL, 'native', MPI_INFO_NULL, ierr)
        
        !Write information to file
        call MPI_FILE_WRITE_ALL(fileid,     Xbuf,     np, MPI_REAL, & 
                                MPI_STATUS_IGNORE, ierr) 

        !-------------Write Y coordinates--------------------
        !Obtain location to write in file
        !If intervals set to zero then full simulation recorded
        disp = disp + globalnp * datasize

!        if (Nvmd_intervals.eq.0) then
!            disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
!                + procdisp  &
!                + globalnp * datasize               !Y Coordinate location
!        else
!            !Otherwise, calculate number of previous intervals
!            disp = recno * nd * globalnp * datasize & !Current iteration
!                + procdisp      &
!                + globalnp * datasize               !Y Coordinate location
!        endif

        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
                               MPI_REAL, 'native', MPI_INFO_NULL, ierr)
        !Write information to file
        call MPI_FILE_WRITE_ALL(fileid,     Ybuf,     np, MPI_REAL, & 
                                MPI_STATUS_IGNORE, ierr) 

        !-------------Write Z coordinates--------------------

        !If intervals set to zero then full simulation recorded
        disp = disp + globalnp * datasize
!        if (Nvmd_intervals.eq.0) then
!            disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
!                + procdisp  &
!                + 2* globalnp * datasize                !Z Coordinate location
!        else
!            !Otherwise, calculate number of previous intervals
!            disp = recno * nd * globalnp * datasize & !Current iteration
!                + procdisp      &
!                + 2* globalnp * datasize                !Z Coordinate location
!        endif

        !Obtain location to write in file
        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
                               MPI_REAL, 'native', MPI_INFO_NULL, ierr)
        
        !Write information to file
        call MPI_FILE_WRITE_ALL(fileid,     Zbuf,     np, MPI_REAL, & 
                                MPI_STATUS_IGNORE, ierr) 

        !-------------- CLOSE -------------------------------   
        call MPI_FILE_CLOSE(fileid, ierr) 
    
    case(1)

        !Allocate buffers
        allocate(Xbufglob(globalnp),Ybufglob(globalnp),Zbufglob(globalnp))

        Xbufglob = 0.0           !Initialise to zero so that global array...
        Ybufglob = 0.0           !...may be found by summation in parallel
        Zbufglob = 0.0           !------------------------------------------

        !Build sparse individual "global" buffers according to global molecular ID of each monomer
        do n=1,np
            globmolno           = monomer(n)%glob_no
            Xbufglob(globmolno) = rtrue(1,n)
            Ybufglob(globmolno) = rtrue(2,n)
            Zbufglob(globmolno) = rtrue(3,n)
        end do

        call globalSum(Xbufglob,globalnp)  !Global summation to complete global buffer
        call globalSum(Ybufglob,globalnp)
        call globalSum(Zbufglob,globalnp)
    
        call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp_true.dcd', & 
                           MPI_MODE_RDWR + MPI_MODE_CREATE,       & 
                           MPI_INFO_NULL, fileid, ierr)

        disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize   !Current iteration

        !Write X positions---------------------------------------------
        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
                               'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(fileid,     Xbufglob,     globalnp, MPI_REAL, &      !Write buffer
                                MPI_STATUS_IGNORE, ierr)

        disp = disp + globalnp*datasize                                          !Update file disp

        !Write Y positions---------------------------------------------
        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
                               'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(fileid,     Ybufglob,     globalnp, MPI_REAL, &      !Write buffer
                                MPI_STATUS_IGNORE, ierr)

        disp = disp + globalnp*datasize

        !Write Z positions---------------------------------------------
        call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
                               'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(fileid,     Zbufglob,     globalnp, MPI_REAL, &      !Write buffer
                                MPI_STATUS_IGNORE, ierr)
            
        call MPI_FILE_CLOSE(fileid, ierr) 
        
    case default
    end select

end subroutine parallel_io_vmd_true

!------------------------------------------------------------------------
!Write solid and liquid positions of molecules to 2 seperate files

subroutine parallel_io_vmd_sl(recno)
    use module_parallel_io
    implicit none

    integer, intent(in)             :: recno

    integer                         :: procdisp
    integer                         :: i,n, datasize
    real,dimension(:),allocatable   :: Xbuf, Ybuf, Zbuf
    integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

    !Allocate buffers
    allocate(Xbuf(np),Ybuf(np),Zbuf(np))

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

    !Determine size of real datatype
    call MPI_type_size(MPI_real,datasize,ierr)

    !Load buffers with single precision r and adjust according
    !to processor topology with r = 0 at centre
    do n = 1, np
        select case(tag(n))
        case(0,4)   !Liquid Molecules
            Xbuf(n) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
            Ybuf(n) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
            Zbuf(n) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
        case default    !Solid molecules
            Xbuf(n) = -halfdomain(1)
            Ybuf(n) = -halfdomain(2)
            Zbuf(n) = -halfdomain(3)
        end select
    enddo

    procdisp = 0
    !Obtain displacement of each processor using all other procs' np
    do i = 1, irank -1
        procdisp = procdisp + procnp(i)*datasize
    enddo

    !================================
    !    Write liquid Molecules =
    !================================
    !Open file on all processors
    call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_liquid_temp.dcd', & 
        MPI_MODE_RDWR + MPI_MODE_CREATE, & 
        MPI_INFO_NULL, fileid, ierr)

    !-------------Write X coordinates--------------------

    !If intervals set to zero then full simulation recorded
    !if (Nvmd_intervals.eq.0) then
    !    disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
    !        + procdisp  
    !else
    !Otherwise, calculate number of previous intervals
    disp = recno * nd * globalnp * datasize & !Current iteration
            + procdisp  
    !endif


    !print*, irank, 'x disp', disp

    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
        MPI_REAL, 'native', MPI_INFO_NULL, ierr)
    
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
        MPI_STATUS_IGNORE, ierr) 

    !-------------Write Y coordinates--------------------

    !If intervals set to zero then full simulation recorded
    disp = disp + globalnp * datasize
!    if (Nvmd_intervals.eq.0) then
!        disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
!            + procdisp  &
!            + globalnp * datasize               !Y Coordinate location
!    else
!        !Otherwise, calculate number of previous intervals
!        disp = recno * nd * globalnp * datasize & !Current iteration
!                + procdisp  &
!                + globalnp * datasize               !Y Coordinate location
!    endif

    !print*, irank, 'y disp', disp
    
    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
        MPI_REAL, 'native', MPI_INFO_NULL, ierr)
        
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
        MPI_STATUS_IGNORE, ierr) 

    !-------------Write Z coordinates--------------------

    !If intervals set to zero then full simulation recorded
    disp = disp + globalnp * datasize
!    if (Nvmd_intervals.eq.0) then
!        disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
!                + procdisp  &
!                + 2 * globalnp * datasize               !Z Coordinate location
!    else
!        !Otherwise, calculate number of previous intervals
!        disp = recno * nd * globalnp * datasize & !Current iteration
!                + procdisp  &
!                + 2 * globalnp * datasize               !Z Coordinate location
!    endif

    !print*, irank, 'z disp', disp

    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
        MPI_REAL, 'native', MPI_INFO_NULL, ierr)
    
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
        MPI_STATUS_IGNORE, ierr) 

    !Close file on all processors
    call MPI_FILE_CLOSE(fileid, ierr) 

    !Load buffers with single precision r and adjust according
    !to processor topology with r = 0 at centre
    do n = 1, np
        select case(tag(n))
        case(0, 4)  !Liquid Molecules
            Xbuf(n) = -halfdomain(1)
            Ybuf(n) = -halfdomain(2)
            Zbuf(n) = -halfdomain(3)
        case default    !Solid molecules
            Xbuf(n) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
            Ybuf(n) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
            Zbuf(n) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
        end select
    enddo

    !================================
    !    Write Solid Molecules  =
    !================================

    !Open file on all processors
    call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_solid_temp.dcd', & 
        MPI_MODE_RDWR + MPI_MODE_CREATE, & 
        MPI_INFO_NULL, fileid, ierr)

    !-------------Write X coordinates--------------------
        
    !If intervals set to zero then full simulation recorded
    if (Nvmd_intervals.eq.0) then
        !Obtain location to write in file
        disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
            + procdisp  
    else
        !Otherwise, calculate number of previous intervals
        disp = recno * nd * globalnp * datasize & !Current iteration
                + procdisp  
    endif

    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
        MPI_REAL, 'native', MPI_INFO_NULL, ierr)
    
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
        MPI_STATUS_IGNORE, ierr) 

    !-------------Write Y coordinates--------------------

    !If intervals set to zero then full simulation recorded
    if (Nvmd_intervals.eq.0) then
        !Obtain location to write in file
        disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
            + procdisp  &
            + globalnp * datasize               !Y Coordinate location
    else
        !Otherwise, calculate number of previous intervals
        disp = recno * nd * globalnp * datasize & !Current iteration
                + procdisp  &
                + globalnp * datasize               !Y Coordinate location
    endif

    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
        MPI_REAL, 'native', MPI_INFO_NULL, ierr)
    
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
        MPI_STATUS_IGNORE, ierr) 

    !-------------Write Z coordinates--------------------

    !If intervals set to zero then full simulation recorded
    if (Nvmd_intervals.eq.0) then
        !Obtain location to write in file
        disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
                + procdisp  &
                + 2 * globalnp * datasize               !Z Coordinate location
    else
        !Otherwise, calculate number of previous intervals
        disp = recno * nd * globalnp * datasize & !Current iteration
                + procdisp  &
                + 2 * globalnp * datasize               !Z Coordinate location
    endif


    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
        MPI_REAL, 'native', MPI_INFO_NULL, ierr)
    
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
        MPI_STATUS_IGNORE, ierr) 

    !Close file on all processors
    call MPI_FILE_CLOSE(fileid, ierr) 

    if (allocated(Xbuf)) deallocate(Xbuf)
    if (allocated(Ybuf)) deallocate(Ybuf)
    if (allocated(Zbuf)) deallocate(Zbuf)

end subroutine parallel_io_vmd_sl

!------------------------------------------------------------------------
!Write positions of molecules to a file using one single write instead of 3

subroutine parallel_io_vmd_optimised(recno)
    use module_parallel_io
    implicit none

    integer, intent(in)             :: recno

    integer                         :: procdisp
    integer                         :: i, datasize
    real,dimension(:),allocatable   :: buf
    integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

    allocate(buf(np*3))

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

    !Determine size of real datatype
    call MPI_type_size(MPI_real,datasize,ierr)

    !Load buffers with single precision r and adjust according
    !to processor topology with r = 0 at centre
    do i=1,np
        buf(3*(i-1)+1) = r(1,i)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
        buf(3*(i-1)+2) = r(2,i)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
        buf(3*(i-1)+3) = r(3,i)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
    enddo

    procdisp = 0
    !Obtain displacement of each processor using all other procs' np
    do i = 1, irank -1
        procdisp = procdisp + procnp(i)*datasize
    enddo

    !Open file on all processors
    call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_temp.dcd', & 
                        MPI_MODE_RDWR + MPI_MODE_CREATE, & 
                        MPI_INFO_NULL, fileid, ierr)

    !-------------Write XYZ coordinates--------------------

    !If intervals set to zero then full simulation recorded
    if (Nvmd_intervals.eq.0) then
        !Obtain location to write in file
        disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1)*nd*globalnp*datasize & !Current iteration
            + procdisp
    else
        !Otherwise, calculate number of previous intervals
        disp = recno * nd * globalnp * datasize & !Current iteration
            + procdisp
    endif

    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
        MPI_REAL, 'native', MPI_INFO_NULL, ierr)
    
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, buf, np*3, MPI_REAL, & 
        MPI_STATUS_IGNORE, ierr) 

    !Close file on all processors
    call MPI_FILE_CLOSE(fileid, ierr) 

end subroutine parallel_io_vmd_optimised

!------------------------------------------------------------------------
!Write positions of molecules in halo to a file

subroutine parallel_io_vmd_halo(recno)
    use module_parallel_io
    implicit none

    integer, intent(in) :: recno

    call error_abort("Cannot print vmd halos in parallel simulation - run in serial")

end subroutine parallel_io_vmd_halo

!-----------------------------------------------------------------------------
! Write a psf file with information in the form
! Note this uses a poor mans MPI where each proc writes its own psf file
! and system cat is called at the end. If this fails, the files can be 
! concatted after. IT IS NOT TRIVIAL TO WRITE ASCII with MPI!!!
!=============================================================================

!PSF CMAP
!       6 !NTITLE
! REMARKS original generated structure x-plor psf file
! REMARKS 2 patches were applied to the molecule.
! REMARKS topology top_all27_prot_lipid.inp
! REMARKS segment U { first NTER; last CTER; auto angles dihedrals }
! REMARKS defaultpatch NTER U:1
! REMARKS defaultpatch CTER U:76

!    1231 !NATOM
!       1 U    1    MET  N    NH3   -0.300000       14.0070           0
!       2 U    1    MET  HT1  HC     0.330000        1.0080           0
!       3 U    1    MET  HT2  HC     0.330000        1.0080           0
!       4 U    1    MET  HT3  HC     0.330000        1.0080           0
!       5 U    1    MET  CA   CT1    0.210000       12.0110           0
!       6 U    1    MET  HA   HB     0.100000        1.0080           0
!       7 U    1    MET  CB   CT2   -0.180000       12.0110           0
! Fields include: 
!   molno, 
!   segment name, 
!   Residue ID, 
!   residue name, 
!   atom name, 
!   atom type, 
!   charge, 
!   mass 
!   A "0" for no reason!
!=============================================================================

subroutine parallel_io_psf()
    use module_parallel_io
    use module_set_parameters, only : moltype_names, mass_lookup, epsilon_lookup
    use librarymod, only : get_new_fileunit
    implicit none

    integer                 :: n, mt, unitno, startmol
    character(2)           :: nprocstr
    character(256)          :: filename, cmd

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp()

    startmol = 0
    do n = 1, irank -1
        startmol = startmol + procnp(n)
    enddo

    unitno = get_new_fileunit()
    if (irank .lt. 10) then
        write(filename,'(a,i1)') trim(prefix_dir)//'results/vmd_out.psf.',irank
    elseif (irank .lt. 100) then
        write(filename,'(a,i2)') trim(prefix_dir)//'results/vmd_out.psf.',irank
    elseif (irank .lt. 1000) then
        write(filename,'(a,i3)') trim(prefix_dir)//'results/vmd_out.psf.',irank
    elseif (irank .lt. 10000) then
        write(filename,'(a,i4)') trim(prefix_dir)//'results/vmd_out.psf.',irank
    else
       stop "Error in parallel_io_psf -- irank > 10000"
    endif
    open(unit=unitno, file=trim(filename), status='replace',action='write')
    if (irank .eq. iroot) then
        write(unitno,'(a)') "PSF CMAP"
        write(unitno,'(a)') ""
        write(unitno,'(a)') "      1 !NTITLE"
        write(unitno,'(a)') " REMARKS psf file written by crompulenceMD code"
        write(unitno,'(a)') ""
        write(unitno,'(i8,a)') globalnp, " !NATOM"
    endif

    do n = 1, np
        mt = moltype(n)
        write(unitno,'(i8,a5,i5,3a5,f12.7,f12.5,i4)') &
                 startmol + n, &            ! molno
                 trim(moltype_names(mt)), & ! segment name
                  1,  &                     ! Residue ID
                 trim(moltype_names(mt)), & ! residue name
                 trim(moltype_names(mt)), & ! atom name
                 trim(moltype_names(mt)), & ! atom type
                 epsilon_lookup(mt,mt), &   ! charge
                 mass_lookup(mt), &         ! mass
                  0                         ! A "0" for no reason!
    enddo
    !Add an empty set of bond information needed
    !for MDTraj 
    write(unitno,'(a)') ""
    write(unitno,'(i8,a)') 0, " !NBOND: bonds"

    close(unitno)


    call MPI_Barrier(MD_COMM,ierr)

    !This is a disgusting hack to concat each processors files on the commandline
    if (irank .eq. iroot) then
        write(cmd,'(3a)') "cat ", trim(prefix_dir)//"results/vmd_out.psf.* > ", & 
                            trim(prefix_dir)//"results/vmd_out.psf"
        call system(cmd)
        if (nproc .gt. 99) then 
            print*, "Warning, manually concat results/vmd_out.psf.* files"
        else if (nproc .gt. 1) then
            write (nprocstr, "(i2)") nproc
            write(cmd,'(4a)') "for i in {2..", trim(nprocstr)//"}; do cat ", &
                                       trim(prefix_dir)//"results/vmd_out.psf.$i >> ", &
                                       trim(prefix_dir)//"results/vmd_out.psf; done;"
            print*, cmd
            call system(cmd)
        else if (nproc .eq. 1) then
            write (nprocstr, "(i2)") nproc
            write(cmd,'(4a)') "cat ", trim(prefix_dir)//"results/vmd_out.psf.1 >> ", &
                                      trim(prefix_dir)//"results/vmd_out.psf"
            print*, cmd
        endif
    endif

end subroutine parallel_io_psf

!subroutine parallel_io_psf()
!    use module_parallel_io
!    use module_set_parameters, only : moltype_names, mass_lookup, epsilon_lookup
!    use librarymod, only : get_new_fileunit
!    implicit none

!    integer         :: n, mt, unitno
!    integer(kind=MPI_OFFSET_KIND) :: headersize, disp, datasize, linesize
!    character(48)  :: headerlines(6)
!    character(256) :: line(np)


!    !Build array of number of particles on neighbouring
!    !processe's subdomain on current proccess
!    call globalGathernp

!    !Write file header on root processor
!    unitno = get_new_fileunit()
!    if (irank .eq. iroot) then
!        headerlines(1) = "PSF CMAP"
!        headerlines(2) = ""
!        headerlines(3) = "      1 !NTITLE"
!        headerlines(4) = " REMARKS psf file written by crompulenceMD code"
!        headerlines(5) = ""
!        write(headerlines(6),'(i8,a)') np, " !NATOM"

!        open(unit=unitno, file=trim(prefix_dir)//'results/vmd_out.psf',status='replace',action='write')
!        headersize = 0
!        do n=1,size(headerlines,1)
!            print*, trim(headerlines(n))
!            write(unitno,'(a)') trim(headerlines(n))
!            headersize = headersize + len_trim(headerlines(n))
!        enddo
!    endif
!    call globalbroadcast(headersize,1,iroot)
!   
!    !Pack data into line character
!    do n = 1, np
!        mt = moltype(n)
!        write(line(n),'(i5,a,i5,3a,f12.7,f12.5,i4,a)') &
!                  n, &                      ! molno
!                 trim(moltype_names(mt)), & ! segment name
!                  1,  &                     ! Residue ID
!                 trim(moltype_names(mt)), & ! residue name
!                 trim(moltype_names(mt)), & ! atom name
!                 trim(moltype_names(mt)), & ! atom type
!                 epsilon_lookup(mt,mt), &   ! charge
!                 mass_lookup(mt), &         ! mass
!                  0, &                      ! A "0" for no reason!
!                 " \n"                       !Newline
!        
!    enddo

!    disp = headersize
!    linesize = len_trim(line(1))
!    datasize = linesize*np
!    !Obtain displacement of each processor using all other procs' np
!    do n = 1, irank -1
!        disp = disp + procnp(n)*datasize
!    enddo


!    print*, irank, headersize, linesize, datasize, disp

!    !Open file on all processors
!    call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_out.psf',      & 
!                       MPI_MODE_RDWR + MPI_MODE_CREATE, & 
!                       MPI_INFO_NULL, fileid,     ierr)

!    !Loop through each line and write
!    do n = 1, np
!        disp = disp + linesize
!        call MPI_FILE_SET_VIEW(fileid, disp, MPI_CHAR,          & 
!                               MPI_CHAR, 'native', MPI_INFO_NULL, ierr)
!        
!        !Write information to file
!        call MPI_FILE_WRITE_ALL(fileid, trim(line(n)), linesize, MPI_CHAR, & 
!                                MPI_STATUS_IGNORE, ierr)
!    enddo

!    !-------------- CLOSE -------------------------------   
!    call MPI_FILE_CLOSE(fileid, ierr)

!end subroutine parallel_io_psf


!-----------------------------------------------------------------------------
! Write cylinder molecules and properties to a file
subroutine parallel_io_write_cylinders()
    use concentric_cylinders, only: r_oo, r_io, r_oi, r_ii, cyl_file
    use module_parallel_io
    use messenger, only: globalise
    implicit none

    integer :: i, datasize
    integer :: n, pos
    integer(selected_int_kind(18))  :: header_pos
    real :: rglob(3)
    real, dimension(:), allocatable :: buf
    integer(kind=MPI_OFFSET_KIND) :: disp, procdisp, filesize

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

    !Determine size of real datatype
    call MPI_type_size(MPI_real,datasize,ierr)

    !Load buffers with single precision r and adjust according
    !to processor topology with r = 0 at centre
    !Allocate buffers
    allocate(buf(3*np))

    pos = 1
    do n = 1,np
        rglob(:) = real(globalise(r(:,n)))
        buf(pos) = rglob(1)
        buf(pos+1) = rglob(2)
        buf(pos+2) = rglob(3)
        pos = pos + 3
    end do

    procdisp = 0
    !Obtain displacement of each processor using all other procs' np
    do i = 1, irank -1
        procdisp = procdisp + procnp(i)*3*datasize
    enddo

    !Open file on all processors
    call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/cylinders', & 
                       MPI_MODE_RDWR + MPI_MODE_CREATE, MPI_INFO_NULL, &
                       fileid, ierr)

    !-------------Write X coordinates--------------------
    !Obtain location to write in file
    disp = procdisp
    call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, MPI_REAL, 'native', & 
                           MPI_INFO_NULL, ierr)
    !Write information to file
    call MPI_FILE_WRITE_ALL(fileid, buf, 3*np, MPI_REAL, MPI_STATUS_IGNORE, &
                            ierr) 
    !-------------- CLOSE -------------------------------   
    call MPI_FILE_CLOSE(fileid, ierr) 
    !Barrier required to get the correct file size for next write
    call MPI_BARRIER(MD_COMM, ierr) 

    !----------------Write header-----------------------
    !Written at the end for performance and simplicity reasons 
    !(See Gropp, lusk & Thakur Using MPI-2)

    !Write the header with one processor only
    if (irank .eq. iroot) then

        call MPI_file_open(MPI_COMM_SELF,trim(prefix_dir)//'results/cylinders', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, fileid, ierr)
                
        if (ierr .ne. 0) then 
            write(0,*) "MD parallel_io: error in MPI_File open"
        endif
        call MPI_FILE_GET_SIZE(fileid,filesize,ierr)
        disp = filesize
        call MPI_FILE_SET_VIEW(fileid,disp,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierr)

        ! Write integers first
        call MPI_FILE_WRITE(fileid,globalnp,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_WRITE(fileid,globaldomain,3,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_WRITE(fileid,r_oo,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_WRITE(fileid,r_io,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_WRITE(fileid,r_oi,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_WRITE(fileid,r_ii,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    
        header_pos = filesize ! just in case offset kind is 32 bit, rather improbable these days  !!!
        call MPI_File_write(fileid,header_pos,1,MPI_INTEGER8,MPI_STATUS_IGNORE,ierr)
        call MPI_file_close(fileid,ierr)

    endif

end subroutine parallel_io_write_cylinders
!-----------------------------------------------------------------------------
! Write value of last output iteration

subroutine update_simulation_progress_file()
    use module_parallel_io
    use librarymod, only: get_new_fileunit
    implicit none

    integer :: fileunit

    if (irank .eq. iroot) then
    	fileunit = get_new_fileunit()
        open (unit=fileunit, file=trim(prefix_dir)//"results/simulation_progress")
        write(fileunit,*) iter
        close(fileunit,status='keep')
    endif

end subroutine update_simulation_progress_file


module field_io

contains

!-----------------------------------------------------------------------------
! Record mass in a slice through the domain
subroutine mass_slice_io(ixyz)
    use module_parallel_io
    use messenger
    use interfaces
    implicit none

    integer                         :: ixyz,jxyz,kxyz
    integer                         :: ijkblock(3)
    integer                         :: slicefileid, int_datasize
    integer(kind=MPI_OFFSET_KIND)   :: disp

    !Get two directions orthogonal to slice direction
    kxyz = mod(ixyz,3)+1
    jxyz = mod(ixyz+1,3)+1

    ijkblock = (/iblock,jblock,kblock/)

    !Sum over all bins using directional sub communicators and gather on {ijk}block=1
    call SubcommSum(slice_mass, nbins(ixyz), jxyz)
    call SubcommSum(slice_mass, nbins(ixyz), kxyz)

    !Only root processor in each directional subcomm writes data
    if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

        !Determine size of datatypes
        call MPI_type_size(MPI_Integer,int_datasize,ierr)

        !Only processors on directional subcomm write
        call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'./results/mslice', & 
                            MPI_MODE_WRONLY+ MPI_MODE_CREATE , & 
                            MPI_INFO_NULL, slicefileid, ierr)

        !Obtain displacement of current record
        disp =   (get_iter()/(tplot*Nmass_ave) - 1)       &       !Current iteration
               * gnbins(ixyz)*int_datasize  &       !Global record size
               + nbins(ixyz)*int_datasize*(ijkblock(ixyz)-1)        !Processor location

        call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_INTEGER, & 
                    MPI_INTEGER, 'native', MPI_INFO_NULL, ierr)

        call MPI_FILE_WRITE_ALL(slicefileid,slice_mass,nbins(ixyz), MPI_INTEGER, & 
                        MPI_STATUS_IGNORE, ierr)

        call MPI_FILE_CLOSE(slicefileid, ierr)

    endif 

end subroutine mass_slice_io

!---------------------------------------------------------------------------------
! Record mass in 3D bins throughout domain
subroutine mass_bin_io(CV_mass_out,io_type)
    use module_parallel_io
    use messenger_bin_handler, only : swaphalos
    !use CV_objects, only : CVcheck_mass, CV_debug
    implicit none

    real(kind(0.d0)), dimension(:,:,:), intent(inout) :: CV_mass_out
    !integer, intent(in)    :: CV_mass_out(nbinso(1),nbinso(2),nbinso(3))
    character(4),intent(in) :: io_type

    real(kind(0.d0)), dimension(:,:,:,:), allocatable :: CVmasscopy
    !integer :: CVmasscopy(nbinso(1),nbinso(2),nbinso(3))
    integer :: m,nresults
    character(30) :: filename, outfile

    !Work out correct filename for i/o type
    write(filename, '(a9,a4)' ) 'results/m', io_type
    outfile = trim(prefix_dir)//filename

    !Swap halo surface fluxes to get correct values for all cells
    nresults = 1

    !Copy CV_mass_out so it is not changed
    allocate(CVmasscopy(nbinso(1),nbinso(2),nbinso(3),1))
    CVmasscopy(:,:,:,1) = CV_mass_out
    call swaphalos(CVmasscopy,nbinso(1),nbinso(2),nbinso(3),nresults)
    CV_mass_out = CVmasscopy(:,:,:,1)
    !deallocate(CVmasscopy)

    !Calculate record number timestep
    if (io_type .eq. 'snap') then
        select case(CV_conserve)
        case(0)
            m = get_iter()/(tplot*Nmflux_ave) + 1 !Initial snapshot taken
        case(1)
            m = get_iter()/(Nmflux_ave) + 1 !Initial snapshot taken
            !if (CV_debug .ne. 0) then
            !    call CVcheck_mass%update_dXdt(CV_mass_out)
            !endif
        case default
            call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')   
        end select
    else
        m = get_iter()/(tplot*Nmass_ave)
    endif

    !Write mass to file
    call write_arrays(CVmasscopy, nresults, outfile, m)
    deallocate(CVmasscopy)

end subroutine mass_bin_io


!---------------------------------------------------------------------------------
! Record momentum in a slice through the domain
subroutine momentum_slice_io(ixyz)
    use module_parallel_io
    use messenger
    implicit none

    integer                         :: ixyz,jxyz,kxyz,ijkblock(3)
    integer                         :: slicefileid, dp_datasize
    integer,dimension(3)            :: idims
    integer(kind=MPI_OFFSET_KIND)   :: disp

    !Write mass
    call mass_slice_io(ixyz)

    !Get two directions orthogonal to slice direction
    kxyz = mod(ixyz,3)+1
    jxyz = mod(ixyz+1,3)+1
    idims(1) = npx; idims(2) = npy; idims(3) = npz

    ijkblock = (/iblock,jblock,kblock/)

    !Sum over all bins using directional sub communicators and gather on root
    call SubcommSum(slice_momentum(:,1), nbins(ixyz), jxyz)
    call SubcommSum(slice_momentum(:,1), nbins(ixyz), kxyz)
    call SubcommSum(slice_momentum(:,2), nbins(ixyz), jxyz)
    call SubcommSum(slice_momentum(:,2), nbins(ixyz), kxyz)
    call SubcommSum(slice_momentum(:,3), nbins(ixyz), jxyz)
    call SubcommSum(slice_momentum(:,3), nbins(ixyz), kxyz)

    !Only root processor in each directional subcomm writes data
    if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

        call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

        !Only processors on directional subcomm write
        call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'./results/vslice', & 
                   MPI_MODE_WRONLY + MPI_MODE_CREATE , & 
                   MPI_INFO_NULL, slicefileid, ierr)

        !Obtain displacement of x record
        disp =   (get_iter()/(tplot*Nmass_ave) - 1)       &   !Current iteration
               * nd*gnbins(ixyz)*dp_datasize    &   !times record size
               + nbins(ixyz)*dp_datasize*(ijkblock(ixyz)-1)     !Processor location

        call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
                    MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,1),nbins(ixyz), & 
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

        !Obtain displacement of y record
        disp =   (get_iter()/(tplot*Nmass_ave) - 1)       &   !Current iteration
               * nd*gnbins(ixyz)*dp_datasize    &   !Record size
               + nbins(ixyz)*dp_datasize*(ijkblock(ixyz)-1) &   !Processor location
               + nbins(ixyz)*dp_datasize*idims(ixyz)        !after x data 

        call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
                    MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,2),nbins(ixyz), & 
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

        !Obtain displacement of z record
        disp =   (get_iter()/(tplot*Nmass_ave) - 1)       &   !Current iteration
               * nd*gnbins(ixyz)*dp_datasize    &   !Record size
               + nbins(ixyz)*dp_datasize*(ijkblock(ixyz)-1) &   !Processor location
               + 2*nbins(ixyz)*dp_datasize*idims(ixyz)      !after x & y data 

        call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
                    MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,3),nbins(ixyz), & 
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

        call MPI_FILE_CLOSE(slicefileid, ierr)

    endif 

end subroutine momentum_slice_io

!------------------------------------------------------------------------
!A routine with each proc writing its own bins in binary

subroutine momentum_bin_io(CV_mass_out,CV_momentum_out,io_type)
    use module_parallel_io
    use CV_objects, only : CVcheck_momentum,CV_debug
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer                 :: m,nresults
    !integer, intent(in)        :: CV_mass_out(nbinso(1),nbinso(2),nbinso(3))
    !real(kind(0.d0)), intent(in) :: CV_momentum_out(nbinso(1),nbinso(2),nbinso(3),nd)
    real(kind(0.d0)), dimension(:,:,:), intent(inout) :: CV_mass_out
    real(kind(0.d0)), dimension(:,:,:,:), intent(inout) :: CV_momentum_out

    character(4)            :: io_type
    character(30)           :: filename,outfile

    !Write mass bins
    if (io_type .ne. 'snap') then
        call mass_bin_io(CV_mass_out,io_type)
    endif

    !Work out correct filename for i/o type
    write(filename, '(a9,a4)' ) 'results/v', io_type
    outfile = trim(prefix_dir)//filename

    ! Swap Halos
    nresults = nd
    call swaphalos(CV_momentum_out,nbinso(1),nbinso(2),nbinso(3),nresults)

    !Setup arrays
    if (io_type .eq. 'snap') then
        !CV_momentum_out = CV_momentum_out / (tplot*Nvflux_ave)
        select case(CV_conserve)
        case(0)
            m = get_iter()/(tplot*Nvflux_ave) + 1 !Initial snapshot taken
        case(1)
            m = get_iter()/(Nvflux_ave) + 1 !Initial snapshot taken
            !Create copy of previous timestep Control Volume mass and calculate time evolution
            if (CV_debug .ne. 0) then
                call CVcheck_momentum%update_dXdt(CV_momentum_out(:,:,:,:))
            endif
        case default
            call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')   
        end select
    else
        !CV_momentum_out = CV_momentum_out / (tplot*Nvel_ave)
        m = get_iter()/(tplot*Nvel_ave)
    endif

    !Write out arrays
    call write_arrays(CV_momentum_out,nresults,outfile,m)

end subroutine momentum_bin_io


!------------------------------------------------------------------------
!A routine with each proc writing its own bins in binary

subroutine centre_of_mass_bin_io(centre_of_mass)
    use module_parallel_io, only : write_arrays, nd, iter, Ncom_ave, & 
                                   initialstep,tplot, nbinso, prefix_dir, &
                                   get_iter
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer                 :: m, nresults
    real(kind(0.d0)), dimension(:,:,:,:), intent(inout) :: centre_of_mass

    ! Swap Halos
    nresults = nd
    call swaphalos(centre_of_mass,nbinso(1),nbinso(2),nbinso(3),nresults)
    m = get_iter()/(tplot*Ncom_ave)

    !Write out arrays
    call write_arrays(centre_of_mass,nresults,trim(prefix_dir)//'results/combin',m)

end subroutine centre_of_mass_bin_io


!---------------------------------------------------------------------------------
! Record velocity PDF in a slice through the domain
subroutine velocity_PDF_slice_io(ixyz,pfdx,pfdy,pfdz)
    use module_parallel_io
    use messenger
    implicit none

    integer,intent(in)                  :: ixyz
    integer,dimension(:,:),intent(in)   :: pfdx,pfdy,pfdz

    integer                             :: jxyz,kxyz,ijkblock(3),nxpz(3),datasize
    integer                             :: xmin,xmax,ymin,ymax,zmin,zmax
    integer                             :: slicefileid, int_datasize
    integer,dimension(3)                :: idims
    integer(kind=MPI_OFFSET_KIND)       :: disp
    integer,dimension(:),allocatable    :: pdfxyz

    !Get two directions orthogonal to slice direction
    kxyz = mod(ixyz,3)+1
    jxyz = mod(ixyz+1,3)+1
    idims(1) = npx; idims(2) = npy; idims(3) = npz

    ijkblock = (/iblock,jblock,kblock/)
    nxpz = (/ npx, npy, npz /)
    !Sum over all bins using directional sub communicators and gather on root
    datasize =size(pfdx)+size(pfdy)+size(pfdz)
    allocate(pdfxyz(datasize)); pdfxyz = 0

    xmin = 1; xmax = size(pfdx)
    ymin = xmax+1; ymax = size(pfdx)+size(pfdy)
    zmin = ymax+1; zmax = size(pfdx)+size(pfdy)+size(pfdz)

    pdfxyz(xmin:xmax) = reshape(pfdx(:,:),(/ size(pfdx) /))
    pdfxyz(ymin:ymax) = reshape(pfdy(:,:),(/ size(pfdy) /))
    pdfxyz(zmin:zmax) = reshape(pfdz(:,:),(/ size(pfdz) /))

    !print'(a,10i10)', 'b4', iblock,jblock,kblock, sum(pdfxyz),sum(pfdx),sum(pfdy),sum(pfdz), sum(pdfxyz(xmin:xmax)),sum(pdfxyz(ymin:ymax)),sum(pdfxyz(zmin:zmax))
    call SubcommSumIntVect(pdfxyz,datasize, jxyz)
    call SubcommSumIntVect(pdfxyz,datasize, kxyz)
    !print'(a,10i10)', 'after', iblock,jblock,kblock, sum(pdfxyz),sum(pfdx),sum(pfdy),sum(pfdz), sum(pdfxyz(xmin:xmax)),sum(pdfxyz(ymin:ymax)),sum(pdfxyz(zmin:zmax))
    !print*, size(pfdx,1),size(pfdx,2),size(pfdx),size(pdfxyz)

    !Only root processor in each directional subcomm writes data
    if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then
        
        call MPI_type_size(MPI_integer,int_datasize,ierr)

        !Only processors on directional subcomm write
        call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'./results/vPDF', & 
                           MPI_MODE_WRONLY + MPI_MODE_CREATE , & 
                           MPI_INFO_NULL, slicefileid, ierr)

        !Obtain displacement of x record
        disp =   (get_iter()/(tplot*Nvpdf_ave) - 1)       &       !Current iteration
                    * datasize*int_datasize                         &       !times record size
                    + datasize*int_datasize*(ijkblock(ixyz)-1)              !Processor location

        !print'(15i10)',iblock, jblock,kblock,ixyz,jxyz,kxyz,nxpz(ixyz),disp,datasize,int_datasize,nbins(ixyz),gnbins(ixyz), & 
        !               ((iter-initialstep+1)/(tplot*Nvpdf_ave)-1),datasize*int_datasize, datasize*int_datasize*(ijkblock(ixyz)-1) 

        call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_INTEGER, & 
                    MPI_INTEGER, 'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(slicefileid,pdfxyz,datasize, & 
                    MPI_INTEGER,MPI_STATUS_IGNORE, ierr)

        call MPI_FILE_CLOSE(slicefileid, ierr)

    endif 

end subroutine velocity_PDF_slice_io

!------------------------------------------------------------------------------
! Cylindrical polar mass bins
subroutine mass_bin_cpol_io(mass_out)
    use module_parallel_io, only: write_zplane, get_iter
    use concentric_cylinders, only: cpol_binso
    use physical_constants_MD, only: nd
    use computational_constants_MD, only: iter, initialstep, tplot, &
                                          Nmass_ave, prefix_dir
    use messenger, only: icomm_grid,iblock,jblock,plane_comm
    implicit none

    !integer, intent(inout) :: mass_out(cpol_binso(1), &
    !                                   cpol_binso(2), &
    !                                   cpol_binso(3))
    
    real(kind(0.d0)), dimension(:,:,:), intent(inout) :: mass_out

    character(200) :: mfile!, vfile
    integer :: m, ierr

    ! Z-Plane global sums of r/theta bins, everyone has a copy of
    ! r/theta bins in their own z-plane
    call ZPlaneReduceMass

    ! Bottom corner z-line of processors write to file  
    if (iblock .eq. 1 .and. jblock .eq. 1) then

        m = get_iter()/(tplot*Nmass_ave)
        mfile = trim(prefix_dir)//'results/mbins'

        call write_zplane(mass_out,1,mfile,m)

    end if

    ! Barrier needed to ensure correct offsets calculated next
    call MPI_Barrier(icomm_grid,ierr)

contains

    subroutine ZPlaneReduceMass
    use messenger_data_exchange, only : PlaneSum
    implicit none

        integer :: nbins
        integer, allocatable :: buf(:)

        nbins = cpol_binso(1) * cpol_binso(2) * cpol_binso(3)

        allocate( buf(nbins) )

        buf = reshape( mass_out, (/nbins/) )
        call PlaneSum(buf,nbins,3)
        mass_out = reshape(buf,(/cpol_binso(1), cpol_binso(2), cpol_binso(3)/))

        deallocate( buf )

    end subroutine ZPlaneReduceMass

end subroutine mass_bin_cpol_io

!------------------------------------------------------------------------------
! Cylindrical polar momentum bins
subroutine momentum_bin_cpol_io(mass_out,mom_out)
    use module_parallel_io, only : write_zplane, get_iter
    use concentric_cylinders, only : cpol_binso
    use physical_constants_MD, only : nd
    use computational_constants_MD, only : iter, initialstep, tplot, &
                                           Nvel_ave, prefix_dir
    use messenger, only: icomm_grid,iblock,jblock,plane_comm
    implicit none

    ! z-plane mass and momentum for all_reduce on zplane subcomms.
    ! The mass and mom arrays on each processor are global in
    ! r and theta, but local in z.  
    !integer, intent(inout) :: mass_out(cpol_binso(1), &
    !                                   cpol_binso(2), &
    !                                   cpol_binso(3))
    !real(kind(0.d0)), intent(inout) :: mom_out (cpol_binso(1), &
    !                                            cpol_binso(2), &
    !                                            cpol_binso(3), &
    !                                            nd)
    real(kind(0.d0)), dimension(:,:,:), intent(inout) :: mass_out
    real(kind(0.d0)), dimension(:,:,:,:), intent(inout) :: mom_out

    character(200) :: vfile
    integer :: m, ierr, nresults

    ! Make sure mass is written out
    call mass_bin_cpol_io(mass_out)

    ! Z-Plane global sums of r/theta bins, everyone has a copy of
    ! r/theta bins in their own z-plane
    call ZPlaneReduceMom

    ! Bottom corner z-line of processors write to file  
    if (iblock .eq. 1 .and. jblock .eq. 1) then

        nresults = 3
        vfile = trim(prefix_dir)//'results/vbins'
        m = get_iter()/(tplot*Nvel_ave)

        call write_zplane(mom_out,nresults,vfile,m)

    end if

    ! Barrier needed to ensure correct offsets calculated next
    call MPI_Barrier(icomm_grid,ierr)
    
contains

    subroutine ZPlaneReduceMom
    use messenger_data_exchange, only : PlaneSum
    implicit none

        integer :: nbins
        real(kind(0.d0)), allocatable :: buf(:)

        nbins = cpol_binso(1) * cpol_binso(2) * cpol_binso(3)

        allocate( buf(nd*nbins) )

        buf = reshape( mom_out, (/nd*nbins/) )
        call PlaneSum(buf, nd*nbins, 3 )
        mom_out = reshape(buf,(/cpol_binso(1),cpol_binso(2),cpol_binso(3),nd/))

        deallocate(buf)

    end subroutine ZPlaneReduceMom

end subroutine momentum_bin_cpol_io

!------------------------------------------------------------------------------
! Cylindrical polar KE bins
subroutine temperature_bin_cpol_io(mass_out,KE_out)
    use module_parallel_io, only: write_zplane, get_iter
    use concentric_cylinders, only: cpol_binso
    use physical_constants_MD, only: nd
    use computational_constants_MD, only: iter, initialstep, tplot, &
                                          NTemp_ave, prefix_dir, &
                                          momentum_outflag
    use messenger, only: icomm_grid,iblock,jblock,plane_comm
    implicit none

    ! z-plane mass and KE for all_reduce on zplane subcomms.
    ! The mass and KE arrays on each processor are global in
    ! r and theta, but local in z.  
    real(kind(0.d0)), dimension(:,:,:), intent(inout) :: mass_out
    real(kind(0.d0)), dimension(:,:,:), intent(inout) :: KE_out

    character(200) :: Tfile
    integer :: m, ierr

    ! Make sure mass is written out
    if (momentum_outflag .ne. 5) then
        call mass_bin_cpol_io(mass_out)
    end if

    ! Z-Plane global sums of r/theta bins, everyone has a copy of
    ! r/theta bins in their own z-plane
    call ZPlaneReduceKE

    ! Bottom corner z-line of processors write to file  
    if (iblock .eq. 1 .and. jblock .eq. 1) then

        m = get_iter()/(tplot*NTemp_ave)
        Tfile = trim(prefix_dir)//'results/Tbins'

        call write_zplane(KE_out,1,Tfile,m)

    end if

    ! Barrier needed to ensure correct offsets calculated next
    call MPI_Barrier(icomm_grid,ierr)
    
contains

    subroutine ZPlaneReduceKE
    use messenger_data_exchange, only : PlaneSum
    implicit none

        integer :: nbins
        real(kind(0.d0)), allocatable :: buf(:)

        nbins = cpol_binso(1) * cpol_binso(2) * cpol_binso(3)

        allocate( buf(nbins) )

        buf = reshape( KE_out, (/nbins/) )
        call PlaneSum(buf, nbins, 3 )
        KE_out = reshape(buf,(/cpol_binso(1),cpol_binso(2),cpol_binso(3)/))

        deallocate(buf)

    end subroutine ZPlaneReduceKE

end subroutine temperature_bin_cpol_io

!------------------------------------------------------------------------------
! Cylindrical polar momentum bins
subroutine VA_stress_cpol_io()
    use module_parallel_io, only: write_zplane, get_iter
    use concentric_cylinders, only: cpol_binso, r_oi, r_io, cpol_bins
    use computational_constants_MD, only: iter, initialstep, tplot, &
                                          Nstress_ave, prefix_dir, domain, &
                                          globaldomain, split_kin_config
    use physical_constants_MD, only: pi
    use messenger, only: icomm_grid,iblock,jblock,plane_comm
    use calculated_properties_MD, only: rfbin, vvbin, Pxybin
    implicit none

    integer :: m, ierr, nresults
    integer :: rbin, tbin, zbin 
    character(200) :: pVAfile
    real(kind(0.d0)), allocatable :: buf9(:,:,:,:)
    real(kind(0.d0)) :: rplus, rminus, dr, dtheta, dz
    real(kind(0.d0)) :: binvolume


    dr     = (r_io - r_oi) / cpol_bins(1)
    dtheta = (2.d0*pi)     / cpol_bins(2)
    dz     = domain(3)     / cpol_bins(3)
    
    !Bins have different volumes
    do rbin = 1, cpol_binso(1)
    do tbin = 1, cpol_binso(2)
    do zbin = 1, cpol_binso(3)

        rplus  = r_oi  +       rbin*dr
        rminus = r_oi  + (rbin - 1)*dr

        binvolume = (1.d0/2.d0)*(rplus**2.d0 - rminus**2.d0)*dtheta*dz

        vvbin(rbin,tbin,zbin,:,:) = vvbin(rbin,tbin,zbin,:,:)/binvolume
        rfbin(rbin,tbin,zbin,:,:) = rfbin(rbin,tbin,zbin,:,:)/binvolume

    end do
    end do
    end do

    ! Average over samples
    vvbin  = vvbin  / Nstress_ave 
    rfbin  = rfbin  / Nstress_ave 
    
    ! Double counted so divide by 2
    rfbin  = rfbin/2.d0

    ! Sum kinetic and configurational contributions
    Pxybin = vvbin + rfbin

    ! The stress arrays on each processor are global in r and theta, 
    ! but local in z. Z-Plane global sums of r/theta bins are performed
    ! so everyone has a copy of global r/theta bins in their own z-plane
    call ZPlaneReduceStress

    ! Only bottom corner z-line of processors write to file 
    if (iblock .eq. 1 .and. jblock .eq. 1) then

        allocate(buf9(cpol_binso(1),cpol_binso(2),cpol_binso(3),9))

        m = get_iter()/(tplot*Nstress_ave)
        nresults = 9
    
        select case (split_kin_config)
        case(0)
            !Together
            pVAfile = trim(prefix_dir)//'results/pVA'
            buf9 = reshape(Pxybin,(/cpol_binso(1),cpol_binso(2),cpol_binso(3),9/))
            call write_zplane(buf9,nresults,pVAfile,m)
        case(1)
            !Separate
            pVAfile = trim(prefix_dir)//'results/pVA_k'
            buf9 = reshape(vvbin,(/cpol_binso(1),cpol_binso(2),cpol_binso(3),9/))
            call write_zplane(buf9,nresults,pVAfile,m)
            pVAfile = trim(prefix_dir)//'results/pVA_c'
            buf9 = reshape(rfbin,(/cpol_binso(1),cpol_binso(2),cpol_binso(3),9/))
            call write_zplane(buf9,nresults,pVAfile,m)
        case default
        end select

        deallocate(buf9)

    end if

    ! Barrier needed to ensure correct offsets calculated next
    call MPI_Barrier(icomm_grid,ierr)
    
contains

    subroutine ZPlaneReduceStress()
    use messenger_data_exchange, only : PlaneSum
    implicit none

        integer :: nbins
        real(kind(0.d0)), allocatable :: buf(:)

        nbins = cpol_binso(1) * cpol_binso(2) * cpol_binso(3)

        select case (split_kin_config)
        case(0)

            !Together
            allocate( buf(9*nbins) )
            buf = reshape(Pxybin,(/9*nbins/))
            call PlaneSum( buf, 9*nbins, 3 )
            Pxybin = reshape(buf,(/cpol_binso(1),cpol_binso(2),cpol_binso(3),3,3/))
            deallocate(buf)

        case(1)

            !Separate
            allocate( buf(9*nbins) )

            buf = reshape(vvbin,(/9*nbins/))
            call PlaneSum( buf, 9*nbins, 3 )
            vvbin = reshape(buf,(/cpol_binso(1),cpol_binso(2),cpol_binso(3),3,3/))

            buf = reshape(rfbin,(/9*nbins/))
            call PlaneSum( buf, 9*nbins, 3 )
            rfbin = reshape(buf,(/cpol_binso(1),cpol_binso(2),cpol_binso(3),3,3/))

            deallocate(buf)

        case default
        end select

    end subroutine ZPlaneReduceStress

end subroutine VA_stress_cpol_io 

!---------------------------------------------------------------------------------
! Record temperature in a slice through the domain

subroutine temperature_slice_io(ixyz)
use module_parallel_io
implicit none

    integer     :: ixyz
    
    ! Remove the following lines when parallelising, they are only here
    ! to satisfy the compiler in pedantic mode
    integer :: dummy
    dummy = ixyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine temperature_slice_io

!---------------------------------------------------------------------------------
! Record temperature in 3D bins throughout domain

subroutine temperature_bin_io(CV_mass_out,CV_temperature_out,io_type)
    use module_parallel_io
    use messenger_bin_handler, only : swaphalos
    implicit none

    real(kind(0.d0)),intent(inout)  :: CV_mass_out(:,:,:)
    real(kind(0.d0)),intent(inout)  :: CV_temperature_out(:,:,:)

    integer             :: m,nresults
    character(4)        :: io_type
    character(30)       :: filename, outfile
    real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: CV_temperature_temp

    !Write mass bins
    if (momentum_outflag .ne. 4) then
        call mass_bin_io(CV_mass_out,io_type)
    endif

    !Work out correct filename for i/o type
    write(filename, '(a9,a4)' ) 'results/T', io_type
    outfile = trim(prefix_dir)//filename

    nresults = 1

    ! Swap Halos
    allocate(CV_temperature_temp(nbinso(1),nbinso(2),nbinso(3),nresults))
    CV_temperature_temp(:,:,:,1) = CV_temperature_out
    call swaphalos(CV_temperature_temp,nbinso(1),nbinso(2),nbinso(3),nresults)
    CV_temperature_out = CV_temperature_temp(:,:,:,1)
    deallocate(CV_temperature_temp)

    if (io_type .eq. 'snap') then
        !CV_temperature_out = CV_temperature_out / (tplot*Nvflux_ave)
        stop "Error in temp_io - Temperature SNAP called and NTflux_ave not defined"
        !m = iter/(NTflux_ave) + 1 !Initial snapshot taken
    else
        !CV_temperature_out = CV_temperature_out / (tplot*Nvel_ave)
        m = get_iter()/(tplot*NTemp_ave)
    endif
    !Write temperature to file
    call write_arrays(CV_temperature_out,nresults,outfile,m)

end subroutine temperature_bin_io

!---------------------------------------------------------------------------------
! Record energy in 3D bins throughout domain

subroutine energy_bin_io(CV_energy_out,io_type)
    use module_parallel_io
    use messenger_bin_handler, only : swaphalos
    use CV_objects, only : CVcheck_energy
    implicit none

    integer                         :: m,nresults
    real(kind(0.d0)),intent(inout)  :: CV_energy_out(:,:,:)

    real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: CV_energy_temp
    character(4)        :: io_type
    character(30)       :: filename, outfile

    !Work out correct filename for i/o type
    write(filename, '(a9,a4)' ) 'results/e', io_type
    outfile = trim(prefix_dir)//filename

    !---------------Correct for surface fluxes on halo cells---------------
    ! Swap Halos
    nresults = 1
    allocate(CV_energy_temp(nbinso(1),nbinso(2),nbinso(3),nresults))
    CV_energy_temp(:,:,:,1) = CV_energy_out
    call swaphalos(CV_energy_temp,nbinso(1),nbinso(2),nbinso(3),nresults)
    CV_energy_out = CV_energy_temp(:,:,:,1)
    deallocate(CV_energy_temp)

    if (io_type .eq. 'snap') then
        !CV_energy_out = CV_energy_out / (tplot*Nvflux_ave)
        select case(CV_conserve)
        case(0)
            m = get_iter()/(tplot*Neflux_ave) + 1 !Initial snapshot taken
        case(1)
            m = get_iter()/(Neflux_ave) + 1 !Initial snapshot taken
            !Create copy of previous timestep Control Volume energy to calculate time evolution
            if (CV_debug .ne. 0) then
                call CVcheck_energy%update_dXdt(CV_energy_out(:,:,:))
            endif
        case default
            call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')   
        end select
    elseif (io_type .eq. 'bins') then
        m = get_iter()/(tplot * Nenergy_ave)
    else
        stop "Energy io error - unknown output type"
    endif

    !Write Energy to file
    call write_arrays(CV_energy_out,nresults,outfile,m)

end subroutine energy_bin_io

!---------------------------------------------------------------------------------

subroutine virial_stress_io()
    use module_parallel_io
    implicit none

    integer     :: m, length

    !call globalAverage(Pxy, 9)

    !Write virial pressure to file
    m = real(get_iter())/real(tplot*Nstress_ave)
    if (m .le. 0) then
        print*, 'WARNING: negative record number in virial_stress_io'
        return
    end if
    
    if (irank .eq. iroot) then
        inquire(iolength=length) Pxy
        open (unit=7, file=trim(prefix_dir)//'results/pvirial',form='unformatted',access='direct',recl=length)
        write(7,rec=m) Pxy
        close(7,status='keep')
    endif

end subroutine virial_stress_io

!---------------------------------------------------------------------------------

subroutine VA_stress_io()
    use module_parallel_io
    implicit none

    integer                                         :: ixyz, jxyz, m, nresults
    real(kind(0.d0))                                :: binvolume
    real(kind(0.d0)),dimension(:,:,:,:),allocatable :: buf

    !Add kinetic and configurational to Pxybin total
    Pxybin(:,:,:,:,:) =     vvbin(:,:,:,:,:)        & 
                          + rfbin(1+nhb(1):nbins(1)+nhb(1),   & 
                                  1+nhb(2):nbins(2)+nhb(2),   & 
                                  1+nhb(3):nbins(3)+nhb(3),:,:)/2.d0

    !Average over samples
    Pxybin = Pxybin / Nstress_ave
    vvbin  = vvbin  / Nstress_ave 
    rfbin  = rfbin  / Nstress_ave 

    !Calculate Virial from sum of Volume Averaged Stresses
    do ixyz = 1,3
    do jxyz = 1,3
        Pxy(ixyz,jxyz) = sum(Pxybin(:,:,:,ixyz,jxyz))
    enddo
    enddo
    Pxy = Pxy / volume

    !Write out Virial stress
    call virial_stress_io

!       vvbin = 0.d0
    !VA pressure per bin
    binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
    Pxybin = Pxybin / binvolume
    vvbin  = vvbin  / binvolume
    rfbin  = rfbin  / (2.d0*binvolume)


!    print'(a,i8,3f20.10)', 'Pressure VA', iter, &
!                            sum(Pxybin(:,:,:,1,1)+Pxybin(:,:,:,2,2)+Pxybin(:,:,:,3,3)), &
!                            sum(vvbin(:,:,:,1,1)+vvbin(:,:,:,2,2)+vvbin(:,:,:,3,3)), &
!                            sum(rfbin(:,:,:,1,1)+rfbin(:,:,:,2,2)+rfbin(:,:,:,3,3))

    !Write VA pressure to file
    nresults = 9
    m = get_iter()/(tplot*Nstress_ave)

    select case (split_kin_config)
    case(0)
        !print'(a,3i8,9f10.5)','Pxybin', iter,irank,size(Pxybin), & 
        !                                           Pxybin(5,5,5,1,1),Pxybin(5,5,5,2,1),Pxybin(5,5,5,3,1), &
        !                                           Pxybin(5,5,5,1,2),Pxybin(5,5,5,2,2),Pxybin(5,5,5,3,2), &
        !                                           Pxybin(5,5,5,1,3),Pxybin(5,5,5,2,3),Pxybin(5,5,5,3,3)

        !Write sum of kinetic and configurational
        !Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
        allocate(buf(nbinso(1), nbinso(2), nbinso(3), nresults))
        buf = 0.d0; 
        buf(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),1:nresults) = &
            reshape(Pxybin,(/nbins(1),nbins(2),nbins(3),nresults/))
        call write_arrays(buf,nresults,trim(prefix_dir)//'results/pVA',m)
        deallocate(buf)
    case(1)
        !Kinetic
        !Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
        allocate(buf(nbinso(1), nbinso(2), nbinso(3), nresults))
        buf = 0.d0 
        buf(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),1:nresults) = &
            reshape(vvbin,(/nbins(1),nbins(2),nbins(3),nresults/))
        call write_arrays(buf,nresults,trim(prefix_dir)//'results/pVA_k',m)
        deallocate(buf)
        !Configurational
        !Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
        allocate(buf(nbinso(1), nbinso(2), nbinso(3), nresults))
        buf = 0.d0; 
        buf(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),1:nresults) = &
            reshape(rfbin(1+nhb(1):nbins(1)+nhb(1), & 
                          1+nhb(2):nbins(2)+nhb(2), & 
                          1+nhb(3):nbins(3)+nhb(3),:,:), & 
                    (/nbins(1),nbins(2),nbins(3),nresults/))
            !reshape(rfbin,(/nbins(1),nbins(2),nbins(3),nresults/))
        call write_arrays(buf,nresults,trim(prefix_dir)//'results/pVA_c',m)
        deallocate(buf)
    case default
        stop 'Error in VA/virial extra flag to split_kinetic_& configuartional parts'
    end select


end subroutine VA_stress_io





!---------------------------------------------------------------------------------

subroutine VA_heatflux_io()
    use module_parallel_io
    implicit none

    integer                                         :: ixyz, jxyz, m, nresults
    real(kind(0.d0))                                :: binvolume
    real(kind(0.d0)), dimension(3)                  :: heatflux
    real(kind(0.d0)),dimension(:,:,:,:),allocatable :: buf

    !Add kinetic and configurational to heatfluxbin total
    heatfluxbin(:,:,:,:) =     evbin(:,:,:,:)        & 
                            + rfvbin(1+nhb(1):nbins(1)+nhb(1),   & 
                                     1+nhb(2):nbins(2)+nhb(2),   & 
                                     1+nhb(3):nbins(3)+nhb(3),:,1)/2.d0

    !Average over samples
    heatfluxbin = heatfluxbin / Nstress_ave
    evbin  = evbin  / Nstress_ave 
    rfvbin  = rfvbin  / Nstress_ave 

	!Calculate total from sum of Volume Averaged heatfluxes
	do ixyz = 1,3
		heatflux(ixyz) = sum(heatfluxbin(:,:,:,ixyz))
	enddo

    !Write out total heatflux
    heatflux = heatflux / volume
    !call total_heatflux_io(heatflux)

    binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
    heatfluxbin = heatfluxbin / binvolume
    evbin  = evbin  / binvolume
    rfvbin = rfvbin  / (2.d0*binvolume)

    !Write VA pressure to file
    nresults = 3
    m = get_iter()/(tplot*Nstress_ave)

    select case (split_hfkin_config)
    case(0)

        !Write sum of kinetic and configurational
        !Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
        allocate(buf(nbinso(1), nbinso(2), nbinso(3), nresults))
        buf = 0.d0; 
        buf(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),1:nresults) = &
            reshape(heatfluxbin,(/nbins(1),nbins(2),nbins(3), nresults/))
        call write_arrays(buf,nresults,trim(prefix_dir)//'results/hfVA',m)
        deallocate(buf)
    case(1)
        !Kinetic
        !Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
        allocate(buf(nbinso(1), nbinso(2), nbinso(3), nresults))
        buf = 0.d0; 
        buf(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),1:nresults) = &
            reshape(evbin,(/nbins(1),nbins(2),nbins(3),nresults/))
        call write_arrays(buf,nresults,trim(prefix_dir)//'results/hfVA_k',m)
        deallocate(buf)
        !Configurational
        !Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
        allocate(buf(nbinso(1), nbinso(2), nbinso(3), nresults))
        buf = 0.d0; 
        buf(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),1:nresults) = &
            reshape(rfvbin(1+nhb(1):nbins(1)+nhb(1), & 
                           1+nhb(2):nbins(2)+nhb(2), & 
                           1+nhb(3):nbins(3)+nhb(3),:,:), & 
                    (/nbins(1),nbins(2),nbins(3),nresults/))
            !reshape(rfbin,(/nbins(1),nbins(2),nbins(3),nresults/))
        call write_arrays(buf,nresults,trim(prefix_dir)//'results/hfVA_c',m)
        deallocate(buf)
    case default
        stop 'Error in Heatflux extra flag to split_kinetic_& configuartional parts'
    end select


end subroutine VA_heatflux_io

subroutine total_heatflux_io(heatflux)
    use module_parallel_io
	use librarymod, only : get_new_fileunit
    implicit none
   
    integer :: fileunit
    real(kind(0.d0)),dimension(:),intent(in)    :: heatflux
    integer     :: m, length

    !Write total_heatflux to file
    m = real(get_iter())/real(tplot*Nheatflux_ave)
    if (m .le. 0) then
        print*, 'WARNING: negative record number in total_heatflux_io'
        return
    end if
    
    if (irank .eq. iroot) then
        inquire(iolength=length) heatflux
        fileunit = get_new_fileunit()
        open (unit=fileunit, file=trim(prefix_dir)//'results/totalheatflux',form='unformatted',access='direct',recl=length)
        write(fileunit,rec=m) heatflux
        close(fileunit,status='keep')
    endif

end subroutine total_heatflux_io






end module field_io

module statistics_io

contains

#if __INTEL_COMPILER > 1200
    subroutine bforce_pdf_write()
        use module_parallel_io, only: get_iter
        use boundary_MD
        use computational_constants_MD, only: separate_outfiles, irank, iroot,&
                                              iter, initialstep, tplot, &
                                              prefix_dir
        use messenger_data_exchange, only : globalSum
        use librarymod, only: get_new_fileunit, get_Timestep_FileName
        implicit none

        integer :: n, ixyz
        integer :: funit, record, length

        integer, dimension(:,:,:), allocatable :: array_out
        integer, dimension(:), allocatable :: buf
        character(30) :: outfile, outfile_t

        allocate(array_out(bforce_pdf_nsubcells,bforce_pdf_nbins,3))

        do n = 1, bforce_pdf_nsubcells  
            do ixyz = 1,3
                array_out(n,:,ixyz) = bforce_pdf(ixyz,n)%hist
            end do
        end do

        outfile = 'bforce_pdf'
        record = get_iter()/(tplot*bforce_pdf_Nave)

        if (separate_outfiles) then
            call get_Timestep_FileName(record-1,outfile,outfile_t)
            record = 1
        else
            outfile_t = outfile 
        endif

        ! Global reduce
        allocate(buf(bforce_pdf_nsubcells*bforce_pdf_nbins*3))
        buf = reshape(array_out, (/bforce_pdf_nsubcells*bforce_pdf_nbins*3/))
        call globalSum(buf,size(buf,1)) 
        array_out = reshape(buf,(/bforce_pdf_nsubcells,bforce_pdf_nbins,3/))

        ! PDFs are a global quantity, so only the root processor writes
        if (irank.eq.iroot) then
        
            funit = get_new_fileunit()
            inquire(iolength=length) array_out
            open(unit=funit,file=trim(prefix_dir)//'results/'//outfile_t,&
                 form='unformatted',access='direct',recl=length)
            write(funit,rec=record) array_out
            close(funit,status='keep') 

        end if

        deallocate(array_out)

    end subroutine bforce_pdf_write

#endif

end module statistics_io

!=================================================================================
! CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV
! CV            Record Fluxes accross surfaces of Control Volumes               CV
! CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV
!=================================================================================

!---------------------------------------------------------------------------------
! Record mass fluxes accross surfaces of Control Volumes

subroutine mass_flux_io()
    use module_parallel_io
    use CV_objects, only : CVcheck_mass, CV_debug
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer             :: m,nresults
    character(30)       :: filename, outfile

    !Work out correct filename for i/o type
    filename='results/mflux'
    outfile = trim(prefix_dir)//filename

    !Include halo surface fluxes to get correct values for all cells
    nresults = 6
    call swaphalos(mass_flux,nbinso(1),nbinso(2),nbinso(3),nresults)

    !Store mass flux value in CV data object
    if (CV_debug .ne. 0) then
        !CVcheck_mass%flux = mass_flux
		call CVcheck_mass%update_flux(mass_flux)
    endif

    !Calculate record number timestep
    select case(CV_conserve)
    case(0)
        m = get_iter()/(tplot*Nmflux_ave)
    case(1)
        m = get_iter()/(Nmflux_ave)
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    !Write mass to file
    call write_arrays(mass_flux,nresults,outfile,m)

end subroutine mass_flux_io

!---------------------------------------------------------------------------------
! Record momentum fluxes accross surfaces of Control Volumes

subroutine momentum_flux_io()
    use module_parallel_io
    use CV_objects, only : CVcheck_momentum, CV_debug
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer                                         :: ixyz,m,nresults
    real(kind(0.d0))                                :: binface
    real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: momentum_flux_temp

    ! Swap Halos
    nresults = 18
    allocate(momentum_flux_temp(size(momentum_flux,1),size(momentum_flux,2),size(momentum_flux,3),nresults))
    momentum_flux_temp = reshape(momentum_flux, & 
                        (/ size(momentum_flux,1), & 
                           size(momentum_flux,2), & 
                           size(momentum_flux,3),nresults /))
    call swaphalos(momentum_flux_temp,nbinso(1),nbinso(2),nbinso(3),nresults)
    momentum_flux = reshape(momentum_flux_temp,(/ size(momentum_flux,1), & 
                                                  size(momentum_flux,2), & 
                                                  size(momentum_flux,3),3,6 /))
    !deallocate(momentum_flux_temp)

    !Divide by size of bin face to give flux per unit area
    do ixyz = 1,3
        binface       = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
                        (domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
        momentum_flux(:,:,:,:,ixyz  )=momentum_flux(:,:,:,:,ixyz  )/(binface) !Bottom
        momentum_flux(:,:,:,:,ixyz+3)=momentum_flux(:,:,:,:,ixyz+3)/(binface) !Top
    enddo

    !Divide momentum flux by averaing period tau=delta_t*Nvflux_ave if CV_conserve=1
    !or Divide momentum flux by sum of the Nvflux_ave times delta_t averaging periods 
    !as sample is taken every tplot steps. The output is then a representaive momentum flux.
    momentum_flux = momentum_flux/(delta_t*Nvflux_ave)

    !Store momentum flux value in CV data object
    if (CV_debug .ne. 0) then
        CVcheck_momentum%flux = 0.d0
        call CVcheck_momentum%update_flux(momentum_flux)
    endif

    !Write momentum to file
    select case(CV_conserve)
    case(0)
        m = get_iter()/(tplot*Nvflux_ave)
    case(1)
        m = get_iter()/(Nvflux_ave)
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    !allocate(momentum_flux_temp(size(momentum_flux,1),size(momentum_flux,2),size(momentum_flux,3),nresults))
    momentum_flux_temp = reshape(momentum_flux,(/ size(momentum_flux,1), & 
												  size(momentum_flux,2), & 
											      size(momentum_flux,3),nresults /))
    call write_arrays(momentum_flux_temp,nresults,trim(prefix_dir)//'results/vflux',m)
    deallocate(momentum_flux_temp)

end subroutine momentum_flux_io

!---------------------------------------------------------------------------------
! Record stress accross plane

subroutine MOP_stress_io(ixyz_in)
    use module_parallel_io
    use messenger
    use field_io, only : mass_slice_io
    implicit none

    integer, intent(in) :: ixyz_in

    integer                         :: m, ixyz, jxyz, kxyz
    integer                         :: slicefileid, dp_datasize
    integer,dimension(3)            :: idims
    integer(kind=MPI_OFFSET_KIND)   :: disp

    !This line added to satisfy pedantic compiler. Remove if you want.
    ixyz = ixyz_in

    ixyz = vflux_outflag

    !Divide by number of samples taken
    Pxy_plane = Pxy_plane/(Nstress_ave)

    !Divide by area of domain and factor of 4 for interactions
    Pxy_plane = Pxy_plane/(4*domain(1)*domain(3))

    !Write plane pressures to file
    m = get_iter()/(tplot*Nvflux_ave)
!************WRITTEN ROUGHLY AND NOT TESTED************************
    !Write mass
    call mass_slice_io(ixyz)

    !Get two directions orthogonal to slice direction
    kxyz = mod(ixyz,3)+1
    jxyz = mod(ixyz+1,3)+1
    idims(1) = npx; idims(2) = npy; idims(3) = npz

    !Sum over all bins using directional sub communicators and gather on root
    call SubcommSum(Pxy_plane(1,:), nplanes, jxyz)
    call SubcommSum(Pxy_plane(1,:), nplanes, kxyz)
    call SubcommSum(Pxy_plane(2,:), nplanes, jxyz)
    call SubcommSum(Pxy_plane(2,:), nplanes, kxyz)
    call SubcommSum(Pxy_plane(3,:), nplanes, jxyz)
    call SubcommSum(Pxy_plane(3,:), nplanes, kxyz)

    !Only root processor in each directional subcomm writes data
    if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

        call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

        !Only processors on directional subcomm write
        call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'/results/pplane', & 
                   MPI_MODE_WRONLY + MPI_MODE_CREATE , & 
                   MPI_INFO_NULL, slicefileid, ierr)

        select case(vflux_outflag)
            case(1)
                !Obtain displacement of x record
                disp =   (get_iter()/(tplot*Nmass_ave) - 1)       &   !Current iteration
                       * nd*gnplanes*dp_datasize            &   !times record size
                       + nplanes*dp_datasize*(iblock-1)         !Processor location
            case(2)
                    !Obtain displacement of x record
                disp =   (get_iter()/(tplot*Nmass_ave) - 1)       &   !Current iteration
                       * nd*gnplanes*dp_datasize            &   !times record size
                       + nplanes*dp_datasize*(jblock-1)         !Processor location
            case(3)
                !Obtain displacement of x record
                disp =   (get_iter()/(tplot*Nmass_ave) - 1)       &   !Current iteration
                       * nd*gnplanes*dp_datasize            &   !times record size
                       + nplanes*dp_datasize*(kblock-1)         !Processor location
            case default
                call error_abort('MOP average output called with incorrect vflux_outflag')
        end select

        call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
                    MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(slicefileid,Pxy_plane,nplanes, & 
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

        call MPI_FILE_CLOSE(slicefileid, ierr)

    endif 
!************WRITTEN ROUGHLY AND NOT TESTED************************

end subroutine MOP_stress_io

!---------------------------------------------------------------------------------
! Record stress accross surfaces of Control Volumes

subroutine surface_stress_io()
    use module_parallel_io
    use CV_objects, only : CVcheck_momentum,CV_debug
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer                                             :: ixyz,m,nresults, i,j,k
    real(kind(0.d0))                                    :: binface
    real(kind(0.d0)),allocatable,dimension(:,:,:,:)     :: Pxyface_temp

    ! Swap Halos
    nresults = 18
    allocate(Pxyface_temp(size(Pxyface,1),size(Pxyface,2),size(Pxyface,3),nresults))
    Pxyface_temp = reshape(Pxyface,(/ size(Pxyface,1),size(Pxyface,2),size(Pxyface,3),nresults /))
    call swaphalos(Pxyface_temp,nbinso(1),nbinso(2),nbinso(3),nresults)
    Pxyface = reshape(Pxyface_temp,(/ size(Pxyface,1),size(Pxyface,2),size(Pxyface,3),3,6 /))
    !deallocate(Pxyface_temp)

    !Divide by size of bin face to give flux per unit area
    do ixyz = 1,3
        binface       = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
                        (domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
        Pxyface(:,:,:,:,ixyz  ) = 0.25d0 * Pxyface(:,:,:,:,ixyz  )/binface !Bottom
        Pxyface(:,:,:,:,ixyz+3) = 0.25d0 * Pxyface(:,:,:,:,ixyz+3)/binface !Top
    enddo

    !Integration of stress using trapizium rule requires multiplication by timestep
    !so delta_t cancels upon division by tau=delta_t*Nvflux_ave resulting in division by Nvflux_ave
    Pxyface = Pxyface/Nvflux_ave

    !Store surface stress value in CV data object
    if (CV_debug .ne. 0) then
        call CVcheck_momentum%update_Pxy(Pxyface)

		! do i =136,137
		! do j =18,19
		! do k =18,19
			! print'(a,i8,3i5, 6f11.4)', "Pxyface ", iter, i, j, k, &
						! Pxyface(i,j,k,1,:)
		! enddo
		! enddo
		! enddo

    endif

    !Write surface pressures to file
    select case(CV_conserve)
    case(0)
        m = get_iter()/(tplot*Nvflux_ave)
    case(1)
        m = get_iter()/(Nvflux_ave)
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    if (m .eq. 0) return

    !Write surface pressures to file
    Pxyface_temp = reshape(Pxyface,(/ size(Pxyface,1),size(Pxyface,2),size(Pxyface,3),nresults /))
    call write_arrays(Pxyface_temp,nresults,trim(prefix_dir)//'results/psurface',m)
    deallocate(Pxyface_temp)

end subroutine surface_stress_io

!---------------------------------------------------------------------------------
! Record external forces applied to molecules inside a volume

subroutine external_force_io()
    use module_parallel_io
    use CV_objects, only : CVcheck_momentum,CV_debug
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer                                             :: m,nresults
    real(kind(0.d0)),dimension(:,:,:,:),allocatable     :: temp

    ! Swap Halos
    nresults = 3
    call swaphalos(F_ext_bin,nbinso(1),nbinso(2),nbinso(3),nresults)

    !Integration of force using trapizium rule requires multiplication by timestep
    !so delta_t cancels upon division by tau=delta_t*Nvflux_ave resulting in division by Nvflux_ave
    F_ext_bin = F_ext_bin/Nvflux_ave

    !Store surface stress value in CV data object
    if (CV_debug .ne. 0) then
        call CVcheck_momentum%update_F_ext(F_ext_bin)
        !CVcheck_momentum%F_ext = F_ext_bin
    endif

    !Write external forces pressures to file
    select case(CV_conserve)
    case(0)
        m = get_iter()/(Nvflux_ave*tplot)
    case(1)
        m = get_iter()/(Nvflux_ave)
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    !Write external forces to file
    !allocate(temp(size(F_ext_bin,1),size(F_ext_bin,2),size(F_ext_bin,3),nresults))
    !temp = reshape(F_ext_bin,(/ size(F_ext_bin,1),size(F_ext_bin,2),size(F_ext_bin,3),nresults /))
    call write_arrays(F_ext_bin,nresults,trim(prefix_dir)//'results/Fext',m)
    !deallocate(temp)

end subroutine external_force_io


!---------------------------------------------------------------------------------
! Record energy fluxes accross surfaces of Control Volumes

subroutine energy_flux_io()
    use module_parallel_io
    use messenger_bin_handler, only : swaphalos
    use CV_objects, only : CVcheck_energy
    implicit none

    integer                 :: ixyz,m,nresults
    real(kind(0.d0))        :: binface

    !Include halo surface fluxes to get correct values for all cells
    nresults = 6
    call swaphalos(energy_flux,nbinso(1),nbinso(2),nbinso(3),nresults)

    do ixyz = 1,3
        binface       = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
                        (domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
        energy_flux(:,:,:,ixyz  )=energy_flux(:,:,:,ixyz  )/binface !Bottom
        energy_flux(:,:,:,ixyz+3)=energy_flux(:,:,:,ixyz+3)/binface !Top
    enddo

    !Divide energy flux by averaing period tau=delta_t*Neflux_ave if CV_conserve=1
    !or Divide sample energy flux by equivalent averaging period delta_t*Neflux_ave
    energy_flux = energy_flux/(delta_t*Neflux_ave)

    !Store momentum flux value in CV data object
    if (CV_debug .ne. 0) then
        CVcheck_energy%flux = 0.d0
        call CVcheck_energy%update_flux(energy_flux)
    endif

    !Write energy flux to file
    select case(CV_conserve)
    case(0)
        m = get_iter()/(Neflux_ave*tplot)
    case(1)
        m = get_iter()/(Neflux_ave)
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    call write_arrays(energy_flux,nresults,trim(prefix_dir)//'results/eflux',m)

end subroutine energy_flux_io

!-----------------------------------
! Density at the surface from 1/velocity at the surface

subroutine surface_density_io()
    use module_parallel_io
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer             :: m,nresults
    character(30)       :: filename, outfile

    !Work out correct filename for i/o type
    filename='results/msurf'
    outfile = trim(prefix_dir)//filename

    !Include halo surface fluxes to get correct values for all cells
    nresults = 6
    call swaphalos(surface_density,nbinso(1),nbinso(2),nbinso(3),nresults)

    !Calculate record number timestep
    select case(CV_conserve)
    case(0)
        m = get_iter()/(tplot*Nsurfm_ave)
    case(1)
        m = get_iter()/(Nsurfm_ave)
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    !Write mass to file
    call write_arrays(surface_density, nresults, outfile, m)

end subroutine surface_density_io

!---------------------------------------------------------------------------------
! Record stress times velocity (power) accross surfaces of Control Volumes

subroutine surface_power_io()
    use module_parallel_io
    use messenger_bin_handler, only : swaphalos
    use CV_objects, only : CVcheck_energy
    implicit none

    integer                 :: ixyz,m,nresults
    real(kind(0.d0))        :: binface

    !Include halo surface stresses to get correct values for all cells
    nresults = 6
    call swaphalos(Pxyvface_integrated,nbinso(1),nbinso(2),nbinso(3),nresults)

    do ixyz = 1,3
        binface = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
                  (domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
        Pxyvface_integrated(:,:,:,ixyz  ) = & 
            0.25d0*Pxyvface_integrated(:,:,:,ixyz  )/binface !Bottom
        Pxyvface_integrated(:,:,:,ixyz+3) = & 
            0.25d0*Pxyvface_integrated(:,:,:,ixyz+3)/binface !Top
    enddo

    !Integration of stress using trapizium rule requires division by number of results
    Pxyvface_integrated = Pxyvface_integrated/Neflux_ave

    !Store surface stress value in CV data object
    if (CV_debug .ne. 0) then
        call CVcheck_energy%update_Pxy(Pxyvface_integrated)
    endif

    !Write surface pressures * velocity to file
    select case(CV_conserve)
    case(0)
        m = get_iter()/(Neflux_ave*tplot)
    case(1)
        m = get_iter()/(Neflux_ave)
    case default
        call error_abort('CV_conserve value used for surface power is incorrectly defined - should be 0=off or 1=on')
    end select

    call write_arrays(Pxyvface_integrated, nresults, & 
                      trim(prefix_dir)//'results/esurface', m)

end subroutine surface_power_io


!---------------------------------------------------------------------------------
! Record external forces times velocity applied to molecules inside a volume

subroutine external_forcev_io()
    use module_parallel_io
    use CV_objects, only : CVcheck_energy
    use messenger_bin_handler, only : swaphalos
    implicit none

    integer                                  :: m,nresults
    real(kind(0.d0)), dimension(:,:,:,:), allocatable :: Fv_ext_copy

    !Copy CV_energy_out so it is not changed
    nresults = 1
    allocate(Fv_ext_copy(nbinso(1),nbinso(2),nbinso(3),nresults))
    Fv_ext_copy(:,:,:,1) = Fv_ext_bin(:,:,:)
    call swaphalos(Fv_ext_copy,nbinso(1),nbinso(2),nbinso(3),nresults)
    Fv_ext_bin = Fv_ext_copy(:,:,:,1)
    deallocate(Fv_ext_copy)

    !Integration of force using trapizium rule requires multiplication by timestep
    !so delta_t cancels upon division by tau=delta_t*Nvflux_ave resulting in division by Nvflux_ave
    Fv_ext_bin = Fv_ext_bin/Neflux_ave

    !Store surface stress value in CV data object
    if (CV_debug .ne. 0) then
        call CVcheck_energy%update_F_ext(Fv_ext_bin)
        !CVcheck_energy%F_ext = Fv_ext_bin
    endif

    !Write external forces pressures to file
    select case(CV_conserve)
    case(0)
        m = get_iter()/(Neflux_ave*tplot)
    case(1)
        m = get_iter()/(Neflux_ave)
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    !Write external forces to file
    call write_arrays(Fv_ext_bin,nresults,trim(prefix_dir)//'results/Fvext',m)

end subroutine external_forcev_io



subroutine surface_evolution_mass_flux_io()
    use module_parallel_io
    use CV_objects, only : CVcheck_mass, CV_debug
    use messenger_bin_handler, only : swaphalos
    use computational_constants_MD, only : II_normal
    implicit none

    integer             :: m,nresults
    character(41)       :: filename, outfile
	real(kind(0.d0)),dimension(:,:,:,:), allocatable	:: surface_flux

    !Work out correct filename for i/o type
    filename='results/dsurf_mflux'
    outfile = trim(prefix_dir)//filename

    !Include halo surface fluxes to get correct values for all cells
    nresults = 6
    call swaphalos(mass_surface_flux,nbinso(1),nbinso(2),nbinso(3),nresults)

    !Store mass flux value in CV data object
    if (CV_debug .ne. 0) then
        !Surface evo stores other components in top and bottom so
        !will ruin conservation, extract just normal for tests
        if (Nsurfevo_outflag .eq. 2) then
            allocate(surface_flux(size(mass_flux,1), size(mass_flux,2), & 
                                  size(mass_flux,3), size(mass_flux,4)))
            surface_flux = 0.d0
            surface_flux(:,:,:,II_normal) = mass_surface_flux(:,:,:,II_normal)
            surface_flux(:,:,:,II_normal+3) = mass_surface_flux(:,:,:,II_normal+3)
            call CVcheck_mass%update_surface(surface_flux)
            deallocate(surface_flux)
        else
            call CVcheck_mass%update_surface(mass_surface_flux)
        endif
    endif

    !Calculate record number timestep
    select case(CV_conserve)
    case(0)
        m = get_iter()/(tplot*Nmflux_ave)+1
    case(1)
        m = get_iter()/(Nmflux_ave)+1
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    !Write mass to file
    call write_arrays(mass_surface_flux, nresults, outfile, m)

end subroutine surface_evolution_mass_flux_io

subroutine surface_evolution_momentum_flux_io()
    use module_parallel_io
    use CV_objects, only : CVcheck_momentum, CV_debug
    use messenger_bin_handler, only : swaphalos
    use computational_constants_MD, only : II_normal
    implicit none

    integer             :: m,nresults, ixyz
    character(41)       :: filename, outfile
    real(kind(0.d0))                                :: binface
    real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: temp
	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable	:: surface_flux

    !Work out correct filename for i/o type
    filename='results/dsurf_vflux'
    outfile = trim(prefix_dir)//filename

    ! Swap Halos
    nresults = 18
    allocate(temp(size(momentum_surface_flux,1),size(momentum_surface_flux,2),size(momentum_surface_flux,3),nresults))
    temp = reshape(momentum_surface_flux, & 
                        (/ size(momentum_surface_flux,1), & 
                           size(momentum_surface_flux,2), & 
                           size(momentum_surface_flux,3),nresults /))
    call swaphalos(temp,nbinso(1),nbinso(2),nbinso(3),nresults)
    momentum_surface_flux = reshape(temp,(/ size(momentum_flux,1), & 
                                            size(momentum_flux,2), & 
                                            size(momentum_flux,3),3,6 /))
											
											
	!Divide by size of bin face to give flux per unit area
    do ixyz = 1,3
        binface       = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
                        (domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
        momentum_surface_flux(:,:,:,:,ixyz  )=momentum_surface_flux(:,:,:,:,ixyz  )/(binface) !Bottom
        momentum_surface_flux(:,:,:,:,ixyz+3)=momentum_surface_flux(:,:,:,:,ixyz+3)/(binface) !Top
    enddo

    !Divide momentum flux by averaing period tau=delta_t*Nvflux_ave if CV_conserve=1
    !or Divide momentum flux by sum of the Nvflux_ave times delta_t averaging periods 
    !as sample is taken every tplot steps. The output is then a representaive momentum flux.
    momentum_surface_flux = momentum_surface_flux/(delta_t*Nvflux_ave)


    !Store momentum flux value in CV data object
    if (CV_debug .ne. 0) then
        !Surface evo stores other components in top and bottom so
        !will ruin conservation, extract just normal for tests
        if (Nsurfevo_outflag .eq. 2) then
            allocate(surface_flux(size(momentum_surface_flux,1), size(momentum_surface_flux,2), & 
                                  size(momentum_surface_flux,3), size(momentum_surface_flux,4), &
                                  size(momentum_surface_flux,5)))
            surface_flux = 0.d0
            surface_flux(:,:,:,:,II_normal) = momentum_surface_flux(:,:,:,:,II_normal)
            surface_flux(:,:,:,:,II_normal+3) = momentum_surface_flux(:,:,:,:,II_normal+3)
            call CVcheck_momentum%update_surface(surface_flux)
            deallocate(surface_flux)
        else
            call CVcheck_momentum%update_surface(momentum_surface_flux)
        endif
    endif

    !Calculate record number timestep
    select case(CV_conserve)
    case(0)
        m = get_iter()/(tplot*Nvflux_ave)+1
    case(1)
        m = get_iter()/(Nvflux_ave)+1
    case default
        call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
    end select

    !Write mass to file
    temp = reshape(momentum_surface_flux,(/ size(momentum_surface_flux,1), & 
											size(momentum_surface_flux,2), & 
											size(momentum_surface_flux,3),nresults /))
    call write_arrays(temp,nresults,trim(prefix_dir)//'results/dsurf_vflux',m)
    deallocate(temp)

end subroutine surface_evolution_momentum_flux_io


!---------------------------------------------------------------------------------
! Record  energy accross plane

subroutine MOP_energy_io(ixyz)
    use module_parallel_io
    implicit none

    integer     :: ixyz

    ! Remove when parallelising, only here to satisfy pedantic compiler
    integer :: dummy
    dummy = ixyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine MOP_energy_io

!-----------------------------------------------------------------------------
! Write macroscopic properties to file
!-----------------------------------------------------------------------------
subroutine macroscopic_properties_header()
    use module_parallel_io
	use librarymod, only : get_new_fileunit
    implicit none

    integer :: fileunit

    if (irank .eq. iroot) then  
	    fileunit = get_new_fileunit()

        open(unit=fileunit,file=trim(prefix_dir)//'results/macroscopic_properties',status='replace')
        if (potential_flag .eq. 0) then
            write(fileunit,'(2a)') &
            ' iter; simtime; VSum; V^2Sum; Temp;', &
            ' KE; PE; TE; Pressure'
        else if (potential_flag.eq.1) then
            write(fileunit,'(2a)') &
            ' iter; simtime; VSum; V^2Sum; Temp;', &
            ' KE; PE (LJ); PE (POLY); PE (Tot); TE; Pressure; Etevtcf; R_g '
        end if
        close(unit=fileunit, status='keep')
    endif
    call macroscopic_properties_record

end subroutine macroscopic_properties_header


subroutine macroscopic_properties_record()
    use module_parallel_io
	use librarymod, only : get_new_fileunit
    implicit none
   
    integer :: fileunit

    if (irank .eq. iroot) then
	    fileunit = get_new_fileunit()
        open(unit=fileunit,file=trim(prefix_dir)//'results/macroscopic_properties',position='append')
        if (potential_flag.eq.0) then   
            write(fileunit,'(1x,i8,a,f15.4,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)') &
            iter,';',simtime,';',vsum,';', mv2sum,';', temperature,';', &
            kinenergy,';',potenergy,';',totenergy,';',pressure
        else if (potential_flag.eq.1) then
            write(fileunit, '(1x,i8,a,f15.4,a,f15.4,a,f15.4,a,f15.4,a,f10.4,a'//&
                      ',f19.15,a,f19.15,a,f19.15,a,f19.15,a,f19.15,a,'//&
                      'f10.4,a,f10.4,a,f10.4)') &
            iter,';',simtime,';',vsum,';', mv2sum,';', temperature,';', &
            kinenergy,';',potenergy_LJ,';',potenergy_POLY,';',potenergy,&
            ';',totenergy,';',pressure,';',etevtcf,';',R_g
        end if
        close(unit=fileunit, status='keep')
    endif

end subroutine macroscopic_properties_record

!===================================================================================
!Integrate virial pressure to get autocorrelations (Green Kubo) viscosity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine viscosity_io
    use module_parallel_io
    use physical_constants_MD
    use librarymod, only : integrate_trap
    implicit none

    integer             :: m, length
    real(kind(0.d0))    :: viscosity

    !call globalAverage(Pxycorrel, Nvisc_ave)

    if (irank .eq. iroot) then
        call integrate_trap(Pxycorrel,tplot*delta_t,Nstress_ave,viscosity)

        viscosity = (viscosity*volume)/(3.0*Nstress_ave*Nvisc_ave*inputtemperature)

        !Write viscosity to file
        m = get_iter()/(tplot*Nstress_ave*Nvisc_ave)
        inquire(iolength=length) viscosity
        open (unit=7, file=trim(prefix_dir)//'results/visc', & 
              form='unformatted',access='direct',recl=length)
        write(7,rec=m) viscosity
        close(7,status='keep')
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine viscosity_io

subroutine viscometrics_io(eta,N1,N2,P)
implicit none

    real(kind(0.d0)), intent(in) :: eta, N1, N2, P

    !Remove all these lines, they are only here to satisfy pedantic compiler    
    real(kind(0.d0)) :: dummy
    dummy = eta + N1 + N2 + P

    print*, 'Viscometrics_io not developed in parallel'

end subroutine viscometrics_io

!-----------------------------------------------------------------------------
! Write end-to-end vector time correlation function
!-----------------------------------------------------------------------------
subroutine etev_io
    use module_parallel_io
    use librarymod, only: get_new_fileunit, get_Timestep_FileName
    implicit none

    integer :: etev_fid, record, length
    character(200) :: outfile_t

    !NB we need the etev at t=0, so it is recorded at the START of the usual
    !averaging interval. If, say, mass_Nave is 1, there will be 1 more record
    !of etevs than mbins (due to etev first record at 0, mbins first record
    !at tplot).
    record = real(get_iter())/real(tplot*1) + 1 !etev Nave must be 1
    if (separate_outfiles) then
        !Either write a separate output files for each timestep
        call get_Timestep_FileName(record-1,'etev',outfile_t)
        record=1
    else
        !or a single file for whole run
        outfile_t = 'etev' 
    endif

    if (irank .eq. iroot) then

        etev_fid = get_new_fileunit()
        inquire(iolength=length) etev 
        open(unit=etev_fid,file=trim(prefix_dir)//'results/'//outfile_t,&
             form='unformatted',access='direct',recl=length)
        write(etev_fid,rec=record) etev 
        close(etev_fid,status='keep') 
    
        !if (iter.eq.etevtcf_iter0) then
        !   call MPI_FILE_DELETE(trim(prefix_dir)//'results/'//outfile_t, MPI_INFO_NULL, ierr)
        !   call MPI_FILE_OPEN(MPI_COMM_SELF, trim(prefix_dir)//'results/etev', &
        !                      MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, etev_fid, ierr)
        !   call MPI_FILE_WRITE(etev_fid,etev,nchains*nd,MPI_DOUBLE_PRECISION, &
        !                       MPI_STATUS_IGNORE, ierr)
        !else if (iter.gt.etevtcf_iter0) then
        !   call MPI_FILE_OPEN(MPI_COMM_SELF,trim(prefix_dir)//'results/'//outfile_t, &
        !                      MPI_MODE_WRONLY + MPI_MODE_APPEND, MPI_INFO_NULL, etev_fid, ierr)
        !   call MPI_FILE_WRITE(etev_fid,etev,nchains*nd,MPI_DOUBLE_PRECISION, &
        !                       MPI_STATUS_IGNORE, ierr)
        !end if
        !
        !call MPI_FILE_CLOSE(etev_fid,ierr)

    endif

end subroutine etev_io

!Write radius of gyration to output file
subroutine r_gyration_io
use module_parallel_io
implicit none

    integer :: Rg_fid

    if (irank .eq. iroot) then

        if (iter .eq. r_gyration_iter0) then

            call MPI_FILE_DELETE    (trim(prefix_dir)//'results/r_gyration', MPI_INFO_NULL, ierr)
            call MPI_FILE_OPEN      (MPI_COMM_SELF, trim(prefix_dir)//'results/r_gyration',       &
                                     MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,            &
                                     Rg_fid, ierr)
            call MPI_FILE_WRITE_ALL (Rg_fid, R_g, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr) 

        else if (iter .gt. r_gyration_iter0) then

            call MPI_FILE_OPEN      (MPI_COMM_SELF, trim(prefix_dir)//'results/r_gyration',       &
                                     MPI_MODE_WRONLY + MPI_MODE_APPEND, MPI_INFO_NULL,            &
                                     Rg_fid, ierr)
            call MPI_FILE_WRITE_ALL (Rg_fid, R_g, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr) 

        end if
        
        call MPI_FILE_CLOSE(Rg_fid, ierr)

    end if
    
end subroutine r_gyration_io 

subroutine rdf3d_io
    use module_parallel_io
	use calculated_properties_MD, only: rdf3d
    implicit none

	integer :: m
	integer :: length

    print*, 'Need to parallelise RDF i/o'

	m = 1
	inquire(iolength=length) rdf3d

	if (iter.eq.0) then
		open(16,file=trim(prefix_dir)//'results/rdf3d', &
		     form = 'unformatted', access = 'direct'       , &
		     status = 'replace',   recl = length              )
		write(16, rec = m) rdf3d
	else
		open(16,file=trim(prefix_dir)//'results/rdf3d', &
		     form = 'unformatted', access = 'direct'       , &
		                           recl = length              )
		write(16, rec = m) rdf3d	
	end if

end subroutine rdf3d_io

subroutine rdf_io
    use module_parallel_io
	use calculated_properties_MD, only: rdf
    implicit none

	integer :: m
	integer :: length
    print*, 'Need to parallelise RDF i/o'

	m = 1
	inquire(iolength=length) rdf

	if (iter.eq.0) then
		open(16,file=trim(prefix_dir)//'results/rdf', &
		     form = 'unformatted', access = 'direct'       , &
		     status = 'replace',   recl = length              )
		write(16, rec = m) rdf
	else
		open(16,file=trim(prefix_dir)//'results/rdf', &
		     form = 'unformatted', access = 'direct'       , &
		                           recl = length              )
		write(16, rec = m) rdf	
	end if

end subroutine rdf_io

subroutine ssf_io 
    use module_parallel_io
    implicit none

    print*, 'Need to parallelise SSF i/o'

end subroutine ssf_io

