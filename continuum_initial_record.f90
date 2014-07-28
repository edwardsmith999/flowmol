!=============================================================================
!                               Setup initial record
! Print initial output for continuum simulation 
!
!-----------------------------------------------------------------------------

module module_continuum_initial_record

	use physical_constants
	use computational_constants
	use grid_arrays
	!use calculated_properties_MD
	!use computational_constants_MD

end module module_continuum_initial_record
!----------------------------------------------------------------------------------

subroutine continuum_initial_record
	use module_continuum_initial_record
	implicit none	

    character(5)           :: file_names_t
    character(5),parameter :: file_names(16) = &
                                (/ "continuum_vbins ",  &
                                   "continuum_tau_xx",  &
                                   "continuum_tau_xy",  &
                                   "continuum_tau_yx",  &
                                   "continuum_tau_yy", /) 


	print*, '=================Continuum Simulation Parameters======================'
	print'(3(a,f10.5))', 'Specified Size in x ', lx , ' & y ', ly, ' Volume of domain ', domain_volume
	print'(2(a,f10.5))', 'Resulting computaional domain size in x ', mx(nx+2)-mx(2),' & y ',my(ny+2)-my(2)
	print'(2(a,i8))', 'Number of points in x ', nx, ' Number of points in y ', ny
	print'(2(a,f10.5))','Cell size in x ', delta_x(1), ' Cell size in y ', delta_y(1)
	print*, 'number of timesteps', continuum_Nsteps
	print*, 'frequency of plots', continuum_tplot
	print*, 'meu', meu
	print*, 'rho', rho

	print*, 'continuum Reynolds number =', Re
	print*, '======================================================================='


    !Delete all files from previous run
    if (irank.eq.iroot) then
        do i=1,size(file_names)
            inquire(file=trim(prefix_dir)//'results/'//file_names(i),exist=file_exist)
            if(file_exist) then
                open (unit=23, file=trim(prefix_dir)//'results/'//file_names(i))
                close(23,status='delete')
            endif
            !Remove indivdual files -- Keep looping until no further increase in number
            do n = 0,9999999
                call get_Timestep_FileName(n,file_names(i),file_names_t)
                inquire(file=trim(prefix_dir)//'results/'//file_names_t,exist=file_exist)
                if(file_exist) then
                    open (unit=23, file=trim(prefix_dir)//'results/'//file_names_t)
                    close(23,status='delete')
                elseif(missing_file_tolerance .eq. 0) then
                    exit !Exit loop if max file reached 
                else
                    missing_file_tolerance = missing_file_tolerance - 1
                endif
            enddo
        enddo
    endif
    call messenger_syncall()

end subroutine continuum_initial_record
