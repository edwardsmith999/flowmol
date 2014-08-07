!-----------------------------------------------------------------------------
!
!                           Continuum Record
!
!-----------------------------------------------------------------------------

module module_continuum_record

	use physical_constants
	!use computational_constants_MD
	use computational_constants
	use grid_arrays
        use continuum_data_export, only : prefix_dir

contains

    !------------------------------------------------------------------------------
    !Pure fortran subroutine to return an updated filename by appending
    !the current timestep to that file
    subroutine get_Timestep_FileName(timestep,basename,filename)
		    implicit none

		    integer,intent(in) 			:: timestep
		    character(*),intent(in) 	:: basename
		    character(*),intent(out)	:: filename

            if(timestep.le.9                         		) &
            write(filename,'(a,a7,i1)') trim(basename),'.000000',timestep
            if(timestep.ge.10      .and. timestep.le.99     ) &
            write(filename,'(a,a6,i2)') trim(basename),'.00000' ,timestep
            if(timestep.ge.100     .and. timestep.le.999    ) &
            write(filename,'(a,a5,i3)') trim(basename),'.0000'  ,timestep
            if(timestep.ge.1000    .and. timestep.le.9999   ) &
            write(filename,'(a,a4,i4)') trim(basename),'.000'   ,timestep
            if(timestep.ge.10000   .and. timestep.le.99999  ) &
            write(filename,'(a,a3,i5)') trim(basename),'.00'    ,timestep
            if(timestep.ge.100000  .and. timestep.le.999999 ) &
            write(filename,'(a,a2,i6)') trim(basename),'.0'     ,timestep
            if(timestep.ge.1000000 .and. timestep.le.9999999) &
            write(filename,'(a,a1,i7)') trim(basename),'.'      ,timestep

		    !Remove any surplus blanks
		    filename = trim(filename)

    end subroutine get_Timestep_FileName


end module module_continuum_record
!----------------------------------------------------------------------------------

!=============================================================================
!                               Setup initial record
! Print initial output for continuum simulation 
!
!-----------------------------------------------------------------------------

subroutine continuum_initial_record
	use module_continuum_record
    use messenger, only : irank, iroot
	implicit none

    character               :: ixyz_char
    integer                 :: i,n,missing_file_tolerance=5
    logical                 :: file_exist
	Character(8)  		    :: the_date
	Character(10)  		    :: the_time
    character(25)           :: file_names_t
    character(16),parameter :: file_names(5) = &
                                (/ "continuum_vbins ",  &
                                   "continuum_tau_xx",  &
                                   "continuum_tau_xy",  &
                                   "continuum_tau_yx",  &
                                   "continuum_tau_yy" /) 

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

	call date_and_time(the_date, the_time)
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
	call continuum_CFL
	select case (continuum_vflag)
	case(0)
		print*, ' ' 
	case(1)
		print*, 'Continuum velocity output slice in x'
	case(2)
		print*, 'Continuum velocity output slice in y'
	case(3)
		print*, 'Continuum velocity 3D bins'
        case(4)
                print*, 'Continuum velocity in text format, gnuplot blocks'
	case default
		stop "Continuum velocity bin input incorrect"
	end select
	print*, '======================================================================='

	!Open and write header for output file
	open(1002,file=trim(prefix_dir)//'results/continuum_header')

	write(1002,*) 'Simulation run on Date;  sim_date ;', the_date
	write(1002,*) 'Simulation start time ;  sim_start_time ;', the_time
	write(1002,*) 'Number of cells in Domain in x ;  nx ;', nx
	write(1002,*) 'Number of cells in Domain in y ;  ny ;', ny
	write(1002,*) 'Number of cells in Domain in z ;  nz ;', nz
	write(1002,*) 'Domain size in x ; lx ;', lx
	write(1002,*) 'Domain size in y ; ly ;', ly
	write(1002,*) 'Domain size in z ; lz ;', lz
	write(1002,*) 'Time Step - delta t ;   continuum_delta_t ;',  	continuum_delta_t
	write(1002,*) 'Total number of steps ; continuum_Nsteps;',  	continuum_Nsteps
	write(1002,*) 'Continuum initial step ; continuum_initialstep;',continuum_initialstep
	write(1002,*) 'Generate output file every steps ;  continuum_tplot ;',  continuum_tplot
	write(1002,*) 'Viscosity ; meu ;',  	meu
	write(1002,*) 'Density ; rho ;',  	rho
	write(1002,*) 'Reynolds Number ; Re ;', Re
	write(1002,*) 'Continuum velocity output flag; continuum_vflag ;', continuum_vflag
	write(1002,*) 'Finite difference of Finite Volume; solver ;', solver

	close(1002,status='keep')

end subroutine continuum_initial_record

!=============================================================================
!                               Setup record
! Print output for continuum simulation 
!
!-----------------------------------------------------------------------------

subroutine continuum_record
	use module_continuum_record
	!use calculated_properties_MD
	implicit none

	integer										    :: i,j, m, length, ixyz
	double precision, allocatable, dimension(:,:)	:: uslice,tauslice
	double precision, allocatable, dimension(:,:,:)	:: buffer

	select case(continuum_vflag)
	case(0)
		!No Output for velocity
	case(1:2)
		!Write Continuum velocities
		m = continuum_iter/continuum_tplot

		!Get directions orthogonal to slice direction
		ixyz = mod(continuum_vflag,2)+1

		if (continuum_vflag .eq. 1) then
			allocate(uslice(nx+2,3))
			uslice(:,1) = sum(uc(1:nx+2,2:ny+1),dim=ixyz)/(ny)
			uslice(:,2) = sum(vc(1:nx+2,2:ny+1),dim=ixyz)/(ny)
			uslice(:,3) = 0.d0
		elseif (continuum_vflag .eq. 2) then
			allocate(uslice(ny+2,3))
			uslice(:,1) = sum(uc(2:nx+1,1:ny+2),dim=ixyz)/(nx)
			uslice(:,2) = sum(vc(2:nx+1,1:ny+2),dim=ixyz)/(nx)
			uslice(:,3) = 0.d0
		endif

		inquire(iolength=length) uslice
		open (unit=1003, file=trim(prefix_dir)//"results/continuum_vslice",form="unformatted",access='direct',recl=length)
		write(1003,rec=m) uslice
		close(1003,status='keep')
		deallocate(uslice)

        if (solver .eq. FV) then
			    tauslice = get_slice(continuum_vflag,tau_xx)
    		    inquire(iolength=length) tauslice
		        open (unit=1004, file=trim(prefix_dir)//"results/continuum_tauslice_xx",form="unformatted",access='direct',recl=length)
		        write(1004,rec=m) tauslice
		        close(1004,status='keep')
    		    deallocate(tauslice)
			    tauslice = get_slice(continuum_vflag,tau_xy)
    		    inquire(iolength=length) tauslice
		        open (unit=1004, file=trim(prefix_dir)//"results/continuum_tauslice_xy",form="unformatted",access='direct',recl=length)
		        write(1004,rec=m) tauslice
		        close(1004,status='keep')
    		    deallocate(tauslice)
			    tauslice = get_slice(continuum_vflag,tau_yx)
    		    inquire(iolength=length) tauslice
		        open (unit=1004, file=trim(prefix_dir)//"results/continuum_tauslice_yx",form="unformatted",access='direct',recl=length)
		        write(1004,rec=m) tauslice
		        close(1004,status='keep')
    		    deallocate(tauslice)
			    tauslice = get_slice(continuum_vflag,tau_yy)
    		    inquire(iolength=length) tauslice
		        open (unit=1004, file=trim(prefix_dir)//"results/continuum_tauslice_yy",form="unformatted",access='direct',recl=length)
		        write(1004,rec=m) tauslice
		        close(1004,status='keep')
    		    deallocate(tauslice)

        endif


	case(3)
		!Write Continuum velocities
		m = continuum_iter/continuum_tplot

		allocate(buffer(nx,ny+2,3))
		buffer(:,:,1) = uc(2:nx+1,:)
		buffer(:,:,2) = vc(2:nx+1,:)
		buffer(:,:,3) = 0.d0
		inquire(iolength=length) buffer
		open (unit=1003, file=trim(prefix_dir)//"results/continuum_vbins",form="unformatted",access='direct',recl=length)
		write(1003,rec=m) buffer
		close(1003,status='keep')

        !Write Stresses
        if (solver .eq. FV) then

		    inquire(iolength=length) tau_xx(2:nx+1,:,:)
    		open (unit=1004, file=trim(prefix_dir)//"results/continuum_tau_xx",form="unformatted",access='direct',recl=length)
	    	write(1004,rec=m) tau_xx(2:nx+1,:,:)
	    	close(1004,status='keep')
    		open (unit=1004, file=trim(prefix_dir)//"results/continuum_tau_xy",form="unformatted",access='direct',recl=length)
	    	write(1004,rec=m) tau_xy(2:nx+1,:,:)
	    	close(1004,status='keep')
    		open (unit=1004, file=trim(prefix_dir)//"results/continuum_tau_yx",form="unformatted",access='direct',recl=length)
	    	write(1004,rec=m) tau_yx(2:nx+1,:,:)
	    	close(1004,status='keep')
    		open (unit=1004, file=trim(prefix_dir)//"results/continuum_tau_yy",form="unformatted",access='direct',recl=length)
	    	write(1004,rec=m) tau_yy(2:nx+1,:,:)
	    	close(1004,status='keep')

        endif

	case default
		stop "Error incorrect input for continuum velocity "
	end select

	!Output simulation progress and residual
	print'(a,i8,2(a,f20.15))', ' elapsed continnum steps ', continuum_iter, &
		 ' residual x, y = ' , sum(xresidual)/nx, ',', sum(yresidual)/ny

contains

    function get_slice(dir,array) result(slice)

	    integer,intent(in)							                :: dir
	    double precision, allocatable, dimension(:,:,:),intent(in)	:: array

	    double precision, allocatable, dimension(:,:)              	:: slice

        allocate(slice(size(array,dir),size(array,3)))

		!Get directions orthogonal to slice direction
		ixyz = mod(dir,2)+1

	    if (dir .eq. 1) then
		    allocate(tauslice(nx+2,4))
		    slice = sum(array(1:nx+2,2:ny+1,:),dim=ixyz)/(ny)
	    elseif (dir .eq. 2) then
		    allocate(tauslice(ny+2,4))
            slice = sum(array(2:nx+1,1:ny+2,:),dim=ixyz)/(nx)
	    endif

    end function get_slice

end subroutine continuum_record


