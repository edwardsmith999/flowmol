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

end module module_continuum_record
!----------------------------------------------------------------------------------

!=============================================================================
!                               Setup initial record
! Print initial output for continuum simulation 
!
!-----------------------------------------------------------------------------

subroutine continuum_initial_record
	use module_continuum_record
	implicit none

	Character(8)  		:: the_date
	Character(10)  		:: the_time

	!Delete existing files
	open (unit=1003, file=trim(prefix_dir)//"results/continuum_vslice")
	close(1003,status='delete')
	open (unit=1004, file=trim(prefix_dir)//"results/continuum_vxbins")
	close(1004,status='delete')

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
	write(1002,*) 'Domain size in x ; lx ;', lx
	write(1002,*) 'Domain size in y ; ly ;', ly
	write(1002,*) 'Time Step - delta t ;   continuum_delta_t ;',  	continuum_delta_t
	write(1002,*) 'Total number of steps ; continuum_Nsteps;',  	continuum_Nsteps
	write(1002,*) 'Continuum initial step ; continuum_initialstep;',continuum_initialstep
	write(1002,*) 'Generate output file every steps ;  continuum_tplot ;',  continuum_tplot
	write(1002,*) 'Viscosity ; meu ;',  	meu
	write(1002,*) 'Density ; rho ;',  	rho
	write(1002,*) 'Reynolds Number ; Re ;', Re
	write(1002,*) 'Continuum velocity output flag; continuum_vflag ;', continuum_vflag

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

	integer						:: m, length, ixyz
	double precision, allocatable, dimension(:)	:: uslice

	select case(continuum_vflag)
	case(0)
		!No Output for velocity
	case(1:2)
		!Write Continuum velocities
		m = continuum_iter/continuum_tplot

		!Get directions orthogonal to slice direction
		ixyz = mod(continuum_vflag,2)+1

		if (continuum_vflag .eq. 1) then
			allocate(uslice(nx+2))
			uslice = sum(uc(1:nx+2,2:ny+1),dim=ixyz)/(ny)
		else 
			allocate(uslice(ny+2))
			uslice = sum(uc(2:nx+1,1:ny+2),dim=ixyz)/(nx)
		endif

		inquire(iolength=length) uslice
		open (unit=1003, file=trim(prefix_dir)//"results/continuum_vslice",form="unformatted",access='direct',recl=length)
		write(1003,rec=m) uslice
		close(1003,status='keep')

		deallocate(uslice)
	case(3)
		!Write Continuum velocities
		m = continuum_iter/continuum_tplot

		inquire(iolength=length) uc(1:nx+2,1:ny+2)
		open (unit=1003, file=trim(prefix_dir)//"results/continuum_vxbins",form="unformatted",access='direct',recl=length)
		write(1003,rec=m) uc(1:nx+2,1:ny+2)
		close(1003,status='keep')

		!inquire(iolength=length) vc(1:nx+2,1:ny+2)
		!open (unit=1004, file="results/continuum_vybins",form="unformatted",access='direct',recl=length)
		!write(1004,rec=m) vc(1:nx+2,1:ny+2)
		!close(1004,status='keep')
        case(4)  ! ASCII output
                call text_output
	case default
		stop "Error input for continuum velocity incorrect"
	end select

	!Output simulation progress and residual
	print'(a,i8,a,f20.15)', ' elapsed continnum steps ', continuum_iter, &
		 ' residual = ' , sum(xresidual)/nx + sum(yresidual)/ny
contains

        subroutine text_output
                implicit none
 

                integer j
                logical, save :: firsttime=.true.
                character(len=32) file_position

                if (firsttime) then
                        firsttime = .false.
                        file_position = "rewind"
                else
                        file_position = "append"
                endif

                open (unit=1003, file=trim(prefix_dir)//"results/continuum_uc.txt",position=file_position)

                write(1003,'(a,I6)')'# step', continuum_iter
                do j=1,ny+2
                        write(1003,'(1000E12.4)') uc(:,j)
                enddo
                write(1003,'(1x/1x)')
                close(1003)
                

       end subroutine text_output

end subroutine continuum_record


