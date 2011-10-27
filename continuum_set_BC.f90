!=============================================================================
!                               Setup Boundary Conditions
! Setup Halo cells which give boundary 
!
!-----------------------------------------------------------------------------

module module_setup_BC

	use physical_constants
	use computational_constants
	use grid_arrays

	double precision, allocatable			:: u_MD(:), v_MD(:)

end module module_setup_BC
!----------------------------------------------------------------------------------

subroutine continuum_set_BC
	use module_setup_BC
        use continuum_coupler_socket_init, only : use_coupling
        use continuum_coupler_socket, only : MD_continuum_BC
	implicit none
	
	integer					:: i, j, n
        logical, save :: first_time = .true.
	!integer					:: averng
	!integer				:: xobs, yobs
	!double precision			:: x_interval,lstsqrsinter,lstsqrsgrad
	!double precision			:: uwall,uwall_LS 

	!Add one to velocity count to prevent division by zero
	!slice_mass = slice_mass + 1	

	!Average from averng above to averng below cell corresponding to halo
	!averng = 0
	!do n = nbins(2)-overlap-averng, nbins(2)-overlap+averng
	!	u_MD = u_MD + volume_momentum(:,n,:,1)
	!enddo
	!u_MD = u_MD / (1.d0 + 2.d0 * averng)

	!Predict using least squares the expect slope based on continuum constraint forces
	!y(:) = slice_momentum((nbins(1)-overlap+1):nbins(1))/slice_mass((nbins(1)-overlap+1):nbins(1))
	!x_interval = 1.d0
	!call least_squares(y,x_interval,overlap,lstsqrsinter,lstsqrsgrad)
	!uwall_LS = lstsqrsgrad * 0.d0 + lstsqrsinter

	!print'(a,i8,a,f10.5,a,f10.5)', '; ; ; ; ; ; ; ; ', iter, '; Least square prediction of uwall ;', uwall_LS, '; uwall MD ;', uwall

	!print'(a,f10.5,a,i3,a,f10.5)', 'Velocity at MD Cell corresponding to Continuum halo location = ', &
	!	 slice_momentum(nbins(1)-overlap)/slice_mass(nbins(1)-overlap), &
	!	' & velocity average over ',averng,' cells above and below = ',  uwall	


        if (first_time) then
                 first_time= .false.
                 allocate(u_MD(nx), v_MD(nx), stat=i)
        endif

	!Loop over all 4 domain boundaries
	do i = 1,4
		select case(BC_flag(i))
		case(0)
			call periodic_BC(i)
		case(1)
			call Dirichlet_BC_set_wall(i,uwall_BC(i),0.d0)
		case(2)
			stop "No Von Neumann BC yet"
		case(3)
                        if (use_coupling) then 
			        call MD_continuum_BC(u_MD,v_MD)
                        else
			!print*, 'u_MD passed to BC', u_MD
                                u_MD = 0.d0; v_MD =0.0
                        endif
                        call Dirichlet_BC_set_halo(i,u_MD,v_MD)
		end select
	enddo

	!Boundaries set up on domain
	!call Dirichlet_BC(1,0.d0,0.d0)
	!call periodic_BC(1)
	!call Dirichlet_BC(2,uwall_BC(2),0.d0)
	!call periodic_BC(2)
	!call Dirichlet_BC(3,0.d0,0.d0)
	!call periodic_BC(3)
	!call Dirichlet_BC(4,uwall_BC(4),0.d0)
	!call periodic_BC(4)

	!Obstacle set up in flow
	!xobs = floor(nx/2.d0)
	!yobs = floor(ny/2.d0)
	
	!call obstacle(xobs,yobs,0)
	!do i=2,nx+1
	!do j=2,ny+1
	!	call obstacle(i, j,solid(i-1,j-1))
	!enddo
	!enddo

end subroutine continuum_set_BC

!----------------------------------------------------------
!Periodic Boundary conditions

subroutine periodic_BC(boundary)
	use module_setup_BC
	implicit none

	integer	:: i, j
	integer	:: boundary


	select case(boundary)
	case(1)
		!Right BC
		i = nx+2
		do j = 1 , ny+2
			uc(i,j)  = uc(i-nx,j)
			vc(i,j)  = vc(i-nx,j)
			flux_xx(i,j) = flux_xx(i-nx,j)
			flux_xy(i,j) = flux_xy(i-nx,j)
			flux_yx(i,j) = flux_yx(i-nx,j)
			flux_yy(i,j) = flux_yy(i-nx,j)
		enddo

	case(2)
		!Top BC
		j = ny+2
		do i = 1 , nx+2
			uc(i,j)  = uc(i,j-ny)
			vc(i,j)  = vc(i,j-ny)
			flux_xx(i,j) = flux_xx(i,j-ny)
			flux_xy(i,j) = flux_xy(i,j-ny)
			flux_yx(i,j) = flux_yx(i,j-ny)
			flux_yy(i,j) = flux_yy(i,j-ny)
		enddo

	case(3)
		!Left BC
		i = 1
		do j = 1 , ny+2
			uc(i,j)  = uc(i+nx,j)
			vc(i,j)  = vc(i+nx,j)
			flux_xx(i,j) = flux_xx(i+nx,j)
			flux_xy(i,j) = flux_xy(i+nx,j)
			flux_yx(i,j) = flux_yx(i+nx,j)
			flux_yy(i,j) = flux_yy(i+nx,j)
		enddo

	case(4)

		!Bottom BC
		j = 1
		do i = 1 , nx+2
			uc(i,j)  = uc(i,j+ny)
			vc(i,j)  = vc(i,j+ny)
			flux_xx(i,j) = flux_xx(i,j+ny)
			flux_xy(i,j) = flux_xy(i,j+ny)
			flux_yx(i,j) = flux_yx(i,j+ny)
			flux_yy(i,j) = flux_yy(i,j+ny)
		enddo

	end select


end subroutine periodic_BC

!----------------------------------------------------------------------------------
!Wall/flow boundary conditions - set halo so that wall is specified value

subroutine Dirichlet_BC_set_wall(boundary,uwall,vwall)
	use module_setup_BC
	implicit none

	integer			:: i, j
	integer			:: boundary
	double precision	:: uwall, vwall

	select case(boundary)
	case(1)
		!Right BC
		i = nx+2
		do j = 2 , ny+1
			uc(i,j) = -uc(i-1,j) + uwall*2.d0
			vc(i,j) = -vc(i-1,j) + vwall*2.d0
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(2)
		!Top BC
		j = ny+2
		do i = 2 , nx+1
			uc(i,j) = -uc(i,j-1) + uwall*2.d0
			vc(i,j) = -vc(i,j-1) + vwall*2.d0
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(3)
		!Left BC
		i = 1
		do j = 2 , ny+1
			uc(i,j) = -uc(i+1,j) + uwall*2.d0
			vc(i,j) = -vc(i+1,j) + vwall*2.d0
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(4)

		!Bottom BC
		j = 1
		do i = 2 , nx+1
			uc(i,j) = -uc(i,j+1) + uwall*2.d0
			vc(i,j) = -vc(i,j+1) + vwall*2.d0
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	end select

end subroutine Dirichlet_BC_set_wall

!----------------------------------------------------------------------------------
!Wall/flow boundary conditions - set halo to specified value

subroutine Dirichlet_BC_set_halo(boundary,uwall,vwall)
	use module_setup_BC
	implicit none

	integer			:: i, j
	integer			:: boundary
	double precision	:: uwall(*), vwall(*)

	!print*, vc(1,2),vc(2,2), vc(nx+2,2),vc(nx+1,2)

	select case(boundary)
	case(1)
		!Right BC
		i = nx+2
		do j = 2 , ny+1
			uc(i,j) = uwall(j-1)
			vc(i,j) = vwall(j-1)
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(2)
		!Top BC
		j = ny+2
		do i = 2 , nx+1
			uc(i,j) = uwall(j-1)
			vc(i,j) = vwall(j-1)
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(3)
		!Left BC
		i = 1
		do j = 2 , ny+1
			uc(i,j) = uwall(j-1)
			vc(i,j) = vwall(j-1)
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(4)

		!Bottom BC
		j = 1
		do i = 2 , nx+1
			uc(i,j) = uwall(i-1)
			vc(i,j) = vwall(i-1)
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	end select

	print*, uc(2:nx+1,1)

end subroutine Dirichlet_BC_set_halo

!---------------------------------------------------------------------

subroutine obstacle(xlocation,ylocation,fix)
	use module_setup_BC
	implicit none
	
	integer	:: xlocation, ylocation,fix

	uc(xlocation,ylocation) = uc(xlocation,ylocation)*fix
	vc(xlocation,ylocation) = vc(xlocation,ylocation)*fix
	flux_xx(xlocation,ylocation) = rho * uc(xlocation,ylocation) * uc(xlocation,ylocation)
	flux_xy(xlocation,ylocation) = rho * uc(xlocation,ylocation) * vc(xlocation,ylocation)
	flux_yx(xlocation,ylocation) = rho * vc(xlocation,ylocation) * uc(xlocation,ylocation)
	flux_yy(xlocation,ylocation) = rho * vc(xlocation,ylocation) * vc(xlocation,ylocation)
	
	return
	
end
