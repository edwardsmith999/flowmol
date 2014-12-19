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
#if USE_COUPLER
    use continuum_coupler_socket, only : socket_coupler_get_md_BC
#endif
	implicit none
	
	integer					:: i, j, n
	logical, save 			:: first_time = .true.

    if (first_time) then
		first_time = .false.
		allocate(u_MD(nx), v_MD(nx))
    endif

	!Loop over all 4 domain boundaries
#if USE_COUPLER
	call socket_coupler_get_md_BC(u_MD,v_MD)
#endif
	do i = 1,4
		select case(BC_flag(i))
		case(0)
			call periodic_BC(i)
		case(1)
			call Dirichlet_BC_set_wall(i,uwall_BC(i),0.d0)
		case(2)
			call Dirichlet_BC_set_halo(i,uwall_BC(i),0.d0)
		case(3)
			call Dirichlet_slip_BC(i)
		case(4)
			call Von_Neumann_BC(i,0.d0,0.d0,0.d0,0.d0)
			!stop "No Von Neumann BC yet"
		case(5)
			!Coupled boundary condition applied
			call Dirichlet_BC_set_halo_cpl(i,u_MD,v_MD)
		case default
			stop "Error in choice of continuum BC"
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


subroutine Dirichlet_BC_set_halo(boundary,uwall,vwall)
	use module_setup_BC
	implicit none

	integer			    :: i, j
	integer			    :: boundary
	double precision	:: uwall, vwall

	!print*, vc(1,2),vc(2,2), vc(nx+2,2),vc(nx+1,2)

	select case(boundary)
	case(1)
		!Right BC
		i = nx+2
		do j = 2 , ny+1
			uc(i,j) = uwall
			vc(i,j) = vwall
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(2)
		!Top BC
		j = ny+2
		do i = 2 , nx+1
			uc(i,j) = uwall
			vc(i,j) = vwall
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(3)
		!Left BC
		i = 1
		do j = 2 , ny+1
			uc(i,j) = uwall
			vc(i,j) = vwall
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(4)

		!Bottom BC
		j = 1
		do i = 2 , nx+1
			uc(i,j) = uwall
			vc(i,j) = vwall
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	end select

	!print*, uc(2:nx+1,1)

end subroutine Dirichlet_BC_set_halo


!----------------------------------------------------------------------------------
!Wall/flow boundary conditions - set halo to specified value

subroutine Dirichlet_BC_set_halo_cpl(boundary,uwall,vwall)
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

	!print*, uc(2:nx+1,1)

end subroutine Dirichlet_BC_set_halo_cpl


!----------------------------------------------------------------------------------
!Set Wall to slip boundary conditions 

subroutine Dirichlet_slip_BC(boundary)
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
			uc(i,j) = uc(i-1,j)
			vc(i,j) = vc(i-1,j)
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(2)
		!Top BC
		j = ny+2
		do i = 2 , nx+1
			uc(i,j) = uc(i,j-1) 
			vc(i,j) = vc(i,j-1)
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(3)
		!Left BC
		i = 1
		do j = 2 , ny+1
			uc(i,j) = uc(i+1,j) 
			vc(i,j) = vc(i+1,j)
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	case(4)

		!Bottom BC
		j = 1
		do i = 2 , nx+1
			uc(i,j) = uc(i,j+1) 
			vc(i,j) = vc(i,j+1) 
			flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
			flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
			flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
			flux_yy(i,j) = rho * vc(i,j) * vc(i,j)
		enddo
	end select

end subroutine Dirichlet_slip_BC

!----------------------------------------------------------------------------------
subroutine Von_Neumann_BC(boundary,	wall_flux_xx,wall_flux_xy, &
								 	wall_flux_yx,wall_flux_yy )
	use module_setup_BC
	implicit none

	integer				:: i, j
	integer				:: boundary
	double precision	:: wall_flux_xx,wall_flux_xy, & 
						   wall_flux_yx,wall_flux_yy

	select case(boundary)
	case(1)
		!Right BC
		i = nx+2
		do j = 2 , ny+1
			flux_xx(i,j) = wall_flux_xx
			flux_xy(i,j) = wall_flux_xy
			flux_yx(i,j) = wall_flux_yx
			flux_yy(i,j) = wall_flux_yy
		enddo
	case(2)
		!Top BC
		j = ny+2
		do i = 2 , nx+1
			flux_xx(i,j) = wall_flux_xx
			flux_xy(i,j) = wall_flux_xy
			flux_yx(i,j) = wall_flux_yx
			flux_yy(i,j) = wall_flux_yy
		enddo
	case(3)
		!Left BC
		i = 1
		do j = 2 , ny+1
			flux_xx(i,j) = wall_flux_xx
			flux_xy(i,j) = wall_flux_xy
			flux_yx(i,j) = wall_flux_yx
			flux_yy(i,j) = wall_flux_yy
		enddo
	case(4)

		!Bottom BC
		j = 1
		do i = 2 , nx+1
			flux_xx(i,j) = wall_flux_xx
			flux_xy(i,j) = wall_flux_xy
			flux_yx(i,j) = wall_flux_yx
			flux_yy(i,j) = wall_flux_yy
		enddo
	end select

end subroutine Von_Neumann_BC

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