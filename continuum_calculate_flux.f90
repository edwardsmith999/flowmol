!=============================================================================
!                                   calculate flux
!
! Calculate the flux on each cells surface and add to give residual
!
!-----------------------------------------------------------------------------

module module_calculate_flux

	use physical_constants
	use computational_constants
	use grid_arrays

end module module_calculate_flux
!----------------------------------------------------------------------------------

subroutine continuum_calculate_flux
use module_calculate_flux
implicit none

	integer	:: i, j, n

	call continuum_calculate_convection
	call continuum_calculate_diffusion

end subroutine continuum_calculate_flux

!----------------------------------------------------------------------------------


subroutine continuum_calculate_convection
use module_calculate_flux
implicit none

	integer			:: i, j, n, m

	do i = 2, nx + 1
	do j = 2, ny + 1

		flux_xx(i,j) = rho * uc(i,j) * uc(i,j)
		flux_xy(i,j) = rho * uc(i,j) * vc(i,j)
		flux_yx(i,j) = rho * vc(i,j) * uc(i,j)
		flux_yy(i,j) = rho * vc(i,j) * vc(i,j)

	enddo
	enddo

    ! set values on the domain boundary
	flux_xx(1,:) = 0.d0 ; flux_xx(nx+2,:) = 0.d0 
   	flux_xy(1,:) = 0.d0 ; flux_xy(nx+2,:) = 0.d0
    flux_yx(1,:) = 0.d0 ; flux_yx(nx+2,:) = 0.d0
    flux_yy(1,:) = 0.d0 ; flux_yy(nx+2,:) = 0.d0

	flux_xx(:,1) = 0.d0 ; flux_xx(:,ny+2) = 0.d0 
   	flux_xy(:,1) = 0.d0 ; flux_xy(:,ny+2) = 0.d0
    flux_yx(:,1) = 0.d0 ; flux_yx(:,ny+2) = 0.d0
    flux_yy(:,1) = 0.d0 ; flux_yy(:,ny+2) = 0.d0
    
	
	!Evaluate fluxes at cell surfaces
	do i = 2, nx + 1
	do j = 2, ny + 1
	
		!Evaluate flux for face [i+1/2,j]

		!Convective flux
		sflux_xx(i,j,1) = (flux_xx(i+1,j) + flux_xx(i,j))/2.d0
		sflux_xy(i,j,1) = (flux_xy(i+1,j) + flux_xy(i,j))/2.d0
		sflux_yx(i,j,1) = (flux_yx(i+1,j) + flux_yx(i,j))/2.d0
		sflux_yy(i,j,1) = (flux_yy(i+1,j) + flux_yy(i,j))/2.d0

		!Evaluate flux for face [i,j+1/2]

		!Convective fluxcontinuum_calculate_diffusion
		sflux_xx(i,j,2) = (flux_xx(i,j+1) + flux_xx(i,j))/2.d0
		sflux_xy(i,j,2) = (flux_xy(i,j+1) + flux_xy(i,j))/2.d0
		sflux_yx(i,j,2) = (flux_yx(i,j+1) + flux_yx(i,j))/2.d0
		sflux_yy(i,j,2) = (flux_yy(i,j+1) + flux_yy(i,j))/2.d0

		!Evaluate flux for face [i-1/2,j]

		!Convective flux
		sflux_xx(i,j,3) = (flux_xx(i-1,j) + flux_xx(i,j))/2.d0
		sflux_xy(i,j,3) = (flux_xy(i-1,j) + flux_xy(i,j))/2.d0
		sflux_yx(i,j,3) = (flux_yx(i-1,j) + flux_yx(i,j))/2.d0
		sflux_yy(i,j,3) = (flux_yy(i-1,j) + flux_yy(i,j))/2.d0
		
		!Evaluate flux for face [i,j-1/2]

		!Convective flux
		sflux_xx(i,j,4) = (flux_xx(i,j-1) + flux_xx(i,j))/2.d0
		sflux_xy(i,j,4) = (flux_xy(i,j-1) + flux_xy(i,j))/2.d0
		sflux_yx(i,j,4) = (flux_yx(i,j-1) + flux_yx(i,j))/2.d0
		sflux_yy(i,j,4) = (flux_yy(i,j-1) + flux_yy(i,j))/2.d0

	enddo
	enddo

end subroutine continuum_calculate_convection

!----------------------------------------------------------------------------------

subroutine continuum_calculate_diffusion
use module_calculate_flux
implicit none

	integer			:: i, j, n

	!Store previous timesteps
	xresidual_t_minus_1 = xresidual

	call continuum_calculate_residual_FV(uc, vc, xresidual, yresidual)

end subroutine continuum_calculate_diffusion


!----------------------------------------------------------------------------------
!Finite difference formulation
subroutine continuum_calculate_residual(uc_t, vc_t, xres_t, yres_t)
use module_calculate_flux
implicit none

	integer							:: i, j
	double precision, dimension(nx+2,ny+2), intent(in)	:: uc_t, vc_t
	double precision, dimension(nx+2,ny+2), intent(out)	:: xres_t, yres_t

	yres_t = 0.d0

	!Based on formulation in Hirsch (2007) using divergence theorm
	do i = 2, nx + 1
	do j = 2, ny + 1

		xres_t(i,j) = (uc_t(i,j+1)+uc_t(i,j-1)-2*uc_t(i,j))/(delta_y(j)*delta_y(j)) &
			    	+ (uc_t(i+1,j)+uc_t(i-1,j)-2*uc_t(i,j))/(delta_x(i)*delta_x(i))

	enddo
	enddo

	!Divide diffusive term by Reynolds number
	xres_t = xres_t/Re
	
end subroutine continuum_calculate_residual

!----------------------------------------------------------------------------------
!Finite Volume formulation
subroutine continuum_calculate_residual_FV(uc_t, vc_t, xres_t, yres_t)
use module_calculate_flux
implicit none

	integer							:: i, j, n
	double precision, dimension(nx+2,ny+2,4)		:: ucc, vcc
	double precision, dimension(nx+2,ny+2), intent(in)	:: uc_t, vc_t
	double precision, dimension(nx+2,ny+2), intent(out)	:: xres_t, yres_t

	!VARIABLE NAMING USED FOR FINITE VOLUME
	!
	!		    		tau(i,j,2)
	!						|
	!						|
	!		ucc(i,j,2)______V_______ucc(i,j,1)
	!				x				x
	!				|				|
	!tau(i,j,3) ___\|   uc(i,j)     |/___ tau(i,j,1) 
 	!      	       /|		x		|\
	!				|				|
	!				x_______________x
	!			ucc(i,j,3)	^	ucc(i,j,4)
	!						|
	!						|
	!		   			 tau(i,j,4)

	!Calculate velocity at cell corners ~> ucc by averaging surrounding cells
	do i = 2, nx + 1
	do j = 2, ny + 1
		ucc(i,j,1) = 0.25d0*(uc_t(i,j)+uc_t(i+1,j)+uc_t(i,j+1)+uc_t(i+1,j+1))
		ucc(i,j,2) = 0.25d0*(uc_t(i,j)+uc_t(i-1,j)+uc_t(i,j+1)+uc_t(i-1,j+1))
		ucc(i,j,3) = 0.25d0*(uc_t(i,j)+uc_t(i-1,j)+uc_t(i,j-1)+uc_t(i-1,j-1))
		ucc(i,j,4) = 0.25d0*(uc_t(i,j)+uc_t(i+1,j)+uc_t(i,j-1)+uc_t(i+1,j-1))
		!print'(2i,4f10.5)', i,j, ucc(i,j,:)
		vcc(i,j,1) = 0.25d0*(vc_t(i,j)+vc_t(i+1,j)+vc_t(i,j+1)+vc_t(i+1,j+1))
		vcc(i,j,2) = 0.25d0*(vc_t(i,j)+vc_t(i-1,j)+vc_t(i,j+1)+vc_t(i-1,j+1))
		vcc(i,j,3) = 0.25d0*(vc_t(i,j)+vc_t(i-1,j)+vc_t(i,j-1)+vc_t(i-1,j-1))
		vcc(i,j,4) = 0.25d0*(vc_t(i,j)+vc_t(i+1,j)+vc_t(i,j-1)+vc_t(i+1,j-1))
	enddo
	enddo

	!Calculate velocity at domain corners ~> ucc by averaging 3 cells
	ucc(nx+1,ny+1,1) = 0.333d0 * (uc_t(nx+1,ny+1)+uc_t(nx+2,ny+1)+uc_t(nx+1,ny+2))
	ucc(2,ny+1,2)    = 0.333d0 * (uc_t(2   ,ny+1)+uc_t(1   ,ny+1)+uc_t(2   ,ny+2))
	ucc(2,2,3)       = 0.333d0 * (uc_t(2   ,   2)+uc_t(1   ,   2)+uc_t(2   ,   1))
	ucc(nx+1,2,4)    = 0.333d0 * (uc_t(nx+1,   2)+uc_t(nx+2,   2)+uc_t(nx+1,   1))

	!Calculate velocity at domain corners ~> ucc by averaging 3 cells
	vcc(nx+1,ny+1,1) = 0.333d0 * (vc_t(nx+1,ny+1)+vc_t(nx+2,ny+1)+vc_t(nx+1,ny+2))
	vcc(2,ny+1,2)    = 0.333d0 * (vc_t(2   ,ny+1)+vc_t(1   ,ny+1)+vc_t(2   ,ny+2))
	vcc(2,2,3)       = 0.333d0 * (vc_t(2   ,   2)+vc_t(1   ,   2)+vc_t(2   ,   1))
	vcc(nx+1,2,4)    = 0.333d0 * (vc_t(nx+1,   2)+vc_t(nx+2,   2)+vc_t(nx+1,   1))


	!Evaluate Diffusive flux or Stresses at cell surfaces
	do i = 2, nx + 1
	do j = 2, ny + 1

		!Evaluate flux for face 1 [i+1/2,j]

		tau_xx(i,j,1) = (uc_t(i+1,j) - uc_t(i,j))/delta_x(i)
		tau_xy(i,j,1) = (ucc(i,j,4) - ucc(i,j,1))/delta_y(j)     
		tau_yx(i,j,1) = (vc_t(i+1,j) - vc_t(i,j))/delta_x(i)		       	
		tau_yy(i,j,1) = (vcc(i,j,4) - vcc(i,j,1))/delta_y(j) 

		!Evaluate flux for face 2 [i,j+1/2]

		tau_xx(i,j,2) = (ucc(i,j,1) - ucc(i,j,2))/delta_x(i)
		tau_xy(i,j,2) = (uc_t(i,j+1) - uc_t(i,j))/delta_y(j)     
		tau_yx(i,j,2) = (vcc(i,j,1) - vcc(i,j,2))/delta_x(i)		       	
		tau_yy(i,j,2) = (vc_t(i,j+1) - vc_t(i,j))/delta_y(j)


		!Evaluate flux for face 3 [i-1/2,j]

		tau_xx(i,j,3) = (uc_t(i,j) - uc_t(i-1,j))/delta_x(i)
		tau_xy(i,j,3) = (ucc(i,j,3) - ucc(i,j,2))/delta_y(j)      
		tau_yx(i,j,3) = (vc_t(i,j) - vc_t(i-1,j))/delta_x(i)		       	
		tau_yy(i,j,3) = (vcc(i,j,3) - vcc(i,j,2))/delta_y(j)

		!Evaluate flux for face 4 [i,j-1/2]

		tau_xx(i,j,4) = (ucc(i,j,4) - ucc(i,j,3))/delta_x(i)
		tau_xy(i,j,4) = (uc_t(i,j) - uc_t(i,j-1))/delta_y(j)     
		tau_yx(i,j,4) = (vcc(i,j,4) - vcc(i,j,3))/delta_x(i)		       	
		tau_yy(i,j,4) = (vc_t(i,j) - vc_t(i,j-1))/delta_y(j)


		!if (tau_xx(i,j,2) .gt. 0.d0) print*,tau_xx(i,j,1)
	enddo
	enddo

	xres_t = 0.d0
	yres_t = 0.d0

	do i = 2, nx + 1
	do j = 2, ny + 1
		do n = 1, 4
			xres_t(i,j) = xres_t(i,j) &
				    + tau_xx(i,j,n)*sx(i,n) &
				    + tau_xy(i,j,n)*sy(j,n)
			yres_t(i,j) = yres_t(i,j) &
				    + tau_yx(i,j,n)*sx(i,n) &
				    + tau_yy(i,j,n)*sy(j,n)  
		enddo

	enddo
	enddo

	!do j = 2, ny + 1
		!print*, j, uc_t(2,j)
		!print'(i,f10.5)',j,xres_t(2,j)
	!enddo

	!Divide diffusive term by Reynolds number
	xres_t = xres_t/Re
	yres_t = yres_t/Re
		
	
end subroutine continuum_calculate_residual_FV


