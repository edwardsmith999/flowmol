
module physical_constants

	double precision		:: meu, rho, Re

end module physical_constants

!-------------------------------------------------------------------------------------
module computational_constants

	integer				:: solver 		!Employ Finite difference or finite volume solver
	integer,parameter	:: FD = 1, FV = 2
	integer				:: continuum_iter, continuum_initialstep
	integer				:: continuum_tplot, continuum_Nsteps
	integer				:: nx, ny  			!Number of points
	integer				:: continuum_vflag
	integer,dimension(4)		:: BC_flag
	double precision,dimension(4)	:: uwall_BC	
	double precision		:: lx, ly, lz  			!Domain extents
	double precision		:: continuum_delta_t  		!Timestep
	double precision		:: domain_volume
	double precision		:: grid_density
	double precision,  dimension(:)  ,  allocatable :: delta_x, delta_y  !grid volume size
	double precision,  dimension(:,:),  allocatable	:: vcell
	double precision,  dimension(:,:),  allocatable	:: sx, sy

end module computational_constants

!-------------------------------------------------------------------------------------
!----------------------------------Grid based Arrays---------------------------------------------

module grid_arrays

	integer*1,        dimension(:,:)  , allocatable	:: solid	!Solid cells in domain
	double precision, dimension(:)    , allocatable	:: mx, my	!Grid locations
	double precision, dimension(:)    , allocatable	:: px, py	!Cell centres
	double precision, dimension(:,:)  , allocatable	:: uc, vc	!Velocity at cell centres
	double precision, dimension(:,:)  , allocatable	:: uc_t_minus_1, vc_t_minus_1 	!Velocity at cell centres
	double precision, dimension(:,:)  , allocatable	:: xresidual, yresidual     !Cell residual
	double precision, dimension(:,:)  , allocatable	:: xresidual_t_minus_1, yresidual_t_minus_1     !Cell residual

	!Cell fluxes - tensors have 4 components to define in 2D
	double precision, dimension(:,:)  , allocatable	:: flux_xx,flux_xy,flux_yx,flux_yy !Cell flux
	double precision, dimension(:,:,:), allocatable	:: sflux_xx,sflux_xy,sflux_yx,sflux_yy !Cell surface flux
	double precision, dimension(:,:,:), allocatable	:: tau_xx,tau_xy,tau_yx,tau_yy 	!Cell diffusive flux
	double precision, dimension(:,:,:), allocatable	:: xfluxdots, yfluxdots


end module grid_arrays
