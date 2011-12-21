module coupler
    implicit none
    save 

    ! dummy or active coupler
    logical, parameter :: COUPLER_IS_ACTIVE = .false.

    ! dummy values
    integer :: COUPLER_CFD=0, COUPLER_MD=0

contains

subroutine coupler_create_comm(realm,realm_comm,ierror)
    implicit none

    integer, intent(in) :: realm ! CFD or MD ?
    integer, intent(out):: realm_comm,ierror

	realm_comm = 0
    ierror=0

end subroutine coupler_create_comm


subroutine coupler_cfd_init(imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,&
   							 kmax,kmaxo,nsteps,x,y,z,dx,dz,npx,npy,npz,icoord,dt,density)
    implicit none

    integer, intent(in) :: imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo,nsteps,&
            npx,npy,npz,icoord(:,:)
    real(kind(0.d0)), intent(in)    :: x(:),y(:),z(:),dx,dz,dt,density

end subroutine coupler_cfd_init


subroutine coupler_md_init(npxin,npyin,npzin,icoordin,dtin)
    implicit none

    integer, intent(in)          :: npxin, npyin, npzin, icoordin(:,:)
    real(kind(0.d0)), intent(in) :: dtin

end subroutine coupler_md_init


subroutine coupler_create_map
    implicit none

end subroutine coupler_create_map


subroutine coupler_md_constrain_forces(np,pressure,r,a)
    implicit none 

    integer, intent(in)                  :: np
    real(kind=kind(0.d0)), intent(inout) :: a(:,:)
    real(kind=kind(0.d0)), intent(in)    :: pressure, r(:,:)

	a = a

end subroutine coupler_md_constrain_forces


subroutine coupler_md_apply_continuum_forces(np,r,v,a,iter)
    implicit none

    real(kind=kind(0.d0)), dimension(:,:), intent(in) :: r,v
    real(kind=kind(0.d0)), dimension(:,:), intent(inout) :: a 
    integer, intent(in) :: np,iter  ! iteration step, it assumes that each MD average

	a = a

end subroutine coupler_md_apply_continuum_forces


subroutine coupler_uc_average_test(np,r,v,lwrite)                 
    implicit none

    integer, intent(in) :: np
    real(kind=kind(0.d0)), intent(in) :: r(:,:),v(:,:)
    logical, intent(in) :: lwrite

end subroutine coupler_uc_average_test


subroutine coupler_md_boundary_cell_average(np,r,v,send_data)
    implicit none

    integer, intent(in) :: np
    real(kind=kind(0.d0)), dimension(:,:), intent(in) :: r,v
    logical, intent(in) :: send_data

end subroutine coupler_md_boundary_cell_average


subroutine coupler_cfd_send_velocity(uc,vc)
    implicit none

    real(kind=kind(0.d0)) uc(:,:),vc(:,:)

end subroutine coupler_cfd_send_velocity


subroutine coupler_cfd_get_velocity(uc,vc,wc)
    implicit none

    real(kind=kind(0.d0)),dimension(:,:,:),intent(out) :: uc, vc, wc

	uc = 0.d0
	vc = 0.d0
	wc = 0.d0

end subroutine coupler_cfd_get_velocity


subroutine coupler_md_get(xL_md,yL_md,zL_md,MD_initial_cellsize,top_dy)
    implicit none

    real(kind(0.d0)), optional, intent(out) :: xL_md, yL_md, zL_md,MD_initial_cellsize,top_dy

	xL_md               = 0.d0
	yL_md               = 0.d0
	zL_md               = 0.d0
        MD_initial_cellsize = 0.d0
        
end subroutine coupler_md_get

subroutine coupler_md_boundary_forces(np,pressure,r,a)
	implicit none 

	integer, intent(in)		  				:: np
	real(kind=kind(0.d0)), intent(inout) 	:: a(:,:)
	real(kind=kind(0.d0)), intent(in)    	:: pressure, r(:,:)

	a = a
       
end subroutine coupler_md_boundary_forces

function coupler_get_save_period() result(p)
    implicit none

    integer p
    p = 1

end function coupler_get_save_period


function coupler_get_average_period() result(p)
    implicit none

    integer p
    p = 1

end function coupler_get_average_period


function coupler_get_nsteps() result(n)
    implicit none
    integer n

    n = 1

end function coupler_get_nsteps


function coupler_md_get_dt_cfd() result(t)
    implicit none
 	real(kind=kind(0.d0)) t

	t = 0.d0

end function coupler_md_get_dt_cfd


subroutine coupler_md_set(zL_md)
	implicit none
        
	real(kind(0.d0)), optional, intent(in) :: zL_md

end subroutine coupler_md_set


function coupler_md_get_density() result(r)
	implicit none

	real(kind(0.d0)) r
	r = 0.d0

end function coupler_md_get_density

end module coupler
