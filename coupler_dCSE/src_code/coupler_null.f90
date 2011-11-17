module coupler
        implicit none
        save 

        ! dummy or active coupler
        logical, parameter :: COUPLER_IS_ACTIVE = .false.

        ! other dummy values
        integer :: COUPLER_CFD=0, COUPLER_MD=0, COUPLER_COMM=0

contains

        subroutine coupler_create_comm(realm,ierror)
                implicit none

                integer, intent(in) :: realm ! CFD or MD ?
                integer, intent(out):: ierror


                ierror=0

        end subroutine coupler_create_comm


        subroutine coupler_get_cfd_info(imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,&
                kmax,kmaxo,nsteps,x,y,z,dx,dz,npx,npy,npz,icoord,dt)

                implicit none

                integer, intent(in) :: imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo,nsteps,&
                        npx,npy,npz,icoord(:,:)
                real(kind(0.d0)), intent(in)    :: x(:),y(:),z(:),dx,dz,dt

        end subroutine coupler_get_cfd_info


        subroutine coupler_get_md_info(npxin,npyin,npzin,icoordin,dtin)
                implicit none

                integer, intent(in)          :: npxin, npyin, npzin, icoordin(:,:)
                real(kind(0.d0)), intent(in) :: dtin

        end subroutine coupler_get_md_info


        subroutine coupler_create_map
                implicit none

        end subroutine coupler_create_map


        subroutine coupler_constrain_forces(np,pressure,r,a)
                implicit none 

                integer, intent(in)                  :: np
                real(kind=kind(0.d0)), intent(inout) :: a(:,:)
                real(kind=kind(0.d0)), intent(in)    :: pressure, r(:,:)

        end subroutine coupler_constrain_forces

        subroutine coupler_apply_continuum_forces(np,r,v,a,iter)
                implicit none

                real(kind=kind(0.d0)), dimension(:,:), intent(in) :: r,v
                real(kind=kind(0.d0)), dimension(:,:), intent(inout) :: a 
                integer, intent(in) :: np,iter  ! iteration step, it assumes that each MD average

        end subroutine coupler_apply_continuum_forces


        subroutine coupler_uc_average_test(np,r,v,lwrite)                 
                implicit none

                integer, intent(in) :: np
                real(kind=kind(0.d0)), intent(in) :: r(:,:),v(:,:)
                logical, intent(in) :: lwrite

        end subroutine coupler_uc_average_test


        subroutine coupler_boundary_cell_average(np,r,v,send_data)
                implicit none

                integer, intent(in) :: np
                real(kind=kind(0.d0)), dimension(:,:), intent(in) :: r,v
                logical, intent(in) :: send_data

        end subroutine coupler_boundary_cell_average


        subroutine coupler_send_CFDvel(uc,vc)
                implicit none

                real(kind=kind(0.d0)) uc(:,:),vc(:,:)

	end subroutine coupler_send_CFDvel


        subroutine coupler_md_vel(uc,vc,wc)
                implicit none

                real(kind=kind(0.d0)),dimension(:,:,:),intent(out) :: uc, vc, wc 

        end subroutine coupler_md_vel

        subroutine coupler_get(xL_md,yL_md,zL_md)
                implicit none

                real(kind(0.d0)), optional, intent(out) :: xL_md, yL_md, zL_md
                
        end subroutine coupler_get


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

                 integer n

                 n = 1
         end function coupler_get_nsteps


end module coupler
