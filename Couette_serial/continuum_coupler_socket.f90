module continuum_coupler_socket
        implicit none
#if USE_COUPLER
contains

        subroutine send_CFDvel
                use grid_arrays, only : uc, vc
                use  coupler, only : coupler_cfd_send_velocity
                implicit none

 		integer i
		real(kind(0.d0)) uc3d(1,size(uc,dim=1),size(uc,dim=2)), & !hack for 2d parallelism, z dimension independent of nx
                                 vc3d(1,size(vc,dim=1),size(vc,dim=2))
 
                do i=1,size(uc3d,dim=1)
                        uc3d(i,:,:) = uc (:,:)
                enddo
                do i=1,size(vc3d,dim=1)
                        vc3d(i,:,:) = vc(:,:)
                enddo

                call coupler_cfd_send_velocity(uc3d,vc3d)

        end subroutine send_CFDvel


        subroutine MD_continuum_BC(u,v)
                use coupler, only : coupler_cfd_get_velocity
                implicit none

!                integer i
!                integer, save :: ncall = 0

                real(kind(0.d0)), intent(out) :: u(:), v(:)
                real(kind(0.d0)) w(1,1,1), u3(1,size(u),1),v3(1,size(v),1)


		!A polymorphic interface is needed here
                call coupler_cfd_get_velocity(u3,v3,w)

                u(:) =  u3(1,:,1)
                v(:) =  v3(1,:,1)
!debug 
!                 u(:) = 0.d0
!                 v(:) = 0.d0

!                ncall = ncall + 1
!                do i = 1,size(u)
!                        write(7000+ncall,*), u(i),v(i)
!                enddo

        end subroutine MD_continuum_BC
#endif
end module continuum_coupler_socket

