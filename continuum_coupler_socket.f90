module continuum_coupler_socket
        implicit none

contains

        subroutine send_CFDvel
                use grid_arrays, only : uc, vc
                use  coupler, only : coupler_send_CFDvel 
                implicit none

!!$                test data transfer
!!$                integer i,j
!!$
!!$                do j=1,size(uc,dim=2)
!!$                 do i=1,size(uc,dim=1)
!!$                  uc(i,j) = 10*j+i
!!$                  vc(i,j) = 100*j+10*i
!!$                 enddo
!!$                enddo

                call coupler_send_CFDvel(uc,vc)
        end subroutine send_CFDvel


        subroutine MD_continuum_BC(u,v)
                use coupler, only : coupler_md_vel
                implicit none

!                integer i
!                integer, save :: ncall = 0

                real(kind(0.d0)), intent(out) :: u(:), v(:)
                real(kind(0.d0)) w(1,1,1), u3(1,size(u),1),v3(1,size(v),1)


!  a polimorfic inteface is needed here
                call coupler_md_vel(u3,v3,w)

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
end module continuum_coupler_socket

