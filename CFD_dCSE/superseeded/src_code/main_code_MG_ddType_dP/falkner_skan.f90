!-------------------------------------------------------------------
!                        f95 -r8 Falkner_Skan.f90
!-------------------------------------------------------------------
!
!                Code to compute Falkner-Skan solutions
!
!-------------------------------------------------------------------
!
!  Differences from (main_tamer/boundaries.blayer.f90)
!        1)  Does not pass (RK4_FSkan_derivs) in (call RK4_FSkan)
!        2)  Added (Beta), the Falkner-Skan angle to all subroutine calls
!        3)  Whenever (F) is used in (RK4_FSkan),
!                it is explicitely replaced by (RK4_FSkan_derivs)
!
!-------------------------------------------------------------------

!=======================================================================
! Runge-Kutta 4th order routine for solving Falkner-Skan flow
!
! RK4_FSkan_derivs(x,y,yp, Beta)
! RK4_FSkan(N,X,Yinit,K,H,F,Re,NY,ygeom,Uin,Uinpp,x0,x1, Beta)
!
!=======================================================================
subroutine RK4_FSkan_derivs(x,y,yp, Beta)
      dimension y(3),yp(3)

      yp(1) = -y(1)*y(3) - Beta*(1.-y(2)*y(2))
!      yp(1) = -y(1)*y(3)/2.
      yp(2) = y(1)
      yp(3) = y(2)

      return
end

!=======================================================================
subroutine RK4_FSkan(N,X,Yinit,K,H,Re,NY,ygeom,Uin,Uinpp,x0,x1, Beta)
        !-----THIS SUBROUTINE ADVANCES THE SOLUTION TO THE SYSTEM OF N ODE'S
        !-----DY(I)/DX = F(I,X,Y(1),Y(2),,.....Y(N)) USING THE FOURTH ORDER
        !-----RUNGE-KUTTA METHOD.  THE OTHER INPUT VARIABLES ARE:
        !-----X = INITIAL VALUE OF THE INDEPENDENT VARIABLE; IT IS INCREASED
        !-----    BY THE SUBROUTINE.
        !-Yinit = INITIAL VALUE OF THE DEPENDENT VARIABLE(S); (deriv at wall)
        !-----    THE ROUTINE COPIES (Yinit) into (Y) TO PREVENT
        !-----    OVERWRITING THIS WITH THE NEW VALUE(S).
        !-----H = Array containing the step sizes (for constant stepsize H
        !-----    should contain the same value).
        !---- ygeom = y(j) geometry information
        !---- Uin = the inlet velocity profile

        ! Need large NY (small H) to get accurate solution
        ! NY = 512 seemed to be fairly good


      DIMENSION Yinit(N),Y(N),YS(N),YSS(N),YSSS(N),T1(N),T2(N), &
                T3(N),T4(N),H(0:K),ygeom(0:K),Uin(0:NY,3),Uinpp(0:NY)
      INTEGER NY,K
      REAL    Beta, mFS                !  Beta = (2m)/(m+1)
      REAL    Uinf                !  Uinf(x) = (x/x0)^m

      ! Initialization
      Uin = 0.
      Uinpp = 0.
      mFS = Beta / (2.-Beta)
      Uinf = (x1/x0)**mFS
      Y = Yinit;

!-----THE MAIN LOOP
      DO I=1,2*K-2
        !-----Independent variable is X=eta=y*sqrt(Re*(mFS+1.)*Uinf/(2.*x1))
        !-----Step by dy/2 for 2 reasons:
        !-----  Increased stability of RK4
        !-----  Need staggered u and v so save each one every other step
        !-----Note: dy(1) = 0. so first time through RK4 gives u|wall

        ! Step size is d(eta), not dy where eta is the Falkner-Skan variable
          HH = H((I-1)/2+1) * sqrt( Re*(mFS+1.)*Uinf/(2.*x1) ) / 2.
        !-----TEMPORARY ARRAYS ARE NEEDED FOR THE FUNCTIONS TO SAVE THEM
        !-----FOR THE FINAL CORRECTOR STEP.

!-----FIRST (HALF STEP) PREDICTOR
                call RK4_FSkan_derivs(X,Y,t1, Beta)
                DO J=1,N
!                       T1(J) = F(J,X,Y)
                        YS(J) = Y(J) + .5 * HH * T1(J)
                ENDDO
!-----SECOND STEP (HALF STEP CORRECTOR)
                X = X + .5 * HH
                call RK4_FSkan_derivs(X,Ys,t2, Beta)
                DO J=1,N
!                       T2(J) = F(J,X,YS)
                        YSS(J) = Y(J) + .5 * HH * T2(J)
                ENDDO
!-----THIRD STEP (FULL STEP MIDPOINT PREDICTOR)
                call RK4_FSkan_derivs(X,Yss,t3, Beta)
                DO J=1,N
!                       T3(J) = F(J,X,YSS)
                        YSSS(J) = Y(J) + HH * T3(J)
                ENDDO
!-----FINAL STEP (SIMPSON'S RULE CORRECTOR)
                X = X + .5 * HH 
                call RK4_FSkan_derivs(X,Ysss,t4, Beta)
                DO J=1,N
!                       T4(J) = F(J,X,YSSS)
                        Y(J)=Y(J) + (HH/6.)*(T1(J) + 2.*(T2(J)+T3(J)) + T4(J))
                ENDDO
        ! f'' = Y1
        ! f'  = Y2
        ! f   = Y3
        ! I is the J index, x is eta

        ! Stability of RK4
        ! lambda*h < -2.8 is unstable where lambda is an eigenvalue of
        ! the matrix y' = Ay
!        print*,'lambda*h = ',-(Y(1)**.33333),HH,-(Y(1)**.33333)*HH
!        write(*,'(8(f13.8))') x,-(Y(1)**.33333),HH,-(Y(1)**.33333)*HH,&
!                y(1),y(2),y(3),-y(1)*y(3)

                !T I changed from [Uin(I/2+2,1) --> Uin(I/2+1,1)]
                !T       and from [Uin(I/2+1,1) --> Uin(I/2  ,1)]
                if (mod(I,2) == 1) then ! save u velocity
                        Uin(I/2+1,1) = Uinf * Y(2)
                        Uinpp(I/2+1) = Uinf * ( -Y(1)*Y(3) - Beta*(1.-Y(2)*Y(2)) )
                        if ( (Uin(I/2+1,1) < Uin(I/2  ,1)) .or. &
                             (Uin(I/2+1,1) > Uinf) ) then
                                Uin(I/2+1,1) = Uinf
                                Uinpp(I/2+1) = 0.
                                vsave = Uin(I/2,2)
                                exit
                        endif
!                       print*,x,x/sqrt(Re*(mFS+1.)*Uinf/(2.*x1)),Uin(I/2+1,1)
                else ! save v velocity (last one when I=2*K-2 isn't used)
                        Uin(I/2+1,2) = sqrt((mFS+1.)*Uinf/(2.*Re*x1)) * (x*Y(2)*(1.-mFS)/(1.+mFS) - Y(3))
!                       print*,x,x/sqrt(Re*(mFS+1.)*Uinf/(2.*x1)),Uin(I/2+1,2)
                endif
!               write(*,'(f5.2,5(f10.5))') x,y,Uin(I/2+1,1),Uin(I/2+1,2)
!                if (Y(2).ge. 0.999999) then
!                        vsave = Uin(I/2+1,2)
!                        exit  ! stepsize may get too large so stop
!                endif
      ENDDO

        ! For no mean pressure gradient and no suction/blowing
        ! u,v,u" = 0 at the wall
        Uin(0,1) = -Uin(1,1)         ! u-no-slip (gets redone later in BC routine)
        Uin(0,2) = -Uin(2,2)         ! v-no-slip (gets redone later in BC routine)
        Uin(1,2) = 0.                 ! v-no-slip (gets redone later in BC routine)
        Uinpp(1) = 0.                 ! u"-for zero mean pressure gradient flow

        if (I .le. 2*k-2) then  ! we broke out of loop
                do j=I,2*k-2
                        if (mod(j,2) == 1) then ! save u velocity
                                Uin(j/2+1,1) = Uinf
                                Uinpp(j/2+1) = 0.
!                               print*,'u: ',x,Uin(j/2+1,1)
                         else
                                Uin(j/2+1,2) = vsave
!                               print*,'v: ',x,Uin(j/2+1,2)
                         endif
                enddo
        endif
      RETURN
      END

