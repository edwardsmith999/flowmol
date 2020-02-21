
module bilnear_intersect

contains

! What is the x,y,z position of a point at params u and v?
subroutine SrfEval(u, v, P, respt)
    implicit none

    real(kind(0.d0)), intent(in) :: u, v
    real(kind(0.d0)), dimension(2,2,3), intent(in) :: P
    real(kind(0.d0)), dimension(3), intent(out) ::  respt

    real(kind(0.d0)), dimension(3) :: P00, P01, P10, P11

    P00(:) = P(1,1,:); P01(:) = P(1,2,:)
    P10(:) = P(2,1,:); P11(:) = P(2,2,:)

    respt(1) = (((1.d0-u) * (1.d0-v) * P00(1) + &
                 (1.d0-u) *       v  * P01(1) + &
                       u  * (1.d0-v) * P10(1) + &
                       u  *       v  * P11(1)))
    respt(2) = (((1.d0-u) * (1.d0-v) * P00(2) + &
                 (1.d0-u) *       v  * P01(2) + &
                       u  * (1.d0-v) * P10(2) + &
                       u  *       v  * P11(2)))
    respt(3) = (((1.d0-u) * (1.d0-v) * P00(3) + &
                 (1.d0-u) *       v  * P01(3) + &
                       u  * (1.d0-v) * P10(3) + &
                       u  *       v  * P11(3)));

end subroutine SrfEval

!choose between the best denominator to avoid singularities
!and to get the most accurate root possible
subroutine getu(v, M1, M2, J1, J2, K1, K2, R1, R2, output)

    real(kind(0.d0)), intent(in) :: v, M1, M2, J1, J2, K1, K2, R1, R2
    real(kind(0.d0)), intent(out) :: output

    denom = (v*(M1-M2)+J1-J2)
    d2 = (v*M1+J1)
    if(abs(denom) > abs(d2)) then ! which denominator is bigger
        output = (v*(K2-K1)+R2-R1)/denom
    else
        output= -(v*K1+R1)/d2
    endif

end subroutine getu

! compute t with the best accuracy by using the component
! of the direction that is largest
subroutine computet(dir, orig, srfpos, output)
    implicit none

    real(kind(0.d0)), dimension(3), intent(in)  :: dir, orig, srfpos
    real(kind(0.d0)), intent(out) :: output

    ! if x is bigger than y and z
    if (abs(dir(1)) >= abs(dir(2)) .and. &
        abs(dir(1)) >= abs(dir(3))) then
        output= (srfpos(1) - orig(1)) / dir(1)
    ! if y is bigger than x and z
    else if (abs(dir(2)) >= abs(dir(3))) then ! .and. abs(dir(2)) >= abs(dir(1)))
        output = (srfpos(2) - orig(2)) / dir(2)
    ! otherwise x isn't bigger than both and y isn't bigger than both
    else  !if(abs(dir(3)) >= abs(dir(1)) .and. abs(dir(3)) >= abs(dir(2)))
        output = (srfpos(3) - orig(3)) / dir(3)
    endif
  
end subroutine computet

subroutine cross_product(a, b, cross)
    implicit none

    real(kind(0.d0)), dimension(3), intent(in)  :: a, b
    real(kind(0.d0)), dimension(3), intent(out) :: cross

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

end subroutine cross_product


subroutine surface_normal(u, v, P, normal)
    implicit none

    real(kind(0.d0)), intent(in)  :: u, v
    real(kind(0.d0)), dimension(2,2,3), intent(in) :: P
    real(kind(0.d0)), dimension(3), intent(out) :: normal

    real(kind(0.d0)), dimension(3) :: P00, P01, P10, P11
    real(kind(0.d0)), dimension(3) :: tanu, tanv

    P00(:) = P(1,1,:); P01(:) = P(1,2,:)
    P10(:) = P(2,1,:); P11(:) = P(2,2,:)

    tanu(1) = ( 1.d0 - v ) * (P10(1) - P00(1)) + v * (P11(1) - P01(1))
    tanu(2) = ( 1.d0 - v ) * (P10(2) - P00(2)) + v * (P11(2) - P01(2))
    tanu(3) = ( 1.d0 - v ) * (P10(3) - P00(3)) + v * (P11(3) - P01(3))

    tanv(1) = ( 1.d0 - u ) * (P01(1) - P00(1)) + u * (P11(1) - P10(1))
    tanv(2) = ( 1.d0 - u ) * (P01(2) - P00(2)) + u * (P11(2) - P10(2))
    tanv(3) = ( 1.d0 - u ) * (P01(3) - P00(3)) + u * (P11(3) - P10(3))

    call cross_product(tanu, tanv, normal)

end subroutine surface_normal

subroutine bicubic_line_intersect(r, q, P, flag, uvsoln, pos)
    implicit none

    real(kind(0.d0)), dimension(3), intent(in) :: r, q
    real(kind(0.d0)), dimension(2,2,3), intent(in) :: P

    integer, intent(out) :: flag
    real(kind(0.d0)), dimension(3,2), intent(out) :: uvsoln, pos

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Equation of the patch:
    ! P(u, v) = (1-u)(1-v)P00 + (1-u)vP01 + u(1-v)P10 + uvP11
    ! Equation of the line:
    ! R(t) = r + tq
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer :: i, num_sol, code ! number of solutions to the quadratic
    real(kind(0.d0)) :: ax, ay, az, bx, by, bz, cx, cy, cz
    real(kind(0.d0)) :: rx, ry, rz, qx, qy, qz
    real(kind(0.d0)) :: dx, dy,dz, A, A1, A2, B, B1, B2, C, C1, C2, D, D1, D2
    real(kind(0.d0)) :: t2, u  ! the t values of the two roots
    real(kind(0.d0)), parameter :: ray_epsilon=1e-12
    real(kind(0.d0)), dimension(2) :: vsol ! the two roots from quadraticroot
    real(kind(0.d0)), dimension(3) :: pos1, pos2, dir, orig !Vector pos = ro + t*rd
    real(kind(0.d0)), dimension(3) :: P00, P01, P10, P11, uv
    complex(kind(0.d0)), dimension(2) :: z

    P00(:) = P(1,1,:); P01(:) = P(1,2,:)
    P10(:) = P(2,1,:); P11(:) = P(2,2,:)

    !Set small numbers to zero, not sure why it's needed 
!    do i=1,3
!        if (P00(i) .lt. 1e-12) P00(i) = 0.d0
!        if (P10(i) .lt. 1e-12) P10(i) = 0.d0
!        if (P01(i) .lt. 1e-12) P01(i) = 0.d0
!        if (P11(i) .lt. 1e-12) P11(i) = 0.d0
!    enddo

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Variables for substitition
    ! a = P11 - P10 - P01 + P00
    ! b = P10 - P00
    ! c = P01 - P00
    ! d = P00 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    ! Find a w.r.t. x, y, z
    ax = P11(1) - P10(1) - P01(1) + P00(1)
    ay = P11(2) - P10(2) - P01(2) + P00(2)
    az = P11(3) - P10(3) - P01(3) + P00(3)

    ! Find b w.r.t. x, y, z
    bx = P10(1) - P00(1)
    by = P10(2) - P00(2)
    bz = P10(3) - P00(3)

    ! Find c w.r.t. x, y, z
    cx = P01(1) - P00(1)
    cy = P01(2) - P00(2)
    cz = P01(3) - P00(3)

    rx = r(1)
    ry = r(2)
    rz = r(3)

    ! Retrieve the xyz of the q part of ray
    qx = q(1)
    qy = q(2)
    qz = q(3)

    ! Find d w.r.t. x, y, z - subtracting r just after  
    dx = P00(1) - r(1)
    dy = P00(2) - r(2)
    dz = P00(3) - r(3)

    ! Find A1 and A2
    A1 = ax*qz - az*qx
    A2 = ay*qz - az*qy

    ! Find B1 and B2
    B1 = bx*qz - bz*qx
    B2 = by*qz - bz*qy

    ! Find C1 and C2
    C1 = cx*qz - cz*qx
    C2 = cy*qz - cz*qy

    ! Find D1 and D2
    D1 = dx*qz - dz*qx
    D2 = dy*qz - dz*qy

    dir = q
    orig = r
    A = A2*C1 - A1*C2
    B = A2*D1 -A1*D2 + B2*C1 -B1*C2
    C = B2*D1 - B1*D2
  
    uv(1) = -2
    uv(2) = -2
    uv(3) = -2

    !Use numerical roots

    call QuadraticRoot_intersect(A, B, C, -ray_epsilon, 1+ray_epsilon, vsol, num_sol)

    !print*, A, B, C, vsol, num_sol
    !Set default null values for position and uv
    uvsoln = -666
    pos = -666

    select case(num_sol)
    case(0)
        flag = 0 ! no solutions found
    case(1)
        uv(2) = vsol(1)
        call getu(vsol(1), A2, A1, B2, B1, C2, C1, D2, D1, uv(1))
        if ((uv(1) < 1+ray_epsilon) .and. & 
            (uv(1) > -ray_epsilon)) then
            call SrfEval(uv(1), uv(2), P, pos1)
            call computet(dir, orig, pos1, uv(3))
            if (uv(3) > 0.d0) then
                !print*, "ONE ROOT, uv, pos1", uv, pos1
                flag = 1
                uvsoln(:,1) = uv
                pos(:,1) = pos1
            else
                flag = 0
            endif
        else
	        flag = 0 ! no other soln - so ret false
        endif 
    case(2) ! two solutions found
        uv(2) = vsol(1)
        call getu(uv(2), A2, A1, B2, B1, C2, C1, D2, D1, uv(1))
        call SrfEval(uv(1),uv(2), P, pos1)
        call computet(dir, orig, pos1, uv(3))
        !print*, "uv, pos1", uv, pos1
        if ((uv(1) < 1+ray_epsilon) .and. &
            (uv(1) > -ray_epsilon) .and. & 
            (uv(3) > 0.d0)) then
            call getu(vsol(2), A2, A1, B2, B1, C2, C1, D2, D1, u)
            flag = 1 ! If u2 is bad, u1 vars are still okay
            uvsoln(:,1) = uv
            pos(:,1) = pos1
            if ((u < 1+ray_epsilon) .and. (u > ray_epsilon)) then
                call SrfEval(u, vsol(2), P, pos2)
                call computet(dir, orig, pos2, t2)
                print*, "TWO ROOTS, uv, pos1, pos2", uv, pos1, pos2, t2
                if ((t2 < 0) .or. (uv(3) < t2)) then
                    flag = 1
                else
                    ! other wise both t2 > 0 and t2 < t1
                    uv(2) = vsol(2)
                    uv(1) = u
                    uv(3) = t2
                    uvsoln(:,2) = uv
                    pos(:,2) = pos2
                    flag = 2
                endif
            endif
        else ! doesn't fit in the root - try other one
            uv(2) = vsol(2)
            call getu(vsol(2), A2, A1, B2, B1, C2, C1, D2, D1, uv(1)) 
            call SrfEval(uv(1), uv(2), P, pos1)
            call computet(dir, orig, pos1, uv(3)) 
            !print*, "NEXT uv, pos1", uv, pos1
            if ((uv(1) < 1+ray_epsilon) .and. (uv(1) > -ray_epsilon) .and. (uv(3) > 0.d0)) then
                flag = 1
                uvsoln(:,1) = uv
                pos(:,1) = pos2
            else
                flag = 0
            endif
        endif
    end select

end subroutine bicubic_line_intersect


! a x ^2 + b x + c = 0
! in this case, the root must be between min and max
! it returns the # of solutions found
! x = [ -b +/- sqrt(b*b - 4 *a*c) ] / 2a
! or x = 2c / [-b +/- sqrt(b*b-4*a*c)]
subroutine QuadraticRoot_intersect(a, b, c, mn, mx, u, num_sol)
    implicit none

    real(kind(0.d0)), intent(in) :: a, b, c, mn, mx

    integer, intent(out) :: num_sol
    real(kind(0.d0)), dimension(2), intent(out) :: u

    real(kind(0.d0)) :: q, d, dummy

    u(1) = mn-mn ! make it lower than min
    u(2) = u(1)
    if(a == 0.0) then ! then its close to 0
        if(b .ne. 0.0) then ! not close to 0
            u(1) = - c / b
            if(u(1) > mn .and. u(1) < mx) then !its in the interval
                num_sol = 1; return  !1 soln found
            else  !its not in the interval
                num_sol = 0; return 
            endif
        else
            num_sol = 0; return 
        endif
    endif

    d = b*b - 4.d0*a*c !discriminant
    if (d <= 0.0) then ! single or no root
        if (d == 0.0) then ! close to 0
            u(1) = -b / a
            if (u(1) > mn .and. u(1) < mx) then ! its in the interval
                num_sol = 1; return 
            else ! its not in the interval
                num_sol = 0; return 
            endif
        else ! no root d must be below 0
            num_sol = 0; return 
        endif
    endif

    q = -0.5d0  * (b + sign(sqrt(d),b))
    u(1) = c / q
    u(2) = q / a

    if((u(1) > mn .and. u(1) < mx) .and. &
       (u(2) > mn .and. u(2) < mx)) then
        num_sol = 2; return 
    else if(u(1) > mn .and. u(1) < mx) then !then one wasn't in interval
        num_sol = 1; return 
    else if(u(2) > mn .and. u(2) < mx) then
      ! make it easier, make u(1) be the valid one always
      dummy = u(1)
      u(1) = u(2)
      u(2) = dummy ! just in case somebody wants to check it
      num_sol = 1; return 
    endif
    num_sol = 0

end subroutine QuadraticRoot_intersect


subroutine line_plane_intersect(ri, rij, P, intersect, normal, flag)
    implicit none


    real(kind(0.d0)), dimension(3), intent(in) :: ri, rij
    real(kind(0.d0)), dimension(2,2,3), intent(in) :: P

    integer, intent(out) :: flag
    real(kind(0.d0)), dimension(3,2), intent(out) :: intersect, normal

    real(kind(0.d0)), dimension(3,2) :: uv

    call bicubic_line_intersect(ri, rij, P, flag, uv, intersect)
    if (flag .eq. 1) then
        !call SrfEval(uv(1), uv(2), P, intersect)
        call surface_normal(uv(1,1), uv(2,1), P, normal(:,1))
    else if (flag .eq. 2) then
        call surface_normal(uv(1,1), uv(2,1), P, normal(:,1))
        call surface_normal(uv(1,2), uv(2,2), P, normal(:,2))
    else
        intersect = -666
        normal = -666
    endif

end subroutine line_plane_intersect


end module bilnear_intersect

    !program test
    !    use line_plane_intersect
    !    implicit none

    !    integer :: unitno=5, flag
    !    real(kind(0.d0)), dimension(3) :: r, q 
    !    real(kind(0.d0)), dimension(3,2) :: uv, intersect, normal
    !    real(kind(0.d0)), dimension(2,2,3) :: P

    !    ! Create 4 points
    !!    P(1,1,:) = (/0.d0, 0.d0, 0.d0 /)
    !!    P(1,2,:) = (/3.d0, 1.d0, 3.d0 /)
    !!    P(2,1,:) = (/1.d0, 3.d0, 1.d0 /)
    !!    P(2,2,:) = (/1.d0,-2.d0, 4.d0 /)

    !    ! Some ray information
    !!    r = (/ 1.d0, 0.3d0, 10.d0 /) !origin of the ray
    !!    q = (/ 0.100499d0, 0.d0, -0.994937d0 /) ! a ray direction

    !	open (unitno,file="./input")
    !	read (unitno,*) P(1,1,:)
    !	read (unitno,*) P(1,2,:)
    !	read (unitno,*) P(2,1,:)
    !	read (unitno,*) P(2,2,:)
    !	read (unitno,*) r
    !	read (unitno,*) q
    !	close(unitno)

    !    print*, P, r, q

    !!    P(1,1,:) = (/0.d0, 0.d0, 0.0d0 /)
    !!    P(1,2,:) = (/0.d0, 1.d0, 1.0d0 /)
    !!    P(2,1,:) = (/1.d0, 0.d0, 1.0d0 /)
    !!    P(2,2,:) = (/1.d0, 1.d0, 0.0d0 /)

    !!    r = (/0.d0, -0.5d0, 0.5d0 /)
    !!    q = (/0.0d0, 1.0d0, 0.5d0 /)

    !    call line_plane_intersect(r, q, P, intersect, normal, flag)
    !    !call bicubic_line_intersect(r, q, P, uv, flag)
    !    !call SrfEval(uv(1), uv(2), P, intersect)
    !    !call surface_normal(uv(1), uv(2), P, normal)

    !    print'(l, 2(a,6f14.5))', flag, " Intersect = ", intersect, " Normal = ", normal

    !end program



