!=======================================================================
! Routines for time-advancing the equations
!
! timeStep_init()
! timeStep_system (system_evaluator, dt)
! timeStep_substep (evaluator)
!
! system_evaluator(dt, iflag)
! evaluator(ipart, dt, iflag)
!
! 3rd order Runge-Kutta time-stepping scheme
!
! Coefficient    Purpose
! -----------    -------
! beta(iRK)      explicit Euler advancement
! zeta(iRK)      RK3 advancement, previous evaluation
! gamma(iRK)     RK3 advancement, current evaluation
! 0.5*beta(iRK)  trapezoidal advancement
!
! iflag = +1: first substep
! iflag = -1: last substep
! iflag =  0: middle substep
!

module RK3
	integer iRK, nRK
	real    deltaT
	real    gamma(3), zeta(3), beta(3)
end module

subroutine timeStep_init()
	use RK3

	! Constants for the RK3 time-stepping algorithm
	nRK = 3
	gamma(1) = 8./15.
	gamma(2) = 5./12.
	gamma(3) = 3./4.
	zeta(1) = 0.
	zeta(2) = -17./60.
	zeta(3) = -5./12.
	beta = gamma + zeta

	return
end

subroutine timeStep_system (system_evaluator, dt)
	use RK3

	deltaT = dt

	do iRK=1,nRK
		iflag = 0
		if (iRK > 1) iflag = iflag - 1
		if (iRK < nRK) iflag = iflag + 1

		call system_evaluator(beta(iRK)*dt, iflag)
	end do

	return
end

subroutine timeStep_substep (evaluator)
	use RK3

	dt = deltaT
	iflag = iRK / nRK  ! signals last substep

	call evaluator (1, beta(iRK)*dt, iflag)      ! Initialize
	if (iRK > 1) &
	call evaluator (2, zeta(iRK)*dt, iflag)      ! Restore part
	call evaluator (3, beta(iRK)*dt, iflag)      ! Explicit Euler part
	call evaluator (4, gamma(iRK)*dt, iflag)     ! RK3 part
	call evaluator (5, beta(iRK)*dt, iflag)      ! Trapezoidal part
	call evaluator (6, beta(iRK)*dt, iflag)      ! Finalize

	return
end

