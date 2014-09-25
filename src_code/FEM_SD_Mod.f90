!-----------------------------------------------------------------------
!	MODULES  FOR  FEM_1D.f90 CODE     
!-----------------------------------------------------------------------

module elements_data
	
	SAVE

	integer, parameter :: NBF_Q = 3
	integer, parameter :: NBF_L = 2

	integer, parameter :: NEQ_Q = 7
     
	integer, parameter :: NEQ_Q_a = 7
	integer, parameter :: NEQ_Q_b = 1
     
	integer, parameter :: NORE = 2

	integer :: NXELa, NNTOLa, NNXa, NODTOLa_QEL, NODTOLa
	integer :: NXELb, NNTOLb, NNXb, NODTOLb_QEL, NODTOLb
	integer :: NXEL,  NNTOL,  NNX,  NODTOL_QEL,  NODTOL

	integer :: IBAND, HBAND, MBAND, IDIM_A

	integer, parameter :: NRHS = 1

end module elements_data

!----------------------------------------------------------------------

module common_arrays
	
	use elements_data, only :  NODTOL, NNX, IBAND, HBAND, NEQ_Q

	save
	
	REAL(8), ALLOCATABLE, DIMENSION(:,:) :: A
	REAL(8), ALLOCATABLE, DIMENSION(:)   :: B, S
	REAL(8), ALLOCATABLE, DIMENSION(:)   :: S_c

	REAL(8) :: A_p, B_p, S_p

	REAL(8), ALLOCATABLE, DIMENSION(:) :: A_ip, A_pi

	INTEGER, ALLOCATABLE, DIMENSION(:) :: IPVT

	real(8), allocatable, dimension(:,:) :: TQa, TQao, TQao1, TQap
	real(8), allocatable, dimension(:,:) :: TQb, TQbo, TQbo1, TQbp
	real(8), allocatable, dimension(:,:) :: TQo_GRAPH

	real(8) :: Xc, Xco, Xco1, Xcp

	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NM_Q_a, NM_QQ_a
	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NM_Q_b, NM_QQ_b

	real(8), allocatable, dimension(:) :: Xa, Xb
						
end module common_arrays

!----------------------------------------------------------------------

module flow_parameters

	SAVE

	real(8) :: eps
	real(8) :: eps2
	real(8) :: As1
	real(8) :: As2
	real(8) :: As12
	real(8) :: d1
	real(8) :: d12
	real(8) :: bslip
	real(8) :: kappa
	integer :: mm
	real(8) :: Thetac
	real(8) :: MASS_SURFACTANT
	real(8) :: Peca
	real(8) :: ba
	real(8) :: ka
	real(8) :: Ra
	LOGICAL :: SOLUBLE_SURFACTANT
	real(8) :: Peb
	real(8) :: Pem
	real(8) :: kb
	real(8) :: nc
	LOGICAL :: WALL_ADSORPTION
	real(8) :: Pecs
	real(8) :: bs
	real(8) :: ks
	real(8) :: Rs 
	real(8) :: kas
	real(8) :: Ras 

	CHARACTER(LEN=5)  :: Equilibrium_Angle
	CHARACTER(LEN=11) :: GEOMETRY
	CHARACTER(LEN=8)  :: ST_MODEL

	real(8) :: Bl1
	real(8) :: Bl2
	real(8) :: Bl12
	real(8) :: Kl1
	real(8) :: Kl2
	real(8) :: Kl12

	real(8) :: Hmaxi = 1.d0

	real(8) :: mi, ci, Cai, Csi
	real(8) :: Thetaa

	REAL(8) :: MASS_FLUIDi, MASS_SURFi

end module flow_parameters

!----------------------------------------------------------------------

module basic_parameters 

	SAVE

	integer, parameter :: IDECI = 1
	integer, parameter :: NITER = 50
	real(8), parameter :: ERROR = 0.5D-8

	LOGICAL, PARAMETER :: FULL_NR = .TRUE.
	CHARACTER(LEN=3) :: FLAG_NR
	LOGICAL :: FORM_JAC

	LOGICAL :: SUPG = .FALSE.

	real(8) :: PO_MIN, PO, PO_MAX, PO_STEP
	real(8) :: PI_MIN, PI, PI_MAX, PI_STEP

	integer :: OUTFE_STEP

	CHARACTER(LEN=1) :: XPACK

	real(8), parameter :: L1 = 1.D0
	real(8), parameter :: L2 = 1.D0

	real(8) :: MPa, MPb

end module basic_parameters

!----------------------------------------------------------------------

module time_integration

	SAVE

	real(8) :: Time
	real(8) :: Time_Limit = 1.0D+4
	real(8) :: Dt, Dtn
	real(8) :: DT_max     = 1.0D-2
	real(8) :: DT_min     = 1.0D-7
	real(8) :: eps_t	= 1.0D-4

end module time_integration

!----------------------------------------------------------------------

module gauss_data

	PRIVATE
	PUBLIC NINT_3P,  GAPT_3P,  WO_3P,  GAUSS_3P,  &
			GAUSS_LINE

	SAVE

	integer, parameter				  :: NINT_3P  = 3
	real(8), dimension ( NINT_3P )    :: GAPT_3P
	real(8), dimension ( NINT_3P )    :: WO_3P

contains

	subroutine GAUSS_LINE
		implicit none

			GAPT_3P(1) = -DSQRT(3.D0/5.D0)
			GAPT_3P(2) =  0.D0
			GAPT_3P(3) = +DSQRT(3.D0/5.D0)

			WO_3P(1)   =  5.D0/9.D0
			WO_3P(2)   =  8.D0/9.D0
			WO_3P(3)   =  5.D0/9.D0

	end subroutine GAUSS_LINE

end module gauss_data

!----------------------------------------------------------------------

module basis_function_1D

	use elements_data, only  : NBF_L,    NBF_Q
	use gauss_data,    only  : GAPT_3P,  NINT_3P

	PRIVATE
	PUBLIC  TF_L, DTF_L, TF_Q, DTF_Q, NF_L, DNF_L, NF_Q, DNF_Q, &
			F_LINEAR, F_QUADRATIC, INITIALIZE_BASIS_FUNCTION_1D

		SAVE

			real(8), dimension (NBF_L,NINT_3P) :: TF_L, DTF_L
			real(8), dimension (NBF_Q,NINT_3P) :: TF_Q, DTF_Q
			real(8), dimension (NBF_L,NBF_Q)   :: NF_L, DNF_L
			real(8), dimension (NBF_Q,NBF_Q)   :: NF_Q, DNF_Q

contains

	subroutine INITIALIZE_BASIS_FUNCTION_1D

		IMPLICIT NONE

		integer :: I
		real(8), dimension(NBF_Q) :: KSI
		real(8), dimension(NBF_Q) :: F1N, DF1N

		DO I = 1, NINT_3P

			CALL  F_LINEAR( F1N(1:NBF_L), DF1N(1:NBF_L), NBF_L, GAPT_3P(I) )
							TF_L(1:NBF_L,I) =  F1N(1:NBF_L)
							DTF_L(1:NBF_L,I) = DF1N(1:NBF_L)

			CALL F_QUADRATIC( F1N(1:NBF_Q), DF1N(1:NBF_Q), NBF_Q, GAPT_3P(I) )
								TF_Q(1:NBF_Q,I) =  F1N(1:NBF_Q)
								DTF_Q(1:NBF_Q,I) = DF1N(1:NBF_Q)

		ENDDO

		KSI(1) = -1.D0
		KSI(2) =  0.D0
		KSI(3) = +1.D0

		DO I = 1, NBF_Q

			CALL  F_LINEAR( F1N(1:NBF_L), DF1N(1:NBF_L), NBF_L, KSI(I) )
			NF_L(1:NBF_L,I) =  F1N(1:NBF_L)
			DNF_L(1:NBF_L,I) = DF1N(1:NBF_L)

			CALL F_QUADRATIC( F1N(1:NBF_Q), DF1N(1:NBF_Q), NBF_Q, KSI(I) )
			NF_Q(1:NBF_Q,I) =  F1N(1:NBF_Q)
			DNF_Q(1:NBF_Q,I) = DF1N(1:NBF_Q)

		ENDDO

	end subroutine INITIALIZE_BASIS_FUNCTION_1D

	subroutine F_LINEAR( F1N, DF1N, IDIM, x )

		IMPLICIT NONE

		integer, intent(in) :: IDIM
		real(8), intent(in) :: x
		real(8), intent(inout), dimension(IDIM) :: F1N
		real(8), intent(inout), dimension(IDIM) :: DF1N

		F1N(1)  =   0.5D0*(1.D0-x)
		F1N(2)  =   0.5D0*(1.D0+x)

		DF1N(1) = - 0.5D0
		DF1N(2) =   0.5D0

	end subroutine F_LINEAR

	subroutine F_QUADRATIC( F3N, DF3N, IDIM, x )

		IMPLICIT NONE

		integer, intent(in) :: IDIM
		real(8), intent(in) :: x
		real(8), intent(inout), dimension(IDIM) :: F3N
		real(8), intent(inout), dimension(IDIM) :: DF3N

		F3N(1)  = -0.5D0*x*(1.D0-x)
		F3N(2)  = (1.D0-x)*(1.D0+x)
		F3N(3)  = +0.5D0*x*(1.D0+x)

		DF3N(1) = + x - 0.5D0
		DF3N(2) = - 2.D0*x
		DF3N(3) = + x + 0.5D0

	end subroutine F_QUADRATIC

end module basis_function_1D

!----------------------------------------------------------------------

module MATRIX_STORAGE_SUBROUTINES

contains

	!----------------------------------------------------------------------
	!					SUBROUTINE MATRIX_STORAGE_RESIDUAL
	!
	!----------------------------------------------------------------------

	SUBROUTINE MATRIX_STORAGE_RESIDUAL( TEMP, NM, IDIM, INEQ, B, IDIM_B )

		!-----------------------------------------------------------------------

		IMPLICIT NONE

		INTEGER					:: I, IEQ, IROW

		INTEGER, INTENT(IN)   :: IDIM, INEQ, IDIM_B

		INTEGER, INTENT(IN),    DIMENSION(IDIM)	:: NM
		REAL(8), INTENT(IN),    DIMENSION(IDIM,INEQ) :: TEMP
		REAL(8), INTENT(INOUT), DIMENSION(IDIM_B)    :: B

		!-----------------------------------------------------------------------
		!			STORE THE RESIDUAL VECTOR IN THE GLOBAL VECTOR B
		!-----------------------------------------------------------------------

		DO I=1,IDIM

			!-----------------------------------------------------------------------
			!			ITERATE OVER IDIM NODES
			!-----------------------------------------------------------------------

			DO IEQ=1,INEQ

				!-----------------------------------------------------------------------
				!			ITERATE OVER INEQ EQUATIONS
				!-----------------------------------------------------------------------

				IROW  = NM(I) + IEQ - 1
				B(IROW) = B(IROW) + TEMP(I,IEQ)

			ENDDO
		ENDDO

	END SUBROUTINE MATRIX_STORAGE_RESIDUAL

	!----------------------------------------------------------------------
	!					SUBROUTINE MATRIX_STORAGE_JACOBIAN_BAND
	!----------------------------------------------------------------------

	SUBROUTINE MATRIX_STORAGE_JACOBIAN_BAND &
		( TEMP_L, NM, NM_L, IDIM, KDIM, INEQ, INEQ_L, &
				A, JDIM_A, IDIM_A, IBAND )

		IMPLICIT NONE

		INTEGER			:: I, J, K, IEQ
		INTEGER			:: IROW, ICOL, JCOL

		INTEGER, INTENT(IN) :: IDIM,  KDIM
		INTEGER, INTENT(IN) :: INEQ,  INEQ_L
		INTEGER, INTENT(IN) :: IBAND, JDIM_A, IDIM_A

		INTEGER, INTENT(IN),    DIMENSION(IDIM)	:: NM
		INTEGER, INTENT(IN),    DIMENSION(KDIM)	:: NM_L
		REAL(8), INTENT(IN),    DIMENSION(IDIM,KDIM,INEQ_L,INEQ)  :: TEMP_L
		REAL(8), INTENT(INOUT), DIMENSION(JDIM_A,IDIM_A)					:: A

		!-----------------------------------------------------------------------
		!			STORE THE ELEMENT INTEGRATION MATRIX AND RESIDUAL VECTOR
		!			IN THE GLOBAL MATRIX A AND VECTOR B
		!-----------------------------------------------------------------------

		DO I = 1, IDIM

			!-----------------------------------------------------------------------
			!			ITERATE OVER IDIM NODES
			!-----------------------------------------------------------------------

			DO IEQ = 1, INEQ

				!-----------------------------------------------------------------------
				!			ITERATE OVER INEQ EQUATIONS
				!-----------------------------------------------------------------------

				IROW  = NM(I) + IEQ - 1

				DO J = 1, KDIM

					!-----------------------------------------------------------------------
					!			ITERATE OVER KDIM BASIS FUNCTIONS
					!-----------------------------------------------------------------------

					DO K = 1, INEQ_L

						!-----------------------------------------------------------------------
						!			ITERATE OVER ALL EQUATIONS
						!-----------------------------------------------------------------------

						ICOL = NM_L(J) + K - 1
						JCOL = ICOL - IROW + IBAND
						A(JCOL,IROW) = A(JCOL,IROW) + TEMP_L(I,J,K,IEQ)

					ENDDO
				ENDDO

			ENDDO
		ENDDO

END SUBROUTINE MATRIX_STORAGE_JACOBIAN_BAND

end module MATRIX_STORAGE_SUBROUTINES

!-----------------------------------------------------------------------

module MACHINE_EPSILON

	PRIVATE
	PUBLIC EPSILON, MACHINE_ERROR

	REAL(8) :: EPSILON

contains

SUBROUTINE MACHINE_ERROR

	IMPLICIT NONE

	REAL(8) :: SREL, TEMP

!-----------------------------------------------------------------------
!     INITIALIZE SREL
!-----------------------------------------------------------------------

	SREL = 1.D0

!-----------------------------------------------------------------------
!    ITERATIONS
!-----------------------------------------------------------------------

	LOOP_UNLIMITED: DO

	SREL = 0.5D0*SREL
	TEMP = SREL + 1.D0

	IF(TEMP .LE. 1.D0)EXIT LOOP_UNLIMITED

	ENDDO LOOP_UNLIMITED

!-----------------------------------------------------------------------
!     DEFINE EPSILON
!-----------------------------------------------------------------------

	EPSILON = SREL*1

	WRITE(42,*)
	WRITE(42,*)'EPSILON PARAMETER =', EPSILON
	WRITE(42,*)

END SUBROUTINE MACHINE_ERROR

end module MACHINE_EPSILON

!-----------------------------------------------------------------------

module functions

	USE FLOW_PARAMETERS

contains

	REAL(8) FUNCTION Tha(Ca,Cs)
		IMPLICIT NONE
		REAL(8) :: st1, st2, st12, Ca, Cs, SR

		st1  = Surface_Tension('SA',Cs) 
		st2  = Surface_Tension('LA',Ca) 
		st12 = Surface_Tension('LS',Cs) 

		SELECT CASE(Equilibrium_Angle)
			CASE('FIXED')
				Tha = Thetac
			CASE('DYNAM')
				SR = (As2+1.D0)*st1/(1.D0+As2*st2)*Thetac**2 + &
				2.D0/(eps**2*(1.D0+As2*st2))*( (1.D0-d1+d12)*(1.D0-st1) + &
				As2*(st2-st1) + As12*d12*(st12-st1) )
				IF ( SR > 0.D0 ) THEN
					Tha = SQRT( SR )
				ELSE
					Tha = 0.D0
				ENDIF
		CASE DEFAULT
			WRITE(42,*) 'WRONG VALUE OF Equilibrium_Angle! / FUNCTION Tha' 
			STOP
		END SELECT

		RETURN
	END FUNCTION Tha

	REAL(8) FUNCTION Surface_Tension(Interface,g)
		IMPLICIT NONE
		REAL(8) :: g, gm, Bl, As, Kl
		CHARACTER(LEN=2) :: Interface

		SELECT CASE(ST_MODEL)
		CASE('LINEAR')

			Surface_Tension = 1.d0 - g

		CASE('LANGMUIR')

			SELECT CASE(Interface)
			CASE ('SA')
				Bl = Bl1
				Kl = Kl1
			CASE ('LA')
				Bl = Bl2
				Kl = Kl2
			CASE ('LS')
				Bl = Bl12
				Kl = Kl12
			CASE DEFAULT
				WRITE(42,*) 'WRONG INTERFACE! / FUNCTION Surface_Tension' 
				STOP
			END SELECT

			gm = 1.d0 - DEXP(-1.D0/Bl)

			IF ( g < gm ) THEN
				Surface_Tension = 1.d0 + Bl*( Log(1.d0-g) + 0.5D0*Kl*g**2 )
			ELSE
				Surface_Tension = 0.D0
			ENDIF 

		CASE('SHELUDKO')

			SELECT CASE(Interface)
			CASE ('SA')
				As = As1
			CASE ('LA')
				As = As2
			CASE ('LS')
				As = As12
			CASE DEFAULT
				WRITE(42,*) 'WRONG INTERFACE! / FUNCTION Surface_Tension' 
				STOP
			END SELECT

			Surface_Tension = - 1.D0/As + (As+1.d0) / &
					( As*(1.d0+((As+1.d0)**(1.d0/3.d0)-1.d0)*g)**3 )

		CASE DEFAULT
			WRITE(42,*) 'WRONG SELECTION OF MODEL! / FUNCTION Surface_Tension' 
			STOP
		END SELECT

		RETURN
	END FUNCTION Surface_Tension

	REAL(8) FUNCTION Surface_Tension_x(Interface,g,gx)
		IMPLICIT NONE
		REAL(8) :: g, gx, gm, As, Bl, Kl
		CHARACTER(LEN=2) :: Interface

		SELECT CASE(ST_MODEL)
		CASE('LINEAR')

			Surface_Tension_x = - gx

		CASE('LANGMUIR')

			SELECT CASE(Interface)
			CASE ('SA')
				Bl = Bl1
				Kl = Kl1
			CASE ('LA')
				Bl = Bl2
				Kl = Kl2
			CASE ('LS')
				Bl = Bl12
				Kl = Kl12
			CASE DEFAULT
				WRITE(42,*) 'WRONG INTERFACE! / FUNCTION Surface_Tension_x' 
				STOP
			END SELECT

			gm = 1.d0 - DEXP(-1.D0/Bl)
			IF ( g < gm ) THEN
				Surface_Tension_x = Bl*( - gx/(1.d0-g) + Kl*g*gx )
			ELSE
				Surface_Tension_x = 0.D0
			ENDIF 

		CASE('SHELUDKO')

			SELECT CASE(Interface)
			CASE ('SA')
				As = As1
			CASE ('LA')
				As = As2
			CASE ('LS')
				As = As12
			CASE DEFAULT
				WRITE(42,*) 'WRONG INTERFACE! / FUNCTION Surface_Tension_x' 
				STOP
			END SELECT

			Surface_Tension_x = &
								-3.d0*(As+1.d0)/As*((As+1.d0)**(1.d0/3.d0)-1.d0)*gx / & 
									(1.d0+((As+1.d0)**(1.d0/3.d0)-1.d0)*g)**4

		CASE DEFAULT
			WRITE(42,*) 'WRONG SELECTION OF MODEL! / FUNCTION Surface_Tension_x' 
			STOP
		END SELECT

		RETURN
	END FUNCTION Surface_Tension_x

	REAL(8) FUNCTION Fm(m)
		IMPLICIT NONE
		REAL(8) :: m
		SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
				Fm = MASS_SURFACTANT - MASS_FLUIDi*(m**(1.d0/nc)+m) - &
				ba*Ra*m**(1.d0/nc)/(Ra*m**(1.d0/nc)+1.d0)
			CASE('CYLINDRICAL')
				Fm = MASS_SURFACTANT - MASS_FLUIDi*(m**(1.d0/nc)+m) - &
				0.5d0*ba*Ra*m**(1.d0/nc)/(Ra*m**(1.d0/nc)+1.d0)
		END SELECT  
		RETURN
		END FUNCTION Fm

		REAL(8) FUNCTION Fm1(m)
		IMPLICIT NONE
		REAL(8) :: m
		SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
					Fm1 = - MASS_FLUIDi*(1.d0/nc*m**((1.d0-nc)/nc)+1.d0) - &
						ba*Ra*1.d0/nc*m**((1.d0-nc)/nc)/(Ra*m**(1.d0/nc)+1.d0)**2
			CASE('CYLINDRICAL')
					Fm1 = - MASS_FLUIDi*(1.d0/nc*m**((1.d0-nc)/nc)+1.d0) - &
								0.5d0*ba*Ra*(1.d0/nc)*m**((1.d0-nc)/nc) / &
								(Ra*m**(1.d0/nc)+1.d0)**2
		END SELECT  
		RETURN
	END FUNCTION Fm1

	REAL(8) FUNCTION Fc(c)
		IMPLICIT NONE
		REAL(8) :: c
		SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
				Fc = MASS_SURFACTANT - MASS_FLUIDi*c - ba*Ra*c/(Ra*c+1.d0)
			CASE('CYLINDRICAL')
				Fc = MASS_SURFACTANT - MASS_FLUIDi*c - 0.5d0*ba*Ra*c/(Ra*c+1.d0)
		END SELECT  
		RETURN
		END FUNCTION Fc

		REAL(8) FUNCTION Fc1(c)
		IMPLICIT NONE
		REAL(8) :: c
		SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
				Fc1 = - MASS_FLUIDi - ba*Ra/(Ra*c+1.d0)**2
			CASE('CYLINDRICAL')
				Fc1 = - MASS_FLUIDi - 0.5d0*ba*Ra/(Ra*c+1.d0)**2
		END SELECT  
		RETURN
	END FUNCTION Fc1

	REAL(8) FUNCTION Fg(g)
		IMPLICIT NONE
		REAL(8) :: g
		SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
				Fg = MASS_SURFACTANT - g
			CASE('CYLINDRICAL')
				Fg = MASS_SURFACTANT - 0.5d0*g
		END SELECT  
		RETURN
	END FUNCTION Fg

	REAL(8) FUNCTION Fg1(g)
		IMPLICIT NONE
		REAL(8) :: g
		SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
				Fg1 = - 1.d0
			CASE('CYLINDRICAL')
				Fg1 = - 0.5d0
		END SELECT  
		RETURN
	END FUNCTION Fg1

    !------------------------------------------------------------------------------
    !Pure fortran subroutine to return an updated filename by appending
    !the current timestep to that file
    subroutine get_Timestep_FileName(timestep,basename,filename)
		implicit none

		integer,intent(in) 			:: timestep
		character(*),intent(in) 	:: basename
		character(*),intent(out)	:: filename

		if(timestep.le.9							) &
		write(filename,'(a,a7,i1)') trim(basename),'.000000',timestep
		if(timestep.ge.10	.and. timestep.le.99     ) &
		write(filename,'(a,a6,i2)') trim(basename),'.00000' ,timestep
		if(timestep.ge.100     .and. timestep.le.999    ) &
		write(filename,'(a,a5,i3)') trim(basename),'.0000'  ,timestep
		if(timestep.ge.1000    .and. timestep.le.9999   ) &
		write(filename,'(a,a4,i4)') trim(basename),'.000'   ,timestep
		if(timestep.ge.10000   .and. timestep.le.99999  ) &
		write(filename,'(a,a3,i5)') trim(basename),'.00'    ,timestep
		if(timestep.ge.100000  .and. timestep.le.999999 ) &
		write(filename,'(a,a2,i6)') trim(basename),'.0'     ,timestep
		if(timestep.ge.1000000 .and. timestep.le.9999999) &
		write(filename,'(a,a1,i7)') trim(basename),'.'	,timestep

		!Remove any surplus blanks
		filename = trim(filename)

    end subroutine get_Timestep_FileName

end module functions
