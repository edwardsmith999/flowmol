!----------------------------------------------------------------------|
!    FEM_1D IS THE MAIN PROGRAM FOR THE SOLUTION OF A SET OF	       |
!    1D EQUATIONS IN A CARTESIAN OR CYLINDRICAL DOMAIN		           |
!															           |
!    PROGRAM WRITTEN BY G. KARAPETSAS, LAST UPDATE 09/12/2009          |
!    CHANGES MADE BY E. R. SMITH, LAST UPDATE 05/10/2014               |
!
!   Functions:
!
!   Variables:
!
!----------------------------------------------------------------------|

PROGRAM FEM_1D

	USE ELEMENTS_DATA
	USE COMMON_ARRAYS
	USE BASIC_PARAMETERS
	USE GAUSS_DATA
	USE BASIS_FUNCTION_1D
	USE TIME_INTEGRATION
	USE FLOW_PARAMETERS
	USE MACHINE_EPSILON
	USE FUNCTIONS
    
	IMPLICIT NONE

	INTEGER :: I, J, ITER, ISIGN, IDIM_Q, IEL
	INTEGER :: JBAND, INOD, INFO, IROW, IEQ, JEQ
	INTEGER :: INCREMENT =  0

	REAL(8) :: TRNS, ER_MAX, SNORM

	INTEGER, DIMENSION(25) :: IERROR

	INTEGER, PARAMETER :: IDIM_ER = NEQ_Q + 1
	REAL(8), DIMENSION(IDIM_ER,3) :: ERIND

	CHARACTER(len=1) :: CCL
	REAL(8) :: RSUM_OLD, RSUM_NEW

	NAMELIST/parameters/ eps2, As1, As2, As12, d1, d12, bslip, kappa, &
	                     mm, Equilibrium_Angle, Thetac, MASS_SURFACTANT, Peca, ba, ka, Ra, &
	                     SOLUBLE_SURFACTANT, Peb, Pem, kb, nc, WALL_ADSORPTION, Pecs, bs,  &
	                     ks, Rs, kas, Ras, ST_MODEL, Bl1, Bl2, Bl12, Kl1, Kl2, Kl12, &
	                     GEOMETRY, NXELa, NXELb, XPACK, MPa, MPb, OUTFE_STEP, NZELa

	!-----------------------------------------------------------------------
	!     READ NAMELIST'S DATA
	!-----------------------------------------------------------------------

	OPEN(32,file='SPARAMETERS',status='OLD', action='READ', &
				position='REWIND')
	READ(32,nml=parameters)
	CLOSE(32)

	!-----------------------------------------------------------------------
	!     OPEN INPUT/OUTPUT FILES
	!-----------------------------------------------------------------------

	OPEN(42, file='results/OUTFE.DAT', status='UNKNOWN', action='WRITE',&
					iostat=IERROR(1))
	IERROR(2) = 0
	OPEN(44, file='results/TIMEDATA1.DAT', status='UNKNOWN', action='WRITE',&
					iostat=IERROR(3))
	OPEN(47, file='results/TIMEDATA2.DAT', status='UNKNOWN', action='WRITE',&
					iostat=IERROR(3))
	OPEN(45, file='results/ERROR.DAT', status='UNKNOWN', action='WRITE',&
					iostat=IERROR(4))
	IERROR(5) = 0

	DO I = 1, 5
		IF ( IERROR(I) /= 0 ) THEN
			WRITE(42,*)'ERROR DURING OPENING PROCEDURE IN file', I
			WRITE(42,*)'PROGRAM STOP !!!'
			STOP
		ENDIF
	ENDDO

	!-----------------------------------------------------------------------
	!     WRITE PARAMATER VALUES TO OUTFE.DAT
	!-----------------------------------------------------------------------

	eps = DSQRT(eps2)

	WRITE(42,nml=parameters)
	WRITE(42,*)
	WRITE(42,*)'-------------------------------------------------------'
	WRITE(42,*)
 
	!-----------------------------------------------------------------------
	!     CHECK PARAMETERS VALUES
	!-----------------------------------------------------------------------

	IF ( MASS_SURFACTANT > 1.d-9 ) THEN

		! Check for Solubility
		IF ( .NOT.SOLUBLE_SURFACTANT ) THEN

			! Insoluble surfactant
			IF ( ka > 0.D0 .OR. ks > 0.D0 .OR. kb > 0.D0 ) THEN
				WRITE(42,*)'ERROR! INSOLUBLE SURFACTANT'
				WRITE(42,*)'ka, ks AND kb SHOULD BE ZERO!'
				STOP
			ENDIF

			SELECT CASE(GEOMETRY)
				CASE('CARTESIAN')
					IF ( MASS_SURFACTANT > 1.D0 ) THEN
						WRITE(42,*)'ERROR! INSOLUBLE SURFACTANT'
						WRITE(42,*)'MASS_SURFACTANT SHOULD BE < 1 !'
						STOP
					ENDIF
				CASE('CYLINDRICAL')
					IF ( MASS_SURFACTANT > 0.5D0 ) THEN
						WRITE(42,*)'ERROR! INSOLUBLE SURFACTANT'
						WRITE(42,*)'MASS_SURFACTANT SHOULD BE < 0.5 !'
						STOP
					ENDIF
			END SELECT

		ELSE

			!	Soluble surfactant

			!	Check for wall adsorption
			IF ( (.NOT.WALL_ADSORPTION) .AND. (ks > 0.D0)  ) THEN
				WRITE(42,*)'ERROR! SURFACTANT IS NOT ADSORBED ON THE WALL'
				WRITE(42,*)'ks SHOULD BE ZERO!'
				STOP
			ELSE IF ( (WALL_ADSORPTION) .AND. (ks < 1.D-9) ) THEN
				WRITE(42,*)'ERROR! ks is equal to ZERO!'
				WRITE(42,*)'WALL_ADSORPTION SHOULD BE .FALSE.'
			ENDIF

			! Check for Critical Micelle Concetration		
			IF ( MASS_SURFACTANT < 1.D0 ) THEN
				IF ( kb > 0.D0 ) THEN
					WRITE(42,*)'ERROR! BELOW CRITICAL MICELLE CONCETRATION'
					WRITE(42,*)'kb SHOULD BE ZERO!'
					STOP
				ENDIF
			ENDIF

		ENDIF

	ENDIF

!-----------------------------------------------------------------------

	CALL GAUSS_LINE
	CALL INITIALIZE_BASIS_FUNCTION_1D
	CALL MACHINE_ERROR

	!-----------------------------------------------------------------------
	!     INITIALIZE FLAG_NR, FORM_JAC
	!-----------------------------------------------------------------------

	FLAG_NR  = 'NRP'
	FORM_JAC = .TRUE.

	!-----------------------------------------------------------------------
	!     ALLOCATE ARRAYS
	!-----------------------------------------------------------------------

	NXEL = NXELa + NXELb
	!Number of elements, total and in a and b ?
	NNTOL = NXEL     ; NNTOLa = NXELa     ; NNTOLb = NXELb
	!Number of nodes, total and in a and b ?
	NNX = 2*NXEL + 1 ; NNXa = 2*NXELa + 1 ; NNXb = 2*NXELb + 1
	NODTOL_QEL = NNX ; NODTOLa_QEL = NNXa ; NODTOLb_QEL = NNXb 

	!Number of quadrature points, total and in a and b ?
	NODTOLa = NEQ_Q_a*NODTOLa_QEL
	NODTOLb = NEQ_Q_b*NODTOLb_QEL
	NODTOL  = NEQ_Q_a*NODTOLa_QEL + NEQ_Q_b*(NODTOLb_QEL-1)

	ALLOCATE( NM_Q_a(NNTOLa,NBF_Q+1),  STAT=IERROR(1) )
	ALLOCATE( NM_Q_b(NNTOLb,NBF_Q+1),  STAT=IERROR(2) )
	ALLOCATE( NM_QQ_a(NNTOLa,NBF_Q+1), STAT=IERROR(3) )
	ALLOCATE( NM_QQ_b(NNTOLb,NBF_Q+1), STAT=IERROR(4) )

	LOOP_CHECK_ALLOCATION: DO I = 1, 4
		IF(IERROR(I) .NE. 0)THEN
			WRITE(42,*)I,'ARRAY IS NOT ALLOCATED / SUB: MAIN_PROGRAM '
			WRITE(42,*)'STAT =', IERROR
			WRITE(42,*)'PROGRAM STOP !!!'
			STOP
		ENDIF
	ENDDO LOOP_CHECK_ALLOCATION

	!-----------------------------------------------------------------------
	!     INITIALIZE CONNECTIVITY NM ARRAYS
	!-----------------------------------------------------------------------

	IBAND  = 0

	NM_Q_a(:,1) = NBF_Q
	IDIM_Q = NM_Q_a(1,1)

    !Elements    1     2     3     4
    !         |-----|-----|-----|-----|
    !         x--x--x--x--x--x--x--x--x
    !Nodes    1  2  3  4  5  6  7  8  9
	DO IEL = 1, NNTOLa
		NM_Q_a(IEL,2) = 2*IEL-1
		NM_Q_a(IEL,3) = NM_Q_a(IEL,2)+1
		NM_Q_a(IEL,4) = NM_Q_a(IEL,3)+1

		NM_QQ_a(IEL,1) = IDIM_Q
		NM_QQ_a(IEL,2) = NM_Q_a(IEL,2)*NEQ_Q_a - NEQ_Q_a + 1
		NM_QQ_a(IEL,3) = NM_QQ_a(IEL,2) + NEQ_Q_a
		NM_QQ_a(IEL,4) = NM_QQ_a(IEL,3) + NEQ_Q_a

		JBAND = 2*( (NM_QQ_a(IEL,4)+NEQ_Q_a-1) - NM_QQ_a(IEL,2) ) + 1
		IBAND = MAX0( IBAND, JBAND )
	ENDDO

	NM_Q_b(:,1) = NBF_Q
	IDIM_Q = NM_Q_b(1,1)

	DO IEL = 1, NNTOLb
		NM_Q_b(IEL,2) = 2*IEL-1
		NM_Q_b(IEL,3) = NM_Q_b(IEL,2)+1
		NM_Q_b(IEL,4) = NM_Q_b(IEL,3)+1

		NM_QQ_b(IEL,1) = IDIM_Q
		NM_QQ_b(IEL,2) = NM_Q_b(IEL,2)*NEQ_Q_b - NEQ_Q_b + 1
		NM_QQ_b(IEL,3) = NM_QQ_b(IEL,2) + NEQ_Q_b
		NM_QQ_b(IEL,4) = NM_QQ_b(IEL,3) + NEQ_Q_b

		NM_QQ_b(IEL,:) = NM_QQ_b(IEL,:) + NODTOLa - NEQ_Q_b 

		JBAND = 2*( (NM_QQ_b(IEL,4)+NEQ_Q_b-1) - NM_QQ_b(IEL,2) ) + 1
		IBAND = MAX0( IBAND, JBAND )
	ENDDO

	MBAND  = ( IBAND + 1 ) / 2
	HBAND  = ( IBAND - 1 ) / 2
	IDIM_A  = IBAND + HBAND

	!-----------------------------------------------------------------------
	!			ALLOCATE ARRAYS
	!-----------------------------------------------------------------------

	ALLOCATE( A(IDIM_A,NODTOL),	          STAT=IERROR(1) )
	ALLOCATE( B(NODTOL),			      STAT=IERROR(2) )
	ALLOCATE( S(NODTOL),			      STAT=IERROR(3) )
	ALLOCATE( A_ip(NODTOL),				  STAT=IERROR(4) )
	ALLOCATE( A_pi(NODTOL),				  STAT=IERROR(5) )
	ALLOCATE( S_c(NODTOL+1),			  STAT=IERROR(6) )
	ALLOCATE( IPVT(NODTOL),				  STAT=IERROR(7) )
	ALLOCATE( TQa(NODTOLa_QEL,NEQ_Q_a),   STAT=IERROR(8) )
	ALLOCATE( TQao(NODTOLa_QEL,NEQ_Q_a),  STAT=IERROR(9) )
	ALLOCATE( TQao1(NODTOLa_QEL,NEQ_Q_a), STAT=IERROR(10) )
	ALLOCATE( TQap(NODTOLa_QEL,NEQ_Q_a),  STAT=IERROR(11) )
	ALLOCATE( TQb(NODTOLb_QEL,NEQ_Q_b),   STAT=IERROR(12) )
	ALLOCATE( TQbo(NODTOLb_QEL,NEQ_Q_b),  STAT=IERROR(13) )
	ALLOCATE( TQbo1(NODTOLb_QEL,NEQ_Q_b), STAT=IERROR(14) )
	ALLOCATE( TQbp(NODTOLb_QEL,NEQ_Q_b),  STAT=IERROR(15) )
	ALLOCATE( TQo_GRAPH(NODTOL_QEL,20),   STAT=IERROR(16) )
	ALLOCATE( Xa(NODTOLa_QEL),		      STAT=IERROR(17) )
	ALLOCATE( Xb(NODTOLb_QEL),		      STAT=IERROR(18) )

	LOOP_CHECK_ALLOCATION2: DO I = 1, 18
		IF(IERROR(I) .NE. 0)THEN
			WRITE(42,*)I,'ARRAY IS NOT ALLOCATED / SUB: MAIN_PROGRAM'
			WRITE(42,*)'STAT =', IERROR
			WRITE(42,*)'PROGRAM STOP !!!'
			STOP
		ENDIF
	ENDDO LOOP_CHECK_ALLOCATION2

	!-----------------------------------------------------------------------
	!		CONSTRUCT MESH
	!-----------------------------------------------------------------------
	CALL MESH
	
	!-----------------------------------------------------------------------
	!		INITIAL STATE
	!-----------------------------------------------------------------------
	MASS_FLUIDi = 2.D0/3.D0

	DO I=1,NODTOLa_QEL
		TQa(I,1) =   1.5D0*MASS_FLUIDi*(1.d0-Xa(I)**2)
		TQa(I,2) = - 3.d0*MASS_FLUIDi*Xa(I)
		TQa(I,3) = - 3.d0*MASS_FLUIDi
	ENDDO

	IF ( MASS_SURFACTANT > 1.d-9 ) THEN
		CALL CALCULATE_INITIAL_CONCETRATIONS
		DO I=1,NODTOLa_QEL
			TQa(I,4) =   cai
			TQa(I,5) =   ci
			TQa(I,6) =   mi
			TQa(I,7) =   csi
		ENDDO
	ELSE
		TQa(:,4:7) =   0.D0
	ENDIF

	TQb(:,:) = 0.D0
	Xc = 1.d0

	TQao  = TQa
	TQbo  = TQb
	Xco   = Xc

	TQao1 = TQao
	TQbo1 = TQbo
	Xco1  = Xco

	TIME = 0.D0
	TRNS = 1.D0
	Dt = Dt_min

	CALL CALCULATE_MASS(TQa, TQb, Xc, MASS_FLUIDi, MASS_SURFi)

	!  WRITE OUTPUT FILES
	CALL WRITE_FILES('Tn', .TRUE., TIME, Dt)

	!-----------------------------------------------------------------------
	!		DEFINE GRAPH
	!-----------------------------------------------------------------------

	PO_MIN  = Time
	PO_MAX  = Time_Limit
	PO_STEP = 5.d-1
	PO	= PO_MIN
   
	PI_MIN  = Time
	PI_MAX  = Time_Limit
	PI_STEP = 5.d-2
	PI	= PI_MIN

	!-----------------------------------------------------------------------
	!-----------------------------------------------------------------------
	!							INTEGRATION ON TIME
	!	This is the top level integration loop
	!-----------------------------------------------------------------------
	!-----------------------------------------------------------------------

	IMPLICIT_EULER: DO WHILE( TIME <= TIME_LIMIT )

		INCREMENT = INCREMENT + 1
		TIME = TIME + DT

		IF(MOD(INCREMENT,OUTFE_STEP) == 0)THEN
			WRITE(42,*)'---------------------------------------------------'
			WRITE(42,*)'---------------------------------------------------'
			WRITE(42,*)' TIME STEP = ',INCREMENT,' TIME = ', TIME
			WRITE(42,*)' Dt = ', Dt, ' FACTOR = ', TRNS
			WRITE(42,*)'---------------------------------------------------'
			WRITE(42,*)
		ENDIF

		!-----------------------------------------------------------------------
		!			NEWTON ITERATION LOOP
		!-----------------------------------------------------------------------

		ITER=0
		RSUM_OLD = 1.D0
		RSUM_NEW = 1.D0

		NEWTON_RAPHSON: DO WHILE( ITER.LT.NITER .AND. RSUM_NEW.GT.ERROR )

			ITER = ITER + 1

			!-----------------------------------------------------------------------
			!					INITIALIZE ARRAYS
			!-----------------------------------------------------------------------

			IF (FORM_JAC) THEN
				A = 0.D0 ; A_p = 0.D0 ; A_ip = 0.D0 ; A_pi = 0.D0
			ENDIF

			B = 0.D0 ; B_p = 0.D0

			!-----------------------------------------------------------------------
			!		COMPUTE AND STORE RESIDUALS AND JACOBIAN
			!-----------------------------------------------------------------------

			CALL DOMI
			CALL BOUNDARY_RESIDUALS

			!-----------------------------------------------------------------------
			!			ESSENTIAL RESIDUAL & ESSENTIAL JACOBIAN
			!-----------------------------------------------------------------------

			INOD = NM_Q_a(1,2)
			IROW = NM_QQ_a(1,2) + 1
			B(IROW) = TQa(INOD,2)	! Hx(x=0) = 0
			IF (FORM_JAC) THEN
				A(:,IROW)     = 0.D0
				A_ip(IROW)    = 0.D0
				A(IBAND,IROW) = 1.D0
			ENDIF

			INOD = NM_Q_a(NNTOLa,4)
			IROW = NM_QQ_a(NNTOLa,4)
			B(IROW) = TQa(INOD,1)	!  H(x=xc) = 0
			IF (FORM_JAC) THEN
				A(:,IROW)     = 0.D0
				A_ip(IROW)    = 0.D0
				A(IBAND,IROW) = 1.D0
			ENDIF

			IF (MASS_SURFACTANT > 1.D-9) THEN
				IF ( MASS_SURFACTANT < 1.D0 ) THEN

					DO IEL = 1, NNTOLa
						DO I = 1, NBF_Q
							INOD = NM_Q_a(IEL,I+1)
							IROW = NM_QQ_a(IEL,I+1) + 5
							B(IROW) = TQa(INOD,6)	! Below CMC M = 0 
							IF (FORM_JAC) THEN
								A(:,IROW)     = 0.D0
								A_ip(IROW)    = 0.D0
								A(IBAND,IROW) = 1.D0
							ENDIF
						ENDDO
					ENDDO

				ENDIF

				IF (.NOT.WALL_ADSORPTION) THEN

					DO IEL = 1, NNTOLa
						DO I = 1, NBF_Q
							INOD = NM_Q_a(IEL,I+1)
							IROW = NM_QQ_a(IEL,I+1) + 6
							B(IROW) = TQa(INOD,7)	! Cs(x<xc) = 0
							IF (FORM_JAC) THEN
								A(:,IROW)     = 0.D0
								A_ip(IROW)    = 0.D0
								A(IBAND,IROW) = 1.D0
							ENDIF
						ENDDO
					ENDDO

					DO IEL = 1, NNTOLb
						DO I = 1, NBF_Q
							INOD = NM_Q_b(IEL,I+1)
							IROW = NM_QQ_b(IEL,I+1)
							B(IROW) = TQb(INOD,1)				! Cs(x>xc) = 0
							IF (FORM_JAC) THEN
								A(:,IROW)     = 0.D0
								A_ip(IROW)    = 0.D0
								A(IBAND,IROW) = 1.D0
							ENDIF
						ENDDO
					ENDDO

				ENDIF

			ELSE

				DO IEL = 1, NNTOLa
				DO I = 1, NBF_Q
				DO JEQ = 4, NEQ_Q_a
					IROW = NM_QQ_a(IEL,I+1) + JEQ - 1
					INOD = NM_Q_a(IEL,I+1)
					B(IROW) = TQa(INOD,JEQ)
					IF (FORM_JAC) THEN
						A(:,IROW)     = 0.D0
						A_ip(IROW)    = 0.D0
						A(IBAND,IROW) = 1.D0
					ENDIF
				ENDDO
				ENDDO
				ENDDO

				DO IEL = 1, NNTOLb
				DO I = 1, NBF_Q
				DO JEQ = 1, NEQ_Q_b
					IROW = NM_QQ_b(IEL,I+1) + JEQ - 1
					INOD = NM_Q_b(IEL,I+1)
					B(IROW) = TQb(INOD,JEQ)
					IF (FORM_JAC) THEN
						A(:,IROW)     = 0.D0
						A_ip(IROW)    = 0.D0
						A(IBAND,IROW) = 1.D0
					ENDIF
				ENDDO
				ENDDO
				ENDDO

			ENDIF

			!-----------------------------------------------------------------------
			!			SOLVE LINEAR SYSTEM
			!-----------------------------------------------------------------------

			IF (FLAG_NR == 'NRP') THEN
				! A Factorization (mkl library)
				CALL DGBTRF(NODTOL, NODTOL, HBAND, & 
                            HBAND, A, IBAND+HBAND, IPVT, INFO)

				! Aip <-- A^(-1).Aip (mkl library)
				CALL DGBTRS('T', NODTOL, HBAND, HBAND, &
                             1, A, IBAND+HBAND, IPVT, A_ip, &
					NODTOL, INFO )
			ENDIF

			! B <-- A^(-1).B (mkl library)
			CALL DGBTRS('T', NODTOL, HBAND, HBAND, 1, A, IBAND+HBAND, IPVT, B, &
					NODTOL, INFO )

			IF (FLAG_NR == 'NRP') THEN
				! Ap =  Ap - Api.A^(-1).Aip
				A_p = A_p - DOT_PRODUCT(A_pi, A_ip)
			ENDIF

			!	S_p = (Ap-Api.A^(-1).Aip)^(-1).(Bp-Api.A^(-1).B)
			S_p = ( B_p - DOT_PRODUCT(A_pi, B) ) / A_p

			!	S = A^(-1).B - A^(-1).Aip.Sp
			S =  B - A_ip*S_p

			!	S_c = (S, Sp)
			S_c(1:NODTOL) = S(1:NODTOL)
			S_c(NODTOL+1) = S_p

			!-----------------------------------------------------------------------
			!							CHECK CONVERGENCE
			!-----------------------------------------------------------------------

			CALL CHECK_CONVERGENCE&
			(INCREMENT, ITER, FLAG_NR, CCL, RSUM_NEW, RSUM_OLD, S_c, NODTOL+1)

			IF (FULL_NR) FLAG_NR = 'NRP'

			IF ( FLAG_NR == 'NRP' ) THEN
				FORM_JAC = .TRUE.
			ELSE
				FORM_JAC = .FALSE.
			ENDIF

			IF ( ISNAN(RSUM_NEW) ) THEN
				WRITE(42,*)
				WRITE(42,*) 'TIME STEP = ',INCREMENT,'	TIME = ', TIME
				WRITE(42,*)' Dt = ', Dt, ' FACTOR = ', TRNS
				WRITE(42,*) 'PROGRAM EXIT: NaN ENCOUNTERED'
				STOP
			ENDIF

			IF( CCL == 'Y' )THEN
				CYCLE NEWTON_RAPHSON
			ENDIF

			!-----------------------------------------------------------------------
			!						UPDATE SOLUTION VECTOR
			!-----------------------------------------------------------------------

			DO IEQ = 1, NEQ_Q_a
				J = IEQ
				TQa(1,IEQ) = TQa(1,IEQ) - S(J)
			ENDDO

			DO IEL = 1, NNTOLa
				IDIM_Q = NM_Q_a(IEL,1)
				DO INOD = 2, IDIM_Q
				DO IEQ = 1, NEQ_Q_a
					I = NM_Q_a(IEL,INOD+1)
					J = NM_QQ_a(IEL,INOD+1) + IEQ - 1
					TQa(I,IEQ) = TQa(I,IEQ) - S(J)
				ENDDO
				ENDDO
			ENDDO

			DO IEQ = 1, NEQ_Q_b
				J = NODTOLa - NEQ_Q_b + IEQ
				TQb(1,IEQ) = TQb(1,IEQ) - S(J)
			ENDDO

			DO IEL = 1, NNTOLb
				IDIM_Q = NM_Q_b(IEL,1)
				DO INOD = 2, IDIM_Q
				DO IEQ = 1, NEQ_Q_b
					I = NM_Q_b(IEL,INOD+1)
					J = NM_QQ_b(IEL,INOD+1) + IEQ - 1
					TQb(I,IEQ) = TQb(I,IEQ) - S(J)
				ENDDO
				ENDDO
			ENDDO

			Xc = Xc - S_p

			!-----------------------------------------------------------------------

		ENDDO  NEWTON_RAPHSON

		!-----------------------------------------------------------------------

		IF(ITER.GE.NITER .AND. RSUM_NEW.GT.ERROR)THEN
			IF ( Dt > Dt_min ) THEN
				TIME = TIME - Dt
				Dtn = Dt/TRNS
				Dt = Dt*0.5D0
				TRNS = Dt/Dtn
				TQap = TQao + (TQao-TQao1)*TRNS
				TQbp = TQbo + (TQbo-TQbo1)*TRNS
				Xcp  = Xco + (Xco-Xco1)*TRNS
				TQa  = TQap
				TQb  = TQbp
				Xc   = Xcp
				CYCLE IMPLICIT_EULER
			ELSE
				WRITE(42,*)'TIME STEP = ',INCREMENT,'	TIME = ', TIME
				WRITE(42,*)'Dt = ', Dt, ' FACTOR = ', TRNS
				WRITE(42,*)'NEWTON METHOD FAILED TO CONVERGE'
				WRITE(42,*)'PROGRAM STOP !!!'
				CALL WRITE_FILES('To', .TRUE., TIME, Dt)
				STOP
			ENDIF
		ENDIF

		!-----------------------------------------------------------------------
		!		WRITE SOLUTION
		!-----------------------------------------------------------------------

		IF( DABS(Time - PO_min) .LE. Dt_min )THEN

			CALL WRITE_FILES('Tn', .TRUE., TIME, Dt)

			IF ( Time > 1.d1-Dt_max ) PO_STEP = 5.0D0
			IF ( Time > 1.d2-Dt_max ) PO_STEP = 5.0D1
			IF ( Time > 1.d3-Dt_max ) PO_STEP = 5.0D2
			IF ( Time > 1.d4-Dt_max ) PO_STEP = 5.0D3
			IF ( Time > 1.d5-Dt_max ) PO_STEP = 5.0D4

			PO = PO_min + PO_STEP

		ENDIF

		PO_min = PO

		IF ( DABS(Time - PI_min) .LE. Dt ) THEN
			CALL WRITE_FILES('Tn',.FALSE.,TIME,Dt)
			PI = 10.d0**(Log10(TIME) + PI_STEP)
		ENDIF

		PI_min = PI

		!-----------------------------------------------------------------------
		!						ADJUST TIME-STEP
		!-----------------------------------------------------------------------

		DO I = 1, NEQ_Q_a 
			ERIND(I,1) = DOT_PRODUCT( TQa(:,I)-TQap(:,I), TQa(:,I)-TQap(:,I) )
			ERIND(I,2) = DOT_PRODUCT( TQa(:,I), TQa(:,I) )
            if (ERIND(I,2) .gt. 1e-7) then
    			ERIND(I,3) = DSQRT( ERIND(I,1)/ERIND(I,2) )
            else
                ERIND(I,3) = 0.d0
            endif
		ENDDO

		ERIND(NEQ_Q_a,1) = ERIND(NEQ_Q_a,1) + &
						DOT_PRODUCT( TQb(:,1)-TQbp(:,1), TQb(:,1)-TQbp(:,1) )
		ERIND(NEQ_Q_a,2) = ERIND(NEQ_Q_a,2) + DOT_PRODUCT( TQb(:,1), TQb(:,1) )
        if (ERIND(NEQ_Q_a,2) .gt. 1e-7) then
		    ERIND(NEQ_Q_a,3) = DSQRT( ERIND(NEQ_Q_a,1)/ERIND(NEQ_Q_a,2) )
        else
            ERIND(NEQ_Q_a,3) = 0.d0
        endif
		ERIND(NEQ_Q_a+1,1) = (Xc-Xcp)**2
		ERIND(NEQ_Q_a+1,2) = Xc**2
        if (ERIND(NEQ_Q_a+1,2) .gt. 1e-7) then
    		ERIND(NEQ_Q_a+1,3) = DSQRT( ERIND(NEQ_Q_a+1,1)/ERIND(NEQ_Q_a+1,2) )
        else
            ERIND(NEQ_Q_a+1,3) = 0.d0
        endif

		IF (MASS_SURFACTANT < 1.d-9) ERIND(4:7,3) = 0.D0

		ER_MAX   = MAXVAL( ERIND(1:IDIM_ER,3) )

		IF(MOD(INCREMENT,OUTFE_STEP) == 0)THEN
			write(42,*)'ERRORS IN TIME PREDICTIONS, USING EXPLICIT EULER'
			write(42,*)
			write(42,*)'H   = ', ERIND(1,3)
			write(42,*)'Hx  = ', ERIND(2,3)
			write(42,*)'Hxx = ', ERIND(3,3)
			write(42,*)'Ca  = ', ERIND(4,3)
			IF (SOLUBLE_SURFACTANT) THEN
				write(42,*)'c   = ', ERIND(5,3)
				write(42,*)'m   = ', ERIND(6,3)
			ENDIF
			IF (WALL_ADSORPTION) write(42,*)'Cs  = ', ERIND(7,3)
			write(42,*)'Xc  = ', ERIND(8,3)
			write(42,*)'MAXIMUM ERROR = ', ER_MAX
		ENDIF
 
		TRNS = eps_t/ER_MAX

		IF (TRNS > 2.d0) TRNS = 2.d0

		!-----------------------------------------------------------------------

		Dtn  = Dt
		Dt   = Dt * TRNS

		IF ( TIME + Dt > PO_min) THEN
			Dt = PO_min - TIME
			TRNS = Dt/Dtn
		ENDIF

		IF( Dt .LT. Dt_min )THEN
			Dt  = Dt_min
			TRNS = Dt/Dtn
		ELSE IF( Dt .GT. Dt_max )THEN
			Dt  = Dt_max
			TRNS = Dt/Dtn
		ENDIF

		!-----------------------------------------------------------------------
		!						UPDATE SOLUTION VECTORS
		!-----------------------------------------------------------------------

		TQap  = TQa + (TQa-TQao)*TRNS
		TQbp  = TQb + (TQb-TQbo)*TRNS
		Xcp   = Xc + (Xc-Xco)*TRNS

		TQao1 = TQao
		TQbo1 = TQbo
		TQao  = TQa

		TQbo  = TQb
		TQa   = TQap
		TQb   = TQbp

		Xco1  = Xco
		Xco   = Xc
		Xc    = Xcp

		!-----------------------------------------------------------------------

	ENDDO IMPLICIT_EULER

	!-----------------------------------------------------------------------
	!			DEALLOCATE ARRAYS
	!-----------------------------------------------------------------------

	DEALLOCATE( A,	    	STAT=IERROR(1)  )
	DEALLOCATE( B,	    	STAT=IERROR(2)  )
	DEALLOCATE( S,  		STAT=IERROR(3)  )
	DEALLOCATE( A_ip,	    STAT=IERROR(4)  )
	DEALLOCATE( A_pi,   	STAT=IERROR(5)  )
	DEALLOCATE( S_c,	    STAT=IERROR(6)  )
	DEALLOCATE( IPVT,   	STAT=IERROR(7)  )
	DEALLOCATE( TQa,	    STAT=IERROR(8)  )
	DEALLOCATE( TQao,   	STAT=IERROR(9)  )
	DEALLOCATE( TQao1,	    STAT=IERROR(10) )
	DEALLOCATE( TQap,	    STAT=IERROR(11) )
	DEALLOCATE( TQb,	    STAT=IERROR(12) )
	DEALLOCATE( TQbo,	    STAT=IERROR(13) )
	DEALLOCATE( TQbo1,	    STAT=IERROR(14) )
	DEALLOCATE( TQbp,	    STAT=IERROR(15) )
	DEALLOCATE( TQo_GRAPH,  STAT=IERROR(16) )
	DEALLOCATE( Xa,			STAT=IERROR(17) )
	DEALLOCATE( Xb,			STAT=IERROR(18) )
	DEALLOCATE( NM_Q_a,		STAT=IERROR(19) )
	DEALLOCATE( NM_Q_b,		STAT=IERROR(20) )
	DEALLOCATE( NM_QQ_a,	STAT=IERROR(21) )
	DEALLOCATE( NM_QQ_b,	STAT=IERROR(22) )

	LOOP_CHECK_DEALLOCATION: DO I = 1, 22
		IF(IERROR(I) .NE. 0)THEN
			WRITE(42,*)I,'ARRAY IS NOT DEALLOCATED / SUB: MAIN_PROGRAM '
			WRITE(42,*)'STAT =', IERROR
			WRITE(42,*)'PROGRAM STOP !!!'
			STOP
		ENDIF
	ENDDO LOOP_CHECK_DEALLOCATION

	!----------------------------------------------------------------------

END PROGRAM FEM_1D

!-----------------------------------------------------------------------
!		SUBROUTINE CALCULATE_INITIAL_CONCETRATIONS
!-----------------------------------------------------------------------

SUBROUTINE CALCULATE_INITIAL_CONCETRATIONS

	USE FLOW_PARAMETERS
	USE FUNCTIONS
	USE BASIC_PARAMETERS, only: ERROR

	implicit none

	INTEGER :: I
	REAL(8) :: RSUM, x_new, xo

	write(42,*) 'Calculate initial concetrations:'
	write(42,*)

	xo = 1.d-5
	RSUM = 1.d0
	I = 0

	LOOP: DO 

		I = I + 1

		IF ( MASS_SURFACTANT .GE. 1.d0) THEN
			x_new = xo - Fm(xo)/Fm1(xo)
		ELSE
			IF (SOLUBLE_SURFACTANT) THEN
				x_new = xo - Fc(xo)/Fc1(xo)
			ELSE
				x_new = xo - Fg(xo)/Fg1(xo)
			ENDIF
		ENDIF

		RSUM = ABS((x_new-xo)/xo)

		write(42,*) 'Iteration ', I
		write(42,*) 'x_new = ', x_new
		write(42,*) 'Error = ', RSUM
    
		IF ( RSUM .LT. ERROR ) EXIT LOOP

		IF ( I .GT. 50 ) THEN
			WRITE(42,*) 'NEWTON RAPHSON FAILED TO CONVERGE', &
					' // SUB: CALCULATE_INITIAL_CONCETRATIONS'
			STOP
		ENDIF 

		xo = x_new

	ENDDO LOOP

	IF ( MASS_SURFACTANT .GE. 1.d0) THEN
		mi  = xo
		ci  = mi**(1.d0/nc)
		Cai = Ra*ci/(Ra*ci+1.d0)
		Csi = 0.D0
	ELSE
		IF (SOLUBLE_SURFACTANT) THEN
			mi  = 0.d0
			ci  = xo
			Cai = Ra*ci/(Ra*ci+1.d0)
			Csi = 0.D0
		ELSE
			mi  = 0.d0
			ci  = 0.d0
			Cai = xo
			Csi = 0.d0
			IF (Cai > 1.d0) THEN
				WRITE(42,*)'ERROR! Cai SHOULD BE LESS OR EQUAL TO 1' 
				WRITE(42,*)'FOR INSOLUBLE SURFACTANT'
				WRITE(42,*)'SUB : CALCULATE_INITIAL_CONCETRATIONS'
				STOP
			ENDIF
		ENDIF
	ENDIF

	write(42,*)
	write(42,'(4A14)') 'Cai', 'ci', 'mi', 'Csi'
	write(42,'(4E14.6)') Cai, ci, mi, Csi
	write(42,*)

END SUBROUTINE CALCULATE_INITIAL_CONCETRATIONS

!-----------------------------------------------------------------------
!		SUBROUTINE DOMI
!-----------------------------------------------------------------------

SUBROUTINE DOMI

	USE ELEMENTS_DATA
	USE COMMON_ARRAYS
	USE GAUSS_DATA
	USE BASIS_FUNCTION_1D
	USE MATRIX_STORAGE_SUBROUTINES
	USE TIME_INTEGRATION, only: Dt, time
	USE FLOW_PARAMETERS
	USE BASIC_PARAMETERS, only: FORM_JAC
	USE MACHINE_EPSILON

	IMPLICIT NONE

	INTEGER, PARAMETER :: IDNM = NBF_Q
	INTEGER, PARAMETER :: IDDFDX = NBF_Q

	INTEGER :: NELEM, I, J, II, IG, IW, JW, JJ, INOD, ICOL, IEQ
	INTEGER :: ROI, INNTOL, INEQ, INDTL
	REAL(8) :: AJACX_Q, WET_Q

	INTEGER, DIMENSION(10)  :: IERROR

	REAL(8), DIMENSION(IDDFDX)  :: DFDX_Q

	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NMQ, NMQQ

	REAL(8), ALLOCATABLE, DIMENSION(:) :: X

	REAL(8), ALLOCATABLE, DIMENSION(:,:) :: TQ, TQo
	REAL(8), ALLOCATABLE, DIMENSION(:,:) :: TERM, TEMP_R, TEMP_X

	REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: TEMP_J

	!-----------------------------------------------------------------------
	!			ITERATE OVER ALL REGIONS
	!-----------------------------------------------------------------------

	LOOP_REGIONS: DO ROI = 1, NORE

		SELECT CASE(ROI)

		CASE(1)

            ! As ROI is droplet, store number of elements NNTOLa, 
            !number of equation NEQ_Q_a and number of nodes NODTOLa_QEL
			INNTOL = NNTOLa
			INEQ   = NEQ_Q_a
			INDTL  = NODTOLa_QEL

			! ALLOCATE ARRAYS

			ALLOCATE( TQ(INDTL,INEQ),		STAT = IERROR(1) )
			ALLOCATE( TQo(INDTL,INEQ),	    STAT = IERROR(2) )
			ALLOCATE( X(INDTL),			    STAT = IERROR(3) )
			ALLOCATE( NMQ(INNTOL,NBF_Q+1),  STAT = IERROR(4) )
			ALLOCATE( NMQQ(INNTOL,NBF_Q+1), STAT = IERROR(5) )

			LOOP_CHECK_ALLOCATION: DO I = 1, 5
				IF(IERROR(I) .NE. 0)THEN
					WRITE(42,*)I,'ARRAY IS NOT ALLOCATED / SUB: DOMI //1'
					WRITE(42,*)'STAT =', IERROR
					WRITE(42,*)'PROGRAM STOP !!!'
					STOP
				ENDIF
			ENDDO LOOP_CHECK_ALLOCATION

			NMQ  = NM_Q_a
			NMQQ = NM_QQ_a

			TQ = TQa ; TQo = TQao
			X  = Xa

		CASE(2)

            ! As ROI is vapour, store number of elements NNTOLb, 
            !number of equation NEQ_Q_b and number of nodes NODTOLb_QEL
			INNTOL = NNTOLb
			INEQ   = NEQ_Q_b
			INDTL  = NODTOLb_QEL

			! ALLOCATE ARRAYS
			ALLOCATE( TQ(INDTL,INEQ),		STAT = IERROR(1) )
			ALLOCATE( TQo(INDTL,INEQ),	    STAT = IERROR(2) )
			ALLOCATE( X(INDTL),			    STAT = IERROR(3) )
			ALLOCATE( NMQ(INNTOL,NBF_Q+1),  STAT = IERROR(4) )
			ALLOCATE( NMQQ(INNTOL,NBF_Q+1), STAT = IERROR(5) )

			LOOP_CHECK_ALLOCATION2: DO I = 1, 5
				IF(IERROR(I) .NE. 0)THEN
					WRITE(42,*)I,'ARRAY IS NOT ALLOCATED / SUB: DOMI //2'
					WRITE(42,*)'STAT =', IERROR
					WRITE(42,*)'PROGRAM STOP !!!'
					STOP
				ENDIF
			ENDDO LOOP_CHECK_ALLOCATION2

			NMQ = NM_Q_b
			NMQQ = NM_QQ_b

			TQ = TQb ; TQo = TQbo
			X  = Xb

		END SELECT

		!-----------------------------------------------------------------------
		!	ALLOCATE TEMPORARY ARRAYS
		!-----------------------------------------------------------------------
		ALLOCATE( TERM(NBF_Q,INEQ),				 STAT = IERROR(1) )
		ALLOCATE( TEMP_R(NBF_Q,INEQ),		     STAT = IERROR(2) )
		ALLOCATE( TEMP_X(NBF_Q,INEQ),		     STAT = IERROR(3) )
		ALLOCATE( TEMP_J(NBF_Q,NBF_Q,INEQ,INEQ), STAT = IERROR(4) )

		LOOP_CHECK_ALLOCATION3: DO I = 1, 4
			IF(IERROR(I) .NE. 0)THEN
				WRITE(42,*)I,'ARRAY IS NOT ALLOCATED / SUB: DOMI //3 '
				WRITE(42,*)'STAT =', IERROR
				WRITE(42,*)'PROGRAM STOP !!!'
				STOP
			ENDIF
		ENDDO LOOP_CHECK_ALLOCATION3

		!-----------------------------------------------------------------------
		!	ITERATE OVER ALL ELEMENTS
		!-----------------------------------------------------------------------
		LOOP_ELEMENTS: DO NELEM = 1, INNTOL

			!-----------------------------------------------------------------------
			!		INITIALIZE ARRAYS
			!-----------------------------------------------------------------------
			TEMP_R = 0.D0 ; TEMP_X = 0.D0 ; TEMP_J = 0.D0

			!-----------------------------------------------------------------------
			!		ITERATE OVER EACH GAUSS POINT IN AN ELEMENT
			!-----------------------------------------------------------------------
			LOOP_GAUSS: DO IG = 1, NINT_3P

				!-----------------------------------------------------------------------
				!		CALCULATE DERIVATIVES OF BASIS FUNCTIONS AND TRANFORMATION
				!		JACOBIAN AT THE GAUSS POINTS
				!-----------------------------------------------------------------------
				AJACX_Q = 0.D0
				DO II = 1, NBF_Q
					INOD = NMQ(NELEM,II+1)
					AJACX_Q = AJACX_Q + X(INOD)*DTF_Q(II,IG)
				ENDDO

				AJACX_Q = DABS(AJACX_Q)
				DFDX_Q(:) = DTF_Q(:,IG) / AJACX_Q
				WET_Q = WO_3P(IG) * AJACX_Q

				!-----------------------------------------------------------------------
				CALL EQUATIONS ( ROI, NELEM, NINT_3P, IG, TF_Q, DFDX_Q, NMQ, INNTOL, &
									TQ, TQo, X, INDTL, INEQ, Xc, TERM )

				!-----------------------------------------------------------------------
				!			FORM THE WORKING RESIDUAL VECTOR IN ELEMENT IEL
				!-----------------------------------------------------------------------
				TEMP_R  = TEMP_R  + TERM * WET_Q

				!-----------------------------------------------------------------------
				!			ITERATE OVER EACH NODE TO FORM JACOBIAN
				!-----------------------------------------------------------------------
				CALL JACOBIAN ( ROI, NELEM, NINT_3P, IG, TF_Q, DFDX_Q, NBF_Q, INEQ, &
								NMQ, INNTOL, TQ, TQo, X, INDTL, TERM, TEMP_X, TEMP_J, &
								WET_Q )

			ENDDO LOOP_GAUSS

			!-----------------------------------------------------------------------
			!		RESIDUAL VECTOR STORAGE
			!-----------------------------------------------------------------------

			CALL MATRIX_STORAGE_RESIDUAL&
				( TEMP_R, NMQQ(NELEM,2:), NBF_Q, INEQ, B, NODTOL )

			!-----------------------------------------------------------------------
			!		JACOBIAN MATRIX STORAGE
			!-----------------------------------------------------------------------

			IF (FORM_JAC) THEN

				CALL MATRIX_STORAGE_JACOBIAN_BAND&
					( TEMP_J, NMQQ(NELEM,2:), NMQQ(NELEM,2:), NBF_Q, NBF_Q, &
						INEQ, INEQ, A, IDIM_A, NODTOL, IBAND )

				! STORE MATRIX A_ip
				DO II = 1, NBF_Q
				DO JJ = 1, INEQ
					ICOL = NMQQ(NELEM,II+1) + JJ - 1
					A_ip(ICOL) = A_ip(ICOL) + TEMP_X(II,JJ)
				ENDDO
				ENDDO

			ENDIF

		ENDDO  LOOP_ELEMENTS

		!-----------------------------------------------------------------------
		!	DEALLOCATE ARRAYS
		!-----------------------------------------------------------------------

		DEALLOCATE( TQ,     STAT = IERROR(1) )
		DEALLOCATE( TQo,    STAT = IERROR(2) )
		DEALLOCATE( X,	    STAT = IERROR(3) )
		DEALLOCATE( NMQ,    STAT = IERROR(4) )
		DEALLOCATE( NMQQ,   STAT = IERROR(5) )
		DEALLOCATE( TERM,   STAT = IERROR(6) )
		DEALLOCATE( TEMP_R, STAT = IERROR(7) )
		DEALLOCATE( TEMP_X, STAT = IERROR(8) )
		DEALLOCATE( TEMP_J, STAT = IERROR(9) )

		LOOP_CHECK_DEALLOCATION: DO I = 1, 9
			IF(IERROR(I) .NE. 0)THEN
				WRITE(42,*)I,'ARRAY IS NOT DEALLOCATED / SUB: DOMI '
				WRITE(42,*)'STAT =', IERROR
				WRITE(42,*)'PROGRAM STOP !!!'
				STOP
			ENDIF
		ENDDO LOOP_CHECK_DEALLOCATION

	ENDDO  LOOP_REGIONS

END SUBROUTINE DOMI

!-----------------------------------------------------------------------
!				SUBROUTINE EQUATIONS
!Inputs:
!        ROI    -- Domain, inside or outside droplet
!        NELEM  -- Current element index
!        ININT  -- Number of Gauss points per element from NINT_3P (3)
!        IG     -- Current Gauss point index 
!        BFN_Q  -- Array of Gauss points (3,3)
!        DFDX_Q -- Array of derivative Gauss points (3,3)
!        NM     -- Array of node indices per element NM_Q in current ROI (NM_Q_a or NM_Q_b)
!        INNTOL -- Total number of elements inside and outside droplet
!        TQ     -- The MAIN array containing current values at all nodes for all equations 
!        TQo    -- MAIN array TQ at the previous timestep
!        X      -- Physical location of nodes (Fixed by MESH)
!        INDTL  -- Number of nodes in current ROI (NODTOLa_QEL or NODTOLb_QEL)
!        INEQ   -- Number of equations (7  i.e. Z,H,Hx,Hxx,Ca,M,Cs)
!        Xc     -- Position of Contact line
!        TERM   -- Working residual vector for element NELEM
!-----------------------------------------------------------------------

SUBROUTINE EQUATIONS ( ROI, NELEM, ININT, IG, BFN_Q, DFDX_Q, NM_Q, &
						INNTOL, TQ, TQo, X, INDTL, INEQ, Xc, TERM )

	USE ELEMENTS_DATA, only: NBF_Q
	USE COMMON_ARRAYS, only: Xco
	USE FLOW_PARAMETERS
	USE BASIC_PARAMETERS, only: SUPG, L1, L2
	USE TIME_INTEGRATION, only: Dt
	USE FUNCTIONS

	IMPLICIT NONE

	INTEGER, PARAMETER :: JDIM_Q = NBF_Q

	REAL(8) :: Z, XX, dZdx, dxdZ
	REAL(8) :: H, Ho, dHdt, dXcdt
	REAL(8) :: Hz, Hxz, Hxxz, Hx, Hxx, Hxxx
	REAL(8) :: Ca, Cao, Caz, Cax, dCadt
	REAL(8) :: C, Co, Cz, Cx, dCdt
	REAL(8) :: M, Mo, Mz, Mx, dMdt
	REAL(8) :: Cs, Cso, Csz, Csx, dCsdt
	REAL(8) :: Q, U, dPdX, S2, S2x, Um, Jcca, Jccs, Jcm

	REAL(8) :: BIFN, DBIFN, SBIFN
	INTEGER :: II, INOD, IW

	INTEGER, INTENT(IN) :: NELEM, IG, ININT, INDTL, INEQ, INNTOL, ROI

	REAL(8), INTENT(IN), DIMENSION (JDIM_Q,ININT)  :: BFN_Q
	REAL(8), INTENT(IN), DIMENSION(JDIM_Q)		   :: DFDX_Q
	INTEGER, INTENT(IN), DIMENSION(INNTOL,NBF_Q+1) :: NM_Q

	REAL(8), INTENT(IN) :: Xc

	REAL(8), INTENT(IN), DIMENSION (INDTL) :: X
	REAL(8), INTENT(IN), DIMENSION (INDTL,INEQ)  :: TQ, TQo

	REAL(8), INTENT(OUT), DIMENSION(JDIM_Q, INEQ) :: TERM

	!-----------------------------------------------------------------------
	!     SELECT REGION 1 OR 2 AND CALCULATE THE CORRESPONDING RESIDUALS
	!-----------------------------------------------------------------------

	SELECT CASE(ROI)

	CASE(1)

		!-----------------------------------------------------------------------
		!    CALCULATE THE DEPENDENT VARIABLES AND THEIR PARTIAL
		!    DERIVATIVES AT THE GAUSS POINTS IN X COORDINATE
		!-----------------------------------------------------------------------

		Z = 0.D0

		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			Z    =  Z  + X(INOD) * BFN_Q(II,IG)
		ENDDO
  
		H  = 0.D0  ; Ho  = 0.D0  ; Hz = 0.D0
		Hx = 0.D0  ; Hxz = 0.D0
		Hxx = 0.D0 ; Hxxz = 0.D0

		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			H  = H  + TQ(INOD,1) * BFN_Q(II,IG)
			Ho = Ho + TQo(INOD,1)* BFN_Q(II,IG)
			Hz = Hz + TQ(INOD,1) * DFDX_Q(II)

			Hx  = Hx  + TQ(INOD,2) * BFN_Q(II,IG)
			Hxz = Hxz + TQ(INOD,2) * DFDX_Q(II)

			Hxx  = Hxx  + TQ(INOD,3) * BFN_Q(II,IG)
			Hxxz = Hxxz + TQ(INOD,3) * DFDX_Q(II)
		ENDDO

		Ca = 0.D0  ; Cao  = 0.D0  ; Caz = 0.D0

		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			Ca  = Ca  + TQ(INOD,4) * BFN_Q(II,IG)
			Cao = Cao + TQo(INOD,4)* BFN_Q(II,IG)
			Caz = Caz + TQ(INOD,4) * DFDX_Q(II)
		ENDDO

		C  = 0.D0  ; Co  = 0.D0  ; Cz = 0.D0

		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			C  = C  + TQ(INOD,5) * BFN_Q(II,IG)
			Co = Co + TQo(INOD,5)* BFN_Q(II,IG)
			Cz = Cz + TQ(INOD,5) * DFDX_Q(II)
		ENDDO

		M  = 0.D0  ; Mo  = 0.D0  ; Mz = 0.D0

		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			M  = M  + TQ(INOD,6) * BFN_Q(II,IG)
			Mo = Mo + TQo(INOD,6)* BFN_Q(II,IG)
			Mz = Mz + TQ(INOD,6) * DFDX_Q(II)
		ENDDO

		Cs  = 0.D0  ; Cso  = 0.D0  ; Csz = 0.D0
				
		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			Cs  = Cs  + TQ(INOD,7) * BFN_Q(II,IG)
			Cso = Cso + TQo(INOD,7)* BFN_Q(II,IG)
			Csz = Csz + TQ(INOD,7) * DFDX_Q(II)
		ENDDO

		!-----------------------------------------------------------------------
		!				CALCULATE TIME DERIVATIVES
		!-----------------------------------------------------------------------

		dxdZ = Xc/L1
		dZdx = L1/Xc

		dXcdt = (Xc-Xco)/Dt
		dHdt  = (H-Ho)/Dt   - dXcdt*Z/Xc*Hz
		dCadt = (Ca-Cao)/Dt - dXcdt*Z/Xc*Caz
		dCdt  = (C-Co)/Dt   - dXcdt*Z/Xc*Cz
		dMdt  = (M-Mo)/Dt   - dXcdt*Z/Xc*Mz
		dCsdt = (Cs-Cso)/Dt - dXcdt*Z/Xc*Csz

		!-----------------------------------------------------------------------

		XX = Z*Xc/L1
		Hxxx = Hxxz*dZdx

		Cax  = Caz*dZdx
		Cx  = Cz*dZdx
		Mx  = Mz*dZdx
		Csx = Csz*dZdx

		S2  =  Surface_Tension('LA',Ca)
		S2x =  Surface_Tension_x('LA',Ca,Cax)

		Jcca = ka*( Ra*C*(1.D0-Ca) - Ca )
		Jccs = ks*( Rs*C*(1.D0-Cs) - Cs )
		Jcm  = kb*( C**nc - M )

		SELECT CASE(GEOMETRY)
		CASE('CARTESIAN')
			dPdX = - eps**2*Hxxx*(1.D0/As2+S2) - eps**2*Hxx*S2x
		CASE('CYLINDRICAL')
			dPdX = - eps**2*( Hxxx + Hxx/XX - Hx/XX**2 )*(1.D0/As2+S2) - &
					 eps**2*( Hxx + Hx/XX )*S2x
		END SELECT

		Q = - dPdX*(H/3.d0+bslip)*H**2 + S2x*H*(H/2.D0+bslip)
		U = - 0.5d0*dPdX*H**2 + H*S2x + bslip*(S2x-dPdX*H)
		Um = Q/H

		!-----------------------------------------------------------------------
		!				ITERATE OVER THE WEIGHT FUNCTIONS
		!-----------------------------------------------------------------------
					
		LOOP_RESIDUALS: DO IW = 1, NBF_Q
 
			BIFN = BFN_Q(IW,IG)
			DBIFN = DFDX_Q(IW)

			IF (SUPG) THEN
				SBIFN = BFN_Q(IW,IG) + Dt*DFDX_Q(IW)*dZdx
			ELSE
				SBIFN = BFN_Q(IW,IG)
			ENDIF

			!-----------------------------------------------------------------------
			!				FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
			!-----------------------------------------------------------------------

			SELECT CASE(GEOMETRY)

			CASE('CARTESIAN')

				TERM(IW,1) = ( dHdt*BIFN - Q*DBIFN*dZdx )*dXdz
				TERM(IW,2) = ( Hx - Hz*dZdx )*BIFN*dXdz
				TERM(IW,3) = ( Hxx - Hxz*dZdx )*BIFN*dXdz
				TERM(IW,4) = ( dCadt*SBIFN - ( U*Ca - Cax/Peca )*DBIFN*dZdx - &
								Jcca*SBIFN )*dXdz
				TERM(IW,5) = ( ( H*dCdt + H*Um*Cx + ba*Jcca + bs*Jccs + &
								H*Jcm )*SBIFN + DBIFN*dZdx*H*Cx/Peb )*dXdz
				TERM(IW,6) = ( ( H*dMdt + H*Um*Mx - H*Jcm )*SBIFN + &
								DBIFN*dZdx*H*Mx/Pem )*dXdz
				TERM(IW,7) = ( dCsdt*SBIFN + Csx/Pecs*DBIFN*dZdx - Jccs*SBIFN )*dXdz

			CASE('CYLINDRICAL')

				TERM(IW,1) = ( dHdt*BIFN - Q*DBIFN*dZdx )*XX*dXdz
				TERM(IW,2) = ( Hx - Hz*dZdx )*BIFN*XX*dXdz
				TERM(IW,3) = ( Hxx - Hxz*dZdx )*BIFN*XX*dXdz
				TERM(IW,4) = ( dCadt*SBIFN - ( U*Ca - Cax/Peca )*DBIFN*dZdx - & 
								Jcca*SBIFN )*XX*dXdz
				TERM(IW,5) = ( ( H*dCdt + H*Um*Cx + ba*Jcca + bs*Jccs + &
								H*Jcm )*SBIFN + DBIFN*dZdx*H*Cx/Peb )*XX*dXdz
				TERM(IW,6) = ( ( H*dMdt + H*Um*Mx - H*Jcm )*SBIFN + &
								DBIFN*dZdx*H*Mx/Pem )*XX*dXdz
				TERM(IW,7) = ( dCsdt*SBIFN + Csx/Pecs*DBIFN*dZdx - & 
								Jccs*SBIFN )*XX*dXdz

			END SELECT

		ENDDO  LOOP_RESIDUALS
	
	CASE(2)

		!-----------------------------------------------------------------------
		!			CALCULATE THE DEPENDENT VARIABLES AND THEIR PARTIAL
		!			DERIVATIVES AT THE GAUSS POINTS IN X COORDINATE
		!-----------------------------------------------------------------------

		Z = 0.D0

		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			Z    =  Z  + X(INOD) * BFN_Q(II,IG)
		ENDDO
  
		Cs  = 0.D0  ; Cso  = 0.D0  ; Csz = 0.D0
 
		DO II=1,NBF_Q
			INOD = NM_Q(NELEM,II+1)
			Cs  = Cs  + TQ(INOD,1) * BFN_Q(II,IG)
			Cso = Cso + TQo(INOD,1)* BFN_Q(II,IG)
			Csz = Csz + TQ(INOD,1) * DFDX_Q(II)
		ENDDO

		!-----------------------------------------------------------------------
		!	CALCULATE TIME DERIVATIVES
		!-----------------------------------------------------------------------

		dxdZ = L2*Xc/(L1+L2-Z)**2
		dZdx = (L1+L2-Z)**2/(L2*Xc)

		dXcdt = (Xc-Xco)/Dt
		dCsdt = (Cs-Cso)/Dt - dXcdt*(L1+L2-Z)/Xc*Csz

		!-----------------------------------------------------------------------

		XX = L2*Xc/(L1+L2-Z)
		Csx = Csz*dZdx

		!-----------------------------------------------------------------------
		!	ITERATE OVER THE WEIGHT FUNCTIONS
		!-----------------------------------------------------------------------

		LOOP_RESIDUALS2: DO IW = 1, NBF_Q

			BIFN = BFN_Q(IW,IG)
			DBIFN = DFDX_Q(IW)

			!-----------------------------------------------------------------------
			!		FORM THE WORKING RESIDUAL VECTOR IN ELEMENT NELEM
			!-----------------------------------------------------------------------

			SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
				TERM(IW,1) = ( dCsdt*BIFN + Csx/Pecs*DBIFN*dZdx )*dXdz
			CASE('CYLINDRICAL')
				TERM(IW,1) = ( dCsdt*BIFN + Csx/Pecs*DBIFN*dZdx )*XX*dXdz
			END SELECT

		ENDDO  LOOP_RESIDUALS2

	END SELECT

END SUBROUTINE EQUATIONS

!-----------------------------------------------------------------
!     SUBROUTINE BOUNDARY_RESIDUALS
!-----------------------------------------------------------------

SUBROUTINE BOUNDARY_RESIDUALS

	USE ELEMENTS_DATA, only: NBF_Q, NNTOLa, NNTOLb, NEQ_Q_a, NEQ_Q_b, &
						NODTOL, IDIM_A, NODTOLa_QEL, NORE
	USE COMMON_ARRAYS
	USE FLOW_PARAMETERS
	USE GAUSS_DATA
	USE BASIS_FUNCTION_1D
	USE MATRIX_STORAGE_SUBROUTINES
	USE TIME_INTEGRATION, only: Dt
	USE BASIC_PARAMETERS, only: FORM_JAC, L1, L2
	USE FUNCTIONS

!-----------------------------------------------------------------

	IMPLICIT NONE

	INTEGER :: IEL, I, IN, IEDGE, IROW, J, ROI
	INTEGER :: JW, KW, IMOD, II, ICOL, IEQ, INOD

	INTEGER, PARAMETER :: IDIM_Q = NBF_Q

	REAL(8) :: Z, XX, AJACX_Q, dXdZ, dZdx
	REAL(8) :: Ca, Caz, Cax, Cax_Xc
	REAL(8) :: C,  Cz,  Cx,  Cx_Xc
	REAL(8) :: M,  Mz,  Mx,  Mx_Xc
	REAL(8) :: Cs, Csz, Csx, Csx_Xc
	REAL(8) :: dZdx_Xc, dXcdt, dXcdt_Xc, XX_Xc
	REAL(8) :: Jcacs, Jcca, Jccs

	REAL(8) :: Theta, Theta_Hx, Thetaa_Ca, Thetaa_Cs, Caj, Csj, epert

	REAL(8) :: BJFN, DBJX

	INTEGER, DIMENSION(2) :: IERROR

	INTEGER, DIMENSION(IDIM_Q) :: NQ, NQ_QL

	REAL(8), DIMENSION(NBF_Q) :: Cax_Ca, Cax_Cs, Csx_Cs
	REAL(8), DIMENSION(NBF_Q) :: Cx_C, Cx_Ca, Cx_Cs, Mx_M
	REAL(8), DIMENSION(NBF_Q) :: Jcacs_Ca, Jcacs_Cs 
	REAL(8), DIMENSION(NBF_Q) :: Jcca_C, Jcca_Ca
	REAL(8), DIMENSION(NBF_Q) :: Jccs_C, Jccs_Cs
	REAL(8), DIMENSION(NBF_Q) :: DFDX_Q

	REAL(8), ALLOCATABLE, DIMENSION(:,:) :: TEMP_X
	REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: TEMP

	!-----------------------------------------------------------------------
	!     ITERATE OVER ALL REGIONS
	!-----------------------------------------------------------------------

	LOOP_REGIONS: DO ROI = 1, NORE

		SELECT CASE(ROI)

		CASE(1)

			!ALLOCATE TEMPORARY ARRAYS
			ALLOCATE( TEMP_X(IDIM_Q,NEQ_Q_a)			,  STAT=IERROR(1) )
			ALLOCATE( TEMP(IDIM_Q,IDIM_Q,NEQ_Q_a,NEQ_Q_a), STAT=IERROR(2) )

			LOOP_CHECK_ALLOCATION: DO I = 1, 2
				IF(IERROR(I) .NE. 0)THEN
					WRITE(42,*)I,'ARRAY IS NOT ALLOCATED / SUB: BOUNDARY_RESIDUALS //1'
					WRITE(42,*)'STAT =', IERROR
					WRITE(42,*)'PROGRAM STOP !!!'
					STOP
				ENDIF
			ENDDO LOOP_CHECK_ALLOCATION

			!-----------------------------------------------------------------------
			!		ITERATE OVER THE 2 EDGES OF REGION 1
			!-----------------------------------------------------------------------
			LOOP_EDGEa: DO IEDGE = 1,2
				SELECT CASE (IEDGE)
				CASE (1)
					IN  = 1
					IEL = 1
					NQ    = NM_Q_a(IEL,2:)
					NQ_QL = NM_QQ_a(IEL,2:)
				CASE (2)
					IN  = 3
					IEL = NNTOLa
					NQ    = NM_Q_a(IEL,2:)
					NQ_QL = NM_QQ_a(IEL,2:)
				END SELECT

				!--------------------------------------------------------------------
				!			INITIALIZE ARRAYS
				!--------------------------------------------------------------------
				TEMP = 0.D0 ; TEMP_X = 0.D0

				!--------------------------------------------------------------------
				!			FIND JACOBIAN OF THE TRANSFORMATION E - X
				!--------------------------------------------------------------------

				AJACX_Q = 0.D0

				DO I = 1, IDIM_Q
					AJACX_Q = AJACX_Q + Xa(NQ(I))*DNF_Q(I,IN)
				ENDDO

				AJACX_Q = DABS(AJACX_Q)

				DFDX_Q(:) = DNF_Q(:,IN) / AJACX_Q

				!--------------------------------------------------------------------
				!			CALCULATE VARIABLES
				!--------------------------------------------------------------------

				Z  = 0.D0
				Ca = 0.D0 ; Caz = 0.d0
				C  = 0.D0 ; Cz  = 0.d0
				M  = 0.D0 ; Mz  = 0.d0
				Cs = 0.D0

                !7 equations of interest (Z,H,Hx,Hxx,Ca,M,Cs)
				DO I = 1, IDIM_Q
					Z   =  Z    +  Xa(NQ(I))*NF_Q(I,IN)
					Ca  =  Ca   +  TQa(NQ(I),4)*NF_Q(I,IN)
					Caz =  Caz  +  TQa(NQ(I),4)*DFDX_Q(I)
					C   =  C    +  TQa(NQ(I),5)*NF_Q(I,IN)
					Cz  =  Cz   +  TQa(NQ(I),5)*DFDX_Q(I)
					M   =  M    +  TQa(NQ(I),6)*NF_Q(I,IN)
					Mz  =  Mz   +  TQa(NQ(I),6)*DFDX_Q(I)
					Cs  =  Cs   +  TQa(NQ(I),7)*NF_Q(I,IN)
				ENDDO

				XX = Z*Xc/L1

				dZdx = L1/Xc

				dXcdt = (Xc-Xco)/Dt

				Cax = Caz*dZdx
				Cx  = Cz*dZdx
				Mx  = Mz*dZdx

				Jcacs = kas*( Ras*Cs*(1.D0-Ca) - Ca*(1.D0-Cs) )
				Jcca  = ka*( Ra*C*(1.D0-Ca) - Ca )
				Jccs  = Ks*( Rs*C*(1.D0-Cs) - Cs )

				XX_Xc = Z/L1

				dZdx_Xc = - L1/Xc**2

				dXcdt_Xc = 1.D0/Dt

				Cax_Ca(:) =  DFDX_Q(:)*dZdx
				Cax_Cs(:) =  0.D0
				Cax_Xc    =  Caz*dZdx_Xc

				Cx_C(:)   =  DFDX_Q(:)*dZdx
				Cx_Ca(:)  =  0.D0
				Cx_Cs(:)  =  0.D0
				Cx_Xc     =  Cz*dZdx_Xc

				Mx_M(:)   =  DFDX_Q(:)*dZdx
				Mx_Xc     =  Mz*dZdx_Xc

				Jcacs_Ca(:) = - kas*Ras*Cs*NF_Q(:,IN) - kas*NF_Q(:,IN)*(1.D0-Cs)
				Jcacs_Cs(:) =   kas*Ras*NF_Q(:,IN)*(1.D0-Ca) + kas*Ca*NF_Q(:,IN)

				Jcca_Ca(:)  = - ka*( Ra*C + 1.D0 )*NF_Q(:,IN)
				Jcca_C(:)   =   ka*Ra*NF_Q(:,IN)*(1.D0-Ca)

				Jccs_Cs(:)  = - ks*( Rs*C + 1.D0 )*NF_Q(:,IN)
				Jccs_C(:)   =   ks*Rs*NF_Q(:,IN)*(1.D0-Cs)

				!-----------------------------------------------------------------------
				!			CALCULATE RESIDUALS
				!-----------------------------------------------------------------------

				SELECT CASE (IEDGE)

				CASE (1)

					IROW = NQ_QL(IN)

					B(IROW+3) = B(IROW+3) ! Cax = 0, u = 0 at x = 0
					B(IROW+4) = B(IROW+4) ! Cx = 0, u = 0 at x = 0
					B(IROW+5) = B(IROW+5) ! Mx = 0, u = 0 at x = 0
					B(IROW+6) = B(IROW+6) ! Csx = 0, u = 0 at x = 0
				
				CASE (2)

					IROW = NQ_QL(IN)

					SELECT CASE (GEOMETRY)

					CASE ('CYLINDRICAL')

						B(IROW+3) = B(IROW+3) + XX*( dXcdt*Ca - Jcacs )
						B(IROW+4) = B(IROW+4)					! Cx = 0 
						B(IROW+5) = B(IROW+5)					! Mx = 0
						B(IROW+6) = B(IROW+6) + XX*ba/bs*Jcacs

						LOOP_JACOBIAN2: DO J = 1, IDIM_Q

							BJFN=NF_Q(J,IN)

							TEMP(IN,J,4,4) =   XX*( dXcdt*NF_Q(J,IN) - Jcacs_Ca(J) ) 
							TEMP(IN,J,7,4) = - XX*Jcacs_Cs(J)
							TEMP(IN,J,4,7) =   XX*ba/bs*Jcacs_Ca(J)
							TEMP(IN,J,7,7) =   XX*ba/bs*Jcacs_Cs(J)

						ENDDO LOOP_JACOBIAN2

						TEMP_X(IN,4) =   XX_Xc*( dXcdt*NF_Q(J,IN) - Jcacs_Ca(J) ) + &
								XX*dXcdt_Xc*Ca
						TEMP_X(IN,7) =   XX_Xc*ba/bs*Jcacs

					CASE ('CARTESIAN')

						B(IROW+3) = B(IROW+3) + dXcdt*Ca - Jcacs
						B(IROW+4) = B(IROW+4)							! Cx = 0 
						B(IROW+5) = B(IROW+5)							! Mx = 0
						B(IROW+6) = B(IROW+6) + ba/bs*Jcacs

						LOOP_JACOBIAN3: DO J = 1, IDIM_Q

							BJFN=NF_Q(J,IN)

							TEMP(IN,J,4,4) =   dXcdt*NF_Q(J,IN) - Jcacs_Ca(J)
							TEMP(IN,J,7,4) = - Jcacs_Cs(J)
							TEMP(IN,J,4,7) =   ba/bs*Jcacs_Ca(J)
							TEMP(IN,J,7,7) =   ba/bs*Jcacs_Cs(J)

						ENDDO LOOP_JACOBIAN3

						TEMP_X(IN,4) =   dXcdt_Xc*Ca
						TEMP_X(IN,7) =   0.D0

					END SELECT

					IF (FORM_JAC) THEN
 
						CALL MATRIX_STORAGE_JACOBIAN_BAND&
							( TEMP, NM_QQ_a(IEL,2:), NM_QQ_a(IEL,2:), NBF_Q, NBF_Q, &
									NEQ_Q_a, NEQ_Q_a, A, IDIM_A, NODTOL, IBAND )

						! STORE MATRIX A_ip

						DO II = 1, NBF_Q
						DO IEQ = 1, NEQ_Q_a
							ICOL = NM_QQ_a(IEL,II+1) + IEQ - 1
							A_ip(ICOL) = A_ip(ICOL) + TEMP_X(II,IEQ)
							ICOL = NM_QQ_a(IEL,II+1) + IEQ
						ENDDO
						ENDDO

					ENDIF

					!-----------------------------------------------------------------------
					!				CONTACT LINE
					!-----------------------------------------------------------------------

					Theta  = - TQa(NODTOLa_QEL,2) 
					Thetaa =   Tha(Ca,Cs)

					Theta_Hx = - 1.D0

					Caj = Ca ; CALL XPERT(Caj, epert)
					Csj = Cs ; CALL XPERT(Csj, epert)

					Thetaa_Ca = (Tha(Caj,Cs)-Thetaa)/epert
					Thetaa_Cs = (Tha(Ca,Csj)-Thetaa)/epert

					B_p = dXcdt - kappa*(Theta-Thetaa)**mm

					IF (FORM_JAC) THEN
						A_p = dXcdt_Xc
	
						IEQ = 2
						INOD = NM_QQ_a(NNTOLa,4) + IEQ - 1
						A_pi(INOD) = - mm*kappa*Theta_Hx*(Theta-Thetaa)**(mm-1)

						IEQ = 4
						INOD = NM_QQ_a(NNTOLa,4) + IEQ - 1
						A_pi(INOD) = mm*kappa*Thetaa_Ca*(Theta-Thetaa)**(mm-1)

						IEQ = 7
						INOD = NM_QQ_a(NNTOLa,4) + IEQ - 1
						A_pi(INOD) = mm*kappa*Thetaa_Cs*(Theta-Thetaa)**(mm-1)
					ENDIF

				END SELECT

			ENDDO LOOP_EDGEa

		CASE(2)

			! ALLOCATE TEMPORARY ARRAYS

			ALLOCATE( TEMP_X(IDIM_Q,NEQ_Q_b)			, STAT=IERROR(1) )
			ALLOCATE( TEMP(IDIM_Q,IDIM_Q,NEQ_Q_b,NEQ_Q_b), STAT=IERROR(2) )

			LOOP_CHECK_ALLOCATION2: DO I = 1, 2
				IF(IERROR(I) .NE. 0)THEN
					WRITE(42,*)I,'ARRAY IS NOT ALLOCATED / SUB: BOUNDARY_RESIDUALS //2'
					WRITE(42,*)'STAT =', IERROR
					WRITE(42,*)'PROGRAM STOP !!!'
					STOP
				ENDIF
			ENDDO LOOP_CHECK_ALLOCATION2

			!-----------------------------------------------------------------------
			!		ITERATE OVER THE 2 EDGES OF REGION 2
			!-----------------------------------------------------------------------

			LOOP_EDGEb: DO IEDGE = 1,2

				SELECT CASE (IEDGE)
 				CASE (1)
					IN  = 1
					IEL = 1
					NQ    = NM_Q_b(IEL,2:)
					NQ_QL = NM_QQ_b(IEL,2:)
				CASE (2)
					IN  = 3
					IEL = NNTOLb
					NQ    = NM_Q_b(IEL,2:)
				    NQ_QL = NM_QQ_b(IEL,2:)
				END SELECT

				!--------------------------------------------------------------------
				!			INITIALIZE ARRAYS
				!--------------------------------------------------------------------

				TEMP = 0.D0 ; TEMP_X = 0.D0

				!--------------------------------------------------------------------
				!			FIND JACOBIAN OF THE TRANSFORMATION E - X
				!--------------------------------------------------------------------

				AJACX_Q = 0.D0

				DO I = 1, IDIM_Q
					AJACX_Q = AJACX_Q + Xb(NQ(I))*DNF_Q(I,IN)
				ENDDO

				AJACX_Q = DABS(AJACX_Q)
				DFDX_Q(:) = DNF_Q(:,IN) / AJACX_Q

				!--------------------------------------------------------------------
				!			CALCULATE VARIABLES
				!--------------------------------------------------------------------

				Z  = 0.D0 ; Cs = 0.D0 ; Csz = 0.D0

				DO I = 1, IDIM_Q
					Z   =  Z    +  Xb(NQ(I))*NF_Q(I,IN)
					Cs  =  Cs   +  TQb(NQ(I),1)*NF_Q(I,IN)
					Caz =  Caz  +  TQb(NQ(I),1)*DFDX_Q(I)
				ENDDO

                if (L1+L2-Z .gt. 1e-7) then
				    XX = L2*Xc/(L1+L2-Z)
                else
                    XX = 0.d0
                endif
				dZdx = (L1+L2-Z)**2/(L2*Xc)
				dXcdt = (Xc-Xco)/Dt
				Csx = Csz*dZdx
                if (L1+L2-Z .gt. 1e-7) then
				    XX_Xc = L2/(L1+L2-Z)
                else
                    XX_Xc = 0.d0
                endif
				dZdx_Xc = - (L1+L2-Z)**2/(L2*Xc**2)
				dXcdt_Xc = 1.D0/Dt
				Csx_Cs(:) = DFDX_Q(:)*dZdx
				Csx_Xc    = Csz*dZdx_Xc

				!-----------------------------------------------------------------------
				!			CALCULATE RESIDUALS
				!-----------------------------------------------------------------------

				SELECT CASE (IEDGE)

				CASE (1)

					IROW = NQ_QL(IN)
					B(IROW) = B(IROW) ! I have added the term ba/bs*Jcacs above
				
				CASE (2)

					IROW = NQ_QL(IN)
					Csx  = 0.D0
					Csx_Cs(:) = 0.D0
					Csx_Xc    = 0.D0

					SELECT CASE (GEOMETRY)

					CASE ('CYLINDRICAL')

!						B(IROW) = B(IROW) - XX*Csx/Pecs
						B(IROW) = B(IROW) ! XX becomes infinite at L1+L2 and Csx(x->oo)=0

						LOOP_JACOBIAN4: DO J = 1, IDIM_Q

							BJFN=NF_Q(J,IN)
					!		TEMP(IN,J,1,1) = - XX*Csx_Cs(J)/Pecs
							TEMP(IN,J,1,1) = 0.D0

						ENDDO LOOP_JACOBIAN4

!						TEMP_X(IN,1) = - (XX_Xc*Csx + XX*Csx_Xc)/Pecs
						TEMP_X(IN,1) = 0.D0

					CASE ('CARTESIAN')

						B(IROW) = B(IROW) - Csx/Pecs

						LOOP_JACOBIAN5: DO J = 1, IDIM_Q

							BJFN=NF_Q(J,IN)
							TEMP(IN,J,1,1) = - Csx_Cs(J)/Pecs

						ENDDO LOOP_JACOBIAN5

						TEMP_X(IN,1) = - Csx_Xc/Pecs

					END SELECT

					IF (FORM_JAC) THEN
 
						CALL MATRIX_STORAGE_JACOBIAN_BAND&
								( TEMP, NM_QQ_b(IEL,2:), NM_QQ_b(IEL,2:), NBF_Q, NBF_Q, &
										NEQ_Q_b, NEQ_Q_b, A, IDIM_A, NODTOL, IBAND )

						! STORE MATRIX A_ip

						DO II = 1, NBF_Q
						DO IEQ = 1, NEQ_Q_b
							ICOL = NM_QQ_b(IEL,II+1) + IEQ - 1
							A_ip(ICOL) = A_ip(ICOL) + TEMP_X(II,IEQ)
							ICOL = NM_QQ_b(IEL,II+1) + IEQ
						ENDDO
						ENDDO

					ENDIF

				END SELECT

			ENDDO LOOP_EDGEb

		END SELECT

!-----------------------------------------------------------------------
!	DEALLOCATE TEMPORARY ARRAYS
!-----------------------------------------------------------------------

		DEALLOCATE( TEMP_X, STAT = IERROR(1) )
		DEALLOCATE( TEMP,   STAT = IERROR(2) )

		LOOP_CHECK_DEALLOCATION: DO I = 1, 2
			IF(IERROR(I) .NE. 0)THEN
				WRITE(42,*)I,'ARRAY IS NOT DEALLOCATED / SUB: BOUNDARY_RESIDUALS '
				WRITE(42,*)'STAT =', IERROR
				WRITE(42,*)'PROGRAM STOP !!!'
				STOP
			ENDIF
		ENDDO LOOP_CHECK_DEALLOCATION

	ENDDO LOOP_REGIONS

!-----------------------------------------------------------------------

END SUBROUTINE BOUNDARY_RESIDUALS

!----------------------------------------------------------------------
!						SUBROUTINE CALCULATE_MASS
!----------------------------------------------------------------------

SUBROUTINE CALCULATE_MASS(TQa, TQb, Xc, MASS_FLUID, MASS_SURFACTANT)

	USE ELEMENTS_DATA, only: NBF_Q, NNTOLa, NNTOLb, NODTOLa_QEL, &
						NODTOLb_QEL, NEQ_Q_a, NEQ_Q_b
	USE COMMON_ARRAYS, only: NM_Q_a, NM_Q_b, Xa, Xb
	USE GAUSS_DATA
	USE BASIS_FUNCTION_1D
	USE FLOW_PARAMETERS, only: GEOMETRY, SOLUBLE_SURFACTANT, ba, bs, &
								WALL_ADSORPTION
	USE BASIC_PARAMETERS, only: L1, L2

	IMPLICIT NONE

	INTEGER :: NELEM, I, IG, KK, IW, JW, INOD
	REAL(8) :: AJACX_Q, WET_Q
	REAL(8) :: Z, H, Ca, C, M, Cs, XX, dxdZ, dZdx

	REAL(8), DIMENSION(NBF_Q) :: DFDX_Q

	REAL(8), INTENT(INOUT) :: MASS_FLUID, MASS_SURFACTANT

    !Arrays which contain all information for all the  
    !nodes in domain (NODTOL*_QEL) and each of the 7 equations (NEQ_Q_*)
	REAL(8), INTENT(IN), DIMENSION (NODTOLa_QEL,NEQ_Q_a) :: TQa
	REAL(8), INTENT(IN), DIMENSION (NODTOLb_QEL,NEQ_Q_b) :: TQb
	REAL(8), INTENT(IN) :: Xc

	!----------------------------------------------------------------------
	!     INITIALIZE MASS_FLUID AND MASS_SURFACTANT
	!----------------------------------------------------------------------

	MASS_FLUID = 0.D0 ; MASS_SURFACTANT = 0.D0

	!----------------------------------------------------------------------
	!     CALCULATE FLUID AND SURFACTANT MASS IN REGION 1
	!----------------------------------------------------------------------

	LOOP_ELEMENTS: DO NELEM = 1, NNTOLa

		!----------------------------------------------------------------------
		!			ITERATE OVER EACH GAUSS POINT
		!----------------------------------------------------------------------

		LOOP_GAUSS: DO IG = 1, NINT_3P

			!--------------------------------------------------------------------
			!     FIND JACOBIAN OF THE TRANSFORMATION E - X
			!--------------------------------------------------------------------

			AJACX_Q = 0.D0

			DO I = 1, NBF_Q
				INOD = NM_Q_a(NELEM,I+1)
				AJACX_Q = AJACX_Q + Xa(INOD)*DTF_Q(I,IG)
			ENDDO

			AJACX_Q = DABS(AJACX_Q)
			WET_Q = WO_3P(IG) * AJACX_Q
			DFDX_Q(:) = DTF_Q(:,IG) / AJACX_Q

			!----------------------------------------------------------------------
			!   CALCULATE THE DEPENDENT VARIABLES AND THEIR PARTIAL
			!   DERIVATIVES AT THE GAUSS POINTS IN X COORDINATE
			!----------------------------------------------------------------------

			Z = 0.D0 ; H = 0.D0 ; Ca = 0.D0 ; C = 0.D0 ; M = 0.D0 ; Cs = 0.D0

			DO I=1,NBF_Q
				INOD = NM_Q_a(NELEM,I+1)
				Z    =  Z  + Xa(INOD)    * TF_Q(I,IG)
				H    =  H  + TQa(INOD,1) * TF_Q(I,IG)
				Ca   =  Ca + TQa(INOD,4) * TF_Q(I,IG)
				C    =  C  + TQa(INOD,5) * TF_Q(I,IG)
				M    =  M  + TQa(INOD,6) * TF_Q(I,IG)
				Cs   =  Cs + TQa(INOD,7) * TF_Q(I,IG)
			ENDDO

			dxdZ = Xc/L1
			dZdx = L1/Xc

			XX = Z*Xc/L1

			!-----------------------------------------------------------------------
			!		CALCULATE TOTAL MASS
			!-----------------------------------------------------------------------

			SELECT CASE(GEOMETRY)

			CASE('CARTESIAN')

				MASS_FLUID = MASS_FLUID + H*dxdZ*WET_Q

				IF (SOLUBLE_SURFACTANT) THEN
					MASS_SURFACTANT = MASS_SURFACTANT + &
								( ba*Ca + bs*Cs + H*C + H*M )*dxdZ*WET_Q
				ELSE
					MASS_SURFACTANT = MASS_SURFACTANT + (Ca+bs/ba*Cs)*dxdZ*WET_Q
				ENDIF					

			CASE('CYLINDRICAL')

				MASS_FLUID = MASS_FLUID + H*XX*dxdZ*WET_Q

				IF (SOLUBLE_SURFACTANT) THEN
					MASS_SURFACTANT = MASS_SURFACTANT + &
					( ba*Ca + bs*Cs + H*C + H*M )*XX*dxdZ*WET_Q
				ELSE
					MASS_SURFACTANT = MASS_SURFACTANT + (Ca+bs/ba*Cs)*XX*dxdZ*WET_Q
				ENDIF					

			END SELECT

		ENDDO LOOP_GAUSS

	ENDDO LOOP_ELEMENTS

	!----------------------------------------------------------------------
	!     FLUID AND SURFACTANT MASS IN REGION 2
	!----------------------------------------------------------------------

	LOOP_ELEMENTS2: DO NELEM = 1, NNTOLb

		!----------------------------------------------------------------------
		!			ITERATE OVER EACH GAUSS POINT
		!----------------------------------------------------------------------

		LOOP_GAUSS2: DO IG = 1, NINT_3P

			!--------------------------------------------------------------------
			!     FIND JACOBIAN OF THE TRANSFORMATION E - X
			!--------------------------------------------------------------------

			AJACX_Q = 0.D0

			DO I = 1, NBF_Q
				INOD = NM_Q_b(NELEM,I+1)
				AJACX_Q = AJACX_Q + Xb(INOD)*DTF_Q(I,IG)
			ENDDO

			AJACX_Q = DABS(AJACX_Q)
			WET_Q = WO_3P(IG) * AJACX_Q
			DFDX_Q(:) = DTF_Q(:,IG) / AJACX_Q

			!----------------------------------------------------------------------
			!  CALCULATE THE DEPENDENT VARIABLES AND THEIR PARTIAL
			!  DERIVATIVES AT THE GAUSS POINTS IN X COORDINATE
			!----------------------------------------------------------------------

			Z = 0.D0 ; Cs = 0.D0

			DO I=1,NBF_Q
				INOD = NM_Q_b(NELEM,I+1)
				Z    =  Z  + Xb(INOD)    * TF_Q(I,IG)
				Cs   =  Cs + TQb(INOD,1) * TF_Q(I,IG)
			ENDDO

			dxdZ = L2*Xc/(L1+L2-Z)**2
			dZdx = (L1+L2-Z)**2/(L2*Xc)

			XX = (L2*Xc)/(L1+L2-Z)

			!-----------------------------------------------------------------------
			!		CALCULATE TOTAL MASS
			!-----------------------------------------------------------------------

			SELECT CASE(GEOMETRY)

			CASE('CARTESIAN')

				IF (SOLUBLE_SURFACTANT) THEN
					MASS_SURFACTANT = MASS_SURFACTANT + bs*Cs*dxdZ*WET_Q
				ELSE
					MASS_SURFACTANT = MASS_SURFACTANT + bs/ba*Cs*dxdZ*WET_Q
				ENDIF

			CASE('CYLINDRICAL')

				IF (SOLUBLE_SURFACTANT) THEN
					MASS_SURFACTANT = MASS_SURFACTANT + bs*Cs*XX*dxdZ*WET_Q
				ELSE
					MASS_SURFACTANT = MASS_SURFACTANT + bs/ba*Cs*XX*dxdZ*WET_Q
				ENDIF

			END SELECT

		ENDDO LOOP_GAUSS2

	ENDDO LOOP_ELEMENTS2

END SUBROUTINE CALCULATE_MASS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!			SUBROUTINE CREATE_GRAPH_ARRAYS
!-----------------------------------------------------------------------

SUBROUTINE CREATE_GRAPH_ARRAYS(TQa,TQb,Xc)

	USE ELEMENTS_DATA,     only: NBF_Q, NNTOLa, NNTOLb, NEQ_Q_a, NEQ_Q_b, &
										NODTOLa_QEL, NODTOLb_QEL, NZELa
	USE COMMON_ARRAYS,     only: Xa, Xb, NM_Q_a, NM_Q_b, TQo_GRAPH
	USE COMMON_ARRAYS,     only: uwarray, Z_loc
	USE BASIS_FUNCTION_1D, only: NF_Q, DNF_Q
	USE FLOW_PARAMETERS
	USE FUNCTIONS
	USE BASIC_PARAMETERS, only: L1, L2

	IMPLICIT NONE

	REAL(8), intent(in)                                 :: Xc
	REAL(8), DIMENSION(NODTOLa_QEL,NEQ_Q_a), intent(in) :: TQa
	REAL(8), DIMENSION(NODTOLb_QEL,NEQ_Q_b), intent(in) :: TQb

	INTEGER :: IEL, I, J, INOD
	REAL(8) :: Z, XX, dZdx, H, Hx, Hxx, Hxxx, Hxxxx
	REAL(8) :: Ca, Caz, Cax, C, Cz, Cx, Caxz, Caxx
	REAL(8) :: M, Mz, Mx, Cs, Csz, Csx, COEFF
	REAL(8) :: S1, S1x, S2, S2x, S12, S12x, S2xx
	REAL(8) :: dPdX, U, Hxxz, AJACX_Q, Hxxxz, Pxx

	REAL(8), DIMENSION(NBF_Q) :: DFDX_Q
	INTEGER, DIMENSION(5) :: IERROR

    integer :: zpoint

	!-----------------------------------------------------------------------
	!     INITIALIZE ARRAY
	!-----------------------------------------------------------------------
	TQo_GRAPH = 0.D0
    if (.not. allocated(Z_loc)) allocate(Z_loc(NODTOLa_QEL,NZELa))
    if (.not. allocated(uwarray)) allocate(uwarray(NODTOLa_QEL,NZELa,2))
    Z_loc = 0.d0; uwarray = 0.d0

	!-----------------------------------------------------------------------
	!     REGION 1
	!-----------------------------------------------------------------------
	LOOP_ELEMENTS1: DO IEL = 1, NNTOLa

		DO I = 1, NBF_Q

			! FIND JACOBIAN OF THE TRANSFORMATION E - X

			AJACX_Q = 0.D0

			DO J = 1, NBF_Q
				INOD = NM_Q_a(IEL,J+1)
				AJACX_Q = AJACX_Q + Xa(INOD)*DNF_Q(J,I)
			ENDDO

			AJACX_Q = DABS(AJACX_Q)

			DFDX_Q(:) = DNF_Q(:,I) / AJACX_Q

			! CALCULATE VARIABLES

			Z    = 0.D0
			H    = 0.D0 ; Hx   = 0.D0
			Hxx  = 0.D0 ; Hxxz = 0.D0
			Ca   = 0.D0 ; Caz  = 0.D0
			C    = 0.D0 ; Cz   = 0.D0
			M    = 0.D0 ; Mz   = 0.D0
			Cs   = 0.D0 ; Csz  = 0.D0

			DO J = 1, NBF_Q
                !Get node index from NM_Q_a, get value at a given node from TQa
                !which contains 7 quantities of interest (Z,H,Hx,Hxx,Ca,M,Cs) 
                !and multiply by Gauss points to integrate
				INOD = NM_Q_a(IEL,J+1)
				Z    = Z    + Xa(INOD)    * NF_Q(J,I)
				H    = H    + TQa(INOD,1) * NF_Q(J,I)
				Hx   = Hx   + TQa(INOD,2) * NF_Q(J,I)
				Hxx  = Hxx  + TQa(INOD,3) * NF_Q(J,I)
				Hxxz = Hxxz + TQa(INOD,3) * DFDX_Q(J)
				Ca   = Ca   + TQa(INOD,4) * NF_Q(J,I)
				Caz  = Caz  + TQa(INOD,4) * DFDX_Q(J)
				C    = C    + TQa(INOD,5) * NF_Q(J,I)
				Cz   = Cz   + TQa(INOD,5) * DFDX_Q(J)
				M    = M    + TQa(INOD,6) * NF_Q(J,I)
				Mz   = Mz   + TQa(INOD,6) * DFDX_Q(J)
				Cs   = Cs   + TQa(INOD,7) * NF_Q(J,I)
				Csz  = Csz  + TQa(INOD,7) * DFDX_Q(J)
			ENDDO

            !x position of node in mapped coordinate system
			XX = Z*Xc/L1

			dZdx = L1/Xc
   
            !Convert from z coordinate derivative to x?
			Hxxx = Hxxz*dZdx
			Cax  = Caz*dZdx
			Cx   = Cz*dZdx
			Mx   = Mz*dZdx
			Csx  = Csz*dZdx

			S12  = Surface_Tension('LS',Cs)
			S12x = Surface_Tension_x('LS',Cs,Csx)

			S2  = Surface_Tension('LA',Ca)
			S2x = Surface_Tension_x('LA',Ca,Cax)

			SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
				dPdX = - eps**2*Hxxx*(1.D0/As2+S2) - eps**2*Hxx*S2x
			CASE('CYLINDRICAL')
				IF ( Z < 1.d-9) THEN
					dPdX = 0.D0
				ELSE
					dPdX = - eps**2*( Hxxx + Hxx/XX - Hx/XX**2 )*(1.D0/As2+S2) - &
							 eps**2*( Hxx + Hx/XX )*S2x
				ENDIF
			END SELECT

			!Note, this is equivalent to Eq 2.32 where 
            !U(z=H) = 0.5d0*dPdX*H**2 + (S2x - dPdX*H)*H + bslip*(S2x-dPdX*H)
			U = - 0.5d0*dPdX*H**2 + H*S2x + bslip*(S2x-dPdX*H)

            !Reinflate the z dimension using equation (2.32) from Karapetsas et al JFM (2011)
            Hxxxz = 0.D0; Caxz = 0.d0
		    DO J=1,NBF_Q
			    Hxxxz = Hxxxz + Hxxx * DFDX_Q(J)
				Caxz  = Caxz  + Cax  * DFDX_Q(J)
		    ENDDO
            Hxxxx = Hxxxz*dZdx
			Caxx  = Caxz*dZdx
            S2xx = Surface_Tension_xx('LA',Ca,Cax,Caxx)
			SELECT CASE(GEOMETRY)
			CASE('CARTESIAN')
                Pxx = -(eps**2) * (Hxxxx*(1.D0/As2+S2) + 2.d0*Hxxx*S2x + Hxx*S2xx )
			CASE('CYLINDRICAL')
                stop "Error in CREATE_GRAPH_ARRAYS -- 2D velocity output not developed for CYLINDRICAL geometry"
			END SELECT

            COEFF = S2xx - Pxx*H - dPdx * Hx
            do zpoint = 1,NZELa
			    INOD = NM_Q_a(IEL,I+1)
                Z_loc(INOD,zpoint)  =  H * zpoint/dble(NZELa)
                uwarray(INOD,zpoint,1) =  0.5d0*dPdX*Z_loc(INOD,zpoint)**2 & 
                                         + Z_loc(INOD,zpoint)*(S2x-dPdX*H) &
                                         + bslip*(S2x-dPdX*H)

                uwarray(INOD,zpoint,2) =-( (1.d0/6.d0)*Pxx * Z_loc(INOD,zpoint)**3 &
                                          + 0.5d0*COEFF*Z_loc(INOD,zpoint)**2 &   
                                          + bslip*COEFF*Z_loc(INOD,zpoint)             )
            enddo

			INOD = NM_Q_a(IEL,I+1)

			TQo_GRAPH(INOD,1)  = Z
			TQo_GRAPH(INOD,2)  = XX
			TQo_GRAPH(INOD,3)  = H
			TQo_GRAPH(INOD,4)  = Hx
			TQo_GRAPH(INOD,5)  = Hxx
			TQo_GRAPH(INOD,6)  = Hxxx
			TQo_GRAPH(INOD,7)  = U
			TQo_GRAPH(INOD,8)  = Ca
			TQo_GRAPH(INOD,9)  = Cax
			TQo_GRAPH(INOD,10) = C
			TQo_GRAPH(INOD,11) = Cx
			TQo_GRAPH(INOD,12) = M
			TQo_GRAPH(INOD,13) = Mx
			TQo_GRAPH(INOD,14) = Cs
			TQo_GRAPH(INOD,15) = Csx
			TQo_GRAPH(INOD,16) = S2
			TQo_GRAPH(INOD,17) = S2x
			TQo_GRAPH(INOD,18) = S12
			TQo_GRAPH(INOD,19) = S12x

		ENDDO

	ENDDO LOOP_ELEMENTS1

!-----------------------------------------------------------------------
!     REGION 2
!-----------------------------------------------------------------------

	LOOP_ELEMENTS2: DO IEL = 1, NNTOLb

		DO I = 2, NBF_Q

			! FIND JACOBIAN OF THE TRANSFORMATION E - X

			AJACX_Q = 0.D0

			DO J = 1, NBF_Q
				INOD = NM_Q_b(IEL,J+1)
				AJACX_Q = AJACX_Q + Xb(INOD)*DNF_Q(J,I)
			ENDDO

			AJACX_Q = DABS(AJACX_Q)

			DFDX_Q(:) = DNF_Q(:,I) / AJACX_Q

			! CALCULATE VARIABLES

			Z = 0.D0 ; Cs = 0.D0 ; Csz = 0.D0

			DO J = 1, NBF_Q
				INOD = NM_Q_b(IEL,J+1)
				Z    = Z    + Xb(INOD)    * NF_Q(J,I)
				Cs   = Cs   + TQb(INOD,1) * NF_Q(J,I)
				Csz  = Csz  + TQb(INOD,1) * DFDX_Q(J)
			ENDDO

            if (L1+L2-Z .gt. 1e-7) then
    			XX = L2*Xc/(L1+L2-Z)
            else
                XX = 0.d0
            endif

			dZdx = (L1+L2-Z)**2/(L2*Xc)

			Csx   = Csz*dZdx

			S1  = Surface_Tension('SA',Cs)
			S1x = Surface_Tension_x('SA',Cs,Csx)

			INOD = NM_Q_b(IEL,I+1) + NODTOLa_QEL - 1

			TQo_GRAPH(INOD,1)     = Z
			TQo_GRAPH(INOD,2)     = XX
			TQo_GRAPH(INOD,3:13)  = 0.D0
			TQo_GRAPH(INOD,14)    = Cs
			TQo_GRAPH(INOD,15)    = Csx
			TQo_GRAPH(INOD,16:17) = 0.D0
			TQo_GRAPH(INOD,18)    = S1
			TQo_GRAPH(INOD,19)    = S1x

		ENDDO

	ENDDO LOOP_ELEMENTS2

END SUBROUTINE CREATE_GRAPH_ARRAYS

!----------------------------------------------------------------------
!     SUBROUTINE WRITE_FILES
!----------------------------------------------------------------------

SUBROUTINE WRITE_FILES(WTS,PRINT_PROFILES,T,Dt)

	USE ELEMENTS_DATA,     only: NBF_Q, NNTOLa, NNTOLb, NEQ_Q_a, NEQ_Q_b, &
										NODTOLa_QEL, NODTOLb_QEL, NODTOL_QEL
	USE COMMON_ARRAYS,     only: TQo_GRAPH, TQa, TQao, TQb, TQbo, Xc, Xco, Xco1
	USE COMMON_ARRAYS,     only: uwarray, Z_loc, xa
	USE FLOW_PARAMETERS,   only: MASS_FLUIDi, MASS_SURFi
    USE FUNCTIONS

	IMPLICIT NONE

	CHARACTER(LEN=2), INTENT(IN) :: WTS
	LOGICAL, INTENT(IN) :: PRINT_PROFILES
	REAL(8), INTENT(IN) :: T, Dt
	
	CHARACTER(200) :: filename
	INTEGER :: I, J, recno, length
	REAL(8) :: MASS_FLUID, MASS_SURF, TIME
	REAL(8) :: MASS_ERROR_F, MASS_ERROR_S
	REAL(8) :: Theta, Hmax, Cac, Cc, Mc, Csc, Cao, Co, Mo, Cso, XXc, XXco

!     SELECT WHICH TIME INSTANT TO PRINT

	SELECT CASE(WTS)

	CASE('Tn')
     
		Time = T

		CALL CREATE_GRAPH_ARRAYS(TQa,TQb,Xc)

		Hmax   =  TQa(1,1)
		Theta  = -TQo_GRAPH(NODTOLa_QEL,4)
		XXc    =  Xc
		XXco   =  Xco
		Cao    =  TQa(1,4)
		Co     =  TQa(1,5)
		Mo     =  TQa(1,6)
		Cso    =  TQa(1,7)
		Cac    =  TQa(NODTOLa_QEL,4)
		Cc     =  TQa(NODTOLa_QEL,5)
		Mc     =  TQa(NODTOLa_QEL,6)
		Csc    =  TQa(NODTOLa_QEL,7)
		Thetaa =  Tha(Cac,Csc)

		CALL CALCULATE_MASS(TQa, TQb, Xc, MASS_FLUID, MASS_SURF)

	CASE('To')
     
		Time = T - Dt

		CALL CREATE_GRAPH_ARRAYS(TQao,TQbo,Xco)

		Hmax   =  TQao(1,1)
		Theta  = -TQo_GRAPH(NODTOLa_QEL,4)
		XXc    =  Xco
		XXco   =  Xco1
		Cao    =  TQao(1,4)
		Co     =  TQao(1,5)
		Mo     =  TQao(1,6)
		Cso    =  TQao(1,7)
		Cac    =  TQao(NODTOLa_QEL,4)
		Cc     =  TQao(NODTOLa_QEL,5)
		Mc     =  TQao(NODTOLa_QEL,6)
		Csc    =  TQao(NODTOLa_QEL,7)
		Thetaa =  Tha(Cac,Csc)

		CALL CALCULATE_MASS(TQao, TQbo, Xco, MASS_FLUID, MASS_SURF)

	END SELECT

!     WRITE OUTPUT FILES

	IF (PRINT_PROFILES) THEN

!     CPROFILE.DAT

		recno = ceiling((Time+Dt)/Dt)
		call get_Timestep_FileName(recno,'results/CPROFILE',filename)
		WRITE(filename,'(a24,a4)') filename,'.DAT'
		OPEN(43, file=filename, status='UNKNOWN', action='WRITE')

		WRITE(43,'(1A24)')'======================'
		WRITE(43,'(1A9,1ES14.6)')'Time = ', Time
		WRITE(43,'(1A24)')'======================'

		WRITE(43,'(10A14)') 'X','Z','Ca','C','M','Cs','Cax', &
							'Cx','Mx','Csx'

		DO I = 1, NODTOLa_QEL
			write(43,'(1X,10ES14.6)') &
			TQo_GRAPH(I,2), TQo_GRAPH(I,1), TQo_GRAPH(I,8), TQo_GRAPH(I,10), &
			TQo_GRAPH(I,12), TQo_GRAPH(I,14), TQo_GRAPH(I,9), TQo_GRAPH(I,11), &
			TQo_GRAPH(I,13), TQo_GRAPH(I,15)
!			| X, Z, Ca, C, M, Cs, Cax, Cx, Mx, Csx |
		ENDDO

		DO I = 2, NODTOLb_QEL
			J = NODTOLa_QEL + I - 1
			write(43,'(1X,2ES14.6,3A14,ES14.6,3A14,ES14.6)') &
			TQo_GRAPH(J,2), TQo_GRAPH(J,1), '-', '-', '-', TQo_GRAPH(J,14), &
			'-', '-', '-', TQo_GRAPH(J,15)
!			| X, Z, - , - , - , Cs, - , - , - , Csx |
		ENDDO
		write(43,*)

		close(43,status='keep')

		!	VPROFILE.DAT

		recno = ceiling((Time+Dt)/Dt)
		call get_Timestep_FileName(recno,'results/VPROFILE',filename)
		WRITE(filename,'(a24,a4)') filename,'.DAT'
		OPEN(46, file=filename, status='UNKNOWN', action='WRITE')

		WRITE(46,'(1A24)')'======================'
		WRITE(46,'(1A9,1ES14.6)')'Time = ', Time
		WRITE(46,'(1A24)')'======================'

		WRITE(46,'(11A14)')'X','Z','H','Hx','Hxx','Hxxx','U', &
						'S2','S2x','S12_or_S1','S12x_or_S1x'

		DO I = 1, NODTOLa_QEL
			write(46,'(1X,11ES14.6)') &
			TQo_GRAPH(I,2), TQo_GRAPH(I,1), TQo_GRAPH(I,3), TQo_GRAPH(I,4), &
			TQo_GRAPH(I,5), TQo_GRAPH(I,6), TQo_GRAPH(I,7), TQo_GRAPH(I,16), &
			TQo_GRAPH(I,17), TQo_GRAPH(I,18), TQo_GRAPH(I,19)
!			| X, Z, H, Hx, Hxx, Hxxx, U, S2, S2x, S12, S12x |
		ENDDO

		DO I = 2, NODTOLb_QEL
			J = NODTOLa_QEL + I - 1
			write(46,'(1X,2ES14.6,7A14,2ES14.6)') &
			TQo_GRAPH(J,2), TQo_GRAPH(J,1), '-', '-', '-', '-', '-', '-', &
			'-', TQo_GRAPH(J,18), TQo_GRAPH(J,19)
!			| X, Z, - , - , - , - , - , - , - , S1, S1x |
		ENDDO
		write(46,*)

		close(46,status='keep')

        !x position of every node
        call get_Timestep_FileName(recno,'results/xgrid',filename)
	    write(filename,'(a21,a4)') filename,'.DAT'
	    inquire(iolength=length) TQo_GRAPH(1:NODTOLa_QEL,2)
        !print*, trim(filename),shape(TQo_GRAPH(1:NODTOLa_QEL,2)),size(TQo_GRAPH(1:NODTOLa_QEL,2)),length
	    open(48, file=filename,action='write',form='unformatted',access='direct',recl=length)
        write(48,rec=1) TQo_GRAPH(1:NODTOLa_QEL,2)
        close(48)

        !2D z location of velocity values at every node
        call get_Timestep_FileName(recno,'results/zgrid',filename)
	    write(filename,'(a21,a4)') filename,'.DAT'
	    inquire(iolength=length) Z_loc
        !print*, trim(filename),shape(z_loc),size(z_loc),length
	    open(48, file=filename,action='write',form='unformatted',access='direct',recl=length)
        write(48,rec=1) Z_loc
        close(48)

        ! 2D domain snapshots of z values of velocity data at every node
        call get_Timestep_FileName(recno,'results/uwgrid',filename)
	    write(filename,'(a21,a4)') filename,'.DAT'
	    inquire(iolength=length) uwarray
        !print*, trim(filename),shape(uwarray),size(uwarray),length
	    open(48, file=filename,action='write',form='unformatted',access='direct',recl=length)
        write(48,rec=1) uwarray
        close(48)

	ENDIF

	!     TIMEDATA1.DAT

	IF (TIME < 1.D-10) &
	WRITE(44,'(14A14)')'Time','Log(Time)','Xc','Log(Xc)','Log(Xc-1)','dXcdt',&
					'H(0,t)','Theta','Thetaa'
					

    if (TIME .gt. 1e-7 .and. XXc .gt. 1e-7 .and. abs(XXc-1.D0) .gt. 1e-7) then
    	WRITE(44,'(1X,14ES14.6)') TIME, log10(TIME), XXc, log10(XXc), &
					log10(XXc-1.D0), (XXc-XXco)/Dt, Hmax, Theta, Thetaa
    else
    	WRITE(44,'(1X,14ES14.6)') TIME, 0.d0, XXc, 0.d0, &
					0.d0, (XXc-XXco)/Dt, Hmax, Theta, Thetaa        
    endif
							
!		|  Time, Log(Time), Xc, Log(Xc), Log(Xc-1), dXcdt, H(0,t), Theta, &
!					Thetaa |

!     TIMEDATA2.DAT

	IF (TIME < 1.D-10) &
	WRITE(47,'(14A14)')'Time','Log(Time)','Ca(0,t)','c(0,t)','m(0,t)',&
					'Cs(0,t)','Ca(xc,t)','c(xc,t)','m(xc,t)','Cs(xc,t)'

    if (TIME .gt. 1e-7 ) then
	    WRITE(47,'(1X,14ES14.6)') TIME, log10(TIME), Cao, Co, Mo, Cso, Cac, Cc, & 
							Mc, Csc
    else
	    WRITE(47,'(1X,14ES14.6)') TIME, 0.d0, Cao, Co, Mo, Cso, Cac, Cc, & 
							Mc, Csc
    endif
!		|  Time, Log(Time), Ca(0,t), c(0,t), m(0,t), Cs(0,t), Ca(xc,t),  
!					c(xc,t), m(xc,t), Cs(xc,t)  |

	!     ERROR.DAT

    if (MASS_FLUIDi .gt. 1e-7) then
        MASS_ERROR_F = 1.d2*(MASS_FLUID - MASS_FLUIDi)/MASS_FLUIDi
    else
        MASS_ERROR_F = 0.d0
    endif
	if (MASS_SURFi .gt. 1e-7) then
        MASS_ERROR_S = 1.d2*(MASS_SURF - MASS_SURFi) /MASS_SURFi
    else
        MASS_ERROR_S = 0.d0
    endif

	IF (TIME < 1.D-10) &
	WRITE(45,'(10A14)')'TIME','total_mass','error_mass', &
					'total_surfact','error_surfact'

	WRITE(45,'(1X,5ES14.6)') TIME, MASS_FLUID, MASS_ERROR_F, MASS_SURF, &
					MASS_ERROR_S
!		| Time, MASS_FLUID, MASS_ERROR, MASS_SURF, MASS_SURF_ERROR |


END SUBROUTINE WRITE_FILES

!----------------------------------------------------------------------
!					SUBROUTINE JACOBIAN
!----------------------------------------------------------------------
!Inputs:
!        ROI    -- Domain, inside or outside droplet
!        NELEM  -- Current element index
!        ININT  -- Number of Gauss points per element from NINT_3P (3)
!        KK     -- Current Gauss point index 
!        BFN_Q  -- Array of Gauss points (3,3)
!        DFDX_Q -- Array of derivative Gauss points (3,3)
!        JDIM   -- Number of Gauss points per element from NBF_Q (3)
!        NEQ    -- Number of equations (7, i.e. Z,H,Hx,Hxx,Ca,M,Cs)
!        NM     -- Array of node indices per element NM_Q in current ROI (NM_Q_a or NM_Q_b)
!        INNTOL -- Total number of elements inside and outside droplet
!        TQ     -- The MAIN array containing current values at all nodes for all equations 
!        TQo    -- MAIN array TQ at the previous timestep
!        X      -- Physical location of nodes (Fixed by MESH)
!        INDTL  -- Number of nodes in current ROI (NODTOLa_QEL or NODTOLb_QEL)
!        TERM   -- Working residual vector for element NELEM
!        TEMP_X -- Cumulative evolution of TERM??
!        TEMP_J -- Returned array for element??
!        WET_Q  -- This seems to employ the GAUSS_LINE and is used to evolve equations??


SUBROUTINE JACOBIAN ( ROI, NELEM, ININT, KK, BFN_Q, DFDX_Q, JDIM, NEQ, &
					  NM, INNTOL, TQ, TQo, X, INDTL, TERM, TEMP_X, TEMP_J, &
					  WET_Q )

	USE ELEMENTS_DATA, only: NBF_Q
	USE COMMON_ARRAYS, only: Xc

	IMPLICIT NONE

	INTEGER, PARAMETER :: IDIM_Q = NBF_Q

	INTEGER :: JW, JEQ, JJ
	REAL(8), save :: epert

	INTEGER, INTENT(IN) :: NELEM, KK, NEQ, JDIM, ININT, INNTOL, INDTL, ROI
	REAL(8), INTENT(IN) :: WET_Q

	INTEGER, INTENT(IN), DIMENSION( INNTOL, JDIM+1 ) :: NM

	REAL(8), INTENT(IN), DIMENSION( IDIM_Q, ININT )  :: BFN_Q
	REAL(8), INTENT(IN), DIMENSION( IDIM_Q )		 :: DFDX_Q

	REAL(8), INTENT(IN), DIMENSION( IDIM_Q, NEQ ) :: TERM

	REAL(8), INTENT(INOUT), DIMENSION( IDIM_Q, NEQ ) :: TEMP_X
	REAL(8), INTENT(INOUT), DIMENSION( IDIM_Q, JDIM, NEQ, NEQ ) :: TEMP_J
	REAL(8), INTENT(IN), DIMENSION( INDTL, NEQ )  :: TQ, TQo

	REAL(8), DIMENSION( IDIM_Q, NEQ ) :: TERM_X_e
	REAL(8), DIMENSION( IDIM_Q, JDIM, NEQ, NEQ ) :: TERM_J_e

	REAL(8) :: Xcj
	REAL(8), DIMENSION( INDTL )         :: X
	REAL(8), DIMENSION( INDTL, NEQ )    :: TQj

!----------------------------------------------------------------------

!	Initialize arrays
	TQj = TQ
	Xcj = Xc

!----------------------------------------------------------------------

    !Use to retain precision when xc only changes by small amounts?
	CALL XPERT(Xcj, epert)
    !Call routine with all of the physical equations 
	CALL EQUATIONS ( ROI, NELEM, ININT, KK, BFN_Q, DFDX_Q, NM, INNTOL, &
			         TQj, TQo, X, INDTL, NEQ, Xcj, TERM_X_e(:,:) )

	Xcj = Xc
	TEMP_X(:,:) = TEMP_X(:,:) + &
					((TERM_X_e(:,:) - TERM(:,:))/epert)*WET_Q

	LOOP_JACOBIAN: DO  JW = 1, JDIM

        !Loop over all equations
		DO JEQ = 1, NEQ
			JJ = NM(NELEM,JW+1)

			CALL XPERT(TQj(JJ,JEQ), epert)
			CALL EQUATIONS ( ROI, NELEM, ININT, KK, BFN_Q, DFDX_Q, NM, INNTOL, &
					         TQj, TQo, X, INDTL, NEQ, Xcj, TERM_J_e(:,JW,JEQ,:) )

			TQj(JJ,JEQ) = TQ(JJ,JEQ)

			TEMP_J(:,JW,JEQ,:) = TEMP_J(:,JW,JEQ,:) + &
			((TERM_J_e(:,JW,JEQ,:) - TERM(:,:))/epert)*WET_Q

		ENDDO

	ENDDO LOOP_JACOBIAN

!----------------------------------------------------------------------

END SUBROUTINE JACOBIAN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!							SUBROUTINE XPERT
!-----------------------------------------------------------------------
! Ensures X is not less than the sqrt of the machine
! precision and return the difference as variable h
! Same as max(X,eps)?

SUBROUTINE XPERT(X, H)
	USE MACHINE_EPSILON
	IMPLICIT NONE

	REAL(8) :: TEMP, EPS
	REAL(8), INTENT(INOUT) ::  X
	REAL(8), INTENT(OUT) ::  H

!-----------------------------------------------------------------------

	TEMP = X
	eps = DSQRT(EPSILON)
	h = eps*DABS(TEMP)
	IF ( h < eps ) h=eps
	X = TEMP + h
	h = X - TEMP

!-----------------------------------------------------------------------

END SUBROUTINE XPERT

!-----------------------------------------------------------------------

SUBROUTINE MESH

	USE ELEMENTS_DATA
	USE COMMON_ARRAYS, only: Xa, Xb, NM_Q_a, NM_Q_b
	USE BASIC_PARAMETERS, only: L1, L2, XPACK, MPa, MPb
	
	IMPLICIT NONE
	INTEGER :: I, J, IEL, I1, I2, I3
	REAL(8) :: XAXIS, STEP, NUMER, DENOM
	REAL(8), DIMENSION(NNXa) :: TMP_MESHXa
	REAL(8), DIMENSION(NNXb) :: TMP_MESHXb

	!------------------------------------------------------------------------
	!     CONSTRUCT UNIFORM MESH
	!------------------------------------------------------------------------

	!     REGION 1
	XAXIS = L1
	STEP  = XAXIS/DBLE(NXELa)
	Xa(1) = 0.D0
	DO I = 2, NNXa
		IF(MOD(I,2).EQ.0)THEN
			Xa(I) = DBLE((I-1)/2)*STEP+STEP/2.D0 + Xa(1)
		ELSE
			Xa(I) = DBLE(I/2)*STEP + Xa(1)
		ENDIF
	ENDDO
	
	!     REGION 2
	XAXIS = L2
	STEP  = XAXIS/DBLE(NXELb)
	Xb(1) = L1
	DO I = 1, NNXb
		IF(MOD(I,2).EQ.0)THEN
			Xb(I) = DBLE((I-1)/2)*STEP+STEP/2.D0 + Xb(1)
		ELSE
			Xb(I) = DBLE(I/2)*STEP + Xb(1)
		ENDIF
	ENDDO

	!------------------------------------------------------------------------
	!		MESH PACKING
	!------------------------------------------------------------------------
	IF ( XPACK == 'Y' ) THEN

		!	REGION 1
		TMP_MESHXa = 0.D0
		DO IEL = 1, NXELa
		DO I = 1, NBF_Q
			J = NM_Q_a(IEL,I+1)
			NUMER = (MPa+1.D0) - &
					(MPa-1.D0)*( (MPa+1.D0)/(MPa-1.D0) )**(1.D0-(L1-Xa(J))/L1)
			DENOM = ( (MPa+1.D0)/(MPa-1.D0) )**(1.D0-(L1-Xa(J))/L1) + 1.D0
			TMP_MESHXa(J) = L1*NUMER/DENOM
		ENDDO
		ENDDO

		Xa = L1 - TMP_MESHXa
		DO IEL = 1, NXELa
			I1 = NM_Q_a(IEL,2) ; I2 = NM_Q_a(IEL,3) ; I3 = NM_Q_a(IEL,4)
			Xa(I2) = (Xa(I1)+Xa(I3))/2.D0
		ENDDO

		!	REGION 2
		TMP_MESHXb = 0.D0
		DO IEL = 1, NXELb
		DO I = 1, NBF_Q
			J = NM_Q_b(IEL,I+1)
			NUMER = (MPb+1.D0) - &
					(MPb-1.D0)*( (MPb+1.D0)/(MPb-1.D0) )**(1.D0-(Xb(J)-L1)/L2)
			DENOM = ( (MPb+1.D0)/(MPb-1.D0) )**(1.D0-(Xb(J)-L1)/L2) + 1.D0
			TMP_MESHXb(J) = L2*NUMER/DENOM
		ENDDO
		ENDDO

		Xb = L1 + TMP_MESHXb
		DO IEL = 1, NXELb
			I1 = NM_Q_b(IEL,2) ; I2 = NM_Q_b(IEL,3) ; I3 = NM_Q_b(IEL,4)
			Xb(I2) = (Xb(I1)+Xb(I3))/2.D0
		ENDDO
	ENDIF

END SUBROUTINE MESH

!-----------------------------------------------------------------------
!							SUBROUTINE CHECK_CONVERGENCE
!-----------------------------------------------------------------------

SUBROUTINE CHECK_CONVERGENCE&
	(INCREMENT, ITER, FLAG_NR, CCL, RSUM_NEW, RSUM_OLD, B, IDIM_B)

	USE COMMON_ARRAYS,    only: TQa, TQap, TQb, TQbp, Xc, Xcp
	USE BASIC_PARAMETERS, only: NITER, ERROR, OUTFE_STEP

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: ITER, INCREMENT, IDIM_B

	CHARACTER(len=3), INTENT(INOUT) :: FLAG_NR
	CHARACTER(len=1), INTENT(OUT)   :: CCL

	REAL(8), INTENT(OUT)   :: RSUM_NEW
	REAL(8), INTENT(INOUT) :: RSUM_OLD

	REAL(8), INTENT(IN), DIMENSION(IDIM_B) :: B

	REAL(8) :: ERROR_NEW

	!-----------------------------------------------------------------------

	CCL = 'N'

	!-----------------------------------------------------------------------

	RSUM_NEW   = DSQRT( DOT_PRODUCT ( B , B ) )
	ERROR_NEW  = MAXVAL( DABS(B) )

	!-----------------------------------------------------------------------

	KIND_OF_NEWTON_RAPHSON: SELECT CASE(FLAG_NR)

	!-----------------------------------------------------------------------

	CASE('NRP')

		IF( RSUM_NEW .LT. 1.0D0 )THEN

			IF( RSUM_NEW .GT. ERROR ) THEN
				FLAG_NR = 'MNR'
			ELSE
				FLAG_NR = 'NRP'
			ENDIF

			IF( ITER .GE. NITER-12 )THEN
				FLAG_NR = 'NRP'
			ENDIF

		ELSE

			FLAG_NR = 'NRP'

		ENDIF

		IF(MOD(INCREMENT,OUTFE_STEP) == 0)THEN
			WRITE(42, 50)ITER, RSUM_NEW
			WRITE(42, 51)FLAG_NR, ERROR_NEW
		ENDIF

!-----------------------------------------------------------------------

	CASE('MNR')

		IF( RSUM_NEW .GE. 1.0D0 )THEN

			FLAG_NR = 'NRP'
			TQa = TQap
			TQb = TQbp
			Xc  = Xcp
			IF(MOD(INCREMENT,OUTFE_STEP) == 0)THEN
				WRITE(42, 52)ITER, RSUM_NEW
				WRITE(42, 51)FLAG_NR, ERROR_NEW
			ENDIF
			CCL = 'Y'
			RETURN
		ELSE

			IF(RSUM_NEW/RSUM_OLD .GT. 1.1D0)THEN
				FLAG_NR = 'NRP'
			ELSE
				FLAG_NR = 'MNR'
			ENDIF

		ENDIF

		IF( ITER .GE. NITER-12 )THEN
			IF( RSUM_NEW .GE. 0.5D-3 )THEN
				FLAG_NR = 'NRP'
				TQa = TQap
				TQb = TQbp
				Xc  = Xcp
				IF(MOD(INCREMENT,OUTFE_STEP) == 0)THEN
					WRITE(42, 52)ITER, RSUM_NEW
					WRITE(42, 51)FLAG_NR, ERROR_NEW
				ENDIF
				CCL = 'Y'
				RETURN
			ELSE
				FLAG_NR = 'NRP'
			ENDIF
		ENDIF

		IF(MOD(INCREMENT,OUTFE_STEP) == 0)THEN
			WRITE(42, 52)ITER, RSUM_NEW
			WRITE(42, 51)FLAG_NR, ERROR_NEW
		ENDIF

	!-----------------------------------------------------------------------
	CASE DEFAULT

		WRITE(42,*)' INCORRECT CHOICE OF NEWTON RAPHSON FLAG '

	!-----------------------------------------------------------------------
	END SELECT KIND_OF_NEWTON_RAPHSON

	!-----------------------------------------------------------------------

	IF( RSUM_NEW .GT. 5.0D+9 )THEN

		WRITE(42,*)'PROGRAM UNABLE TO CONVERGE'
		WRITE(42,*)'TOO LARGE NORMA !!!'
		WRITE(42,*)'PROGRAM STOP !!!'
		STOP
	ENDIF

	!-----------------------------------------------------------------------

	RSUM_OLD = RSUM_NEW

	!-----------------------------------------------------------------------
	!			FORMAT STATEMENTS
	!-----------------------------------------------------------------------

50    FORMAT(I3,'. ITERATION. FULL NEWTON RAPSHON, NORM= ', E15.6)
51    FORMAT(5x,' FLAG= ', A4, '. MAX RESIDUAL= ', E15.6)
52    FORMAT(I3,'. ITERATION. MODIFIED NEWTON RAPSHON, NORM= ', E15.6)

	!-----------------------------------------------------------------------

END SUBROUTINE CHECK_CONVERGENCE
