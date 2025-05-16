MODULE SCATTERING_NN_VARIATIONAL
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: PHASE_SHIFT_RESULT
    DOUBLE PRECISION :: delta1_BB
    DOUBLE PRECISION :: delta2_BB
    DOUBLE PRECISION :: epsilon_BB
    DOUBLE PRECISION :: delta1_S
    DOUBLE PRECISION :: delta2_S
    DOUBLE PRECISION :: epsilon_S
  END TYPE PHASE_SHIFT_RESULT

  INTEGER, PARAMETER :: NNE = 80
  INTEGER, PARAMETER :: NCH_MAX = 2
  INTEGER, PARAMETER :: NNN = NNE * 2
  INTEGER, PARAMETER :: NNR = 100000
  INTEGER :: NCH
  INTEGER :: NDIM ! GIVEN BY CORE_CORE_MATRIX_ELEMENTS
  INTEGER :: NEQ  ! GIVEN BY ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS

  DOUBLE PRECISION, PARAMETER :: HC = 197.327053D0
  DOUBLE PRECISION, PARAMETER :: MP = 938.272029D0
  DOUBLE PRECISION, PARAMETER :: MN = 939.565630D0
  DOUBLE PRECISION, PARAMETER :: MR = MP * MN / (MP + MN)
  DOUBLE PRECISION, PARAMETER :: PI = 4*DATAN(1.0D0)
  COMPLEX*16, PARAMETER ::       IM = (0.0D0, 1.0D0)
  DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0, ZERO = 0.0D0

  DOUBLE PRECISION :: K, HTM, M
  INTEGER :: LC(NCH_MAX), T, T1Z, T2Z

  TYPE, PUBLIC :: VARIATIONAL_PARAMETERS
    INTEGER :: J, L, S, TZ, IPOT, ILB, LEMP
    DOUBLE PRECISION :: E
    DOUBLE PRECISION :: HR1, H, RANGE, GAMMA, EPS, AF
    INTEGER :: NNL
  END TYPE VARIATIONAL_PARAMETERS

  TYPE(VARIATIONAL_PARAMETERS), PRIVATE :: VARIATIONAL_PARAMS

  PUBLIC :: SET_VARIATIONAL_PARAMETERS
  PUBLIC :: NN_SCATTERING_VARIATIONAL
  PRIVATE:: PRINT_DIVIDER

CONTAINS
  SUBROUTINE SET_VARIATIONAL_PARAMETERS(E, J, L, S, TZ, IPOT, ILB, LEMP)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: E
    INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP

    VARIATIONAL_PARAMS%J = J
    VARIATIONAL_PARAMS%L = L
    VARIATIONAL_PARAMS%S = S
    VARIATIONAL_PARAMS%TZ = TZ
    VARIATIONAL_PARAMS%IPOT = IPOT
    VARIATIONAL_PARAMS%ILB = ILB
    VARIATIONAL_PARAMS%LEMP = LEMP
    VARIATIONAL_PARAMS%E = E
    VARIATIONAL_PARAMS%HR1 = 0.01D0
    VARIATIONAL_PARAMS%H   = 0.02D0
    VARIATIONAL_PARAMS%AF  = 1.02D0
    VARIATIONAL_PARAMS%RANGE = 40.0D0
    VARIATIONAL_PARAMS%GAMMA = 4.D0
    VARIATIONAL_PARAMS%EPS = 0.25D0
    VARIATIONAL_PARAMS%NNL = 32
    
    SELECT CASE (TZ)
      CASE (1)
        M = MP
        T1Z = 1
        T2Z = 1
      CASE (0)
        M = MP * MN / (MP + MN)
        T1Z = 1
        T2Z =-1
      CASE (-1)
        M = MN
        T1Z =-1
        T2Z =-1
      CASE DEFAULT
        PRINT *, "Invalid TZ value"
        STOP
    END SELECT

    HTM = HC**2 / (2 * M)

  ! Ensure T is set such that T + L + S is odd
    IF (MOD(L + S, 2) == 0) THEN
      T = 1
    ELSE
      T = 0
    END IF

    K = DSQRT(2*E*MR) / HC

    ! PRINT 5
    ! PRINT 15, "L", "S", "T", "TZ", "J"
    ! PRINT 5
    LC(1) = L
    NCH = 1
    ! PRINT 20, LC(1), S, T, TZ, J
    IF (J-L.EQ.1) THEN 
      LC(2) = L + 2
      ! PRINT 20, LC(2), S, T, TZ, J
      NCH = 2
    ENDIF
    ! PRINT 5


  5 FORMAT(30("-"))
 15 FORMAT(" ", A5, A5, A5, A5, A5, A5)
 20 FORMAT(" ", I5, I5, I5, I5, I5, I5)
  END SUBROUTINE SET_VARIATIONAL_PARAMETERS
  

  ! This subroutine calculates the NN scattering using a variational method
  ! INPUT PARAMETERS:
  ! E     : Energy in MeV
  ! J     : Total angular momentum
  ! L     : Orbital angular momentum
  ! S     : Spin
  ! TZ    : Isospin projection
  ! IPOT  : Potential type
  ! ILB   : Interaction type
  ! LEMP  : Set EM potential, if 0 then only pure Coulomb
FUNCTION NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP, PRINT_COEFFICIENTS) RESULT(PHASE_SHIFT)
  IMPLICIT NONE
! INPUT PARAMETERS
  DOUBLE PRECISION, INTENT(IN) :: E
  INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP
  LOGICAL, INTENT(IN), OPTIONAL :: PRINT_COEFFICIENTS
  TYPE(PHASE_SHIFT_RESULT) :: PHASE_SHIFT

! VARIABLES AND PARAMETERS FOR DGESV
  DOUBLE PRECISION, DIMENSION(NNN, NNN) :: C, CC, CCC
  DOUBLE PRECISION, DIMENSION(NNN, NCH_MAX) :: CAR, CAI
  DOUBLE PRECISION, DIMENSION(NNN)      :: CARR, CAII
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NNN) :: XRCOEFF, XICOEFF
  DOUBLE PRECISION, DIMENSION(NNN) :: IPIV
  INTEGER :: INFO, IAK, IAB

! MATRICES FOR THE VARIATIONAL METHOD
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: ARI, AIR, ARR, AII
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: BD1, BD2, BD3, BD4
  DOUBLE PRECISION, DIMENSION(NCH_MAX,NCH_MAX) :: AM, AN, AMM, ANN

! COEFFICIENT R2
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: RD

! PHASE-SHIFTS AND MIXING ANGLES
  DOUBLE PRECISION :: AMIXR, AMIXG, DELTA1, DELTA2, DELTA1G, DELTA2G
  DOUBLE PRECISION :: AMIXGS, DELTA1S, DELTA2S

! S-MATRIX
  DOUBLE COMPLEX :: SMAT(NCH_MAX, NCH_MAX)
  
! INITIALIZE THE VARIATIONAL PARAMETERS
  CALL SET_VARIATIONAL_PARAMETERS(E, J, L, S, TZ, IPOT, ILB, LEMP)

  ! CALL PRINT_INFO()

! Evaluating the matrix elements
  CALL PRINT_DIVIDER
  CALL ASYMPTOTIC_CORE_MATRIX_ELEMENTS      (CAR, CAI)
  CALL PRINT_DIVIDER
  CALL CORE_CORE_MATRIX_ELEMENTS            (C)
  CALL PRINT_DIVIDER

  DO IAK=1, NEQ
! Preparing the matrix elements for the diagonalization
    CARR =-CAR(:,IAK)
    CAII =-CAI(:,IAK)
    CC  = C
    CCC = C

  ! Solving the eigenvalue problem using LAPACK dgesv
  ! Evaluating for the "c_{n, alpha}" coefficients
    CALL DGESV(NDIM, 1, CC , NNN, IPIV, CARR, NNN, INFO)
    CALL HANDLE_INFO_ERROR()  ! Handle the error after the first DGESV call
    !  WRITE(*,*) "INFO: ", INFO
    CALL DGESV(NDIM, 1, CCC, NNN, IPIV, CAII, NNN, INFO)
    CALL HANDLE_INFO_ERROR()  ! Handle the error after the second DGESV call
    !  WRITE(*,*) "INFO: ", INFO

    XRCOEFF(IAK,:) = CARR
    XICOEFF(IAK,:) = CAII
  ENDDO

! Calculating R coefficients
  !  WRITE(*,*) neq, NDIM
  ! WRITE(900,'(F20.16)') XICOEFF
  ! WRITE(901,'(F20.16)') XRCOEFF
  ! WRITE(902,'(F20.16)') CAR
  ! WRITE(903,'(F20.16)') CAI
  ! WRITE(904,'(F50.16)') C


  CALL MATRIX_MULTIPLY(BD1,XICOEFF,CAI)
  CALL MATRIX_MULTIPLY(BD2,XRCOEFF,CAI)
  CALL MATRIX_MULTIPLY(BD3,XICOEFF,CAR)
  CALL MATRIX_MULTIPLY(BD4,XRCOEFF,CAR)

  CALL PRINT_DIVIDER
  CALL ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS(ARI, AIR, ARR, AII)
  CALL PRINT_DIVIDER

  DO IAB=1,NEQ       
  DO IAK=1,NEQ
    AM(IAB,IAK)=BD1(IAK,IAB)+BD1(IAB,IAK)+AII(IAB,IAK)+AII(IAK,IAB)
    AN(IAB,IAK)=BD2(IAK,IAB)+BD3(IAB,IAK)+2.D0*AIR(IAB,IAK)   
  ENDDO
  ENDDO

  

  AMM = AM
  ANN =-AN
  CALL PRINT_DIVIDER
  !  WRITE(*,*) "AMM = ", AMM
  !  WRITE(*,*) "ANN = ", ANN
  CALL PRINT_DIVIDER


! Evaluating the "R_{alpha, beta}" matrix elements
  CALL DGESV(NEQ, NEQ, AMM, NCH_MAX, IPIV, ANN, NCH_MAX, INFO)
  CALL HANDLE_INFO_ERROR()  ! Handle the error after the third DGESV call
  !  WRITE(*,*) "INFO: ", INFO

  DO IAB=1,NEQ
  DO IAK=1,NEQ
    !  WRITE(*,*)"COEFF R",ANN(IAB,IAK)
  ENDDO 
  ENDDO 
  
! Evaluating the "R_{alpha, beta}" matrix elements to the second order
  CALL R_SECOND_ORDER()
  !  WRITE(*,*)
  DO IAB=1,NEQ
  DO IAK=1,NEQ
    !  WRITE(*,*)"COEFF R2 NORMALIZZATO", -RD(IAB,IAK)
  ENDDO 
  ENDDO 

! Writing the coefficients to a file in order torecreate the wave function
  IF (PRESENT(PRINT_COEFFICIENTS)) THEN
    IF (PRINT_COEFFICIENTS) THEN
      CALL WRITE_COEFFICIENTS_TO_RECREATE_THE_WAVE_FUNCTION()
    ENDIF
  ENDIF

! Calculating the phase shifts and mixing angles in the Blatt-Biedenharn convention
  CALL CALCULATE_PHASE_SHIFTS_BLATT
  !  WRITE(*,*) 
  !  WRITE(*,*)"BLATT-BIEDENHARN"
  !  WRITE(*,*)"MIXING ANGLE=",AMIXG
  !  WRITE(*,*)"SFASAMENTO1=",DELTA1G
  !  WRITE(*,*)"SFASAMENTO2=",DELTA2G

  PHASE_SHIFT%delta1_BB = DELTA1G
  PHASE_SHIFT%delta2_BB = DELTA2G
  PHASE_SHIFT%epsilon_BB = AMIXG


! Calculating the S-matrix
  CALL CALCULATE_S_MATRIX
  !  WRITE(*,*)
  !  WRITE(*,*) "S(1,1)=" , SMAT(1,1)
  !  WRITE(*,*) "S(1,2)=" , SMAT(1,2)
  !  WRITE(*,*) "S(2,2)=" , SMAT(2,2)


! Calculating the phase shifts and mixing angles in the Stapp convention
  CALL CALCULATE_PHASE_SHIFTS_STAPP
  !  WRITE(*,*)
  !  WRITE(*,*)"STAPP"
  !  WRITE(*,*)"MIXING ANGLE=",AMIXGS
  !  WRITE(*,*)"SFASAMENTO1=",DELTA1S
  !  WRITE(*,*)"SFASAMENTO2=",DELTA2S

  PHASE_SHIFT%delta1_S = DELTA1S
  PHASE_SHIFT%delta2_S = DELTA2S
  PHASE_SHIFT%epsilon_S = AMIXGS

  !  WRITE(*,*) DELTA1S, DELTA2S, AMIXGS

  RETURN 

  CONTAINS
    SUBROUTINE HANDLE_INFO_ERROR()
      IMPLICIT NONE
      IF (INFO.NE.0) THEN
        PRINT *, "Error in LAPACK dgesv: INFO = ", INFO
        STOP
      ENDIF
    END SUBROUTINE HANDLE_INFO_ERROR

    SUBROUTINE PRINT_INFO()
      IMPLICIT NONE
      PRINT *, "E:    ",                  VARIATIONAL_PARAMS%E, " MeV"
      PRINT *, "HTM:  ",                  HTM, " MeV fm^2"
      PRINT *, "k:    ",                  K, " fm^-1"
      PRINT 10, "J:     ",                VARIATIONAL_PARAMS%J
      PRINT 10, "NCH:   ",                NCH
      PRINT 10, "L0:    ",                LC(1)
      IF (NCH.EQ.2) PRINT 10, "L1:    ",  LC(2)
      PRINT 10, "S:     ",                VARIATIONAL_PARAMS%S
      PRINT 10, "T:     ",                T
      PRINT 10, "TZ:    ",                VARIATIONAL_PARAMS%TZ

      10 FORMAT(" ",A, I2)
    END SUBROUTINE PRINT_INFO

    SUBROUTINE MATRIX_MULTIPLY(M1,M2,M3)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(NCH_MAX,NCH_MAX) :: M1
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NCH_MAX, NNN) :: M2
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NNN, NCH_MAX) :: M3
      INTEGER :: II, JJ, LL

      M1 = 0
      DO II = 1, NEQ
      DO JJ = 1, NEQ
        DO LL = 1, NDIM
          M1(II,JJ) = M1(II,JJ) + M2(II,LL)*M3(LL,JJ)
        ENDDO
      ENDDO
      ENDDO
    END SUBROUTINE MATRIX_MULTIPLY

    SUBROUTINE R_SECOND_ORDER()
      IMPLICIT NONE
      INTEGER :: I, IK, IB
      DOUBLE PRECISION :: SOMMA, SOMMA1, SOMMA2, SOMMA3, SOMMA4, SOMMA5, SOMMA6
      DOUBLE PRECISION, DIMENSION(NNN, NCH_MAX) :: XRCOEFV, XICOEFV
      DOUBLE PRECISION, DIMENSION(NCH_MAX, NNN) :: BD5, BD6
      DOUBLE PRECISION, DIMENSION(NCH_MAX, NNN) :: RCI
      DOUBLE PRECISION, DIMENSION(NNN, NCH_MAX) :: RCIV
      DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: CD0, CD1, CD2, CD3
      DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: CD4, CD5, CD6
      DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: CD7, CD8
      DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: RD0, RD1, RD2, RD3
      DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: ASS, ANNS

      DO IAK=1,NEQ
        DO IK=1,NDIM
          XRCOEFV(IK,IAK) = XRCOEFF(IAK,IK)
          XICOEFV(IK,IAK) = XICOEFF(IAK,IK)
        ENDDO
      ENDDO


      DO IAB=1,NEQ
        DO IK=1,NDIM
          SOMMA  = 0.D0
          SOMMA1 = 0.D0
          DO IB=1,NDIM
            SOMMA  = SOMMA  + XRCOEFF(IAB,IB)*C(IB,IK) 
            SOMMA1 = SOMMA1 + XICOEFF(IAB,IB)*C(IB,IK) 
          ENDDO
          BD5(IAB,IK) = SOMMA
          BD6(IAB,IK) = SOMMA1
        ENDDO
      ENDDO


      DO IAB=1,NEQ
        DO IB=1,NDIM
          SOMMA=0.D0
          DO IAK=1,NEQ
            SOMMA = SOMMA + ANN(IAB,IAK)*XICOEFF(IAK,IB)
          ENDDO
          RCI(IAB,IB)  = SOMMA
          RCIV(IB,IAB) = RCI(IAB,IB)
        ENDDO
      ENDDO


      DO IAB=1,NEQ
      DO IAK=1,NEQ
        SOMMA  = 0.D0
        SOMMA1 = 0.D0
        SOMMA2 = 0.D0
        DO IK=1,NDIM
          SOMMA  = SOMMA  + BD5(IAB,IK)*XRCOEFV(IK,IAK)
          SOMMA1 = SOMMA1 + BD5(IAB,IK)*RCIV(IK,IAK)
          SOMMA2 = SOMMA2 + BD6(IAB,IK)*XRCOEFV(IK,IAK)
        ENDDO
        CD0(IAB,IAK) = SOMMA
        CD1(IAB,IAK) = SOMMA1
        RD0(IAB,IAK) = SOMMA2
      ENDDO
      ENDDO

      DO I=1,NEQ
      DO IAB=1,NEQ
        SOMMA  = 0.D0
        SOMMA1 = 0.D0
        SOMMA2 = 0.D0
        SOMMA3 = 0.D0
        SOMMA4 = 0.D0
        SOMMA5 = 0.D0
        SOMMA6 = 0.D0
        DO IAK=1,NEQ
          SOMMA =SOMMA  + BD2(I,IAK) * ANN(IAK,IAB)
          SOMMA1=SOMMA1 + RD0(I,IAK) * ANN(IAK,IAB)
          SOMMA2=SOMMA2 + BD3(I,IAK) * ANN(IAK,IAB) 
          SOMMA3=SOMMA3 + ARI(I,IAK) * ANN(IAK,IAB) 
          SOMMA4=SOMMA4 + AIR(I,IAK) * ANN(IAK,IAB) 
          SOMMA5=SOMMA5 + BD1(I,IAK) * ANN(IAK,IAB)
          SOMMA6=SOMMA6 + AII(I,IAK) * ANN(IAK,IAB)
        ENDDO
        CD2(I,IAB) = SOMMA
        CD3(I,IAB) = SOMMA1
        CD4(I,IAB) = SOMMA2
        CD5(I,IAB) = SOMMA3
        CD6(I,IAB) = SOMMA4
        RD2(I,IAB) = SOMMA5
        RD3(I,IAB) = SOMMA6
      ENDDO
      ENDDO
      
      DO I=1,NEQ
      DO IAB=1,NEQ
        SOMMA  = 0.D0
        SOMMA1 = 0.D0
        DO IAK=1,NEQ
          SOMMA  = SOMMA  + RD2(I,IAK) * ANN(IAK,IAB)
          SOMMA1 = SOMMA1 + RD3(I,IAK) * ANN(IAK,IAB)
        ENDDO
        CD7(I,IAB) = SOMMA
        CD8(I,IAB) = SOMMA1
      ENDDO
      ENDDO


      !  WRITE(*,*)
      DO IAB=1,NEQ
      DO IAK=1,NEQ
        ASS(IAB,IAK)=CD0(IAB,IAK)+2*BD4(IAB,IAK)+CD4(IAB,IAK)    &
                    +ARR(IAB,IAK)+CD5(IAB,IAK)+CD2(IAB,IAK)    &
                    +CD7(IAB,IAK)+CD6(IAB,IAK)    &
                    +CD8(IAB,IAK)
  
        ANNS(IAB,IAK) = ANN(IAB,IAK) + ASS(IAB,IAK)
        !  WRITE(*,*)"COEFF R2",ANNS(IAB,IAK), ANN(IAB,IAK), ASS(IAB,IAK)
      ENDDO
      ENDDO

      !SYMMETRIZATION
      DO IAB=1,NEQ
      DO IAK=1,NEQ
        RD(IAB,IAK) =-0.5D0*( ANNS(IAB,IAK) + ANNS(IAK,IAB) )
      ENDDO 
      ENDDO 
    END SUBROUTINE R_SECOND_ORDER

    SUBROUTINE WRITE_COEFFICIENTS_TO_RECREATE_THE_WAVE_FUNCTION(FILE)
      IMPLICIT NONE
      INTEGER :: JP, NNL, I
      INTEGER :: LLA(NCH_MAX), LSA(NCH_MAX), LTA(NCH_MAX), LJA(NCH_MAX)
      DOUBLE PRECISION :: EPS, GAMMA
      CHARACTER(LEN=*), OPTIONAL :: FILE

      IF (PRESENT(FILE)) THEN
        OPEN(UNIT=19, FILE=TRIM(FILE), STATUS='UNKNOWN', ACTION="WRITE")
      ELSE
        OPEN(UNIT=19, FILE='fort.19', STATUS='UNKNOWN', ACTION="WRITE")
      ENDIF

      NNL = VARIATIONAL_PARAMS%NNL
      EPS = VARIATIONAL_PARAMS%EPS
      GAMMA = VARIATIONAL_PARAMS%GAMMA

      JP = 0
      IF ( MOD(L,2) .EQ. 0) JP = 1

      IF(J.EQ.0.AND.JP.EQ.1)THEN ! 1S0
         LLA(1)=0
         LSA(1)=0
         LTA(1)=1
         LJA(1)=0
      ENDIF
      IF(J.EQ.0.AND.JP.EQ.0)THEN ! 3P0
         LLA(1)=1
         LSA(1)=1
         LTA(1)=1
         LJA(1)=0
      ENDIF
      IF(J.EQ.1.AND.JP.EQ.0)THEN ! 3P1
         LLA(1)=1
         LSA(1)=1
         LTA(1)=1
         LJA(1)=1
      ENDIF
      IF(J.EQ.2.AND.JP.EQ.1)THEN ! 1D2
         LLA(1)=2
         LSA(1)=0
         LTA(1)=1
         LJA(1)=2
      ENDIF
      IF(J.EQ.2.AND.JP.EQ.0)THEN ! 3P2-3F2
         LLA(1)=1
         LSA(1)=1
         LTA(1)=1
         LJA(1)=2
         
         LLA(2)=3
         LSA(2)=1
         LTA(2)=1
         LJA(2)=2
      ENDIF

      WRITE(19,*)NEQ,GAMMA,NNL
      WRITE(19,*)EPS

      DO IAK=1,NEQ
        WRITE(19,*)LLA(IAK),LSA(IAK),LTA(IAK),LJA(IAK)

        DO I=1,NDIM
          WRITE(19,*) XICOEFF(IAK,I)
        ENDDO

        DO I=1,NDIM
          WRITE(19,*) XRCOEFF(IAK,I)
        ENDDO

      ENDDO

      DO IAK=1,NEQ
      DO IAB=1,NEQ
        WRITE(19,*)-RD(IAK,IAB)
      ENDDO
      ENDDO

    END SUBROUTINE WRITE_COEFFICIENTS_TO_RECREATE_THE_WAVE_FUNCTION

    SUBROUTINE CALCULATE_PHASE_SHIFTS_BLATT()
      IMPLICIT NONE
      AMIXR=0.5D0*ATAN(2.*RD(1,2)/(RD(1,1)-RD(2,2)))
      AMIXG=(AMIXR*180.D0)/PI     
  
      DELTA1=ATAN((COS(AMIXR)*COS(AMIXR)*RD(1,1)  &
                   +SIN(AMIXR)*SIN(AMIXR)*RD(2,2)  &
                   +2*COS(AMIXR)*SIN(AMIXR)*RD(1,2)))

      DELTA2=ATAN((SIN(AMIXR)*SIN(AMIXR)*RD(1,1)  &
                   +COS(AMIXR)*COS(AMIXR)*RD(2,2)  &
                   -2*COS(AMIXR)*SIN(AMIXR)*RD(1,2)))
      DELTA1G=(DELTA1*180.D0)/PI 
      DELTA2G=(DELTA2*180.D0)/PI

      !  WRITE(*,*)"BLATT-BIEDENHARN"
      !  WRITE(*,*)"MIXING ANGLE=",AMIXG
      !  WRITE(*,*)"SFASAMENTO1=",DELTA1G
      !  WRITE(*,*)"SFASAMENTO2=",DELTA2G

    END SUBROUTINE CALCULATE_PHASE_SHIFTS_BLATT

    SUBROUTINE CALCULATE_S_MATRIX()
      IMPLICIT NONE
      DOUBLE PRECISION :: COS1, SIN1
      DOUBLE COMPLEX :: SM1, SM2
      SM1  = CDEXP(2.D0*IM*DELTA1)
      SM2  = CDEXP(2.D0*IM*DELTA2)
      COS1 = COS(AMIXR)
      SIN1 = SIN(AMIXR)
      SMAT(1,1) = COS1*COS1*SM1 + SIN1*SIN1*SM2
      SMAT(2,2) = COS1*COS1*SM2 + SIN1*SIN1*SM1
      SMAT(1,2) = COS1*SIN1*(SM1 - SM2)
    END SUBROUTINE CALCULATE_S_MATRIX
    
    SUBROUTINE CALCULATE_PHASE_SHIFTS_STAPP()
      IMPLICIT NONE
      DOUBLE COMPLEX :: SM1, SM2
      DOUBLE PRECISION :: SI2E, CI2E, CX, SX
      SM1  = SMAT(1,1)*SMAT(2,2)-SMAT(1,2)*SMAT(1,2)
      SI2E =-SMAT(1,2)*SMAT(1,2)/SM1
      CI2E = 1.D0 - SI2E
      SI2E = DSQRT(SI2E)
      CI2E = DSQRT(CI2E)
      
    ! I MIXING ANGLES DELLE ONDE DISPARI (JP=0) VENGONO COL SEGNO SBAGLIATO!
      SM1 = CDSQRT(SMAT(1,1)/CI2E)
      SM2 = CDSQRT(SMAT(2,2)/CI2E)
      CX  = DREAL(SM1)
      SX  = DIMAG(SM1)
      DELTA1S = DACOS(CX)*180.D0/PI
      IF(SX.LT.0.D0) DELTA1S = -DELTA1S
      
      CX = DREAL(SM2)
      SX = DIMAG(SM2)
      DELTA2S = DACOS(CX)*180.D0/PI
      IF(SX.LT.0.D0) DELTA2S = -DELTA2S

      AMIXGS = (0.5D0*DASIN(SI2E))*180.D0/PI

    END SUBROUTINE CALCULATE_PHASE_SHIFTS_STAPP

  END FUNCTION NN_SCATTERING_VARIATIONAL



  SUBROUTINE CORE_CORE_MATRIX_ELEMENTS(AM)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NNN, NNN), INTENT(OUT) :: AM
    INTEGER, PARAMETER :: NNRCC=200
    DOUBLE PRECISION :: XPNT(NNRCC), PWEIGHT(NNRCC)
    DOUBLE PRECISION :: XX(NNRCC), WG(NNRCC), YY(NNRCC)
    DOUBLE PRECISION :: U0(0:NNE,NNRCC), U1(0:NNE,NNRCC), U2(0:NNE,NNRCC)
    DOUBLE PRECISION :: V0(NNE,NNRCC), V1(NNE,NNRCC), V2(NNE,NNRCC)
    DOUBLE PRECISION :: GAMMA, APF, XG, ANL, R
    INTEGER :: I, IX, NX, NMX
    INTEGER :: L, S, J, NNL
    DOUBLE PRECISION :: V(NNRCC, NCH_MAX, NCH), VPW(NCH_MAX, NCH_MAX)
    INTEGER :: ICONT(NCH_MAX,NNE), LIK, IAB, IAK, IL, IR, IB, IK
    DOUBLE PRECISION :: SUM, AKEM(NNN,NNN), APEM(NNN,NNN), FUN(NNRCC)

    DOUBLE PRECISION, SAVE :: HM(NNN,NNN), AXX(NNN,NNN)
    LOGICAL, SAVE :: FIRST_CALL = .TRUE.

    GAMMA = VARIATIONAL_PARAMS%GAMMA
    NNL = VARIATIONAL_PARAMS%NNL

    NX = 100
    IF (FIRST_CALL) THEN
      CALL GAULAG(NX, XPNT,PWEIGHT)

      DO I=1,NX                                                    
        XX(I)=XPNT(I)/GAMMA
        WG(I)=PWEIGHT(I)                                   
        YY(I)=XPNT(I)     ! grid for evaluating Laguerre
      ENDDO  

      !  WRITE(*,*)'PRIMO E ULTIMO PUNTO =',XX(1),XX(NX)
      
      NMX = VARIATIONAL_PARAMS%NNL-1                                   
      APF=2.D0
      CALL LAGUERRE_POLYNOMIAL(YY, NX, APF, NMX, U0, U1, U2)
      DO IX=1, NX
        DO I=0, NMX
          ANL = DSQRT(DGAMMA(I+1.D0)*GAMMA**3/DGAMMA(I+3.D0))

          V0(I+1,IX) = ANL * U0(I,IX)
          V1(I+1,IX) = ANL * GAMMA * (U1(I,IX) - 0.5D0*U0(I,IX))
          V2(I+1,IX) = ANL * GAMMA**2 * (U2(I,IX) - 0.5D0*U1(I,IX)) &
                        - 0.5D0*GAMMA*V1(I+1,IX)
        ENDDO
      ENDDO

      L = LC(1)
      S = VARIATIONAL_PARAMS%S
      J = VARIATIONAL_PARAMS%J

      DO I=1, NX
        R = XX(I)
        CALL AV18PW90(1, L, S, J, T, T1Z, T2Z, R, VPW, VARIATIONAL_PARAMS%LEMP)
        V(I,1,1) = VPW(1,1)
        V(I,1,2) = VPW(1,2)
        V(I,2,1) = VPW(2,1)
        V(I,2,2) = VPW(2,2)
      ENDDO

      CALL PREPARE_INDECES

      DO IAB=1,NEQ          
      DO IAK=1,NEQ
            
        LIK=LC(IAK)*(LC(IAK)+1)

        DO IL=1,NNL            
        DO IR=1,NNL            
                
          IB=ICONT(IAB,IL)       
          IK=ICONT(IAK,IR)       

    ! SI CALCOLA LA NORMA
          AXX(IB,IK)=0.D0
          IF(IB.EQ.IK) AXX(IB,IK) = 1.D0
    ! SI CALCOLA ENERGIA CINETICA
          AKEM(IB,IK)=0.D0              
          IF(IAB.EQ.IAK)THEN
            SUM=0.D0
            DO I=1,NX
              FUN(I)=V0(IL,I)*(V2(IR,I)+2.D0*V1(IR,I)/XX(I) &  
                              -LIK*V0(IR,I)/XX(I)**2 )
              SUM=SUM + XX(I)**2*FUN(I)*WG(I) 
            ENDDO                                                                 
            AKEM(IB,IK)=-HTM*SUM/GAMMA                                       
          ENDIF

          IF(IB.EQ.1.AND.IK.EQ.1)THEN
            !  WRITE(*,*)
            !  WRITE(*,*)'C-C MATRIX'
            !  WRITE(*,*)'KINETIC',AKEM(1,1)
          ENDIF

    ! SI CALCOLA ENERGIA POTENZIALE
          SUM=0.D0                                         
          DO I=1,NX
            FUN(I)=V0(IL,I)*V0(IR,I)*V(I,IAB,IAK)     
            SUM=SUM+XX(I)*XX(I)*FUN(I)*WG(I)
          ENDDO
          APEM(IB,IK)=1./GAMMA*SUM                                    

    ! SI CALCOLA HAMILTONIANA    
          ! AM(IB,IK)=1./HTM*(AKEM(IB,IK)+APEM(IB,IK)-AXX(IB,IK))                                  
          IF(IB.EQ.1.AND.IK.EQ.1)THEN
            !  WRITE(*,*)'POTENTIAL',APEM(1,1)
            !  WRITE(*,*)'C-C MATRIX',HTM*AM(IB,IK)
          ENDIF

        ENDDO ! IR
        ENDDO ! IL

      ENDDO ! IAK
      ENDDO ! IAB
      FIRST_CALL = .FALSE.
      HM = ( AKEM + APEM )
    ENDIF

    AM = (HM - AXX * VARIATIONAL_PARAMS%E )/ HTM
    !  WRITE(*,*)'C-C MATRIX',HTM*AM(1,1)


  CONTAINS 
    SUBROUTINE PREPARE_INDECES()
      IMPLICIT NONE
      INTEGER II, JJ

      II = 0
      DO I=1, NEQ
        DO JJ=1, VARIATIONAL_PARAMS%NNL
          II = II + 1
          ICONT(I, JJ) = II
        ENDDO
      ENDDO
      NDIM = II
    END SUBROUTINE PREPARE_INDECES
  END SUBROUTINE CORE_CORE_MATRIX_ELEMENTS

  SUBROUTINE ASYMPTOTIC_CORE_MATRIX_ELEMENTS(AM, AM1)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NNN, NCH_MAX), INTENT(OUT) :: AM(NNN, NCH_MAX), AM1(NNN, NCH_MAX)

    DOUBLE PRECISION :: H5, HR ! Step size in r
    INTEGER :: I, IX, NX, NMX ! Number evenly spaced points
    DOUBLE PRECISION :: XX(NNR), AJ(NNR), YYB(NNR), YYL(NNR), A(NNR)
    DOUBLE PRECISION :: U0(0:NNE, NNR), U1(0:NNE, NNR), U2(0:NNE, NNR) 
    DOUBLE PRECISION :: V0(NNE, NNR), V1(NNE, NNR), V2(NNE, NNR) 
    DOUBLE PRECISION :: FBES(NCH_MAX, NNR), GBES(NCH_MAX, NNR)
    DOUBLE PRECISION :: APF, GAMMA, ANL, XG, FEXP, RR
    INTEGER :: NEQC, L, S, J
    DOUBLE PRECISION :: VPW(2, 2), VV(NNR, NCH_MAX, NCH_MAX)
    INTEGER :: ICONT(NCH_MAX, NNE)
    INTEGER :: IAB, IAK, LIK, IL, IB
    DOUBLE PRECISION :: FUN(NNR), FUN1(NNR)
    DOUBLE PRECISION :: AXXM1(NNN, NCH_MAX), AXX1
    DOUBLE PRECISION :: AKEM(NNN, NCH_MAX), AKE1, AKEM1(NNN, NCH_MAX)
    DOUBLE PRECISION :: APE, APE1, APEM(NNN, NCH_MAX), APEM1(NNN, NCH_MAX)
    DOUBLE PRECISION, EXTERNAL :: B5

    GAMMA = VARIATIONAL_PARAMS%GAMMA

    HR = VARIATIONAL_PARAMS%HR1
    H5 = HR/22.5D0
    NX = INT(VARIATIONAL_PARAMS%RANGE/HR) + 10

    IF (VARIATIONAL_PARAMS%RANGE.LT.H5 .OR. VARIATIONAL_PARAMS%RANGE.GT.200.D0) THEN
      PRINT *, "Error: RANGE out of bounds"
      STOP
    ENDIF


    IF (ABS(NX).GT.NNR) THEN
      PRINT *, "Error: NX exceeds NNR"
      STOP
    ENDIF

  ! Initialize grid with r values
    !  WRITE(*,*)'NX =',NX
    DO I=1, NX
      XX(I)   = HR*I
      AJ(I)   = XX(I)**2
      YYB(I)  = K*XX(I)
      YYL(I)  = GAMMA*XX(I)
      A(I)    = 1.D0 - DEXP(-VARIATIONAL_PARAMS%EPS*XX(I))
    ENDDO
    !  WRITE(*,*)'PRIMO E ULTIMO PUNTO =',XX(1),XX(NX),NX

  ! Laguerre polynomial and their first two derivatives grid
    NMX = VARIATIONAL_PARAMS%NNL - 1
    APF = 2.D0

    CALL LAGUERRE_POLYNOMIAL(YYL, NX, APF, NMX, U0, U1, U2)
    DO IX = 1, NX
      XG = YYL(IX)
      FEXP = DEXP(-XG/2.D0)
      DO I = 0, NMX
        ANL = DSQRT(DGAMMA(I+1.D0)*GAMMA**3/DGAMMA(I+3.D0))*FEXP

        V0(I+1, IX) = ANL * U0(I,IX)
        V1(I+1, IX) = ANL * GAMMA * ( U1(I,IX) -0.5D0*U0(I,IX) )
        V2(I+1, IX) = ANL * GAMMA * ( GAMMA * ( U2(I,IX) -0.5D0*U1(I,IX) ) ) &
                    - 0.5D0*GAMMA*V1(I+1,IX) 
      ENDDO
    ENDDO
    
    NEQ = NCH
    NEQC= NCH
    
  ! Prepare the Bessel functions
    CALL SPHERICAL_BESSEL_FUNCTIONS()
    ! do m=1, NX
    !   do i=0, nmx
    !     Write(1000,*) yyl(M), v0(i+1,m), v1(i+1,m), v2(i+1,m), aj(m)
    !   enddo
    ! enddo


    ! DO I = 1, NCH
    !   DO M = 1, NX
    !     WRITE(20+I, *) YYB(M), FBES(I,M), GBES(I,M)
    !   ENDDO
    ! ENDDO
    
  ! Prepare the potential 
    L = LC(1)
    S = VARIATIONAL_PARAMS%S
    J = VARIATIONAL_PARAMS%J

    DO I = 1, NX
      RR = XX(I)
      CALL AV18PW90(1, LC(1), S, J, T, T1Z, T2Z, RR, VPW, VARIATIONAL_PARAMS%LEMP)
      VV(I, 1, 1) = VPW(1, 1)
      VV(I, 1, 2) = VPW(1, 2)
      VV(I, 2, 1) = VPW(2, 1)
      VV(I, 2, 2) = VPW(2, 2)
      ! WRITE(23, *) XX(I), VV(I, 1, 1), VV(I, 1, 2), VV(I, 2, 1), VV(I, 2, 2)  
    ENDDO

  ! Prepare the indeces for the matrix elements
    CALL PREPARE_INDECES()


    ! Evaluate the matrix elements
    DO IAB = 1, NEQC
      DO IAK = 1, NEQ
        LIK = LC(IAB)*(LC(IAB)+1)

        DO IL = 1, VARIATIONAL_PARAMS%NNL
          IB = ICONT(IAB, IL)

        ! Evaluate the normalization core-irregular (axx1)
          AXXM1(IB,IAK)=0.D0
          IF(IAB.EQ.IAK)THEN 
            FUN1(1)=0.D0                                         
            DO I=1,NX  
            FUN1(I+1)=AJ(I)*V0(IL,I)*GBES(IAK,I)
            ! WRITE(150,*) DSQRT(AJ(I)), V0(IL,I), GBES(IAK,I)
            ENDDO
            ! STOP
            !  WRITE(111,*) FUN1

            AXX1=VARIATIONAL_PARAMS%E * B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)
            AXXM1(IB,IAK)=AXX1
            ! write(111,*) iab, iak, il, IB, axx1
          ENDIF

        ! Evaluate the kinetic energy core-irregular (ake1)
          AKE1=0.D0
          AKEM1(IB,IAK)=0.D0     
          IF(IAB.EQ.IAK)THEN
             FUN1(1)=0.D0
             DO I=1,NX
              FUN1(I+1)=AJ(I)*GBES(IAK,I)*(V2(IL,I) + 2.D0*V1(IL,I)/XX(I)-LIK*V0(IL,I)/XX(I)**2)
             ENDDO
             AKE1= -HTM * B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)                                
             AKEM1(IB,IAK)=AKE1
            ! write(112,*) iab, iak, il, IB, ake1
          ENDIF
          IF(IB.EQ.1.AND.IAK.EQ.1)THEN
            !  WRITE(*,*)
            !  WRITE(*,*)'C-A MATRIX'
            !  WRITE(*,*)'IRREGULAR A'
            !  WRITE(*,*)'NORM ',AXXM1(1,1)
            !  WRITE(*,*)'KINETIC',AKEM1(1,1)
          ENDIF
        
        
        ! Evaluate the potential energy core-regular (ape), core-irregular (ape1)
          FUN(1)=0.D0
          FUN1(1)=0.D0                                 
          DO I=1,NX   
            FUN(I+1)=AJ(I)*V0(IL,I)*FBES(IAK,I)*VV(I,IAB,IAK)     
            FUN1(I+1)=AJ(I)*V0(IL,I)*GBES(IAK,I)*VV(I,IAB,IAK)                             
          ENDDO
          APE=B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)
          APEM(IB,IAK)=APE
          APE1=B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)
          APEM1(IB,IAK)=APE1
          ! write(113,*) iab, iak, il, IB, ape, ape1
          IF(IB.EQ.1.AND.IAK.EQ.1)THEN
            !  WRITE(*,*)
            !  WRITE(*,*)'C-A MATRIX'
            !  WRITE(*,*)'IRREGULAR A'
            !  WRITE(*,*)'POTENTIAL ',APEM1(1,1)
            !  WRITE(*,*)'REGULAR A'
            !  WRITE(*,*)'POTENTIAL ',APEM(1,1)
          ENDIF

        ! Evaluate the Hamiltonian: core-regular (am), core-irregular (am1)

          AM(IB,IAK) = 1.D0/HTM * APEM(IB,IAK)
          AM1(IB,IAK)= 1.D0/HTM * (AKEM1(IB,IAK)+APEM1(IB,IAK)-AXXM1(IB,IAK))
          
          IF(IB.EQ.1.AND.IAK.EQ.1)THEN
            !  WRITE(*,*)
            !  WRITE(*,*)'C-A MATRIX'     ,IB,IAK
            !  WRITE(*,*)"CORE-REGULAR=  ",AM(IB,IAK),IB,IAK
            !  WRITE(*,*)"CORE-IRREGULAR=",AM1(IB,IAK),IB,IAK
          ENDIF
        ENDDO
      ENDDO
    ENDDO 
    RETURN

  CONTAINS
    SUBROUTINE SPHERICAL_BESSEL_FUNCTIONS()
      use gsl_bessel
      IMPLICIT NONE
      
      DOUBLE PRECISION :: AG, GBSS, FBSS
      INTEGER :: LL
      
      DO I=1,NCH
        DO IX=1, NX                                 
          XG=YYB(IX)
          AG=A(IX)
          LL=LC(I)
          IF(K.LE.1.D-8)THEN                                   !(K->0)
            FBES(I,IX)=XX(IX)**LL
            GBES(I,IX)=-1./((2*LL+1.)*XX(IX)**(LL+1.D0))*AG**(2*LL+1.D0)
          ELSE  
            FBSS = SPHERICAL_J(LL, XG)
            GBSS = SPHERICAL_Y(LL, XG)
            FBES(I,IX)=K**(LL+0.5D0)*FBSS/(K**LL)
            GBES(I,IX)=-(GBSS*K**(LL+1.D0)*AG**(2*LL+1.D0))/(K**(LL+0.5D0))
            ! write(500+l,*) xg, FBSS, GBSS, fbes(i,IX), gbes(i,IX)
          ENDIF
        ENDDO
      ENDDO
    END SUBROUTINE SPHERICAL_BESSEL_FUNCTIONS

    SUBROUTINE PREPARE_INDECES()
      IMPLICIT NONE
      INTEGER II, JJ

      II = 0
      DO I=1, NEQC
        DO JJ=1, VARIATIONAL_PARAMS%NNL
          II = II + 1
          ICONT(I, JJ) = II
        ENDDO
      ENDDO
    END SUBROUTINE PREPARE_INDECES

  END SUBROUTINE ASYMPTOTIC_CORE_MATRIX_ELEMENTS




  SUBROUTINE ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS(AM, AM1, AM2, AM3)
    USE gsl_bessel
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NCH, NCH), INTENT(OUT) :: AM, AM1, AM2, AM3
    INTEGER, PARAMETER :: NNRAA = 200
    DOUBLE PRECISION :: H, H5, RANGE, R
    DOUBLE PRECISION :: AF, EPS, YY(NNRAA)
    DOUBLE PRECISION, DIMENSION(NCH_MAX,NNRAA) :: GBES, GBES0, GBES1, GBES2, FBES, HNOR
    INTEGER :: I, IX, L, IAB, IAK
    DOUBLE PRECISION :: XG, AG, BG
    DOUBLE PRECISION :: GBSS, GBSS1, GBSS2, FBSS
    INTEGER :: S, J
    DOUBLE PRECISION :: VPW(NCH_MAX,NCH_MAX)
    DOUBLE PRECISION, DIMENSION(NNRAA) :: FUN, FUN1, FUN2, FUN3
    DOUBLE PRECISION :: AXX, AXX3, AKE, AKE3, APE, APE1, APE2, APE3
    DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: AXXM, AXXM3, AKEM, AKEM3, APEM, APEM1, APEM2, APEM3, VER

    INTEGER, SAVE :: NX
    DOUBLE PRECISION, SAVE :: XX(NNRAA), AJ(NNRAA), A(NNRAA), B(NNRAA)
    DOUBLE PRECISION, SAVE :: VV(NNRAA,NCH_MAX,NCH_MAX)
    LOGICAL, SAVE :: FIRST_CALL = .TRUE.
    DOUBLE PRECISION, EXTERNAL :: B5

    H = VARIATIONAL_PARAMS%H
    H5= H/22.5D0
    RANGE = VARIATIONAL_PARAMS%RANGE
    AF = VARIATIONAL_PARAMS%AF
    IF (FIRST_CALL) CALL EXPONENTIALLY_GROWING_GRID(H, AF, RANGE, XX, AJ, NNR, NX)
    VARIATIONAL_PARAMS%RANGE = RANGE
    !  WRITE(*,*)'PRIMO E ULTIMO PUNTO =',XX(1),XX(NX)

    EPS = VARIATIONAL_PARAMS%EPS
    IF (FIRST_CALL) THEN
      DO I=1, NX
        R = XX(I)
        YY(I) = K*R
        A(I)  = 1.0 - DEXP(-EPS*R)
        B(I)  = EPS * DEXP(-EPS*R)
      ENDDO
    ELSE
      DO I=1, NX
        YY(I) = K*XX(I)
      ENDDO
    ENDIF

  !definisco funzioni bessel regolarizzate e con giuste dimensioni (K**(l+0.5d0) ed andamenti asintotici  
    DO I=1, NEQ
      L=LC(I)
      !  WRITE(*,*) "Preparing Bessels for L = ", L, " with k = ", K
      DO IX=1, NX                       
        XG=YY(IX)
        AG=A(IX)
        BG=B(IX)
        IF(K.LE.1.D-8) THEN
          FBES(I,IX)=XX(IX)**L
          GBES(I,IX)=-1./((2*L+1.)*XX(IX)**(L+1.D0))*AG**(2*L+1.D0)
          GBES0(I,IX)=1./((2*L+1.)*XX(IX)**(L+1.D0))*(EPS*BG*(2*L+1.)*((2*L+1.)*(BG/EPS)-1.) &
                    +2*(2*L+1.)*BG*(AG/XX(IX))-L*(L+1.)*(AG/XX(IX))**2)
          GBES1(I,IX)=-2.*AG*((2*L+1.)*BG+AG/XX(IX))*(L+1.)/((2*L+1.)*XX(IX)**(L+2.))
          GBES2(I,IX)=AG**2*(L+1.)*(L+2.)/((2*L+1.)*XX(IX)**(L+3.))
          HNOR(I,IX)=AG**(2*L-1.)
          WRITE(200+L,*) XX(IX), FBES(I,IX), GBES(I,IX), GBES1(I,IX), GBES2(I,IX)
        ELSE
          FBSS = SPHERICAL_J(L, XG)
          GBSS = SPHERICAL_Y(L, XG)
          GBSS1= SPHERICAL_YP(L, XG)
          GBSS2= SPHERICAL_YPP(L, XG)
          
          FBES(I,IX)=K**(L+0.5D0)*FBSS/(K**L)
          GBES(I,IX)=-(GBSS*K**(L+1.D0)*AG**(2*L+1.D0))/(K**(L+0.5D0))
          GBES0(I,IX)=GBSS*(EPS*BG*(2*L+1.)*((2*L+1.)*(BG/EPS)-1.) &
                    +2*(2*L+1.)*BG*(AG/XX(IX))-L*(L+1.)*(AG/XX(IX))**2)
          GBES1(I,IX)=GBSS1*2.*K*AG*((2*L+1.)*BG + AG/XX(IX))                   
          GBES2(I,IX)=(K**2)*(AG**2)*GBSS2
          HNOR(I,IX)=(K**(L+1.D0))*(AG**(2*L-1.))/(K**(L+0.5D0))
          ! WRITE(200+L,*) XX(IX), FBES(I,IX), GBES(I,IX), GBES1(I,IX), GBES2(I,IX), HNOR(I,IX)
        ENDIF
      ENDDO
    ENDDO

  ! si tabula il potenziale per i vari canali
    L = LC(1)
    S = VARIATIONAL_PARAMS%S
    J = VARIATIONAL_PARAMS%J
    IF (FIRST_CALL) THEN
      DO I=1, NX
        R = XX(I)
        CALL AV18PW90(1,L,S,J,T,T1Z,T2Z,R,VPW,VARIATIONAL_PARAMS%LEMP)
        VV(I,1,1) = VPW(1,1)
        VV(I,1,2) = VPW(1,2)
        VV(I,2,1) = VPW(2,1)
        VV(I,2,2) = VPW(2,2)
      ENDDO
      FIRST_CALL = .FALSE.
    ENDIF

  !SI CALCOLANO ELEMENTI MATRICE  
    DO IAB=1,NEQ          
    DO IAK=1,NEQ
  ! SI CALCOLA NORMA DEL CASO REGOLARE-IRREGOLARE (AXX), CASO IRREGOLARE-IRREGOLARE (AXX3)  
      AXX=0.D0
      AXX3=0.D0
      AXXM(IAB,IAK)=0.D0
      AXXM3(IAB,IAK)=0.D0
      IF(IAB.EQ.IAK)THEN 
          FUN(1)=0.D0
          FUN3(1)=0.D0                                     
        DO I=1,NX
          FUN(I+1)= AJ(I)*FBES(IAB,I)*GBES(IAK,I)
          FUN3(I+1)=AJ(I)*GBES(IAB,I)*GBES(IAK,I)                                             
        ENDDO
        AXX=  VARIATIONAL_PARAMS%E * B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)
        AXXM(IAB,IAK)=AXX
        AXX3= VARIATIONAL_PARAMS%E * B5(1,NX,0,0,H5,0.D0,0.D0,FUN3,1)
        AXXM3(IAB,IAK)=AXX3
      ENDIF    
      IF (IAB==1 .AND. IAK==1) THEN
        !  WRITE(*,*) 'NORM RI(1,1)', AXX
        !  WRITE(*,*) 'NORM II(1,1)', AXX3
      ENDIF

  ! SI CALCOLA ENERGIA CINETICA DEL CASO REGOLARE-IRREGOLARE (AKE) E CASO IRREGOLARE-IRREGOLARE (AKE3)
      AKE=0.D0
      AKE3=0.D0
      AKEM(IAB,IAK)=0.D0 
      AKEM3(IAB,IAK)=0.D0       
      IF(IAB.EQ.IAK)THEN
        FUN(1)=0.D0
        FUN3(1)=0.D0
        DO I=1,NX
          FUN(I+1) =AJ(I)*FBES(IAB,I)*HNOR(IAK,I)*(GBES2(IAK,I)+GBES1(IAK,I)+GBES0(IAK,I))
          FUN3(I+1)=AJ(I)*GBES(IAB,I)*HNOR(IAK,I)*(GBES2(IAK,I)+GBES1(IAK,I)+GBES0(IAK,I))
        ENDDO
        AKE=  HTM * B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)                                         
        AKEM(IAB,IAK)=AKE
        AKE3= HTM * B5(1,NX,0,0,H5,0.D0,0.D0,FUN3,1)
        AKEM3(IAB,IAK)=AKE3
      ENDIF 
      IF (IAB==1 .AND. IAK==1) THEN
        !  WRITE(*,*) 'KINETIC RI(1,1)', AKE
        !  WRITE(*,*) 'KINETIC II(1,1)', AKE3
      ENDIF

  ! SI CALCOLA ENERGIA POTENZIALE DEL CASO REGOLARE-IRREGOLARE (APE), 
  !      IRREGOLARE-REGOLARE (APE1), REGOLARE-REGOLARE (APE2), IRREGOLARE-IRREGOLARE (APE3)
      APE=0.D0
      APE1=0.D0
      APE2=0.D0
      APE3=0.D0
      FUN(1)=0.D0 
      FUN1(1)=0.D0
      FUN2(1)=0.D0 
      FUN3(1)=0.D0                                          
      DO I=1,NX
        FUN(I+1) =AJ(I)*FBES(IAB,I)*GBES(IAK,I)*VV(I,IAB,IAK) 
        FUN1(I+1)=AJ(I)*GBES(IAB,I)*FBES(IAK,I)*VV(I,IAB,IAK)  
        FUN2(I+1)=AJ(I)*FBES(IAB,I)*FBES(IAK,I)*VV(I,IAB,IAK)
        FUN3(I+1)=AJ(I)*GBES(IAB,I)*GBES(IAK,I)*VV(I,IAB,IAK)
      ENDDO
      APE=  B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)
      APEM(IAB,IAK)=APE
      APE1= B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)
      APEM1(IAB,IAK)=APE1
      APE2= B5(1,NX,0,0,H5,0.D0,0.D0,FUN2,1)
      APEM2(IAB,IAK)=APE2
      APE3= B5(1,NX,0,0,H5,0.D0,0.D0,FUN3,1)
      APEM3(IAB,IAK)=APE3
      IF (IAB==1 .AND. IAK==1) THEN
        !  WRITE(*,*) 'POTENTIAL RR(1,1)', APE2
        !  WRITE(*,*) 'POTENTIAL RI(1,1)', APE
        !  WRITE(*,*) 'POTENTIAL IR(1,1)', APE1
        !  WRITE(*,*) 'POTENTIAL II(1,1)', APE3
      ENDIF



  ! SI CALCOLA HAMILTONIANA PER I VARI CASI:REGOLARE-IRREGOLARE(AM), IRREGOLARE-REGOLARE(AM1),
  !         REGOLARE-REGOLARE(AM2),IRREGOLARE-IRREGOLARE(AM3)
      AM(IAB,IAK)= 1.D0/HTM*(AKEM(IAB,IAK)+APEM(IAB,IAK)-AXXM(IAB,IAK)) 
      AM1(IAB,IAK)=1.D0/HTM*APEM1(IAB,IAK)
      AM2(IAB,IAK)=1.D0/HTM*APEM2(IAB,IAK)
      AM3(IAB,IAK)=1.D0/HTM*(AKEM3(IAB,IAK)+APEM3(IAB,IAK)-AXXM3(IAB,IAK))

      VER(IAB,IAK)=(AM1(IAB,IAK)-AM(IAB,IAK)) 
    ENDDO 
    ENDDO 

    !  WRITE(*,*)
    !  WRITE(*,*)'A-A MATRIX'      
    DO IAB=1,NEQ
    DO IAK=1,NEQ
      !  WRITE(*,'(2I6,4D17.7)') IAB,IAK,AM2(IAB,IAK),AM(IAB,IAK),AM1(IAB,IAK),AM3(IAB,IAK)
      !  WRITE(*,*)"1=",VER(IAB,IAK)
    END DO
    END DO
  END SUBROUTINE ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS

  SUBROUTINE PRINT_DIVIDER()
    IMPLICIT NONE
    !  WRITE(*,*) '====================================================================================='
  END SUBROUTINE PRINT_DIVIDER

END MODULE SCATTERING_NN_VARIATIONAL