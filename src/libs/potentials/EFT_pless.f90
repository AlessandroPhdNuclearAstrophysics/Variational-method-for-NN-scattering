!   CLO( 0)  => c10
!   CLO( 1)  => c01
!   CNLO(1)  => c2c00
!   CNLO(2)  => c2c10
!   CNLO(3)  => c2c01
!   CNLO(4)  => c2c11
!   CNLO(5)  => c2t0
!   CNLO(6)  => c2t1
!   CNLO(7)  => c2b
!   CIT( 0)  => c0ct

MODULE EFT_PLESS
  USE REALLOCATE_UTILS
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  PRIVATE

  DOUBLE PRECISION, PARAMETER :: HTC = 197.32697D0

  TYPE, PUBLIC :: LECS_EFT_PLESS
    INTEGER :: ILB = -1
    INTEGER :: ORDER = -1
    DOUBLE PRECISION :: RC(0:1,0:1) = 0.D0
    DOUBLE PRECISION :: CLO(0:1)   = 0.D0
    DOUBLE PRECISION :: CNLO(7) = 0.D0
    DOUBLE PRECISION :: CN3LO(11)= 0.D0
    DOUBLE PRECISION :: CIT(0:4)= 0.D0
  END TYPE LECS_EFT_PLESS

  LOGICAL, PRIVATE :: FIRST_CALL = .TRUE., LECS_SET = .FALSE.
  INTEGER :: NMODELS = -1
  TYPE(LECS_EFT_PLESS), ALLOCATABLE :: LECS_ALL(:)
  TYPE(LECS_EFT_PLESS) :: LECS
  TYPE(SCATTERING_CHANNEL) :: CHANNEL
  INTEGER :: ORDER =-1

  INTEGER          :: I2   (2,2)
  INTEGER          :: LS   (2,2)
  INTEGER          :: L2   (2,2)
  DOUBLE PRECISION :: S12  (2,2)
  INTEGER          :: TZ_OP

  PUBLIC :: EFT_PLESS_PW, LECS_TO_ST_LECS, ST_LECTS_TO_LECS
  PUBLIC :: GET_LECS, PRINT_LECS
  PRIVATE:: SET_ALL_LECS, SET_LECS, SET_OPERATORS
  PRIVATE:: PREPARE, S12_OPERATOR, LS_OPERATOR, L2_OPERATOR
  PRIVATE:: CD_OP, CR, CR_1, CR_2, CR_3, CR_4, CR_5, CR_6, CR_7

CONTAINS 
  SUBROUTINE SET_ALL_LECS
    IMPLICIT NONE
    INTEGER :: I, IOS
    CHARACTER(LEN=256) :: FILENAME
    CHARACTER(LEN=1024) :: LINE
    INTEGER :: UNIT
    INTEGER :: ILB
    DOUBLE PRECISION :: RC(0:1,0:1), CLO(0:1), CNLO(7), CN3LO(11), CIT(0:4)

    IF (LECS_SET) RETURN

    ! Set the filename (you can change this as needed)
    FILENAME = 'src/libs/potentials/lecs_eft.dat'

    ! Open the file
    UNIT = 100
    OPEN(NEWUNIT=UNIT, FILE=FILENAME, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (ios /= 0) THEN
      PRINT *, 'Error opening file: ', TRIM(FILENAME),"   ->   ", ios
      STOP
    END IF

    ! Read the number of lines (number of LECS)
    READ(UNIT, *, IOSTAT=IOS) NMODELS
    IF (IOS /= 0) THEN
      PRINT *, 'Error reading number of lines from file.'
      CLOSE(UNIT)
      STOP
    END IF

    ! ALLOCATE THE ARRAY
    IF (ALLOCATED(LECS_ALL)) DEALLOCATE(LECS_ALL)
    ALLOCATE(LECS_ALL(NMODELS))

    ! Read each LECS_EFT_PLESS entry
    DO I = 1, NMODELS
      READ(UNIT, *, IOSTAT=IOS) ILB, RC(0,0), RC(1,0), RC(0,1), RC(1,1), CLO(1), CLO(0), CNLO, CN3LO, CIT
      LECS_ALL(I)%ILB = ILB
      LECS_ALL(I)%RC  = RC
      LECS_ALL(I)%CLO = CLO
      LECS_ALL(I)%CNLO = CNLO
      LECS_ALL(I)%CN3LO = CN3LO
      LECS_ALL(I)%CIT = CIT
      CALL SET_LECS_ORDER(LECS_ALL(I))
      IF (IOS /= 0) THEN
        PRINT *, 'Error reading LECS entry at line ', i
        CLOSE(UNIT)
        STOP
      END IF
    END DO
    CLOSE(UNIT)
    LECS_SET = .TRUE.
  END SUBROUTINE

  SUBROUTINE SET_LECS_ORDER(LECS_IN)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(INOUT) :: LECS_IN
    INTEGER :: I, J, K
    INTEGER :: NLO = SIZE(LECS_IN%CNLO)
    INTEGER :: N3LO = SIZE(LECS_IN%CN3LO)

    LECS_IN%ORDER = 0
    DO I=1, NLO
      IF (LECS_IN%CNLO(I) /= 0.D0) LECS_IN%ORDER = 1
    ENDDO

    DO J=1, N3LO
      IF (LECS_IN%CN3LO(J) /= 0.D0) LECS_IN%ORDER = 3
    ENDDO

  END SUBROUTINE SET_LECS_ORDER

  SUBROUTINE SET_LECS(ILB)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB
    IF (ILB <= 0 .AND. ILB > NMODELS) THEN
      WRITE(*,*) "THIS MODEL (ILB) IS NOT RECOGNIZED"
      STOP
    ENDIF
    
    LECS = LECS_ALL(ILB)
  END SUBROUTINE SET_LECS

  FUNCTION S12_OPERATOR(L, S, J) RESULT(S12_RES)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L, S, J
    DOUBLE PRECISION :: S12_RES(2,2)

    S12_RES = 0
    IF (S==0) RETURN
    IF (L==J) THEN
      S12_RES(1,1) = 2.0
      RETURN
    ELSEIF (ABS(J-L)==1) THEN
      IF (J==0) THEN
        IF (J==L+1) THEN
          S12_RES(1,1) =-2*(J-1.D0)/(2*J+1.D0)
        ELSEIF (J==L-1) THEN
          S12_RES(1,1) =-2*(J+2.D0)/(2*J+1.D0)
        ELSE
          WRITE(*,*) "Wrong LSJ combination"
          STOP
        ENDIF
      ELSE
        S12_RES(1,1) =-2*(J-1.D0)/(2*J+1.D0)
        S12_RES(1,2) = 6*SQRT(J*(J+1.D0))/(2*J+1.D0)
        S12_RES(2,1) = 6*SQRT(J*(J+1.D0))/(2*J+1.D0)
        S12_RES(2,2) =-2*(J+2.D0)/(2*J+1.D0)
      ENDIF
    ELSE
      WRITE(*,*) "Wrong LSJ combination"
      STOP
    ENDIF
  END FUNCTION S12_OPERATOR

  FUNCTION LS_OPERATOR(LMIN, S, J) RESULT(LS_OP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LMIN, S, J
    INTEGER :: LS_OP(2,2)
    
    LS_OP = 0
    IF (LMIN==J .OR. J==0) THEN
      LS_OP(1,1) = (J*(J+1)-S*(S+1)-LMIN*(LMIN+1))/2
    ELSEIF ((J-LMIN)==1 .AND. S==1) THEN
      LS_OP(1,1) = (J*(J+1)-S*(S+1)-LMIN*(LMIN+1))/2
      LS_OP(2,2) = (J*(J+1)-S*(S+1)-(LMIN+2)*(LMIN+3))/2
    ELSE
      WRITE(*,*) "Wrong LSJ combination"
      STOP
    ENDIF
  END FUNCTION LS_OPERATOR

  FUNCTION L2_OPERATOR(L) RESULT(L2_MAT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    INTEGER :: L2_MAT(2,2)
    L2_MAT = 0
    L2_MAT(1,1) = L*(L+1)
    IF (IS_CHANNEL_COUPLED(CHANNEL)) L2_MAT(2,2) = (L+2)*(L+3)
  END FUNCTION L2_OPERATOR

  FUNCTION CD_OP(T, TZ) RESULT(CD)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: T, TZ
    INTEGER :: CD

    CD = 0
    IF (T==0) RETURN
    IF (ABS(TZ)==1) THEN
      CD = 2
    ELSEIF (TZ==0) THEN
      CD =-4
    ELSE
      WRITE(*,*) "Wrong T Tz combination"

    ENDIF
  END FUNCTION CD_OP


  SUBROUTINE SET_OPERATORS()
    IMPLICIT NONE
    INTEGER :: L, S, J, T, TZ, I

    I2 = 0
    DO I=1, 2
      I2(I,I) = 1
    ENDDO

    L = GET_CHANNEL_L (CHANNEL, 1)
    S = GET_CHANNEL_S (CHANNEL, 1)
    T = GET_CHANNEL_T (CHANNEL, 1)
    TZ= GET_CHANNEL_TZ(CHANNEL)
    J = GET_CHANNEL_J (CHANNEL)
    
    S12   = S12_OPERATOR(L, S, J)
    LS    = LS_OPERATOR (L, S, J)
    L2    = L2_OPERATOR (L)
    TZ_OP = CD_OP       (T,TZ)

  END SUBROUTINE SET_OPERATORS



  SUBROUTINE PREPARE(L, S, J, TZ, ILB)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: MINR = 1.D-1, ZERO = 1.D-10
    INTEGER, INTENT(IN) :: L, S, J, TZ, ILB
    LECS%ILB = ILB
    CALL SET_LECS(ILB)
    CALL SET_OPERATORS
    IF (SUM(ABS(LECS%RC)) < MINR) STOP "RC not set"
  END SUBROUTINE PREPARE

  SUBROUTINE EFT_PLESS_PW(ILB, L, S, J, T1Z, T2Z, R, VPW, LEMP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB, L, S, J, T1Z, T2Z, LEMP
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT):: VPW(2,2)
    INTEGER :: T, TZ
    TYPE(SCATTERING_CHANNEL) :: CHANNEL_NEW
    LOGICAL :: IS_NEW_CHANNEL
    DOUBLE PRECISION :: RC
    INTEGER :: LS2(2,2)
    
    VPW = 0
    TZ = (T1Z+T2Z)/2
    T = MOD(MOD((L+S),2)+1,2)
    CALL SET_CHANNEL(CHANNEL_NEW, J, L, S, TZ)
    IS_NEW_CHANNEL = .NOT.IS_SAME_CHANNEL(CHANNEL, CHANNEL_NEW)
    CHANNEL = CHANNEL_NEW

    IF (FIRST_CALL) THEN
      CALL SET_ALL_LECS
      CALL PREPARE(L, S, J, TZ, ILB)
      FIRST_CALL = .FALSE.
    ELSEIF (IS_NEW_CHANNEL) THEN
      CALL PREPARE(L, S, J, TZ, ILB)
    ENDIF
    
    RC = LECS%RC(S,T)
    ORDER = LECS%ORDER
    
    IF ( S==0 .AND. T==0 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RETURN
      CASE (1)
        VPW = LECS%CNLO(1) * CR_1(R, RC) * I2
      CASE (3)
        VPW = LECS%CNLO(1) * CR_1(R, RC) * I2
        VPW = VPW + LECS%CN3LO(1)  * CR_4(R, RC) *I2 &
                  + LECS%CN3LO(10) * CR_7(R, RC) *L2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=0 AND T=0"
      END SELECT
    ENDIF

    IF ( S==1 .AND. T==0 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        VPW = LECS%CLO(T) * I2
      CASE (1)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(2) * CR_1(R, RC) * I2 &
                  + LECS%CNLO(5) * CR_2(R, RC) * S12 &
                  + LECS%CNLO(7) * CR_3(R, RC) * LS
      CASE (3)
        LS2 = MATMUL(LS, LS)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(2) * CR_1(R, RC) * I2 &
                  + LECS%CNLO(5) * CR_2(R, RC) * S12 &
                  + LECS%CNLO(7) * CR_3(R, RC) * LS
        VPW = VPW + LECS%CN3LO(2) * CR_4(R, RC) * I2 &
                  + LECS%CN3LO(5) * CR_5(R, RC) * S12 &
                  + LECS%CN3LO(7) * CR_6(R, RC) * LS &
                  + LECS%CN3LO(9) * CR_7(R, RC) * LS2 &
                  + LECS%CN3LO(11)* CR_7(R, RC) * L2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=1 AND T=0"
      END SELECT
    ENDIF

    IF ( S==0 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        VPW = LECS%CLO(T) * I2
      CASE (1)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(3) * CR_1(R, RC) * I2 &
                  + LECS%CIT(0)                * TZ_OP * I2
      CASE (3)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(3)  * CR_1(R, RC) * I2 &
                  + LECS%CIT(0)                 * TZ_OP * I2
        VPW = VPW + LECS%CN3LO(3) * CR_4(R, RC)     * I2 &
                  + LECS%CN3LO(10)* CR_7(R, RC)     * L2 &
                  + LECS%CIT(1)   * CR_1(R, RC)     * TZ_OP * I2 
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=0 AND T=1"
      END SELECT
    ENDIF

    IF ( S==1 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RETURN
      CASE (1)
        VPW =   LECS%CNLO(4) * CR_1(R, RC) *  I2 &
              + LECS%CNLO(6) * CR_2(R, RC) *  S12 &
              + LECS%CNLO(7) * CR_3(R, RC) *  LS &
              + LECS%CIT(0)                *  TZ_OP * I2
      CASE (3)
        LS2 = MATMUL(LS, LS)
        VPW =   LECS%CNLO(4) * CR_1(R, RC) *  I2 &
              + LECS%CNLO(6) * CR_2(R, RC) *  S12 &
              + LECS%CNLO(7) * CR_3(R, RC) *  LS &
              + LECS%CIT(0)                *  TZ_OP * I2
        VPW =   VPW &
              + LECS%CN3LO(4) * CR_4(R, RC)    * I2 &
              + LECS%CN3LO(6) * CR_5(R, RC)    * S12 &
              + LECS%CN3LO(8) * CR_6(R, RC)    * LS &
              + LECS%CN3LO(9) * CR_7(R, RC)    * LS2 &
              + LECS%CN3LO(11)* CR_7(R, RC)    * L2 &
              +(  LECS%CIT(2) * CR_1(R, RC)    * I2  + &
                  LECS%CIT(3) * CR_2(R, RC)    * S12 + &
                  LECS%CIT(4) * CR_3(R, RC)    * LS    & 
                                                        ) * TZ_OP
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=1 AND T=1"
      END SELECT
    ENDIF
    
    VPW = VPW * CR(R, RC)
    VPW = VPW * HTC

    IF (GET_CHANNEL_NCH(CHANNEL) == 1 ) THEN
      VPW(1,2) = 0.D0
      VPW(2,1) = 0.D0
      VPW(2,2) = 0.D0
    ENDIF
    RETURN
  END SUBROUTINE EFT_PLESS_PW


  PURE ELEMENTAL FUNCTION CR(R, RC) RESULT(CR_OUT)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: CR_OUT
    CR_OUT = DEXP(-R**2/RC**2)/(PI**(3.D0/2.D0)*RC**3)
  END FUNCTION CR

  PURE ELEMENTAL FUNCTION CR_1(R, RC) RESULT(F1)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F1
    F1 = (6*RC**2 - 4*R**2)/(RC**4)
  END FUNCTION CR_1

  PURE ELEMENTAL FUNCTION CR_2(R, RC) RESULT(F2)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F2
    F2 =-4*R**2/RC**4
  END FUNCTION CR_2

  PURE ELEMENTAL FUNCTION CR_3(R, RC) RESULT(F3)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F3
    F3 = 2.D0/RC**2
  END FUNCTION CR_3

  PURE ELEMENTAL FUNCTION CR_4(R, RC) RESULT(F4)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F4
    F4 = 4.D0*( 4*R**4 - 20*RC**2*R**2 + 15*RC**4 )/RC**8
  END FUNCTION CR_4

  PURE ELEMENTAL FUNCTION CR_5(R, RC) RESULT(F5)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F5
    F5 = 8*R**2*( 2*R**2 - 7*RC**2 )/RC**8
  END FUNCTION CR_5

  PURE ELEMENTAL FUNCTION CR_6(R, RC) RESULT(F6)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F6
    F6 = ( 20*RC**2 - 8*R**2 )/RC**6
  END FUNCTION CR_6

  PURE ELEMENTAL FUNCTION CR_7(R, RC) RESULT(F7)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F7
    F7 =-CR_3(R, RC)**2
  END FUNCTION CR_7


  PURE FUNCTION LECS_TO_ST_LECS(LECS_OLD) RESULT(ST_LECS)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_OLD
    TYPE(LECS_EFT_PLESS) :: ST_LECS
    DOUBLE PRECISION :: C(7)
    DOUBLE PRECISION :: D(11)
    DOUBLE PRECISION :: CIT(0:4)

    C   = LECS_OLD%CNLO
    D   = LECS_OLD%CN3LO
    CIT = LECS_OLD%CIT

    ST_LECS%ILB       = LECS_OLD%ILB
    ST_LECS%RC        = LECS_OLD%RC  
    ST_LECS%CLO       = LECS_OLD%CLO
    ST_LECS%CNLO(1)   = C(1) - 3*C(2) - 3*C(3) + 9*C(4)
    ST_LECS%CNLO(2)   = C(1) - 3*C(2) + C(3) - 3*C(4)
    ST_LECS%CNLO(3)   = C(1) + C(2) - 3*C(3) - 3*C(4)
    ST_LECS%CNLO(4)   = C(1) + C(2) + C(3) + C(4)
    ST_LECS%CNLO(5)   = C(5) - 3*C(6)
    ST_LECS%CNLO(6)   = C(5) + C(6)
    ST_LECS%CNLO(7)   = C(7)
    ST_LECS%CN3LO( 1) = D(1) - 3*D(2) - 3*D(3) + 9*D(4)
    ST_LECS%CN3LO( 2) = D(1) - 3*D(2) + D(3) - 3*D(4)
    ST_LECS%CN3LO( 3) = D(1) + D(2) - 3*D(3) - 3*D(4)
    ST_LECS%CN3LO( 4) = D(1) + D(2) + D(3) + D(4)
    ST_LECS%CN3LO( 5) = D(5) - 3*D(6)
    ST_LECS%CN3LO( 6) = D(5) + D(6)
    ST_LECS%CN3LO( 7) = D(7) - 3*D(8)
    ST_LECS%CN3LO( 8) = D(7) + D(8)
    ST_LECS%CN3LO( 9) = D(9)
    ST_LECS%CN3LO(10) = D(10) - 3*D(11)
    ST_LECS%CN3LO(11) = D(10) + D(11)
    ST_LECS%CIT(0)    = CIT(0)
    ST_LECS%CIT(1)    = CIT(1) - 3*CIT(2)
    ST_LECS%CIT(2)    = CIT(1) + CIT(2)
    ST_LECS%CIT(3)    = CIT(3)
    ST_LECS%CIT(4)    = CIT(4)
  END FUNCTION LECS_TO_ST_LECS

  FUNCTION ST_LECTS_TO_LECS(ST_LECS) RESULT(OP_LECS)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: ST_LECS
    TYPE(LECS_EFT_PLESS) :: OP_LECS
    DOUBLE PRECISION :: C1(7)
    DOUBLE PRECISION :: D1(11)
    DOUBLE PRECISION :: CIT(0:4)

    C1   = ST_LECS%CNLO
    D1   = ST_LECS%CN3LO
    CIT  = ST_LECS%CIT
    
    OP_LECS%ILB= ST_LECS%ILB
    OP_LECS%RC = ST_LECS%RC
    OP_LECS%CLO = ST_LECS%CLO
    OP_LECS%CNLO(1)   = (C1(1) + 3*C1(2) + 3*C1(3) + 9*C1(4))/16.D0
    OP_LECS%CNLO(2)   = (-C1(1) - 3*C1(2) + C1(3) + 3*C1(4))/16.D0
    OP_LECS%CNLO(3)   = (-C1(1) + C1(2) - 3*C1(3) + 3*C1(4))/16.D0
    OP_LECS%CNLO(4)   = (C1(1) - C1(2) - C1(3) + C1(4))/16.D0
    OP_LECS%CNLO(5)   = (C1(5) + 3*C1(6))/4.D0
    OP_LECS%CNLO(6)   = (-C1(5) + C1(6))/4.D0
    OP_LECS%CNLO(7)   = C1(7)
    OP_LECS%CN3LO(1)  = (D1(1) + 3*(D1(2) + D1(3) + 3*D1(4)))/16.D0
    OP_LECS%CN3LO(2)  = (-D1(1) - 3*D1(2) + D1(3) + 3*D1(4))/16.D0
    OP_LECS%CN3LO(3)  = (-D1(1) + D1(2) - 3*D1(3) + 3*D1(4))/16.D0
    OP_LECS%CN3LO(4)  = (D1(1) - D1(2) - D1(3) + D1(4))/16.D0
    OP_LECS%CN3LO(5)  = (D1(5) + 3*D1(6))/4.D0
    OP_LECS%CN3LO(6)  = (-D1(5) + D1(6))/4.D0
    OP_LECS%CN3LO(7)  = (D1(7) + 3*D1(8))/4.D0
    OP_LECS%CN3LO(8)  = (-D1(7) + D1(8))/4.D0
    OP_LECS%CN3LO(9)  = D1(9)
    OP_LECS%CN3LO(10) = (D1(10) + 3*D1(11))/4.D0
    OP_LECS%CN3LO(11) = (-D1(10) + D1(11))/4.D0
    OP_LECS%CIT(0)    = CIT(0)
    OP_LECS%CIT(1)    = (CIT(1) + 3*CIT(2))/4.D0
    OP_LECS%CIT(2)    = (-CIT(1) + CIT(2))/4.D0
    OP_LECS%CIT(3)    = CIT(3)
    OP_LECS%CIT(4)    = CIT(4)
  END FUNCTION ST_LECTS_TO_LECS

  FUNCTION GET_LECS(ILB) RESULT(LECS_OUT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB
    TYPE(LECS_EFT_PLESS) :: LECS_OUT

    IF (.NOT.LECS_SET) CALL SET_ALL_LECS

    IF (ILB <= 0 .AND. ILB > NMODELS) THEN
      WRITE(*,*) "THIS MODEL (ILB) IS NOT RECOGNIZED"
      STOP
    ENDIF

    LECS_OUT = LECS_ALL(ILB)
  END FUNCTION GET_LECS

  SUBROUTINE PRINT_LECS(LECS_IN)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_IN
    INTEGER :: I, J
    CHARACTER(LEN=16) :: ORD

    SELECT CASE (LECS_IN%ORDER)
    CASE (0)
      ORD = "LO"
    CASE (1)
      ORD = "NLO"
    CASE (3)
      ORD = "N3LO"
    CASE DEFAULT
      ORD = "UNKNOWN"
    END SELECT

    PRINT *, "LECS_EFT_PLESS:"
    PRINT *, "  ILB   = ", LECS_IN%ILB
    PRINT *, "  ORDER =          ", ORD
    PRINT *, "  RC:"
    DO I = 0, 1
      WRITE(*,'(A,2F12.8)') "    ", LECS_IN%RC(I,0), LECS_IN%RC(I,1)
    END DO
    PRINT *, "  CLO:"
    WRITE(*,'(A,2ES18.8E2)') "    ", LECS_IN%CLO(1), LECS_IN%CLO(0)
    PRINT *, "  CNLO:"
    WRITE(*,'(A,7ES18.8E2)') "    ", (LECS_IN%CNLO(J), J=1,7)
    PRINT *, "  CN3LO:"
    WRITE(*,'(A,7ES18.8E2)') "    ", (LECS_IN%CN3LO(J), J=1,7)
    WRITE(*,'(A,4ES18.8E2)') "    ", (LECS_IN%CN3LO(J), J=8,11)
    PRINT *, "  CIT:"
    WRITE(*,'(A,5ES18.8E2)') "    ", (LECS_IN%CIT(J), J=0,4)
  END SUBROUTINE PRINT_LECS

END MODULE EFT_PLESS

