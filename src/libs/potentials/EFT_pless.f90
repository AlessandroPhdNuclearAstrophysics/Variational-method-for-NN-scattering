!   CLO( 1)  => c01
!   CLO( 2)  => c10
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

  TYPE, PUBLIC :: LECS_EFT_PLESS
    INTEGER :: ILB = -1
    DOUBLE PRECISION :: RC(2,2) = 0.D0
    DOUBLE PRECISION :: CLO(2)   = 0.D0
    DOUBLE PRECISION :: CNLO(7) = 0.D0
    DOUBLE PRECISION :: CN3LO(11)= 0.D0
    DOUBLE PRECISION :: CIT(0:4)= 0.D0
  END TYPE LECS_EFT_PLESS

  LOGICAL, PRIVATE :: FIRST_CALL = .TRUE.
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

CONTAINS 
  SUBROUTINE SET_ALL_LECS
    IMPLICIT NONE
    INTEGER :: I, IOS
    CHARACTER(LEN=256) :: FILENAME
    CHARACTER(LEN=1024) :: LINE
    INTEGER :: UNIT

    ! Set the filename (you can change this as needed)
    FILENAME = 'src/potentials/lecs_eft.dat'

    ! Open the file
    UNIT = 100
    OPEN(NEWUNIT=UNIT, FILE=FILENAME, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (ios /= 0) THEN
      PRINT *, 'Error opening file: ', FILENAME
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
      READ(UNIT, *, IOSTAT=IOS) LECS_ALL(I)
      IF (IOS /= 0) THEN
      PRINT *, 'Error reading LECS entry at line ', i
      CLOSE(UNIT)
      STOP
      END IF
    END DO
    CLOSE(UNIT)
  END SUBROUTINE

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
      IF (J==0 .AND. L<J) THEN
        S12_RES(1,1) =-2*(J-1.D0)/(2*J+1.D0)
      ELSEIF (J==0 .AND. L>J) THEN
        S12_RES(1,1) =-2*(J+2.D0)/(2*J+1.D0)
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
      LS_OP(1,1) = (J*(J+1)-S*(S+1)-LMIN*(LMIN+1))
    ELSEIF ((J-LMIN)==1 .AND. S==1) THEN
      LS_OP(1,1) = (J*(J+1)-S*(S+1)-LMIN*(LMIN+1))
      LS_OP(2,2) = (J*(J+1)-S*(S+1)-(LMIN+2)*(LMIN+3))
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


  SUBROUTINE PREPARE_OPERATORS()
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

  END SUBROUTINE PREPARE_OPERATORS!
    
    
  

  SUBROUTINE PREPARE(L, S, J, TZ, ILB)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: MINR = 1.D-1, ZERO = 1.D-10
    INTEGER, INTENT(IN) :: L, S, J, TZ, ILB
    LECS%ILB = ILB
    CALL SET_LECS(ILB)
    CALL SET_OPERATORS
    IF (SUM(ABS(LECS%RC)) < MINR) STOP "RC not set"
    IF (IS_ARRAY_ZERO(LECS%CLO, ZERO)) STOP "LO is not set"
    ORDER = 0
    IF (IS_ARRAY_ZERO(LECS%CNLO, ZERO) .AND. LECS%CIT(0) < ZERO) RETURN
    ORDER = 1
    IF (IS_ARRAY_ZERO(LECS%CN3LO, ZERO) .AND. IS_ARRAY_ZERO(LECS%CIT(1:), ZERO)) RETURN
    ORDER = 3
  CONTAINS
    FUNCTION IS_ARRAY_ZERO(ARRAY, MIN_VAL) RESULT(RES)
      IMPLICIT NONE 
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(:), MIN_VAL
      LOGICAL :: RES
      IF (SIZE(ARRAY) == 0) STOP "EFT_PLESS::IS_ARRAY_ZERO -> NO ARRAY"
      RES = SUM(ABS(ARRAY)) < MIN_VAL
    END FUNCTION IS_ARRAY_ZERO
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
    
    VPW = 0
    TZ = (T1Z+T2Z)/2
    T = MOD(MOD((L+S),2)+1,2)
    CALL SET_CHANNEL(CHANNEL_NEW, J, L, S, TZ)
    IS_NEW_CHANNEL = IS_SAME_CHANNEL(CHANNEL, CHANNEL_NEW)

    IF (FIRST_CALL) THEN
      CALL SET_ALL_LECS
      CALL PREPARE(L, S, J, TZ, ILB)
      FIRST_CALL = .FALSE.
    ELSEIF (IS_NEW_CHANNEL) THEN
      CALL PREPARE(L, S, J, TZ, ILB)
    ENDIF
    
    RC = LECS%RC(S,T)
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
      VPW = VPW * CR_FUN(R, RC)
    ENDIF

    IF ( S==1 .AND. T==0 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        VPW = LECS%CLO(2) * I2
      CASE (1)
        VPW = LECS%CLO(2) *I2
        VPW = VPW + LECS%CNLO(2) * CR_1(R, RC) * I2 &
                  + LECS%CNLO(5) * CR_2(R, RC) * S12 &
                  + LECS%CNLO(7) * CR_3(R, RC) * LS
      CASE (3)
        VPW = LECS%CLO(2) *I2
        VPW = VPW + LECS%CNLO(2) * CR_1(R, RC) * I2 &
                  + LECS%CNLO(5) * CR_2(R, RC) * S12 &
                  + LECS%CNLO(7) * CR_3(R, RC) * LS
        VPW = VPW + LECS%CN3LO(2) * CR_4(R, RC) * I2 &
                  + LECS%CN3LO(5) * CR_5(R, RC) * S12 &
                  + LECS%CN3LO(7) * CR_6(R, RC) * LS &
                  + LECS%CN3LO(9) * CR_7(R, RC) * MATMUL(LS,LS) &
                  + LECS%CN3LO(11)* CR_7(R, RC) *L2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=0 AND T=0"
      END SELECT
    ENDIF

    IF ( S==0 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        VPW = LECS%CLO(1) * I2
      CASE (1)
        VPW = LECS%CLO(1) *I2
        VPW = VPW + LECS%CNLO(3) * CR_1(R, RC) * I2 &
                  + LECS%CIT(0)                * TZ_OP * I2
      CASE (3)
        VPW = LECS%CLO(1) *I2
        VPW = VPW + LECS%CNLO(3) * CR_1(R, RC) * I2 &
                  + LECS%CIT(0)                * TZ_OP * I2
        VPW = VPW + LECS%CN3LO(3) * CR_4(R, RC)     * I2 &
                  + LECS%CN3LO(10)* CR_7(R, RC)     * L2 &
                  + LECS%CIT(1)                     * TZ_OP * I2 
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=0 AND T=0"
      END SELECT
    ENDIF

    IF ( S==1 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RETURN
      CASE (1)
        VPW =   LECS%CNLO(4) * I2 &
              + LECS%CNLO(6) * CR_5(R, RC) *  S12 &
              + LECS%CNLO(7) * CR_3(R, RC) *  LS &
              + LECS%CIT(0)                *  TZ_OP * I2
      CASE (3)
        VPW =   LECS%CNLO(4) * CR_1(R, RC) *  I2 &
              + LECS%CNLO(6) * CR_5(R, RC) *  S12 &
              + LECS%CNLO(7) * CR_3(R, RC) *  LS &
              + LECS%CIT(0)                *  TZ_OP * I2
        VPW =   VPW &
              + LECS%CN3LO(4) * CR_4(R, RC)    * I2 &
              + LECS%CN3LO(6) * CR_5(R, RC)    * S12 &
              + LECS%CN3LO(8) * CR_6(R, RC)    * LS &
              + LECS%CN3LO(9) * CR_7(R, RC)    *MATMUL(LS,LS) &
              + LECS%CN3LO(11)* CR_7(R, RC)    * L2 &
              +(LECS%CIT(2)*CR_1(R,RC) - LECS%CIT(3)*CR_2(R, RC) + LECS%CIT(4)*CR_3(R, RC) ) * TZ_OP * I2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=0 AND T=0"
      END SELECT
    ENDIF
    VPW = VPW/CR_FUN(R, RC)
    RETURN
  END SUBROUTINE EFT_PLESS_PW


  PURE ELEMENTAL FUNCTION CR_FUN(R, RC) RESULT(CR_OUT)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: CR_OUT
    CR_OUT = DEXP(-R**2/RC**2)/(PI**(3/2)*RC**3)
  END FUNCTION CR_FUN

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
    F4 = 4.D0*( 4*R**2 - 20*RC*R + 15*RC**4 )/RC**8
  END FUNCTION CR_4

  PURE ELEMENTAL FUNCTION CR_5(R, RC) RESULT(F5)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F5
    F5 = 8*R**2*( 2*R**2 - 7*RC**2 )/RC**6
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



END MODULE EFT_PLESS

