PROGRAM FIT_LOW_ENERGY_OBSERVABLES
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS
  USE OPERATING_SYSTEM_LINUX
  USE NUMBER_DIFFERENCES
  IMPLICIT NONE
  INTEGER, PARAMETER :: LMAX = 2
  INTEGER, PARAMETER :: JMAX = 2
  INTEGER, PARAMETER :: TZ = 0
  INTEGER, PARAMETER :: FIT_ORDER = 2
  CHARACTER(LEN=*), PARAMETER :: AV18_DIR = "output/AV18/"
  CHARACTER(LEN=*), PARAMETER :: AV18_KCOTD_DIR = "output/test_AV18_kcotd/"
  CHARACTER(LEN=*), PARAMETER :: EFT_PLESS_10_DIR = "output/EFT_pless_10_3MeV/"
  CHARACTER(LEN=*), PARAMETER :: EFT_PLESS_10_KCOTD_DIR = "output/test_EFT_pless_10_kcotd/"
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
  INTEGER :: NCHANNELS, ICH
  CHARACTER(LEN=8) :: CHANNEL_NAME
  INTEGER :: IE, NDATA, IOS, IEQ
  CHARACTER(LEN=256) :: LINE
  REAL(8), ALLOCATABLE :: ENERGIES(:)
  TYPE(PHASE_SHIFT_RESULT), ALLOCATABLE :: PHASE_SHIFT_ARRAY(:)
  DOUBLE PRECISION :: CONST_RESULT(2,FIT_ORDER+1)
  DOUBLE PRECISION, ALLOCATABLE :: K2(:), KCOTD(:,:)
  LOGICAL :: FITTED


  DOUBLE PRECISION, ALLOCATABLE :: AEXP(:,:), BEXP(:,:)

  
  ! AV18 TEST
  CALL CREATE_DIRECTORY(TRIM(AV18_KCOTD_DIR))
  WRITE(*,*) "========================================================"
  WRITE(*,*) "========================================================"
  PRINT *, "Testing AV18 potential..."
  WRITE(*,*) "========================================================"
  WRITE(*,*) "========================================================"

  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  NCHANNELS = SIZE(CHANNELS)

  CALL REALLOCATE(AEXP, NCHANNELS, 2)
  CALL REALLOCATE(BEXP, NCHANNELS, 2)
  AEXP(1,1) =  4.331D-2  ! 1S0
  AEXP(2,1) = -1.851D-1  ! 3S1
  AEXP(3,1) = -3.592D-1  ! 1P1
  AEXP(4,1) =  3.971D-1  ! 3P0
  AEXP(5,1) = -6.541D-1  ! 3P1
  AEXP(6,1) =  3.413D+0  ! 3P2
  AEXP(7,1) =  7.211D-1  ! 1D2
  AEXP(2,2) = -1.553D-1  ! 3D1
  AEXP(8,1) =  1.351D-1  ! 3D2
  AEXP(6,2) =  1.025D+0  ! 3F2

  BEXP(1,1) =  1.356D0   ! 1S0
  BEXP(2,1) =  8.521D-1  ! 3S1
  BEXP(3,1) = -3.270D0   ! 1P1
  BEXP(4,1) =  1.935D0   ! 3P0
  BEXP(5,1) = -4.373D0   ! 3P1
  BEXP(6,1) = -4.090D0   ! 3P2
  BEXP(7,1) =  7.384D0   ! 1D2
  BEXP(2,2) = -2.010D0   ! 3D1
  BEXP(8,1) =  1.419D0   ! 3D2
  BEXP(6,2) =  1.339D1   ! 3F2
  
  DO ICH = 1, NCHANNELS
    WRITE(*,*) "========================================================"
    CHANNEL_NAME = TRIM(GET_CHANNEL_NAME(CHANNELS(ICH)))
    PRINT *, "  Channel: ", CHANNEL_NAME
    OPEN(UNIT=10, FILE=TRIM(AV18_DIR)//"delta_"//TRIM(CHANNEL_NAME)//".dat", STATUS='OLD', ACTION='READ')

    ! First, count the number of data lines (not starting with #)
    NDATA = 0
    DO
      READ(10, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      IF (LEN_TRIM(LINE) == 0) CYCLE
      LINE = ADJUSTL(LINE)
      IF (LINE(1:1) == "#") CYCLE
      NDATA = NDATA + 1
    END DO
    REWIND(10)

    IF (ALLOCATED(ENERGIES)) THEN
      DEALLOCATE(ENERGIES)
    END IF
    ALLOCATE(ENERGIES(NDATA))
    IF (ALLOCATED(PHASE_SHIFT_ARRAY)) THEN
      DEALLOCATE(PHASE_SHIFT_ARRAY)
    END IF
    ALLOCATE(PHASE_SHIFT_ARRAY(NDATA))

    WRITE(*,*) "  Number of data points: ", NDATA
    IE = 0
    DO
      READ(10, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      IF (LEN_TRIM(LINE) == 0) CYCLE
      LINE = ADJUSTL(LINE)
      IF (LINE(1:1) == "#") CYCLE
      IE = IE + 1
      READ(LINE, *) ENERGIES(IE), PHASE_SHIFT_ARRAY(IE)%delta1_S, &
            PHASE_SHIFT_ARRAY(IE)%delta2_S, PHASE_SHIFT_ARRAY(IE)%epsilon_S
    END DO
    CLOSE(UNIT=10)


    IF (ALLOCATED(K2)) THEN
      DEALLOCATE(K2)
    END IF
    ALLOCATE(K2(NDATA))
    IF (ALLOCATED(KCOTD)) THEN
      DEALLOCATE(KCOTD)
    END IF
    ALLOCATE(KCOTD(2,NDATA))

    OPEN(UNIT=10, FILE=TRIM(AV18_KCOTD_DIR)//"k2_"//TRIM(CHANNEL_NAME)//".dat", STATUS='UNKNOWN', ACTION='WRITE')

    CONST_RESULT = 0.0D0
    FITTED = FIT_CHANNEL_LOW_ENERGY(CHANNELS(ICH), ENERGIES, PHASE_SHIFT_ARRAY, CONST_RESULT, FIT_ORDER, K2, KCOTD)
    IF (.NOT. FITTED) THEN
      WRITE(*,*) "  WARNING: Fitting failed for channel ", CHANNEL_NAME
      WRITE(*,*) "  This may be due to insufficient data points or a poor fit order."
    END IF

    DO IE = 1, NDATA
      WRITE(10, *) K2(IE), KCOTD(1,IE), KCOTD(2,IE)
    END DO


    WRITE(*,*) "  Fitted constants for channel ", CHANNEL_NAME, ":"
    DO IEQ = 1, GET_CHANNEL_NCH(CHANNELS(ICH))
      IF (IEQ==2) WRITE(*,*)
      WRITE(*,'(A)', ADVANCE='NO') "    Fit: "
      DO IE = 1, FIT_ORDER+1
        IF (IE > 1) WRITE(*,'(A)', ADVANCE='NO') " + "
        IF (IE == 1) THEN
          WRITE(*,'(F12.6)', ADVANCE='NO') CONST_RESULT(IEQ,IE)
        ELSE
          WRITE(*,'(F12.6,A,I0)', ADVANCE='NO') CONST_RESULT(IEQ, IE), " * X"
          IF (IE > 2) WRITE(*,'(A,I0)', ADVANCE='NO') "**", IE-1
        END IF
      END DO 
      WRITE(*,*)
      CALL WRITE_OK_OR_FAIL(CONST_RESULT(IEQ,1), AEXP(ICH, IEQ), "The c constant for channel "//TRIM(CHANNEL_NAME))
      CALL WRITE_OK_OR_FAIL(CONST_RESULT(IEQ,2), BEXP(ICH, IEQ), "The b constant for channel "//TRIM(CHANNEL_NAME))
    ENDDO ! IEQ
    WRITE(*,*)
    CLOSE(UNIT=10)
  ENDDO ! ICH



  ! EFT_pless 10 TEST
  CALL CREATE_DIRECTORY(TRIM(EFT_PLESS_10_KCOTD_DIR))
  WRITE(*,*) "========================================================"
  WRITE(*,*) "========================================================"
  PRINT *, "   Testing EFT_pless ilb=10 potential..."
  WRITE(*,*) "========================================================"
  WRITE(*,*) "========================================================"

  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  NCHANNELS = SIZE(CHANNELS)

  CALL REALLOCATE(AEXP, NCHANNELS, 2)
  CALL REALLOCATE(BEXP, NCHANNELS, 2)
  AEXP(1,1) =  4.358D-2  ! 1S0
  AEXP(2,1) = -1.849D-1  ! 3S1
  AEXP(3,1) = -5.830D-1  ! 1P1
  AEXP(4,1) =  4.665D-1  ! 3P0
  AEXP(5,1) = -7.495D-1  ! 3P1
  AEXP(6,1) =  2.827D+0  ! 3P2
  AEXP(7,1) =  2.664D+0  ! 1D2
  AEXP(2,2) = -6.348D-1  ! 3D1
  AEXP(8,1) =  1.124D+0  ! 3D2
  AEXP(6,2) =  1.628D+1  ! 3F2

  BEXP(1,1) =  1.381D+0  ! 1S0
  BEXP(2,1) =  8.656D-1  ! 3S1
  BEXP(3,1) = -1.949D+0  ! 1P1
  BEXP(4,1) =  7.174D-1  ! 3P0
  BEXP(5,1) = -2.270D+0  ! 3P1
  BEXP(6,1) =  3.366D+0  ! 3P2
  BEXP(7,1) =  4.319D+0  ! 1D2
  BEXP(2,2) = -5.554D+0  ! 3D1
  BEXP(8,1) =  1.469D+0  ! 3D2
  BEXP(6,2) =  8.903D+0  ! 3F2

  DO ICH = 1, NCHANNELS
    WRITE(*,*) "========================================================"
    CHANNEL_NAME = TRIM(GET_CHANNEL_NAME(CHANNELS(ICH)))
    PRINT *, "  Channel: ", CHANNEL_NAME
    OPEN(UNIT=10, FILE=TRIM(EFT_PLESS_10_DIR)//"delta_"//TRIM(CHANNEL_NAME)//".dat", STATUS='OLD', ACTION='READ')

    ! First, count the number of data lines (not starting with #)
    NDATA = 0
    DO
      READ(10, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      IF (LEN_TRIM(LINE) == 0) CYCLE
      LINE = ADJUSTL(LINE)
      IF (LINE(1:1) == "#") CYCLE
      NDATA = NDATA + 1
    END DO
    REWIND(10)

    IF (ALLOCATED(ENERGIES)) THEN
      DEALLOCATE(ENERGIES)
    END IF
    ALLOCATE(ENERGIES(NDATA))
    IF (ALLOCATED(PHASE_SHIFT_ARRAY)) THEN
      DEALLOCATE(PHASE_SHIFT_ARRAY)
    END IF
    ALLOCATE(PHASE_SHIFT_ARRAY(NDATA))

    WRITE(*,*) "  Number of data points: ", NDATA
    IE = 0
    DO
      READ(10, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      IF (LEN_TRIM(LINE) == 0) CYCLE
      LINE = ADJUSTL(LINE)
      IF (LINE(1:1) == "#") CYCLE
      IE = IE + 1
      READ(LINE, *) ENERGIES(IE), PHASE_SHIFT_ARRAY(IE)%delta1_S, &
            PHASE_SHIFT_ARRAY(IE)%delta2_S, PHASE_SHIFT_ARRAY(IE)%epsilon_S
    END DO
    CLOSE(UNIT=10)


    IF (ALLOCATED(K2)) THEN
      DEALLOCATE(K2)
    END IF
    ALLOCATE(K2(NDATA))
    IF (ALLOCATED(KCOTD)) THEN
      DEALLOCATE(KCOTD)
    END IF
    ALLOCATE(KCOTD(2,NDATA))

    OPEN(UNIT=10, FILE=TRIM(EFT_PLESS_10_KCOTD_DIR)//"k2_"//TRIM(CHANNEL_NAME)//".dat", STATUS='UNKNOWN', ACTION='WRITE')

    CONST_RESULT = 0.0D0
    FITTED = FIT_CHANNEL_LOW_ENERGY(CHANNELS(ICH), ENERGIES, PHASE_SHIFT_ARRAY, CONST_RESULT, FIT_ORDER, K2, KCOTD)

    IF (.NOT. FITTED) THEN
      WRITE(*,*) "  WARNING: Fitting failed for channel ", CHANNEL_NAME
      WRITE(*,*) "  This may be due to insufficient data points or a poor fit order."
    END IF

    DO IE = 1, NDATA
      WRITE(10, *) K2(IE), KCOTD(1,IE), KCOTD(2,IE)
    END DO


    WRITE(*,*) "  Fitted constants for channel ", CHANNEL_NAME, ":"
    DO IEQ = 1, GET_CHANNEL_NCH(CHANNELS(ICH))
      IF (IEQ==2) WRITE(*,*)
      WRITE(*,'(A)', ADVANCE='NO') "    Fit: "
      DO IE = 1, FIT_ORDER+1
        IF (IE > 1) WRITE(*,'(A)', ADVANCE='NO') " + "
        IF (IE == 1) THEN
          WRITE(*,'(F12.6)', ADVANCE='NO') CONST_RESULT(IEQ,IE)
        ELSE
          WRITE(*,'(F12.6,A,I0)', ADVANCE='NO') CONST_RESULT(IEQ, IE), " * X"
          IF (IE > 2) WRITE(*,'(A,I0)', ADVANCE='NO') "**", IE-1
        END IF
      END DO 
      WRITE(*,*)
      CALL WRITE_OK_OR_FAIL(CONST_RESULT(IEQ,1), AEXP(ICH, IEQ), "The c constant for channel "//TRIM(CHANNEL_NAME))
      CALL WRITE_OK_OR_FAIL(CONST_RESULT(IEQ,2), BEXP(ICH, IEQ), "The b constant for channel "//TRIM(CHANNEL_NAME))
    ENDDO ! IEQ
    WRITE(*,*)
    CLOSE(UNIT=10)
  ENDDO ! ICH

CONTAINS

  SUBROUTINE WRITE_OK_OR_FAIL(RESULT, EXPECTED, LABEL)
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    DOUBLE PRECISION, INTENT(IN) :: RESULT, EXPECTED
    IF (ABS_DIFF_PROCENT(RESULT, EXPECTED) < 5.D0) THEN
      WRITE(*,'(A,A, F9.6, A, F9.6, A, F8.3,A)') CHAR(27)//"[32m[OK]"//CHAR(27)//"[0m "//LABEL, &
      " is within 5% of the expected value: ", &
            RESULT, "  in accord to ", EXPECTED, "   (", ABS_DIFF_PROCENT(RESULT, EXPECTED), " %)"
    ELSE
      WRITE(*,'(A,A, F9.6, A, F9.6, A, F8.3,A)') CHAR(27)//"[31m[FAIL]"//CHAR(27)//"[0m "//LABEL, &
            " is NOT within 5% of the expected value: ", &
            RESULT, "  and not ", EXPECTED, "(", ABS_DIFF_PROCENT(RESULT, EXPECTED), " %)"
    END IF
  END SUBROUTINE WRITE_OK_OR_FAIL


END PROGRAM FIT_LOW_ENERGY_OBSERVABLES