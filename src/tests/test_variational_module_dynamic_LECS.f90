PROGRAM VARIATIONAL_WITH_DYNAMIC_LECS
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS
  USE EFT_PLESS
  USE REALLOCATE_UTILS
  IMPLICIT NONE
  TYPE RESULT_FROM_FILE
    DOUBLE PRECISION :: E
    DOUBLE PRECISION :: PHASES(3)
  END TYPE RESULT_FROM_FILE

  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
  TYPE(CHARACTER(LEN=256)), ALLOCATABLE :: FILE_NAMES(:)
  TYPE(LECS_EFT_PLESS) :: LEC_15
  INTEGER :: NCHANNELS, LMAX = 2, JMAX = 2, TZ = 0
  INTEGER :: I
  TYPE(RESULT_FROM_FILE), ALLOCATABLE :: RESULTS(:)
  DOUBLE PRECISION, ALLOCATABLE :: ENERGIES(:)
  INTEGER :: ICH, IE, NE, IPOT = -1, ILB = -1, LEMP = 0
  DOUBLE PRECISION :: E
  TYPE(PHASE_SHIFT_RESULT) :: PHASE_SHIFTS
  INTEGER :: J, L, S, T
  LOGICAL :: ENERGIES_SET = .FALSE.


  CHARACTER(LEN=*), PARAMETER :: RED = CHAR(27)//'[31m'
  CHARACTER(LEN=*), PARAMETER :: RESET = CHAR(27)//'[0m'
  CHARACTER(LEN=*), PARAMETER :: GREEN = CHAR(27)//'[32m'
  
  TYPE(LECS_EFT_PLESS), EXTERNAL :: READ_LECS_EFT_PLESS
  LOGICAL, EXTERNAL :: FILE_EXISTS

  INTEGER, ALLOCATABLE :: NLINES(:)

  !---------------------------------------------------------------------------
  ! Prepare the list of channels based on the specified maximum angular momentum
  ! (LMAX), maximum total angular momentum (JMAX), and isospin projection (TZ).
  ! The resulting CHANNELS array is populated accordingly.
  ! The number of channels (NCHANNELS) is determined by the size of the CHANNELS array.
  ! Allocate the FILE_NAMES array to store file names for each channel.
  ! Print the total number of channels to standard output for verification.
  !---------------------------------------------------------------------------
  IPOT = 19
  ILB = 15
  
  ! LMAX = 1  !if not 0, doesn't work
  LMAX = 2
  JMAX = 2
  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  NCHANNELS = SIZE(CHANNELS)
  ALLOCATE(FILE_NAMES(NCHANNELS))
  ALLOCATE(NLINES(NCHANNELS))
  PRINT *, 'Number of channels:', NCHANNELS 
  
  !========================================================================================
  
  
  
  
  ! Loop over all scattering channels, construct file names, print info, and check file existence
  DO I=1, NCHANNELS
    FILE_NAMES(I) = "bin/tests/test_files/scattering_phase_shifts/EFT_pless_15/delta_"&
                        //TRIM(GET_CHANNEL_NAME(CHANNELS(I))) // ".dat"
    WRITE(*,*) 'Channel ', I, ':', TRIM(GET_CHANNEL_NAME(CHANNELS(I)))
    IF (FILE_EXISTS(TRIM(FILE_NAMES(I)))) THEN
      WRITE(*,*) '            > File exists : ', TRIM(FILE_NAMES(I))
    ELSE
      WRITE(*,*) '            > File does not exist: ', TRIM(FILE_NAMES(I))
    END IF
  END DO

  LEC_15 = READ_LECS_EFT_PLESS(ILB)
  CALL PRINT_LECS(LEC_15)

  CALL SET_DYNAMIC(.TRUE.)
  
  PRINT *, REPEAT("=", 80)
  DO ICH = 1, NCHANNELS
    WRITE(*,*) 'Processing channel: ', TRIM(GET_CHANNEL_NAME(CHANNELS(ICH)))
    CALL READ_RESULTS_FROM_FILE(TRIM(FILE_NAMES(ICH)), RESULTS)
    IF (.NOT.ENERGIES_SET) THEN
      NE = SIZE(RESULTS)
      CALL REALLOCATE(ENERGIES, NE)
      DO IE = 1, NE
        ENERGIES(IE) = RESULTS(IE)%E
      END DO
      CALL SET_ENERGIES(ENERGIES)
      ENERGIES_SET = .TRUE.
      CALL SET_NEW_LECS(LEC_15)
      CALL SET_CHANNELS(CHANNELS)
    END IF

    S = GET_CHANNEL_S(CHANNELS(ICH),1)
    T = GET_CHANNEL_T(CHANNELS(ICH),1)
    IF ( (S==0 .AND. T==0 .OR. S==1 .AND. T==1) .AND. ILB < 6 ) CYCLE
    ! DO I = 1, NE
    !   WRITE(*, '(F10.4, 3F10.4)') RESULTS(I)%E, RESULTS(I)%PHASES(1), RESULTS(I)%PHASES(2), RESULTS(I)%PHASES(3)
    ! END DO
    CALL PRINT_SCATTERING_CHANNEL(CHANNELS(ICH))

    DO I = 1, NE
      E =  ENERGIES(I)
      J = GET_CHANNEL_J(CHANNELS(ICH))
      L = GET_CHANNEL_L(CHANNELS(ICH),1)
      S = GET_CHANNEL_S(CHANNELS(ICH),1)
      TZ = GET_CHANNEL_TZ(CHANNELS(ICH))
      CALL NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, -1, -1, LEMP, PHASE_SHIFTS, PRINT_COEFFICIENTS=.FALSE.)
      ! Print energy and absolute percentage differences for all three on the same line
      IF (I <= SIZE(RESULTS)) THEN
        IF (print_diff_pct_line(E,  PHASE_SHIFTS%delta1_S , RESULTS(I)%PHASES(1), &
                                    PHASE_SHIFTS%delta2_S , RESULTS(I)%PHASES(2), &
                                    PHASE_SHIFTS%epsilon_S, RESULTS(I)%PHASES(3))) THEN
          NLINES(ICH) = NLINES(ICH) + 1
        ENDIF
      END IF
    ENDDO
  ENDDO

  DO ICH = 1, NCHANNELS
    IF (NLINES(ICH) == 0) THEN
      WRITE(*,*) GREEN, 'Channel ', GET_CHANNEL_NAME_FROM_OBJECT(CHANNELS(ICH)), ': ', "PASSED.", RESET
    ELSE
      WRITE(*,*) RED, 'Channel ', GET_CHANNEL_NAME_FROM_OBJECT(CHANNELS(ICH)), ': ', "FAILED. ", &
          RESET, TRIM(GET_CHANNEL_NAME(CHANNELS(ICH))), ' has significant differences.', RESET
    ENDIF
  ENDDO

  DEALLOCATE(NLINES)
  DEALLOCATE(ENERGIES)
  DO I = 1, NCHANNELS
    CALL RESET_CHANNEL(CHANNELS(I))
  ENDDO
  DEALLOCATE(CHANNELS)
  DEALLOCATE(FILE_NAMES)
  DEALLOCATE(RESULTS)
  CALL RESET_SCATTERING_NN_VARIATIONAL


CONTAINS
  FUNCTION print_diff_pct_line(ENERGY, val1, ref1, val2, ref2, val3, ref3) RESULT(DIFFERENT_LINES)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: ENERGY, val1, ref1, val2, ref2, val3, ref3
    DOUBLE PRECISION :: pct1, pct2, pct3
    DOUBLE PRECISION, PARAMETER :: SMALL = 1.D-5
    DOUBLE PRECISION, PARAMETER :: THRESH = 1.D-2  ! Threshold for percentage difference
    DOUBLE PRECISION, PARAMETER :: THRESH2= 1.D-1  ! Threshold for percentage difference
    LOGICAL :: DIFFERENT_LINES

    IF (ABS(val1) < SMALL .AND. ABS(ref1) < SMALL) THEN
      pct1 = 0.0D0
    ELSE
      pct1 = ABS(val1 - ref1) / MAX(1D-12, ABS(ref1)) * 100.0D0
    END IF
    IF (ABS(val2) < SMALL .AND. ABS(ref2) < SMALL) THEN
      pct2 = 0.0D0
    ELSE
      pct2 = ABS(val2 - ref2) / MAX(1D-12, ABS(ref2)) * 100.0D0
    END IF
    IF (ABS(val3) < SMALL .AND. ABS(ref3) < SMALL) THEN
      pct3 = 0.0D0
    ELSE
      pct3 = ABS(val3 - ref3) / MAX(1D-12, ABS(ref3)) * 100.0D0
    END IF

    DIFFERENT_LINES = .FALSE.
    IF (pct1 > THRESH .OR. pct2 > THRESH2 .OR. pct3 > THRESH) THEN
      DIFFERENT_LINES = .TRUE.
      IF (pct1 > THRESH) THEN
        WRITE(*,'(F10.3,1X,A,F20.12,A,1X)', ADVANCE='NO') ENERGY, RED, pct1, RESET
      ELSE
        WRITE(*,'(F10.3,1X,F20.12,1X)', ADVANCE='NO') ENERGY, pct1
      END IF
      IF (pct2 > THRESH) THEN
        WRITE(*,'(A,F20.12,A,1X)', ADVANCE='NO') RED, pct2, RESET
      ELSE
        WRITE(*,'(F20.12,1X)', ADVANCE='NO') pct2
      END IF
      IF (pct3 > THRESH) THEN
        WRITE(*,'(A,F20.12,A)') RED, pct3, RESET
      ELSE
        WRITE(*,'(F20.12)') pct3
      END IF
    END IF
  END FUNCTION print_diff_pct_line

  SUBROUTINE READ_RESULTS_FROM_FILE(FILENAME, RES)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    TYPE(RESULT_FROM_FILE), ALLOCATABLE, INTENT(OUT) :: RES(:)
    CHARACTER(LEN=512) :: LINE
    INTEGER :: UNIT, IOS, COUNT, II
    DOUBLE PRECISION :: ENERGY, DELTA1, DELTA2, MIXING_ANGLE
    CHARACTER :: FIRST_CHAR

    ! First, count the number of valid data lines
    COUNT = 0
    OPEN(NEWUNIT=UNIT, FILE=TRIM(FILENAME), STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (IOS /= 0) THEN
      ALLOCATE(RES(0))
      RETURN
    END IF
    DO
      READ(UNIT, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      LINE = ADJUSTL(LINE)
      IF (LEN_TRIM(LINE) == 0) CYCLE
      FIRST_CHAR = LINE(1:1)
      IF (FIRST_CHAR == '#') CYCLE
      COUNT = COUNT + 1
    END DO
    REWIND(UNIT)
    ALLOCATE(RES(COUNT))

    II = 0
    DO
      READ(UNIT, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      LINE = ADJUSTL(LINE)
      IF (LEN_TRIM(LINE) == 0) CYCLE
      FIRST_CHAR = LINE(1:1)
      IF (FIRST_CHAR == '#') CYCLE
      READ(LINE, *, IOSTAT=IOS) ENERGY, DELTA1, DELTA2, MIXING_ANGLE
      IF (IOS /= 0) CYCLE
      II = II + 1
      RES(II)%E = ENERGY
      RES(II)%PHASES(1) = DELTA1
      RES(II)%PHASES(2) = DELTA2
      RES(II)%PHASES(3) = MIXING_ANGLE
      ! J, L, S, TZ are left uninitialized (set as needed)
    END DO
    CLOSE(UNIT)
  END SUBROUTINE READ_RESULTS_FROM_FILE
END PROGRAM VARIATIONAL_WITH_DYNAMIC_LECS
















LOGICAL FUNCTION FILE_EXISTS(FILENAME)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: FILENAME
  INQUIRE(FILE=FILENAME, EXIST=FILE_EXISTS)
END FUNCTION FILE_EXISTS

! Reads the LECS from src/libs/physics/potentials/lecs_eft.dat for the row with first column = 15
FUNCTION READ_LECS_EFT_PLESS(ID) RESULT(LECS)
  USE EFT_PLESS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ID
  TYPE(LECS_EFT_PLESS) :: LECS
  CHARACTER(LEN=256) :: FILENAME
  INTEGER :: UNIT, IOS

  FILENAME = 'src/libs/physics/potentials/lecs_eft.dat'
  LECS = LECS_EFT_PLESS() ! Initialize to default

  OPEN(NEWUNIT=UNIT, FILE=TRIM(FILENAME), STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS /= 0) THEN
    PRINT *, 'Error opening file: ', TRIM(FILENAME)
    RETURN
  END IF

  ! Skip the first line (number of models)
  READ(UNIT, *, IOSTAT=IOS)

  DO
    READ(UNIT, *, IOSTAT=IOS) LECS%ILB, LECS%RC(0,0), LECS%RC(1,0), LECS%RC(0,1), LECS%RC(1,1), LECS%CLO(1), LECS%CLO(0),&
                 LECS%CNLO, LECS%CN3LO, LECS%CIT
    IF (IOS /= 0) EXIT
    IF (INT(LECS%ILB) == ID) EXIT
  END DO
  CLOSE(UNIT)
  IF (ID <= 5) LECS%ORDER = 0
  IF (ID > 5 .AND. ID <= 10) LECS%ORDER = 1
  IF (ID > 10 .AND. ID <= 15) LECS%ORDER = 3
END FUNCTION READ_LECS_EFT_PLESS