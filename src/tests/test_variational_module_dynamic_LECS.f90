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

  TYPE(RESULT_FROM_FILE) :: RESULT
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
  INTEGER :: J, L, S
  LOGICAL :: PRINT_COEFFICIENTS = .FALSE., ENERGIES_SET = .FALSE.
  
  TYPE(LECS_EFT_PLESS), EXTERNAL :: READ_LECS_EFT_PLESS
  LOGICAL, EXTERNAL :: FILE_EXISTS

  !---------------------------------------------------------------------------
  ! Prepare the list of channels based on the specified maximum angular momentum
  ! (LMAX), maximum total angular momentum (JMAX), and isospin projection (TZ).
  ! The resulting CHANNELS array is populated accordingly.
  ! The number of channels (NCHANNELS) is determined by the size of the CHANNELS array.
  ! Allocate the FILE_NAMES array to store file names for each channel.
  ! Print the total number of channels to standard output for verification.
  !---------------------------------------------------------------------------
  LMAX = 0
  JMAX = 0
  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  NCHANNELS = SIZE(CHANNELS)
  ALLOCATE(FILE_NAMES(NCHANNELS))
  PRINT *, 'Number of channels:', NCHANNELS 

    !========================================================================================


  
  ! Loop over all scattering channels, construct file names, print info, and check file existence
  DO I=1, NCHANNELS
    FILE_NAMES(I) = "output/EFT_pless_15/delta_"//TRIM(GET_CHANNEL_NAME(CHANNELS(I))) // ".dat"
    WRITE(*,*) 'Channel ', I, ':', TRIM(GET_CHANNEL_NAME(CHANNELS(I)))
    IF (FILE_EXISTS(TRIM(FILE_NAMES(I)))) THEN
      WRITE(*,*) '            > File exists : ', TRIM(FILE_NAMES(I))
    ELSE
      WRITE(*,*) '            > File does not exist: ', TRIM(FILE_NAMES(I))
    END IF
  END DO

  PRINT *, REPEAT("=", 80)
  DO ICH = 1, NCHANNELS
    WRITE(*,*) 'Processing channel: ', TRIM(GET_CHANNEL_NAME(CHANNELS(ICH)))
    CALL READ_RESULTS_FROM_FILE(TRIM(FILE_NAMES(ICH)), RESULTS)
    IF (.NOT.ENERGIES_SET) THEN
      NE = SIZE(RESULTS)
      CALL REALLOCATE_1D_1(ENERGIES, NE)
      DO IE = 1, NE
        ENERGIES(IE) = RESULTS(IE)%E
      END DO
      CALL SET_ENERGIES(ENERGIES)
      CALL SET_CHANNELS(CHANNELS)
      ENERGIES_SET = .TRUE.
    END IF
    CALL PRINT_SCATTERING_CHANNEL(CHANNELS(ICH))

    DO I = 1, NE
      IF (I /= 1) CYCLE
      E =  ENERGIES(I)
      J = GET_CHANNEL_J(CHANNELS(ICH))
      L = GET_CHANNEL_L(CHANNELS(ICH),1)
      S = GET_CHANNEL_S(CHANNELS(ICH),1)
      TZ = GET_CHANNEL_TZ(CHANNELS(ICH))
      CALL NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, 19, 15, LEMP, PHASE_SHIFTS, PRINT_COEFFICIENTS=.FALSE.)
      WRITE(20, *) E, PHASE_SHIFTS%delta1_S, PHASE_SHIFTS%delta2_S, PHASE_SHIFTS%epsilon_S
    ENDDO
  ENDDO

    call RESET_SCATTERING_NN_VARIATIONAL
    ENERGIES_SET = .FALSE.

    !========================================================================================

  LEC_15 = READ_LECS_EFT_PLESS(15)
  CALL PRINT_LECS(LEC_15)

  CALL SET_DYNAMIC(.TRUE.)
  
  PRINT *, REPEAT("=", 80)
  DO ICH = 1, NCHANNELS
    WRITE(*,*) 'Processing channel: ', TRIM(GET_CHANNEL_NAME(CHANNELS(ICH)))
    CALL READ_RESULTS_FROM_FILE(TRIM(FILE_NAMES(ICH)), RESULTS)
    IF (.NOT.ENERGIES_SET) THEN
      NE = SIZE(RESULTS)
      CALL REALLOCATE_1D_1(ENERGIES, NE)
      DO IE = 1, NE
        ENERGIES(IE) = RESULTS(IE)%E
      END DO
      CALL SET_ENERGIES(ENERGIES)
      ENERGIES_SET = .TRUE.
      CALL SET_NEW_LECS(LEC_15)
      CALL SET_CHANNELS(CHANNELS)
    END IF
    ! DO I = 1, NE
    !   WRITE(*, '(F10.4, 3F10.4)') RESULTS(I)%E, RESULTS(I)%PHASES(1), RESULTS(I)%PHASES(2), RESULTS(I)%PHASES(3)
    ! END DO
    CALL PRINT_SCATTERING_CHANNEL(CHANNELS(ICH))

    DO I = 1, NE
      IF (I /= 1) CYCLE
      E =  ENERGIES(I)
      J = GET_CHANNEL_J(CHANNELS(ICH))
      L = GET_CHANNEL_L(CHANNELS(ICH),1)
      S = GET_CHANNEL_S(CHANNELS(ICH),1)
      TZ = GET_CHANNEL_TZ(CHANNELS(ICH))
      CALL NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP, PHASE_SHIFTS, PRINT_COEFFICIENTS=.FALSE.)
      WRITE(21, *) E, PHASE_SHIFTS%delta1_S, PHASE_SHIFTS%delta2_S, PHASE_SHIFTS%epsilon_S
    ENDDO
  ENDDO


  

CONTAINS
  SUBROUTINE READ_RESULTS_FROM_FILE(FILENAME, RES)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    TYPE(RESULT_FROM_FILE), ALLOCATABLE, INTENT(OUT) :: RES(:)
    CHARACTER(LEN=512) :: LINE
    INTEGER :: UNIT, IOS, COUNT, I
    DOUBLE PRECISION :: E, DELTA1, DELTA2, MIXING_ANGLE
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


    I = 0
    DO
      READ(UNIT, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      LINE = ADJUSTL(LINE)
      IF (LEN_TRIM(LINE) == 0) CYCLE
      FIRST_CHAR = LINE(1:1)
      IF (FIRST_CHAR == '#') CYCLE
      READ(LINE, *, IOSTAT=IOS) E, DELTA1, DELTA2, MIXING_ANGLE
      IF (IOS /= 0) CYCLE
      I = I + 1
      RES(I)%E = E
      RES(I)%PHASES(1) = DELTA1
      RES(I)%PHASES(2) = DELTA2
      RES(I)%PHASES(3) = MIXING_ANGLE
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
  DOUBLE PRECISION :: TMP(40) ! Large enough for all columns
  INTEGER :: I

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
    READ(UNIT, *, IOSTAT=IOS) (TMP(I), I=1, 30)
    IF (IOS /= 0) EXIT
    IF (INT(TMP(1)) == ID) THEN
      LECS%ILB = INT(TMP(1))
      LECS%RC(0,0) = TMP(2)
      LECS%RC(1,0) = TMP(3)
      LECS%RC(0,1) = TMP(4)
      LECS%RC(1,1) = TMP(5)
      LECS%CLO(1)  = TMP(6)
      LECS%CLO(0)  = TMP(7)
      LECS%CNLO(1) = TMP(8)
      LECS%CNLO(2) = TMP(9)
      LECS%CNLO(3) = TMP(10)
      LECS%CNLO(4) = TMP(11)
      LECS%CNLO(5) = TMP(12)
      LECS%CNLO(6) = TMP(13)
      LECS%CNLO(7) = TMP(14)
      LECS%CN3LO(1) = TMP(15)
      LECS%CN3LO(2) = TMP(16)
      LECS%CN3LO(3) = TMP(17)
      LECS%CN3LO(4) = TMP(18)
      LECS%CN3LO(5) = TMP(19)
      LECS%CN3LO(6) = TMP(20)
      LECS%CN3LO(7) = TMP(21)
      LECS%CN3LO(8) = TMP(22)
      LECS%CN3LO(9) = TMP(23)
      LECS%CN3LO(10) = TMP(24)
      LECS%CN3LO(11) = TMP(25)
      LECS%CIT(0) = TMP(26)
      LECS%CIT(1) = TMP(27)
      LECS%CIT(2) = TMP(28)
      LECS%CIT(3) = TMP(29)
      LECS%CIT(4) = TMP(30)
      EXIT
    END IF
  END DO
  CLOSE(UNIT)
  IF (ID <= 5) LECS%ORDER = 0
  IF (ID > 5 .AND. ID <= 10) LECS%ORDER = 1
  IF (ID > 10 .AND. ID <= 15) LECS%ORDER = 3
END FUNCTION READ_LECS_EFT_PLESS