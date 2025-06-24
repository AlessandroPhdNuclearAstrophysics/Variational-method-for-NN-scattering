!> \file quantum_numbers.f90
!! \brief Quantum number utilities and SCATTERING_CHANNEL type for nuclear/particle physics.
!!
!! This module defines the SCATTERING_CHANNEL type and provides helper routines
!! for initializing, setting, and querying quantum numbers (L, S, J, T, TZ) for
!! scattering channels. It also includes utilities for naming channels and checking
!! physical validity and coupling.
!!
!! \author Alessandro
!! \date 2025
MODULE QUANTUM_NUMBERS
  USE REALLOCATE_UTILS
  IMPLICIT NONE

  TYPE, PUBLIC :: SCATTERING_CHANNEL
    INTEGER, PRIVATE :: J
    INTEGER, ALLOCATABLE, PRIVATE :: L(:)
    INTEGER, ALLOCATABLE, PRIVATE :: S(:)
    INTEGER, ALLOCATABLE, PRIVATE :: T(:)
    INTEGER, PRIVATE :: TZ
    INTEGER, PRIVATE :: NCH
    LOGICAL, PRIVATE :: COUPLED = .FALSE.
  ENDTYPE SCATTERING_CHANNEL

  INTERFACE GET_CHANNEL_NAME
    MODULE PROCEDURE GET_CHANNEL_NAME_LSJ
    MODULE PROCEDURE GET_CHANNEL_NAME_FROM_OBJECT
  END INTERFACE GET_CHANNEL_NAME

  PUBLIC :: init_scattering_channel
  PUBLIC :: SET_CHANNEL
  PUBLIC :: GET_CHANNEL_NAME
  PUBLIC :: GET_CHANNEL_NCH
  PUBLIC :: GET_CHANNEL_L
  PUBLIC :: GET_CHANNEL_S
  PUBLIC :: GET_CHANNEL_T
  PUBLIC :: GET_CHANNEL_J
  PUBLIC :: GET_CHANNEL_TZ
  PUBLIC :: IS_CHANNEL_COUPLED
  PUBLIC :: IS_PHYSICAL_CHANNEL
  PUBLIC :: GET_CHANNEL_FROM_NAME
  PUBLIC :: PREPARE_CHANNELS
  PUBLIC :: PRINT_SCATTERING_CHANNEL
  PUBLIC :: RESET_CHANNEL
  PUBLIC :: TZ_TO_T1Z_T2Z
  PUBLIC :: T1Z_T2Z_TO_TZ
  PUBLIC :: T_FROM_L_S
  
CONTAINS
  !> \brief Constructor for SCATTERING_CHANNEL.
  !! \param[in] J Total angular momentum
  !! \param[in] IS_EVEN Logical for parity
  !! \param[in] TZ Isospin projection
  !! \return Initialized SCATTERING_CHANNEL object
  FUNCTION init_scattering_channel(J, IS_EVEN, TZ) RESULT(channel)
    INTEGER, INTENT(IN) :: J, TZ
    LOGICAL, INTENT(IN) :: IS_EVEN
    TYPE(SCATTERING_CHANNEL) :: CHANNEL

    INTEGER :: ICH

    IF (J == 0) THEN
      CHANNEL%NCH = 1
    ELSE
      CHANNEL%NCH = 2
    ENDIF

    ! ALLOCATE ARRAYS FOR L AND S
    CALL REALLOCATE_1D_1_INT(CHANNEL%L, CHANNEL%NCH)
    CALL REALLOCATE_1D_1_INT(CHANNEL%S, CHANNEL%NCH)
    CALL REALLOCATE_1D_1_INT(CHANNEL%T, CHANNEL%NCH)
    IF (J == 0) THEN
      IF (IS_EVEN) THEN
        CHANNEL%L(1) = 0
        CHANNEL%S(1) = 0
      ELSE
        CHANNEL%L(1) = 1
        CHANNEL%S(1) = 1
      ENDIF
    ELSE
      IF (MOD(J, 2) == 0) THEN
        IF (IS_EVEN) THEN
          CHANNEL%L(1) = J
          CHANNEL%S(1) = 0
          CHANNEL%L(2) = J
          CHANNEL%S(2) = 1
        ELSE
          CHANNEL%L(1) = J-1
          CHANNEL%S(1) = 1
          CHANNEL%L(2) = J+1
          CHANNEL%S(2) = 1
          CHANNEL%COUPLED = .TRUE.
        ENDIF
      ELSE
        IF (IS_EVEN) THEN
          CHANNEL%L(1) = J - 1
          CHANNEL%S(1) = 1
          CHANNEL%L(2) = J + 1
          CHANNEL%S(2) = 1
          CHANNEL%COUPLED = .TRUE.
        ELSE
          CHANNEL%L(1) = J
          CHANNEL%S(1) = 0
          CHANNEL%L(2) = J
          CHANNEL%S(2) = 1
        ENDIF
      ENDIF
    ENDIF

    CHANNEL%J = J
    DO ICH = 1, CHANNEL%NCH
      CHANNEL%T(ICH) = EVALUATE_T(CHANNEL%L(ICH), CHANNEL%S(ICH))
    ENDDO
    CHANNEL%TZ = TZ
  END FUNCTION init_scattering_channel

  !> \brief Evaluate isospin T for given L, S, and TZ.
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !!
  !! \return T The isospin quantum number
  FUNCTION EVALUATE_T(L, S) RESULT(T)
    INTEGER, INTENT(IN) :: L, S
    INTEGER :: T

    ! Calculate T based on the values of J, L, S, and TZ
    T = MOD( MOD(L+S, 2) + 1 , 2)
  END FUNCTION EVALUATE_T

  !> \brief Check if a channel is physical.
  !! \param[in] CHANNEL The SCATTERING_CHANNEL object
  !! \return .TRUE. if physical, .FALSE. otherwise
  FUNCTION IS_PHYSICAL_CHANNEL(CHANNEL) RESULT(IS_PHYSICAL)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    LOGICAL :: IS_PHYSICAL
    INTEGER :: L, S, J, T, ICH, TZ

    J = CHANNEL%J
    IS_PHYSICAL = .TRUE.

    DO ICH=1, CHANNEL%NCH
      L = CHANNEL%L(ICH)
      S = CHANNEL%S(ICH)
      T = CHANNEL%T(ICH)
      TZ = CHANNEL%TZ
      IF (ABS(TZ) > T) IS_PHYSICAL = .FALSE.

      ! CHECK IF THE SCATTERING CHANNEL IS PHYSICAL
      IS_PHYSICAL = IS_PHYSICAL .AND. (J >= 0 .AND. L >= 0 .AND. S >= 0 .AND. T >= 0)
      IF (IS_PHYSICAL) THEN
        IS_PHYSICAL = IS_PHYSICAL .AND. IS_LSJ_PHYSICAL(L, S, J)
      ENDIF
    ENDDO
  END FUNCTION IS_PHYSICAL_CHANNEL

  !> \brief Check if LSJ quantum numbers are physical.
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] J Total angular momentum
  !! \return .TRUE. if physical, .FALSE. otherwise
  FUNCTION IS_LSJ_PHYSICAL(L, S, J) RESULT(IS_LSJ)
    INTEGER, INTENT(IN) :: L, S, J
    LOGICAL :: IS_LSJ
    IS_LSJ = ABS(L-S) <= J .AND. J <= (L+S)
  END FUNCTION IS_LSJ_PHYSICAL

  !> \brief Set the quantum numbers of a SCATTERING_CHANNEL.
  !! \brief Set the quantum numbers for a scattering channel.
  !! Sets the J, L, S, TZ quantum numbers and determines if the channel is coupled or not.
  !! Allocates and fills the L, S, T arrays for the channel.
  !! \param[inout] CHANNEL Scattering channel object to set
  !! \param[in] J Total angular momentum
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] TZ Isospin projection
  SUBROUTINE SET_CHANNEL(CHANNEL, J, L, S, TZ)
    TYPE(SCATTERING_CHANNEL), INTENT(INOUT) :: CHANNEL
    INTEGER, INTENT(IN) :: J, L, S, TZ
    LOGICAL :: IS_EVEN

    IS_EVEN = MOD(L, 2) == 0

    CHANNEL%J = J
    CHANNEL%TZ = TZ
    CHANNEL%COUPLED = .FALSE.

    IF (J == 0) THEN
      CHANNEL%NCH = 1
      CALL REALLOCATE_1D_1_INT(CHANNEL%L, 1)
      CALL REALLOCATE_1D_1_INT(CHANNEL%S, 1)
      CALL REALLOCATE_1D_1_INT(CHANNEL%T, 1)
      CHANNEL%L(1) = L
      CHANNEL%S(1) = S
      CHANNEL%T(1) = EVALUATE_T(L, S)
    ELSE
      IF ((L==(J-1) .OR. L==(J+1))) THEN
        CHANNEL%NCH = 2
        CALL REALLOCATE_1D_1_INT(CHANNEL%L, 2)
        CALL REALLOCATE_1D_1_INT(CHANNEL%S, 2)
        CALL REALLOCATE_1D_1_INT(CHANNEL%T, 2)
        CHANNEL%L(1) = J-1
        CHANNEL%S(1) = S
        CHANNEL%L(2) = J+1
        CHANNEL%S(2) = S
        CHANNEL%T(1) = EVALUATE_T(J-1, S)
        CHANNEL%T(2) = EVALUATE_T(J+1, S)
        CHANNEL%COUPLED = .TRUE.
      ELSE
        CHANNEL%NCH = 1
        CALL REALLOCATE_1D_1_INT(CHANNEL%L, 1)
        CALL REALLOCATE_1D_1_INT(CHANNEL%S, 1)
        CALL REALLOCATE_1D_1_INT(CHANNEL%T, 1)
        CHANNEL%L(1) = J
        CHANNEL%S(1) = S
        CHANNEL%T(1) = EVALUATE_T(J, S)
      ENDIF
    ENDIF
  END SUBROUTINE SET_CHANNEL

  !> \brief Get the spectroscopic name for a channel from L, S, J.
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] J Total angular momentum
  !! \return Name as CHARACTER(LEN=3)
  FUNCTION GET_CHANNEL_NAME_LSJ(L, S, J) RESULT(NAME)
    INTEGER, INTENT(IN) :: L, S, J
    CHARACTER(LEN=3) :: NAME

    ! Generate a name for the scattering channel based on L, S, and J
    ! Format: '<2S+1><L_letter><J>', e.g., '3SD' for S=1, L=2, J=2

    NAME = '   '
    SELECT CASE (L)
      CASE (0)
        NAME(2:2) = 'S'
      CASE (1)
        NAME(2:2) = 'P'
      CASE (2)
        NAME(2:2) = 'D'
      CASE (3)
        NAME(2:2) = 'F'
      CASE DEFAULT
        ! Map L >= 4 to corresponding spectroscopic letter (G=4, H=5, I=6, etc.)
        NAME(2:2) = CHAR(71 + (L-4))  ! 71 is ASCII for 'G'
    END SELECT

    WRITE(NAME(1:1), '(I1)') 2*S+1
    WRITE(NAME(3:3), '(I1)') J
  END FUNCTION GET_CHANNEL_NAME_LSJ

  !> \brief Get the spectroscopic name for a SCATTERING_CHANNEL object.
  !! \param[in] CHANNEL The channel object
  !! \return Name as CHARACTER(LEN=16)
  FUNCTION GET_CHANNEL_NAME_FROM_OBJECT(CHANNEL) RESULT(NAME)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    CHARACTER(LEN=16) :: NAME
    INTEGER :: ICH

    IF (.NOT.ALLOCATED(CHANNEL%L) .OR. .NOT.ALLOCATED(CHANNEL%S)) THEN
      NAME=""
      RETURN
    ENDIF
    ICH = 1
    NAME = GET_CHANNEL_NAME_LSJ(CHANNEL%L(ICH), CHANNEL%S(ICH), CHANNEL%J)
    DO ICH = 2, CHANNEL%NCH
      NAME = TRIM(NAME) // "-" // GET_CHANNEL_NAME_LSJ(CHANNEL%L(ICH), CHANNEL%S(ICH), CHANNEL%J)
    ENDDO
  END FUNCTION GET_CHANNEL_NAME_FROM_OBJECT

  !> \brief Get the number of channels (NCH).
  !! \param[in] CHANNEL The channel object
  !! \return Number of channels
  FUNCTION GET_CHANNEL_NCH(CHANNEL) RESULT (NCH)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: NCH
    NCH = CHANNEL%NCH
  END FUNCTION GET_CHANNEL_NCH

  !> \brief Get the L quantum number for a given channel index.
  !! \param[in] CHANNEL The channel object
  !! \param[in] I Channel index
  !! \return L quantum number
  FUNCTION GET_CHANNEL_L(CHANNEL, I) RESULT(L)
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER, INTENT(IN) :: I
    INTEGER :: L
    IF (I < 1 .OR. I > CHANNEL%NCH) THEN
      PRINT *, "Error: Index out of bounds in GET_CHANNEL_L"
      STOP
    ENDIF
    L = CHANNEL%L(I)
  END FUNCTION GET_CHANNEL_L

  !> \brief Get the S quantum number for a given channel index.
  !! \param[in] CHANNEL The channel object
  !! \param[in] I Channel index
  !! \return S quantum number
  FUNCTION GET_CHANNEL_S(CHANNEL, I) RESULT(S)
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER, INTENT(IN) :: I
    INTEGER :: S
    IF (I < 1 .OR. I > CHANNEL%NCH) THEN
      PRINT *, "Error: Index out of bounds in GET_CHANNEL_S"
      STOP
    ENDIF
    S = CHANNEL%S(I)
  END FUNCTION GET_CHANNEL_S

  !> \brief Get the T quantum number for a given channel index.
  !! \param[in] CHANNEL The channel object
  !! \param[in] I Channel index
  !! \return T quantum number
  FUNCTION GET_CHANNEL_T(CHANNEL, I) RESULT(T)
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER, INTENT(IN) :: I
    INTEGER :: T
    IF (I < 1 .OR. I > CHANNEL%NCH) THEN
      PRINT *, "Error: Index out of bounds in GET_CHANNEL_T"
      STOP
    ENDIF
    T = CHANNEL%T(I)
  END FUNCTION GET_CHANNEL_T

  !> \brief Get the J quantum number for a channel.
  !! \param[in] CHANNEL The channel object
  !! \return J quantum number
  FUNCTION GET_CHANNEL_J(CHANNEL) RESULT(J)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: J
    J = CHANNEL%J
  END FUNCTION GET_CHANNEL_J

  !> \brief Get the TZ quantum number for a channel.
  !! \param[in] CHANNEL The channel object
  !! \return TZ quantum number
  FUNCTION GET_CHANNEL_TZ(CHANNEL) RESULT(TZ)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: TZ
    TZ = CHANNEL%TZ
  END FUNCTION GET_CHANNEL_TZ

  !> \brief Check if the channel is coupled.
  !! \param[in] CHANNEL The channel object
  !! \return .TRUE. if coupled, .FALSE. otherwise
  FUNCTION IS_CHANNEL_COUPLED(CHANNEL) RESULT(COUPLED)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    LOGICAL :: COUPLED
    COUPLED = CHANNEL%COUPLED
  END FUNCTION IS_CHANNEL_COUPLED

  !> \brief Check if two channels are the same.
  !! \param[in] CHANNEL1 First channel
  !! \param[in] CHANNEL2 Second channel
  !! \return .TRUE. if the channels are the same, .FALSE. otherwise
  FUNCTION IS_SAME_CHANNEL(CHANNEL1, CHANNEL2) RESULT(SAME)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL1, CHANNEL2
    LOGICAL :: SAME
    INTEGER :: ICH

    SAME = .TRUE.
    IF (CHANNEL1%NCH /= CHANNEL2%NCH .OR. CHANNEL1%J /= CHANNEL2%J .OR. CHANNEL1%TZ /= CHANNEL2%TZ) THEN
      SAME = .FALSE.
      RETURN
    ENDIF

    DO ICH = 1, CHANNEL1%NCH
      IF (CHANNEL1%L(ICH) /= CHANNEL2%L(ICH) .OR. &
          CHANNEL1%S(ICH) /= CHANNEL2%S(ICH) .OR. &
          CHANNEL1%T(ICH) /= CHANNEL2%T(ICH)) THEN
        SAME = .FALSE.
        RETURN
      ENDIF
    ENDDO
  END FUNCTION IS_SAME_CHANNEL

  FUNCTION GET_CHANNEL_FROM_NAME(NAME) RESULT(CHANNEL)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    TYPE(SCATTERING_CHANNEL) :: CHANNEL
    CHARACTER(LEN=16) :: TMPNAME
    INTEGER :: L(2), S(2), J, TZ, NCH, POS, LEN1, LEN2

    IF (LEN_TRIM(NAME) < 3) THEN
      PRINT *, "Error: NAME must be at least 3 characters long.", NAME
      STOP
    ENDIF
    IF (LEN_TRIM(NAME) > 7) THEN
      PRINT *, "Error: NAME must not exceed 7 characters.", NAME
      STOP
    ENDIF
    ! Default TZ
    TZ = 0

    ! Check if the name contains a '-' (coupled channel)
    POS = INDEX(NAME, '-')
    IF (POS == 0) THEN
      IF (INDEX('SPDFGHIJKLMNOPQRSTUVWXYZ', NAME(2:2)) == 0) THEN
        PRINT *, "Error: Invalid spectroscopic letter in NAME."
        STOP
      ENDIF
      IF (NAME(1:1) < '1' .OR. NAME(1:1) > '9') THEN
        PRINT *, "Error: Invalid spin multiplicity in NAME."
        STOP
      ENDIF
      IF (NAME(3:3) < '0' .OR. NAME(3:3) > '9') THEN
        PRINT *, "Error: Invalid total angular momentum in NAME."
        STOP
      ENDIF

      READ(NAME(1:1), '(I1)') S(1)
      S(1) = (S(1) - 1) / 2  ! Convert from 2S+1 to S
      SELECT CASE (NAME(2:2))
        CASE ('S')
          L(1) = 0
        CASE ('P')
          L(1) = 1
        CASE ('D')
          L(1) = 2
        CASE ('F')
          L(1) = 3
        CASE DEFAULT
          L(1) = INDEX('GHIJKLMNOPQRSTUVWXYZ', NAME(2:2)) + 4 - 1
      END SELECT
      READ(NAME(3:3), '(I1)') J

      NCH = 1
      CHANNEL%COUPLED = .FALSE.
      CHANNEL%NCH = 1
      ALLOCATE(CHANNEL%L(1))
      ALLOCATE(CHANNEL%S(1))
      ALLOCATE(CHANNEL%T(1))
      CHANNEL%L(1) = L(1)
      CHANNEL%S(1) = S(1)
      CHANNEL%T(1) = EVALUATE_T(L(1), S(1))
      CHANNEL%J = J
      CHANNEL%TZ = TZ
    ELSE
      ! Coupled channel: split at '-'
      LEN1 = POS - 1
      LEN2 = LEN_TRIM(NAME) - POS
      TMPNAME = '                '
      TMPNAME(1:LEN1) = NAME(1:LEN1)

      ! First part
      IF (LEN1 < 3) THEN
        PRINT *, "Error: First channel name too short."
        STOP
      ENDIF
      READ(TMPNAME(1:1), '(I1)') S(1)
      S(1) = (S(1) - 1) / 2  ! Convert from 2S+1 to S
      SELECT CASE (TMPNAME(2:2))
        CASE ('S')
          L(1) = 0
        CASE ('P')
          L(1) = 1
        CASE ('D')
          L(1) = 2
        CASE ('F')
          L(1) = 3
        CASE DEFAULT
          L(1) = INDEX('GHIJKLMNOPQRSTUVWXYZ', TMPNAME(2:2)) + 4 - 1
      END SELECT
      READ(TMPNAME(3:3), '(I1)') J

      ! Second part
      TMPNAME = '                '
      TMPNAME(1:LEN2) = NAME(POS+1:POS+LEN2)
      IF (LEN2 < 3) THEN
        PRINT *, "Error: Second channel name too short."
        STOP
      ENDIF
      READ(TMPNAME(1:1), '(I1)') S(2)
      S(2) = (S(2) - 1) / 2  ! Convert from 2S+1 to S
      SELECT CASE (TMPNAME(2:2))
        CASE ('S')
          L(2) = 0
        CASE ('P')
          L(2) = 1
        CASE ('D')
          L(2) = 2
        CASE ('F')
          L(2) = 3
        CASE DEFAULT
          L(2) = INDEX('GHIJKLMNOPQRSTUVWXYZ', TMPNAME(2:2)) + 4 - 1
      END SELECT
      ! J must be the same for both parts, so skip reading again

      NCH = 2
      CHANNEL = init_scattering_channel(J, MOD(L(1), 2) == 0, TZ)

      ! Set both channels explicitly
      CALL REALLOCATE_1D_1_INT(CHANNEL%L, 2)
      CALL REALLOCATE_1D_1_INT(CHANNEL%S, 2)
      CALL REALLOCATE_1D_1_INT(CHANNEL%T, 2)
      CHANNEL%NCH = 2
      CHANNEL%L(1) = L(1)
      CHANNEL%S(1) = S(1)
      CHANNEL%T(1) = EVALUATE_T(L(1), S(1))
      CHANNEL%L(2) = L(2)
      CHANNEL%S(2) = S(2)
      CHANNEL%T(2) = EVALUATE_T(L(2), S(2))
      CHANNEL%COUPLED = .TRUE.
      CHANNEL%J = J
      CHANNEL%TZ = TZ
    ENDIF
  END FUNCTION GET_CHANNEL_FROM_NAME

  !> @brief Prepares the list of physical scattering channels up to given quantum number limits.
  !!
  !! This subroutine generates all possible physical scattering channels for given maximum
  !! orbital angular momentum (LMAX), total angular momentum (JMAX), and isospin projection (TZ).
  !! It allocates and fills the output array CHANNELS with all valid channels, filtering out
  !! unphysical combinations according to selection rules.
  !!
  !! \param[in] LMAX Maximum orbital angular momentum
  !! \param[in] JMAX Maximum total angular momentum
  !! \param[in] TZ Isospin projection
  !! \param[out] CHANNELS Array of physical scattering channels
  SUBROUTINE PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LMAX, JMAX, TZ
    TYPE(SCATTERING_CHANNEL), ALLOCATABLE, INTENT(OUT) :: CHANNELS(:)
    INTEGER :: ICH, NCH, L, S, J
    TYPE(SCATTERING_CHANNEL) :: CHANNEL

    NCH = 0
    DO L = 0, LMAX
      DO S = 0, 1
        DO J = ABS(L-S), MIN(L+S, JMAX)
          CALL SET_CHANNEL(CHANNEL, J, L, S, TZ)
          IF (.NOT.IS_PHYSICAL_CHANNEL(CHANNEL)) CYCLE
          IF ( J /= 0 .AND. L > J .AND. S == 1) CYCLE
          NCH = NCH + 1
        ENDDO
      ENDDO
    ENDDO

    ALLOCATE(CHANNELS(NCH))
    ICH = 1
    DO L = 0, LMAX
      DO S = 0, 1
        DO J = ABS(L-S), MIN(L+S, JMAX)
          CALL SET_CHANNEL(CHANNEL, J, L, S, TZ)
          IF (.NOT.IS_PHYSICAL_CHANNEL(CHANNEL)) CYCLE
          IF ( J /= 0 .AND. L > J .AND. S == 1) CYCLE
          CHANNELS(ICH) = CHANNEL
          ICH = ICH + 1
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE PREPARE_CHANNELS

  ! FUNCTION EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(FILENAME) RESULT(CHANNELS)
  !   IMPLICIT NONE
  !   CHARACTER(LEN=*), INTENT(IN) :: FILENAME
  !   TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
  !   TYPE(SCATTERING_CHANNEL) :: TMP

  !   CHARACTER(LEN=16) :: CH_NAME
  !   INTEGER :: I, N, LENF, START, ENDCH, COUNT
  !   LOGICAL :: FIRST_NUMBER_IS_READ = .FALSE.
  !   LOGICAL :: LAST_NUMBER_IS_READ = .FALSE.


  !   ! Find and print all channels in the filename

  !   LENF = LEN_TRIM(FILENAME)
  !   COUNT = 0

  !   I = 1
  !   DO WHILE (I <= LENF - 2)
  !     IF (.NOT.FIRST_NUMBER_IS_READ) THEN
  !       ! Look for the first digit
  !       IF (FILENAME(I:I) >= '0' .AND. FILENAME(I:I) <= '9') THEN
  !         START = I
  !         FIRST_NUMBER_IS_READ = .TRUE.
  !       ENDIF
  !       CH_NAME(1:1) = FILENAME(I:I)
  !     ELSE
  !       IF (FILENAME(I:I) < 'A' .AND. FILENAME(I:I) > 'Z') THEN
  !         CH_NAME(2:2) = FILENAME(I:I)
  !       ELSE
  !         IF (FILENAME(I:I) >= '0' .AND. FILENAME(I:I) <= '9') THEN
  !           CH_NAME(3:3) = FILENAME(I:I)
  !           LAST_NUMBER_IS_READ = .TRUE.
  !         ELSE
  !           FIRST_NUMBER_IS_READ = .FALSE.
  !         ENDIF
  !       ENDIF
  !     ENDIF
  !     I = I + 1
  !     IF (LAST_NUMBER_IS_READ) THEN
  !       TMP = GET_CHANNEL_FROM_NAME(TRIM(CH_NAME))
  !       =>>>>> TO FINISH
  !   ENDDO

  !   IF (COUNT > 0) THEN
  !     ALLOCATE(CHANNELS(COUNT))
  !     CHANNELS = TMP
  !   ELSE
  !     ALLOCATE(CHANNELS(0))
  !   END IF
  ! END FUNCTION EXTRACT_CHANNELS_FROM_WHOLE_FILENAME

  !> \brief Print all quantum numbers and info for a SCATTERING_CHANNEL object.
  !! \brief Print information about a scattering channel.
  !! Prints all quantum numbers and spectroscopic name for the given channel.
  !! \param[in] CHANNEL Scattering channel to print
  SUBROUTINE PRINT_SCATTERING_CHANNEL(CHANNEL)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: ICH
    PRINT *, '--- SCATTERING_CHANNEL INFO ---'
    PRINT *, '  J  =', CHANNEL%J
    PRINT *, '  TZ =', CHANNEL%TZ
    PRINT *, '  NCH=', CHANNEL%NCH
    PRINT *, '  COUPLED =', CHANNEL%COUPLED
    DO ICH = 1, CHANNEL%NCH
      PRINT *, '    Channel index:', ICH
      PRINT *, '      L =', CHANNEL%L(ICH)
      PRINT *, '      S =', CHANNEL%S(ICH)
      PRINT *, '      T =', CHANNEL%T(ICH)
    END DO
    PRINT *, '  Spectroscopic name: ', TRIM(GET_CHANNEL_NAME_FROM_OBJECT(CHANNEL))
    PRINT *, '------------------------------'
  END SUBROUTINE PRINT_SCATTERING_CHANNEL

  !> \brief Reset a scattering channel to default/uninitialized state.
  !! Deallocates arrays and resets all quantum numbers.
  !! \param[inout] CHANNEL Scattering channel to reset
  SUBROUTINE RESET_CHANNEL(CHANNEL)
    TYPE(SCATTERING_CHANNEL), INTENT(INOUT) :: CHANNEL

    CHANNEL%J = 0
    CHANNEL%TZ = 0
    CHANNEL%NCH = 0
    CHANNEL%COUPLED = .FALSE.
    IF (ALLOCATED(CHANNEL%L)) DEALLOCATE(CHANNEL%L)
    IF (ALLOCATED(CHANNEL%S)) DEALLOCATE(CHANNEL%S)
    IF (ALLOCATED(CHANNEL%T)) DEALLOCATE(CHANNEL%T)
  END SUBROUTINE RESET_CHANNEL


  !> \brief Compute all unique (L1, L2) combinations for a set of channels.
  !! Fills LEFT_RIGHT_L_COMBINATIONS with all unique pairs of L quantum numbers from CHANNELS.
  !! \param[in] CHANNELS Array of scattering channels
  !! \param[inout] LEFT_RIGHT_L_COMBINATIONS Output array of unique (L1, L2) pairs
  SUBROUTINE L_COMBINATIONS(CHANNELS, LEFT_RIGHT_L_COMBINATIONS)
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNELS(:)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: LEFT_RIGHT_L_COMBINATIONS(:,:)
    
    INTEGER :: N_COMB
    INTEGER, ALLOCATABLE :: TMP(:,:)
    INTEGER :: ICH, I, J, NCH, L1, L2, K, NCHANNELS
    LOGICAL :: FOUND

    NCHANNELS = SIZE(CHANNELS)
    ALLOCATE(TMP(4*NCHANNELS, 2))
    TMP = -1  ! Initialize TMP to zero
    N_COMB = 0

    ! Loop over all channels
    DO ICH = 1, NCHANNELS
      NCH = CHANNELS(ICH)%NCH
      ! Loop over all pairs (i, j) of L values
      DO I = 1, NCH
        L1 = CHANNELS(ICH)%L(I)
        DO J = 1, NCH
          L2 = CHANNELS(ICH)%L(J)
          ! Check if this combination already exists
          FOUND = .FALSE.
          DO K = 1, N_COMB
            IF (TMP(K,1) == L1 .AND. TMP(K,2) == L2) THEN
              FOUND = .TRUE.
              EXIT
            ENDIF
          ENDDO
          IF (FOUND) CYCLE  ! Skip if found
          ! If not found, add it
          N_COMB = N_COMB + 1
          TMP(N_COMB,1) = L1
          TMP(N_COMB,2) = L2
        ENDDO
      ENDDO
    ENDDO
    IF (ALLOCATED(LEFT_RIGHT_L_COMBINATIONS)) DEALLOCATE(LEFT_RIGHT_L_COMBINATIONS)
    ALLOCATE(LEFT_RIGHT_L_COMBINATIONS(N_COMB, 2))
    LEFT_RIGHT_L_COMBINATIONS = TMP(1:N_COMB, :)
  END SUBROUTINE L_COMBINATIONS

  !> \brief Convert isospin projection TZ to individual nucleon projections T1Z, T2Z.
  !! \param[in] TZ Total isospin projection
  !! \param[out] T1Z Isospin projection of nucleon 1
  !! \param[out] T2Z Isospin projection of nucleon 2
  SUBROUTINE TZ_TO_T1Z_T2Z(TZ, T1Z, T2Z)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: TZ
    INTEGER, INTENT(OUT) :: T1Z, T2Z

    IF (ABS(TZ) == 1) THEN
      T1Z = TZ
      T2Z = TZ
    ELSE
      T1Z = 1
      T2Z = -1
    END IF
  END SUBROUTINE TZ_TO_T1Z_T2Z

  !> \brief Convert individual nucleon isospin projections to total TZ.
  !! \param[in] T1Z Isospin projection of nucleon 1
  !! \param[in] T2Z Isospin projection of nucleon 2
  !! \param[out] TZ Total isospin projection
  SUBROUTINE T1Z_T2Z_TO_TZ(T1Z, T2Z, TZ)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: T1Z, T2Z
    INTEGER, INTENT(OUT) :: TZ

    IF (T1Z == T2Z) THEN
      TZ = T1Z
    ELSE
      TZ = 0
    END IF
  END SUBROUTINE T1Z_T2Z_TO_TZ

  !> @brief Evaluates T from L and S quantum numbers assuming (-1)^(L+S+T)==-1.
  !>
  !> @param[in] L Orbital angular momentum quantum number
  !> @param[in] S Spin quantum number
  !> 
  !> @return T Total angular momentum quantum number
  FUNCTION T_FROM_L_S(L, S) RESULT(T)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L, S
    INTEGER :: T

    T = MOD(MOD(L + S, 2) + 1, 2)
  END FUNCTION T_FROM_L_S

END MODULE QUANTUM_NUMBERS