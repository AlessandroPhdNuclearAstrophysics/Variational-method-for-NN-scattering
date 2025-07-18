!> \defgroup physics Physics Utilities
!! \brief All physics-related utilities.
!> \file quantum_numbers.f90
!! \defgroup quantum_numbers Quantum Numbers
!! \ingroup physics
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

  !> \brief Data type representing a scattering channel in quantum physics.
  !> \ingroup quantum_numbers
  !>
  !> This type stores quantum numbers and properties associated with a scattering channel,
  !> including total angular momentum, orbital angular momentum, spin, isospin, and coupling information.
  TYPE, PUBLIC :: SCATTERING_CHANNEL
  PRIVATE
    INTEGER :: J_ !< Total angular momentum quantum number
    INTEGER, ALLOCATABLE :: L_(:) !< Orbital angular momentum quantum numbers (allocatable array)
    INTEGER, ALLOCATABLE :: S_(:) !< Spin quantum numbers (allocatable array)
    INTEGER, ALLOCATABLE :: T_(:) !< Isospin quantum numbers (allocatable array)
    INTEGER :: TZ_ !< Isospin projection quantum number
    INTEGER :: NCH_ !< Number of channels
    LOGICAL :: COUPLED_ = .FALSE. !< Logical flag indicating if the channel is coupled (default: .FALSE.)
  CONTAINS 
    PROCEDURE :: SET => SET_CHANNEL
    PROCEDURE :: ASSIGN_CHANNEL
    PROCEDURE :: IS_SAME_CHANNEL
    PROCEDURE :: IS_NOT_SAME_CHANNEL
    PROCEDURE :: RESET => RESET_CHANNEL
    PROCEDURE :: PRINT => PRINT_SCATTERING_CHANNEL
    PROCEDURE :: TO_STRING => CHANNEL_TO_STRING
    PROCEDURE :: IS_COUPLED => IS_CHANNEL_COUPLED
    PROCEDURE :: NAME => GET_CHANNEL_NAME_FROM_OBJECT
    PROCEDURE :: J => GET_CHANNEL_J
    PROCEDURE :: TZ => GET_CHANNEL_TZ
    PROCEDURE :: NCH => GET_CHANNEL_NCH
    PROCEDURE :: L => GET_CHANNEL_L
    PROCEDURE :: S => GET_CHANNEL_S
    PROCEDURE :: T => GET_CHANNEL_T

    GENERIC :: ASSIGNMENT(=) => ASSIGN_CHANNEL
    GENERIC :: OPERATOR(==) => IS_SAME_CHANNEL
    GENERIC :: OPERATOR(/=) => IS_NOT_SAME_CHANNEL
  ENDTYPE SCATTERING_CHANNEL


  INTERFACE GET_CHANNEL_NAME
    MODULE PROCEDURE GET_CHANNEL_NAME_LSJ
    MODULE PROCEDURE GET_CHANNEL_NAME_FROM_OBJECT
  END INTERFACE GET_CHANNEL_NAME

  PUBLIC :: NEW_SCATTERING_CHANNEL
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

  !> \brief Assignment procedure for SCATTERING_CHANNEL type
  !! \param[out] LHS Left-hand side (target)
  !! \param[in] RHS Right-hand side (source)
  SUBROUTINE ASSIGN_CHANNEL(LHS, RHS)
    CLASS(SCATTERING_CHANNEL), INTENT(OUT) :: LHS
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: RHS

    ! Simple assignment - Fortran handles allocatable copying automatically
    LHS%J_ = RHS%J_
    LHS%TZ_ = RHS%TZ_
    LHS%NCH_ = RHS%NCH_
    LHS%COUPLED_ = RHS%COUPLED_
    LHS%L_ = RHS%L_     ! Automatic deep copy
    LHS%S_ = RHS%S_     ! Automatic deep copy
    LHS%T_ = RHS%T_     ! Automatic deep copy
  END SUBROUTINE ASSIGN_CHANNEL


  !> \brief Constructor for SCATTERING_CHANNEL.
  !! \ingroup quantum_numbers
  !! \param[in] J Total angular momentum
  !! \param[in] IS_EVEN Logical for parity
  !! \param[in] TZ Isospin projection
  !! \return Initialized SCATTERING_CHANNEL object
  FUNCTION NEW_SCATTERING_CHANNEL(J, IS_EVEN, TZ) RESULT(channel)
    INTEGER, INTENT(IN) :: J, TZ
    LOGICAL, INTENT(IN) :: IS_EVEN
    TYPE(SCATTERING_CHANNEL) :: CHANNEL

    INTEGER :: ICH

    IF (J == 0) THEN
      CHANNEL%NCH_ = 1
    ELSE
      CHANNEL%NCH_ = 2
    ENDIF

    ! ALLOCATE ARRAYS FOR L AND S
    CALL REALLOCATE(CHANNEL%L_, CHANNEL%NCH_)
    CALL REALLOCATE(CHANNEL%S_, CHANNEL%NCH_)
    CALL REALLOCATE(CHANNEL%T_, CHANNEL%NCH_)
    IF (J == 0) THEN
      IF (IS_EVEN) THEN
        CHANNEL%L_(1) = 0
        CHANNEL%S_(1) = 0
      ELSE
        CHANNEL%L_(1) = 1
        CHANNEL%S_(1) = 1
      ENDIF
    ELSE
      IF (MOD(J, 2) == 0) THEN
        IF (IS_EVEN) THEN
          CHANNEL%L_(1) = J
          CHANNEL%S_(1) = 0
          CHANNEL%L_(2) = J
          CHANNEL%S_(2) = 1
        ELSE
          CHANNEL%L_(1) = J-1
          CHANNEL%S_(1) = 1
          CHANNEL%L_(2) = J+1
          CHANNEL%S_(2) = 1
          CHANNEL%COUPLED_ = .TRUE.
        ENDIF
      ELSE
        IF (IS_EVEN) THEN
          CHANNEL%L_(1) = J - 1
          CHANNEL%S_(1) = 1
          CHANNEL%L_(2) = J + 1
          CHANNEL%S_(2) = 1
          CHANNEL%COUPLED_ = .TRUE.
        ELSE
          CHANNEL%L_(1) = J
          CHANNEL%S_(1) = 0
          CHANNEL%L_(2) = J
          CHANNEL%S_(2) = 1
        ENDIF
      ENDIF
    ENDIF

    CHANNEL%J_ = J
    DO ICH = 1, CHANNEL%NCH_
      CHANNEL%T_(ICH) = EVALUATE_T(CHANNEL%L_(ICH), CHANNEL%S_(ICH))
    ENDDO
    CHANNEL%TZ_ = TZ
  END FUNCTION NEW_SCATTERING_CHANNEL

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

    J = CHANNEL%J_
    IS_PHYSICAL = .TRUE.

    DO ICH=1, CHANNEL%NCH_
      L = CHANNEL%L_(ICH)
      S = CHANNEL%S_(ICH)
      T = CHANNEL%T_(ICH)
      TZ = CHANNEL%TZ_
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
  !! \ingroup quantum_numbers
  !! Sets the J, L, S, TZ quantum numbers and determines if the channel is coupled or not.
  !! Allocates and fills the L, S, T arrays for the channel.
  !! \param[inout] CHANNEL Scattering channel object to set
  !! \param[in] J Total angular momentum
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] TZ Isospin projection
  SUBROUTINE SET_CHANNEL(CHANNEL, J, L, S, TZ)
    CLASS(SCATTERING_CHANNEL), INTENT(INOUT) :: CHANNEL
    INTEGER, INTENT(IN) :: J, L, S, TZ
    LOGICAL :: IS_EVEN

    IF (.NOT.IS_LSJ_PHYSICAL(L, S, J)) THEN
      PRINT *, "QUANTUM_NUMBERS::SET_CHANNEL: Invalid quantum numbers L=", L, ", S=", S, ", J=", J
      STOP
    ENDIF

    IS_EVEN = MOD(L, 2) == 0

    CHANNEL%J_ = J
    CHANNEL%TZ_ = TZ
    CHANNEL%COUPLED_ = .FALSE.

    IF (J == 0) THEN
      CHANNEL%NCH_ = 1
      CALL REALLOCATE(CHANNEL%L_, 1)
      CALL REALLOCATE(CHANNEL%S_, 1)
      CALL REALLOCATE(CHANNEL%T_, 1)
      CHANNEL%L_(1) = L
      CHANNEL%S_(1) = S
      CHANNEL%T_(1) = EVALUATE_T(L, S)
    ELSE
      IF ((L==(J-1) .OR. L==(J+1))) THEN
        CHANNEL%NCH_ = 2
        CALL REALLOCATE(CHANNEL%L_, 2)
        CALL REALLOCATE(CHANNEL%S_, 2)
        CALL REALLOCATE(CHANNEL%T_, 2)
        CHANNEL%L_(1) = J-1
        CHANNEL%S_(1) = S
        CHANNEL%L_(2) = J+1
        CHANNEL%S_(2) = S
        CHANNEL%T_(1) = EVALUATE_T(J-1, S)
        CHANNEL%T_(2) = EVALUATE_T(J+1, S)
        CHANNEL%COUPLED_ = .TRUE.
      ELSE
        CHANNEL%NCH_ = 1
        CALL REALLOCATE(CHANNEL%L_, 1)
        CALL REALLOCATE(CHANNEL%S_, 1)
        CALL REALLOCATE(CHANNEL%T_, 1)
        CHANNEL%L_(1) = J
        CHANNEL%S_(1) = S
        CHANNEL%T_(1) = EVALUATE_T(J, S)
      ENDIF
    ENDIF
  END SUBROUTINE SET_CHANNEL

  !> \brief Get the spectroscopic name for a channel from L, S, J.
  !! \ingroup quantum_numbers
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] J Total angular momentum
  !! \return Name as CHARACTER(LEN=3)
  FUNCTION GET_CHANNEL_NAME_LSJ(L, S, J) RESULT(NAME)
    INTEGER, INTENT(IN) :: L, S, J
    CHARACTER(LEN=3) :: NAME

    ! Generate a name for the scattering channel based on L, S, and J
    ! Format: '<2S+1><L_letter><J>', e.g., '3SD' for S=1, L=2, J=2

    IF (IS_LSJ_PHYSICAL(L, S, J) .EQV. .FALSE.) THEN
      PRINT *, "QUANTUM_NUMBERS::GET_CHANNEL_NAME_LSJ: Invalid quantum numbers L=", L, ", S=", S, ", J=", J
      STOP
    ENDIF

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
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \return Name as CHARACTER(LEN=16)
  FUNCTION GET_CHANNEL_NAME_FROM_OBJECT(CHANNEL) RESULT(NAME)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    CHARACTER(LEN=16) :: NAME
    INTEGER :: ICH

    IF (.NOT.ALLOCATED(CHANNEL%L_) .OR. .NOT.ALLOCATED(CHANNEL%S_)) THEN
      NAME=""
      RETURN
    ENDIF
    ICH = 1
    NAME = GET_CHANNEL_NAME_LSJ(CHANNEL%L_(ICH), CHANNEL%S_(ICH), CHANNEL%J_)
    DO ICH = 2, CHANNEL%NCH_
      NAME = TRIM(NAME) // "-" // GET_CHANNEL_NAME_LSJ(CHANNEL%L_(ICH), CHANNEL%S_(ICH), CHANNEL%J_)
    ENDDO
  END FUNCTION GET_CHANNEL_NAME_FROM_OBJECT

  !> \brief Get the number of channels (NCH).
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \return Number of channels
  FUNCTION GET_CHANNEL_NCH(CHANNEL) RESULT (NCH)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: NCH
    NCH = CHANNEL%NCH_
  END FUNCTION GET_CHANNEL_NCH

  !> \brief Get the L quantum number for a given channel index.
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \param[in] I Channel index
  !! \return L quantum number
  FUNCTION GET_CHANNEL_L(CHANNEL, I) RESULT(L)
    IMPLICIT NONE
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER, OPTIONAL, INTENT(IN) :: I
    INTEGER :: L
    IF (PRESENT(I)) THEN
      IF (I < 1 .OR. I > CHANNEL%NCH_) THEN
        PRINT *, "Error: Index out of bounds in GET_CHANNEL_L"
        STOP
      ENDIF
      L = CHANNEL%L_(I)
      RETURN
    ENDIF
    L = CHANNEL%L_(1)  ! Default to first channel if no index provided
  END FUNCTION GET_CHANNEL_L

  !> \brief Get the S quantum number for a given channel index.
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \param[in] I Channel index
  !! \return S quantum number
  FUNCTION GET_CHANNEL_S(CHANNEL, I) RESULT(S)
    IMPLICIT NONE
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER, OPTIONAL, INTENT(IN) :: I
    INTEGER :: S
    IF (PRESENT(I)) THEN
      IF (I < 1 .OR. I > CHANNEL%NCH_) THEN
        PRINT *, "Error: Index out of bounds in GET_CHANNEL_S"
        STOP
      ENDIF
      S = CHANNEL%S_(I)
      RETURN
    ENDIF
    S = CHANNEL%S_(1)  ! Default to first channel if no index provided
  END FUNCTION GET_CHANNEL_S

  !> \brief Get the T quantum number for a given channel index.
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \param[in] I Channel index
  !! \return T quantum number
  FUNCTION GET_CHANNEL_T(CHANNEL, I) RESULT(T)
    IMPLICIT NONE
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER, OPTIONAL, INTENT(IN) :: I
    INTEGER :: T
    IF (PRESENT(I)) THEN
      IF (I < 1 .OR. I > CHANNEL%NCH_) THEN
        PRINT *, "Error: Index out of bounds in GET_CHANNEL_T"
        STOP
      ENDIF
      T = CHANNEL%T_(I)
      RETURN
    ENDIF
    T = CHANNEL%T_(1)  ! Default to first channel if no index provided
  END FUNCTION GET_CHANNEL_T

  !> \brief Get the J quantum number for a channel.
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \return J quantum number
  FUNCTION GET_CHANNEL_J(CHANNEL) RESULT(J)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: J
    J = CHANNEL%J_
  END FUNCTION GET_CHANNEL_J

  !> \brief Get the TZ quantum number for a channel.
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \return TZ quantum number
  FUNCTION GET_CHANNEL_TZ(CHANNEL) RESULT(TZ)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: TZ
    TZ = CHANNEL%TZ_
  END FUNCTION GET_CHANNEL_TZ

  !> \brief Check if the channel is coupled.
  !! \ingroup quantum_numbers
  !! \param[in] CHANNEL The channel object
  !! \return .TRUE. if coupled, .FALSE. otherwise
  FUNCTION IS_CHANNEL_COUPLED(CHANNEL) RESULT(COUPLED)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    LOGICAL :: COUPLED
    COUPLED = CHANNEL%COUPLED_
  END FUNCTION IS_CHANNEL_COUPLED

  !> \brief Check if two channels are the same.
  !! \ingroup quantum_numbers 
  !! \param[in] CHANNEL1 First channel
  !! \param[in] CHANNEL2 Second channel
  !! \return .TRUE. if the channels are the same, .FALSE. otherwise
  FUNCTION IS_SAME_CHANNEL(CHANNEL1, CHANNEL2) RESULT(SAME)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL1
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL2
    LOGICAL :: SAME
    INTEGER :: ICH

    SAME = .TRUE.
    IF (CHANNEL1%NCH_ /= CHANNEL2%NCH_ .OR. CHANNEL1%J_ /= CHANNEL2%J_ .OR. CHANNEL1%TZ_ /= CHANNEL2%TZ_) THEN
      SAME = .FALSE.
      RETURN
    ENDIF

    DO ICH = 1, CHANNEL1%NCH_
      IF (CHANNEL1%L_(ICH) /= CHANNEL2%L_(ICH) .OR. &
          CHANNEL1%S_(ICH) /= CHANNEL2%S_(ICH) .OR. &
          CHANNEL1%T_(ICH) /= CHANNEL2%T_(ICH)) THEN
        SAME = .FALSE.
        RETURN
      ENDIF
    ENDDO
  END FUNCTION IS_SAME_CHANNEL

  !> \brief Check if two channels are not the same.
  !! \ingroup quantum_numbers 
  !! \param[in] CHANNEL1 First channel
  !! \param[in] CHANNEL2 Second channel
  !! \return .TRUE. if the channels are different, .FALSE. otherwise
  FUNCTION IS_NOT_SAME_CHANNEL(CHANNEL1, CHANNEL2) RESULT(NOT_SAME)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL1
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL2
    LOGICAL :: NOT_SAME

    NOT_SAME = .NOT. IS_SAME_CHANNEL(CHANNEL1, CHANNEL2)
  END FUNCTION IS_NOT_SAME_CHANNEL

  !> \brief Get a SCATTERING_CHANNEL object from its spectroscopic name.
  !! \ingroup quantum_numbers
  !! \param[in] NAME Spectroscopic name (e.g., '3SD-1SD' or '3F2')
  !! \return SCATTERING_CHANNEL object
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
      CHANNEL%COUPLED_ = .FALSE.
      CHANNEL%NCH_ = 1
      ALLOCATE(CHANNEL%L_(1))
      ALLOCATE(CHANNEL%S_(1))
      ALLOCATE(CHANNEL%T_(1))
      CHANNEL%L_(1) = L(1)
      CHANNEL%S_(1) = S(1)
      CHANNEL%T_(1) = EVALUATE_T(L(1), S(1))
      CHANNEL%J_ = J
      CHANNEL%TZ_ = TZ
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
      CHANNEL = NEW_SCATTERING_CHANNEL(J, MOD(L(1), 2) == 0, TZ)

      ! Set both channels explicitly
      CALL REALLOCATE(CHANNEL%L_, 2)
      CALL REALLOCATE(CHANNEL%S_, 2)
      CALL REALLOCATE(CHANNEL%T_, 2)
      CHANNEL%NCH_ = 2
      CHANNEL%L_(1) = L(1)
      CHANNEL%S_(1) = S(1)
      CHANNEL%T_(1) = EVALUATE_T(L(1), S(1))
      CHANNEL%L_(2) = L(2)
      CHANNEL%S_(2) = S(2)
      CHANNEL%T_(2) = EVALUATE_T(L(2), S(2))
      CHANNEL%COUPLED_ = .TRUE.
      CHANNEL%J_ = J
      CHANNEL%TZ_ = TZ
    ENDIF
  END FUNCTION GET_CHANNEL_FROM_NAME

  !> @brief Prepares the list of physical scattering channels up to given quantum number limits.
  !! \ingroup quantum_numbers
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
          CALL CHANNEL%SET(J, L, S, TZ)
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
          CALL CHANNEL%SET(J, L, S, TZ)
          IF (.NOT.IS_PHYSICAL_CHANNEL(CHANNEL)) CYCLE
          IF ( J /= 0 .AND. L > J .AND. S == 1) CYCLE
          CHANNELS(ICH) = CHANNEL
          ICH = ICH + 1
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE PREPARE_CHANNELS

  !> \brief Print all quantum numbers and info for a SCATTERING_CHANNEL object.
  !! \ingroup quantum_numbers
  !! \brief Print information about a scattering channel.
  !! Prints all quantum numbers and spectroscopic name for the given channel.
  !! \param[in] CHANNEL Scattering channel to print
  SUBROUTINE PRINT_SCATTERING_CHANNEL(CHANNEL, UNIT)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    INTEGER :: UNIT_
    INTEGER :: ICH
    IF (PRESENT(UNIT)) THEN
      UNIT_ = UNIT
    ELSE
      UNIT_ = 6  ! Default to standard output
    ENDIF
    WRITE(UNIT_, *) '--- SCATTERING_CHANNEL INFO ---'
    WRITE(UNIT_, *) '  J  =', CHANNEL%J_
    WRITE(UNIT_, *) '  TZ =', CHANNEL%TZ_
    WRITE(UNIT_, *) '  NCH=', CHANNEL%NCH_
    WRITE(UNIT_, *) '  COUPLED =', CHANNEL%COUPLED_
    DO ICH = 1, CHANNEL%NCH_
      WRITE(UNIT_, *) '    Channel index:', ICH
      WRITE(UNIT_, *) '      L =', CHANNEL%L_(ICH)
      WRITE(UNIT_, *) '      S =', CHANNEL%S_(ICH)
      WRITE(UNIT_, *) '      T =', CHANNEL%T_(ICH)
    END DO
    WRITE(UNIT_, *) '  Spectroscopic name: ', TRIM(GET_CHANNEL_NAME_FROM_OBJECT(CHANNEL))
    WRITE(UNIT_, *) '------------------------------'
  END SUBROUTINE PRINT_SCATTERING_CHANNEL

  SUBROUTINE CHANNEL_TO_STRING(CHANNEL, STRING)
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    CHARACTER(LEN=*), INTENT(OUT) :: STRING
    CHARACTER(LEN=256) :: TEMP_STRING
    INTEGER :: ICH

    WRITE(TEMP_STRING, *) '--- SCATTERING_CHANNEL INFO ---'
    STRING = TEMP_STRING
    WRITE(TEMP_STRING, '(A,I0)') '  J  = ', CHANNEL%J_
    STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
    WRITE(TEMP_STRING, '(A,I0)') '  TZ = ', CHANNEL%TZ_
    STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
    WRITE(TEMP_STRING, '(A,I0)') '  NCH=', CHANNEL%NCH_
    STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
    WRITE(TEMP_STRING, '(A,L1)') '  COUPLED = ', CHANNEL%COUPLED_
    STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
    DO ICH = 1, CHANNEL%NCH_
      WRITE(TEMP_STRING, '(A,I0)') '    Channel index: ', ICH
      STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
      WRITE(TEMP_STRING, '(A,I0)') '      L = ', CHANNEL%L_(ICH)
      STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
      WRITE(TEMP_STRING, '(A,I0)') '      S = ', CHANNEL%S_(ICH)
      STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
      WRITE(TEMP_STRING, '(A,I0)') '      T = ', CHANNEL%T_(ICH)
      STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
    END DO
    WRITE(TEMP_STRING, '(A,A)') '  Spectroscopic name: ', TRIM(GET_CHANNEL_NAME_FROM_OBJECT(CHANNEL))
    STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
    WRITE(TEMP_STRING, '(A)') '------------------------------'
    STRING = TRIM(STRING) // NEW_LINE('A') // TEMP_STRING
    
    STRING = TRIM(STRING)  ! Remove trailing spaces
  END SUBROUTINE CHANNEL_TO_STRING

  !> \brief Reset a scattering channel to default/uninitialized state.
  !! \ingroup quantum_numbers
  !! Deallocates arrays and resets all quantum numbers.
  !! \param[inout] CHANNEL Scattering channel to reset
  SUBROUTINE RESET_CHANNEL(CHANNEL)
    CLASS(SCATTERING_CHANNEL), INTENT(INOUT) :: CHANNEL

    CHANNEL%J_ = 0
    CHANNEL%TZ_ = 0
    CHANNEL%NCH_ = 0
    CHANNEL%COUPLED_ = .FALSE.
    IF (ALLOCATED(CHANNEL%L_)) DEALLOCATE(CHANNEL%L_)
    IF (ALLOCATED(CHANNEL%S_)) DEALLOCATE(CHANNEL%S_)
    IF (ALLOCATED(CHANNEL%T_)) DEALLOCATE(CHANNEL%T_)
  END SUBROUTINE RESET_CHANNEL


  !> \brief Compute all unique (L1, L2) combinations for a set of channels.
  !! \ingroup quantum_numbers
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
      NCH = CHANNELS(ICH)%NCH_
      ! Loop over all pairs (i, j) of L values
      DO I = 1, NCH
        L1 = CHANNELS(ICH)%L_(I)
        DO J = 1, NCH
          L2 = CHANNELS(ICH)%L_(J)
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
  !! \ingroup quantum_numbers
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
  !! \ingroup quantum_numbers
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
  !! \ingroup quantum_numbers
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