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

  PUBLIC :: &
    init_scattering_channel, &
    SET_CHANNEL, &
    GET_CHANNEL_NAME, &
    GET_CHANNEL_NCH, &
    GET_CHANNEL_L, &
    GET_CHANNEL_S, &
    GET_CHANNEL_T, &
    GET_CHANNEL_J, &
    GET_CHANNEL_TZ, &
    IS_CHANNEL_COUPLED, &
    IS_PHYSICAL_CHANNEL, &
    GET_CHANNEL_FROM_NAME

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
      CHANNEL%T(ICH) = EVALUATE_T(CHANNEL%L(ICH), CHANNEL%S(ICH), TZ)
    ENDDO
    CHANNEL%TZ = TZ
  END FUNCTION init_scattering_channel

  !> \brief Evaluate isospin T for given L, S, and TZ.
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] TZ Isospin projection
  !! \return Isospin T
  FUNCTION EVALUATE_T(L, S, TZ) RESULT(T)
    INTEGER, INTENT(IN) :: L, S, TZ
    INTEGER :: T

    ! Calculate T based on the values of J, L, S, and TZ
    T = MOD( MOD(L+S, 2) + 1 , 2)
    IF (ABS(TZ) > T) THEN
      PRINT *, "Error: Invalid value for TZ. It must be less than or equal to T."
      STOP
    ENDIF
  END FUNCTION EVALUATE_T

  !> \brief Check if a channel is physical.
  !! \param[in] CHANNEL The SCATTERING_CHANNEL object
  !! \return .TRUE. if physical, .FALSE. otherwise
  FUNCTION IS_PHYSICAL_CHANNEL(CHANNEL) RESULT(IS_PHYSICAL)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    LOGICAL :: IS_PHYSICAL
    INTEGER :: L, S, J, T, ICH
    
    J = CHANNEL%J
    IS_PHYSICAL = .TRUE.

    DO ICH=1, CHANNEL%NCH
      L = CHANNEL%L(ICH)
      S = CHANNEL%S(ICH)
      T = CHANNEL%T(ICH)

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
  !! \param[inout] CHANNEL The channel to set
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
      CHANNEL%T(1) = EVALUATE_T(L, S, TZ)
    ELSE
      IF ((L.EQ.(J-1) .OR. L.EQ.(J+1))) THEN
        CHANNEL%NCH = 2
        CALL REALLOCATE_1D_1_INT(CHANNEL%L, 2)
        CALL REALLOCATE_1D_1_INT(CHANNEL%S, 2)
        CALL REALLOCATE_1D_1_INT(CHANNEL%T, 2)
        CHANNEL%L(1) = J-1
        CHANNEL%S(1) = S
        CHANNEL%L(2) = J+1
        CHANNEL%S(2) = S
        CHANNEL%T(1) = EVALUATE_T(J-1, S, TZ)
        CHANNEL%T(2) = EVALUATE_T(J+1, S, TZ)
        CHANNEL%COUPLED = .TRUE.
      ELSE
        CHANNEL%NCH = 1
        CALL REALLOCATE_1D_1_INT(CHANNEL%L, 1)
        CALL REALLOCATE_1D_1_INT(CHANNEL%S, 1)
        CALL REALLOCATE_1D_1_INT(CHANNEL%T, 1)
        CHANNEL%L(1) = J
        CHANNEL%S(1) = S
        CHANNEL%T(1) = EVALUATE_T(J, S, TZ)
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
    CHARACTER(LEN=:), ALLOCATABLE :: PARTS(:)
    CHARACTER(LEN=16) :: TMPNAME
    INTEGER :: L(2), S(2), J, TZ, NCH, ICH, POS, LEN1, LEN2

    ! Default TZ
    TZ = 0

    ! Check if the name contains a '-' (coupled channel)
    POS = INDEX(NAME, '-')
    IF (POS == 0) THEN
      ! Uncoupled channel
      IF (LEN_TRIM(NAME) < 3) THEN
        PRINT *, "Error: NAME must be at least 3 characters long."
        STOP
      ENDIF
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
      CHANNEL = init_scattering_channel(J, MOD(L(1), 2) == 0, TZ)
      CALL SET_CHANNEL(CHANNEL, J, L(1), S(1), TZ)
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
      CHANNEL%T(1) = EVALUATE_T(L(1), S(1), TZ)
      CHANNEL%L(2) = L(2)
      CHANNEL%S(2) = S(2)
      CHANNEL%T(2) = EVALUATE_T(L(2), S(2), TZ)
      CHANNEL%COUPLED = .TRUE.
      CHANNEL%J = J
      CHANNEL%TZ = TZ
    ENDIF
  END FUNCTION GET_CHANNEL_FROM_NAME
END MODULE QUANTUM_NUMBERS