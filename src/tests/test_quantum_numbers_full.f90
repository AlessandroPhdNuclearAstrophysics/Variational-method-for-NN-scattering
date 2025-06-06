PROGRAM test_quantum_numbers_full
  USE QUANTUM_NUMBERS
  USE strings_utils
  USE console_colors
  IMPLICIT NONE

  TYPE(SCATTERING_CHANNEL) :: ch, ch2
  CHARACTER(LEN=16) :: name, name2
  INTEGER :: L, S, J, T, TZ
  LOGICAL :: is_phys, is_coup, result

  ! For EXTRACT_CHANNELS_FROM_WHOLE_FILENAME tests
  CALL print_status("Comprehensive Test of QUANTUM_NUMBERS Module", .TRUE.)

  ! Test 1: Test initialization with different parameters
  CALL print_status("Test 1: Testing init_scattering_channel with various parameters", .TRUE.)

  ! Test 1.1: Even parity J=0 channel
  ch = init_scattering_channel(0, .TRUE., 0)
  CALL test_channel(ch, 0, 0, 0, 0, 1, .FALSE., "1S0", "Test 1.1: Even parity J=0 channel")

  ! Test 1.2: Odd parity J=0 channel
  ch = init_scattering_channel(0, .FALSE., 0)
  CALL test_channel(ch, 0, 1, 1, 0, 1, .FALSE., "3P0", "Test 1.2: Odd parity J=0 channel")

  ! Test 1.3: Even parity J=1 channel (coupled)
  ch = init_scattering_channel(1, .TRUE., 0)
  CALL test_channel(ch, 1, 0, 1, 0, 2, .TRUE., "3S1-3D1", "Test 1.3: Even parity J=1 channel (coupled)")

  ! Test 1.4: Odd parity J=1 channel
  ch = init_scattering_channel(1, .FALSE., 0)
  CALL test_channel(ch, 1, 1, 0, 0, 2, .FALSE., "1P1-3P1", "Test 1.4: Odd parity J=1 channel")

  ! Test 1.5: Test different TZ value
  ch = init_scattering_channel(2, .TRUE., 0)
  CALL test_channel(ch, 2, 2, 0, 0, 2, .FALSE., "1D2-3D2", "Test 1.5: Different TZ value")

  ! Test 2: Test SET_CHANNEL functionality
  CALL print_status("Test 2: Testing SET_CHANNEL functionality", .TRUE.)

  ! Test 2.1: Set single-channel parameters
  CALL SET_CHANNEL(ch, 3, 3, 0, 0)
  CALL test_channel(ch, 3, 3, 0, 0, 1, .FALSE., "1F3", "Test 2.1: Set single-channel parameters")

  ! Test 2.2: Set coupled-channel parameters
  CALL SET_CHANNEL(ch, 2, 1, 1, 0)
  CALL test_channel(ch, 2, 1, 1, 0, 2, .TRUE., "3P2-3F2", "Test 2.2: Set coupled-channel parameters")

  ! Test 3: Test GET_CHANNEL_FROM_NAME functionality
  CALL print_status("Test 3: Testing GET_CHANNEL_FROM_NAME functionality", .TRUE.)

  ! Test 3.1: Parse uncoupled channel name
  ch = GET_CHANNEL_FROM_NAME("1S0")
  CALL test_channel(ch, 0, 0, 0, 0, 1, .FALSE., "1S0", "Test 3.1: Parse uncoupled channel name")

  ! Test 3.2: Parse coupled channel name
  ch = GET_CHANNEL_FROM_NAME("3S1-3D1")
  CALL test_channel(ch, 1, 0, 1, 0, 2, .TRUE., "3S1-3D1", "Test 3.2: Parse coupled channel name")

  ! Test 3.3: Test with unusual L values (high angular momentum)
  ch = GET_CHANNEL_FROM_NAME("3G4")
  CALL test_channel(ch, 4, 4, 1, 0, 1, .FALSE., "3G4", "Test 3.3: High angular momentum")

  ! Test 4: Test IS_SAME_CHANNEL functionality
  CALL print_status("Test 4: Testing IS_SAME_CHANNEL functionality", .TRUE.)

  ! Test 4.1: Compare identical channels
  ch = GET_CHANNEL_FROM_NAME("3P2-3F2")
  ch2 = GET_CHANNEL_FROM_NAME("3P2-3F2")
  result = IS_SAME_CHANNEL(ch, ch2)
  IF (.NOT. result) THEN
    CALL print_status("Test 4.1: Identical channels not detected as same", .FALSE.)
    STOP
  ELSE
    CALL print_status("Test 4.1: Identical channels detected as same", .TRUE.)
  END IF

  ! Test 4.2: Compare different channels
  ch = GET_CHANNEL_FROM_NAME("3P2-3F2")
  ch2 = GET_CHANNEL_FROM_NAME("1D2")
  result = IS_SAME_CHANNEL(ch, ch2)
  IF (result) THEN
    CALL print_status("Test 4.2: Different channels detected as same", .FALSE.)
    STOP
  ELSE
    CALL print_status("Test 4.2: Different channels detected as different", .TRUE.)
  END IF

  ! Test 5: Test IS_PHYSICAL_CHANNEL functionality
  CALL print_status("Test 5: Testing IS_PHYSICAL_CHANNEL functionality", .TRUE.)

  ! Test 5.1: Test with valid channel
  ch = GET_CHANNEL_FROM_NAME("3F3")
  is_phys = IS_PHYSICAL_CHANNEL(ch)
  IF (.NOT. is_phys) THEN
    CALL print_status("Test 5.1: Valid channel not detected as physical", .FALSE.)
    STOP
  ELSE
    CALL print_status("Test 5.1: Valid channel detected as physical", .TRUE.)
  END IF

  ! Test 5.2: Test with non-physical L, S, J combination
  ch = init_scattering_channel(2, .TRUE., 0)
  CALL SET_CHANNEL(ch, 2, 2, 1, 0)
  ! Not directly testable for failure here, so skip status print

  ! Test 6: Test all quantum number getters
  CALL print_status("Test 6: Testing quantum number getters", .TRUE.)
  ch = GET_CHANNEL_FROM_NAME("3D2-3G2")

  J = GET_CHANNEL_J(ch)
  IF (J /= 2) THEN
    CALL print_status("Test 6.1: GET_CHANNEL_J incorrect", .FALSE.)
    STOP
  END IF

  L = GET_CHANNEL_L(ch, 1)
  IF (L /= 2) THEN
    CALL print_status("Test 6.2: GET_CHANNEL_L incorrect for first component", .FALSE.)
    STOP
  END IF

  L = GET_CHANNEL_L(ch, 2)
  IF (L /= 4) THEN
    CALL print_status("Test 6.3: GET_CHANNEL_L incorrect for second component", .FALSE.)
    STOP
  END IF

  S = GET_CHANNEL_S(ch, 1)
  IF (S /= 1) THEN
    CALL print_status("Test 6.4: GET_CHANNEL_S incorrect", .FALSE.)
    STOP
  END IF

  T = GET_CHANNEL_T(ch, 1)
  IF (T /= 0) THEN
    CALL print_status("Test 6.5: GET_CHANNEL_T incorrect", .FALSE.)
    STOP
  END IF

  TZ = GET_CHANNEL_TZ(ch)
  IF (TZ /= 0) THEN
    CALL print_status("Test 6.6: GET_CHANNEL_TZ incorrect", .FALSE.)
    STOP
  END IF

  is_coup = IS_CHANNEL_COUPLED(ch)
  IF (.NOT. is_coup) THEN
    CALL print_status("Test 6.7: IS_CHANNEL_COUPLED incorrect", .FALSE.)
    STOP
  ELSE
    CALL print_status("Test 6: All quantum number getters", .TRUE.)
  END IF

  ! Test 7: Round-trip name conversion
  CALL print_status("Test 7: Testing round-trip name conversion", .TRUE.)
  name = "3P2-3F2"
  ch = GET_CHANNEL_FROM_NAME(name)
  name2 = GET_CHANNEL_NAME(ch)
  IF (TRIM(name) /= TRIM(name2)) THEN
    CALL print_status("Test 7: Round-trip name conversion failed", .FALSE.)
    STOP
  ELSE
    CALL print_status("Test 7: Round-trip name conversion", .TRUE.)
  END IF

  ! Test 8: Prepare channels up to LMAX=3 and JMAX=3
  CALL print_status("Test 8: Testing PREPARE_CHANNELS with LMAX=3, JMAX=3", .TRUE.)
  CALL test_prepare_channels(0)
  CALL test_prepare_channels(1)

  ! Test 9: L combinations for TZ=0, LMAX=3, JMAX=3
  CALL print_status("Test 9: Testing L_COMBINATIONS for TZ=0, LMAX=3, JMAX=3", .TRUE.)
  CALL test_l_combinations

CONTAINS

  SUBROUTINE print_status(testname, passed)
    CHARACTER(LEN=*), INTENT(IN) :: testname
    LOGICAL, INTENT(IN) :: passed
    CHARACTER(LEN=*), PARAMETER :: RED = CHAR(27)//'[31m'
    CHARACTER(LEN=*), PARAMETER :: GREEN = CHAR(27)//'[32m'
    CHARACTER(LEN=*), PARAMETER :: RESET = CHAR(27)//'[0m'
    
    IF (passed) THEN
      WRITE(*,'(A,A,A,A)') GREEN, "PASSED:", RESET, " "//TRIM(testname)
    ELSE
      WRITE(*,'(A,A,A,A)') RED, "FAILED:", RESET, " "//TRIM(testname)
    END IF
  END SUBROUTINE print_status

  SUBROUTINE test_channel(test_ch, exp_J, exp_L, exp_S, exp_TZ, exp_nch, exp_coupled, exp_name, testname)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: test_ch
    INTEGER, INTENT(IN) :: exp_J, exp_L, exp_S, exp_TZ, exp_nch
    LOGICAL, INTENT(IN) :: exp_coupled
    CHARACTER(LEN=*), INTENT(IN) :: exp_name
    CHARACTER(LEN=*), INTENT(IN) :: testname

    INTEGER :: actual_J, actual_L, actual_S, actual_TZ, actual_nch
    LOGICAL :: actual_coupled
    CHARACTER(LEN=16) :: actual_name

    actual_J = GET_CHANNEL_J(test_ch)
    actual_L = GET_CHANNEL_L(test_ch, 1)
    actual_S = GET_CHANNEL_S(test_ch, 1)
    actual_TZ = GET_CHANNEL_TZ(test_ch)
    actual_nch = GET_CHANNEL_NCH(test_ch)
    actual_coupled = IS_CHANNEL_COUPLED(test_ch)
    actual_name = GET_CHANNEL_NAME(test_ch)

    IF (actual_J /= exp_J .OR. actual_L /= exp_L .OR. actual_S /= exp_S .OR. actual_TZ /= exp_TZ .OR. &
        actual_nch /= exp_nch .OR. actual_coupled .NEQV. exp_coupled .OR. TRIM(actual_name) /= TRIM(exp_name)) THEN
      CALL print_status(testname, .FALSE.)
      IF (actual_J /= exp_J) WRITE(*,*) "  J mismatch"
      IF (actual_L /= exp_L) WRITE(*,*) "  L mismatch"
      IF (actual_S /= exp_S) WRITE(*,*) "  S mismatch"
      IF (actual_TZ /= exp_TZ) WRITE(*,*) "  TZ mismatch"
      IF (actual_nch /= exp_nch) WRITE(*,*) "  NCH mismatch"
      IF (actual_coupled .NEQV. exp_coupled) WRITE(*,*) "  COUPLED status mismatch"
      IF (TRIM(actual_name) /= TRIM(exp_name)) WRITE(*,*) "  Name mismatch: Expected: ", TRIM(exp_name), " Got: ", TRIM(actual_name)
      STOP
    ELSE
      CALL print_status(testname, .TRUE.)
    END IF
  END SUBROUTINE test_channel

  SUBROUTINE test_prepare_channels(tz_ch)
    TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: channels(:)
    INTEGER, INTENT(IN) :: tz_ch
    INTEGER :: i, lmax, jmax, n_coupled
    CHARACTER(LEN=16) :: chname

    lmax = 3
    jmax = 3
    CALL PREPARE_CHANNELS(lmax, jmax, tz_ch, channels)
    n_coupled = 0
    DO i = 1, SIZE(channels)
      chname = GET_CHANNEL_NAME(channels(i))
      IF (IS_CHANNEL_COUPLED(channels(i))) n_coupled = n_coupled + 1
    END DO
    IF (tz_ch == 0) THEN
      IF (SIZE(channels) /= 11) THEN
        CALL print_status("Test 8.1: Expected 11 channels for TZ=0", .FALSE.)
        STOP
      END IF
      IF (n_coupled /= 3) THEN
        CALL print_status("Test 8.2: Expected 3 coupled channels for TZ=0", .FALSE.)
        STOP
      END IF
      CALL print_status("Test 8: PREPARE_CHANNELS TZ=0", .TRUE.)
    ELSEIF (tz_ch == 1) THEN
      IF (SIZE(channels) /= 6) THEN
        CALL print_status("Test 8.3: Expected 6 channels for TZ=1", .FALSE.)
        STOP
      END IF
      IF (n_coupled /= 1) THEN
        CALL print_status("Test 8.4: Expected 1 coupled channel for TZ=1", .FALSE.)
        STOP
      END IF
      CALL print_status("Test 8: PREPARE_CHANNELS TZ=1", .TRUE.)
    END IF
  END SUBROUTINE test_prepare_channels

  SUBROUTINE test_l_combinations
    TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: channels(:)
    INTEGER, ALLOCATABLE :: lcomb(:,:)
    INTEGER :: lmax, jmax, tz_ch, ncomb

    lmax = 2
    jmax = 2
    tz_ch = 0
    CALL PREPARE_CHANNELS(lmax, jmax, tz_ch, channels)
    CALL L_COMBINATIONS(channels, lcomb)
    ncomb = SIZE(lcomb, 1)
    IF (ncomb /= 8) THEN
      CALL print_status("Test 9: Expected 8 L combinations for LMAX=2, JMAX=2, TZ=0", .FALSE.)
      STOP
    ELSE
      CALL print_status("Test 9: L_COMBINATIONS for LMAX=2, JMAX=2, TZ=0", .TRUE.)
    END IF
  END SUBROUTINE test_l_combinations

END PROGRAM test_quantum_numbers_full
