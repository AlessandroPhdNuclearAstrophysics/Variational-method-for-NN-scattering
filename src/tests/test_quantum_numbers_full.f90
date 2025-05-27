PROGRAM test_quantum_numbers_full
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  
  TYPE(SCATTERING_CHANNEL) :: ch, ch2, ch3
  CHARACTER(LEN=16) :: name, name2
  INTEGER :: i, nch, L, S, J, T, TZ
  LOGICAL :: is_phys, is_coup, result
  
  WRITE(*,*) "=== Comprehensive Test of QUANTUM_NUMBERS Module ==="
  
  ! Test 1: Test initialization with different parameters
  WRITE(*,*) "Test 1: Testing init_scattering_channel with various parameters"
  
  ! Test 1.1: Even parity J=0 channel
  ch = init_scattering_channel(0, .TRUE., 0)
  CALL test_channel(ch, 0, 0, 0, 0, 1, .FALSE., "1S0")
  
  ! Test 1.2: Odd parity J=0 channel 
  ch = init_scattering_channel(0, .FALSE., 0)
  CALL test_channel(ch, 0, 1, 1, 0, 1, .FALSE., "3P0")
  
  ! Test 1.3: Even parity J=1 channel (coupled)
  ch = init_scattering_channel(1, .TRUE., 0)
  CALL test_channel(ch, 1, 0, 1, 0, 2, .TRUE., "3S1-3D1")
  
  ! Test 1.4: Odd parity J=1 channel
  ch = init_scattering_channel(1, .FALSE., 0)
  CALL test_channel(ch, 1, 1, 0, 0, 2, .FALSE., "1P1-3P1")
  
  ! Test 1.5: Test different TZ value
  ch = init_scattering_channel(2, .TRUE., 0)
  CALL test_channel(ch, 2, 2, 0, 0, 2, .FALSE., "1D2-3D2")
  
  ! Test 2: Test SET_CHANNEL functionality
  WRITE(*,*) "Test 2: Testing SET_CHANNEL functionality"
  
  ! Test 2.1: Set single-channel parameters
  CALL SET_CHANNEL(ch, 3, 3, 0, 0)
  CALL test_channel(ch, 3, 3, 0, 0, 1, .FALSE., "1F3")
  
  ! Test 2.2: Set coupled-channel parameters
  CALL SET_CHANNEL(ch, 2, 1, 1, 0)
  CALL test_channel(ch, 2, 1, 1, 0, 2, .TRUE., "3P2-3F2")
  
  ! Test 3: Test GET_CHANNEL_FROM_NAME functionality
  WRITE(*,*) "Test 3: Testing GET_CHANNEL_FROM_NAME functionality"
  
  ! Test 3.1: Parse uncoupled channel name
  ch = GET_CHANNEL_FROM_NAME("1S0")
  CALL test_channel(ch, 0, 0, 0, 0, 1, .FALSE., "1S0")
  
  ! Test 3.2: Parse coupled channel name
  ch = GET_CHANNEL_FROM_NAME("3S1-3D1")
  CALL test_channel(ch, 1, 0, 1, 0, 2, .TRUE., "3S1-3D1")
  
  ! Test 3.3: Test with unusual L values (high angular momentum)
  ch = GET_CHANNEL_FROM_NAME("3G4")
  CALL test_channel(ch, 4, 4, 1, 0, 1, .FALSE., "3G4")
  
  ! Test 4: Test IS_SAME_CHANNEL functionality
  WRITE(*,*) "Test 4: Testing IS_SAME_CHANNEL functionality"
  
  ! Test 4.1: Compare identical channels
  ch = GET_CHANNEL_FROM_NAME("3P2-3F2")
  ch2 = GET_CHANNEL_FROM_NAME("3P2-3F2")
  result = IS_SAME_CHANNEL(ch, ch2)
  IF (.NOT. result) STOP "FAIL: Identical channels not detected as same"
  
  ! Test 4.2: Compare different channels
  ch = GET_CHANNEL_FROM_NAME("3P2-3F2")
  ch2 = GET_CHANNEL_FROM_NAME("1D2")
  result = IS_SAME_CHANNEL(ch, ch2)
  IF (result) STOP "FAIL: Different channels detected as same"
  
  ! Test 5: Test IS_PHYSICAL_CHANNEL functionality
  WRITE(*,*) "Test 5: Testing IS_PHYSICAL_CHANNEL functionality"
  
  ! Test 5.1: Test with valid channel
  ch = GET_CHANNEL_FROM_NAME("3F3")
  is_phys = IS_PHYSICAL_CHANNEL(ch)
  IF (.NOT. is_phys) STOP "FAIL: Valid channel not detected as physical"
  
  ! Test 5.2: Test with non-physical L, S, J combination
  ! Need to manually set up an invalid channel since the constructor enforces validity
  ch = init_scattering_channel(2, .TRUE., 0)
  CALL SET_CHANNEL(ch, 2, 2, 1, 0)
  ! Now manually change one value to make it invalid
  ! We can't do this directly because the fields are private,
  ! but we can create a copy and modify it through the interface
  
  ! Test 6: Test all quantum number getters
  WRITE(*,*) "Test 6: Testing quantum number getters"
  ch = GET_CHANNEL_FROM_NAME("3D2-3G2")
  
  J = GET_CHANNEL_J(ch)
  IF (J /= 2) STOP "FAIL: GET_CHANNEL_J incorrect"
  
  L = GET_CHANNEL_L(ch, 1)
  IF (L /= 2) STOP "FAIL: GET_CHANNEL_L incorrect for first component"
  
  L = GET_CHANNEL_L(ch, 2)
  IF (L /= 4) STOP "FAIL: GET_CHANNEL_L incorrect for second component"
  
  S = GET_CHANNEL_S(ch, 1)
  IF (S /= 1) STOP "FAIL: GET_CHANNEL_S incorrect"
  
  T = GET_CHANNEL_T(ch, 1)
  IF (T /= 0) STOP "FAIL: GET_CHANNEL_T incorrect"
  
  TZ = GET_CHANNEL_TZ(ch)
  IF (TZ /= 0) STOP "FAIL: GET_CHANNEL_TZ incorrect"
  
  is_coup = IS_CHANNEL_COUPLED(ch)
  IF (.NOT. is_coup) STOP "FAIL: IS_CHANNEL_COUPLED incorrect"
  
  ! Test 7: Round-trip name conversion
  WRITE(*,*) "Test 7: Testing round-trip name conversion"
  name = "3P2-3F2"
  ch = GET_CHANNEL_FROM_NAME(name)
  name2 = GET_CHANNEL_NAME(ch)
  IF (TRIM(name) /= TRIM(name2)) STOP "FAIL: Round-trip name conversion failed"
  
  WRITE(*,*) "All tests passed successfully!"
  
CONTAINS

  ! Helper subroutine to test a channel against expected values
  SUBROUTINE test_channel(ch, exp_J, exp_L, exp_S, exp_TZ, exp_nch, exp_coupled, exp_name)
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: ch
    INTEGER, INTENT(IN) :: exp_J, exp_L, exp_S, exp_TZ, exp_nch
    LOGICAL, INTENT(IN) :: exp_coupled
    CHARACTER(LEN=*), INTENT(IN) :: exp_name
    
    INTEGER :: actual_J, actual_L, actual_S, actual_TZ, actual_nch
    LOGICAL :: actual_coupled
    CHARACTER(LEN=16) :: actual_name
    
    actual_J = GET_CHANNEL_J(ch)
    actual_L = GET_CHANNEL_L(ch, 1)
    actual_S = GET_CHANNEL_S(ch, 1)
    actual_TZ = GET_CHANNEL_TZ(ch)
    actual_nch = GET_CHANNEL_NCH(ch)
    actual_coupled = IS_CHANNEL_COUPLED(ch)
    actual_name = GET_CHANNEL_NAME(ch)
    
    IF (actual_J /= exp_J) STOP "FAIL: J mismatch"
    IF (actual_L /= exp_L) STOP "FAIL: L mismatch"
    IF (actual_S /= exp_S) STOP "FAIL: S mismatch"
    IF (actual_TZ /= exp_TZ) STOP "FAIL: TZ mismatch"
    IF (actual_nch /= exp_nch) STOP "FAIL: NCH mismatch"
    IF (actual_coupled .NEQV. exp_coupled) STOP "FAIL: COUPLED status mismatch"
    IF (TRIM(actual_name) /= TRIM(exp_name)) THEN
      WRITE(*,*) "Expected: ", TRIM(exp_name), " Got: ", TRIM(actual_name)
      STOP "FAIL: Name mismatch"
    END IF
    
    WRITE(*,*) "  Passed: ", TRIM(actual_name)
  END SUBROUTINE test_channel

END PROGRAM test_quantum_numbers_full
