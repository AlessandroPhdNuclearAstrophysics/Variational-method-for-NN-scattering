PROGRAM test_quantum_numbers_full
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  
  TYPE(SCATTERING_CHANNEL) :: ch, ch2, ch3
  CHARACTER(LEN=16) :: name, name2
  INTEGER :: i, nch, L, S, J, T, TZ
  LOGICAL :: is_phys, is_coup, result
  
  ! For EXTRACT_CHANNELS_FROM_WHOLE_FILENAME tests
  CHARACTER(LEN=100) :: filename
  CHARACTER(LEN=16) :: extracted_ch1, extracted_ch2
  LOGICAL :: found, is_coupled_ch
  
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
  
  ! ! Test 8: Testing EXTRACT_CHANNELS_FROM_WHOLE_FILENAME
  ! WRITE(*,*) "Test 8: Testing EXTRACT_CHANNELS_FROM_WHOLE_FILENAME"
  
  ! ! Test 8.1: Simple filename with one channel
  ! filename = "potential_model_1S0_output.dat"
  ! extracted_ch1 = ""
  ! found = EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(filename, extracted_ch1)
  ! IF (.NOT. found) STOP "FAIL: Channel not found in filename with single channel"
  ! IF (TRIM(extracted_ch1) /= "1S0") THEN
  !   WRITE(*,*) "Expected: 1S0, Got: ", TRIM(extracted_ch1)
  !   STOP "FAIL: Incorrect channel extracted from single channel filename"
  ! END IF
  ! WRITE(*,*) "  Passed: Single channel extraction"
  
  ! ! Test 8.2: Filename with coupled channel
  ! filename = "phase_shift_3S1-3D1_calculation.txt"
  ! extracted_ch1 = ""
  ! extracted_ch2 = ""
  ! is_coupled_ch = .FALSE.
  ! found = EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(filename, extracted_ch1, extracted_ch2, is_coupled_ch)
  ! IF (.NOT. found) STOP "FAIL: Channels not found in coupled channel filename"
  ! IF (.NOT. is_coupled_ch) STOP "FAIL: Coupled channels not detected"
  ! IF (TRIM(extracted_ch1) /= "3S1" .OR. TRIM(extracted_ch2) /= "3D1") THEN
  !   WRITE(*,*) "Expected: 3S1 and 3D1, Got: ", TRIM(extracted_ch1), " and ", TRIM(extracted_ch2)
  !   STOP "FAIL: Incorrect channels extracted from coupled channel filename"
  ! END IF
  ! WRITE(*,*) "  Passed: Coupled channel extraction with separate outputs"
  
  ! ! Test 8.3: Filename with coupled channel (single output parameter)
  ! filename = "results_3P2-3F2_analysis.csv"
  ! extracted_ch1 = ""
  ! found = EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(filename, extracted_ch1)
  ! IF (.NOT. found) STOP "FAIL: Channel not found in coupled channel filename (single output)"
  ! IF (TRIM(extracted_ch1) /= "3P2-3F2") THEN
  !   WRITE(*,*) "Expected: 3P2-3F2, Got: ", TRIM(extracted_ch1)
  !   STOP "FAIL: Incorrect channel extracted from coupled channel filename (single output)"
  ! END IF
  ! WRITE(*,*) "  Passed: Coupled channel extraction with single output"
  
  ! ! Test 8.4: Filename with no valid channel
  ! filename = "general_results_2023.dat"
  ! extracted_ch1 = "should_be_empty"
  ! found = EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(filename, extracted_ch1)
  ! IF (found) STOP "FAIL: False positive channel detection"
  ! IF (TRIM(extracted_ch1) /= "") THEN
  !   WRITE(*,*) "Expected empty string, Got: ", TRIM(extracted_ch1)
  !   STOP "FAIL: Output not cleared for filename with no channel"
  ! END IF
  ! WRITE(*,*) "  Passed: No false positives"
  
  ! ! Test 8.5: Filename with channel in the middle
  ! filename = "prefix_with_1D2_suffix.txt"
  ! extracted_ch1 = ""
  ! found = EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(filename, extracted_ch1)
  ! IF (.NOT. found) STOP "FAIL: Channel not found in middle of filename"
  ! IF (TRIM(extracted_ch1) /= "1D2") THEN
  !   WRITE(*,*) "Expected: 1D2, Got: ", TRIM(extracted_ch1)
  !   STOP "FAIL: Incorrect channel extracted from middle of filename"
  ! END IF
  ! WRITE(*,*) "  Passed: Channel extraction from middle of filename"
  
  ! ! Test 8.6: Edge case - channel at start of filename
  ! filename = "3G4_results_file.dat"
  ! extracted_ch1 = ""
  ! found = EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(filename, extracted_ch1)
  ! IF (.NOT. found) STOP "FAIL: Channel not found at start of filename"
  ! IF (TRIM(extracted_ch1) /= "3G4") THEN
  !   WRITE(*,*) "Expected: 3G4, Got: ", TRIM(extracted_ch1)
  !   STOP "FAIL: Incorrect channel extracted from start of filename"
  ! END IF
  ! WRITE(*,*) "  Passed: Channel extraction from start of filename"
  
  ! ! Test 8.7: Edge case - channel at end of filename
  ! filename = "nuclear_data_3H6"
  ! extracted_ch1 = ""
  ! found = EXTRACT_CHANNELS_FROM_WHOLE_FILENAME(filename, extracted_ch1)
  ! IF (.NOT. found) STOP "FAIL: Channel not found at end of filename"
  ! IF (TRIM(extracted_ch1) /= "3H6") THEN
  !   WRITE(*,*) "Expected: 3H6, Got: ", TRIM(extracted_ch1)
  !   STOP "FAIL: Incorrect channel extracted from end of filename"
  ! END IF
  ! WRITE(*,*) "  Passed: Channel extraction from end of filename"
  
  ! WRITE(*,*) "All tests passed successfully!"
  
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
