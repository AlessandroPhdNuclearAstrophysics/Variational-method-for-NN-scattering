PROGRAM test_quantum_numbers
  USE QUANTUM_NUMBERS
  IMPLICIT NONE

  TYPE(SCATTERING_CHANNEL) :: ch, ch2
  CHARACTER(LEN=16) :: name
  INTEGER :: i, nch
  LOGICAL :: ok

  PRINT *, "Testing init_scattering_channel (J=0, even, TZ=0)..."
  ch = init_scattering_channel(0, .TRUE., 0)
  IF (GET_CHANNEL_NCH(ch) /= 1) STOP "FAIL: NCH /= 1 for J=0"
  IF (GET_CHANNEL_L(ch,1) /= 0 .OR. GET_CHANNEL_S(ch,1) /= 0) STOP "FAIL: L/S wrong for J=0, even"
  IF (.NOT. IS_PHYSICAL_CHANNEL(ch)) STOP "FAIL: Channel not physical"

  PRINT *, "Testing init_scattering_channel (J=1, odd, TZ=1)..."
  ch = init_scattering_channel(1, .FALSE., 1)
  IF (GET_CHANNEL_NCH(ch) /= 2) STOP "FAIL: NCH /= 2 for J=1, odd"
  IF (.NOT. IS_CHANNEL_COUPLED(ch)) STOP "FAIL: Channel not coupled for J=1, odd"
  IF (.NOT. IS_PHYSICAL_CHANNEL(ch)) STOP "FAIL: Channel not physical"

  PRINT *, "Testing SET_CHANNEL..."
  CALL SET_CHANNEL(ch, 2, 2, 1, 0)
  IF (GET_CHANNEL_J(ch) /= 2 .OR. GET_CHANNEL_L(ch,1) /= 2 .OR. GET_CHANNEL_S(ch,1) /= 1) STOP "FAIL: SET_CHANNEL wrong"
  IF (.NOT. IS_PHYSICAL_CHANNEL(ch)) STOP "FAIL: Channel not physical after SET_CHANNEL"

  PRINT *, "Testing GET_CHANNEL_NAME..."
  name = GET_CHANNEL_NAME(ch)
  IF (TRIM(name) /= '3D2') STOP "FAIL: GET_CHANNEL_NAME wrong, got "//TRIM(name)

  PRINT *, "Testing GET_CHANNEL_FROM_NAME and IS_SAME_CHANNEL..."
  ch2 = GET_CHANNEL_FROM_NAME('3D2')
  IF (.NOT. IS_SAME_CHANNEL(ch, ch2)) STOP "FAIL: IS_SAME_CHANNEL failed for 3D2"

  PRINT *, "Testing coupled channel name and parsing..."
  ch = GET_CHANNEL_FROM_NAME('3P1-1P1')
  name = GET_CHANNEL_NAME(ch)
  IF (TRIM(name) /= '3P1-1P1') STOP "FAIL: Coupled channel name wrong, got "//TRIM(name)
  IF (.NOT. IS_CHANNEL_COUPLED(ch)) STOP "FAIL: Coupled channel not detected"

  PRINT *, "Testing quantum number getters..."
  nch = GET_CHANNEL_NCH(ch)
  IF (nch /= 2) STOP "FAIL: GET_CHANNEL_NCH wrong"
  IF (GET_CHANNEL_L(ch,1) /= 1 .OR. GET_CHANNEL_S(ch,1) /= 1) STOP "FAIL: GET_CHANNEL_L/S wrong"
  IF (GET_CHANNEL_L(ch,2) /= 1 .OR. GET_CHANNEL_S(ch,2) /= 0) STOP "FAIL: GET_CHANNEL_L/S wrong (2nd)"

  PRINT *, "All tests passed."

END PROGRAM test_quantum_numbers