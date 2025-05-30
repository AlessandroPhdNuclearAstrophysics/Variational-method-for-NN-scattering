!> \file test_quantum_operators.f90
!! \brief Unit tests for the QUANTUM_OPERATORS module

PROGRAM TEST_QUANTUM_OPERATORS
  USE QUANTUM_OPERATORS
  USE QUANTUM_NUMBERS
  IMPLICIT NONE

  INTEGER :: i, j
  DOUBLE PRECISION, ALLOCATABLE :: id2(:,:)
  DOUBLE PRECISION :: s12_mat(2,2), s12_chan(2,2)
  INTEGER :: ls_mat(2,2), ls_chan(2,2)
  INTEGER :: l2_mat(2,2), l2_chan(2,2)
  INTEGER :: t12_val, t12_chan
  TYPE(SCATTERING_CHANNEL) :: ch

  CALL SET_CHANNEL(ch, 1, 0, 1, 0)

  PRINT *, "Testing IDENTITY_MATRIX..."
  id2 = IDENTITY_MATRIX()
  IF (ALL(id2 == RESHAPE([1.0D0,0.0D0,0.0D0,1.0D0],[2,2]))) THEN
    PRINT *, "  IDENTITY_MATRIX (default 2x2) OK"
  ELSE
    PRINT *, "  IDENTITY_MATRIX (default 2x2) FAILED"
  END IF

  id2 = IDENTITY_MATRIX(3)
  IF (ALL(id2 == RESHAPE([1.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,1.0D0],[3,3]))) THEN
    PRINT *, "  IDENTITY_MATRIX (3x3) OK"
  ELSE
    PRINT *, "  IDENTITY_MATRIX (3x3) FAILED"
  END IF

  PRINT *, "Testing S12_OPERATOR..."
  s12_mat = S12_OPERATOR(1,1,1)
  IF (ABS(s12_mat(1,1)-2.0D0)<1E-12 .AND. ALL(s12_mat(2,:) == 0.0D0) .AND. ALL(s12_mat(:,2) == 0.0D0)) THEN
    PRINT *, "  S12_OPERATOR uncoupled OK"
  ELSE
    PRINT *, "  S12_OPERATOR uncoupled FAILED"
  END IF

  s12_mat = S12_OPERATOR(1,1,2)
  IF (ABS(s12_mat(1,1)+2.0D0/5.0D0)<1E-12 .AND. ABS(s12_mat(2,2)+8.0D0/5.0D0)<1E-12) THEN
    PRINT *, "  S12_OPERATOR coupled OK"
  ELSE
    PRINT *, "  S12_OPERATOR coupled FAILED"
  END IF

  ! Test S12_OPERATOR
  CALL SET_CHANNEL(ch, 2, 1, 1, 0)  ! Set channel to L=1, S=1, J=2
  s12_chan = S12_OPERATOR(ch)
  s12_mat = S12_OPERATOR(1,1,2)
  IF (ALL(ABS(s12_chan-s12_mat)<1E-12)) THEN
    PRINT *, "  S12_OPERATOR matches quantum numbers version OK"
  ELSE
    PRINT *, "  S12_OPERATOR matches quantum numbers version FAILED"
  END IF

  PRINT *, "Testing LS_OPERATOR..."
  ls_mat = LS_OPERATOR(1,1,1)
  IF (ls_mat(1,1) == -1) THEN
    PRINT *, "  LS_OPERATOR uncoupled OK"
  ELSE
    PRINT *, "  LS_OPERATOR uncoupled FAILED"
  END IF

  ls_mat = LS_OPERATOR(1,1,2)
  IF (ls_mat(1,1) == 1 .AND. ls_mat(2,2) == -4) THEN
    PRINT *, "  LS_OPERATOR coupled OK"
  ELSE
    PRINT *, "  LS_OPERATOR coupled FAILED"
  END IF

  ! Test LS_OPERATOR
  CALL SET_CHANNEL(ch, 2, 1, 1, 0)
  ls_chan = LS_OPERATOR(ch)
  ls_mat = LS_OPERATOR(1,1,2)
  IF (ALL(ls_chan == ls_mat)) THEN
    PRINT *, "  LS_OPERATOR matches quantum numbers version OK"
  ELSE
    PRINT *, "  LS_OPERATOR matches quantum numbers version FAILED"
  END IF

  PRINT *, "Testing L2_OPERATOR..."
  l2_mat = L2_OPERATOR(1)
  IF (l2_mat(1,1) == 2 .AND. l2_mat(2,2) == 0) THEN
    PRINT *, "  L2_OPERATOR single OK"
  ELSE
    PRINT *, "  L2_OPERATOR single FAILED"
  END IF

  l2_mat = L2_OPERATOR(1, .TRUE.)
  IF (l2_mat(1,1) == 2 .AND. l2_mat(2,2) == 12) THEN
    PRINT *, "  L2_OPERATOR coupled OK"
  ELSE
    PRINT *, "  L2_OPERATOR coupled FAILED"
  END IF

  ! Test L2_OPERATOR
  CALL SET_CHANNEL(ch, 1, 1, 0, 0)
  l2_chan = L2_OPERATOR(ch)
  l2_mat = L2_OPERATOR(1)
  IF (ALL(l2_chan == l2_mat)) THEN
    PRINT *, "  L2_OPERATOR matches quantum numbers version (single) OK"
  ELSE
    PRINT *, "  L2_OPERATOR matches quantum numbers version (single) FAILED"
  END IF

  CALL SET_CHANNEL(ch, 2, 1, 1, 0)
  l2_chan = L2_OPERATOR(ch)
  l2_mat = L2_OPERATOR(1, .TRUE.)
  IF (ALL(l2_chan == l2_mat)) THEN
    PRINT *, "  L2_OPERATOR matches quantum numbers version (coupled) OK"
  ELSE
    PRINT *, "  L2_OPERATOR matches quantum numbers version (coupled) FAILED"
  END IF

  PRINT *, "Testing T12_OPERATOR..."
  t12_val = T12_OPERATOR(1,1)
  IF (t12_val == 2) THEN
    PRINT *, "  T12_OPERATOR Tz=1 OK"
  ELSE
    PRINT *, "  T12_OPERATOR Tz=1 FAILED"
  END IF

  t12_val = T12_OPERATOR(1,0)
  IF (t12_val == -4) THEN
    PRINT *, "  T12_OPERATOR Tz=0 OK"
  ELSE
    PRINT *, "  T12_OPERATOR Tz=0 FAILED"
  END IF

  ! Test T12_OPERATOR
  CALL SET_CHANNEL(ch, 0, 0, 0, 1)
  t12_chan = T12_OPERATOR(ch)
  t12_val = T12_OPERATOR(1,1)
  IF (t12_chan == t12_val) THEN
    PRINT *, "  T12_OPERATOR matches quantum numbers version (Tz=1) OK"
  ELSE
    PRINT *, "  T12_OPERATOR matches quantum numbers version (Tz=1) FAILED"
  END IF

  CALL SET_CHANNEL(ch, 1, 0, 0, 0)
  t12_chan = T12_OPERATOR(ch)
  t12_val = T12_OPERATOR(1,0)
  IF (t12_chan == t12_val) THEN
    PRINT *, "  T12_OPERATOR matches quantum numbers version (Tz=0) OK"
  ELSE
    PRINT *, "  T12_OPERATOR matches quantum numbers version (Tz=0) FAILED"
  END IF

  ! Additional edge case tests
  PRINT *, "Testing S12_OPERATOR with S=0 (should be zero)..."
  s12_mat = S12_OPERATOR(1,0,1)
  IF (ALL(s12_mat == 0.0D0)) THEN
    PRINT *, "  S12_OPERATOR S=0 OK"
  ELSE
    PRINT *, "  S12_OPERATOR S=0 FAILED"
  END IF

  PRINT *, "Testing T12_OPERATOR with T=0 (should be zero)..."
  t12_val = T12_OPERATOR(0,1)
  IF (t12_val == 0) THEN
    PRINT *, "  T12_OPERATOR T=0 OK"
  ELSE
    PRINT *, "  T12_OPERATOR T=0 FAILED"
  END IF

  PRINT *, "All extended tests completed."

END PROGRAM TEST_QUANTUM_OPERATORS