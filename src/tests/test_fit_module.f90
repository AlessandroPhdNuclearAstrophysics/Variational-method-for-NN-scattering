PROGRAM test_fit_module
  USE FIT_MODULE
  IMPLICIT NONE

  ! Arrays for test data
  INTEGER, PARAMETER :: N = 100
  DOUBLE PRECISION :: X(N), Y_LINEAR(N), Y_QUADRATIC(N), Y_CUBIC(N)

  ! Expected coefficients (ground truth)
  DOUBLE PRECISION :: EXP_M = 2.5D0, EXP_Q = -1.2D0
  DOUBLE PRECISION :: EXP_A = 0.5D0, EXP_B = -1.5D0, EXP_C = 2.0D0
  DOUBLE PRECISION :: EXP_A3 = 0.3D0, EXP_B3 = -0.8D0, EXP_C3 = 1.2D0, EXP_D3 = 0.5D0

  ! Fitted coefficients
  DOUBLE PRECISION :: FIT_M, FIT_Q
  DOUBLE PRECISION :: FIT_A, FIT_B, FIT_C
  DOUBLE PRECISION :: FIT_A3, FIT_B3, FIT_C3, FIT_D3

  ! Tolerances for comparisons
  DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-3  ! Changed from 1.0D-10 to be compatible with noise level

  ! Loop variables
  INTEGER :: STEP

  WRITE(*,*) "=== Testing FIT_MODULE ==="

  ! Generate test data
  CALL GENERATE_TEST_DATA()

  ! Test 1: Basic linear regression
  WRITE(*,*) "Test 1: Basic linear regression"
  CALL LINEAR_REGRESSION(Y_LINEAR, X, 1, N, FIT_M, FIT_Q)
  CALL CHECK_LINEAR_RESULTS(FIT_M, FIT_Q)

  ! Test 2: Linear regression with step parameter
  WRITE(*,*) "Test 2: Linear regression with step parameter"
  STEP = 2
  CALL LINEAR_REGRESSION(Y_LINEAR, X, 1, N, FIT_M, FIT_Q, STEP)
  CALL CHECK_LINEAR_RESULTS(FIT_M, FIT_Q)

  ! Test 3: Linear regression with subset of data
  WRITE(*,*) "Test 3: Linear regression with subset of data"
  CALL LINEAR_REGRESSION(Y_LINEAR, X, N/2, N, FIT_M, FIT_Q)
  CALL CHECK_LINEAR_RESULTS(FIT_M, FIT_Q)

  ! Test 4: Basic quadratic regression
  WRITE(*,*) "Test 4: Basic quadratic regression"
  CALL QUADRATIC_REGRESSION(Y_QUADRATIC, X, 1, N, FIT_A, FIT_B, FIT_C)
  CALL CHECK_QUADRATIC_RESULTS(FIT_A, FIT_B, FIT_C)

  ! Test 5: Quadratic regression with step parameter
  WRITE(*,*) "Test 5: Quadratic regression with step parameter"
  STEP = 2
  CALL QUADRATIC_REGRESSION(Y_QUADRATIC, X, 1, N, FIT_A, FIT_B, FIT_C, STEP)
  CALL CHECK_QUADRATIC_RESULTS(FIT_A, FIT_B, FIT_C)

  ! Test 6: Quadratic regression with subset of data
  WRITE(*,*) "Test 6: Quadratic regression with subset of data"
  CALL QUADRATIC_REGRESSION(Y_QUADRATIC, X, N/2, N, FIT_A, FIT_B, FIT_C)
  CALL CHECK_QUADRATIC_RESULTS(FIT_A, FIT_B, FIT_C)

  ! Test 7: Basic cubic regression
  WRITE(*,*) "Test 7: Basic cubic regression"
  CALL CUBIC_REGRESSION(Y_CUBIC, X, 1, N, FIT_A3, FIT_B3, FIT_C3, FIT_D3)
  CALL CHECK_CUBIC_RESULTS(FIT_A3, FIT_B3, FIT_C3, FIT_D3)

  ! Test 8: Cubic regression with step parameter
  WRITE(*,*) "Test 8: Cubic regression with step parameter"
  STEP = 2
  CALL CUBIC_REGRESSION(Y_CUBIC, X, 1, N, FIT_A3, FIT_B3, FIT_C3, FIT_D3, STEP)
  CALL CHECK_CUBIC_RESULTS(FIT_A3, FIT_B3, FIT_C3, FIT_D3)

  ! Test 9: Cubic regression with subset of data
  WRITE(*,*) "Test 9: Cubic regression with subset of data"
  CALL CUBIC_REGRESSION(Y_CUBIC, X, N/2, N, FIT_A3, FIT_B3, FIT_C3, FIT_D3)
  CALL CHECK_CUBIC_RESULTS(FIT_A3, FIT_B3, FIT_C3, FIT_D3)

  WRITE(*,*) "All tests passed successfully!"

CONTAINS

  ! Generate test data with known coefficients
  SUBROUTINE GENERATE_TEST_DATA()
    IMPLICIT NONE
    INTEGER :: I
    REAL :: R

    ! Initialize random number generator for reproducibility
    CALL RANDOM_SEED()

    ! Generate x values from -5 to 5
    DO I = 1, N
      X(I) = -5.0D0 + 10.0D0 * (I - 1.0D0) / (N - 1.0D0)

      ! Linear function: y = m*x + q
      Y_LINEAR(I) = EXP_M * X(I) + EXP_Q

      ! Quadratic function: y = a*x^2 + b*x + c
      Y_QUADRATIC(I) = EXP_A * X(I)**2 + EXP_B * X(I) + EXP_C

      ! Cubic function: y = a*x^3 + b*x^2 + c*x + d
      Y_CUBIC(I) = EXP_A3 * X(I)**3 + EXP_B3 * X(I)**2 + EXP_C3 * X(I) + EXP_D3
    END DO

    ! Add small random noise to make it more realistic
    DO I = 1, N
      CALL RANDOM_NUMBER(R)
      Y_LINEAR(I) = Y_LINEAR(I) + (2.0D0 * R - 1.0D0) * 0.001D0

      CALL RANDOM_NUMBER(R)
      Y_QUADRATIC(I) = Y_QUADRATIC(I) + (2.0D0 * R - 1.0D0) * 0.001D0

      CALL RANDOM_NUMBER(R)
      Y_CUBIC(I) = Y_CUBIC(I) + (2.0D0 * R - 1.0D0) * 0.001D0
    END DO
  END SUBROUTINE GENERATE_TEST_DATA

  ! Check if linear regression results are correct
  SUBROUTINE CHECK_LINEAR_RESULTS(M, Q)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: M, Q

    WRITE(*,*) "  Expected: M =", EXP_M, "Q =", EXP_Q
    WRITE(*,*) "  Fitted:   M =", M, "Q =", Q

    IF (ABS(M - EXP_M) > TOL) THEN
      WRITE(*,*) "FAIL: Slope (M) does not match expected value"
      WRITE(*,*) "Expected:", EXP_M, "Got:", M, "Diff:", ABS(M - EXP_M)
      STOP
    END IF

    IF (ABS(Q - EXP_Q) > TOL) THEN
      WRITE(*,*) "FAIL: Intercept (Q) does not match expected value"
      WRITE(*,*) "Expected:", EXP_Q, "Got:", Q, "Diff:", ABS(Q - EXP_Q)
      STOP
    END IF

    WRITE(*,*) "  PASSED"
  END SUBROUTINE CHECK_LINEAR_RESULTS

  ! Check if quadratic regression results are correct
  SUBROUTINE CHECK_QUADRATIC_RESULTS(A, B, C)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: A, B, C

    WRITE(*,*) "  Expected: A =", EXP_A, "B =", EXP_B, "C =", EXP_C
    WRITE(*,*) "  Fitted:   A =", A, "B =", B, "C =", C

    IF (ABS(A - EXP_A) > TOL) THEN
      WRITE(*,*) "FAIL: Quadratic coefficient (A) does not match expected value"
      WRITE(*,*) "Expected:", EXP_A, "Got:", A, "Diff:", ABS(A - EXP_A)
      STOP
    END IF

    IF (ABS(B - EXP_B) > TOL) THEN
      WRITE(*,*) "FAIL: Linear coefficient (B) does not match expected value"
      WRITE(*,*) "Expected:", EXP_B, "Got:", B, "Diff:", ABS(B - EXP_B)
      STOP
    END IF

    IF (ABS(C - EXP_C) > TOL) THEN
      WRITE(*,*) "FAIL: Constant term (C) does not match expected value"
      WRITE(*,*) "Expected:", EXP_C, "Got:", C, "Diff:", ABS(C - EXP_C)
      STOP
    END IF

    WRITE(*,*) "  PASSED"
  END SUBROUTINE CHECK_QUADRATIC_RESULTS

  ! Check if cubic regression results are correct
  SUBROUTINE CHECK_CUBIC_RESULTS(A, B, C, D)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: A, B, C, D

    WRITE(*,*) "  Expected: A =", EXP_A3, "B =", EXP_B3, "C =", EXP_C3, "D =", EXP_D3
    WRITE(*,*) "  Fitted:   A =", A, "B =", B, "C =", C, "D =", D

    IF (ABS(A - EXP_A3) > TOL) THEN
      WRITE(*,*) "FAIL: Cubic coefficient (A) does not match expected value"
      WRITE(*,*) "Expected:", EXP_A3, "Got:", A, "Diff:", ABS(A - EXP_A3)
      STOP
    END IF

    IF (ABS(B - EXP_B3) > TOL) THEN
      WRITE(*,*) "FAIL: Quadratic coefficient (B) does not match expected value"
      WRITE(*,*) "Expected:", EXP_B3, "Got:", B, "Diff:", ABS(B - EXP_B3)
      STOP
    END IF

    IF (ABS(C - EXP_C3) > TOL) THEN
      WRITE(*,*) "FAIL: Linear coefficient (C) does not match expected value"
      WRITE(*,*) "Expected:", EXP_C3, "Got:", C, "Diff:", ABS(C - EXP_C3)
      STOP
    END IF

    IF (ABS(D - EXP_D3) > TOL) THEN
      WRITE(*,*) "FAIL: Constant term (D) does not match expected value"
      WRITE(*,*) "Expected:", EXP_D3, "Got:", D, "Diff:", ABS(D - EXP_D3)
      STOP
    END IF

    WRITE(*,*) "  PASSED"
  END SUBROUTINE CHECK_CUBIC_RESULTS

END PROGRAM test_fit_module
