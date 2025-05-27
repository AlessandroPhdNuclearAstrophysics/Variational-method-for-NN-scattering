MODULE FIT_MODULE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LINEAR_REGRESSION, QUADRATIC_REGRESSION


CONTAINS
  SUBROUTINE LINEAR_REGRESSION(Y, X, NMIN, NMAX, M, Q, STEP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NMIN, NMAX
    DOUBLE PRECISION, INTENT(IN) :: Y(NMAX), X(NMAX)
    DOUBLE PRECISION, INTENT(OUT) :: M, Q
    INTEGER, OPTIONAL, INTENT(IN) :: STEP
    INTEGER :: I, N, STEP_I
    DOUBLE PRECISION :: SX, SX2, SY, SXY, DENOM

    IF (PRESENT(STEP)) THEN
      STEP_I = STEP
    ELSE
      STEP_I = 1
    END IF

    ! Compute number of points
    N = NMAX - NMIN + 1
    IF (N < 2*STEP_I) THEN
      PRINT *, "ERROR: NOT ENOUGH POINTS FOR LINEAR_REGRESSION"
      STOP
    END IF

    ! Initialize sums
    SX  = 0.0D0
    SX2 = 0.0D0
    SY  = 0.0D0
    SXY = 0.0D0

    ! Compute sums
    N = 0
    DO I = NMIN, NMAX, STEP_I
      SX  = SX  + X(I)
      SX2 = SX2 + X(I) * X(I)
      SY  = SY  + Y(I)
      SXY = SXY + X(I) * Y(I)
      N = N + 1
    END DO

    ! Compute denominator (check for zero)
    DENOM = N * SX2 - SX * SX
    IF (ABS(DENOM) < 1.0D-12) THEN
      PRINT *, "ERROR: LINEAR REGRESSION SINGULARITY (DIVISION BY ZERO)"
      STOP
    END IF

    ! Compute slope and intercept
    M = (N * SXY - SX * SY) / DENOM
    Q = (SY - M * SX) / N

    RETURN
  END SUBROUTINE LINEAR_REGRESSION



  SUBROUTINE QUADRATIC_REGRESSION(Y, X, NMIN, NMAX, A, B, C, STEP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NMIN, NMAX
    DOUBLE PRECISION, INTENT(IN) :: Y(NMAX), X(NMAX)
    DOUBLE PRECISION, INTENT(OUT) :: A, B, C
    INTEGER, OPTIONAL, INTENT(IN) :: STEP
    INTEGER :: I, N, STEP_I
    DOUBLE PRECISION :: S4, S3, S2, S1, S2Y, S1Y, S0Y, DET

    IF (PRESENT(STEP)) THEN
      STEP_I = STEP
    ELSE
      STEP_I = 1
    END IF

    ! Compute number of points
    N = NMAX - NMIN + 1
    IF (N < 3*STEP_I) THEN
      PRINT *, "ERROR: NOT ENOUGH POINTS FOR LINEAR_REGRESSION"
      STOP
    END IF

    ! Initialize sums
    S4  = 0.D0  
    S3  = 0.D0
    S2  = 0.D0
    S1  = 0.D0
    S2Y = 0.D0
    S1Y = 0.D0
    S0Y = 0.D0

    ! Compute sums
    N = 0
    DO I = NMIN, NMAX, STEP_I
      S4  = S4  + X(I)**4
      S3  = S3  + X(I)**3
      S2  = S2  + X(I)**2
      S1  = S1  + X(I)
      S2Y = S2Y + X(I)**2*Y(I)
      S1Y = S1Y + X(I)*Y(I)
      S0Y = S0Y + Y(I)
      N = N + 1
    END DO

    DET = -S2**3 + 2*S1*S2*S3 - N*S3**2 - S1**2*S4 + N*S2*S4
    IF (ABS(DET).LT.1.D-15) STOP "DETERMINANT IS ZERO"

    A = S1*S1Y*S2 - S0Y*S2**2 - S1**2*S2Y + N*S2*S2Y + S0Y*S1*S3 - N*S1Y*S3
    B = -S1Y*S2**2 + S1*S2*S2Y + S0Y*S2*S3 - N*S2Y*S3 - S0Y*S1*S4 + N*S1Y*S4
    C = -S2**2*S2Y + S1Y*S2*S3 + S1*S2Y*S3 - S0Y*S3**2 - S1*S1Y*S4 + S0Y*S2*S4
    
    A = A/DET
    B = B/DET
    C = C/DET

    RETURN
  END SUBROUTINE QUADRATIC_REGRESSION




END MODULE FIT_MODULE