MODULE FIT_MODULE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LINEAR_REGRESSION, QUADRATIC_REGRESSION, CUBIC_REGRESSION, POLYNOMIAL_REGRESSION


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

  SUBROUTINE CUBIC_REGRESSION(Y, X, NMIN, NMAX, A, B, C, D, STEP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NMIN, NMAX
    DOUBLE PRECISION, INTENT(IN) :: Y(NMAX), X(NMAX)
    DOUBLE PRECISION, INTENT(OUT) :: A, B, C, D
    INTEGER, OPTIONAL, INTENT(IN) :: STEP
    INTEGER :: I, N, STEP_I
    DOUBLE PRECISION :: S6, S5, S4, S3, S2, S1, S0
    DOUBLE PRECISION :: S3Y, S2Y, S1Y, S0Y
    DOUBLE PRECISION :: M(4,4), RHS(4), DET
    DOUBLE PRECISION :: M11, M12, M13, M14, M21, M22, M23, M24
    DOUBLE PRECISION :: M31, M32, M33, M34, M41, M42, M43, M44
    DOUBLE PRECISION :: V1, V2, V3, V4

    IF (PRESENT(STEP)) THEN
      STEP_I = STEP
    ELSE
      STEP_I = 1
    END IF

    ! Compute number of points
    N = NMAX - NMIN + 1
    IF (N < 4*STEP_I) THEN
      PRINT *, "ERROR: NOT ENOUGH POINTS FOR CUBIC_REGRESSION"
      STOP
    END IF

    ! Initialize sums
    S6  = 0.D0
    S5  = 0.D0  
    S4  = 0.D0
    S3  = 0.D0
    S2  = 0.D0
    S1  = 0.D0
    S0  = 0.D0
    S3Y = 0.D0
    S2Y = 0.D0
    S1Y = 0.D0
    S0Y = 0.D0

    ! Compute sums
    N = 0
    DO I = NMIN, NMAX, STEP_I
      S6  = S6  + X(I)**6
      S5  = S5  + X(I)**5
      S4  = S4  + X(I)**4
      S3  = S3  + X(I)**3
      S2  = S2  + X(I)**2
      S1  = S1  + X(I)
      S3Y = S3Y + X(I)**3*Y(I)
      S2Y = S2Y + X(I)**2*Y(I)
      S1Y = S1Y + X(I)*Y(I)
      S0Y = S0Y + Y(I)
      N = N + 1
    END DO
    S0 = N

    ! Setup the normal equations matrix
    ! [S6 S5 S4 S3] [A]   [S3Y]
    ! [S5 S4 S3 S2] [B] = [S2Y]
    ! [S4 S3 S2 S1] [C]   [S1Y]
    ! [S3 S2 S1 S0] [D]   [S0Y]
    
    M11 = S6; M12 = S5; M13 = S4; M14 = S3
    M21 = S5; M22 = S4; M23 = S3; M24 = S2
    M31 = S4; M32 = S3; M33 = S2; M34 = S1
    M41 = S3; M42 = S2; M43 = S1; M44 = S0
    
    V1 = S3Y
    V2 = S2Y
    V3 = S1Y
    V4 = S0Y
    
    ! Calculate determinant of the coefficient matrix
    DET = M11*(M22*(M33*M44 - M34*M43) - M23*(M32*M44 - M34*M42) + M24*(M32*M43 - M33*M42)) - &
          M12*(M21*(M33*M44 - M34*M43) - M23*(M31*M44 - M34*M41) + M24*(M31*M43 - M33*M41)) + &
          M13*(M21*(M32*M44 - M34*M42) - M22*(M31*M44 - M34*M41) + M24*(M31*M42 - M32*M41)) - &
          M14*(M21*(M32*M43 - M33*M42) - M22*(M31*M43 - M33*M41) + M23*(M31*M42 - M32*M41))
    
    IF (ABS(DET) < 1.D-15) THEN
      PRINT *, "ERROR: DETERMINANT IS ZERO IN CUBIC_REGRESSION"
      STOP
    END IF
    
    ! Solve for A using Cramer's rule
    A = (V1*(M22*(M33*M44 - M34*M43) - M23*(M32*M44 - M34*M42) + M24*(M32*M43 - M33*M42)) - &
         M12*(V2*(M33*M44 - M34*M43) - M23*(V3*M44 - M34*V4) + M24*(V3*M43 - M33*V4)) + &
         M13*(V2*(M32*M44 - M34*M42) - M22*(V3*M44 - M34*V4) + M24*(V3*M42 - M32*V4)) - &
         M14*(V2*(M32*M43 - M33*M42) - M22*(V3*M43 - M33*V4) + M23*(V3*M42 - M32*V4))) / DET
    
    ! Solve for B using Cramer's rule
    B = (M11*(V2*(M33*M44 - M34*M43) - M23*(V3*M44 - M34*V4) + M24*(V3*M43 - M33*V4)) - &
         V1*(M21*(M33*M44 - M34*M43) - M23*(M31*M44 - M34*M41) + M24*(M31*M43 - M33*M41)) + &
         M13*(M21*(V3*M44 - M34*V4) - V2*(M31*M44 - M34*M41) + M24*(M31*V4 - V3*M41)) - &
         M14*(M21*(V3*M43 - M33*V4) - V2*(M31*M43 - M33*M41) + M23*(M31*V4 - V3*M41))) / DET
    
    ! Solve for C using Cramer's rule
    C = (M11*(M22*(V3*M44 - M34*V4) - V2*(M32*M44 - M34*M42) + M24*(M32*V4 - V3*M42)) - &
         M12*(M21*(V3*M44 - M34*V4) - V2*(M31*M44 - M34*M41) + M24*(M31*V4 - V3*M41)) + &
         V1*(M21*(M32*M44 - M34*M42) - M22*(M31*M44 - M34*M41) + M24*(M31*M42 - M32*M41)) - &
         M14*(M21*(M32*V4 - V3*M42) - M22*(M31*V4 - V3*M41) + V2*(M31*M42 - M32*M41))) / DET
    
    ! Solve for D using Cramer's rule
    D = (M11*(M22*(M33*V4 - V3*M43) - M23*(M32*V4 - V3*M42) + V2*(M32*M43 - M33*M42)) - &
         M12*(M21*(M33*V4 - V3*M43) - M23*(M31*V4 - V3*M41) + V2*(M31*M43 - M33*M41)) + &
         M13*(M21*(M32*V4 - V3*M42) - M22*(M31*V4 - V3*M41) + V2*(M31*M42 - M32*M41)) - &
         V1*(M21*(M32*M43 - M33*M42) - M22*(M31*M43 - M33*M41) + M23*(M31*M42 - M32*M41))) / DET

    RETURN
  END SUBROUTINE CUBIC_REGRESSION


  SUBROUTINE POLYNOMIAL_REGRESSION(y, x, degree, n_points, coeffs)
    ! Performs polynomial regression of specified degree using LAPACK
    ! y = coeffs(1) + coeffs(2)*x + coeffs(3)*x^2 + ... + coeffs(degree+1)*x^degree
    
    INTEGER, INTENT(IN) :: degree, n_points
    DOUBLE PRECISION, INTENT(IN) :: x(n_points), y(n_points)
    DOUBLE PRECISION, INTENT(OUT) :: coeffs(degree+1)
    
    ! Local variables
    DOUBLE PRECISION, ALLOCATABLE :: work(:), A(:,:), b(:)
    INTEGER :: i, j, lwork, info, lda, ldb, nrhs
    CHARACTER(1) :: trans
    
    ! Set up the Vandermonde matrix A
    ALLOCATE(A(n_points, degree+1), b(n_points))
    
    ! Fill the Vandermonde matrix
    DO i = 1, n_points
      A(i, 1) = 1.0d0  ! Constant term
      DO j = 2, degree+1
        A(i, j) = A(i, j-1) * x(i)  ! x^(j-1)
      END DO
    END DO
    
    ! Copy y values to b (DGELS modifies the right-hand side)
    b = y
    
    ! Prepare for DGELS call
    lda = n_points
    ldb = n_points
    nrhs = 1
    trans = 'N'  ! No transpose
    
    ! Query optimal workspace size
    ALLOCATE(work(1))
    CALL DGELS(trans, n_points, degree+1, nrhs, A, lda, b, ldb, work, -1, info)
    lwork = INT(work(1))
    DEALLOCATE(work)
    ALLOCATE(work(lwork))
    
    ! Solve the least squares problem
    CALL DGELS(trans, n_points, degree+1, nrhs, A, lda, b, ldb, work, lwork, info)
    
    IF (info /= 0) THEN
      WRITE(*,*) "DGELS failed with error code:", info
    ELSE
      ! The solution is in the first degree+1 elements of b
      coeffs = b(1:degree+1)
    END IF
    
    DEALLOCATE(A, b, work)
  END SUBROUTINE POLYNOMIAL_REGRESSION


END MODULE FIT_MODULE