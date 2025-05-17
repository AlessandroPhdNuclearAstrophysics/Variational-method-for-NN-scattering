MODULE LAGUERRE_POLYNOMIAL_MOD
CONTAINS
  SUBROUTINE LAGUERRE_POLYNOMIAL(XX, APF, U, U1, U2)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: XX(:), APF
    DOUBLE PRECISION, INTENT(OUT) :: U(0:, :), U1(0:, :), U2(0:, :)
    INTEGER :: NX, I, N, NMAX
    DOUBLE PRECISION :: X, D, A0, A1, A2

    NX = SIZE(XX)
    NMAX = UBOUND(U,1) - 1

    DO I = 1, NX
      X = XX(I)
      U (0,I) = 1.D0
      U (1,I) = APF + 1.D0 - X
      U1(0,I) = 0.D0
      U1(1,I) = -1.D0
      U2(0,I) = 0.D0
      U2(1,I) = 0.D0
    END DO

    DO N = 1, NMAX
      D = 1.D0 / (N + 1.D0)
      A1 = 2.D0 * N + APF + 1.D0
      A0 = N + APF
      A2 = N + 1 + APF
      DO I = 1, NX
        X = XX(I)
        U (N+1,I) = ((A1 - X) * U(N,I) - A0 * U(N-1,I)) * D
        U1(N+1,I) = ((N+1) * U(N+1,I) - A2 * U(N,I)) / X
        U2(N+1,I) = (-U(1,I) * U1(N+1,I) - (N+1) * U(N+1,I)) / X
      END DO
    END DO
  END SUBROUTINE LAGUERRE_POLYNOMIAL
END MODULE LAGUERRE_POLYNOMIAL_MOD