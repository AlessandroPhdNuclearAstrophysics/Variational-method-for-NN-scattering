!> \file Laguerre_polynomials.f90
!! \brief Laguerre polynomial evaluation utilities for variational calculations.
!! \defgroup laguerre Laguerre Polynomials
!! \ingroup math
!!
!! This module provides routines to compute Laguerre polynomials and their first and
!! second derivatives for use in basis expansions and quadrature.
!!
!! \author Alessandro
!! \date 2025
MODULE LAGUERRE_POLYNOMIAL_MOD
CONTAINS

  !> \brief Compute Laguerre polynomials and their derivatives.
  !! \ingroup laguerre
  !!
  !! This subroutine evaluates the generalized Laguerre polynomials L_n^{(α)}(x)
  !! and their first and second derivatives for a set of points.
  !!
  !! \param[in]  XX   Array of points at which to evaluate the polynomials
  !! \param[in]  APF  The α parameter of the generalized Laguerre polynomial
  !! \param[out] U    Array (0:NMAX, NX) of polynomial values
  !! \param[out] U1   Array (0:NMAX, NX) of first derivatives
  !! \param[out] U2   Array (0:NMAX, NX) of second derivatives
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