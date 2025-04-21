PROGRAM main_test_coulomb
  USE, intrinsic :: iso_c_binding
  USE gsl_coulomb
  IMPLICIT NONE

  ! Declare variables
  TYPE(gsl_sf_result) :: F, G, FP, GP
  REAL(c_double) :: eta, x, F_exp, G_exp
  INTEGER :: L, status

  ! Initialize parameters
  eta = 1.0d0
  x = 2.0d0
  L = 1

  ! Call the Fortran wrapper for GSL's coulomb_wave_FG_e
  CALL coulomb_wave_FG(eta, x, L, F, FP, G, GP, F_exp, G_exp, status)

  ! Print results
  PRINT *, "Input Parameters:"
  PRINT *, "eta: ", eta
  PRINT *, "x: ", x
  PRINT *, "L: ", L
  IF (status /= 0) THEN
    PRINT *, "Error: GSL function call failed with status ", status
    STOP
  END IF
  PRINT *, "Coulomb Wave Function Results:"
  PRINT *, "F: ", F
  PRINT *, "G: ", G
  PRINT *, "FP: ", FP
  PRINT *, "GP: ", GP
  PRINT *, "F_exp: ", F_exp
  PRINT *, "G_exp: ", G_exp

END PROGRAM main_test_coulomb