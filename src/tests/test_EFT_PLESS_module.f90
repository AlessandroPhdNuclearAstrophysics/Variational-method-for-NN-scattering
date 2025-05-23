! filepath: /home/alessandro/Dropbox/Variazionale_mine/tests/test_eft_pless.f90
PROGRAM test_eft_pless
  USE EFT_PLESS
  IMPLICIT NONE

  INTEGER :: ilb, l, s, j, t1z, t2z, lemp
  DOUBLE PRECISION :: r
  DOUBLE PRECISION :: vpw(2,2)
  TYPE(LECS_EFT_PLESS) :: lecs

  ! Example values (adjust as needed for your data)
  ilb = 1
  l = 0
  s = 1
  j = 1
  t1z = 1
  t2z = -1
  r = 1.0d0
  lemp = 0

  ! Test: get LECS
  lecs = GET_LECS(ilb)
  PRINT *, "LECS%ILB: ", lecs%ILB
  PRINT *, "LECS%RC: ", lecs%RC

  ! Test: call EFT_PLESS_PW
  CALL EFT_PLESS_PW(ilb, l, s, j, t1z, t2z, r, vpw, lemp)
  PRINT *, "VPW:"
  PRINT *, vpw

END PROGRAM test_eft_pless