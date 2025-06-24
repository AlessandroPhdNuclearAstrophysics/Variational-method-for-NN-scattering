!> \file random_numbers.f90
!! \brief Utilities for generating random numbers in a given interval for various types.
!! \defgroup random_numbers Random Number Generators
!! \ingroup utils
!!
!! This module provides functions to generate random numbers of type integer, real,
!! double precision, and complex (real/imaginary parts) in a specified interval.
!!
!! \author Alessandro Grassi
!! \date 2025

MODULE RANDOM_NUMBERS
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: RANDOM_INT, RANDOM_REAL, RANDOM_DOUBLE, RANDOM_COMPLEX

CONTAINS

  !> \brief Generate a random integer in [MIN, MAX] (inclusive).
  !! \ingroup random_numbers
  !! \param[in] MIN Minimum integer value
  !! \param[in] MAX Maximum integer value
  !! \return Random integer in [MIN, MAX]
  FUNCTION RANDOM_INT(MIN, MAX) RESULT(VAL)
    INTEGER, INTENT(IN) :: MIN, MAX
    INTEGER :: VAL
    REAL :: R
    IF (MIN > MAX) THEN
      VAL = MIN
      RETURN
    END IF
    CALL RANDOM_NUMBER(R)
    VAL = MIN + INT((MAX - MIN + 1) * R)
    IF (VAL > MAX) VAL = MAX
  END FUNCTION RANDOM_INT

  !> \brief Generate a random real in [MIN, MAX].
  !! \ingroup random_numbers
  !! \param[in] MIN Minimum real value
  !! \param[in] MAX Maximum real value
  !! \return Random real in [MIN, MAX]
  FUNCTION RANDOM_REAL(MIN, MAX) RESULT(VAL)
    REAL, INTENT(IN) :: MIN, MAX
    REAL :: VAL, R
    CALL RANDOM_NUMBER(R)
    VAL = MIN + (MAX - MIN) * R
  END FUNCTION RANDOM_REAL

  !> \brief Generate a random double precision number in [MIN, MAX].
  !! \ingroup random_numbers
  !! \param[in] MIN Minimum double value
  !! \param[in] MAX Maximum double value
  !! \return Random double in [MIN, MAX]
  FUNCTION RANDOM_DOUBLE(MIN, MAX) RESULT(VAL)
    DOUBLE PRECISION, INTENT(IN) :: MIN, MAX
    DOUBLE PRECISION :: VAL, R
    CALL RANDOM_NUMBER(R)
    VAL = MIN + (MAX - MIN) * R
  END FUNCTION RANDOM_DOUBLE

  !> \brief Generate a random complex number in [MIN, MAX] for both real and imaginary parts.
  !! \ingroup random_numbers
  !! \param[in] MIN Minimum value for real and imaginary parts
  !! \param[in] MAX Maximum value for real and imaginary parts
  !! \return Random complex number with real and imaginary parts in [MIN, MAX]
  FUNCTION RANDOM_COMPLEX(MIN, MAX) RESULT(VAL)
    DOUBLE PRECISION, INTENT(IN) :: MIN, MAX
    COMPLEX(KIND=KIND(1.0D0)) :: VAL
    DOUBLE PRECISION :: R1, R2
    CALL RANDOM_NUMBER(R1)
    CALL RANDOM_NUMBER(R2)
    VAL = CMPLX(MIN + (MAX - MIN) * R1, MIN + (MAX - MIN) * R2, KIND=KIND(1.0D0))
  END FUNCTION RANDOM_COMPLEX

END MODULE RANDOM_NUMBERS
