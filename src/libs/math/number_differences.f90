!> \file number_differences.f90
!! \brief Provides functions for numerical difference and comparison operations.
!! \defgroup num_diff Numerical Comparison Utilities
!! \ingroup math
!!
!! This module defines utility functions for comparing floating-point numbers, including
!! absolute, relative, and percent differences, with optional tolerance and threshold parameters.
!!
!! \author Alessandro Grassi
!! \date 2025
MODULE NUMBER_DIFFERENCES
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ABS_DIFF_PROCENT
  PUBLIC :: ABS_DIFF_RELATIVE
  PUBLIC :: ABS_DIFF

CONTAINS
  !> \brief Computes the percent difference between two numbers.
  !! \ingroup num_diff
  !! \param[in] NUM The number to compare.
  !! \param[in] REF The reference value.
  !! \param[in] TOLERANCE Optional tolerance for comparison (default: 1e-10).
  !! \param[in] THRESHOLD Optional threshold below which values are considered zero (default: 0).
  !! \return Percent difference (in percent units).
  FUNCTION ABS_DIFF_PROCENT(NUM, REF, TOLERANCE, THRESHOLD) RESULT(DIFF)
    DOUBLE PRECISION, INTENT(IN) :: NUM, REF
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: TOLERANCE, THRESHOLD
    DOUBLE PRECISION :: DIFF, TOLERANCE_, THRESHOLD_

    TOLERANCE_ = 1.D-10
    THRESHOLD_ = 0.D0
    IF (PRESENT(TOLERANCE)) THEN
      TOLERANCE_ = TOLERANCE
    END IF
    IF (PRESENT(THRESHOLD)) THEN
      THRESHOLD_ = THRESHOLD
    END IF
    
    DIFF = ABS_DIFF_RELATIVE(NUM, REF, TOLERANCE_, THRESHOLD_)
    DIFF = DIFF * 100.0D0
  END FUNCTION ABS_DIFF_PROCENT

  !> \ingroup num_diff
  !! \brief Computes the relative difference between two numbers.
  !! \param[in] NUM The number to compare.
  !! \param[in] REF The reference value.
  !! \param[in] TOLERANCE Optional tolerance for comparison (default: 1e-10).
  !! \param[in] THRESHOLD Optional threshold below which values are considered zero (default: 0).
  !! \return Relative difference (0 to 1).
  FUNCTION ABS_DIFF_RELATIVE(NUM, REF, TOLERANCE, THRESHOLD) RESULT(DIFF)
    DOUBLE PRECISION, INTENT(IN) :: NUM, REF
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: TOLERANCE, THRESHOLD
    DOUBLE PRECISION :: DIFF, TOLERANCE_, THRESHOLD_

    TOLERANCE_ = 1.D-10
    THRESHOLD_ = 0.D0
    IF (PRESENT(TOLERANCE)) THEN
      TOLERANCE_ = TOLERANCE
    END IF
    IF (PRESENT(THRESHOLD)) THEN
      THRESHOLD_ = THRESHOLD
    END IF

    IF (ABS(NUM) < THRESHOLD_ .AND. ABS(REF) < THRESHOLD_) THEN
      DIFF = 0.D0
    ELSEIF (ABS(NUM) < THRESHOLD_ .AND. ABS(REF) >= THRESHOLD_ .OR. &
             ABS(NUM) >= THRESHOLD_ .AND. ABS(REF) < THRESHOLD_) THEN
      DIFF = 1.D0
    ELSE
      DIFF = ABS_DIFF(NUM - REF, THRESHOLD_) / MAX(MIN(ABS(REF),ABS(NUM)), TOLERANCE_)
    END IF

  END FUNCTION ABS_DIFF_RELATIVE
    

  !> \ingroup num_diff
  !! \brief Computes the absolute difference between two numbers, with optional threshold.
  !! \param[in] A First number.
  !! \param[in] B Second number.
  !! \param[in] THRESHOLD Optional threshold below which values are considered zero (default: 0).
  !! \return Absolute difference.
  FUNCTION ABS_DIFF(A, B, THRESHOLD) RESULT(DIFF)
    DOUBLE PRECISION, INTENT(IN) :: A, B
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: THRESHOLD
    DOUBLE PRECISION :: DIFF, THRESHOLD_

    THRESHOLD_ = 0.D0
    IF (PRESENT(THRESHOLD)) THEN
      THRESHOLD_ = THRESHOLD
    END IF

    IF (ABS(A) < THRESHOLD_ .AND. ABS(B) < THRESHOLD_) THEN
      DIFF = 0.D0
    ELSE
      DIFF = ABS(A - B)
    ENDIF
  END FUNCTION ABS_DIFF

END MODULE NUMBER_DIFFERENCES