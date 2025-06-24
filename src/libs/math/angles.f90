!> \defgroup math Mathematical Utilities
!! \file angles.f90
!! \defgroup angles Angles and Conversions
!! \ingroup math
!! \brief Provides constants and functions for angle conversions between degrees and radians.
!! \author Alessandro Grassi
!! \date 2025
!!
!! This module defines mathematical constants for angle conversions and provides utility functions
!! to convert between degrees and radians. All constants and functions are documented and grouped
!! for Doxygen output.

MODULE ANGLES
  IMPLICIT NONE
  PRIVATE

  !> \ingroup angles
  !! \brief The mathematical constant \f$\pi\f$ (pi).
  DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)

  !> \ingroup angles
  !! \brief Conversion factor from radians to degrees (180/\f$\pi\f$).
  DOUBLE PRECISION, PARAMETER :: ONEEIGHTY_OVER_PI = 180.D0 / PI

  !> \ingroup angles
  !! \brief Conversion factor from degrees to radians (\f$\pi\f$/180).
  DOUBLE PRECISION, PARAMETER :: PI_OVER_ONEEIGHTY = PI / 180.D0

  PUBLIC :: DEG_TO_RAD, RAD_TO_DEG
  PUBLIC :: PI, ONEEIGHTY_OVER_PI, PI_OVER_ONEEIGHTY

CONTAINS

  !> \ingroup angles
  !! \brief Converts an angle in degrees to radians.
  !! \param[in] DEG Angle in degrees.
  !! \return Angle in radians.
  FUNCTION DEG_TO_RAD(DEG) RESULT(RAD)
    DOUBLE PRECISION, INTENT(IN) :: DEG
    DOUBLE PRECISION :: RAD

    RAD = DEG * PI_OVER_ONEEIGHTY
  END FUNCTION DEG_TO_RAD

  !> \ingroup angles
  !! \brief Converts an angle in radians to degrees.
  !! \param[in] RAD Angle in radians.
  !! \return Angle in degrees.
  FUNCTION RAD_TO_DEG(RAD) RESULT(DEG)
    DOUBLE PRECISION, INTENT(IN) :: RAD
    DOUBLE PRECISION :: DEG

    DEG = RAD * ONEEIGHTY_OVER_PI
  END FUNCTION RAD_TO_DEG

END MODULE ANGLES