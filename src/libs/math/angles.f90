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

  TYPE, PUBLIC :: ANGLE
    DOUBLE PRECISION :: VALUE ! Angle value in radians
    LOGICAL :: IS_DEGREES = .FALSE. ! Flag indicating if the value is in degrees
  CONTAINS
    PROCEDURE :: CONVERT_TO_RAD
    PROCEDURE :: CONVERT_TO_DEG
    PROCEDURE :: TO_RAD
    PROCEDURE :: TO_DEG
    PROCEDURE :: TO_STRING
    PROCEDURE :: ASSIGN_REAL_TO_ANGLE
    PROCEDURE :: ADD_ANGLES
    PROCEDURE :: SUBTRACT_ANGLES
    PROCEDURE :: MULTIPLY_ANGLE_BY_REAL
    PROCEDURE :: DIVIDE_ANGLE_BY_REAL
    PROCEDURE :: COS => COSINE
    PROCEDURE :: SIN => SINE
    PROCEDURE :: TAN => TANGENT
    PROCEDURE :: COTAN => COTANGENT
    GENERIC :: ASSIGNMENT(=) => ASSIGN_REAL_TO_ANGLE
    GENERIC :: OPERATOR(+) => ADD_ANGLES
    GENERIC :: OPERATOR(-) => SUBTRACT_ANGLES
    GENERIC :: OPERATOR(*) => MULTIPLY_ANGLE_BY_REAL
    GENERIC :: OPERATOR(/) => DIVIDE_ANGLE_BY_REAL
  END TYPE ANGLE

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
  PUBLIC :: NEW_ANGLE_DEG, NEW_ANGLE_RAD
  PUBLIC :: ASSIGN_REAL_TO_ANGLE
  PUBLIC :: COSINE, SINE, TANGENT, COTANGENT
  PUBLIC :: ADD_ANGLES, SUBTRACT_ANGLES
  PUBLIC :: MULTIPLY_ANGLE_BY_REAL, DIVIDE_ANGLE_BY_REAL

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

  !> \ingroup angles
  !! \brief Converts an angle to radians if it is in degrees.
  !! \param[in] ANGLE_IN Input angle (TYPE(ANGLE)).
  SUBROUTINE CONVERT_TO_RAD(THIS)
    CLASS(ANGLE), INTENT(INOUT) :: THIS
    IF (THIS%IS_DEGREES) THEN
      THIS%VALUE = DEG_TO_RAD(THIS%VALUE)
      THIS%IS_DEGREES = .FALSE.
    END IF
  END SUBROUTINE CONVERT_TO_RAD

  !> \ingroup angles
  !! \brief Converts an angle to degrees if it is in radians.
  !! \param[in] ANGLE_IN Input angle (TYPE(ANGLE)).
  SUBROUTINE CONVERT_TO_DEG(THIS)
    CLASS(ANGLE), INTENT(INOUT) :: THIS
    IF (.NOT. THIS%IS_DEGREES) THEN
      THIS%VALUE = RAD_TO_DEG(THIS%VALUE)
      THIS%IS_DEGREES = .TRUE.
    END IF
  END SUBROUTINE CONVERT_TO_DEG

  !> \ingroup angles
  !! \brief Converts an angle to radians.
  !! \param[in] THIS Input angle (TYPE(ANGLE)).
  !! \return Angle in radians.
  FUNCTION TO_RAD(THIS) RESULT(RAD)
    CLASS(ANGLE), INTENT(IN) :: THIS
    DOUBLE PRECISION :: RAD
    IF (THIS%IS_DEGREES) THEN
      RAD = DEG_TO_RAD(THIS%VALUE)
    ELSE
      RAD = THIS%VALUE
    END IF
  END FUNCTION TO_RAD

  !> \ingroup angles
  !! \brief Converts an angle to degrees.
  !! \param[in] THIS Input angle (TYPE(ANGLE)).
  !! \return Angle in degrees.
  FUNCTION TO_DEG(THIS) RESULT(DEG)
    CLASS(ANGLE), INTENT(IN) :: THIS
    DOUBLE PRECISION :: DEG
    IF (THIS%IS_DEGREES) THEN
      DEG = THIS%VALUE
    ELSE
      DEG = RAD_TO_DEG(THIS%VALUE)
    END IF
  END FUNCTION TO_DEG

  !> \ingroup angles
  !! \brief Converts an angle to a string representation.
  !! \param[in] THIS Input angle (TYPE(ANGLE)).
  !! \return String representation of the angle.
  FUNCTION TO_STRING(THIS) RESULT(STR)
    CLASS(ANGLE), INTENT(IN) :: THIS
    CHARACTER(LEN=32) :: STR
    IF (THIS%IS_DEGREES) THEN
      WRITE(STR, *) THIS%VALUE, " deg"
    ELSE
      WRITE(STR, *) THIS%VALUE, " rad"
    END IF
  END FUNCTION TO_STRING

  !> \ingroup angles
  !! \brief Assigns a real number to an angle.
  !! \param[in] REAL_VAL Real number to assign.
  !! \param[out] ANGLE_OUT Output angle (TYPE(ANGLE)).
  SUBROUTINE ASSIGN_REAL_TO_ANGLE(ANGLE_OUT, REAL_VAL)
    CLASS(ANGLE), INTENT(OUT) :: ANGLE_OUT
    DOUBLE PRECISION, INTENT(IN) :: REAL_VAL
    ANGLE_OUT%VALUE = REAL_VAL
    ANGLE_OUT%IS_DEGREES = .FALSE. ! Default to radians
  END SUBROUTINE ASSIGN_REAL_TO_ANGLE

  !> \ingroup angles
  !! \brief Constructor for ANGLE from degrees
  FUNCTION NEW_ANGLE_DEG(VALUE_DEG) RESULT(ANG)
    DOUBLE PRECISION, INTENT(IN) :: VALUE_DEG
    TYPE(ANGLE) :: ANG
    
    ANG%VALUE = VALUE_DEG
    ANG%IS_DEGREES = .TRUE.
  END FUNCTION NEW_ANGLE_DEG

  !> \ingroup angles
  !! \brief Constructor for ANGLE from radians
  FUNCTION NEW_ANGLE_RAD(VALUE_RAD) RESULT(ANG)
    DOUBLE PRECISION, INTENT(IN) :: VALUE_RAD
    TYPE(ANGLE) :: ANG
    
    ANG%VALUE = VALUE_RAD
    ANG%IS_DEGREES = .FALSE.
  END FUNCTION NEW_ANGLE_RAD

  !> \ingroup angles
  !! \brief Computes the cosine of an angle.
  !! \param[in] THIS Input angle (TYPE(ANGLE)).
  !! \return Cosine of the angle.
  FUNCTION COSINE(THIS) RESULT(COS_VAL)
    CLASS(ANGLE), INTENT(IN) :: THIS
    DOUBLE PRECISION :: COS_VAL
    COS_VAL = DCOS(THIS%TO_RAD())
  END FUNCTION COSINE
  
  !> \ingroup angles
  !! \brief Computes the sine of an angle.
  !! \param[in] THIS Input angle (TYPE(ANGLE)).
  !! \return Sine of the angle.
  FUNCTION SINE(THIS) RESULT(SIN_VAL)
    CLASS(ANGLE), INTENT(IN) :: THIS
    DOUBLE PRECISION :: SIN_VAL
    SIN_VAL = DSIN(THIS%TO_RAD())
  END FUNCTION SINE

  !> \ingroup angles
  !! \brief Computes the tangent of an angle.
  !! \param[in] THIS Input angle (TYPE(ANGLE)).
  !! \return Tangent of the angle.
  FUNCTION TANGENT(THIS) RESULT(TAN_VAL)
    CLASS(ANGLE), INTENT(IN) :: THIS
    DOUBLE PRECISION :: TAN_VAL
    TAN_VAL = DTAN(THIS%TO_RAD())
  END FUNCTION TANGENT

  !> \ingroup angles
  !! \brief Computes the cotangent of an angle.
  !! \param[in] THIS Input angle (TYPE(ANGLE)).
  !! \return Cotangent of the angle.
  FUNCTION COTANGENT(THIS) RESULT(COTAN_VAL)
    CLASS(ANGLE), INTENT(IN) :: THIS
    DOUBLE PRECISION :: COTAN_VAL
    IF (DTAN(THIS%TO_RAD()) < 1.D-15) THEN
      STOP "COTAN: Angle is zero, cotangent is undefined."
    ELSE
      COTAN_VAL = 1.D0 / DTAN(THIS%TO_RAD())
    END IF
  END FUNCTION COTANGENT

  !> \ingroup angles
  !! \brief Adds two angles.
  !! \param[in] ANGLE1 First angle (TYPE(ANGLE)).
  !! \param[in] ANGLE2 Second angle (TYPE(ANGLE)).
  !! \return Sum of the two angles (TYPE(ANGLE)), if they aren't in the same unit, return in rad.
  FUNCTION ADD_ANGLES(ANGLE1, ANGLE2) RESULT(RESULT_ANGLE)
    CLASS(ANGLE), INTENT(IN) :: ANGLE1, ANGLE2
    TYPE(ANGLE) :: RESULT_ANGLE

    IF (ANGLE1%IS_DEGREES .AND. ANGLE2%IS_DEGREES) THEN
      RESULT_ANGLE%VALUE = ANGLE1%VALUE + ANGLE2%VALUE
      RESULT_ANGLE%IS_DEGREES = .TRUE.
    ELSE
      RESULT_ANGLE%VALUE = ANGLE1%TO_RAD() + ANGLE2%TO_RAD()
      RESULT_ANGLE%IS_DEGREES = .FALSE.
    END IF
  END FUNCTION ADD_ANGLES

  !> \ingroup angles
  !! \brief Subtracts two angles.
  !! \param[in] ANGLE1 First angle (TYPE(ANGLE)).
  !! \param[in] ANGLE2 Second angle (TYPE(ANGLE)).
  !! \return Difference of the two angles (TYPE(ANGLE)), if they aren't in the same unit, return in rad.
  FUNCTION SUBTRACT_ANGLES(ANGLE1, ANGLE2) RESULT(RESULT_ANGLE)
    CLASS(ANGLE), INTENT(IN) :: ANGLE1, ANGLE2
    TYPE(ANGLE) :: RESULT_ANGLE 

    IF (ANGLE1%IS_DEGREES .AND. ANGLE2%IS_DEGREES) THEN
      RESULT_ANGLE%VALUE = ANGLE1%VALUE - ANGLE2%VALUE
      RESULT_ANGLE%IS_DEGREES = .TRUE.
    ELSE
      RESULT_ANGLE%VALUE = ANGLE1%TO_RAD() - ANGLE2%TO_RAD()
      RESULT_ANGLE%IS_DEGREES = .FALSE.
    END IF
  END FUNCTION SUBTRACT_ANGLES

  !> \ingroup angles
  !! \brief Multiplies an angle by a real number.
  !! \param[in] ANGLE_IN Input angle (TYPE(ANGLE)).
  !! \param[in] REAL_VAL Real number to multiply by.
  !! \return Resulting angle (TYPE(ANGLE)).
  FUNCTION MULTIPLY_ANGLE_BY_REAL(ANGLE_IN, REAL_VAL) RESULT(RESULT_ANGLE)
    CLASS(ANGLE), INTENT(IN) :: ANGLE_IN
    DOUBLE PRECISION, INTENT(IN) :: REAL_VAL
    TYPE(ANGLE) :: RESULT_ANGLE

    RESULT_ANGLE%VALUE = ANGLE_IN%VALUE * REAL_VAL
    RESULT_ANGLE%IS_DEGREES = ANGLE_IN%IS_DEGREES
  END FUNCTION MULTIPLY_ANGLE_BY_REAL


  !> \ingroup angles
  !! \brief Divides an angle by a real number.
  !! \param[in] ANGLE_IN Input angle (TYPE(ANGLE)).
  !! \param[in] REAL_VAL Real number to divide by.
  !! \return Resulting angle (TYPE(ANGLE)).
  FUNCTION DIVIDE_ANGLE_BY_REAL(ANGLE_IN, REAL_VAL) RESULT(RESULT_ANGLE)
    CLASS(ANGLE), INTENT(IN) :: ANGLE_IN
    DOUBLE PRECISION, INTENT(IN) :: REAL_VAL
    TYPE(ANGLE) :: RESULT_ANGLE 

    RESULT_ANGLE%VALUE = ANGLE_IN%VALUE / REAL_VAL
    RESULT_ANGLE%IS_DEGREES = ANGLE_IN%IS_DEGREES
  END FUNCTION DIVIDE_ANGLE_BY_REAL

END MODULE ANGLES