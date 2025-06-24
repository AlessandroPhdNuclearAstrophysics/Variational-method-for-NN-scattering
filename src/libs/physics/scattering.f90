!> \file scattering.f90
!! \defgroup nn_scattering Scattering NN
!! \defgroup nn_scattering_utils Scattering NN Utilities
!! \ingroup nn_scattering
!! \brief Utilities for phase shift and scattering calculations in nuclear/particle physics.
!!
!! This module provides types and routines for phase shift calculations, S-matrix construction,
!! and zero-energy observables in both Blatt-Biedenharn and Stapp conventions.
!!
!! \author Alessandro
!! \date 2025
MODULE SCATTERING
  USE ANGLES
  IMPLICIT NONE
  PRIVATE

  !> \brief Structure to store the results of phase shift calculations.
  !! \ingroup nn_scattering
  !! Contains phase shifts and mixing angles in both Blatt-Biedenharn and Stapp conventions.
  TYPE, PUBLIC :: PHASE_SHIFTS_STRUCT
    DOUBLE PRECISION :: DELTA1 = 0.D0        !< Phase shift 1 [deg]
    DOUBLE PRECISION :: DELTA2 = 0.D0        !< Phase shift 2 [deg]
    DOUBLE PRECISION :: MIXING = 0.D0        !< Mixing angle [deg]
    LOGICAL :: DEGREES = .TRUE.              !< Flag to indicate if angles are in degrees
    LOGICAL :: STAPP = .TRUE.                !< Flag to indicate if angles are in Stapp convention
  END TYPE PHASE_SHIFTS_STRUCT

  !> \brief Structure to store zero-energy scattering observables.
  !! \ingroup nn_scattering
  !! Contains scattering lengths and mixing angle in both conventions.
  TYPE, PUBLIC :: ZERO_ENERGY_OBSERVABLES
    DOUBLE PRECISION :: a1 = 0.D0         !< Scattering length a_1 [fm^(2l_1+1)]
    DOUBLE PRECISION :: a2 = 0.D0         !< Scattering length a_2 [fm^(2l_2+1)]
    DOUBLE PRECISION :: e  = 0.D0         !< Mixing angle [fm^(l_1+l_2+1) (Stapp) or fm^-2 (BB)]
    LOGICAL :: STAPP = .TRUE.             !< Flag to indicate if the values are in the Stapp convention
  END TYPE ZERO_ENERGY_OBSERVABLES

  DOUBLE COMPLEX, PARAMETER :: IM = (0.0D0, 1.0D0)  !< Imaginary unit

  PUBLIC :: CALCULATE_PHASE_SHIFTS_BLATT_RAD
  PUBLIC :: CALCULATE_PHASE_SHIFTS_BLATT_DEG
  PUBLIC :: CALCULATE_S_MATRIX_FROM_BLATT
  PUBLIC :: CALCULATE_PHASE_SHIFTS_STAPP_RAD
  PUBLIC :: CALCULATE_PHASE_SHIFTS_STAPP_DEG
  PUBLIC :: EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB
  PUBLIC :: EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP
  PUBLIC :: TRANSFORM_OBSERVABLES_TO_STAPP
  PUBLIC :: TRANSFORM_OBSERVABLES_TO_BB

CONTAINS

  !> \brief Calculate phase shifts in Blatt-Biedenharn convention (radians).
  !! \ingroup nn_scattering
  !! Given an R-matrix and its dimension, computes the phase shifts and mixing angle
  !! in the Blatt-Biedenharn convention (in radians).
  !! \param[in] R_MATRIX The R-matrix (1x1 or 2x2)
  !! \param[in] DIM Dimension of the matrix (1 or 2)
  !! \return PHASE_SHIFTS_STRUCT with phase shifts and mixing angle (radians)
  FUNCTION CALCULATE_PHASE_SHIFTS_BLATT_RAD(R_MATRIX, DIM) RESULT(PS)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: DIM
    DOUBLE PRECISION, INTENT(IN) :: R_MATRIX(:,:)
    TYPE(PHASE_SHIFTS_STRUCT) :: PS
    DOUBLE PRECISION :: AMIX_BB, DELTA_1BB, DELTA_2BB
    PS%DEGREES = .FALSE.
    PS%STAPP   = .FALSE.

    IF (DIM == 1) THEN
      PS%DELTA1 = ATAN((R_MATRIX(1,1)))
      PS%DELTA2 = 0.D0
      PS%MIXING = 0.D0
      RETURN
    ENDIF

    IF (DIM > 2) THEN
      PRINT *, "ERROR: R_MATRIX dimension is greater than 2 in CALCULATE_PHASE_SHIFTS_BLATT"
      STOP
    ENDIF
    IF (R_MATRIX(1,2) /= 0.D0 .AND. R_MATRIX(1,1) == R_MATRIX(2,2)) THEN
      PRINT *, "ERROR: R_MATRIX is singular in CALCULATE_PHASE_SHIFTS_BLATT"
      STOP
    ENDIF
    AMIX_BB=0.5D0*ATAN(2.*R_MATRIX(1,2)/(R_MATRIX(1,1)-R_MATRIX(2,2)))

    DELTA_1BB=ATAN((COS(AMIX_BB)*COS(AMIX_BB)*R_MATRIX(1,1)  &
                  +SIN(AMIX_BB)*SIN(AMIX_BB)*R_MATRIX(2,2)  &
                  +2*COS(AMIX_BB)*SIN(AMIX_BB)*R_MATRIX(1,2)))

    DELTA_2BB=ATAN((SIN(AMIX_BB)*SIN(AMIX_BB)*R_MATRIX(1,1)  &
                  +COS(AMIX_BB)*COS(AMIX_BB)*R_MATRIX(2,2)  &
                  -2*COS(AMIX_BB)*SIN(AMIX_BB)*R_MATRIX(1,2)))

    PS%DELTA1 = DELTA_1BB
    PS%DELTA2 = DELTA_2BB
    PS%MIXING = AMIX_BB
  END FUNCTION CALCULATE_PHASE_SHIFTS_BLATT_RAD

  !> \brief Calculate phase shifts in Blatt-Biedenharn convention (degrees).
  !! \ingroup nn_scattering
  !! Converts the output of CALCULATE_PHASE_SHIFTS_BLATT_RAD to degrees.
  !! \param[in] R_MATRIX The R-matrix (1x1 or 2x2)
  !! \param[in] DIM Dimension of the matrix (1 or 2)
  !! \return PHASE_SHIFTS_STRUCT with phase shifts and mixing angle (degrees)
  FUNCTION CALCULATE_PHASE_SHIFTS_BLATT_DEG(R_MATRIX, DIM) RESULT(PS)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: DIM
    DOUBLE PRECISION, INTENT(IN) :: R_MATRIX(:,:)
    TYPE(PHASE_SHIFTS_STRUCT) :: PS
    PS = CALCULATE_PHASE_SHIFTS_BLATT_RAD(R_MATRIX, DIM)
    PS%DELTA1 = RAD_TO_DEG(PS%DELTA1)
    PS%DELTA2 = RAD_TO_DEG(PS%DELTA2)
    PS%MIXING = RAD_TO_DEG(PS%MIXING)
    
    PS%DEGREES = .TRUE.
    PS%STAPP   = .FALSE.
  END FUNCTION CALCULATE_PHASE_SHIFTS_BLATT_DEG

  !> \brief Construct S-matrix from Blatt-Biedenharn phase shifts.
  !! \ingroup nn_scattering
  !! Given phase shifts and mixing angle in Blatt-Biedenharn convention, constructs the S-matrix.
  !! \param[in] PS_BB Phase shifts structure (Blatt-Biedenharn)
  !! \param[in] DIM Dimension (1 or 2)
  !! \param[out] S_MATRIX The resulting S-matrix (complex)
  SUBROUTINE CALCULATE_S_MATRIX_FROM_BLATT(PS_BB, DIM, S_MATRIX)
    IMPLICIT NONE
    TYPE(PHASE_SHIFTS_STRUCT), INTENT(IN) :: PS_BB
    INTEGER, INTENT(IN) :: DIM
    DOUBLE COMPLEX, INTENT(OUT) :: S_MATRIX(:,:)
    DOUBLE COMPLEX :: SM1, SM2
    DOUBLE PRECISION :: COS1, SIN1
    DOUBLE PRECISION :: DELTA_1, DELTA_2, MIXING
    IF (DIM < 1 .OR. DIM > 2) THEN
      PRINT *, "Error: Dimension must be 1 or 2 in CALCULATE_S_MATRIX_FROM_BLATT"
      STOP
    ENDIF
    IF (PS_BB%STAPP) THEN
      PRINT *, "Error: Stapp convention is not supported in CALCULATE_S_MATRIX_FROM_BLATT"
      STOP
    ENDIF

    IF (PS_BB%DEGREES) THEN
      DELTA_1 = RAD_TO_DEG(PS_BB%DELTA1)
      DELTA_2 = RAD_TO_DEG(PS_BB%DELTA2)
      MIXING  = RAD_TO_DEG(PS_BB%MIXING)
    ELSE
      DELTA_1 = PS_BB%DELTA1
      DELTA_2 = PS_BB%DELTA2
      MIXING  = PS_BB%MIXING
    ENDIF

    SM1  = CDEXP(2.D0*IM*PS_BB%DELTA1)
    SM2  = CDEXP(2.D0*IM*PS_BB%DELTA2)
    COS1 = COS(PS_BB%MIXING)
    SIN1 = SIN(PS_BB%MIXING)
    S_MATRIX(1,1) = COS1*COS1*SM1 + SIN1*SIN1*SM2
    IF ( DIM == 1 ) RETURN
    S_MATRIX(2,2) = COS1*COS1*SM2 + SIN1*SIN1*SM1
    S_MATRIX(1,2) = COS1*SIN1*(SM1 - SM2)
    S_MATRIX(2,1) = S_MATRIX(1,2)
  END SUBROUTINE CALCULATE_S_MATRIX_FROM_BLATT

  !> \brief Calculate phase shifts in Stapp convention (radians).
  !! \ingroup nn_scattering
  !! Computes phase shifts and mixing angle in the Stapp convention (radians) from
  !! a given R-matrix and S-matrix.
  !! \param[in] R_MATRIX_BB The R-matrix (Blatt-Biedenharn)
  !! \param[in] S_MATRIX The S-matrix
  !! \param[in] DIM Dimension (1 or 2)
  !! \return PHASE_SHIFTS_STRUCT with phase shifts and mixing angle (radians)
  FUNCTION CALCULATE_PHASE_SHIFTS_STAPP_RAD(R_MATRIX_BB, S_MATRIX, DIM) RESULT(PS_S)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R_MATRIX_BB(:,:)
    DOUBLE COMPLEX, INTENT(IN) :: S_MATRIX(:,:)
    INTEGER, INTENT(IN) :: DIM
    TYPE(PHASE_SHIFTS_STRUCT) :: PS_S

    DOUBLE PRECISION, PARAMETER :: TOL = 1.D-14
    DOUBLE COMPLEX :: SM1, SM2
    DOUBLE PRECISION :: SI2E, CI2E, CX, SX

    PS_S%DEGREES = .FALSE.
    PS_S%STAPP   = .TRUE.

    IF (DIM == 1) THEN
      PS_S%DELTA1 = DATAN(R_MATRIX_BB(1,1))
      PS_S%DELTA2 = 0.D0
      PS_S%MIXING  = 0.D0
      RETURN
    ENDIF

    SM1  = REAL(S_MATRIX(1,1)*S_MATRIX(2,2)-S_MATRIX(1,2)*S_MATRIX(1,2))
    SI2E =-REAL(S_MATRIX(1,2)*S_MATRIX(1,2)/SM1)
    CI2E = 1.D0 - SI2E
    IF (CI2E < 0.D0) THEN
      PRINT *, "Error: CI2E is negative in CALCULATE_PHASE_SHIFTS_STAPP", CI2E
      STOP
    ENDIF
    IF (SI2E < 0.D0) THEN
      PRINT *, "Error: SI2E is negative in CALCULATE_PHASE_SHIFTS_STAPP", SI2E
      STOP
    ENDIF
    SI2E = DSQRT(SI2E)
    CI2E = DSQRT(CI2E)

  ! I MIXING ANGLES DELLE ONDE DISPARI (JP=0) VENGONO COL SEGNO SBAGLIATO!
    SM1 = CDSQRT(S_MATRIX(1,1)/CI2E)
    SM2 = CDSQRT(S_MATRIX(2,2)/CI2E)
    CX  = DREAL(SM1)
    SX  = DIMAG(SM1)
    IF (ABS(CX) > 1.D0) THEN
      IF (ABS(CX-1) < TOL) THEN
        CX = CX/DABS(CX)
      ELSE
        PRINT *, "Error: CX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", CX, " CX - 1", CX-1
        STOP
      ENDIF 
    ENDIF
    IF (ABS(SX) > 1.D0) THEN
      IF (ABS(SX-1) < TOL) THEN
        SX = SX/DABS(SX)
      ELSE
        PRINT *, "Error: SX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", SX, " SX - 1", SX-1
        STOP
      ENDIF
    ENDIF
    PS_S%DELTA1 = DACOS(CX)
    IF(SX < 0.D0) PS_S%DELTA1 = -PS_S%DELTA1

    CX = DREAL(SM2)
    SX = DIMAG(SM2)
    IF (ABS(CX) > 1.D0) THEN
      IF (ABS(CX-1) < TOL) THEN
        CX = CX/DABS(CX)
      ELSE
        PRINT *, "Error: CX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", CX, " CX - 1", CX-1
        STOP
      ENDIF
    ENDIF
    IF (ABS(SX) > 1.D0) THEN
      IF (ABS(SX-1) < TOL) THEN
        SX = SX/DABS(SX)
      ELSE
        PRINT *, "Error: SX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", SX, " SX - 1", SX-1
        STOP
      ENDIF
    ENDIF
    PS_S%DELTA2 = DACOS(CX)
    IF(SX < 0.D0) PS_S%DELTA2 = -PS_S%DELTA2

    IF (ABS(SI2E) > 1.D0) THEN
      IF (ABS(SI2E-1) < TOL) THEN
        SI2E = SI2E/DABS(SI2E)
      ELSE
        PRINT *, "Error: SI2E is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", SI2E, " SI2E - 1", SI2E-1
        STOP
      ENDIF
    ENDIF
    PS_S%MIXING = 0.5D0*DASIN(SI2E)
  END FUNCTION CALCULATE_PHASE_SHIFTS_STAPP_RAD

  !> \brief Calculate phase shifts in Stapp convention (degrees).
  !! \ingroup nn_scattering
  !! Converts the output of CALCULATE_PHASE_SHIFTS_STAPP_RAD to degrees.
  !! \param[in] R_MATRIX_BB The R-matrix (Blatt-Biedenharn)
  !! \param[in] S_MATRIX The S-matrix
  !! \param[in] DIM Dimension (1 or 2)
  !! \return PHASE_SHIFTS_STRUCT with phase shifts and mixing angle (degrees)
  FUNCTION CALCULATE_PHASE_SHIFTS_STAPP_DEG(R_MATRIX_BB, S_MATRIX, DIM) RESULT(PS_S)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R_MATRIX_BB(:,:)
    DOUBLE COMPLEX, INTENT(IN) :: S_MATRIX(:,:)
    INTEGER, INTENT(IN) :: DIM
    TYPE(PHASE_SHIFTS_STRUCT) :: PS_S

    PS_S = CALCULATE_PHASE_SHIFTS_STAPP_RAD(R_MATRIX_BB, S_MATRIX, DIM)
    PS_S%DELTA1 = RAD_TO_DEG(PS_S%DELTA1)
    PS_S%DELTA2 = RAD_TO_DEG(PS_S%DELTA2)
    PS_S%MIXING = RAD_TO_DEG(PS_S%MIXING)

    PS_S%DEGREES = .TRUE.
    PS_S%STAPP   = .TRUE.
  END FUNCTION CALCULATE_PHASE_SHIFTS_STAPP_DEG

  !> \brief Evaluate zero-energy observables in Blatt-Biedenharn convention.
  !! \ingroup nn_scattering
  !! Computes zero-energy scattering observables (scattering lengths, mixing) in the Blatt-Biedenharn convention.
  !! \param[in] R_MATRIX The R-matrix (1x1 or 2x2)
  !! \param[in] DIM Dimension of the matrix (1 or 2)
  !! \param[in] LMIN Minimum orbital angular momentum
  !! \return ZERO_ENERGY_OBSERVABLES structure
  FUNCTION EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R_MATRIX, DIM, LMIN) RESULT(OBS)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R_MATRIX(:,:)
    INTEGER, INTENT(IN) :: DIM, LMIN
    TYPE(ZERO_ENERGY_OBSERVABLES) :: OBS
    DOUBLE PRECISION :: A1, A2, E
    INTEGER, EXTERNAL :: DOUBLE_FACTORIAL

    IF (DIM == 1) THEN
      OBS%a1 =-R_MATRIX(1,1)/DOUBLE_FACTORIAL(2*LMIN+1)**2
      OBS%e  = 0.D0
      OBS%a2 = 0.D0
    ELSEIF (DIM == 2) THEN
      OBS%a1 =-R_MATRIX(1,1)/DOUBLE_FACTORIAL(2*LMIN+1)**2
      OBS%e  = R_MATRIX(1,2)/R_MATRIX(1,1) * DOUBLE_FACTORIAL(2*LMIN+1)/DOUBLE_FACTORIAL(2*LMIN+5) 
      OBS%a2 =(R_MATRIX(1,2)**2 - R_MATRIX(1,1)*R_MATRIX(2,2)) / (DOUBLE_FACTORIAL(2*LMIN+5)**2 * R_MATRIX(1,1))
    ELSE
      PRINT *, "Error: Invalid dimension in EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB"
      STOP
    ENDIF
    OBS%STAPP = .FALSE.
  END FUNCTION EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB

  !> \brief Evaluate zero-energy observables in Stapp convention.
  !! \ingroup nn_scattering
  !! Computes zero-energy scattering observables (scattering lengths, mixing) in the Stapp convention.
  !! \param[in] R_MATRIX The R-matrix (1x1 or 2x2)
  !! \param[in] DIM Dimension (1 or 2)
  !! \param[in] LMIN Minimum orbital angular momentum
  !! \return ZERO_ENERGY_OBSERVABLES structure
  FUNCTION EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R_MATRIX, DIM, LMIN) RESULT(OBS)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R_MATRIX(:,:)
    INTEGER, INTENT(IN) :: DIM, LMIN
    TYPE(ZERO_ENERGY_OBSERVABLES) :: OBS, OBS_BB
    DOUBLE PRECISION :: A1, A2, E

    OBS_BB = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R_MATRIX, DIM, LMIN)
    OBS = TRANSFORM_OBSERVABLES_TO_STAPP(OBS_BB)
  END FUNCTION EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP

  !> \brief Transform zero-energy observables to Stapp convention.
  !! \ingroup nn_scattering
  !! Converts zero-energy observables from Blatt-Biedenharn to Stapp convention.
  !! \param[in] OBS_BB Observables in Blatt-Biedenharn convention
  !! \return Observables in Stapp convention
  FUNCTION TRANSFORM_OBSERVABLES_TO_STAPP(OBS_BB) RESULT(OBS_S)
    IMPLICIT NONE
    TYPE(ZERO_ENERGY_OBSERVABLES), INTENT(IN) :: OBS_BB
    TYPE(ZERO_ENERGY_OBSERVABLES) :: OBS_S

    IF (OBS_BB%STAPP) THEN
      PRINT *, "ERROR: OBS_BB is already in Stapp convention in TRANSFORM_OBSERVABLES_TO_STAPP"
      STOP
    ENDIF

    OBS_S%a1 = OBS_BB%a1
    OBS_S%a2 = OBS_BB%a2 + OBS_BB%a1 * OBS_BB%e**2
    OBS_S%e  =-OBS_BB%a1 * OBS_BB%e
    OBS_S%STAPP = .TRUE.
  END FUNCTION TRANSFORM_OBSERVABLES_TO_STAPP

  !> \brief Transform zero-energy observables to Blatt-Biedenharn convention.
  !! \ingroup nn_scattering
  !! Converts zero-energy observables from Stapp to Blatt-Biedenharn convention.
  !! \param[in] OBS_S Observables in Stapp convention
  !! \return Observables in Blatt-Biedenharn convention
  FUNCTION TRANSFORM_OBSERVABLES_TO_BB(OBS_S) RESULT(OBS_BB)
    IMPLICIT NONE
    TYPE(ZERO_ENERGY_OBSERVABLES), INTENT(IN) :: OBS_S
    TYPE(ZERO_ENERGY_OBSERVABLES) :: OBS_BB

    IF (.NOT. OBS_S%STAPP) THEN
      PRINT *, "Error: OBS_S is not in Stapp convention in TRANSFORM_OBSERVABLES_TO_BLATT"
      STOP
    ENDIF

    OBS_BB%a1 = OBS_S%a1
    OBS_BB%a2 = OBS_S%a2 - OBS_S%e**2 / OBS_S%a1
    OBS_BB%e  =-OBS_S%e / OBS_S%a1
    OBS_BB%STAPP = .FALSE.
  END FUNCTION TRANSFORM_OBSERVABLES_TO_BB

END MODULE SCATTERING