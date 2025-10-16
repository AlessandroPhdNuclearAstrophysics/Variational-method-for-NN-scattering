!> \file EFT_pless.f90
!! \brief Effective Field Theory (EFT) potential module for NN scattering.
!! \defgroup eft_pless EFT-pionless
!! \ingroup nn_potentials
!!
!! This module defines the low-energy constants (LECs) and radial functions for
!! the pionless EFT potential, and provides routines to evaluate the
!! potential matrix elements for use in variational and Numerov solvers.
!!
!! \details
!! - Supports LO, NLO, and N3LO orders.
!! - Provides routines to set LECs, print them, and combine operator matrix elements.
!! - Used by the main variational scattering code.
!!
!! \author Alessandro
!! \date 2025

MODULE EFT_PLESS
  USE REALLOCATE_UTILS
  USE QUANTUM_NUMBERS
  USE QUANTUM_OPERATORS
  IMPLICIT NONE
  PRIVATE

  DOUBLE PRECISION, PARAMETER :: HTC = 197.32697D0

  !> \brief Structure for storing all LECs for a given EFT model.
  TYPE, PUBLIC :: LECS_EFT_PLESS
    INTEGER :: ILB = -1
    INTEGER :: ORDER = -1
    DOUBLE PRECISION :: RC(0:1,0:1) = 0.D0
    DOUBLE PRECISION :: CLO(0:1)   = 0.D0
    DOUBLE PRECISION :: CNLO(7) = 0.D0
    DOUBLE PRECISION :: CN3LO(11)= 0.D0
    DOUBLE PRECISION :: CIT(0:4)= 0.D0
    LOGICAL :: ADIMENSIONAL = .FALSE.
  CONTAINS 
    PROCEDURE :: PRINT => PRINT_LECS
    PROCEDURE :: CONVERT_TO_ADIMENSIONAL => CONVERT_LECS_TO_ADIMENSIONAL
    PROCEDURE :: CONVERT_TO_DIMENSIONAL => CONVERT_LECS_TO_DIMENSIONAL
  END TYPE LECS_EFT_PLESS

  !> \brief Structure for storing all LECs for a given EFT model with R depending on the order.
  TYPE, PUBLIC :: LECS_EFT_PLESS_R_ORDER
    INTEGER :: ILB = -1
    INTEGER :: ORDER = -1
    DOUBLE PRECISION :: R0(0:1,0:1) = 0.D0
    DOUBLE PRECISION :: R1(0:1,0:1) = 0.D0
    DOUBLE PRECISION :: R3(0:1,0:1) = 0.D0
    DOUBLE PRECISION :: CLO(0:1)   = 0.D0
    DOUBLE PRECISION :: CNLO(7) = 0.D0
    DOUBLE PRECISION :: CN3LO(11)= 0.D0
    DOUBLE PRECISION :: CIT(0:4)= 0.D0
  END TYPE LECS_EFT_PLESS_R_ORDER

  !> \brief Structure for storing radial functions for a given cutoff and order.
  TYPE, PUBLIC :: EFT_RADIAL_FUNCTIONS
    INTEGER :: ORDER = -1
    DOUBLE PRECISION :: RC(0:1,0:1) = 0.D0
    DOUBLE PRECISION, ALLOCATABLE :: FR_I(:,:,:,:)  ! FR_I(S, T, I, R)
  END TYPE EFT_RADIAL_FUNCTIONS

  !> \brief Generic interface to set the Low-Energy Constants (LECs) for the EFT pionless potential.
  !! \ingroup eft_pless
  !!
  !! This interface allows setting the LECs using different argument types:
  !! - From another LECS_EFT_PLESS structure (SET_LECS_FROM_LECS)
  !! - By specifying individual LEC values (SET_LECS_SINGLE)
  !! - By model index (ILB) (SET_LECS_ILB)
  INTERFACE SET_LECS
    MODULE PROCEDURE SET_LECS_FROM_LECS
    MODULE PROCEDURE SET_LECS_SINGLE
    MODULE PROCEDURE SET_LECS_ILB
  END INTERFACE

  LOGICAL, PRIVATE :: FIRST_CALL = .TRUE., LECS_SET = .FALSE., LECS_ALL_SET = .FALSE.
  INTEGER :: NMODELS = -1
  TYPE(LECS_EFT_PLESS), ALLOCATABLE :: LECS_ALL(:)
  TYPE(LECS_EFT_PLESS) :: LECS
  TYPE(SCATTERING_CHANNEL) :: CHANNEL
  INTEGER :: ORDER =-1

  INTEGER, PARAMETER :: I2   (2,2) = RESHAPE([1,0,0,1], [2,2])
  INTEGER            :: LS   (2,2)
  INTEGER            :: L2   (2,2)
  DOUBLE PRECISION   :: S12  (2,2)
  INTEGER            :: T12

  PUBLIC :: EFT_PLESS_PW, EFT_PLESS_PW_LECS, EFT_PLESS_PW_FITTED, LECS_TO_ST_LECS, ST_LECTS_TO_LECS
  PUBLIC :: SET_LECS, GET_LECS, PRINT_LECS
  PUBLIC :: GET_EFT_RADIAL_FUNCTIONS, COMBINE_POTENTIAL
  PUBLIC :: ADIMENSIONAL_LECS_FROM_DIMENSIONAL, DIMENSIONAL_LECS_FROM_ADIMENSIONAL
  PUBLIC :: CONVERT_LECS_TO_ADIMENSIONAL, CONVERT_LECS_TO_DIMENSIONAL

  PRIVATE:: CR, EFT_RADIAL_1, EFT_RADIAL_2, EFT_RADIAL_3, EFT_RADIAL_4, EFT_RADIAL_5, EFT_RADIAL_6, EFT_RADIAL_7
  PRIVATE:: SET_LECS_FROM_LECS, SET_LECS_SINGLE, SET_ALL_LECS, SET_LECS_ILB
  PRIVATE:: SET_OPERATORS
  PRIVATE:: PREPARE

CONTAINS
  !> \brief Set LECs from another LECS_EFT_PLESS structure.
  !! \param[in] LECS_IN Structure containing all LECs for the EFT potential
  SUBROUTINE SET_LECS_FROM_LECS(LECS_IN)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_IN
    LECS = LECS_IN
    LECS_SET = .TRUE.
  END SUBROUTINE SET_LECS_FROM_LECS

  !> @brief Set the low-energy constants (LECs) for the EFT pionless potential using individual values.
  !>
  !> @details
  !> Sets the values of the LECs in the global LECS structure. Each argument corresponds to a specific LEC or radial cutoff.
  !> All arguments are optional; at least one R* and one C* parameter must be provided. If a parameter is not present, its value is not updated.
  !>
  !> R00, R10, R01, R11 are the radial cutoff parameters for different spin/isospin channels.
  !> C10, C01, C200, ..., C4Q are the LECs for various orders and operator structures.
  !> CIT0, CIT1, CIT2, CIT3, CIT4 are the isospin-dependent LECs.
  !>
  !> @param[in] R00   (Optional) Radial cutoff for S=0, T=0 channel.
  !> @param[in] R10   (Optional) Radial cutoff for S=1, T=0 channel.
  !> @param[in] R01   (Optional) Radial cutoff for S=0, T=1 channel.
  !> @param[in] R11   (Optional) Radial cutoff for S=1, T=1 channel.
  !> @param[in] C10   (Optional) Leading-order LEC for S=1, T=0 channel.
  !> @param[in] C01   (Optional) Leading-order LEC for S=1, T=1 channel.
  !> @param[in] C200, C201, C210, C211, C2T0, C2T1, C2B
  !>                (Optional) Next-to-leading-order LECs.
  !> @param[in] C400, C410, C401, C411, C4T0, C4T1, C4B0, C4B1, C4BB0, C4BB1, C4Q
  !>                (Optional) N3LO LECs.
  !> @param[in] CIT0, CIT1, CIT2, CIT3, CIT4
  !>                (Optional) Isospin-dependent LECs.
  !>
  !> @note
  !> All arguments are optional, but at least one R* and one C* must be present.
  !> The LECs are stored in the global LECS structure for use by other routines.
  SUBROUTINE SET_LECS_SINGLE(R00, R10, R01, R11, C10, C01, C200, C201, C210, C211, C2T0, C2T1, C2B, &
    C400, C410, C401, C411, C4T0, C4T1, C4B0, C4B1, C4BB0, C4BB1, C4Q, &
    CIT0, CIT1, CIT2, CIT3, CIT4)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: NUMERICAL_ZERO = 1.D-16
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: R00, R10, R01, R11
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: C10, C01
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: C200, C201, C210, C211, C2T0, C2T1, C2B
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: C400, C410, C401, C411, C4T0, C4T1, C4B0, C4B1, C4BB0, C4BB1, C4Q
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: CIT0, CIT1, CIT2, CIT3, CIT4
    LOGICAL :: AT_LEAST_ONE_LEC_IS_NOT_NULL

    IF (.NOT.(PRESENT(R00) .OR. PRESENT(R10) .OR. PRESENT(R01) .OR. PRESENT(R11))) THEN
      STOP "AT LEAST ONE OF R MUST BE PRESENT"
    END IF
    IF (.NOT.(PRESENT(C10) .OR. PRESENT(C01) .OR. PRESENT(C200) .OR. PRESENT(C201) .OR. PRESENT(C210) .OR. &
          PRESENT(C211) .OR. PRESENT(C2T0) .OR. PRESENT(C2T1) .OR. PRESENT(C2B) .OR. PRESENT(C400) .OR. &
          PRESENT(C410) .OR. PRESENT(C401) .OR. PRESENT(C411) .OR. PRESENT(C4T0) .OR. PRESENT(C4T1) .OR. &
          PRESENT(C4B0) .OR. PRESENT(C4B1) .OR. PRESENT(C4BB0) .OR. PRESENT(C4BB1) .OR. PRESENT(C4Q) .OR. &
          PRESENT(CIT0) .OR. PRESENT(CIT1) .OR. PRESENT(CIT2) .OR. PRESENT(CIT3) .OR. PRESENT(CIT4))) THEN
      STOP "AT LEAST ONE C MUST BE PRESENT"
    END IF
    AT_LEAST_ONE_LEC_IS_NOT_NULL = .FALSE.
    IF (PRESENT(C10))  THEN
      IF (DABS(C10) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C01))  THEN
      IF (DABS(C01) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C200)) THEN
      IF (DABS(C200) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C201)) THEN
      IF (DABS(C201) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C210)) THEN
      IF (DABS(C210) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C211)) THEN
      IF (DABS(C211) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C2T0)) THEN
      IF (DABS(C2T0) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C2T1)) THEN
      IF (DABS(C2T1) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C2B))  THEN
      IF (DABS(C2B) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C400)) THEN
      IF (DABS(C400) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C410)) THEN
      IF (DABS(C410) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C401)) THEN
      IF (DABS(C401) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C411)) THEN
      IF (DABS(C411) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C4T0)) THEN
      IF (DABS(C4T0) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C4T1)) THEN
      IF (DABS(C4T1) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C4B0)) THEN
      IF (DABS(C4B0) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C4B1)) THEN
      IF (DABS(C4B1) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C4BB0)) THEN
      IF (DABS(C4BB0) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C4BB1)) THEN
      IF (DABS(C4BB1) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(C4Q)) THEN
      IF (DABS(C4Q) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(CIT0)) THEN
      IF (DABS(CIT0) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(CIT1)) THEN
      IF (DABS(CIT1) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(CIT2)) THEN
      IF (DABS(CIT2) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(CIT3)) THEN
      IF (DABS(CIT3) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF
    IF (PRESENT(CIT4)) THEN
      IF (DABS(CIT4) < NUMERICAL_ZERO) AT_LEAST_ONE_LEC_IS_NOT_NULL = .TRUE.
    END IF

    IF (.NOT. AT_LEAST_ONE_LEC_IS_NOT_NULL) THEN
      STOP "AT LEAST ONE OF C MUST BE PRESENT"
    END IF

    IF (PRESENT(R00)) LECS%RC(0,0) = R00
    IF (PRESENT(R10)) LECS%RC(1,0) = R10
    IF (PRESENT(R01)) LECS%RC(0,1) = R01
    IF (PRESENT(R11)) LECS%RC(1,1) = R11

    IF (PRESENT(C10)) LECS%CLO(0) = C10
    IF (PRESENT(C01)) LECS%CLO(1) = C01

    IF (PRESENT(C200)) LECS%CNLO(1) = C200
    IF (PRESENT(C201)) LECS%CNLO(2) = C201
    IF (PRESENT(C210)) LECS%CNLO(3) = C210
    IF (PRESENT(C211)) LECS%CNLO(4) = C211
    IF (PRESENT(C2T0)) LECS%CNLO(5) = C2T0
    IF (PRESENT(C2T1)) LECS%CNLO(6) = C2T1
    IF (PRESENT(C2B))  LECS%CNLO(7) = C2B

    IF (PRESENT(C400))  LECS%CN3LO(1)  = C400
    IF (PRESENT(C410))  LECS%CN3LO(2)  = C410
    IF (PRESENT(C401))  LECS%CN3LO(3)  = C401
    IF (PRESENT(C411))  LECS%CN3LO(4)  = C411
    IF (PRESENT(C4T0))  LECS%CN3LO(5)  = C4T0
    IF (PRESENT(C4T1))  LECS%CN3LO(6)  = C4T1
    IF (PRESENT(C4B0))  LECS%CN3LO(7)  = C4B0
    IF (PRESENT(C4B1))  LECS%CN3LO(8)  = C4B1
    IF (PRESENT(C4BB0)) LECS%CN3LO(9)  = C4BB0
    IF (PRESENT(C4BB1)) LECS%CN3LO(10) = C4BB1
    IF (PRESENT(C4Q))   LECS%CN3LO(11) = C4Q

    IF (PRESENT(CIT0)) LECS%CIT(0) = CIT0
    IF (PRESENT(CIT1)) LECS%CIT(1) = CIT1
    IF (PRESENT(CIT2)) LECS%CIT(2) = CIT2
    IF (PRESENT(CIT3)) LECS%CIT(3) = CIT3
    IF (PRESENT(CIT4)) LECS%CIT(4) = CIT4

    CALL SET_LECS_ORDER(LECS)
    LECS_SET = .TRUE.
  END SUBROUTINE SET_LECS_SINGLE

  !> \brief Return the path to LECS_EFT.DAT depending on working directory.
  FUNCTION GET_LECS_FILE_PATH() RESULT(lecs_file_path)
    USE OPERATING_SYSTEM_LINUX
    IMPLICIT NONE
    CHARACTER(LEN=255) :: lecs_file_path
    CHARACTER(LEN=256) :: current_dir

    current_dir = GET_CURRENT_WORKING_DIRECTORY()
    IF (INDEX(TRIM(current_dir), 'build') > 0) THEN
      lecs_file_path = FIND_FILE('lecs_eft.dat', '../')
    ELSE
      lecs_file_path = FIND_FILE('lecs_eft.dat', '.')
    END IF
  END FUNCTION GET_LECS_FILE_PATH

  !> \brief Set all LECs from file.
  SUBROUTINE SET_ALL_LECS
    IMPLICIT NONE
    INTEGER :: I, IOS
    CHARACTER(LEN=256) :: FILENAME
    INTEGER :: UNIT
    INTEGER :: ILB
    DOUBLE PRECISION :: RC(0:1,0:1), CLO(0:1), CNLO(7), CN3LO(11), CIT(0:4)

    IF (LECS_ALL_SET) RETURN

    ! Set the filename (you can change this as needed)
    
    FILENAME = GET_LECS_FILE_PATH()

    ! Open the file
    UNIT = 100
    OPEN(NEWUNIT=UNIT, FILE=FILENAME, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
    IF (ios /= 0) THEN
      PRINT *, 'Error opening file: ', TRIM(FILENAME),"   ->   ", ios
      STOP
    END IF

    ! Read the number of lines (number of LECS)
    READ(UNIT, *, IOSTAT=IOS) NMODELS
    IF (IOS /= 0) THEN
      PRINT *, 'Error reading number of lines from file.'
      CLOSE(UNIT)
      STOP
    END IF

    ! ALLOCATE THE ARRAY
    IF (ALLOCATED(LECS_ALL)) DEALLOCATE(LECS_ALL)
    ALLOCATE(LECS_ALL(NMODELS))

    ! Read each LECS_EFT_PLESS entry
    DO I = 1, NMODELS
      READ(UNIT, *, IOSTAT=IOS) ILB, RC(0,0), RC(1,0), RC(0,1), RC(1,1), CLO(1), CLO(0), CNLO, CN3LO, CIT
      LECS_ALL(I)%ILB = ILB
      LECS_ALL(I)%RC  = RC
      LECS_ALL(I)%CLO = CLO
      LECS_ALL(I)%CNLO = CNLO
      LECS_ALL(I)%CN3LO = CN3LO
      LECS_ALL(I)%CIT = CIT
      CALL SET_LECS_ORDER(LECS_ALL(I))
      IF (IOS /= 0) THEN
        PRINT *, 'Error reading LECS entry at line ', i
        CLOSE(UNIT)
        STOP
      END IF
    END DO
    CLOSE(UNIT)
    LECS_ALL_SET = .TRUE.
  END SUBROUTINE

  !> \brief Set the order of the LECs structure based on nonzero entries.
  SUBROUTINE SET_LECS_ORDER(LECS_IN)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: NUMERICAL_ZERO = 1.D-16
    TYPE(LECS_EFT_PLESS), INTENT(INOUT) :: LECS_IN
    INTEGER :: I, J
    INTEGER :: NLO = SIZE(LECS_IN%CNLO)
    INTEGER :: N3LO = SIZE(LECS_IN%CN3LO)

    LECS_IN%ORDER = 0
    DO I=1, NLO
      IF (DABS(LECS_IN%CNLO(I)) > NUMERICAL_ZERO) LECS_IN%ORDER = 1
    ENDDO

    DO J=1, N3LO
      IF (DABS(LECS_IN%CN3LO(J)) > NUMERICAL_ZERO) LECS_IN%ORDER = 3
    ENDDO

  END SUBROUTINE SET_LECS_ORDER

  !> \brief Set LECs by model index (ILB).
  !! \param[in] ILB Model index
  SUBROUTINE SET_LECS_ILB(ILB)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB
    IF (ILB <= 0 .AND. ILB > NMODELS) THEN
      WRITE(*,*) "THIS MODEL (ILB) IS NOT RECOGNIZED"
      STOP
    ENDIF
    LECS = LECS_ALL(ILB)
    LECS_SET = .TRUE.
  END SUBROUTINE SET_LECS_ILB

  !> \brief Set up operator matrices for the current channel.
  SUBROUTINE SET_OPERATORS()
    IMPLICIT NONE
    
    S12   = S12_OPERATOR(CHANNEL)
    LS    = LS_OPERATOR (CHANNEL)
    L2    = L2_OPERATOR (CHANNEL)
    T12   = T12_OPERATOR(CHANNEL)

  END SUBROUTINE SET_OPERATORS



  !> \brief Prepare the LECs and operators for a given model.
  !! \param[in] ILB Model index
  SUBROUTINE PREPARE(ILB)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: MINR = 1.D-1
    INTEGER, INTENT(IN) :: ILB
    LECS%ILB = ILB
    CALL SET_LECS_ILB(ILB)
    CALL SET_OPERATORS
    IF (SUM(ABS(LECS%RC)) < MINR) STOP "RC not set"
  END SUBROUTINE PREPARE

  SUBROUTINE EFT_PLESS_PW_LECS(L, S, J, T1Z, T2Z, R, VPW, LEMP, ILB, LOW_ENERGY_CONSTANTS)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L, S, J, T1Z, T2Z, LEMP
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VPW(2,2)
    INTEGER, OPTIONAL :: ILB
    TYPE(LECS_EFT_PLESS), OPTIONAL :: LOW_ENERGY_CONSTANTS

    INTEGER :: TZ, T, LS2(2,2)
    DOUBLE PRECISION :: RC
    TYPE(LECS_EFT_PLESS) :: LECS_

    IF (PRESENT(LOW_ENERGY_CONSTANTS) .AND. PRESENT(ILB)) THEN
      WRITE(*,*) "ERROR: BOTH LOW_ENERGY_CONSTANTS AND ILB ARE PRESENT. CHOOSE ONE."
      STOP
    ELSEIF (PRESENT(ILB)) THEN
      CALL EFT_PLESS_PW(ILB, L, S, J, T1Z, T2Z, R, VPW, LEMP)
    ELSEIF (.NOT.PRESENT(LOW_ENERGY_CONSTANTS)) THEN
        WRITE(*,*) "ERROR: LOW_ENERGY_CONSTANTS MUST BE PRESENT IF ILB IS NOT SET."
        STOP
    ELSE
      
      VPW = 0
      TZ = (T1Z+T2Z)/2
      T = T_FROM_L_S(L,S)

      BLOCK 
        TYPE(SCATTERING_CHANNEL) :: CHANNEL_NEW
        LOGICAL :: IS_NEW_CHANNEL
        CALL CHANNEL_NEW%SET(J, L, S, TZ)
        IS_NEW_CHANNEL = (CHANNEL_NEW /= CHANNEL)
        CHANNEL = CHANNEL_NEW
        
        IF (FIRST_CALL .OR. IS_NEW_CHANNEL) THEN
          CALL SET_OPERATORS
        ENDIF
      END BLOCK
      
      LECS_ = DIMENSIONAL_LECS_FROM_ADIMENSIONAL(LOW_ENERGY_CONSTANTS,S,T)
      RC = LECS_%RC(S,T)
      CALL SET_LECS_ORDER(LECS_)
      ORDER = LECS_%ORDER

      IF ( S==0 .AND. T==0 ) THEN
        SELECT CASE (ORDER)
        CASE (0)
          RETURN
        CASE (1)
          VPW = LECS_%CNLO(1) * EFT_RADIAL_1(R, RC) * I2
        CASE (3)
          VPW = LECS_%CNLO(1) * EFT_RADIAL_1(R, RC) * I2
          VPW = VPW + LECS_%CN3LO(1)  * EFT_RADIAL_4(R, RC) *I2 &
                    + LECS_%CN3LO(10) * EFT_RADIAL_7(RC) *L2
        CASE DEFAULT
          STOP "ERROR IN EFT_PLESS_PW_LECS:: S=0 AND T=0"
        END SELECT
      ENDIF

      IF ( S==1 .AND. T==0 ) THEN
        SELECT CASE (ORDER)
        CASE (0)
          VPW = LECS_%CLO(T) * I2
        CASE (1)
          VPW = LECS_%CLO(T) * I2
          VPW = VPW + LECS_%CNLO(2) * EFT_RADIAL_1(R, RC) * I2 &
                    + LECS_%CNLO(5) * EFT_RADIAL_2(R, RC) * S12 &
                    + LECS_%CNLO(7) * EFT_RADIAL_3(RC) * LS
        CASE (3)
          LS2 = MATMUL(LS, LS)
          VPW = LECS_%CLO(T) * I2
          VPW = VPW + LECS_%CNLO(2) * EFT_RADIAL_1(R, RC) * I2 &
                    + LECS_%CNLO(5) * EFT_RADIAL_2(R, RC) * S12 &
                    + LECS_%CNLO(7) * EFT_RADIAL_3(RC) * LS
          VPW = VPW + LECS_%CN3LO(2) * EFT_RADIAL_4(R, RC) * I2 &
                    + LECS_%CN3LO(5) * EFT_RADIAL_5(R, RC) * S12 &
                    + LECS_%CN3LO(7) * EFT_RADIAL_6(R, RC) * LS &
                    + LECS_%CN3LO(9) * EFT_RADIAL_7(RC) * LS2 &
                    + LECS_%CN3LO(11)* EFT_RADIAL_7(RC) * L2
        CASE DEFAULT
          STOP "ERROR IN EFT_PLESS_PW_LECS:: S=1 AND T=0"
        END SELECT
      ENDIF

      IF ( S==0 .AND. T==1 ) THEN
        SELECT CASE (ORDER)
        CASE (0)
          VPW = LECS_%CLO(T) * I2
        CASE (1)
          VPW = LECS_%CLO(T) * I2
          VPW = VPW + LECS_%CNLO(3) * EFT_RADIAL_1(R, RC) * I2 &
                    + LECS_%CIT(0)                * T12 * I2
        CASE (3)
          VPW = LECS_%CLO(T) * I2
          VPW = VPW + LECS_%CNLO(3)  * EFT_RADIAL_1(R, RC) * I2 &
                    + LECS_%CIT(0)                 * T12 * I2
          VPW = VPW + LECS_%CN3LO(3) * EFT_RADIAL_4(R, RC)     * I2 &
                    + LECS_%CN3LO(10)* EFT_RADIAL_7(RC)     * L2 &
                    + LECS_%CIT(1)   * EFT_RADIAL_1(R, RC)     * T12 * I2
        CASE DEFAULT
          STOP "ERROR IN EFT_PLESS_PW_LECS:: S=0 AND T=1"
        END SELECT
      ENDIF

      IF ( S==1 .AND. T==1 ) THEN
        SELECT CASE (ORDER)
        CASE (0)
          RETURN
        CASE (1)
          VPW =   LECS_%CNLO(4) * EFT_RADIAL_1(R, RC) *  I2 &
                + LECS_%CNLO(6) * EFT_RADIAL_2(R, RC) *  S12 &
                + LECS_%CNLO(7) * EFT_RADIAL_3(RC) *  LS &
                + LECS_%CIT(0)                *  T12 * I2
        CASE (3)
          LS2 = MATMUL(LS, LS)
          VPW =   LECS_%CNLO(4) * EFT_RADIAL_1(R, RC) *  I2 &
                + LECS_%CNLO(6) * EFT_RADIAL_2(R, RC) *  S12 &
                + LECS_%CNLO(7) * EFT_RADIAL_3(RC) *  LS &
                + LECS_%CIT(0)                *  T12 * I2
          VPW =   VPW &
                + LECS_%CN3LO(4) * EFT_RADIAL_4(R, RC)    * I2 &
                + LECS_%CN3LO(6) * EFT_RADIAL_5(R, RC)    * S12 &
                + LECS_%CN3LO(8) * EFT_RADIAL_6(R, RC)    * LS &
                + LECS_%CN3LO(9) * EFT_RADIAL_7(RC)    * LS2 &
                + LECS_%CN3LO(11)* EFT_RADIAL_7(RC)    * L2 &
                +(  LECS_%CIT(2) * EFT_RADIAL_1(R, RC)    * I2  + &
                    LECS_%CIT(3) * EFT_RADIAL_2(R, RC)    * S12 + &
                    LECS_%CIT(4) * EFT_RADIAL_3(RC)    * LS    &
                                                          ) * T12
        CASE DEFAULT
          STOP "ERROR IN EFT_PLESS_PW_LECS:: S=1 AND T=1"
        END SELECT
      ENDIF

      VPW = VPW * CR(R, RC)
      VPW = VPW * HTC

      IF (GET_CHANNEL_NCH(CHANNEL) == 1 ) THEN
        VPW(1,2) = 0.D0
        VPW(2,1) = 0.D0
        VPW(2,2) = 0.D0
      ENDIF
    ENDIF
  END SUBROUTINE EFT_PLESS_PW_LECS



  !> \brief Evaluate the EFT pionless potential matrix elements for given quantum numbers and radius.
  !! \ingroup eft_pless
  !! Computes the 2x2 potential matrix VPW for the specified channel and radius R using the current LECs.
  !! Handles all spin/isospin channels and EFT orders (LO, NLO, N3LO).
  !! \param[in] ILB Interaction label
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] J Total angular momentum
  !! \param[in] T1Z Isospin projection of nucleon 1
  !! \param[in] T2Z Isospin projection of nucleon 2
  !! \param[in] R Radial distance [fm]
  !! \param[out] VPW 2x2 potential matrix
  !! \param[in] LEMP Electromagnetic flag
  SUBROUTINE EFT_PLESS_PW(ILB, L, S, J, T1Z, T2Z, R, VPW, LEMP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB, L, S, J, T1Z, T2Z, LEMP
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT):: VPW(2,2)
    INTEGER :: T, TZ
    DOUBLE PRECISION :: RC
    INTEGER :: LS2(2,2)

    ! FITTED POTENTIAL
    IF     (ILB >= 35) THEN
      CALL EFT_PLESS_PW_FITTED(ILB, L, S, J, T1Z, T2Z, R, VPW, LEMP)
      RETURN
    ENDIF

    ! NORMAL POTENTIAL
    VPW = 0
    TZ = (T1Z+T2Z)/2
    T = MOD(MOD((L+S),2)+1,2)
    
    BLOCK
      TYPE(SCATTERING_CHANNEL) :: CHANNEL_NEW
      LOGICAL :: IS_NEW_CHANNEL
      CALL CHANNEL_NEW%SET(J, L, S, TZ)
      IS_NEW_CHANNEL = (CHANNEL /= CHANNEL_NEW)
      CHANNEL = CHANNEL_NEW
      
      IF (FIRST_CALL) THEN
        CALL SET_ALL_LECS
        CALL PREPARE(ILB)
        FIRST_CALL = .FALSE.
      ELSEIF (IS_NEW_CHANNEL) THEN
        CALL PREPARE(ILB)
      ENDIF
    END BLOCK

    RC = LECS%RC(S,T)
    ORDER = LECS%ORDER

    IF ( S==0 .AND. T==0 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RETURN
      CASE (1)
        VPW = LECS%CNLO(1) * EFT_RADIAL_1(R, RC) * I2
      CASE (3)
        VPW = LECS%CNLO(1) * EFT_RADIAL_1(R, RC) * I2
        VPW = VPW + LECS%CN3LO(1)  * EFT_RADIAL_4(R, RC) *I2 &
                  + LECS%CN3LO(10) * EFT_RADIAL_7(RC) *L2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=0 AND T=0"
      END SELECT
    ENDIF

    IF ( S==1 .AND. T==0 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        VPW = LECS%CLO(T) * I2
      CASE (1)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(2) * EFT_RADIAL_1(R, RC) * I2 &
                  + LECS%CNLO(5) * EFT_RADIAL_2(R, RC) * S12 &
                  + LECS%CNLO(7) * EFT_RADIAL_3(RC) * LS
      CASE (3)
        LS2 = MATMUL(LS, LS)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(2) * EFT_RADIAL_1(R, RC) * I2 &
                  + LECS%CNLO(5) * EFT_RADIAL_2(R, RC) * S12 &
                  + LECS%CNLO(7) * EFT_RADIAL_3(RC) * LS
        VPW = VPW + LECS%CN3LO(2) * EFT_RADIAL_4(R, RC) * I2 &
                  + LECS%CN3LO(5) * EFT_RADIAL_5(R, RC) * S12 &
                  + LECS%CN3LO(7) * EFT_RADIAL_6(R, RC) * LS &
                  + LECS%CN3LO(9) * EFT_RADIAL_7(RC) * LS2 &
                  + LECS%CN3LO(11)* EFT_RADIAL_7(RC) * L2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=1 AND T=0"
      END SELECT
    ENDIF

    IF ( S==0 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        VPW = LECS%CLO(T) * I2
      CASE (1)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(3) * EFT_RADIAL_1(R, RC) * I2 &
                  + LECS%CIT(0)                * T12 * I2
      CASE (3)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(3)  * EFT_RADIAL_1(R, RC) * I2 &
                  + LECS%CIT(0)                 * T12 * I2
        VPW = VPW + LECS%CN3LO(3) * EFT_RADIAL_4(R, RC)     * I2 &
                  + LECS%CN3LO(10)* EFT_RADIAL_7(RC)     * L2 &
                  + LECS%CIT(1)   * EFT_RADIAL_1(R, RC)     * T12 * I2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=0 AND T=1"
      END SELECT
    ENDIF

    IF ( S==1 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RETURN
      CASE (1)
        VPW =   LECS%CNLO(4) * EFT_RADIAL_1(R, RC) *  I2 &
              + LECS%CNLO(6) * EFT_RADIAL_2(R, RC) *  S12 &
              + LECS%CNLO(7) * EFT_RADIAL_3(RC) *  LS &
              + LECS%CIT(0)                *  T12 * I2
      CASE (3)
        LS2 = MATMUL(LS, LS)
        VPW =   LECS%CNLO(4) * EFT_RADIAL_1(R, RC) *  I2 &
              + LECS%CNLO(6) * EFT_RADIAL_2(R, RC) *  S12 &
              + LECS%CNLO(7) * EFT_RADIAL_3(RC) *  LS &
              + LECS%CIT(0)                *  T12 * I2
        VPW =   VPW &
              + LECS%CN3LO(4) * EFT_RADIAL_4(R, RC)    * I2 &
              + LECS%CN3LO(6) * EFT_RADIAL_5(R, RC)    * S12 &
              + LECS%CN3LO(8) * EFT_RADIAL_6(R, RC)    * LS &
              + LECS%CN3LO(9) * EFT_RADIAL_7(RC)    * LS2 &
              + LECS%CN3LO(11)* EFT_RADIAL_7(RC)    * L2 &
              +(  LECS%CIT(2) * EFT_RADIAL_1(R, RC)    * I2  + &
                  LECS%CIT(3) * EFT_RADIAL_2(R, RC)    * S12 + &
                  LECS%CIT(4) * EFT_RADIAL_3(RC)    * LS    &
                                                        ) * T12
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW:: S=1 AND T=1"
      END SELECT
    ENDIF

    VPW = VPW * CR(R, RC)
    VPW = VPW * HTC

    IF (GET_CHANNEL_NCH(CHANNEL) == 1 ) THEN
      VPW(1,2) = 0.D0
      VPW(2,1) = 0.D0
      VPW(2,2) = 0.D0
    ENDIF
    RETURN
  END SUBROUTINE EFT_PLESS_PW


  PURE ELEMENTAL FUNCTION CR(R, RC) RESULT(CR_OUT)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: CR_OUT
    CR_OUT = DEXP(-R**2/RC**2)/(PI**(3.D0/2.D0)*RC**3)
  END FUNCTION CR

  PURE ELEMENTAL FUNCTION EFT_RADIAL_1(R, RC) RESULT(F1)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F1
    F1 = (6*RC**2 - 4*R**2)/(RC**4)
  END FUNCTION EFT_RADIAL_1

  PURE ELEMENTAL FUNCTION EFT_RADIAL_2(R, RC) RESULT(F2)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F2
    F2 =-4*R**2/RC**4
  END FUNCTION EFT_RADIAL_2

  PURE ELEMENTAL FUNCTION EFT_RADIAL_3(RC) RESULT(F3)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RC
    DOUBLE PRECISION :: F3
    F3 = 2.D0/RC**2
  END FUNCTION EFT_RADIAL_3

  PURE ELEMENTAL FUNCTION EFT_RADIAL_4(R, RC) RESULT(F4)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F4
    F4 = 4.D0*( 4*R**4 - 20*RC**2*R**2 + 15*RC**4 )/RC**8
  END FUNCTION EFT_RADIAL_4

  PURE ELEMENTAL FUNCTION EFT_RADIAL_5(R, RC) RESULT(F5)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F5
    F5 = 8*R**2*( 2*R**2 - 7*RC**2 )/RC**8
  END FUNCTION EFT_RADIAL_5

  PURE ELEMENTAL FUNCTION EFT_RADIAL_6(R, RC) RESULT(F6)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: F6
    F6 = ( 20*RC**2 - 8*R**2 )/RC**6
  END FUNCTION EFT_RADIAL_6

  PURE ELEMENTAL FUNCTION EFT_RADIAL_7(RC) RESULT(F7)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RC
    DOUBLE PRECISION :: F7
    F7 =-EFT_RADIAL_3(RC)**2
  END FUNCTION EFT_RADIAL_7


  !> \brief Transform LECs from the operator basis to the S,T-coupled basis.
  !! \ingroup eft_pless
  !!
  !! This function converts a LECS_EFT_PLESS structure containing Low-Energy Constants (LECs)
  !! in the operator basis (used for constructing the EFT pionless potential) into the S,T-coupled basis
  !! (used for fitting or storage).
  !!
  !! The transformation applies the appropriate linear combinations to the CNLO, CN3LO,
  !! and CIT arrays, following the mapping between operator and S,T representations.
  !!
  !! \param[in]  LECS_OLD  LECs in the operator basis
  !! \return     ST_LECS   LECs in the S,T-coupled basis, ready for fitting or storage
  FUNCTION LECS_TO_ST_LECS(LECS_OLD) RESULT(ST_LECS)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_OLD
    TYPE(LECS_EFT_PLESS) :: ST_LECS
    DOUBLE PRECISION :: C(7)
    DOUBLE PRECISION :: D(11)
    DOUBLE PRECISION :: CIT(0:4)

    C   = LECS_OLD%CNLO
    D   = LECS_OLD%CN3LO
    CIT = LECS_OLD%CIT

    ST_LECS%ILB       = LECS_OLD%ILB
    ST_LECS%RC        = LECS_OLD%RC
    ST_LECS%CLO       = LECS_OLD%CLO
    ST_LECS%CNLO(1)   = C(1) - 3*C(2) - 3*C(3) + 9*C(4)
    ST_LECS%CNLO(2)   = C(1) - 3*C(2) + C(3) - 3*C(4)
    ST_LECS%CNLO(3)   = C(1) + C(2) - 3*C(3) - 3*C(4)
    ST_LECS%CNLO(4)   = C(1) + C(2) + C(3) + C(4)
    ST_LECS%CNLO(5)   = C(5) - 3*C(6)
    ST_LECS%CNLO(6)   = C(5) + C(6)
    ST_LECS%CNLO(7)   = C(7)
    ST_LECS%CN3LO( 1) = D(1) - 3*D(2) - 3*D(3) + 9*D(4)
    ST_LECS%CN3LO( 2) = D(1) - 3*D(2) + D(3) - 3*D(4)
    ST_LECS%CN3LO( 3) = D(1) + D(2) - 3*D(3) - 3*D(4)
    ST_LECS%CN3LO( 4) = D(1) + D(2) + D(3) + D(4)
    ST_LECS%CN3LO( 5) = D(5) - 3*D(6)
    ST_LECS%CN3LO( 6) = D(5) + D(6)
    ST_LECS%CN3LO( 7) = D(7) - 3*D(8)
    ST_LECS%CN3LO( 8) = D(7) + D(8)
    ST_LECS%CN3LO( 9) = D(9)
    ST_LECS%CN3LO(10) = D(10) - 3*D(11)
    ST_LECS%CN3LO(11) = D(10) + D(11)
    ST_LECS%CIT(0)    = CIT(0)
    ST_LECS%CIT(1)    = CIT(1) - 3*CIT(2)
    ST_LECS%CIT(2)    = CIT(1) + CIT(2)
    ST_LECS%CIT(3)    = CIT(3)
    ST_LECS%CIT(4)    = CIT(4)
  END FUNCTION LECS_TO_ST_LECS

  !> \brief Transform LECs from the S,T-coupled basis to the operator basis.
  !! \ingroup eft_pless
  !!
  !! This function converts a LECS_EFT_PLESS structure containing Low-Energy Constants (LECs)
  !! in the S,T-coupled basis (used for fitting or storage) into the operator basis
  !! used for constructing the EFT pionless potential.
  !!
  !! The transformation applies the appropriate linear combinations to the CNLO, CN3LO,
  !! and CIT arrays, following the mapping between S,T and operator representations.
  !!
  !! \param[in]  ST_LECS  LECs in the S,T-coupled basis
  !! \return     OP_LECS  LECs in the operator basis, ready for use in the potential
  FUNCTION ST_LECTS_TO_LECS(ST_LECS) RESULT(OP_LECS)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: ST_LECS
    TYPE(LECS_EFT_PLESS) :: OP_LECS
    DOUBLE PRECISION :: C1(7)
    DOUBLE PRECISION :: D1(11)
    DOUBLE PRECISION :: CIT(0:4)

    C1   = ST_LECS%CNLO
    D1   = ST_LECS%CN3LO
    CIT  = ST_LECS%CIT

    OP_LECS%ILB= ST_LECS%ILB
    OP_LECS%RC = ST_LECS%RC
    OP_LECS%CLO = ST_LECS%CLO
    OP_LECS%CNLO(1)   = (C1(1) + 3*C1(2) + 3*C1(3) + 9*C1(4))/16.D0
    OP_LECS%CNLO(2)   = (-C1(1) - 3*C1(2) + C1(3) + 3*C1(4))/16.D0
    OP_LECS%CNLO(3)   = (-C1(1) + C1(2) - 3*C1(3) + 3*C1(4))/16.D0
    OP_LECS%CNLO(4)   = (C1(1) - C1(2) - C1(3) + C1(4))/16.D0
    OP_LECS%CNLO(5)   = (C1(5) + 3*C1(6))/4.D0
    OP_LECS%CNLO(6)   = (-C1(5) + C1(6))/4.D0
    OP_LECS%CNLO(7)   = C1(7)
    OP_LECS%CN3LO(1)  = (D1(1) + 3*(D1(2) + D1(3) + 3*D1(4)))/16.D0
    OP_LECS%CN3LO(2)  = (-D1(1) - 3*D1(2) + D1(3) + 3*D1(4))/16.D0
    OP_LECS%CN3LO(3)  = (-D1(1) + D1(2) - 3*D1(3) + 3*D1(4))/16.D0
    OP_LECS%CN3LO(4)  = (D1(1) - D1(2) - D1(3) + D1(4))/16.D0
    OP_LECS%CN3LO(5)  = (D1(5) + 3*D1(6))/4.D0
    OP_LECS%CN3LO(6)  = (-D1(5) + D1(6))/4.D0
    OP_LECS%CN3LO(7)  = (D1(7) + 3*D1(8))/4.D0
    OP_LECS%CN3LO(8)  = (-D1(7) + D1(8))/4.D0
    OP_LECS%CN3LO(9)  = D1(9)
    OP_LECS%CN3LO(10) = (D1(10) + 3*D1(11))/4.D0
    OP_LECS%CN3LO(11) = (-D1(10) + D1(11))/4.D0
    OP_LECS%CIT(0)    = CIT(0)
    OP_LECS%CIT(1)    = (CIT(1) + 3*CIT(2))/4.D0
    OP_LECS%CIT(2)    = (-CIT(1) + CIT(2))/4.D0
    OP_LECS%CIT(3)    = CIT(3)
    OP_LECS%CIT(4)    = CIT(4)
  END FUNCTION ST_LECTS_TO_LECS

  !> \brief Retrieve the Low-Energy Constants (LECs) for a given model index.
  !! \ingroup eft_pless
  !!
  !! Returns the LECS_EFT_PLESS structure corresponding to the specified model index (ILB).
  !! If the LECs have not been set, they are loaded from file.
  !! If the index is invalid, the routine stops with an error.
  !!
  !! \param[in] ILB Model index for which to retrieve the LECs
  !! \return LECS_OUT Structure containing the LECs for the specified model
  FUNCTION GET_LECS(ILB) RESULT(LECS_OUT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB
    TYPE(LECS_EFT_PLESS) :: LECS_OUT

    IF (.NOT.LECS_SET) CALL SET_ALL_LECS

    IF (ILB <= 0 .AND. ILB > NMODELS) THEN
      WRITE(*,*) "THIS MODEL (ILB) IS NOT RECOGNIZED"
      STOP
    ENDIF
    LECS_OUT = LECS_ALL(ILB)
  END FUNCTION GET_LECS

  !> @brief Prints the contents of a LECS_EFT_PLESS structure.
  !! \ingroup eft_pless
  !>
  !> @details It prints the potential submodel (ILB), order of the potential (ORDER),
  !> the cutoff parameters (RC), and the coefficients for the LO, NLO, and N3LO potentials.
  !>
  !> @param[in] LECS_IN Input LECS_EFT_PLESS structure.
  SUBROUTINE PRINT_LECS(LECS_IN)
    IMPLICIT NONE
    CLASS(LECS_EFT_PLESS), INTENT(IN) :: LECS_IN
    INTEGER :: I, J
    CHARACTER(LEN=16) :: ORD

    SELECT CASE (LECS_IN%ORDER)
    CASE (0)
      ORD = "LO"
    CASE (1)
      ORD = "NLO"
    CASE (3)
      ORD = "N3LO"
    CASE DEFAULT
      ORD = "UNKNOWN"
    END SELECT

    PRINT *, "LECS_EFT_PLESS:"
    PRINT *, "  ILB   = ", LECS_IN%ILB
    PRINT *, "  ORDER =          ", ORD
    PRINT *, "  DIMENSIONAL: ", .NOT.(LECS_IN%ADIMENSIONAL)
    PRINT *, "  RC:"
    DO I = 0, 1
      WRITE(*,'(A,2F12.8)') "    ", LECS_IN%RC(I,0), LECS_IN%RC(I,1)
    END DO
    PRINT *, "  CLO:"
    WRITE(*,'(A,2ES18.8E2)') "    ", LECS_IN%CLO(1), LECS_IN%CLO(0)
    PRINT *, "  CNLO:"
    WRITE(*,'(A,7ES18.8E2)') "    ", (LECS_IN%CNLO(J), J=1,7)
    PRINT *, "  CN3LO:"
    WRITE(*,'(A,7ES18.8E2)') "    ", (LECS_IN%CN3LO(J), J=1,7)
    WRITE(*,'(A,4ES18.8E2)') "    ", (LECS_IN%CN3LO(J), J=8,11)
    PRINT *, "  CIT:"
    WRITE(*,'(A,5ES18.8E2)') "    ", (LECS_IN%CIT(J), J=0,4)
  END SUBROUTINE PRINT_LECS


  !> @brief It returns all the possible radial function multiplying the LECS in 
  !> the EFT pionless potential.
  !! \ingroup eft_pless
  !>
  !> @details The possible radial functions are of 8 kind. At LO there is only one radial function (I=0),
  !> at NLO there are 5 radial functions (I=0,1,2,3,4), and at N3LO there are 8 radial functions (I=0,1,2,3,4,5,6,7).
  !>
  !> @param[in] R_ARRAY Array of radial points (in fm).
  !> @param[in] RC Cutoff parameters for the radial functions (2x2 array).
  !> @param[out] FUNCTIONS Structure containing the radial functions.
  !> @param[in] ORDER_POTENTIAL (optional) Order of the potential (0 for LO, 1 for NLO, 3 for N3LO).
  !> If not specified, defaults to 3 (N3LO).
  SUBROUTINE GET_EFT_RADIAL_FUNCTIONS(R_ARRAY, RC, FUNCTIONS, ORDER_POTENTIAL)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RC(0:1,0:1)
    DOUBLE PRECISION, INTENT(IN) :: R_ARRAY(:)
    TYPE(EFT_RADIAL_FUNCTIONS), INTENT(OUT) :: FUNCTIONS
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER_POTENTIAL

    INTEGER :: NR, I, IMAX, S, T

    IF (PRESENT(ORDER_POTENTIAL)) THEN
      FUNCTIONS%ORDER = ORDER_POTENTIAL
    ELSE
      FUNCTIONS%ORDER = 3  ! Default to N3LO if not specified
    END IF

    NR = SIZE(R_ARRAY)
    IF (NR < 1) THEN
      PRINT *, "ERROR: No radial points provided for EFT radial functions."
      RETURN
    END IF

    IMAX = 0
    IF (FUNCTIONS%ORDER == 1) THEN
      IMAX = 4
    ELSEIF (FUNCTIONS%ORDER == 3) THEN
      IMAX = 7
    END IF

    IF (ALLOCATED(FUNCTIONS%FR_I)) DEALLOCATE(FUNCTIONS%FR_I)
    ALLOCATE(FUNCTIONS%FR_I(0:1, 0:1, 0:IMAX, NR))
    FUNCTIONS%RC = RC
    FUNCTIONs%FR_I = 0.D0
    DO S = 0, 1
    DO T = 0, 1
      IF (RC(S,T) <= 0.D0) CYCLE
      DO I=0, IMAX
        FUNCTIONS%FR_I(S,T,I,:) = EFT_RADIAL(I, R_ARRAY, RC(S,T))
      END DO
    ENDDO
    ENDDO
  END SUBROUTINE GET_EFT_RADIAL_FUNCTIONS

  PURE ELEMENTAL FUNCTION EFT_RADIAL(I, R, RC) RESULT(FI)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I
    DOUBLE PRECISION, INTENT(IN) :: R, RC
    DOUBLE PRECISION :: FI

    SELECT CASE (I)
    CASE (0)
      FI = 1.d0 ! The CR function for I=0
    CASE (1)
      FI = EFT_RADIAL_1(R, RC)
    CASE (2)
      FI = EFT_RADIAL_2(R, RC)
    CASE (3)
      FI = EFT_RADIAL_3(RC)
    CASE (4)
      FI = EFT_RADIAL_4(R, RC)
    CASE (5)
      FI = EFT_RADIAL_5(R, RC)
    CASE (6)
      FI = EFT_RADIAL_6(R, RC)
    CASE (7)
      FI = EFT_RADIAL_7(RC)
    CASE DEFAULT
      FI = 0.D0
    END SELECT
    FI = FI * CR(R, RC)
  END FUNCTION EFT_RADIAL

  !> \brief Combine operator radial matrix elements with LECs to form the full potential.
  !! \ingroup eft_pless
  !! \param[in] CHANNELS Array of channels
  !! \param[in] FR_MATRIX_EL Operator radial matrix elements with dimensions (NOPERATOR, NCHANNELS, NR, NALPHAL, NALPHAR)
  !!                        <left| FR_i(r) |right> for each operator i, channel, subchannel
  !! \param[in] LECS_IN LECs structure
  !! \param[out] POTENTIAL_OUT Output potential matrix elements with dimensions (NCHANNELS, NOPERATOR, NALPHAL, NALPHAR)
  SUBROUTINE COMBINE_POTENTIAL(CHANNELS, FR_MATRIX_EL, LECS_IN, POTENTIAL_OUT)
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNELS(:)
    DOUBLE PRECISION, INTENT(IN) :: FR_MATRIX_EL(:,:,:,:,:)
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_IN
    DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT) :: POTENTIAL_OUT(:,:,:,:)
    INTEGER :: ORDER_
    INTEGER :: S, T, LSC(2,2), LS2C(2,2), L2C(2,2), T12C
    DOUBLE PRECISION :: S12C(2,2)
    INTEGER :: ICH, NCHANNELS, NALPHAL, NALPHAR, NEQ, NEQ_MAX, NNLL, NNLR
    INTEGER :: IL, IR, AR, AL, NL, NR
    INTEGER :: SL, SR, TL, TR

    IF (LECS_IN%ADIMENSIONAL) THEN
      WRITE(*,*) "ERROR IN COMBINE_POTENTIAL:: LECs are in dimensional form, not suitable for this routine."
      STOP
    END IF

    ORDER_ = LECS_IN%ORDER
    IF (ORDER_ /= 0 .AND. ORDER_ /= 1 .AND. ORDER_ /= 3) THEN
      STOP "ERROR IN COMBINE_POTENTIAL:: ORDER must be 0, 1, or 3"
    END IF
    CALL REALLOCATE(POTENTIAL_OUT, SIZE(CHANNELS), SIZE(FR_MATRIX_EL, 3), SIZE(FR_MATRIX_EL, 4), SIZE(FR_MATRIX_EL, 5))

    NCHANNELS = SIZE(CHANNELS)
    NALPHAL = SIZE(FR_MATRIX_EL, 4)
    NALPHAR = SIZE(FR_MATRIX_EL, 5)
    NEQ_MAX = 2
    NNLL = NALPHAL / NEQ_MAX
    NNLR = NALPHAR / NEQ_MAX
    DO ICH = 1, NCHANNELS
      S12C = S12_OPERATOR(CHANNELS(ICH))
      LSC = LS_OPERATOR(CHANNELS(ICH))
      LS2C = MATMUL(LSC, LSC)
      L2C  = L2_OPERATOR(CHANNELS(ICH))
      T12C = T12_OPERATOR(CHANNELS(ICH))
      NEQ = GET_CHANNEL_NCH(CHANNELS(ICH))

      IR = 0
      IL = 0
      DO AL = 1, NEQ
      DO AR = 1, NEQ
      DO NL = 1, NNLL
      DO NR = 1, NNLR
        IL = (AL-1)*NNLL + NL
        IR = (AR-1)*NNLR + NR

        SL = GET_CHANNEL_S(CHANNELS(ICH), AL)
        TL = GET_CHANNEL_T(CHANNELS(ICH), AL)
        SR = GET_CHANNEL_S(CHANNELS(ICH), AR)
        TR = GET_CHANNEL_T(CHANNELS(ICH), AR)
        IF ( SL /= SR .OR. TL /= TR ) CYCLE
        S = SL
        T = TL

        IF ( S==0 .AND. T==0 ) THEN
          SELECT CASE (ORDER_)
          CASE (0)
            POTENTIAL_OUT(ICH,:,IL,IR) = 0.D0
          CASE (1)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CNLO(1) * FR_MATRIX_EL(2,ICH,:,IL,IR) * I2(AL,AR)
          CASE (3)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CNLO(1) * FR_MATRIX_EL(2,ICH,:,IL,IR) * I2(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) = POTENTIAL_OUT(ICH,:,IL,IR) &
                      + LECS_IN%CN3LO(1)  * FR_MATRIX_EL(5,ICH,:,IL,IR) * I2(AL,AR) &
                      + LECS_IN%CN3LO(10) * FR_MATRIX_EL(8,ICH,:,IL,IR) * L2C(AL,AR)
          CASE DEFAULT
            WRITE(*,*) "ERROR IN EFT_PLESS_PW:: S=0 AND T=0, ORDER: ", ORDER_
            STOP
          END SELECT
        ENDIF

        IF ( S==1 .AND. T==0 ) THEN
          SELECT CASE (ORDER_)
          CASE (0)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CLO(T) * FR_MATRIX_EL(1,ICH,:,IL,IR) * I2(AL,AR)
          CASE (1)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CLO(T) * FR_MATRIX_EL(1,ICH,:,IL,IR) * I2(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) = POTENTIAL_OUT(ICH,:,IL,IR) &
                      + LECS_IN%CNLO(2) * FR_MATRIX_EL(2,ICH,:,IL,IR) * I2(AL,AR) &
                      + LECS_IN%CNLO(5) * FR_MATRIX_EL(3,ICH,:,IL,IR) * S12C(AL,AR) &
                      + LECS_IN%CNLO(7) * FR_MATRIX_EL(4,ICH,:,IL,IR) * LSC(AL,AR)
          CASE (3)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CLO(T) * FR_MATRIX_EL(1,ICH,:,IL,IR) * I2(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) = POTENTIAL_OUT(ICH,:,IL,IR) &
                      + LECS_IN%CNLO(2) * FR_MATRIX_EL(2,ICH,:,IL,IR) * I2(AL,AR) &
                      + LECS_IN%CNLO(5) * FR_MATRIX_EL(3,ICH,:,IL,IR) * S12C(AL,AR) &
                      + LECS_IN%CNLO(7) * FR_MATRIX_EL(4,ICH,:,IL,IR) * LSC(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) = POTENTIAL_OUT(ICH,:,IL,IR) &
                      + LECS_IN%CN3LO(2) * FR_MATRIX_EL(5,ICH,:,IL,IR) * I2(AL,AR) &
                      + LECS_IN%CN3LO(5) * FR_MATRIX_EL(6,ICH,:,IL,IR) * S12C(AL,AR) &
                      + LECS_IN%CN3LO(7) * FR_MATRIX_EL(7,ICH,:,IL,IR) * LSC(AL,AR) &
                      + LECS_IN%CN3LO(9) * FR_MATRIX_EL(8,ICH,:,IL,IR) * LS2C(AL,AR) &
                      + LECS_IN%CN3LO(11)* FR_MATRIX_EL(8,ICH,:,IL,IR) * L2C(AL,AR)
          CASE DEFAULT
            WRITE(*,*) "ERROR IN EFT_PLESS_PW:: S=1 AND T=0, ORDER: ", ORDER_
            STOP
          END SELECT
        ENDIF

        IF ( S==0 .AND. T==1 ) THEN
          SELECT CASE (ORDER_)
          CASE (0)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CLO(T) * FR_MATRIX_EL(1,ICH,:,IL,IR) * I2(AL,AR)
          CASE (1)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CLO(T) * FR_MATRIX_EL(1,ICH,:,IL,IR) * I2(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) = POTENTIAL_OUT(ICH,:,IL,IR) &
                      + LECS_IN%CNLO(3)    * FR_MATRIX_EL(2,ICH,:,IL,IR) * I2(AL,AR) &
                      + LECS_IN%CIT(0)     * FR_MATRIX_EL(1,ICH,:,IL,IR) * T12C * I2(AL,AR)
          CASE (3)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CLO(T) * FR_MATRIX_EL(1,ICH,:,IL,IR) * I2(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) = POTENTIAL_OUT(ICH,:,IL,IR) &
                      + LECS_IN%CNLO(3)     * FR_MATRIX_EL(2,ICH,:,IL,IR) * I2(AL,AR) &
                      + LECS_IN%CIT(0)      * FR_MATRIX_EL(1,ICH,:,IL,IR) * T12C * I2(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) = POTENTIAL_OUT(ICH,:,IL,IR) &
                      + LECS_IN%CN3LO(3) * FR_MATRIX_EL(5,ICH,:,IL,IR)     * I2(AL,AR) &
                      + LECS_IN%CN3LO(10)* FR_MATRIX_EL(8,ICH,:,IL,IR)     * L2C(AL,AR) &
                      + LECS_IN%CIT(1)   * FR_MATRIX_EL(2,ICH,:,IL,IR)     * T12C * I2(AL,AR)
          CASE DEFAULT
            WRITE(*,*) "ERROR IN EFT_PLESS_PW:: S=0 AND T=1, ORDER: ", ORDER_
            STOP
          END SELECT
        ENDIF

        IF ( S==1 .AND. T==1 ) THEN
          SELECT CASE (ORDER_)
          CASE (0)
            POTENTIAL_OUT(ICH,:,IL,IR) = 0.D0
          CASE (1)
            POTENTIAL_OUT(ICH,:,IL,IR) = &
                  + LECS_IN%CNLO(4)         * FR_MATRIX_EL(2,ICH,:,IL,IR) *  I2(AL,AR) &
                  + LECS_IN%CNLO(6)         * FR_MATRIX_EL(3,ICH,:,IL,IR) *  S12C(AL,AR) &
                  + LECS_IN%CNLO(7)         * FR_MATRIX_EL(4,ICH,:,IL,IR) *  LSC(AL,AR) &
                  + LECS_IN%CIT(0)          * FR_MATRIX_EL(1,ICH,:,IL,IR) *  T12C * I2(AL,AR)
          CASE (3)
            POTENTIAL_OUT(ICH,:,IL,IR) = LECS_IN%CNLO(4) * FR_MATRIX_EL(2,ICH,:,IL,IR) *  I2(AL,AR) &
                  + LECS_IN%CNLO(6)         * FR_MATRIX_EL(3,ICH,:,IL,IR) *  S12C(AL,AR) &
                  + LECS_IN%CNLO(7)         * FR_MATRIX_EL(4,ICH,:,IL,IR) *  LSC(AL,AR) &
                  + LECS_IN%CIT(0)          * FR_MATRIX_EL(1,ICH,:,IL,IR) *  T12C * I2(AL,AR)
            POTENTIAL_OUT(ICH,:,IL,IR) =   POTENTIAL_OUT(ICH,:,IL,IR) &
                  + LECS_IN%CN3LO(4) * FR_MATRIX_EL(5,ICH,:,IL,IR)    * I2(AL,AR) &
                  + LECS_IN%CN3LO(6) * FR_MATRIX_EL(6,ICH,:,IL,IR)    * S12C(AL,AR) &
                  + LECS_IN%CN3LO(8) * FR_MATRIX_EL(7,ICH,:,IL,IR)    * LSC(AL,AR) &
                  + LECS_IN%CN3LO(9) * FR_MATRIX_EL(8,ICH,:,IL,IR)    * LS2C(AL,AR) &
                  + LECS_IN%CN3LO(11)* FR_MATRIX_EL(8,ICH,:,IL,IR)    * L2C(AL,AR) &
                  +(  LECS_IN%CIT(2) * FR_MATRIX_EL(2,ICH,:,IL,IR)    * I2(AL,AR)  + &
                      LECS_IN%CIT(3) * FR_MATRIX_EL(3,ICH,:,IL,IR)    * S12C(AL,AR) + &
                      LECS_IN%CIT(4) * FR_MATRIX_EL(4,ICH,:,IL,IR)    * LSC(AL,AR) &
                                                            ) * T12C
          CASE DEFAULT
            WRITE(*,*) "ERROR IN EFT_PLESS_PW:: S=1 AND T=1, ORDER: ", ORDER_
            STOP
          END SELECT
        ENDIF
      ENDDO ! IR
      ENDDO ! IL
      ENDDO ! AR
      ENDDO ! AL
    ENDDO ! ICH
    POTENTIAL_OUT = POTENTIAL_OUT * HTC
  END SUBROUTINE COMBINE_POTENTIAL


  FUNCTION ADIMENSIONAL_LECS_FROM_DIMENSIONAL(LECS_IN, S, T) RESULT(LECS_OUT)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_IN
    INTEGER, INTENT(IN) :: S, T
    TYPE(LECS_EFT_PLESS) :: LECS_OUT

    DOUBLE PRECISION :: RC

    LECS_OUT = LECS_IN
    IF (LECS_IN%ADIMENSIONAL) RETURN
    RC = LECS_IN%RC(S,T)

    LECS_OUT%CLO = LECS_IN%CLO / RC**2
    LECS_OUT%CNLO = LECS_IN%CNLO / RC**4
    LECS_OUT%CN3LO = LECS_IN%CN3LO / RC**6
    LECS_OUT%CIT(0) = LECS_IN%CIT(0) / RC**2
    LECS_OUT%CIT(1:4) = LECS_IN%CIT(1:4) / RC**4
    LECS_OUT%ADIMENSIONAL = .TRUE.
  END FUNCTION ADIMENSIONAL_LECS_FROM_DIMENSIONAL

  FUNCTION DIMENSIONAL_LECS_FROM_ADIMENSIONAL(LECS_IN, S, T) RESULT(LECS_OUT)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_IN
    INTEGER, INTENT(IN) :: S, T
    TYPE(LECS_EFT_PLESS) :: LECS_OUT

    DOUBLE PRECISION :: RC

    LECS_OUT = LECS_IN
    IF (.NOT.LECS_IN%ADIMENSIONAL) RETURN
    RC = LECS_IN%RC(S,T)

    LECS_OUT%CLO = LECS_IN%CLO * RC**2
    LECS_OUT%CNLO = LECS_IN%CNLO * RC**4
    LECS_OUT%CN3LO = LECS_IN%CN3LO * RC**6
    LECS_OUT%CIT(0) = LECS_IN%CIT(0) * RC**2
    LECS_OUT%CIT(1:4) = LECS_IN%CIT(1:4) * RC**4
    LECS_OUT%ADIMENSIONAL = .FALSE.
  END FUNCTION DIMENSIONAL_LECS_FROM_ADIMENSIONAL

  SUBROUTINE CONVERT_LECS_TO_ADIMENSIONAL(THIS, S, T)
    IMPLICIT NONE
    CLASS(LECS_EFT_PLESS), INTENT(INOUT) :: THIS
    INTEGER, INTENT(IN) :: S, T

    IF (THIS%ADIMENSIONAL) RETURN
    THIS%CLO = THIS%CLO / THIS%RC(S,T)**2
    THIS%CNLO = THIS%CNLO / THIS%RC(S,T)**4
    THIS%CN3LO = THIS%CN3LO / THIS%RC(S,T)**6
    THIS%CIT(0) = THIS%CIT(0) / THIS%RC(S,T)**2
    THIS%CIT(1:4) = THIS%CIT(1:4) / THIS%RC(S,T)**4
    THIS%ADIMENSIONAL = .TRUE.
  END SUBROUTINE CONVERT_LECS_TO_ADIMENSIONAL

  SUBROUTINE CONVERT_LECS_TO_DIMENSIONAL(THIS, S, T)
    IMPLICIT NONE
    CLASS(LECS_EFT_PLESS), INTENT(INOUT) :: THIS
    INTEGER, INTENT(IN) :: S, T

    IF (.NOT.THIS%ADIMENSIONAL) RETURN
    THIS%CLO = THIS%CLO * THIS%RC(S,T)**2
    THIS%CNLO = THIS%CNLO * THIS%RC(S,T)**4
    THIS%CN3LO = THIS%CN3LO * THIS%RC(S,T)**6
    THIS%CIT(0) = THIS%CIT(0) * THIS%RC(S,T)**2
    THIS%CIT(1:4) = THIS%CIT(1:4) * THIS%RC(S,T)**4
    THIS%ADIMENSIONAL = .FALSE.
  END SUBROUTINE CONVERT_LECS_TO_DIMENSIONAL




  SUBROUTINE EFT_PLESS_PW_FITTED(ILB, L, S, J, T1Z, T2Z, R, VPW, LEMP)
    USE AV18
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB, L, S, J, T1Z, T2Z, LEMP
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VPW(2,2)
    
    INTEGER :: TZ, T, LS2(2,2)
    DOUBLE PRECISION :: VEM(14)
    DOUBLE PRECISION :: RC

    VPW = 0
    TZ = (T1Z+T2Z)/2
    T = MOD(MOD((L+S),2)+1,2)
    
    BLOCK
      TYPE(SCATTERING_CHANNEL) :: CHANNEL_NEW
      LOGICAL :: IS_NEW_CHANNEL
      CALL CHANNEL_NEW%SET(J, L, S, TZ)
      IS_NEW_CHANNEL = (CHANNEL /= CHANNEL_NEW)
      CHANNEL = CHANNEL_NEW
      
      CALL PREPARE_LECS_FITTED(ILB, IS_NEW_CHANNEL)
    END BLOCK


    VEM = 0
    IF (LEMP.EQ.0) THEN
      IF (T1Z+T2Z .EQ. 2) VEM(1) = ALPHA*HTC/R
    ELSE
      CALL EMPOT(LEMP,R,VEM)
    ENDIF

    ! IMPLEMENT LEMP != 0

    RC = LECS%RC(S,T)

    IF ( S==0 .AND. T==0 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RETURN
      CASE (1)
        VPW = LECS%CNLO(1) * EFT_RADIAL_1(R, RC) * I2
      CASE (3)
        VPW = LECS%CNLO(1) * EFT_RADIAL_1(R, RC) * I2
        VPW = VPW + LECS%CN3LO(1)  * EFT_RADIAL_4(R, RC) *I2 &
                  + LECS%CN3LO(10) * EFT_RADIAL_7(RC) *L2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW_FITTED:: S=0 AND T=0"
      END SELECT
      VPW = VPW * CR(R, RC)
    ENDIF

    BLOCK 
      USE ANGLES, ONLY: PI
      DOUBLE PRECISION :: RHO, R10
      DOUBLE PRECISION :: C10, D10, E10
    IF ( S==1 .AND. T==0 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RHO = 1.56640D0
        C10 =-0.52342D0 * PI**(1.5D0) * RHO**2
        IF ( L == 0 ) VPW(1,1) = C10 * CR(R, RHO)
      CASE (1)
        RHO = 1.97840D0
        C10 =-0.26039D0 * PI**(1.5D0) * RHO**2
        IF ( L == 0 ) VPW(1,1) = C10 * CR(R, RHO)

        R10 = 1.8D0
        D10 =-0.17194D0  * PI**(1.5D0) * R10**4 / 3.D0
        E10 = 0.36566D-1 * PI**(1.5D0) * R10**4

        VPW(1,1) = VPW(1,1) + (D10 * ( 6.D0*R10**2 - 4.D0*R**2)/R10**4 - E10 * 4.D0 * (R/R10**2)**2 * S12(1,1))* CR(R,R10)
        VPW(2,2) =            (D10 * ( 6.D0*R10**2 - 4.D0*R**2)/R10**4 - E10 * 4.D0 * (R/R10**2)**2 * S12(2,2))* CR(R,R10)
        VPW(1,2) =            (                                        - E10 * 4.D0 * (R/R10**2)**2 * S12(2,1))* CR(R,R10)
        VPW(2,1) =            (                                        - E10 * 4.D0 * (R/R10**2)**2 * S12(1,2))* CR(R,R10)
      CASE (3)
        STOP "NOT IMPLEMENTED POTENTIAL"
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW_FITTED:: S=1 AND T=0"
      END SELECT
    ENDIF
    END BLOCK

    IF ( S==0 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        VPW = LECS%CLO(T) * I2
      CASE (1)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(3) * EFT_RADIAL_1(R, RC) * I2 &
                  + LECS%CIT(0)                * T12 * I2
      CASE (3)
        VPW = LECS%CLO(T) * I2
        VPW = VPW + LECS%CNLO(3)  * EFT_RADIAL_1(R, RC) * I2 &
                  + LECS%CIT(0)                 * T12 * I2
        VPW = VPW + LECS%CN3LO(3) * EFT_RADIAL_4(R, RC)     * I2 &
                  + LECS%CN3LO(10)* EFT_RADIAL_7(RC)     * L2 &
                  + LECS%CIT(1)   * EFT_RADIAL_1(R, RC)     * T12 * I2
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW_FITTED:: S=0 AND T=1"
      END SELECT
      VPW = VPW * CR(R, RC)
    ENDIF

    IF ( S==1 .AND. T==1 ) THEN
      SELECT CASE (ORDER)
      CASE (0)
        RETURN
      CASE (1)
        VPW =   LECS%CNLO(4) * EFT_RADIAL_1(R, RC) *  I2 &
              + LECS%CNLO(6) * EFT_RADIAL_2(R, RC) *  S12 &
              + LECS%CNLO(7) * EFT_RADIAL_3(RC) *  LS &
              + LECS%CIT(0)                *  T12 * I2
      CASE (3)
        LS2 = MATMUL(LS, LS)
        VPW =   LECS%CNLO(4) * EFT_RADIAL_1(R, RC) *  I2 &
              + LECS%CNLO(6) * EFT_RADIAL_2(R, RC) *  S12 &
              + LECS%CNLO(7) * EFT_RADIAL_3(RC) *  LS &
              + LECS%CIT(0)                *  T12 * I2
        VPW =   VPW &
              + LECS%CN3LO(4) * EFT_RADIAL_4(R, RC)    * I2 &
              + LECS%CN3LO(6) * EFT_RADIAL_5(R, RC)    * S12 &
              + LECS%CN3LO(8) * EFT_RADIAL_6(R, RC)    * LS &
              + LECS%CN3LO(9) * EFT_RADIAL_7(RC)    * LS2 &
              + LECS%CN3LO(11)* EFT_RADIAL_7(RC)    * L2 &
              +(  LECS%CIT(2) * EFT_RADIAL_1(R, RC)    * I2  + &
                  LECS%CIT(3) * EFT_RADIAL_2(R, RC)    * S12 + &
                  LECS%CIT(4) * EFT_RADIAL_3(RC)    * LS    &
                                                        ) * T12
      CASE DEFAULT
        STOP "ERROR IN EFT_PLESS_PW_FITTED:: S=1 AND T=1"
      END SELECT
      VPW = VPW * CR(R, RC)
    ENDIF

    VPW = VPW * HTC
    ! IF (R > 1.D-30) VPW = VPW + VEM(1) * I2

    ! Check for NaN in VPW and stop if found
    IF (ANY(VPW /= VPW)) THEN
      PRINT *, "ERROR: VPW contains NaN."
      STOP
    END IF

    IF (GET_CHANNEL_NCH(CHANNEL) == 1 ) THEN
      VPW(1,2) = 0.D0
      VPW(2,1) = 0.D0
      VPW(2,2) = 0.D0
    ENDIF
    RETURN
  END SUBROUTINE EFT_PLESS_PW_FITTED

  SUBROUTINE PREPARE_LECS_FITTED(ILB, NEW_CHANNEL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ILB
    LOGICAL, INTENT(INOUT) :: NEW_CHANNEL
    INTEGER :: REL_ILB
    
    IF (FIRST_CALL) THEN
        CALL SET_ALL_LECS
        FIRST_CALL = .FALSE.
    ENDIF
    IF (NEW_CHANNEL .OR. FIRST_CALL) THEN
      IF (ILB >= 35 .AND. ILB < 40) THEN            ! LO, only 01 and 10 contribute
          CALL PREPARE(5)
          ORDER = 0
        ELSEIF (ILB >= 40 .AND. ILB < 50) THEN        ! NLO
          CALL PREPARE(10)
          ORDER = 1
          SELECT CASE (ILB)
            CASE (40)  ! v11 old fit
              LECS%RC(0,0) = 0.77D0
              LECS%CNLO(1) =-2.213333D0

              LECS%RC(1,1) = 2.895676896183044D0
              LECS%CNLO(4) = 0.722440350712463D0
              LECS%CNLO(6) =-1.145046034106316D0
              ! LECS%CNLO(7) -> standard from old ILB=10

            CASE (41)  ! v11 aJ best fit
              LECS%RC(0,0) = 0.77D0
              LECS%CNLO(1) =-2.213333D0

              LECS%RC(1,1) = 2.8D0
              LECS%CNLO(4) =-0.607942908915985D0
              LECS%CNLO(6) =-1.075304038853378D0
              LECS%CNLO(7) =-0.837789899746530D0

            CASE (42)  ! v11 rJ best fit
              LECS%RC(0,0) = 0.77D0
              LECS%CNLO(1) =-2.213333D0

              LECS%RC(1,1) = 2.43D0
              LECS%CNLO(4) =-0.708536842938800D0
              LECS%CNLO(6) =-0.758503669329703D0
              LECS%CNLO(7) =-0.252811651876925D0

            CASE (43)  ! v11 aJ-rJ best fit
              LECS%RC(0,0) = 0.77D0
              LECS%CNLO(1) =-2.213333D0

              LECS%RC(1,1) = 2.74D0
              LECS%CNLO(4) =-0.786306173096756D0
              LECS%CNLO(6) =-1.089791521668032D0
              LECS%CNLO(7) =-1.136047421482726D0

            CASE (44)  ! v11 old fit
              LECS%RC(0,0) = 0.761781945701797D0
              LECS%CNLO(1) =-2.144536513263479D0

              LECS%RC(1,1) = 2.895676896183044D0
              LECS%CNLO(4) = 0.722440350712463D0
              LECS%CNLO(6) =-1.145046034106316D0
              ! LECS%CNLO(7) -> standard from old ILB=10

            CASE (45)  ! v11 aJ best fit
              LECS%RC(0,0) = 0.761781945701797D0
              LECS%CNLO(1) =-2.144536513263479D0

              LECS%RC(1,1) = 2.8D0
              LECS%CNLO(4) =-0.607942908915985D0
              LECS%CNLO(6) =-1.075304038853378D0
              LECS%CNLO(7) =-0.837789899746530D0

            CASE (46)  ! v11 rJ best fit
              LECS%RC(0,0) = 0.761781945701797D0
              LECS%CNLO(1) =-2.144536513263479D0

              LECS%RC(1,1) = 2.43D0
              LECS%CNLO(4) =-0.708536842938800D0
              LECS%CNLO(6) =-0.758503669329703D0
              LECS%CNLO(7) =-0.252811651876925D0

            CASE (47)  ! v11 aJ-rJ best fit
              LECS%RC(0,0) = 0.761781945701797D0
              LECS%CNLO(1) =-2.144536513263479D0

              LECS%RC(1,1) = 2.74D0
              LECS%CNLO(4) =-0.786306173096756D0
              LECS%CNLO(6) =-1.089791521668032D0
              LECS%CNLO(7) =-1.136047421482726D0

            CASE DEFAULT
              STOP "ERROR IN EFT_PLESS::PREPARE_LECS_FITTED:: Choose the right ILB"
          END SELECT
        ELSEIF (ILB >= 50) THEN                       ! N3LO
          CALL PREPARE(15)
          ORDER = 3
          REL_ILB = ILB - 50 + 1
          SELECT CASE (REL_ILB)

            CASE DEFAULT
              STOP "ERROR IN EFT_PLESS::PREPARE_LECS_FITTED:: Choose the right ILB"
          END SELECT
        ELSE
          STOP "ERROR IN EFT_PLESS::PREPARE_LECS_FITTED:: Choose the right ILB"
        ENDIF
      ENDIF
  END SUBROUTINE PREPARE_LECS_FITTED


END MODULE EFT_PLESS