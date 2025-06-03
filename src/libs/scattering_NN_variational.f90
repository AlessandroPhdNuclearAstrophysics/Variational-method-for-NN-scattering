!> \file scattering_NN_variational.f90
!! \defgroup scattering_nn_variational_mod Scattering NN Variational Module
!! \brief Module for evaluating nucleon-nucleon (NN) scattering wave functions and phase shifts
!!        using the variational and Kohn principles.
!! \ingroup scattering_nn_variational_mod
!! This module provides routines to compute the scattering wave function and phase shifts
!! for a given set of energies and channels (partial waves) in nucleon-nucleon scattering.
!! The calculations are based on the variational principle and the Kohn variational method.
!! The module supports coupled and uncoupled channels, and allows for flexible configuration
!! of the variational basis and grid parameters.
!!
!! \author Alessandro
!! \date 2025
MODULE SCATTERING_NN_VARIATIONAL
  USE REALLOCATE_UTILS
  USE QUANTUM_NUMBERS
  USE EFT_PLESS
  IMPLICIT NONE
  PRIVATE

  !> \ingroup scattering_nn_variational_mod
  !> \brief Structure to store the results of phase shift calculations.
  !! Contains phase shifts and mixing angles in both Blatt-Biedenharn and Stapp conventions.
  TYPE, PUBLIC :: PHASE_SHIFT_RESULT
    DOUBLE PRECISION :: delta1_BB   !< Phase shift 1 (Blatt-Biedenharn convention) [deg]
    DOUBLE PRECISION :: delta2_BB   !< Phase shift 2 (Blatt-Biedenharn convention) [deg]
    DOUBLE PRECISION :: epsilon_BB  !< Mixing angle (Blatt-Biedenharn convention) [deg]
    DOUBLE PRECISION :: delta1_S    !< Phase shift 1 (Stapp convention) [deg]
    DOUBLE PRECISION :: delta2_S    !< Phase shift 2 (Stapp convention) [deg]
    DOUBLE PRECISION :: epsilon_S   !< Mixing angle (Stapp convention) [deg]
  END TYPE PHASE_SHIFT_RESULT

  !> \ingroup scattering_nn_variational_mod
  !> \brief Structure containing all variational parameters for the scattering calculation.
  !! This includes quantum numbers, grid parameters, and basis size.
  TYPE, PUBLIC :: VARIATIONAL_PARAMETERS
    INTEGER :: J      !< Total angular momentum
    INTEGER :: L      !< Orbital angular momentum
    INTEGER :: S      !< Spin
    INTEGER :: T      !< Isospin
    INTEGER :: TZ     !< Isospin projection
    INTEGER :: T1Z    !< Isospin projection of particle 1
    INTEGER :: T2Z    !< Isospin projection of particle 2
    INTEGER :: IPOT = 0  !< Potential type index
    INTEGER :: ILB = 1    !< Interaction type (default 1)
    INTEGER :: LEMP = 0   !< Electromagnetic potential flag (default 0)
    DOUBLE PRECISION :: E     !< Scattering energy [MeV]
    DOUBLE PRECISION :: K     !< Scattering momentum [fm^-1]
    DOUBLE PRECISION :: HR1 = 0.01D0   !< Step size for core grid [fm]
    DOUBLE PRECISION :: H   = 0.02D0   !< Step size for asymptotic grid [fm]
    DOUBLE PRECISION :: RANGE = 40.0D0 !< Maximum radial range [fm]
    DOUBLE PRECISION :: GAMMA = 4.D0   !< Scaling parameter for Laguerre basis
    DOUBLE PRECISION :: EPS = 0.25D0   !< Exponential grid parameter
    DOUBLE PRECISION :: AF = 1.02D0    !< Exponential grid parameter
    INTEGER :: NX_AA         !< Number of points in asymptotic-asymptotic grid
    INTEGER :: NX_AC         !< Number of points in asymptotic-core grid
    INTEGER :: NX_CC = 100   !< Number of points in core-core grid (default 100)
    INTEGER :: NNL = 32      !< Number of Laguerre basis functions (default 32)
  END TYPE VARIATIONAL_PARAMETERS

  INTEGER, PARAMETER :: NNE = 80
  INTEGER, PARAMETER :: NCH_MAX = 2
  INTEGER :: NCH
  INTEGER :: NNN_MAX
  LOGICAL :: USE_DYNAMIC = .FALSE.

  DOUBLE PRECISION, PARAMETER :: HC = 197.327053D0
  DOUBLE PRECISION, PARAMETER :: MP = 938.272029D0
  DOUBLE PRECISION, PARAMETER :: MN = 939.565630D0
  DOUBLE PRECISION, PARAMETER :: MR = MP * MN / (MP + MN)
  DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0, ZERO = 0.0D0
  DOUBLE PRECISION, PARAMETER :: PI = 4*DATAN(ONE)
  DOUBLE COMPLEX,   PARAMETER :: IM = (ZERO, ONE)

  DOUBLE PRECISION :: HTM, M
  LOGICAL :: HTM_SET = .FALSE.
  INTEGER :: LC(NCH_MAX)
  LOGICAL :: PRINT_I

  INTEGER :: LMAX=-1
  LOGICAL :: TZ_SET = .FALSE., LMAX_SET = .FALSE.

  ! Channels analyzed
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS_(:)
  TYPE(SCATTERING_CHANNEL) :: CHANNEL
  INTEGER :: CH_INDEX, NCHANNELS
  LOGICAL :: CHANNELS_SET = .FALSE.

  ! Commonly used r grids
  DOUBLE PRECISION, ALLOCATABLE :: XX_CC(:), YY_CC(:), XX_AC(:), XX_AA(:)
  DOUBLE PRECISION, ALLOCATABLE :: A_CC(:), A_AC(:), A_AA(:), B_AA(:), AJ_AA(:), AJ_AC(:), YYL_AC(:)
  LOGICAL :: GRID_SET = .FALSE.

  ! Potential for all the channels (first index) V_XX(CHANNEL, IR, ALPHA, ALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: V_CC(:,:,:,:), V_AC(:,:,:,:), V_AA(:,:,:,:)
  LOGICAL :: POTENTIAL_SET = .FALSE.
  LOGICAL :: IPOT_SET = .FALSE.

  ! Commonly used variables depending on the channels and/or the energy
  DOUBLE PRECISION, ALLOCATABLE :: ENERGIES_(:), KK(:), K2(:)
  INTEGER :: NE = -1
  LOGICAL :: ENERGIES_SET = .FALSE.

  ! Bessel functions F_L (E, L, IX)
  DOUBLE PRECISION, ALLOCATABLE :: FBES_AA (:,:,:), FBES_AC (:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GBES_AA (:,:,:), GBES_AC (:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GBES0_AA(:,:,:), GBES1_AA(:,:,:), GBES2_AA(:,:,:), HNOR_AA(:,:,:)
  LOGICAL :: BESSELS_SET = .FALSE.

  ! Laguerre functions modified and their derivatives
  DOUBLE PRECISION, ALLOCATABLE :: V0_CC(:,:), V1_CC(:,:), V2_CC(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: V0_AC(:,:), V1_AC(:,:), V2_AC(:,:)
  LOGICAL :: LAGUERRE_SET = .FALSE.


  ! Matrix elements for all energies and channels, evaluated once in case USE_DYNAMIC is .FALSE.
  ! H_MINUS_E_CC    = < n alpha | H - E | n' alpha' >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! H_MINUS_E_AC_R  = < F_ALPHA | H - E | n' alpha' >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! H_MINUS_E_AC_I  = < G_ALPHA | H - E | n' alpha' >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! H_MINUS_E_AA_RR = < F_ALPHA | H - E | F_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! H_MINUS_E_AA_RI = < F_ALPHA | H - E | G_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! H_MINUS_E_AA_IR = < G_ALPHA | H - E | F_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! H_MINUS_E_AA_II = < G_ALPHA | H - E | G_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  DOUBLE PRECISION, ALLOCATABLE :: H_MINUS_E_CC   (:,:,:,:)  ! H_MINUS_E_CC   (CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: H_MINUS_E_AC_R (:,:,:,:)  ! H_MINUS_E_AC_R (CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: H_MINUS_E_AC_I (:,:,:,:)  ! H_MINUS_E_AC_I (CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: H_MINUS_E_AA_RR(:,:,:,:)  ! H_MINUS_E_AA_RR(CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: H_MINUS_E_AA_RI(:,:,:,:)  ! H_MINUS_E_AA_RI(CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: H_MINUS_E_AA_IR(:,:,:,:)  ! H_MINUS_E_AA_IR(CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: H_MINUS_E_AA_II(:,:,:,:)  ! H_MINUS_E_AA_II(CH_INDEX, E, NALPHA, NALPHA')

  ! K_MINUS_E_CC    = < n alpha | K - E | n' alpha' >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! K_MINUS_E_AC_R  = < F_ALPHA | K - E | n' alpha' >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! K_MINUS_E_AC_I  = < G_ALPHA | K - E | n' alpha' >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! K_MINUS_E_AA_RR = < F_ALPHA | K - E | F_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! K_MINUS_E_AA_RI = < F_ALPHA | K - E | G_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! K_MINUS_E_AA_IR = < G_ALPHA | K - E | F_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  ! K_MINUS_E_AA_II = < G_ALPHA | K - E | G_ALPHA'  >, depends on (n alpha), (n' alpha'), E and the channel (reverse order)
  TYPE(EFT_RADIAL_FUNCTIONS) :: EFT_RADIAL_CC, EFT_RADIAL_AC, EFT_RADIAL_AA
  TYPE(LECS_EFT_PLESS) :: LECS
  LOGICAL :: LECS_SET = .FALSE.
  LOGICAL :: NEW_LECS = .TRUE.
  DOUBLE PRECISION, ALLOCATABLE :: K_MINUS_E_CC   (:,:,:,:)  ! K_MINUS_E_CC   (CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: K_MINUS_E_AC_R (:,:,:,:)  ! K_MINUS_E_AC_R (CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: K_MINUS_E_AC_I (:,:,:,:)  ! K_MINUS_E_AC_I (CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: K_MINUS_E_AA_RR(:,:,:,:)  ! K_MINUS_E_AA_RR(CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: K_MINUS_E_AA_RI(:,:,:,:)  ! K_MINUS_E_AA_RI(CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: K_MINUS_E_AA_IR(:,:,:,:)  ! K_MINUS_E_AA_IR(CH_INDEX, E, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: K_MINUS_E_AA_II(:,:,:,:)  ! K_MINUS_E_AA_II(CH_INDEX, E, NALPHA, NALPHA')

  ! VM_CC    = < n alpha | V | n' alpha' >, depends on (n alpha), (n' alpha') and the channel (reverse order)
  ! VM_AC_R  = < F_ALPHA | V | n' alpha' >, depends on (n alpha), (n' alpha') and the channel (reverse order)
  ! VM_AC_I  = < G_ALPHA | V | n' alpha' >, depends on (n alpha), (n' alpha') and the channel (reverse order)
  ! VM_AA_RR = < F_ALPHA | V | F_ALPHA'  >, depends on (n alpha), (n' alpha') and the channel (reverse order)
  ! VM_AA_RI = < F_ALPHA | V | G_ALPHA'  >, depends on (n alpha), (n' alpha') and the channel (reverse order)
  ! VM_AA_IR = < G_ALPHA | V | F_ALPHA'  >, depends on (n alpha), (n' alpha') and the channel (reverse order)
  ! VM_AA_II = < G_ALPHA | V | G_ALPHA'  >, depends on (n alpha), (n' alpha') and the channel (reverse order)
  DOUBLE PRECISION, ALLOCATABLE :: VM_CC   (:,:,:,:)  ! VM_CC   (CH_INDEX, NE, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: VM_AC_R (:,:,:,:)  ! VM_AC_R (CH_INDEX, NE, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: VM_AC_I (:,:,:,:)  ! VM_AC_I (CH_INDEX, NE, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: VM_AA_RR(:,:,:,:)  ! VM_AA_RR(CH_INDEX, NE, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: VM_AA_RI(:,:,:,:)  ! VM_AA_RI(CH_INDEX, NE, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: VM_AA_IR(:,:,:,:)  ! VM_AA_IR(CH_INDEX, NE, NALPHA, NALPHA')
  DOUBLE PRECISION, ALLOCATABLE :: VM_AA_II(:,:,:,:)  ! VM_AA_II(CH_INDEX, NE, NALPHA, NALPHA')

  ! Matrix elements for the variational method, these are
  ! < n alpha | V | n' alpha' > = SUM_i C_i (int_0^Range  dr r^2 f_n(r) f_n'(r) < n alpha | V_i(r) | n' alpha' >)
  DOUBLE PRECISION, ALLOCATABLE :: FMAT_CC   (:,:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: FMAT_AC_R (:,:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: FMAT_AC_I (:,:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: FMAT_AA_RR(:,:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: FMAT_AA_RI(:,:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: FMAT_AA_IR(:,:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: FMAT_AA_II(:,:,:,:,:)

  TYPE(VARIATIONAL_PARAMETERS), PRIVATE :: VAR_P

  PUBLIC :: SET_NNL
  PUBLIC :: SET_AF
  PUBLIC :: SET_IPOT
  PUBLIC :: SET_CHANNELS
  PUBLIC :: SET_ENERGIES
  PUBLIC :: NN_SCATTERING_VARIATIONAL
  PUBLIC :: SET_VARIATIONAL_PARAMETERS
  PUBLIC :: GET_HTM
  PUBLIC :: SET_DYNAMIC
  PUBLIC :: SET_NEW_LECS
  PUBLIC :: RESET_SCATTERING_NN_VARIATIONAL

  PRIVATE:: PRINT_DIVIDER
  PRIVATE:: SET_M_T1Z_T2Z_HTM
  PRIVATE:: PREPARE_GRID
  PRIVATE:: PREPARE_LAGUERRE
  PRIVATE:: PREPARE_POTENTIAL
  PRIVATE:: SET_VARIATIONAL_PARAMETERS_
  PRIVATE:: PREPARE_ASYMPTOTIC_FUNCTIONS
  PRIVATE:: PREPARE_ASYMPTOTIC_CORE_MATRIX_ELEMENTS
  PRIVATE:: PREPARE_CORE_CORE_MATRIX_ELEMENTS
  PRIVATE:: PREPARE_ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS

CONTAINS
!> \ingroup scattering_nn_variational_mod
!> \brief Set all variational parameters at once.
  !! \param[in] J    Total angular momentum
  !! \param[in] L    Orbital angular momentum
  !! \param[in] S    Spin
  !! \param[in] T    Isospin
  !! \param[in] TZ   Isospin projection
  !! \param[in] IPOT Potential type index
  !! \param[in] ILB  Interaction type (optional)
  !! \param[in] LEMP Electromagnetic potential flag (optional)
  !! \param[in] HR1  Step size for core grid [fm] (optional)
  !! \param[in] H    Step size for asymptotic grid [fm] (optional)
  !! \param[in] RANGE Maximum radial range [fm] (optional)
  !! \param[in] GAMMA Scaling parameter for Laguerre basis (optional)
  !! \param[in] EPS  Exponential grid parameter (optional)
  !! \param[in] AF   Exponential grid parameter (optional)
  !! \param[in] NX_AA Number of points in asymptotic-asymptotic grid (optional)
  !! \param[in] NX_AC Number of points in asymptotic-core grid (optional)
  !! \param[in] NX_CC Number of points in core-core grid (optional)
  !! \param[in] NNL   Number of Laguerre basis functions (optional)
  SUBROUTINE SET_VARIATIONAL_PARAMETERS(J, L, S, TZ, IPOT, T, ILB, LEMP, HR1, H, &
                                            RANGE, GAMMA, EPS, AF, NX_AA, NX_AC, NX_CC, NNL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT
    INTEGER, OPTIONAL, INTENT(IN) :: T, NX_AA, NX_CC, NX_AC, NNL, ILB, LEMP
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: HR1, H, RANGE, GAMMA, EPS, AF

    VAR_P%J     = J
    VAR_P%L     = L
    VAR_P%S     = S
    IF (PRESENT(T)) THEN
      VAR_P%T = T
    ELSE
      VAR_P%T = MOD(MOD((L+S),2)+1,2)
    ENDIF
    VAR_P%TZ    = TZ
    IF (ABS(TZ)>VAR_P%T) THEN
      PRINT *, "Error: |Tz| > T!"
      STOP
    ENDIF
    VAR_P%IPOT  = IPOT
    IF (PRESENT(ILB))   VAR_P%ILB   = ILB
    IF (PRESENT(LEMP))  VAR_P%LEMP  = LEMP
    IF (PRESENT(HR1))   VAR_P%HR1   = HR1
    IF (PRESENT(H))     VAR_P%H     = H
    IF (PRESENT(RANGE)) VAR_P%RANGE = RANGE
    IF (PRESENT(GAMMA)) VAR_P%GAMMA = GAMMA
    IF (PRESENT(EPS))   VAR_P%EPS   = EPS
    IF (PRESENT(AF))    VAR_P%AF    = AF
    IF (PRESENT(NX_AA)) VAR_P%NX_AA = NX_AA
    IF (PRESENT(NX_AC)) VAR_P%NX_AC = NX_AC
    IF (PRESENT(NX_CC)) VAR_P%NX_CC = NX_CC
    IF (PRESENT(NNL))   VAR_P%NNL   = NNL

    SELECT CASE (TZ)
    CASE (1)
      VAR_P%T1Z = 1
      VAR_P%T2Z = 1
    CASE (0)
      VAR_P%T1Z =-1
      VAR_P%T2Z = 1
    CASE (-1)
      VAR_P%T1Z =-1
      VAR_P%T2Z =-1
    END SELECT
    CALL SET_M_T1Z_T2Z_HTM(TZ)

  END SUBROUTINE SET_VARIATIONAL_PARAMETERS

  !> \ingroup scattering_nn_variational_mod
  !> \brief Set the energies at which to compute phase shifts and wave functions.
  !! \param[in] ENERGIES Array of energies [MeV]
  SUBROUTINE SET_ENERGIES(ENERGIES)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: ENERGIES(:)

    NE = SIZE(ENERGIES)
    CALL REALLOCATE_1D_1(ENERGIES_, NE)
    ENERGIES_ = ENERGIES
    ENERGIES_SET = .TRUE.
    IF (CHANNELS_SET .AND. .NOT.BESSELS_SET .AND. HTM_SET) CALL PREPARE_ASYMPTOTIC_FUNCTIONS
  END SUBROUTINE SET_ENERGIES

  SUBROUTINE SET_IPOT(IPOT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IPOT
    VAR_P%IPOT = IPOT
    IPOT_SET = .TRUE.
  END SUBROUTINE SET_IPOT

  !> \ingroup scattering_nn_variational_mod
  !> \brief Set the scattering channels to be analyzed.
  !! \param[in] CHANNELS Array of SCATTERING_CHANNEL structures
  SUBROUTINE SET_CHANNELS(CHANNELS)
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNELS(:)
    INTEGER :: I1, I2

    NCHANNELS = SIZE(CHANNELS)
    ALLOCATE(CHANNELS_(NCHANNELS))
    CHANNELS_ = CHANNELS
    CHANNELS_SET = .TRUE.
    if (.NOT.POTENTIAL_SET .AND. IPOT_SET) CALL PREPARE_POTENTIAL(CHANNELS_)
    DO I1 = 1, SIZE(CHANNELS_)
      DO I2 = 1, GET_CHANNEL_NCH(CHANNELS_(I1))
        IF (GET_CHANNEL_L(CHANNELS_(I1), I2) > LMAX) LMAX = GET_CHANNEL_L(CHANNELS_(I1), I2)
      ENDDO
    ENDDO
    VAR_P%TZ = GET_CHANNEL_TZ(CHANNELS_(1))
    TZ_SET = .TRUE.
    LMAX_SET = .TRUE.
    IF (ENERGIES_SET .AND. .NOT.BESSELS_SET) CALL PREPARE_ASYMPTOTIC_FUNCTIONS
  END SUBROUTINE SET_CHANNELS

  !> \ingroup scattering_nn_variational_mod
  !> \brief Set the number of Laguerre basis functions for the variational calculation.
  !! \param[in] NNL Number of Laguerre basis functions
  SUBROUTINE SET_NNL(NNL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NNL
    VAR_P%NNL = NNL
  END SUBROUTINE SET_NNL

  !> \ingroup scattering_nn_variational_mod
  !> \brief Set the exponential grid parameter AF.
  !! \param[in] AF Exponential grid parameter
  SUBROUTINE SET_AF(AF)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: AF
    VAR_P%AF = AF
  END SUBROUTINE SET_AF

  !> \ingroup scattering_nn_variational_mod
  !> \brief Set all variational parameters for a single calculation (internal use).
  !! \param[in] E, J, L, S, TZ, IPOT, ILB, LEMP
  SUBROUTINE SET_VARIATIONAL_PARAMETERS_(E, J, L, S, TZ, IPOT, ILB, LEMP)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: E
    INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP

    WRITE(*,*) "Setting variational parameters..."

    VAR_P%J = J
    VAR_P%L = L
    VAR_P%S = S
    VAR_P%TZ = TZ
    VAR_P%IPOT = IPOT
    VAR_P%ILB = ILB
    VAR_P%LEMP = LEMP
    VAR_P%E = E

    VAR_P%NX_AC = INT(VAR_P%RANGE/VAR_P%HR1) + 10

    SELECT CASE (TZ)
    CASE (1)
      VAR_P%T1Z = 1
      VAR_P%T2Z = 1
    CASE (0)
      VAR_P%T1Z =-1
      VAR_P%T2Z = 1
    CASE (-1)
      VAR_P%T1Z =-1
      VAR_P%T2Z =-1
    END SELECT

    CALL SET_M_T1Z_T2Z_HTM(TZ)

  ! Ensure T is set such that T + L + S is odd
    IF (MOD(L + S, 2) == 0) THEN
      VAR_P%T = 1
    ELSE
      VAR_P%T = 0
    END IF

    VAR_P%K = DSQRT(2*E*MR) / HC

    IF (PRINT_I) THEN
      PRINT 5
      PRINT 15, "L", "S", "T", "TZ", "J"
      PRINT 5
    ENDIF
    LC(1) = L
    NCH = 1
    IF (PRINT_I) PRINT 20, LC(1), S, VAR_P%T, VAR_P%TZ, J
    IF (J-L.EQ.1) THEN
      LC(2) = L + 2
      IF (PRINT_I) PRINT 20, LC(2), S, VAR_P%T, VAR_P%TZ, J
      NCH = 2
    ENDIF
    IF (PRINT_I) PRINT 5


  5 FORMAT(30("-"))
 15 FORMAT(" ", A5, A5, A5, A5, A5, A5)
 20 FORMAT(" ", I5, I5, I5, I5, I5, I5)
  END SUBROUTINE SET_VARIATIONAL_PARAMETERS_


  !> \ingroup scattering_nn_variational_mod
  !> \brief Main routine to perform NN scattering calculation using the variational method.
  !! Computes phase shifts and mixing angles for given quantum numbers and energy.
  !!
  !! \param[in]  E      Scattering energy [MeV]
  !! \param[in]  J      Total angular momentum
  !! \param[in]  L      Orbital angular momentum
  !! \param[in]  S      Spin
  !! \param[in]  TZ     Isospin projection
  !! \param[in]  IPOT   Potential type index
  !! \param[in]  ILB    Interaction type
  !! \param[in]  LEMP   Electromagnetic potential flag
  !! \param[out] PHASE_SHIFT  Structure with phase shifts and mixing angles
  !! \param[in]  PRINT_COEFFICIENTS (optional) Print wave function coefficients
  !! \param[in]  PRINT_INFORMATIONS (optional) Print detailed calculation info
  SUBROUTINE NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP, PHASE_SHIFT, &
   PRINT_COEFFICIENTS, PRINT_INFORMATIONS)
  IMPLICIT NONE
! INPUT PARAMETERS
  DOUBLE PRECISION, INTENT(IN) :: E
  INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP
  LOGICAL, INTENT(IN), OPTIONAL :: PRINT_COEFFICIENTS, PRINT_INFORMATIONS
  TYPE(PHASE_SHIFT_RESULT), INTENT(OUT) :: PHASE_SHIFT

  LOGICAL :: PRINT_C, FIRST_CALL = .TRUE., PREPARE = .TRUE.
! VARIABLES AND PARAMETERS FOR DGESV
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: CAR(:,:), CAI(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: CARR(:), CAII(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: XRCOEFF(:,:), XICOEFF(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: IPIV(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: C(:,:), CC(:,:), CCC(:,:)
  INTEGER :: INFO, IAK, IAB, IE
  INTEGER, SAVE :: NNN

  ! MATRICES FOR THE VARIATIONAL METHOD
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: ARI, AIR, ARR, AII
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: BD1(:,:), BD2(:,:), BD3(:,:), BD4(:,:)
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: AM, AN, AMM, RMAT

  ! COEFFICIENT RMAT2
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: RMAT2

  ! PHASE-SHIFTS AND MIXING ANGLES
  DOUBLE PRECISION :: AMIXR, AMIXG, DELTA1, DELTA2, DELTA1G, DELTA2G
  DOUBLE PRECISION :: AMIXGS, DELTA1S, DELTA2S

! S-MATRIX
  DOUBLE COMPLEX :: SMAT(NCH_MAX, NCH_MAX)

  FIRST_CALL = FIRST_CALL .OR. IS_FIRST_CALL(J, L, S, TZ, IPOT, ILB, LEMP)
  IF (.NOT.POTENTIAL_SET) THEN
    IF (.NOT.CHANNELS_SET) THEN
      PRINT *, "Error: Potential not set and channels not set!"
      STOP
    ENDIF
    IF (ENERGIES_SET) THEN
      CALL SET_VARIATIONAL_PARAMETERS(J, L, S, TZ, IPOT, ILB=ILB, LEMP=LEMP)
      CALL PREPARE_POTENTIAL(CHANNELS_)
      POTENTIAL_SET = .TRUE.
    ELSE
      PRINT *, "Error: Potential not set and energies not set!"
      STOP
    ENDIF
  ENDIF


  IF (PRESENT(PRINT_COEFFICIENTS)) THEN
    PRINT_C = PRINT_COEFFICIENTS
  ELSE
    PRINT_C = .FALSE.
  ENDIF

  IF (PRESENT(PRINT_INFORMATIONS)) THEN
    PRINT_I = PRINT_INFORMATIONS
  ELSE
    PRINT_I = .FALSE.
  ENDIF


  ! INITIALIZE THE VARIATIONAL PARAMETERS
  IF (FIRST_CALL) THEN
    CALL SET_VARIATIONAL_PARAMETERS_(E, J, L, S, TZ, IPOT, ILB, LEMP)
    IF (.NOT.GRID_SET) CALL PREPARE_GRID
    NNN     = VAR_P%NNL * NCH
    NNN_MAX = VAR_P%NNL * NCH_MAX
    CALL REALLOCATE_2D_3(C, CC, CCC, NNN, NNN)
    CALL REALLOCATE_2D_2(CAR, CAI, NNN, NCH)
    CALL REALLOCATE_1D_1(CARR, NNN)
    CALL REALLOCATE_1D_1(CAII, NNN)
    CALL REALLOCATE_2D_2(XRCOEFF, XICOEFF, NCH, NNN)
    CALL REALLOCATE_1D_1(IPIV, NNN)
    CALL REALLOCATE_2D_4(BD1, BD2, BD3, BD4, NCH, NCH)


    CALL SET_CHANNEL(CHANNEL, J, L, S, TZ)
    CH_INDEX = FIND_CHANNEL_INDEX()
    WRITE(*,*) "CH_INDEX", CH_INDEX
  ELSE
    VAR_P%E = E
    VAR_P%K = DSQRT(2*E*MR) / HC
  ENDIF

  IF (.NOT.GRID_SET .OR. .NOT.BESSELS_SET) THEN
    STOP "Grid not ready or Bessels not ready"
  ENDIF

  IF (PRINT_I) CALL PRINT_INFO()
  IF (PREPARE) THEN
    CALL PREPARE_CORE_CORE_MATRIX_ELEMENTS
    IF (PRINT_I) CALL PRINT_DIVIDER
    CALL PREPARE_ASYMPTOTIC_CORE_MATRIX_ELEMENTS
  ENDIF

  IF (PRINT_I) CALL PRINT_DIVIDER

  IE = FIND_ENERGY_INDEX(E)

  IF (USE_DYNAMIC .AND. NEW_LECS) THEN
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_CC,    LECS, VM_CC    )
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AC_R,  LECS, VM_AC_R  )
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AC_I,  LECS, VM_AC_I  )
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_RR, LECS, VM_AA_RR )
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_RI, LECS, VM_AA_RI )
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_IR, LECS, VM_AA_IR )
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_II, LECS, VM_AA_II )
    
    
    H_MINUS_E_CC(CH_INDEX, IE, 1:NNN, 1:NNN)   = K_MINUS_E_CC(CH_INDEX, IE, 1:NNN, 1:NNN) &
                                                + VM_CC(CH_INDEX, IE, 1:NNN, 1:NNN)
    H_MINUS_E_AC_R(CH_INDEX, IE, 1:NNN, 1:NCH) = K_MINUS_E_AC_R(CH_INDEX, IE, 1:NNN, 1:NCH) &
                                                + VM_AC_R(CH_INDEX, IE, 1:NNN, 1:NCH)
    H_MINUS_E_AC_I(CH_INDEX, IE, 1:NNN, 1:NCH) = K_MINUS_E_AC_I(CH_INDEX, IE, 1:NNN, 1:NCH) &
                                                + VM_AC_I(CH_INDEX, IE, 1:NNN, 1:NCH)
    H_MINUS_E_AA_RR(CH_INDEX, IE, :, :)       = K_MINUS_E_AA_RR(CH_INDEX, IE, :, :) &
                                              + VM_AA_RR(CH_INDEX, IE, :, :)
    H_MINUS_E_AA_RI(CH_INDEX, IE, :, :)       = K_MINUS_E_AA_RI(CH_INDEX, IE, :, :) &
                                              + VM_AA_RI(CH_INDEX, IE, :, :)
    H_MINUS_E_AA_IR(CH_INDEX, IE, :, :)       = K_MINUS_E_AA_IR(CH_INDEX, IE, :, :) &
                                              + VM_AA_IR(CH_INDEX, IE, :, :)
    H_MINUS_E_AA_II(CH_INDEX, IE, :, :)       = K_MINUS_E_AA_II(CH_INDEX, IE, :, :) &
                                              + VM_AA_II(CH_INDEX, IE, :, :)

    NEW_LECS = .FALSE.
  ENDIF

  C   = H_MINUS_E_CC  (CH_INDEX, IE, 1:NNN, 1:NNN)  ! H - E
  CAR = H_MINUS_E_AC_R(CH_INDEX, IE, 1:NNN, 1:NCH)  ! H - E
  CAI = H_MINUS_E_AC_I(CH_INDEX, IE, 1:NNN, 1:NCH)  ! H - E
  IF (PRINT_I) CALL PRINT_DIVIDER

  DO IAK=1, NCH
! Preparing the matrix elements for the diagonalization
    CARR =-CAR(:,IAK)
    CAII =-CAI(:,IAK)
    CC  = C
    CCC = C

  ! Solving the eigenvalue problem using LAPACK dgesv
  ! Evaluating for the "c_{n, alpha}" coefficients
    CALL DGESV(NNN, 1, CC , NNN, IPIV, CARR, NNN, INFO)
    CALL HANDLE_INFO_ERROR()  ! Handle the error after the first DGESV call
    IF (PRINT_I) WRITE(*,*) "INFO: ", INFO
    CALL DGESV(NNN, 1, CCC, NNN, IPIV, CAII, NNN, INFO)
    CALL HANDLE_INFO_ERROR()  ! Handle the error after the second DGESV call
    IF (PRINT_I) WRITE(*,*) "INFO: ", INFO

    XRCOEFF(IAK,:) = CARR
    XICOEFF(IAK,:) = CAII
  ENDDO

! Calculating R coefficients
  IF (PRINT_I) WRITE(*,*) NCH, NNN

  ! This performs matrix multiplication using DGEMM -> BD# = 1.d0*(X#COEFF * CA#) + 0.d0*BD#
  CALL DGEMM('N', 'N', NCH, NCH, NNN, 1.0D0, XICOEFF, SIZE(XICOEFF,1), CAI, SIZE(CAI,1), 0.0D0, BD1, SIZE(BD1,1))
  CALL DGEMM('N', 'N', NCH, NCH, NNN, 1.0D0, XRCOEFF, SIZE(XRCOEFF,1), CAI, SIZE(CAI,1), 0.0D0, BD2, SIZE(BD2,1))
  CALL DGEMM('N', 'N', NCH, NCH, NNN, 1.0D0, XICOEFF, SIZE(XICOEFF,1), CAR, SIZE(CAR,1), 0.0D0, BD3, SIZE(BD3,1))
  CALL DGEMM('N', 'N', NCH, NCH, NNN, 1.0D0, XRCOEFF, SIZE(XRCOEFF,1), CAR, SIZE(CAR,1), 0.0D0, BD4, SIZE(BD4,1))

  IF (PREPARE) THEN
    IF (PRINT_I) CALL PRINT_DIVIDER
    CALL PREPARE_ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS
    IF (PRINT_I) CALL PRINT_DIVIDER
    PREPARE = .FALSE.
  ENDIF

  ARR = H_MINUS_E_AA_RR(CH_INDEX, IE, :, :)
  ARI = H_MINUS_E_AA_RI(CH_INDEX, IE, :, :)
  AIR = H_MINUS_E_AA_IR(CH_INDEX, IE, :, :)
  AII = H_MINUS_E_AA_II(CH_INDEX, IE, :, :)

  DO IAB = 1, NCH
  DO IAK = 1, NCH
    AM(IAB,IAK)=BD1(IAK,IAB)+BD1(IAB,IAK)+AII(IAB,IAK)+AII(IAK,IAB)
    AN(IAB,IAK)=BD2(IAK,IAB)+BD3(IAB,IAK)+2.D0*AIR(IAB,IAK)
  ENDDO
  ENDDO



  AMM = AM
  RMAT =-AN
  IF (PRINT_I) THEN
    CALL PRINT_DIVIDER
    WRITE(*,*) "AMM = ", AMM
    WRITE(*,*) "RMAT = ", RMAT
    CALL PRINT_DIVIDER
  ENDIF


! Evaluating the "R_{alpha, beta}" matrix elements
  CALL DGESV(NCH, NCH, AMM, NCH_MAX, IPIV, RMAT, NCH_MAX, INFO)
  CALL HANDLE_INFO_ERROR()  ! Handle the error after the third DGESV call
  IF (PRINT_I)  WRITE(*,*) "INFO: ", INFO

  IF (PRINT_I) THEN
    DO IAB = 1, NCH
    DO IAK = 1, NCH
      WRITE(*,*)"COEFF R",RMAT(IAB,IAK)
    ENDDO
    ENDDO
  ENDIF

! Evaluating the "R_{alpha, beta}" matrix elements to the second order
  CALL R_SECOND_ORDER()
  IF (PRINT_I) THEN
    WRITE(*,*)
    DO IAB = 1, NCH
    DO IAK = 1, NCH
      WRITE(*,*)"COEFF RMAT2 NORMALIZZATO", -RMAT2(IAB,IAK)
    ENDDO
    ENDDO
  ENDIF

! Writing the coefficients to a file in order torecreate the wave function
  IF (PRESENT(PRINT_COEFFICIENTS)) THEN
    IF (PRINT_COEFFICIENTS) THEN
      CALL WRITE_COEFFICIENTS_TO_RECREATE_THE_WAVE_FUNCTION()
    ENDIF
  ENDIF

! Calculating the phase shifts and mixing angles in the Blatt-Biedenharn convention
  CALL CALCULATE_PHASE_SHIFTS_BLATT

  IF (PRINT_I) THEN
    WRITE(*,*)
    WRITE(*,*)"BLATT-BIEDENHARN"
    WRITE(*,*)"MIXING ANGLE=",AMIXG
    WRITE(*,*)"SFASAMENTO1=",DELTA1G
    WRITE(*,*)"SFASAMENTO2=",DELTA2G
  ENDIF

  PHASE_SHIFT%delta1_BB = DELTA1G
  PHASE_SHIFT%delta2_BB = DELTA2G
  PHASE_SHIFT%epsilon_BB = AMIXG


! Calculating the S-matrix
  CALL CALCULATE_S_MATRIX

  IF (PRINT_I) THEN
    WRITE(*,*)
    WRITE(*,*)"S-MATRIX"
    WRITE(*,*) "S(1,1)=" , SMAT(1,1)
    WRITE(*,*) "S(1,2)=" , SMAT(1,2)
    WRITE(*,*) "S(2,2)=" , SMAT(2,2)
  ENDIF


! Calculating the phase shifts and mixing angles in the Stapp convention
  CALL CALCULATE_PHASE_SHIFTS_STAPP

  IF (PRINT_I) THEN
    WRITE(*,*)
    WRITE(*,*)"STAPP"
    WRITE(*,*) "MIXING ANGLE=",AMIXGS
    WRITE(*,*) "SFASAMENTO1=",DELTA1S
    WRITE(*,*) "SFASAMENTO2=",DELTA2S
  ENDIF

  PHASE_SHIFT%delta1_S = DELTA1S
  PHASE_SHIFT%delta2_S = DELTA2S
  PHASE_SHIFT%epsilon_S = AMIXGS

  IF (PRINT_I) WRITE(*,*) DELTA1S, DELTA2S, AMIXGS

  RETURN

  CONTAINS
    !> \brief Handle LAPACK DGESV INFO error.
    SUBROUTINE HANDLE_INFO_ERROR()
      IMPLICIT NONE
      IF (INFO.NE.0) THEN
        PRINT *, "Error in LAPACK dgesv: INFO = ", INFO
        STOP
      ENDIF
    END SUBROUTINE HANDLE_INFO_ERROR

    !> \brief Print information about the current calculation.
    SUBROUTINE PRINT_INFO()
      IMPLICIT NONE
      PRINT *, "E:    ",                  VAR_P%E, " MeV"
      PRINT *, "HTM:  ",                  HTM, " MeV fm^2"
      PRINT *, "k:    ",                  VAR_P%K, " fm^-1"
      PRINT 10, "J:     ",                VAR_P%J
      PRINT 10, "NCH:   ",                NCH
      PRINT 10, "L0:    ",                LC(1)
      IF (NCH.EQ.2) PRINT 10, "L1:    ",  LC(2)
      PRINT 10, "S:     ",                VAR_P%S
      PRINT 10, "T:     ",                VAR_P%T
      PRINT 10, "TZ:    ",                VAR_P%TZ

      10 FORMAT(" ",A, I2)
    END SUBROUTINE PRINT_INFO

    !> \brief Compute the second order R-matrix.
    SUBROUTINE R_SECOND_ORDER()
      IMPLICIT NONE
      INTEGER :: I, IK, IB
      DOUBLE PRECISION :: SOMMA, SOMMA1, SOMMA2, SOMMA3, SOMMA4, SOMMA5, SOMMA6
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: XRCOEFV(:,:), XICOEFV(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: BD5(:,:), BD6(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: RCI(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: RCIV(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: CD0(:,:), CD1(:,:), CD2(:,:), CD3(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: CD4(:,:), CD5(:,:), CD6(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: CD7(:,:), CD8(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: RD0(:,:), RD2(:,:), RD3(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: ASS(:,:), RMAT2_ASYM(:,:)

      IF (FIRST_CALL) THEN
        CALL REALLOCATE_2D_2(XRCOEFV, XICOEFV, NNN, NCH)
        CALL REALLOCATE_2D_2(BD5, BD6, NCH, NNN)
        CALL REALLOCATE_2D_9(CD0, CD1, CD2, CD3, CD4, CD5, CD6, CD7, CD8, NCH, NCH)
        CALL REALLOCATE_2D_2(ASS, RMAT2_ASYM, NCH, NCH)
        CALL REALLOCATE_2D_3(RD0, RD2, RD3, NCH, NCH)
        CALL REALLOCATE_2D_1(RCI,  NCH, NNN)
        CALL REALLOCATE_2D_1(RCIV, NNN, NCH)
        FIRST_CALL = .FALSE.
      ENDIF

      DO IAK=1,NCH
        DO IK=1,NNN
          XRCOEFV(IK,IAK) = XRCOEFF(IAK,IK)
          XICOEFV(IK,IAK) = XICOEFF(IAK,IK)
        ENDDO
      ENDDO


      DO IAB=1,NCH
        DO IK=1,NNN
          SOMMA  = ZERO
          SOMMA1 = ZERO
          DO IB=1,NNN
            SOMMA  = SOMMA  + XRCOEFF(IAB,IB)*C(IB,IK)
            SOMMA1 = SOMMA1 + XICOEFF(IAB,IB)*C(IB,IK)
          ENDDO
          BD5(IAB,IK) = SOMMA
          BD6(IAB,IK) = SOMMA1
        ENDDO
      ENDDO


      DO IAB=1,NCH
        DO IB=1,NNN
          SOMMA = ZERO
          DO IAK=1,NCH
            SOMMA = SOMMA + RMAT(IAB,IAK)*XICOEFF(IAK,IB)
          ENDDO
          RCI (IAB,IB) = SOMMA
          RCIV(IB,IAB) = RCI(IAB,IB)
        ENDDO
      ENDDO


      DO IAB=1,NCH
      DO IAK=1,NCH
        SOMMA  = ZERO
        SOMMA1 = ZERO
        SOMMA2 = ZERO
        DO IK=1,NNN
          SOMMA  = SOMMA  + BD5(IAB,IK)*XRCOEFV(IK,IAK)
          SOMMA1 = SOMMA1 + BD5(IAB,IK)*RCIV(IK,IAK)
          SOMMA2 = SOMMA2 + BD6(IAB,IK)*XRCOEFV(IK,IAK)
        ENDDO
        CD0(IAB,IAK) = SOMMA
        CD1(IAB,IAK) = SOMMA1
        RD0(IAB,IAK) = SOMMA2
      ENDDO
      ENDDO

      DO I=1,NCH
      DO IAB=1,NCH
        SOMMA  = ZERO
        SOMMA1 = ZERO
        SOMMA2 = ZERO
        SOMMA3 = ZERO
        SOMMA4 = ZERO
        SOMMA5 = ZERO
        SOMMA6 = ZERO
        DO IAK=1,NCH
          SOMMA  = SOMMA  + BD2(I,IAK) * RMAT(IAK,IAB)
          SOMMA1 = SOMMA1 + RD0(I,IAK) * RMAT(IAK,IAB)
          SOMMA2 = SOMMA2 + BD3(I,IAK) * RMAT(IAK,IAB)
          SOMMA3 = SOMMA3 + ARI(I,IAK) * RMAT(IAK,IAB)
          SOMMA4 = SOMMA4 + AIR(I,IAK) * RMAT(IAK,IAB)
          SOMMA5 = SOMMA5 + BD1(I,IAK) * RMAT(IAK,IAB)
          SOMMA6 = SOMMA6 + AII(I,IAK) * RMAT(IAK,IAB)
        ENDDO
        CD2(I,IAB) = SOMMA
        CD3(I,IAB) = SOMMA1
        CD4(I,IAB) = SOMMA2
        CD5(I,IAB) = SOMMA3
        CD6(I,IAB) = SOMMA4
        RD2(I,IAB) = SOMMA5
        RD3(I,IAB) = SOMMA6
      ENDDO
      ENDDO

      DO I=1,NCH
      DO IAB=1,NCH
        SOMMA  = ZERO
        SOMMA1 = ZERO
        DO IAK=1,NCH
          SOMMA  = SOMMA  + RD2(I,IAK) * RMAT(IAK,IAB)
          SOMMA1 = SOMMA1 + RD3(I,IAK) * RMAT(IAK,IAB)
        ENDDO
        CD7(I,IAB) = SOMMA
        CD8(I,IAB) = SOMMA1
      ENDDO
      ENDDO

      IF (PRINT_I) WRITE(*,*)
      DO IAB=1,NCH
      DO IAK=1,NCH
        ASS(IAB,IAK)=CD0(IAB,IAK)+2*BD4(IAB,IAK)+CD4(IAB,IAK)    &
                    +ARR(IAB,IAK)+CD5(IAB,IAK)+CD2(IAB,IAK)    &
                    +CD7(IAB,IAK)+CD6(IAB,IAK)    &
                    +CD8(IAB,IAK)

        RMAT2_ASYM(IAB,IAK) = RMAT(IAB,IAK) + ASS(IAB,IAK)
        IF (PRINT_I) WRITE(*,*)"COEFF RMAT2",RMAT2_ASYM(IAB,IAK), RMAT(IAB,IAK), ASS(IAB,IAK)
      ENDDO
      ENDDO

      !SYMMETRIZATION
      RMAT2 = ZERO
      DO IAB=1,NCH
      DO IAK=1,NCH
        RMAT2(IAB,IAK) =-0.5D0*( RMAT2_ASYM(IAB,IAK) + RMAT2_ASYM(IAK,IAB) )
      ENDDO
      ENDDO
    END SUBROUTINE R_SECOND_ORDER

    !> \brief Write coefficients to file to recreate the wave function.
    !! \param[in] FILE (optional) Output file name
    SUBROUTINE WRITE_COEFFICIENTS_TO_RECREATE_THE_WAVE_FUNCTION(FILE)
      IMPLICIT NONE
      INTEGER :: NNL, I
      INTEGER :: LLA(NCH_MAX), LSA(NCH_MAX), LTA(NCH_MAX), LJA(NCH_MAX)
      DOUBLE PRECISION :: EPS, GAMMA
      CHARACTER(LEN=*), OPTIONAL :: FILE

      IF (PRESENT(FILE)) THEN
        OPEN(UNIT=19, FILE=TRIM(FILE), STATUS='UNKNOWN', ACTION="WRITE")
      ELSE
        OPEN(UNIT=19, FILE='fort.19', STATUS='UNKNOWN', ACTION="WRITE")
      ENDIF

      NNL   = VAR_P%NNL
      EPS   = VAR_P%EPS
      GAMMA = VAR_P%GAMMA

      DO I = 1, NCH
        LLA(1) = GET_CHANNEL_L(CHANNEL, I)
        LSA(1) = GET_CHANNEL_S(CHANNEL, I)
        LTA(1) = GET_CHANNEL_T(CHANNEL, I)
        LJA(1) = GET_CHANNEL_J(CHANNEL)
      ENDDO

      WRITE(19,*)NCH,GAMMA,NNL
      WRITE(19,*)EPS

      DO IAK=1,NCH
        WRITE(19,*)LLA(IAK),LSA(IAK),LTA(IAK),LJA(IAK)

        DO I=1,NNN
          WRITE(19,*) XICOEFF(IAK,I)
        ENDDO

        DO I=1,NNN
          WRITE(19,*) XRCOEFF(IAK,I)
        ENDDO

      ENDDO

      DO IAK=1,NCH
      DO IAB=1,NCH
        WRITE(19,*)-RMAT2(IAK,IAB)
      ENDDO
      ENDDO
      CLOSE(19)
    END SUBROUTINE WRITE_COEFFICIENTS_TO_RECREATE_THE_WAVE_FUNCTION

    !> \brief Calculate phase shifts and mixing angle in the Blatt-Biedenharn convention.
    SUBROUTINE CALCULATE_PHASE_SHIFTS_BLATT()
      IMPLICIT NONE
      AMIXR=0.5D0*ATAN(2.*RMAT2(1,2)/(RMAT2(1,1)-RMAT2(2,2)))
      AMIXG=(AMIXR*180.D0)/PI

      DELTA1=ATAN((COS(AMIXR)*COS(AMIXR)*RMAT2(1,1)  &
                   +SIN(AMIXR)*SIN(AMIXR)*RMAT2(2,2)  &
                   +2*COS(AMIXR)*SIN(AMIXR)*RMAT2(1,2)))

      DELTA2=ATAN((SIN(AMIXR)*SIN(AMIXR)*RMAT2(1,1)  &
                   +COS(AMIXR)*COS(AMIXR)*RMAT2(2,2)  &
                   -2*COS(AMIXR)*SIN(AMIXR)*RMAT2(1,2)))
      DELTA1G=(DELTA1*180.D0)/PI
      DELTA2G=(DELTA2*180.D0)/PI

      IF (PRINT_I) WRITE(*,*)"BLATT-BIEDENHARN"
      IF (PRINT_I) WRITE(*,*)"MIXING ANGLE=",AMIXG
      IF (PRINT_I) WRITE(*,*)"SFASAMENTO1=",DELTA1G
      IF (PRINT_I) WRITE(*,*)"SFASAMENTO2=",DELTA2G

    END SUBROUTINE CALCULATE_PHASE_SHIFTS_BLATT

    !> \brief Calculate the S-matrix from the R-matrix.
    SUBROUTINE CALCULATE_S_MATRIX()
      IMPLICIT NONE
      DOUBLE PRECISION :: COS1, SIN1
      DOUBLE COMPLEX :: SM1, SM2
      SM1  = CDEXP(2.D0*IM*DELTA1)
      SM2  = CDEXP(2.D0*IM*DELTA2)
      COS1 = COS(AMIXR)
      SIN1 = SIN(AMIXR)
      SMAT(1,1) = COS1*COS1*SM1 + SIN1*SIN1*SM2
      SMAT(2,2) = COS1*COS1*SM2 + SIN1*SIN1*SM1
      SMAT(1,2) = COS1*SIN1*(SM1 - SM2)
    END SUBROUTINE CALCULATE_S_MATRIX

    !> \brief Calculate phase shifts and mixing angle in the Stapp convention.
    SUBROUTINE CALCULATE_PHASE_SHIFTS_STAPP()
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: TOL = 1.D-14
      DOUBLE COMPLEX :: SM1, SM2
      DOUBLE PRECISION :: SI2E, CI2E, CX, SX
      SM1  = REAL(SMAT(1,1)*SMAT(2,2)-SMAT(1,2)*SMAT(1,2))
      SI2E =-REAL(SMAT(1,2)*SMAT(1,2)/SM1)
      CI2E = ONE - SI2E
      IF (CI2E.LT.ZERO) THEN
        PRINT *, "Error: CI2E is negative in CALCULATE_PHASE_SHIFTS_STAPP", CI2E
        STOP
      ENDIF
      IF (SI2E.LT.ZERO) THEN
        PRINT *, "Error: SI2E is negative in CALCULATE_PHASE_SHIFTS_STAPP", SI2E
        STOP
      ENDIF
      SI2E = DSQRT(SI2E)
      CI2E = DSQRT(CI2E)

    ! I MIXING ANGLES DELLE ONDE DISPARI (JP=0) VENGONO COL SEGNO SBAGLIATO!
      SM1 = CDSQRT(SMAT(1,1)/CI2E)
      SM2 = CDSQRT(SMAT(2,2)/CI2E)
      CX  = DREAL(SM1)
      SX  = DIMAG(SM1)
      IF (ABS(CX).GT.ONE) THEN
        IF (ABS(CX-1).LT.TOL) THEN
          CX = CX/DABS(CX)
        ELSE
          PRINT *, "Error: CX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", CX, " CX - 1", CX-1
          STOP
        ENDIF 
      ENDIF
      IF (ABS(SX).GT.ONE) THEN
        IF (ABS(SX-1).LT.TOL) THEN
          SX = SX/DABS(SX)
        ELSE
          PRINT *, "Error: SX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", SX, " SX - 1", SX-1
          STOP
        ENDIF
      ENDIF
      DELTA1S = DACOS(CX)*180.D0/PI
      IF(SX.LT.ZERO) DELTA1S = -DELTA1S

      CX = DREAL(SM2)
      SX = DIMAG(SM2)
      IF (ABS(CX).GT.ONE) THEN
        IF (ABS(CX-1).LT.TOL) THEN
          CX = CX/DABS(CX)
        ELSE
          PRINT *, "Error: CX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", CX, " CX - 1", CX-1
          STOP
        ENDIF
      ENDIF
      IF (ABS(SX).GT.ONE) THEN
        IF (ABS(SX-1).LT.TOL) THEN
          SX = SX/DABS(SX)
        ELSE
          PRINT *, "Error: SX is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", SX, " SX - 1", SX-1
          STOP
        ENDIF
      ENDIF
      DELTA2S = DACOS(CX)*180.D0/PI
      IF(SX.LT.ZERO) DELTA2S = -DELTA2S

      IF (ABS(SI2E).GT.ONE) THEN
        IF (ABS(SI2E-1).LT.TOL) THEN
          SI2E = SI2E/DABS(SI2E)
        ELSE
          PRINT *, "Error: SI2E is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", SI2E, " SI2E - 1", SI2E-1
          STOP
        ENDIF
      ENDIF
      AMIXGS = (0.5D0*DASIN(SI2E))*180.D0/PI

    END SUBROUTINE CALCULATE_PHASE_SHIFTS_STAPP

  END SUBROUTINE NN_SCATTERING_VARIATIONAL


  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare core-core matrix elements for the variational calculation.
  SUBROUTINE PREPARE_CORE_CORE_MATRIX_ELEMENTS()
    USE LAGUERRE_POLYNOMIAL_MOD
    USE INTEGRATION_MOD
    USE EFT_PLESS

    IMPLICIT NONE
    INTEGER :: NNL, NNN, NEQ
    DOUBLE PRECISION :: SUM
    INTEGER :: I, NX, IE, ICH, II, JJ, IPOT
    DOUBLE PRECISION :: GAMMA
    INTEGER :: COMMON_INDEX(NCH_MAX, NNE), LIK, IAB, IAK, IL, IR, IB, IK, LR
    DOUBLE PRECISION, ALLOCATABLE :: KIN_MATRIX(:,:), POT_MATRIX(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: HCC(:,:,:), ENCC(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: INTEGRAND(:)

    GAMMA = VAR_P%GAMMA
    NNL = VAR_P%NNL

    NX = VAR_P%NX_CC
    CALL REALLOCATE_1D_1(INTEGRAND, NX)
    CALL REALLOCATE_2D_1(KIN_MATRIX, NNN_MAX, NNN_MAX)
    CALL REALLOCATE_3D_1(POT_MATRIX, NCHANNELS, NNN_MAX, NNN_MAX)
    CALL REALLOCATE_3D_1(HCC, NCHANNELS, NNN_MAX, NNN_MAX)
    CALL REALLOCATE_3D_1(ENCC, NE, NNN_MAX, NNN_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_CC, NCHANNELS, NE, NNN_MAX, NNN_MAX)

    IF (USE_DYNAMIC) THEN
      CALL REALLOCATE_4D_1(K_MINUS_E_CC, NCHANNELS, NE, NNN_MAX, NNN_MAX)
      IF (.NOT.ALLOCATED(FMAT_CC)) THEN
        ALLOCATE(FMAT_CC(0:7, NCHANNELS, NE, NNN_MAX, NNN_MAX))
      ENDIF
    ENDIF

    ENCC = ZERO
    DO IE =1, NE
      DO I = 1, NNN_MAX
        ENCC(IE,I,I) = ENERGIES_(IE)
      ENDDO
    ENDDO

    DO ICH = 1, NCHANNELS
      NEQ = GET_CHANNEL_NCH(CHANNELS_(ICH))

      II = 0
      DO I=1, NEQ
        DO JJ=1, VAR_P%NNL
          II = II + 1
          COMMON_INDEX(I, JJ) = II
        ENDDO
      ENDDO

      NNN = NEQ*NNL
      DO IAB=1, NEQ
      DO IAK=1, NEQ
        LR = GET_CHANNEL_L(CHANNELS_(ICH), IAK)
        LIK= LR*(LR+1)

        DO IL=1, NNL
        DO IR=1, NNL

          IB=COMMON_INDEX(IAB,IL)
          IK=COMMON_INDEX(IAK,IR)

    ! CALCULATING KINETIC ENERGY
          KIN_MATRIX(IB,IK) = ZERO
          IF(IAB.EQ.IAK)THEN
            INTEGRAND = V0_CC(IL,:)*( V2_CC(IR,:) + 2.D0*V1_CC(IR,:)/XX_CC &
                  -LIK*V0_CC(IR,:)/XX_CC**2 )
            SUM = ZERO
            DO I=1,NX
              SUM = SUM + XX_CC(I)**2*INTEGRAND(I)*A_CC(I)
            ENDDO
            KIN_MATRIX(IB,IK)=-HTM*SUM/GAMMA
          ENDIF

          IF (PRINT_I .AND. IB.EQ.1.AND.IK.EQ.1)THEN
            WRITE(*,*)
            WRITE(*,*)'C-C MATRIX'
            WRITE(*,*)'KINETIC',KIN_MATRIX(1,1)
          ENDIF

    ! CALCULATING POTENTIAL ENERGY
          IF (.NOT.USE_DYNAMIC) THEN
            SUM = ZERO
            INTEGRAND = V0_CC(IL,:)*V0_CC(IR,:)*V_CC(ICH,:,IAB,IAK)
            DO I=1,NX
              SUM = SUM + XX_CC(I)**2 *INTEGRAND(I) * A_CC(I)
            ENDDO
            POT_MATRIX(ICH,IB,IK) = SUM/GAMMA
          ELSE
            ! FILL FOR CALCULATING DYNAMIC POTENTIAL ENERGY
            DO IPOT=0, 7
              SUM = ZERO
              INTEGRAND = V0_CC(IL,:)*V0_CC(IR,:)*EFT_RADIAL_CC%FR_I(IPOT,:)
              DO I=1,NX
                SUM = SUM + XX_CC(I)**2 * INTEGRAND(I) * A_CC(I)
              ENDDO
              FMAT_CC(IPOT,ICH,:,IB,IK) = SUM/GAMMA
            ENDDO
          ENDIF

        ENDDO ! IR
        ENDDO ! IL

      ENDDO ! IAK
      ENDDO ! IAB
      IE = FIND_ENERGY_INDEX(VAR_P%E)
      IF (.NOT. USE_DYNAMIC) THEN
        HCC(ICH,1:NNN,1:NNN) = ( KIN_MATRIX(1:NNN,1:NNN) + POT_MATRIX(ICH,1:NNN,1:NNN) )
        DO IE = 1, NE
          H_MINUS_E_CC(ICH, IE, 1:NNN, 1:NNN) = (HCC(ICH,1:NNN,1:NNN) - ENCC(IE,1:NNN,1:NNN))/ HTM
        ENDDO
      ELSE
        DO IE = 1, NE
          K_MINUS_E_CC(ICH, IE, 1:NNN, 1:NNN) = (KIN_MATRIX(1:NNN,1:NNN) - ENCC(IE,1:NNN,1:NNN))/ HTM
        ENDDO ! IE
      ENDIF 
    ENDDO ! ICH

    ! IF (PRINT_I) WRITE(*,*)'C-C MATRIX',HTM*AM(1,1)
    DEALLOCATE(INTEGRAND)
    DEALLOCATE(KIN_MATRIX, POT_MATRIX)
    DEALLOCATE(HCC, ENCC)
  END SUBROUTINE PREPARE_CORE_CORE_MATRIX_ELEMENTS

  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare asymptotic-core matrix elements for the variational calculation.
  SUBROUTINE PREPARE_ASYMPTOTIC_CORE_MATRIX_ELEMENTS()
    USE LAGUERRE_POLYNOMIAL_MOD
    USE INTEGRATION_MOD
    IMPLICIT NONE

    DOUBLE PRECISION :: H5, HR ! Step size in r
    INTEGER :: NX, NEQ, NNN
    INTEGER :: IE, ICH, LR, II, JJ, I
    INTEGER :: IAB, IAK, LIK, IL, IB, LL, IPOT
    DOUBLE PRECISION :: AXX1, AKE1, APE, APE1
    DOUBLE PRECISION, ALLOCATABLE, SAVE :: AXXM1(:,:)
    DOUBLE PRECISION, ALLOCATABLE, SAVE :: AKEM1(:,:)
    DOUBLE PRECISION, ALLOCATABLE, SAVE :: APEM(:,:), APEM1(:,:)

    DOUBLE PRECISION, ALLOCATABLE :: AM(:,:), AM1(:,:)
    DOUBLE PRECISION, ALLOCATABLE, SAVE :: FUN(:), FUN1(:)
    INTEGER, SAVE :: COMMON_INDEX(NCH_MAX, NNE)

    CALL REALLOCATE_4D_1(H_MINUS_E_AC_R, NCHANNELS, NE, NNN_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AC_I, NCHANNELS, NE, NNN_MAX, NCH_MAX)
    CALL REALLOCATE_2D_1(AM, NNN_MAX, NCH_MAX)
    CALL REALLOCATE_2D_1(AM1, NNN_MAX, NCH_MAX)

    IF (USE_DYNAMIC) THEN
      CALL REALLOCATE_4D_1(K_MINUS_E_AC_R, NCHANNELS, NE, NNN_MAX, NCH_MAX)
      CALL REALLOCATE_4D_1(K_MINUS_E_AC_I, NCHANNELS, NE, NNN_MAX, NCH_MAX)
      IF (.NOT.ALLOCATED(FMAT_AC_R)) THEN
        ALLOCATE(FMAT_AC_R(0:7, NCHANNELS, NE, NNN_MAX, NCH_MAX))
        ALLOCATE(FMAT_AC_I(0:7, NCHANNELS, NE, NNN_MAX, NCH_MAX))
      ENDIF
    ENDIF

    HR = VAR_P%HR1
    H5 = HR/22.5D0
    NX = VAR_P%NX_AC

    IF (VAR_P%RANGE.LT.H5 .OR. VAR_P%RANGE.GT.200.D0) THEN
      PRINT *, "Error: RANGE out of bounds"
      STOP
    ENDIF

  ! Initialize grid with r values
    IF (PRINT_I) WRITE(*,*)'NX =',NX
    CALL REALLOCATE_2D_4(AXXM1, AKEM1, APEM, APEM1, NNN_MAX, NCH_MAX)
    IF (PRINT_I) WRITE(*,*)'FIRST AND LAST POINT =',XX_AC(1),XX_AC(NX),NX

      ! Prepare the indeces for the matrix elements
    CALL REALLOCATE_1D_2(FUN, FUN1, NX+1)


    ! Evaluate the matrix elements
    H_MINUS_E_AC_R = ZERO
    H_MINUS_E_AC_I = ZERO
    DO IE = 1, NE
      DO ICH = 1, NCHANNELS
        NEQ = GET_CHANNEL_NCH(CHANNELS_(ICH))
        NNN = NEQ * VAR_P%NNL

        II = 0
        DO I=1, NEQ
          DO JJ=1, VAR_P%NNL
            II = II + 1
            COMMON_INDEX(I, JJ) = II
          ENDDO
        ENDDO

        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(IAB,IAK,IL,IB,LL,LR,LIK,FUN,FUN1) SCHEDULE(static)
        DO IAB = 1, NEQ
          DO IAK = 1, NEQ
            LL = GET_CHANNEL_L(CHANNELS_(ICH), IAK)
            LR = GET_CHANNEL_L(CHANNELS_(ICH), IAB)
            LIK = LR*(LR+1)

            DO IL = 1, VAR_P%NNL
              IB = COMMON_INDEX(IAB, IL)

            ! Evaluate the normalization core-irregular (axx1)
              AXXM1(IB,IAK) = ZERO
              IF(IAB.EQ.IAK)THEN
                FUN1(1)  = ZERO
                FUN1(2:) = AJ_AC*V0_AC(IL,:)*GBES_AC(IE,LL,:)

                AXX1= ENERGIES_(IE) * B5_SINGLE(NX,H5,FUN1,1)
                AXXM1(IB,IAK)=AXX1
                ! write(111,*) iab, iak, il, IB, axx1
              ENDIF

            ! Evaluate the kinetic energy core-irregular (ake1)
              AKE1 = ZERO
              AKEM1(IB,IAK) = ZERO
              IF(IAB.EQ.IAK)THEN
                FUN1(1) = ZERO
                FUN1(2:) = AJ_AC*GBES_AC(IE,LL,:)*( V2_AC(IL,:) + 2.D0*V1_AC(IL,:)/XX_AC &
                                                            - LIK*V0_AC(IL,:)/XX_AC**2)

                AKE1 = -HTM * B5_SINGLE(NX,H5,FUN1,1)
                AKEM1(IB,IAK) = AKE1
                ! write(112,*) iab, iak, il, IB, ake1
              ENDIF
              IF(PRINT_I .AND. IB.EQ.1.AND.IAK.EQ.1)THEN
                WRITE(*,*)
                WRITE(*,*)'C-A MATRIX'
                WRITE(*,*)'IRREGULAR A'
                WRITE(*,*)'NORM ',AXXM1(1,1)
                WRITE(*,*)'KINETIC',AKEM1(1,1)
              ENDIF


            ! Evaluate the potential energy core-regular (ape), core-irregular (ape1)
              IF (.NOT. USE_DYNAMIC) THEN
                FUN (1) = ZERO
                FUN1(1) = ZERO
                FUN (2:) = AJ_AC*V0_AC(IL,:)*FBES_AC(IE,LL,:)*V_AC(ICH,:,IAB,IAK)
                FUN1(2:) = AJ_AC*V0_AC(IL,:)*GBES_AC(IE,LL,:)*V_AC(ICH,:,IAB,IAK)
                
                APE  = B5_SINGLE(NX,H5,FUN,1)
                APEM(IB,IAK) = APE

                APE1 = B5_SINGLE(NX,H5,FUN1,1)
                APEM1(IB,IAK) = APE1
              ELSE
                DO IPOT = 0, 7
                  FUN (1) = ZERO
                  FUN1(1) = ZERO
                  FUN(2:)  = AJ_AC*V0_AC(IL,:)*FBES_AC(IE,LL,:)*EFT_RADIAL_AC%FR_I(IPOT,:)
                  FUN1(2:) = AJ_AC*V0_AC(IL,:)*GBES_AC(IE,LL,:)*EFT_RADIAL_AC%FR_I(IPOT,:)
                  FMAT_AC_R(IPOT,ICH,IE,IB,IAK) = B5_SINGLE(NX,H5,FUN,1)
                  FMAT_AC_I(IPOT,ICH,IE,IB,IAK) = B5_SINGLE(NX,H5,FUN1,1)
                ENDDO
              ENDIF

              ! write(113,*) iab, iak, il, IB, ape, ape1
              IF(PRINT_I .AND. IB.EQ.1.AND.IAK.EQ.1)THEN
                WRITE(*,*)
                WRITE(*,*)'C-A MATRIX'
                WRITE(*,*)'IRREGULAR A'
                WRITE(*,*)'POTENTIAL ',APEM1(1,1)
                WRITE(*,*)'REGULAR A'
                WRITE(*,*)'POTENTIAL ',APEM(1,1)
              ENDIF

            ! Evaluate the Hamiltonian: core-regular (am), core-irregular (am1)

              AM(IB,IAK) = APEM(IB,IAK) / HTM
              AM1(IB,IAK)= (AKEM1(IB,IAK)+APEM1(IB,IAK)-AXXM1(IB,IAK)) / HTM

              IF(PRINT_I .AND. IB.EQ.1.AND.IAK.EQ.1)THEN
                WRITE(*,*)
                WRITE(*,*)'C-A MATRIX'     ,IB,IAK
                WRITE(*,*)"CORE-REGULAR=  ",AM(IB,IAK),IB,IAK
                WRITE(*,*)"CORE-IRREGULAR=",AM1(IB,IAK),IB,IAK
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
        IF (.NOT. USE_DYNAMIC) THEN ! DO NOT USE_DYNAMIC
          ! Store the matrix elements in the global arrays
          H_MINUS_E_AC_R(ICH, IE, 1:NNN, 1:NEQ) = AM(1:NNN, 1:NEQ)
          H_MINUS_E_AC_I(ICH, IE, 1:NNN, 1:NEQ) = AM1(1:NNN, 1:NEQ)
        ELSE ! USE_DYNAMIC
          K_MINUS_E_AC_R(ICH, IE, 1:NNN, 1:NEQ) = AKEM1(1:NNN, 1:NEQ) / HTM
          K_MINUS_E_AC_I(ICH, IE, 1:NNN, 1:NEQ) = AXXM1(1:NNN, 1:NEQ) / HTM
        ENDIF ! END USE_DYNAMIC
      ENDDO
    ENDDO

    DEALLOCATE(AM, AM1)
    DEALLOCATE(AXXM1, AKEM1, APEM, APEM1)
    DEALLOCATE(FUN, FUN1)
    RETURN
  END SUBROUTINE PREPARE_ASYMPTOTIC_CORE_MATRIX_ELEMENTS




  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare asymptotic-asymptotic matrix elements for the variational calculation.
  SUBROUTINE PREPARE_ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS
    USE gsl_bessel
    USE INTEGRATION_MOD
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: K_SMALL = 1.D-8 ! fm^-1

    DOUBLE PRECISION :: H, H5
    INTEGER :: IAB, IAK, IE, LL, LR, ICH, NEQ, IPOT
    DOUBLE PRECISION :: AXX, AXX3, AKE, AKE3, APE, APE1, APE2, APE3
    DOUBLE PRECISION, ALLOCATABLE :: AM(:,:), AM1(:,:), AM2(:,:), AM3(:,:)
    DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: AXXM, AXXM3, AKEM, AKEM3, APEM, APEM1, APEM2, APEM3, CHECK

    INTEGER, SAVE :: NX
    DOUBLE PRECISION, ALLOCATABLE, SAVE :: FUN(:), FUN1(:), FUN2(:), FUN3(:)

    CALL REALLOCATE_4D_1(H_MINUS_E_AA_RR, NCHANNELS, NE, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AA_IR, NCHANNELS, NE, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AA_RI, NCHANNELS, NE, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AA_II, NCHANNELS, NE, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_2D_1(AM , NCH_MAX, NCH_MAX)
    CALL REALLOCATE_2D_1(AM1, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_2D_1(AM2, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_2D_1(AM3, NCH_MAX, NCH_MAX)

    IF (USE_DYNAMIC) THEN
      CALL REALLOCATE_4D_1(K_MINUS_E_AA_RR, NCHANNELS, NE, NCH_MAX, NCH_MAX)
      CALL REALLOCATE_4D_1(K_MINUS_E_AA_IR, NCHANNELS, NE, NCH_MAX, NCH_MAX)
      CALL REALLOCATE_4D_1(K_MINUS_E_AA_RI, NCHANNELS, NE, NCH_MAX, NCH_MAX)
      CALL REALLOCATE_4D_1(K_MINUS_E_AA_II, NCHANNELS, NE, NCH_MAX, NCH_MAX)
      IF (.NOT.ALLOCATED(FMAT_AA_RR)) THEN
        ALLOCATE(FMAT_AA_RR(0:7, NCHANNELS, NE, NCH_MAX, NCH_MAX))
        ALLOCATE(FMAT_AA_IR(0:7, NCHANNELS, NE, NCH_MAX, NCH_MAX))
        ALLOCATE(FMAT_AA_RI(0:7, NCHANNELS, NE, NCH_MAX, NCH_MAX))
        ALLOCATE(FMAT_AA_II(0:7, NCHANNELS, NE, NCH_MAX, NCH_MAX))
      ENDIF
    ENDIF

    H = VAR_P%H
    H5= H/22.5D0
    NX = VAR_P%NX_AA

    CALL REALLOCATE_1D_4(FUN, FUN1, FUN2, FUN3, NX+1)

    IF (PRINT_I) WRITE(*,*)'FIRST AND LAST POINT =',XX_AA(1),XX_AA(NX)

    !SI CALCOLANO ELEMENTI MATRICE
    DO ICH = 1, NCHANNELS
      DO IE = 1, NE
        AXXM  = ZERO
        AXXM3 = ZERO
        AKEM  = ZERO
        AKEM3 = ZERO

        NEQ = GET_CHANNEL_NCH(CHANNELS_(ICH))
        DO IAB=1,NEQ
        DO IAK=1,NEQ
          LL = GET_CHANNEL_L(CHANNELS_(ICH),IAB)
          LR = GET_CHANNEL_L(CHANNELS_(ICH),IAK)
      ! SI CALCOLA NORMA DEL CASO REGOLARE-IRREGOLARE (AXX), CASO IRREGOLARE-IRREGOLARE (AXX3)
          AXX = ZERO
          AXX3= ZERO
          IF(IAB.EQ.IAK)THEN
            FUN (1) = ZERO
            FUN3(1) = ZERO
            FUN (2:) = AJ_AA*FBES_AA(IE, LL,:)*GBES_AA(IE, LR,:)
            FUN3(2:) = AJ_AA*GBES_AA(IE, LL,:)*GBES_AA(IE, LR,:)

            AXX=  ENERGIES_(IE) * B5_SINGLE(NX,H5,FUN,1)
            AXXM(IAB,IAK)=AXX
            AXX3= ENERGIES_(IE) * B5_SINGLE(NX,H5,FUN3,1)
            AXXM3(IAB,IAK)=AXX3
          ENDIF
          IF (PRINT_I .AND. IAB==1 .AND. IAK==1) THEN
            WRITE(*,*) 'NORM RI(1,1)', AXX
            WRITE(*,*) 'NORM II(1,1)', AXX3
          ENDIF

      ! SI CALCOLA ENERGIA CINETICA DEL CASO REGOLARE-IRREGOLARE (AKE) E CASO IRREGOLARE-IRREGOLARE (AKE3)
          AKE  = ZERO
          AKE3 = ZERO
          IF(IAB.EQ.IAK)THEN
            FUN(1)  = ZERO
            FUN3(1) = ZERO
            FUN (2:) = AJ_AA*FBES_AA(IE,LL,:)*HNOR_AA(IE,LR,:)*(GBES2_AA(IE,LR,:)+GBES1_AA(IE,LR,:)+GBES0_AA(IE,LR,:))
            FUN3(2:) = AJ_AA*GBES_AA(IE,LL,:)*HNOR_AA(IE,LR,:)*(GBES2_AA(IE,LR,:)+GBES1_AA(IE,LR,:)+GBES0_AA(IE,LR,:))

            AKE=  HTM * B5_SINGLE(NX,H5,FUN,1)
            AKEM(IAB,IAK)=AKE
            AKE3= HTM * B5_SINGLE(NX,H5,FUN3,1)
            AKEM3(IAB,IAK)=AKE3
          ENDIF
          IF (PRINT_I .AND. IAB==1 .AND. IAK==1) THEN
            WRITE(*,*) 'KINETIC RI(1,1)', AKE
            WRITE(*,*) 'KINETIC II(1,1)', AKE3
          ENDIF

      ! SI CALCOLA ENERGIA POTENZIALE DEL CASO REGOLARE-IRREGOLARE (APE),
      !      IRREGOLARE-REGOLARE (APE1), REGOLARE-REGOLARE (APE2), IRREGOLARE-IRREGOLARE (APE3)
          IF (.NOT.USE_DYNAMIC) THEN
            APE  = ZERO
            APE1 = ZERO
            APE2 = ZERO
            APE3 = ZERO
            FUN (1) = ZERO
            FUN1(1) = ZERO
            FUN2(1) = ZERO
            FUN3(1) = ZERO
            FUN (2:) = AJ_AA*FBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)
            FUN1(2:) = AJ_AA*GBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)
            FUN2(2:) = AJ_AA*FBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)
            FUN3(2:) = AJ_AA*GBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)

            APE=  B5_SINGLE(NX,H5,FUN,1)
            APEM(IAB,IAK)=APE

            APE1= B5_SINGLE(NX,H5,FUN1,1)
            APEM1(IAB,IAK)=APE1

            APE2= B5_SINGLE(NX,H5,FUN2,1)
            APEM2(IAB,IAK)=APE2

            APE3= B5_SINGLE(NX,H5,FUN3,1)
            APEM3(IAB,IAK)=APE3
            IF (PRINT_I .AND. IAB==1 .AND. IAK==1) THEN
              WRITE(*,*) 'POTENTIAL RR(1,1)', APE2
              WRITE(*,*) 'POTENTIAL RI(1,1)', APE
              WRITE(*,*) 'POTENTIAL IR(1,1)', APE1
              WRITE(*,*) 'POTENTIAL II(1,1)', APE3
            ENDIF
          ELSE
            ! FILL HERE TO CALCULATE DYNAMIC POTENTIAL ENERGY
            DO IPOT = 0, 7
              FUN (1) = ZERO
              FUN1(1) = ZERO
              FUN2(1) = ZERO
              FUN3(1) = ZERO
              FUN (2:) = AJ_AA*FBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(IPOT,:)
              FUN1(2:) = AJ_AA*GBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(IPOT,:)
              FUN2(2:) = AJ_AA*FBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(IPOT,:)
              FUN3(2:) = AJ_AA*GBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(IPOT,:)
              
              FMAT_AA_RR(IPOT,ICH,IE,IAB,IAK) = B5_SINGLE(NX,H5,FUN2,1)
              FMAT_AA_RI(IPOT,ICH,IE,IAB,IAK) = B5_SINGLE(NX,H5,FUN,1)
              FMAT_AA_IR(IPOT,ICH,IE,IAB,IAK) = B5_SINGLE(NX,H5,FUN1,1)
              FMAT_AA_II(IPOT,ICH,IE,IAB,IAK) = B5_SINGLE(NX,H5,FUN3,1)
            ENDDO
          ENDIF

      ! SI CALCOLA HAMILTONIANA PER I VARI CASI:REGOLARE-IRREGOLARE(AM), IRREGOLARE-REGOLARE(AM1),
      !         REGOLARE-REGOLARE(AM2),IRREGOLARE-IRREGOLARE(AM3)
          AM   (IAB,IAK) = (AKEM(IAB,IAK)+APEM(IAB,IAK)-AXXM(IAB,IAK)) / HTM
          AM1  (IAB,IAK) =  APEM1(IAB,IAK) / HTM
          AM2  (IAB,IAK) =  APEM2(IAB,IAK) / HTM
          AM3  (IAB,IAK) = (AKEM3(IAB,IAK)+APEM3(IAB,IAK)-AXXM3(IAB,IAK)) / HTM

          CHECK(IAB,IAK) = (AM1(IAB,IAK)-AM(IAB,IAK))
        ENDDO !IAB
        ENDDO !IAK
        IF (.NOT.USE_DYNAMIC) THEN ! NOT USE_DYNAMIC
          H_MINUS_E_AA_RR(ICH, IE, :, :) = AM2
          H_MINUS_E_AA_RI(ICH, IE, :, :) = AM
          H_MINUS_E_AA_IR(ICH, IE, :, :) = AM1
          H_MINUS_E_AA_II(ICH, IE, :, :) = AM3
        ELSE ! USE_DYNAMIC
          K_MINUS_E_AA_RR(ICH, IE, :, :) = (AKEM - AXXM) / HTM
          K_MINUS_E_AA_RI(ICH, IE, :, :) =  0.D0
          K_MINUS_E_AA_IR(ICH, IE, :, :) =  0.D0
          K_MINUS_E_AA_II(ICH, IE, :, :) = (AKEM3 - AXXM3) / HTM
        ENDIF ! ENDIF USE_DYNAMIC
      ENDDO ! IE
    ENDDO ! ICH

    IF (PRINT_I) WRITE(*,*)
    IF (PRINT_I) WRITE(*,*)'A-A MATRIX'
    DO IAB=1,NEQ
    DO IAK=1,NEQ
      IF (PRINT_I) WRITE(*,'(2I6,4D17.7)') IAB,IAK,AM2(IAB,IAK),AM(IAB,IAK),AM1(IAB,IAK),AM3(IAB,IAK)
      IF (PRINT_I) WRITE(*,*)"1=",CHECK(IAB,IAK)
    END DO
    END DO

    DEALLOCATE(FUN, FUN1, FUN2, FUN3)
    DEALLOCATE(AM, AM1, AM2, AM3)

  END SUBROUTINE PREPARE_ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS

  !> \ingroup scattering_nn_variational_mod
  !> \brief Print a divider line to the output.
  SUBROUTINE PRINT_DIVIDER()
    IMPLICIT NONE
      WRITE(*,*) '====================================================================================='
  END SUBROUTINE PRINT_DIVIDER

  !> \ingroup scattering_nn_variational_mod
  !> \brief Check if this is the first call with a given set of quantum numbers and parameters.
  !! \param[in] J, L, S, TZ, IPOT, ILB, LEMP
  !! \return .TRUE. if first call, .FALSE. otherwise
  FUNCTION IS_FIRST_CALL(J, L, S, TZ, IPOT, ILB, LEMP) RESULT(FIRST_CALL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP
    INTEGER :: JOLD = -1, LOLD = -1, SOLD = -1, TZOLD = -1, IPOTOLD = -1, ILBOLD = -1, LEMPOOLD = -1
    LOGICAL :: FIRST_CALL

    ! Check if the current call is the first call based on the input parameters
    FIRST_CALL = ( J /= JOLD .OR. L /= LOLD .OR. S /= SOLD .OR. TZ /= TZOLD .OR. IPOT /= IPOTOLD .OR. &
                  ILB /= ILBOLD .OR. LEMP /= LEMPOOLD )
    IF (FIRST_CALL) THEN
      JOLD = J
      LOLD = L
      SOLD = S
      TZOLD = TZ
      IPOTOLD = IPOT
      ILBOLD = ILB
      LEMPOOLD = LEMP
    ENDIF
  END FUNCTION IS_FIRST_CALL

  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare the radial grids for the calculation.
  SUBROUTINE PREPARE_GRID()
    USE INTEGRATION_MOD
    IMPLICIT NONE
    INTEGER :: NX, I
    DOUBLE PRECISION :: RANGE, EPS, GAMMA
    DOUBLE PRECISION :: HR

    ! A-A
    RANGE = VAR_P%RANGE
    CALL EXPONENTIALLY_GROWING_GRID(VAR_P%H, VAR_P%AF, RANGE, XX_AA, AJ_AA, NX)
    ALLOCATE(A_AA(NX), B_AA(NX))
    EPS = VAR_P%EPS
    A_AA  = ONE - DEXP(-EPS*XX_AA)
    B_AA  = EPS * DEXP(-EPS*XX_AA)
    VAR_P%NX_AA = NX

    ! A-C
    HR = VAR_P%HR1
    NX = INT(VAR_P%RANGE/VAR_P%HR1) + 10
    VAR_P%NX_AC = NX

    CALL REALLOCATE_1D_2(XX_AC, A_AC, NX)
    GAMMA = VAR_P%GAMMA
    XX_AC = HR * [(I, I=1,NX)]
    AJ_AC  = XX_AC**2
    YYL_AC = GAMMA*XX_AC
    A_AC   = ONE - DEXP(-EPS*XX_AC)
    VAR_P%NX_AC = NX

    ! C-C
    NX = VAR_P%NX_CC
    ALLOCATE(XX_CC(NX), YY_CC(NX), A_CC(NX))
    CALL GAULAG(NX, YY_CC,A_CC)
    XX_CC = YY_CC/GAMMA

    GRID_SET = .TRUE.

    CALL PREPARE_LAGUERRE

  END SUBROUTINE PREPARE_GRID

  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare Bessel functions for the asymptotic region.
  SUBROUTINE PREPARE_ASYMPTOTIC_FUNCTIONS
    USE gsl_bessel
    USE INTEGRATION_MOD
    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: K_SMALL = 1.D-8 ! fm^-1
    INTEGER :: NX, IE, IX, L
    DOUBLE PRECISION :: K, FBSS, GBSS, GBSS1, GBSS2
    DOUBLE PRECISION :: AG, BG, XG, EPS

    IF (.NOT.TZ_SET .OR. .NOT.LMAX_SET) THEN
      PRINT *, "Error: TZ or LMAX not set"
      STOP
    ENDIF

    IF (.NOT.ENERGIES_SET) THEN
      PRINT *, "Error: energies not set"
      STOP
    ENDIF

    CALL SET_M_T1Z_T2Z_HTM(VAR_P%TZ)
    IF (.NOT.GRID_SET) CALL PREPARE_GRID

    WRITE(*,*)'PREPARING BESSEL FUNCTIONS'

    KK = DSQRT(ENERGIES_/HTM)
    K2 = KK**2

    NX = VAR_P%NX_AA
    EPS= VAR_P%EPS
    ALLOCATE(FBES_AA(NE, 0:LMAX, NX), GBES_AA(NE, 0:LMAX, NX))
    ALLOCATE(GBES0_AA(NE, 0:LMAX, NX), GBES1_AA(NE, 0:LMAX, NX), GBES2_AA(NE, 0:LMAX, NX), HNOR_AA(NE, 0:LMAX, NX))

    DO IE = 1, NE
      K = KK(IE)
      DO IX = 1, NX
        XG = XX_AA(IX)*K
        AG = A_AA(IX)
        BG = B_AA(IX)
        DO L = 0, LMAX
          IF(K.LE.K_SMALL) THEN
            FBES_AA (IE,L,IX) = XX_AA(IX)**L
            GBES_AA (IE,L,IX) =-ONE/((2*L+ONE)*XX_AA(IX)**(L+ONE))*AG**(2*L+ONE)
            GBES0_AA(IE,L,IX) = ONE/((2*L+ONE)*XX_AA(IX)**(L+ONE))*(EPS*BG*(2*L+ONE)*((2*L+ONE)*(BG/EPS)-ONE) &
                                    +2*(2*L+ONE)*BG*(AG/XX_AA(IX))-L*(L+ONE)*(AG/XX_AA(IX))**2)
            GBES1_AA(IE,L,IX) =-2.*AG*((2*L+ONE)*BG+AG/XX_AA(IX))*(L+ONE)/((2*L+ONE)*XX_AA(IX)**(L+2.))
            GBES2_AA(IE,L,IX) = AG**2*(L+ONE)*(L+2.)/((2*L+ONE)*XX_AA(IX)**(L+3.))
            HNOR_AA (IE,L,IX) = AG**(2*L-ONE)
          ELSE
            FBSS = SPHERICAL_J  (L, XG)
            GBSS = SPHERICAL_Y  (L, XG)
            GBSS1= SPHERICAL_YP (L, XG)
            GBSS2= SPHERICAL_YPP(L, XG)

            FBES_AA (IE,L,IX) = K**(L+0.5D0)*FBSS/(K**L)
            GBES_AA (IE,L,IX) =-(GBSS*K**(L+ONE)*AG**(2*L+ONE))/(K**(L+0.5D0))
            GBES0_AA(IE,L,IX) = GBSS*(EPS*BG*(2*L+ONE)*((2*L+ONE)*(BG/EPS)-ONE) &
                                      +2*(2*L+ONE)*BG*(AG/XX_AA(IX))-L*(L+ONE)*(AG/XX_AA(IX))**2)
            GBES1_AA(IE,L,IX) = GBSS1*2.*K*AG*((2*L+ONE)*BG + AG/XX_AA(IX))
            GBES2_AA(IE,L,IX) = (K**2)*(AG**2)*GBSS2
            HNOR_AA (IE,L,IX) = (K**(L+ONE))*(AG**(2*L-ONE))/(K**(L+0.5D0))
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    NX = VAR_P%NX_AC

    ALLOCATE(FBES_AC(NE, 0:LMAX, NX), GBES_AC(NE, 0:LMAX, NX))
    DO IE = 1, NE
      K = KK(IE)
      DO L = 0, LMAX
        !$OMP PARALLEL DO PRIVATE(IX, XG, AG, FBSS, GBSS) SHARED(FBES_AC, GBES_AC, XX_AC, A_AC, K, L, IE)
        DO IX = 1, NX
          XG = XX_AC(IX)*K
          AG = A_AC(IX)
          IF(K.LE.K_SMALL)THEN                                   !(K->0)
            FBES_AC(IE,L,IX) = XX_AC(IX)**L
            GBES_AC(IE,L,IX) =-ONE/((2*L+ONE)*XX_AC(IX)**(L+ONE))*AG**(2*L+ONE)
          ELSE
            FBSS = SPHERICAL_J(L, XG)
            GBSS = SPHERICAL_Y(L, XG)
            FBES_AC(IE,L,IX) = K**(L+0.5D0)*FBSS/(K**L)
            GBES_AC(IE,L,IX) =-(GBSS*K**(L+ONE)*AG**(2*L+ONE))/(K**(L+0.5D0))
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ENDDO
      BESSELS_SET = .TRUE.
    ENDDO

    WRITE(*,*)'BESSEL FUNCTIONS PREPARED'
  END SUBROUTINE PREPARE_ASYMPTOTIC_FUNCTIONS

  !> \ingroup scattering_nn_variational_mod
  !> \brief Find the index of a given energy in the ENERGIES array.
  !! \param[in] E Energy value
  !! \return Index in ENERGIES_ array
  FUNCTION FIND_ENERGY_INDEX(E) RESULT(IE)
    IMPLICIT NONE
    INTEGER :: IE
    DOUBLE PRECISION, INTENT(IN) :: E
    INTEGER :: I

    ! Find the index of the energy in the ENERGIES array
    DO I = 1, SIZE(ENERGIES_)
      IF (ABS(E - ENERGIES_(I)) < 1.D-10) THEN
        IE = I
        RETURN
      ENDIF
    ENDDO

    PRINT *, "Error: Energy not found in ENERGIES array"
    STOP
  END FUNCTION FIND_ENERGY_INDEX

  !> \ingroup scattering_nn_variational_mod
  !> \brief Set mass and isospin projections for the calculation.
  !! \param[in] TZ Isospin projection
  SUBROUTINE SET_M_T1Z_T2Z_HTM(TZ)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: TZ
    SELECT CASE (TZ)
      CASE (1)
        M = MP
        VAR_P%T1Z = 1
        VAR_P%T2Z = 1
      CASE (0)
        M = MP * MN / (MP + MN)
        VAR_P%T1Z = 1
        VAR_P%T2Z =-1
      CASE (-1)
        M = MN
        VAR_P%T1Z =-1
        VAR_P%T2Z =-1
      CASE DEFAULT
        PRINT *, "Invalid TZ value"
        STOP
    END SELECT

    HTM = HC**2 / (2 * M)
    HTM_SET = .TRUE.
  END SUBROUTINE SET_M_T1Z_T2Z_HTM

  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare the potential matrices for all channels.
  !! \param[in] CHANNELS Array of SCATTERING_CHANNEL structures
  SUBROUTINE PREPARE_POTENTIAL(CHANNELS)
    USE POTENTIALS
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNELS(:)
    INTEGER :: NC, NEQ_C, ICH, IR, L, S, J, T, TZ, T1Z, T2Z
    LOGICAL :: COUPLED
    DOUBLE PRECISION :: R, V2(2,2)

    IF (VAR_P%IPOT.EQ.0) THEN
      PRINT *, "Error: IPOT not set"
      RETURN
    ENDIF

    IF (.NOT.GRID_SET) CALL PREPARE_GRID

    NC = SIZE(CHANNELS)
    
    IF (.NOT.USE_DYNAMIC) THEN
      ALLOCATE(V_CC(NC, VAR_P%NX_CC, 2, 2), V_AC(NC, VAR_P%NX_AC, 2, 2), V_AA(NC, VAR_P%NX_AA, 2, 2))


      DO ICH = 1, NC
        NEQ_C = GET_CHANNEL_NCH(CHANNELS(ICH))
        J   = GET_CHANNEL_J(CHANNELS(ICH))
        TZ  = GET_CHANNEL_TZ(CHANNELS(ICH))
        COUPLED = IS_CHANNEL_COUPLED(CHANNELS(ICH))
        L   = GET_CHANNEL_L(CHANNELS(ICH), 1)
        S   = GET_CHANNEL_S(CHANNELS(ICH), 1)
        T   = GET_CHANNEL_T(CHANNELS(ICH), 1)
        CALL EVAL_T1Z_T2Z

        DO IR = 1, VAR_P%NX_CC
          R = XX_CC(IR)
          CALL POT_PW(VAR_P%IPOT, VAR_P%ILB, VAR_P%LEMP, L, S, J, T1Z, T2Z, R, V2)
          IF (COUPLED) THEN
            V_CC(ICH,IR,:,:) = V2
          ELSEIF (NEQ_C == 2) THEN
            V_CC(ICH,IR,:,:) = ZERO
            V_CC(ICH,IR,1,1) = V2(1,1)
            L   = GET_CHANNEL_L(CHANNELS(ICH), 2)
            S   = GET_CHANNEL_S(CHANNELS(ICH), 2)
            T   = GET_CHANNEL_T(CHANNELS(ICH), 2)
            CALL POT_PW(VAR_P%IPOT, VAR_P%ILB, VAR_P%LEMP, L, S, J, T1Z, T2Z, R, V2)
            V_CC(ICH,IR,2,2) = V2(1,1)
          ELSE
            V_CC(ICH,IR,:,:) = ZERO
            V_CC(ICH,IR,1,1) = V2(1,1)
          ENDIF
        ENDDO

        DO IR = 1, VAR_P%NX_AC
          R = XX_AC(IR)
          CALL POT_PW(VAR_P%IPOT, VAR_P%ILB, VAR_P%LEMP, L, S, J, T1Z, T2Z, R, V2)
          IF (COUPLED) THEN
            V_AC(ICH,IR,:,:) = V2
            ELSEIF (NEQ_C == 2) THEN
              V_AC(ICH,IR,:,:) = ZERO
              V_AC(ICH,IR,1,1) = V2(1,1)
              L   = GET_CHANNEL_L(CHANNELS(ICH), 2)
              S   = GET_CHANNEL_S(CHANNELS(ICH), 2)
              T   = GET_CHANNEL_T(CHANNELS(ICH), 2)
              CALL POT_PW(VAR_P%IPOT, VAR_P%ILB, VAR_P%LEMP, L, S, J, T1Z, T2Z, R, V2)
              V_AC(ICH,IR,2,2) = V2(1,1)
            ELSE
              V_AC(ICH,IR,:,:) = ZERO
              V_AC(ICH,IR,1,1) = V2(1,1)
            ENDIF
            ! if (r<=20) WRITE(1000+(2*S+1)*100+L*10+J,*) R, V2(1,1), V2(1,2), V2(2,1), V2(2,2)
        ENDDO

        DO IR = 1, VAR_P%NX_AA
          R = XX_AA(IR)
          CALL POT_PW(VAR_P%IPOT, VAR_P%ILB, VAR_P%LEMP, L, S, J, T1Z, T2Z, R, V2)
          IF (COUPLED) THEN
            V_AA(ICH,IR,:,:) = V2
          ELSEIF (NEQ_C == 2) THEN
            V_AA(ICH,IR,:,:) = ZERO
            V_AA(ICH,IR,1,1) = V2(1,1)
            L   = GET_CHANNEL_L(CHANNELS(ICH), 2)
            S   = GET_CHANNEL_S(CHANNELS(ICH), 2)
            T   = GET_CHANNEL_T(CHANNELS(ICH), 2)
            CALL POT_PW(VAR_P%IPOT, VAR_P%ILB, VAR_P%LEMP, L, S, J, T1Z, T2Z, R, V2)
            V_AA(ICH,IR,2,2) = V2(1,1)
          ELSE
            V_AA(ICH,IR,:,:) = ZERO
            V_AA(ICH,IR,1,1) = V2(1,1)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      ! FILL TO PREPARE POTENTIALS FOR DYNAMIC CASE
      IF (.NOT.LECS_SET) THEN
        PRINT *, "Error: LECS_NOT_SET"
        STOP
      ENDIF
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_CC, LECS%RC(VAR_P%S,VAR_P%T), EFT_RADIAL_CC, ORDER_POTENTIAL = LECS%ORDER)
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_AC, LECS%RC(VAR_P%S,VAR_P%T), EFT_RADIAL_AC, ORDER_POTENTIAL = LECS%ORDER)
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_AA, LECS%RC(VAR_P%S,VAR_P%T), EFT_RADIAL_AA, ORDER_POTENTIAL = LECS%ORDER)
    ENDIF

    POTENTIAL_SET = .TRUE.
    ! stop

  CONTAINS
    !> \brief Evaluate T1Z and T2Z from TZ.
    SUBROUTINE EVAL_T1Z_T2Z()
      IMPLICIT NONE
      SELECT CASE (TZ)
        CASE (1)
          T1Z = 1
          T2Z = 1
        CASE (0)
          T1Z = 1
          T2Z =-1
        CASE (-1)
          T1Z =-1
          T2Z =-1
        CASE DEFAULT
          PRINT *, "Invalid TZ value"
          STOP
      END SELECT
    END SUBROUTINE EVAL_T1Z_T2Z
  END SUBROUTINE PREPARE_POTENTIAL

  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare Laguerre basis functions and their derivatives.
  SUBROUTINE PREPARE_LAGUERRE()
    USE LAGUERRE_POLYNOMIAL_MOD
    IMPLICIT NONE
    INTEGER :: I, J, NMX, NX
    DOUBLE PRECISION :: APF, GAMMA, XG, FEXP, ANJ
    DOUBLE PRECISION, ALLOCATABLE :: U0(:,:), U1(:,:), U2(:,:)

    IF (.NOT.GRID_SET) CALL PREPARE_GRID

    NMX = VAR_P%NNL - 1
    APF = 2.D0

    NX = VAR_P%NX_CC
    ALLOCATE(U0(0:NMX, NX), U1(0:NMX, NX), U2(0:NMX, NX))
    ALLOCATE(V0_CC(NMX+1, NX), V1_CC(NMX+1, NX), V2_CC(NMX+1, NX))
    CALL LAGUERRE_POLYNOMIAL(YY_CC, APF, U0, U1, U2)
    GAMMA = VAR_P%GAMMA
    DO I = 1, NX
      XG = YY_CC(I)
      DO J = 0, NMX
        ANJ = DSQRT(DGAMMA(J+ONE)*GAMMA**3/DGAMMA(J+3.D0))  ! Here no need for the exp, it is in the weight

        V0_CC(J+1,I) = ANJ * U0(J,I)
        V1_CC(J+1,I) = ANJ * GAMMA * (U1(J,I) - 0.5D0*U0(J,I))
        V2_CC(J+1,I) = ANJ * GAMMA * (GAMMA * (U2(J,I) - 0.5D0*U1(J,I))) &
                  - 0.5D0*GAMMA*V1_CC(J+1,I)
      ENDDO
    ENDDO
    DEALLOCATE(U0, U1, U2)

    NX = VAR_P%NX_AC
    ALLOCATE(U0(0:NMX, NX), U1(0:NMX, NX), U2(0:NMX, NX))
    ALLOCATE(V0_AC(NMX+1, NX), V1_AC(NMX+1, NX), V2_AC(NMX+1, NX))
    CALL LAGUERRE_POLYNOMIAL(YYL_AC, APF, U0, U1, U2)
    DO I = 1, NX
      XG = YYL_AC(I)
      FEXP = DEXP(-XG/2.D0)
      DO J = 0, NMX
        ANJ = DSQRT(DGAMMA(J+ONE)*GAMMA**3/DGAMMA(J+3.D0))*FEXP

        V0_AC(J+1,I) = ANJ * U0(J,I)
        V1_AC(J+1,I) = ANJ * GAMMA * (U1(J,I) - 0.5D0*U0(J,I))
        V2_AC(J+1,I) = ANJ * GAMMA * (GAMMA * (U2(J,I) - 0.5D0*U1(J,I))) &
                  - 0.5D0*GAMMA*V1_AC(J+1,I)
      ENDDO
    ENDDO
    DEALLOCATE(U0, U1, U2)

    LAGUERRE_SET = .TRUE.
  END SUBROUTINE PREPARE_LAGUERRE


  !> \ingroup scattering_nn_variational_mod
  !> \brief Find the index of the current channel in the CHANNELS array.
  !! \return Index in CHANNELS_ array
  FUNCTION FIND_CHANNEL_INDEX() RESULT(INDX)
    USE QUANTUM_NUMBERS
    IMPLICIT NONE
    INTEGER :: INDX

    ! Find the index of the channel in the CHANNELS array
    INDX = 0
    DO INDX = 1, SIZE(CHANNELS_)
      IF (IS_SAME_CHANNEL(CHANNELS_(INDX),CHANNEL)) THEN
        RETURN
      ENDIF
    ENDDO

    PRINT *, "Error: Channel not found in CHANNELS array"
    STOP
  END FUNCTION FIND_CHANNEL_INDEX

  !> @brief Returns the Hamiltonian time constant (HTM).
  !>
  !> This function retrieves the value of the Hamiltonian time constant (HTM).
  !> If the HTM value has not been set, an error message is printed and the program stops.
  !>
  !> @return The value of the Hamiltonian time constant (HTM).
  FUNCTION GET_HTM() RESULT(HTM_VALUE)
    IMPLICIT NONE
    DOUBLE PRECISION :: HTM_VALUE

    ! Return the Hamiltonian time constant
    IF (.NOT.HTM_SET) THEN
      PRINT *, "Error: HTM not set"
      STOP
    ENDIF

    HTM_VALUE = HTM
  END FUNCTION GET_HTM

  !> @brief Sets the dynamic flag for the calculation.
  !>
  !> This subroutine assigns the value of the input logical variable FLAG to the module variable USE_DYNAMIC,
  !> controlling whether dynamic behavior is enabled or disabled in the calculation.
  !>
  !> @param[in] FLAG Logical flag to enable (TRUE) or disable (FALSE) dynamic behavior.
  SUBROUTINE SET_DYNAMIC(FLAG)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: FLAG
    USE_DYNAMIC = FLAG
  END SUBROUTINE SET_DYNAMIC

  !> @brief Sets new Low Energy Constants (LECs) for the EFT potential.
  !>
  !> This subroutine updates the global LECs variable with the provided values.
  !> If the regulator cutoffs have changed, it recalculates the potential functions,
  !> their integrals, and prepares all necessary matrix elements and asymptotic functions.
  !>
  !> @param[in] LECS_NEW  New set of LECs (type LECS_EFT_PLESS) to be used.
  !>
  !> The subroutine also sets flags to indicate that the LECs and potential are up-to-date.
  !> If required, it prepares Bessel and Laguerre functions and matrix elements for the new potential.
  SUBROUTINE SET_NEW_LECS(LECS_NEW)
    IMPLICIT NONE
    TYPE(LECS_EFT_PLESS), INTENT(IN) :: LECS_NEW
    LOGICAL :: NEW_CUTOFFS
    IF (ANY(LECS_NEW%RC /= LECS%RC)) THEN
      NEW_CUTOFFS = .TRUE.
    ELSE
      NEW_CUTOFFS = .FALSE.
    ENDIF
    LECS = LECS_NEW
    LECS_SET = .TRUE.
    IF (NEW_CUTOFFS) THEN
      WRITE(*,*) "EVALUATING THE NEW POTENTIAL FUNCTIONS AND THEIR INTEGRALS"
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_CC, LECS%RC(VAR_P%S,VAR_P%T), EFT_RADIAL_CC, ORDER_POTENTIAL = LECS%ORDER)
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_AC, LECS%RC(VAR_P%S,VAR_P%T), EFT_RADIAL_AC, ORDER_POTENTIAL = LECS%ORDER)
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_AA, LECS%RC(VAR_P%S,VAR_P%T), EFT_RADIAL_AA, ORDER_POTENTIAL = LECS%ORDER)
      POTENTIAL_SET = .TRUE.
      IF (.NOT.BESSELS_SET) THEN
        CALL PREPARE_ASYMPTOTIC_FUNCTIONS
      ENDIF
      IF (.NOT.LAGUERRE_SET) THEN
        CALL PREPARE_LAGUERRE
      ENDIF
      CALL PREPARE_CORE_CORE_MATRIX_ELEMENTS
      CALL PREPARE_ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS
      CALL PREPARE_ASYMPTOTIC_FUNCTIONS
      WRITE(*,*) "NEW POTENTIAL FUNCTIONS AND THEIR INTEGRALS EVALUATED"
    ENDIF
    NEW_LECS = .TRUE.
  END SUBROUTINE SET_NEW_LECS

  !> \brief Reset the SCATTERING_NN_VARIATIONAL module to its initial state.
  !! This subroutine deallocates all allocatable arrays, resets logical flags,
  !! and restores all module variables to their default values.
  SUBROUTINE RESET_SCATTERING_NN_VARIATIONAL()
    IMPLICIT NONE

    WRITE(*,*) 'RESETTING SCATTERING_NN_VARIATIONAL MODULE'

    ! Deallocate all allocatable arrays if allocated
    IF (ALLOCATED(CHANNELS_))      DEALLOCATE(CHANNELS_)
    IF (ALLOCATED(XX_CC))          DEALLOCATE(XX_CC)
    IF (ALLOCATED(YY_CC))          DEALLOCATE(YY_CC)
    IF (ALLOCATED(XX_AC))          DEALLOCATE(XX_AC)
    IF (ALLOCATED(XX_AA))          DEALLOCATE(XX_AA)
    IF (ALLOCATED(A_CC))           DEALLOCATE(A_CC)
    IF (ALLOCATED(A_AC))           DEALLOCATE(A_AC)
    IF (ALLOCATED(A_AA))           DEALLOCATE(A_AA)
    IF (ALLOCATED(B_AA))           DEALLOCATE(B_AA)
    IF (ALLOCATED(AJ_AA))          DEALLOCATE(AJ_AA)
    IF (ALLOCATED(AJ_AC))          DEALLOCATE(AJ_AC)
    IF (ALLOCATED(YYL_AC))         DEALLOCATE(YYL_AC)
    IF (ALLOCATED(V_CC))           DEALLOCATE(V_CC)
    IF (ALLOCATED(V_AC))           DEALLOCATE(V_AC)
    IF (ALLOCATED(V_AA))           DEALLOCATE(V_AA)
    IF (ALLOCATED(ENERGIES_))      DEALLOCATE(ENERGIES_)
    IF (ALLOCATED(KK))             DEALLOCATE(KK)
    IF (ALLOCATED(K2))             DEALLOCATE(K2)
    IF (ALLOCATED(FBES_AA))        DEALLOCATE(FBES_AA)
    IF (ALLOCATED(FBES_AC))        DEALLOCATE(FBES_AC)
    IF (ALLOCATED(GBES_AA))        DEALLOCATE(GBES_AA)
    IF (ALLOCATED(GBES_AC))        DEALLOCATE(GBES_AC)
    IF (ALLOCATED(GBES0_AA))       DEALLOCATE(GBES0_AA)
    IF (ALLOCATED(GBES1_AA))       DEALLOCATE(GBES1_AA)
    IF (ALLOCATED(GBES2_AA))       DEALLOCATE(GBES2_AA)
    IF (ALLOCATED(HNOR_AA))        DEALLOCATE(HNOR_AA)
    IF (ALLOCATED(V0_CC))          DEALLOCATE(V0_CC)
    IF (ALLOCATED(V1_CC))          DEALLOCATE(V1_CC)
    IF (ALLOCATED(V2_CC))          DEALLOCATE(V2_CC)
    IF (ALLOCATED(V0_AC))          DEALLOCATE(V0_AC)
    IF (ALLOCATED(V1_AC))          DEALLOCATE(V1_AC)
    IF (ALLOCATED(V2_AC))          DEALLOCATE(V2_AC)
    IF (ALLOCATED(H_MINUS_E_CC))   DEALLOCATE(H_MINUS_E_CC)
    IF (ALLOCATED(H_MINUS_E_AC_R)) DEALLOCATE(H_MINUS_E_AC_R)
    IF (ALLOCATED(H_MINUS_E_AC_I)) DEALLOCATE(H_MINUS_E_AC_I)
    IF (ALLOCATED(H_MINUS_E_AA_RR))DEALLOCATE(H_MINUS_E_AA_RR)
    IF (ALLOCATED(H_MINUS_E_AA_RI))DEALLOCATE(H_MINUS_E_AA_RI)
    IF (ALLOCATED(H_MINUS_E_AA_IR))DEALLOCATE(H_MINUS_E_AA_IR)
    IF (ALLOCATED(H_MINUS_E_AA_II))DEALLOCATE(H_MINUS_E_AA_II)
    IF (ALLOCATED(K_MINUS_E_CC))   DEALLOCATE(K_MINUS_E_CC)
    IF (ALLOCATED(K_MINUS_E_AC_R)) DEALLOCATE(K_MINUS_E_AC_R)
    IF (ALLOCATED(K_MINUS_E_AC_I)) DEALLOCATE(K_MINUS_E_AC_I)
    IF (ALLOCATED(K_MINUS_E_AA_RR))DEALLOCATE(K_MINUS_E_AA_RR)
    IF (ALLOCATED(K_MINUS_E_AA_RI))DEALLOCATE(K_MINUS_E_AA_RI)
    IF (ALLOCATED(K_MINUS_E_AA_IR))DEALLOCATE(K_MINUS_E_AA_IR)
    IF (ALLOCATED(K_MINUS_E_AA_II))DEALLOCATE(K_MINUS_E_AA_II)
    IF (ALLOCATED(VM_CC))          DEALLOCATE(VM_CC)
    IF (ALLOCATED(VM_AC_R))        DEALLOCATE(VM_AC_R)
    IF (ALLOCATED(VM_AC_I))        DEALLOCATE(VM_AC_I)
    IF (ALLOCATED(VM_AA_RR))       DEALLOCATE(VM_AA_RR)
    IF (ALLOCATED(VM_AA_RI))       DEALLOCATE(VM_AA_RI)
    IF (ALLOCATED(VM_AA_IR))       DEALLOCATE(VM_AA_IR)
    IF (ALLOCATED(VM_AA_II))       DEALLOCATE(VM_AA_II)
    IF (ALLOCATED(FMAT_CC))        DEALLOCATE(FMAT_CC)
    IF (ALLOCATED(FMAT_AC_R))      DEALLOCATE(FMAT_AC_R)
    IF (ALLOCATED(FMAT_AC_I))      DEALLOCATE(FMAT_AC_I)
    IF (ALLOCATED(FMAT_AA_RR))     DEALLOCATE(FMAT_AA_RR)
    IF (ALLOCATED(FMAT_AA_RI))     DEALLOCATE(FMAT_AA_RI)
    IF (ALLOCATED(FMAT_AA_IR))     DEALLOCATE(FMAT_AA_IR)
    IF (ALLOCATED(FMAT_AA_II))     DEALLOCATE(FMAT_AA_II)

    ! Reset logical flags
    USE_DYNAMIC     = .FALSE.
    HTM_SET         = .FALSE.
    PRINT_I         = .FALSE.
    GRID_SET        = .FALSE.
    POTENTIAL_SET   = .FALSE.
    IPOT_SET        = .FALSE.
    ENERGIES_SET    = .FALSE.
    BESSELS_SET     = .FALSE.
    LAGUERRE_SET    = .FALSE.
    CHANNELS_SET    = .FALSE.
    LECS_SET        = .FALSE.
    NEW_LECS        = .TRUE.
    TZ_SET          = .FALSE.
    LMAX_SET        = .FALSE.

    ! Reset integer and real variables
    NCH         = 0
    NNN_MAX     = 0
    LMAX        = -1
    NE          = -1
    NCHANNELS   = 0
    CH_INDEX    = 0
    HTM         = 0.0D0
    M           = 0.0D0

    ! Reset arrays to default values
    LC          = 0

    ! Reset types to default
    VAR_P = VARIATIONAL_PARAMETERS(0,0,0,0,0,0,0,0,1,0,0.0D0,0.0D0,0.01D0,0.02D0,40.0D0,4.0D0,0.25D0,1.02D0,0,0,100,32)
    EFT_RADIAL_CC%ORDER = -1
    EFT_RADIAL_AC%ORDER = -1
    EFT_RADIAL_AA%ORDER = -1
    EFT_RADIAL_CC%RC = 0.0D0
    EFT_RADIAL_AC%RC = 0.0D0
    EFT_RADIAL_AA%RC = 0.0D0
    IF (ALLOCATED(EFT_RADIAL_CC%FR_I)) DEALLOCATE(EFT_RADIAL_CC%FR_I)
    IF (ALLOCATED(EFT_RADIAL_AC%FR_I)) DEALLOCATE(EFT_RADIAL_AC%FR_I)
    IF (ALLOCATED(EFT_RADIAL_AA%FR_I)) DEALLOCATE(EFT_RADIAL_AA%FR_I)
    LECS = LECS_EFT_PLESS()

    WRITE(*,*) 'SCATTERING_NN_VARIATIONAL MODULE RESET COMPLETED'

  END SUBROUTINE RESET_SCATTERING_NN_VARIATIONAL

END MODULE SCATTERING_NN_VARIATIONAL