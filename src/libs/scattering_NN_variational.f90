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

  DOUBLE PRECISION :: HTM = 0, M = 0
  LOGICAL :: HTM_SET = .FALSE.
  INTEGER :: LC(NCH_MAX)
  LOGICAL :: PRINT_I = .FALSE.

  INTEGER :: LMAX=-1
  LOGICAL :: TZ_SET = .FALSE., LMAX_SET = .FALSE., PREPARE = .TRUE.

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
  PUBLIC :: NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS
  PUBLIC :: SET_VARIATIONAL_PARAMETERS
  PUBLIC :: GET_HTM
  PUBLIC :: SET_DYNAMIC
  PUBLIC :: SET_NEW_LECS
  PUBLIC :: RESET_SCATTERING_NN_VARIATIONAL
  PUBLIC :: PRINT_ALL_DATA_IN_MODULE
  PUBLIC :: FIT_CHANNEL_LOW_ENERGY
  PUBLIC :: FIT_CHANNELS_LOW_ENERGY

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

  SUBROUTINE PRINT_ALL_DATA_IN_MODULE(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER :: i
    CHARACTER(LEN=256) :: fname
    
    WRITE(fname, '(A,I0,A)') 'scattering_nn_variational_dump_unit', unit, '.txt'
    OPEN(unit, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')

    WRITE(unit,*) '==== Module SCATTERING_NN_VARIATIONAL Dump ===='

    ! Scalars
    WRITE(unit,*) 'NCH =', NCH
    WRITE(unit,*) 'NNN_MAX =', NNN_MAX
    WRITE(unit,*) 'USE_DYNAMIC =', USE_DYNAMIC
    WRITE(unit,*) 'HTM =', HTM
    WRITE(unit,*) 'M =', M
    WRITE(unit,*) 'HTM_SET =', HTM_SET
    WRITE(unit,*) 'PRINT_I =', PRINT_I
    WRITE(unit,*) 'LMAX =', LMAX
    WRITE(unit,*) 'TZ_SET =', TZ_SET
    WRITE(unit,*) 'LMAX_SET =', LMAX_SET
    WRITE(unit,*) 'PREPARE =', PREPARE
    WRITE(unit,*) 'NCHANNELS =', NCHANNELS
    WRITE(unit,*) 'CHANNELS_SET =', CHANNELS_SET
    WRITE(unit,*) 'GRID_SET =', GRID_SET
    WRITE(unit,*) 'POTENTIAL_SET =', POTENTIAL_SET
    WRITE(unit,*) 'IPOT_SET =', IPOT_SET
    WRITE(unit,*) 'NE =', NE
    WRITE(unit,*) 'ENERGIES_SET =', ENERGIES_SET
    WRITE(unit,*) 'BESSELS_SET =', BESSELS_SET
    WRITE(unit,*) 'LAGUERRE_SET =', LAGUERRE_SET
    WRITE(unit,*) 'LECS_SET =', LECS_SET
    WRITE(unit,*) 'NEW_LECS =', NEW_LECS
    WRITE(unit,*) 'CH_INDEX =', CH_INDEX

    ! Arrays and allocatables: print only sizes
    WRITE(unit,*) 'LC: size =', SIZE(LC)

    IF (ALLOCATED(CHANNELS_)) THEN
      WRITE(unit,*) 'CHANNELS_ size =', SIZE(CHANNELS_)
      DO i=1, SIZE(CHANNELS_)
        WRITE(unit,*) '  CHANNELS_(',i,') = ', GET_CHANNEL_NAME(CHANNELS_(i))
      END DO
    END IF

    WRITE(unit,*) 'CHANNEL: ', GET_CHANNEL_NAME(CHANNEL)

    IF (ALLOCATED(XX_CC)) WRITE(unit,*) 'XX_CC size =', SIZE(XX_CC)
    IF (ALLOCATED(YY_CC)) WRITE(unit,*) 'YY_CC size =', SIZE(YY_CC)
    IF (ALLOCATED(XX_AC)) WRITE(unit,*) 'XX_AC size =', SIZE(XX_AC)
    IF (ALLOCATED(XX_AA)) WRITE(unit,*) 'XX_AA size =', SIZE(XX_AA)
    IF (ALLOCATED(A_CC))  WRITE(unit,*) 'A_CC size =', SIZE(A_CC)
    IF (ALLOCATED(A_AC))  WRITE(unit,*) 'A_AC size =', SIZE(A_AC)
    IF (ALLOCATED(A_AA))  WRITE(unit,*) 'A_AA size =', SIZE(A_AA)
    IF (ALLOCATED(B_AA))  WRITE(unit,*) 'B_AA size =', SIZE(B_AA)
    IF (ALLOCATED(AJ_AA)) WRITE(unit,*) 'AJ_AA size =', SIZE(AJ_AA)
    IF (ALLOCATED(AJ_AC)) WRITE(unit,*) 'AJ_AC size =', SIZE(AJ_AC)
    IF (ALLOCATED(YYL_AC))WRITE(unit,*) 'YYL_AC size =', SIZE(YYL_AC)

    IF (ALLOCATED(V_CC))  WRITE(unit,*) 'V_CC size =', SHAPE(V_CC)
    IF (ALLOCATED(V_AC))  WRITE(unit,*) 'V_AC size =', SHAPE(V_AC)
    IF (ALLOCATED(V_AA))  WRITE(unit,*) 'V_AA size =', SHAPE(V_AA)

    IF (ALLOCATED(ENERGIES_)) WRITE(unit,*) 'ENERGIES_ size =', SIZE(ENERGIES_)
    IF (ALLOCATED(KK))       WRITE(unit,*) 'KK size =', SIZE(KK)
    IF (ALLOCATED(K2))       WRITE(unit,*) 'K2 size =', SIZE(K2)

    IF (ALLOCATED(FBES_AA))  WRITE(unit,*) 'FBES_AA size =', SHAPE(FBES_AA)
    IF (ALLOCATED(FBES_AC))  WRITE(unit,*) 'FBES_AC size =', SHAPE(FBES_AC)
    IF (ALLOCATED(GBES_AA))  WRITE(unit,*) 'GBES_AA size =', SHAPE(GBES_AA)
    IF (ALLOCATED(GBES_AC))  WRITE(unit,*) 'GBES_AC size =', SHAPE(GBES_AC)
    IF (ALLOCATED(GBES0_AA)) WRITE(unit,*) 'GBES0_AA size =', SHAPE(GBES0_AA)
    IF (ALLOCATED(GBES1_AA)) WRITE(unit,*) 'GBES1_AA size =', SHAPE(GBES1_AA)
    IF (ALLOCATED(GBES2_AA)) WRITE(unit,*) 'GBES2_AA size =', SHAPE(GBES2_AA)
    IF (ALLOCATED(HNOR_AA))  WRITE(unit,*) 'HNOR_AA size =', SHAPE(HNOR_AA)

    IF (ALLOCATED(V0_CC)) WRITE(unit,*) 'V0_CC size =', SHAPE(V0_CC)
    IF (ALLOCATED(V1_CC)) WRITE(unit,*) 'V1_CC size =', SHAPE(V1_CC)
    IF (ALLOCATED(V2_CC)) WRITE(unit,*) 'V2_CC size =', SHAPE(V2_CC)
    IF (ALLOCATED(V0_AC)) WRITE(unit,*) 'V0_AC size =', SHAPE(V0_AC)
    IF (ALLOCATED(V1_AC)) WRITE(unit,*) 'V1_AC size =', SHAPE(V1_AC)
    IF (ALLOCATED(V2_AC)) WRITE(unit,*) 'V2_AC size =', SHAPE(V2_AC)

    IF (ALLOCATED(H_MINUS_E_CC))    WRITE(unit,*) 'H_MINUS_E_CC size =', SHAPE(H_MINUS_E_CC)
    IF (ALLOCATED(H_MINUS_E_AC_R))  WRITE(unit,*) 'H_MINUS_E_AC_R size =', SHAPE(H_MINUS_E_AC_R)
    IF (ALLOCATED(H_MINUS_E_AC_I))  WRITE(unit,*) 'H_MINUS_E_AC_I size =', SHAPE(H_MINUS_E_AC_I)
    IF (ALLOCATED(H_MINUS_E_AA_RR)) WRITE(unit,*) 'H_MINUS_E_AA_RR size =', SHAPE(H_MINUS_E_AA_RR)
    IF (ALLOCATED(H_MINUS_E_AA_RI)) WRITE(unit,*) 'H_MINUS_E_AA_RI size =', SHAPE(H_MINUS_E_AA_RI)
    IF (ALLOCATED(H_MINUS_E_AA_IR)) WRITE(unit,*) 'H_MINUS_E_AA_IR size =', SHAPE(H_MINUS_E_AA_IR)
    IF (ALLOCATED(H_MINUS_E_AA_II)) WRITE(unit,*) 'H_MINUS_E_AA_II size =', SHAPE(H_MINUS_E_AA_II)

    IF (ALLOCATED(K_MINUS_E_CC))    WRITE(unit,*) 'K_MINUS_E_CC size =', SHAPE(K_MINUS_E_CC)
    IF (ALLOCATED(K_MINUS_E_AC_R))  WRITE(unit,*) 'K_MINUS_E_AC_R size =', SHAPE(K_MINUS_E_AC_R)
    IF (ALLOCATED(K_MINUS_E_AC_I))  WRITE(unit,*) 'K_MINUS_E_AC_I size =', SHAPE(K_MINUS_E_AC_I)
    IF (ALLOCATED(K_MINUS_E_AA_RR)) WRITE(unit,*) 'K_MINUS_E_AA_RR size =', SHAPE(K_MINUS_E_AA_RR)
    IF (ALLOCATED(K_MINUS_E_AA_RI)) WRITE(unit,*) 'K_MINUS_E_AA_RI size =', SHAPE(K_MINUS_E_AA_RI)
    IF (ALLOCATED(K_MINUS_E_AA_IR)) WRITE(unit,*) 'K_MINUS_E_AA_IR size =', SHAPE(K_MINUS_E_AA_IR)
    IF (ALLOCATED(K_MINUS_E_AA_II)) WRITE(unit,*) 'K_MINUS_E_AA_II size =', SHAPE(K_MINUS_E_AA_II)

    IF (ALLOCATED(VM_CC))    WRITE(unit,*) 'VM_CC size =', SHAPE(VM_CC)
    IF (ALLOCATED(VM_AC_R))  WRITE(unit,*) 'VM_AC_R size =', SHAPE(VM_AC_R)
    IF (ALLOCATED(VM_AC_I))  WRITE(unit,*) 'VM_AC_I size =', SHAPE(VM_AC_I)
    IF (ALLOCATED(VM_AA_RR)) WRITE(unit,*) 'VM_AA_RR size =', SHAPE(VM_AA_RR)
    IF (ALLOCATED(VM_AA_RI)) WRITE(unit,*) 'VM_AA_RI size =', SHAPE(VM_AA_RI)
    IF (ALLOCATED(VM_AA_IR)) WRITE(unit,*) 'VM_AA_IR size =', SHAPE(VM_AA_IR)
    IF (ALLOCATED(VM_AA_II)) WRITE(unit,*) 'VM_AA_II size =', SHAPE(VM_AA_II)

    IF (ALLOCATED(FMAT_CC))    WRITE(unit,*) 'FMAT_CC size =', SHAPE(FMAT_CC)
    IF (ALLOCATED(FMAT_AC_R))  WRITE(unit,*) 'FMAT_AC_R size =', SHAPE(FMAT_AC_R)
    IF (ALLOCATED(FMAT_AC_I))  WRITE(unit,*) 'FMAT_AC_I size =', SHAPE(FMAT_AC_I)
    IF (ALLOCATED(FMAT_AA_RR)) WRITE(unit,*) 'FMAT_AA_RR size =', SHAPE(FMAT_AA_RR)
    IF (ALLOCATED(FMAT_AA_RI)) WRITE(unit,*) 'FMAT_AA_RI size =', SHAPE(FMAT_AA_RI)
    IF (ALLOCATED(FMAT_AA_IR)) WRITE(unit,*) 'FMAT_AA_IR size =', SHAPE(FMAT_AA_IR)
    IF (ALLOCATED(FMAT_AA_II)) WRITE(unit,*) 'FMAT_AA_II size =', SHAPE(FMAT_AA_II)

    ! Dump types
    WRITE(unit,*) 'VAR_P:'
    WRITE(unit,*) '  J=', VAR_P%J, ' L=', VAR_P%L, ' S=', VAR_P%S, ' T=', VAR_P%T, ' TZ=', VAR_P%TZ
    WRITE(unit,*) '  T1Z=', VAR_P%T1Z, ' T2Z=', VAR_P%T2Z, ' IPOT=', VAR_P%IPOT, ' ILB=', VAR_P%ILB, ' LEMP=', VAR_P%LEMP
    WRITE(unit,*) '  E=', VAR_P%E, ' K=', VAR_P%K, ' HR1=', VAR_P%HR1, ' H=', VAR_P%H, ' RANGE=', VAR_P%RANGE
    WRITE(unit,*) '  GAMMA=', VAR_P%GAMMA, ' EPS=', VAR_P%EPS, ' AF=', VAR_P%AF
    WRITE(unit,*) '  NX_AA=', VAR_P%NX_AA, ' NX_AC=', VAR_P%NX_AC, ' NX_CC=', VAR_P%NX_CC, ' NNL=', VAR_P%NNL

    WRITE(unit,*) 'LECS:'
    WRITE(unit,*) '  ILB=', LECS%ILB, ' ORDER=', LECS%ORDER
    WRITE(unit,*) '  RC(0:1,0:1)'
    WRITE(unit,*) '  CLO(0:1)'
    WRITE(unit,*) '  CNLO(7)'
    WRITE(unit,*) '  CN3LO(11)'
    WRITE(unit,*) '  CIT(0:4)'

    WRITE(unit,*) 'EFT_RADIAL_CC:'
    WRITE(unit,*) '  ORDER=', EFT_RADIAL_CC%ORDER
    WRITE(unit,*) '  RC(0:1,0:1)'
    IF (ALLOCATED(EFT_RADIAL_CC%FR_I)) WRITE(unit,*) '  FR_I size =', SHAPE(EFT_RADIAL_CC%FR_I)

    WRITE(unit,*) 'EFT_RADIAL_AC:'
    WRITE(unit,*) '  ORDER=', EFT_RADIAL_AC%ORDER
    WRITE(unit,*) '  RC(0:1,0:1)'
    IF (ALLOCATED(EFT_RADIAL_AC%FR_I)) WRITE(unit,*) '  FR_I size =', SHAPE(EFT_RADIAL_AC%FR_I)

    WRITE(unit,*) 'EFT_RADIAL_AA:'
    WRITE(unit,*) '  ORDER=', EFT_RADIAL_AA%ORDER
    WRITE(unit,*) '  RC(0:1,0:1)'
    IF (ALLOCATED(EFT_RADIAL_AA%FR_I)) WRITE(unit,*) '  FR_I size =', SHAPE(EFT_RADIAL_AA%FR_I)

    CLOSE(unit)
  END SUBROUTINE PRINT_ALL_DATA_IN_MODULE


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
    IF (ALLOCATED(CHANNELS_)) THEN
      DEALLOCATE(CHANNELS_)
    ENDIF
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
    IF (J-L==1) THEN
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
   PRINT_COEFFICIENTS, PRINT_INFORMATIONS, RESET)
  IMPLICIT NONE
! INPUT PARAMETERS
  DOUBLE PRECISION, INTENT(IN) :: E
  INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP
  TYPE(PHASE_SHIFT_RESULT), INTENT(OUT) :: PHASE_SHIFT
  LOGICAL, INTENT(IN), OPTIONAL :: PRINT_COEFFICIENTS, PRINT_INFORMATIONS, RESET

  LOGICAL :: PRINT_C, FIRST_CALL = .TRUE.
! VARIABLES AND PARAMETERS FOR DGESV
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: CAR(:,:), CAI(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: CARR(:), CAII(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: XRCOEFF(:,:), XICOEFF(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: IPIV(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: C(:,:), CC(:,:), CCC(:,:)
  INTEGER :: INFO, IAK, IAB, IE
  INTEGER, SAVE :: NNN

  ! MATRICES FOR THE VARIATIONAL METHOD
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: ARI(:,:), AIR(:,:), ARR(:,:), AII(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: BD1(:,:), BD2(:,:), BD3(:,:), BD4(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: AM(:,:), AN(:,:), AMM(:,:), RMAT(:,:)

  ! COEFFICIENT RMAT2
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: RMAT2(:,:)

  ! PHASE-SHIFTS AND MIXING ANGLES
  DOUBLE PRECISION :: AMIXR, AMIXG, DELTA1, DELTA2, DELTA1G, DELTA2G
  DOUBLE PRECISION :: AMIXGS, DELTA1S, DELTA2S

! S-MATRIX
  DOUBLE COMPLEX, ALLOCATABLE, SAVE :: SMAT(:,:)

  IF (PRESENT(RESET)) THEN
    IF (.NOT. RESET) RETURN
    NNN = 0
    IF (ALLOCATED(CAR)) DEALLOCATE(CAR)
    IF (ALLOCATED(CAI)) DEALLOCATE(CAI)
    IF (ALLOCATED(CARR)) DEALLOCATE(CARR)
    IF (ALLOCATED(CAII)) DEALLOCATE(CAII)
    IF (ALLOCATED(XRCOEFF)) DEALLOCATE(XRCOEFF)
    IF (ALLOCATED(XICOEFF)) DEALLOCATE(XICOEFF)
    IF (ALLOCATED(IPIV)) DEALLOCATE(IPIV)
    IF (ALLOCATED(C)) DEALLOCATE(C)
    IF (ALLOCATED(CC)) DEALLOCATE(CC)
    IF (ALLOCATED(CCC)) DEALLOCATE(CCC)
    IF (ALLOCATED(BD1)) DEALLOCATE(BD1)
    IF (ALLOCATED(BD2)) DEALLOCATE(BD2)
    IF (ALLOCATED(BD3)) DEALLOCATE(BD3)
    IF (ALLOCATED(BD4)) DEALLOCATE(BD4)
    IF (ALLOCATED(ARI)) DEALLOCATE(ARI)
    IF (ALLOCATED(AIR)) DEALLOCATE(AIR)
    IF (ALLOCATED(ARR)) DEALLOCATE(ARR)
    IF (ALLOCATED(AII)) DEALLOCATE(AII)
    IF (ALLOCATED(AM)) DEALLOCATE(AM)
    IF (ALLOCATED(AN)) DEALLOCATE(AN)
    IF (ALLOCATED(AMM)) DEALLOCATE(AMM)
    IF (ALLOCATED(RMAT)) DEALLOCATE(RMAT)
    IF (ALLOCATED(RMAT2)) DEALLOCATE(RMAT2)
    IF (ALLOCATED(SMAT)) DEALLOCATE(SMAT)
    CALL R_SECOND_ORDER
    RETURN
  ENDIF
  
  IF (NEW_LECS .AND. USE_DYNAMIC) FIRST_CALL = .TRUE.

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
    CALL REALLOCATE_2D_9(ARI, AIR, ARR, AII, AM, AN, AMM, RMAT, RMAT2, NCH, NCH)
    IF (ALLOCATED(SMAT)) DEALLOCATE(SMAT)
    ALLOCATE(SMAT(NCH, NCH))


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
    WRITE(*,*) "Combining the CC potential"
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_CC,    LECS, VM_CC    )
    WRITE(*,*) "Combining the AC potential (real part)"
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AC_R,  LECS, VM_AC_R  )
    WRITE(*,*) "Combining the AC potential (imaginary part)"
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AC_I,  LECS, VM_AC_I  )
    VM_CC = VM_CC / HTM
    VM_AC_R = VM_AC_R / HTM
    VM_AC_I = VM_AC_I / HTM
    
    H_MINUS_E_CC   = K_MINUS_E_CC + VM_CC
    H_MINUS_E_AC_R = K_MINUS_E_AC_R + VM_AC_R
    H_MINUS_E_AC_I = K_MINUS_E_AC_I + VM_AC_I
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

  IF (USE_DYNAMIC .AND. NEW_LECS) THEN
    WRITE(*,*) "Combining the AA potential (real part)"
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_RR, LECS, VM_AA_RR )
    WRITE(*,*) "Combining the AA potential (real-imaginary part)"
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_RI, LECS, VM_AA_RI )
    WRITE(*,*) "Combining the AA potential (imaginary-real part)"
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_IR, LECS, VM_AA_IR )
    WRITE(*,*) "Combining the AA potential (imaginary part)"
    CALL COMBINE_POTENTIAL( CHANNELS_, FMAT_AA_II, LECS, VM_AA_II )
    VM_AA_RR = VM_AA_RR / HTM
    VM_AA_RI = VM_AA_RI / HTM
    VM_AA_IR = VM_AA_IR / HTM
    VM_AA_II = VM_AA_II / HTM
    WRITE(*,*) "Finished combining potentials"    
    H_MINUS_E_AA_RR = K_MINUS_E_AA_RR + VM_AA_RR
    H_MINUS_E_AA_RI = K_MINUS_E_AA_RI + VM_AA_RI
    H_MINUS_E_AA_IR = K_MINUS_E_AA_IR + VM_AA_IR
    H_MINUS_E_AA_II = K_MINUS_E_AA_II + VM_AA_II

    NEW_LECS = .FALSE.
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
  CALL CALCULATE_PHASE_SHIFTS_BLATT(RMAT2, NCH, DELTA1, DELTA2, AMIXR)
  DELTA1G = TO_DEGREES(DELTA1)
  DELTA2G = TO_DEGREES(DELTA2)
  AMIXG = TO_DEGREES(AMIXR)

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
  CALL CALCULATE_S_MATRIX(SMAT, NCH, DELTA1, DELTA2, AMIXR)

  IF (PRINT_I) THEN
    WRITE(*,*)
    WRITE(*,*)"S-MATRIX"
    WRITE(*,*) "S(1,1)=" , SMAT(1,1)
    WRITE(*,*) "S(1,2)=" , SMAT(1,2)
    WRITE(*,*) "S(2,2)=" , SMAT(2,2)
  ENDIF

  ! Calculating the phase shifts and mixing angles in the Stapp convention
  CALL CALCULATE_PHASE_SHIFTS_STAPP(RMAT2, SMAT, DELTA1S, DELTA2S, AMIXGS)

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
      IF (INFO/=0) THEN
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
      IF (NCH==2) PRINT 10, "L1:    ",  LC(2)
      PRINT 10, "S:     ",                VAR_P%S
      PRINT 10, "T:     ",                VAR_P%T
      PRINT 10, "TZ:    ",                VAR_P%TZ

      10 FORMAT(" ",A, I2)
    END SUBROUTINE PRINT_INFO

    !> \brief Compute the second order R-matrix.
    SUBROUTINE R_SECOND_ORDER()
      IMPLICIT NONE
      INTEGER :: IB
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: XRCOEFV(:,:), XICOEFV(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: BD5(:,:), BD6(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: RCI(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: RCIV(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: CD0(:,:), CD1(:,:), CD2(:,:), CD3(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: CD4(:,:), CD5(:,:), CD6(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: CD7(:,:), CD8(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: RD0(:,:), RD2(:,:), RD3(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: ASS(:,:), RMAT2_ASYM(:,:)

      IF (PRESENT(RESET)) THEN
        IF (RESET) THEN
          IF (ALLOCATED(XRCOEFV)) DEALLOCATE(XRCOEFV)
          IF (ALLOCATED(XICOEFV)) DEALLOCATE(XICOEFV)
          IF (ALLOCATED(BD5)) DEALLOCATE(BD5)
          IF (ALLOCATED(BD6)) DEALLOCATE(BD6)
          IF (ALLOCATED(RCI)) DEALLOCATE(RCI)
          IF (ALLOCATED(RCIV)) DEALLOCATE(RCIV)
          IF (ALLOCATED(CD0)) DEALLOCATE(CD0)
          IF (ALLOCATED(CD1)) DEALLOCATE(CD1)
          IF (ALLOCATED(CD2)) DEALLOCATE(CD2)
          IF (ALLOCATED(CD3)) DEALLOCATE(CD3)
          IF (ALLOCATED(CD4)) DEALLOCATE(CD4)
          IF (ALLOCATED(CD5)) DEALLOCATE(CD5)
          IF (ALLOCATED(CD6)) DEALLOCATE(CD6)
          IF (ALLOCATED(CD7)) DEALLOCATE(CD7)
          IF (ALLOCATED(CD8)) DEALLOCATE(CD8)
          IF (ALLOCATED(RD0)) DEALLOCATE(RD0)
          IF (ALLOCATED(RD2)) DEALLOCATE(RD2)
          IF (ALLOCATED(RD3)) DEALLOCATE(RD3)
          IF (ALLOCATED(ASS)) DEALLOCATE(ASS)
          IF (ALLOCATED(RMAT2_ASYM)) DEALLOCATE(RMAT2_ASYM)
          RETURN
        ENDIF
      ENDIF
      
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
      ! Transpose coefficient matrices for more efficient matrix operations
      ! XRCOEFF/XICOEFF = coefficients for regular/irregular solutions
      XRCOEFV = TRANSPOSE(XRCOEFF)  ! Transpose real coefficients matrix
      XICOEFV = TRANSPOSE(XICOEFF)  ! Transpose imaginary coefficients matrix

      ! Calculate BD5 = XRCOEFF * C (matrix product of real coefficients with core Hamiltonian matrix)
      ! Parameters: (matrix_order, transpose_A, transpose_B, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
      CALL DGEMM('N', 'N', NCH, NNN, NNN, 1.0D0, XRCOEFF, SIZE(XRCOEFF,1), C, SIZE(C,1), 0.0D0, BD5, SIZE(BD5,1))
      ! Calculate BD6 = XICOEFF * C (matrix product of imaginary coefficients with core Hamiltonian matrix)
      CALL DGEMM('N', 'N', NCH, NNN, NNN, 1.0D0, XICOEFF, SIZE(XICOEFF,1), C, SIZE(C,1), 0.0D0, BD6, SIZE(BD6,1))


      ! First block: Calculate product of R-matrix with imaginary coefficients
      ! RCI = RMAT * XICOEFF (NCH×NCH matrix multiplied by NCH×NNN matrix)
      CALL DGEMM('N', 'N', NCH, NNN, NCH, 1.0D0, RMAT, SIZE(RMAT,1), XICOEFF, SIZE(XICOEFF,1), 0.0D0, RCI, SIZE(RCI,1))

      ! Transpose RCI to RCIV for subsequent calculations
      ! This changes from (channel×basis) format to (basis×channel) format
      DO IAB=1,NCH
        DO IB=1,NNN
          RCIV(IB,IAB) = RCI(IAB,IB)
        ENDDO
      ENDDO

      ! Second block: Calculate three matrix products needed for second-order terms
      ! CD0 = BD5 * XRCOEFV (matrix product giving channel-channel interactions)
      CALL DGEMM('N', 'N', NCH, NCH, NNN, 1.0D0, BD5, SIZE(BD5,1), XRCOEFV, SIZE(XRCOEFV,1), 0.0D0, CD0, SIZE(CD0,1))

      ! CD1 = BD5 * RCIV (combines real coefficients with transposed R-matrix product)
      CALL DGEMM('N', 'N', NCH, NCH, NNN, 1.0D0, BD5, SIZE(BD5,1), RCIV, SIZE(RCIV,1), 0.0D0, CD1, SIZE(CD1,1))

      ! RD0 = BD6 * XRCOEFV (combines imaginary coefficients with real coefficients)
      CALL DGEMM('N', 'N', NCH, NCH, NNN, 1.0D0, BD6, SIZE(BD6,1), XRCOEFV, SIZE(XRCOEFV,1), 0.0D0, RD0, SIZE(RD0,1))

      ! Third block: Calculate seven matrix products with the R-matrix (RMAT)
      ! Each product represents different components of the second-order correction
      ! CD2 = BD2 * RMAT (BD2 contains core-asymptotic coupling terms)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, BD2, SIZE(BD2,1), RMAT, SIZE(RMAT,1), 0.0D0, CD2, SIZE(CD2,1))

      ! CD3 = RD0 * RMAT (combines previously calculated RD0 with R-matrix)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, RD0, SIZE(RD0,1), RMAT, SIZE(RMAT,1), 0.0D0, CD3, SIZE(CD3,1))

      ! CD4 = BD3 * RMAT (BD3 contains additional core-asymptotic terms)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, BD3, SIZE(BD3,1), RMAT, SIZE(RMAT,1), 0.0D0, CD4, SIZE(CD4,1))

      ! CD5 = ARI * RMAT (ARI contains asymptotic regular-irregular matrix elements)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, ARI, SIZE(ARI,1), RMAT, SIZE(RMAT,1), 0.0D0, CD5, SIZE(CD5,1))

      ! CD6 = AIR * RMAT (AIR contains asymptotic irregular-regular matrix elements)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, AIR, SIZE(AIR,1), RMAT, SIZE(RMAT,1), 0.0D0, CD6, SIZE(CD6,1))

      ! RD2 = BD1 * RMAT (BD1 contains another set of coupling terms)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, BD1, SIZE(BD1,1), RMAT, SIZE(RMAT,1), 0.0D0, RD2, SIZE(RD2,1))

      ! RD3 = AII * RMAT (AII contains asymptotic irregular-irregular matrix elements)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, AII, SIZE(AII,1), RMAT, SIZE(RMAT,1), 0.0D0, RD3, SIZE(RD3,1))

      ! Fourth block: Calculate two more second-order matrix products
      ! CD7 = RD2 * RMAT (Squares the effect of RD2 by multiplying it with R-matrix)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, RD2, SIZE(RD2,1), RMAT, SIZE(RMAT,1), 0.0D0, CD7, SIZE(CD7,1))

      ! CD8 = RD3 * RMAT (Squares the effect of RD3 by multiplying it with R-matrix)
      CALL DGEMM('N', 'N', NCH, NCH, NCH, 1.0D0, RD3, SIZE(RD3,1), RMAT, SIZE(RMAT,1), 0.0D0, CD8, SIZE(CD8,1))

      ! Print a blank line before R-matrix output if debug printing is enabled
      IF (PRINT_I) WRITE(*,*)

      ! Combine all matrix products to form the second-order correction to the R-matrix
      DO IAB=1,NCH
        DO IAK=1,NCH
          ! Calculate the second-order correction ASS by summing all contributing terms
          ASS(IAB,IAK)=CD0(IAB,IAK)+2*BD4(IAB,IAK)+CD4(IAB,IAK)    &  ! Terms involving real solutions
                      +ARR(IAB,IAK)+CD5(IAB,IAK)+CD2(IAB,IAK)    &    ! Mixed coupling terms
                      +CD7(IAB,IAK)+CD6(IAB,IAK)    &                 ! Second-order terms
                      +CD8(IAB,IAK)                                   ! Final contribution

          ! Add correction to the original R-matrix to get the asymmetric second-order R-matrix
          RMAT2_ASYM(IAB,IAK) = RMAT(IAB,IAK) + ASS(IAB,IAK)
          
          ! Print R-matrix components if debug mode is enabled
          IF (PRINT_I) WRITE(*,*)"COEFF RMAT2",RMAT2_ASYM(IAB,IAK), RMAT(IAB,IAK), ASS(IAB,IAK)
        ENDDO
      ENDDO

      ! Symmetrize the R-matrix (make it symmetric) for physical consistency
      ! The negative factor relates to the convention used in the Kohn variational principle
      RMAT2 = ZERO  ! Initialize the symmetric R-matrix to zero
      DO IAB=1,NCH
        DO IAK=1,NCH
          ! Average the asymmetric R-matrix with its transpose and multiply by -0.5
          RMAT2(IAB,IAK) =-0.5D0*( RMAT2_ASYM(IAB,IAK) + RMAT2_ASYM(IAK,IAB) )
        ENDDO
      ENDDO
    END SUBROUTINE R_SECOND_ORDER

    PURE FUNCTION TO_DEGREES(ANGLE_RADIANS) RESULT(ANGLE_DEGREES)
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: PI_C = 4.D0*DATAN(1.D0)
      DOUBLE PRECISION, INTENT(IN) :: ANGLE_RADIANS
      DOUBLE PRECISION :: ANGLE_DEGREES
      ANGLE_DEGREES = (ANGLE_RADIANS*180.D0)/PI_C
    END FUNCTION TO_DEGREES

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
    SUBROUTINE CALCULATE_PHASE_SHIFTS_BLATT(R, DIM, DELTA_1BB, DELTA_2BB, AMIX_BB)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: DIM
      DOUBLE PRECISION, INTENT(IN) :: R(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: DELTA_1BB, DELTA_2BB, AMIX_BB
      IF (DIM == 1) THEN
        DELTA_1BB = ATAN((R(1,1)))
        DELTA_2BB = 0.D0
        AMIX_BB  = 0.D0
        RETURN
      ENDIF

      IF (DIM > 2) THEN
        PRINT *, "Error: R-matrix dimension is greater than 2 in CALCULATE_PHASE_SHIFTS_BLATT"
        STOP
      ENDIF
      IF (R(1,2) /= ZERO .AND. R(1,1) == R(2,2)) THEN
        PRINT *, "Error: R is singular in CALCULATE_PHASE_SHIFTS_BLATT"
        STOP
      ENDIF
      AMIX_BB=0.5D0*ATAN(2.*R(1,2)/(R(1,1)-R(2,2)))

      DELTA_1BB=ATAN((COS(AMIX_BB)*COS(AMIX_BB)*R(1,1)  &
                   +SIN(AMIX_BB)*SIN(AMIX_BB)*R(2,2)  &
                   +2*COS(AMIX_BB)*SIN(AMIX_BB)*R(1,2)))

      DELTA_2BB=ATAN((SIN(AMIX_BB)*SIN(AMIX_BB)*R(1,1)  &
                   +COS(AMIX_BB)*COS(AMIX_BB)*R(2,2)  &
                   -2*COS(AMIX_BB)*SIN(AMIX_BB)*R(1,2)))

    END SUBROUTINE CALCULATE_PHASE_SHIFTS_BLATT

    !> \brief Calculate the S-matrix from the R-matrix.
    SUBROUTINE CALCULATE_S_MATRIX(SM, DIM, DELTA_BB1, DELTA_BB2, MIXING_ANGLE_BB)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: DELTA_BB1, DELTA_BB2, MIXING_ANGLE_BB
      INTEGER, INTENT(IN) :: DIM
      DOUBLE COMPLEX, INTENT(OUT) :: SM(:,:)

      DOUBLE PRECISION :: COS1, SIN1
      DOUBLE COMPLEX :: SM1, SM2
      SM1  = CDEXP(2.D0*IM*DELTA_BB1)
      SM2  = CDEXP(2.D0*IM*DELTA_BB2)
      COS1 = COS(MIXING_ANGLE_BB)
      SIN1 = SIN(MIXING_ANGLE_BB)
      SM(1,1) = COS1*COS1*SM1 + SIN1*SIN1*SM2
      IF ( DIM == 1 ) RETURN
      SM(2,2) = COS1*COS1*SM2 + SIN1*SIN1*SM1
      SM(1,2) = COS1*SIN1*(SM1 - SM2)
      SM(2,1) = SM(1,2)
    END SUBROUTINE CALCULATE_S_MATRIX

    !> \brief Calculate phase shifts and mixing angle in the Stapp convention.
    SUBROUTINE CALCULATE_PHASE_SHIFTS_STAPP(R, SM, DELTA1_S, DELTA2_S, AMIXG_S)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: R(:,:)
      DOUBLE COMPLEX, INTENT(IN) :: SM(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: DELTA1_S, DELTA2_S, AMIXG_S

      DOUBLE PRECISION, PARAMETER :: TOL = 1.D-14
      DOUBLE COMPLEX :: SM1, SM2
      DOUBLE PRECISION :: SI2E, CI2E, CX, SX
      IF (NCH == 1) THEN
        DELTA1_S = ATAN((R(1,1)))
        DELTA1_S = (DELTA1_S*180.D0)/PI
        DELTA2_S = 0.D0
        AMIXG_S  = 0.D0
        RETURN
      ENDIF
      SM1  = REAL(SM(1,1)*SM(2,2)-SM(1,2)*SM(1,2))
      SI2E =-REAL(SM(1,2)*SM(1,2)/SM1)
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
      SM1 = CDSQRT(SM(1,1)/CI2E)
      SM2 = CDSQRT(SM(2,2)/CI2E)
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
      DELTA1_S = DACOS(CX)*180.D0/PI
      IF(SX.LT.ZERO) DELTA1_S = -DELTA1_S

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
      DELTA2_S = DACOS(CX)*180.D0/PI
      IF(SX.LT.ZERO) DELTA2_S = -DELTA2_S

      IF (ABS(SI2E).GT.ONE) THEN
        IF (ABS(SI2E-1).LT.TOL) THEN
          SI2E = SI2E/DABS(SI2E)
        ELSE
          PRINT *, "Error: SI2E is out of bounds in CALCULATE_PHASE_SHIFTS_STAPP", SI2E, " SI2E - 1", SI2E-1
          STOP
        ENDIF
      ENDIF
      AMIXG_S = (0.5D0*DASIN(SI2E))*180.D0/PI
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
    INTEGER :: LIK, IAB, IAK, IL, IR, IB, IK, LR
    INTEGER, ALLOCATABLE :: COMMON_INDEX(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: KIN_MATRIX(:,:), POT_MATRIX(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: HCC(:,:,:), ENCC(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: INTEGRAND(:)
    INTEGER :: SL, SR, TL, TR, T, S

    IF (.NOT.ENERGIES_SET .OR. .NOT.GRID_SET .OR. .NOT.LAGUERRE_SET) THEN
      PRINT *, "Error: Energies, grid, Bessels or Laguerre polynomials not set"
      STOP
    ENDIF
    
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

    DO I = 1, NNN_MAX
      ENCC(:,I,I) = ENERGIES_
    ENDDO

    WRITE(*,*) "Preparing core-core matrix elements for channels: ", NCHANNELS
    IF (USE_DYNAMIC) THEN
      WRITE(*,*) "Using dynamic potential for core-core matrix elements"
    ELSE
      WRITE(*,*) "Using static potential for core-core matrix elements"
    ENDIF
    
    DO ICH = 1, NCHANNELS
      NEQ = GET_CHANNEL_NCH(CHANNELS_(ICH))

      IF (ALLOCATED(COMMON_INDEX)) DEALLOCATE(COMMON_INDEX)
      ALLOCATE(COMMON_INDEX(NEQ, NNL))
      COMMON_INDEX = ZERO
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
        SL = GET_CHANNEL_S(CHANNELS_(ICH), IAB)
        SR = GET_CHANNEL_S(CHANNELS_(ICH), IAK)
        IF (SL/=SR) CYCLE ! Only diagonal elements are calculated
        TL = GET_CHANNEL_T(CHANNELS_(ICH), IAB)
        TR = GET_CHANNEL_T(CHANNELS_(ICH), IAK)
        IF (TL/=TR) CYCLE ! Only diagonal elements are calculated
        T = TL
        S = SL

        LR = GET_CHANNEL_L(CHANNELS_(ICH), IAK)
        LIK= LR*(LR+1)

        DO IL=1, NNL
        DO IR=1, NNL

          IB=COMMON_INDEX(IAB,IL)
          IK=COMMON_INDEX(IAK,IR)

    ! CALCULATING KINETIC ENERGY
          KIN_MATRIX(IB,IK) = ZERO
          IF(IAB==IAK)THEN
            INTEGRAND = V0_CC(IL,:)*( V2_CC(IR,:) + 2.D0*V1_CC(IR,:)/XX_CC &
                  -LIK*V0_CC(IR,:)/XX_CC**2 )
            SUM = ZERO
            DO I=1,NX
              SUM = SUM + XX_CC(I)**2*INTEGRAND(I)*A_CC(I)
            ENDDO
            KIN_MATRIX(IB,IK)=-HTM*SUM/GAMMA
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
            DO IPOT=0, ORDER_TO_NMAX(EFT_RADIAL_CC%ORDER)
              SUM = ZERO
              INTEGRAND = V0_CC(IL,:)*V0_CC(IR,:)*EFT_RADIAL_CC%FR_I(S,T,IPOT,:)
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

  END SUBROUTINE PREPARE_CORE_CORE_MATRIX_ELEMENTS

  !> \ingroup scattering_nn_variational_mod
  !> \brief Prepare asymptotic-core matrix elements for the variational calculation.
  ! filepath: /home/alessandro/Dropbox/Variazionale_mine/src/libs/scattering_NN_variational.f90
  SUBROUTINE PREPARE_ASYMPTOTIC_CORE_MATRIX_ELEMENTS()
    USE LAGUERRE_POLYNOMIAL_MOD
    USE INTEGRATION_MOD
    USE OMP_LIB
    IMPLICIT NONE

    DOUBLE PRECISION :: H5, HR ! Step size in r
    INTEGER :: NX, NEQ, NNN
    INTEGER :: IE, ICH, II, JJ, I, LIK
    INTEGER :: IAB, IAK, IL, IB, LL, IPOT
    INTEGER :: SL, SR, TL, TR, S, T
    INTEGER :: NUM_THREADS
    LOGICAL :: PARALLEL_ENABLED = .TRUE.

    DOUBLE PRECISION, ALLOCATABLE :: INTEGRAND(:)
    INTEGER, ALLOCATABLE :: COMMON_INDEX(:,:)

    INTEGER :: NCOMB
    INTEGER, ALLOCATABLE :: LCOMBINATIONS(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: KIN_MIN_E(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: POT_AC_R(:,:,:,:,:,:), POT_AC_I(:,:,:,:,:,:)

    ! Check if prerequisites are set
    IF (.NOT.ENERGIES_SET .OR. .NOT.GRID_SET .OR. .NOT.BESSELS_SET .OR. .NOT.LAGUERRE_SET) THEN
      PRINT *, "Error: Energies, grid, Bessels or Laguerre polynomials not set"
      STOP
    ENDIF

    WRITE(*,*) "Preparing asymptotic-core matrix elements for channels: ", NCHANNELS
    IF (USE_DYNAMIC) THEN
      WRITE(*,*) "Using dynamic potential for asymptotic-core matrix elements"
    ELSE
      WRITE(*,*) "Using static potential for asymptotic-core matrix elements"
    ENDIF

    ! Allocate output arrays
    CALL REALLOCATE_4D_1(H_MINUS_E_AC_R, NCHANNELS, NE, NNN_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AC_I, NCHANNELS, NE, NNN_MAX, NCH_MAX)

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

    ! Initialize integration array with padding for index 1
    CALL REALLOCATE_1D_1(INTEGRAND, NX+1)

    ! Calculate channel combinations
    CALL L_COMBINATIONS(CHANNELS_, LCOMBINATIONS)
    NCOMB = SIZE(LCOMBINATIONS, DIM=1)
    
    ! Initialize kinetic energy matrix and potential matrices if dynamic
    ALLOCATE(KIN_MIN_E(NE, VAR_P%NNL, 0:LMAX))
    KIN_MIN_E = ZERO
    
    IF (USE_DYNAMIC) THEN
      ALLOCATE(POT_AC_R(0:ORDER_TO_NMAX(EFT_RADIAL_AC%ORDER), NE, 0:1, 0:1, VAR_P%NNL, 0:LMAX))
      ALLOCATE(POT_AC_I(0:ORDER_TO_NMAX(EFT_RADIAL_AC%ORDER), NE, 0:1, 0:1, VAR_P%NNL, 0:LMAX))
      POT_AC_R = ZERO
      POT_AC_I = ZERO
    ENDIF

    ! Set OpenMP parameters - determine number of threads to use
    NUM_THREADS = 1
    !$OMP PARALLEL
    !$OMP MASTER
    IF (PARALLEL_ENABLED) NUM_THREADS = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
    !$OMP END PARALLEL
    
    WRITE(*,*) "Computing matrix elements using", NUM_THREADS, "threads"

    ! Compute kinetic energy matrix elements and potential terms - this is the computationally intensive part
    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(IE, IL, LL, LIK, INTEGRAND, S, T, IPOT) &
    !$OMP SHARED(KIN_MIN_E, POT_AC_R, POT_AC_I, NE, LMAX, NX, H5) &
    !$OMP SCHEDULE(dynamic) IF(PARALLEL_ENABLED)
    DO IE = 1, NE
      DO IL = 1, VAR_P%NNL
        DO LL = 0, LMAX
          LIK = LL*(LL+1)
          
          ! Compute kinetic energy terms - extracted to a helper routine
          CALL COMPUTE_KIN_MIN_E(IE, IL, LL, LIK, NX, H5, KIN_MIN_E(IE,IL,LL))
          
          ! Compute potential terms if using dynamic potential
          IF (USE_DYNAMIC) THEN
            CALL COMPUTE_POTENTIAL_TERMS(IE, IL, LL, NX, H5, POT_AC_R, POT_AC_I)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Evaluate the matrix elements for each channel and energy
    DO IE = 1, NE
      DO ICH = 1, NCHANNELS
        NEQ = GET_CHANNEL_NCH(CHANNELS_(ICH))
        NNN = NEQ * VAR_P%NNL

        ! Create a mapping from (channel, basis function) to linear index
        IF (ALLOCATED(COMMON_INDEX)) DEALLOCATE(COMMON_INDEX)
        ALLOCATE(COMMON_INDEX(NEQ, VAR_P%NNL))
        COMMON_INDEX = 0
        II = 0
        DO I=1, NEQ
          DO JJ=1, VAR_P%NNL
            II = II + 1
            COMMON_INDEX(I, JJ) = II
          ENDDO
        ENDDO

        !$OMP PARALLEL DO COLLAPSE(2) &
        !$OMP PRIVATE(IAB, IAK, IL, IB, LL, SL, SR, TL, TR, S, T, IPOT, INTEGRAND) &
        !$OMP SHARED(H_MINUS_E_AC_R, H_MINUS_E_AC_I, K_MINUS_E_AC_I, FMAT_AC_R, FMAT_AC_I) &
        !$OMP SCHEDULE(dynamic) IF(PARALLEL_ENABLED)
        DO IAB = 1, NEQ
          DO IAK = 1, NEQ
            ! Get quantum numbers for current pair of channels
            SL = GET_CHANNEL_S(CHANNELS_(ICH), IAK)
            SR = GET_CHANNEL_S(CHANNELS_(ICH), IAB)
            IF (SL/=SR) CYCLE ! Only diagonal elements are calculated
            
            TL = GET_CHANNEL_T(CHANNELS_(ICH), IAK)
            TR = GET_CHANNEL_T(CHANNELS_(ICH), IAB)
            IF (TL/=TR) CYCLE ! Only diagonal elements are calculated
            
            S = SL
            T = TL

            LL = GET_CHANNEL_L(CHANNELS_(ICH), IAK)

            DO IL = 1, VAR_P%NNL
              IB = COMMON_INDEX(IAB, IL)

              ! Handle dynamic and static potential differently
              IF (USE_DYNAMIC) THEN
                ! Copy precomputed matrix elements to output arrays
                DO IPOT = 0, ORDER_TO_NMAX(EFT_RADIAL_AC%ORDER)
                  FMAT_AC_R(IPOT,ICH,IE,IB,IAK) = POT_AC_R(IPOT,IE,S,T,IL,LL)
                  FMAT_AC_I(IPOT,ICH,IE,IB,IAK) = POT_AC_I(IPOT,IE,S,T,IL,LL)
                ENDDO
                
                ! Set kinetic matrix elements for matching diagonal IAK==IAB elements
                IF (IAK == IAB) THEN
                  K_MINUS_E_AC_I(ICH, IE, IB, IAK) = KIN_MIN_E(IE,IL,LL) / HTM
                ENDIF
              ELSE
                ! Handle static potential case
                CALL COMPUTE_STATIC_POTENTIAL_MATRIX_ELEMENTS(ICH, IE, IB, IAK, IL, LL, NX, H5)
              ENDIF
            ENDDO ! IL
          ENDDO ! IAK
        ENDDO ! IAB
        !$OMP END PARALLEL DO
      ENDDO ! ICH
    ENDDO ! IE
    
    ! Clean up temporary arrays
    IF (ALLOCATED(KIN_MIN_E)) DEALLOCATE(KIN_MIN_E)
    IF (ALLOCATED(POT_AC_R)) DEALLOCATE(POT_AC_R)
    IF (ALLOCATED(POT_AC_I)) DEALLOCATE(POT_AC_I)
    IF (ALLOCATED(LCOMBINATIONS)) DEALLOCATE(LCOMBINATIONS)
    
  CONTAINS

    !> \brief Compute kinetic energy minus potential energy matrix elements
    SUBROUTINE COMPUTE_KIN_MIN_E(IE_LOCAL, IL_LOCAL, LL_LOCAL, LIK_LOCAL, NX_LOCAL, H5_LOCAL, RESULT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IE_LOCAL, IL_LOCAL, LL_LOCAL, NX_LOCAL
      INTEGER, INTENT(IN) :: LIK_LOCAL
      DOUBLE PRECISION, INTENT(IN) :: H5_LOCAL
      DOUBLE PRECISION, INTENT(OUT) :: RESULT
      DOUBLE PRECISION :: AXX1_LOCAL, AKE1_LOCAL
      DOUBLE PRECISION, ALLOCATABLE :: INTEGRAND_LOCAL(:)
      
      ALLOCATE(INTEGRAND_LOCAL(NX_LOCAL+1))
      
      ! NORMALIZATIONS CORE-IRREGULAR
      INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
      INTEGRAND_LOCAL(2:) = AJ_AC*V0_AC(IL_LOCAL,:)*GBES_AC(IE_LOCAL,LL_LOCAL,:)
      AXX1_LOCAL = B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)

      ! KINETIC ENERGY CORE-IRREGULAR
      INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
      INTEGRAND_LOCAL(2:) = AJ_AC*GBES_AC(IE_LOCAL,LL_LOCAL,:)*( &
                      V2_AC(IL_LOCAL,:) + 2.D0*V1_AC(IL_LOCAL,:)/XX_AC &
                      - LIK_LOCAL*V0_AC(IL_LOCAL,:)/XX_AC**2)
      AKE1_LOCAL = -HTM * B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)
      
      ! KINETIC ENERGY - ENERGY CORE-IRREGULAR
      RESULT = AKE1_LOCAL - ENERGIES_(IE_LOCAL) * AXX1_LOCAL
      
      DEALLOCATE(INTEGRAND_LOCAL)
    END SUBROUTINE COMPUTE_KIN_MIN_E

    !> \brief Compute potential terms for dynamic case
    SUBROUTINE COMPUTE_POTENTIAL_TERMS(IE_LOCAL, IL_LOCAL, LL_LOCAL, NX_LOCAL, H5_LOCAL, &
                                      POT_R_LOCAL, POT_I_LOCAL)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IE_LOCAL, IL_LOCAL, LL_LOCAL, NX_LOCAL
      DOUBLE PRECISION, INTENT(IN) :: H5_LOCAL
      DOUBLE PRECISION, INTENT(INOUT) :: POT_R_LOCAL(0:,:,0:,0:,:,0:), POT_I_LOCAL(0:,:,0:,0:,:,0:)
      INTEGER :: T_LOCAL, IPOT_LOCAL
      DOUBLE PRECISION, ALLOCATABLE :: INTEGRAND_LOCAL(:)
      
      ALLOCATE(INTEGRAND_LOCAL(NX_LOCAL+1))
      
      DO T_LOCAL = 0, 1
        DO IPOT_LOCAL = 0, ORDER_TO_NMAX(EFT_RADIAL_AC%ORDER)
          ! Regular component
          INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
          INTEGRAND_LOCAL(2:) = AJ_AC*V0_AC(IL_LOCAL,:)*FBES_AC(IE_LOCAL,LL_LOCAL,:)* &
                                EFT_RADIAL_AC%FR_I(0,T_LOCAL,IPOT_LOCAL,:)
          POT_R_LOCAL(IPOT_LOCAL,IE_LOCAL,0,T_LOCAL,IL_LOCAL,LL_LOCAL) = &
                      B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)

          ! Irregular component
          INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
          INTEGRAND_LOCAL(2:) = AJ_AC*V0_AC(IL_LOCAL,:)*GBES_AC(IE_LOCAL,LL_LOCAL,:)* &
                                EFT_RADIAL_AC%FR_I(0,T_LOCAL,IPOT_LOCAL,:)
          POT_I_LOCAL(IPOT_LOCAL,IE_LOCAL,0,T_LOCAL,IL_LOCAL,LL_LOCAL) = &
                      B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)
          
          ! Handle spin=1 case, check if cutoffs differ
          IF (DABS(LECS%RC(0,T_LOCAL)-LECS%RC(1,T_LOCAL)) > 1.D-10) THEN
            ! Regular component S=1
            INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
            INTEGRAND_LOCAL(2:) = AJ_AC*V0_AC(IL_LOCAL,:)*FBES_AC(IE_LOCAL,LL_LOCAL,:)* &
                                  EFT_RADIAL_AC%FR_I(1,T_LOCAL,IPOT_LOCAL,:)
            POT_R_LOCAL(IPOT_LOCAL,IE_LOCAL,1,T_LOCAL,IL_LOCAL,LL_LOCAL) = &
                        B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)

            ! Irregular component S=1
            INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
            INTEGRAND_LOCAL(2:) = AJ_AC*V0_AC(IL_LOCAL,:)*GBES_AC(IE_LOCAL,LL_LOCAL,:)* &
                                  EFT_RADIAL_AC%FR_I(1,T_LOCAL,IPOT_LOCAL,:)
            POT_I_LOCAL(IPOT_LOCAL,IE_LOCAL,1,T_LOCAL,IL_LOCAL,LL_LOCAL) = &
                        B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)
          ELSE
            ! Use same values as S=0 case if cutoffs are identical
            POT_R_LOCAL(IPOT_LOCAL,IE_LOCAL,1,T_LOCAL,IL_LOCAL,LL_LOCAL) = &
                        POT_R_LOCAL(IPOT_LOCAL,IE_LOCAL,0,T_LOCAL,IL_LOCAL,LL_LOCAL)
            POT_I_LOCAL(IPOT_LOCAL,IE_LOCAL,1,T_LOCAL,IL_LOCAL,LL_LOCAL) = &
                        POT_I_LOCAL(IPOT_LOCAL,IE_LOCAL,0,T_LOCAL,IL_LOCAL,LL_LOCAL)
          ENDIF
        ENDDO ! IPOT_LOCAL
      ENDDO ! T_LOCAL
      
      DEALLOCATE(INTEGRAND_LOCAL)
    END SUBROUTINE COMPUTE_POTENTIAL_TERMS

    !> \brief Compute static potential matrix elements
    SUBROUTINE COMPUTE_STATIC_POTENTIAL_MATRIX_ELEMENTS(ICH_LOCAL, IE_LOCAL, IB_LOCAL, IAK_LOCAL, &
                                                      IL_LOCAL, LL_LOCAL, NX_LOCAL, H5_LOCAL)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ICH_LOCAL, IE_LOCAL, IB_LOCAL, IAK_LOCAL, IL_LOCAL, LL_LOCAL, NX_LOCAL
      DOUBLE PRECISION, INTENT(IN) :: H5_LOCAL
      DOUBLE PRECISION :: APE_LOCAL, APE1_LOCAL
      DOUBLE PRECISION, ALLOCATABLE :: INTEGRAND_LOCAL(:)
      INTEGER :: IAB_LOCAL

      ALLOCATE(INTEGRAND_LOCAL(NX_LOCAL+1))
      
      ! Get IAB value from IB
      IAB_LOCAL = (IB_LOCAL-1) / VAR_P%NNL + 1
      
      ! Evaluate the potential energy core-regular (ape_local)
      INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
      INTEGRAND_LOCAL(2:) = AJ_AC*V0_AC(IL_LOCAL,:)*FBES_AC(IE_LOCAL,LL_LOCAL,:)* &
                          V_AC(ICH_LOCAL,:,IAB_LOCAL,IAK_LOCAL)
      APE_LOCAL = B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)

      ! Evaluate the potential energy core-irregular (ape1_local)
      INTEGRAND_LOCAL(1) = ZERO  ! Pad first element
      INTEGRAND_LOCAL(2:) = AJ_AC*V0_AC(IL_LOCAL,:)*GBES_AC(IE_LOCAL,LL_LOCAL,:)* &
                          V_AC(ICH_LOCAL,:,IAB_LOCAL,IAK_LOCAL)
      APE1_LOCAL = B5_SINGLE(NX_LOCAL, H5_LOCAL, INTEGRAND_LOCAL, 1)

      ! Store results in output arrays
      H_MINUS_E_AC_R(ICH_LOCAL, IE_LOCAL, IB_LOCAL, IAK_LOCAL) = APE_LOCAL / HTM
      
      IF (IAK_LOCAL == IAB_LOCAL) THEN
        H_MINUS_E_AC_I(ICH_LOCAL, IE_LOCAL, IB_LOCAL, IAK_LOCAL) = &
                        (KIN_MIN_E(IE_LOCAL,IL_LOCAL,LL_LOCAL) + APE1_LOCAL) / HTM
      ELSE
        H_MINUS_E_AC_I(ICH_LOCAL, IE_LOCAL, IB_LOCAL, IAK_LOCAL) = APE1_LOCAL / HTM
      ENDIF
      
      DEALLOCATE(INTEGRAND_LOCAL)
    END SUBROUTINE COMPUTE_STATIC_POTENTIAL_MATRIX_ELEMENTS

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

    INTEGER :: SL, SR, TL, TR, S, T
    INTEGER :: NX
    DOUBLE PRECISION, ALLOCATABLE :: INTEGRAND(:)

    INTEGER :: NCOMB
    INTEGER, ALLOCATABLE :: LCOMBINATIONS(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: KIN_MIN_E_RI(:,:,:), KIN_MIN_E_II(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: POT_AA_RR(:,:,:,:,:,:), POT_AA_IR(:,:,:,:,:,:), &
                                     POT_AA_RI(:,:,:,:,:,:), POT_AA_II(:,:,:,:,:,:)

    IF (.NOT.ENERGIES_SET .OR. .NOT.GRID_SET .OR. .NOT.BESSELS_SET) THEN
      PRINT *, "Error: Energies, grid or Bessel functions not set in PREPARE_ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS"
      STOP
    ENDIF

    CALL REALLOCATE_4D_1(H_MINUS_E_AA_RR, NCHANNELS, NE, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AA_IR, NCHANNELS, NE, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AA_RI, NCHANNELS, NE, NCH_MAX, NCH_MAX)
    CALL REALLOCATE_4D_1(H_MINUS_E_AA_II, NCHANNELS, NE, NCH_MAX, NCH_MAX)

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
    CALL REALLOCATE_1D_1(INTEGRAND, NX+1)


    WRITE(*,*) "Preparing asymptotic-asymptotic matrix elements for channels: ", NCHANNELS
    IF (USE_DYNAMIC) THEN
      WRITE(*,*) "Using dynamic potential for asymptotic-asymptotic matrix elements"
    ELSE
      WRITE(*,*) "Using static potential for asymptotic-asymptotic matrix elements"
    ENDIF

    ! MATRIX ELEMENTS DEPENDING ONLY ON (LL, LR)
    CALL L_COMBINATIONS(CHANNELS_, LCOMBINATIONS)
    NCOMB = SIZE(LCOMBINATIONS, 1)
    ALLOCATE(KIN_MIN_E_RI(NE, 0:LMAX, 0:LMAX))
    ALLOCATE(KIN_MIN_E_II(NE, 0:LMAX, 0:LMAX))
    KIN_MIN_E_RI = ZERO
    KIN_MIN_E_II = ZERO
    IF (USE_DYNAMIC) THEN
      ALLOCATE(POT_AA_RR(0:7, NE, 0:1, 0:1, 0:LMAX, 0:LMAX))
      ALLOCATE(POT_AA_IR(0:7, NE, 0:1, 0:1, 0:LMAX, 0:LMAX))
      ALLOCATE(POT_AA_RI(0:7, NE, 0:1, 0:1, 0:LMAX, 0:LMAX))
      ALLOCATE(POT_AA_II(0:7, NE, 0:1, 0:1, 0:LMAX, 0:LMAX))
      POT_AA_RR = ZERO
      POT_AA_IR = ZERO
      POT_AA_RI = ZERO
      POT_AA_II = ZERO
    ENDIF

    DO ICH =1, NCOMB
      LL = LCOMBINATIONS(ICH,1)
      LR = LCOMBINATIONS(ICH,2)
      DO IE=1, NE
        IF (LL==LR) THEN ! Only diagonal elements are calculated
          ! NORM REGULAR-IRREGULAR (AXX) AND IRREGULAR-IRREGULAR(AXX3)
          INTEGRAND(2:) = AJ_AA*FBES_AA(IE, LL,:)*GBES_AA(IE, LR,:)
          AXX=  B5_SINGLE(NX,H5,INTEGRAND,1)
          
          INTEGRAND(2:) = AJ_AA*GBES_AA(IE, LL,:)*GBES_AA(IE, LR,:)
          AXX3= B5_SINGLE(NX,H5,INTEGRAND,1)

          ! KINETIC ENERGY
          INTEGRAND(2:) = AJ_AA*FBES_AA(IE,LL,:)*HNOR_AA(IE,LR,:)*(GBES2_AA(IE,LR,:)+GBES1_AA(IE,LR,:)+GBES0_AA(IE,LR,:))
          AKE=  HTM * B5_SINGLE(NX,H5,INTEGRAND,1)

          INTEGRAND(2:) = AJ_AA*GBES_AA(IE,LL,:)*HNOR_AA(IE,LR,:)*(GBES2_AA(IE,LR,:)+GBES1_AA(IE,LR,:)+GBES0_AA(IE,LR,:))
          AKE3= HTM * B5_SINGLE(NX,H5,INTEGRAND,1)

          KIN_MIN_E_RI(IE,LL,LR) = AKE  - ENERGIES_(IE) * AXX
          KIN_MIN_E_II(IE,LL,LR) = AKE3 - ENERGIES_(IE) * AXX3
        ENDIF

        ! POTENTIAL RADIAL ENERGY (DYNAMIC CASE)
        IF (USE_DYNAMIC) THEN
          DO S = 0, 1
          DO T = 0, 1
          DO IPOT = 0, ORDER_TO_NMAX(EFT_RADIAL_AA%ORDER)
            INTEGRAND(2:) = AJ_AA*FBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(S,T,IPOT,:)
            POT_AA_RI(IPOT,IE,S,T,LL,LR) = B5_SINGLE(NX,H5,INTEGRAND,1)

            INTEGRAND(2:) = AJ_AA*GBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(S,T,IPOT,:)
            POT_AA_IR(IPOT,IE,S,T,LL,LR) = B5_SINGLE(NX,H5,INTEGRAND,1)

            INTEGRAND(2:) = AJ_AA*FBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(S,T,IPOT,:)
            POT_AA_RR(IPOT,IE,S,T,LL,LR) = B5_SINGLE(NX,H5,INTEGRAND,1)

            INTEGRAND(2:) = AJ_AA*GBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*EFT_RADIAL_AA%FR_I(S,T,IPOT,:)
            POT_AA_II(IPOT,IE,S,T,LL,LR) = B5_SINGLE(NX,H5,INTEGRAND,1)
          ENDDO
          ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    ! END PREPARATION OF MATRIX ELEMENTS DEPENDING ON (LL, LR)

    ! EVALUATION OF HAMILTONIAN AND KINETIC MINUS ENERGY MATRIX ELEMENTS
    DO ICH = 1, NCHANNELS
      DO IE = 1, NE
        NEQ = GET_CHANNEL_NCH(CHANNELS_(ICH))
        DO IAB=1,NEQ
        DO IAK=1,NEQ

          SL = GET_CHANNEL_S(CHANNELS_(ICH), IAB)
          SR = GET_CHANNEL_S(CHANNELS_(ICH), IAK)
          IF (SL/=SR) CYCLE ! Only diagonal elements are calculated
          TL = GET_CHANNEL_T(CHANNELS_(ICH), IAB)
          TR = GET_CHANNEL_T(CHANNELS_(ICH), IAK)
          IF (TL/=TR) CYCLE ! Only diagonal elements are calculated
          S = SL
          T = TL

          LL = GET_CHANNEL_L(CHANNELS_(ICH),IAB)
          LR = GET_CHANNEL_L(CHANNELS_(ICH),IAK)
      

      ! SI CALCOLA ENERGIA POTENZIALE DEL CASO REGOLARE-IRREGOLARE (APE),
      !      IRREGOLARE-REGOLARE (APE1), REGOLARE-REGOLARE (APE2), IRREGOLARE-IRREGOLARE (APE3)
          IF (.NOT.USE_DYNAMIC) THEN
            INTEGRAND(2:) = AJ_AA*FBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)
            APE=  B5_SINGLE(NX,H5,INTEGRAND,1)

            INTEGRAND(2:)  = AJ_AA*GBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)
            APE1= B5_SINGLE(NX,H5,INTEGRAND,1)

            INTEGRAND(2:) = AJ_AA*FBES_AA(IE,LL,:)*FBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)
            APE2= B5_SINGLE(NX,H5,INTEGRAND,1)

            INTEGRAND(2:) = AJ_AA*GBES_AA(IE,LL,:)*GBES_AA(IE,LR,1:NX)*V_AA(ICH,:,IAB,IAK)
            APE3= B5_SINGLE(NX,H5,INTEGRAND,1)
          ELSE
            ! FILL HERE TO CALCULATE DYNAMIC POTENTIAL ENERGY
            DO IPOT = 0, ORDER_TO_NMAX(EFT_RADIAL_AA%ORDER)
              FMAT_AA_RI(IPOT,ICH,IE,IAB,IAK) = POT_AA_RI(IPOT,IE,S,T,LL,LR)
              FMAT_AA_IR(IPOT,ICH,IE,IAB,IAK) = POT_AA_IR(IPOT,IE,S,T,LL,LR)
              FMAT_AA_RR(IPOT,ICH,IE,IAB,IAK) = POT_AA_RR(IPOT,IE,S,T,LL,LR)
              FMAT_AA_II(IPOT,ICH,IE,IAB,IAK) = POT_AA_II(IPOT,IE,S,T,LL,LR)
            ENDDO
          ENDIF

      ! SI CALCOLA HAMILTONIANA PER I VARI CASI:REGOLARE-IRREGOLARE(AM), IRREGOLARE-REGOLARE(AM1),
      !         REGOLARE-REGOLARE(AM2),IRREGOLARE-IRREGOLARE(AM3)
          IF (.NOT.USE_DYNAMIC) THEN 
            H_MINUS_E_AA_RI(ICH, IE, IAB, IAK) = (KIN_MIN_E_RI(IE,LL,LR)+APE) / HTM    ! RI
            H_MINUS_E_AA_IR(ICH, IE, IAB, IAK) =  APE1 / HTM                           ! IR
            H_MINUS_E_AA_RR(ICH, IE, IAB, IAK) =  APE2 / HTM                           ! RR
            H_MINUS_E_AA_II(ICH, IE, IAB, IAK) = (KIN_MIN_E_II(IE,LL,LR)+APE3) / HTM   ! II
          ELSE ! USE_DYNAMIC
            K_MINUS_E_AA_RR(ICH, IE, IAB, IAK) =  0.D0
            K_MINUS_E_AA_RI(ICH, IE, IAB, IAK) = KIN_MIN_E_RI(IE,LL,LR) / HTM
            K_MINUS_E_AA_IR(ICH, IE, IAB, IAK) =  0.D0
            K_MINUS_E_AA_II(ICH, IE, IAB, IAK) = KIN_MIN_E_II(IE,LL,LR) / HTM
          ENDIF ! ENDIF USE_DYNAMIC
          
        ENDDO !IAB
        ENDDO !IAK
      ENDDO ! IE
    ENDDO ! ICH
  END SUBROUTINE PREPARE_ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS
 
  FUNCTION ORDER_TO_NMAX(ORDER) RESULT(N)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ORDER
    INTEGER :: N

    SELECT CASE (ORDER)
      CASE (0)
        N = 0
      CASE (1)
        N = 4
      CASE (3)
        N = 7
      CASE DEFAULT
        PRINT *, "Error: Invalid order for EFT radial function"
        STOP
    END SELECT
  END FUNCTION ORDER_TO_NMAX

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
  FUNCTION IS_FIRST_CALL(J, L, S, TZ, IPOT, ILB, LEMP, RESET) RESULT(FIRST_CALL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP
    LOGICAL, INTENT(IN), OPTIONAL :: RESET
    INTEGER :: JOLD = -1, LOLD = -1, SOLD = -1, TZOLD = -1, IPOTOLD = -1, ILBOLD = -1, LEMPOOLD = -1
    LOGICAL :: FIRST_CALL

    IF (PRESENT(RESET)) THEN
      IF (RESET) THEN
        JOLD = -1
        LOLD = -1
        SOLD = -1
        TZOLD = -1
        IPOTOLD = -1
        ILBOLD = -1
        LEMPOOLD = -1
        FIRST_CALL = .TRUE.
        RETURN
      ENDIF
    ENDIF
    
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

    IF (BESSELS_SET) RETURN

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
    ENDDO
    BESSELS_SET = .TRUE.
    
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

    IF (VAR_P%IPOT==0) THEN
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
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_CC, LECS%RC, EFT_RADIAL_CC, ORDER_POTENTIAL = LECS%ORDER)
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_AC, LECS%RC, EFT_RADIAL_AC, ORDER_POTENTIAL = LECS%ORDER)
      CALL GET_EFT_RADIAL_FUNCTIONS(XX_AA, LECS%RC, EFT_RADIAL_AA, ORDER_POTENTIAL = LECS%ORDER)
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
    IF (.NOT.ENERGIES_SET) THEN
      PRINT *, "Error: ENERGIES not set, first set them before calling SET_NEW_LECS"
      STOP
    ENDIF
        IF (ANY(LECS_NEW%RC /= LECS%RC)) THEN
      NEW_CUTOFFS = .TRUE.
      IF (USE_DYNAMIC) PREPARE = .TRUE.
    ELSE
      NEW_CUTOFFS = .FALSE.
    ENDIF
    LECS = LECS_NEW
    LECS_SET = .TRUE.
    NEW_LECS = .TRUE.
  END SUBROUTINE SET_NEW_LECS
  

  !> \brief Reset the SCATTERING_NN_VARIATIONAL module to its initial state.
  !! This subroutine deallocates all allocatable arrays, resets logical flags,
  !! and restores all module variables to their default values.
  SUBROUTINE RESET_SCATTERING_NN_VARIATIONAL()
    IMPLICIT NONE
    TYPE(PHASE_SHIFT_RESULT) :: PSR
    LOGICAL :: TMP

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
    CALL RESET_CHANNEL(CHANNEL)

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
    PREPARE         = .TRUE. 
    NEW_LECS        = .TRUE.

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

    CALL NN_SCATTERING_VARIATIONAL(0.D0, 0, 0, 0, 0, 0, 0, 0, PSR, .FALSE., .FALSE., .TRUE.)
    TMP = IS_FIRST_CALL(0,0,0,0,0,0,0,.TRUE.)

  END SUBROUTINE RESET_SCATTERING_NN_VARIATIONAL

  SUBROUTINE NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PHASE_SHIFTS, IPOT, ILB, LECS_FOR_PLESS, &
        FIT_CONSTANTS, ORDER_OF_THE_FIT)
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNELS(:)
    DOUBLE PRECISION, INTENT(IN) :: ENERGIES(:)
    INTEGER, INTENT(IN) :: LEMP
    INTEGER, INTENT(IN), OPTIONAL :: IPOT, ILB
    TYPE(PHASE_SHIFT_RESULT), INTENT(OUT) :: PHASE_SHIFTS(:,:)
    TYPE(LECS_EFT_PLESS), INTENT(IN), OPTIONAL :: LECS_FOR_PLESS
    DOUBLE PRECISION, INTENT(IN), OPTIONAL, ALLOCATABLE :: FIT_CONSTANTS(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: ORDER_OF_THE_FIT

    DOUBLE PRECISION :: E
    INTEGER :: IPOT_, ILB_, IE, ICH
    INTEGER :: L, S, J, T, TZ
    INTEGER :: ORDER_OF_THE_FIT_
    DOUBLE PRECISION, ALLOCATABLE :: FIT_CONSTANTS_(:,:)

    ! Set the energies and channels
    NE = SIZE(ENERGIES)
    IF (NE <= 0) THEN
      PRINT *, "Error: No energies provided"
      STOP
    ENDIF
    CALL SET_ENERGIES(ENERGIES)
    NCHANNELS = SIZE(CHANNELS)
    IF (NCHANNELS <= 0) THEN
      PRINT *, "Error: No channels provided"
      STOP
    ENDIF
    CALL SET_CHANNELS(CHANNELS)
    
    IF (PRESENT(LECS_FOR_PLESS)) THEN
      IF (PRESENT(IPOT) .OR. PRESENT(ILB)) THEN
        PRINT *, "Error: Cannot set LECS and IPOT/ILB at the same time"
        STOP
      ENDIF
      IPOT_ = -1
      ILB_ = -1
      CALL SET_NEW_LECS(LECS_FOR_PLESS)
    ELSEIF (.NOT.PRESENT(IPOT)) THEN
      PRINT *, "Error: LECS and IPOT not set, set one of them"
      STOP
    ELSEIF (PRESENT(ILB)) THEN
      IPOT_ = IPOT
      IF (IPOT /= 19) THEN
        IF (ILB /= 1) THEN
          WRITE(*,*) "ILB for this IPOT could be only 1, setting it to 1"
        ENDIF
        ILB_ = 1
      ELSE
        ILB_ = ILB
      ENDIF
    ELSE
      IF (IPOT /= 19) THEN
        WRITE(*,*) "Setting ILB to 1 for IPOT=", IPOT
      ELSE
        WRITE(*,*) "ILB for IPOT=19 is not set, using default value 15"
        ILB_ = 15
      ENDIF
    ENDIF

    DO ICH = 1, NCHANNELS
      J = GET_CHANNEL_J(CHANNELS(ICH))
      TZ = GET_CHANNEL_TZ(CHANNELS(ICH))
      L = GET_CHANNEL_L(CHANNELS(ICH), 1)
      S = GET_CHANNEL_S(CHANNELS(ICH), 1)
      T = GET_CHANNEL_T(CHANNELS(ICH), 1)
      DO IE = 1, NE
        E = ENERGIES(IE)
        CALL NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT_, ILB_, LEMP, PHASE_SHIFTS(ICH,IE))
      ENDDO
    ENDDO
    CALL RESET_SCATTERING_NN_VARIATIONAL

    IF (PRESENT(FIT_CONSTANTS) .OR. PRESENT(ORDER_OF_THE_FIT)) THEN
      IF (.NOT.PRESENT(ORDER_OF_THE_FIT)) THEN
        WRITE(*,*) "FIT_CONSTATS: using default ORDER_OF_THE_FIT=2"
        ORDER_OF_THE_FIT_ = 2
      ELSE
        ORDER_OF_THE_FIT_ = ORDER_OF_THE_FIT
      ENDIF
      CALL FIT_CHANNELS_LOW_ENERGY(CHANNELS, K2, PHASE_SHIFTS, FIT_CONSTANTS_, ORDER_OF_THE_FIT_)
    ENDIF
  END SUBROUTINE NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS

  SUBROUTINE FIT_CHANNELS_LOW_ENERGY(CHANNELS_TO_FIT, ENERGIES, PHASE_SHIFTS, FIT_CONSTANTS, ORDER_OF_THE_FIT)
    USE FIT_MODULE
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNELS_TO_FIT(:)
    DOUBLE PRECISION, INTENT(IN) :: ENERGIES(:)
    TYPE(PHASE_SHIFT_RESULT), INTENT(IN) :: PHASE_SHIFTS(:,:)
    DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT) :: FIT_CONSTANTS(:,:)
    INTEGER, INTENT(IN) :: ORDER_OF_THE_FIT

    TYPE(SCATTERING_CHANNEL) :: CHANNEL_TO_FIT
    INTEGER :: NK2, NCH_TO_FIT, ICH

    NK2 = SIZE(ENERGIES)
    NCH_TO_FIT = SIZE(PHASE_SHIFTS, 1)
    CALL REALLOCATE_2D_1(FIT_CONSTANTS, NCH_TO_FIT, ORDER_OF_THE_FIT + 1)
    
    DO ICH = 1, NCH_TO_FIT
      CHANNEL_TO_FIT = CHANNELS_TO_FIT(ICH)
      CALL FIT_CHANNEL_LOW_ENERGY(CHANNEL_TO_FIT, ENERGIES, PHASE_SHIFTS(ICH,:), FIT_CONSTANTS, ORDER_OF_THE_FIT)
    ENDDO
  END SUBROUTINE FIT_CHANNELS_LOW_ENERGY

  SUBROUTINE FIT_CHANNEL_LOW_ENERGY(CHANNEL_TO_FIT, ENERGIES, PHASE_SHIFTS, FIT_CONSTANTS, ORDER_OF_THE_FIT, KSQUARED, K2L1COTD)
    USE FIT_MODULE
    IMPLICIT NONE
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL_TO_FIT
    DOUBLE PRECISION, INTENT(IN) :: ENERGIES(:)
    TYPE(PHASE_SHIFT_RESULT), INTENT(IN) :: PHASE_SHIFTS(:)
    DOUBLE PRECISION, INTENT(OUT) :: FIT_CONSTANTS(:,:)
    INTEGER, INTENT(IN) :: ORDER_OF_THE_FIT
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: KSQUARED(:)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: K2L1COTD(:,:)

    INTEGER :: I, IMIN, NK2, IEQ, NEQ, L
    DOUBLE PRECISION, ALLOCATABLE :: X(:), Y(:)
    DOUBLE PRECISION :: DELTA
    DOUBLE PRECISION :: FIT_CONSTANTS_(ORDER_OF_THE_FIT +1)

    NK2 = SIZE(ENERGIES)
    NEQ = GET_CHANNEL_NCH(CHANNEL_TO_FIT)
    CALL REALLOCATE_1D_2(X, Y, NK2)
    CALL SET_M_T1Z_T2Z_HTM(GET_CHANNEL_TZ(CHANNEL_TO_FIT))
    X = ENERGIES / HTM
    IF (PRESENT(KSQUARED)) THEN
      IF (SIZE(KSQUARED) /= NK2) THEN
        PRINT *, "Error: KSQUARED size does not match ENERGIES size"
        STOP
      ENDIF
      KSQUARED = X
    ENDIF
    
    DO IEQ = 1, NEQ
      L = GET_CHANNEL_L(CHANNEL_TO_FIT, IEQ)
      IMIN = -1
      DO I = 1, NK2
        IF (IEQ==1) DELTA = PHASE_SHIFTS(I)%delta1_BB
        IF (IEQ==2) DELTA = PHASE_SHIFTS(I)%delta2_BB
        IF (ABS(DELTA) > 1.D-6) THEN
          IF (PRESENT(K2L1COTD)) THEN
            K2L1COTD(IEQ,I) = X(I)**((2.D0*L + 1.D0)/2.D0) / DTAN(DELTA*PI/180.D0)
          ENDIF
          IF (IMIN == -1) IMIN = I
          Y(I) = X(I)**((2.D0*L + 1.D0)/2.D0) / DTAN(DELTA*PI/180.D0)
        ELSE
          Y(I) = 0.D0
          IMIN =-1
        ENDIF
      ENDDO
      FIT_CONSTANTS_ = 0.D0
      CALL POLYNOMIAL_REGRESSION(Y(IMIN:), X(IMIN:), ORDER_OF_THE_FIT, NK2-IMIN, FIT_CONSTANTS_)
      FIT_CONSTANTS(IEQ,:) = FIT_CONSTANTS_
    ENDDO
  END SUBROUTINE FIT_CHANNEL_LOW_ENERGY

END MODULE SCATTERING_NN_VARIATIONAL