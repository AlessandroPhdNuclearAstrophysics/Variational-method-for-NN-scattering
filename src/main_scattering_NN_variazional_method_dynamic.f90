!> \file main_scattering_NN_variazional_method_dynamic.f90
!! \brief Main program for computing NN scattering phase shifts using the variational method.
!!
!! This program sets up quantum numbers, energies, and channels for nucleon-nucleon
!! scattering, then computes phase shifts and mixing angles for all physical channels
!! using the variational approach. Results are written to output files for each channel.
!!
!! \details
!! - Reads input parameters from a namelist or prompts the user.
!! - Allocates energy grid and prepares all physical channels.
!! - Loops over all L, S, J combinations and computes phase shifts.
!! - Output: Files named delta_XXX.dat for each channel.
!!
!! \author Alessandro
!! \date 2025

PROGRAM SCATTERING_NN_VARIATIONAL_METHOD_DYNAMIC
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS, ONLY: SCATTERING_CHANNEL, GET_CHANNEL_NAME, GET_CHANNEL_NCH, PREPARE_CHANNELS, GET_CHANNEL_FROM_NAME
  USE EFT_PLESS, ONLY: LECS_EFT_PLESS, GET_LECS
  IMPLICIT NONE

  INTEGER, PARAMETER :: JMAX = 2, LMAX = 2
  INTEGER, PARAMETER :: LEMP = 0
  DOUBLE PRECISION, PARAMETER :: PI = 4.D0 * ATAN(1.D0)
  INTEGER :: NE = 200
  INTEGER :: NCH = 0
  INTEGER :: NEQ
  DOUBLE PRECISION :: EMAX = 1.D0, HTM
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)

  INTEGER :: I, TZ
  DOUBLE PRECISION, ALLOCATABLE :: ENERGIES(:), K2(:), K3COT(:)
  DOUBLE PRECISION :: E, HE, RANGE =40., H = 0.01
  CHARACTER(LEN=256) :: OUT_DIR
  NAMELIST /IN/ EMAX, NE, TZ, OUT_DIR, RANGE, H

  CHARACTER(LEN=256) :: ARGUMENT
  INTEGER :: NARGS
  TYPE(LECS_EFT_PLESS) :: LECS
  TYPE(PHASE_SHIFT_RESULT), ALLOCATABLE :: PS(:)


  CALL SETUP_FROM_ARGS

  !> \brief Allocate and fill the energy grid.
  ALLOCATE(ENERGIES(NE), K2(NE), K3COT(NE), PS(NE))
  HE = EMAX / NE
  ENERGIES = (/ (I * HE, I = 1, NE) /)

  !> \brief Prepare the list of all physical channels.
  CALL SET_VARIATIONAL_PARAMETERS( RANGE = RANGE, H = H )
  CALL SET_MAX_LOG_LEVEL(1)
  CALL SET_DYNAMIC(.TRUE.)
  
  LECS = GET_LECS(10)
  LECS%RC(1,1) = 2.6D0
  CALL LECS%CONVERT_TO_ADIMENSIONAL(1,1)
  LECS%CNLO(4) = 4.D0
  LECS%CNLO(6) = 0.5D0
  LECS%CNLO(7) =-0.5D0
  CALL LECS%CONVERT_TO_DIMENSIONAL(1,1)

  CALL SET_ENERGIES(ENERGIES)
  CALL SET_NEW_LECS(LECS)
  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  NCH = SIZE(CHANNELS)
  CALL SET_CHANNELS(CHANNELS)

  !> \brief Loop over all possible L, S, J combinations and compute phase shifts for each channel.
  BLOCK ! MAIN LOOP
    USE SCATTERING_NN_VARIATIONAL, ONLY: NN_SCATTERING_VARIATIONAL, GET_HTM, FIT => FIT_CHANNEL_LOW_ENERGY
    !> \brief Structure to hold phase shift results.
    CHARACTER(LEN=16) :: CHANNEL_NAME
    INTEGER :: ICH, LMIN, S, J

    DO ICH = 1, NCH
      LMIN = CHANNELS(ICH)%L()
      S = CHANNELS(ICH)%S()
      J = CHANNELS(ICH)%J()
      TZ = CHANNELS(ICH)%TZ()
      HTM = GET_HTM()

      CHANNEL_NAME = GET_CHANNEL_NAME(CHANNELS(ICH))
      CALL SET_VARIATIONAL_PARAMETERS(J=J, L=LMIN, S=S, T=TZ, LEMP=LEMP)
      PRINT *, "Scattering channel name: ", TRIM(CHANNEL_NAME)
      NEQ = GET_CHANNEL_NCH(CHANNELS(ICH))
      OPEN(21, FILE=TRIM(OUT_DIR)//'delta_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
      IF (NEQ == 2) OPEN(22, FILE=TRIM(OUT_DIR)//'delta_BB_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
      DO I = 1, NE
        E =  ENERGIES(I)
        K2(I) = E / HTM
        
        CALL NN_SCATTERING_VARIATIONAL(E, J, LMIN, S, TZ, -1, -1, LEMP, PS(I), PRINT_COEFFICIENTS=.FALSE.)
        
        K3COT(I) = K2(I)**((2*LMIN+1.D0)/2) / TAN(PS(I)%delta1_S * PI / 180.D0)

        WRITE(21, *) E, PS(I)%delta1_S, PS(I)%delta2_S, PS(I)%epsilon_S
        IF (NEQ == 2) WRITE(22, *) E, PS(I)%delta1_BB, PS(I)%delta2_BB, PS(I)%epsilon_BB
      ENDDO

      CLOSE(21)
      IF (NEQ == 2) CLOSE(22)

    ENDDO

  END BLOCK ! MAIN LOOP

  CALL CLEANUP

CONTAINS
  SUBROUTINE SETUP_FROM_ARGS
    !> \brief Initialize default values for quantum numbers and options.
    RANGE = 80.D0
    TZ = 0
    OUT_DIR = 'output/EFT_pless_23/'

    WRITE(*,IN)
    ! READ(*,NML=IN)
    WRITE(*,IN)

    !> \brief Create output directory if it does not exist.
    CALL SYSTEM('mkdir -p "' // TRIM(OUT_DIR) // '"')

    !> \brief Print the final values of the namelist for confirmation.
    PRINT *, "Final values of the namelist:"
    WRITE(*, NML=IN)
  END SUBROUTINE SETUP_FROM_ARGS



  SUBROUTINE CLEANUP
    USE SCATTERING_NN_VARIATIONAL, ONLY: RESET_SCATTERING_NN_VARIATIONAL
    USE QUANTUM_NUMBERS, ONLY: RESET_CHANNEL
    INTEGER :: ICH
    !> \brief Cleanup routine to deallocate resources.
    IF (ALLOCATED(ENERGIES)) DEALLOCATE(ENERGIES)
    DO ICH=1, NCH
      CALL RESET_CHANNEL(CHANNELS(ICH))
    END DO
    IF (ALLOCATED(CHANNELS)) DEALLOCATE(CHANNELS)
    CALL RESET_SCATTERING_NN_VARIATIONAL
  END SUBROUTINE CLEANUP

END PROGRAM SCATTERING_NN_VARIATIONAL_METHOD_DYNAMIC