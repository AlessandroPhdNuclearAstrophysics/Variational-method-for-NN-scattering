!> \file main_scattering_NN_variazional_method.f90
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

PROGRAM SCATTERING_NN_VARIATIONAL_METHOD
  USE SCATTERING_NN_VARIATIONAL, ONLY: SET_VARIATIONAL_PARAMETERS, SET_ENERGIES, SET_CHANNELS, DUMP_MODULE_DATA
  USE QUANTUM_NUMBERS, ONLY: SCATTERING_CHANNEL, GET_CHANNEL_NAME, GET_CHANNEL_NCH, PREPARE_CHANNELS
  IMPLICIT NONE

  INTEGER, PARAMETER :: JMAX = 2, LMAX = 2

  !> \brief Number of energy points and number of channels.
  INTEGER :: NE = 20
  INTEGER :: NCH = 0
  INTEGER :: NEQ
  !> \brief Maximum energy [MeV].
  DOUBLE PRECISION :: EMAX = 1.D0
  !> \brief Array of all physical channels.
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)

  INTEGER :: I, TZ
  !> \brief Array of energy values.
  DOUBLE PRECISION, ALLOCATABLE :: ENERGIES(:)
  !> \brief Current energy and energy step.
  DOUBLE PRECISION :: E, HE
  !> \brief Potential, interaction, and EM flag.
  INTEGER :: IPOT, ILB, LEMP
  !> \brief Print coefficients flag.
  LOGICAL :: PRINT_COEFFICIENTS = .FALSE.
  !> \brief Namelist for input parameters.
  CHARACTER(LEN=256) :: OUT_DIR
  NAMELIST /IN/ EMAX, NE, TZ, IPOT, ILB, LEMP, PRINT_COEFFICIENTS, OUT_DIR

  !> \brief Input file name for namelist (if provided).
  CHARACTER(LEN=256) :: INPUT_FILE
  !> \brief Flag: are there command-line arguments?
  LOGICAL :: HAS_ARGUMENTS

  CALL SETUP_FROM_ARGS

  !> \brief Allocate and fill the energy grid.
  ALLOCATE(ENERGIES(NE))
  HE = EMAX / NE
  ENERGIES = (/ (I * HE, I = 1, NE) /)

  !> \brief Prepare the list of all physical channels.
  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  NCH = SIZE(CHANNELS)
  CALL SET_ENERGIES(ENERGIES)
  CALL SET_CHANNELS(CHANNELS)

  !> \brief Loop over all possible L, S, J combinations and compute phase shifts for each channel.
  BLOCK ! MAIN LOOP
    USE SCATTERING_NN_VARIATIONAL, ONLY: PHASE_SHIFT_RESULT, NN_SCATTERING_VARIATIONAL
    !> \brief Structure to hold phase shift results.
    TYPE(PHASE_SHIFT_RESULT) :: PHASE_SHIFTS
    CHARACTER(LEN=16) :: CHANNEL_NAME
    INTEGER :: ICH, LMIN, S, J

    DO ICH = 1, NCH

      LMIN = CHANNELS(ICH)%L()
      S = CHANNELS(ICH)%S()
      J = CHANNELS(ICH)%J()
      TZ = CHANNELS(ICH)%TZ()

      CHANNEL_NAME = GET_CHANNEL_NAME(CHANNELS(ICH))
      CALL SET_VARIATIONAL_PARAMETERS(J, LMIN, S, TZ, IPOT, ILB, LEMP)
      PRINT *, "Scattering channel name: ", TRIM(CHANNEL_NAME)
      NEQ = GET_CHANNEL_NCH(CHANNELS(ICH))
      OPEN(21, FILE=TRIM(OUT_DIR)//'delta_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
      IF (NEQ == 2) OPEN(22, FILE=TRIM(OUT_DIR)//'delta_BB_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
      DO I = 1, NE
        E =  ENERGIES(I)
        CALL NN_SCATTERING_VARIATIONAL(E, J, LMIN, S, TZ, IPOT, ILB, LEMP, PHASE_SHIFTS, PRINT_COEFFICIENTS=.FALSE.)
        WRITE(21, *) E, PHASE_SHIFTS%delta1_S, PHASE_SHIFTS%delta2_S, PHASE_SHIFTS%epsilon_S
        IF (NEQ == 2) WRITE(22, *) E, PHASE_SHIFTS%delta1_BB, PHASE_SHIFTS%delta2_BB, PHASE_SHIFTS%epsilon_BB
      ENDDO
      CLOSE(21)
      IF (NEQ == 2) CLOSE(22)
    ENDDO
  END BLOCK ! MAIN LOOP

  CALL CLEANUP

CONTAINS
  SUBROUTINE SETUP_FROM_ARGS
    !> \brief Check if there are command-line arguments.
    HAS_ARGUMENTS = COMMAND_ARGUMENT_COUNT() > 0

    !> \brief Initialize default values for quantum numbers and options.
    TZ = 0
    IPOT = 18
    ILB = 1
    LEMP = 0
    OUT_DIR = 'output/AV18/'


    !> \brief Read namelist from file if provided, otherwise prompt user.
    IF (HAS_ARGUMENTS) THEN
      CALL GET_COMMAND_ARGUMENT(1, INPUT_FILE)
      OPEN(UNIT=10, FILE=TRIM(INPUT_FILE), STATUS='OLD', ACTION='READ')
      READ(10, NML=IN)
      CLOSE(10)
    ELSE
      ! If no arguments, write the namelist and allow user to modify it
      WRITE(*, NML=IN)
      PRINT *, "Modify the above values if needed and press Enter to continue."
      READ(*, NML=IN)
    END IF

    !> \brief Ensure OUT_DIR ends with a slash.
    IF (LEN_TRIM(OUT_DIR) > 0) THEN
      IF (OUT_DIR(LEN_TRIM(OUT_DIR):LEN_TRIM(OUT_DIR)) /= '/') THEN
        OUT_DIR = TRIM(OUT_DIR) // '/'
      END IF
    END IF

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

END PROGRAM SCATTERING_NN_VARIATIONAL_METHOD