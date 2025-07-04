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
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS
  IMPLICIT NONE

  !> \brief Number of energy points and number of channels.
  INTEGER :: NE = 200, NCH = 0
  INTEGER :: NEQ
  !> \brief Maximum energy [MeV].
  DOUBLE PRECISION :: EMAX = 1.D0
  !> \brief Name of the scattering channel (for output).
  CHARACTER(LEN=16) :: CHANNEL_NAME
  !> \brief Temporary channel object.
  TYPE(SCATTERING_CHANNEL) :: CHANNEL
  !> \brief Array of all physical channels.
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)

  !> \brief Quantum numbers and loop variables.
  INTEGER :: J, L, S, TZ, I
  !> \brief Array of energy values.
  DOUBLE PRECISION, ALLOCATABLE :: ENERGIES(:)
  !> \brief Current energy and energy step.
  DOUBLE PRECISION :: E, HE
  !> \brief Potential, interaction, and EM flag.
  INTEGER :: IPOT, ILB, LEMP
  !> \brief Print coefficients flag.
  LOGICAL :: PRINT_COEFFICIENTS = .FALSE.
  !> \brief Structure to hold phase shift results.
  TYPE(PHASE_SHIFT_RESULT) :: PHASE_SHIFTS
  !> \brief Namelist for input parameters.
  CHARACTER(LEN=256) :: OUT_DIR
  NAMELIST /IN/ EMAX, NE, TZ, IPOT, ILB, LEMP, PRINT_COEFFICIENTS, OUT_DIR

  !> \brief Input file name for namelist (if provided).
  CHARACTER(LEN=256) :: input_file
  !> \brief Flag: are there command-line arguments?
  LOGICAL :: has_arguments

  !> \brief Check if there are command-line arguments.
  has_arguments = COMMAND_ARGUMENT_COUNT() > 0

  !> \brief Initialize default values for quantum numbers and options.
  TZ = 0
  IPOT = 18
  ILB = 1
  LEMP = 0
  OUT_DIR = 'output/AV18/'


  !> \brief Read namelist from file if provided, otherwise prompt user.
  IF (has_arguments) THEN
    CALL GET_COMMAND_ARGUMENT(1, input_file)
    OPEN(UNIT=10, FILE=TRIM(input_file), STATUS='OLD', ACTION='READ')
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

  !> \brief Allocate and fill the energy grid.
  ALLOCATE(ENERGIES(NE))
  HE = EMAX / NE
  ENERGIES = (/ (I * HE, I = 1, NE) /)

  !> \brief Prepare the list of all physical channels.
  CALL PREPARE_CHANNELS(2, 2, 0, CHANNELS)
  NCH = SIZE(CHANNELS)
  CALL SET_ENERGIES(ENERGIES)
  CALL SET_CHANNELS(CHANNELS)

  !> \brief Loop over all possible L, S, J combinations and compute phase shifts for each channel.
  DO L = 0, 2
    DO S = 0, 1
      DO J = ABS(L-S), MIN(L+S, 2)
        CALL SET_CHANNEL(CHANNEL, J, L, S, TZ)
        IF (.NOT.IS_PHYSICAL_CHANNEL(CHANNEL)) CYCLE
        IF ( J .NE. 0 .AND. L > J .AND. S == 1) CYCLE
        CHANNEL_NAME = get_channel_name(CHANNEL)
        CALL SET_VARIATIONAL_PARAMETERS(J, L, S, TZ, IPOT, ILB, LEMP)
        PRINT *, "Scattering channel name: ", TRIM(CHANNEL_NAME)
        NEQ = GET_CHANNEL_NCH(CHANNEL)
        OPEN(21, FILE=TRIM(OUT_DIR)//'delta_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
        IF (NEQ == 2) OPEN(22, FILE=TRIM(OUT_DIR)//'delta_BB_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
        DO I = 1, NE
          E =  ENERGIES(I)
          CALL NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP, PHASE_SHIFTS, PRINT_COEFFICIENTS=.FALSE.)
          WRITE(21, *) E, PHASE_SHIFTS%delta1_S, PHASE_SHIFTS%delta2_S, PHASE_SHIFTS%epsilon_S
          IF (NEQ == 2) WRITE(22, *) E, PHASE_SHIFTS%delta1_BB, PHASE_SHIFTS%delta2_BB, PHASE_SHIFTS%epsilon_BB
        ENDDO
        CLOSE(21)
        IF (NEQ == 2) CLOSE(22)
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE(ENERGIES)
  DO I = 1, NCH
    CALL RESET_CHANNEL(CHANNELS(I))
  ENDDO
  DEALLOCATE(CHANNELS)
  CALL RESET_CHANNEL(CHANNEL)

  CALL RESET_SCATTERING_NN_VARIATIONAL

END PROGRAM SCATTERING_NN_VARIATIONAL_METHOD