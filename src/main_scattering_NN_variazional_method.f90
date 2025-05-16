PROGRAM SCATTERING_NN_VARIATIONAL_METHOD
  USE SCATTERING_NN_VARIATIONAL
  IMPLICIT NONE

  INTEGER :: NE = 200
  DOUBLE PRECISION :: EMAX = 1.D0

  ! Variable declarations
  INTEGER :: J, L, S, TZ, I
  DOUBLE PRECISION, ALLOCATABLE :: ENERGIES(:)
  DOUBLE PRECISION :: E, HE
  INTEGER :: IPOT, ILB, LEMP
  LOGICAL :: PRINT_COEFFICIENTS = .FALSE.
  TYPE(PHASE_SHIFT_RESULT) :: PHASE_SHIFTS
  NAMELIST /IN/ EMAX, NE, J, L, S, TZ, IPOT, ILB, LEMP, PRINT_COEFFICIENTS

  ! Local variables
  CHARACTER(LEN=256) :: input_file
  LOGICAL :: has_arguments

  ! Check if there are command-line arguments
  has_arguments = COMMAND_ARGUMENT_COUNT() > 0

  ! Initialize default values
  J = 1
  L = 0
  S = 1
  TZ = 0
  E = 1.0D0
  IPOT = 18
  ILB = 1
  LEMP = 0

  ! If there are arguments, treat the first as the namelist input file
  IF (has_arguments) THEN
    CALL GET_COMMAND_ARGUMENT(1, input_file)
    OPEN(UNIT=10, FILE=TRIM(input_file), STATUS='OLD', ACTION='READ')
    READ(10, NML=IN)
    CLOSE(10)
  ELSE
    ! If no arguments, write the namelist and allow user to modify it
    WRITE(*, NML=IN)
    PRINT *, "Modify the above values if needed and press Enter to continue."
    ! READ(*, NML=IN)
  END IF

  ! Print the final values of the namelist
  PRINT *, "Final values of the namelist:"
  WRITE(*, NML=IN)

  ALLOCATE(ENERGIES(NE))
  HE = EMAX / NE
  ENERGIES = (/ (I * HE, I = 1, NE) /)

  DO I = 1, NE
    E =  ENERGIES(I)
    PHASE_SHIFTS = NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP, PRINT_COEFFICIENTS=.FALSE.)
    WRITE(20, *) E, PHASE_SHIFTS%delta1_BB, PHASE_SHIFTS%delta2_BB, PHASE_SHIFTS%epsilon_BB
    WRITE(21, *) E, PHASE_SHIFTS%delta1_S, PHASE_SHIFTS%delta2_S, PHASE_SHIFTS%epsilon_S
  ENDDO

END PROGRAM SCATTERING_NN_VARIATIONAL_METHOD