PROGRAM SCATTERING_NN_VARIATIONAL_METHOD
  USE SCATTERING_NN_VARIATIONAL
  IMPLICIT NONE

  LOGICAL, EXTERNAL :: CHECK_QUANTUM_NUMBERS

  ! Variable declarations
  INTEGER :: J, L, S, TZ
  DOUBLE PRECISION :: E
  INTEGER :: IPOT, ILB, LEMP
  NAMELIST /IN/ E, J, L, S, TZ, IPOT, ILB, LEMP

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

  IF (.NOT. CHECK_QUANTUM_NUMBERS(J, L, S, TZ)) THEN
    PRINT *, "Invalid quantum numbers. Please check the input."
    STOP
  END IF


  CALL NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP)

END PROGRAM SCATTERING_NN_VARIATIONAL_METHOD