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
  INTEGER :: NE = 10
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
  LECS%RC(1,1) = 1.D0
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

  OPEN(UNIT=1100, FILE="convergences_1S0.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1301, FILE="convergences_3S1.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1321, FILE="convergences_3D1.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1111, FILE="convergences_1P1.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1310, FILE="convergences_3P0.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1311, FILE="convergences_3P1.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1312, FILE="convergences_3P2.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1122, FILE="convergences_1D2.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')
  OPEN(UNIT=1322, FILE="convergences_3D2.dat", POSITION='APPEND', STATUS='UNKNOWN', ACTION='WRITE')

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
      OPEN(22, FILE=TRIM(OUT_DIR)//'kcotd_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
      IF (NEQ == 2) OPEN(22, FILE=TRIM(OUT_DIR)//'delta_BB_'//TRIM(CHANNEL_NAME)//'.dat', STATUS='UNKNOWN', ACTION='WRITE')
      DO I = 1, NE
        E =  ENERGIES(I)
        K2(I) = E / HTM
        
        CALL NN_SCATTERING_VARIATIONAL(E, J, LMIN, S, TZ, -1, -1, LEMP, PS(I), PRINT_COEFFICIENTS=.FALSE.)
        
        K3COT(I) = K2(I)**((2*LMIN+1.D0)/2) / TAN(PS(I)%delta1_S * PI / 180.D0)

        WRITE(21, *) E, PS(I)%delta1_S, PS(I)%delta2_S, PS(I)%epsilon_S
        WRITE(22, *) K2(I), K3COT(I)
        IF (NEQ == 2) WRITE(22, *) E, PS(I)%delta1_BB, PS(I)%delta2_BB, PS(I)%epsilon_BB
      ENDDO

      BLOCK
        USE SCATTERING_NN_VARIATIONAL, ONLY: PHASE_SHIFT_RESULT
        DOUBLE PRECISION :: COEFFS(2,3)
        LOGICAL :: FITTED
        IF (CHANNELS(ICH) == GET_CHANNEL_FROM_NAME('3P2-3F2')) THEN
          CHANNELS(ICH) = GET_CHANNEL_FROM_NAME('3P2')
        ENDIF

        FITTED = FIT(CHANNELS(ICH), ENERGIES, PS, COEFFS, 2)
        
        IF (FITTED) THEN 
          BLOCK 
            INTEGER :: FILE_ID   
            INTEGER :: FILE_ID_2 
            FILE_ID = 1000 + (2*S+1)*100 + LMIN*10 + J
            FILE_ID_2 = 1000 + (2*S+1)*100 + (LMIN+2)*10 + J
            write(*,*) "oadosaidjaosid ", file_id
            WRITE(*,*) -1/COEFFS(1,1)
            WRITE(*,*) 2*COEFFS(1,2)
            WRITE(*,*) COEFFS(1,3)
            WRITE(FILE_ID,*) RANGE, H, -1/COEFFS(1,1), 2*COEFFS(1,2),  COEFFS(1,3)
            IF (CHANNELS(ICH)%NCH() == 2) THEN
              WRITE(*,*) -1/COEFFS(2,1)
              WRITE(*,*) 2*COEFFS(2,2)
              WRITE(*,*) COEFFS(2,3)
              WRITE(FILE_ID_2,*) RANGE, H, -1/COEFFS(2,1), 2*COEFFS(2,2), COEFFS(2,3)
            END IF
            END BLOCK
          ELSE
            WRITE(*,*) "Fitting failed for channel: ", TRIM(CHANNEL_NAME)
          END IF
      END BLOCK

      CLOSE(21)
      CLOSE(22)

      IF (NEQ == 2) CLOSE(22)

    ENDDO
  END BLOCK ! MAIN LOOP
  CLOSE(1100)
  CLOSE(1301)
  CLOSE(1321)
  CLOSE(1111)
  CLOSE(1310)
  CLOSE(1311)
  CLOSE(1312)
  CLOSE(1122)
  CLOSE(1322)

  CALL CLEANUP

CONTAINS
  SUBROUTINE SETUP_FROM_ARGS
    !> \brief Check if there are command-line arguments.
    NARGS = COMMAND_ARGUMENT_COUNT()

    !> \brief Initialize default values for quantum numbers and options.
    TZ = 0
    OUT_DIR = 'output/EFT_pless_dynamic/'


    !> \brief Read namelist from file if provided, otherwise prompt user.
    IF (NARGS==2) THEN
      CALL GET_COMMAND_ARGUMENT(1, ARGUMENT)
      ! If the first argument is provided, treat it as RANGE
      READ(ARGUMENT, *) RANGE
      PRINT *, "RANGE set from command-line argument: ", RANGE
      CALL GET_COMMAND_ARGUMENT(2, ARGUMENT)
      ! If the second argument is provided, treat it as H
      READ(ARGUMENT, *) H
      PRINT *, "H set from command-line argument: ", H
    ELSEIF (NARGS==0) THEN
      ! If no arguments, write the namelist and allow user to modify it
      WRITE(*, NML=IN)
      PRINT *, "Modify the above values if needed and press Enter to continue."
      READ(*, NML=IN)
    ELSE
      STOP "Invalid number of command-line arguments. Expected 0 or 2 arguments."
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

END PROGRAM SCATTERING_NN_VARIATIONAL_METHOD_DYNAMIC