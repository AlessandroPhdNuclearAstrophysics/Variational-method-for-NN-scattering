! Here I try first to find the depth of the Gaussian assuming no polynomial
! behaviour, so the potential looks like
!           V(r) = C exp(-r^2/R^2)
! To do this I use C_7 which is already multiplied by only a constant

PROGRAM FITTER_3PX_CHANNELS_ONLY_DEPTH_CONSTANT
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS
  USE SCATTERING
  USE NUMBER_DIFFERENCES
  USE EFT_PLESS, ONLY: LECS_EFT_PLESS, GET_LECS
  IMPLICIT NONE

  ! Parameters for my case of interest
  INTEGER, PARAMETER :: LMAX = 0, JMAX = 0, TZ = 0
  INTEGER, PARAMETER :: LEMP = 0
  INTEGER, PARAMETER :: L = 1, S = 1

  ! Fit to
  DOUBLE PRECISION, PARAMETER :: A0_18(3) = -1.D0 / (/ 0.397417, -0.654595, 3.415  /)
  DOUBLE PRECISION, PARAMETER :: R0_18(3) =  2.D0 * (/ 1.8846 , -4.2900, -4.8 /)
  

  ! Size parameters and some physical constants
  INTEGER, PARAMETER :: NE = 2, NCHANNELS = 3, RDIM = 2, FIT_ORDER = 2
  DOUBLE PRECISION, PARAMETER :: K_SMALL = 1.d-8
  DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)

  ! Variables for logging and messages
  CHARACTER(LEN=1024) :: MESSAGE

  ! Inputs and outputs for the scattering modules
  TYPE(SCATTERING_CHANNEL) :: CHANNELS(NCHANNELS)
  DOUBLE PRECISION :: ENERGIES(NE)
  TYPE(PHASE_SHIFT_RESULT) :: PS(NCHANNELS, NE)
  TYPE(ZERO_ENERGY_OBSERVABLES) :: ZEOBS
  TYPE(LECS_EFT_PLESS) :: LECS

  ! R_BB matrices
  DOUBLE PRECISION RBB0(RDIM,RDIM), RBB(RDIM,RDIM), K2, HTM, R0, A0, E

  INTEGER :: NEQ, ICH
  INTEGER, EXTERNAL :: DOUBLE_FACTORIAL

  ! Read E from first argument if provided
  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, MESSAGE)
    READ(MESSAGE, *) E
  ELSE
    STOP "No energy provided. Please provide an energy value as the first argument."
  END IF
  ENERGIES(1) = 0
  ENERGIES(2) = E

  CALL SET_MAX_LOG_LEVEL(0)

  CALL SET_CHANNEL(CHANNELS(1), 0, L, S, TZ)
  CALL SET_CHANNEL(CHANNELS(2), 1, L, S, TZ)
  CALL SET_CHANNEL(CHANNELS(3), 2, L, S, TZ)

  CALL SET_DYNAMIC(.TRUE.)
  LECS = GET_LECS(10)

  BLOCK ! Loop over parameters
  DOUBLE PRECISION, PARAMETER :: RMIN = 1.D0, HR = 0.4D0, RMAX = 2.6D0
  DOUBLE PRECISION, PARAMETER :: C7MAX = 15.D0
  INTEGER, PARAMETER :: NR = INT((RMAX-RMIN)/HR) + 1
  INTEGER, PARAMETER :: N7 = 30000
  DOUBLE PRECISION, PARAMETER :: H7 = 2*C7MAX / (N7 - 1)
  INTEGER :: I7, IR

  LECS%CNLO = 0.D0
  DO IR = 1, NR
    LECS%RC(1,1) = RMIN + (IR-1) * HR
    CALL SET_DYNAMIC(.TRUE.)
    DO I7 = 1, N7
      ! Setting the LECs for the current iteration
      LECS%CNLO(7) =-C7MAX   + (I7  - 1) * H7
    
      CALL NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PS, LECS_FOR_PLESS=LECS)
      HTM = GET_HTM()

      K2 = ENERGIES(2)/HTM
      DO ICH = 1, NCHANNELS
        NEQ = CHANNELS(ICH)%NCH()

        ! Saving the R_BB matrix at zero energy and saving the scattering lengths
        RBB0 = PS(ICH,1)%R_BB
        ZEOBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(RBB0, NEQ, L)
        A0 = ZEOBS%a1
        
        ! Evaluating r_e using the R_BB matrix at non-zero energy
        RBB = PS(ICH,2)%R_BB
        IF ( SQRT(K2) > K_SMALL) RBB(1,1) =   RBB(1,1) * DOUBLE_FACTORIAL(2*L+1)**2 / SQRT(K2)**(2*L+1)
        R0 = -2 * DOUBLE_FACTORIAL(2*L+1)**2 * (RBB(1,1) - RBB0(1,1))/(K2 * RBB0(1,1)**2);

        IF ( ABS_DIFF_PROCENT(A0_18(ICH), A0) < 50.D0 .AND. ABS_DIFF_PROCENT(R0_18(ICH), R0) < 50.D0 ) THEN
          WRITE(200+ICH-1,'(F6.2, 3F10.4, 2F20.5)') LECS%RC(1,1), LECS%CNLO(4), LECS%CNLO(6), LECS%CNLO(7), A0, R0
        ENDIF
      ENDDO
    ENDDO
    CALL RESET_SCATTERING_NN_VARIATIONAL
  ENDDO
  END BLOCK


END PROGRAM FITTER_3PX_CHANNELS_ONLY_DEPTH_CONSTANT