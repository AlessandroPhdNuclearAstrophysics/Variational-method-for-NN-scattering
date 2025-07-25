PROGRAM LECS_FROM_GRID
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  INTEGER, PARAMETER :: N_CHANNELS = 3
  TYPE(SCATTERING_CHANNEL) :: CHANNELS(N_CHANNELS)

  CHANNELS(1) = GET_CHANNEL_FROM_NAME('3P0')
  CHANNELS(2) = GET_CHANNEL_FROM_NAME('3P1')
  CHANNELS(3) = GET_CHANNEL_FROM_NAME('3P2-3F2')

  BLOCK ! Print channels
    INTEGER :: I
    DO I = 1, N_CHANNELS
      CALL CHANNELS(I)%PRINT()
    END DO
  END BLOCK

  BLOCK ! GRID
    USE SCATTERING_NN_VARIATIONAL, ONLY:  PHASE_SHIFT_RESULT, SET_ENERGIES, SET_CHANNELS, SET_NEW_LECS, &
                                          NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS, SET_DYNAMIC
    USE EFT_PLESS, ONLY: LECS_EFT_PLESS
    TYPE(LECS_EFT_PLESS) :: LECS
    INTEGER, PARAMETER :: NE = 10
    INTEGER, PARAMETER :: LEMP = 0
    INTEGER, PARAMETER :: FIT_ORDER = 2
    DOUBLE PRECISION, DIMENSION(NE), PARAMETER :: ENERGIES = (/ 0.1D0, 0.2D0, 0.3D0, 0.4D0, 0.5D0, 0.6D0, 0.7D0, 0.8D0, 0.9D0, 1.0D0 /)
    LOGICAL :: PRINT_ = .FALSE.
    LOGICAL :: FITTED
    TYPE(PHASE_SHIFT_RESULT) :: PHASES(N_CHANNELS, NE)
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ABC

    CALL SET_DYNAMIC(.TRUE.)

    LECS%ADIMENSIONAL = .FALSE.
    LECS%ORDER = 1
    LECS%CIT(0) = 0.219960910D-1 ! From article

    BLOCK ! LOOPS
      INTEGER, PARAMETER :: NR = 5, N11 = 40, N6 = 8, N7 = 8
      DOUBLE PRECISION, PARAMETER :: C11MAX = 10.D0, C6MAX = 10.D0, C7MAX = 10.D0
      DOUBLE PRECISION, PARAMETER :: RMIN = 1.D0, RMAX = 2.6D0, HR = (RMAX - RMIN) / DBLE(NR-1)
      DOUBLE PRECISION :: H11 = 2*C11MAX / DBLE(N11-1)
      DOUBLE PRECISION :: H6 = 2*C6MAX / DBLE(N6-1)
      DOUBLE PRECISION :: H7 = 2*C7MAX / DBLE(N7-1)
      DOUBLE PRECISION :: R
      INTEGER :: IR, I11, I6, I7, I


      DO IR = 1, NR
        R = RMIN + (IR-1)*HR
      DO I11= 1, N11
      DO I6 = 1, N6
      DO I7 = 1, N7
        LECS%RC(1,1) = R
        LECS%CNLO(4) = -C11MAX + H11*(I11-1)
        LECS%CNLO(6) = -C6MAX + H6*(I6-1)
        LECS%CNLO(7) = -C7MAX + H7*(I7-1)

        CALL NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PHASES, LECS_FOR_PLESS= LECS, &
                                                          FIT_CONSTANTS=ABC, ORDER_OF_THE_FIT = FIT_ORDER, &
                                                          FITTED=FITTED)
        IF (FITTED) THEN
          PRINT *, 'Fitted LECs for R = ', R, ' C11 = ', LECS%CNLO(4), ' C6 = ', LECS%CNLO(6), ' C7 = ', LECS%CNLO(7)
          DO I = 1, N_CHANNELS
            PRINT *, 'Channel: ', CHANNELS(I)%NAME(), ' a ,b ,c: ', ABC(I,1,1), ABC(I,1,2), ABC(I,1,3)
            IF (CHANNELS(I)%IS_COUPLED()) THEN
              PRINT *, 'Cahannel: ', CHANNELS(I)%NAME(), ' a ,b ,c: ', ABC(I,2,1), ABC(I,2,2), ABC(I,2,3)
            ENDIF
          END DO
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    END BLOCK ! LOOPS
  END BLOCK ! GRID


END PROGRAM LECS_FROM_GRID
