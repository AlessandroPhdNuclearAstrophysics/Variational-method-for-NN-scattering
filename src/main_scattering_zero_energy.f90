PROGRAM SCATTERING_NN_ZERO_ENERGY
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS
  USE SCATTERING
  IMPLICIT NONE
  INTEGER, PARAMETER :: LMAX=2, JMAX=2
  INTEGER, PARAMETER :: TZ=0
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
  TYPE(PHASE_SHIFT_RESULT), ALLOCATABLE :: PS(:,:)
  TYPE(ZERO_ENERGY_OBSERVABLES) :: OBS
  DOUBLE PRECISION :: ENERGIES(1), R(2,2)
  INTEGER :: IPOT, ILB, LEMP=0
  INTEGER :: L, DIM

  IPOT = 18
  ILB = 1
  
  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  ALLOCATE(PS(SIZE(CHANNELS), 1))

  ENERGIES(1) = 0.D0
  CALL SET_CHANNEL(CHANNELS(1), 0, 0, 0, 0)
  CALL SET_CHANNEL(CHANNELS(2), 1, 0, 1, 0)
  CALL SET_CHANNEL(CHANNELS(3), 1, 1, 0, 0)
  CALL SET_CHANNEL(CHANNELS(4), 0, 1, 1, 0)

  CALL SET_ENERGIES(ENERGIES)
  CALL SET_CHANNELS(CHANNELS)
  
  CALL NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PS, IPOT=IPOT, ILB=ILB)
  
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"SINGLE CHANNEL RESULTS"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(1))
  R = PS(1,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(1), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(1))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.0433

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(3))
  R = PS(3,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(3), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(3))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.359170

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(4))
  R = PS(4,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(4), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(4))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.397155

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(5))
  R = PS(5,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(5), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(5))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.654212


  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(7))
  R = PS(7,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(7), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(7))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.721184

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(8))
  R = PS(8,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(8), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(8))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.135085



  

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"BLATT-BIEDERMANN PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(2), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(2))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.185057
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 0.3004
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.16822

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(6), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(6))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/3.413315
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", -5.512
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", 1/0.126009

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"STAPP PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(2), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(2))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.185057
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 1.62084
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.155404

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = GET_CHANNEL_L(CHANNELS(6), 1)
  DIM = GET_CHANNEL_NCH(CHANNELS(6))
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/3.413315
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 1.61436 
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", -1/1.03225

  DEALLOCATE(CHANNELS)
  DEALLOCATE(PS)

END PROGRAM SCATTERING_NN_ZERO_ENERGY