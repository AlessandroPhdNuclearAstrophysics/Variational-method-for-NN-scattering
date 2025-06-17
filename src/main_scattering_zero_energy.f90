PROGRAM SCATTERING_NN_ZERO_ENERGY
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  INTEGER, PARAMETER :: LMAX=2, JMAX=2
  INTEGER, PARAMETER :: TZ=0
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
  TYPE(PHASE_SHIFT_RESULT), ALLOCATABLE :: PS(:,:)
  DOUBLE PRECISION :: ENERGIES(1), HTM, K
  INTEGER :: IPOT, ILB, LEMP=0

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
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(1))
  WRITE(*,*) "a_1", PS(1,1)%R_BB(1,1)
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.0433

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  WRITE(*,*) "a_1", PS(2,1)%R_BB(1,1)
  WRITE(*,*) "e  ", PS(2,1)%R_BB(1,2)
  WRITE(*,*) "a_2", PS(2,1)%R_BB(2,2)
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.185057
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.166542

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(3))
  WRITE(*,*) "a_1", PS(3,1)%R_BB(1,1)
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.359170

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(4))
  WRITE(*,*) "a_1", PS(4,1)%R_BB(1,1)
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.397155

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(5))
  WRITE(*,*) "a_1", PS(5,1)%R_BB(1,1)
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.654212

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  WRITE(*,*) "a_1", PS(6,1)%R_BB(1,1)
  WRITE(*,*) "e  ", PS(6,1)%R_BB(1,2)
  WRITE(*,*) "a_2", PS(6,1)%R_BB(2,2)
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/3.413315
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", 1/0.112842

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(7))
  WRITE(*,*) "a_1", PS(7,1)%R_BB(1,1)
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.721184

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(8))
  WRITE(*,*) "a_1", PS(8,1)%R_BB(1,1)
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.135085




  DEALLOCATE(CHANNELS)
  DEALLOCATE(PS)

END PROGRAM SCATTERING_NN_ZERO_ENERGY