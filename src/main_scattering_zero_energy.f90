PROGRAM SCATTERING_NN_ZERO_ENERGY
  USE SCATTERING_NN_VARIATIONAL
  USE QUANTUM_NUMBERS
  USE SCATTERING
  USE EFT_PLESS
  IMPLICIT NONE
  INTEGER, PARAMETER :: LMAX=2, JMAX=2
  INTEGER, PARAMETER :: TZ=0
  TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
  TYPE(PHASE_SHIFT_RESULT), ALLOCATABLE :: PS(:,:)
  TYPE(ZERO_ENERGY_OBSERVABLES) :: OBS
  TYPE(LECS_EFT_PLESS) :: LECS
  DOUBLE PRECISION :: ENERGIES(1), R(2,2)
  INTEGER :: IPOT, ILB, LEMP=0
  INTEGER :: L, DIM

  CALL PREPARE_CHANNELS(LMAX, JMAX, TZ, CHANNELS)
  ALLOCATE(PS(SIZE(CHANNELS), 1))

  WRITE(*,*) 
  WRITE(*,*) "===================== AV18 =========================="
  WRITE(*,*) 
  IPOT = 18
  ILB = 1
  
  ENERGIES(1) = 0.D0
  CALL SET_ENERGIES(ENERGIES)
  CALL SET_CHANNELS(CHANNELS)
  
  CALL NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PS, IPOT=IPOT, ILB=ILB)
  
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"SINGLE CHANNELS RESULTS"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(1))
  R = PS(1,1)%R_BB
  L = CHANNELS(1)%L(1)
  DIM = CHANNELS(1)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.0433

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(3))
  R = PS(3,1)%R_BB
  L = CHANNELS(3)%L(1)
  DIM = CHANNELS(3)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.359170

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(4))
  R = PS(4,1)%R_BB
  L = CHANNELS(4)%L(1)
  DIM = CHANNELS(4)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.397155

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(5))
  R = PS(5,1)%R_BB
  L = CHANNELS(5)%L(1)
  DIM = CHANNELS(5)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.654212


  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(7))
  R = PS(7,1)%R_BB
  L = CHANNELS(7)%L(1)
  DIM = CHANNELS(7)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.721184

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(8))
  R = PS(8,1)%R_BB
  L = CHANNELS(8)%L(1)
  DIM = CHANNELS(8)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.135085
  

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"BLATT-BIEDERMANN PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
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
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
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
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
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
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/3.413315
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 1.61436 
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", -1/1.03225


  CALL RESET_SCATTERING_NN_VARIATIONAL
  WRITE(*,*) 
  WRITE(*,*) "===================== EFT_pless_15 =========================="
  WRITE(*,*) 
  IPOT = 19
  ILB = 15
  
  CALL SET_ENERGIES(ENERGIES)
  CALL SET_CHANNELS(CHANNELS)
  
  CALL NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PS, IPOT=IPOT, ILB=ILB)
  
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"SINGLE CHANNELS RESULTS"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(1))
  R = PS(1,1)%R_BB
  L = CHANNELS(1)%L(1)
  DIM = CHANNELS(1)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.043617

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(3))
  R = PS(3,1)%R_BB
  L = CHANNELS(3)%L(1)
  DIM = CHANNELS(3)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.581388

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(4))
  R = PS(4,1)%R_BB
  L = CHANNELS(4)%L(1)
  DIM = CHANNELS(4)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.438454

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(5))
  R = PS(5,1)%R_BB
  L = CHANNELS(5)%L(1)
  DIM = CHANNELS(5)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.717147


  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(7))
  R = PS(7,1)%R_BB
  L = CHANNELS(7)%L(1)
  DIM = CHANNELS(7)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/2.795050

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(8))
  R = PS(8,1)%R_BB
  L = CHANNELS(8)%L(1)
  DIM = CHANNELS(8)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/1.032862

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"BLATT-BIEDERMANN PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.185051
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 0.284509
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.829435
  
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/2.528886
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", -1.63131
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", 1/1.02352

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"STAPP PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.185051
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 1.53556
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.608997

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/2.528886
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 0.644972
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", -1/14.548




  CALL RESET_SCATTERING_NN_VARIATIONAL
  WRITE(*,*) 
  WRITE(*,*) "===================== EFT_pless_10 =========================="
  WRITE(*,*) 
  IPOT = 19
  ILB = 10
  
  CALL SET_ENERGIES(ENERGIES)
  CALL SET_CHANNELS(CHANNELS)
  
  CALL NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PS, IPOT=IPOT, ILB=ILB)
  
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"SINGLE CHANNELS RESULTS"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(1))
  R = PS(1,1)%R_BB
  L = CHANNELS(1)%L(1)
  DIM = CHANNELS(1)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.043585

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(3))
  R = PS(3,1)%R_BB
  L = CHANNELS(3)%L(1)
  DIM = CHANNELS(3)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.583100

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(4))
  R = PS(4,1)%R_BB
  L = CHANNELS(4)%L(1)
  DIM = CHANNELS(4)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.466579

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(5))
  R = PS(5,1)%R_BB
  L = CHANNELS(5)%L(1)
  DIM = CHANNELS(5)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.749634


  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(7))
  R = PS(7,1)%R_BB
  L = CHANNELS(7)%L(1)
  DIM = CHANNELS(7)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/2.664854

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(8))
  R = PS(8,1)%R_BB
  L = CHANNELS(8)%L(1)
  DIM = CHANNELS(8)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/1.124504


  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"BLATT-BIEDERMANN PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.1849
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 0.285121
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.878179

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/2.825603
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", -1.66873
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", 1/1.08302

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"STAPP PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.1849
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 1.5402 
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.633951

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/2.825603
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 0.590358
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", -1/14.1517




  CALL RESET_SCATTERING_NN_VARIATIONAL
  WRITE(*,*) 
  WRITE(*,*) "===================== EFT_pless_10 DYNAMIC ==================="
  WRITE(*,*) 
  
  LECS = GET_LECS(10)
  CALL SET_DYNAMIC(.TRUE.)
  CALL SET_ENERGIES(ENERGIES)
  CALL SET_CHANNELS(CHANNELS)
  
  CALL NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS(ENERGIES, CHANNELS, LEMP, PS, LECS_FOR_PLESS=LECS)
  
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"SINGLE CHANNELS RESULTS"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(1))
  R = PS(1,1)%R_BB
  L = CHANNELS(1)%L(1)
  DIM = CHANNELS(1)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.043585

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(3))
  R = PS(3,1)%R_BB
  L = CHANNELS(3)%L(1)
  DIM = CHANNELS(3)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.583100

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(4))
  R = PS(4,1)%R_BB
  L = CHANNELS(4)%L(1)
  DIM = CHANNELS(4)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/0.466579

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(5))
  R = PS(5,1)%R_BB
  L = CHANNELS(5)%L(1)
  DIM = CHANNELS(5)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", 1/0.749634


  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(7))
  R = PS(7,1)%R_BB
  L = CHANNELS(7)%L(1)
  DIM = CHANNELS(7)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/2.664854

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(8))
  R = PS(8,1)%R_BB
  L = CHANNELS(8)%L(1)
  DIM = CHANNELS(8)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,'(A,F10.3)') " Expected scattering length: ", -1/1.124504


  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"BLATT-BIEDERMANN PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.1849
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 0.285121
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.878179

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_BB(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/2.825603
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", -1.66873
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", 1/1.08302

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) CHAR(27)//'[34m'//"STAPP PARAMETRIZATION"//CHAR(27)//'[0m'
  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(2))
  R = PS(2,1)%R_BB
  L = CHANNELS(2)%L(1)
  DIM = CHANNELS(2)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=0: ", 1/0.1849
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 1.5402 
  WRITE(*,'(A,F10.3)') " Expected scattering length L=2: ", 1/0.633951

  WRITE(*,*)
  WRITE(*,*) "Scattering lengths at zero energy for channel: ", GET_CHANNEL_NAME(CHANNELS(6))
  R = PS(6,1)%R_BB
  L = CHANNELS(6)%L(1)
  DIM = CHANNELS(6)%NCH()
  OBS = EVALUATE_ZERO_ENERGIES_OBSERVABLES_STAPP(R, DIM, L)
  WRITE(*,*) "a_1", OBS%a1
  WRITE(*,*) "e  ", OBS%e
  WRITE(*,*) "a_2", OBS%a2
  WRITE(*,'(A,F10.3)') " Expected scattering length L=1: ", -1/2.825603
  WRITE(*,'(A,F10.3)') " Expected mixing angle constant: ", 0.590358
  WRITE(*,'(A,F10.3)') " Expected scattering length L=3: ", -1/14.1517


  DEALLOCATE(CHANNELS)
  DEALLOCATE(PS)

END PROGRAM SCATTERING_NN_ZERO_ENERGY