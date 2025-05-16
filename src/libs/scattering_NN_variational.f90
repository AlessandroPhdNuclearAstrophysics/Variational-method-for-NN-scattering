MODULE SCATTERING_NN_VARIATIONAL
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: NNE = 80
  INTEGER, PARAMETER :: NCH_MAX = 2
  INTEGER, PARAMETER :: NNN = NNE * 2
  INTEGER, PARAMETER :: NNR = 100000
  INTEGER :: NCH
  INTEGER :: NDIM ! GIVEN BY CORE_CORE_MATRIX_ELEMENTS
  INTEGER :: NEQ  ! GIVEN BY ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS

  DOUBLE PRECISION, PARAMETER :: HC = 197.327053D0
  DOUBLE PRECISION, PARAMETER :: MP = 938.272029D0
  DOUBLE PRECISION, PARAMETER :: MN = 939.565630D0
  DOUBLE PRECISION, PARAMETER :: MR = MP * MN / (MP + MN)
  DOUBLE PRECISION, PARAMETER :: PI = 4*DATAN(1.0D0)
  COMPLEX*16, PARAMETER ::       IM = (0.0D0, 1.0D0)
  DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0, ZERO = 0.0D0

  DOUBLE PRECISION :: K, HTM, M
  INTEGER :: LC(NCH_MAX), T, T1Z, T2Z

  TYPE, PUBLIC :: VARIATIONAL_PARAMETERS
    INTEGER :: J, L, S, TZ, IPOT, ILB, LEMP
    LOGICAL :: VCOUL
    DOUBLE PRECISION :: E
    DOUBLE PRECISION :: HR1, H, RANGE, GAMMA, EPS, AF
    INTEGER :: NNL
  END TYPE VARIATIONAL_PARAMETERS

  TYPE(VARIATIONAL_PARAMETERS), PRIVATE :: VARIATIONAL_PARAMS

  PUBLIC :: SET_VARIATIONAL_PARAMETERS
  PUBLIC :: NN_SCATTERING_VARIATIONAL
  PRIVATE:: PRINT_DIVIDER

CONTAINS
  SUBROUTINE SET_VARIATIONAL_PARAMETERS(E, J, L, S, TZ, IPOT, ILB, LEMP, VCOUL)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: E
    INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP
    LOGICAL, INTENT(IN) :: VCOUL

    VARIATIONAL_PARAMS%J = J
    VARIATIONAL_PARAMS%L = L
    VARIATIONAL_PARAMS%S = S
    VARIATIONAL_PARAMS%TZ = TZ
    VARIATIONAL_PARAMS%IPOT = IPOT
    VARIATIONAL_PARAMS%ILB = ILB
    VARIATIONAL_PARAMS%LEMP = LEMP
    VARIATIONAL_PARAMS%E = E
    VARIATIONAL_PARAMS%VCOUL = VCOUL
    VARIATIONAL_PARAMS%HR1 = 0.01D0
    VARIATIONAL_PARAMS%H   = 0.02D0
    VARIATIONAL_PARAMS%AF  = 1.02D0
    VARIATIONAL_PARAMS%RANGE = 40.0D0
    VARIATIONAL_PARAMS%GAMMA = 4.D0
    VARIATIONAL_PARAMS%EPS = 0.25D0
    VARIATIONAL_PARAMS%NNL = 32
    
    SELECT CASE (TZ)
      CASE (1)
        M = MP
        T1Z = 1
        T2Z = 1
      CASE (0)
        M = MP * MN / (MP + MN)
        T1Z = 1
        T2Z =-1
      CASE (-1)
        M = MN
        T1Z =-1
        T2Z =-1
      CASE DEFAULT
        PRINT *, "Invalid TZ value"
        STOP
    END SELECT

    HTM = HC**2 / (2 * M)

  ! Ensure T is set such that T + L + S is odd
    IF (MOD(L + S, 2) == 0) THEN
      T = 1
    ELSE
      T = 0
    END IF

    K = DSQRT(2*E*MR) / HC

    PRINT 5
    PRINT 15, "L", "S", "T", "TZ", "J"
    PRINT 5
    LC(1) = L
    NCH = 1
    PRINT 20, LC(1), S, T, TZ, J
    IF (J-L.EQ.1) THEN 
      LC(2) = L + 2
      PRINT 20, LC(2), S, T, TZ, J
      NCH = 2
    ENDIF
    PRINT 5


  5 FORMAT(30("-"))
 15 FORMAT(" ", A5, A5, A5, A5, A5, A5)
 20 FORMAT(" ", I5, I5, I5, I5, I5, I5)
  END SUBROUTINE SET_VARIATIONAL_PARAMETERS
  

  ! This subroutine calculates the NN scattering using a variational method
  ! INPUT PARAMETERS:
  ! E     : Energy in MeV
  ! J     : Total angular momentum
  ! L     : Orbital angular momentum
  ! S     : Spin
  ! TZ    : Isospin projection
  ! IPOT  : Potential type
  ! ILB   : Interaction type
  ! LEMP  : Orbital angular momentum of the last particle
  ! VCOUL: Coulomb potential flag
SUBROUTINE NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP, VCOUL)
  IMPLICIT NONE
! INPUT PARAMETERS
  DOUBLE PRECISION, INTENT(IN) :: E
  INTEGER, INTENT(IN) :: J, L, S, TZ, IPOT, ILB, LEMP
  LOGICAL, INTENT(IN) :: VCOUL

! VARIABLES AND PARAMETERS FOR DGESV
  DOUBLE PRECISION, DIMENSION(NNN, NNN) :: C, CC, CCC
  DOUBLE PRECISION, DIMENSION(NNN, NCH_MAX) :: CAR, CAI
  DOUBLE PRECISION, DIMENSION(NNN)      :: CARR, CAII
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NNN) :: XRCOEFF, XICOEFF
  DOUBLE PRECISION, DIMENSION(NNN) :: IPIV
  INTEGER :: INFO, IAK, IAB

! MATRICES FOR THE VARIATIONAL METHOD
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: ARI, AIR, ARR, AII
  DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: BD1, BD2, BD3, BD4
  DOUBLE PRECISION, DIMENSION(NCH_MAX,NCH_MAX) :: AM, AN, AMM, ANN
  
! INITIALIZE THE VARIATIONAL PARAMETERS
  CALL SET_VARIATIONAL_PARAMETERS(E, J, L, S, TZ, IPOT, ILB, LEMP, VCOUL)

  CALL PRINT_INFO()

! Evaluating the matrix elements
  CALL PRINT_DIVIDER
  CALL ASYMPTOTIC_CORE_MATRIX_ELEMENTS      (CAR, CAI)
  CALL PRINT_DIVIDER
  CALL CORE_CORE_MATRIX_ELEMENTS            (C)
  CALL PRINT_DIVIDER

  DO IAK=1, NEQ
! Preparing the matrix elements for the diagonalization
    CARR =-CAR(:,IAK)
    CAII =-CAI(:,IAK)
    CC  = C
    CCC = C

  ! Solving the eigenvalue problem using LAPACK dgesv
  ! Evaluating for the "c_{n, alpha}" coefficients
    CALL DGESV(NDIM, 1, CC , NNN, IPIV, CARR, NNN, INFO)
    CALL HANDLE_INFO_ERROR()  ! Handle the error after the first DGESV call
    WRITE(*,*) "INFO: ", INFO
    CALL DGESV(NDIM, 1, CCC, NNN, IPIV, CAII, NNN, INFO)
    CALL HANDLE_INFO_ERROR()  ! Handle the error after the second DGESV call
    WRITE(*,*) "INFO: ", INFO

    XRCOEFF(IAK,:) = CARR
    XICOEFF(IAK,:) = CAII
  ENDDO

! Calculating R coefficients
  write(*,*) neq, NDIM
  WRITE(900,'(F20.16)') XICOEFF
  WRITE(901,'(F20.16)') XRCOEFF
  WRITE(902,'(F20.16)') CAR
  WRITE(903,'(F20.16)') CAI
  WRITE(904,'(F50.16)') C


  CALL MATRIX_MULTIPLY(BD1,XICOEFF,CAI)
  CALL MATRIX_MULTIPLY(BD2,XRCOEFF,CAI)
  CALL MATRIX_MULTIPLY(BD3,XICOEFF,CAR)
  CALL MATRIX_MULTIPLY(BD4,XRCOEFF,CAR)

  CALL PRINT_DIVIDER
  CALL ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS(ARI, AIR, ARR, AII)
  CALL PRINT_DIVIDER

  WRITE(*,*) "ARI = ", ARI
  WRITE(*,*) "AIR = ", AIR
  WRITE(*,*) "ARR = ", ARR
  WRITE(*,*) "AII = ", AII
  WRITE(*,*) "BD1 = ", BD1
  WRITE(*,*) "BD2 = ", BD2
  WRITE(*,*) "BD3 = ", BD3
  WRITE(*,*) "BD4 = ", BD4

  DO IAB=1,NEQ       
  DO IAK=1,NEQ
    AM(IAB,IAK)=BD1(IAK,IAB)+BD1(IAB,IAK)+AII(IAB,IAK)+AII(IAK,IAB)
    AN(IAB,IAK)=BD2(IAK,IAB)+BD3(IAB,IAK)+2.D0*AIR(IAB,IAK)   
    WRITE(*,*) "---> ", BD1(IAK,IAB), BD1(IAB,IAK), AII(IAB,IAK), AII(IAK,IAB)
    WRITE(*,*) "---> ", BD2(IAK,IAB), BD3(IAB,IAK), AIR(IAB,IAK)
  ENDDO
  ENDDO

  

  AMM = AM
  ANN =-AN
  CALL PRINT_DIVIDER
  WRITE(*,*) "AMM = ", AMM
  WRITE(*,*) "ANN = ", ANN
  CALL PRINT_DIVIDER


! Evaluating the "R_{alpha, beta}" matrix elements
  CALL DGESV(NEQ, NEQ, AMM, NCH_MAX, IPIV, ANN, NCH_MAX, INFO)
  CALL HANDLE_INFO_ERROR()  ! Handle the error after the third DGESV call
  WRITE(*,*) "INFO: ", INFO

  DO IAB=1,NEQ
  DO IAK=1,NEQ
    WRITE(*,*)"COEFF R",ANN(IAB,IAK)
  ENDDO 
  ENDDO 
  
  STOP


! Evaluating the "R_{alpha, beta}" matrix elements to the second order



! Calculating the phase shifts and mixing angles







  CONTAINS
    SUBROUTINE HANDLE_INFO_ERROR()
      IMPLICIT NONE
      IF (INFO.NE.0) THEN
        PRINT *, "Error in LAPACK dgesv: INFO = ", INFO
        STOP
      ENDIF
    END SUBROUTINE HANDLE_INFO_ERROR

    SUBROUTINE PRINT_INFO()
      IMPLICIT NONE
      PRINT *, "E:    ",                  VARIATIONAL_PARAMS%E, " MeV"
      PRINT *, "HTM:  ",                  HTM, " MeV fm^2"
      PRINT *, "k:    ",                  K, " fm^-1"
      PRINT 10, "J:     ",                VARIATIONAL_PARAMS%J
      PRINT 10, "NCH:   ",                NCH
      PRINT 10, "L0:    ",                LC(1)
      IF (NCH.EQ.2) PRINT 10, "L1:    ",  LC(2)
      PRINT 10, "S:     ",                VARIATIONAL_PARAMS%S
      PRINT 10, "T:     ",                T
      PRINT 10, "TZ:    ",                VARIATIONAL_PARAMS%TZ

      10 FORMAT(" ",A, I2)
    END SUBROUTINE PRINT_INFO

    SUBROUTINE MATRIX_MULTIPLY(M1,M2,M3)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(NCH_MAX,NCH_MAX) :: M1
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NCH_MAX, NNN) :: M2
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NNN, NCH_MAX) :: M3
      INTEGER :: II, JJ, LL

      M1 = 0
      DO II = 1, NEQ
      DO JJ = 1, NEQ
        DO LL = 1, NDIM
          M1(II,JJ) = M1(II,JJ) + M2(II,LL)*M3(LL,JJ)
        ENDDO
      ENDDO
      ENDDO
    END SUBROUTINE MATRIX_MULTIPLY

  END SUBROUTINE NN_SCATTERING_VARIATIONAL



  SUBROUTINE CORE_CORE_MATRIX_ELEMENTS(AM)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NNN, NNN), INTENT(OUT) :: AM
    INTEGER, PARAMETER :: NNR=200
    DOUBLE PRECISION :: XPNT(NNR), PWEIGHT(NNR)
    DOUBLE PRECISION :: XX(NNR), WG(NNR), YY(NNR)
    DOUBLE PRECISION :: U0(0:NNE,NNR), U1(0:NNE,NNR), U2(0:NNE,NNR)
    DOUBLE PRECISION :: V0(NNE,NNR), V1(NNE,NNR), V2(NNE,NNR)
    DOUBLE PRECISION :: GAMMA, APF, XG, ANL, R
    INTEGER :: I, M, NX, NMX
    INTEGER :: L, S, J, NNL
    DOUBLE PRECISION :: V(NNR, NCH_MAX, NCH), VPW(NCH_MAX, NCH_MAX)
    INTEGER :: ICONT(NCH_MAX,NNE), LIK, IAB, IAK, IL, IR, IB, IK
    DOUBLE PRECISION :: SUM, AKEM(NNN,NNN), APEM(NNN,NNN), AXX(NNN,NNN),FUN(NNR)

    GAMMA = VARIATIONAL_PARAMS%GAMMA
    NNL = VARIATIONAL_PARAMS%NNL

    NX = 100
    CALL GAULAG(NX, XPNT,PWEIGHT)

    DO I=1,NX                                                    
      XX(I)=XPNT(I)/GAMMA
      WG(I)=PWEIGHT(I)                                   
      YY(I)=XPNT(I)     ! grid for evaluating Laguerre
    ENDDO  

    WRITE(*,*)'PRIMO E ULTIMO PUNTO =',XX(1),XX(NX)
     
    NMX = VARIATIONAL_PARAMS%NNL-1                                   
    APF=2.D0
    CALL LAGUERRE_POLYNOMIAL(YY, NX, APF, NMX, U0, U1, U2)
    DO M=1, NX
      DO I=0, NMX
        ANL = DSQRT(DGAMMA(I+1.D0)*GAMMA**3/DGAMMA(I+3.D0))

        V0(I+1,M ) = ANL * U0(I,M)
        V1(I+1,M ) = ANL * GAMMA * (U1(I,M) - 0.5D0*U0(I,M))
        V2(I+1,M ) = ANL * GAMMA**2 * (U2(I,M) - 0.5D0*U1(I,M)) &
                          - 0.5D0*GAMMA*V1(I+1,M )
      ENDDO
    ENDDO

    L = LC(1)
    S = VARIATIONAL_PARAMS%S
    J = VARIATIONAL_PARAMS%J

    DO I=1, NX
      R = XX(I)
      CALL AV18PW90(1, L, S, J, T, T1Z, T2Z, R, VPW, VARIATIONAL_PARAMS%LEMP)
      V(I,1,1) = VPW(1,1)
      V(I,1,2) = VPW(1,2)
      V(I,2,1) = VPW(2,1)
      V(I,2,2) = VPW(2,2)
    ENDDO

    CALL PREPARE_INDECES

    DO IAB=1,NEQ          
    DO IAK=1,NEQ
           
      LIK=LC(IAK)*(LC(IAK)+1)

      DO IL=1,NNL            
      DO IR=1,NNL            
              
        IB=ICONT(IAB,IL)       
        IK=ICONT(IAK,IR)       

  ! SI CALCOLA LA NORMA
        AXX(IB,IK)=0.D0
        IF(IB.EQ.IK) AXX(IB,IK) = VARIATIONAL_PARAMS%E
  ! SI CALCOLA ENERGIA CINETICA
        AKEM(IB,IK)=0.D0              
        IF(IAB.EQ.IAK)THEN
          SUM=0.D0
          DO I=1,NX
            FUN(I)=V0(IL,I)*(V2(IR,I)+2.D0*V1(IR,I)/XX(I) &  
                            -LIK*V0(IR,I)/XX(I)**2 )
            SUM=SUM + XX(I)**2*FUN(I)*WG(I) 
          ENDDO                                                                 
          AKEM(IB,IK)=-HTM*SUM/GAMMA                                       
        ENDIF

        IF(IB.EQ.1.AND.IK.EQ.1)THEN
          WRITE(*,*)
          WRITE(*,*)'C-C MATRIX'
          WRITE(*,*)'KINETIC',AKEM(1,1)
        ENDIF

  ! SI CALCOLA ENERGIA POTENZIALE
        SUM=0.D0                                         
        DO I=1,NX
          FUN(I)=V0(IL,I)*V0(IR,I)*V(I,IAB,IAK)     
          SUM=SUM+XX(I)*XX(I)*FUN(I)*WG(I)
        ENDDO
        APEM(IB,IK)=1./GAMMA*SUM                                    

  ! SI CALCOLA HAMILTONIANA    
        AM(IB,IK)=1./HTM*(AKEM(IB,IK)+APEM(IB,IK)-AXX(IB,IK))                                  
        IF(IB.EQ.1.AND.IK.EQ.1)THEN
          WRITE(*,*)'POTENTIAL',APEM(1,1)
          WRITE(*,*)'C-C MATRIX',HTM*AM(IB,IK)
        ENDIF

      ENDDO ! IR
      ENDDO ! IL

    ENDDO ! IAK
    ENDDO ! IAB

  CONTAINS 
    SUBROUTINE PREPARE_INDECES()
      IMPLICIT NONE
      INTEGER II, J

      II = 0
      DO I=1, NEQ
        DO J=1, VARIATIONAL_PARAMS%NNL
          II = II + 1
          ICONT(I, J) = II
        ENDDO
      ENDDO
      NDIM = II
    END SUBROUTINE PREPARE_INDECES
  END SUBROUTINE CORE_CORE_MATRIX_ELEMENTS

!
  SUBROUTINE ASYMPTOTIC_CORE_MATRIX_ELEMENTS(AM, AM1)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NNN, NCH_MAX), INTENT(OUT) :: AM(NNN, NCH_MAX), AM1(NNN, NCH_MAX)

    DOUBLE PRECISION :: H5, HR ! Step size in r
    INTEGER :: I, M, NX, NMX ! Number evenly spaced points
    DOUBLE PRECISION :: XX(NNR), AJ(NNR), YYB(NNR), YYL(NNR), A(NNR)
    DOUBLE PRECISION :: U0(0:NNE, NNR), U1(0:NNE, NNR), U2(0:NNE, NNR) 
    DOUBLE PRECISION :: V0(NNE, NNR), V1(NNE, NNR), V2(NNE, NNR) 
    DOUBLE PRECISION :: FBES(NCH_MAX, NNR), GBES(NCH_MAX, NNR)
    DOUBLE PRECISION :: APF, GAMMA, ANL, XG, FEXP, RR
    INTEGER :: NEQC, L, S, J
    DOUBLE PRECISION :: VPW(2, 2), VV(NNR, NCH_MAX, NCH_MAX)
    INTEGER :: ICONT(NCH_MAX, NNE)
    INTEGER :: IAB, IAK, LIK, IL, IB
    DOUBLE PRECISION :: FUN(NNR), FUN1(NNR)
    DOUBLE PRECISION :: AXXM1(NNN, NCH_MAX), AXX1
    DOUBLE PRECISION :: AKEM(NNN, NCH_MAX), AKE1, AKEM1(NNN, NCH_MAX)
    DOUBLE PRECISION :: APE, APE1, APEM(NNN, NCH_MAX), APEM1(NNN, NCH_MAX)
    DOUBLE PRECISION, EXTERNAL :: B5

    GAMMA = VARIATIONAL_PARAMS%GAMMA

    HR = VARIATIONAL_PARAMS%HR1
    H5 = HR/22.5D0
    NX = INT(VARIATIONAL_PARAMS%RANGE/HR) + 10

    IF (VARIATIONAL_PARAMS%RANGE.LT.H5 .OR. VARIATIONAL_PARAMS%RANGE.GT.200.D0) THEN
      PRINT *, "Error: RANGE out of bounds"
      STOP
    ENDIF


    IF (ABS(NX).GT.NNR) THEN
      PRINT *, "Error: NX exceeds NNR"
      STOP
    ENDIF

  ! Initialize grid with r values
    WRITE(*,*)'NX =',NX
    DO I=1, NX
      XX(I)   = HR*I
      AJ(I)   = XX(I)**2
      YYB(I)  = K*XX(I)
      YYL(I)  = GAMMA*XX(I)
      A(I)    = 1.D0 - DEXP(-VARIATIONAL_PARAMS%EPS*XX(I))
    ENDDO
    WRITE(*,*)'PRIMO E ULTIMO PUNTO =',XX(1),XX(NX),NX

  ! Laguerre polynomial and their first two derivatives grid
    NMX = VARIATIONAL_PARAMS%NNL - 1
    APF = 2.D0

    CALL LAGUERRE_POLYNOMIAL(YYL, NX, APF, NMX, U0, U1, U2)
    DO M = 1, NX
      XG = YYL(M)
      FEXP = DEXP(-XG/2.D0)
      DO I = 0, NMX
        ANL = DSQRT(DGAMMA(I+1.D0)*GAMMA**3/DGAMMA(I+3.D0))*FEXP

        V0(I+1, M) = ANL * U0(I,M)
        V1(I+1, M) = ANL * GAMMA * ( U1(I,M) -0.5D0*U0(I,M) )
        V2(I+1, M) = ANL * GAMMA * ( GAMMA * ( U2(I,M) -0.5D0*U1(I,M) ) ) &
                    - 0.5D0*GAMMA*V1(I+1,M) 
      ENDDO
    ENDDO
    
    NEQ = NCH
    NEQC= NCH
    
  ! Prepare the Bessel functions
    CALL SPHERICAL_BESSEL_FUNCTIONS()
    do m=1, NX
      do i=0, nmx
        Write(1000,*) yyl(M), v0(i+1,m), v1(i+1,m), v2(i+1,m), aj(m)
      enddo
    enddo


    ! DO I = 1, NCH
    !   DO M = 1, NX
    !     WRITE(20+I, *) YYB(M), FBES(I,M), GBES(I,M)
    !   ENDDO
    ! ENDDO
    
  ! Prepare the potential 
    L = LC(1)
    S = VARIATIONAL_PARAMS%S
    J = VARIATIONAL_PARAMS%J

    DO I = 1, NX
      RR = XX(I)
      CALL AV18PW90(1, LC(1), S, J, T, T1Z, T2Z, RR, VPW, VARIATIONAL_PARAMS%LEMP)
      VV(I, 1, 1) = VPW(1, 1)
      VV(I, 1, 2) = VPW(1, 2)
      VV(I, 2, 1) = VPW(2, 1)
      VV(I, 2, 2) = VPW(2, 2)
      WRITE(23, *) XX(I), VV(I, 1, 1), VV(I, 1, 2), VV(I, 2, 1), VV(I, 2, 2)  
    ENDDO

  ! Prepare the indeces for the matrix elements
    CALL PREPARE_INDECES()


    ! Evaluate the matrix elements
    DO IAB = 1, NEQC
      DO IAK = 1, NEQ
        LIK = LC(IAB)*(LC(IAB)+1)

        DO IL = 1, VARIATIONAL_PARAMS%NNL
          IB = ICONT(IAB, IL)

        ! Evaluate the normalization core-irregular (axx1)
          AXXM1(IB,IAK)=0.D0
          IF(IAB.EQ.IAK)THEN 
            FUN1(1)=0.D0                                         
            DO I=1,NX  
            FUN1(I+1)=AJ(I)*V0(IL,I)*GBES(IAK,I)
            ! WRITE(150,*) DSQRT(AJ(I)), V0(IL,I), GBES(IAK,I)
            ENDDO
            ! STOP
            !  WRITE(111,*) FUN1

            AXX1=VARIATIONAL_PARAMS%E * B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)
            AXXM1(IB,IAK)=AXX1
            write(111,*) iab, iak, il, IB, axx1
          ENDIF

        ! Evaluate the kinetic energy core-irregular (ake1)
          AKE1=0.D0
          AKEM1(IB,IAK)=0.D0     
          IF(IAB.EQ.IAK)THEN
             FUN1(1)=0.D0
             DO I=1,NX
              FUN1(I+1)=AJ(I)*GBES(IAK,I)*(V2(IL,I) + 2.D0*V1(IL,I)/XX(I)-LIK*V0(IL,I)/XX(I)**2)
             ENDDO
             AKE1= -HTM * B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)                                
             AKEM1(IB,IAK)=AKE1
            write(112,*) iab, iak, il, IB, ake1
          ENDIF
          IF(IB.EQ.1.AND.IAK.EQ.1)THEN
            WRITE(*,*)
            WRITE(*,*)'C-A MATRIX'
            WRITE(*,*)'IRREGULAR A'
            WRITE(*,*)'NORM ',AXXM1(1,1)
            WRITE(*,*)'KINETIC',AKEM1(1,1)
          ENDIF
        
        
        ! Evaluate the potential energy core-regular (ape), core-irregular (ape1)
          FUN(1)=0.D0
          FUN1(1)=0.D0                                 
          DO I=1,NX   
            FUN(I+1)=AJ(I)*V0(IL,I)*FBES(IAK,I)*VV(I,IAB,IAK)     
            FUN1(I+1)=AJ(I)*V0(IL,I)*GBES(IAK,I)*VV(I,IAB,IAK)                             
          ENDDO
          APE=B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)
          APEM(IB,IAK)=APE
          APE1=B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)
          APEM1(IB,IAK)=APE1
          write(113,*) iab, iak, il, IB, ape, ape1
          IF(IB.EQ.1.AND.IAK.EQ.1)THEN
            WRITE(*,*)
            WRITE(*,*)'C-A MATRIX'
            WRITE(*,*)'IRREGULAR A'
            WRITE(*,*)'POTENTIAL ',APEM1(1,1)
            WRITE(*,*)'REGULAR A'
            WRITE(*,*)'POTENTIAL ',APEM(1,1)
          ENDIF

        ! Evaluate the Hamiltonian: core-regular (am), core-irregular (am1)

          AM(IB,IAK) = 1.D0/HTM * APEM(IB,IAK)
          AM1(IB,IAK)= 1.D0/HTM * (AKEM1(IB,IAK)+APEM1(IB,IAK)-AXXM1(IB,IAK))
          
          IF(IB.EQ.1.AND.IAK.EQ.1)THEN
            WRITE(*,*)
            WRITE(6,*)'C-A MATRIX'     ,IB,IAK
            WRITE(6,*)"CORE-REGULAR=  ",AM(IB,IAK),IB,IAK
            WRITE(6,*)"CORE-IRREGULAR=",AM1(IB,IAK),IB,IAK
          ENDIF
        ENDDO
      ENDDO
    ENDDO 
    RETURN

  CONTAINS
    SUBROUTINE SPHERICAL_BESSEL_FUNCTIONS()
      use gsl_bessel
      IMPLICIT NONE
      
      DOUBLE PRECISION :: AG, GBSS, FBSS
      INTEGER :: L
      
      DO I=1,NCH
        DO M=1, NX                                 
          XG=YYB(M)
          AG=A(M)
          L=LC(I)
          IF(K.LE.1.D-8)THEN                                   !(K->0)
            FBES(I,M)=XX(M)**L
            GBES(I,M)=-1./((2*L+1.)*XX(M)**(L+1.D0))*AG**(2*L+1.D0)
          ELSE  
            FBSS = SPHERICAL_J(L, XG)
            GBSS = SPHERICAL_Y(L, XG)
            FBES(I,M)=K**(L+0.5D0)*FBSS/(K**L)
            GBES(I,M)=-(GBSS*K**(L+1.D0)*AG**(2*L+1.D0))/(K**(L+0.5D0))
            write(500+l,*) xg, FBSS, GBSS, fbes(i,m), gbes(i,m)
          ENDIF
        ENDDO
      ENDDO
    END SUBROUTINE SPHERICAL_BESSEL_FUNCTIONS

    SUBROUTINE PREPARE_INDECES()
      IMPLICIT NONE
      INTEGER II, JJ

      II = 0
      DO I=1, NEQC
        DO JJ=1, VARIATIONAL_PARAMS%NNL
          II = II + 1
          ICONT(I, JJ) = II
        ENDDO
      ENDDO
    END SUBROUTINE PREPARE_INDECES

  END SUBROUTINE ASYMPTOTIC_CORE_MATRIX_ELEMENTS




  SUBROUTINE ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS(AM, AM1, AM2, AM3)
    USE gsl_bessel
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NCH, NCH), INTENT(OUT) :: AM, AM1, AM2, AM3
    INTEGER, PARAMETER :: NNR = 200
    DOUBLE PRECISION :: H, H5, RANGE, R
    DOUBLE PRECISION :: AF, EPS, XX(NNR), AJ(NNR), YY(NNR), A(NNR), B(NNR)
    DOUBLE PRECISION, DIMENSION(NCH_MAX,NNR) :: GBES, GBES0, GBES1, GBES2, FBES, HNOR
    INTEGER :: I, NX, M, L, IAB, IAK
    DOUBLE PRECISION :: XG, AG, BG
    DOUBLE PRECISION :: GBSS, GBSS1, GBSS2, FBSS
    INTEGER :: S, J
    DOUBLE PRECISION :: VPW(NCH_MAX,NCH_MAX), VV(NNR,NCH_MAX,NCH_MAX)
    DOUBLE PRECISION, DIMENSION(NNR) :: FUN, FUN1, FUN2, FUN3
    DOUBLE PRECISION :: AXX, AXX3, AKE, AKE3, APE, APE1, APE2, APE3
    DOUBLE PRECISION, DIMENSION(NCH_MAX, NCH_MAX) :: AXXM, AXXM3, AKEM, AKEM3, APEM, APEM1, APEM2, APEM3, VER

    DOUBLE PRECISION, EXTERNAL :: B5

    H = VARIATIONAL_PARAMS%H
    H5= H/22.5D0
    RANGE = VARIATIONAL_PARAMS%RANGE
    AF = VARIATIONAL_PARAMS%AF
    CALL EXPONENTIALLY_GROWING_GRID(H, AF, RANGE, XX, AJ, NNR, NX)
    VARIATIONAL_PARAMS%RANGE = RANGE
    WRITE(*,*)'PRIMO E ULTIMO PUNTO =',XX(1),XX(NX)

    EPS = VARIATIONAL_PARAMS%EPS
    DO I=1, NX
      R = XX(I)
      YY(I) = K*R
      A(I)  = 1.0 - DEXP(-EPS*R)
      B(I)  = EPS * DEXP(-EPS*R)
    ENDDO

  !definisco funzioni bessel regolarizzate e con giuste dimensioni (K**(l+0.5d0) ed andamenti asintotici  
    DO I=1, NEQ
      L=LC(I)
      write(*,*) "Preparing Bessels for L = ", L, " with k = ", K
      DO M=1, NX                       
        XG=YY(M)
        AG=A(M)
        BG=B(M)
        IF(K.LE.1.D-8) THEN
          FBES(I,M)=XX(M)**L
          GBES(I,M)=-1./((2*L+1.)*XX(M)**(L+1.D0))*AG**(2*L+1.D0)
          GBES0(I,M)=1./((2*L+1.)*XX(M)**(L+1.D0))*(EPS*BG*(2*L+1.)*((2*L+1.)*(BG/EPS)-1.) &
                    +2*(2*L+1.)*BG*(AG/XX(M))-L*(L+1.)*(AG/XX(M))**2)
          GBES1(I,M)=-2.*AG*((2*L+1.)*BG+AG/XX(M))*(L+1.)/((2*L+1.)*XX(M)**(L+2.))
          GBES2(I,M)=AG**2*(L+1.)*(L+2.)/((2*L+1.)*XX(M)**(L+3.))
          HNOR(I,M)=AG**(2*L-1.)
          WRITE(200+L,*) XX(M), FBES(I,M), GBES(I,M), GBES1(I,M), GBES2(I,M)
        ELSE
          FBSS = SPHERICAL_J(L, XG)
          GBSS = SPHERICAL_Y(L, XG)
          GBSS1= SPHERICAL_YP(L, XG)
          GBSS2= SPHERICAL_YPP(L, XG)
          
          FBES(I,M)=K**(L+0.5D0)*FBSS/(K**L)
          GBES(I,M)=-(GBSS*K**(L+1.D0)*AG**(2*L+1.D0))/(K**(L+0.5D0))
          GBES0(I,M)=GBSS*(EPS*BG*(2*L+1.)*((2*L+1.)*(BG/EPS)-1.) &
                    +2*(2*L+1.)*BG*(AG/XX(M))-L*(L+1.)*(AG/XX(M))**2)
          GBES1(I,M)=GBSS1*2.*K*AG*((2*L+1.)*BG + AG/XX(M))                   
          GBES2(I,M)=(K**2)*(AG**2)*GBSS2
          HNOR(I,M)=(K**(L+1.D0))*(AG**(2*L-1.))/(K**(L+0.5D0))
          ! WRITE(200+L,*) XX(M), FBES(I,M), GBES(I,M), GBES1(I,M), GBES2(I,M), HNOR(I,M)
        ENDIF
      ENDDO
    ENDDO

  ! si tabula il potenziale per i vari canali
    L = LC(1)
    S = VARIATIONAL_PARAMS%S
    J = VARIATIONAL_PARAMS%J
    DO I=1, NX
      R = XX(I)
      CALL AV18PW90(1,L,S,J,T,T1Z,T2Z,R,VPW,VARIATIONAL_PARAMS%LEMP)
      VV(I,1,1) = VPW(1,1)
      VV(I,1,2) = VPW(1,2)
      VV(I,2,1) = VPW(2,1)
      VV(I,2,2) = VPW(2,2)
    ENDDO

  !SI CALCOLANO ELEMENTI MATRICE  
    DO IAB=1,NEQ          
    DO IAK=1,NEQ
  ! SI CALCOLA NORMA DEL CASO REGOLARE-IRREGOLARE (AXX), CASO IRREGOLARE-IRREGOLARE (AXX3)  
      AXX=0.D0
      AXX3=0.D0
      AXXM(IAB,IAK)=0.D0
      AXXM3(IAB,IAK)=0.D0
      IF(IAB.EQ.IAK)THEN 
          FUN(1)=0.D0
          FUN3(1)=0.D0                                     
        DO I=1,NX
          FUN(I+1)= AJ(I)*FBES(IAB,I)*GBES(IAK,I)
          FUN3(I+1)=AJ(I)*GBES(IAB,I)*GBES(IAK,I)                                             
        ENDDO
        AXX=  VARIATIONAL_PARAMS%E * B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)
        AXXM(IAB,IAK)=AXX
        AXX3= VARIATIONAL_PARAMS%E * B5(1,NX,0,0,H5,0.D0,0.D0,FUN3,1)
        AXXM3(IAB,IAK)=AXX3
      ENDIF    
      IF (IAB==1 .AND. IAK==1) THEN
        WRITE(*,*) 'NORM RI(1,1)', AXX
        WRITE(*,*) 'NORM II(1,1)', AXX3
      ENDIF

  ! SI CALCOLA ENERGIA CINETICA DEL CASO REGOLARE-IRREGOLARE (AKE) E CASO IRREGOLARE-IRREGOLARE (AKE3)
      AKE=0.D0
      AKE3=0.D0
      AKEM(IAB,IAK)=0.D0 
      AKEM3(IAB,IAK)=0.D0       
      IF(IAB.EQ.IAK)THEN
        FUN(1)=0.D0
        FUN3(1)=0.D0
        DO I=1,NX
          FUN(I+1) =AJ(I)*FBES(IAB,I)*HNOR(IAK,I)*(GBES2(IAK,I)+GBES1(IAK,I)+GBES0(IAK,I))
          FUN3(I+1)=AJ(I)*GBES(IAB,I)*HNOR(IAK,I)*(GBES2(IAK,I)+GBES1(IAK,I)+GBES0(IAK,I))
        ENDDO
        AKE=  HTM * B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)                                         
        AKEM(IAB,IAK)=AKE
        AKE3= HTM * B5(1,NX,0,0,H5,0.D0,0.D0,FUN3,1)
        AKEM3(IAB,IAK)=AKE3
      ENDIF 
      IF (IAB==1 .AND. IAK==1) THEN
        WRITE(*,*) 'KINETIC RI(1,1)', AKE
        WRITE(*,*) 'KINETIC II(1,1)', AKE3
      ENDIF

  ! SI CALCOLA ENERGIA POTENZIALE DEL CASO REGOLARE-IRREGOLARE (APE), 
  !      IRREGOLARE-REGOLARE (APE1), REGOLARE-REGOLARE (APE2), IRREGOLARE-IRREGOLARE (APE3)
      APE=0.D0
      APE1=0.D0
      APE2=0.D0
      APE3=0.D0
      FUN(1)=0.D0 
      FUN1(1)=0.D0
      FUN2(1)=0.D0 
      FUN3(1)=0.D0                                          
      DO I=1,NX
        FUN(I+1) =AJ(I)*FBES(IAB,I)*GBES(IAK,I)*VV(I,IAB,IAK) 
        FUN1(I+1)=AJ(I)*GBES(IAB,I)*FBES(IAK,I)*VV(I,IAB,IAK)  
        FUN2(I+1)=AJ(I)*FBES(IAB,I)*FBES(IAK,I)*VV(I,IAB,IAK)
        FUN3(I+1)=AJ(I)*GBES(IAB,I)*GBES(IAK,I)*VV(I,IAB,IAK)
      ENDDO
      APE=  B5(1,NX,0,0,H5,0.D0,0.D0,FUN,1)
      APEM(IAB,IAK)=APE
      APE1= B5(1,NX,0,0,H5,0.D0,0.D0,FUN1,1)
      APEM1(IAB,IAK)=APE1
      APE2= B5(1,NX,0,0,H5,0.D0,0.D0,FUN2,1)
      APEM2(IAB,IAK)=APE2
      APE3= B5(1,NX,0,0,H5,0.D0,0.D0,FUN3,1)
      APEM3(IAB,IAK)=APE3
      IF (IAB==1 .AND. IAK==1) THEN
        WRITE(*,*) 'POTENTIAL RR(1,1)', APE2
        WRITE(*,*) 'POTENTIAL RI(1,1)', APE
        WRITE(*,*) 'POTENTIAL IR(1,1)', APE1
        WRITE(*,*) 'POTENTIAL II(1,1)', APE3
      ENDIF



  ! SI CALCOLA HAMILTONIANA PER I VARI CASI:REGOLARE-IRREGOLARE(AM), IRREGOLARE-REGOLARE(AM1),
  !         REGOLARE-REGOLARE(AM2),IRREGOLARE-IRREGOLARE(AM3)
      AM(IAB,IAK)= 1.D0/HTM*(AKEM(IAB,IAK)+APEM(IAB,IAK)-AXXM(IAB,IAK)) 
      AM1(IAB,IAK)=1.D0/HTM*APEM1(IAB,IAK)
      AM2(IAB,IAK)=1.D0/HTM*APEM2(IAB,IAK)
      AM3(IAB,IAK)=1.D0/HTM*(AKEM3(IAB,IAK)+APEM3(IAB,IAK)-AXXM3(IAB,IAK))

      VER(IAB,IAK)=(AM1(IAB,IAK)-AM(IAB,IAK)) 
    ENDDO 
    ENDDO 

    write(*,*)
    WRITE(*,*)'A-A MATRIX'      
    DO IAB=1,NEQ
    DO IAK=1,NEQ
      WRITE(6,'(2I6,4D17.7)') IAB,IAK,AM2(IAB,IAK),AM(IAB,IAK),AM1(IAB,IAK),AM3(IAB,IAK)
      WRITE(*,*)"1=",VER(IAB,IAK)
    END DO
    END DO
  END SUBROUTINE ASYMPTOTIC_ASYMPTOTIC_MATRIX_ELEMENTS

  SUBROUTINE PRINT_DIVIDER()
    IMPLICIT NONE
    WRITE(*,*) '====================================================================================='
  END SUBROUTINE PRINT_DIVIDER

END MODULE SCATTERING_NN_VARIATIONAL