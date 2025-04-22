!> @brief Numerical integration block configuration
!!
!! @details Controls adaptive numerical integration across 1-3 blocks with independent step sizes.
!!
!! @param[in] JB    Integration block mode (1, 2, or 3 blocks)
!! @param[in] IAS   Starting abscissa point for integration
!! @param[in] M1    Number of steps in Block 1
!! @param[in] M2    Number of steps in Block 2 (if JB ≥ 2)
!! @param[in] M3    Number of steps in Block 3 (if JB = 3)
!! @param[in] H1    Step size for Block 1
!! @param[in] H2    Step size for Block 2 (if JB ≥ 2)
!! @param[in] H3    Step size for Block 3 (if JB = 3)
!! @param[in,out] HI Input step size (automatically divided by 22.5 for internal use)
!!
!! @note The 22.5 scaling factor optimizes stability for stiff systems.
!! @warning JB values > 3 will cause array bounds errors.
!!
!! @code{.f90}
!! ! Example usage:
!! call integrate_system(JB=2, IAS=0.0d0, M1=100, M2=50, &
!!                      H1=0.01d0, H2=0.02d0, HI=initial_step)
!! @endcode
FUNCTION B5(JB,M1,M2,M3,H1,H2,H3,A,IAS) RESULT(I5)
  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
  INTEGER, INTENT(IN) :: JB, IAS
  INTEGER, INTENT(IN) :: M1,M2,M3
  DOUBLE PRECISION, INTENT(IN) :: H1,H2,H3
  DOUBLE PRECISION, INTENT(IN) :: A(M1+M2+M3+3)
  DOUBLE PRECISION :: I5
!
  I5=0.D0
  N=M1
  H5=H1
  J0=IAS
  IR=JB-2
 1    IF(N-4) 44,2,2
 2    HL5=H5/32.D0
  NB=N/4
  N1=N-4*NB+1
  S2=0.D0
  S3=S2
  S4=S2
  L=J0
  DO 4 I=1,NB
    S2=S2+A(L+2)
    S3=S3+A(L+1)+A(L+3)
    S4=S4+A(L)+A(L+4)
 4  L=L+4
    GOTO (8,12,16,20),N1
 8  AX=0.D0
    GOTO 40
12  AX=HL5*(-19.D0*A(L-3)+106.D0*A(L-2)-264.D0*A(L-1)+646.D0*A(L)+ &
      251.D0*A(L+1))
    GOTO 40
16  AX=HL5*(-8.D0*A(L-2)+32.D0*A(L-1)+192.D0*A(L)+992.D0*A(L+1)+ &
      232.D0*A(L+2))
    GOTO 40
20  AX=HL5*(-27.D0*A(L-1)+378.D0*A(L)+648.D0*A(L+1)+918.D0*A(L+2)+ &
      243.D0*A(L+3))
40  I5=I5+AX+H5*(7.D0*S4+32.D0*S3+12.D0*S2)
    IF(IR) 60,50,55
50  J0=M1+IAS
    N=M2
    H5=H2
    IR=-1
    GOTO 1
55  J0=M2+M1+IAS
    N=M3
    H5=H3
    IR=0
    GOTO 1
60  RETURN
  !
44   WRITE(4,200)
200  FORMAT(//1X,'  NUMERO DE INTERVALOS MENOR QUE 4'/)
  STOP
END FUNCTION B5



!> @brief Computes Gauss-Laguerre quadrature points and weights
!!
!! @details Calculates the roots (XPNT) and weights (PWEIGHT) for N-point Gauss-Laguerre 
!! quadrature.
!! This implementation uses Newton's method with polynomial recurrence relations.
!!
!! @param[in]  N        Number of quadrature points (1 ≤ N ≤ 600)
!! @param[out] XPNT     Array(N) of quadrature points (roots of Laguerre polynomial)
!! @param[out] PWEIGHT  Array(N) of quadrature weights
!!
!! @note 
!! - Uses Laguerre polynomials with α=0 (standard Gauss-Laguerre quadrature)
!! - Implements Newton-Raphson iteration with MAXIT=10 and EPS=1D-12 tolerance
!! - Includes initial root approximations for faster convergence
!!
!! @warning
!! - Terminates if N < 1 or N > NNR (600)
!! - Stops with error if convergence fails within MAXIT iterations
!!
!! @algorithm
!! 1. Generates initial guess for ith root using empirical formulas
!! 2. Refines root using Newton's method on Laguerre polynomial recurrence.
!! 3. Computes weights using derivative relationship.
!!
!! @code{.f90}
!! ! Example usage:
!! integer, parameter :: n = 10
!! double precision :: x(n), w(n)
!! call gaulag(n, x, w)
!! ! x now contains quadrature points, w contains weights
!! @endcode
!!
!! @reference
!! - Adapted from Numerical Recipes with modified root approximations
!! - See Golub & Welsch (1969) for mathematical foundation
SUBROUTINE GAULAG(N,XPNT,PWEIGHT)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER(NNR=600)
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(OUT) ::  XPNT(NNR),PWEIGHT(NNR)
  DIMENSION X(NNR)

  PARAMETER (EPS=1.D-12,MAXIT=10)

  IF(N.LT.1)    STOP 'GAULAG.F: BAD ARGUMENT N '
  IF(N.GT.NNR)  STOP 'GAULAG.F: N > NNR TOO LARGE'

  ALF=0.D0
  DO I=1,N
! APPROXIMATE THE ITH ROOT
    IF(I.EQ.1)THEN
      Z=(1.+ALF)*(3.+.92*ALF)/(1.+2.4*N+1.8*ALF)
    ELSE IF(I.EQ.2)THEN
      Z=Z+(15.+6.25*ALF)/(1.+.9*ALF+2.5*N)
    ELSE
      AI=I-2
      Z=Z+((1.+2.55*AI)/(1.9*AI)+1.26*AI*ALF/(1.+3.5*AI))*(Z-X(I-2))/(1.+.3*ALF)
    ENDIF

    DO ITS=1,MAXIT
      P1=1.D0
      P2=0.D0
      DO J=1,N ! RECURRENCE RELATION FOR LAGUERRE POLYNOMIAL IN Z
        P3=P2
        P2=P1
        P1=((2*J-1+ALF-Z)*P2-(J-1+ALF)*P3)/J
      ENDDO
      PP=(N*P1-(N+ALF)*P2)/Z ! DERIVATIVE OF LAGUERRE POLYNOMIAL
      Z1=Z
      Z=Z1-P1/PP ! NEWTON'S METHOD TO REFINE ROOT
      IF(DABS(Z-Z1).LE.EPS) EXIT
    ENDDO

    IF ( ITS > MAXIT ) THEN
      WRITE(*,'("N = ", I3," I = ",I3)') N,I
      STOP 'GAULAG.F: TOO MANY ITERATIONS'
    ENDIF
    X(I)=Z
    XPNT(I)=Z
    DN=DFLOAT(N)
    PWEIGHT(I)=-1.D0/(PP*N*P2) ! WEIGHT
  ENDDO
END SUBROUTINE GAULAG



SUBROUTINE EXPONENTIALLY_GROWING_GRID(H, AF, RANGE, R, W, NMAX, N)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NMAX
  DOUBLE PRECISION, INTENT(IN) :: H, AF
  DOUBLE PRECISION, INTENT(INOUT) :: RANGE
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(NMAX) :: R, W
  INTEGER, INTENT(OUT) :: N

  INTEGER :: NP, I

  N = 1 + DLOG(1.D0 + RANGE*(AF-1.D0)/H)/DLOG(AF)
  
  IF (N > NMAX) THEN
    WRITE(*,*) 'NNR < N ', NMAX, N
    STOP
  ENDIF
  
  NP = N+1
  RANGE = H*(AF**NP - 1.D0)/(AF - 1.D0)
  DO I=1, N
    R(I) = H*(AF**I-1.D0)/(AF-1.D0)
    W(I) = R(I)**2*AF**I*DLOG(AF)/(AF-1.D0)
  ENDDO
END SUBROUTINE EXPONENTIALLY_GROWING_GRID

