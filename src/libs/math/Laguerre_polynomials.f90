SUBROUTINE LAGUERRE_POLYNOMIAL(XX,NX,APF,N1,U,U1,U2)
!*calcula el polinomio de laguerre en las absisas xx(1,...nx) con numeros
!*cuanticos apf,n con n=0,....n1
  IMPLICIT NONE
  INTEGER, PARAMETER :: NNR=200,NNE=80
  INTEGER, INTENT(IN) :: NX
  INTEGER, INTENT(INOUT) :: N1
  DOUBLE PRECISION, INTENT(IN) :: XX(NNR), APF
  DOUBLE PRECISION, INTENT(OUT) :: U(0:NNE,NX), U1(0:NNE,NX), U2(0:NNE,NX)

  INTEGER :: I, N, NMAX
  DOUBLE PRECISION :: X, D, A0, A1, A2
  NMAX = N1 - 1
  IF (NMAX.GT.NNE-1) THEN
      WRITE(*,*) 'LAGUERRE POLYNOMIAL: N1 TOO LARGE', N1, NNE
      STOP
  ENDIF
  DO I=1,NX
    X=XX(I)
    U (0,I)= 1.D0
    U (1,I)= APF+1.D0-X
    U1(0,I)= 0.D0
    U1(1,I)=-1.D0
    U2(0,I)= 0.D0
    U2(1,I)= 0.D0
  ENDDO
  DO N=1,NMAX
    N1=N+1
    D=1.D0/(N+1.D0)
    A1=2.D0*N+APF+1.D0
    A0=N+APF
    A2=N1+APF
    DO I=1,NX
      X=XX(I)
      U (N1,I)=((A1-X)*U(N,I)-A0*U(N-1,I))*D
      U1(N1,I)=(N1*U(N1,I)-A2*U(N,I))/X
      U2(N1,I)=(-U(1,I)*U1(N1,I)-N1*U(N1,I))/X
    ENDDO
  ENDDO
END SUBROUTINE LAGUERRE_POLYNOMIAL