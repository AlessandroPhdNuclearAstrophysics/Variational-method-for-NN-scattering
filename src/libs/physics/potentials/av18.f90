!> **********************************************************************
!! \file av18.f90
!! Argonne v18 and vn' and Super-Soft Core (C) potential package
!! \defgroup av18 Argonne v18
!! \ingroup nn_potentials
!! \brief Module for the Argonne v18 nuclear potential.
!!
!! Prepared 1 Sep 1994 by R. B. Wiringa, Physics Division,  
!! Argonne National Laboratory, Argonne, IL 60439  
!! e-mail: wiringa@theory.phy.anl.gov
!!
!! This modular implementation has been prepared by A. Grassi (May 2025)   
!! e-mail: alessandro.grassi@df.unipi.it
!!
!! \par References:
!! - "Accurate nucleon-nucleon potential with charge-independence breaking",  R. B. Wiringa, V. G. J. Stoks, and R. Schiavilla,  Physical Review C51, 38 (1995) - WSS95.
!! - "Quantum Monte Carlo calculations of nuclei with A<=7",  B. S. Pudliner, V. R. Pandharipande, J. Carlson, Steven C. Pieper,  and R. B. Wiringa, Physical Review C56, 1720 (1997) - PPCPW97
!! - "Evolution of Nuclear Spectra with Nuclear Forces",  R. B. Wiringa and Steven C. Pieper, Physical Review Letters 89, 182501 (2002) - WP02
!! - "Construction d'un potentiel nucleon-nucleon a coeur tres mou (SSC)",  R. de Tourreil and D. W. L. Sprung, 
!>   Nuclear Physics A201, 193 (1973) - TS73
!!
!! \par Options:
!! - v8' reprojection of v18 potential (added 10 Jan 2001)
!! - v6', v4', vx', v2', v1' potentials (added 16 Jul 2002)
!! - Super-Soft Core (C) (added 14 Feb 2007)
!! - Modified Argonne v8' and modified SSCC v8' (added 14 Feb 2007, corrected 4 Apr 2007)
!!
!! \par Contents:
!! This module contains 4 subroutines:
!!   - AV18PW : full potential in a particular partial wave
!!   - AV18OP : strong interaction part in operator format
!!   - EMPOT  : electromagnetic part in operator format
!!   - CONSTS : values of fundamental constants and masses used
!!
!! The variable \c lpot selects between v18, v8' and other options.
!!
!! \par Notes:
!! - av18pw90 includes full EM interaction for lpot=1; for lpot>1 it includes only C1(pp), i.e., Coulomb with a form factor for pp channels.
!! - empot90 does not include the energy-dependence of the Coulomb interaction used in eq.(4) of WSS95, i.e., it uses alpha, not alpha'.
!! - The vacuum polarization in empot90 is a short-range approximation to eq.(7) suitable for bound states, but not for scattering. It is taken from eq.(3.13) of Rev. Mod. Phys. 44, 48 (1972).  
!!    (8/28/97: error in this formula detected and corrected: should be -(gamma+5/6) instead of printed (-gamma+5/6))
!! - These subroutines should be compiled with a compiler option that forces all floating point constants to be evaluated at DOUBLE PRECISION significance.  
!!    For example: on IBM RS6000 use xlf option qdpc=e; on SGI use -r8; on Cray no action is needed.  
!!    If such an option is not available and the default precision is real*4 (32 bits), then all constants should be explicitly converted to double precision by appending a D0.
!! - consts90 now (14 Feb 2007) depends upon potential: need to call to generate appropriate hbar**2/M.
MODULE AV18
  IMPLICIT NONE
  PUBLIC :: AV18PW, EMPOT
  PRIVATE:: AV18OP, CONSTS
  DOUBLE PRECISION :: ALPHA = 1.D0 / 137.035989D0
  DOUBLE PRECISION :: HC = 197.327053D0
  PUBLIC :: ALPHA, HC
  CONTAINS
  !> \brief Partial-wave projection of the Argonne v18 (or related) potential.
  !! \ingroup av18
  !!
  !! Computes the nucleon-nucleon potential matrix for a given partial wave,
  !! including strong and electromagnetic terms, for Argonne v18, v8', v6', v4', vx', v2', v1', and Super-Soft Core (C) models.
  !!
  !! Calls subroutines AV18OP and EMPOT.
  !!
  !! \par Potential choice (\a lpot):
  !! \code
  !!   -----------------------------------------------
  !!       Argonne                Super-Soft Core (C)
  !!     = 1 : av18              = 101 : sscc v14
  !!     = 2 : av8'              = 102 : sscc v8'
  !!     = 3 : av6'
  !!     = 4 : av4'
  !!     = 5 : avx'
  !!     = 6 : av2'
  !!     = 7 : av1'
  !!     = 8 : modified av8'     = 108 : modified sscc v8'
  !!   -----------------------------------------------
  !! \endcode
  !!
  !! \par Matrix structure:
  !! \li Single channel: v(1,1) = V( l, s, j, t, t1z, t2z ); others zero
  !! \li Coupled channel (l=j-1, s=1): v(1,1), v(2,2), v(1,2), v(2,1) as described in documentation
  !!
  !! \param[in]  lpot  Switch for potential choice
  !! \param[in]  l     Orbital angular momentum of pair (0,1,2,...)
  !! \param[in]  s     Total spin of pair (0 or 1)
  !! \param[in]  j     Total angular momentum of pair (0,1,2,...)
  !! \param[in]  t     Total isospin of pair (0 or 1)
  !! \param[in]  t1z   Isospin of particle 1 (1 for p, -1 for n)
  !! \param[in]  t2z   Isospin of particle 2 (1 for p, -1 for n)
  !! \param[in]  r     Separation in fm
  !! \param[out] vpw   Returned potential in MeV (2x2 array, includes all strong and EM terms)
  !! \param[in]  lemp  Electromagnetic potential flag
  SUBROUTINE AV18PW(LPOT,L,S,J,T,T1Z,T2Z,R,VPW,LEMP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LPOT, L, S, J, T, T1Z, T2Z, LEMP
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VPW(2,2)
  !
    INTEGER S1DS2, T1DT2, T12, NCC
    DOUBLE PRECISION :: VNN(18),VEM(14)
    DOUBLE PRECISION :: VC, VT, VLS, VL2, VLS2
    DOUBLE PRECISION :: S12, S12M, S12P, LS, LSM, LSP
  ! ------------------------
  ! strong interaction terms
  ! ------------------------
    CALL AV18OP(LPOT,R,VNN)
    S1DS2=4*S-3
    T1DT2=4*T-3
    T12=3*T1Z*T2Z-T1DT2
    VC=VNN(1)+T1DT2*VNN(2)+S1DS2*VNN(3)+S1DS2*T1DT2*VNN(4) &
          +T12*VNN(15)+S1DS2*T12*VNN(16)+(T1Z+T2Z)*VNN(18)
    VT=VNN(5)+T1DT2*VNN(6)+T12*VNN(17)
    VLS=VNN(7)+T1DT2*VNN(8)
    VL2=VNN(9)+T1DT2*VNN(10)+S1DS2*VNN(11)+S1DS2*T1DT2*VNN(12)
    VLS2=VNN(13)+T1DT2*VNN(14)
  ! ---------------------
  ! electromagnetic terms
  ! ---------------------
    VEM = 0
    IF (LEMP.EQ.0) THEN
      IF (T1Z+T2Z .EQ. 2) VEM(1) = ALPHA*HC/R
    ELSE
      CALL EMPOT(LPOT,R,VEM)
    ENDIF
    SELECT CASE (T1Z+T2Z)
    CASE (-2)
      VC=VC+S1DS2*VEM(7)
      VT=VT+VEM(10)
    CASE (0)
      VC=VC+VEM(5)+S1DS2*VEM(8)
      VT=VT+VEM(11)
      VLS=VLS+VEM(14)
    CASE (2)
      VC=VC+VEM(1)+VEM(2)+VEM(3)+VEM(4)+S1DS2*VEM(6)
      VT=VT+VEM(9)
      VLS=VLS+VEM(12)
    END SELECT

  ! ---------------------
    NCC=1
    IF (S.EQ.1.AND.J.GT.L) NCC=2
    IF (NCC.EQ.1) THEN
      S12=0.
      IF (S.EQ.1.AND.L.EQ.J) S12=2.
      IF (L.EQ.(J+1)) S12=-2.*(J+2.)/(2.*J+1.)
      LS=(J*(J+1)-L*(L+1)-S*(S+1))/2
      VPW(1,1)=VC+S12*VT+LS*VLS+L*(L+1)*VL2+LS**2*VLS2
      VPW(2,1)=0
      VPW(1,2)=0
      VPW(2,2)=0
    ELSE IF (NCC.EQ.2) THEN
      S12M=-2.*(J-1.)/(2.*J+1.)
      S12=SQRT(36.*J*(J+1))/(2.*J+1.)
      S12P=-2.*(J+2.)/(2.*J+1.)
      LSM=J-1
      LSP=-(J+2)
      VPW(1,1)=VC+S12M*VT+LSM*VLS+L*(L+1)*VL2+LSM**2*VLS2
      VPW(2,1)=S12*VT
      VPW(1,2)=S12*VT
      VPW(2,2)=VC+S12P*VT+LSP*VLS+(L+2)*(L+3)*VL2+LSP**2*VLS2
    END IF
  END SUBROUTINE AV18PW


  ! *id* av18op90 **********************************************************
  ! subroutine for strong interaction part of argonne v18 potential
  ! or super-soft core (C) v14 potential
  ! or reprojected vn' potential in operator format
  ! calls subroutine consts90
  ! ----------------------------------------------------------------------
  ! arguments for av18pot
  ! lpot: switch for potential choice
  !     -----------------------------------------------
  !         Argonne                Super-Soft Core (C)
  !       = 1 : av18              = 101 : sscc v14
  !       = 2 : av8'              = 102 : sscc v8'
  !       = 3 : av6'
  !       = 4 : av4'
  !       = 5 : avx'
  !       = 6 : av2'
  !       = 7 : av1'
  !       = 8 : modified av8'     = 108 : modified sscc v8'
  !     -----------------------------------------------
  ! r:    separation in fm
  ! vnn:  output potential in MeV (18 component array)
  ! ----------------------------------------------------------------------
  ! order of operators l in vnn(l):
  ! l:    1=1                              2=t1.t2
  !       3=s1.s2                          4=(s1.s2)(t1.t2)
  !       5=S12 [=3(s1.r)(s2.r)-s1.s2]     6=S12(t1.t2)
  !       7=L.S                            8=L.S(t1.t2)
  !       9=L**2                          10=L**2(t1.t2)
  !      11=L**2(s1.s2)                   12=L**2(s1.s2)(t1.t2)
  !      13=(L.S)**2                      14=(L.S)**2(t1.t2)
  !      15=T12 [=3*t1z*t2z-t1.t2]        16=(s1.s2)T12
  !      17=S12*T12                       18=t1z+t2z
  ! where s1=sigma_1, t1=tau_1, t1z=tau_1(z), etc.
  ! ----------------------------------------------------------------------
  SUBROUTINE AV18OP(LPOT,R,VNN)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LPOT
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VNN(18)

    DOUBLE PRECISION, PARAMETER :: SMALL=1E-4, VSMALL=1E-10
    DOUBLE PRECISION :: MPI0,MPIC,MP,MN,MUP,MUN
    DOUBLE PRECISION :: MPI,MU0,MUC,MU
    DOUBLE PRECISION :: X,X0,XC
    DOUBLE PRECISION :: P11PP,P11NP,P11NN,P11CD,P11CS,P10,P01,P00
    DOUBLE PRECISION :: P11,PT1,PT1CD,PT1CS, PT1PP, PT1NP,PT1NN
    DOUBLE PRECISION :: P01CD,P01CS,PT0,PLS1,PLS0,PL211,PL210
    DOUBLE PRECISION :: PL201,PL200,PLS21,PLS20,PQ0
    DOUBLE PRECISION :: P01PP,P01NP,P01NN
    DOUBLE PRECISION :: TPI,TPI0,TPIC,TPI2
    DOUBLE PRECISION :: YPI,YPI0,YPIC
    DOUBLE PRECISION :: YPI0P,YPICP
    DOUBLE PRECISION :: RCUT
    DOUBLE PRECISION :: FSQ, CPI,RWS,AIWS
    DOUBLE PRECISION :: WS,WS0,WSP,WSX,WSX2
    DOUBLE PRECISION :: DYPI00,DYPIC0
    DOUBLE PRECISION :: RR4, RC4, RC6
    DOUBLE PRECISION :: HR, RR
  ! -------------------
    VNN = 0
    RR = R
  ! ---------------------------------
  ! Argonne potential and derivatives
  ! ---------------------------------
    IF (LPOT.LT.100) THEN
    CALL CONSTS(LPOT,MPI0,MPIC,MP,MN,MUP,MUN)
    MPI=(MPI0+2.*MPIC)/3.
    MU0=MPI0/HC
    MUC=MPIC/HC
    MU=MPI/HC
    FSQ=.075
    CPI=2.1
    RWS=.5
    AIWS=5.
    X=MU*RR
    X0=MU0*RR
    XC=MUC*RR
    IF (RR.LE.SMALL) THEN
      TPI=3*CPI**2*RR/MU**3
      YPI0=(MPI0/MPIC)**2*(MPI0/3)*CPI*RR/MU0
      TPI0=3*CPI*YPI0/MU0**2
      YPIC=(MPIC/3)*CPI*RR/MUC
      TPIC=3*CPI*YPIC/MUC**2
    ELSE
      RCUT=1-EXP(-CPI*RR*RR)
      YPI=EXP(-X)*RCUT/X
      TPI=(1+(3+3/X)/X)*YPI*RCUT
      YPI0=(MPI0/MPIC)**2*(MPI0/3)*EXP(-X0)*RCUT/X0
      TPI0=(1+(3+3/X0)/X0)*YPI0*RCUT
      YPIC=(MPIC/3)*EXP(-XC)*RCUT/XC
      TPIC=(1+(3+3/XC)/XC)*YPIC*RCUT
    END IF
    YPI0=FSQ*YPI0
    YPIC=FSQ*YPIC
    TPI0=FSQ*TPI0
    TPIC=FSQ*TPIC
    TPI2=TPI*TPI
    WS=1/(1+EXP((RR-RWS)*AIWS))
    WS0=1/(1+EXP(-RWS*AIWS))
    WSP=WS*(1+AIWS*EXP(-RWS*AIWS)*WS0*RR)
    WSX=WS*X
    WSX2=WSX*X
    DYPI00=(MPI0/MPIC)**2*(MPI0/3)*CPI/MU0
    DYPIC0=(MPIC/3)*CPI/MUC
    YPI0P=YPI0-FSQ*DYPI00*WS*RR/WS0
    YPICP=YPIC-FSQ*DYPIC0*WS*RR/WS0
    YPI=(YPI0+2*YPIC)/3
    TPI=(TPI0+2*TPIC)/3
    P11PP=  -7.62701*TPI2+1815.4920*WSP+1847.8059*WSX2+YPI0P
    P11NP=  -7.62701*TPI2+1813.5315*WSP+1847.8059*WSX2-YPI0P+2*YPICP
    P11NN=  -7.62701*TPI2+1811.5710*WSP+1847.8059*WSX2+YPI0P
    PT1PP=   1.07985*TPI2 -190.0949*WSX -811.2040*WSX2+TPI0
    PT1NP=   1.07985*TPI2 -190.0949*WSX -811.2040*WSX2-TPI0+2*TPIC
    PT1NN=   1.07985*TPI2 -190.0949*WSX -811.2040*WSX2+TPI0
    PLS1=    -.62697*TPI2 -570.5571*WSP +819.1222*WSX2
    PL211=    .06709*TPI2 +342.0669*WSP -615.2339*WSX2
    PLS21=    .74129*TPI2   +9.3418*WSP -376.4384*WSX2
    P10=    -8.62770*TPI2+2605.2682*WSP +441.9733*WSX2-YPI0P-2*YPICP
    PT0=    1.485601*TPI2-1126.8359*WSX +370.1324*WSX2-TPI0-2*TPIC
    PLS0=     .10180*TPI2  +86.0658*WSP -356.5175*WSX2
    PL210=   -.13201*TPI2 +253.4350*WSP   -1.0076*WSX2
    PLS20=    .07357*TPI2 -217.5791*WSP  +18.3935*WSX2
    P01PP= -11.27028*TPI2+3346.6874*WSP-3*YPI0P
    P01NP= -10.66788*TPI2+3126.5542*WSP-3*(-YPI0P+2*YPICP)
    P01NN= -11.27028*TPI2+3342.7664*WSP-3*YPI0P
    PL201=    .12472*TPI2  +16.7780*WSP
    P00=    -2.09971*TPI2+1204.4301*WSP-3*(-YPI0P-2*YPICP)
    PL200=   -.31452*TPI2 +217.4559*WSP
    P11=(P11PP+P11NN+P11NP)/3
    P11CD=(.5*(P11PP+P11NN)-P11NP)/6
    P11CS=(P11PP-P11NN)/4
    PT1=(PT1PP+PT1NN+PT1NP)/3
    PT1CD=(.5*(PT1PP+PT1NN)-PT1NP)/6
    PT1CS=(PT1PP-PT1NN)/4
    P01=(P01PP+P01NN+P01NP)/3
    P01CD=(.5*(P01PP+P01NN)-P01NP)/6
    P01CS=(P01PP-P01NN)/4
  ! ------------------------
  ! option for v8' reduction
  ! ------------------------
    IF (LPOT.GE.2) THEN
      P00=P00+2*PL200
      P11=P11+2*PL211+4*PLS21/3
      PT1=PT1-5*PLS21/12
      PLS1=PLS1-.5*PLS21
      PLS0=PLS0-2*PL210-3*PLS20
    END IF
  ! ------------------------
  ! option for v6' redutcion
  ! ------------------------
    IF (LPOT.GE.3 .AND. LPOT.LE.7) P10=P10-.3*PLS0
  ! ------------------------
  ! option for v4' reduction
  ! ------------------------
    IF (LPOT.GE.4 .AND. LPOT.LE.7) P10=P10+.8735*PT0
  ! ------------------------
  ! option for vx' reduction
  ! ------------------------
    IF (LPOT.EQ.5) THEN
      VNN(1)=.0625*(9*P11+3*P10+3*P01+P00)
      VNN(2)=.0125*(9*P11-5*P10-5*P01+P00)
      VNN(3)=VNN(2)
      VNN(4)=VNN(2)
      RETURN
  ! ------------------------
  ! option for v2' reduction
  ! ------------------------
    ELSE IF (LPOT.EQ.6) THEN
      VNN(1)=.25*(3*P01+P10)
      VNN(2)=.25*(  P01-P10)
      RETURN
  ! ------------------------
  ! option for v1' reduction
  ! ------------------------
    ELSE IF (LPOT.EQ.7) THEN
      VNN(1)=.5*(P01+P10)
      RETURN
    END IF
  ! ------------------------
  ! option for modified v8'
  ! ------------------------
    IF (LPOT.EQ.8) P11=P11-.37*TPI2
  ! ------------------------
    VNN(1)=.0625*(9*P11+3*P10+3*P01+P00)
    VNN(2)=.0625*(3*P11-3*P10  +P01-P00)
    VNN(3)=.0625*(3*P11  +P10-3*P01-P00)
    VNN(4)=.0625*(  P11  -P10  -P01+P00)
  ! ------------------------
    IF (LPOT.EQ.4) RETURN
  ! ------------------------
    VNN(5)=.25*(3*PT1+PT0)
    VNN(6)=.25*(  PT1-PT0)
  ! ------------------------
    IF (LPOT.EQ.3) RETURN
  ! ------------------------
    VNN(7)=.25*(3*PLS1+PLS0)
    VNN(8)=.25*(  PLS1-PLS0)
  ! ------------------------
    IF (LPOT.EQ.2 .OR. LPOT.EQ.8) RETURN
  ! ------------------------
    VNN(9)= .0625*(9*PL211+3*PL210+3*PL201+PL200)
    VNN(10)=.0625*(3*PL211-3*PL210+  PL201-PL200)
    VNN(11)=.0625*(3*PL211+  PL210-3*PL201-PL200)
    VNN(12)=.0625*(  PL211-  PL210-  PL201+PL200)
    VNN(13)=.25*(3*PLS21+PLS20)
    VNN(14)=.25*(  PLS21-PLS20)
    VNN(15)=.25*(3*P11CD+P01CD)
    VNN(16)=.25*(  P11CD-P01CD)
    VNN(17)=PT1CD
    VNN(18)=P01CS
    RETURN
  ! ---------------------------------------------
  ! super-soft core (C) potential and derivatives
  ! ---------------------------------------------
    ELSE IF (LPOT.GT.100) THEN
    IF (RR.LE.VSMALL) RR=VSMALL
    X=.7*RR
    RR4=RR**4
    RC4=1-EXP(-RR4)
    RC6=1-EXP(-RR**6)
    HR=10.463
    P11=144.83*EXP(-RR4/.88787**2)&
        +(-241.34*YC(3.3788*X)+(HR/3)*YC(X))*RC4
    P10=215.32*EXP(-RR4/.85807**2)&
        +(-883.6*YC(3.5042*X)-HR*YC(X))*RC4
    P01=375.*EXP(-RR4/.47552**2)&
        +(-1001.6*YC(3.6071*X)-HR*YC(X))*RC4
    P00=75.653*EXP(-RR4/3.0000**2)&
        +(-286.26*YC(2.0254*X)+3*HR*YC(X))*RC4
    PT1=36.*EXP(-RR4/1.0805**2)&
        +(-110.*YT(3.9529*X)+(HR/3)*YT(X))*RC6
    PT0=-58.951*EXP(-RR4/1.3171**2)&
        +(395.18*YT(4.3098*X)-HR*YT(X))*RC6
    PLS1=(520.*YLS(5.661*X)-54.85*YLS(4.0141*X))*RC6
    PLS0=(-40.466*YLS(5.768*X)-40.408*YLS(4.0676*X))*RC6
    PL211=(6.65*YL2(1.965*X)-.959*YL2(X))*RC6
    PL210=(17.626*YL2(2.6463*X)-.35261*YL2(X))*RC6
    PL201=(14.*YL2(2.5*X)-.35*YL2(X))*RC6
    PL200=(15.633*YL2(2.01*X)+.72581*YL2(X))*RC6
    PQ0=-3.9904*YL2(2.4583*X)*RC6
  ! ------------------------
  ! option for v8' reduction
  ! ------------------------
    IF (LPOT.GE.102) THEN
        P00=P00+2*PL200
        PLS0=PLS0-2*PL210-10*PQ0
        P11=P11+2*PL211
    END IF
  ! ------------------------
  ! option for modified v8'
  ! ------------------------
    IF (LPOT.EQ.108) P11=P11-111*YC(3.3788*X)*RC4
  ! ------------------------
    VNN(1)=.0625*(9*P11+3*P10+3*P01+P00)
    VNN(2)=.0625*(3*P11-3*P10+  P01-P00)
    VNN(3)=.0625*(3*P11+  P10-3*P01-P00)
    VNN(4)=.0625*(  P11-  P10-  P01+P00)
    VNN(5)=.25*(3*PT1+PT0)
    VNN(6)=.25*(  PT1-PT0)
    VNN(7)=.25*(3*PLS1+PLS0)+.75*PQ0
    VNN(8)=.25*(  PLS1-PLS0)-.75*PQ0
  ! ------------------------
    IF (LPOT.GE.102) RETURN
  ! ------------------------
    VNN(9)= .0625*(9*PL211+3*PL210+3*PL201+PL200)-.75*PQ0
    VNN(10)=.0625*(3*PL211-3*PL210+  PL201-PL200)+.75*PQ0
    VNN(11)=.0625*(3*PL211+  PL210-3*PL201-PL200)-.25*PQ0
    VNN(12)=.0625*(  PL211-  PL210-  PL201+PL200)+.25*PQ0
    VNN(13)=1.5*PQ0
    VNN(14)=-1.5*PQ0
    END IF

  ! -------------------
  ! statement functions
  ! -------------------
    CONTAINS
      PURE DOUBLE PRECISION FUNCTION YC(T)
        DOUBLE PRECISION, INTENT(IN) :: T
        YC = EXP(-T) / X
      END FUNCTION YC

      PURE DOUBLE PRECISION FUNCTION YT(T)
        DOUBLE PRECISION, INTENT(IN) :: T
        YT = (1 + 3/T + 3/T**2) * EXP(-T) / X
      END FUNCTION YT

      PURE DOUBLE PRECISION FUNCTION YLS(T)
        DOUBLE PRECISION, INTENT(IN) :: T
        YLS = -(1 + T) * EXP(-T) / X**3
      END FUNCTION YLS

      PURE DOUBLE PRECISION FUNCTION YL2(T)
        DOUBLE PRECISION, INTENT(IN) :: T
        YL2 = (1 + 2/T) * EXP(-T) / X**3
      END FUNCTION YL2
  END SUBROUTINE AV18OP

  !> \brief Electromagnetic part of the Argonne v18 potential.
  !! \ingroup av18
  !!
  !! Computes the electromagnetic part of the Argonne v18 potential.
  !! For avn' models, returns pp Coulomb only.
  !! Calls subroutine CONSTS.
  !!
  !! \par Potential choice (\a lpot):
  !!   - \c lpot = 1 : full EM potential
  !!   - \c lpot > 1 : C1(pp) only
  !!
  !! \param[in]  lpot  Switch for potential choice
  !! \param[in]  r     Input separation in fm
  !! \param[out] vem   Output potential in MeV (14 component array)
  !!
  !! \par Order of operators in \a vem(l):
  !! \code
  !! l:  1=C1    (pp)          2=DF    (pp)           3=C2      (pp)
  !!     4=VP    (pp)                             5=C1      (np)
  !!     6=s1.s2 (pp)       7=s1.s2 (nn)        8=s1.s2   (np)
  !!     9=S12   (pp)         10=S12   (nn)         11=S12     (np)
  !!    12=L.S   (pp)         13=L.S   (nn)         14=L.S     (np)
  !! \endcode
  !!
  !! - C1 = one-photon-exchange Coulomb with form factor
  !! - C2 = two-photon-exchange Coulomb
  !! - DF = Darwin-Foldy
  !! - VP = vacuum polarization (short-range approximation)
  !! - All other terms from magnetic moment (MM) interactions
  SUBROUTINE EMPOT(LPOT, R, VEM)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LPOT
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VEM(14)

    DOUBLE PRECISION :: MPI0, MPIC, MP, MN, MUP, MUN
    DOUBLE PRECISION :: B, BR, PI, ME, MR, GAMMA_VAL, BETA
    DOUBLE PRECISION :: FCOULR, FTR3, FLSR3, KR, FIVP, FDELTA, FNPR
    DOUBLE PRECISION, PARAMETER :: SMALL = 1.0E-5
    DOUBLE PRECISION, PARAMETER :: GAMMA_EULER = 0.577216

    VEM = 0.0D0
    CALL CONSTS(LPOT, MPI0, MPIC, MP, MN, MUP, MUN)

    B = 4.27
    BR = B * R
    PI = ACOS(-1.0)
    ME = 0.510999
    MR = MP * MN / (MP + MN)
    GAMMA_VAL = GAMMA_EULER
    BETA = 0.0189

    IF (R < SMALL) THEN
      FCOULR = 5.0 * B / 16.0
      FTR3 = B**3 * BR**2 / 720.0
      FLSR3 = B**3 / 48.0
      KR = ME * SMALL / HC
    ELSE
      FCOULR = (1.0 - (1.0 + 11.0*BR/16.0 + 3.0*BR**2/16.0 + BR**3/48.0) * EXP(-BR)) / R
      FTR3 = (1.0 - (1.0 + BR + BR**2/2.0 + BR**3/6.0 + BR**4/24.0 + BR**5/144.0) * EXP(-BR)) / R**3
      FLSR3 = (1.0 - (1.0 + BR + BR**2/2.0 + 7.0*BR**3/48.0 + BR**4/48.0) * EXP(-BR)) / R**3
      KR = ME * R / HC
    END IF

    FIVP = -GAMMA_VAL - 5.0/6.0 + ABS(LOG(KR)) + 6.0*PI*KR/8.0
    FDELTA = B**3 * (1.0 + BR + BR**2/3.0) * EXP(-BR) / 16.0
    FNPR = B**3 * (15.0 + 15.0*BR + 6.0*BR**2 + BR**3) * EXP(-BR) / 384.0

    VEM(1) = ALPHA * HC * FCOULR

  ! ------------------------
    IF (LPOT >= 2) RETURN
  ! ------------------------

    VEM(2) = -ALPHA * HC**3 * FDELTA / (4.0 * MP**2)
    VEM(3) = -VEM(1)**2 / MP
    VEM(4) = 2.0 * ALPHA * VEM(1) * FIVP / (3.0 * PI)
    VEM(5) = ALPHA * HC * BETA * FNPR
    VEM(6) = -ALPHA * HC**3 * MUP**2 * FDELTA / (6.0 * MP**2)
    VEM(7) = -ALPHA * HC**3 * MUN**2 * FDELTA / (6.0 * MN**2)
    VEM(8) = -ALPHA * HC**3 * MUP * MUN * FDELTA / (6.0 * MN * MP)
    VEM(9) = -ALPHA * HC**3 * MUP**2 * FTR3 / (4.0 * MP**2)
    VEM(10) = -ALPHA * HC**3 * MUN**2 * FTR3 / (4.0 * MN**2)
    VEM(11) = -ALPHA * HC**3 * MUP * MUN * FTR3 / (4.0 * MP * MN)
    VEM(12) = -ALPHA * HC**3 * (4.0 * MUP - 1.0) * FLSR3 / (2.0 * MP**2)
    VEM(13) = 0.0
    VEM(14) = -ALPHA * HC**3 * MUN * FLSR3 / (2.0 * MN * MR)
  END SUBROUTINE EMPOT


  !> *id* consts90 **********************************************************
  !! subroutine for constants in av18 and sscc potentials
  !! ----------------------------------------------------------------------
  !! \param[in] lpot  input potential
  !! \param[out] hc    output value for hbar*c (MeV-fm)
  !! \param[out] mpi0  neutral pion mass (MeV)
  !! \param[out] mpic  charged pion mass (MeV)
  !! \param[out] mp    proton mass (MeV)
  !! \param[out] mn    neutron mass (MeV)
  !! \param[out] alpha electromagnetic constant alpha
  !! \param[out] mup   proton magnetic moment (nm)
  !! \param[out] mun   neutron magnetic moment (nm)
  !! arguments for consts90
  !! lpot:  input potential
  !! hc:    output value for hbar*c (MeV-fm)
  !! mpi0:    "      "    "  neutral pion mass (MeV)
  !! mpic:    "      "    "  charged pion mass (MeV)
  !! mp:      "      "    "  proton mass (MeV)
  !! mn:      "      "    "  neutron mass (MeV)
  !! alpha:   "      "    "  electromagnetic constant alpha
  !! mup:     "      "    "  proton magnetic moment (nm)
  !! mun:     "      "    "  neutron magnetic moment (nm)
  !! ----------------------------------------------------------------------
  SUBROUTINE CONSTS(LPOT, MPI0, MPIC, MP, MN, MUP, MUN)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LPOT
    DOUBLE PRECISION, INTENT(OUT) :: MPI0, MPIC, MP, MN, MUP, MUN

    HC = 197.327053

    IF (LPOT < 100) THEN
      MPI0 = 134.9739
      MPIC = 139.5675
      MP = 938.27231
      MN = 939.56563
    ELSE
      MPI0 = 0.7 * HC
      MPIC = 0.7 * HC
      MP = HC**2 / 41.47
      MN = HC**2 / 41.47
    END IF

    ALPHA = 1.0 / 137.035989
    MUP = 2.7928474
    MUN = -1.9130427

  END SUBROUTINE CONSTS
END MODULE AV18