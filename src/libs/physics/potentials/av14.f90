!> \file av14.f90
!! \defgroup av14 Argonne v14
!! \ingroup nn_potentials
!! \brief Module for the Argonne v14 nuclear potential.
!!
!! This module provides the AV14 nuclear potential function
!! for the calculation of nucleon-nucleon interactions in different partial-wave channels,
!! including isospin effects and electromagnetic terms.
!!
!! \author Alessandro Grassi
!! \date 2025
MODULE AV14
  USE AV18
  IMPLICIT NONE
  PRIVATE

  !> Total angular momentum (J) of the nucleon-nucleon system.
  INTEGER, PRIVATE :: J_
  !> Type of electromagnetic interaction (0=none, >0=present).
  INTEGER, PRIVATE :: LEMP_
  !> Flag to include the electromagnetic tensor term.
  INTEGER, PRIVATE :: IPTE
  !> Flag to include the spin-orbit term.
  INTEGER, PRIVATE :: IPLS
  !> Flag to include the quadrupole term.
  INTEGER, PRIVATE :: IPQQ
  !> Flag to include the Breit term.
  INTEGER, PRIVATE :: IPBB
  !> Flag to include the spin-bremsstrahlung term.
  INTEGER, PRIVATE :: IPSB
  !> Interparticle distance (fm).
  DOUBLE PRECISION, PRIVATE :: R_
  !> Inverse interparticle distance (1/fm).
  DOUBLE PRECISION, PRIVATE :: UR

  !> Public subroutine for AV14 potential calculation in partial-wave channels.
  !! \ingroup av14
  PUBLIC :: AV14PW
  !> Private support subroutines.
  PRIVATE:: POTL, POT
CONTAINS

  !> \brief Compute the AV14 potential matrix for given quantum numbers and distance.
  !! \ingroup av14
  !! \param[in] LEMP Type of electromagnetic interaction (0=Coulomb, >0=submodel from AV18)
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Total spin (0 or 1)
  !! \param[in] J Total angular momentum
  !! \param[in] T1Z Isospin projection nucleon 1 (+1=proton, -1=neutron)
  !! \param[in] T2Z Isospin projection nucleon 2 (+1=proton, -1=neutron)
  !! \param[in] R Interparticle distance (fm)
  !! \param[out] VPW 2x2 potential matrix in coupled basis (mixed channels)
  !!
  !! The VPW matrix contains the potential values for the different spin and angular momentum coupling channels.
  SUBROUTINE AV14PW(LEMP, L, S, J, T1Z, T2Z, R, VPW)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LEMP
    INTEGER, INTENT(IN) :: T1Z
    INTEGER, INTENT(IN) :: T2Z
    INTEGER, INTENT(IN) :: S
    INTEGER, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: J
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VPW(2, 2)
    INTEGER :: INN
    DOUBLE PRECISION :: VX1, VX2, VX3, VX4, VX5
    DOUBLE PRECISION :: VPP1, VPP2, VPP3, VPP4, VPP5
    DOUBLE PRECISION :: VPN1, VPN2, VPN3, VPN4, VPN5
    DOUBLE PRECISION :: VNN1, VNN2, VNN3, VNN4, VNN5

    J_ = J
    LEMP_ = LEMP
    R_ = R
    UR = 1.D0/R_

    IPTE=1
    IPLS=1
    IPQQ=1
    IPBB=1
    IPSB=1

    VPW = 0.D0

    SELECT CASE (T1Z + T2Z)
      CASE (-2)
        INN = 1 ! pp
      CASE (0)
        INN = 2 ! pn
      CASE (2)
        INN = 3 ! nn
      CASE DEFAULT
        WRITE(*, *) "ERROR IN AV14PW: Invalid T1Z + T2Z = ", T1Z + T2Z
        STOP
    END SELECT

    CALL POTL(VPP1, VPP2, VPP3, VPP4, VPP5, &
              VPN1, VPN2, VPN3, VPN4, VPN5, &
              VNN1, VNN2, VNN3, VNN4, VNN5)

    SELECT CASE (INN)
      CASE (1)
        VX1 = VPP1; VX2 = VPP2; VX3 = VPP3; VX4 = VPP4; VX5 = VPP5
      CASE (2)
        VX1 = VPN1; VX2 = VPN2; VX3 = VPN3; VX4 = VPN4; VX5 = VPN5
      CASE (3)
        VX1 = VNN1; VX2 = VNN2; VX3 = VNN3; VX4 = VNN4; VX5 = VNN5
    END SELECT

    IF (J == L .AND. S == 0) THEN
      VPW(1, 1) = VX1
    ELSE IF (J == L .AND. S == 1) THEN
      VPW(1, 1) = VX2
    ELSE IF (J == 0) THEN
      VPW(1, 1) = VX4
    ELSE
      VPW(1, 1) = VX3
      VPW(1, 2) = VX5
      VPW(2, 1) = VX5
      VPW(2, 2) = VX4
    END IF
  END SUBROUTINE AV14PW

  !> \brief Compute the AV14 potential terms for all isospin combinations.
  !! \ingroup av14
  !! \param[out] VPP1,VPP2,VPP3,VPP4,VPP5 Terms for proton-proton (pp)
  !! \param[out] VPN1,VPN2,VPN3,VPN4,VPN5 Terms for proton-neutron (pn)
  !! \param[out] VNN1,VNN2,VNN3,VNN4,VNN5 Terms for neutron-neutron (nn)
  !!
  !! Each term corresponds to a different spin and angular momentum coupling channel.
  SUBROUTINE POTL(VPP1, VPP2, VPP3, VPP4, VPP5, VPN1, VPN2, VPN3, VPN4, VPN5, &
                  VNN1, VNN2, VNN3, VNN4, VNN5)
    IMPLICIT NONE
    !> Conversion factor for the Coulomb term (MeV*fm)
    DOUBLE PRECISION, PARAMETER :: EQMF=197.327053D0/137.035989D0
    !> Output terms for the different channels and interactions
    DOUBLE PRECISION, INTENT(OUT) :: VPP1, VPP2, VPP3, VPP4, VPP5, VPN1, VPN2, VPN3, VPN4, VPN5, &
                  VNN1, VNN2, VNN3, VNN4, VNN5
    !> Basic potential terms
    DOUBLE PRECISION :: POTC, POTT, POTS, POTM, PTEC, PTET
    DOUBLE PRECISION :: PLSC, PLST, PQQC, PQQT, PQQS, PQQM, PBBC, PBBT
    DOUBLE PRECISION :: PSB1, PSB2, PSB3, PSB4
    !> Electromagnetic terms and temporary variables
    DOUBLE PRECISION :: VCOUL, VEM(14), P1, P2, P3, P4
    DOUBLE PRECISION :: VEMPP, VEMPN, VEMNN
    DOUBLE PRECISION :: XTE, XLS, XLL, XLS2, UJ, XJJ
    DOUBLE PRECISION :: ZERO, ONE, U2, U3, U4, U6
    INTEGER :: JP, LTI, ICONT
    INTEGER :: LTS0(0:1), LTS1(0:1), LTIN(0:1)
    DATA ZERO /0.D0/, ONE /1.D0/, U2 /2.D0/, U3 /3.D0/, U4 /4.D0/, U6 /6.D0/
    DATA LTS0/1,0/     !LTS0(MOD(L))=VALUE OF T: CASE S=0
    DATA LTS1/0,1/     !LTS1(MOD(L))=VALUE OF T: CASE S=1
    DATA LTIN/1,0/     !LTIN(T=0)=1 AND LTIN(T=1)=0
    DATA ICONT /0/

    JP = MOD(J_, 2)
    UJ = ONE / DBLE(2 * J_ + 1)
    XJJ = DBLE(J_ * (J_ + 1))

    VPP1 = ZERO
    VPP2 = ZERO
    VPP3 = ZERO
    VPP4 = ZERO
    VPP5 = ZERO
    VPN1 = ZERO
    VPN2 = ZERO
    VPN3 = ZERO
    VPN4 = ZERO
    VPN5 = ZERO
    VNN1 = ZERO
    VNN2 = ZERO
    VNN3 = ZERO
    VNN4 = ZERO
    VNN5 = ZERO

    CALL POT(POTC, POTT, POTS, POTM, PTEC, PTET, PLSC, PLST, &
      PQQC, PQQT, PQQS, PQQM, PBBC, PBBT, PSB1, PSB2, PSB3, PSB4)

    IF (IPTE == 0) THEN
      PTEC = ZERO
      PTET = ZERO
    END IF
    IF (IPLS == 0) THEN
      PLSC = ZERO
      PLST = ZERO
    END IF
    IF (IPQQ == 0) THEN
      PQQC = ZERO
      PQQT = ZERO
      PQQS = ZERO
      PQQM = ZERO
    END IF
    IF (IPBB == 0) THEN
      PBBC = ZERO
      PBBT = ZERO
    END IF
    IF (IPSB == 0) THEN
      PSB1 = ZERO
      PSB2 = ZERO
      PSB3 = ZERO
      PSB4 = ZERO
    END IF

    ICONT = ICONT + 1
    IF (ICONT == 1) THEN
      WRITE(6, *) IPTE, IPLS, IPQQ, IPBB, IPSB
      WRITE(6, *) POTC, POTT, POTS, POTM
      WRITE(6, *) PTEC, PTET, PLSC, PLST
      WRITE(6, *) PQQC, PQQT, PQQS, PQQM
      WRITE(6, *) PBBC, PBBT
      WRITE(6, *) PSB1, PSB2, PSB3, PSB4
    END IF

    VCOUL = ZERO
    VEM   = ZERO
    IF(LEMP_.EQ.0) VCOUL=EQMF*UR
    IF(LEMP_.GT.0) CALL EMPOT(LEMP_,R_,VEM)

    ! CASE 1: SI=SJ=0, LI=LJ=J
    LTI = LTS0(JP)
    P1 = POTC - U3 * POTS + XJJ * (PQQC - U3 * PQQS)
    P2 = POTT - U3 * POTM + XJJ * (PQQT - U3 * PQQM)
    P3 = PSB1 - U3 * PSB2
    P4 = PSB4

    IF (LTI == 1) THEN
      VEMPP = VCOUL + VEM(1) + VEM(2) + VEM(3) + VEM(4) - U3 * VEM(6)
      VEMPN = VEM(5) - U3 * VEM(8)
      VEMNN = -U3 * VEM(7)
      VPP1 = P1 + P2 + U2 * (P3 + P4) + VEMPP
      VPN1 = P1 + P2 - U4 * P3 + VEMPN
      VNN1 = P1 + P2 + U2 * (P3 - P4) + VEMNN
    ELSE
      VEMPN = VEM(5) - U3 * VEM(8)
      VPN1 = P1 - U3 * P2 + VEMPN
      VNN1 = 0.0D0
      VPP1 = 0.0D0
    END IF

    ! CASE 2: SI=SJ=1, LI=LJ=J
    LTI = LTS1(JP)
    IF (J_ > 0) THEN
      P1 = POTC + POTS + U2 * PTEC - PLSC + XJJ * (PQQC + PQQS) + PBBC
      P2 = POTT + POTM + U2 * PTET - PLST + XJJ * (PQQT + PQQM) + PBBT
      P3 = PSB1 + PSB2 + U2 * PSB3
      P4 = PSB4

      IF (LTI == 1) THEN
        VEMPP = VCOUL + VEM(1) + VEM(2) + VEM(3) + VEM(4) + VEM(6) + U2 * VEM(9) - VEM(12)
        VEMPN = VEM(5) + VEM(8) + U2 * VEM(11) - VEM(14)
        VEMNN = VEM(7) + U2 * VEM(10) - VEM(13)
        VPP2 = P1 + P2 + U2 * (P3 + P4) + VEMPP
        VPN2 = P1 + P2 - U4 * P3 + VEMPN
        VNN2 = P1 + P2 + U2 * (P3 - P4) + VEMNN
      ELSE
        VEMPN = VEM(5) + VEM(8) + U2 * VEM(11) - VEM(14)
        VPP2 = 0.0D0
        VPN2 = P1 - U3 * P2 + VEMPN
        VNN2 = 0.0D0
      END IF
    END IF

    ! CASE 3: SI=SJ=1, LI=LJ=J-1
    LTI = LTIN(LTI)
    IF (J_ > 0) THEN
      XTE = -U2 * (J_ - 1) * UJ
      XLS = J_ - 1
      XLL = J_ * (J_ - 1)
      XLS2 = XLS * XLS
      P1 = POTC + POTS + XTE * PTEC + XLS * PLSC + XLL * (PQQC + PQQS) + XLS2 * PBBC
      P2 = POTT + POTM + XTE * PTET + XLS * PLST + XLL * (PQQT + PQQM) + XLS2 * PBBT
      P3 = PSB1 + PSB2 + XTE * PSB3
      P4 = PSB4

      IF (LTI == 1) THEN
        VEMPP = VCOUL + VEM(1) + VEM(2) + VEM(3) + VEM(4) + VEM(6) + XTE * VEM(9) + XLS * VEM(12)
        VEMPN = VEM(5) + VEM(8) + XTE * VEM(11) + XLS * VEM(14)
        VEMNN = VEM(7) + XTE * VEM(10) + XLS * VEM(13)
        VPP3 = P1 + P2 + U2 * (P3 + P4) + VEMPP
        VPN3 = P1 + P2 - U4 * P3 + VEMPN
        VNN3 = P1 + P2 + U2 * (P3 - P4) + VEMNN
      ELSE
        VEMPN = VEM(5) + VEM(8) + XTE * VEM(11) + XLS * VEM(14)
        VPN3 = P1 - U3 * P2 + VEMPN
        VNN3 = 0.0D0
        VPP3 = 0.0D0
      END IF
    END IF

    ! CASE 4: SI=SJ=1, LI=LJ=J+1
    XTE = -U2 * (J_ + 2) * UJ
    XLS = -J_ - 2
    XLL = (J_ + 1) * (J_ + 2)
    XLS2 = XLS * XLS
    P1 = POTC + POTS + XTE * PTEC + XLS * PLSC + XLL * (PQQC + PQQS) + XLS2 * PBBC
    P2 = POTT + POTM + XTE * PTET + XLS * PLST + XLL * (PQQT + PQQM) + XLS2 * PBBT
    P3 = PSB1 + PSB2 + XTE * PSB3
    P4 = PSB4

    IF (LTI == 1) THEN
      VEMPP = VCOUL + VEM(1) + VEM(2) + VEM(3) + VEM(4) + VEM(6) + XTE * VEM(9) + XLS * VEM(12)
      VEMPN = VEM(5) + VEM(8) + XTE * VEM(11) + XLS * VEM(14)
      VEMNN = VEM(7) + XTE * VEM(10) + XLS * VEM(13)
      VPP4 = P1 + P2 + U2 * (P3 + P4) + VEMPP
      VPN4 = P1 + P2 - U4 * P3 + VEMPN
      VNN4 = P1 + P2 + U2 * (P3 - P4) + VEMNN
    ELSE
      VEMPN = VEM(5) + VEM(8) + XTE * VEM(11) + XLS * VEM(14)
      VPN4 = P1 - U3 * P2 + VEMPN
      VNN4 = 0.0D0
      VPP4 = 0.0D0
    END IF

    ! CASE 5: SI=SJ=1, LI=J-1, LJ=J+1
    IF (J_ > 0) THEN
      XTE = U6 * DSQRT(XJJ) * UJ
      P1 = XTE * PTEC
      P2 = XTE * PTET
      P3 = XTE * PSB3
      P4 = 0.0D0

      IF (LTI == 1) THEN
        VEMPP = XTE * VEM(9)
        VEMPN = XTE * VEM(11)
        VEMNN = XTE * VEM(10)
        VPP5 = P1 + P2 + U2 * (P3 + P4) + VEMPP
        VPN5 = P1 + P2 - U4 * P3 + VEMPN
        VNN5 = P1 + P2 + U2 * (P3 - P4) + VEMNN
      ELSE
        VEMPN = XTE * VEM(11)
        VPN5 = P1 - U3 * P2 + VEMPN
        VNN5 = 0.0D0
        VPP5 = 0.0D0
      END IF
    END IF
    ! write(6,*) "vpp1,vpp2,vpp3,vpp4,vpp5"
    ! write(6,*) vpp1,vpp2,vpp3,vpp4,vpp5
    ! write(6,*) "vpn1,vpn2,vpn3,vpn4,vpn5"
    ! write(6,*) vpn1,vpn2,vpn3,vpn4,vpn5
    ! write(6,*) "vnn1,vnn2,vnn3,vnn4,vnn5"
    ! write(6,*) vnn1,vnn2,vnn3,vnn4,vnn5
    ! stop
    RETURN
  END SUBROUTINE POTL

  !> \brief Compute the basic AV14 potential terms for a given distance.
  !! \ingroup av14
  !! \param[out] POTC Central term
  !! \param[out] POTT Tensor term
  !! \param[out] POTS Spin-spin term
  !! \param[out] POTM Spin-orbit term
  !! \param[out] PTEC Electromagnetic tensor term
  !! \param[out] PTET Transverse tensor term
  !! \param[out] PLSC Central spin-orbit term
  !! \param[out] PLST Transverse spin-orbit term
  !! \param[out] PQQC Central quadrupole term
  !! \param[out] PQQT Tensor quadrupole term
  !! \param[out] PQQS Spin-spin quadrupole term
  !! \param[out] PQQM Mixed quadrupole term
  !! \param[out] PBBC Central Breit term
  !! \param[out] PBBT Tensor Breit term
  !! \param[out] PSB1, PSB2, PSB3, PSB4 Spin-bremsstrahlung terms
  !!
  !! Each term represents a specific physical contribution to the nuclear potential.
  SUBROUTINE POT(POTC, POTT, POTS, POTM, PTEC, PTET, PLSC, PLST, &
          PQQC, PQQT, PQQS, PQQM, PBBC, PBBT, PSB1, PSB2, PSB3, PSB4)
      IMPLICIT NONE
      !> Physical constants and model parameters
      DOUBLE PRECISION, PARAMETER :: DMU2  = 0.699488167D0
      DOUBLE PRECISION, PARAMETER :: UDMU2 = 1.D0 / DMU2
      !> Output terms
      DOUBLE PRECISION, INTENT(OUT) :: POTC, POTT, POTS, POTM, PTEC, PTET
      DOUBLE PRECISION, INTENT(OUT) :: PLSC, PLST, PQQC, PQQT, PQQS, PQQM, PBBC, PBBT
      DOUBLE PRECISION, INTENT(OUT) :: PSB1, PSB2, PSB3, PSB4
      !> Temporary variables for the calculation of the terms
      DOUBLE PRECISION :: X, UX, X2, CUTOFF, YX, YP, TP, TP2, WW


      POTC = 0.D0
      POTS = 0.D0
      POTT = 0.D0
      POTM = 0.D0
      PTEC = 0.D0
      PTET = 0.D0
      PLSC = 0.D0
      PLST = 0.D0
      PQQC = 0.D0
      PQQT = 0.D0
      PQQS = 0.D0
      PQQM = 0.D0
      PBBC = 0.D0
      PBBT = 0.D0
      PSB1 = 0.D0
      PSB2 = 0.D0
      PSB3 = 0.D0
      PSB4 = 0.D0

      IF (R_ .GT. 90.D0) RETURN

      X = R_ * DMU2
      UX = UR * UDMU2
      X2 = 2.D0 * R_ * R_
      IF (X2 .GT. 40.D0) THEN
        CUTOFF = 1.D0
      ELSE
        CUTOFF = 1.D0 - DEXP(-X2)
      END IF

      YX = DEXP(-X) * UX
      YP = YX * CUTOFF
      TP = (1.D0 + 3.D0 * UX + 3.D0 * UX * UX) * YP * CUTOFF
      TP2 = TP * TP
      WW = 1.D0 / (1.D0 + DEXP((R_ - 0.5D0) * 5.D0))

      POTC = -4.801125D0 * TP2 + 2061.5625D0 * WW
      POTT =  0.798925D0 * TP2 -  477.3125D0 * WW
      POTS =  1.189325D0 * TP2 -  502.3125D0 * WW
      POTM =  0.182875D0 * TP2 +   97.0625D0 * WW + 3.72681D0 * YP
      PTEC = -0.1575D0   * TP2 +  108.75D0   * WW
      PTET = -0.7525D0   * TP2 +  297.25D0   * WW + 3.72681D0 * TP

      IF (IPLS .EQ. 0) RETURN
      PLSC =  0.5625D0   * TP2 -  719.75D0   * WW
      PLST =  0.0475D0   * TP2 -  159.25D0   * WW

      IF (IPQQ .EQ. 0) RETURN
      PQQC =  0.070625D0 * TP2 +    8.625D0  * WW
      PQQT = -0.148125D0 * TP2 +    5.625D0  * WW
      PQQS = -0.040625D0 * TP2 +   17.375D0  * WW
      PQQM = -0.001875D0 * TP2 -   33.625D0  * WW

      IF (IPBB .EQ. 0) RETURN
      PBBC = -0.5425D0   * TP2 +  391.0D0    * WW
      PBBT =  0.0025D0   * TP2 +  145.0D0    * WW

      RETURN
    END SUBROUTINE POT
END MODULE AV14