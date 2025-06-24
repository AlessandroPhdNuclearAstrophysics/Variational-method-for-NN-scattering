!> \file physical_constants.f90
!! \defgroup phys_const Physical Constants
!! \ingroup physics
!! \brief Physical constants and utility functions for nuclear and particle physics calculations.
!!
!! This module defines fundamental physical constants (masses, hbar*c, fine-structure constant, etc.)
!! and provides utility functions for reduced mass and kinetic energy factors for nucleons and leptons.
!!
MODULE PHYSICAL_CONSTANTS
  USE ANGLES
  IMPLICIT NONE

  !> \brief Planck constant times speed of light (\f$\hbar c\f$) in MeV fm
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: HBARC = 197.327053D0  ! MeV fm
  !> \brief Proton mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_PROTON = 938.272029D0  ! MeV
  !> \brief Neutron mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_NEUTRON = 939.565630D0  ! MeV
  !> \brief Average nucleon mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_NUCLEON = (MASS_PROTON + MASS_NEUTRON) / 2.D0  ! MeV
  !> \brief Reduced mass of neutron-proton system in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_REDUCED_NP = MASS_PROTON*MASS_NEUTRON/(MASS_PROTON+MASS_NEUTRON)  ! MeV
  !> \brief Deuteron mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_DEUTERON = 1875.612D0  ! MeV
  !> \brief Triton mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_TRITON = 2808.921D0  ! MeV
  !> \brief Helium-3 mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_HELIUM3 = 2809.415D0  ! MeV
  !> \brief Helium-4 mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_HELIUM4 = 3727.379D0  ! MeV
  !> \brief Charged pion mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_PION = 139.57039D0  ! MeV

  !> \brief \f$\hbar^2/(2m_p)\f$ in MeV fm^2
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: HBARM_P = HBARC**2 / (2.0D0 * MASS_PROTON)  ! MeV fm^2
  !> \brief \f$\hbar^2/(2m_n)\f$ in MeV fm^2
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: HBARM_N = HBARC**2 / (2.0D0 * MASS_NEUTRON)  ! MeV fm^2
  !> \brief \f$\hbar^2/(2\mu_{np})\f$ in MeV fm^2
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: HBARM_NP = HBARC**2 / (2.0D0 * MASS_REDUCED_NP)  ! MeV fm^2
  
  !> \brief Electron mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_ELECTRON = 0.510998946D0  ! MeV
  !> \brief Muon mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_MUON = 105.6583745D0  ! MeV
  !> \brief Tau mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_TAU = 1776.86D0  ! MeV
  !> \brief Z boson mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_Z_BOSON = 9.11876D3  ! MeV
  !> \brief W boson mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_W_BOSON = 8.0379D3  ! MeV
  !> \brief Higgs boson mass in MeV
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: MASS_HIGGS_BOSON = 125.1D3  ! MeV

  !> \brief Fine-structure constant (dimensionless)
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: ALPHA_FINE_STRUCTURE = 1.0D0 / 137.035999084D0  ! Fine-structure constant
  !> \brief Fine-structure constant times hbar*c (MeV fm)
  !! \ingroup phys_const
  DOUBLE PRECISION, PARAMETER :: ALPHA_HBARC = ALPHA_FINE_STRUCTURE * HBARC  ! Fine-structure constant in MeV fm

CONTAINS
  !> \brief Compute the reduced mass of two particles.
  !! \ingroup phys_const
  !! \param[in] M1 Mass of first particle (MeV)
  !! \param[in] M2 Mass of second particle (MeV)
  !! \return Reduced mass (MeV)
  PURE FUNCTION REDUCED_MASS(M1, M2) RESULT(MU)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: M1, M2
    DOUBLE PRECISION :: MU

    MU = (M1 * M2) / (M1 + M2)
  END FUNCTION REDUCED_MASS

  !> \brief Compute \f$\hbar^2/(2m)\f$ for a given mass.
  !! \ingroup phys_const
  !! \param[in] MASS Mass in MeV
  !! \return Value of \f$\hbar^2/(2m)\f$ in MeV fm^2
  PURE FUNCTION HBARM(MASS)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: MASS
    DOUBLE PRECISION :: HBARM

    HBARM = HBARC**2 / (2*MASS)  ! MeV fm^2
  END FUNCTION HBARM

  !> \brief Return \f$\hbar^2/(2\mu)\f$ for nucleon-nucleon system by isospin projection.
  !! \ingroup phys_const
  !! \param[in] TZ Isospin projection (0: np, 1: pp, -1: nn)
  !! \return Value of \f$\hbar^2/(2\mu)\f$ in MeV fm^2
  FUNCTION HBARM_NUCLEON_NUCLEON(TZ) RESULT(HBARM_NN)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: TZ
    DOUBLE PRECISION :: HBARM_NN

    IF (TZ == 0) THEN
      HBARM_NN = HBARM_NP  ! Neutron-proton reduced mass
    ELSE IF (TZ == 1) THEN
      HBARM_NN = HBARM_P  ! Proton-proton reduced mass
    ELSE IF (TZ == -1) THEN
      HBARM_NN = HBARM_N  ! Neutron-neutron reduced mass
    ELSE
      STOP "Invalid TZ value in HBARM_NUCLEON_NUCLEON"
    END IF
  END FUNCTION HBARM_NUCLEON_NUCLEON

  !> \brief Set reduced mass and kinetic factor for a nucleon-nucleon system.
  !! \ingroup phys_const
  !! \param[in] TZ Isospin projection (0: np, 1: pp, -1: nn)
  !! \param[out] MASS Mass to use (MeV)
  !! \param[out] HTM Value of \f$\hbar^2/(2m)\f$ (MeV fm^2)
  SUBROUTINE SET_REDUCED_MASS_AND_HTM(TZ, MASS, HTM)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: TZ  ! 0 for np, 1 for pp, -1 for nn
    DOUBLE PRECISION, INTENT(OUT) :: MASS  ! Mass of the particle in MeV
    DOUBLE PRECISION, INTENT(OUT) :: HTM  ! HTM value in MeV fm^2

    SELECT CASE (TZ)
      CASE (0)  ! Neutron-proton
        MASS = MASS_REDUCED_NP  ! Use reduced mass for np
        HTM = HBARM_NP
      CASE (1)  ! Proton-proton
        MASS = MASS_PROTON  ! Use proton mass for pp
        HTM = HBARM_P
      CASE (-1)  ! Neutron-neutron
        MASS = MASS_NEUTRON  ! Use neutron mass for nn
        HTM = HBARM_N
      CASE DEFAULT
        STOP "Invalid TZ value in SET_REDUCED_MASS_AND_HTM"
    END SELECT
  END SUBROUTINE SET_REDUCED_MASS_AND_HTM

END MODULE PHYSICAL_CONSTANTS