!> \file potential_pw.f90
!! \brief Interface for the calculation of partial-wave nuclear potentials.
!!
!! This module provides interfaces and routines to compute the nucleon-nucleon
!! potential matrix for various models (AV14, AV18, EFT_PLESS, ...).
!!
!! \author Alessandro Grassi
!! \date 2025

MODULE POTENTIALS
  USE QUANTUM_NUMBERS
  IMPLICIT NONE

  PUBLIC :: POT_PW
  PUBLIC :: POT_PW_PARAMS_ALL
  PUBLIC :: POT_PW_PARAMS
  PUBLIC :: POT_PW_PARAMS_CHANNEL
  PUBLIC :: POTENTIAL_PARAMETERS

  INTERFACE POT_PW
    MODULE PROCEDURE POT_PW_PARAMS_ALL
    MODULE PROCEDURE POT_PW_PARAMS
    MODULE PROCEDURE POT_PW_PARAMS_CHANNEL
  END INTERFACE

  !> Structure containing all parameters needed to specify a nuclear potential.
  TYPE:: POTENTIAL_PARAMETERS
    INTEGER :: POT_MODEL      = 0   !< Potential identifier (14=AV14, 18=AV18, 19=EFT_PLESS, 21=EFT_PLESS_FITTED)
    INTEGER :: POT_SUBMODEL   = 0   !< Channel index (used by some potentials)
    INTEGER :: EM_INTERACTION = -1  !< Type of electromagnetic interaction (0=Coulomb, >0=Check AV18 EMPOT)
  END TYPE POTENTIAL_PARAMETERS

CONTAINS

  !> Computes the partial-wave nuclear potential matrix for different models.
  !! Selects the desired potential model based on IPOT and calls the corresponding routine.
  !! \param[in] IPOT Potential identifier (14=AV14, 18=AV18, 19=EFT_PLESS, 21=EFT_PLESS_FITTED)
  !! \param[in] ILB Channel index (used by some potentials)
  !! \param[in] LEMP Type of electromagnetic interaction (0=none, >0=present)
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Total spin
  !! \param[in] J Total angular momentum
  !! \param[in] T1Z Isospin projection nucleon 1 (+1=proton, -1=neutron)
  !! \param[in] T2Z Isospin projection nucleon 2 (+1=proton, -1=neutron)
  !! \param[in] R Interparticle distance (fm)
  !! \param[out] VPW 2x2 potential matrix in coupled basis (mixed channels)
  SUBROUTINE POT_PW_PARAMS_ALL(IPOT, ILB, LEMP, L, S, J, T1Z, T2Z, R, VPW)
    USE AV18
    USE AV14
    USE EFT_PLESS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IPOT, LEMP, L, S, J, T1Z, T2Z
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT):: VPW(2,2)
    INTEGER, INTENT(IN) :: ILB
  !
    INTEGER :: T
  !
    T = T_FROM_L_S(L, S)
    VPW = 0.D0
  !
    SELECT CASE (IPOT)
    CASE (14)
      CALL AV14PW(LEMP, L, S, J, T1Z, T2Z, R, VPW)
      RETURN
    CASE (18)
      CALL AV18PW(ILB, L, S, J, T, T1Z, T2Z, R, VPW, LEMP)
      RETURN
    CASE (19)
      CALL EFT_PLESS_PW(ILB, L, S, J, T1Z, T2Z, R, VPW, LEMP)
      RETURN
    CASE (21)
      CALL EFT_PLESS_PW_FITTED(ILB, R, L, S, J, T1Z + T2Z, VPW)
    CASE DEFAULT
      STOP "ERROR, POTENTIAL NOT FOUND"
    END SELECT
    RETURN
  END SUBROUTINE POT_PW_PARAMS_ALL

  !> Computes the partial-wave potential matrix using a parameter structure.
  !! \param[in] POTENTIAL_PARAMS Structure containing all potential parameters
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Total spin
  !! \param[in] J Total angular momentum
  !! \param[in] T1Z Isospin projection nucleon 1
  !! \param[in] T2Z Isospin projection nucleon 2
  !! \param[in] R Interparticle distance (fm)
  !! \param[out] VPW 2x2 potential matrix
  SUBROUTINE POT_PW_PARAMS(POTENTIAL_PARAMS, L, S, J, T1Z, T2Z, R, VPW)
    IMPLICIT NONE
    TYPE(POTENTIAL_PARAMETERS), INTENT(IN) :: POTENTIAL_PARAMS
    INTEGER, INTENT(IN) :: L, S, J, T1Z, T2Z
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VPW(2,2)

    CALL POT_PW_PARAMS_ALL(POTENTIAL_PARAMS%POT_MODEL, POTENTIAL_PARAMS%POT_SUBMODEL, POTENTIAL_PARAMS%EM_INTERACTION, &
             L, S, J, T1Z, T2Z, R, VPW)
  END SUBROUTINE POT_PW_PARAMS

  !> Computes the partial-wave potential matrix for a given scattering channel.
  !! \param[in] POTENTIAL_PARAMS Structure containing all potential parameters
  !! \param[in] CHANNEL Scattering channel structure
  !! \param[in] R Interparticle distance (fm)
  !! \param[out] VPW 2x2 potential matrix
  SUBROUTINE POT_PW_PARAMS_CHANNEL(POTENTIAL_PARAMS, CHANNEL, R, VPW)
    IMPLICIT NONE
    TYPE(POTENTIAL_PARAMETERS), INTENT(IN) :: POTENTIAL_PARAMS
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT) :: VPW(2,2)
    INTEGER :: L, S, J, T1Z, T2Z, TZ, ICH
    DOUBLE PRECISION :: VPW_(2,2)
    
    IF (IS_CHANNEL_COUPLED(CHANNEL)) THEN
      J = GET_CHANNEL_J(CHANNEL)
      TZ= GET_CHANNEL_TZ(CHANNEL)
      CALL TZ_TO_T1Z_T2Z(TZ, T1Z, T2Z)
      L = GET_CHANNEL_L(CHANNEL,1)
      S = GET_CHANNEL_S(CHANNEL,1)
      CALL POT_PW_PARAMS_ALL(POTENTIAL_PARAMS%POT_MODEL, POTENTIAL_PARAMS%POT_SUBMODEL, POTENTIAL_PARAMS%EM_INTERACTION, &
                  L, S, J, T1Z, T2Z, R, VPW)
    ELSE
      VPW_ = 0.D0
      DO ICH = 1, GET_CHANNEL_NCH(CHANNEL)
        J = GET_CHANNEL_J(CHANNEL)
        TZ = GET_CHANNEL_TZ(CHANNEL)
        CALL TZ_TO_T1Z_T2Z(TZ, T1Z, T2Z)
        L = GET_CHANNEL_L(CHANNEL, ICH)
        S = GET_CHANNEL_S(CHANNEL, ICH)
        CALL POT_PW_PARAMS_ALL(POTENTIAL_PARAMS%POT_MODEL, POTENTIAL_PARAMS%POT_SUBMODEL, POTENTIAL_PARAMS%EM_INTERACTION, &
                 L, S, J, T1Z, T2Z, R, VPW_)
        VPW(ICH, ICH) = VPW_(1,1)
      END DO
    END IF
  END SUBROUTINE POT_PW_PARAMS_CHANNEL
END MODULE POTENTIALS