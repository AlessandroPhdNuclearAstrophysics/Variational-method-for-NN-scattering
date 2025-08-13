!> \file potential_partial_wave_caller.f90
!! \defgroup nn_potentials N-N Potentials
!! \ingroup physics
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
  PUBLIC :: SET_POTENTIAL_PARAMETERS
  PUBLIC :: POT_PW_RVALUES


  !> \brief Generic interface for partial-wave potential calculation.
  !! \ingroup nn_potentials
  !! This interface allows the user to call POT_PW with different argument lists:
  !! - All parameters specified individually (see POT_PW_PARAMS_ALL)
  !! - Using a parameter structure (see POT_PW_PARAMS)
  !! - Using a parameter structure and a channel object (see POT_PW_PARAMS_CHANNEL)
  INTERFACE POT_PW
    MODULE PROCEDURE POT_PW_PARAMS_ALL
    MODULE PROCEDURE POT_PW_PARAMS
    MODULE PROCEDURE POT_PW_PARAMS_CHANNEL
  END INTERFACE

  !> Structure containing all parameters needed to specify a nuclear potential.
  !! \ingroup nn_potentials
  TYPE:: POTENTIAL_PARAMETERS
    INTEGER :: POT_MODEL      = 0   !< Potential identifier (14=AV14, 18=AV18, 19=EFT_PLESS, 21=EFT_PLESS_FITTED)
    INTEGER :: POT_SUBMODEL   = 0   !< Channel index (used by some potentials)
    INTEGER :: EM_INTERACTION = -1  !< Type of electromagnetic interaction (0=Coulomb, >0=Check AV18 EMPOT)
  END TYPE POTENTIAL_PARAMETERS

CONTAINS

  !> Computes the partial-wave nuclear potential matrix for different models.
  !! \ingroup nn_potentials
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
    CASE DEFAULT
      STOP "ERROR, POTENTIAL NOT FOUND"
    END SELECT
    RETURN
  END SUBROUTINE POT_PW_PARAMS_ALL

  !> Computes the partial-wave potential matrix using a parameter structure.
  !! \ingroup nn_potentials
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
  !! \ingroup nn_potentials
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

  !> Computes the potential matrix for a set of R values for a given channel.
  !! \ingroup nn_potentials
  !! \param[in] POTENTIAL_PARAMS Structure containing all potential parameters
  !! \param[in] CHANNEL Scattering channel structure
  !! \param[in] RVALS Array of interparticle distances (fm)
  !! \param[out] VPW Array of potential matrices, shape (size(RVALS), NCH, NCH)
  SUBROUTINE POT_PW_RVALUES(POTENTIAL_PARAMS, CHANNEL, RVALS, VPW)
    IMPLICIT NONE
    TYPE(POTENTIAL_PARAMETERS), INTENT(IN) :: POTENTIAL_PARAMS
    TYPE(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    DOUBLE PRECISION, INTENT(IN) :: RVALS(:)
    DOUBLE PRECISION, INTENT(OUT) :: VPW(:,:,:)

    INTEGER :: I
    INTEGER :: NCH
    DOUBLE PRECISION :: VPW_(2,2)

    ! Check that the first dimension of VPW matches the size of RVALS
    IF (SIZE(VPW,1) /= SIZE(RVALS)) THEN
      STOP "ERROR: First dimension of VPW must match size of RVALS"
    END IF

    ! Check that the last two dimensions of VPW are equal and either 1 or 2
    NCH = SIZE(VPW,2)
    IF (NCH /= SIZE(VPW,3)) THEN
      STOP "ERROR: Last two dimensions of VPW must be equal"
    END IF
    IF ((NCH /= 1) .AND. (NCH /= 2)) THEN
      STOP "ERROR: Last two dimensions of VPW must be size 1 or 2"
    END IF

    DO I = 1, SIZE(RVALS)
      CALL POT_PW_PARAMS_CHANNEL(POTENTIAL_PARAMS, CHANNEL, RVALS(I), VPW_)
      VPW(I,:,:) = VPW_
    END DO
  END SUBROUTINE POT_PW_RVALUES

  !> Sets the potential parameters structure.
  !! \ingroup nn_potentials
  !! \param[in] IPOT Potential identifier (14=AV14, 18=AV18, 19=EFT_PLESS)
  !! \param[in] ILB Channel index (used by some potentials)
  !! \param[in] LEMP Type of electromagnetic interaction (0=Coulomb, >0=see AV18 EMPOT)
  !! \param[out] POTENTIAL_PARAMS Structure to be filled with potential parameters
  !! This subroutine initializes the potential parameters structure with the given values.
  !! It sets the potential model, submodel, and electromagnetic interaction type.
  SUBROUTINE SET_POTENTIAL_PARAMETERS(IPOT, ILB, LEMP, POTENTIAL_PARAMS)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IPOT, ILB, LEMP
    TYPE(POTENTIAL_PARAMETERS), INTENT(OUT) :: POTENTIAL_PARAMS

    POTENTIAL_PARAMS%POT_MODEL = IPOT
    POTENTIAL_PARAMS%POT_SUBMODEL = ILB
    POTENTIAL_PARAMS%EM_INTERACTION = LEMP
  END SUBROUTINE SET_POTENTIAL_PARAMETERS

END MODULE POTENTIALS