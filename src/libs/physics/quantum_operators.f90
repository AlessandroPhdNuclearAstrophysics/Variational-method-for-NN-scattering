!> \file quantum_operators.f90
!! \defgroup quantum_op Quantum Operators
!! \ingroup physics
!! \brief Module providing quantum mechanical operators for nuclear physics calculations
!!
!! This module implements various quantum mechanical operators commonly used in nuclear physics,
!! including tensor (S12), spin-orbit (LS), orbital angular momentum squared (L²),
!! and isospin (T12) operators.
!!
!! \author Alessandro
!! \date 2025

MODULE QUANTUM_OPERATORS
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  PRIVATE

  !> \brief Interface for the S12 operator.
  !! \ingroup quantum_op
  !! This interface provides overloaded procedures for the S12 operator,
  !! allowing usage with either quantum numbers or channel representations.
  !!
  !! \par Overloads:
  !! - S12_OPERATOR_QUANTUM_NUMBERS(L, S, J): 
  !!     \arg[in] L Orbital angular momentum (integer)
  !!     \arg[in] S Spin (integer)
  !!     \arg[in] J Total angular momentum (integer)
  !! - S12_OPERATOR_CHANNEL(CHANNEL): 
  !!     \arg[in] CHANNEL Scattering channel object (TYPE(SCATTERING_CHANNEL))
  INTERFACE S12_OPERATOR
    MODULE PROCEDURE S12_OPERATOR_QUANTUM_NUMBERS
    MODULE PROCEDURE S12_OPERATOR_CHANNEL
  END INTERFACE

  !> \brief Interface for the LS operator.
  !! \ingroup quantum_op
  !! This interface provides overloaded procedures for the LS operator,
  !! supporting both quantum numbers and channel representations.
  !!
  !! \par Overloads:
  !! - LS_OPERATOR_QUANTUM_NUMBERS(LMIN, S, J): 
  !!     \arg[in] LMIN Minimum orbital angular momentum (integer)
  !!     \arg[in] S Spin (integer)
  !!     \arg[in] J Total angular momentum (integer)
  !! - LS_OPERATOR_CHANNEL(CHANNEL): 
  !!     \arg[in] CHANNEL Scattering channel object (TYPE(SCATTERING_CHANNEL))
  INTERFACE LS_OPERATOR
    MODULE PROCEDURE LS_OPERATOR_QUANTUM_NUMBERS
    MODULE PROCEDURE LS_OPERATOR_CHANNEL
  END INTERFACE

  !> \brief Interface for the L2 operator.
  !! \ingroup quantum_op
  !! This interface provides overloaded procedures for the L2 operator,
  !! supporting both quantum numbers and channel representations.
  !!
  !! \par Overloads:
  !! - L2_OPERATOR_QUANTUM_NUMBERS(L, COUPLED): 
  !!     \arg[in] L Orbital angular momentum (integer)
  !!     \arg[in] COUPLED Logical flag for coupled channels (optional, logical)
  !! - L2_OPERATOR_CHANNEL(CHANNEL): 
  !!     \arg[in] CHANNEL Scattering channel object (TYPE(SCATTERING_CHANNEL))
  INTERFACE L2_OPERATOR
    MODULE PROCEDURE L2_OPERATOR_QUANTUM_NUMBERS
    MODULE PROCEDURE L2_OPERATOR_CHANNEL
  END INTERFACE

  !> \brief Interface for the T12 operator.
  !! \ingroup quantum_op
  !! This interface provides overloaded procedures for the T12 operator,
  !! supporting both quantum numbers and channel representations.
  !!
  !! \par Overloads:
  !! - T12_OPERATOR_QUANTUM_NUMBERS(T, TZ): 
  !!     \arg[in] T Total isospin (integer)
  !!     \arg[in] TZ Isospin projection (integer)
  !! - T12_OPERATOR_CHANNEL(CHANNEL): 
  !!     \arg[in] CHANNEL Scattering channel object (TYPE(SCATTERING_CHANNEL))
  INTERFACE T12_OPERATOR
    MODULE PROCEDURE T12_OPERATOR_QUANTUM_NUMBERS
    MODULE PROCEDURE T12_OPERATOR_CHANNEL
  END INTERFACE

  ! Public interfaces for the operators
  PUBLIC :: IDENTITY_MATRIX
  PUBLIC :: S12_OPERATOR, LS_OPERATOR, L2_OPERATOR, T12_OPERATOR
  PRIVATE:: S12_OPERATOR_QUANTUM_NUMBERS, S12_OPERATOR_CHANNEL
  PRIVATE:: LS_OPERATOR_QUANTUM_NUMBERS, LS_OPERATOR_CHANNEL
  PRIVATE:: L2_OPERATOR_QUANTUM_NUMBERS, L2_OPERATOR_CHANNEL
  PRIVATE:: T12_OPERATOR_QUANTUM_NUMBERS, T12_OPERATOR_CHANNEL

CONTAINS

  !> \brief Identity operator
  !! \ingroup quantum_op
  !! \param[in] N Dimension of the identity matrix (default 2)
  !! \return Identity matrix of dimension N×N
  FUNCTION IDENTITY_MATRIX(N) RESULT(I_MAT)
    IMPLICIT NONE
    INTEGER, INTENT(IN), OPTIONAL :: N
    INTEGER :: DIM, I
    DOUBLE PRECISION, ALLOCATABLE :: I_MAT(:,:)
    
    IF (PRESENT(N)) THEN
      DIM = N
    ELSE
      DIM = 2  ! Default dimension
    END IF
    ALLOCATE(I_MAT(DIM,DIM))
    I_MAT = 0.0D0
    DO I = 1, DIM
      I_MAT(I,I) = 1.0D0
    END DO
  END FUNCTION IDENTITY_MATRIX

  !> \brief Tensor operator S₁₂ = 3(σ₁·r̂)(σ₂·r̂) - σ₁·σ₂
  !! \ingroup quantum_op
  !! \param[in] L Orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] J Total angular momentum
  !! \return S₁₂ matrix representation in coupled |LSJ> basis
  FUNCTION S12_OPERATOR_QUANTUM_NUMBERS(L, S, J) RESULT(S12_RES)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L, S, J
    DOUBLE PRECISION :: S12_RES(2,2)

    S12_RES = 0.0D0
    
    ! If S=0, tensor force doesn't contribute
    IF (S==0) RETURN
    
    IF (L==J) THEN
      ! Uncoupled channel
      S12_RES(1,1) = 2.0D0
      RETURN
    ELSEIF (ABS(J-L)==1) THEN
      IF (J==0) THEN
        ! Special case for J=0
        IF (J==L+1) THEN
          S12_RES(1,1) = -2.0D0*(J-1.0D0)/(2.0D0*J+1.0D0)
        ELSEIF (J==L-1) THEN
          S12_RES(1,1) = -2.0D0*(J+2.0D0)/(2.0D0*J+1.0D0)
        ELSE
          WRITE(*,*) "Error in S12_OPERATOR_QUANTUM_NUMBERS: Invalid LSJ combination"
          STOP
        ENDIF
      ELSE
        ! Coupled channels case
        S12_RES(1,1) = -2.0D0*(J-1.0D0)/(2.0D0*J+1.0D0)
        S12_RES(1,2) = 6.0D0*SQRT(J*(J+1.0D0))/(2.0D0*J+1.0D0)
        S12_RES(2,1) = 6.0D0*SQRT(J*(J+1.0D0))/(2.0D0*J+1.0D0)
        S12_RES(2,2) = -2.0D0*(J+2.0D0)/(2.0D0*J+1.0D0)
      ENDIF
    ELSE
      WRITE(*,*) "Error in S12_OPERATOR_QUANTUM_NUMBERS: Invalid LSJ combination"
      STOP
    ENDIF
  END FUNCTION S12_OPERATOR_QUANTUM_NUMBERS

  !> \brief Tensor operator S₁₂ for a specific scattering channel
  !! \ingroup quantum_op
  !! \param[in] CHANNEL Scattering channel object
  !! \return S₁₂ matrix representation for the given channel
  FUNCTION S12_OPERATOR_CHANNEL(CHANNEL) RESULT(S12_RES)
    IMPLICIT NONE
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    DOUBLE PRECISION :: S12_RES(2,2)
    INTEGER :: L, S, J

    L = CHANNEL%L(1)
    S = CHANNEL%S(1)
    J = CHANNEL%J()
    
    S12_RES = S12_OPERATOR_QUANTUM_NUMBERS(L, S, J)
  END FUNCTION S12_OPERATOR_CHANNEL

  !> \brief Spin-orbit operator L·S
  !! \ingroup quantum_op
  !! \param[in] LMIN Minimum orbital angular momentum
  !! \param[in] S Spin
  !! \param[in] J Total angular momentum
  !! \return L·S matrix representation in coupled |LSJ> basis
  FUNCTION LS_OPERATOR_QUANTUM_NUMBERS(LMIN, S, J) RESULT(LS_OP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LMIN, S, J
    INTEGER :: LS_OP(2,2)

    LS_OP = 0
    
    IF (LMIN==J .OR. J==0) THEN
      ! Uncoupled channel or special case
      LS_OP(1,1) = (J*(J+1)-S*(S+1)-LMIN*(LMIN+1))/2
    ELSEIF ((J-LMIN)==1 .AND. S==1) THEN
      ! Coupled channels case with S=1
      LS_OP(1,1) = (J*(J+1)-S*(S+1)-LMIN*(LMIN+1))/2
      LS_OP(2,2) = (J*(J+1)-S*(S+1)-(LMIN+2)*(LMIN+3))/2
    ELSE
      WRITE(*,*) "Error in LS_OPERATOR_QUANTUM_NUMBERS: Invalid LSJ combination"
      STOP
    ENDIF
  END FUNCTION LS_OPERATOR_QUANTUM_NUMBERS

  !> \brief Spin-orbit operator L·S for a specific scattering channel
  !! \ingroup quantum_op
  !! \param[in] CHANNEL Scattering channel object
  !! \return L·S matrix representation for the given channel
  FUNCTION LS_OPERATOR_CHANNEL(CHANNEL) RESULT(LS_OP)
    IMPLICIT NONE
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: LS_OP(2,2)
    INTEGER :: L, S, J
    
    L = CHANNEL%L(1)
    S = CHANNEL%S(1)
    J = CHANNEL%J()
    
    LS_OP = LS_OPERATOR_QUANTUM_NUMBERS(L, S, J)
  END FUNCTION LS_OPERATOR_CHANNEL

  !> \brief Orbital angular momentum squared operator L²
  !! \ingroup quantum_op
  !! \param[in] L Orbital angular momentum
  !! \param[in] COUPLED Flag indicating if channels are coupled
  !! \return L² matrix representation in |L> or coupled |LSJ> basis
  FUNCTION L2_OPERATOR_QUANTUM_NUMBERS(L, COUPLED) RESULT(L2_MAT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    LOGICAL, INTENT(IN), OPTIONAL :: COUPLED
    INTEGER :: L2_MAT(2,2)
    LOGICAL :: IS_COUPLED
    
    L2_MAT = 0
    
    IF (PRESENT(COUPLED)) THEN
      IS_COUPLED = COUPLED
    ELSE
      IS_COUPLED = .FALSE.
    END IF
    
    L2_MAT(1,1) = L*(L+1)
    
    IF (IS_COUPLED) L2_MAT(2,2) = (L+2)*(L+3)
  END FUNCTION L2_OPERATOR_QUANTUM_NUMBERS

  !> \brief L² operator for a specific scattering channel
  !! \ingroup quantum_op
  !! \param[in] CHANNEL Scattering channel object
  !! \return L² matrix representation for the given channel
  FUNCTION L2_OPERATOR_CHANNEL(CHANNEL) RESULT(L2_MAT)
    IMPLICIT NONE
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: L2_MAT(2,2)
    INTEGER :: L
    LOGICAL :: COUPLED

    L = CHANNEL%L(1)
    COUPLED = CHANNEL%IS_COUPLED()
    L2_MAT = L2_OPERATOR_QUANTUM_NUMBERS(L, COUPLED)
  END FUNCTION L2_OPERATOR_CHANNEL

  !> \brief Isospin operator T₁₂ = τ₁·τ₂
  !! \ingroup quantum_op
  !! \param[in] T Total isospin
  !! \param[in] TZ Isospin projection
  !! \return T₁₂ value (scalar for two-nucleon system)
  FUNCTION T12_OPERATOR_QUANTUM_NUMBERS(T, TZ) RESULT(CD)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: T, TZ
    INTEGER :: CD

    CD = 0
    
    IF (T==0) RETURN
    
    IF (ABS(TZ)==1) THEN
      ! Tz=1 gives T₁₂ = 2
      CD = 2
      ELSEIF (TZ==0) THEN
      ! Tz=0 gives T₁₂ = -4
      CD = -4
    ELSE
      WRITE(*,*) "Error in T12_OPERATOR_QUANTUM_NUMBERS: Invalid T,TZ combination"
    ENDIF
  END FUNCTION T12_OPERATOR_QUANTUM_NUMBERS

  !> \brief Isospin operator T₁₂ for a specific scattering channel
  !! \ingroup quantum_op
  !! \param[in] CHANNEL Scattering channel object
  !! \return T₁₂ value for the given channel
  FUNCTION T12_OPERATOR_CHANNEL(CHANNEL) RESULT(CD)
    IMPLICIT NONE
    CLASS(SCATTERING_CHANNEL), INTENT(IN) :: CHANNEL
    INTEGER :: CD
    INTEGER :: T, TZ

    T = CHANNEL%T(1)
    TZ = CHANNEL%TZ()
    
    CD = T12_OPERATOR_QUANTUM_NUMBERS(T, TZ)
  END FUNCTION T12_OPERATOR_CHANNEL

END MODULE QUANTUM_OPERATORS
