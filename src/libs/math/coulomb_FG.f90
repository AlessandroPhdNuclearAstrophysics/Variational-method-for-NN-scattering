!> \file coulomb_FG.f90
!! \brief Fortran interface to GSL Coulomb wave functions F and G.
!!
!! This module provides a Fortran interface to the GNU Scientific Library (GSL)
!! routines for computing regular (F) and irregular (G) Coulomb wave functions,
!! as well as their derivatives, using the C API via iso_c_binding.
!!
!! \author Alessandro
!! \date 2025

module gsl_coulomb
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> \brief Expose main types and procedures.
  public :: coulomb_wave_FG, gsl_sf_result

  !> \brief Fortran type matching GSL's gsl_sf_result struct.
  !! Holds a value and an absolute error estimate.
  type, bind(C), public :: gsl_sf_result
    real(c_double) :: val  !< Function value
    real(c_double) :: err  !< Absolute error estimate
  end type gsl_sf_result

  !> \brief Interface to GSL's coulomb_wave_FG_e C function.
  interface
    function gsl_sf_coulomb_wave_FG_e(eta, x, L, k, F, FP, G, GP, F_exp, G_exp) &
        bind(C, name='gsl_sf_coulomb_wave_FG_e')
      import :: c_int, c_double, gsl_sf_result
      integer(c_int) :: gsl_sf_coulomb_wave_FG_e  !< Return status (0 = success)
      real(c_double), value :: eta, x, L          !< Input parameters
      integer(c_int), value :: k                  !< Boundary condition flag
      type(gsl_sf_result) :: F, FP, G, GP         !< Output structures
      real(c_double) :: F_exp, G_exp              !< Scaling exponents
    end function
  end interface

  private :: gsl_sf_coulomb_wave_FG_e  !< C function to compute Coulomb wave functions
contains

  !> \brief Compute regular and irregular Coulomb wave functions and their derivatives.
  !!
  !! This is a Fortran wrapper for GSL's gsl_sf_coulomb_wave_FG_e.
  !! Computes F_L(η,ρ), F'_L(η,ρ), G_L(η,ρ), G'_L(η,ρ) and their scaling exponents.
  !!
  !! \param[in]  eta    Sommerfeld parameter (η)
  !! \param[in]  x      Radial coordinate (ρ = kr)
  !! \param[in]  L      Angular momentum quantum number
  !! \param[out] F      Regular Coulomb function F_L(η,ρ) and error
  !! \param[out] FP     Derivative F'_L(η,ρ) and error
  !! \param[out] G      Irregular Coulomb function G_L(η,ρ) and error
  !! \param[out] GP     Derivative G'_L(η,ρ) and error
  !! \param[out] F_exp  Exponent for F underflow/overflow
  !! \param[out] G_exp  Exponent for G underflow/overflow
  !! \param[out] status (optional) GSL status code (0=success)
  subroutine coulomb_wave_FG(eta, x, L, F, FP, G, GP, F_exp, G_exp, status)
    real(c_double), intent(in) :: eta            ! [in] Sommerfeld parameter (η)
    real(c_double), intent(in) :: x              ! [in] Radial coordinate (ρ = kr)
    integer, intent(in) :: L                     ! [in] Angular momentum quantum number
    type(gsl_sf_result), intent(out) :: F        ! [out] Regular solution F_L(η,ρ) and error
    type(gsl_sf_result), intent(out) :: FP       ! [out] Derivative F'_L(η,ρ) and error
    type(gsl_sf_result), intent(out) :: G        ! [out] Irregular solution G_L(η,ρ) and error
    type(gsl_sf_result), intent(out) :: GP       ! [out] Derivative G'_L(η,ρ) and error
    real(c_double), intent(out) :: F_exp         ! [out] Exponent for F underflow/overflow
    real(c_double), intent(out) :: G_exp         ! [out] Exponent for G underflow/overflow
    integer, intent(out), optional :: status     ! [out] Optional error status

    integer(c_int) :: k, stat

    k = 0  ! Default boundary condition (GSL recommended)

    stat = gsl_sf_coulomb_wave_FG_e( &
          eta,               &  ! Sommerfeld parameter
          x,                 &  ! Radial coordinate
          real(L, c_double), &  ! Angular momentum (converted to double)
          k,                 &  ! Boundary condition flag
          F,                 &  ! Regular Coulomb function F
          FP,                &  ! Derivative of F
          G,                 &  ! Irregular Coulomb function G
          GP,                &  ! Derivative of G
          F_exp,             &  ! Exponent for F scaling
          G_exp              &  ! Exponent for G scaling
    )

    if (present(status)) status = stat
  end subroutine coulomb_wave_FG

end module gsl_coulomb