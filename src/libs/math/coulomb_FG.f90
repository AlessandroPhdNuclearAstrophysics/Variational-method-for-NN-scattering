module gsl_coulomb
  use, intrinsic :: iso_c_binding
  implicit none
  private

  ! Exposed types and procedures
  public :: coulomb_wave_FG, gsl_sf_result

  ! GSL result type (matches C struct)
  type, bind(C), public :: gsl_sf_result
    real(c_double) :: val  ! Function value
    real(c_double) :: err  ! Absolute error estimate
  end type gsl_sf_result

  ! Interface to GSL's coulomb_wave_FG_e
  interface
    function gsl_sf_coulomb_wave_FG_e(eta, x, L, k, F, FP, G, GP, F_exp, G_exp) &
        bind(C, name='gsl_sf_coulomb_wave_FG_e')
      import :: c_int, c_double, gsl_sf_result
      integer(c_int) :: gsl_sf_coulomb_wave_FG_e  ! Return status (0 = success)
      real(c_double), value :: eta, x, L          ! Input parameters
      integer(c_int), value :: k                  ! Boundary condition flag
      type(gsl_sf_result) :: F, FP, G, GP        ! Output structures
      real(c_double) :: F_exp, G_exp              ! Scaling exponents
    end function
  end interface

contains

  ! Wrapper subroutine that returns complete result structures
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