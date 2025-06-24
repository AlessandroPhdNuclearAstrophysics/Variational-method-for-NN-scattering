!> \file spherical_Bessel_functions.f90
!! \brief Fortran interface and helpers for spherical Bessel functions using GSL.
!! \defgroup gsl_bessel Spherical Bessel Functions
!! \ingroup math
!!
!! This module provides Fortran interfaces to the GNU Scientific Library (GSL) routines
!! for regular and irregular spherical Bessel functions, as well as convenient Fortran
!! wrappers for their values and derivatives.
!!
!! \author Alessandro
!! \date 2025
!!
!! \note This module requires the GSL library to be installed and linked with your Fortran project.
module gsl_bessel
  use iso_c_binding
  implicit none
  private

  !> \brief Interface to GSL's spherical Bessel functions.
  interface
    !> \brief Regular spherical Bessel function j_l(x) from GSL.
    function gsl_sf_bessel_jl(l, x) bind(C, name="gsl_sf_bessel_jl") result(gsl_sf_bessel_jl)
      import c_int, c_double
      integer(c_int), value :: l
      real(c_double), value :: x
      real(c_double) :: gsl_sf_bessel_jl
    end function

    !> \brief Irregular spherical Bessel function y_l(x) from GSL.
    function gsl_sf_bessel_yl(l, x) bind(C, name="gsl_sf_bessel_yl") result(gsl_sf_bessel_yl)
      import c_int, c_double
      integer(c_int), value :: l
      real(c_double), value :: x
      real(c_double) :: gsl_sf_bessel_yl
    end function
  end interface

  public :: spherical_j, spherical_y
  public :: spherical_jp, spherical_jpp
  public :: spherical_yp, spherical_ypp

contains

  !> \brief Compute the regular spherical Bessel function j_l(x).
  !! \ingroup gsl_bessel
  !! \param[in] l Angular momentum quantum number
  !! \param[in] x Argument
  !! \return j Value of j_l(x)
  function spherical_j(l, x) result(j)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: j
    j = gsl_sf_bessel_jl(int(l, c_int), real(x, c_double))
  end function

  !> \brief Compute the irregular spherical Bessel function y_l(x).
  !! \ingroup gsl_bessel
  !! \param[in] l Angular momentum quantum number
  !! \param[in] x Argument
  !! \return y Value of y_l(x)
  function spherical_y(l, x) result(y)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: y
    y = -gsl_sf_bessel_yl(int(l, c_int), real(x, c_double))
  end function

  !> \brief Compute the derivative of the regular spherical Bessel function j_l'(x).
  !! \ingroup gsl_bessel
  !! \param[in] l Angular momentum quantum number
  !! \param[in] x Argument
  !! \return dj Value of j_l'(x)
  function spherical_jp(l, x) result(dj)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: dj
    real(8) :: jl, jlp1

    jl = spherical_j(l, x)
    jlp1 = spherical_j(l + 1, x)
    dj = (l / x) * jl - jlp1
  end function

  !> \brief Compute the second derivative of the regular spherical Bessel function j_l''(x).
  !! \ingroup gsl_bessel
  !! \param[in] l Angular momentum quantum number
  !! \param[in] x Argument
  !! \return d2j Value of j_l''(x)
  function spherical_jpp(l, x) result(d2j)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: d2j
    real(8) :: jl, jlp1, jlp2
    real(8) :: x2, l1

    x2 = x * x
    jl = spherical_j(l, x)
    jlp1 = spherical_j(l + 1, x)
    jlp2 = spherical_j(l + 2, x)
    l1 = real(2 * l + 1, 8)

    d2j = ((l * (l + 1)) / x2 - 1.0d0) * jl - (l1 / x) * jlp1 + jlp2
  end function

  !> \brief Compute the derivative of the irregular spherical Bessel function y_l'(x).
  !! \ingroup gsl_bessel
  !! \param[in] l Angular momentum quantum number
  !! \param[in] x Argument
  !! \return dy Value of y_l'(x)
  function spherical_yp(l, x) result(dy)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: dy
    real(8) :: yl, ylp1, ylm1

    if (l==0) then
      dy = -(dcos(x) + x*dsin(x))/x**2
    else
      ylm1 = spherical_y(l-1, x)
      yl = spherical_y(l, x)
      ylp1 = spherical_y(l + 1, x)
      dy = 0.5d0*( ylm1 - ylp1 ) - yl/(2*x)
    endif
  end function

  !> \brief Compute the second derivative of the irregular spherical Bessel function y_l''(x).
  !! \ingroup gsl_bessel
  !! \param[in] l Angular momentum quantum number
  !! \param[in] x Argument
  !! \return d2y Value of y_l''(x)
  function spherical_ypp(l, x) result(d2y)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: d2y
    real(8) :: yl, ylp1
    real(8) :: x2

    x2 = x * x
    yl = spherical_y(l, x)
    ylp1 = spherical_y(l + 1, x)

    d2y = ( 2*x*ylp1 + ( l*(l-1) - x2 )*yl )/x2
  end function

end module gsl_bessel