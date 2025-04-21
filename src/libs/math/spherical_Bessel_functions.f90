module gsl_bessel
  use iso_c_binding
  implicit none

  ! Interface for GSL's spherical Bessel functions
  interface
    ! Regular spherical Bessel function j_l(x)
    function gsl_sf_bessel_jl(l, x) bind(C, name="gsl_sf_bessel_jl")
      import c_int, c_double
      integer(c_int), value :: l
      real(c_double), value :: x
      real(c_double) :: gsl_sf_bessel_jl
    end function

    ! Irregular spherical Bessel function y_l(x)
    function gsl_sf_bessel_yl(l, x) bind(C, name="gsl_sf_bessel_yl")
      import c_int, c_double
      integer(c_int), value :: l
      real(c_double), value :: x
      real(c_double) :: gsl_sf_bessel_yl
    end function
  end interface

contains

  ! Helper function to compute j_l(x) (optional, for convenience)
  function spherical_j(l, x) result(j)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: j
    j = gsl_sf_bessel_jl(int(l, c_int), real(x, c_double))
  end function

  ! Helper function to compute y_l(x) (optional, for convenience)
  function spherical_y(l, x) result(y)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: y
    y = gsl_sf_bessel_yl(int(l, c_int), real(x, c_double))
  end function

end module gsl_bessel