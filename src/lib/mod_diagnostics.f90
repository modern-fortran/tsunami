module mod_diagnostics

  ! Provides various diagnostic functions.

  use mod_kinds, only: ik, rk

  implicit none

  private
  public :: ke, mean

  interface mean
    module procedure :: mean_1d, mean_2d
  end interface mean

contains

  pure elemental real(rk) function ke(u, v)
    ! Computes the kinetic energy as 1/2 (u^2 + v^2)
    real(rk), intent(in) :: u, v
    ke = 0.5_rk * sqrt(u**2 + v**2)
  end function ke

  pure real(rk) function mean_1d(x) result(mean)
    real(rk), intent(in) :: x(:)
    mean = sum(x) / size(x)
  end function mean_1d

  pure real(rk) function mean_2d(x) result(mean)
    real(rk), intent(in) :: x(:, :)
    mean = sum(x) / size(x)
  end function mean_2d

end module mod_diagnostics
