module mod_diff

  use iso_fortran_env, only: int32, real32
  implicit none

  private
  public :: diff

contains

  pure function diff(x) result(dx)
    ! Returns a 1st-order upstream finite difference of a 1-d array.
    real(real32), intent(in) :: x(:)
    real(real32) :: dx(size(x))
    integer(int32) :: im
    im = size(x)
    dx(1) = x(1) - x(im)
    dx(2:im) = x(2:im) - x(1:im-1)
  end function diff

end module mod_diff
