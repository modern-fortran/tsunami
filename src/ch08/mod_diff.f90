module mod_diff

  ! Finite difference functions.

  use iso_fortran_env, only: int32, real32
  implicit none

  private
  public :: diffx, diffy

contains

  pure function diffx(x) result(dx)
    ! Centered finite difference in x.
    real(real32), intent(in) :: x(:,:)
    real(real32) :: dx(size(x, dim=1), size(x, dim=2))
    integer(int32) :: i, im
    im = size(x, dim=1)
    dx = 0
    dx(2:im-1,:) = 0.5 * (x(3:im,:) - x(1:im-2,:))
  end function diffx

  pure function diffy(x) result(dx)
    ! Centered finite difference in y.
    real(real32), intent(in) :: x(:,:)
    real(real32) :: dx(size(x, dim=1), size(x, dim=2))
    integer(int32) :: j, jm
    jm = size(x, dim=2)
    dx = 0
    dx(:,2:jm-1) = 0.5 * (x(:,3:jm) - x(:,1:jm-2))
  end function diffy

end module mod_diff
