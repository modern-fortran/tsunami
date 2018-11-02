module mod_diff

  ! Finite difference functions.

  use mod_kinds, only: ik, rk

  implicit none

  private
  public :: diffx, diffy

contains

  pure function diffx(x) result(dx)
    ! Centered finite difference in x.
    real(rk), intent(in) :: x(:,:)
    real(rk) :: dx(size(x, dim=1), size(x, dim=2))
    integer(ik) :: i, im
    im = size(x, dim=1)
    dx = 0
    do concurrent(i = 2:im-1)
      dx(i,:) = 0.5 * (x(i+1,:) - x(i-1,:))
    end do
  end function diffx

  pure function diffy(x) result(dx)
    ! Centered finite difference in y.
    real(rk), intent(in) :: x(:,:)
    real(rk) :: dx(size(x, dim=1), size(x, dim=2))
    integer(ik) :: j, jm
    jm = size(x, dim=2)
    dx = 0
    do concurrent(j = 2:jm-1)
      dx(:,j) = 0.5 * (x(:,j+1) - x(:,j-1))
    end do
  end function diffy

end module mod_diff
