module mod_diff

  ! Module that provides finite difference functions.

  use mod_kinds, only: ik, rk

  implicit none

  private
  public :: diffc, diffc_periodic, diffu
  public :: diffc_2d_x, diffc_2d_y

contains

  pure function diffc_2d_x(x) result(dx)
    real(rk), intent(in) :: x(:,:)
    real(rk) :: dx(size(x, dim=1), size(x, dim=2))
    integer(ik) :: i, im
    im = size(x, dim=1)
    dx = 0
    dx(1,:) = 0.5 * (x(2,:) - x(im,:))
    dx(im,:) = 0.5 * (x(1,:) - x(im-1,:))
    do concurrent(i = 2:im-1)
      dx(i,:) = 0.5 * (x(i+1,:) - x(i-1,:))
    end do
  end function diffc_2d_x

  pure function diffc_2d_y(x) result(dx)
    real(rk), intent(in) :: x(:,:)
    real(rk) :: dx(size(x, dim=1), size(x, dim=2))
    integer(ik) :: j, jm
    jm = size(x, dim=2)
    dx = 0
    dx(:,1) = 0.5 * (x(:,2) - x(:,jm))
    dx(:,jm) = 0.5 * (x(:,1) - x(:,jm-1))
    do concurrent(j = 2:jm-1)
      dx(:,j) = 0.5 * (x(:,j+1) - x(:,j-1))
    end do
  end function diffc_2d_y

  pure function diffc(x) result(dx)
    ! Returns a centered difference of a 1-d array,
    ! without any boundary conditions applied.
    real(rk), intent(in) :: x(:)
    real(rk) :: dx(size(x))
    integer(ik) :: i, idm

    idm = size(x)
    dx = 0

    do concurrent(i = 2:idm-1)
      dx(i) = 0.5 * (x(i+1) - x(i-1))
    end do

  end function diffc

  pure function diffc_periodic(x) result(dx)
    ! Returns a centered difference of a 1-d array,
    ! with periodic boundary condition.
    real(rk), intent(in) :: x(:)
    real(rk) :: dx(size(x))
    integer(ik) :: i, idm

    idm = size(x)

    dx(1) = 0.5 * (x(2) - x(idm))
    dx(idm) = 0.5 * (x(1) - x(idm-1))

    do concurrent(i = 2:idm-1)
      dx(i) = 0.5 * (x(i+1) - x(i-1))
    end do

  end function diffc_periodic

  pure function diffu(x) result(dx)
    ! Returns an upstream difference of a 1-d array,
    ! with periodic boundary condition.
    real(rk), intent(in) :: x(:)
    real(rk) :: dx(size(x))
    integer(ik) :: i, idm

    idm = size(x)

    dx(1) = x(1) - x(idm)

    do concurrent(i = 2:idm)
      dx(2:idm) = x(2:idm) - x(1:idm-1)
    end do

  end function diffu

end module mod_diff
