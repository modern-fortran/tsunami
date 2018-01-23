module mod_diff

  ! Module that provides 

  use iso_fortran_env, only: real32

  implicit none

  private
  public :: diff

contains

  pure function diff(x) result(dx)
    ! Returns an upstream difference of a 1-d array,
    ! with periodic boundary condition.
    real(kind=real32), dimension(:), intent(in) :: x
    real(kind=real32), dimension(:), allocatable :: dx
    integer :: i, idm

    idm = size(x)
    allocate(dx(idm))

    dx(1) = x(1) - x(idm)

    do concurrent(i = 2:idm)
      dx(2:idm) = x(2:idm) - x(1:idm-1)
    end do

  end function diff

end module mod_diff
