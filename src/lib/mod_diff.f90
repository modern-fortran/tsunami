module mod_diff
use iso_fortran_env, only: real32

implicit none

private

public :: diffu, diffc

contains


pure function diffu(x, periodic) result(dx)
  !! Returns a upstream difference of a 1-d array. Applies periodic boundary
  !! conditions if `periodic .eqv. .true.`.
  real(kind=real32), dimension(:), intent(in) :: x
    !! Input array
  logical, intent(in), optional :: periodic
    !! If true, periodic boundary conditions area applied. Default false.
  real(kind=real32), dimension(:), allocatable :: dx
  logical :: periodic_
  integer :: idm

  if (present(periodic)) then
    periodic_ = periodic
  else
    periodic_ = .false.
  endif

  idm = size(x)
  allocate(dx(idm))

  if(idm == 0)then
    return
  elseif(idm == 1)then
    dx = 0
    return
  endif

  dx(2:idm) = x(2:idm)-x(1:idm-1)

  if (periodic_) then
    dx(1) = x(1)-x(idm)
  else
    dx(1) = dx(2)
  endif

end function diffu


pure function diffc(x, periodic) result(dx)
  !! Returns a centered difference of a 1-d array. Applies periodic boundary
  !! conditions if `periodic .eqv. .true.`.
  real(kind=real32), dimension(:), intent(in) :: x
    !! Input array
  logical, intent(in), optional :: periodic
    !! If true, periodic boundary conditions area applied. Default false.
  real(kind=real32), dimension(:), allocatable :: dx
  logical :: periodic_
  integer :: idm

  if (present(periodic)) then
    periodic_ = periodic
  else
    periodic_ = .false.
  endif

  idm = size(x)
  allocate(dx(idm))

  if(idm == 0)then
    return
  elseif(idm == 1)then
    dx = 0
    return
  endif

  dx(2:idm-1) = 0.5*(x(3:idm)-x(1:idm-2))

  if (periodic_) then
    dx(1) = 0.5*(x(2)-x(idm))
    dx(idm) = 0.5*(x(1)-x(idm-1))
  else
    dx(1) = x(2)-x(1)
    dx(idm) = x(idm)-x(idm-1)
  endif

end function diffc


end module mod_diff
