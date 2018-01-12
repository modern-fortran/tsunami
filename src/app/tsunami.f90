program tsunami

! Tsunami simulator
!
! This version solves the non-linear 1-d shallow water equation:
!
!     du/dt + u du/dx + g dh/dx= 0
!
!     dh/dt + h_mean du/dx = 0

use iso_fortran_env, only: wip=>int32, wp=>real32, output_unit
use mod_diff, only: diffc, diffu

implicit none

integer(kind=wip) :: i, n

integer(kind=wip), parameter :: im = 100 ! grid size in x
integer(kind=wip), parameter :: nm = 1000 ! number of time steps

real(kind=wp), parameter :: dt = 0.1 ! time step [s]
real(kind=wp), parameter :: dx = 1 ! grid spacing [m]
real(kind=wp), parameter :: c = 1 ! phase speed [m/s]

real(kind=wp), parameter :: g = 9.8 ! gravitational acceleration [m/s]

real(kind=wp), dimension(im) :: du, u, h
real(kind=wp), dimension(im) :: h_mean

! initialize a gaussian blob centered at i = 25
do i = 1, im
  h(i) = exp(-0.02*(i-25)**2)
end do

u = 0

h_mean = h + 10

! write initial state to screen
write(unit=output_unit, fmt=*)0, h

time_loop: do n = 1, nm

  du = u * diffc(u, periodic=.true.) + g * diffc(h, periodic=.true.)
  u = u - du / dx * dt

  h = h - h_mean * diffc(u, periodic=.true.) / dx * dt

  ! write current state to screen
  write(unit=output_unit, fmt=*)n, h

end do time_loop

end program tsunami
