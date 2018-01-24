program tsunami

! Tsunami simulator
!
! This version solves the non-linear 1-d shallow water equation:
!
!     du/dt + u du/dx + g dh/dx= 0
!
!     dh/dt + d(hu)/dx = 0

use iso_fortran_env, only: int32, real32, output_unit
use mod_diff, only: diff => diffc

implicit none

integer(kind=int32) :: i, n

integer(kind=int32), parameter :: im = 100 ! grid size in x
integer(kind=int32), parameter :: nm = 5000 ! number of time steps

real(kind=real32), parameter :: dt = 0.02 ! time step [s]
real(kind=real32), parameter :: dx = 1 ! grid spacing [m]

real(kind=real32), parameter :: g = 9.8 ! gravitational acceleration [m/s]

real(kind=real32), dimension(im) :: h, hmean, u

integer(kind=int32), parameter :: ipos = 25
real(kind=real32), parameter :: decay = 0.02

! initialize a gaussian blob centered at i = 25
do i = 1, im
  h(i) = exp(-decay * (i - ipos)**2)
end do

! set initial velocity to zero
u = 0

! set mean water depth to 10 m
hmean = 10

! write initial state to screen
write(unit=output_unit, fmt=*) 0, h

time_loop: do n = 1, nm

  ! compute u at next time step
  u = u - (u * diff(u) + g * diff(h)) / dx * dt

  ! compute h at next time step
  h = h - diff(u * (hmean + h)) / dx * dt

  ! write current state to screen
  write(unit=output_unit, fmt=*) n, h

end do time_loop

end program tsunami
