program tsunami

  ! This version solves the linearized 1-d advection equation:
  !
  !     du/dt + c du/dx = 0
  !
  ! The initial conditions and the finite difference calculation 
  ! are abstracted as external procedures.

  use iso_fortran_env, only: int32, real32
  use mod_diff, only: diff => diff_centered
  use mod_initial, only: set_gaussian

  implicit none

  integer(int32) :: n

  integer(int32), parameter :: im = 100 ! grid size in x
  integer(int32), parameter :: nm = 5000 ! number of time steps

  real(real32), parameter :: dt = 0.02 ! time step [s]
  real(real32), parameter :: dx = 1 ! grid spacing [m]
  real(real32), parameter :: g = 9.8 ! gravitational acceleration [m/s^2]

  real(real32) :: h(im), hmean(im), u(im)

  integer(int32), parameter :: icenter = 25
  real(real32), parameter :: decay = 0.02

  ! initialize water height to a Gaussian blob
  call set_gaussian(h, icenter, decay)

  ! initialize water velocity to zero
  u = 0

  ! initialize mean water depth to 10 m
  hmean = 10

  ! write initial state to screen
  print *, 0, h

  time_loop: do n = 1, nm

    ! compute u at next time step
    u = u - (u * diff(u) + g * diff(h)) / dx * dt

    ! compute h at next time step
    h = h - diff(u * (hmean + h)) / dx * dt

    ! write current state to screen
    print *, n, h

  end do time_loop

end program tsunami
