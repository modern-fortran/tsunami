program tsunami

  ! This version solves the non-linear 1-d shallow water equations:
  !
  !     du/dt + u du/dx = - g dh/dx
  ! 
  !     du/dt + u du/dx = - g dh/dx
  !
  ! The initial conditions and the finite difference calculation 
  ! are abstracted as external procedures.

  use iso_fortran_env, only: int32, real32
  use mod_diff, only: diff => diff_centered
  use mod_initial, only: set_gaussian

  implicit none

  integer(int32) :: n

  integer(int32), parameter :: grid_size = 100 ! grid size in x
  integer(int32), parameter :: num_time_steps = 5000 ! number of time steps

  real(real32), parameter :: dt = 0.02 ! time step [s]
  real(real32), parameter :: dx = 1 ! grid spacing [m]
  real(real32), parameter :: g = 9.8 ! gravitational acceleration [m/s^2]
  real(real32), parameter :: hmean = 10 ! mean water depth [m]

  real(real32) :: h(grid_size), u(grid_size)

  integer(int32), parameter :: icenter = 25
  real(real32), parameter :: decay = 0.02

  character(*), parameter :: fmt = '(i0,*(1x,es15.8e2))'

  ! check input parameter values
  if (grid_size <= 0) stop 'grid_size must be > 0'
  if (dt <= 0) stop 'time step dt must be > 0'
  if (dx <= 0) stop 'grid spacing dx must be > 0'

  ! initialize water height to a Gaussian blob
  call set_gaussian(h, icenter, decay)

  ! initialize water velocity to zero
  u = 0

  ! write initial state to screen
  print fmt, 0, h

  time_loop: do n = 1, num_time_steps

    ! compute u at next time step
    u = u - (u * diff(u) + g * diff(h)) / dx * dt

    ! compute h at next time step
    h = h - diff(u * (hmean + h)) / dx * dt

    ! write current state to screen
    print fmt, n, h

  end do time_loop

end program tsunami
