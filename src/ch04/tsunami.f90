program tsunami

  ! This version solves the linearized 1-d advection equation:
  !
  !     du/dt + c du/dx = 0
  !
  ! The initial conditions and the finite difference calculation 
  ! are abstracted as external procedures.

  use iso_fortran_env, only: int32, real32
  use mod_diff, only: diff
  use mod_initial, only: set_gaussian_blob

  implicit none

  integer(int32) :: i, n

  integer(int32), parameter :: im = 100 ! grid size in x
  integer(int32), parameter :: nm = 100 ! number of time steps

  real(real32), parameter :: dt = 1 ! time step [s]
  real(real32), parameter :: dx = 1 ! grid spacing [m]
  real(real32), parameter :: c = 1 ! phase speed [m/s]

  real(real32) :: u(im)

  integer(int32), parameter :: icenter = 25
  real(real32), parameter :: decay = 0.02

  ! initialize to a Gaussian blob
  call set_gaussian_blob(u, icenter, decay)

  ! write initial state to screen
  print *, 0, u

  time_loop: do n = 1, nm

    ! compute u at next time step
    u = u - c * diff(u) / dx * dt

    ! write current state to screen
    print *, n, u

  end do time_loop

end program tsunami
