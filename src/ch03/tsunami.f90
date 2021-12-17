program tsunami

  ! This version solves the linearized 1-d advection equation:
  !
  !     du/dt + c du/dx = 0
  !
  ! The initial conditions and the finite difference calculation 
  ! are abstracted as external procedures.

  implicit none

  integer :: n

  integer, parameter :: grid_size = 100 ! grid size in x
  integer, parameter :: num_time_steps = 100 ! number of time steps

  real, parameter :: dt = 1 ! time step [s]
  real, parameter :: dx = 1 ! grid spacing [m]
  real, parameter :: c = 1 ! phase speed [m/s]

  real :: h(grid_size)

  integer, parameter :: icenter = 25
  real, parameter :: decay = 0.02

  character(*), parameter :: fmt = '(i0,*(1x,es15.8e2))'

  ! check input parameter values
  if (grid_size <= 0) stop 'grid_size must be > 0'
  if (dt <= 0) stop 'time step dt must be > 0'
  if (dx <= 0) stop 'grid spacing dx must be > 0'
  if (c <= 0) stop 'background flow speed c must be > 0'

  ! initialize to a Gaussian shape
  call set_gaussian(h, icenter, decay)

  ! write initial state to screen
  print fmt, 0, h

  time_loop: do n = 1, num_time_steps

    ! compute u at next time step
    h = h - c * diff(h) / dx * dt

    ! write current state to screen
    print fmt, n, h

  end do time_loop

contains

  pure function diff(x) result(dx)
    ! Returns a 1st-order upstream finite difference of a 1-d array.
    real, intent(in) :: x(:)
    real :: dx(size(x))
    integer :: im
    im = size(x)
    dx(1) = x(1) - x(im)
    dx(2:im) = x(2:im) - x(1:im-1)
  end function diff

  pure subroutine set_gaussian(x, icenter, decay)
    ! Sets the values of x to a Gaussian shape centered on icenter
    ! that decays with the given input decay.
    real, intent(in out) :: x(:)
    integer, intent(in) :: icenter
    real, intent(in) :: decay
    integer :: i
    do concurrent(i = 1:size(x))
      x(i) = exp(-decay * (i - icenter)**2)
    end do
  end subroutine set_gaussian

end program tsunami
