program tsunami

  ! This version solves the linearized 1-d advection equation:
  !
  !     du/dt + c du/dx = 0
  !
  ! The finite difference calculation is abstracted in the function diff().

  implicit none

  integer :: i, n

  integer, parameter :: im = 100 ! grid size in x
  integer, parameter :: nm = 100 ! number of time steps

  real, parameter :: dt = 1 ! time step [s]
  real, parameter :: dx = 1 ! grid spacing [m]
  real, parameter :: c = 1 ! phase speed [m/s]

  real :: u(im)

  integer, parameter :: ipos = 25
  real, parameter :: decay = 0.02

  ! initialize a gaussian blob centered at i = 25
  do concurrent(i = 1:im)
    u(i) = exp(-decay * (i - ipos)**2)
  end do

  ! write initial state to screen
  print *, 0, u

  time_loop: do n = 1, nm

    ! compute u at next time step
    u = u - c * diff(u) / dx * dt

    ! write current state to screen
    print *, n, u

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

end program tsunami
