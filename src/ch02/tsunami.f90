program tsunami

  ! This version solves the linearized 1-d advection equation:
  !
  !     du/dt + c du/dx = 0

  implicit none

  integer :: i, n

  integer, parameter :: im = 100 ! grid size in x
  integer, parameter :: nm = 100 ! number of time steps

  real, parameter :: dt = 1 ! time step [s]
  real, parameter :: dx = 1 ! grid spacing [m]
  real, parameter :: c = 1 ! phase speed [m/s]

  real :: du(im), u(im)

  integer, parameter :: ipos = 25
  real, parameter :: decay = 0.02

  ! initialize a gaussian blob centered at i = 25
  do i = 1, im
    u(i) = exp(-decay * (i - ipos)**2)
  end do

  ! write initial state to screen
  print *, 0, u

  time_loop: do n = 1, nm

    ! apply the periodic boundary condition
    du(1) = u(1) - u(im)

    ! calculate the difference of u in space
    do concurrent (i = 2:im)
      du(i) = u(i) - u(i-1)
    end do

    ! compute u at next time step
    do concurrent (i = 1:im)
      u(i) = u(i) - c * du(i) / dx * dt
    end do

    ! write current state to screen
    print *, n, u

  end do time_loop

end program tsunami
