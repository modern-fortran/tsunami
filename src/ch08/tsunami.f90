program tsunami

  ! Tsunami simulator.
  !
  ! Solves the non-linear 2-d shallow water equation system:
  !
  !     du/dt + u du/dx + v du/dy + g dh/dx = 0
  !     dv/dt + u dv/dx + v dv/dy + g dh/dy = 0
  !     dh/dt + d(hu)/dx + d(hv)/dy = 0
  !
  ! This version is parallelized, and uses derived types.

  use iso_fortran_env, only: int32, real32
  use mod_field, only: Field, diffx, diffy

  implicit none

  integer(int32) :: n

  integer(int32), parameter :: im = 101 ! grid size in x
  integer(int32), parameter :: jm = 101 ! grid size in y
  integer(int32), parameter :: num_time_steps = 1000 ! number of time steps

  real(real32), parameter :: dt = 0.02 ! time step [s]
  real(real32), parameter :: dx = 1 ! grid spacing [m]
  real(real32), parameter :: dy = 1 ! grid spacing [m]
  real(real32), parameter :: g = 9.8 ! gravitational acceleration [m/s]

  integer(int32), parameter :: ic = 51, jc = 51
  real(real32), parameter :: decay = 0.02

  type(Field) :: h, u, v, hm

  if (this_image() == 1) print *, 'Tsunami started'

  u = Field('u', [im, jm])
  v = Field('v', [im, jm])
  h = Field('h', [im, jm])
  hm = Field('hm', [im, jm])

  ! initialize a gaussian blob in the center
  call h % set_gaussian(decay, ic, jc)
  call h % sync_edges()

  ! set mean water depth
  hm = 10.

  call h % write(0)

  time_loop: do n = 1, num_time_steps

    if (this_image() == 1) then
      print *, 'Computing time step', n, '/', num_time_steps
    end if

    ! compute u at next time step
    u = u - (u * diffx(u) / dx + v * diffy(u) / dy &
      + g * diffx(h) / dx) * dt
    call u % sync_edges()

    ! compute v at next time step
    v = v - (u * diffx(v) / dx + v * diffy(v) / dy &
      + g * diffy(h) / dy) * dt
    call v % sync_edges()

    ! compute h at next time step
    h = h - (diffx(u * (hm + h)) / dx + diffy(v * (hm + h)) / dy) * dt
    call h % sync_edges()

    call h % write(n)

  end do time_loop

end program tsunami
