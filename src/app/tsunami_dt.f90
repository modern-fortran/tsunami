program tsunami_dt

  ! Tsunami simulator.
  !
  ! Solves the non-linear 2-d shallow water equation system:
  !
  !     du/dt + u du/dx + v du/dy + g dh/dx = 0
  !     dv/dt + u dv/dx + v dv/dy + g dh/dy = 0
  !     dh/dt + d(hu)/dx + d(hv)/dy = 0
  !
  ! This version is parallelized, and uses derived types.

  use mod_diagnostics, only: mean
  use mod_field, only: Field, diffx, diffy
  use mod_io, only: write_field
  use mod_kinds, only: ik, rk
  use mod_parallel, only: sync_edges

  implicit none

  integer(ik) :: i, j, n

  integer(ik), parameter :: im = 101 ! grid size in x
  integer(ik), parameter :: jm = 101 ! grid size in y
  integer(ik), parameter :: nm = 1000 ! number of time steps

  real(rk), parameter :: dt = 0.02 ! time step [s]
  real(rk), parameter :: dx = 1 ! grid spacing [m]
  real(rk), parameter :: dy = 1 ! grid spacing [m]
  real(rk), parameter :: g = 9.8 ! gravitational acceleration [m/s]

  real(rk), allocatable :: gather(:,:)

  integer(ik), parameter :: ic = 51, jc = 51
  real(rk), parameter :: decay = 0.02

  type(Field) :: h, u, v, hm

  if (this_image() == 1) print *, 'Tsunami started'

  u = Field('x-component of velocity', [im, jm])
  v = Field('y-component of velocity', [im, jm])
  h = Field('Water height displacement', [im, jm])
  hm = Field('Mean water height', [im, jm])

  ! initialize a gaussian blob centered at i = 25
  call h % init_gaussian(decay, ic, jc)
  call h % sync_edges()

  ! set initial velocity and mean water depth
  u = 0.
  v = 0.
  hm = 10.

  gather = h % gather(1)
  n = 0
  if (this_image() == 1) then
    print *, n, mean(gather)
    call write_field(gather, 'h', n)
  end if

  time_loop: do n = 1, nm

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

    gather = h % gather(1)
    if (this_image() == 1) then
      print *, n, mean(gather)
      call write_field(gather, 'h', n)
    end if

  end do time_loop

end program tsunami_dt
