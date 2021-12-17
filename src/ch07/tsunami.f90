program tsunami

  ! Tsunami simulator
  !
  ! Solves the non-linear 1-d shallow water equation:
  !
  !     du/dt + u du/dx + g dh/dx = 0
  !     dh/dt + d(hu)/dx = 0
  !
  ! This version is parallelized.

  use iso_fortran_env, only: int32, real32
  use mod_diff, only: diff => diff_centered
  use mod_parallel, only: tile_indices, tile_neighbors

  implicit none

  integer(int32) :: i, n

  integer(int32), parameter :: grid_size = 100 ! grid size in x
  integer(int32), parameter :: num_time_steps = 5000 ! number of time steps

  real(real32), parameter :: dt = 0.02 ! time step [s]
  real(real32), parameter :: dx = 1 ! grid spacing [m]
  real(real32), parameter :: g = 9.8 ! gravitational acceleration [m/s]

  real(real32), allocatable :: h(:)[:], u(:)[:]
  real(real32), allocatable :: gather(:)[:]
  real(real32), allocatable :: hmean(:)

  integer(int32), parameter :: ipos = 25
  real(real32), parameter :: decay = 0.02

  integer(int32) :: indices(2), neighbors(2)
  integer(int32) :: left, right
  integer(int32) :: is, ie ! global start and end indices
  integer(int32) :: ils, ile ! local start and end computational indices
  integer(int32) :: ims, ime ! local start and end memory indices
  integer(int32) :: tile_size

  character(*), parameter :: fmt = '(i0,*(1x,es15.8e2))'

  if (mod(grid_size, num_images()) > 0) then
    error stop 'Error: grid_size must be divisible by number of images'
  end if

  neighbors = tile_neighbors()
  left = neighbors(1)
  right = neighbors(2)

  indices = tile_indices(grid_size)
  is = indices(1)
  ie = indices(2)

  tile_size = grid_size / num_images()
  ils = 1
  ile = tile_size
  ims = ils - 1
  ime = ile + 1

  allocate(h(ims:ime)[*])
  allocate(u(ims:ime)[*])
  allocate(hmean(ims:ime))

  allocate(gather(grid_size)[*])

  ! initialize a gaussian blob centered at i = 25
  do i = is - 1, ie + 1
    h(i-is+1) = exp(-decay * (i - ipos)**2)
  end do

  ! set initial velocity and mean water depth
  u = 0
  hmean = 10

  ! gather to image 1 and write current state to screen
  gather(is:ie)[1] = h(ils:ile)
  sync all
  if (this_image() == 1) print fmt, 0, gather

  time_loop: do n = 1, num_time_steps

    ! update halo for h
    h(ime)[left] = h(ils)
    h(ims)[right] = h(ile)
    sync all

    ! compute u at next time step
    u = u - (u * diff(u) + g * diff(h)) / dx * dt

    sync all

    ! update halo for u
    u(ime)[left] = u(ils)
    u(ims)[right] = u(ile)
    sync all

    ! compute h at next time step
    h = h - diff(u * (hmean + h)) / dx * dt

    ! gather to image 1 and write current state to screen
    gather(is:ie)[1] = h(ils:ile)
    sync all
    if (this_image() == 1) print fmt, n, gather

  end do time_loop

end program tsunami
