program tsunami

  ! Tsunami simulator
  !
  ! Solves the non-linear 1-d shallow water equation:
  !
  !     du/dt + u du/dx + g dh/dx = 0
  !
  !     dh/dt + d(hu)/dx = 0
  !
  ! This version is parallelized.

  use mod_diff, only: diff => diff_centered
  use mod_kinds, only: ik, rk
  use mod_parallel, only: tile_indices, tile_neighbors

  implicit none

  integer(ik) :: i, n

  integer(ik), parameter :: grid_size = 100 ! grid size in x
  integer(ik), parameter :: num_time_steps = 5000 ! number of time steps

  real(rk), parameter :: dt = 0.02 ! time step [s]
  real(rk), parameter :: dx = 1 ! grid spacing [m]
  real(rk), parameter :: g = 9.8 ! gravitational acceleration [m/s]

  real(rk), allocatable :: h(:)[:], u(:)[:]
  real(rk), allocatable :: gather(:)[:]
  real(rk), allocatable :: hmean(:)

  integer(ik), parameter :: ipos = 25
  real(rk), parameter :: decay = 0.02

  integer(ik) :: indices(2), neighbors(2)
  integer(ik) :: left, right
  integer(ik) :: is, ie ! global start and end indices
  integer(ik) :: ils, ile ! local start and end computational indices
  integer(ik) :: ims, ime ! local start and end memory indices
  integer(ik) :: tile_size

  if (mod(im, num_images()) > 0) then
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
  if (this_image() == 1) print *, 0, gather

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
    if (this_image() == 1) print *, n, gather

  end do time_loop

end program tsunami
