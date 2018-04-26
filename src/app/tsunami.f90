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

  use iso_fortran_env, only: output_unit
  use mod_diff, only: diff => diffc
  use mod_kinds, only: ik, rk
  use mod_parallel, only: tile_indices, tile_neighbors

  implicit none

  integer(kind=ik) :: i, n

  integer(kind=ik), parameter :: im = 100 ! grid size in x
  integer(kind=ik), parameter :: nm = 5000 ! number of time steps

  real(kind=rk), parameter :: dt = 0.02 ! time step [s]
  real(kind=rk), parameter :: dx = 1 ! grid spacing [m]

  real(kind=rk), parameter :: g = 9.8 ! gravitational acceleration [m/s]

  real(kind=rk), allocatable :: h(:)[:], u(:)[:]
  real(kind=rk), allocatable :: gather(:)[:]
  real(kind=rk), allocatable :: hmean(:)

  integer(kind=ik), parameter :: ipos = 25
  real(kind=rk), parameter :: decay = 0.02

  integer(kind=ik), dimension(2) :: indices, neighbors
  integer(kind=ik) :: left, right
  integer(kind=ik) :: is, ie ! global start and end indices
  integer(kind=ik) :: ils, ile ! local start and end computational indices
  integer(kind=ik) :: ims, ime ! local start and end memory indices
  integer(kind=ik) :: tile_size

  if (mod(im, num_images()) > 0) then
    error stop 'Error: im must be divisible by number of images'
  end if

  neighbors = tile_neighbors()
  left = neighbors(1)
  right = neighbors(2)

  indices = tile_indices(im)
  is = indices(1)
  ie = indices(2)

  tile_size = im / num_images()
  ils = 1
  ile = tile_size
  ims = ils - 1
  ime = ile + 1

  allocate(h(ims:ime)[*])
  allocate(u(ims:ime)[*])
  allocate(hmean(ims:ime))

  allocate(gather(im)[*])

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
  if (this_image() == 1) write(unit=output_unit, fmt=*) 0, gather

  time_loop: do n = 1, nm

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
    if (this_image() == 1) write(unit=output_unit, fmt=*) n, gather

  end do time_loop

end program tsunami
