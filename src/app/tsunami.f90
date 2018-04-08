program tsunami

  ! Tsunami simulator
  !
  ! Solves the non-linear 1-d shallow water equation:
  !
  !     du/dt + u du/dx + g dh/dx= 0
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

  real(kind=rk), dimension(:), codimension[:], allocatable :: h, u
  real(kind=rk), dimension(:), codimension[:], allocatable :: h_gather
  real(kind=rk), dimension(:), allocatable :: hmean

  integer(kind=ik), parameter :: ipos = 25
  real(kind=rk), parameter :: decay = 0.02

  integer(kind=ik), dimension(2) :: indices, neighbors
  integer(kind=ik) :: left, right
  integer(kind=ik) :: its, ite, is, ie ! start and end tile indices
  integer(kind=ik) :: ils, ile
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
  its = 0
  ite = tile_size + 1

  ils = 1
  ile = tile_size

  allocate(h(its:ite)[*])
  allocate(u(its:ite)[*])
  allocate(hmean(its:ite))

  allocate(h_gather(im)[*])
  h_gather = 0

  ! initialize a gaussian blob centered at i = 25
  do i = is - 1, ie + 1
    h(i-is+1) = exp(-decay * (i - ipos)**2)
  end do

  ! set initial velocity to zero
  u = 0

  ! set mean water depth to 10 m
  hmean = 10

  ! gather to image 1 and write current state to screen
  h_gather(is:ie)[1] = h(ils:ile)
  sync all
  if (this_image() == 1)  write(unit=output_unit, fmt=*) 0, h_gather

  time_loop: do n = 1, nm

    ! update halo for h
    h(ite)[left] = h(ils)
    h(its)[right] = h(ile)
    sync all

    ! compute u at next time step
    u = u - (u * diff(u) + g * diff(h)) / dx * dt

    sync all

    ! update halo for u
    u(ite)[left] = u(ils)
    u(its)[right] = u(ile)
    sync all

    ! compute h at next time step
    h = h - diff(u * (hmean + h)) / dx * dt

    ! gather to image 1 and write current state to screen
    h_gather(is:ie)[1] = h(ils:ile)
    sync all
    if (this_image() == 1)  write(unit=output_unit, fmt=*) n, h_gather

  end do time_loop

end program tsunami
