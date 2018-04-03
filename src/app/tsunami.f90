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
  real(kind=rk), dimension(:), allocatable :: hmean

  integer(kind=ik), parameter :: ipos = 25
  real(kind=rk), parameter :: decay = 0.02

  integer(kind=ik), dimension(2) :: indices, neighbors
  integer(kind=ik) :: image_left, image_right
  integer(kind=ik) :: its, ite, is, ie ! start and end tile indices

  neighbors = tile_neighbors()
  image_left = neighbors(1)
  image_right = neighbors(2)

  indices = tile_indices(im)
  is = indices(1)
  ie = indices(2)
 
  its = is - 1
  ite = ie + 1

  allocate(h(its:ite)[*])
  allocate(u(its:ite)[*])
  allocate(hmean(its:ite))

  ! initialize a gaussian blob centered at i = 25
  do i = its, ite
    h(i) = exp(-decay * (i - ipos)**2)
  end do

  ! set initial velocity to zero
  u = 0

  ! set mean water depth to 10 m
  hmean = 10

  ! write current state to screen
  do i = 1, num_images()
    sync all
    if (this_image() == i) then
      write(unit=output_unit, fmt=*) 0, this_image(), h(is:ie)
    end if
  end do

  time_loop: do n = 1, nm

    ! update halo for h
    h(its) = h(ie)[image_left]
    h(ite) = h(is)[image_right]
    sync all

    ! compute u at next time step
    u = u - (u * diff(u) + g * diff(h)) / dx * dt

    sync all

    ! update halo for u
    u(its) = u(ie)[image_left]
    u(ite) = u(is)[image_right]
    sync all

    ! compute h at next time step
    h = h - diff(u * (hmean + h)) / dx * dt

    ! write current state to screen
    do i = 1, num_images()
      sync all
      if (this_image() == i) then
        write(unit=output_unit, fmt=*) n, this_image(), h(is:ie)
      end if
    end do

  end do time_loop

end program tsunami
