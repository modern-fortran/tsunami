program tsunami

  ! Tsunami simulator
  !
  ! This version solves the non-linear 1-d shallow water equation:
  !
  !     du/dt + u du/dx + g dh/dx= 0
  !
  !     dh/dt + d(hu)/dx = 0

  use iso_fortran_env, only: output_unit
  use mod_diff, only: diff => diffc
  use mod_kinds, only: ik, rk
  use mod_parallel, only: tile_indices

  implicit none

  integer(kind=ik) :: i, n

  integer(kind=ik), parameter :: im = 100 ! grid size in x
  integer(kind=ik), parameter :: nm = 5000 ! number of time steps

  real(kind=rk), parameter :: dt = 0.02 ! time step [s]
  real(kind=rk), parameter :: dx = 1 ! grid spacing [m]

  real(kind=rk), parameter :: g = 9.8 ! gravitational acceleration [m/s]

  !real(kind=rk), dimension(im) :: h, hmean, u
  real(kind=rk), dimension(:), codimension[:], allocatable :: h, hmean, u

  integer(kind=ik), parameter :: ipos = 25
  real(kind=rk), parameter :: decay = 0.02

  integer(kind=ik) :: image_left, image_right
  integer(kind=ik) :: its, ite, is, ie ! start and end tile indices
  integer(kind=ik), dimension(2) :: indices

  logical :: parallel

  parallel = num_images() > 1

  if (parallel) then
    if (this_image() == 1) then
      image_left = num_images()
      image_right = 2
    else if (this_image() > 1 .and. this_image() < num_images()) then
      image_left = this_image() - 1
      image_right = this_image() + 1
    else
      image_left = num_images() - 1
      image_right = 1
    end if
  else
    image_left = 1
    image_right = 1
  end if

  indices = tile_indices(im)
  is = indices(1)
  ie = indices(2)
 
  its = is - 1
  ite = ie + 1

  allocate(h(its:ite)[*])
  allocate(u(its:ite)[*])
  allocate(hmean(its:ite)[*])

  ! TODO Solution is reproducible only if tiles are even

  !write(*,*)this_image(), image_left, image_right
  !write(*,*)this_image(), is, ie, its, ite
  !stop

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
