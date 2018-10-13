program tsunami

  ! Tsunami simulator
  !
  ! Solves the non-linear 2-d shallow water equation system:
  !
  !     du/dt + u du/dx + v du/dy + g dh/dx = 0
  !     dv/dt + u dv/dx + v dv/dy + g dh/dy = 0
  !     dh/dt + d(hu)/dx + d(hv)/dy = 0
  !
  ! This version is parallelized.

  use iso_fortran_env, only: output_unit

  use mod_diagnostics, only: ke, mean
  use mod_diff, only: diffx => diffc_2d_x, diffy => diffc_2d_y
  use mod_io, only: write_field
  use mod_kinds, only: ik, rk
  use mod_parallel, only: tile_indices, update_halo

  implicit none

  integer(ik) :: i, j, n

  integer(ik), parameter :: im = 101 ! grid size in x
  integer(ik), parameter :: jm = 101 ! grid size in y
  integer(ik), parameter :: nm = 1000 ! number of time steps

  real(rk), parameter :: dt = 0.02 ! time step [s]
  real(rk), parameter :: dx = 1 ! grid spacing [m]
  real(rk), parameter :: dy = 1 ! grid spacing [m]

  real(rk), parameter :: g = 9.8 ! gravitational acceleration [m/s]

  real(rk), allocatable :: h(:,:), u(:,:), v(:,:)
  real(rk), allocatable :: gather(:,:)[:]
  real(rk), allocatable :: hmean(:,:)

  integer(ik), parameter :: ipos = 51, jpos = 51
  real(rk), parameter :: decay = 0.02

  integer(ik) :: is, ie, js, je ! global start and end indices
  integer(ik) :: indices(4)

  if (this_image() == 1) print *, 'Tsunami started'

  indices = tile_indices([im, jm])
  is = indices(1)
  ie = indices(2)
  js = indices(3)
  je = indices(4)

  allocate(h(is-1:ie+1, js-1:je+1))
  allocate(u(is-1:ie+1, js-1:je+1))
  allocate(v(is-1:ie+1, js-1:je+1))
  allocate(hmean(is-1:ie+1, js-1:je+1))

  allocate(gather(im, jm)[*])

  ! initialize a gaussian blob centered at i = 25
  do concurrent(i = is-1:ie+1, j = js-1:je+1)
    h(i, j) = exp(-decay * ((i - ipos)**2 + (j - jpos)**2))
  end do

  ! set initial velocity and mean water depth
  u = 0
  v = 0
  hmean = 10

  time_loop: do n = 1, nm

    call update_halo(h, indices)

    ! compute u at next time step
    u = u - (u * diffx(u) / dx + v * diffy(u) / dy &
      + g * diffx(h) / dx) * dt

    ! compute v at next time step
    v = v - (u * diffx(v) / dx + v * diffy(v) / dy &
      + g * diffy(h) / dy) * dt

    sync all
    call update_halo(u, indices)
    call update_halo(v, indices)

    ! compute h at next time step
    h = h - diffx(u * (hmean + h)) / dx * dt&
          - diffy(v * (hmean + h)) / dy * dt
    sync all

    ! gather to image 1 and write water height to file
    gather(is:ie, js:je)[1] = h(is:ie, js:je)
    sync all
    if (this_image() == 1) then
      print *, n, gather(3 * im / 4, jm / 2)
      call write_field(gather, 'h', n)
    end if

  end do time_loop

end program tsunami
