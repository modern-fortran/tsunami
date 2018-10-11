program test_halo_update

  use mod_kinds
  use mod_parallel

  implicit none

  integer(ik) :: i, j, n

  integer(ik), parameter :: im = 10, jm = 10 ! grid size in x and y

  integer(ik) :: indices(4), tiles(2), is, ie, js, je
  integer(ik) :: imax, jmax

  real(rk), allocatable :: u(:,:)

  indices = tile_indices([im, jm])
  is = indices(1)
  ie = indices(2)
  js = indices(3)
  je = indices(4)

  allocate(u(is-1:ie+1, js-1:je+1))
  u = this_image()

  if (this_image() == 1) then
    print *, 'before halo, west,', u(is-1, js:je)
    print *, 'before halo, east,', u(ie+1, js:je)
    print *, 'before halo, south,', u(is:ie, js-1)
    print *, 'before halo, north,', u(is:ie, je+1)
  end if

  call update_halo(u, indices)

  if (this_image() == 1) then
    print *, 'after halo, west,', u(is-1, js:je)
    print *, 'after halo, east,', u(ie+1, js:je)
    print *, 'after halo, south,', u(is:ie, js-1)
    print *, 'after halo, north,', u(is:ie, je+1)
  end if

end program test_halo_update
