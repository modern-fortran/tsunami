program test_tiling

  use mod_kinds
  use mod_parallel

  implicit none

  integer(ik) :: i, j, n

  integer(ik), parameter :: im = 100 ! grid size in x
  integer(ik), parameter :: jm = 100 ! grid size in y

  integer(ik) :: neighbors(4)

  integer(ik) :: ix(2), iy(2), tiles(2), ij(2)
  integer(ik) :: imax, jmax

  print *, 'this_image = ', this_image(), 'of', num_images()

  sync all
  do n = 1, num_images()
    if (this_image() == n) print *, 'tile_ij2n(tile_n2ij(n)):', this_image(), tile_ij2n(tile_n2ij(this_image()))
    sync all
  end do

  if (this_image() == 1) print *, 'neighbors, periodic'
  sync all
  do n = 1, num_images()
    if (this_image() == n) print *, this_image(), tile_neighbors_2d(periodic=.true.)
    sync all
  end do

  if (this_image() == 1) print *, 'neighbors, closed'
  sync all
  do n = 1, num_images()
    if (this_image() == n) print *, this_image(), tile_neighbors_2d(periodic=.false.)
    sync all
  end do

end program test_tiling
