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

  ! get number of tiles in x and y
  tiles = num_tiles(num_images())
  imax = tiles(1)
  jmax = tiles(2)

  ij = tile_index_2d_from_1d(this_image())
  i = ij(1)
  j = ij(2)

  neighbors = tile_neighbors_2d(periodic=.true.) 

  sync all
  do n = 1, num_images()
    if (this_image() == n) print *, '2d layout:', this_image(), i, j
    sync all
  end do

  sync all
  do n = 1, num_images()
    if (this_image() == n) print *, '1d layout:', this_image(), tile_index_1d_from_2d(i, j)
    sync all
  end do

  if (this_image() == 1) print *, 'neighbors'
  sync all
  do n = 1, num_images()
    if (this_image() == n) print *, this_image(), neighbors
    sync all
  end do

end program test_tiling
