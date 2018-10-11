program test_halo_update

  use mod_kinds
  use mod_parallel

  implicit none

  integer(ik) :: i, j, n

  integer(ik), parameter :: im = 3 ! grid size in x
  integer(ik), parameter :: jm = 3 ! grid size in y

  integer(ik) :: neighbors(4)

  integer(ik) :: ix(2), iy(2), tiles(2), ij(2)
  integer(ik) :: imax, jmax

  print *, 'this_image = ', this_image(), 'of', num_images()

end program test_halo_update
