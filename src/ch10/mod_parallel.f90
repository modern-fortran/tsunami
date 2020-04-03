module mod_parallel

  ! A module to provide parallel facilities
  ! to the shallow water solver.

  use iso_fortran_env, only: int32, real32

  implicit none

  private
  public :: num_tiles, tile_indices, &
            tile_neighbors_1d, tile_neighbors_2d

  interface tile_indices
    module procedure :: tile_indices_1d, tile_indices_2d
  end interface tile_indices

contains

  pure function denominators(n)
    ! Returns all common denominators of n.
    integer(int32), intent(in) :: n
    integer(int32), allocatable :: denominators(:)
    integer(int32) :: i
    denominators = [integer(int32) ::]
    do i = 1, n
      if (mod(n, i) == 0) denominators = [denominators, i]
    end do
  end function denominators

  pure function num_tiles(n)
    ! Returns the optimal number of tiles in 2 dimensions
    ! given total number of tiles n.
    !
    ! Examples:
    !   * num_tiles(1) = [1, 1]
    !   * num_tiles(2) = [2, 1]
    !   * num_tiles(3) = [3, 1]
    !   * num_tiles(4) = [2, 2]
    !   * num_tiles(5) = [5, 1]
    !   * num_tiles(6) = [3, 2]
    !
    integer(int32), intent(in) :: n
    integer(int32) :: num_tiles(2)
    integer(int32), allocatable :: denoms(:)
    integer(int32), allocatable :: dim1(:), dim2(:)
    integer(int32) :: i, j, n1, n2

    ! find all common denominators of the total number of images
    denoms = denominators(n)

    ! find all combinations of common denominators
    ! whose product equals the total number of images
    dim1 = [integer(int32) ::]
    dim2 = [integer(int32) ::]
    do j = 1, size(denoms)
      do i = 1, size(denoms)
        if (denoms(i) * denoms(j) == n) then
          dim1 = [dim1, denoms(i)]
          dim2 = [dim2, denoms(j)]
        end if
      end do
    end do

    ! pick the set of common denominators with the minimal norm
    ! between two elements -- rectangle closest to a square
    num_tiles = [dim1(1), dim2(1)]
    do i = 2, size(dim1)
      n1 = norm2([dim1(i), dim2(i)] - sqrt(real(n)))
      n2 = norm2(num_tiles - sqrt(real(n)))
      if (n1 < n2) num_tiles = [dim1(i), dim2(i)]
    end do

  end function num_tiles

  pure function tile_indices_1d(dims, i, n) result(indices)
    ! Given input global array size, return start and end index
    ! of a parallel 1-d tile that correspond to this image.
    integer(int32), intent(in) :: dims, i, n
    integer(int32) :: indices(2)
    integer(int32) :: offset, tile_size

    tile_size = dims / n

    ! start and end indices assuming equal tile sizes
    indices(1) = (i - 1) * tile_size + 1
    indices(2) = indices(1) + tile_size - 1

    ! if we have any remainder, distribute it to the tiles at the end
    offset = n - mod(dims, n)
    if (i > offset) then
      indices(1) = indices(1) + i - offset - 1
      indices(2) = indices(2) + i - offset
    end if

  end function tile_indices_1d


  pure function tile_indices_2d(dims) result(indices)
    ! Given an input x- and y- dimensions of the total computational domain [im, jm].
    ! returns an array of start and end indices in x- and y-, [is, ie, js, je].
    integer(int32), intent(in) :: dims(2)
    integer(int32) :: indices(4)
    integer(int32) :: tiles(2), tiles_ij(2)
    tiles = num_tiles(num_images())
    tiles_ij = tile_n2ij(this_image())
    indices(1:2) = tile_indices_1d(dims(1), tiles_ij(1), tiles(1))
    indices(3:4) = tile_indices_1d(dims(2), tiles_ij(2), tiles(2))
  end function tile_indices_2d


  pure function tile_neighbors_1d() result(neighbors)
    ! Returns the image indices corresponding
    ! to left and right neighbor tiles.
    integer(int32) :: neighbors(2)
    integer(int32) :: left, right
    if (num_images() > 1) then
      left = this_image() - 1
      right = this_image() + 1
      if (this_image() == 1) then
        left = num_images()
      else if (this_image() == num_images()) then
        right = 1
      end if
    else
      left = 1
      right = 1
    end if
    neighbors = [left, right]
  end function tile_neighbors_1d


  pure function tile_n2ij(n) result(ij)
    ! Given tile index in a 1-d layout, returns the
    ! corresponding tile indices in a 2-d layout.
    !
    !    +---+---+---+
    !  2 | 4 | 5 | 6 |
    !    +---+---+---+
    !  1 | 1 | 2 | 3 |
    !  j +---+---+---+
    !    i 1   2   3
    !
    ! Examples:
    !   * tile_n2ij(2) = [2, 1]
    !   * tile_n2ij(4) = [1, 2]
    !   * tile_n2ij(6) = [3, 2]
    !
    integer(int32), intent(in) :: n
    integer(int32) :: ij(2), i, j, tiles(2)
    if (n == 0) then
      ij = 0
    else
      tiles = num_tiles(num_images())
      j = (n - 1) / tiles(1) + 1
      i = n - (j - 1) * tiles(1)
      ij = [i, j]
    end if
  end function tile_n2ij


  pure function tile_ij2n(ij) result(n)
    ! Given tile indices in a 2-d layout, returns the
    ! corresponding tile index in a 1-d layout:
    !
    !    +---+---+---+
    !  2 | 4 | 5 | 6 |
    !    +---+---+---+
    !  1 | 1 | 2 | 3 |
    !  j +---+---+---+
    !    i 1   2   3
    !
    ! Examples:
    !   * tile_ij2n([2, 1]) = 2
    !   * tile_ij2n([1, 2]) = 4
    !   * tile_ij2n([3, 2]) = 6
    !
    integer(int32), intent(in) :: ij(2)
    integer(int32) :: n, tiles(2)
    if (any(ij == 0)) then
      n = 0
    else
      tiles = num_tiles(num_images())
      n = (ij(2) - 1) * tiles(1) + ij(1)
    end if
  end function tile_ij2n


  pure function tile_neighbors_2d(periodic) result(neighbors)
    ! Returns the neighbor image indices given.
    logical, intent(in) :: periodic
    integer(int32) :: neighbors(4)
    integer(int32) :: tiles(2), tiles_ij(2), itile, jtile
    integer(int32) :: left, right, down, up
    integer(int32) :: ij_left(2), ij_right(2), ij_down(2), ij_up(2)

    tiles = num_tiles(num_images())
    tiles_ij = tile_n2ij(this_image())
    itile = tiles_ij(1)
    jtile = tiles_ij(2)

    ! i, j tile indices for each of the neighbors
    ij_left = [itile - 1, jtile]
    ij_right = [itile + 1, jtile]
    ij_down = [itile, jtile - 1]
    ij_up = [itile, jtile + 1]

    if (periodic) then
      ! set neighbor to wrap around the edge
      if (ij_left(1) < 1) ij_left(1) = tiles(1)
      if (ij_right(1) > tiles(1)) ij_right(1) = 1
      if (ij_down(2) < 1) ij_down(2) = tiles(2)
      if (ij_up(2) > tiles(2)) ij_up(2) = 1
    else
      ! set neighbor to 0 -- no neighbor
      if (ij_left(1) < 1) ij_left = 0
      if (ij_right(1) > tiles(1)) ij_right = 0
      if (ij_down(2) < 1) ij_down = 0
      if (ij_up(2) > tiles(2)) ij_up = 0
    end if

    left = tile_ij2n(ij_left)
    right = tile_ij2n(ij_right)
    down = tile_ij2n(ij_down)
    up = tile_ij2n(ij_up)

    neighbors = [left, right, down, up]

  end function tile_neighbors_2d

end module mod_parallel
