module mod_parallel

  ! A module to provide parallel facilities
  ! to the shallow water solver.

  use mod_kinds, only: ik, rk

  implicit none

  private
  public :: num_tiles, tile_indices, tile_neighbors_1d, tile_neighbors_2d

contains

  pure function denominators(n)
    ! Returns all common denominators of n.
    integer(ik), intent(in) :: n
    integer(ik), allocatable :: denominators(:)
    integer(ik) :: i
    denominators = [integer(ik) ::]
    do i = 1, n
      if (mod(n, i) == 0) denominators = [denominators, i]
    end do
  end function denominators

  pure function num_tiles(n)
    ! Returns the optimal number of tiles in 2 dimensions
    ! given total number of tiles n.
    integer(ik), intent(in) :: n
    integer(ik) :: num_tiles(2)
    integer(ik), allocatable :: denoms(:)
    integer(ik), allocatable :: dim1(:), dim2(:)
    integer(ik) :: i, j, n1, n2

    ! find all common denominators of the total number of images
    denoms = denominators(n)

    ! find all combinations of common denominators
    ! whose product equals the total number of images
    dim1 = [integer(ik) ::]
    dim2 = [integer(ik) ::]
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

  pure function tile_indices(dims, i, n)

    ! Given input global array size, return start and end index
    ! of a parallel 1-d tile that correspond to this image.

    integer(ik), intent(in) :: dims, i, n
    integer(ik) :: tile_indices(2)
    integer(ik) :: offset, tile_size

    tile_size = dims / n

    ! start and end indices assuming equal tile sizes
    tile_indices(1) = (i - 1) * tile_size + 1
    tile_indices(2) = tile_indices(1) + tile_size - 1

    ! if we have any remainder, distribute it to the tiles at the end
    offset = n - mod(dims, n)
    if (i > offset) then
      tile_indices(1) = tile_indices(1) + i - offset - 1
      tile_indices(2) = tile_indices(2) + i - offset
    end if

  end function tile_indices

  pure function tile_neighbors_1d() result(neighbors)
    ! Returns the image indices corresponding
    ! to left and right neighbor tiles.
    integer(ik) :: neighbors(2)
    integer(ik) :: left, right
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
    !tile_neighbors(1) = left
    !tile_neighbors(2) = right
    neighbors = [left, right]
  end function tile_neighbors_1d

  pure function tile_neighbors_2d(periodic) result(neighbors)
    ! Returns the neighbor image indices given
    logical, intent(in) :: periodic
    integer(ik) :: neighbors(4)
    integer(ik) :: tiles(2), itile, jtile
    integer(ik) :: left, right, down, up

    tiles = num_tiles(num_images())
    jtile = (this_image() - 1) / tiles(1) + 1
    itile = this_image() - (jtile - 1) * tiles(1)

    left = itile - 1
    right = itile + 1
    down = jtile - 1
    up = jtile + 1

    if (periodic) then
      ! set neighbor to wrap around the edge
      if (left < 1) left = tiles(1)
      if (right > tiles(1)) right = 1
      if (down < 1) down = tiles(2)
      if (up > tiles(2)) up = 1
    else
      ! set neighbor to 0 -- no neighbor
      if (left < 1) left = 0
      if (right > tiles(1)) right = 0
      if (down < 1) down = 0
      if (up > tiles(2)) up = 0
    end if

    neighbors = [left, right, down, up]

  end function tile_neighbors_2d

end module mod_parallel
