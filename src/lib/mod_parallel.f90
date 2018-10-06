module mod_parallel

  ! A module to provide parallel facilities
  ! to the shallow water solver.

  use mod_kinds, only: ik, rk

  implicit none

  private
  public :: num_tiles, tile_indices, tile_neighbors_1d, &
            tile_neighbors_2d, update_halo, allocate_coarray, &
            tile_index_2d_from_1d, tile_index_1d_from_2d

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
    neighbors = [left, right]
  end function tile_neighbors_1d


  pure function tile_index_2d_from_1d(n) result(ij)
    ! Given tile index in a 1-d layout, returns the 
    ! corresponding tile indices in a 2-d layout.
    integer(ik), intent(in) :: n
    integer(ik) :: ij(2), i, j, tiles(2)
    tiles = num_tiles(num_images())
    j = (n - 1) / tiles(1) + 1
    i = n - (j - 1) * tiles(1)
    ij = [i, j]
  end function tile_index_2d_from_1d
    
    
  pure function tile_index_1d_from_2d(i, j) result(n)
    ! Given tile indices in a 2-d layout, returns the 
    ! corresponding tile index in a 1-d layout.
    integer(ik), intent(in) :: i, j
    integer(ik) :: n, tiles(2)
    tiles = num_tiles(num_images())
    n = (j - 1) * tiles(1) + i
  end function tile_index_1d_from_2d


  pure function tile_neighbors_2d(periodic) result(neighbors)
    ! Returns the neighbor image indices given
    logical, intent(in) :: periodic
    integer(ik) :: neighbors(4)
    integer(ik) :: tiles(2), tiles_ij(2), itile, jtile
    integer(ik) :: left, right, down, up
    integer(ik) :: ij_left(2), ij_right(2), ij_down(2), ij_up(2)

    tiles = num_tiles(num_images())
    tiles_ij = tile_index_2d_from_1d(this_image())
    itile = tiles_ij(1)
    jtile = tiles_ij(2)

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
      if (left < 1) left = 0
      if (right > tiles(1)) right = 0
      if (down < 1) down = 0
      if (up > tiles(2)) up = 0
    end if

    left = tile_index_1d_from_2d(ij_left(1), ij_left(2))
    right = tile_index_1d_from_2d(ij_right(1), ij_right(2))
    down = tile_index_1d_from_2d(ij_down(1), ij_down(2))
    up = tile_index_1d_from_2d(ij_up(1), ij_up(2))

    neighbors = [left, right, down, up]

  end function tile_neighbors_2d

  subroutine update_halo(a)
    real(rk), allocatable, intent(in out) :: a(:,:)
    real(rk), allocatable :: halo(:,:)[:]
    integer(ik) :: ix(2), iy(2), tiles(2)
    integer(ik) :: itile, jtile, is, ie, js, je, im, jm
    integer(ik) :: neighbors(4)
    integer(ik) :: stat
    character(len=100) :: errmsg

    if (.not. allocated(a)) then
      stop 'Error in update_halo: input array not allocated.'
    end if

    im = size(a, dim=1)
    jm = size(a, dim=2)
    allocate(halo(100, 4)[*])
    !if (.not. allocated(halo)) allocate(halo(100, 4)[*], stat=stat, errmsg=errmsg)

    ! tile layout in 2-d
    tiles = num_tiles(num_images())
    jtile = (this_image() - 1) / tiles(1) + 1
    itile = this_image() - (jtile - 1) * tiles(1)

    neighbors = tile_neighbors_2d(periodic=.true.)
    print *, 'update_halo, this_image, neighbors:', this_image(), neighbors

    ix = tile_indices(im, itile, tiles(1)) ! start and end index in x
    iy = tile_indices(jm, jtile, tiles(2)) ! start and end index in y

    is = ix(1)
    ie = ix(2)
    js = iy(1)
    je = iy(2)

    ! send to neighbors
    halo(1:je-js+1,1)[neighbors(1)] = a(is,js:je) ! send left
    halo(1:je-js+1,2)[neighbors(2)] = a(ie,js:je) ! send right
    halo(1:ie-is+1,3)[neighbors(3)] = a(is:ie,js) ! send down
    halo(1:ie-is+1,4)[neighbors(4)] = a(is:ie,je) ! send up

    sync images([this_image(), neighbors])

    ! copy from halo buffer into array
    a(is-1,js:je) = halo(1:je-js+1,2) ! from left
    a(ie+1,js:je) = halo(1:je-js+1,1) ! from right
    a(is:ie,js-1) = halo(1:ie-is+1,4) ! from down
    a(is:ie,je+1) = halo(1:ie-is+1,3) ! from up

    deallocate(halo)

  end subroutine update_halo

  subroutine allocate_coarray()
    real(rk), allocatable :: halo(:,:)[:]
    allocate(halo(1000, 4)[*])
    print *, 'allocated halo array'
    deallocate(halo)
  end subroutine allocate_coarray

end module mod_parallel
