module mod_field

  ! Provides the Field class and its methods.

  use mod_kinds, only: ik, rk
  use mod_parallel, only: tile_indices, tile_neighbors_2d
  implicit none

  private
  public :: Field

  type :: Field
    character(:), allocatable :: name
    integer(ik) :: lb(2), ub(2)
    integer(ik) :: neighbors(4)
    integer(ik) :: edge_size
    real(rk), allocatable :: data(:,:)
  end type Field

  interface Field
    module procedure :: field_constructor
  end interface Field

contains

  type(Field) function field_constructor(name, dims) result(self)
    character(*), intent(in) :: name ! field name
    integer(ik), intent(in) :: dims(2) ! domain size in x and y
    integer(ik) :: edge_size, indices(4)
    self % name = name
    indices = tile_indices(dims)
    self % lb = indices([1, 3])
    self % ub = indices([2, 4])
    allocate(self % data(self % lb(1)-1:self % ub(1)+1,&
                         self % lb(2)-1:self % ub(2)+1))
    self % data = 0
    self % neighbors = tile_neighbors_2d(periodic=.true.)
    self % edge_size = max(self % ub(1)-self % lb(1)+1,&
                           self % ub(2)-self % lb(2)+1)
    call co_max(self % edge_size)
  end function field_constructor

end module mod_field
