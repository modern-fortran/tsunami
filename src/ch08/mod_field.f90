module mod_field

  ! Provides the Field class and its methods.

  use iso_fortran_env, only: int32, real32
  use mod_diff, only: diffx_real => diffx, diffy_real => diffy
  use mod_io, only: write_field
  use mod_parallel, only: tile_indices, tile_neighbors_2d

  implicit none

  private
  public :: Field, diffx, diffy

  type :: Field

    character(:), allocatable :: name
    integer(int32) :: lb(2), ub(2)
    integer(int32) :: dims(2)
    integer(int32) :: neighbors(4)
    integer(int32) :: edge_size
    real(real32), allocatable :: data(:,:)

  contains

    procedure, private, pass(self) :: assign_array, assign_const_int, assign_const_real
    procedure, private, pass(self) :: array_mult_field, field_mult_array, field_mult_real, field_mult_field
    procedure, private, pass(self) :: field_div_real
    procedure, private, pass(self) :: field_add_field, field_add_real, real_add_field
    procedure, private, pass(self) :: field_sub_field, field_sub_array
    procedure, public, pass(self) :: gather
    procedure, public, pass(self) :: set_gaussian
    procedure, public, pass(self) :: sync_edges
    procedure, public, pass(self) :: write

    generic :: assignment(=) => assign_array, assign_const_int, assign_const_real
    generic :: operator(+) => field_add_field, field_add_real, real_add_field
    generic :: operator(-) => field_sub_field, field_sub_array
    generic :: operator(*) => field_mult_array, array_mult_field, field_mult_real, field_mult_field
    generic :: operator(/) => field_div_real

  end type Field

  interface Field
    module procedure :: field_constructor
  end interface Field

contains

  type(Field) function field_constructor(name, dims) result(self)
    character(*), intent(in) :: name ! field name
    integer(int32), intent(in) :: dims(2) ! domain size in x and y
    integer(int32) :: edge_size, indices(4)
    self % name = name
    self % dims = dims
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

  pure subroutine assign_array(self, a)
    class(Field), intent(in out) :: self
    real(real32), intent(in) :: a(:,:)
    self % data = a
  end subroutine assign_array

  pure subroutine assign_const_int(self, a)
    class(Field), intent(in out) :: self
    integer(real32), intent(in) :: a
    self % data = a
  end subroutine assign_const_int

  pure subroutine assign_const_real(self, a)
    class(Field), intent(in out) :: self
    real(real32), intent(in) :: a
    self % data = a
  end subroutine assign_const_real

  pure function diffx(input_field)
    ! Returns the finite difference in x of input_field as a 2-d array.
    class(Field), intent(in) :: input_field
    real(real32), allocatable :: diffx(:,:)
    diffx = diffx_real(input_field % data)
  end function diffx

  pure function diffy(input_field)
    ! Returns the finite difference in y of input_field as a 2-d array.
    class(Field), intent(in) :: input_field
    real(real32), allocatable :: diffy(:,:)
    diffy = diffy_real(input_field % data)
  end function diffy

  function gather(self, image)
    ! Performs a gather of field data to image.
    class(Field), intent(in) :: self
    integer(int32), intent(in) :: image
    real(real32), allocatable :: gather_coarray(:,:)[:]
    real(real32) :: gather(self % dims(1), self % dims(2))
    allocate(gather_coarray(self % dims(1), self % dims(2))[*])
    associate(is => self % lb(1), ie => self % ub(1),&
              js => self % lb(2), je => self % ub(2))
      gather_coarray(is:ie, js:je)[image] = self % data(is:ie, js:je)
      sync all
      if (this_image() == image) gather = gather_coarray
    end associate
    deallocate(gather_coarray)
  end function gather

  pure subroutine set_gaussian(self, decay, ic, jc)
    class(Field), intent(in out) :: self
    real(real32), intent(in) :: decay ! the rate of decay of gaussian
    integer(int32), intent(in) :: ic, jc ! center indices of the gaussian blob
    integer(int32) :: i, j
    do concurrent(i = self % lb(1)-1:self % ub(1)+1,&
                  j = self % lb(2)-1:self % ub(2)+1)
      self % data(i, j) = exp(-decay * ((i - ic)**2 + (j - jc)**2))
    end do
  end subroutine set_gaussian

  pure type(Field) function field_add_field(self, f) result(res)
    class(Field), intent(in) :: self, f
    res = self
    res = self % data + f % data
  end function field_add_field

  pure type(Field) function field_add_real(self, x) result(res)
    class(Field), intent(in) :: self
    real(real32), intent(in) :: x(:,:)
    res = self
    res = self % data + x
  end function field_add_real

  pure type(Field) function real_add_field(x, self) result(res)
    class(Field), intent(in) :: self
    real(real32), intent(in) :: x(:,:)
    res = self
    res = self % data + x
  end function real_add_field

  pure type(Field) function field_div_real(self, x) result(res)
    class(Field), intent(in) :: self
    real(real32), intent(in) :: x
    res = self
    res = self % data / x
  end function field_div_real

  pure type(Field) function field_mult_array(self, x) result(res)
    class(Field), intent(in) :: self
    real(real32), intent(in) :: x(:,:)
    res = self
    res = self % data * x
  end function field_mult_array

  pure type(Field) function field_mult_real(self, x) result(res)
    class(Field), intent(in) :: self
    real(real32), intent(in) :: x
    res = self
    res = self % data * x
  end function field_mult_real

  pure type(Field) function field_mult_field(self, f) result(res)
    class(Field), intent(in) :: self, f
    res = self
    res = self % data * f % data
  end function field_mult_field

  pure type(Field) function array_mult_field(x, self) result(res)
    real(real32), intent(in) :: x(:,:)
    class(Field), intent(in) :: self
    res = self
    res = self % data * x
  end function array_mult_field

  pure type(Field) function field_sub_array(self, x) result(res)
    class(Field), intent(in) :: self
    real(real32), intent(in) :: x(:,:)
    res = self
    res = self % data - x
  end function field_sub_array

  pure type(Field) function field_sub_field(self, f) result(res)
    class(Field), intent(in) :: self, f
    res = self
    res = self % data - f % data
  end function field_sub_field

  subroutine sync_edges(self)
    class(Field), intent(in out) :: self
    real(real32), allocatable :: edge(:,:)[:]
    integer(int32) :: is, ie, js, je

    is = self % lb(1)
    ie = self % ub(1)
    js = self % lb(2)
    je = self % ub(2)

    if (.not. allocated(edge)) allocate(edge(self % edge_size, 4)[*])
    edge = 0

    sync all

    ! copy data into coarray buffer
    edge(1:je-js+1,1)[self % neighbors(1)] = self % data(is,js:je) ! send left
    edge(1:je-js+1,2)[self % neighbors(2)] = self % data(ie,js:je) ! send right
    edge(1:ie-is+1,3)[self % neighbors(3)] = self % data(is:ie,js) ! send down
    edge(1:ie-is+1,4)[self % neighbors(4)] = self % data(is:ie,je) ! send up

    sync all

    ! copy from halo buffer into array
    self % data(is-1,js:je) = edge(1:je-js+1,2) ! from left
    self % data(ie+1,js:je) = edge(1:je-js+1,1) ! from right
    self % data(is:ie,js-1) = edge(1:ie-is+1,4) ! from down
    self % data(is:ie,je+1) = edge(1:ie-is+1,3) ! from up

    deallocate(edge)

  end subroutine sync_edges

  subroutine write(self, n)
    class(Field), intent(in) :: self
    integer(int32), intent(in) :: n
    real(real32), allocatable :: gather(:,:)
    gather = self % gather(1)
    if (this_image() == 1) then
      call write_field(gather, self % name, n)
    end if
  end subroutine write

end module mod_field
