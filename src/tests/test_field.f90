program test_field

  use mod_kinds, only: ik, rk
  use mod_field, only: Field

  implicit none

  integer(ik), parameter :: im = 10, jm = 10
  type(Field) :: u

  ! creating the field
  u = Field('x-component of velocity', [im, jm])

  print *, 'Field name: ', u % name
  print *, 'Lower bounds: ', u % lb
  print *, 'Upper bounds: ', u % ub
  print *, 'Neighbors (E, W, S, N): ', u % neighbors
  print *, 'Edge size: ', u % edge_size
  print *, 'Is allocated: ', allocated(u % data)

  ! test assignment operator
  print *, 'Field values before assignment: ', u % data
  u = 10._rk
  print *, 'Field values after real assignment: ', u % data

  call u % init_gaussian(0.07, 3, 3)
  print *, 'Field values after gaussian init: ', u % data

end program test_field
