program test_field

  use mod_kinds, only: ik, rk
  use mod_field, only: Field

  implicit none

  integer(ik), parameter :: im = 100, jm = 100
  type(Field) :: u

  u = Field('x-component of velocity', [im, jm])

  print *, 'Field name: ', u % name
  print *, 'Lower bounds: ', u % lb
  print *, 'Upper bounds: ', u % ub
  print *, 'Neighbors (E, W, S, N): ', u % neighbors
  print *, 'Is allocated: ', allocated(u % data)

end program test_field
