program test_sync_on_assign

  use mod_kinds, only: ik, rk
  use mod_field, only: Field, from_field

  implicit none

  integer(ik), parameter :: im = 4, jm = 4 ! grid size in x and y

  type(Field) :: a, b 

  a = Field('Source field', [im, jm])
  b = Field('Target field', [im, jm])
  a % data(a % lb(1):a % ub(1), a % lb(2):a % ub(2)) = this_image()
  b = a
  print *, 'after', this_image(), b % data

end program test_sync_on_assign
