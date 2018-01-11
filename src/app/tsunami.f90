program tsunami

use iso_fortran_env, only: int32, real32

implicit none

integer(kind=int32) :: i, n

integer(kind=int32), parameter :: im = 100 ! grid size in space
integer(kind=int32), parameter :: nm = 100 ! number of time steps

real(kind=real32), parameter :: dt = 1 ! time step [s]
real(kind=real32), parameter :: dx = 1 ! grid spacing [m]
real(kind=real32), parameter :: c = 1 ! phase speed [m/s]

real(kind=real32), dimension(im) :: du, u

end program tsunami
