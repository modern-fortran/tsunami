program tsunami

implicit none

integer :: i, n ! indices in space and time
integer :: im ! grid size in space
integer :: nm ! number of time steps

real :: dt ! time step [s]
real :: dx ! grid spacing [m]
real :: c ! phase speed [m/s]

im = 100
nm = 100

dt = 1.
dx = 1.
c = 1.

end program tsunami
