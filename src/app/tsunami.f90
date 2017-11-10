program tsunami

! Tsunami simulator
!
! This version solves the linearized 1-d advection equation:
!
!     du/dt + c du/dx = 0
!
! Discretized as:
!
!     u^{n+1}_i = u^n_i - c * dt * (u^n_{i} - u^n_{i-1}) / dx

implicit none

integer, parameter :: idm = 100 ! grid size in x
integer, parameter :: ndm = 100 ! number of time steps

real, parameter :: dt = 1 ! time step [s]
real, parameter :: dx = 1 ! grid spacing [m]
real, parameter :: c = 1 ! phase speed [m/s]

real, parameter :: pi = 3.14159256

integer :: i,n

real, dimension(idm) :: du, u, u_init

! initialize a gaussian blob centered at i = 25
do i = 1, idm
  u(i) = exp(-0.02*(i-25)**2)
end do

! store initial state
u_init = u

! write initial state to screen
write(*,*)0, u

! time loop
do n = 1, ndm

  ! calculate the upstream difference of h in x
  do i = 2, idm
    du(i) = u(i) - u(i-1)
  end do

  ! apply periodic boundary condition on the left 
  du(1) = u(1) - u(idm)

  ! compute u at next time step
  do i = 1, idm
    u(i) = u(i) - c * du(i) / dx * dt
  end do

  ! write current state to screen
  write(*,*)n, u

enddo

end program tsunami
