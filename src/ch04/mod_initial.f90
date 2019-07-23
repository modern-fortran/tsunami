module mod_initial

  use iso_fortran_env, only: int32, real32
  implicit none

  private
  public :: set_gaussian

contains

  pure subroutine set_gaussian(x, icenter, decay)
    ! Sets the values of x to a Gaussian shape centered on icenter
    ! that decays with the given input decay.
    real(real32), intent(in out) :: x(:)
    integer(int32), intent(in) :: icenter
    real(real32), intent(in) :: decay
    integer(int32) :: i
    do concurrent(i = 1:size(x))
      x(i) = exp(-decay * (i - icenter)**2)
    end do
  end subroutine set_gaussian

end module mod_initial
