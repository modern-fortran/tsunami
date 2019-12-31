module mod_boundary

  ! Provides subroutines to deal with boundary conditions.

  use mod_kinds, only: ik, rk

  implicit none

  private
  public :: reflective_boundary

contains

  subroutine reflective_boundary(u, v, h)
    real(rk), allocatable, intent(in out) :: u(:,:), v(:,:), h(:,:)
    integer(ik) :: ims, ime, jms, jme
    associate(ims => lbound(u, dim=1), ime => ubound(u, dim=1),&
              jms => lbound(u, dim=2), jme => ubound(u, dim=2))
      u([ims, ime],:) = 0
      u(:,[jms, jme]) = 0
      v([ims, ime],:) = 0
      v(:,[jms, jme]) = 0
      h([ims, ime],:) = h([ims+1, ime-1],:)
      h(:,[jms, jme]) = h(:,[jms+1, jme-1])
    end associate
  end subroutine reflective_boundary

end module mod_boundary
