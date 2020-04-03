program tsunami

  ! Tsunami simulator.
  !
  ! Solves the non-linear 2-d shallow water equation system:
  !
  !     du/dt + u du/dx + v du/dy + g dh/dx = 0
  !     dv/dt + u dv/dx + v dv/dy + g dh/dy = 0
  !     dh/dt + d(hu)/dx + d(hv)/dy = 0
  !
  ! This version is parallelized and uses derived types.

  use iso_fortran_env, only: int32, real32
  use mod_field, only: Field, diffx, diffy

  implicit none

  integer(int32) :: n

  integer(int32), parameter :: im = 201, jm = 201
  integer(int32), parameter :: num_time_steps = 1000

  real(real32), parameter :: dt = 0.02
  real(real32), parameter :: dx = 1, dy = 1
  real(real32), parameter :: g = 9.8

  integer(int32), parameter :: ic = im / 2 + 1
  integer(int32), parameter :: jc = jm / 2 + 1
  real(real32), parameter :: decay = 0.02

  type(Field) :: h, hm, u, v

  real(real32) :: hmin, hmax, hmean

  u = Field('u', [im, jm])
  v = Field('v', [im, jm])
  h = Field('h', [im, jm])
  hm = Field('h_mean', [im, jm])

  ! initialize a gaussian blob in the center
  call h % init_gaussian(decay, ic, jc)

  hm = 10.

  call h % write(0)

  time_loop: do n = 1, num_time_steps

    ! compute u at next time step
    u = u - (u * diffx(u) / dx + v * diffy(u) / dy &
      + g * diffx(h) / dx) * dt

    ! compute v at next time step
    v = v - (u * diffx(v) / dx + v * diffy(v) / dy &
      + g * diffy(h) / dy) * dt

    ! compute h at next time step
    h = h - (diffx(u * (hm + h)) / dx &
           + diffy(v * (hm + h)) / dy) * dt

    hmin = minval(h % data)
    call co_min(hmin, 1)

    hmax = maxval(h % data)
    call co_max(hmax, 1)

    hmean = sum(h % data(h % lb(1):h % ub(1),h % lb(2):h % ub(2))) &
          / size(h % data(h % lb(1):h % ub(1),h % lb(2):h % ub(2)))
    call co_sum(hmean, 1)
    hmean = hmean / num_images()

    if (this_image() == 1) &
      print '(a, i5, 3(f10.6))', &
        'step, min(h), max(h), mean(h):', &
         n, hmin, hmax, hmean

    call h % write(n)

  end do time_loop

end program tsunami
