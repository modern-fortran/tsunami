program tsunami_dt

  ! Tsunami simulator.
  !
  ! Solves the non-linear 2-d shallow water equation system:
  !
  !     du/dt + u du/dx + v du/dy + g dh/dx = 0
  !     dv/dt + u dv/dx + v dv/dy + g dh/dy = 0
  !     dh/dt + d(hu)/dx + d(hv)/dy = 0
  !
  ! This version is parallelized, and uses derived types.

  use mod_config, only: Config, config_from_namelist
  use mod_field, only: Field, diffx, diffy
  use mod_kinds, only: ik, rk

  implicit none

  integer(ik) :: n

  real(rk), parameter :: g = 9.8 ! gravitational acceleration [m/s^2]
  integer(ik), parameter :: ic = 51, jc = 51
  real(rk), parameter :: decay = 0.02

  type(Config) :: conf
  type(Field) :: h, u, v, hm

  if (this_image() == 1) print *, 'Tsunami started'

  conf = config_from_namelist('../tsunami.nml')

  associate (im => conf % grid_size_x, jm => conf % grid_size_y, &
             nm => conf % num_time_steps, dt => conf % time_step, &
             dx => conf % grid_spacing_x, dy => conf % grid_spacing_y)

  u = Field('u', [im, jm])
  v = Field('v', [im, jm])
  h = Field('h', [im, jm])
  hm = Field('hm', [im, jm])

  ! initialize a gaussian blob in the center
  call h % init_gaussian(decay, ic, jc)

  ! set mean water depth
  hm = 10.0

  call h % write(0)

  time_loop: do n = 1, nm

    if (this_image() == 1) then
      print *, 'Computing time step', n, '/', nm
    end if

    ! compute u at next time step
    u = u - (u * diffx(u) / dx + v * diffy(u) / dy &
      + g * diffx(h) / dx) * dt

    ! compute v at next time step
    v = v - (u * diffx(v) / dx + v * diffy(v) / dy &
      + g * diffy(h) / dy) * dt

    ! compute h at next time step
    h = h - (diffx(u * (hm + h)) / dx + diffy(v * (hm + h)) / dy) * dt

    call h % write(n)

  end do time_loop

  end associate

end program tsunami_dt
