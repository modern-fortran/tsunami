module mod_config

  ! Provides input/output procedures.

  use mod_kinds, only: ik, rk

  implicit none

  private
  public :: Config, config_from_namelist

  type :: Config
    integer :: grid_size_x
    integer :: grid_size_y
    integer :: num_time_steps
    real :: time_step
    real :: grid_spacing_x
    real :: grid_spacing_y
  end type Config

contains

  function config_from_namelist(filename) result(conf)
    character(len=*), intent(in) :: filename
    type(Config) :: conf
    integer :: fileunit
    integer :: grid_size_x, grid_size_y, num_time_steps
    real :: time_step, grid_spacing_x, grid_spacing_y
    namelist /domain/ grid_size_x, grid_size_y, num_time_steps, time_step, &
                      grid_spacing_x, grid_spacing_y
    open(newunit=fileunit, file=filename, status='old', action='read')
    read(fileunit, nml=domain)
    close(fileunit)
    conf = Config(grid_size_x, grid_size_y, num_time_steps, time_step, &
                  grid_spacing_x, grid_spacing_y)
  end function config_from_namelist

end module mod_config
