module mod_io

  ! Provides input/output procedures.

  use iso_fortran_env, only: int32, real32

  implicit none

  private
  public :: write_field

contains

  subroutine write_field(field, fieldname, time)
    ! Writes a field into a binary file.
    real(real32), intent(in) :: field(:,:)
    character(*), intent(in) :: fieldname
    integer(int32), intent(in) :: time
    integer(int32) :: fileunit, record_length
    character(100) :: filename, timestr
    write(timestr, '(i4.4)') time
    filename = 'tsunami_' // fieldname // '_' // trim(timestr) // '.dat'
    record_length = storage_size(field) / 8 * size(field)
    open(newunit=fileunit, file=filename, access='direct', recl=record_length)
    write(unit=fileunit, rec=1) field
    close(fileunit)
  end subroutine write_field

end module mod_io
