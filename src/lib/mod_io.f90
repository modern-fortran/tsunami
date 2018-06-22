module mod_io

  implicit none

  private
  public :: write_field
  
contains

  subroutine write_field(field, fieldname, time)
    ! Writes a field into a binary file.
    real, allocatable, intent(in) :: field(:,:)
    character(len=*), intent(in) :: fieldname
    integer, intent(in) :: time
    integer :: fileunit, record_length
    character(len=100) :: filename, timestr
    write(timestr, '(i4.4)') time
    filename = 'tsunami_' // fieldname // '_' // trim(timestr) // '.dat'
    record_length = storage_size(field) / 8 * size(field)
    write(*,*) trim(filename), record_length, storage_size(field), size(field)
    open(newunit=fileunit, file=filename, access='direct', recl=record_length)
    write(unit=fileunit, rec=1) field
    close(fileunit)
  end subroutine write_field

end module mod_io
