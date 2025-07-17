module pw2intw_utils

  implicit none

  ! functions and subroutines
  public :: write_tag

  private

contains

  subroutine write_tag(string,i,tag)
    !-----------------------------------------------
    ! This subroutine creates a character string of
    ! the form "string"integer, where the integer
    ! will be immediately after the end of "string",
    ! without blank spaces.
    !-----------------------------------------------
    implicit none

    integer :: i
    character(*) :: string
    character(256) :: integer_part, tag


    if (i < 10) then
      write(integer_part,100) i
    elseif (i < 100 ) then
      write(integer_part,200) i
    elseif (i < 1000 ) then
      write(integer_part,300) i
    elseif (i < 10000 ) then
      write(integer_part,400) i
    elseif (i < 100000 ) then
      write(integer_part,500) i
    end if

    tag = trim(string)//trim(integer_part)

    100 format(I1)
    200 format(I2)
    300 format(I3)
    400 format(I4)
    500 format(I5)

  end subroutine write_tag

end module pw2intw_utils