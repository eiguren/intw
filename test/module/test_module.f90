module intw_test_module

  use kinds, only: dp

  implicit none

  public :: assert

  private

  interface assert
    module procedure assert_logical, assert_integer, assert_integer_vec, &
                     assert_real, assert_real_vec
  end interface assert


contains

  logical function assert_logical(a,b,prec) result(r)

    use kinds, only: dp

    logical, intent(in) :: a, b
    real(kind=dp), intent(in) :: prec

    r = .true.
    if (a .neqv. b) then
      r = .false.
    endif

  end function


  logical function assert_integer(a,b,prec) result(r)

    use kinds, only: dp

    integer, intent(in) :: a, b
    real(kind=dp), intent(in) :: prec

    r = .true.
    if (a /= b) then
      r = .false.
    endif

  end function

  logical function assert_integer_vec(a,b,prec) result(r)

    use kinds, only: dp

    integer, intent(in), dimension(:) :: a, b
    real(kind=dp), intent(in) :: prec

    integer :: i

    if (size(a) /= size(b)) then
      r = .false.
      write(6,*) "WARNING: assert_integer_vec: different size vectors"
      return
    endif

    r = .true.
    do i=1,size(a)
      if (a(i) /= b(i)) then
        r = .false.
        return
      endif
    enddo

  end function

  logical function assert_real(a,b,prec) result(r)

    use kinds, only: dp

    real(kind=dp), intent(in) :: a, b
    real(kind=dp), intent(in) :: prec

    r = .true.
    if (abs(a-b) > prec) then
      r = .false.
    endif

  end function

  logical function assert_real_vec(a,b,prec) result(r)

    use kinds, only: dp

    real(kind=dp), intent(in), dimension(:) :: a, b
    real(kind=dp), intent(in) :: prec

    integer :: i

    if (size(a) /= size(b)) then
      r = .false.
      write(6,*) "WARNING: assert_real_vec: different size vectors"
      return
    endif

    r = .true.
    do i=1,size(a)
      if (abs(a(i)-b(i)) > prec) then
        r = .false.
        return
      endif
    enddo

  end function

end module intw_test_module
