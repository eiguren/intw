integer function simple_test_test() result(r)

  use kinds, only: dp
  use intw_test_module, only: assert

  implicit none

  real(kind=dp), parameter :: precision = 0.0001_dp
  !
  integer, parameter :: reference_calculation_result = 8
  !
  integer :: num1, num2, result

  ! Execute some program
  num1 = 5
  num2 = 3
  result = num1 + num2

  ! Compare the results with previous reference calculations
  if ( assert(result, reference_calculation_result, precision) ) then
    ! Exit with error code 0 if results are correct
    r = 0
  else
    ! Exit with error code 1 if results are wrong
    r = 1
  end if

  return

end function
