!----------------------------------------------------------------------------------
subroutine gaussian(e0,e_min,e_max,n,sigma,gauss)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!==========================================================================================!
! We create a list containing values for Delta[e-e0] simulated by the gaussian function    !
! with e defined by us from e_min to e_max, and with n points, each on separated from its  !
! neighbours by abs[e_max-e_min]/n                                                         !
!==========================================================================================!

  use kinds, only: dp
  use intw_useful_constants, only: pi

  implicit none

  !I/O variables

  integer,intent(in) :: n
  real(dp),intent(in) :: e0,e_min,e_max,sigma
  real(dp),intent(inout) :: gauss(n)

  !local variables

  integer :: i
  real(dp) :: e

  gauss(:)=0.d0
  !
  do i=1,n
     !
     e=e_min+(e_max-e_min)*(i-1)/(n-1)
     !
     if (abs(e-e0).gt.(4.d0*sigma)) then
        gauss(i)=0.d0
     else
        gauss(i)=exp(-((e-e0)/sigma)**2.d0)/(sigma*sqrt(pi))
     endif
     !
  enddo
  !
  return

end subroutine gaussian
