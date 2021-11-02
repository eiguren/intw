!----------------------------------------------------------------------------
subroutine spectral_function_im_re_selfen(epts,e_min,e_max,d_e,ekk,ef,imselfen_w,reselfen_w,a_w)
!----------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 31/03/2017
!
!=========================================================================!
! We calculate the Real part of the Selfenergy at T=0 K once we know the  !
! Imaginary part of SE calculated using the direct many-body formula      !
! (NO ELIASHBERG FORMALISM) by means of Asier Eiguren thought "trick" of  !
! resolving the Kramers-Kronig relation                                   !
!=========================================================================!
!
!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables

  integer,intent(in) :: epts
  real(dp),intent(in) :: e_min,e_max,d_e,ekk,ef
  real(dp),intent(in) :: imselfen_w(epts)
  real(dp),intent(in) :: reselfen_w(epts)
  real(dp),intent(inout) :: a_w(epts)

  !local variables

  integer :: ipts
  real(dp) :: omega

  do ipts=1,epts
     !
     omega=e_min+d_e*(ipts-1)
     !
     a_w(ipts)=(1/pi)*imselfen_w(ipts)/((omega-(ekk-ef)-reselfen_w(ipts))**2.d0+imselfen_w(ipts)**2.d0)
     !
  enddo !ipts
  !
  return

end subroutine spectral_function_im_re_selfen
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
!---------------------------------------------------------------------------------------------
