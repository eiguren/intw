!----------------------------------------------------------------------------
subroutine imselfen2reselfen_t0(epts,e_min,e_max,d_e,imselfen_w,reselfen_w)
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
  real(dp),intent(in) :: e_min,e_max,d_e
  real(dp),intent(in) :: imselfen_w(epts)
  real(dp),intent(inout) :: reselfen_w(epts)

  !local variables

  real(dp) :: e,dimselfen_w(epts),integral,tita,modz
  integer :: ipts,jpts,i
  complex(dp) :: z
  real(dp) :: omega,oomega

  ! Derivative of the Imaginary part
  !
  dimselfen_w(:)=0.d0
  !
  do ipts=2,epts-1
     !
     dimselfen_w(ipts)=(imselfen_w(ipts+1)-imselfen_w(ipts-1))/(2*d_e)
     !
  enddo !ipts
  !
  dimselfen_w(1)=dimselfen_w(2)
  dimselfen_w(epts)=dimselfen_w(epts-1)
  !
  ! We make Kramers-Kronig relation by means of Asier's trick
  !
  !$omp parallel default(none) &
  !$omp shared(epts,e_min,d_e,dimselfen_w,reselfen_w,pi,cmplx_0) &
  !$omp private(ipts,i,omega,integral,jpts,oomega,z,modz,tita)
  !
  !$omp do
  !
  do ipts=1,epts
     !
     i=0
     !
     omega=e_min+d_e*(ipts-1)
     integral=cmplx_0
     !
     do jpts=1,epts
        !
        i=i+1
        !
        oomega=e_min+d_e*(jpts-1)
        !
        z=cmplx(omega-oomega,d_e)
        modz=abs(z)
        tita=aimag(z)/real(z)
        !
        if (jpts.lt.2.or.jpts.gt.epts-1) then
           !
           integral=integral+dimselfen_w(jpts)*real(log(z))
           !
        else
           !
           if ((i/2)*2==i) then
              !
              integral=integral+2*dimselfen_w(jpts)*real(log(z))
              !
           else
              !
              integral=integral+4*dimselfen_w(jpts)*real(log(z))
              !
           endif
           !
        endif
        !
     enddo !jpts
     !
     reselfen_w(ipts)=integral*d_e/(pi*3.d0)
!     reselfen_w(ipts)=dimselfen_w(ipts)
     !
  enddo !ipts
  !
  !$omp end parallel
  !
  reselfen_w(:)=reselfen_w(:)-reselfen_w(((epts-1)/2)+1)
  !
  return

end subroutine imselfen2reselfen_t0
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
!---------------------------------------------------------------------------------------------
