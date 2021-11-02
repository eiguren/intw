!---------------------------------------------------------------------------------------------
subroutine realselfe_eliash_qe_t0(wpts,w_min,w_max,d_w,epts,e_min,e_max,d_e,eliash_f,realse)
!---------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 17/03/2016
!
!========================================================================================!
! We calculate the real part of the self-energy for the state (k,n) from the Eliashberg  !
! EPC Quasi-Elastic Function we know from previous calculation. For Quasi-Elastic case:  !
!                                                                                        !
! Re[SE(k,n,e)]=Integrate[a^2F(k,n,w)*Re[L(w,e)],{w,0,w_max}]                            !
!                                                                                        !
! Re[L(w,e)]=0.5*Log[((e-w)^2+sig^2)/((e+w)^2+sig^2)]                                    !
!                                                                                        !
! We make this integration using the Composite Simpson's rule:                           !
!                                                                                        !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                      !
!                                                                                        !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                      !
!                                                                                        !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                      !
!========================================================================================!

  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables
  
  integer,intent(in) :: wpts,epts
  real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
  real(dp),intent(in) :: eliash_f(wpts)
  real(dp),intent(inout) :: realse(epts)

  !local variables

  integer :: ifs,i,j,ipts,jpts,ipol,jpol
  real(dp) :: omega,eps,sigma_e,rel

  realse(:)=0.0d0
  !
  sigma_e=20.0d0*d_e
  !
  do ipts=1,epts
     !
     eps=e_min+d_e*(ipts-1)
     !
     j=0
     !
     do jpts=2,wpts
        !
        j=j+1
        !
        omega=w_min+d_w*(jpts-1)
        !
        rel=0.5d0*log(((eps-omega)**2+sigma_e**2)/((eps+omega)**2+sigma_e**2))
        !
        if (jpts.gt.wpts-1) then
           !
           realse(ipts)=realse(ipts)+eliash_f(jpts)*rel
           exit
           !
        endif
        !
        if ((j/2)*2==j) then
           !
           realse(ipts)=realse(ipts)+2*eliash_f(jpts)*rel
           !
        else
           !
           realse(ipts)=realse(ipts)+4*eliash_f(jpts)*rel
           !
        endif
        !
     enddo !jpts
     !
     realse(ipts)=(1.d0/3.d0)*realse(ipts)*d_w
     !
  enddo !ipts
  !
  return

end subroutine realselfe_eliash_qe_t0
!----------------------------------------------------------------------------------------------
!**********************************************************************************************
!----------------------------------------------------------------------------------------------
subroutine imagselfe_eliash_qe_t0(wpts,w_min,w_max,d_w,epts,e_min,e_max,d_e,eliash_f,imagse)
!----------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 17/03/2016
!
!============================================================================================!
! We calculate the imaginary part of the self-energy for the state (k,n) from the Eliashberg !
! EPC Quasi-Elastic Function we know from previous calculation. For Quasi-Elastic case:      !
!                                                                                            !
! Im[SE(k,n,e)]=Integrate[a^2F(k,n,w)*Im[L(w,e)],{w,0,w_max}]                                !
!                                                                                            !
! Im[L(w,e)]=Pi*(Tita[e-w]+Tita[-(e+w)])                                                              !
!                                                                                            !
! We make this integration using the Composite Simpson's rule:                               !
!                                                                                            !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                          !
!                                                                                            !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                          !
!                                                                                            !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                          !
!============================================================================================!

  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables
  
  integer,intent(in) :: wpts,epts
  real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
  real(dp),intent(in) :: eliash_f(wpts)
  real(dp),intent(inout) :: imagse(epts)

  !local variables

  integer :: ifs,i,j,ipts,jpts,ipol,jpol
  real(dp) :: omega,eps,sigma_e,iml

  imagse(:)=0.0d0
  !
  sigma_e=20.0d0*d_e
  !
  do ipts=1,epts
     !
     eps=e_min+d_e*(ipts-1)
     !
     j=0
     !
     do jpts=2,wpts
        !
        j=j+1
        !
        omega=w_min+d_w*(jpts-1)
        !
        iml=0.0d0
        !
        if (eps.ge.omega) then
           iml=iml+pi
        elseif (eps.lt.-1.d0*omega) then
           iml=iml+pi
        endif
        !
        if (jpts.gt.wpts-1) then
           !
           imagse(ipts)=imagse(ipts)+eliash_f(jpts)*iml
           exit
           !
        endif
        !
        if ((j/2)*2==j) then
           !
           imagse(ipts)=imagse(ipts)+2*eliash_f(jpts)*iml
           !
        else
           !
           imagse(ipts)=imagse(ipts)+4*eliash_f(jpts)*iml
           !
        endif
        !
     enddo !jpts
     !
     imagse(ipts)=(1.d0/3.d0)*imagse(ipts)*d_w
     !
  enddo !ipts
  !
  return

end subroutine imagselfe_eliash_qe_t0
!-----------------------------------------------------------------------------------------------------------------
!*****************************************************************************************************************
!-----------------------------------------------------------------------------------------------------------------
subroutine realselfe_eliash_full_t0(wpts,w_min,w_max,d_w,epts,e_min,e_max,d_e,eliash_fplus,eliash_fminus,realse)
!-----------------------------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 17/03/2016
!
!========================================================================================!
! We calculate the real part of the self-energy for the state (k,n) from the Eliashberg  !
! EPC Function we know from previous calculation.                                        !
!                                                                                        !
! Re[SE(k,n,e)]=0.5Integrate[a^2F-(k,n,w)*Log[(e-w)^2+sig^2]-                            !
!                            a^2F+(k,n,w)*Log[(e+w)^2+sig^2],{w,0,w_max}]                !
!                                                                                        !
! We make this integration using the Composite Simpson's rule:                           !
!                                                                                        !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                      !
!                                                                                        !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                      !
!                                                                                        !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                      !
!========================================================================================!

  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables
  
  integer,intent(in) :: wpts,epts
  real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
  real(dp),intent(in) :: eliash_fplus(wpts),eliash_fminus(wpts)
  real(dp),intent(inout) :: realse(epts)

  !local variables

  integer :: ifs,i,j,ipts,jpts,ipol,jpol
  real(dp) :: omega,eps,sigma_e,relminus,relplus,relminust,relplust

  realse(:)=0.0d0
  !
  sigma_e=100.0d0*d_e
  !
  do ipts=1,epts
     !
     eps=e_min+d_e*(ipts-1)
     !
     j=0
     !
     do jpts=2,wpts
        !
        j=j+1
        !
        omega=w_min+d_w*(jpts-1)
        !
        if (eps.gt.0.d0) then
           !
           relminus=0.5d0*log((abs(eps-omega)**2+sigma_e**2)/(abs(eps+omega)**2+sigma_e**2))
           relplus=0.0d0
           !
        else
           !
           relminus=0.0d0
           relplus=0.5d0*log((abs(eps+omega)**2+sigma_e**2)/(abs(eps-omega)**2+sigma_e**2))
           !
        endif
!        relplust=0.5d0*log(abs(eps-e_min)**2+sigma_e**2)
!        relminust=0.5d0*log(abs(eps-e_max)**2+sigma_e**2)
        relplust=0.0d0
        relminust=0.0d0
        !
        if (jpts.gt.wpts-1) then
           !
           realse(ipts)=realse(ipts)+&
           (eliash_fminus(jpts)*relminus-eliash_fplus(jpts)*relplus)+&
           (eliash_fplus(jpts)*relplust-eliash_fminus(jpts)*relminust)
           exit
           !
        endif
        !
        if ((j/2)*2==j) then
           !
           realse(ipts)=realse(ipts)+&
           2*(eliash_fminus(jpts)*relminus-eliash_fplus(jpts)*relplus)+&
           2*(eliash_fplus(jpts)*relplust-eliash_fminus(jpts)*relminust)
           !
        else
           !
           realse(ipts)=realse(ipts)+&
           4*(eliash_fminus(jpts)*relminus-eliash_fplus(jpts)*relplus)+&
           4*(eliash_fplus(jpts)*relplust-eliash_fminus(jpts)*relminust)
           !
        endif
        !
     enddo !jpts
     !
     realse(ipts)=(1.d0/3.d0)*realse(ipts)*d_w
     !
  enddo !ipts
  !
  return

end subroutine realselfe_eliash_full_t0
!------------------------------------------------------------------------------------------------------------------
!******************************************************************************************************************
!------------------------------------------------------------------------------------------------------------------
subroutine imagselfe_eliash_full_t0(wpts,w_min,w_max,d_w,epts,e_min,e_max,d_e,eliash_fplus,eliash_fminus,imagse)
!------------------------------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 17/03/2016
!
!============================================================================================!
! We calculate the imaginary part of the self-energy for the state (k,n) from the Eliashberg !
! EPC Function we know from previous calculation.                                            !
!                                                                                            !
! Im[SE(k,n,e)]=Pi*Integrate[a^2F+(k,n,w)*Tita[-(e+w)]+                                      !
!                            a^2F-(k,n,w)*Tita[e-w],{w,0,w_max}]                             !
!                                                                                            !
! We make this integration using the Composite Simpson's rule:                               !
!                                                                                            !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                          !
!                                                                                            !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                          !
!                                                                                            !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                          !
!============================================================================================!

  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables
  
  integer,intent(in) :: wpts,epts
  real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
  real(dp),intent(in) :: eliash_fplus(wpts),eliash_fminus(wpts)
  real(dp),intent(inout) :: imagse(epts)

  !local variables

  integer :: ifs,i,j,ipts,jpts,ipol,jpol
  real(dp) :: omega,eps,sigma_e,imlplus,imlminus

  imagse(:)=0.0d0
  !
  sigma_e=40.0d0*d_e
  !
  do ipts=1,epts
     !
     eps=e_min+d_e*(ipts-1)
     !
     j=0
     !
     do jpts=2,wpts
        !
        j=j+1
        !
        omega=w_min+d_w*(jpts-1)
        !
        imlplus=0.0d0
        imlminus=0.0d0
        !
        if (eps.ge.omega) then
           imlminus=imlminus+pi
        elseif (eps.lt.-1.d0*omega) then
           imlplus=imlplus+pi
        endif
        !
        if (jpts.gt.wpts-1) then
           !
           imagse(ipts)=imagse(ipts)+(eliash_fminus(jpts)*imlminus+eliash_fplus(jpts)*imlplus)
           exit
           !
        endif
        !
        if ((j/2)*2==j) then
           !
           imagse(ipts)=imagse(ipts)+2*(eliash_fminus(jpts)*imlminus+eliash_fplus(jpts)*imlplus)
           !
        else
           !
           imagse(ipts)=imagse(ipts)+4*(eliash_fminus(jpts)*imlminus+eliash_fplus(jpts)*imlplus)
           !
        endif
        !
     enddo !jpts
     !
     imagse(ipts)=(1.d0/3.d0)*imagse(ipts)*d_w
     !
  enddo !ipts
  !
  return

end subroutine imagselfe_eliash_full_t0
!-------------------------------------------------------------------------------------------------
