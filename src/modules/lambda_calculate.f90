!---------------------------------------------------------------------------------------
subroutine lambda_calculate_qe_t0_so(wpts,w_min,w_max,d_w,eliash_f_so,lambda_so)
!---------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 15/03/2016
!
!===================================================================================!
! We calculate the EPC parameter lambda for the state (k,n) which Eliashberg EPC    !
! Quasi-Elastic Function we know from previous calculation. For Quasi-Elastic case: !
!                                                                                   !
! l(k,n)=2*Integrate[a^2F(k,n,w)/w,{w,0,w_max}]                                     !
!                                                                                   !
! We make this integration using the Composite Simpson's rule:                      !
!                                                                                   !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                 !
!                                                                                   !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                 !
!                                                                                   !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                 !
!                                                                                   !
! So that:                                                                          !
!                                                                                   !
! l(k,n)=(2*d_w/3)*(4*I+2*P+a^2F(k,n,w_max)/w_max)                                  !
!                                                                                   !
! Where we have taken into account that lim_w->0[a^2F(k,n,w)/w]=0                   !
!===================================================================================!

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables

  integer,intent(in) :: wpts
  real(dp),intent(in) :: w_min,w_max,d_w
  real(dp),intent(in) :: eliash_f_so(wpts,nspin,nspin)
  real(dp),intent(inout) :: lambda_so(nspin*nspin)

  !local variables

  real(dp) :: eli(wpts)
  integer :: ifs,i,ipts,nfs
  real(dp) :: omega

  nfs=nspin*nspin
  !
  do ifs=1,nfs
     !
     if (ifs==1) then
        !
        eli(:)=eliash_f_so(:,1,1)
        !
     elseif (ifs==2) then
        !
        eli(:)=eliash_f_so(:,1,2)
        !
     elseif (ifs==3) then
        !
        eli(:)=eliash_f_so(:,2,1)
        !
     elseif (ifs==4) then
        !
        eli(:)=eliash_f_so(:,2,2)
        !
     endif
     !
     i=0
     !
     lambda_so(ifs)=0.0d0
     !
     do ipts=2,wpts
        !
        i=i+1
        !
        omega=w_min+d_w*(ipts-1)
        !
        if (ipts.gt.wpts-1) then
           !
           lambda_so(ifs)=lambda_so(ifs)+(eli(ipts)/omega)
           exit
           !
        endif
        !
        if ((i/2)*2==i) then
           !
           lambda_so(ifs)=lambda_so(ifs)+2*(eli(ipts)/omega)
           !
        else
           !
           lambda_so(ifs)=lambda_so(ifs)+4*(eli(ipts)/omega)
           !
        endif
        !
     enddo !ipts
  enddo !ifs
  !
  lambda_so(:)=(2.d0/3.d0)*lambda_so(:)*d_w
  !
  return

end subroutine lambda_calculate_qe_t0_so
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
!---------------------------------------------------------------------------------------------
subroutine lambda_calculate_qe_t0(wpts,w_min,w_max,d_w,eliash_f,lambda)
!---------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 15/03/2016
!
!===================================================================================!
! We calculate the EPC parameter lambda for the state (k,n) which Eliashberg EPC    !
! Quasi-Elastic Function we know from previous calculation. For Quasi-Elastic case: !
!                                                                                   !
! l(k,n)=2*Integrate[a^2F(k,n,w)/w,{w,0,w_max}]                                     !
!                                                                                   !
! We make this integration using the Composite Simpson's rule:                      !
!                                                                                   !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                 !
!                                                                                   !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                 !
!                                                                                   !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                 !
!                                                                                   !
! So that:                                                                          !
!                                                                                   !
! l(k,n)=(2*d_w/3)*(4*I+2*P+a^2F(k,n,w_max)/w_max)                                  !
!                                                                                   !
! Where we have taken into account that lim_w->0[a^2F(k,n,w)/w]=0                   !
!===================================================================================!

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables

  integer,intent(in) :: wpts
  real(dp),intent(in) :: w_min,w_max,d_w
  real(dp),intent(in) :: eliash_f(wpts)
  real(dp),intent(inout) :: lambda

  !local variables

  real(dp) :: eli(wpts)
  integer :: ifs,i,ipts,nfs
  real(dp) :: omega

  eli(:)=eliash_f(:)
  !
  i=0
  !
  lambda=0.0d0
  !
  do ipts=2,wpts
     !
     i=i+1
     !
     omega=w_min+d_w*(ipts-1)
     !
     if (ipts.gt.wpts-1) then
        !
        lambda=lambda+(eli(ipts)/omega)
        exit
        !
     endif
     !
     if ((i/2)*2==i) then
        !
        lambda=lambda+2*(eli(ipts)/omega)
        !
     else
        !
        lambda=lambda+4*(eli(ipts)/omega)
        !
     endif
     !
  enddo !ipts
  !
  lambda=(2.d0/3.d0)*lambda*d_w
  !
  return

end subroutine lambda_calculate_qe_t0
!--------------------------------------------------------------------------------------------------------
!********************************************************************************************************
!--------------------------------------------------------------------------------------------------------
subroutine lambda_calculate_full_t0_so(wpts,w_min,w_max,d_w,eliash_fplus_so,eliash_fminus_so,lambda_so)
!--------------------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 15/03/2016
!
!===================================================================================!
! We calculate the EPC parameter lambda for the state (k,n) which Eliashberg EPC    !
! Full Function we know from previous calculation. For Quasi-Elastic case:          !
!                                                                                   !
! l(k,n)=Integrate[(1/w)*(a^2F+(k,n,w)+a^2F-(k,n,w)),{w,0,w_max}]                                     !
!                                                                                   !
! We make this integration using the Composite Simpson's rule:                      !
!                                                                                   !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                 !
!                                                                                   !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                 !
!                                                                                   !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                 !
!                                                                                   !
! So that:                                                                          !
!                                                                                   !
! l(k,n)=(d_w/3)*(4*I+2*P+a^2F(k,n,w_max)/w_max)                                    !
!                                                                                   !
! Where we have taken into account that lim_w->0[a^2F+/-(k,n,w)/w]=0                !
!===================================================================================!

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables

  integer,intent(in) :: wpts
  real(dp),intent(in) :: w_min,w_max,d_w
  real(dp),intent(in) :: eliash_fplus_so(wpts,nspin,nspin)
  real(dp),intent(in) :: eliash_fminus_so(wpts,nspin,nspin)
  real(dp),intent(inout) :: lambda_so(nspin*nspin)

  !local variables

  real(dp) :: eliplus_so(wpts),eliminus_so(wpts)
  integer :: ifs,i,ipts,nfs
  real(dp) :: omega

  nfs=nspin*nspin
  !
  do ifs=1,nfs
     !
     if (ifs==1) then
        !
        eliplus_so(:)=eliash_fplus_so(:,1,1)
        eliminus_so(:)=eliash_fminus_so(:,1,1)
        !
     elseif (ifs==2) then
        !
        eliplus_so(:)=eliash_fplus_so(:,1,2)
        eliminus_so(:)=eliash_fminus_so(:,1,2)
        !
     elseif (ifs==3) then
        !
        eliplus_so(:)=eliash_fplus_so(:,2,1)
        eliminus_so(:)=eliash_fminus_so(:,2,1)
        !
     elseif (ifs==4) then
        !
        eliplus_so(:)=eliash_fplus_so(:,2,2)
        eliminus_so(:)=eliash_fminus_so(:,2,2)
        !
     endif
     !
     i=0
     !
     lambda_so(ifs)=0.0d0
     !
     do ipts=2,wpts
        !
        i=i+1
        !
        omega=w_min+d_w*(ipts-1)
        !
        if (ipts.gt.wpts-1) then
           !
           lambda_so(ifs)=lambda_so(ifs)+(eliplus_so(ipts)+eliminus_so(ipts))/omega
           exit
           !
        endif
        !
        if ((i/2)*2==i) then
           !
           lambda_so(ifs)=lambda_so(ifs)+2*(eliplus_so(ipts)+eliminus_so(ipts))/omega
           !
        else
           !
           lambda_so(ifs)=lambda_so(ifs)+4*(eliplus_so(ipts)+eliminus_so(ipts))/omega
           !
        endif
        !
     enddo !ipts
  enddo !ifs
  !
  lambda_so(:)=(1.d0/3.d0)*lambda_so(:)*d_w
  !
  return

end subroutine lambda_calculate_full_t0_so
!----------------------------------------------------------------------------------------------
!**********************************************************************************************
!----------------------------------------------------------------------------------------------
subroutine lambda_calculate_full_t0(wpts,w_min,w_max,d_w,eliash_fplus,eliash_fminus,lambda)
!----------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 15/03/2016
!
!===================================================================================!
! We calculate the EPC parameter lambda for the state (k,n) which Eliashberg EPC    !
! Full Function we know from previous calculation. For Quasi-Elastic case:          !
!                                                                                   !
! l(k,n)=Integrate[(1/w)*(a^2F+(k,n,w)+a^2F-(k,n,w)),{w,0,w_max}]                                     !
!                                                                                   !
! We make this integration using the Composite Simpson's rule:                      !
!                                                                                   !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                 !
!                                                                                   !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                 !
!                                                                                   !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                 !
!                                                                                   !
! So that:                                                                          !
!                                                                                   !
! l(k,n)=(d_w/3)*(4*I+2*P+a^2F(k,n,w_max)/w_max)                                    !
!                                                                                   !
! Where we have taken into account that lim_w->0[a^2F+/-(k,n,w)/w]=0                !
!===================================================================================!

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables

  integer,intent(in) :: wpts
  real(dp),intent(in) :: w_min,w_max,d_w
  real(dp),intent(in) :: eliash_fplus(wpts)
  real(dp),intent(in) :: eliash_fminus(wpts)
  real(dp),intent(inout) :: lambda

  !local variables

  real(dp) :: eliplus(wpts),eliminus(wpts)
  integer :: ifs,i,ipts,nfs
  real(dp) :: omega

  eliplus(:)=eliash_fplus(:)
  eliminus(:)=eliash_fminus(:)
  !
  i=0
  !
  lambda=0.0d0
  !
  do ipts=2,wpts
     !
     i=i+1
     !
     omega=w_min+d_w*(ipts-1)
     !
     if (ipts.gt.wpts-1) then
        !
        lambda=lambda+(eliplus(ipts)+eliminus(ipts))/omega
        exit
        !
     endif
     !
     if ((i/2)*2==i) then
        !
        lambda=lambda+2*(eliplus(ipts)+eliminus(ipts))/omega
        !
     else
        !
        lambda=lambda+4*(eliplus(ipts)+eliminus(ipts))/omega
        !
     endif
     !
  enddo !ipts
  !
  lambda=(1.d0/3.d0)*lambda*d_w
  !
  return

end subroutine lambda_calculate_full_t0
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
subroutine hlambda_calculate_full_t0(wpts,w_min,w_max,d_w,eliash_fplus,lambda)
!----------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 15/03/2016
!
!===================================================================================!
! We calculate the EPC parameter lambda for the state (k,n) which Eliashberg EPC    !
! Full Function we know from previous calculation. For Quasi-Elastic case:          !
!                                                                                   !
! l(k,n)=Integrate[(1/w)*(a^2F+(k,n,w)+a^2F-(k,n,w)),{w,0,w_max}]                                     !
!                                                                                   !
! We make this integration using the Composite Simpson's rule:                      !
!                                                                                   !
! Integrate[f(x),{x,a,b}]=(h/3)*(f(a)+4*I+2*P+f(b))                                 !
!                                                                                   !
! I=Sum_i_odd_1,N-1_[f(xi)]=f(x1)+f(x3)+...+f(xN-1)                                 !
!                                                                                   !
! P=Sum_i_odd_2,N-2_[f(xi)]=f(x2)+f(x4)+...+f(xN-2)                                 !
!                                                                                   !
! So that:                                                                          !
!                                                                                   !
! l(k,n)=(d_w/3)*(4*I+2*P+a^2F(k,n,w_max)/w_max)                                    !
!                                                                                   !
! Where we have taken into account that lim_w->0[a^2F+/-(k,n,w)/w]=0                !
!===================================================================================!

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants

  implicit none

  !I/O variables

  integer,intent(in) :: wpts
  real(dp),intent(in) :: w_min,w_max,d_w
  real(dp),intent(in) :: eliash_fplus(wpts)
  real(dp),intent(inout) :: lambda

  !local variables

  real(dp) :: eliplus(wpts)
  integer :: ifs,i,ipts,nfs
  real(dp) :: omega

  eliplus(:)=eliash_fplus(:)
  !
  i=0
  !
  lambda=0.0d0
  !
  do ipts=2,wpts
     !
     i=i+1
     !
     omega=w_min+d_w*(ipts-1)
     !
     if (ipts.gt.wpts-1) then
        !
        lambda=lambda+(eliplus(ipts))/omega
        exit
        !
     endif
     !
     if ((i/2)*2==i) then
        !
        lambda=lambda+2*(eliplus(ipts))/omega
        !
     else
        !
        lambda=lambda+4*(eliplus(ipts))/omega
        !
     endif
     !
  enddo !ipts
  !
  lambda=(1.d0/3.d0)*lambda*d_w
  !
  return

end subroutine hlambda_calculate_full_t0
!----------------------------------------------------------------------------------------------
