module intw_selfen

  implicit none

  public :: realselfe_eliash_qe_t0, &
            imagselfe_eliash_qe_t0, &
            realselfe_eliash_full_t0, &
            imagselfe_eliash_full_t0, &
            imselfen_calculate_full_t0, &
            reselfen_calculate_full_t0, &
            add_q_m_nu_imag, &
            add_q_m_nu_imag_minus, &
            add_q_m_nu_imag_plus, &
            add_q_m_nu_real, &
            add_q_m_nu_real_minus, &
            add_q_m_nu_real_plus, &
            imselfen2reselfen_t0

  private


contains

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

    use kinds, only: dp

    implicit none

    !I/O variables

    integer,intent(in) :: wpts,epts
    real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
    real(dp),intent(in) :: eliash_f(wpts)
    real(dp),intent(inout) :: realse(epts)

    !local variables

    integer :: j,ipts,jpts
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

    use kinds, only: dp
    use intw_useful_constants, only: pi

    implicit none

    !I/O variables

    integer,intent(in) :: wpts,epts
    real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
    real(dp),intent(in) :: eliash_f(wpts)
    real(dp),intent(inout) :: imagse(epts)

    !local variables

    integer :: j,ipts,jpts
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

    use kinds, only: dp

    implicit none

    !I/O variables

    integer,intent(in) :: wpts,epts
    real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
    real(dp),intent(in) :: eliash_fplus(wpts),eliash_fminus(wpts)
    real(dp),intent(inout) :: realse(epts)

    !local variables

    integer :: j,ipts,jpts
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
  !           relminus=0.5d0*log((abs(eps-omega)**2+sigma_e**2)/(abs(omega)**2+sigma_e**2))
  !           relplus=0.0d0
            !
          else
            !
  !           relminus=0.0d0
            relplus=0.5d0*log((abs(eps+omega)**2+sigma_e**2)/(abs(eps-omega)**2+sigma_e**2))
  !           relplus=0.5d0*log((abs(eps+omega)**2+sigma_e**2)/(abs(omega)**2+sigma_e**2))
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

    use kinds, only: dp
    use intw_useful_constants, only: pi

    implicit none

    !I/O variables

    integer,intent(in) :: wpts,epts
    real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e
    real(dp),intent(in) :: eliash_fplus(wpts),eliash_fminus(wpts)
    real(dp),intent(inout) :: imagse(epts)

    !local variables

    integer :: j,ipts,jpts
    real(dp) :: omega,eps,sigma_e,imlplus,imlminus

    imagse(:)=0.0d0
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
  !---------------------------------------------------------------------------------------------------------------------
  !*********************************************************************************************************************
  !---------------------------------------------------------------------------------------------------------------------
  subroutine imselfen_calculate_full_t0(wpts,w_min,w_max,d_w,epts,e_min,e_max,d_e,eliash_fplus,eliash_fminus,imagse,ef)
  !---------------------------------------------------------------------------------------------------------------------
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

    use kinds, only: dp
    use intw_useful_constants, only: pi

    implicit none

    !I/O variables

    integer,intent(in) :: wpts,epts
    real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e,ef
    real(dp),intent(in) :: eliash_fplus(wpts),eliash_fminus(wpts)
    real(dp),intent(inout) :: imagse(epts)

    !local variables

    integer :: j,ipts,jpts
    real(dp) :: omega,eps,imlplus,imlminus

    imagse(:)=0.0d0
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
  !        if (eps+ef-omega.gt.0.d0) then
          if (eps+ef.gt.omega) then
            imlminus=imlminus+pi
          elseif (eps+ef.lt.-1.d0*omega) then
  !        elseif (eps+ef+omega.lt.0.d0) then
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

  end subroutine imselfen_calculate_full_t0
  !----------------------------------------------------------------------------------------------------------------------
  !**********************************************************************************************************************
  !----------------------------------------------------------------------------------------------------------------------
  subroutine reselfen_calculate_full_t0(wpts,w_min,w_max,d_w,epts,e_min,e_max,d_e,eliash_fplus,eliash_fminus,realse,ef)
  !----------------------------------------------------------------------------------------------------------------------
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

    use kinds, only: dp

    implicit none

    !I/O variables

    integer,intent(in) :: wpts,epts
    real(dp),intent(in) :: w_min,w_max,d_w,e_min,e_max,d_e,ef
    real(dp),intent(in) :: eliash_fplus(wpts),eliash_fminus(wpts)
    real(dp),intent(inout) :: realse(epts)

    !local variables

    integer :: j,ipts,jpts
    real(dp) :: omega,eps,sigma_e,relminus,relplus

    realse(:)=0.0d0
    !
  !  sigma_e=15.0d0*d_e
    sigma_e=0.050d0
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
          if (eps+ef.gt.omega) then
            !
            relminus=0.5d0*log((abs(eps+ef-omega)**2+sigma_e**2)/(abs(eps+ef+omega)**2+sigma_e**2))
            !
          elseif (eps+ef.lt.-1.d0*omega) then
            !
            relplus=-0.5d0*log((abs(eps+ef+omega)**2+sigma_e**2)/(abs(eps+ef-omega)**2+sigma_e**2))
            !
          endif
          !
          if (jpts.gt.wpts-1) then
            !
            realse(ipts)=realse(ipts)+&
            (eliash_fminus(jpts)*relminus+eliash_fplus(jpts)*relplus)
            exit
            !
          endif
          !
          if ((j/2)*2==j) then
            !
            realse(ipts)=realse(ipts)+&
            2*(eliash_fminus(jpts)*relminus+eliash_fplus(jpts)*relplus)
            !
          else
            !
            realse(ipts)=realse(ipts)+&
            4*(eliash_fminus(jpts)*relminus+eliash_fplus(jpts)*relplus)
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

  end subroutine reselfen_calculate_full_t0
  !------------------------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------------------------
  subroutine add_q_m_nu_imag(sigma_w,sigma_e,ekk,ekq,wqv,g_me,termplus,termminus)
  !---------------------------------------------------------------------------------------------------------------------
  !
  ! Description not related with the subroutine
  ! Created by Peio G. Goiricelaya 14/03/2016
  !
  !==================================================================================================================!
  ! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
  ! associated to a (q,m) state we have selected before                                                              !
  !                                                                                                                  !
  ! (What we want at the end)                                                                                        !
  ! a^2F+/-(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
  !                                                                                                                  !
  ! (What we do here)                                                                                                !
  ! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
  !==================================================================================================================!

    use kinds, only: dp
    use intw_reading, only: nspin
    use intw_useful_constants, only: eps_5, pi
    use intw_utility, only: weight_ph

    implicit none

    !I/O variables

    real(dp),intent(in) :: sigma_w,sigma_e
    real(dp),intent(in) :: ekk,ekq
    real(dp),intent(in) :: wqv
    complex(dp),intent(in) :: g_me(nspin,nspin)
    real(dp),intent(inout) :: termplus
    real(dp),intent(inout) :: termminus

    !local variables

    real(dp) :: g2
    real(dp) :: delta_eplus,delta_eminus
    complex(dp) :: g

    if (wqv.lt.eps_5) then
      !
      g2=0.0d0
      !
    else
      !
      if (nspin.gt.1) then
          !
          g=g_me(1,1)+g_me(1,2)+g_me(2,1)+g_me(2,2)
          !
          g2=weight_ph(wqv)*(abs(g)**2.d0)/(2.d0*abs(wqv))
          !
      else
          !
          g2=weight_ph(wqv)*(abs(g_me(1,1))**2.d0)/(2.d0*abs(wqv))
          !
      endif
      !
    endif !wqv > 0
    !
    if (abs(ekk-ekq+wqv).gt.(4.d0*sigma_e)) then
      delta_eplus=0.d0
    else
      delta_eplus=exp(-((ekk-ekq+wqv)**2.d0)/(sigma_e**2.d0))/(sigma_e*sqrt(pi))
    endif
    !
    if (abs(ekk-ekq-wqv).gt.(4.d0*sigma_e)) then
      delta_eminus=0.d0
    else
      delta_eminus=exp(-((ekk-ekq-wqv)**2.d0)/(sigma_e**2.d0))/(sigma_e*sqrt(pi))
    endif
    !
    termplus=g2*delta_eplus
    termminus=g2*delta_eminus
    !
    return

  end subroutine add_q_m_nu_imag
  !-------------------------------------------------------------------------------------------------
  !*************************************************************************************************
  !-------------------------------------------------------------------------------------------------
  subroutine add_q_m_nu_imag_minus(sigma_w,sigma_e,ekk,ekq,wqv,g_me,termminus)
  !---------------------------------------------------------------------------------------------------------------------
  !
  ! Description not related with the subroutine
  ! Created by Peio G. Goiricelaya 14/03/2016
  !
  !==================================================================================================================!
  ! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
  ! associated to a (q,m) state we have selected before                                                              !
  !                                                                                                                  !
  ! (What we want at the end)                                                                                        !
  ! a^2F+/-(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
  !                                                                                                                  !
  ! (What we do here)                                                                                                !
  ! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
  !==================================================================================================================!

    use kinds, only: dp
    use intw_reading, only: nspin
    use intw_useful_constants, only: eps_5, pi
    use intw_utility, only: weight_ph

    implicit none

    !I/O variables

    real(dp),intent(in) :: sigma_w,sigma_e
    real(dp),intent(in) :: ekk,ekq
    real(dp),intent(in) :: wqv
    complex(dp),intent(in) :: g_me(nspin,nspin)
    real(dp),intent(inout) :: termminus

    !local variables

    real(dp) :: g2
    real(dp) :: delta_eminus
    complex(dp) :: g

    if (wqv.lt.eps_5) then
      !
      g2=0.0d0
      !
    else
      !
      if (nspin.gt.1) then
          !
          g=g_me(1,1)+g_me(1,2)+g_me(2,1)+g_me(2,2)
          !
          g2=weight_ph(wqv)*(abs(g)**2.d0)/(2.d0*abs(wqv))
          !
      else
          !
          g2=weight_ph(wqv)*(abs(g_me(1,1))**2.d0)/(2.d0*abs(wqv))
          !
      endif
      !
    endif !wqv > 0
    !
    if (abs(ekk-ekq-wqv).gt.(4.d0*sigma_e)) then
      delta_eminus=0.d0
    else
      delta_eminus=exp(-((ekk-ekq-wqv)**2.d0)/(sigma_e**2.d0))/(sigma_e*sqrt(pi))
    endif
    !
    termminus=g2*delta_eminus
    !
    return

  end subroutine add_q_m_nu_imag_minus
  !------------------------------------------------------------
  !************************************************************
  !---------------------------------------------------------------------------------------------------------------------
  subroutine add_q_m_nu_imag_plus(sigma_w,sigma_e,ekk,ekq,wqv,g_me,termplus)
  !---------------------------------------------------------------------------------------------------------------------
  !
  ! Description not related with the subroutine
  ! Created by Peio G. Goiricelaya 14/03/2016
  !
  !==================================================================================================================!
  ! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
  ! associated to a (q,m) state we have selected before                                                              !
  !                                                                                                                  !
  ! (What we want at the end)                                                                                        !
  ! a^2F+/-(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
  !                                                                                                                  !
  ! (What we do here)                                                                                                !
  ! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
  !==================================================================================================================!

    use kinds, only: dp
    use intw_reading, only: nspin
    use intw_useful_constants, only: eps_5, pi
    use intw_utility, only: weight_ph

    implicit none

    !I/O variables

    real(dp),intent(in) :: sigma_w,sigma_e
    real(dp),intent(in) :: ekk,ekq
    real(dp),intent(in) :: wqv
    complex(dp),intent(in) :: g_me(nspin,nspin)
    real(dp),intent(inout) :: termplus

    !local variables

    real(dp) :: g2
    real(dp) :: delta_eplus
    complex(dp) :: g

    if (wqv.lt.eps_5) then
      !
      g2=0.0d0
      !
    else
      !
      if (nspin.gt.1) then
          !
          g=g_me(1,1)+g_me(1,2)+g_me(2,1)+g_me(2,2)
          !
          g2=weight_ph(wqv)*(abs(g)**2.d0)/(2.d0*abs(wqv))
          !
      else
          !
          g2=weight_ph(wqv)*(abs(g_me(1,1))**2.d0)/(2.d0*abs(wqv))
          !
      endif
      !
    endif !wqv > 0
    !
    if (abs(ekk-ekq+wqv).gt.(4.d0*sigma_e)) then
      delta_eplus=0.d0
    else
      delta_eplus=exp(-((ekk-ekq+wqv)**2.d0)/(sigma_e**2.d0))/(sigma_e*sqrt(pi))
    endif
    !
    termplus=g2*delta_eplus
    !
    return

  end subroutine add_q_m_nu_imag_plus
  !----------------------------------------------------
  !****************************************************
  !----------------------------------------------------
  subroutine add_q_m_nu_real(sigma_w,sigma_e,ekk,ekq,wqv,g_me,termplus,termminus)
  !---------------------------------------------------------------------------------------------------------------------
  !
  ! Description not related with the subroutine
  ! Created by Peio G. Goiricelaya 14/03/2016
  !
  !==================================================================================================================!
  ! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
  ! associated to a (q,m) state we have selected before                                                              !
  !                                                                                                                  !
  ! (What we want at the end)                                                                                        !
  ! a^2F+/-(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
  !                                                                                                                  !
  ! (What we do here)                                                                                                !
  ! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
  !==================================================================================================================!

    use kinds, only: dp
    use intw_reading, only: nspin
    use intw_useful_constants, only: eps_5
    use intw_utility, only: weight_ph

    implicit none

    !I/O variables

    real(dp),intent(in) :: sigma_w,sigma_e
    real(dp),intent(in) :: ekk,ekq
    real(dp),intent(in) :: wqv
    complex(dp),intent(in) :: g_me(nspin,nspin)
    real(dp),intent(inout) :: termplus
    real(dp),intent(inout) :: termminus

    !local variables

    real(dp) :: g2
    real(dp) :: term_eplus,term_eminus
    complex(dp) :: g

    if (wqv.lt.eps_5) then
      !
      g2=0.0d0
      !
    else
      !
      if (nspin.gt.1) then
          !
          g=g_me(1,1)+g_me(1,2)+g_me(2,1)+g_me(2,2)
          !
          g2=weight_ph(wqv)*(abs(g)**2.d0)/(2.d0*abs(wqv))
          !
      else
          !
          g2=weight_ph(wqv)*(abs(g_me(1,1))**2.d0)/(2.d0*abs(wqv))
          !
      endif
      !
    endif !wqv > 0
    !
    term_eplus=(ekk-ekq+wqv)/((ekk-ekq+wqv)**2.d0+sigma_e**2.d0)
    !
    term_eminus=(ekk-ekq-wqv)/((ekk-ekq-wqv)**2.d0+sigma_e**2.d0)
    !
    termplus=g2*term_eplus
    termminus=g2*term_eminus
    !
    return

  end subroutine add_q_m_nu_real
  !-------------------------------------------------------------------------------------------------
  !*************************************************************************************************
  !-------------------------------------------------------------------------------------------------
  subroutine add_q_m_nu_real_minus(sigma_w,sigma_e,ekk,ekq,wqv,g_me,termminus)
  !---------------------------------------------------------------------------------------------------------------------
  !
  ! Description not related with the subroutine
  ! Created by Peio G. Goiricelaya 14/03/2016
  !
  !==================================================================================================================!
  ! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
  ! associated to a (q,m) state we have selected before                                                              !
  !                                                                                                                  !
  ! (What we want at the end)                                                                                        !
  ! a^2F+/-(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
  !                                                                                                                  !
  ! (What we do here)                                                                                                !
  ! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
  !==================================================================================================================!

    use kinds, only: dp
    use intw_reading, only: nspin
    use intw_useful_constants, only: eps_5
    use intw_utility, only: weight_ph

    implicit none

    !I/O variables

    real(dp),intent(in) :: sigma_w,sigma_e
    real(dp),intent(in) :: ekk,ekq
    real(dp),intent(in) :: wqv
    complex(dp),intent(in) :: g_me(nspin,nspin)
    real(dp),intent(inout) :: termminus

    !local variables

    real(dp) :: g2
    real(dp) :: term_eminus
    complex(dp) :: g

    if (wqv.lt.eps_5) then
      !
      g2=0.0d0
      !
    else
      !
      if (nspin.gt.1) then
          !
          g=g_me(1,1)+g_me(1,2)+g_me(2,1)+g_me(2,2)
          !
          g2=weight_ph(wqv)*(abs(g)**2.d0)/(2.d0*abs(wqv))
          !
      else
          !
          g2=weight_ph(wqv)*(abs(g_me(1,1))**2.d0)/(2.d0*abs(wqv))
          !
      endif
      !
    endif !wqv > 0
    !
    term_eminus=(ekk-ekq-wqv)/((ekk-ekq-wqv)**2.d0+sigma_e**2.d0)
    !
    termminus=g2*term_eminus
    !
    return

  end subroutine add_q_m_nu_real_minus
  !------------------------------------------------------------
  !************************************************************
  !---------------------------------------------------------------------------------------------------------------------
  subroutine add_q_m_nu_real_plus(sigma_w,sigma_e,ekk,ekq,wqv,g_me,termplus)
  !---------------------------------------------------------------------------------------------------------------------
  !
  ! Description not related with the subroutine
  ! Created by Peio G. Goiricelaya 14/03/2016
  !
  !==================================================================================================================!
  ! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
  ! associated to a (q,m) state we have selected before                                                              !
  !                                                                                                                  !
  ! (What we want at the end)                                                                                        !
  ! a^2F+/-(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
  !                                                                                                                  !
  ! (What we do here)                                                                                                !
  ! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
  !==================================================================================================================!

    use kinds, only: dp
    use intw_reading, only: nspin
    use intw_useful_constants, only: eps_5
    use intw_utility, only: weight_ph

    implicit none

    !I/O variables

    real(dp),intent(in) :: sigma_w,sigma_e
    real(dp),intent(in) :: ekk,ekq
    real(dp),intent(in) :: wqv
    complex(dp),intent(in) :: g_me(nspin,nspin)
    real(dp),intent(inout) :: termplus

    !local variables

    real(dp) :: g2
    real(dp) :: term_eplus
    complex(dp) :: g

    if (wqv.lt.eps_5) then
      !
      g2=0.0d0
      !
    else
      !
      if (nspin.gt.1) then
          !
          g=g_me(1,1)+g_me(1,2)+g_me(2,1)+g_me(2,2)
          !
          g2=weight_ph(wqv)*(abs(g)**2.d0)/(2.d0*abs(wqv))
          !
      else
          !
          g2=weight_ph(wqv)*(abs(g_me(1,1))**2.d0)/(2.d0*abs(wqv))
          !
      endif
      !
    endif !wqv > 0
    !
    term_eplus=(ekk-ekq+wqv)/((ekk-ekq+wqv)**2.d0+sigma_e**2.d0)
    !
    termplus=g2*term_eplus
    !
    return

  end subroutine add_q_m_nu_real_plus
  !---------------------------------------------------------------------------------------------
  !*********************************************************************************************
  !---------------------------------------------------------------------------------------------
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
    use kinds, only: dp
    use intw_useful_constants, only: pi, cmplx_0

    implicit none

    !I/O variables

    integer,intent(in) :: epts
    real(dp),intent(in) :: e_min,e_max,d_e
    real(dp),intent(in) :: imselfen_w(epts)
    real(dp),intent(inout) :: reselfen_w(epts)

    !local variables

    real(dp) :: dimselfen_w(epts),integral,tita,modz
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

end module intw_selfen