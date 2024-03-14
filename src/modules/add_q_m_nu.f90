!-----------------------------------------------------------------------------------------------------------
subroutine add_q_m_nu_term_qe_so(w_min,w_max,wpts,sigma_w,delta_energ,omega_phon,g_me,eliash_f_so)
!-----------------------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================================!
! For a given state (k,n) we add to its respective eliashberg function the different nu terms               !
! associated to a (q,m) state we have selected before                                                       !
!                                                                                                           !
! (What we want at the end)                                                                                 !
! a^2F(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)]] !
!                                                                                                           !
! (What we do here)                                                                                         !
! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)]
!===========================================================================================================!

  use kinds, only: dp
  use intw_reading, only: nat, nspin
  use intw_utility, only: gaussian

  implicit none

  !I/O variables

  integer :: wpts
  real(dp), intent(in) :: w_min, w_max, sigma_w
  real(dp), intent(in) :: delta_energ
  real(dp), intent(in) :: omega_phon(3*nat)
  complex(dp), intent(in) :: g_me(nspin,nspin,3*nat)
  real(dp), intent(inout) :: eliash_f_so(wpts,nspin,nspin)

  !local variables

  integer :: imode, ipts, ispin, jspin
  real(dp) :: wqv, delta_omega(wpts), g2

  do imode=1,3*nat-4
     !
     if (imode.eq.4.or.imode.eq.9.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24) cycle
     !
     wqv=omega_phon(imode)
     call gaussian(wqv,w_min,w_max,wpts,sigma_w,delta_omega)
     !
     ! We finaly add the (q,m,nu) term to the eliashber function: a^2F(k,n,w)
     !
     do ispin=1,nspin
        do jspin=1,nspin
           !
           if (abs(wqv).lt.0.00001) then
              !
              g2=0.0d0
              !
           else
              !
              g2=(abs(g_me(ispin,jspin,imode))**2.d0)/(2.d0*abs(wqv))
              !
           endif
           !
           do ipts=1,wpts
              !
              eliash_f_so(ipts,ispin,jspin)=eliash_f_so(ipts,ispin,jspin)+g2*delta_omega(ipts)*delta_energ
              !
           enddo !ipts
           !
        enddo !jspin
     enddo !ispin
     !
  enddo !imode
  !
  return

end subroutine add_q_m_nu_term_qe_so
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
!---------------------------------------------------------------------------------------------
subroutine add_q_m_nu_term_qe(w_min,w_max,wpts,sigma_w,delta_energ,omega_phon,g_me,eliash_f)
!---------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================================!
! For a given state (k,n) we add to its respective eliashberg function the different nu terms               !
! associated to a (q,m) state we have selected before                                                       !
!                                                                                                           !
! (What we want at the end)                                                                                 !
! a^2F(k,n,si,sj,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)]] !
!                                                                                                           !
! (What we do here)                                                                                         !
! Sum_nu[|g(q,k,m,n,si,sj,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)]
!===========================================================================================================!

  use kinds, only: dp
  use intw_reading, only: nat, nspin
  use intw_utility, only: gaussian

  implicit none

  !I/O variables

  integer :: wpts
  real(dp), intent(in) :: w_min, w_max, sigma_w
  real(dp), intent(in) :: delta_energ
  real(dp), intent(in) :: omega_phon(3*nat)
  complex(dp), intent(in) :: g_me(nspin,nspin,3*nat)
  real(dp), intent(inout) :: eliash_f(wpts)

  !local variables

  integer :: imode, ipts
  real(dp) :: wqv, delta_omega(wpts), g2, g

  do imode=1,3*nat-4
!  do imode=1,3
     !
     if (imode.eq.4.or.imode.eq.9.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24) cycle
     !
     wqv=omega_phon(imode)
     call gaussian(wqv,w_min,w_max,wpts,sigma_w,delta_omega)
     !
     ! We finaly add the (q,m,nu) term to the eliashber function: a^2F(k,n,w)
     !
     if (abs(wqv).lt.0.00001) then
        !
        g2=0.0d0
        !
     else
        !
        if (nspin.gt.1) then
           !
           g=abs(g_me(1,1,imode))**2.d0+abs(g_me(1,2,imode))**2.d0 + &
             abs(g_me(2,1,imode))**2.d0+abs(g_me(2,2,imode))**2.d0
           g2=g/(2.d0*abs(wqv))
           !
        else
           !
           g2=(abs(g_me(1,1,imode))**2.d0)/(2.d0*abs(wqv))
           !
        endif
        !
     endif !wqv > 0
     !
     do ipts=1,wpts
        !
        eliash_f(ipts)=eliash_f(ipts)+g2*delta_omega(ipts)*delta_energ
        !
     enddo !ipts
     !
  enddo !imode
  !
  return

end subroutine add_q_m_nu_term_qe
!------------------------------------------------------------------------------------------------------------------------------------
!************************************************************************************************************************************
!------------------------------------------------------------------------------------------------------------------------------------
subroutine add_q_m_nu_term_full_so(w_min,w_max,wpts,sigma_w,sigma_e,ekk,ekq,omega_phon,g_me,eliash_fplus_so,eliash_fminus_so,nmode)
!------------------------------------------------------------------------------------------------------------------------------------
!
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
  use intw_reading, only: nat, nspin
  use intw_useful_constants, only: eps_5
  use intw_utility, only: weight_ph, gaussian

  implicit none

  !I/O variables

  integer,intent(in) :: wpts,nmode
  real(dp),intent(in) :: w_min,w_max,sigma_w,sigma_e
  real(dp),intent(in) :: ekk,ekq
  real(dp),intent(in) :: omega_phon(3*nat)
  complex(dp),intent(in) :: g_me(nspin,nspin,3*nat)
  real(dp),intent(inout) :: eliash_fplus_so(wpts,nspin,nspin)
  real(dp),intent(inout) :: eliash_fminus_so(wpts,nspin,nspin)

  !local variables

  integer :: imode,ipts,ispin,jspin
  real(dp) :: wqv,delta_omega(wpts),g2
  real(dp) :: delta_eplus(wpts),delta_eminus(wpts)

  do imode=1,nmode
     !
     !We discard modes related to the termination we are not interested in
     !
     if (imode.eq.4.or.imode.eq.9.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24) cycle
     !
     wqv=omega_phon(imode)
     call gaussian(wqv,w_min,w_max,wpts,sigma_w,delta_omega)
     call gaussian(-(ekk-ekq),w_min,w_max,wpts,sigma_e,delta_eplus)
     call gaussian((ekk-ekq),w_min,w_max,wpts,sigma_e,delta_eminus)
     !
     ! We finaly add the (q,m,nu) term to the eliashber function: a^2F+/-(k,n,w)
     !
     do ispin=1,nspin
        do jspin=1,nspin
           !
           if (wqv.lt.eps_5) then
              !
              g2=0.0d0
              !
           else
              !
              g2=weight_ph(wqv)*(abs(g_me(ispin,jspin,imode))**2.d0)/(2.d0*abs(wqv))
              !
           endif
           !
           do ipts=1,wpts
              !
              eliash_fplus_so(ipts,ispin,jspin)=eliash_fplus_so(ipts,ispin,jspin)+ &
                                          g2*delta_omega(ipts)*delta_eplus(ipts)
              !
              eliash_fminus_so(ipts,ispin,jspin)=eliash_fminus_so(ipts,ispin,jspin)+ &
                                           g2*delta_omega(ipts)*delta_eminus(ipts)
              !
           enddo !ipts
           !
        enddo !jspin
     enddo !ispin
     !
  enddo !imode
  !
  return

end subroutine add_q_m_nu_term_full_so
!---------------------------------------------------------------------------------------------------------------------
!*********************************************************************************************************************
!---------------------------------------------------------------------------------------------------------------------
subroutine add_q_m_nu_term_full(w_min,w_max,wpts,sigma_w,sigma_e,ekk,ekq,omega_phon,g_me,eliash_fplus,eliash_fminus,nmode)
!---------------------------------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 14/03/2016
!
!==================================================================================================================!
! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
! associated to a (q,m) state we have selected before                                                              !
!                                                                                                                  !
! (What we want at the end)                                                                                        !
! a^2F+/-(k,n,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
!                                                                                                                  !
! (What we do here)                                                                                                !
! Sum_nu[|g(q,k,m,n,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
!==================================================================================================================!

  use kinds, only: dp
  use intw_reading, only: nat, nspin
  use intw_useful_constants, only: eps_5
  use intw_utility, only: weight_ph, gaussian

  implicit none

  !I/O variables

  integer, intent(in) :: wpts,nmode
  real(dp), intent(in) :: w_min, w_max, sigma_w, sigma_e
  real(dp), intent(in) :: ekk,ekq
  real(dp), intent(in) :: omega_phon(3*nat)
  complex(dp), intent(in) :: g_me(nspin,nspin,3*nat)
  real(dp), intent(inout) :: eliash_fplus(wpts)
  real(dp), intent(inout) :: eliash_fminus(wpts)

  !local variables

  integer :: imode, ipts
  real(dp) :: wqv, delta_omega(wpts), g2
  real(dp) :: delta_eplus(wpts), delta_eminus(wpts)
  complex(dp) :: g

  do imode=1,nmode
     !
     !We discard modes related to the termination we are not interested in
     !
     if (imode.eq.4.or.imode.eq.9.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24) cycle
     !
     wqv=omega_phon(imode)
     call gaussian(wqv,w_min,w_max,wpts,sigma_w,delta_omega)
     call gaussian(-(ekk-ekq),w_min,w_max,wpts,sigma_e,delta_eplus)
     call gaussian((ekk-ekq),w_min,w_max,wpts,sigma_e,delta_eminus)
     !
     ! We finaly add the (q,m,nu) term to the eliashber function: a^2F+/-(k,n,w)
     !
     if (wqv.lt.eps_5) then
        !
        g2=0.0d0
        !
     else
        !
        if (nspin.gt.1) then
           !
           g=g_me(1,1,imode)+g_me(1,2,imode)+g_me(2,1,imode)+g_me(2,2,imode)
           !
           g2=weight_ph(wqv)*(abs(g)**2.d0)/(2.d0*abs(wqv))
           !
        else
           !
           g2=weight_ph(wqv)*(abs(g_me(1,1,imode))**2.d0)/(2.d0*abs(wqv))
           !
        endif
        !
     endif !wqv > 0
     !
     do ipts=1,wpts
        !
        eliash_fplus(ipts)=eliash_fplus(ipts)+ &
                           g2*delta_omega(ipts)*delta_eplus(ipts)
        !
        eliash_fminus(ipts)=eliash_fminus(ipts)+ &
                            g2*delta_omega(ipts)*delta_eminus(ipts)
        !
     enddo !ipts
     !
  enddo !imode
  !
  return

end subroutine add_q_m_nu_term_full
!-----------------------------------------------------------------------------------------------------------------------------------
!***********************************************************************************************************************************
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine add_q_m_nu_term_full_no_spin(w_min,w_max,wpts,sigma_w,sigma_e,ekk,ekq,omega_phon,g_me,eliash_fplus,eliash_fminus,nmode)
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 14/03/2016
!
!==================================================================================================================!
! For a given state (k,n) we add to its respective eliashberg function the different nu terms                      !
! associated to a (q,m) state we have selected before                                                              !
!                                                                                                                  !
! (What we want at the end)                                                                                        !
! a^2F+/-(k,n,w)=(1/nqfmesh)*Sum_q_m_nu[|g(q,k,m,n,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]] !
!                                                                                                                  !
! (What we do here)                                                                                                !
! Sum_nu[|g(q,k,m,n,nu|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)+/-w]                                        !
!==================================================================================================================!

  use kinds, only: dp
  use intw_reading, only: nat
  use intw_useful_constants, only: eps_5, eps_6
  use intw_utility, only: weight_ph, gaussian

  implicit none

  !I/O variables

  integer, intent(in) :: wpts,nmode
  real(dp), intent(in) :: w_min,w_max,sigma_w,sigma_e
  real(dp), intent(in) :: ekk,ekq
  real(dp), intent(in) :: omega_phon(3*nat)
  complex(dp), intent(in) :: g_me(3*nat)
  real(dp), intent(inout) :: eliash_fplus(wpts)
  real(dp) ,intent(inout) :: eliash_fminus(wpts)

  !local variables

  integer :: imode, ipts
  real(dp) :: wqv, delta_omega(wpts), g2
  real(dp) :: delta_eplus(wpts), delta_eminus(wpts)

  !$omp parallel default(none) &
  !$omp shared(nmode,w_min,w_max,wpts,sigma_w,ekk,ekq,sigma_e,eps_5,eps_6,g_me,eliash_fplus,eliash_fminus,omega_phon) &
  !$omp private(imode,wqv,delta_omega,delta_eplus,delta_eminus,g2)
  !
  !$omp do
  !
  do imode=1,nmode
     !
     !We discard modes related to the termination we are not interested in
     !
!     if (imode.eq.4.or.imode.eq.9.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24) cycle
     if (imode.eq.4.or.imode.eq.6.or.imode.eq.8.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24.or.imode.eq.27 &
                   .or.imode.eq.28.or.imode.eq.30.or.imode.eq.34.or.imode.eq.35.or.imode.eq.36) cycle
     !
     wqv=abs(omega_phon(imode))
     call gaussian(wqv,w_min,w_max,wpts,sigma_w,delta_omega)
     call gaussian((ekq-ekk),w_min,w_max,wpts,sigma_e,delta_eplus)
     call gaussian(-(ekq-ekk),w_min,w_max,wpts,sigma_e,delta_eminus)
     !
     ! We finaly add the (q,m,nu) term to the eliashber function: a^2F+/-(k,n,w)
     !
     if (wqv.lt.eps_5) then
        !
        g2=0.0d0
        !
     else
        !
        g2=weight_ph(wqv)*(abs(g_me(imode))**2.d0)/(2.d0*wqv)
        !
     endif !wqv > 0
     !
     do ipts=1,wpts
        !
        eliash_fplus(ipts)=eliash_fplus(ipts)+ &
                           g2*delta_omega(ipts)*delta_eplus(ipts)
        !
        eliash_fminus(ipts)=eliash_fminus(ipts)+ &
                            g2*delta_omega(ipts)*delta_eminus(ipts)
        !
     enddo !ipts
     !
  enddo !imode
  !
  !$omp end parallel
  !
  return

end subroutine add_q_m_nu_term_full_no_spin
!-------------------------------------------------------------------------------------------------
