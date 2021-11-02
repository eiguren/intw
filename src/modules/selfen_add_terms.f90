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
  use intw_input_parameters
  use intw_reading
  use w90_parameters, only: num_wann
  use intw_useful_constants, only: eps_5
  use intw_utility

  implicit none

  !I/O variables

  real(dp),intent(in) :: sigma_w,sigma_e
  real(dp),intent(in) :: ekk,ekq
  real(dp),intent(in) :: wqv
  complex(dp),intent(in) :: g_me(nspin,nspin)
  real(dp),intent(inout) :: termplus
  real(dp),intent(inout) :: termminus

  !local variables

  integer :: imode,ipts,is,js
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
  use intw_input_parameters
  use intw_reading
  use w90_parameters, only: num_wann
  use intw_useful_constants, only: eps_5
  use intw_utility

  implicit none

  !I/O variables

  real(dp),intent(in) :: sigma_w,sigma_e
  real(dp),intent(in) :: ekk,ekq
  real(dp),intent(in) :: wqv
  complex(dp),intent(in) :: g_me(nspin,nspin)
  real(dp),intent(inout) :: termminus

  !local variables

  integer :: imode,ipts,is,js
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
  use intw_input_parameters
  use intw_reading
  use w90_parameters, only: num_wann
  use intw_useful_constants, only: eps_5
  use intw_utility

  implicit none

  !I/O variables

  real(dp),intent(in) :: sigma_w,sigma_e
  real(dp),intent(in) :: ekk,ekq
  real(dp),intent(in) :: wqv
  complex(dp),intent(in) :: g_me(nspin,nspin)
  real(dp),intent(inout) :: termplus

  !local variables

  integer :: imode,ipts,is,js
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
  use intw_input_parameters
  use intw_reading
  use w90_parameters, only: num_wann
  use intw_useful_constants, only: eps_5
  use intw_utility

  implicit none

  !I/O variables

  real(dp),intent(in) :: sigma_w,sigma_e
  real(dp),intent(in) :: ekk,ekq
  real(dp),intent(in) :: wqv
  complex(dp),intent(in) :: g_me(nspin,nspin)
  real(dp),intent(inout) :: termplus
  real(dp),intent(inout) :: termminus

  !local variables

  integer :: imode,ipts,is,js
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
  use intw_input_parameters
  use intw_reading
  use w90_parameters, only: num_wann
  use intw_useful_constants, only: eps_5
  use intw_utility

  implicit none

  !I/O variables

  real(dp),intent(in) :: sigma_w,sigma_e
  real(dp),intent(in) :: ekk,ekq
  real(dp),intent(in) :: wqv
  complex(dp),intent(in) :: g_me(nspin,nspin)
  real(dp),intent(inout) :: termminus

  !local variables

  integer :: imode,ipts,is,js
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
  use intw_input_parameters
  use intw_reading
  use w90_parameters, only: num_wann
  use intw_useful_constants, only: eps_5
  use intw_utility

  implicit none

  !I/O variables

  real(dp),intent(in) :: sigma_w,sigma_e
  real(dp),intent(in) :: ekk,ekq
  real(dp),intent(in) :: wqv
  complex(dp),intent(in) :: g_me(nspin,nspin)
  real(dp),intent(inout) :: termplus

  !local variables

  integer :: imode,ipts,is,js
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
