!-----------------------------------------------------------------------------------------------------------------------------------
subroutine add_m_nu_lambda_text_no_spin(sigma_e,ekk,ekq,omega_phon,g_me,lambda_kq,nmode)
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 15/09/2017

  use kinds, only: dp
  use intw_reading, only: nat
  use intw_useful_constants, only: pi, eps_5, eps_6
  use intw_utility, only: weight_ph

  implicit none

  !I/O variables

  integer, intent(in) :: nmode
  real(dp), intent(in) :: sigma_e
  real(dp), intent(in) :: ekk,ekq
  real(dp), intent(in) :: omega_phon(3*nat)
  complex(dp), intent(in) :: g_me(3*nat)
  real(dp), intent(inout) :: lambda_kq

  !local variables

  integer :: imode
  real(dp) :: wqv, g2
  real(dp) :: delta_eplus, delta_eminus

  !$omp parallel default(none) &
  !$omp shared(nmode,ekk,ekq,sigma_e,eps_5,eps_6,g_me,lambda_kq,omega_phon,pi) &
  !$omp private(imode,wqv,delta_eplus,delta_eminus,g2)
  !
  !$omp do
  !
  do imode=1,nmode
     !
!     if (imode.eq.4.or.imode.eq.9.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24) cycle
     if (imode.eq.4.or.imode.eq.6.or.imode.eq.8.or.imode.eq.14.or.imode.eq.16.or.imode.eq.24.or.imode.eq.27 &
                   .or.imode.eq.28.or.imode.eq.30.or.imode.eq.34.or.imode.eq.35.or.imode.eq.36) cycle
     !
     wqv=abs(omega_phon(imode))
     !
     delta_eplus=0.d0
     delta_eminus=0.d0
     !
     if (abs(ekk-ekq+wqv).le.4.d0*sigma_e) delta_eplus=exp(-((ekk-ekq+wqv)/sigma_e)**2.d0)/(sigma_e*sqrt(pi))
     if (abs(ekk-ekq-wqv).le.4.d0*sigma_e) delta_eminus=exp(-((ekk-ekq-wqv)/sigma_e)**2.d0)/(sigma_e*sqrt(pi))
     !
     if (wqv.lt.eps_5) then
        !
     else
        !
        g2=weight_ph(wqv)*(abs(g_me(imode))**2.d0)/(2.d0*wqv)
        !
        lambda_kq=lambda_kq+(g2*(delta_eplus+delta_eminus)/wqv)
        !
     endif !wqv > 0
     !
  enddo !imode
  !
  !$omp end parallel
  !
  return

end subroutine add_m_nu_lambda_text_no_spin
!-------------------------------------------------------------------------------------------------
