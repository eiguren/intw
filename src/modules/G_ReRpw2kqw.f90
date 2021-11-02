!----------------------------------------------------------------------------------
subroutine G_ReRpw2kRpw(nrr_k,nrr_q,irr_k,ndegen_k,kpoint,g_ReRp,g_Rp)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(Rp,Re,...) from Real Wannier representation
! (Re) to Reciprocal Wannier representation associated to electrons (k), so that:
!
!  g(Rp,k,...) = Sum_in_Re_WS[g(Rp,Re,...)*Exp[-i*k*Re]]/(NRe)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k,nrr_q
  integer,intent(in) :: irr_k(3,nrr_k),ndegen_k(nrr_k)
  real(dp),intent(in) :: kpoint(3)
  complex(dp),intent(in) :: g_ReRp(nrr_q,nrr_k,num_wann,num_wann,npol,npol,3*nat)
  complex(dp),intent(inout) :: g_Rp(nrr_q,num_wann,num_wann,npol,npol,3*nat)
  
  !local variables

  integer :: irq,irk,imode,ipol,jpol,ibnd,jbnd
  real(dp) :: arg
  complex(dp) :: fac

  do irq=1,nrr_q
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              do ipol=1,npol
                 do jpol=1,npol
                    !
                    do irk=1,nrr_k
                       arg=tpi*dot_product(kpoint(:),dble(irr_k(:,irk)))
                       fac=exp(cmplx_i*arg)/dble(ndegen_k(irk))
                       g_Rp(irq,ibnd,jbnd,ipol,jpol,imode)=g_Rp(irq,ibnd,jbnd,ipol,jpol,imode)+ &
                                                         fac*g_ReRp(irq,irk,ibnd,jbnd,ipol,jpol,imode)
                    enddo !irk
                    !
                 enddo !jpol
              enddo !ipol
           enddo !jbnd
        enddo !ibnd
     enddo !imode
  enddo !irq
  !
  return

end subroutine G_ReRpw2kRpw
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_ReRpw2kRpw_no_spin(nrr_k,nrr_q,irr_k,ndegen_k,kpoint,g_ReRp,g_Rp)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(Rp,Re,...) from Real Wannier representation
! (Re) to Reciprocal Wannier representation associated to electrons (k), so that:
!
!  g(Rp,k,...) = Sum_in_Re_WS[g(Rp,Re,...)*Exp[-i*k*Re]]/(NRe)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k,nrr_q
  integer,intent(in) :: irr_k(3,nrr_k),ndegen_k(nrr_k)
  real(dp),intent(in) :: kpoint(3)
  complex(dp),intent(in) :: g_ReRp(nrr_q,nrr_k,num_wann,num_wann,3*nat)
  complex(dp),intent(inout) :: g_Rp(nrr_q,num_wann,num_wann,3*nat)
  
  !local variables

  integer :: irq,irk,imode,ibnd,jbnd
  real(dp) :: arg
  complex(dp) :: fac

  do irq=1,nrr_q
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              !
              do irk=1,nrr_k
                 arg=tpi*dot_product(kpoint(:),dble(irr_k(:,irk)))
                 fac=exp(cmplx_i*arg)/dble(ndegen_k(irk))
                 g_Rp(irq,ibnd,jbnd,imode)=g_Rp(irq,ibnd,jbnd,imode)+ &
                                           fac*g_ReRp(irq,irk,ibnd,jbnd,imode)
              enddo !irk
              !
           enddo !jbnd
        enddo !ibnd
     enddo !imode
  enddo !irq
  !
  return

end subroutine G_ReRpw2kRpw_no_spin
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_ReRpw2Reqw(nrr_k,nrr_q,irr_q,ndegen_q,qpoint,g_ReRp,g_Re)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(Rp,Re,...) from Real Wannier representation
! (Rp) to Reciprocal Wannier representation associated to phonons (q), so that:
!
!  g(q,Re,...) = Sum_in_Rp_WS[g(Rp,Re,...)*Exp[-i*q*Rp]]/(NRp)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k,nrr_q
  integer,intent(in) :: irr_q(3,nrr_q),ndegen_q(nrr_q)
  real(dp),intent(in) :: qpoint(3)
  complex(dp),intent(in) :: g_ReRp(nrr_q,nrr_k,num_wann,num_wann,npol,npol,3*nat)
  complex(dp),intent(inout) :: g_Re(nrr_k,num_wann,num_wann,npol,npol,3*nat)
  
  !local variables

  integer :: irq,irk,imode,ipol,jpol,ibnd,jbnd
  real(dp) :: arg
  complex(dp) :: fac

  do irk=1,nrr_k
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              do ipol=1,npol
                 do jpol=1,npol
                    !
                    do irq=1,nrr_q
                       arg=tpi*dot_product(qpoint(:),dble(irr_q(:,irq)))
                       fac=exp(cmplx_i*arg)/dble(ndegen_q(irq))
                       g_Re(irk,ibnd,jbnd,ipol,jpol,imode)=g_Re(irk,ibnd,jbnd,ipol,jpol,imode)+ &
                                                         fac*g_ReRp(irq,irk,ibnd,jbnd,ipol,jpol,imode)
                    enddo !irq
                    !
                 enddo !jpol
              enddo !ipol
           enddo !jbnd
        enddo !ibnd
     enddo !imode
  enddo !irk
  !
  return

end subroutine G_ReRpw2Reqw
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_ReRpw2Reqw_no_spin(nrr_k,nrr_q,irr_q,ndegen_q,qpoint,g_ReRp,g_Re)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(Rp,Re,...) from Real Wannier representation
! (Rp) to Reciprocal Wannier representation associated to phonons (q), so that:
!
!  g(q,Re,...) = Sum_in_Rp_WS[g(Rp,Re,...)*Exp[-i*q*Rp]]/(NRp)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k,nrr_q
  integer,intent(in) :: irr_q(3,nrr_q),ndegen_q(nrr_q)
  real(dp),intent(in) :: qpoint(3)
  complex(dp),intent(in) :: g_ReRp(nrr_q,nrr_k,num_wann,num_wann,3*nat)
  complex(dp),intent(inout) :: g_Re(nrr_k,num_wann,num_wann,3*nat)
  
  !local variables

  integer :: irq,irk,imode,ipol,jpol,ibnd,jbnd
  real(dp) :: arg
  complex(dp) :: fac

  do irk=1,nrr_k
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              !
              do irq=1,nrr_q
                 arg=tpi*dot_product(qpoint(:),dble(irr_q(:,irq)))
                 fac=exp(cmplx_i*arg)/dble(ndegen_q(irq))
                 g_Re(irk,ibnd,jbnd,imode)=g_Re(irk,ibnd,jbnd,imode)+ &
                                                   fac*g_ReRp(irq,irk,ibnd,jbnd,imode)
              enddo !irq
              !
           enddo !jbnd
        enddo !ibnd
     enddo !imode
  enddo !irk
  !
  return

end subroutine G_ReRpw2Reqw_no_spin
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_kRpw2kqw(nrr_q,irr_q,ndegen_q,qpoint,g_Rp,g_wannier)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(Rp,k,...) from Real Wannier representation
! (Rp) to Reciprocal Wannier representation associated to phonons (q), so that:
!
!  g(q,k,...) = Sum_in_Rp_WS[g(Rp,k,...)*Exp[-i*q*Rp]]/(NRp)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_q
  integer,intent(in) :: irr_q(3,nrr_q),ndegen_q(nrr_q)
  real(dp),intent(in) :: qpoint(3)
  complex(dp),intent(in) :: g_Rp(nrr_q,num_wann,num_wann,npol,npol,3*nat)
  complex(dp),intent(inout) :: g_wannier(num_wann,num_wann,npol,npol,3*nat)

  !local variables

  integer :: irq,imode,ibnd,jbnd,ipol,jpol
  real(dp) :: arg
  complex(dp) :: fac

  do imode=1,3*nat
     do ibnd=1,num_wann
        do jbnd=1,num_wann
           do ipol=1,npol
              do jpol=1,npol
                 !
                 do irq=1,nrr_q
                    arg=tpi*dot_product(qpoint(:),dble(irr_q(:,irq)))
                    fac=exp(cmplx_i*arg)/dble(ndegen_q(irq))
                    g_wannier(ibnd,jbnd,ipol,jpol,imode)=g_wannier(ibnd,jbnd,ipol,jpol,imode)+ &
                                                     fac*g_Rp(irq,ibnd,jbnd,ipol,jpol,imode)
                 enddo !irq
                 !
              enddo !jpol
           enddo !ipol
        enddo !jbnd
     enddo !ibnd
  enddo !imode
  !
  return

end subroutine G_kRpw2kqw
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_kRpw2kqw_no_spin(nrr_q,irr_q,ndegen_q,qpoint,g_Rp,g_wannier)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(Rp,k,...) from Real Wannier representation
! (Rp) to Reciprocal Wannier representation associated to phonons (q), so that:
!
!  g(q,k,...) = Sum_in_Rp_WS[g(Rp,k,...)*Exp[-i*q*Rp]]/(NRp)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_q
  integer,intent(in) :: irr_q(3,nrr_q),ndegen_q(nrr_q)
  real(dp),intent(in) :: qpoint(3)
  complex(dp),intent(in) :: g_Rp(nrr_q,num_wann,num_wann,3*nat)
  complex(dp),intent(inout) :: g_wannier(num_wann,num_wann,3*nat)

  !local variables

  integer :: irq,imode,ibnd,jbnd
  real(dp) :: arg
  complex(dp) :: fac

  do imode=1,3*nat
     do ibnd=1,num_wann
        do jbnd=1,num_wann
           !
           do irq=1,nrr_q
              arg=tpi*dot_product(qpoint(:),dble(irr_q(:,irq)))
              fac=exp(cmplx_i*arg)/dble(ndegen_q(irq))
              g_wannier(ibnd,jbnd,imode)=g_wannier(ibnd,jbnd,imode)+ &
                                         fac*g_Rp(irq,ibnd,jbnd,imode)
           enddo !irq
           !
        enddo !jbnd
     enddo !ibnd
  enddo !imode
  !
  return

end subroutine G_kRpw2kqw_no_spin
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_Reqw2kqw(nrr_k,irr_k,ndegen_k,kpoint,g_Re,g_wannier)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(q,Re,...) from Real Wannier representation
! (Re) to Reciprocal Wannier representation associated to electrons (k), so that:
!
!  g(q,k,...) = Sum_in_Re_WS[g(q,Re,...)*Exp[-i*q*Re]]/(NRe)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k
  integer,intent(in) :: irr_k(3,nrr_k),ndegen_k(nrr_k)
  real(dp),intent(in) :: kpoint(3)
  complex(dp),intent(in) :: g_Re(nrr_k,num_wann,num_wann,npol,npol,3*nat)
  complex(dp),intent(inout) :: g_wannier(num_wann,num_wann,npol,npol,3*nat)

  !local variables

  integer :: irk,imode,ibnd,jbnd,ipol,jpol
  real(dp) :: arg
  complex(dp) :: fac

  do imode=1,3*nat
     do ibnd=1,num_wann
        do jbnd=1,num_wann
           do ipol=1,npol
              do jpol=1,npol
                 !
                 do irk=1,nrr_k
                    arg=tpi*dot_product(kpoint(:),dble(irr_k(:,irk)))
                    fac=exp(cmplx_i*arg)/dble(ndegen_k(irk))
                    g_wannier(ibnd,jbnd,ipol,jpol,imode)=g_wannier(ibnd,jbnd,ipol,jpol,imode)+ &
                                                      fac*g_Re(irk,ibnd,jbnd,ipol,jpol,imode)
                 enddo !irk
                 !
              enddo !jpol
           enddo !ipol
        enddo !jbnd
     enddo !ibnd
  enddo !imode
  !
  return

end subroutine G_Reqw2kqw
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_Reqw2kqw_no_spin(nrr_k,irr_k,ndegen_k,kpoint,g_Re,g_wannier)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 11/03/2016
!
!===========================================================================================
! Fourier Antitransform the matrix elements g(q,Re,...) from Real Wannier representation
! (Re) to Reciprocal Wannier representation associated to electrons (k), so that:
!
!  g(q,k,...) = Sum_in_Re_WS[g(q,Re,...)*Exp[-i*q*Re]]/(NRe)
!===========================================================================================

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k
  integer,intent(in) :: irr_k(3,nrr_k),ndegen_k(nrr_k)
  real(dp),intent(in) :: kpoint(3)
  complex(dp),intent(in) :: g_Re(nrr_k,num_wann,num_wann,3*nat)
  complex(dp),intent(inout) :: g_wannier(num_wann,num_wann,3*nat)

  !local variables

  integer :: irk,imode,ibnd,jbnd
  real(dp) :: arg
  complex(dp) :: fac

  do imode=1,3*nat
     do ibnd=1,num_wann
        do jbnd=1,num_wann
           !
           do irk=1,nrr_k
              arg=tpi*dot_product(kpoint(:),dble(irr_k(:,irk)))
              fac=exp(cmplx_i*arg)/dble(ndegen_k(irk))
              g_wannier(ibnd,jbnd,imode)=g_wannier(ibnd,jbnd,imode)+ &
                                         fac*g_Re(irk,ibnd,jbnd,imode)
           enddo !irk
           !
        enddo !jbnd
     enddo !ibnd
  enddo !imode
  !
  return

end subroutine G_Reqw2kqw_no_spin
