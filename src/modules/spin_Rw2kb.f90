!----------------------------------------------------------------------------------
subroutine spin_Rw2kb(nrr_k,irr_k,ndegen_k,kpoint,P_diag,spin_Rw,spin_kb)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 18/07/2016
!
!=======================================================================================!
! 1) Fourier Antitransform of the Real Wannier spin components spin_Rw(R,n_w,3) to      !
! Reciprocal Wannier representation spin_kw(k,n_w,3), so that:                          !
!                                                                                       !
! spin_kw(k,m,3)=Sum_in_R[spin_Rw(R,m,3)*Exp[i*k*R]]/ndeg(R)                            !
!                                                                                       !
! 2) Unrotation of the Reciprocal Wannier spin components spin_kw(k,n_w,3) in order to  !
! obtain Reciprocal Bloch spin components                                               !
!                                                                                       !
! spin_kb(k,n,3)=Sum_m[P+(m,n)*spin_kw(k,m,3)*P(m,n)]                                   !
!=======================================================================================!

  use kinds, only: dp
  use intw_reading
  use intw_useful_constants
  use intw_W90
  use w90_hamiltonian, only: irvec
  use w90_parameters, only: num_wann,num_bands

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k
  integer,intent(in) :: irr_k(3,nrr_k),ndegen_k(nrr_k)
  real(dp),intent(in) :: kpoint(3)
  complex(dp),intent(in) :: spin_Rw(nrr_k,num_wann,num_wann,3)
  complex(dp),intent(in) :: P_diag(num_wann,num_wann)
  complex(dp),intent(out) :: spin_kb(num_wann,3)

  !local variables

  integer :: irk,ibnd,jbnd,ipol,i
  real(dp) :: arg
  complex(dp) :: fac,spin_kw(num_wann,num_wann,3)

  ! 1) Inverse Fourier transform of the Real Wannier spin components
  ! to Reciprocal Wannier representation
  !
  spin_kw(:,:,:)=cmplx_0
  !
  do ipol=1,3
     do ibnd=1,num_wann
        do jbnd=1,num_wann
           !
           do irk=1,nrr_k
              !
              arg=tpi*dot_product(kpoint(:),dble(irr_k(:,irk)))
              fac=exp(cmplx_i*arg)/dble(ndegen_k(irk))
              spin_kw(ibnd,jbnd,ipol)=spin_kw(ibnd,jbnd,ipol)+fac*spin_Rw(irk,ibnd,jbnd,ipol)
              !
           enddo !irk
           !
        enddo !jbnd
     enddo !ibnd
  enddo !ipol
  !
  ! 2) We wannier-unrotate the Reciprocal Wannier spin components in order to obtain Reciprocal
  ! Bloch spin components at each k point of the fine mesh.
  !
  do ipol=1,3
     do i=1,num_wann
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              !
              spin_kb(i,ipol)=spin_kb(i,ipol)+ &
                              conjg(P_diag(ibnd,i))*spin_kw(ibnd,jbnd,ipol)*P_diag(jbnd,i)
              !
           enddo !jbnd
        enddo !ibnd
     enddo ! i
  enddo !ipol
  !
  return

end subroutine spin_Rw2kb
