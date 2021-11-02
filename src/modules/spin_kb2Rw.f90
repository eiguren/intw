!----------------------------------------------------------------------------------
subroutine spin_kb2Rw(nrr_k,irr_k,kmesh,spin_kb,spin_Rw)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 18/07/2016
!
!===================================================================================!
! 1) Wannier Rotation of the spin components sb_k(k,n_b,3) in Reciprocal Bloch to   !
! Reciprocal Wannier representation sw_k(k,n_w,3), so that:                         !
!                                                                                   !
! sw_k(k,m)=<wfc_w(m,k)|S|wfc_w(m,k)>                                               !
! sw_k(k,m)=Sum_n[U+(n,m,k)*<wfc(n,k)|S|wfc(n,k)>*U(n,m,k)]                         !
!                                                                                   !
! 2) Fourier Transform of the Reciprocal Wannier spin components sw_k(k,n_w,3) to   !
! Real Wannier representation sw_r(R,n_w,3), so that:                               !
!                                                                                   !
! sw_r(R,m,3)=Sum_in_k_1BZ[sw_k(k,m,3)*Exp[-i*k*R]]/nkmesh                          !
!===================================================================================!

  use intw_reading
  use intw_useful_constants
  use intw_W90
  use intw_allwfcs
  use w90_parameters, only: num_wann,num_bands

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_k
  integer,intent(in) :: irr_k(3,nrr_k)
  real(dp),intent(in) :: kmesh(3,nk1*nk2*nk3)
  complex(dp),intent(in) :: spin_kb(nk1*nk2*nk3,num_bands,3)
  complex(dp),intent(out) :: spin_Rw(nrr_k,num_wann,num_wann,3)

  !local variables

  integer :: ik,irk,ibnd,jbnd,i,nkmesh,ipol
  real(dp) :: arg
  complex(dp) :: spin_kw(nk1*nk2*nk3,num_wann,num_wann,3)
  complex(dp) :: fac

  nkmesh=nk1*nk2*nk3
  spin_kw=cmplx_0
  !
  ! 1) Calculate the spin components in Reciprocal Wannier Representation for
  ! each k point in electrons coarse kmesh
  !
  write(*,*)'1) We calculate Reciprocal Wannier spin components in coarse mesh'
  !
  do ik=1,nkmesh
     !
     write(*,*)'ik = ',ik,'/',nkmesh
     !
     do ipol=1,3
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              do i=1,num_bands
                 !
                 spin_kw(ik,ibnd,jbnd,ipol)=spin_kw(ik,ibnd,jbnd,ipol) &
                                           +conjg(u_mesh(i,ibnd,ik))*spin_kb(ik,i,ipol)*u_mesh(i,jbnd,ik)
                 !
              enddo !i
           enddo !jbnd
        enddo !ibnd
     enddo !ipol
     !
  enddo !ik
  !
  write(*,*)' '
  !
  ! 2) Fourier Transform of the Reciprocal Wannier spin components to
  ! to the Real Wannier representation
  !
  write(*,*)'2) Fourier transform: Reciprocal Wannier spin comp. -> Real Wannier spin comp.'
  !
  do irk=1,nrr_k
     !
     write(*,*)'irk = ',irk,'/',nrr_k
     !
     do ibnd=1,num_wann
        do jbnd=1,num_wann
           do ipol=1,3
              !
              do ik=1,nkmesh
                 !
                 arg=tpi*dot_product(kmesh(:,ik),dble(irr_k(:,irk)))
                 fac=exp(-cmplx_i*arg)/dble(nkmesh)
                 spin_Rw(irk,ibnd,jbnd,ipol)=spin_Rw(irk,ibnd,jbnd,ipol) &
                                            +fac*spin_kw(ik,ibnd,jbnd,ipol)
                 !
              enddo !ik
              !
           enddo !ipol
        enddo !jbnd
     enddo !ibnd
     !
  enddo !irk
  !
  return

end subroutine spin_kb2Rw
!----------------------------------------------------------------------------------
