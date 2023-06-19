!----------------------------------------------------------------------------------
subroutine H_kb2Rew(nrr_k,irr_k,kmesh,ham_rw)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!============================================================================================!
! 1) Wannier Rotation of the Hamiltonians H_k(n_b,n_b,k) in Reciprocal Bloch to Reciprocal   !
! Wannier representation Hw_k(n_w,n_w,k), so that:                                           !
!                                                                                            !
! Hw_k(m,n,k)=<wfc_w(m,k)|H|wfc_w(n,k)>                                                      !
! Hw_k(m,n,k)=Sum_m',n'[U+(m',m,k)*<wfc(m',k)|H|wfc(n',k)>*U(n',n,k)]                        !
! Hw_k(m,n,k)=Sum_m',n'[U+(m',m,k)*eig(n',k)*d(m'-n')*U(n',n,k)]                             !
! Hw_k(m,n,k)=Sum_m',n'[U+(n',m,k)*eig(n',k)*U(n',n,k)]                                      !
!                                                                                            !
! 2) Fourier Transform of the Reciprocal Wannier Hamiltonians Hw_k(n_w,n_w,k) to             !
! Real Wannier representation Hw_r(n_w,n_w,R), so that:                                      !
!                                                                                            !
! Hw_r(m,n,R)=Sum_in_k_1BZ[Hw_k(m,n,k)*Exp[-i*k*R]]/nkmesh                                   !
!=============================================================================================

  use kinds, only: dp
  use intw_input_parameters, only: nk1, nk2, nk3
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
  complex(dp),intent(inout) :: ham_rw(num_wann,num_wann,nrr_k)

  !local variables

  integer :: ik,irk,ibnd,jbnd,i,j,k,ikpt_k,nkmesh,G_plusk(3),GKQ_bz(3)
  real(dp) :: kpoint(3),kpoint_in_1bz(3),arg,term,QE_eig_k(num_bands)
  real(dp) :: rvec(3),len
  integer :: list_iGk(nG_max)
  complex(dp) :: wfc_k(nG_max,num_bands,nspin)
  complex(dp) :: ham_kw(num_wann,num_wann,nk1*nk2*nk3)
  complex(dp) :: fac

  nkmesh=nk1*nk2*nk3
  ham_kw=cmplx_0
  !
  ! 1) Calculate the KS Hamiltonians in Reciprocal Wannier Representation for
  ! each k point in electrons coarse kmesh
  !
  write(*,*)'1) We calculate Reciprocal Wannier Hamiltonians in coarse mesh'
  !
  do ik=1,nkmesh
     write(*,*)'ik = ',ik,'/',nkmesh
     !
     kpoint(:)=kmesh(:,ik)
     call get_psi_general_k_all_wfc(kpoint,list_iGk,wfc_k,QE_eig_k)
     !
     do jbnd=1,num_wann
        do ibnd=1,jbnd
           do i=1,num_bands
              !
              ham_kw(ibnd,jbnd,ik)=ham_kw(ibnd,jbnd,ik)&
                       +conjg(u_mesh(i,ibnd,ik))*QE_eig_k(i)*u_mesh(i,jbnd,ik)
              !
           enddo !i
           !
           ! we enforce hermiticity here (based in Wannier90 knowhow)
           !
           ham_kw(jbnd,ibnd,ik)=conjg(ham_kw(ibnd,jbnd,ik))
           !
        enddo !ibnd
     enddo !jbnd
  enddo !ik
  !
  write(*,*)' '
  !
  ! 2) Fourier Transform of the Reciprocal Wannier Hamiltonians to
  ! to the Real Wannier representation
  ! + check of the spatial decay of the Hamiltonian in Real Wannier basis
  !
  write(*,*)'2) Fourier transform: Reciprocal Wannier Hamiltonians -> Real Wannier Hamiltonians'
  !
!!  open(unit=222,file='h_r_decay',status='replace',action='readwrite')
  !
  do irk=1,nrr_k
     write(*,*)'irk = ',irk,'/',nrr_k
     !
     do ik=1,nkmesh
        !
        arg=tpi*dot_product(kmesh(:,ik),dble(irr_k(:,irk)))
        fac=exp(-cmplx_i*arg)/dble(nkmesh)
        ham_rw(:,:,irk)=ham_rw(:,:,irk)+fac*ham_kw(:,:,ik)
        !
     enddo !ik
     !
!!     rvec(:)=dble(irr_k(1,irk))*at(:,1)+dble(irr_k(2,irk))*at(:,2)+dble(irr_k(3,irk))*at(:,3)
!!     len=sqrt((rvec(1)**2.d0)+(rvec(2)**2.d0)+(rvec(3)**2.d0))
!!     !
!!     term=maxval(abs(ham_rw(:,:,irk)))
!!     write(222,'(3f16.8)') len*alat*bohr,term
     !
  enddo !irk
  !
!!  close(222)
  !
  return

end subroutine H_kb2Rew
!----------------------------------------------------------------------------------
