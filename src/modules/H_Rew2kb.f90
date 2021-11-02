!----------------------------------------------------------------------------------
subroutine H_Rew2kb(nrr_k,irr_k,ndegen_k,kpoint,ham_rw,P_diag,eig)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!===========================================================================================
! 1) Fourier Antitransform of the Real Wannier Hamiltonians Hw_r(n_w,n_w,R) to Reciprocal  !
! Wannier representation Hw_k(n_w,n_w,k), so that:                                         !
!                                                                                          !
! Hw_k(m,n,k)=Sum_in_R[Hw_R(m,n,R)*Exp[i*k*R]]/ndeg(R)                                     !
!                                                                                          !
! 2) Diagonalization of the Reciprocal Wannier Hamiltonians Hw_k(m_w,n_w,k) in order to    !
! obtain Reciprocal Bloch Hamiltonians (eigenvalues) and its respective unitary change     !
! matrix (eigenvalues)                                                                     !
!                                                                                          !
! H_k(m,n)=Sum_m'_n'[P+(m',m)*Hw_k(m',n')*P(n',n)]                                           !
!                                                                                          !
! H_k(m,n)=eig(n)*d(m-n)                                                                   !
!===========================================================================================

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
  complex(dp),intent(in) :: ham_rw(num_wann,num_wann,nrr_k)
  complex(dp),intent(inout) :: P_diag(num_wann,num_wann)
  real(dp),intent(inout) :: eig(num_wann)

  !local variables

  integer :: irk,ibnd,jbnd,ifail(num_wann),nfound,info
  integer :: iwork(5*num_wann)
  real(dp) :: arg,rwork(7*num_wann)
  complex(dp) :: fac,ham_kw_f(num_wann,num_wann)
  complex(dp) :: ham_pack(num_wann*(num_wann+1)/2)
  complex(dp) :: cwork(2*num_wann)

  ! 1) Inverse Fourier transform of the Real Wannier Hamiltonian
  ! to Reciprocal Wannier representation
  !
  ham_kw_f(:,:)=cmplx_0
  !
  do ibnd=1,num_wann
     do jbnd=1,num_wann
        !
        do irk=1,nrr_k
           !
           arg=tpi*dot_product(kpoint(:),dble(irr_k(:,irk)))
           fac=exp(cmplx_i*arg)/dble(ndegen_k(irk))
           ham_kw_f(ibnd,jbnd)=ham_kw_f(ibnd,jbnd)+fac*ham_rw(ibnd,jbnd,irk)
           !
        enddo !irk
        !
     enddo !jbnd
  enddo !ibnd
  !
  ! 2) We diagonalize the Reciprocal Wannier Bloch Hamiltonians in order to obtain eigenstates 
  ! and eigenvalues at each k point of the electron fine mesh (Reciprocal Bloch KS representation)
  ! We also keep the diagonalizing matrix in order to unrotate afterwards the matrix elements,
  ! that is Wannier unrotation of the Reciprocal Wannier matrix elements to Reciprocal Bloch
  ! matrix elements
  !
  ! Diagonalise Hw_k of the fine grid (->basis of eigenstates) (based on Wannier90 knowhow)
  !
  ham_pack(:)=cmplx_0
  P_diag(:,:)=cmplx_0
  eig(:)=0.0d0  ! in eV
  !
  ! We pack the complex hamiltonian (upper triangular part for zhpevx)
  !
  do jbnd=1,num_wann
     do ibnd=1,jbnd
        !
        ham_pack(ibnd+((jbnd-1)*jbnd)/2)=ham_kw_f(ibnd,jbnd)
        !
     enddo !ibnd
  enddo !jbnd
  !
  ! And we diagonalize
  !
  call zhpevx('V','A','U',num_wann,ham_pack,0.0_dp,0.0_dp,0,0,-1.0_dp,nfound,eig,P_diag,&
                                           num_wann,cwork,rwork,iwork,ifail,info)
  !
  return

end subroutine H_Rew2kb
