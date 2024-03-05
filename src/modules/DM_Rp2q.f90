!----------------------------------------------------------------------------------
subroutine DM_Rp2q(nrr_q,irr_q,ndegen_q,qpoint,dyn_r,P_diag,omega)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!===========================================================================================
! 1) Fourier Antitransform of the Real Dynamical Matrix DM_r(3*nat,3*nat,R) to Reciprocal  !
! representation DM_q(3*nat,3*nat,q), so that:                                             !
!                                                                                          !
! DM_q(m,n,q)=Sum_in_R[DM_R(m,n,R)*Exp[i*q*R]]/ndeg(R)                                     !
!                                                                                          !
! 2) Diagonalization of the Reciprocal Dynamical Matrix DM_q(3*nat,3*nat,q) in order to    !
! obtain the eigenmodes and eigenvalues:                                                   !
!                                                                                          !
! W_q(m,n)=Sum_m',n'[P+(m',m)*DM_q(m,n)*P(n',n)]                                           !
!                                                                                          !
! W_k(m,n)=w2(n)*d(m-n)                                                                    !
!                                                                                          !
! But these eigenmodes have not physical meaning and have to be divide by sqrt(mass)       !
! in order to get the real displacements of the ions of these normal modes:                !
!                                                                                          !
! P(m,n)=P(m,n)/Sqrt[M(m)]                                                                 !
!                                                                                          !
! We use these displacements in order to recover physical meaning matrix elements from     !
! our canonical representation matrix elements                                             !
!===========================================================================================

  use kinds, only: dp
  use intw_reading, only: nat, amass, ityp
  use intw_useful_constants, only: tpi, cmplx_0, cmplx_i

  implicit none

  !I/O variables

  integer,intent(in) :: nrr_q
  integer,intent(in) :: irr_q(3,nrr_q)
  integer,intent(in) :: ndegen_q(nrr_q)
  real(dp),intent(in) :: qpoint(3)
  complex(dp),intent(in) :: dyn_r(3*nat,3*nat,nrr_q)
  complex(dp),intent(inout) :: P_diag(3*nat,3*nat)
  real(dp),intent(inout) :: omega(3*nat)

  !local variables

  integer :: irq,imode,jmode,na,ifail(3*nat),nfound,info
  integer :: iwork(5*3*nat)
  real(dp) :: arg,w2(3*nat),rwork(7*3*nat)
  complex(dp) :: fac,dyn_q_f(3*nat,3*nat)
  complex(dp) :: dyn_pack(3*nat*(3*nat+1)/2)
  complex(dp) :: cwork(2*3*nat)

  ! 1) Inverse Fourier transform of the Real Dynamical Matrix
  ! to Reciprocal representation
  !
  dyn_q_f(:,:)=cmplx_0
  !
  do imode=1,3*nat
     do jmode=1,3*nat
        !
        do irq=1,nrr_q
           !
           arg=tpi*dot_product(qpoint(:),dble(irr_q(:,irq)))
           fac=exp(cmplx_i*arg)/dble(ndegen_q(irq))
           dyn_q_f(imode,jmode)=dyn_q_f(imode,jmode)+fac*dyn_r(imode,jmode,irq)
           !
        enddo !irq
        !
     enddo !jmode
  enddo !imode
  !
  ! 2) We diagonalize the Reciprocal Dynamical Matrix in order to obtain the eigenmodes
  ! and squared eigenfrequencies at each q point of the phonon fine mesh (real vibrational
  ! phonon representation). We use the diagonalizing matrix to recover the physical meaning
  ! vibration modes
  !
  ! Diagonalise DM_k of the fine grid (->basis of eigenmodes) (based on Wannier90 knowhow)
  !
  dyn_pack(:)=cmplx_0
  P_diag(:,:)=cmplx_0
  w2(:)=0.0d0
  omega(:)=0.0d0
  !
  ! We pack the complex hamiltonian (upper triangular part for zhpevx)
  !
  do jmode=1,3*nat
     do imode=1,jmode
        !
        dyn_pack(imode+((jmode-1)*jmode)/2)=dyn_q_f(imode,jmode)
        !
     enddo !imode
  enddo !jmode
  !
  ! And we diagonalize
  !
  call zhpevx('V','A','U',3*nat,dyn_pack,0.0,0.0,0,0,-1.0,nfound,w2,P_diag,&
                                          3*nat,cwork,rwork,iwork,ifail,info)
  !
  do imode=1,3*nat
     !
     if (w2(imode).gt.0.d0) then
        omega(imode)=sqrt(abs(w2(imode)))
     else
        omega(imode)=-sqrt(abs(w2(imode)))
     endif
     omega(imode)=omega(imode)*0.001d0  ! from meV to eV
     !
     ! At this step, some folks calculate the real displacements (polarizations) of the atoms
     ! that are the vibrational eigenvectors divided by sqrt(amass)
     !
     do jmode=1,3*nat
        !
        na=((jmode-1)/3)+1
        P_diag(jmode,imode)=P_diag(jmode,imode)/sqrt(amass(ityp(na)))
        !
     enddo !jmode
     !
     ! But, if we calculate the displacements, the matrix made by them is not the change matrix P
     ! that diagonalize A (Dynamical Matrix) through P^-1*A*P, because this one is made by the
     ! eigenvectors (and not the displacements) of A
     ! But don't worry, because this displacements allow you to recover the physical meaning
     ! of the vibraional mode, so let's work with them
     !
  enddo !imode
  !
  return

end subroutine DM_Rp2q
