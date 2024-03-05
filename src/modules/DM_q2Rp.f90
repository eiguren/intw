!----------------------------------------------------------------------------------
subroutine DM_q2Rp(nrr_q,irr_q,dyn_r)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!====================================================================================
! Fourier Transform to the Reciprocal Dynamical Matrix DM_q(3*nat,3*nat,q) to Real  !
! representation DM_r(3*nat,3*nat,R), so that:                                      !
!                                                                                   !
! DM_r(m,n,R)=Sum_in_q_1BZ[DM_q(m,n,k)*Exp[-i*q*R]]/nqmesh                          !
!====================================================================================
  use kinds, only: dp
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_reading, only: nat
  use intw_useful_constants, only: tpi, cmplx_0, cmplx_i
  use intw_ph, only: nqmesh, qmesh, frc

  implicit none

  !I/O vaiables

  integer,intent(in) :: nrr_q
  integer,intent(in) :: irr_q(3,nrr_q)
  complex(dp),intent(inout) :: dyn_r(3*nat,3*nat,nrr_q)

  !local variables

  integer :: iq, irq
  real(dp) :: qpoint(3), arg
  complex(dp) :: dyn_q(3*nat,3*nat,nqmesh), dynm(3*nat,3*nat), fac

  nqmesh=nq1*nq2*nq3
  dyn_q(:,:,:)=cmplx_0
  !
  ! 1) We have to calculate the Dynamical Matrix
  ! for each q point in phonons coarse qmesh
  !
  write(*,*)'1) We calculate the Reciprocal Dynamical Matrices in coarse mesh'
  !
  do iq=1,nqmesh
     write(*,*)'iq = ',iq,'/',nqmesh
     !
     qpoint(:)=qmesh(:,iq)
     !
     call mat_inv_four_t(qpoint,nq1,nq2,nq3,3*nat,frc,dynm)
     !
     dyn_q(:,:,iq)=dynm(:,:)
     !
  enddo !iq
  !
  write(*,*)' '
  !
  ! 2) Fourier transform of the Reciprocal Dynamical Matrix to
  ! to the Real Representation
  ! + check of the spatial decay of Dynamical Matrix in Real space
  !
  write(*,*)'2) Fourier transform: Reciprocal Dynamical Matrix -> Real Dynamical Matrix'
  !
!!  open(unit=222,file='dm_r_decay',status='replace',action='readwrite')
  !
  do irq=1,nrr_q
     write(*,*)'irq = ',irq,'/',nrr_q
     !
     do iq=1,nqmesh
        !
        arg=tpi*dot_product(qmesh(:,iq),dble(irr_q(:,irq)))
        fac=exp(-cmplx_i*arg)/dble(nqmesh)
        dyn_r(:,:,irq)=dyn_r(:,:,irq)+fac*dyn_q(:,:,iq)
        !
     enddo !iq
     !
!!     rvec(:)=dble(irr_q(1,irq))*at(:,1)+dble(irr_q(2,irq))*at(:,2)+dble(irr_q(3,irq))*at(:,3)
!!     len=sqrt((rvec(1)**2.d0)+(rvec(2)**2.d0)+(rvec(3)**2.d0))
!!     !
!!     term=maxval(abs(dyn_r(:,:,irq)))
!!     write(222,'(3f16.8)') len*alat*bohr, term
     !
  enddo !irq
  !
!!  close(222)
  !
  return

end subroutine DM_q2Rp
!----------------------------------------------------------------------------------
