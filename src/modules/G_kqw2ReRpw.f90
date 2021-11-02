!----------------------------------------------------------------------------------
subroutine G_kqw2Reqw(nrr_k,irr_k,kmesh,g_kq,g_Req)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!==========================================================================================
! Fourier transform the matrix elements g(q,k,...) from Reciprocal Wannier representation !
! (k) to Real Wannier representation associated to electrons (Re), so that:               !
!                                                                                         !
! g(q,Re,...) = Sum_in_k_1BZ[g(q,k,...)*Exp[-i*k*Re]]/nkmesh                              !
!==========================================================================================

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
  integer,intent(in) :: irr_k(3,nrr_k)
  real(dp),intent(in) :: kmesh(3,nk1*nk2*nk3)
  complex(dp),intent(inout) :: g_Req(nq1*nq2*nq3,nrr_k,num_wann,num_wann,npol,npol,3*nat)
  complex(dp),intent(in) :: g_kq(nq1*nq2*nq3,nk1*nk2*nk3,num_wann,num_wann,npol,npol,3*nat)

  !local variables

  integer :: irk,ik,imode,ibnd,jbnd,ipol,jpol,iq,nkmesh,nqmesh
  real(dp) :: arg,term,rvec(3),len
  complex(dp) :: fac
  character(len=4) :: iq_loc
  character(len=256) :: file_name

  nkmesh=nk1*nk2*nk3
  nqmesh=nq1*nq2*nq3
  !
  do iq=1,nqmesh
     !
!!     if (iq<10)        write(iq_loc,"(i1)")iq
!!     if (10<=iq<100)   write(iq_loc,"(i2)")iq
!!     if (100<=iq<1000) write(iq_loc,"(i3)")iq
!!     !
!!     file_name=trim('./greq_decay_'//adjustl(iq_loc))
!!     open(unit=222,file=file_name,status='replace',action='readwrite')
     !
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              do ipol=1,npol
                 do jpol=1,npol
                    !
                    !$omp parallel default(none) &
                    !$omp shared(nrr_k,nkmesh,kmesh,irr_k,g_kq,g_Req,iq,imode,ibnd,jbnd,cmplx_i,tpi,ipol,jpol) &
                    !$omp private(irk,ik,arg,fac)
                    !
                    !$omp do
                    !
                    do irk=1,nrr_k
                       do ik=1,nkmesh
                          arg=tpi*dot_product(kmesh(:,ik),dble(irr_k(:,irk)))
                          fac=exp(-cmplx_i*arg)/dble(nkmesh)
                          g_Req(iq,irk,ibnd,jbnd,ipol,jpol,imode)=g_Req(iq,irk,ibnd,jbnd,ipol,jpol,imode)&
                                                                +fac*g_kq(iq,ik,ibnd,jbnd,ipol,jpol,imode)
                       enddo !ik
                    enddo !irk
                    !
                    !$omp end parallel
                    !
                 enddo !jpol
              enddo !ipol
           enddo !jbnd
        enddo !ibnd
     enddo !imode
     !
!!     do irk=1,nrr_k
!!        !
!!        rvec(:)=dble(irr_k(1,irk))*at(:,1)+dble(irr_k(2,irk))*at(:,2)+dble(irr_k(3,irk))*at(:,3)
!!        len=sqrt((rvec(1)**2.d0)+(rvec(2)**2.d0)+(rvec(3)**2.d0))
!!        !
!!        term=maxval(abs(g_ReRp(iq,irk,:,:,:,:,:)))
!!        write(222,'(3f16.8)') len*alat*bohr,term
!!        !
!!     enddo
!!     !
!!     close(222)
     !
  enddo !iq
  !
  return

end subroutine G_kqw2Reqw
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_kqw2Reqw_no_spin(nrr_k,irr_k,kmesh,g_kq,g_Req)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!==========================================================================================
! Fourier transform the matrix elements g(q,k,...) from Reciprocal Wannier representation !
! (k) to Real Wannier representation associated to electrons (Re), so that:               !
!                                                                                         !
! g(q,Re,...) = Sum_in_k_1BZ[g(q,k,...)*Exp[-i*k*Re]]/nkmesh                              !
!==========================================================================================

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
  integer,intent(in) :: irr_k(3,nrr_k)
  real(dp),intent(in) :: kmesh(3,nk1*nk2*nk3)
  complex(dp),intent(inout) :: g_Req(nq1*nq2*nq3,nrr_k,num_wann,num_wann,3*nat)
  complex(dp),intent(in) :: g_kq(nq1*nq2*nq3,nk1*nk2*nk3,num_wann,num_wann,3*nat)

  !local variables

  integer :: irk,ik,imode,ibnd,jbnd,ipol,jpol,iq,nkmesh,nqmesh
  real(dp) :: arg,term,rvec(3),len
  complex(dp) :: fac
  character(len=4) :: iq_loc
  character(len=256) :: file_name

  nkmesh=nk1*nk2*nk3
  nqmesh=nq1*nq2*nq3
  !
  do iq=1,nqmesh
     !
!!     if (iq<10)        write(iq_loc,"(i1)")iq
!!     if (10<=iq<100)   write(iq_loc,"(i2)")iq
!!     if (100<=iq<1000) write(iq_loc,"(i3)")iq
!!     !
!!     file_name=trim('./greq_decay_'//adjustl(iq_loc))
!!     open(unit=222,file=file_name,status='replace',action='readwrite')
     !
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              !
              !$omp parallel default(none) &
              !$omp shared(nrr_k,nkmesh,kmesh,irr_k,g_kq,g_Req,iq,imode,ibnd,jbnd,cmplx_i,tpi) &
              !$omp private(irk,ik,arg,fac)
              !
              !$omp do
              !
              do irk=1,nrr_k
                 do ik=1,nkmesh
                    arg=tpi*dot_product(kmesh(:,ik),dble(irr_k(:,irk)))
                    fac=exp(-cmplx_i*arg)/dble(nkmesh)
                    g_Req(iq,irk,ibnd,jbnd,imode)=g_Req(iq,irk,ibnd,jbnd,imode)&
                                                 +fac*g_kq(iq,ik,ibnd,jbnd,imode)
                 enddo !ik
              enddo !irk
              !
              !$omp end parallel
              !
           enddo !jbnd
        enddo !ibnd
     enddo !imode
     !
!!     do irk=1,nrr_k
!!        !
!!        rvec(:)=dble(irr_k(1,irk))*at(:,1)+dble(irr_k(2,irk))*at(:,2)+dble(irr_k(3,irk))*at(:,3)
!!        len=sqrt((rvec(1)**2.d0)+(rvec(2)**2.d0)+(rvec(3)**2.d0))
!!        !
!!        term=maxval(abs(g_ReRp(iq,irk,:,:,:)))
!!        write(222,'(3f16.8)') len*alat*bohr,term
!!        !
!!     enddo
!!     !
!!     close(222)
     !
  enddo !iq
  !
  return

end subroutine G_kqw2Reqw_no_spin
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_Reqw2ReRpw(nrr_q,nrr_k,irr_q,irr_k,qmesh,g_Req,g_ReRp)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!===========================================================================================
! Fourier transform the matrix elements g(q,Re,...) from Reciprocal Wannier representation !
! (q) to Real Wannier representation associated to phonons (Rp), so that:                  !
!                                                                                          !
! g(Rp,Re,...) = Sum_in_q_1BZ[g(q,Re,...)*Exp[-i*q*Rp]]/nqmesh                             !
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

  integer,intent(in) :: nrr_q,nrr_k
  integer,intent(in) :: irr_q(3,nrr_q),irr_k(3,nrr_k)
  real(dp),intent(in) :: qmesh(3,nq1*nq2*nq3)
  complex(dp),intent(in) :: g_Req(nq1*nq2*nq3,nrr_k,num_wann,num_wann,npol,npol,3*nat)
  complex(dp),intent(inout) :: g_ReRp(nrr_q,nrr_k,num_wann,num_wann,npol,npol,3*nat)

  !local variables

  integer :: irq,irk,iq,nqmesh,ibnd,jbnd,imode,ipol,jpol
  real(dp) :: arg,term,len1,len2,rvec1(3),rvec2(3)
  complex(dp) :: fac

  nqmesh=nq1*nq2*nq3
  !
!!  open(unit=222,file='grerp_decay',status='replace',action='readwrite')
  !
  do irk=1,nrr_k
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              do ipol=1,npol
                 do jpol=1,npol
                    !
                    !$omp parallel default(none) &
                    !$omp shared(nrr_q,nqmesh,qmesh,irr_q,g_ReRp,g_Req,irk,imode,ibnd,jbnd,cmplx_i,tpi,ipol,jpol) &
                    !$omp private(irq,iq,arg,fac)
                    !
                    !$omp do
                    !
                    do irq=1,nrr_q
                       do iq=1,nqmesh
                          arg=tpi*dot_product(qmesh(:,iq),dble(irr_q(:,irq)))
                          fac=exp(-cmplx_i*arg)/dble(nqmesh)
                          g_ReRp(irq,irk,ibnd,jbnd,ipol,jpol,imode)=g_ReRp(irq,irk,ibnd,jbnd,ipol,jpol,imode)&
                                                                  +fac*g_Req(iq,irk,ibnd,jbnd,ipol,jpol,imode)
                       enddo !iq
                    enddo !irq
                    !
                    !$omp end parallel
                    !
                 enddo !jpol
              enddo !ipol
           enddo !jbnd
        enddo !ibnd
     enddo !imode
     !
!!     rvec1(:)=dble(irr_k(1,irk))*at(:,1)+dble(irr_k(2,irk))*at(:,2)+dble(irr_k(3,irk))*at(:,3)
!!     len1=sqrt((rvec1(1)**2.d0)+(rvec1(2)**2.d0)+(rvec1(3)**2.d0))
!!     !
!!     do irq=1,nrr_q
!!        !
!!        rvec2(:)=dble(irr_q(1,irq))*at(:,1)+dble(irr_q(2,irq))*at(:,2)+dble(irr_q(3,irq))*at(:,3)
!!        len2=sqrt((rvec2(1)**2.d0)+(rvec2(2)**2.d0)+(rvec2(3)**2.d0))
!!        !
!!        term=maxval(abs(g_ReRp(irq,irk,:,:,:,:,:)))
!!        write(222,'(3f16.8)') len1*alat*bohr,len2*alat*bohr,term
!!     enddo
  !
  enddo !irk
  !
!!  close(222)
  !
  return

end subroutine G_Reqw2ReRpw
!----------------------------------------------------------------------------------
!**********************************************************************************
!----------------------------------------------------------------------------------
subroutine G_Reqw2ReRpw_no_spin(nrr_q,nrr_k,irr_q,irr_k,qmesh,g_Req,g_ReRp)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!===========================================================================================
! Fourier transform the matrix elements g(q,Re,...) from Reciprocal Wannier representation !
! (q) to Real Wannier representation associated to phonons (Rp), so that:                  !
!                                                                                          !
! g(Rp,Re,...) = Sum_in_q_1BZ[g(q,Re,...)*Exp[-i*q*Rp]]/nqmesh                             !
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

  integer,intent(in) :: nrr_q,nrr_k
  integer,intent(in) :: irr_q(3,nrr_q),irr_k(3,nrr_k)
  real(dp),intent(in) :: qmesh(3,nq1*nq2*nq3)
  complex(dp),intent(in) :: g_Req(nq1*nq2*nq3,nrr_k,num_wann,num_wann,3*nat)
  complex(dp),intent(inout) :: g_ReRp(nrr_q,nrr_k,num_wann,num_wann,3*nat)

  !local variables

  integer :: irq,irk,iq,nqmesh,ibnd,jbnd,imode,ipol,jpol
  real(dp) :: arg,term,len1,len2,rvec1(3),rvec2(3)
  complex(dp) :: fac

  nqmesh=nq1*nq2*nq3
  !
!!  open(unit=222,file='grerp_decay',status='replace',action='readwrite')
  !
  do irk=1,nrr_k
     do imode=1,3*nat
        do ibnd=1,num_wann
           do jbnd=1,num_wann
              !
              !$omp parallel default(none) &
              !$omp shared(nrr_q,nqmesh,qmesh,irr_q,g_ReRp,g_Req,irk,imode,ibnd,jbnd,cmplx_i,tpi) &
              !$omp private(irq,iq,arg,fac)
              !
              !$omp do
              !
              do irq=1,nrr_q
                 do iq=1,nqmesh
                    arg=tpi*dot_product(qmesh(:,iq),dble(irr_q(:,irq)))
                    fac=exp(-cmplx_i*arg)/dble(nqmesh)
                    g_ReRp(irq,irk,ibnd,jbnd,imode)=g_ReRp(irq,irk,ibnd,jbnd,imode)&
                                                   +fac*g_Req(iq,irk,ibnd,jbnd,imode)
                 enddo !iq
              enddo !irq
              !
              !$omp end parallel
              !
           enddo !jbnd
        enddo !ibnd
     enddo !imode
     !
!!     rvec1(:)=dble(irr_k(1,irk))*at(:,1)+dble(irr_k(2,irk))*at(:,2)+dble(irr_k(3,irk))*at(:,3)
!!     len1=sqrt((rvec1(1)**2.d0)+(rvec1(2)**2.d0)+(rvec1(3)**2.d0))
!!     !
!!     do irq=1,nrr_q
!!        !
!!        rvec2(:)=dble(irr_q(1,irq))*at(:,1)+dble(irr_q(2,irq))*at(:,2)+dble(irr_q(3,irq))*at(:,3)
!!        len2=sqrt((rvec2(1)**2.d0)+(rvec2(2)**2.d0)+(rvec2(3)**2.d0))
!!        !
!!        term=maxval(abs(g_ReRp(irq,irk,:,:,:)))
!!        write(222,'(3f16.8)') len1*alat*bohr,len2*alat*bohr,term
!!     enddo
  !
  enddo !irk
  !
!!  close(222)
  !
  return

end subroutine G_Reqw2ReRpw_no_spin

