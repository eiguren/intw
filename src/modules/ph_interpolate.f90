!----------------------------------------------------------------------------!
!       INTW project
!
! MBR 24/01/24
! intw_ph_interpolate contains routines to generate Wigner-Seitz real meshes,
! Fourier transforms of the dynamical matrix, interpolation at a given q-point.
! So far, it contains the necessary tools to interpolate phonons. 
!
! Those tools are similar to some routines in w90_setup, which have been 
! adapted here for phonon meshes.
!
! This is a mid result for EP interpolation. The goal is to have two ways of
! doing this: Wannier and Asier's method.
!----------------------------------------------------------------------------!
!
module intw_ph_interpolate

  use kinds, only: dp
  implicit none
  !
  ! variables
  public :: n_wss_q, n_ws_search_q     !!! search space for WS vectors
  public :: irvec_q, nrpts_q, ndegen_q   !!! these will substitute w90_hamiltonian: irvec, nrpts, ndegen
  public :: dyn_q, dyn_r        !!! dynamical matrix in qmesh and real space (WS vectors)
  public :: w2_q, u_q            !!! omega^q and eigenvectors in qmesh
  !
  ! subroutines
  public :: allocate_and_build_ws_irvec_q, allocate_and_build_dyn_qmesh, dyn_q_to_dyn_r, dyn_interp_1q, &
          dyn_diagonalize_1q , allocate_and_build_dyn_qmesh2
  !
  private
  !
  save
  !
  integer :: nrpts_q
  integer, parameter :: n_wss_q=27  !! TODO give somewhere as input
  integer, dimension(3) , parameter :: n_ws_search_q =(/ 1,1,1 /) !! TODO give somewhere as input
  integer, allocatable :: ndegen_q(:), irvec_q(:,:)
  real(kind=dp), allocatable :: w2_q(:,:)
  complex(kind=dp), allocatable :: dyn_r(:,:,:), dyn_q(:,:,:), u_q(:,:,:)
  !
  contains
  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ws_irvec_q ()
  !----------------------------------------------------------------------------!
  !  Calculate real-space Wigner-Seitz lattice vectors for phonon grid.
  !  Similar to w90_setup allocate_and_build_ws_irvec routine.
  !----------------------------------------------------------------------------!
  !
  use intw_reading, only: alat, at
  use intw_useful_constants, only: eps_6
  use intw_utility, only: generate_kmesh, cryst_to_cart, HPSORT_real
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_ph, only: nqmesh, qmesh
  !
  implicit none
  !
  integer :: iq, nws, i,j,k,l
  integer :: permu(n_wss_q), rdeg_ws_max(nq1*nq2*nq3)
  integer :: r_ws_max(3,n_wss_q,nq1*nq2*nq3)
  real(kind=dp) :: r_cryst(3,n_wss_q), r_length(n_wss_q), r_cart(3)
  !
  nws = 0  ! total number of WS vectors
  do iq = 1,nq1*nq2*nq3
      l = 0
      do i = -n_ws_search_q(1),n_ws_search_q(1)
      do j = -n_ws_search_q(2),n_ws_search_q(2)
      do k = -n_ws_search_q(3),n_ws_search_q(3)
         l = l + 1
         ! generate equivalents to point kmesh(:,ik)
         r_cryst(1,l) =  qmesh(1,iq) + real(i,dp)
         r_cryst(2,l) =  qmesh(2,iq) + real(j,dp)
         r_cryst(3,l) =  qmesh(3,iq) + real(k,dp)
         r_cart = r_cryst(:,l) * (/ nq1, nq2, nq3 /)
         !R-vector from crystallographic to cartesian
         call cryst_to_cart (1, r_cart, at, 1)
         r_cart = r_cart * alat  ! bohr units
         r_length(l) = sqrt ( sum(r_cart*r_cart) )
     end do
     end do
     end do
     ! order by ascending length
     call HPSORT_real(n_wss_q,r_length,permu)
     ! store first vector (shortest)
     r_ws_max(:,1,iq) = INT(r_cryst(:,permu(1)) * (/ nq1, nq2, nq3 /) )
     rdeg_ws_max(iq) = 1
     nws = nws + 1
     ! detect degeneracies and store vectors if degenerate
     do l = 2,n_wss_q
        if ( abs(r_length(l) - r_length(l-1))<eps_6) then
           r_ws_max(:,l,iq) = INT(r_cryst(:,permu(l)) * (/ nq1, nq2, nq3 /) )
           rdeg_ws_max(iq) = rdeg_ws_max(iq) + 1
           nws = nws + 1
        else
           exit
        end if
     end do
  end do
  !
  ! build WS kpoint list and degeneracies
  nrpts_q = nws
  allocate ( irvec_q(3,nrpts_q) )
  allocate ( ndegen_q(nrpts_q) )
  ! allocate and build WS kpoint list and degeneracies (irvec_q, ndegen_q)
  i = 0
  do iq = 1,nq1*nq2*nq3
     do l = 1, rdeg_ws_max(iq)
         i = i + 1
         irvec_q(:,i) = r_ws_max(:,l,iq)
         ndegen_q(i) = rdeg_ws_max(iq)
     end do
  end do
  !        
  return
  end subroutine allocate_and_build_ws_irvec_q
  !
  !----------------------------------------------------------------------------
  subroutine allocate_and_build_dyn_qmesh (fcfile)
  ! calculate the dynamical matrix, eiegnvectors and omega^2 on the full qmesh
  !----------------------------------------------------------------------------
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_reading, only:  alat, nat, ntyp, amass
  use intw_ph, only: nqmesh, qmesh, frc, mat_inv_four_t, readfc
  !
  implicit none
  !
  character(256) , intent(in) :: fcfile
  integer :: ir, iq, i,j,k
  real(dp) :: at_frc(3,3), qpoint(3)
  !
  ! Read the force constant matrix from the QE directory
  call readfc(fcfile,i,j,k, nat, alat, at_frc, ntyp, amass)
  !
  ! allocate and calculate dynmat at qmesh
  allocate( dyn_q(3*nat,3*nat,nqmesh), w2_q(3*nat,nqmesh), u_q(3*nat,3*nat,nqmesh) )
  !
  ! transform to dyn(q) and diagonalize
  do iq=1,nqmesh
     qpoint = qmesh(:,iq)
     ! mat_inv_four_t returns data in meV units
     call mat_inv_four_t(qpoint, nq1, nq2, nq3, 3*nat, frc, dyn_q(:,:,iq))
     call dyn_diagonalize_1q(3*nat, dyn_q(:,:,iq), u_q(:,:,iq), w2_q(:,iq))
  end do
  !
  return 
  end subroutine allocate_and_build_dyn_qmesh
  !
  !----------------------------------------------------------------------------
  subroutine allocate_and_build_dyn_qmesh2()
  ! Same as allocate_and_build_dyn_qmesh, but reading directly the prefix.dyn_q* 
  ! files generated by QE
  !----------------------------------------------------------------------------
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_reading, only:  alat, nat, ntyp, amass
  use intw_ph, only: nqmesh, qmesh, read_dynq
  !
  implicit none
  !
  integer :: ir, iq, i,j,k
  real(dp) :: at_frc(3,3), qpoint(3)
  !
  ! allocate and calculate dynmat at qmesh
  allocate( dyn_q(3*nat,3*nat,nqmesh), w2_q(3*nat,nqmesh), u_q(3*nat,3*nat,nqmesh) )
  !
  ! read dynmat
  call read_dynq(dyn_q)
  !
  ! diagonalize
  do iq=1,nqmesh
     qpoint = qmesh(:,iq)
     call dyn_diagonalize_1q(3*nat, dyn_q(:,:,iq), u_q(:,:,iq), w2_q(:,iq))
  end do
  !
  return
  end subroutine allocate_and_build_dyn_qmesh2
  !
  !----------------------------------------------------------------------------
  subroutine dyn_q_to_dyn_r ()
  ! allocate_and_build_ws_irvec_q and allocate_and_build_dyn_qmesh must have been previously run
  !----------------------------------------------------------------------------
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only:  nat
  use intw_ph, only: nqmesh, qmesh
  !
  implicit none
  !
  integer :: ir, iq
  complex(kind=dp) :: fac
  !
  allocate (dyn_r(3*nat,3*nat,nrpts_q))
  dyn_r = cmplx_0
  !
  do ir = 1,nrpts_q
     do iq = 1, nqmesh
        fac = exp(-cmplx_i*tpi*dot_product(qmesh(:,iq),irvec_q(:,ir)))
        dyn_r(:,:,ir) = dyn_r(:,:,ir) + fac*dyn_q(:,:,iq)
     end do
  end do
  dyn_r = dyn_r / real(nqmesh,dp)
  !
  return
  end subroutine dyn_q_to_dyn_r
  !
  !----------------------------------------------------------------------------
  subroutine dyn_interp_1q (qpoint, dyn_qint)
  ! allocate_and_build_ws_irvec_q iand  dyn_q_to_dyn_r ust have been previously run,
  ! as this uses irvec_q and dyn_r variables.
  !----------------------------------------------------------------------------
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only:  nat
  !
  implicit none
  !
  real(dp) , intent(in) :: qpoint(3)
  !
  complex(dp) , intent(out) :: dyn_qint(3*nat,3*nat)
  !
  integer :: ir, iq
  complex(kind=dp) :: fac
  !
  dyn_qint = cmplx_0
  do ir = 1,nrpts_q
     fac = exp(cmplx_i*tpi*dot_product(qpoint(:),irvec_q(:,ir)))/real(ndegen_q(ir),dp)
     dyn_qint = dyn_qint + fac * dyn_r(:,:,ir)
  end do
  !
  return
  end subroutine dyn_interp_1q
  !
  !----------------------------------------------------------------------------
  subroutine dyn_diagonalize_1q (n, dynq, uq, w2q)
  ! For a given dynq at a certain qpoint, with n=3*nat,
  ! this subroutine is a driver for zhpevx, which returns omega^2 and eigenvectors
  !
  ! TODO this is c+p Haritz's diagonalize_cmat, which should be in intw_utilities, since it is useful
  !----------------------------------------------------------------------------
  !
  implicit none
  !
  integer, intent(in)  :: n
  complex(dp) , intent(in) :: dynq(n,n)
  complex(dp) , intent(out) :: uq(n,n)
  real(dp) , intent(out) :: w2q(n)
  !
  integer :: i,j, nfound
  complex(dp) :: a_pack(n*(n+1)/2)
  complex(dp) :: cwork(2*n)
  real   (dp) :: rwork(7*n)
  integer     :: iwork(5*n), ifail(n), info
  !
  ! force hermiticity
  uq = 0.5_dp * ( dynq + transpose(conjg(dynq)) )
  !
  do j=1,n
    do i=1,j
      a_pack(i+((j-1)*j)/2)=uq(i,j)
    enddo
  enddo
  !
  call zhpevx('V', 'A', 'U', n, a_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp,  &
         nfound, w2q, uq, n, cwork, rwork, iwork, ifail, info)
  !
  return
  end subroutine dyn_diagonalize_1q

end module intw_ph_interpolate
