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
  !
  public :: irvec_qtau, nrpts_qtau, ndegen_qtau, nrpts_qtau12   !!! as above, but corrected with atom positions
  public :: dyn_rtau
  !
  public :: w2_q, u_q            !!! omega^q and eigenvectors in qmesh
  !
  ! subroutines
  public :: allocate_and_build_ws_irvec_q, allocate_and_build_dyn_qmesh, dyn_q_to_dyn_r, dyn_interp_1q, &
          dyn_diagonalize_1q , allocate_and_build_dyn_qmesh2
  !
  public :: allocate_and_build_ws_irvec_qtau, dyn_q_to_dyn_rtau, dyn_interp_1q_tau !!corrected with atom positions
  !
  private
  !
  save
  !
  integer :: nrpts_q
  integer, parameter :: n_wss_q=27  !! TODO give somewhere as input
  integer, dimension(3) , parameter :: n_ws_search_q = (/ 1,1,1 /) !! TODO give somewhere as input
  integer, allocatable :: ndegen_q(:), irvec_q(:,:)
  real(kind=dp), allocatable :: w2_q(:,:)
  complex(kind=dp), allocatable :: dyn_r(:,:,:), dyn_q(:,:,:), u_q(:,:,:) ! in a.u. (without the mass factor)
  !
  integer :: nrpts_qtau
  integer, allocatable :: ndegen_qtau(:,:,:), irvec_qtau(:,:,:,:), nrpts_qtau12(:,:)
  complex(kind=dp), allocatable :: dyn_rtau(:,:,:)
  !
  contains
  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ws_irvec_q_hold ()
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
  integer :: r_ws_max(3,n_wss_q,nq1*nq2*nq3), r_cryst_int(3,n_wss_q)
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
         r_cryst_int(:,l) = nint(r_cart)
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
     r_ws_max(:,1,iq) = r_cryst_int(:,permu(1))
     rdeg_ws_max(iq) = 1
     nws = nws + 1
     ! detect degeneracies and store vectors if degenerate
     do l = 2,n_wss_q
       if ( abs(r_length(l) - r_length(l-1))<eps_6) then
           r_ws_max(:,l,iq) = r_cryst_int(:,permu(l))
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
  print *,' #nrpts_q= ', nrpts_q
  !
  return
  !
  end subroutine allocate_and_build_ws_irvec_q_hold
  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ws_irvec_qtau_hold()
  !----------------------------------------------------------------------------!
  !  Calculate real-space Wigner-Seitz lattice vectors for phonon grid.
  !  Similar to w90_setup allocate_and_build_ws_irvec routine.
  !  Unlike the allocate_and_build_ws_irvec_q version, here we choose as WS
  !  criterion that, for each tau1,tau2 atom pair, the WS cell is centered at tau1.
  !  So, we impose the truncation using R+tau2-tau1,
  !  as in S. Poncé et al, Phys. Rev. Research 3, 043022 (2021)  (appendix D)
  !  and  G. Pizzi et al, J. Phys.: Condens. Matter 32 165902 (2020) (section 4.2).
  !----------------------------------------------------------------------------!
  !
  !
  use intw_reading, only: alat, at, nat, tau_cryst
  use intw_useful_constants, only: eps_6
  use intw_utility, only: cryst_to_cart, HPSORT_real
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_ph, only: nqmesh, qmesh
  !
  implicit none
  !
  integer :: iq, nws, i,j,k,l, iat1, iat2
  integer :: permu(n_wss_q), rdeg_ws_max(nq1*nq2*nq3,nat,nat)
  integer :: r_ws_max(3,n_wss_q,nq1*nq2*nq3,nat,nat), r_cryst_int(3,n_wss_q)
  real(kind=dp) :: r_cryst(3,n_wss_q), r_length(n_wss_q), r_cart(3)

  nrpts_qtau = 0 !max. number of WS vectors among all pairs
  allocate( nrpts_qtau12(nat,nat) )
  do iat1=1,nat
  do iat2=1,nat
  !
  nws = 0  ! total number of WS vectors for this atom pair
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
         r_cryst_int(:,l) = nint(  r_cryst(:,l) * (/ nq1, nq2, nq3 /) )
         r_cart = r_cryst(:,l) * (/ nq1, nq2, nq3 /) + tau_cryst(:,iat2) - tau_cryst(:,iat1)
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
     r_ws_max(:,1,iq,iat1,iat2) = r_cryst_int(:,permu(1))
     rdeg_ws_max(iq,iat1,iat2) = 1
     nws = nws + 1
     ! detect degeneracies and store vectors if degenerate
     do l = 2,n_wss_q
        if ( abs(r_length(l) - r_length(l-1))<eps_6) then
           r_ws_max(:,l,iq,iat1,iat2) = r_cryst_int(:,permu(l))
           rdeg_ws_max(iq,iat1,iat2) = rdeg_ws_max(iq,iat1,iat2) + 1
           nws = nws + 1
        else
           exit
        end if
     end do
  end do
  !
  nrpts_qtau12(iat1,iat2) = nws
  if (nws > nrpts_qtau) nrpts_qtau = nws
  !
  end do !iat2
  end do !iat1
  !
  ! build WS kpoint list and degeneracies
  allocate ( irvec_qtau(3,nrpts_qtau,nat,nat) )
  allocate ( ndegen_qtau(nrpts_qtau,nat,nat) )
  ! allocate and build WS kpoint list and degeneracies (irvec_q, ndegen_q)
  do iat1 = 1,nat
  do iat2 = 1,nat
  !
  i = 0
  do iq = 1,nq1*nq2*nq3
     do l = 1, rdeg_ws_max(iq,iat1,iat2)
         i = i + 1
         irvec_qtau(:,i,iat1,iat2) = r_ws_max(:,l,iq,iat1,iat2)
         ndegen_qtau(i,iat1,iat2) = rdeg_ws_max(iq,iat1,iat2)
     end do
  end do
  !
  end do !iat2
  end do !iat1
  !
  ! print
  print *, '#max nrpts_q = ',  nrpts_qtau12
  !do iat1=1,nat
  !do iat2=1,nat
  !   print *, '#at1, at2, nrpts_q = ', iat1, iat2, nrpts_qtau12(iat1,iat2)
  !   print *, '# irvec, degen:'
  !   do i=1,nrpts_qtau12(iat1,iat2)
  !      write(*,'(3i5,3x,i5)') irvec_qtau(:,i,iat1,iat2), ndegen_qtau(i,iat1,iat2)
  !   end do
  !   print *, '  '
  !   print *, '  '
  !end do
  !end do
  !
  return
  end subroutine allocate_and_build_ws_irvec_qtau_hold
  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ws_irvec_q ()  !NEW VERSION
  !----------------------------------------------------------------------------!
  !  Calculate real-space Wigner-Seitz lattice vectors for phonon grid.
  !  Similar to w90_setup allocate_and_build_ws_irvec routine.
  !----------------------------------------------------------------------------!
  !
  use intw_reading, only: at, alat
  use intw_useful_constants, only: eps_8
  use intw_utility, only:  cryst_to_cart
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_ph, only: nqmesh, qmesh
  !
  implicit none
  !
  integer :: iq, i,j,k,l, l0,l1, nboundary
  logical :: in_ws
  integer :: Rs(3,n_wss_q), r_cryst_int(3), ndegen_ws(nq1*nq2*nq3*n_wss_q), irvec_ws(3,nq1*nq2*nq3*n_wss_q)
  real(kind=dp) :: r_cryst(3), r_length_l, r_length_l1, r_cart(3)
  !
    ! generate superlattice replica vectors search mesh
  l = 0
  do i = -n_ws_search_q(1),n_ws_search_q(1)
  do j = -n_ws_search_q(2),n_ws_search_q(2)
  do k = -n_ws_search_q(3),n_ws_search_q(3)
         l = l + 1
         Rs(:,l) = (/ i,j,k /)
         if (i == 0 .and. j == 0 .and. k == 0) l0=l ! Origin O
  end do
  end do
  end do
  !
  nrpts_q = 0  ! total number of WS vectors
  do iq = 1,nq1*nq2*nq3
     !
     do l = 1, n_wss_q
        ! r-R(l), where for r-supercell-vector I use a conventional cell mesh of size nq1, nq2, nq3
        ! and R(l) runs over replicas
        r_cryst = ( qmesh(:,iq) - real(Rs(:,l),dp) ) * real( (/ nq1, nq2, nq3 /), dp)
        r_cryst_int = nint(r_cryst)
        !R-vector from crystallographic to cartesian
        r_cart = r_cryst
        call cryst_to_cart (1, r_cart, at, 1)
        r_length_l = alat * sqrt ( sum(r_cart*r_cart) )  ! distance of r-R(l) to O (cartesian, bohr)
        !
        ! r-R(l) is in the WS if its distance to O is shorter than its
        ! distance to any other O' origin.
        ! If it is equidistant, it lies on the boundary and is degenerate.
        in_ws = .true.
        nboundary = 1
        !
        ! Loop over origins O' given by R(l1)
        do l1 = 1, n_wss_q
          ! r-R(l)-R(l1)
          r_cryst = ( qmesh(:,iq) - real(Rs(:,l)+Rs(:,l1),dp) ) * real( (/ nq1, nq2, nq3 /), dp)
          r_cart = r_cryst
          call cryst_to_cart (1, r_cart, at, 1)
          r_length_l1 = alat * sqrt ( sum(r_cart*r_cart) )  ! distance of r-R(l) to O' (cartesian, bohr)
          ! compare distances leaving a gap eps_8
             ! TODO !!! put tolerance as parameter. It depends a lot on this!
             ! I guess that we need a smaller one the less nq...
          if ( r_length_l > r_length_l1 + eps_8*1000. .and. l1/=l0) then ! not in the WS => remove vector from list
                  in_ws = .false.
                  exit
          else if ( abs(r_length_l-r_length_l1)<=eps_8*1000. .and. l1/=l0) then ! on the boundary => add degeneracy
                  nboundary = nboundary + 1
          end if
        end do
        !
        ! store r-R(l) and its degeneracy if it is inside WS
        if (in_ws) then
             nrpts_q=nrpts_q+1
             irvec_ws(:,nrpts_q) = r_cryst_int
             ndegen_ws(nrpts_q) = nboundary
        end if
        !
     end do
  end do ! iq
  !
  ! Data for Wannier: WS kpoint list and degeneracies.
  ! Simply dismiss the array at sites >nrpts, which have not been used
  allocate ( irvec_q(3,nrpts_q) )
  allocate ( ndegen_q(nrpts_q) )
  ndegen_q = ndegen_ws(1:nrpts_q)
  do i=1,3
       irvec_q(i,:) = irvec_ws(i,1:nrpts_q)
  end do
  print *, '#max nrpts_q = ',  nrpts_q
  !
  return
  !
  end subroutine allocate_and_build_ws_irvec_q
  !
    !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ws_irvec_qtau ()  !NEW VERSION
  !----------------------------------------------------------------------------!
  !  Calculate real-space Wigner-Seitz lattice vectors for phonon grid.
  !  Similar to w90_setup allocate_and_build_ws_irvec routine.
  !  Unlike the allocate_and_build_ws_irvec_q version, here we choose as WS
  !  criterion that, for each tau1,tau2 atom pair, the WS cell is centered at tau1.
  !  So, we impose the truncation using R+tau2-tau1,
  !  as in S. Poncé et al, Phys. Rev. Research 3, 043022 (2021)  (appendix D)
  !  and  G. Pizzi et al, J. Phys.: Condens. Matter 32 165902 (2020) (section 4.2).
  !----------------------------------------------------------------------------!
  !
  !
  use intw_reading, only: at, alat, nat, tau_cryst
  use intw_useful_constants, only: eps_8
  use intw_utility, only: cryst_to_cart, HPSORT_real
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_ph, only: nqmesh, qmesh
  !
  implicit none
  !
  integer :: iq, i,j,k,l, l0,l1, nboundary, iat1, iat2, nws
  logical :: in_ws
  integer :: Rs(3,n_wss_q), r_cryst_int(3), ndegen_ws(nq1*nq2*nq3*n_wss_q,nat,nat), irvec_ws(3,nq1*nq2*nq3*n_wss_q,nat,nat)
  real(kind=dp) :: r_cryst(3), r_length_l, r_length_l1, r_cart(3)
  !
  ! generate superlattice replica vectors search mesh
  l = 0
  do i = -n_ws_search_q(1),n_ws_search_q(1)
  do j = -n_ws_search_q(2),n_ws_search_q(2)
  do k = -n_ws_search_q(3),n_ws_search_q(3)
         l = l + 1
         Rs(:,l) = (/ i,j,k /)
         if (i == 0 .and. j == 0 .and. k == 0) l0=l ! Origin O
  end do
  end do
  end do
  !
  nrpts_q = 0
  !
  allocate( nrpts_qtau12(nat,nat) )
  !
  do iat1=1,nat
  do iat2=1,nat
     nws = 0  ! total number of WS vectors for this atoms pair
     do iq = 1,nq1*nq2*nq3
     !
        do l = 1, n_wss_q
           ! r-R(l), where for r-supercell-vector I use a conventional cell mesh of size nq1, nq2, nq3
           ! and R(l) runs over replicas
           r_cryst = ( qmesh(:,iq) - real(Rs(:,l),dp) ) * real( (/ nq1, nq2, nq3 /), dp)
           r_cryst_int = nint(r_cryst)
           r_cart = r_cryst + tau_cryst(:,iat2) - tau_cryst(:,iat1)  ! add interatomic distance criterion (O is at iat1 position)
           !R-vector from crystallographic to cartesian
           call cryst_to_cart (1, r_cart, at, 1)
           r_length_l = alat * sqrt ( sum(r_cart*r_cart) )  ! distance of r-R(l) to O  (cartesian, bohr)
           !
           ! r-R(l) is in the WS if its distance to O is shorter than its
           ! distance to any other O' origin.
           ! If it is equidistant, it lies on the boundary and is degenerate.
           in_ws = .true.
           nboundary = 1
           !
           ! Loop over origins O' given by R(l1)
           do l1 = 1, n_wss_q
             ! r-R(l)-R(l1)
             r_cryst = ( qmesh(:,iq) - real(Rs(:,l)+Rs(:,l1),dp) ) * real( (/ nq1, nq2, nq3 /), dp)
             r_cart = r_cryst + tau_cryst(:,iat2) - tau_cryst(:,iat1)  ! add interatomic distance to distance criterion (O' is at l1+iat1 position)
             call cryst_to_cart (1, r_cart, at, 1)
             r_length_l1 =  alat * sqrt ( sum(r_cart*r_cart) )  ! distance of r-R(l) to O' (cartesian, bohr)
             ! compare distances leaving a gap eps_8*1000.
             ! TODO !!! put tolerance as parameter. It depends a lot on this!
             ! I guess that we need a smaller one the less nq...
             if ( r_length_l > r_length_l1 + eps_8*1000. .and. l1/=l0) then ! not in the WS => remove vector from list
                  in_ws = .false.
                  exit
             else if ( abs(r_length_l-r_length_l1)<=eps_8*1000. .and. l1/=l0) then ! on the boundary => add degeneracy
                  nboundary = nboundary + 1
             end if
           end do
           !
           ! store r-R(l) and its degeneracy if it is inside WS
           if (in_ws) then
             nws=nws+1
             irvec_ws(:,nws,iat1,iat2) = r_cryst_int
             ndegen_ws(nws,iat1,iat2) = nboundary
           end if
           !
        end do ! l
      end do ! iq
      !
      nrpts_qtau12(iat1,iat2) = nws
      !
  end do !iat2
  end do !iat1
  !
  ! max. number of WS vectors among all pairs is used as irvec_q array dimension
  nrpts_qtau = maxval(nrpts_qtau12)
  !
  ! Data for Wannier: WS kpoint list and degeneracies.
  ! Simply dismiss the array at sites >nrpts, which have not been used
  allocate ( irvec_qtau(3,nrpts_qtau,nat,nat) )
  allocate ( ndegen_qtau(nrpts_qtau,nat,nat) )
  ndegen_qtau = 0
  irvec_qtau = 0
  do iat1 = 1,nat
  do iat2 = 1,nat
     nws = nrpts_qtau12(iat1,iat2)
     ndegen_qtau(1:nws,iat1,iat2) = ndegen_ws(1:nws,iat1,iat2)
     do i=1,3
       irvec_qtau(i,1:nws,iat1,iat2) = irvec_ws(i,1:nws,iat1,iat2)
     end do
  end do
  end do
  !
  print *, '#max nrpts_q = ',  nrpts_qtau12
  !
  return
  !
  end subroutine allocate_and_build_ws_irvec_qtau
  !----------------------------------------------------------------------------
  subroutine allocate_and_build_dyn_qmesh (fcfile)
  ! calculate the dynamical matrix, eiegnvectors and omega^2 on the full qmesh
  !----------------------------------------------------------------------------
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_input_parameters, only: nq1, nq2, nq3
  use intw_reading, only: nat
  use intw_ph, only: nqmesh, qmesh, mat_inv_four_t, readfc
  !
  implicit none
  !
  character(256) , intent(in) :: fcfile
  integer :: ir, iq, i,j,k
  real(dp) :: at_frc(3,3), qpoint(3)
  complex(dp) :: frc(nq1,nq2,nq3,3,3,nat,nat)
  !
  ! Read the force constant matrix from the QE directory
  call readfc(fcfile, frc)
  !
  ! allocate and calculate dynmat at qmesh
  allocate( dyn_q(3*nat,3*nat,nqmesh), w2_q(3*nat,nqmesh), u_q(3*nat,3*nat,nqmesh) )
  !
  ! transform to dyn(q) and diagonalize
  do iq=1,nqmesh
     qpoint = qmesh(:,iq)
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
  subroutine dyn_q_to_dyn_rtau ()
  ! allocate_and_build_ws_irvec_q and allocate_and_build_dyn_qmesh must have been previously run
  !----------------------------------------------------------------------------
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only:  nat
  use intw_ph, only: nqmesh, qmesh
  !
  implicit none
  !
  integer :: ir, iq, iat1, iat2
  complex(kind=dp) :: fac
  !
  allocate (dyn_rtau(3*nat,3*nat,nrpts_qtau))
  dyn_rtau = cmplx_0
  !
  do iat1=1,nat
  do iat2=1,nat
     do ir = 1,nrpts_qtau12(iat1,iat2)
        do iq = 1, nqmesh
           fac = exp(-cmplx_i*tpi*dot_product(qmesh(:,iq),irvec_qtau(:,ir,iat1,iat2)))
           dyn_rtau((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, ir ) = &
                   dyn_rtau((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, ir ) + &
                   fac*dyn_q((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, iq)
        end do
     end do
  end do
  end do
  dyn_rtau = dyn_rtau / real(nqmesh,dp)
  !
  return
  end subroutine dyn_q_to_dyn_rtau
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
  subroutine dyn_interp_1q_tau (qpoint, dyn_qint)
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
  integer :: ir, iq, iat1, iat2
  complex(kind=dp) :: fac
  !
  dyn_qint = cmplx_0
  do iat1 = 1,nat
  do iat2 = 1,nat
     do ir = 1,nrpts_qtau12(iat1,iat2)
       fac = exp(cmplx_i*tpi*dot_product(qpoint(:),irvec_qtau(:,ir,iat1,iat2)))/real(ndegen_qtau(ir,iat1,iat2),dp)
       dyn_qint((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3) = &
           dyn_qint((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3) + &
           fac * dyn_rtau((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3,ir)
     end do
  end do
  end do
  !
  return
  end subroutine dyn_interp_1q_tau
  !
  !----------------------------------------------------------------------------
  subroutine dyn_diagonalize_1q (n, dynq, uq, w2q)
    ! For a given dynq at a certain qpoint, with n=3*nat,
    ! this subroutine is a driver for zhpevx, which returns omega^2 and eigenvectors
    !
    ! TODO this is c+p Haritz's diagonalize_cmat, which should be in intw_utilities, since it is useful
    !----------------------------------------------------------------------------
    !
    use intw_reading, only: nat, amass, ityp
    use intw_utility, only: diagonalize_cmat
    !
    implicit none
    !
    integer, intent(in)  :: n
    complex(dp), intent(in) :: dynq(n,n) ! Dynamical matrix in a.u. (without the mass factor)
    complex(dp), intent(out) :: uq(n,n) ! Phonon polarization vectors in a.u.
    real(dp), intent(out) :: w2q(n) ! Phonon frequencies^2 in a.u.
    !
    integer :: iat1, iat2
    real(dp), parameter :: pmass = 1822.88848426_dp
    real(dp), parameter :: aumev = 27211.396132_dp
    !
    !
    ! Add mass factor
    do iat1 = 1,nat
      do iat2 = 1,nat
        !
        uq( (iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3 ) = dynq( (iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3 ) &
                / sqrt( amass(ityp(iat1)) * amass(ityp(iat2)) ) / pmass ! in a.u. (Hartree/Bohr^2*/me = Hartree^2/hbar^2)
        !
      end do
    end do !atoms
    !
    ! force hermiticity
    uq = 0.5_dp * ( uq + transpose(conjg(uq)) )
    !
    call diagonalize_cmat(uq, w2q)

  end subroutine dyn_diagonalize_1q

end module intw_ph_interpolate
