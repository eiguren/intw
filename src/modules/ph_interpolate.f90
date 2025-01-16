!
! Copyright (C) 2024 INTW group
!
! This file is part of INTW.
!
! INTW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! INTW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
module intw_ph_interpolate

  !-------------------------------------------------------------------------------!
  ! MBR 24/01/24                                                                  !
  ! intw_ph_interpolate contains routines to generate Wigner-Seitz real meshes,   !
  ! Fourier transforms of the dynamical matrix, interpolation at a given q-point. !
  ! So far, it contains the necessary tools to interpolate phonons.               !
  !                                                                               !
  ! Those tools are similar to some routines in w90_setup, which have been        !
  ! adapted here for phonon meshes.                                               !
  !                                                                               !
  ! This is a mid result for EP interpolation. The goal is to have two ways of    !
  ! doing this: Wannier and Asier's method.                                       !
  !-------------------------------------------------------------------------------!

  use kinds, only: dp

  implicit none
  !
  ! variables
  public :: n_wss_q, n_ws_search_q
  public :: irvec_q, nrpts_q, ndegen_q
  public :: irvec_qtau, nrpts_qtau, ndegen_qtau, nrpts_qtau12
  public :: dyn_q, dyn_r
  public :: w2_q, u_q
  !
  ! subroutines
  public :: allocate_and_build_ws_irvec_q, allocate_and_build_ws_irvec_qtau, &
            allocate_and_build_dyn_qmesh, allocate_and_build_dyn_qmesh_from_fc, &
            dyn_q_to_dyn_r, dyn_interp_1q, dyn_diagonalize_1q
  !
  private
  !
  save
  !
  !!! these will substitute w90_hamiltonian: irvec, nrpts, ndegen
  integer :: nrpts_q
  integer, allocatable :: ndegen_q(:), irvec_q(:,:)

  !!! as above, but corrected with atom positions
  integer :: nrpts_qtau
  integer, allocatable :: ndegen_qtau(:,:,:), irvec_qtau(:,:,:,:), nrpts_qtau12(:,:)

  !!! search space for WS vectors
  integer, parameter :: n_wss_q=27  !! TODO give somewhere as input
  integer, dimension(3) , parameter :: n_ws_search_q = (/ 1,1,1 /) !! TODO give somewhere as input

  !!! dynamical matrix in qmesh and real space (WS vectors)
  complex(kind=dp), allocatable :: dyn_q(:,:,:), dyn_r(:,:,:)

  !!! omega^q and eigenvectors in qmesh
  real(kind=dp), allocatable :: w2_q(:,:)
  complex(kind=dp), allocatable :: u_q(:,:,:) ! in a.u. (without the mass factor)
  !
  !
  contains

  subroutine allocate_and_build_ws_irvec_q() !NEW VERSION
    !  Calculate real-space Wigner-Seitz lattice vectors for phonon grid.
    !  Similar to w90_setup allocate_and_build_ws_irvec routine.
    !----------------------------------------------------------------------------!
    use intw_reading, only: at, alat
    use intw_useful_constants, only: eps_8
    use intw_utility, only:  cryst_to_cart
    use intw_input_parameters, only: nq1, nq2, nq3
    use intw_ph, only: qmesh
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
        call cryst_to_cart(1, r_cart, at, 1)
        r_length_l = alat * sqrt( sum(r_cart*r_cart) )  ! distance of r-R(l) to O (cartesian, bohr)
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
          call cryst_to_cart(1, r_cart, at, 1)
          r_length_l1 = alat * sqrt( sum(r_cart*r_cart) )  ! distance of r-R(l) to O' (cartesian, bohr)
          ! compare distances leaving a gap eps_8
            ! TODO !!! put tolerance as parameter. It depends a lot on this!
            ! I guess that we need a smaller one the less nq...
          if ( r_length_l > r_length_l1 + eps_8*1000. .and. l1/=l0 ) then ! not in the WS => remove vector from list
            in_ws = .false.
            exit
          else if ( abs(r_length_l-r_length_l1)<=eps_8*1000. .and. l1/=l0 ) then ! on the boundary => add degeneracy
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
    allocate( irvec_q(3,nrpts_q) )
    allocate( ndegen_q(nrpts_q) )
    ndegen_q = ndegen_ws(1:nrpts_q)
    do i=1,3
        irvec_q(i,:) = irvec_ws(i,1:nrpts_q)
    end do
    !
    print *, '#max nrpts_q = ',  nrpts_q

  end subroutine allocate_and_build_ws_irvec_q


  subroutine allocate_and_build_ws_irvec_qtau() !NEW VERSION
    !  Calculate real-space Wigner-Seitz lattice vectors for phonon grid.
    !  Similar to w90_setup allocate_and_build_ws_irvec routine.
    !  Unlike the allocate_and_build_ws_irvec_q version, here we choose as WS
    !  criterion that, for each tau1,tau2 atom pair, the WS cell is centered at tau1.
    !  So, we impose the truncation using R+tau2-tau1,
    !  as in S. Poncé et al, Phys. Rev. Research 3, 043022 (2021)  (appendix D)
    !  and  G. Pizzi et al, J. Phys.: Condens. Matter 32 165902 (2020) (section 4.2).
    !----------------------------------------------------------------------------!
    use intw_reading, only: at, alat, nat, tau_cryst
    use intw_useful_constants, only: eps_8
    use intw_utility, only: cryst_to_cart
    use intw_input_parameters, only: nq1, nq2, nq3
    use intw_ph, only: qmesh
    !
    implicit none
    !
    integer :: iq, i,j,k,l, l0,l1, nboundary, iat1, iat2, nws
    logical :: in_ws
    integer :: Rs(3,n_wss_q), r_cryst_int(3), ndegen_ws(nq1*nq2*nq3*n_wss_q,nat,nat), irvec_ws(3,nq1*nq2*nq3*n_wss_q,nat,nat)
    real(kind=dp) :: r_cryst(3), r_length_l, r_length_l1, r_cart(3)

    irvec_ws = 0
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
    nrpts_qtau12 = 0
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
            call cryst_to_cart(1, r_cart, at, 1)
            r_length_l = alat * sqrt( sum(r_cart*r_cart) )  ! distance of r-R(l) to O  (cartesian, bohr)
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
              call cryst_to_cart(1, r_cart, at, 1)
              r_length_l1 =  alat * sqrt( sum(r_cart*r_cart) )  ! distance of r-R(l) to O' (cartesian, bohr)
              ! compare distances leaving a gap eps_8*1000.
              ! TODO !!! put tolerance as parameter. It depends a lot on this!
              ! I guess that we need a smaller one the less nq...
              if ( r_length_l > r_length_l1 + eps_8*1000. .and. l1/=l0 ) then ! not in the WS => remove vector from list
                in_ws = .false.
                exit
              else if ( abs(r_length_l-r_length_l1)<=eps_8*1000. .and. l1/=l0 ) then ! on the boundary => add degeneracy
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
    allocate( irvec_qtau(3,nrpts_qtau,nat,nat) )
    allocate( ndegen_qtau(nrpts_qtau,nat,nat) )
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
    print*, '#max nrpts_q = ',  nrpts_qtau12

  end subroutine allocate_and_build_ws_irvec_qtau


  subroutine allocate_and_build_dyn_qmesh()
    ! Read the dynamical matrices and compute eiegnvectors and omega^2 on the full qmesh
    !----------------------------------------------------------------------------
    use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
    use intw_reading, only: nat
    use intw_ph, only: nqmesh, qmesh, read_dynq, read_dynq_sym
    !
    implicit none
    !
    integer :: iq
    real(dp) :: qpoint(3)
    !
    !
    ! allocate and calculate dynmat at qmesh
    allocate( dyn_q(3*nat,3*nat,nqmesh), w2_q(3*nat,nqmesh), u_q(3*nat,3*nat,nqmesh) )
    !
    ! read dynmat
    call read_dynq_sym(dyn_q)
    !
    ! diagonalize
    do iq=1,nqmesh
      qpoint = qmesh(:,iq)
      call dyn_diagonalize_1q(3*nat, dyn_q(:,:,iq), u_q(:,:,iq), w2_q(:,iq))
    end do

  end subroutine allocate_and_build_dyn_qmesh


  subroutine allocate_and_build_dyn_qmesh_from_fc(fcfile)
    ! Same as but read force constants and then compute dynamical matrices
    !----------------------------------------------------------------------------
    use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
    use intw_input_parameters, only: nq1, nq2, nq3
    use intw_reading, only: nat
    use intw_ph, only: nqmesh, qmesh, mat_inv_four_t, readfc
    !
    implicit none
    !
    character(256) , intent(in) :: fcfile
    integer :: iq
    real(dp) :: qpoint(3)
    complex(dp) :: frc(nq1,nq2,nq3,3,3,nat,nat)
    !
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

  end subroutine allocate_and_build_dyn_qmesh_from_fc


  subroutine dyn_q_to_dyn_r()
    ! allocate_and_build_ws_irvec_qtau and allocate_and_build_dyn_qmesh must have been previously run
    !----------------------------------------------------------------------------
    use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
    use intw_reading, only: nat
    use intw_ph, only: nqmesh, qmesh
    !
    implicit none
    !
    integer :: ir, iq, iat1, iat2
    complex(kind=dp) :: fac
    !
    allocate (dyn_r(3*nat,3*nat,nrpts_qtau))
    !
    !
    dyn_r = cmplx_0
    !
    do iat1=1,nat
    do iat2=1,nat
      do ir = 1,nrpts_qtau12(iat1,iat2)
          do iq = 1, nqmesh
            fac = exp(-cmplx_i*tpi*dot_product(qmesh(:,iq), irvec_qtau(:,ir,iat1,iat2)))
            dyn_r((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, ir ) = &
                    dyn_r((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, ir ) + &
                    fac*dyn_q((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, iq)
          end do
      end do
    end do
    end do
    dyn_r = dyn_r / real(nqmesh,dp)

  end subroutine dyn_q_to_dyn_r


  subroutine dyn_interp_1q(qpoint, dyn_qint)
    ! allocate_and_build_ws_irvec_qtau and dyn_q_to_dyn_r must have been previously run,
    ! as this uses irvec_qtau and dyn_r variables.
    !----------------------------------------------------------------------------
    use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
    use intw_reading, only: nat
    !
    implicit none
    !
    real(dp) , intent(in) :: qpoint(3)
    !
    complex(dp) , intent(out) :: dyn_qint(3*nat,3*nat)
    !
    integer :: ir, iat1, iat2
    complex(kind=dp) :: fac
    !
    !
    dyn_qint = cmplx_0
    do iat1 = 1,nat
      do iat2 = 1,nat
        do ir = 1,nrpts_qtau12(iat1,iat2)
          fac = exp(cmplx_i*tpi*dot_product(qpoint(:), irvec_qtau(:,ir,iat1,iat2)))/real(ndegen_qtau(ir,iat1,iat2),dp)
          dyn_qint((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3) = &
              dyn_qint((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3) + &
              fac * dyn_r((iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3,ir)
        end do
      end do
    end do

  end subroutine dyn_interp_1q


  subroutine dyn_diagonalize_1q(n, dynq, uq, w2q)
    ! For a given dynq at a certain qpoint, with n=3*nat,
    ! this subroutine is a driver for zhpevx, which returns omega^2 and eigenvectors
    !
    ! TODO this is c+p Haritz's diagonalize_cmat, which should be in intw_utilities, since it is useful
    !----------------------------------------------------------------------------
    use intw_reading, only: nat, amass, ityp
    use intw_utility, only: diagonalize_cmat
    use intw_useful_constants, only: pmass
    !
    implicit none
    !
    integer, intent(in)  :: n
    complex(dp), intent(in) :: dynq(n,n) ! Dynamical matrix in a.u. (without the mass factor)
    complex(dp), intent(out) :: uq(n,n) ! Phonon polarization vectors in a.u.
    real(dp), intent(out) :: w2q(n) ! Phonon frequencies^2 in a.u.
    !
    integer :: iat1, iat2
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
