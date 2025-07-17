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
module intw_pseudo_non_local

  !------------------------------------------------------------------------!
  ! This module contains variables and subroutines to obtain the non-local !
  ! part of the pseudo-potentials.                                         !
  !------------------------------------------------------------------------!

  use kinds, only: dp

  implicit none

  ! variables
  public :: nkb, DKB

  ! subroutines
  public :: init_KB_PP, init_KB_projectors, &
            multiply_psi_by_vKB, multiply_psi_by_dvKB

  private


  integer, parameter :: l_kb_max = 3 ! Max non local angular momentum (l=0 to l_kb_max)

  real(kind=dp), parameter   :: dq = 0.01d0    ! Space between points in the pseudopotential tab
  integer                    :: nqx            ! Number of interpolation points
  real(kind=dp), allocatable :: tab(:,:,:)     ! Interpolation table for PPs
  real(kind=dp), allocatable :: tab_d2y(:,:,:) ! For cubic splines

  integer                       :: lmaxkb          ! Max angular momentum of beta functions
  integer                       :: nbetam          ! Max number of beta functions per atomic type
  integer,          allocatable :: nh(:)           ! Number of beta(lm) functions per atomic type
  integer                       :: nhm             ! Max number of beta(lm) functions per atomic type
  integer,          allocatable :: nhtonbeta(:,:)  ! Link between index of beta(lm) function in the atomic type -> index of beta function in the atomic type
  integer,          allocatable :: nhtol(:,:)      ! Link between index of beta(lm) function in the atomic type -> angular momentum l
  integer,          allocatable :: nhtolm(:,:)     ! Link between index of beta(lm) function in the atomic type -> combined lm angular momentum index l*l+m
  real(kind=dp),    allocatable :: nhtoj(:,:)      ! Link between index of beta(lm) function in the atomic type -> total angular momentum j
  complex(kind=dp), allocatable :: Dion(:,:,:,:,:) ! D_{mu,nu} matrix for beta(lm) functions for each atomic type

  integer                       :: nkb          ! Total number of beta(lm) functions in the solid
  integer,          allocatable :: nkbtona(:)   ! Link between index of beta(lm) function in the solid -> index of the atom
  complex(kind=dp), allocatable :: DKB(:,:,:,:) ! D_{mu,nu} matrix for beta(lm) functions in the solid


contains

  subroutine init_KB_PP()
    !
    ! This subroutine is based on the init_us_1 subroutine distributed as part of
    ! the Quantum Espresso project:
    !   Copyright (C) 2001-2007 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    !  Modifications by INTW group, 2024:
    !   - Split the subroutine in small parts that make certain tasks.
    !   - Remove parts related to ultrasoft pseudo potentials.
    !   - Calculate Dij matrix for the KB projectors of the solid.
    !
    use intw_reading, only: nat, ntyp, ityp
    use intw_pseudo, only: upf
    use intw_useful_constants, only: cmplx_0

    implicit none

    integer :: nt, nb, na


    !
    ! Calculate the total number of beta(lm) projectors for each atomic type and maximum angular momentum
    !
    if (allocated(nh)) deallocate(nh)
    allocate(nh(ntyp))
    lmaxkb = - 1
    do nt = 1, ntyp
      !
      nh(nt) = 0
      !
      do nb = 1, upf(nt)%nbeta
        nh(nt) = nh(nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = max(lmaxkb, upf(nt)%lll(nb))
      end do
      !
    end do
    !
    ! Calculate the maximum number of beta(lm) and beta functions
    !
    nhm = maxval(nh)
    nbetam = maxval(upf(:)%nbeta)
    !
    ! Calculate the number of beta(lm) functions in the solid
    !
    nkb = 0
    do na = 1, nat
      nt = ityp(na)
      nkb = nkb + nh(nt)
    end do
    !
    ! Calculate the arrays to link indices of the beta functions
    !
    call init_KB_link_indices()
    !
    ! Calculate Dij matrix for each atomic type
    !
    call init_Dion()
    !
    ! Calculate Dij matrix for the solid
    !
    call init_DKB()
    !
    ! Calculate interpolation table
    !
    call init_interpolation_table()

  end subroutine init_KB_PP


  subroutine init_KB_link_indices()
    ! Initialize the arrays used to link the index of the beta(lm) projector for
    ! each atomic type with the index of the projetor of the PP file (nhtonbeta),
    ! with the angular momentum of the projetor (nhtol, nhtolm and nhtoj), and the
    ! index of the beta(lm) projector in the solid with the index of the atom (nkbtona).
    !
    ! This subroutine is based on the init_us_1 subroutine distributed as part of
    ! the Quantum Espresso project:
    !   Copyright (C) 2001-2007 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use intw_reading, only: ntyp, nat, ityp
    use intw_pseudo, only: upf

    implicit none

    !local variables
    integer :: nt, ih, nb, l, m, na, ikb
    real(kind=dp) :: j


    if(allocated(nhtonbeta)) deallocate(nhtonbeta)
    allocate(nhtonbeta(nhm,ntyp))
    !
    if(allocated(nhtol)) deallocate(nhtol)
    allocate(nhtol(nhm,ntyp))
    !
    if(allocated(nhtolm)) deallocate(nhtolm)
    allocate(nhtolm(nhm,ntyp))
    !
    if(allocated(nhtoj)) deallocate(nhtoj)
    allocate(nhtoj(nhm,ntyp))
    !
    do nt = 1, ntyp
      !
      ih = 1
      do nb = 1, upf(nt)%nbeta
        !
        l = upf(nt)%lll(nb)
        !
        do m = 1, 2 * l + 1
          !
          nhtol(ih,nt) = l
          nhtolm(ih,nt) = l*l+m
          nhtonbeta(ih,nt) = nb
          ih = ih + 1
          !
        end do !m
        !
      end do !nb
      !
      if ( upf(nt)%has_so ) then
        !
        ih = 1
        do nb = 1, upf(nt)%nbeta
          !
          l = upf(nt)%lll(nb)
          j = upf(nt)%jjj(nb)
          !
          do m = 1, 2 * l + 1
            !
            nhtoj(ih, nt) = j
            ih = ih + 1
            !
          end do !m
          !
        end do !nb
        !
      end if
      !
    end do !nt

    if(allocated(nkbtona)) deallocate(nkbtona)
    allocate(nkbtona(nkb))
    ikb = 0
    do na = 1, nat
      nt = ityp(na)
      do ih = 1, nh(nt)
        ikb = ikb + 1
        nkbtona(ikb) = na
      end do !ih
    end do !na

  end subroutine init_KB_link_indices


  subroutine init_Dion()
    ! Initialize the Dion matrix for the beta(lm) projetors for all atomic types
    ! Dion includes the spin indices in the case of a SO calculation
    !
    ! This subroutine is based on the init_us_1 subroutine distributed as part of
    ! the Quantum Espresso project:
    !   Copyright (C) 2001-2007 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use intw_useful_constants, only: sqrt2, cmplx_0, cmplx_1, cmplx_i
    use intw_reading, only: ntyp, lspinorb, nspin
    use intw_pseudo, only: upf

    implicit none

    complex(kind=dp) :: rot_ylm(2*l_kb_max+1,2*l_kb_max+1) ! to transform real spherical harmonics into complex ones
    complex(kind=dp), allocatable :: fcoef(:,:,:,:,:) ! needed to account for spinors
    integer :: n1, n, l, m
    integer :: nt, ih, jh, kh, ir, is
    integer :: m0, m1, li, lk, mi, mk, vi, vj, ispin, jspin
    complex(kind=dp) :: coeff
    real(kind=dp) :: ji, jk


    !
    ! Initialization of the variables
    !
    if (lspinorb) then
      allocate(Dion(nhm,nhm,nspin,nspin,ntyp))
      allocate(fcoef(nhm,nhm,2,2,ntyp))
    else
      allocate(Dion(nhm,nhm,1,1,ntyp))
    end if
    !
    Dion = cmplx_0
    !
    ! the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
    ! but in some versions of the PP files lmax is not set to the maximum
    ! l of the beta functions but includes the l of the local potential
    !
    if (lspinorb) then
      !
      ! In the spin-orbit case we need the unitary matrix u which rotates the
      ! real spherical harmonics and yields the complex ones.
      !
      rot_ylm = cmplx_0
      l = l_kb_max
      rot_ylm(l+1, 1) = cmplx_1
      do n1 = 2, 2*l+1, 2
        m = n1/2
        n = l + 1 - m
        rot_ylm(n, n1) = cmplx_1/sqrt2 * (-1)**m
        rot_ylm(n, n1+1) = -cmplx_i/sqrt2 * (-1)**m
        n = l + 1 + m
        rot_ylm(n, n1) = cmplx_1/sqrt2
        rot_ylm(n, n1+1) = cmplx_i/sqrt2
      end do
      fcoef = cmplx_0
    end if
    !
    ! For each pseudopotential we initialize the indices nhtol, nhtolm,
    ! nhtoj, nhtonbeta, and if the pseudopotential is of KB type we initialize the
    ! atomic D terms
    !
    !
    ! Here we initialize the D of the solid
    !
    do nt = 1, ntyp
      !
      if (upf(nt)%has_so .and. lspinorb) then
        !
        ! first calculate the fcoef coefficients
        !
        do ih = 1, nh(nt)
          li = nhtol(ih,nt)
          ji = nhtoj(ih,nt)
          mi = nhtolm(ih,nt) - li*li
          do kh = 1, nh(nt)
            lk = nhtol(kh,nt)
            jk = nhtoj(kh,nt)
            mk = nhtolm(kh,nt) - lk*lk
            if (li == lk .and. abs(ji-jk) < 1.d-7) then
              do ispin = 1, 2
                do jspin = 1, 2
                  coeff = cmplx_0
                  do m = -li-1, li
                    m0 = sph_ind(li,ji,m,ispin) + l_kb_max + 1
                    m1 = sph_ind(lk,jk,m,jspin) + l_kb_max + 1
                    coeff = coeff +         rot_ylm(m0,mi)  * spinor(li,ji,m,ispin) &
                                    * conjg(rot_ylm(m1,mk)) * spinor(lk,jk,m,jspin)
                  end do
                  fcoef(ih,kh,ispin,jspin,nt) = coeff
                end do
              end do
            end if
          end do
        end do
        !
        ! and calculate the bare coefficients
        !
        do ih = 1, nh(nt)
          vi = nhtonbeta(ih,nt)
          do jh = 1, nh(nt)
            vj = nhtonbeta(jh,nt)
            do ispin = 1, 2
              do jspin = 1, 2
                Dion(ih,jh,ispin,jspin,nt) = upf(nt)%dion(vi,vj) * fcoef(ih,jh,ispin,jspin,nt)
                if (vi.ne.vj) fcoef(ih,jh,ispin,jspin,nt) = cmplx_0
              end do
            end do
          end do
        end do
      else
        do ih = 1, nh(nt)
          do jh = 1, nh(nt)
            if ( nhtol(ih,nt) == nhtol(jh,nt) .and. nhtolm(ih,nt) == nhtolm(jh,nt) ) then
              ir = nhtonbeta(ih,nt)
              is = nhtonbeta(jh,nt)
              if (lspinorb) then
                Dion(ih,jh,1,1,nt) = upf(nt)%dion(ir,is)
                Dion(ih,jh,2,2,nt) = upf(nt)%dion(ir,is)
              else
                Dion(ih,jh,1,1,nt) = upf(nt)%dion(ir,is)
              end if
            end if
          end do
        end do
      end if
    end do

  end subroutine init_Dion


  subroutine init_DKB()
    ! Initilize the DKB matrix for all beta(lm) projectors in the solid

    use intw_useful_constants, only: cmplx_0
    use intw_reading, only: nat, ntyp, nspin, lspinorb, ityp

    implicit none

    integer :: ikb, jkb, nt, na, ih, jh, ntj, naj, ispin


    allocate(DKB(nkb,nkb,nspin,nspin))
    DKB = cmplx_0
    ikb = 0
    do nt = 1, ntyp
      do na = 1, nat
        !
        if (ityp(na)==nt) then
          !
          do ih = 1, nh(ityp(na))
            ikb = ikb + 1
            jkb = 0
            do ntj = 1, ntyp
              do naj = 1, nat
                if (ityp(naj)==ntj) then
                  do jh = 1, nh(ityp(naj))
                    jkb = jkb + 1
                    if (na==naj) then
                      if (lspinorb) then
                        DKB(ikb,jkb,:,:) = Dion(ih,jh,:,:,ityp(na))
                      else
                        do ispin=1,nspin
                          DKB(ikb,jkb,ispin,ispin) = Dion(ih,jh,1,1,ityp(na))
                        end do
                      end if
                    end if
                  end do !jh
                end if
              end do !naj
            end do !ntj
          end do !ih
          !
        end if
        !
      end do !na
    end do !nt

  end subroutine init_DKB


  subroutine init_interpolation_table()
    ! Initialize the interpoation table used to calculate the KB projectors in reciprocal space
    !
    ! This subroutine is based on the init_us_1 subroutine distributed as part of
    ! the Quantum Espresso project:
    !   Copyright (C) 2001-2007 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use intw_reading, only: ecutwfc, ntyp, volume0
    use intw_pseudo, only: upf
    use intw_useful_constants, only: fpi, ZERO
    use intw_utility, only: sphb, simpson, spline

    implicit none

    real(kind=dp), allocatable :: aux(:), xdata(:)
    real(kind=dp) ::  vqint
    integer :: iq, ndm, nt, nb, l


    !
    ! Calculate dimensions for array tab (including a possible factor
    ! coming from cell contraction during variable cell relaxation/MD)
    !
    nqx = int( sqrt(2*ecutwfc) / dq + 4 ) ! x2 zeren Ry -> Hartree egin behar
    !
    ! q-point grid for interpolation
    !
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    end do
    !
    ! Fill the interpolation table
    !
    allocate(tab(nqx,nbetam,ntyp))
    allocate(tab_d2y(nqx,nbetam,ntyp))
    !
    ndm = maxval(upf(:)%kkbeta)
    allocate(aux(ndm))
    !
    tab = ZERO
    do nt = 1, ntyp
      do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll(nb)
        do iq = 1, nqx

          aux = upf(nt)%beta(:,nb) * sphb(l, xdata(iq)*upf(nt)%r) * upf(nt)%r
          call simpson(upf(nt)%kkbeta, aux, upf(nt)%rab, vqint)

          ! ASIER 29/07/2021
          ! Integrating by spline + gauss 2. order
          ! vqint = intgr_spline_gaussq( upf(nt)%r(1:upf(nt)%kkbeta), aux )

          tab(iq,nb,nt) = vqint * fpi / sqrt(volume0)

        end do
      end do
    end do
    !
    deallocate(aux)
    !
    ! Initialize spline interpolation
    !
    do nt = 1, ntyp
      do nb = 1, upf(nt)%nbeta
        call spline(xdata, tab(:,nb,nt), 0.0_dp, 0.0_dp, tab_d2y(:,nb,nt))
      end do
    end do

  end subroutine init_interpolation_table


  subroutine init_KB_projectors(npw, igk, kpoint_cryst, kb_projectors)
    !----------------------------------------------------------------------
    !
    ! Calculates beta functions (Kleinman-Bylander projectors), with
    ! structure factor, for all atoms, in reciprocal space
    !
    ! This subroutine is based on the init_us_2 subroutine distributed as part of
    ! the Quantum Espresso project:
    !   Copyright (C) 2001 PWSCF group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use intw_reading, only: nat, nGk_max, ntyp, ityp, tau, bg, tpiba
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i
    use intw_fft, only: gvec_cart
    use intw_pseudo, only: upf
    use intw_utility, only: real_ylmr2, splint

    implicit none

    !I/O variables

    integer, intent(in) :: npw !number of PW's in k-point
    integer, intent(in) :: igk(nGk_max) !G list of k-point
    real(kind=dp), intent(in) :: kpoint_cryst(3) !k-point vector
    complex(kind=dp), intent(out) :: kb_projectors(nGk_max,nkb) !beta functions

    !local variables

    integer :: ig, l, lm, na, nt, nb, ih, jkb, iq
    real(kind=dp) :: kpoint_cart(3), kg_cart(3,npw), vkb1(npw,nhm), xdata(nqx), ylm(npw, (lmaxkb+1)**2)
    real(kind=dp), dimension(npw) :: kg(npw), vk(npw)
    complex(kind=dp) :: pref, sk(npw)


    kb_projectors = cmplx_0
    !
    kpoint_cart = matmul(bg, kpoint_cryst)
    !
    if (lmaxkb.lt.0) return
    !
    do ig = 1, npw
      !
      kg_cart(1,ig) = kpoint_cart(1) + gvec_cart(1, igk(ig))
      kg_cart(2,ig) = kpoint_cart(2) + gvec_cart(2, igk(ig))
      kg_cart(3,ig) = kpoint_cart(3) + gvec_cart(3, igk(ig))
      !
      kg(ig) = kg_cart(1,ig)**2 + kg_cart(2,ig)**2 + kg_cart(3,ig)**2
      !
    end do !ig
    !
    call real_ylmr2(lmaxkb, npw, kg_cart, kg, ylm)
    !
    ! set now kg=|k+G| in atomic units
    !
    kg(:) = tpiba * sqrt(kg(:))
    !
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    end do !iq
    !
    ! |beta_lm(kg)> = (4pi/omega) * (-i^l) * Y_lm(kg/|kg|) * f_l(|kg|) * S(kg)
    !
    jkb = 0
    vk = 0.d0
    !
    do nt = 1, ntyp
      !
      ! calculate beta in G-space using an interpolation table f_l(|kg|)=\int _0 ^\infty dr r^2 f_l(r) j_l(|kg|.r)
      !
      do nb = 1, upf(nt)%nbeta
        !
        do ig = 1, npw
          !
          vk(ig) = splint(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), kg(ig))
          !
        end do !ig
        !
        ! add spherical harmonic part: Y_lm(kg/|kg|) * f_l(|kg|)
        !
        do ih = 1, nh(nt)
          !
          if (nb /= nhtonbeta(ih, nt)) cycle
          !
          l  = nhtol(ih,nt)
          lm = nhtolm(ih,nt)
          !
          vkb1(:,ih) = ylm(:,lm) * vk(:)
          !
        end do !ih
        !
      end do !nb
      !
      ! vkb1 contains all betas including angular part for type nt
      ! now add the structure factor and factor (-i)^l
      !
      do na = 1, nat
        !
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        !
        if (ityp(na) /= nt) cycle
        !
        do ig = 1, npw
          !
          sk(ig) = exp(-tpi*cmplx_i*(  kg_cart(1,ig)*tau(1,na) &
                                     + kg_cart(2,ig)*tau(2,na) &
                                     + kg_cart(3,ig)*tau(3,na) ) )
          !
        end do !ig
        !
        do ih = 1, nh(nt)
          !
          jkb = jkb + 1
          pref = (-cmplx_i)**nhtol(ih,nt)
          !
          do ig = 1, npw
            !
            kb_projectors(ig, jkb) = pref * vkb1(ig,ih) * sk(ig)
            !
          end do !ig
          !
        end do !ih
        !
      end do !na
      !
    end do !ntyp

  end subroutine init_KB_projectors


  subroutine multiply_psi_by_vKB(k_cryst, list_iGk, num_bands, psi_k, vnl_psi)
    !INTW project: KB projection by wave functions.
    !

    use intw_useful_constants, only: cmplx_0
    use intw_reading, only: nspin, nGk_max

    implicit none

    real(kind=dp), intent(in)       :: k_cryst(3)
    integer, intent(in)             :: list_iGk(nGk_max)
    integer, intent(in)             :: num_bands
    complex(kind=dp), intent(in)    :: psi_k(nGk_max,num_bands,nspin)
    complex(kind=dp), intent(inout) :: vnl_psi(nGk_max,num_bands,nspin)

    complex(kind=dp)                :: vkb_k(nGk_max,nkb)
    complex(kind=dp)                :: projec_d(nkb,nspin), DKB_projec_d(nkb,nspin)
    integer                         :: iband, ikb, ispin, jspin
    integer                         :: iG, nGk


    ! Compute non local |beta> projectors of KB PP for k
    nGk = 0
    do iG=1,nGk_max
      if (list_iGk(iG)==0) exit
      nGk = nGk + 1
    enddo
    call init_KB_projectors(nGk, list_iGk, k_cryst, vkb_k)

    do iband = 1, num_bands
      !
      projec_d = cmplx_0
      !
      !
      do ikb = 1, nkb
        !
        do ispin = 1, nspin
          !
          do iG = 1, nGk
            projec_d(ikb,ispin) = projec_d(ikb,ispin) + conjg(vkb_k(iG,ikb)) * psi_k(iG,iband,ispin)
          end do !iG
          !
        end do !ispin
        !
      end do !ikb

      ! multiplay the projections <\beta_j|\psi_n> by the matrix DKB
      DKB_projec_d = cmplx_0
      do ispin = 1, nspin
        do jspin = 1, nspin
          !
          DKB_projec_d(:,ispin) = DKB_projec_d(:,ispin) + matmul( DKB(:,:,ispin,jspin), projec_d(:,jspin) )
          !
        end do !jspin
      end do !ispin

      !
      do ikb = 1, nkb
        !
        do ispin = 1, nspin
          !
          vnl_psi(:,iband,ispin) = DKB_projec_d(ikb,ispin) * vkb_k(:,ikb)
          !
        end do !ispin
        !
      end do !ikb
      !
    end do !iband

  end subroutine multiply_psi_by_vKB


  subroutine multiply_psi_by_dvKB(k_cryst, q_cryst, list_iGk, list_iGkq, num_bands, psi_k, dvnl_psi)

    use intw_useful_constants, only: cmplx_0, cmplx_i
    use intw_reading, only: nat, nspin, nGk_max, tpiba, bg
    use intw_fft, only: gvec_cart

    implicit none

    real(kind=dp), intent(in)       :: k_cryst(3), q_cryst(3)
    integer, intent(in)             :: list_iGk(nGk_max), list_iGkq(nGk_max)
    integer, intent(in)             :: num_bands
    complex(kind=dp), intent(in)    :: psi_k(nGk_max,num_bands,nspin)
    complex(kind=dp), intent(inout) :: dvnl_psi(nGk_max,num_bands,nspin,nspin,3*nat)

    complex(kind=dp)                :: vkb_k(nGk_max,nkb), vkb_kq(nGk_max,nkb)
    complex(kind=dp)                :: projec_1(nkb,3,nspin), projec_2(nkb,3,nspin)
    complex(kind=dp)                :: DKB_projec_1(nkb,3,nspin,nspin), DKB_projec_2(nkb,3,nspin,nspin)

    real(kind=dp)                   :: k_cart(3), q_cart(3)
    integer                         :: iG, iGk, iGkq, nGk, nGkq
    integer                         :: iband, ipol, ikb, ispin, jspin, na, imode


    k_cart = matmul(bg, k_cryst)
    q_cart = matmul(bg, q_cryst)

    ! Compute non local |beta> projectors of KB PP for k and k+q
    nGk = 0
    do iG=1,nGk_max
      if (list_iGk(iG)==0) exit
      nGk = nGk + 1
    enddo
    call init_KB_projectors(nGk, list_iGk, k_cryst, vkb_k)
    !
    nGkq = 0
    do iG=1,nGk_max
      if (list_iGkq(iG)==0) exit
      nGkq = nGkq + 1
    enddo
    call init_KB_projectors(nGkq, list_iGkq, k_cryst+q_cryst, vkb_kq)

    do iband = 1, num_bands

      projec_1 = cmplx_0
      projec_2 = cmplx_0

      ! Asier: KB potentziala hurrengo eran emanik dago: sum_l |b(l)> <b(l)|
      !               beraz deribatuak bi gai dauzka:
      !               sum_l d|b(l,r)> <b(l,r)| + |b(l,r)> d<b(l,r)| ~ Fourier ~
      !               sum_l i(k+G)|b(l,G)> <b(l,G)| + |b(l,G)> <b(l,G)|
      !

      do ipol = 1, 3 ! Cart. coord.
        !
        do ikb = 1, nkb
          !
          do ispin = 1, nspin
            !
            do iG = 1, nGk
              !
              iGk = list_iGk(iG)
              !
              projec_1(ikb,ipol,ispin) = projec_1(ikb,ipol,ispin) &
                  + conjg(vkb_k(iG,ikb)) * psi_k(iG,iband,ispin)
              !
              projec_2(ikb,ipol,ispin) = projec_2(ikb,ipol,ispin) &
                  + conjg(vkb_k(iG,ikb)) * psi_k(iG,iband,ispin) * &
                    tpiba * cmplx_i * ( k_cart(ipol) + gvec_cart(ipol,iGk) )
              !
            end do !iG
            !
          end do !ispin
          !
        end do !ikb
        !
      end do !ipol


      ! multiplay the projections <\beta_j|\psi_n> by the matrix DKB
      DKB_projec_1 = cmplx_0
      DKB_projec_2 = cmplx_0
      do ipol = 1, 3 ! Cart. coord.
        !
        do ispin = 1, nspin
          do jspin = 1, nspin
            !
            DKB_projec_1(:,ipol,ispin,jspin) = matmul( DKB(:,:,ispin,jspin), projec_1(:,ipol,jspin) )
            DKB_projec_2(:,ipol,ispin,jspin) = matmul( DKB(:,:,ispin,jspin), projec_2(:,ipol,jspin) )
            !
          end do !jspin
        end do !ispin
        !
      end do !ipol


      do ipol = 1, 3 ! Cart. coord.
        !
        do ikb = 1, nkb
          !
          na = nkbtona(ikb)
          !
          imode = (na-1)*3 + ipol
          !
          do iG = 1, nGkq
            !
            iGkq = list_iGkq(iG)
            !
            do ispin = 1, nspin
              do jspin = 1, nspin
                !
                dvnl_psi(iG,iband,ispin,jspin,imode) = dvnl_psi(iG,iband,ispin,jspin,imode) &
                    + DKB_projec_2(ikb,ipol,ispin,jspin) * vkb_kq(iG,ikb)
                !
                dvnl_psi(iG,iband,ispin,jspin,imode) = dvnl_psi(iG,iband,ispin,jspin,imode) &
                    - DKB_projec_1(ikb,ipol,ispin,jspin) * vkb_kq(iG,ikb) * &
                      tpiba * cmplx_i * ( k_cart(ipol) + q_cart(ipol) + gvec_cart(ipol,iGkq) )
                !
              end do !jspin
            end do !ispin
            !
          end do !iG
          !
        end do !ikb
        !
      end do !ipol

    end do !iband

  end subroutine multiply_psi_by_dvKB


  function sph_ind(l, j, m, spin)
    ! This function calculates the m index of the spherical harmonic
    ! in a spinor with orbital angular momentum l, total angular
    ! momentum j, projection along z of the total angular momentum m+-1/2.
    ! Spin selects the up (spin=1) or down (spin=2) coefficient.
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso project:
    !   Copyright (C) 2004 PWSCF group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    !  Modifications by INTW group, 2024:
    !   - Intent(in) variables specified explicitly.
    !
    use kinds, only: dp
    use intw_utility, ONLY : errore

    implicit none

    integer, intent(in)       :: l, & ! orbital angular momentum
                                 m, & ! projection of the total angular momentum+-1/2
                                 spin ! 1 or 2 select the component
    real(kind=dp), intent(in) :: j    ! total angular momentum
    integer :: sph_ind

    if (spin.ne.1.and.spin.ne.2) call errore('sph_ind','spin direction unknown',1)
    if (m.lt.-l-1.or.m.gt.l) call errore('sph_ind','m not allowed',1)

    if (abs(j-l-0.5d0).lt.1.d-8) then
      if (spin.eq.1) sph_ind= m
      if (spin.eq.2) sph_ind= m+1
    elseif (abs(j-l+0.5d0).lt.1.d-8) then
      if (m.lt.-l+1) then
          sph_ind=0
      else
          if (spin.eq.1) sph_ind= m-1
          if (spin.eq.2) sph_ind= m
      endif
    else
      write(6,*) l, j
      call errore('sph_ind','l and j not compatible',1)
    endif
    if (sph_ind.lt.-l.or.sph_ind.gt.l) sph_ind=0

    return
  end function sph_ind

  function spinor(l, j, m, spin)
    ! This function calculates the numerical coefficient of a spinor
    ! with orbital angular momentum l, total angular momentum j,
    ! projection along z of the total angular momentum m+-1/2. Spin selects
    ! the up (spin=1) or down (spin=2) coefficient.
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso project:
    !   Copyright (C) 2004 PWSCF group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    !  Modifications by INTW group, 2024:
    !   - Intent(in) variables specified explicitly.
    !
    use kinds, only: dp
    use intw_utility, ONLY : errore

    implicit none

    integer, intent(in)       :: l, & ! orbital angular momentum
                                 m, & ! projection of the total angular momentum+-1/2
                                 spin ! 1 or 2 select the component
    real(kind=dp), intent(in) :: j    ! total angular momentum
    real(kind=dp) :: spinor

    real(kind=dp) :: denom ! denominator

    if (spin.ne.1.and.spin.ne.2) call errore('spinor','spin direction unknown',1)
    if (m.lt.-l-1.or.m.gt.l) call errore('spinor','m not allowed',1)

    denom=1.d0/(2.d0*l+1.d0)
    if (abs(j-l-0.5d0).lt.1.d-8) then
      if (spin.eq.1) spinor= sqrt((l+m+1.d0)*denom)
      if (spin.eq.2) spinor= sqrt((l-m)*denom)
    elseif (abs(j-l+0.5d0).lt.1.d-8) then
      if (m.lt.-l+1) then
          spinor=0.d0
      else
          if (spin.eq.1) spinor= sqrt((l-m+1.d0)*denom)
          if (spin.eq.2) spinor= -sqrt((l+m)*denom)
      endif
    else
      call errore('spinor','j and l not compatible',1)
    endif

    return
  end function spinor

end module intw_pseudo_non_local