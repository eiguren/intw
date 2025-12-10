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
module intw_pseudo_local

  !! display: none
  !!
  !! This module contains variables and subroutines for obtaining the local
  !! part of the pseudo-potentials.
  !!

  use kinds, only: dp

  implicit none

  ! variables
  public :: vloc

  ! subroutines
  public :: init_local_PP, deallocate_local_PP, init_vlocq, &
            calculate_local_part_v, calculate_local_part_dv, &
            dvqpsi_local

  private


  real(kind=dp), allocatable :: vloc(:,:) ! Local potential in reciprocal space for each atomic type


contains

  subroutine init_local_PP()
    !

    use intw_reading, only: ntyp, nG

    implicit none

    real(kind=dp) :: q_cryst(3)

    if (.not. allocated(vloc)) allocate(vloc(nG, ntyp))
    !
    q_cryst = (/0.0_dp, 0.0_dp, 0.0_dp/)
    !
    call init_vlocq(q_cryst, vloc)

  end subroutine init_local_PP


  subroutine deallocate_local_PP()

    implicit none

    deallocate(vloc)

  end subroutine deallocate_local_PP


  subroutine init_vlocq(q_cryst, vlocq)
    !
    ! This subroutine is based on the phq_init subroutine distributed as part of
    ! the Quantum Espresso project:
    !   Copyright (C) 2001-2008 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use intw_reading, only: ntyp, tpiba2, nG, bg, volume0
    use intw_fft, only: gvec_cart
    use intw_pseudo, only: upf
    !
    implicit none
    !
    real(kind=dp), intent(in) :: q_cryst(3)
    real(kind=dp), intent(out) :: vlocq(nG, ntyp)
    !
    ! local variables
    real(kind=dp) :: q_cart(3)
    integer :: nt


    q_cart = matmul(bg, q_cryst)
    !
    do nt = 1, ntyp
      call setlocq( q_cart, upf(nt)%mesh, upf(nt)%mesh, upf(nt)%rab, upf(nt)%r,&
                    upf(nt)%vloc, upf(nt)%zp, tpiba2, nG, gvec_cart, volume0, &
                    vlocq(:,nt) )
    end do

  end subroutine init_vlocq


  subroutine calculate_local_part_v(v_local)
    !======================================================================
    ! Add the local part of the PP (V_loc) to v_local                     !
    !======================================================================

    use intw_reading, only: nr1, nr2, nr3, nG, ntyp
    use intw_fft, only: nl, strf
    use intw_useful_constants, only: cmplx_0


    implicit none

    external :: cfftnd

    !I/O variables

    complex(kind=dp), intent(inout) :: v_local(nr1*nr2*nr3)

    !local variables
    integer :: nt, ig
    complex(kind=dp) :: aux(nr1*nr2*nr3)


    aux = cmplx_0
    do nt = 1, ntyp
      !
      do ig = 1, nG
        !
        aux(nl(ig)) = aux(nl(ig)) + strf(ig,nt) * vloc(ig,nt)
        !
      enddo ! ig
      !
    end do ! nt
    !
    call cfftnd(3, (/nr1,nr2,nr3/), 1, aux)
    !
    v_local = v_local + aux

  end subroutine calculate_local_part_v


  subroutine calculate_local_part_dv(q_cryst, dvq_local)
    !======================================================================
    ! We have dV_scf as input and we add to it the derivative of the PP   !
    !======================================================================

    use intw_reading, only: nat, nspin, nr1, nr2, nr3, nG, tpiba, ityp, bg, tau, ntyp
    use intw_fft, only: nl, gvec_cart, phase
    use intw_useful_constants, only: cmplx_i, cmplx_0, tpi

    implicit none

    external :: cfftnd

    !I/O variables

    real(kind=dp), intent(in) :: q_cryst(3)
    complex(kind=dp), intent(inout) :: dvq_local(nr1*nr2*nr3,3*nat,nspin,nspin)

    !local variables
    complex(kind=dp) :: eigqts(nat)
    integer :: imode, na, nt, ipol, ig, ispin, ir
    complex(kind=dp) :: aux(nr1*nr2*nr3)
    real(kind=dp) :: q_cart(3), vlocq(nG,ntyp)


    q_cart = matmul(bg, q_cryst)

    ! Compute local potential of the KB PP for q
    call init_vlocq(q_cryst, vlocq)

    do na = 1, nat
      !
      eigqts(na) = exp(-cmplx_i*tpi*dot_product(q_cart, tau(:,na)))
      !
    end do
    !
    do imode = 1, 3*nat
      !
      na = (imode-1)/3 + 1
      ipol = modulo(imode-1, 3) + 1
      nt = ityp(na)
      !
      aux = cmplx_0
      !
      do ig = 1, nG
        !
        aux(nl(ig)) = aux(nl(ig)) - cmplx_i * tpiba * ( q_cart(ipol)+gvec_cart(ipol,ig) ) * &
                                    eigqts(na) * phase(ig,na) * vlocq(ig,nt)
        !
      enddo !ig
      !
      call cfftnd(3, (/nr1,nr2,nr3/), 1, aux)
      !
      do ispin = 1, nspin
        do ir = 1, nr1*nr2*nr3
          !
          dvq_local(ir,imode,ispin,ispin) = dvq_local(ir,imode,ispin,ispin) + aux(ir)
          !
        enddo !ir
      enddo !ispin
      !
    enddo !imode

  end subroutine calculate_local_part_dv


  subroutine dvqpsi_local(nbands, list_iGk, list_iGkq, wfc_k, dvq_local, dvpsi_local)

    use intw_useful_constants, only: cmplx_0
    use intw_reading, only: nat, nspin, nGk_max, nr1, nr2, nr3
    use intw_fft, only: nl

    implicit none

    external :: cfftnd

    !I/O variables

    integer, intent(in) :: nbands, list_iGk(nGk_max), list_iGkq(nGk_max)
    complex(kind=dp), intent(in) :: dvq_local(nr1*nr2*nr3,3*nat,nspin,nspin), wfc_k(nGk_max,nbands,nspin)
    complex(kind=dp), intent(out) :: dvpsi_local(nGk_max,nbands,nspin,nspin,3*nat)

    !local variables

    integer :: ibnd, ispin, jspin, ig, imode, ir
    complex(kind=dp) :: wfc_r(nr1*nr2*nr3,nspin,nspin), wfc_r1(nr1*nr2*nr3,nspin)


    dvpsi_local = cmplx_0
    !
    do imode = 1, 3*nat
      !
      do ibnd= 1, nbands
        !
        ! Fourier transform the wave function to real space
        wfc_r1 = cmplx_0
        do ispin = 1, nspin
          do ig = 1, nGk_max
            !
            if (list_iGk(ig)==0) exit
            wfc_r1(nl(list_iGk(ig)),ispin) = wfc_k(ig,ibnd,ispin)
            !
          enddo !ig
          !
          call cfftnd(3, (/nr1,nr2,nr3/), 1, wfc_r1(:,ispin))
          !
        enddo !ispin
        !
        ! Multiply the wave function and the potential in real space
        wfc_r = cmplx_0
        do ispin = 1, nspin
          do jspin = 1, nspin
            !
            do ir = 1, nr1*nr2*nr3
              !
              wfc_r(ir,ispin,jspin) = dvq_local(ir,imode,ispin,jspin) * wfc_r1(ir,jspin)
              !
            enddo !ir
            !
          enddo !jspin
        enddo !ispin
        !
        !
        do ispin = 1, nspin
          do jspin = 1, nspin
            !
            call cfftnd(3, (/nr1,nr2,nr3/), -1, wfc_r(:,ispin,jspin))
            !
            do ig = 1, nGk_max
              !
              if (list_iGkq(ig)==0) exit
              !
              dvpsi_local(ig,ibnd,ispin,jspin,imode) = dvpsi_local(ig,ibnd,ispin,jspin,imode) &
                  + wfc_r(nl(list_iGkq(ig)),ispin,jspin)
              !
            enddo !ig
            !
          enddo !jspin
        enddo !ispin
        !
      enddo !ibnd
      !
    enddo !imode

  end subroutine dvqpsi_local


  subroutine setlocq(q_cart, mesh, msh, rab, r, vloc_at, zp, tpiba2, nG, g_cart, omega, vlocq)
    !----------------------------------------------------------------------
    !
    ! This routine computes the Fourier transform of the local
    ! part of the pseudopotential in the q+G vectors.
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2001-2009 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use kinds, only: dp
    use intw_useful_constants, only: TWO, fpi, pi
    use intw_utility, only: simpson, qe_erf, qe_erfc
    !
    implicit none
    !
    ! I/O variables
    !
    integer, intent(in) :: nG, mesh, msh
    ! input: the number of G vectors
    ! input: the dimensions of the mesh
    ! input: mesh points for radial integration

    real(kind=dp), intent(in) :: q_cart(3), zp, rab(mesh), r(mesh), vloc_at(mesh), tpiba2, omega, g_cart(3,nG)
    ! input: the q point
    ! input: valence pseudocharge
    ! input: the derivative of mesh points
    ! input: the mesh points
    ! input: the pseudo on the radial
    ! input: 2 pi / alat
    ! input: the volume of the unit cell
    ! input: the g_cart vectors
    real(kind=dp), intent(out) :: vlocq(nG)
    ! output: the fourier transform of the potential
    !
    ! local variables
    !
    real(kind=dp), parameter :: eps = 1.d-8
    real(kind=dp) :: vlcp, vloc0, fac, g2a, aux(mesh), aux1(mesh), gx
    ! auxiliary variables
    ! gx = modulus of g_cart vectors
    integer :: ig, ir
    ! counters
    !
    ! Pseudopotentials in numerical form (Vnl(lloc) contain the local part)
    ! in order to perform the Fourier transform, a term erf(r)/r is
    ! subtracted in real space and added again in G space
    !
    ! first the G=0 term
    !
    do ir = 1, msh
      aux(ir) = r(ir) * (r(ir) * vloc_at(ir) + TWO * zp)
    enddo
    call simpson(msh, aux, rab, vloc0)
    !
    !   here the G<>0 terms, we first compute the part of the integrand func
    !   indipendent of |G| in real space
    !
    do ir = 1, msh
      aux1(ir) = r(ir) * vloc_at(ir) + TWO * zp * qe_erf(r(ir))
    enddo
    fac = TWO * zp / tpiba2
    !
    !    and here we perform the integral, after multiplying for the |G|
    !    dependent  part
    !
    do ig = 1, nG
      g2a = (q_cart(1) + g_cart(1,ig))**2 + (q_cart(2) + g_cart(2,ig))**2 + (q_cart(3) + g_cart(3,ig))**2
      if (g2a < eps) then
        vlocq(ig) = vloc0
      else
        gx = sqrt(g2a * tpiba2)
        do ir = 1, msh
          aux(ir) = aux1(ir) * sin(gx * r(ir)) / gx
        enddo
        call simpson(msh, aux, rab, vlcp)
        !
        !     here we add the analytic fourier transform of the erf function
        !
        vlocq(ig) = vlcp - fac * exp( - g2a * tpiba2 * 0.25d0) / g2a
      endif
    enddo

    vlocq(:) = vlocq(:) * fpi / omega

  end subroutine setlocq


  subroutine setlocq_coul(q_cart, zp, tpiba2, nG, g_cart, omega, vlocq)
    !----------------------------------------------------------------------
    !
    ! Fourier transform of the Coulomb potential - For all-electron
    ! calculations, in specific cases only, for testing purposes
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2001-2009 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use kinds, only: dp
    use intw_useful_constants, only: fpi, TWO, eps_8
    implicit none
    !
    integer, intent(in) :: nG
    real(kind=dp), intent(in) :: q_cart(3), zp, tpiba2, omega, g_cart(3,nG)
    real(kind=dp), intent(out) :: vlocq(nG)
    !
    real(kind=dp) :: g2a
    integer :: ig

    do ig = 1, nG
      g2a = (q_cart(1) + g_cart(1,ig))**2 + (q_cart(2) + g_cart(2,ig))**2 + (q_cart(3) + g_cart(3,ig))**2
      if (g2a < eps_8) then
          vlocq(ig) = 0.d0
      else
          vlocq(ig) = - fpi * TWO * zp / omega / tpiba2 / g2a
      endif
    enddo

  end subroutine setlocq_coul

end module intw_pseudo_local
