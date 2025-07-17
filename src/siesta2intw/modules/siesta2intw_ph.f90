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
module siesta2intw_ph

  use precision, only: dp

  implicit none

  public :: allocate_vloc, deallocate_vloc, phq_init, &
            calculate_local_part_dv

  private

  real(kind=dp), allocatable, dimension(:,:) :: vlocq


contains

  subroutine calculate_local_part_dv(tpiba, bg, qcryst, imode, dv_local)
    !
    ! We have dV_scf as input and we add to it the derivative of the PP
    !

    ! variables
    use siesta_geom, only: isa, na_u, xa
    use m_ntm, only: ntm
    use siesta2intw_fft, only: nG, gvec_cart, gamma_only, nl, nlm
    use siesta2intw_utils, only: cmplx_0, cmplx_i

    implicit none

    ! I/O variables
    real(kind=dp), intent(in) :: tpiba, bg(3, 3), qcryst(3)
    integer, intent(in) :: imode
    complex(dp), intent(out) :: dv_local(ntm(1)*ntm(2)*ntm(3))

    ! local variables
    real(kind=dp) :: qcart(3)
    integer :: ia, id, nt, ig
    complex(dp) :: qtau, gtau
    external :: cfftnd

    !
    ia = (imode-1)/3+1
    if (ia >na_u) stop "ERROR: calculate_local_part_dv: ia is larger than na_u"
    id = modulo(imode-1, 3)+1
    nt = isa(ia)
    !
    dv_local = cmplx_0
    !
    qcart = tpiba * matmul(bg, qcryst)
    qtau = exp(-cmplx_i*sum(qcart(:)*xa(:,ia)))
    do ig=1,nG
      gtau = exp(-cmplx_i*sum(gvec_cart(:,ig)*xa(:,ia)))
      dv_local(nl(ig)) = dv_local(nl(ig)) -cmplx_i*(qcart(id)+gvec_cart(id, ig))*qtau*gtau*vlocq(ig, nt)
      if (gamma_only) dv_local(nlm(ig)) = dv_local(nlm(ig)) -conjg(cmplx_i*(qcart(id)+gvec_cart(id, ig))*qtau*gtau*vlocq(ig, nt))
    enddo !ig
    !
    !
    call cfftnd(3, (/ntm(1), ntm(2), ntm(3)/), 1, dv_local)
    !

  end subroutine calculate_local_part_dv


  subroutine allocate_vloc()

    use atm_types, only: nspecies
    use siesta2intw_utils, only: cmplx_0
    use siesta2intw_fft, only: nG

    allocate(vlocq(nG, nspecies))
    vlocq = cmplx_0

  end subroutine allocate_vloc


  subroutine deallocate_vloc()

    deallocate(vlocq)

  end subroutine deallocate_vloc


  subroutine phq_init(tpiba, bg, qcryst)
    !
    ! This subroutine is based on the phq_init subroutine distributed as part of
    ! the Quantum Espresso project:
    !   Copyright (C) 2001-2008 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !

    ! variables
    use atm_types, only: nspecies
    use siesta2intw_pp, only: intwPPs

    implicit none

    ! I/O variables
    real(kind=dp), intent(in) :: tpiba, bg(3, 3), qcryst(3)

    ! local variables
    real(kind=dp) :: qcart(3)
    integer :: nt

    !
    qcart = tpiba * matmul(bg, qcryst)
    !
    DO nt = 1, nspecies
      CALL setlocq( qcart, intwPPs(nt)%mesh, intwPPs(nt)%r, intwPPs(nt)%rab, &
                    intwPPs(nt)%vpsloc, intwPPs(nt)%Zval, &
                    vlocq(:, nt) )
    END DO

  END SUBROUTINE phq_init


  subroutine setlocq(q_cart, mesh, r, rab, vloc_r, Zval, vloc_g)
    !
    !    This routine computes the Fourier transform of the local
    !    part of the pseudopotential in the q+G vectors.
    !    In order to perform the Fourier transform, a term erf(r)/r
    !    is subtracted in real space and added again in G space
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2001-2009 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !

    ! variables
    USE units, ONLY: pi
    use siesta_geom, only: volume_of_some_cell
    use siesta2intw_utils, only: cmplx_0
    use siesta2intw_fft, only: nG, gvec_cart
    ! functions and subroutines
    use siesta2intw_utils, only: simpson

    implicit none

    real(kind=dp), intent(in) :: q_cart(3) ! the q point in Cartesian units
    integer, intent(in) :: mesh ! the dimensions of the mesh
    real(kind=dp), intent(in) :: r(mesh) ! the mesh points
    real(kind=dp), intent(in) :: rab(mesh) ! the derivative of mesh points
    real(kind=dp), intent(in) :: vloc_r(mesh) ! the pseudo on the radial
    real(kind=dp), intent(in) :: Zval ! valence pseudocharge
    real(kind=dp), intent(out) :: vloc_g(nG) ! the fourier transform of the potential

    real(kind=dp), parameter :: eps = 1.d-8
    !
    real(kind=dp) :: vlcp, vloc0, g2a, aux(mesh), aux1(mesh), gx
    integer :: ig, ir

    !
    vloc_g = cmplx_0

    !
    ! the G == 0 term
    aux(:) = r(:) * (r(:) * vloc_r(:) + 2.0_dp*Zval)
    call simpson(mesh, rab, aux, vloc0)
    !
    ! the G /= 0 terms, we first compute the part of the integrand func
    ! indipendent of |G| in real space
    do ir = 1, mesh
      aux1(ir) = r(ir) * vloc_r(ir) + 2.0_dp*Zval * erf(r(ir))
    enddo
    !
    ! perform the integral, after multiplying for the |G| dependent part
    !
    do ig = 1, nG
      g2a = sum( (q_cart(:)+gvec_cart(:, ig))**2 )
      if (g2a < eps) then
          vloc_g(ig) = vloc0
      else
          gx = sqrt( g2a )
          do ir = 1, mesh
            aux(ir) = aux1(ir) * sin( gx*r(ir) ) / gx
          enddo
          call simpson(mesh, rab, aux, vlcp)
          !
          ! add the analytic fourier transform of the erf function
          vloc_g(ig) = vlcp - 2.0_dp*Zval * exp( -g2a*0.25d0 ) / g2a
      endif
    enddo

    vloc_g = 4.0_dp*pi * vloc_g / volume_of_some_cell

  end subroutine setlocq

end module siesta2intw_ph