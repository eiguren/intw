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
module siesta2intw_utils

  use precision, only: dp

  implicit none

  public :: cmplx_0, cmplx_1, cmplx_i

  public :: atomic_magnetic_moment, ainv, simpson, dosplineint, spline, splint

  private

  complex(kind=dp), parameter :: cmplx_0 = (0.0_dp, 0.0_dp)
  complex(kind=dp), parameter :: cmplx_1 = (1.0_dp, 0.0_dp)
  complex(kind=dp), parameter :: cmplx_i = (0.0_dp, 1.0_dp)

contains

  subroutine atomic_magnetic_moment(Svec, Stot)
    !
    ! Compute the spin magnetic moment (vector and magnitude) for each atom
    !

    ! variables
    use precision, only: dp
    use siesta_geom, only: na_u
    use m_spin, only: spin
    use atomlist, only: no_l, iaorb
    use sparse_matrices, only: listhptr, numh
    use sparse_matrices, only: S, Dscf   ! Dscf could have 1, 2, 4, or 8 components
    use parallel, only: ionode, Node, Nodes
    use siesta2intw_io, only: stdout
    ! functions and subroutines
    use parallelsubs, only: LocalToGlobalOrb
#ifdef MPI
    use m_mpi_utils, only: globalize_sum
#endif

    implicit none

    real(kind=dp), intent(out) :: Svec(3, na_u) ! Spin magnetic moment vector in Cartesian coordinates
    real(kind=dp), intent(out) :: Stot(na_u) ! Spin magnetic moment magnitude

    real(kind=dp) :: qspin(spin%Grid, na_u) ! Spin density matrix
    real(kind=dp) :: qtot ! Total charge
#ifdef MPI
    real(kind=dp) :: qtmp(spin%Grid)
#endif
    integer :: io_l, io_g, ia, j, ind


    !
    ! Compute the charge density matrix for each atom
    qspin = 0.0_dp
    do io_l = 1, no_l
      call LocalToGlobalOrb(io_l, Node, Nodes, io_g)
      ia = iaorb(io_g)
      do j = 1, numh(io_l)
        ind = listhptr(io_l) + j
        if ( spin%SO ) then
          ! In the SOC case, hermitify Dscf
          qspin(1:2, ia) = qspin(1:2, ia) + Dscf(ind,1:2) * S(ind)
          qspin(3, ia) = qspin(3, ia) + 0.5_dp*(Dscf(ind,3)+Dscf(ind,7)) * S(ind)
          qspin(4, ia) = qspin(4, ia) + 0.5_dp*(Dscf(ind,4)+Dscf(ind,8)) * S(ind)
        else
          qspin(:, ia) = qspin(:, ia) + Dscf(ind,:) * S(ind)
        endif
      end do
    end do
    !
#ifdef MPI
    ! Global reduction of spin components
    do ia = 1, na_u
      call globalize_sum(qspin(1:spin%Grid, ia), qtmp(1:spin%Grid))
      qspin(1:spin%Grid, ia) = qtmp(1:spin%Grid)
    enddo
#endif
    !
    ! Get magnetic moment from charge density
    do ia = 1, na_u
      call spnvec(spin%Grid, qspin(:,ia), qtot, Stot(ia), Svec(:, ia))
    enddo

#ifdef DEBUG
    if (ionode) write(stdout, "(a4,a30,a10)") "ia", "Svec          ", "Stot"
    do ia = 1, na_u
      if (ionode) write(stdout, "(i4,3f10.6,f10.6)") ia, Svec(:, ia), Stot(ia)
    enddo
#endif

  end subroutine atomic_magnetic_moment


  function ainv(a)
    !
    ! Compute the inverse of a 3x3 matrix.

    ! variables
    use precision, only: dp
    !
    implicit none
    !
    real(kind=dp), dimension(3,3), intent(in) :: a
    real(kind=dp), dimension(3,3)             :: ainv
    !
    real(kind=dp), parameter :: eps = 1.0d-13
    real(kind=dp) :: det
    real(kind=dp), dimension(3,3) :: cofactor
    !
    !
    det =   a(1,1)*a(2,2)*a(3,3)  &
          - a(1,1)*a(2,3)*a(3,2)  &
          - a(1,2)*a(2,1)*a(3,3)  &
          + a(1,2)*a(2,3)*a(3,1)  &
          + a(1,3)*a(2,1)*a(3,2)  &
          - a(1,3)*a(2,2)*a(3,1)
    !
    if (abs(det) .le. eps) then
      ainv = 0.0d0
      write(*,*)"ERROR IN CALCULATING MATRIX INVERSE"
      stop
    end if
    !
    cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
    cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
    cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
    cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
    cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
    cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
    cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))
    !
    ainv = transpose(cofactor) / det

  end function ainv


  subroutine simpson(nr, dr, f, integral)
    !
    ! simpson's rule integration.
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2001 PWSCF group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    implicit none

    integer, intent(in) :: nr ! number of grid points (should be odd)
    real(DP), intent(in) :: dr (nr) ! dr(i)/di
    real(DP), intent(in) :: f (nr) ! function to be integrated
    real(DP), intent(out):: integral ! sum_i c_i f(i)*dr(i) = int_0^infty f(r) dr
    !
    real(DP) :: f1, f2, f3, r12
    integer :: i

    if ( mod(nr, 2) /= 1 ) then
      stop "ERROR: simpson: nr must be odd."
    endif

    integral = 0.0d0
    r12 = 1.0d0 / 12.0d0
    f3 = f(1) * dr(1) * r12

    do i = 2, nr - 1, 2
      f1 = f3
      f2 = f(i) * dr(i) * r12
      f3 = f(i + 1) * dr(i + 1) * r12
      integral = integral + 4.0d0 * f1 + 16.0d0 * f2 + 4.0d0 * f3
    enddo

  end subroutine simpson


  SUBROUTINE dosplineint( old_mesh, old_vec, new_mesh, new_vec)
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2004-2006 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    IMPLICIT NONE
    !
    REAL (DP), INTENT(IN)  :: old_mesh(:), new_mesh(:)
    REAL (DP), INTENT(IN)  :: old_vec(:)
    REAL (DP), INTENT(OUT) :: new_vec(:)
    !
    REAL (DP), ALLOCATABLE :: d2y(:)
    INTEGER                :: i
    INTEGER                :: old_dim, new_dim
    !
    !
    old_dim = SIZE( old_vec )
    new_dim = SIZE( new_vec )
    !
    ! IF ( old_dim /= SIZE( old_mesh ) ) CALL errore( 'dosplineint', 'dimensions of old_mesh and old_vec do not match', 1 )
    IF ( old_dim /= SIZE( old_mesh ) ) stop 'ERROR: dosplineint: dimensions of old_mesh and old_vec do not match'
    !
    ! IF ( new_dim /= SIZE( new_mesh ) ) CALL errore( 'dosplineint', 'dimensions of new_mesh and new_vec do not match', 1 )
    IF ( new_dim /= SIZE( new_mesh ) ) stop 'ERROR: dosplineint: dimensions of new_mesh and new_vec do not match'
    !
    ALLOCATE( d2y( old_dim ) )
    !
    d2y = 0
    !
    CALL spline( old_mesh , old_vec(:), 0.0_DP, 0.0_DP, d2y  )
    !
    DO i = 1, new_dim
        !
        new_vec(i) = splint( old_mesh, old_vec(:), d2y, new_mesh(i) )
        !
    END DO
    !
    DEALLOCATE( d2y )
    !
  END SUBROUTINE dosplineint


  SUBROUTINE spline( xdata, ydata, startu, startd, d2y )
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2004-2006 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: xdata(:), ydata(:), startu, startd
    REAL(DP), INTENT(OUT) :: d2y(:)
    !
    INTEGER               :: i, k, ydim
    REAL(DP)              :: p, sig
    REAL(DP), ALLOCATABLE :: u(:)
    !
    !
    ydim = SIZE( ydata )
    !
    ALLOCATE( u( ydim ) )
    !
    u(1)   = startu
    d2y(1) = startd
    !
    DO  i = 2, ydim - 1
        !
        sig    = ( xdata(i) - xdata(i-1) ) / ( xdata(i+1) - xdata(i-1) )
        p      = sig * d2y(i- 1) + 2.0_DP
        d2y(i) = ( sig - 1.0_DP ) / p
        u(i)   = ( 6.0_DP * ( ( ydata(i+1) - ydata(i) ) / &
                  ( xdata(i+1) - xdata(i) ) - ( ydata(i) - ydata(i-1) ) / &
                  ( xdata(i) - xdata(i-1) ) ) / &
                  ( xdata(i+1) - xdata(i-1) ) - sig * u(i-1) ) / p
        !
    END DO
    !
    d2y(ydim) = 0
    !
    DO  k = ydim - 1, 1, -1
        !
        d2y(k) = d2y(k) * d2y(k+1) + u(k)
        !
    END DO
    !
    DEALLOCATE( u )
    !
  END SUBROUTINE spline


  FUNCTION splint( xdata, ydata, d2y, x )
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2004-2006 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xdata(:), ydata(:), d2y(:)
    REAL(DP), INTENT(IN) :: x
    !
    REAL(DP) :: splint
    INTEGER  :: khi, klo, xdim
    REAL(DP) :: a, b, h
    !
    !
    xdim = SIZE( xdata )
    !
    klo = 1
    khi = xdim
    !
    klo = MAX( MIN( locate( xdata, x ), ( xdim - 1 ) ), 1 )
    !
    khi = klo + 1
    !
    h = xdata(khi) - xdata(klo)
    !
    a = ( xdata(khi) - x ) / h
    b = ( x - xdata(klo) ) / h
    !
    splint = a * ydata(klo) + b * ydata(khi) + &
              ( ( a**3 - a ) * d2y(klo) + ( b**3 - b ) * d2y(khi) ) * &
              ( h**2 ) / 6.0_DP

  END FUNCTION splint


  FUNCTION locate( xx, x )
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2004-2006 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xx(:)
    REAL(DP), INTENT(IN) :: x
    !
    INTEGER :: locate
    INTEGER :: n, jl, jm, ju
    LOGICAL :: ascnd
    !
    !
    n     = SIZE( xx )
    ascnd = ( xx(n) >= xx(1) )
    jl    = 0
    ju    = n + 1
    !
    main_loop: DO
      !
      IF ( ( ju - jl ) <= 1 ) EXIT main_loop
      !
      jm = ( ju + jl ) / 2
      !
      IF ( ascnd .EQV. ( x >= xx(jm) ) ) THEN
          !
          jl = jm
          !
      ELSE
          !
          ju = jm
          !
      END IF
      !
    END DO main_loop
    !
    IF ( x  ==  xx(1) ) THEN
      !
      locate = 1
      !
    ELSE IF ( x  ==  xx(n) ) THEN
      !
      locate = n - 1
      !
    ELSE
      !
      locate = jl
      !
    END IF
    !
  END FUNCTION locate

end module siesta2intw_utils