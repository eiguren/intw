module siesta2intw_utils

  use precision, only: dp

  implicit none

  public :: cmplx_0, cmplx_1, cmplx_i

  public :: atomic_magnetic_moment, ainv, simpson

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

end module siesta2intw_utils