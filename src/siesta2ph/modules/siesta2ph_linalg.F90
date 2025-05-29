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
module siesta2ph_linalg

  implicit none

  public :: ainv, &
            diagonalize_cmat, &
            rank, &
            pseudoinverse

  private

contains

  function ainv(a)
    !
    ! Compute the inverse of a 3x3 matrix.
    !
    use precision, only: dp
    !
    implicit none
    !
    real(kind=dp), dimension(3,3), intent(in) :: a
    real(kind=dp), dimension(3,3) :: ainv
    !
    real(kind=dp), parameter :: eps = 1.0d-13
    real(kind=dp) :: det
    real(kind=dp), dimension(3,3) :: cofactor


    ! Compute determinant
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


  subroutine diagonalize_cmat(A,w)
    !
    ! Diagonalize a complex Hermitian matrix
    !
    use precision, only: dp
    !
    implicit none
    !
    external :: zheev
    !
    ! input variables
    complex(kind=dp), intent(inout) :: A(:,:) ! matrix to diagonalize on input, eigenvector on output
    real(kind=dp), intent(out) :: w(:) ! eigenvalues
    !
    complex(kind=dp), allocatable, dimension(:) :: WORK
    real(kind=dp), allocatable, dimension(:) :: RWORK
    integer :: N, LWORK, INFO


    N = size((A),DIM=1)
    if ( size((A),DIM=2) .ne. N ) stop "diagonalize_cmat: ERROR"
    if ( size((w),DIM=1) .ne. N ) stop "diagonalize_cmat: ERROR"
    !
    LWORK = -1
    allocate(WORK(1))
    allocate(RWORK(max(1, 3*N-2)))
    !
    call zheev( "V", "U", N, A, N, w, WORK, LWORK, RWORK, INFO )
    !
    LWORK = int(WORK(1))
    !
    deallocate(WORK)
    allocate(WORK(LWORK))
    !
    call zheev( "V", "U", N, A, N, w, WORK, LWORK, RWORK, INFO )

  end subroutine diagonalize_cmat


  subroutine rank(A,r)
    !
    ! Compute the rank of a matrix using column pivoting QR factorization.
    !
    use precision, only: dp
    !
    implicit none
    !
    external :: dgeqp3
    !
    real(kind=dp), dimension(:,:), intent(in) :: A
    integer, intent(out) :: r
    !
    integer :: M ! rows of A
    integer :: N ! columns of A
    integer :: LDA ! rows of A
    integer :: LWORK
    integer :: INFO
    integer, allocatable, dimension(:) :: JPVT ! Permutations
    real(kind=dp), allocatable, dimension(:) :: TAU
    real(kind=dp), allocatable, dimension(:) :: WORK
    !
    real(kind=dp), allocatable, dimension(:,:) :: A_tmp
    !
    integer :: i
    !
#ifdef DEBUG
    external :: DORG2R
    integer :: LDT
    integer, allocatable, dimension(:,:) :: Pmat
    real(kind=dp), allocatable, dimension(:,:) :: Qmat, Rmat, II
    integer :: r_check
#endif


    M = size((A),DIM=1)
    N = size((A),DIM=2)
    allocate(A_tmp(M,N))
    A_tmp = A
    LDA = M
    allocate(JPVT(N))
    JPVT = 0
    allocate(TAU(min(M,N)))
    LWORK = -1
    allocate(WORK(1))
    !
    call dgeqp3( M, N, A_tmp, LDA, JPVT, TAU, WORK, LWORK, INFO )
    !
    LWORK = int(WORK(1))
    !
    deallocate(WORK)
    !
    allocate(WORK(LWORK))
    call dgeqp3( M, N, A_tmp, LDA, JPVT, TAU, WORK, LWORK, INFO )
    !
    r = 0
    do i=1,N
      if (abs(A_tmp(i,i)) > 0.00001_dp ) r = r + 1
    enddo
    !
#ifdef DEBUG

    !getting R matrix
    allocate(Rmat(M,N))
    Rmat = 0.0_dp
    do i=1,N
      Rmat(i,i:N) = A_tmp(i,i:N)
    enddo

    !getting Q matrix
    allocate(Qmat(M,M), II(M,M))
    Qmat = A_tmp
    call DORG2R(M, M, LDT, A_tmp, M, TAU, WORK, INFO)
    Qmat = A_tmp(1:3,1:3)

    !getting the P matrix
    allocate(Pmat(N,N))
    Pmat = 0
    do i=1,N
      Pmat(JPVT(i),i) = 1
    enddo

    print*, "A="
    print*, A(1,:)
    print*, A(2,:)
    print*, A(3,:)

     ! print*, "A_tmp="
     ! print*, A_tmp(1,:)
     ! print*, A_tmp(2,:)
     ! print*, A_tmp(3,:)

    print*, "Q="
    print*, Qmat(1,:)
    print*, Qmat(2,:)
    print*, Qmat(3,:)

    print*, "Q.Q**T="
    A_tmp = matmul(Qmat, transpose(Qmat))
    print*, A_tmp(1,:)
    print*, A_tmp(2,:)
    print*, A_tmp(3,:)

    print*, "R="
    print*, Rmat(1,:)
    print*, Rmat(2,:)
    print*, Rmat(3,:)

    print*, "P="
    print*, Pmat(1,:)
    print*, Pmat(2,:)
    print*, Pmat(3,:)

    ! print*, "Q.R="
    !A_tmp = matmul(Qmat, (Rmat))
    ! print*, A_tmp(1,:)
    ! print*, A_tmp(2,:)
    ! print*, A_tmp(3,:)

    print*, "Q.R.P="
    A_tmp = matmul(matmul(Qmat, Rmat), Pmat)
    print*, A_tmp(1,:)
    print*, A_tmp(2,:)
    print*, A_tmp(3,:)

    ! print*, "JPVT="
    ! print*, JPVT

    ! print*, "TAU="
    ! print*, TAU

    call rank_old(A, r_check)
    print*, "rank:", r
    print*, "rank_old:", r_check
    if (r.ne.r_check) stop "rank: ERROR"

#endif

  end subroutine rank


  subroutine rank_old(A,r)
    !
    ! Compute the rank of a matrix using SVD.
    !
    use precision, only: dp
    !
    implicit none
    !
    external :: dgesvd
    !
    real(kind=dp), dimension(:,:), intent(in) :: A
    integer, intent(out) :: r
    !
    integer :: M ! rows of A
    integer :: N ! columns of A
    integer :: LDA
    integer :: LDU
    integer :: LDVT
    integer :: LWORK
    integer :: INFO
    real(kind=dp), allocatable, dimension(:) :: S
    real(kind=dp), allocatable, dimension(:,:) :: U
    real(kind=dp), allocatable, dimension(:,:) :: VT
    real(kind=dp), allocatable, dimension(:) :: WORK
    !
    real(kind=dp), allocatable, dimension(:,:) :: A_tmp


    M = size((A),DIM=1)
    N = size((A),DIM=2)
    allocate(A_tmp(M,N))
    A_tmp=A
    LDA = M
    LDU = M
    LDVT = min(M,N)
    allocate(S(LDVT))
    allocate(U(LDU,LDVT))
    allocate(VT(LDVT,N))
    LWORK = -1
    allocate(WORK(1))
    !
    call dgesvd( "N", "N", M, N, A_tmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    !
    LWORK = int(WORK(1))
    !
    deallocate(WORK)
    !
    allocate(WORK(LWORK))
    !
    call dgesvd( "N", "N", M, N, A_tmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    !
    r = count(abs(S) .gt. 0.0000000001_dp)

  end subroutine rank_old


  subroutine pseudoinverse(A)
    !
    ! Compute the pseudo-inverse of a matrix.
    ! See https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse
    !
    use precision, only: dp
    !
    implicit none
    !
    external :: dgesvd, dscal, dgemm
    !
    real(kind=dp), dimension(:,:), intent(inout) :: A
    !
    integer :: M ! rows of A
    integer :: N ! columns of A
    integer :: LDA
    integer :: LDU
    integer :: LDVT
    integer :: LWORK
    integer :: INFO
    real(kind=dp), allocatable, dimension(:) :: S
    real(kind=dp), allocatable, dimension(:,:) :: U
    real(kind=dp), allocatable, dimension(:,:) :: VT
    real(kind=dp), allocatable, dimension(:) :: WORK
    !
    integer :: i
    real(kind=dp) :: ss
    real(kind=dp), allocatable, dimension(:) :: tmp


    M = size((A),DIM=1)
    N = size((A),DIM=2)
    LDA = M
    LDU = M
    LDVT = min(M,N)
    allocate(S(LDVT))
    allocate(U(LDU,LDVT))
    allocate(VT(LDVT,N))
    LWORK = -1
    allocate(WORK(1))
    !
    call dgesvd( "S", "S", M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    !
    LWORK = int(WORK(1))
    !
    deallocate(WORK)
    !
    allocate(WORK(LWORK))
    !
    call dgesvd( "S", "S", M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    !
    do i=1,LDVT
    if(S(i) > 1.0e-9) then
      ss=1.0_dp/S(i)
    else
      ss=S(i)
    endif
      call dscal(M, ss, U(:,i), 1)
    enddo
    !
    call dgemm("T", "T", N, M, LDVT, 1.0_dp, VT, LDVT, U, LDU, 0.0_dp, A, N)
    !
    allocate(tmp(M*N))
    tmp = reshape(A,shape(tmp))
    do i=1,M
      A(i,:) = tmp(N*(i-1)+1:i*N)
    enddo

  end subroutine pseudoinverse

end module siesta2ph_linalg
