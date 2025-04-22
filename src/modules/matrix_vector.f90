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
module intw_matrix_vector

  !------------------------------------------------------------------!
  ! Module that contains all the necessary functions and subroutines !
  ! to perform matrix and vector operations.                         !
  !------------------------------------------------------------------!

  use kinds, only: dp

  public :: norma, cross, area_vec
  public :: ainv, det, minor_matrix
  public :: cmplx_ainv, cmplx_det, cmplx_minor_matrix, cmplx_trace

  private

contains

  function norma(v)

    implicit none

    real(kind=dp), intent(in) :: v(:)
    real(kind=dp) :: norma

    norma = norm2(v)
    ! norma = dot_product(v, v)

  end function norma


  function cross(a, b)

    implicit none

    real(kind=dp), dimension(3) :: cross
    real(kind=dp), dimension(3), intent(in) :: a, b

    CROSS(1) = A(2) * B(3) - A(3) * B(2)
    CROSS(2) = A(3) * B(1) - A(1) * B(3)
    CROSS(3) = A(1) * B(2) - A(2) * B(1)

  end function cross


  function area_vec(v1,v2)

    implicit none

    real(kind=dp), dimension(3), intent(in) :: v1,v2
    real(kind=dp) :: area_vec


    area_vec = (v1(2)*v2(3) - v1(3)*v2(2))**2 + &
               (v1(3)*v2(1) - v1(1)*v2(3))**2 + &
               (v1(1)*v2(2) - v1(2)*v2(1))**2
    area_vec = sqrt(area_vec)*0.5_dp

  end function area_vec


  function ainv(A)

    use intw_useful_constants, only: eps_14
    implicit none

    real(kind=dp), intent(in) :: A(:,:)
    real(kind=dp) :: ainv(size(A, dim=1),size(A, dim=2))

    real(kind=dp) :: determinant
    real(kind=dp) :: cofactor(size(A, dim=1),size(A, dim=2))
    integer :: N, M, i, j


    N = size(A, dim=1)
    M = size(A, dim=2)

    if (N/=M) stop "ERROR: ainv: A must be a square matrix"

    determinant = det(A)

    if ( abs(determinant) < eps_14 ) stop "ERROR IN CALCULATING MATRIX INVERSE"

    do i = 1, N
      do j = 1, M
        cofactor(i,j) = (-1)**(i+j)*det(minor_matrix(A, i, j))
      enddo
    end do

    ainv = transpose(cofactor) / determinant

  end function ainv


  recursive function det(A) result(d)

    implicit none

    real(kind=dp), intent(in) :: A(:,:)
    real(kind=dp) :: d

    integer :: N, M, i


    N = size(A, dim=1)
    M = size(A, dim=2)

    if (N/=M) stop "ERROR: det: A must be a square matrix"

    if (N==1) then
      d = A(1,1)
    else if (N>1) then
      d = 0.0_dp
      do i=1,N
        d = d + (-1)**(i+1) * A(1,i) * det(minor_matrix(A, 1, i))
      enddo
    else
      stop "ERROR: det: Wrong size"
    end if

  end function det


  function minor_matrix(A, ROW, COL)

    implicit none

    real(kind=dp), intent(in) :: A(:,:)
    integer, intent(in) :: ROW, COL

    real(kind=dp) :: minor_matrix(size(A, dim=1)-1, size(A, dim=2)-1)

    integer :: N, M, i, ii, j, jj


    N = size(A, dim=1)
    M = size(A, dim=2)

    ii = 0
    do i = 1, N
      if (i==ROW) cycle
      ii = ii + 1
      jj = 0
      do j = 1, N
        if (j==COL) cycle
        jj = jj + 1
        minor_matrix(ii, jj) = A(i, j)
      end do
    end do

  end function minor_matrix


  function cmplx_ainv(A)

    use intw_useful_constants, only: cmplx_0, eps_14

    implicit none

    complex(kind=dp), intent(in) :: A(:,:)
    complex(kind=dp) :: cmplx_ainv(size(A, dim=1),size(A, dim=2))

    complex(kind=dp) :: determinant
    complex(kind=dp) :: cofactor(size(A, dim=1),size(A, dim=2))
    integer :: N, M, i, j


    N = size(A, dim=1)
    M = size(A, dim=2)

    if (N/=M) stop "ERROR: cmplx_ainv: A must be a square matrix"

    determinant = cmplx_det(A)

    if ( abs(determinant) < eps_14 ) stop "ERROR IN CALCULATING MATRIX INVERSE"

    do i = 1, N
      do j = 1, M
        cofactor(i,j) = (-1)**(i+j)*cmplx_det(cmplx_minor_matrix(A, i, j))
      enddo
    end do

    cmplx_ainv = transpose(cofactor) / determinant

  end function cmplx_ainv


  recursive function cmplx_det(A) result(d)

    use intw_useful_constants, only: cmplx_0

    implicit none

    complex(kind=dp), intent(in) :: A(:,:)
    complex(kind=dp) :: d

    integer :: N, M, i


    N = size(A, dim=1)
    M = size(A, dim=2)

    if (N/=M) stop "ERROR: cmplx_det: A must be a square matrix"

    if (N==1) then
      d = A(1,1)
    else if (N>1) then
      d = cmplx_0
      do i=1,N
        d = d + (-1)**(1+i) * A(1,i) * cmplx_det(cmplx_minor_matrix(A, 1, i))
      enddo
    else
      stop "ERROR: cmplx_det: Wrong size"
    end if

  end function cmplx_det


  function cmplx_minor_matrix(A, ROW, COL)

    implicit none

    complex(kind=dp), intent(in) :: A(:,:)
    integer, intent(in) :: ROW, COL

    complex(kind=dp) :: cmplx_minor_matrix(size(A, dim=1)-1, size(A, dim=2)-1)

    integer :: N, M, i, ii, j, jj


    N = size(A, dim=1)
    M = size(A, dim=2)

    ii = 0
    do i = 1, N
      if (i==ROW) cycle
      ii = ii + 1
      jj = 0
      do j = 1, N
        if (j==COL) cycle
        jj = jj + 1
        cmplx_minor_matrix(ii, jj) = A(i, j)
      end do
    end do

  end function cmplx_minor_matrix


  function cmplx_trace(A)

    use intw_useful_constants, only: cmplx_0

    implicit none

    complex(kind=dp), intent(in) :: A(:,:)
    complex(kind=dp) :: cmplx_trace
    integer :: N, M, i

    N = size(A, dim=1)
    M = size(A, dim=2)

    if (N/=M) stop "ERROR: cmplx_trace: A must be a square matrix"

    cmplx_trace = cmplx_0
    do i = 1, N
      cmplx_trace = cmplx_trace + A(i,i)
    enddo

  end function cmplx_trace

end module intw_matrix_vector
