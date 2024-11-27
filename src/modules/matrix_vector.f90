module intw_matrix_vector

!------------------------------------------------------------------!
! Module that contains all the necessary functions and subroutines !
! to perform matrix and vector operations.                         !
!------------------------------------------------------------------!

  use kinds, only: dp

  public :: norma, cross, area_vec
  public :: ainv, det, cof, minor_matrix
  public :: cmplx_ainv, cmplx_det, cmplx_cof, cmplx_minor_matrix, cmplx_trace

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

    real(kind=dp), dimension(3)  :: cross
    real(kind=dp), dimension(3), intent(in) :: a, b

    CROSS(1) = A(2) * B(3) - A(3) * B(2)
    CROSS(2) = A(3) * B(1) - A(1) * B(3)
    CROSS(3) = A(1) * B(2) - A(2) * B(1)

  end function cross


  function  area_vec(v1,v2)

    implicit none

    real(kind=dp), dimension(3), intent(in) :: v1,v2
    real(kind=dp) :: area_vec


    area_vec =  (v1(2)*v2(3) - v1(3)*v2(2))**2 + &
                (v1(3)*v2(1) - v1(1)*v2(3))**2 + &
                (v1(1)*v2(2) - v1(2)*v2(1))**2
    area_vec = sqrt(area_vec)*0.5_dp

  end function area_vec


  function ainv(A)

    implicit none

    real(kind=dp), intent(in) :: A(:,:)
    real(kind=dp) :: ainv(size(A, dim=1),size(A, dim=2))

    real(kind=dp), parameter :: eps = 1.0d-13
    real(kind=dp) :: determinant
    real(kind=dp) :: cofactor(size(A, dim=1),size(A, dim=2))
    integer :: N, M, i, j


    N = size(A, dim=1)
    M = size(A, dim=2)

    if ( N /= M ) stop "ERROR: ainv: A must be a square matrix"

    determinant = det(A)

    if (abs(determinant) < eps) then
      ainv = 0.0d0
      write(*,*)"ERROR IN CALCULATING MATRIX INVERSE"
      stop
    end if

    do i = 1, N
      do j = 1, M
       cofactor(i,j) = cof(A, i, j)
      enddo
    end do

    ainv = transpose(cofactor) / determinant

  end function ainv


  recursive function det(A)

    implicit none

    real(kind=dp), intent(in) :: A(:,:)
    real(kind=dp) :: det

    integer :: N, M, i


    N = size(A, dim=1)
    M = size(A, dim=2)

    if  ( N /= M ) stop "ERROR: det: A must be a square matrix"

    if (N == 1) then
      det = A(1,1)
    else if (N > 1) then
      det = 0.0_dp
      do i=1,N
        det = det + A(1,i) * cof(A, 1, i)
      enddo
    else
      stop "ERROR: det: Wrong size"
    end if

  end function det


  function cof(A, ROW, COL)

    implicit none

    real(kind=dp), intent(in) :: A(:,:)
    integer, intent(in) :: ROW, COL

    real(kind=dp) :: cof


    cof = (-1)**(ROW+COL)*det(minor_matrix(A, ROW, COL))

  end function cof


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

    use intw_useful_constants, only: cmplx_0

    implicit none

    complex(kind=dp), intent(in) :: A(:,:)
    complex(kind=dp) :: cmplx_ainv(size(A, dim=1),size(A, dim=2))

    real(kind=dp), parameter :: eps = 1.0d-13
    complex(kind=dp) :: determinant
    complex(kind=dp) :: cofactor(size(A, dim=1),size(A, dim=2))
    integer :: N, M, i, j


    N = size(A, dim=1)
    M = size(A, dim=2)

    if  ( N /= M ) stop "ERROR: cmplx_ainv: A must be a square matrix"

    determinant =  cmplx_det(A)

    if (abs(determinant) < eps) then
      cmplx_ainv = cmplx_0
      write(*,*)"ERROR IN CALCULATING MATRIX INVERSE"
      stop
    end if

    do i = 1, N
      do j = 1, M
       cofactor(i,j) = cmplx_cof(A, i, j)
      enddo
    end do

    cmplx_ainv = transpose(cofactor) / determinant

  end function cmplx_ainv


  recursive function cmplx_det(A)

    use intw_useful_constants, only: cmplx_0

    implicit none

    complex(kind=dp), intent(in) :: A(:,:)
    complex(kind=dp) :: cmplx_det

    integer :: N, M, i


    N = size(A, dim=1)
    M = size(A, dim=2)

    if  ( N /= M ) stop "ERROR: cmplx_det: A must be a square metrix"

    if (N == 1) then
      cmplx_det = A(1,1)
    else if (N > 1) then
      cmplx_det = cmplx_0
      do i=1,N
        cmplx_det = cmplx_det + A(1,i) * cmplx_cof(A, 1, i)
      enddo
    else
      stop "ERROR: cmplx_det: Wrong size"
    end if

  end function cmplx_det


  function cmplx_cof(A, ROW, COL)

    implicit none

    complex(kind=dp), intent(in) :: A(:,:)
    integer, intent(in) :: ROW, COL

    real(kind=dp) :: cmplx_cof


    cmplx_cof = (-1)**(ROW+COL)*cmplx_det(cmplx_minor_matrix(A, ROW, COL))

  end function cmplx_cof


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


  function cmplx_trace (mat)

    use intw_useful_constants, only: cmplx_0

    implicit none

    complex(kind=dp), intent(in) :: mat(:,:)
    complex(kind=dp) :: cmplx_trace
    integer :: d1, d2, i

    d1 = size(mat(:,1))
    d2 = size(mat(1,:))

    if (d1.ne.d2) then
      write(*,*)"Errorea cmplx_trace"
      stop
    end if

    cmplx_trace = cmplx_0
    do i=1,d1
      cmplx_trace = cmplx_trace + mat(i,i)
    enddo

  end function  cmplx_trace

end module intw_matrix_vector
