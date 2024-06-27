!----------------------------------------------------------------------------!
! intw project.
!
!----------------------------------------------------------------------------!
module intw_utility
!----------------------------------------------------------------------------!
!
!       This module contains useful functions which implement
!       common tasks. It will be VERY useful to call these functions
!       instead of reimplementing them every time they are needed,
!       especially to insure CONSISTENCY.
!
!----------------------------------------------------------------------------!
use kinds, only: dp

  implicit none
  !
  ! subroutines
  public :: get_timing, &
            joint_to_triple_index_g, triple_to_joint_index_g, &
            joint_to_triple_index_r, triple_to_joint_index_r, &
            generate_kmesh, find_neighbor, find_maximum_index_int, test_qpt_on_fine_mesh, &
            find_k_1BZ_and_G, cryst_to_cart, &
            HPSORT, HPSORT_real, hpsort_eps, &
            find_r_in_WS_cell, errore, simpson, gaussian, sphb, &
            real_ylmr2
  !
  ! functions
  public :: intgr_spline_gaussq, ainv, det, cmplx_ainv, cmplx_det, cmplx_trace, multiple, weight_ph, &
            qe_erf, qe_erfc, find_free_unit, conmesurate_and_coarser, &
            smeared_delta, fermi_dirac, area_vec
  !
  private
  !

contains

function intgr_spline_gaussq(xdata,ydata) result(batura)
! This is original from MBGROUP
! This function integrates exactly (numerical)
! the cubic spline interpolation of a given data in 1D
! this is so because 2nd order gauss cuadrature
! is exact of pol. of order 3.

use mcf_spline, only: spline_mcf, splint_mcf

real(kind=dp),intent(in) :: xdata(:),ydata(:)
!out
real(kind=dp):: batura
!local
real(kind=dp), dimension(size(xdata)) :: ys2
real(kind=dp) :: x1, x2, y1, y2
real(kind=dp),parameter :: onesqrt3=0.57735026919_dp
integer :: i


call SPLINE_MCF (xdata, ydata,  size(xdata), ys2)

batura=0.0_dp
do i=1,size(xdata)-1

 x1=(xdata(i+1)+xdata(i))/2.0_dp - (xdata(i+1)-xdata(i))/2.0_dp/onesqrt3
 x2=(xdata(i+1)+xdata(i))/2.0_dp + (xdata(i+1)-xdata(i))/2.0_dp/onesqrt3

 call splint_MCF (xdata, ydata, ys2, size(xdata), x1, y1)
 call splint_MCF (xdata, ydata, ys2, size(xdata), x2, y2)

 batura=batura + (y2 + y1) * (xdata(i+1)-xdata(i))/2.0_dp
enddo
end function intgr_spline_gaussq


  subroutine get_timing(time)
    !----------------------------------------------------------------------------!
    ! Timing the code is very important for optimized performance.
    ! However, some experimence indicates that "timing" is not so easy:
    ! time can overflow the buffer in which it is held, or it can give
    ! crazy numbers for (so far) unknown reasons. Thus, in order to be
    ! able to quickly modify how time is measured, a single subroutine
    ! will be called throuhout the code to measure time. This subroutine
    ! can then easily be modified as the code writer becomes
    ! aware of "better" timing algorithms.
    !----------------------------------------------------------------------------!
    implicit none

    real(dp), external :: dsecnd

    real(dp) :: time

    time = dsecnd()

  end subroutine get_timing


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

    if  ( N /= M ) stop "ERROR: ainv: A must be a square matrix"

    determinant =  det(A)

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
      cmplx_ainv = 0.0d0
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
      cmplx_det = 0.0_dp
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


  function conmesurate_and_coarser (nk1, nk2, nk3, nq1, nq2, nq3)
    !Asier&&Idoia 23 06 2014
    logical :: conmesurate_and_coarser
    integer, intent(in) :: nk1, nk2, nk3, nq1, nq2, nq3

    conmesurate_and_coarser = .true.

    !check if electron k grid does not contain the phonon q grid.
    if ((nk1<nq1).or.(nk2<nq2).or.(nk3<nq3)) conmesurate_and_coarser = .false.

    !check if q is contained in k

    if ( .not.multiple(nk1,nq1) ) conmesurate_and_coarser = .false.
    if ( .not.multiple(nk2,nq2) ) conmesurate_and_coarser = .false.
    if ( .not.multiple(nk3,nq3) ) conmesurate_and_coarser = .false.

  end function conmesurate_and_coarser


  function multiple (i1, i2)
    !Asier&&Idoia 23 06 2014
    !Check if two integers are multiple of each other.
    logical :: multiple
    integer, intent(in) :: i1, i2

    multiple = .false.

    if ( (i1/i2*i2==i1) .or. (i2/i1*i1==i2) ) then
      multiple = .true.
    endif

  end function  multiple


  subroutine joint_to_triple_index_g(n1, n2, n3, i_joint, i, j, k)
    !----------------------------------------------------------------------------!
    ! Compute the triple indices i, j, k of a point in a n1, n2, n3 mesh on
    ! reciprocal space from the joint index i_joint.
    !
    ! i_joint is a single integer from 1 to n1*n2*n3, which uniquely
    ! labels all the points of the mesh.
    !
    ! i, j, k are the triple indices, which uniquely label all the points
    ! of the mesh according to their coordinates in the mesh.
    !
    ! The 3rd triple index loops fastest:
    !   i_joint = k + (j-1)*n3 + (i-1)*n2*n3
    !----------------------------------------------------------------------------!

    implicit none

    integer, intent(in) :: n1, n2, n3
    integer, intent(in) :: i_joint
    integer, intent(out) :: i, j, k


    i = (i_joint - 1)/(n2*n3) + 1
    j = (i_joint - 1)/n3 - (i-1)*n2 + 1
    k = i_joint - int((i_joint-1)/n3) * n3 ! k = mod(i_joint-1, n3) + 1


  end subroutine joint_to_triple_index_g


  subroutine triple_to_joint_index_g(n1, n2, n3, i_joint, i, j, k)
    !----------------------------------------------------------------------------!
    ! Compute the joint index i_joint of a point in a n1, n2, n3 mesh on
    ! reciprocal space from the triple indices i, j, k.
    !
    ! i_joint is a single integer from 1 to n1*n2*n3, which uniquely
    ! labels all the points of the mesh.
    !
    ! i, j, k are the triple indices, which uniquely label all the points
    ! of the mesh according to their coordinates in the mesh.
    !
    ! The 3rd triple index loops fastest:
    !   i_joint = k + (j-1)*n3 + (i-1)*n2*n3
    !----------------------------------------------------------------------------!

    implicit none

    integer, intent(in) :: n1, n2, n3
    integer, intent(out) :: i_joint
    integer, intent(in) :: i, j, k


    i_joint = k + (j-1 + (i-1)*n2)*n3

  end subroutine triple_to_joint_index_g


  subroutine joint_to_triple_index_r(n1, n2, n3, i_joint, i, j, k)
    !----------------------------------------------------------------------------!
    ! Compute the triple indices i, j, k of a point in a n1, n2, n3 mesh on
    ! real space from the joint index i_joint.
    !
    ! i_joint is a single integer from 1 to n1*n2*n3, which uniquely
    ! labels all the points of the mesh.
    !
    ! i, j, k are the triple indices, which uniquely label all the points
    ! of the mesh according to their coordinates in the mesh.
    !
    ! The 1st triple index loops fastest:
    !   i_joint = i + (j-1)*n2 + (k-1)*n1*n2
    !----------------------------------------------------------------------------!

    implicit none

    integer, intent(in) :: n1, n2, n3
    integer, intent(in) :: i_joint
    integer, intent(out) :: i, j, k


    k = (i_joint - 1)/(n1*n2) + 1
    j = (i_joint - 1)/n1 - (k-1)*n2 + 1
    i = i_joint - int((i_joint-1)/n1) * n1 ! k = mod(i_joint-1, n1) + 1

  end subroutine joint_to_triple_index_r


  subroutine triple_to_joint_index_r(n1, n2, n3, i_joint, i, j, k)
    !----------------------------------------------------------------------------!
    ! Compute the triple indices i, j, k of a point in a n1, n2, n3 mesh on
    ! real space from the joint index i_joint.
    !
    ! i_joint is a single integer from 1 to n1*n2*n3, which uniquely
    ! labels all the points of the mesh.
    !
    ! i, j, k are the triple indices, which uniquely label all the points
    ! of the mesh according to their coordinates in the mesh.
    !
    ! The 1st triple index loops fastest:
    !   i_joint = i + (j-1)*n2 + (k-1)*n1*n2
    !----------------------------------------------------------------------------!

    implicit none

    integer, intent(in) :: n1, n2, n3
    integer, intent(out) :: i_joint
    integer, intent(in) :: i, j, k


    i_joint = i + (j-1 + (k-1)*n2)*n1

  end subroutine triple_to_joint_index_r


  subroutine generate_kmesh(kmesh,nk_1,nk_2,nk_3)
    !----------------------------------------------------------------------------!
    !     This subroutine builds the array of k vectors corresponding
    !     to MP indices nk_1,nk_2,nk_3, ordered according to their joint
    !     index.
    !----------------------------------------------------------------------------!
    implicit none

    integer, intent(in) :: nk_1, nk_2, nk_3
    real(dp), intent(out) :: kmesh(3,nk_1*nk_2*nk_3)

    integer :: i, j, k
    integer :: ikpt, nkmesh


    ! Build kmesh
    nkmesh = nk_1*nk_2*nk_3

    do ikpt = 1, nkmesh
      call joint_to_triple_index_g(nk_1,nk_2,nk_3,ikpt,i,j,k)

      kmesh(1, ikpt) = dble(i-1)/nk_1
      kmesh(2, ikpt) = dble(j-1)/nk_2
      kmesh(3, ikpt) = dble(k-1)/nk_3
    end do


  end subroutine generate_kmesh

  subroutine find_neighbor(kpoint,nk_1,nk_2,nk_3,i_k,j_k,k_k)
    !----------------------------------------------------------------------------!
    !     Given a kpoint = (kx, ky, kz) in crystal coordinates, which lies
    !     inside the 1BZ, namely 0 <= kx,ky,kz < 1, this subroutine returns
    !     the coordinates of the neighbor which lies at the origin of the tricubic
    !     4x4x4 interpolation grid, in the format (i_k,j_k,k_k).
    !     The crystal coordinates of the neigbhor are given by
    !                     ((i_k-1)/nk_1, (j_k-1)/nk_2, (k_k-1)/nk_3).
    !
    !----------------------------------------------------------------------------!
    implicit none

    real(dp), intent(in) :: kpoint(3)
    integer, intent(in) :: nk_1, nk_2, nk_3

    integer, intent(out) :: i_k, j_k, k_k

    i_k = 1 + nint(kpoint(1)*nk_1)
    j_k = 1 + nint(kpoint(2)*nk_2)
    k_k = 1 + nint(kpoint(3)*nk_3)

  end subroutine find_neighbor

  subroutine find_maximum_index_int(array,sze,i_max)
    !----------------------------------------------------------------------------!
    !     Given an array of integers array(sze), this subroutine returns
    !     the largest value in the array.
    !----------------------------------------------------------------------------!
    implicit none

    integer, intent(in) :: sze, array(sze)
    integer, intent(out) :: i_max

    integer :: i, maximum

    i_max   = 1
    maximum = array(1)

    do i=2,sze
      if ( array(i) > maximum) then
        maximum = array(i)
        i_max   = i
      end if
    end do

  end subroutine find_maximum_index_int

  subroutine test_qpt_on_fine_mesh(qpt,nk1s,nk2s,nk3s,test_qpt,i_qpt1,i_qpt2,i_qpt3)
    !----------------------------------------------------------------------------!
    !
    ! This subroutine tests whether the q-point, qpt(3), is of the form
    !
    !     qpt(:) = [ i_qpt1-1   i_qpt2-1   i_qpt3-1 ] + (I,J,K)
    !              [ -------- , -------- , -------- ]
    !              [   nk1s       nk2s       nk3s   ]
    !
    !
    ! where I,J,K are integers, and i_qpt(1,2,3) are integers between 1 and nk(123)s.
    !
    ! It is important to perform this task in a defensive way.
    ! The fortran internal functions nint, floor, and modulo are DANGEROUS.
    ! IT IS CRUCIAL TO DO THIS RIGHT ONCE AND FOR ALL.
    !
    !----------------------------------------------------------------------------!
    use intw_useful_constants, only: eps_8

    implicit none

    ! input
    integer, intent(in) :: nk1s, nk2s, nk3s ! the fine mesh parameters
    real(dp), intent(in) :: qpt(3) ! the q-point

    ! output
    logical, intent(out):: test_qpt ! is the q-point on the fine mesh?
    integer, intent(out):: i_qpt1, i_qpt2, i_qpt3

    ! internal variables
    real(dp) :: nqpt_1, nqpt_2, nqpt_3
    real(dp) :: err_1, err_2, err_3 ! the q-point


    ! find err
    nqpt_1 = dble(nk1s)*qpt(1)
    nqpt_2 = dble(nk2s)*qpt(2)
    nqpt_3 = dble(nk3s)*qpt(3)

    err_1 = (nqpt_1-nint(nqpt_1))/dble(nk1s)
    err_2 = (nqpt_2-nint(nqpt_2))/dble(nk2s)
    err_3 = (nqpt_3-nint(nqpt_3))/dble(nk3s)


    test_qpt =  err_1 <  eps_8 .and.  &
                err_2 <  eps_8 .and.  &
                err_3 <  eps_8

    if (test_qpt) then
      ! q is on mesh; find its coordinates.

      i_qpt1 = modulo(nint(nqpt_1),nk1s)+1
      i_qpt2 = modulo(nint(nqpt_2),nk2s)+1
      i_qpt3 = modulo(nint(nqpt_3),nk3s)+1

    else
      ! q is not on mesh! throw an error outside the subroutine
      i_qpt1 = 0
      i_qpt2 = 0
      i_qpt3 = 0
    end if

  end subroutine test_qpt_on_fine_mesh

  subroutine find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpt_in_1BZ,G)
    !----------------------------------------------------------------------------!
    ! Given a kpoint(3) in crystal coordinates, this subroutine
    ! generates:
    !            - k_1BZ(3),   with 0 <= k_1BZ(i) < 1
    !            - G(3),       G(i) an integer
    !            - i,j,k       the triple coordinates of the point in the mesh.
    !
    ! It it EXTREMELY important to have a subroutine which performs this
    ! task in a defensive way. The fortran internal functions nint, floor,
    ! and modulo are DANGEROUS. IT IS CRUCIAL TO DO THIS RIGHT ONCE AND FOR
    ! ALL.
    !
    ! The basic relationship is
    !    kpoint = k_1BZ + G
    !
    ! Assume
    !    kpoint(l) =  G(l) + (i_l-1)/nk_l+ epsilon ,
    ! with
    !    l    = 1, 2, 3
    !    i_l  = i, j ,k
    !    nk_l = nk1, nk2, nk3
    !    G(l) an integer
    ! epsilon numerical noise
    !----------------------------------------------------------------------------!
    implicit none

    ! input
    integer, intent(in) :: nk1, nk2, nk3             ! the mesh parameters
    real(dp), intent(in) :: kpoint(3)                ! the kpoint of interest

    ! output
    integer, intent(out) :: i, j, k                  ! the triple coordinates of k_1BZ
    real(dp), intent(out) :: kpt_in_1BZ(3)           ! the k point in the 1BZ
    integer, intent(out) :: G(3)                     ! the translation vector

    ! internal variables
    integer :: nG_im1, nG_jm1, nG_km1
    integer :: im1, jm1, km1


    ! this step kills epsilon
    nG_im1 = nint(kpoint(1)*dble(nk1))
    nG_jm1 = nint(kpoint(2)*dble(nk2))
    nG_km1 = nint(kpoint(3)*dble(nk3))

    ! this step gets rid of G
    im1 = modulo(nG_im1,nk1)
    jm1 = modulo(nG_jm1,nk2)
    km1 = modulo(nG_km1,nk3)

    ! G can be extracted. This division must be exact.
    G(1) = (nG_im1 - im1)/nk1
    G(2) = (nG_jm1 - jm1)/nk2
    G(3) = (nG_km1 - km1)/nk3

    ! finally we have the triple coordinates
    i = im1 + 1
    j = jm1 + 1
    k = km1 + 1

    ! compute the k point in the 1BZ
    kpt_in_1BZ(1) = dble(i-1)/dble(nk1)
    kpt_in_1BZ(2) = dble(j-1)/dble(nk2)
    kpt_in_1BZ(3) = dble(k-1)/dble(nk3)

  end subroutine find_k_1BZ_and_G


  function find_free_unit()
    !--------------------------------------------------------------------------
    ! This function finds a free input/output unit. It is copied from
    ! the QE distribution and slightly modified to suit our needs.
    !--------------------------------------------------------------------------
    !
    implicit none
    !
    integer :: find_free_unit
    integer :: io_unit
    logical :: opnd

    do io_unit = 99, 1, -1
      !
      inquire( unit = io_unit, opened = opnd )
      if ( .not. opnd ) then
        find_free_unit = io_unit
        return
      end if
      !
    end do
    !
    write(*,*) 'NO free units available!?!'
    !
    return
    !
  end function find_free_unit


  subroutine cryst_to_cart (nvec, vec, trmat, iflag)
    !-----------------------------------------------------------------------
    !
    ! This routine transforms the atomic positions or the k-point
    ! components from crystallographic to cartesian coordinates
    ! ( iflag=1 ) and viceversa ( iflag=-1 ).
    ! Output cartesian coordinates are stored in the input ('vec') array
    !
    !-----------------------------------------------------------------------
    !
    implicit none
    !
    integer, intent(in) :: nvec, iflag
    ! nvec:  number of vectors (atomic positions or k-points)
    !        to be transformed from crystal to cartesian and vice versa
    ! iflag: gives the direction of the transformation
    real(8), intent(in) :: trmat (3, 3)
    ! trmat: transformation matrix
    ! if iflag=1:
    !    trmat = at ,  basis of the real-space lattice,       for atoms   or
    !          = bg ,  basis of the reciprocal-space lattice, for k-points
    ! if iflag=-1: the opposite
    real(8), intent(inout) :: vec (3, nvec)
    ! coordinates of the vector (atomic positions or k-points) to be
    ! transformed - overwritten on output
    !
    !    local variables
    !
    integer :: nv, kpol
    ! counter on vectors
    ! counter on polarizations
    real(8) :: vau (3)
    ! workspace
    !
    !     Compute the cartesian coordinates of each vectors
    !     (atomic positions or k-points components)
    !
    do nv = 1, nvec
      if (iflag.eq.1) then
        do kpol = 1, 3
          vau (kpol) = trmat (kpol, 1) * vec (1, nv) + &
                       trmat (kpol, 2) * vec (2, nv) + &
                       trmat (kpol, 3) * vec (3, nv)
        enddo
      else
        do kpol = 1, 3
          vau (kpol) = trmat (1, kpol) * vec (1, nv) + &
                       trmat (2, kpol) * vec (2, nv) + &
                       trmat (3, kpol) * vec (3, nv)
        enddo
      endif
      do kpol = 1, 3
        vec (kpol, nv) = vau (kpol)
      enddo
    enddo
    !
    return
    !
  end subroutine cryst_to_cart


    SUBROUTINE HPSORT(N,RA,P)
    !------------------------------------------------------------
    ! subroutine which performs heap sort on a list of integers
    ! and also returns an array identifying the permutation
    ! which sorted the array.
    !
    ! adapted from Numerical Recipes pg. 329 (new edition)
    !
    !*****************************************************
    !*  Sorts an array RA of length N in ascending order *
    !*                by the Heapsort method             *
    !* ------------------------------------------------- *
    !* INPUTS:                                           *
    !*      N   size of table RA                         *
    !*      RA  table to be sorted                       *
    !* OUTPUT:                                           *
    !*      RA  table sorted in ascending order          *
    !*      P   table of indices showing transform       *
    !*                                                   *
    !* NOTE: The Heapsort method is a N Log N routine,   *
    !*       and can be used for very large arrays.      *
    !* ------------------------------------------------- *
    !* REFERENCE:                                        *
    !*  "NUMERICAL RECIPES by W.H. Press, B.P. Flannery, *
    !*   S.A. Teukolsky and W.T. Vetterling, Cambridge   *
    !*   University Press, 1986".                        *
    !*****************************************************

    implicit none

    integer, intent(in) :: N
    integer, intent(inout) :: RA(N)
    integer, intent(out) :: P(N)

    integer :: I, L, IR, RRA, PP, J

    ! nothing to order
    if (N < 2) return

    do I=1,N
      P(I) = I
    end do

    L=N/2+1
    IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
  10 continue
    if(L > 1)then
      L=L-1

      RRA=RA(L)
      PP =P(L)

    else
      RRA=RA(IR)
      PP =P (IR)

      RA(IR)=RA(1)
      P(IR) =P(1)

      IR=IR-1
      if(IR.eq.1)then

        RA(1)=RRA
        P (1)=PP

        return
      end if
    end if
    I=L
    J=L+L
  20 if(J.le.IR)then
    if(J < IR)then
      if(RA(J) < RA(J+1))  J=J+1
    end if
    if(RRA < RA(J))then
      RA(I)=RA(J)
      P (I)=P (J)

      I=J; J=J+J
    else
      J=IR+1
    end if
    goto 20
    end if
    RA(I)=RRA
    P (I)=PP
    goto 10

  END SUBROUTINE HPSORT


  SUBROUTINE HPSORT_real(N,RA,P)
    !------------------------------------------------------------
    ! subroutine which performs heap sort on a list of real numbers
    ! and also returns an array identifying the permutation
    ! which sorted the array.
    !
    ! adapted from Numerical Recipes pg. 329 (new edition)
    !
    !*****************************************************
    !*  Sorts an array RA of length N in ascending order *
    !*                by the Heapsort method             *
    !* ------------------------------------------------- *
    !* INPUTS:                                           *
    !*      N   size of table RA                         *
    !*      RA  table to be sorted                       *
    !* OUTPUT:                                           *
    !*      RA  table sorted in ascending order          *
    !*      P   table of indices showing transform       *
    !*                                                   *
    !* NOTE: The Heapsort method is a N Log N routine,   *
    !*       and can be used for very large arrays.      *
    !* ------------------------------------------------- *
    !* REFERENCE:                                        *
    !*  "NUMERICAL RECIPES by W.H. Press, B.P. Flannery, *
    !*   S.A. Teukolsky and W.T. Vetterling, Cambridge   *
    !*   University Press, 1986".                        *
    !*****************************************************

    implicit none

    integer, intent(in) :: N
    real(dp), intent(inout) :: RA(N)
    integer, intent(out) :: P(N)

    integer :: I, L, IR, RRA, PP, J


    ! nothing to order
    if (N < 2) return

    do I=1,N
      P(I) = I
    end do

    L=N/2+1
    IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
    10 continue
    if(L > 1)then
      L=L-1

      RRA=RA(L)
      PP =P(L)

    else
      RRA=RA(IR)
      PP =P (IR)

      RA(IR)=RA(1)
      P(IR) =P(1)

      IR=IR-1
      if(IR.eq.1)then

        RA(1)=RRA
        P (1)=PP

        return
      end if
    end if
    I=L
    J=L+L
  20 if(J.le.IR)then
    if(J < IR)then
      if(RA(J) < RA(J+1))  J=J+1
    end if
    if(RRA < RA(J))then
      RA(I)=RA(J)
      P (I)=P (J)

      I=J; J=J+J
    else
      J=IR+1
    end if
    goto 20
    end if
    RA(I)=RRA
    P (I)=PP
    goto 10
  END SUBROUTINE HPSORT_real


  subroutine hpsort_eps(n, ra, ind, eps)
    !---------------------------------------------------------------------
    ! sort an array ra(1:n) into ascending order using heapsort algorithm,
    ! and considering two elements being equal if their values differ
    ! for less than "eps".
    ! n is input, ra is replaced on output by its sorted rearrangement.
    ! create an index table (ind) by making an exchange in the index array
    ! whenever an exchange is made on the sorted data array (ra).
    ! in case of equal values in the data array (ra) the values in the
    ! index array (ind) are used to order the entries.
    ! if on input ind(1)  = 0 then indices are initialized in the routine,
    ! if on input ind(1) != 0 then indices are assumed to have been
    !                initialized before entering the routine and these
    !                indices are carried around during the sorting process
    !
    ! no work space needed !
    ! free us from machine-dependent sorting-routines !
    !
    ! adapted from Numerical Recipes pg. 329 (new edition)
    !

    implicit none

    !-input/output variables
    integer, intent(in) :: n
    integer, intent(inout) :: ind(*)
    real(dp), intent(inout) :: ra(*)
    real(dp), intent(in) :: eps

    !-local variables
    integer :: i, ir, j, l, iind
    real(dp) :: rra


    ! nothing to order
    if (n < 2) return

    ! initialize index array
    if (ind(1) == 0) then
      do i = 1, n
          ind(i) = i
      enddo
    endif

    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1

    ir = n

    sorting: do

      ! still in hiring phase
      if ( l > 1 ) then
        l    = l - 1
        rra  = ra(l)
        iind = ind(l)
        ! in retirement-promotion phase.
      else
        ! clear a space at the end of the array
        rra  = ra(ir)
        !
        iind = ind(ir)
        ! retire the top of the heap into it
        ra(ir) = ra(1)
        !
        ind(ir) = ind(1)
        ! decrease the size of the corporation
        ir = ir - 1
        ! done with the last promotion
        if ( ir == 1 ) then
            ! the least competent worker at all !
            ra(1)  = rra
            !
            ind(1) = iind
            exit sorting
        endif
      endif
      ! wheter in hiring or promotion phase, we
      i = l
      ! set up to place rra in its proper level
      j = l + l
      !
      do while ( j <= ir )
        if ( j < ir ) then
            ! compare to better underling
            if ( abs(ra(j)-ra(j+1)) >= eps ) then
              if (ra(j) < ra(j+1)) j = j + 1
            else
              ! this means ra(j) == ra(j+1) within tolerance
              if (ind(j) < ind(j + 1) ) j = j + 1
            endif
        endif
        ! demote rra
        if ( abs(rra - ra(j)) >= eps ) then
            if (rra < ra(j)) then
              ra(i) = ra(j)
              ind(i) = ind(j)
              i = j
              j = j + j
            else
              ! set j to terminate do-while loop
              j = ir + 1
            end if
        else
            !this means rra == ra(j) within tolerance
            ! demote rra
            if (iind < ind(j) ) then
              ra(i) = ra(j)
              ind(i) = ind(j)
              i = j
              j = j + j
            else
              ! set j to terminate do-while loop
              j = ir + 1
            endif
        end if
      enddo
      ra(i) = rra
      ind(i) = iind

    end do sorting
    !
  end subroutine hpsort_eps


  subroutine find_r_in_WS_cell(at,rvec_cryst,nr1,nr2,nr3,rvec_WS_cryst)
    !----------------------------------------------------------------------------!
    ! Given a real space vector in crystal coordinates rvec_cryst(3),
    ! which has coordinates on a nr1 x nr2 x nr3 grid,
    ! this subroutine finds the coordinates of the corresponding vector
    ! in the Wigner Seitz cell, namely the corresponding vector which is
    ! closest to the origin.
    !----------------------------------------------------------------------------!
    use intw_useful_constants, only: one, zero

    implicit none

    ! input
    real(dp), intent(in) :: at(3,3)                   ! the real space basis vectors
    real(dp), intent(in) :: rvec_cryst(3)             ! the real space vector
    integer, intent(in) :: nr1, nr2, nr3              ! the real mesh parameters

    ! output
    real(dp), intent(out) :: rvec_WS_cryst(3)         ! vector in the WS cell

    ! internal variables
    integer :: i, j, k                   ! the triple coordinates of k_1BZ
    integer :: Rlat(3)                   ! the translation vector


    real(dp) :: r_UC(3)                   ! temporary UC vector
    real(dp) :: r_WS(3)                   ! temporary WS vector
    real(dp) :: list_T(3)                 ! translation coefficients
    real(dp) :: T1, T2, T3
    integer :: it1, it2, it3

    real(dp) :: ai_dot_aj(3,3)            ! inner product of basis vectors
    real(dp) :: square_norm
    real(dp) :: minimum_square_norm

    ! First, find the coordinates of the vector in the unit cell coordinates
    ! Note that the subroutine below was initially designed for reciprocal
    ! space, but the relevant algorithm is the same.
    call find_k_1BZ_and_G(rvec_cryst,nr1,nr2,nr3,i,j,k,r_UC,Rlat)

    ai_dot_aj = matmul(at(:,:),transpose(at(:,:)))

    list_T(:) = (/ -one, zero, one /)

    ! initialize
    minimum_square_norm = dot_product(r_UC,matmul(ai_dot_aj(:,:),r_UC))
    rvec_WS_cryst(:)    = r_UC(:)

    do it1 = 1, 3
      T1  = list_T(it1)
      do it2 = 1, 3
        T2  = list_T(it2)
        do it3 = 1, 3
          T3  = list_T(it3)

          ! Define a possible vector in cartesian coordinates as
          ! r_WS = sum_{i=1}^3 rvec_US(i) a_i
          !
          ! its norm^2 is given by r_WS*r_WS

          r_WS(:)  =  r_UC(:) + (/T1,T2,T3/)

          square_norm = dot_product(r_WS,matmul(ai_dot_aj(:,:),r_WS))

          if ( square_norm < minimum_square_norm ) then
            minimum_square_norm  = square_norm
            rvec_WS_cryst(:)     = r_WS(:)

          end if

        end do !it3
      end do !it2
    end do !it1

  end subroutine find_r_in_WS_cell

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

  function weight_ph(x)

    implicit none

    real(dp), intent(in) :: x
    real(dp) :: weight_ph
    real(dp),parameter :: x0=0.002d0 !x0=0.0001d0 !x0=0.0004d0 !x0=0.0006d0 !x0=0.003
    real(dp),parameter :: sigma=0.00001d0 !sigma=0.000001d0 !sigma=0.0001d0 !sigma=0.00002d0 !sigma=0.0003

    weight_ph = 1/(1+exp((x0-x)/sigma))

  end function weight_ph

  subroutine errore (a,b,i)
    character(len=*), intent(in) :: a,b
    integer :: i
    write(*,*)"ERROR:", a,b,i
    stop
  end subroutine errore

  subroutine simpson (mesh, func, rab, asum)
    !-----------------------------------------------------------------------
    !
    !     simpson's rule integration. On input:
    !       mesh = mhe number of grid points (should be odd)
    !       func(i)= function to be integrated
    !       rab(i) = r(i) * dr(i)/di * di
    !     For the logarithmic grid not including r=0 :
    !       r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
    !     For the logarithmic grid including r=0 :
    !       r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
    !     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
    !     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    !
    implicit none

    integer, intent(in) :: mesh
    real(DP), intent(in) :: rab (mesh), func (mesh)
    real(DP), intent(out):: asum
    !
    real(DP) :: f1, f2, f3, r12
    integer :: i
    !
    !     routine assumes that mesh is an odd number so run check
    !     if ( mesh+1 - ( (mesh+1) / 2 ) * 2 .ne. 1 ) then
    !       write(*,*) '***error in subroutine radlg'
    !       write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1
    !       stop
    !     endif
    asum = 0.0d0
    r12 = 1.0d0 / 12.0d0
    f3 = func (1) * rab (1) * r12

    do i = 2, mesh - 1, 2
      f1 = f3
      f2 = func (i) * rab (i) * r12
      f3 = func (i + 1) * rab (i + 1) * r12
      asum = asum + 4.0d0 * f1 + 16.0d0 * f2 + 4.0d0 * f3
    enddo

    return

  end subroutine simpson

  function qe_erf (x)
    !---------------------------------------------------------------------
    !
    !     Error function - computed from the rational approximations of
    !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
    !
    !     for abs(x) le 0.47 erf is calculated directly
    !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
    !
    use kinds, only: DP

    implicit none

    real(DP), intent(in) :: x
    real(DP) :: x2, p1 (4), q1 (4)
    !real(DP), external :: qe_erfc
    real(DP) :: qe_erf
    data p1 / 2.426679552305318E2_DP, 2.197926161829415E1_DP, &
              6.996383488619136_DP,  -3.560984370181538E-2_DP /
    data q1 / 2.150588758698612E2_DP, 9.116490540451490E1_DP, &
              1.508279763040779E1_DP, 1.000000000000000_DP /
    !
    if (abs (x) > 6.0_DP) then
      !
      !  erf(6)=1-10^(-17) cannot be distinguished from 1
      !
      qe_erf = sign (1.0_DP, x)
    else
      if (abs (x)  <= 0.47_DP) then
        x2 = x**2
        qe_erf = x * (p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 (4) ) ) ) &
                   / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 (4) ) ) )
      else
        qe_erf = 1.0_DP - qe_erfc (x)
      endif
    endif
    !
    return

  end function qe_erf

  function qe_erfc (x)
    !---------------------------------------------------------------------
    !
    !     erfc(x) = 1-erf(x)  - See comments in erf
    !
    use kinds, only: DP

    implicit none

    real(DP),intent(in) :: x
    real(DP) :: qe_erfc
    real(DP) :: ax, x2, xm2, p2 (8), q2 (8), p3 (5), q3 (5), pim1
    !real(DP), external :: qe_erf
    data p2 / 3.004592610201616E2_DP,  4.519189537118719E2_DP, &
              3.393208167343437E2_DP,  1.529892850469404E2_DP, &
              4.316222722205674E1_DP,  7.211758250883094_DP,   &
              5.641955174789740E-1_DP,-1.368648573827167E-7_DP /
    data q2 / 3.004592609569833E2_DP,  7.909509253278980E2_DP, &
              9.313540948506096E2_DP,  6.389802644656312E2_DP, &
              2.775854447439876E2_DP,  7.700015293522947E1_DP, &
              1.278272731962942E1_DP,  1.000000000000000_DP /
    data p3 /-2.996107077035422E-3_DP,-4.947309106232507E-2_DP, &
             -2.269565935396869E-1_DP,-2.786613086096478E-1_DP, &
             -2.231924597341847E-2_DP /
    data q3 / 1.062092305284679E-2_DP, 1.913089261078298E-1_DP, &
              1.051675107067932_DP,    1.987332018171353_DP,    &
              1.000000000000000_DP /

    data pim1 / 0.56418958354775629_DP /
    !        ( pim1= sqrt(1/pi) )
    ax = abs (x)
    if (ax > 26.0_DP) then
      !
      !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
      !
      qe_erfc = 0.0_DP
    elseif (ax > 4.0_DP) then
      x2 = x**2
      xm2 = (1.0_DP / ax) **2
      qe_erfc = (1.0_DP / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
            + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
            ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
            (q3 (4) + xm2 * q3 (5) ) ) ) ) )
    elseif (ax > 0.47_DP) then
      x2 = x**2
      qe_erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
            + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
            + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
            (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
            (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
    else
      qe_erfc = 1.0_DP - qe_erf (ax)
    endif
    !
    ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
    !
    if (x < 0.0_DP) qe_erfc = 2.0_DP - qe_erfc
    !
    return

  end function qe_erfc


  subroutine gaussian(e0, e_min, e_max, n, sigma, gauss)
    !
    ! Created by Peio G. Goiricelaya 11/03/2016
    !
    !==========================================================================================!
    ! We create a list containing values for Delta[e-e0] simulated by the gaussian function    !
    ! with e defined by us from e_min to e_max, and with n points, each on separated from its  !
    ! neighbours by abs[e_max-e_min]/n                                                         !
    !==========================================================================================!

    use kinds, only: dp
    use intw_useful_constants, only: ZERO, pi

    implicit none

    !I/O variables

    integer, intent(in) :: n
    real(kind=dp), intent(in) :: e0, e_min, e_max, sigma
    real(kind=dp), intent(inout) :: gauss(n)

    !local variables

    integer :: i
    real(kind=dp) :: e

    gauss = ZERO
    !
    do i = 1, n
      !
      e = e_min+(e_max-e_min)*(i-1)/(n-1)
      !
      if ( abs(e-e0) > 4.d0*sigma ) then
        gauss(i) = ZERO
      else
        gauss(i) = exp(-((e-e0)/sigma)**2)/(sigma*sqrt(pi))
      endif
      !
    enddo

  end subroutine gaussian


  recursive function sphb(n,x) result(sphb_)
    !
    ! Spherical bessel functions
    !

    use intw_useful_constants, only: ZERO, ONE

    implicit none

    real(kind=dp), intent(in) :: x(:)
    integer, intent(in) :: n
    real(kind=dp) :: sphb_(size(x))

    integer :: i0


    if (any(x < ZERO)) stop "ERROR: sphb"

    do i0 = 1, size(x)
      if (abs(x(i0)) > epsilon(x(1))) exit
    enddo

    if (n==0) then
      sphb_(:i0-1) = ONE
      sphb_(i0:) = sin(x(i0:))/x(i0:)
    else if (n==1) then
      sphb_(:i0-1) = ZERO
      sphb_(i0:) = ( sin(x(i0:))/x(i0:) - cos(x(i0:)) ) / x(i0:)
    else if (n==2) then
      sphb_(:i0-1) = ZERO
      sphb_(i0:) = (3.0_dp/x(i0:)**2 - 1) * sin(x(i0:))/x(i0:) - 3.0_dp*cos(x(i0:))/x(i0:)**2
    else if (n==3) then
      sphb_(:i0-1) = ZERO
      sphb_(i0:) = (15.0_dp/x(i0:)**2 - 6.0_dp) * sin(x(i0:))/x(i0:)**2 - (15.0_dp/x(i0:)**2 - 1.0_dp) * cos(x(i0:))/x(i0:)
    else if (n>=4) then
      sphb_ = - sphb(n-2,x) + (2*n-1) * sphb(n-1,x)/x
    end if

  end function sphb


  subroutine real_ylmr2(lmax, ng, g, gg, ylm)
    !-----------------------------------------------------------------------
    ! Real spherical harmonics ylm(G) up to l=lmax
    ! lmax2 = (lmax+1)^2 is the total number of spherical harmonics
    ! Numerical recursive algorithm based on the one given in Numerical
    ! Recipes but avoiding the calculation of factorials that generate
    ! overflow for lmax > 11
    !-----------------------------------------------------------------------
    use kinds, only: dp
    use intw_useful_constants, only: pi, tpi, fpi
    !
    implicit none
    !
    integer, intent(in) :: lmax, ng
    real(kind=dp), intent(in) :: g(3,ng), gg(ng)
    real(kind=dp), intent(out) :: ylm(ng,(lmax+1)**2)
    !
    ! local variables
    real(kind=dp), parameter :: eps = 1.0d-9
    real(kind=dp) :: cost(ng), sent(ng), phi(ng)
    real(kind=dp) :: Q(ng,0:lmax,0:lmax)
    real(kind=dp) :: c, gmod
    integer :: lmax2, ig, l, m, lm


    if (ng < 1 .or. lmax < 1) return

    lmax2 = (lmax+1)**2

  10 continue

    if (lmax == 0) then
      ylm(:,1) =  sqrt(1.d0 / fpi)
      return
    end if
    !
    !  theta and phi are polar angles, cost = cos(theta)
    !
!$omp parallel default(none) &
!$omp shared(ng,g,gg,cost,phi,sent) &
!$omp private(ig,gmod)
!$omp do
    do ig = 1, ng
      gmod = sqrt(gg(ig))
      if (gmod < eps) then
        cost(ig) = 0.d0
      else
        cost(ig) = g(3,ig)/gmod
      endif
      !
      !  beware the arc tan, it is defined modulo pi
      !
      if (g(1,ig) > eps) then
        phi(ig) = atan( g(2,ig)/g(1,ig) )
      else if (g(1,ig) < -eps) then
        phi(ig) = atan( g(2,ig)/g(1,ig) ) + pi
      else
        phi(ig) = sign( pi/2.d0,g(2,ig) )
      end if
      sent(ig) = sqrt(max(0d0,1.d0-cost(ig)**2))
    enddo
    !
!$omp end parallel
    !
    !  Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
    !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
    !
    lm = 0
    do l = 0, lmax
      !
      c = sqrt(dble(2*l+1) / fpi)
      if ( l == 0 ) then
        Q(:,0,0) = 1.d0
      else if ( l == 1 ) then
        Q(:,1,0) = cost(:)
        Q(:,1,1) =-sent(:)/sqrt(2.d0)
      else
        !
        !  recursion on l for Q(:,l,m)
        !
        do m = 0, l - 2
          Q(:,l,m) = cost(:)*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(:,l-1,m) &
                    - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(:,l-2,m)
        end do
        Q(:,l,l-1) = cost(:) * sqrt(DBLE(2*l-1)) * Q(:,l-1,l-1)
        Q(:,l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent(:)*Q(:,l-1,l-1)
      end if
      !
      ! Y_lm, m = 0
      !
      lm = lm + 1
      ylm(:, lm) = c * Q(:,l,0)
      !
      do m = 1, l
        !
        ! Y_lm, m > 0
        !
        lm = lm + 1
        ylm(:, lm) = c * sqrt(2.d0) * Q(:,l,m) * cos (m*phi(:))
        !
        ! Y_lm, m < 0
        !
        lm = lm + 1
        ylm(:, lm) = c * sqrt(2.d0) * Q(:,l,m) * sin (m*phi(:))
      end do
      !
    end do

  end subroutine real_ylmr2


  ! subroutine complex_ylmr2(lmax, ng, g, gg, ylm)
  !   !-----------------------------------------------------------------------
  !   !     INTW project
  !   !     Complex spherical harmonics ylm(G) up to l=lmax
  !   !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !   !-----------------------------------------------------------------------
  !
  !   use kinds, only: dp
  !   use intw_useful_constants, only: pi, tpi, fpi
  !
  !   implicit none
  !
  !   integer, intent(in) :: lmax, ng
  !   real(kind=dp), intent(in) :: g(3,ng), gg(ng)
  !   complex(kind=dp), intent(out) :: ylm(ng,(lmax+1)**2)
  !
  !   real(kind=dp), parameter :: eps = 1.0d-9
  !   integer :: lmax2
  !
  !
  !   lmax2 =(lmax+1)**2
  !
  !   ylm=(0.0_dp,0.0_dp)
  !
  ! end subroutine complex_ylmr2



  ! MBR 24/04/24
  ! smearing functions for integrals

  function smeared_delta(x,s)
    use kinds, only: dp
    use intw_useful_constants, only: tpi
    implicit none
    real(dp) :: x,s, smeared_delta
    !gaussian
    smeared_delta = exp(-0.5_dp*(x/s)**2 ) / (s*sqrt(tpi))
  return
  end function smeared_delta

  function fermi_dirac(x, kt)
          ! x = energy - e_fermi, kT = Boltzmann * temp, same units
    use kinds, only: dp
    use intw_useful_constants, only: tpi
    implicit none
    real(dp) :: x, kt, fermi_dirac
    ! With finite T
    fermi_dirac = 1.0_dp/(exp(x/kt) + 1.0_dp)
    !  T=0 limit
    ! fermi_dirac = 0.5_dp* (-sign(1.0_dp,x) + 1.0_dp)
  return
  end function fermi_dirac


  ! MBR 24/04/24
  ! triangle area (taken from FSH/modules/geometry.f90)
  function  area_vec(v1,v2)
    use kinds, only: dp
    implicit none
    real(kind=dp), dimension(3), intent(in) :: v1,v2
    real(kind=dp) :: area_vec
    area_vec =  (v1(2)*v2(3) - v1(3)*v2(2))**2 + &
                (v1(3)*v2(1) - v1(1)*v2(3))**2 + &
                (v1(1)*v2(2) - v1(2)*v2(1))**2
    area_vec = sqrt(area_vec)*0.5_dp
  end function area_vec

end module intw_utility
