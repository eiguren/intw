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
module intw_fft_interp
  ! Haritz 29/11/2024: This module is not used anywhere and it is not tested, so it could be removed.

  ! TODO: Add a description.

  use kinds, only: dp
  !
  implicit none
  !
  ! subroutines
  public :: fft_interp_3d_real
  !
  private

contains

  subroutine fft_interp_3d_real(nr1, nr2, nr3, nr1s, nr2s, nr3s, ecut, fr, frs)
    !
    ! This subroutine interpolates a real function (charge density, ...) from a
    ! (nr1,nr2,nr3) coarse FFT grid to a (nr1s,nr2s,nr3s) finer grid using
    ! Fourier interpolation.
    !
    ! module variables
    use intw_useful_constants, only: cmplx_0
    use intw_reading, only: bg, tpiba2
    !
    implicit none
    !
    external :: cfftnd
    !
    ! input variables
    integer, intent(in)  :: nr1, nr2, nr3,&
                            nr1s, nr2s, nr3s ! real space discretization
    real(kind=dp), intent(in) :: fr(nr1*nr2*nr3) ! original function
    real(kind=dp), intent(in) :: ecut
    !
    ! output variables
    real(kind=dp), intent(out) :: frs(nr1s*nr2s*nr3s) ! interpolated function
    !
    ! local variables
    integer, allocatable :: nl(:), nls(:) ! map between origin centered G's and FFT G's
    complex(kind=dp), allocatable :: fg(:)
    complex(kind=dp), allocatable :: fft_dummy(:)
    integer  :: nri1, nrj2, nrk3
    integer :: i, j, k, ig
    real(kind=dp) :: G(3)
    integer, allocatable :: gvec(:,:) ! origin centered G vectors
    integer  :: n1, n2, n3
    integer :: ind(nr1,nr2,nr3), inds(nr1s,nr2s,nr3s) ! (i,j,k) 3D --> it 1D mappings
    integer :: nG, cnt
    !
    !
    ! Initialize 3D --> 1D mappings
    !
    cnt = 0
    do k = 1, nr3
      do j = 1, nr2
        do i = 1, nr1
          cnt = cnt + 1
          ind(i,j,k) = cnt ! 3D --> 1D
        enddo
      enddo
    enddo
    !
    cnt = 0
    do k = 1, nr3s
      do j = 1, nr2s
        do i = 1, nr1s
          cnt = cnt + 1
          inds(i,j,k) = cnt ! 3D --> 1D
        enddo
      enddo
    enddo
    !
    !
    ! Count G vectors inside Ecut
    !
    nri1 = int(nr1/2.0_dp)
    nrj2 = int(nr2/2.0_dp)
    nrk3 = int(nr3/2.0_dp)
    !
    nG = 0
    do  i = -nri1, nri1
      do  j = -nrj2, nrj2
        do  k = -nrk3, nrk3
          !
          G = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
          if ( G(1)**2 + G(2)**2 + G(3)**2 <= ecut/tpiba2 ) then
            nG = nG + 1
          end if
          !
        enddo
      enddo
    enddo
    !
    !
    ! Initialize G vectors
    !
    allocate( gvec(3,nG) )
    !
    ig = 0
    do i = -nri1, nri1
      do j = -nrj2, nrj2
        do k = -nrk3, nrk3
          !
          G = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
          if ( G(1)**2 + G(2)**2 + G(3)**2 <= ecut/tpiba2 ) then
            ig = ig + 1
            gvec(:,ig) = (/i, j, k/)
          endif
          !
        enddo
      enddo
    enddo
    !
    !
    ! Initialize G vector mappings
    !
    allocate( nl(nG), nls(nG) )
    allocate( fg(nG) )
    !
    do ig = 1, nG
      !
      n1      = modulo((gvec(1,ig)), nr1) + 1
      n2      = modulo((gvec(2,ig)), nr2) + 1
      n3      = modulo((gvec(3,ig)), nr3) + 1
      nl (ig) = ind (n1,n2,n3)
      !
      n1      = modulo((gvec(1,ig)), nr1s) + 1
      n2      = modulo((gvec(2,ig)), nr2s) + 1
      n3      = modulo((gvec(3,ig)), nr3s) + 1
      nls(ig) = inds(n1,n2,n3)
      !
    enddo
    !
    !
    ! Fourier transform fr to G space
    !
    allocate( fft_dummy(nr1*nr2*nr3) )
    !
    fft_dummy = cmplx(fr, 0.0_dp, kind=dp)
    !
    call cfftnd(3, (/nr1, nr2, nr3/), -1, fft_dummy)
    !
    do ig = 1, nG
      fg(ig) = fft_dummy(nl(ig))
    enddo
    !
    deallocate( fft_dummy )
    !
    !
    ! Obtain Fourier transform of frs from fg
    !
    allocate( fft_dummy(nr1s*nr2s*nr3s) )
    !
    fft_dummy = cmplx_0
    !
    do ig = 1, nG
      fft_dummy(nls(ig)) = fg(ig)
    enddo
    !
    if ( mod(nr1, 2) == 0 .and. nr1 /= nr1s ) then
      do j = 1, nr2s
        do k = 1, nr3s
          fft_dummy(inds(nr1/2+1,j,k)) = 0.0_dp
        enddo
      enddo
    endif
    !
    if ( mod(nr2, 2) == 0 .and. nr2 /= nr2s ) then
      do i = 1, nr1s
        do k = 1, nr3s
          fft_dummy(inds(i,nr2/2+1,k)) = 0.0_dp
        enddo
      enddo
    endif
    !
    if ( mod(nr3, 2) == 0 .and. nr3 /= nr3s ) then
      do i = 1, nr1s
        do j = 1, nr2s
          fft_dummy(inds(i,j,nr3/2+1)) = 0.0_dp
        enddo
      enddo
    endif
    !
    !
    ! Fourier transform frs to real space
    !
    call cfftnd(3, (/nr1s, nr2s, nr3s/), +1, fft_dummy)
    !
    frs = real(fft_dummy)
    !
  end subroutine fft_interp_3d_real

end module intw_fft_interp
