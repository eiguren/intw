module intw_fft_interp
  !
  !
  use kinds, only: dp
  !
  implicit none
  !
  ! subroutines
  public :: fft_interp_3d_real, nufft_interp_3d_cmplx
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
    use intw_reading, only: bg, tpiba2
    !
    implicit none

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
    integer, allocatable, dimension(:)     :: nl, nls      ! map entre g's centrados en origen y gs de FFT
    complex(kind=dp), allocatable, dimension(:) :: fg
    complex(kind=dp), allocatable, dimension(:) :: fft_dummy(:)
    integer  :: nri1, nrj2, nrk3, nri1s, nrj2s, nrk3s
    integer :: i, j, k, ig
    real(kind=dp), dimension(3) :: G
    integer, allocatable, dimension(:,:) :: gvec !vectores g centrados en origen
    integer  :: n1, n2, n3
    ! integer  :: ng
    integer :: ind(nr1,nr2,nr3), inds(nr1s,nr2s,nr3s)     ! indice de (i,j,k) 3D --> it 1D
    integer :: ngm, cnt
    !

    ! QE Gamma point ngm: 30227
    !
    !
    cnt=0
    do k=1,nr3
      do j=1,nr2
        do i=1,nr1
          cnt=cnt+1
          ind(i,j,k) = cnt ! 3D --> 1D
        enddo
      enddo
    enddo


    cnt=0
    do k=1,nr3s
      do j=1,nr2s
        do i=1,nr1s
          cnt=cnt+1
          inds(i,j,k) = cnt ! 3D --> 1D
        enddo
      enddo
    enddo
    !
    !
    nri1=int(nr1/2.0)
    nrj2=int(nr2/2.0)
    nrk3=int(nr3/2.0)
    !
    ngm=0
    do  i=-nri1,nri1! contar cuantos gs caen en Ecut
      do  j=-nrj2,nrj2
        do  k=-nrk3,nrk3
          !
          G = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
          if ( G(1)**2 + G(2)**2 + G(3)**2 <= ecut/tpiba2 ) then
            ngm=ngm+1
          end if
          !
        enddo
      enddo
    enddo
    !
    !
    nri1s=int(nr1s/2.0)
    nrj2s=int(nr2s/2.0)
    nrk3s=int(nr3s/2.0)
    !
    !
    !
    allocate( nl(ngm), nls(ngm) )
    allocate( fg (ngm) )
    allocate( gvec(3,ngm) )
    !
    ig=0
    do i=-nri1,nri1
      do j=-nrj2,nrj2
        do k=-nrk3,nrk3
          !
          G = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
          if ( G(1)**2 + G(2)**2 + G(3)**2 <= ecut/tpiba2 ) then
            ig=ig+1
            gvec(:,ig) = (/i,j,k/)
          endif
          !
        enddo
      enddo
    enddo
    !
    !
    !
    !
    do ig = 1, ngm
      !
      n1       = modulo((gvec(1,ig)), nr1) +1
      n2       = modulo((gvec(2,ig)), nr2) +1
      n3       = modulo((gvec(3,ig)), nr3) +1
      nl (ig) = ind ( n1, n2, n3)
      !
      n1      = modulo((gvec(1,ig)), nr1s) +1
      n2      = modulo((gvec(2,ig)), nr2s) +1
      n3      = modulo((gvec(3,ig)), nr3s) +1
      nls(ig) = inds(n1,n2,n3)
      !
    enddo
    !
    !
    allocate( fft_dummy(nr1*nr2*nr3) )
    fft_dummy = cmplx(fr,0.0)
    call cfftnd(3,(/nr1,nr2,nr3/),-1,fft_dummy)
    !
    do ig=1,ngm
      fg(ig)=fft_dummy(nl(ig))
    enddo
    deallocate( fft_dummy )
    !
    !
    ! open(unit=123,file="FFTfr.dat",status="replace")
    ! write(123,*) '#          g              FFT[f(x)]'
    ! do ig=1,ngm
    !   write(123,"(3i7,2x,2f15.10)") gvec(:,ig),fg(ig)
    ! enddo
    ! close(123)
    !
    !
    allocate( fft_dummy(nr1s*nr2s*nr3s) )
    fft_dummy  =  (0.0_dp,0.0_dp)
    do ig=1,ngm
         fft_dummy(nls(ig))=fg(ig)
    enddo
    !
    if ( mod(nr1,2)==0 .and. nr1/=nr1s ) then
      do j=1,nr2s
        do k=1,nr3s
          fft_dummy(inds(nr1/2+1,j,k))=0.0
        enddo
      enddo
    endif
    !
    if ( mod(nr2,2)==0 .and. nr2/=nr2s ) then
      do i=1,nr1s
        do k=1,nr3s
          fft_dummy(inds(i,nr2/2+1,k))=0.0
        enddo
      enddo
    endif
    !
    if ( mod(nr3,2)==0 .and. nr3/=nr3s ) then
      do i=1,nr1s
        do j=1,nr2s
          fft_dummy(inds(i,j,nr3/2+1))=0.0
        enddo
      enddo
    endif
    !
    !
    call cfftnd(3,(/nr1s,nr2s,nr3s/),+1,fft_dummy)
    !
    frs = real(fft_dummy)
    deallocate(fft_dummy)
    !
    !
    !
  end subroutine fft_interp_3d_real

    ! MBR 29/02/2024
    subroutine nufft_interp_3d_cmplx(nk1, nk2, nk3, kpoint, cfr, cfkG)

    ! module variables
    use intw_useful_constants, only : tpi
    use intw_reading, only: ngm, gvec,  nr1, nr2, nr3
    !
    implicit none

    external :: nufft3d2f90
    !
    ! input variables
    integer, intent(in)  :: nk1, nk2, nk3
    real(kind=dp), intent(in) :: kpoint(3) ! cryst coordinates in the conventional 1BZ
    complex(kind=dp), intent(in) :: cfr(nr1*nk1,nr2*nk2,nr3*nk3) ! original function in equispaced fractional r+R grid
    ! output variables
    complex(kind=dp), intent(out) :: cfkG(ngm) ! inverse FT at k+G for all G in gvec
    ! internal
    integer :: ierr, iG, G(3)
    real(kind=dp) :: kpG(3), kG1(ngm), kG2(ngm), kG3(ngm)
    !
    ! Prepare list of k+G vectors
    do iG = 1,ngm
      G = gvec(:,iG)
      !Fold G vector from -nr/2:nr/2 to 0:nr
      if (G(1)<0) G(1)=G(1)+nr1
      if (G(2)<0) G(2)=G(2)+nr2
      if (G(3)<0) G(3)=G(3)+nr3
      ! kpoint in 0:1 already. So kpG = (kpoint+G)/nr1 is in 0:1. Finally,
      ! k+G must be in -pi:pi range
      kpG = kpoint(:) + real(G(:),dp)
      kpG(1) =  kpG(1) / real(nr1,dp)
      kpG(2) =  kpG(2) / real(nr2,dp)
      kpG(3) =  kpG(3) / real(nr3,dp)
      if ( kpG(1) > 0.5_dp) kpG(1) = kpG(1) - 1.0_dp
      if ( kpG(2) > 0.5_dp) kpG(2) = kpG(2) - 1.0_dp
      if ( kpG(3) > 0.5_dp) kpG(3) = kpG(3) - 1.0_dp
      kpG = kpG * tpi
      ! Add to list of reciprocal 1,2,3 coordinates
      kG1(iG) = kpG(1)
      kG2(iG) = kpG(2)
      kG3(iG) = kpG(3)
    end do
    !
    call nufft3d2f90 (ngm, kG1, kG2, kG3, cfkG, -1, 1.0d-12, nr1*nk1,nr2*nk2,nr3*nk3, cfr, ierr)
    if (ierr .gt. 0) then
            print *, 'Error ', ierr, 'in nufft3d2f90. Stopping.'
            stop
    end if
    !
  return
  end  subroutine nufft_interp_3d_cmplx


end module intw_fft_interp
