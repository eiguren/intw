module intw_fft_interp
  !
  !
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


end module intw_fft_interp
