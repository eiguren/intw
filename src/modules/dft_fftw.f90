!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
!
module intw_dft_fftw
  !
  !----------------------------------------------------------------------------!
  !
  !       The subroutines in this module handle the g -> r and
  !       r -> g tasks using FFTW.
  !
  !----------------------------------------------------------------------------!
  !
  use kinds, only : dp
  !
  implicit none

  include 'fftw3.f'

  ! variables
  public :: fftw_forward_3D_plan, fftw_backward_3D_plan, fftw_3D_in, fftw_3D_out
  !
  ! subroutines
  public :: setup_fftw_3D, cleanup_fftw_3D, &
            func_from_r_to_g_3D_FFTW, &
            func_from_g_to_r_3D_FFTW, &
            func_from_1D_to_3D_FFTW
  !
  private

  integer*8     ::  fftw_forward_3D_plan
  integer*8     ::  fftw_backward_3D_plan

  complex(dp), allocatable :: fftw_3D_in  (:,:,:)
  complex(dp), allocatable :: fftw_3D_out (:,:,:)
  !

contains

  subroutine setup_fftw_3D(nr1,nr2,nr3)
    !--------------------------------------------------------
    !  This subroutine is a driver which uses the 3D-FFT
    !--------------------------------------------------------

    implicit none

    external :: DFFTW_PLAN_DFT_3D

    integer	::	nr1, nr2, nr3


    allocate(fftw_3D_in(nr1,nr2,nr3))
    allocate(fftw_3D_out(nr1,nr2,nr3))

    ! give me  plans!
    call DFFTW_PLAN_DFT_3D(fftw_forward_3D_plan,nr1,nr2,nr3,     &
  			fftw_3D_in,fftw_3D_out, &
  			FFTW_FORWARD, FFTW_MEASURE)

    call DFFTW_PLAN_DFT_3D(fftw_backward_3D_plan,nr1,nr2,nr3,     &
  			fftw_3D_in,fftw_3D_out, &
  			FFTW_BACKWARD, FFTW_MEASURE)


  end subroutine setup_fftw_3D


  subroutine cleanup_fftw_3D()
    !--------------------------------------------------------
    !  This subroutine is a driver which uses the 3D-FFT
    !--------------------------------------------------------

    implicit none

    external :: dfftw_destroy_plan

    deallocate(fftw_3D_in)
    deallocate(fftw_3D_out)

    call dfftw_destroy_plan(fftw_forward_3D_plan)
    call dfftw_destroy_plan(fftw_backward_3D_plan)

  end subroutine cleanup_fftw_3D


  subroutine func_from_r_to_g_3D_FFTW(nr1,nr2,nr3,fr,fg)
    !----------------------------------------------------------------------------!
    !     This subroutine uses FFTW to compute the Fourier transform
    !     of a function of space in 3D.
    !----------------------------------------------------------------------------!

    implicit none

    external :: DFFTW_EXECUTE

    ! input variables
    integer     :: nr1, nr2, nr3
    complex(dp) :: fr(nr1,nr2,nr3)

    ! output variables
    complex(dp) :: fg(nr1,nr2,nr3)


    fftw_3D_in(:,:,:) = fr(:,:,:)

    CALL DFFTW_EXECUTE(fftw_forward_3D_plan)

    fg(:,:,:) =  fftw_3D_out(:,:,:)/dble(nr1*nr2*nr3)

  end subroutine func_from_r_to_g_3D_FFTW

  subroutine func_from_g_to_r_3D_FFTW(nr1,nr2,nr3,fr,fg)
    !----------------------------------------------------------------------------!
    !     This subroutine uses FFTW to compute the Fourier transform
    !     of a function of space in 3D.
    !----------------------------------------------------------------------------!

    implicit none

    external :: DFFTW_EXECUTE

    ! input variables
    integer     :: nr1, nr2, nr3
    complex(dp) :: fg(nr1,nr2,nr3)

    ! output variables
    complex(dp) :: fr(nr1,nr2,nr3)


    fftw_3D_in(:,:,:) = fg(:,:,:)

    CALL DFFTW_EXECUTE(fftw_backward_3D_plan)

    fr(:,:,:) =  fftw_3D_out(:,:,:)

  end subroutine func_from_g_to_r_3D_FFTW


  subroutine func_from_1D_to_3D_FFTW(nr1,nr2,nr3,f_1D,f_3D)
    !----------------------------------------------------------------------------!
    !     This subroutine takes a function in 1D and puts in in a 3D array,
    !     using switch_indices to store the data.
    !----------------------------------------------------------------------------!
    use intw_utility, only: switch_indices

    implicit none

    ! input variables
    integer     :: nr1, nr2, nr3
    complex(dp) :: f_1D(nr1*nr2*nr3)

    ! output variables
    complex(dp) :: f_3D(nr1,nr2,nr3)

    ! computation variables

    integer     :: i_singlet, i1,i2,i3

    integer     :: switch_singlet_to_triplet


    switch_singlet_to_triplet = -1

    do i_singlet = 1, nr1*nr2*nr3
  	  call switch_indices(nr1,nr2,nr3,i_singlet,i1,i2,i3,switch_singlet_to_triplet)
      f_3D(i1,i2,i3) = f_1D(i_singlet)
    end do

  end subroutine func_from_1D_to_3D_FFTW



!----------------------------------------------------------------------------!
!
!
  end module intw_dft_fftw
!
!
!----------------------------------------------------------------------------!
