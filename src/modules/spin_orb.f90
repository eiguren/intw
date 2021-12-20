MODULE intw_spin_orb

  USE kinds, ONLY: DP
  USE intw_pseudo, ONLY : lmaxx

  SAVE
  !
  ! variables
  public :: starting_spin_angle, domag, rot_ylm, fcoef
  !
  private
  !
  LOGICAL :: &
!      lspinorb,  &       ! if .TRUE. this is a spin-orbit calculation
      starting_spin_angle, & ! if .TRUE. the initial wavefunctions are
                             ! spin-angle functions.
      domag              ! if .TRUE. magnetization is computed


  COMPLEX (DP) :: rot_ylm(2*lmaxx+1,2*lmaxx+1)  ! transform real
                         ! spherical harmonics into complex ones
  COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:) ! function needed to
                         ! account for spinors.
END MODULE intw_spin_orb
