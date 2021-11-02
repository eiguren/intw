MODULE intw_vlocal
  !
  ! ... The variables needed for the local potential in reciprocal space
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE :: &
       strf(:,:)              ! the structure factor
  REAL(DP), ALLOCATABLE :: &
       vloc(:,:)              ! the local potential for each atom type
  !
END MODULE intw_vlocal

