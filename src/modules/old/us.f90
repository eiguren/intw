MODULE intw_us
  !
  ! ... These parameters are needed with the US pseudopotentials
  !
  USE kinds,      ONLY : DP
  !
  SAVE
  !
  INTEGER :: &
       nqxq,             &! size of interpolation table
       nqx                ! number of interpolation points
  REAL(DP), PARAMETER:: &
       dq = 0.01D0           ! space between points in the pseudopotential tab.
  REAL(DP), ALLOCATABLE :: &
       qrad(:,:,:,:),         &! radial FT of Q functions
       tab(:,:,:),            &! interpolation table for PPs
       tab_at(:,:,:)           ! interpolation table for atomic wfc
  LOGICAL :: spline_ps = .false.
  REAL(DP), ALLOCATABLE :: &
       tab_d2y(:,:,:)            ! for cubic splines
  !
END MODULE intw_us

