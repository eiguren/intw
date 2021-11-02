!
! Adapted from Quantum ESPRESSO group
!
MODULE intw_uspp_param
  !
  USE kinds,        ONLY : DP

  USE intw_pseudo2, ONLY : intw_pseudo
  !
  SAVE
  !
  TYPE (intw_pseudo),  ALLOCATABLE, TARGET :: upf(:)

  INTEGER, PARAMETER :: npsx=10

  INTEGER :: &
       nh(npsx),             &! number of beta functions per atomic type
       nhm,                  &! max number of different beta functions per atom
       nbetam                 ! max number of beta functions
  INTEGER :: &
       lmaxkb                ! max angular momentum
  !
END MODULE intw_uspp_param

MODULE intw_uspp
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: indv, nhtol, nhtolm, nkb,  &
       vkb, vkqb, dvan,  nhtoj, ijtoh, beta, becsum 
  PUBLIC ::  dvan_so
  INTEGER, PARAMETER :: &
       lmaxx  = 3      ! max non local angular momentum (l=0 to lmaxx)   
  !
  INTEGER :: nkb        ! total number of beta functions, with struct.fact.
  !
  INTEGER, ALLOCATABLE ::&
       indv(:,:),        &! indes linking  atomic beta's to beta's in the solid
       nhtol(:,:),       &! correspondence n <-> angular momentum l
       nhtolm(:,:),      &! correspondence n <-> combined lm index for (l,m)
       ijtoh(:,:,:)       ! correspondence beta indexes ih,jh -> composite index ijh
  !
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       vkb(:,:),   vkqb(:,:)             ! all beta functions in reciprocal space
  REAL(DP), ALLOCATABLE :: &
       becsum(:,:,:)           ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  REAL(DP), ALLOCATABLE :: &
       dvan(:,:,:),           &! the D functions of the solid
       nhtoj(:,:)              ! correspondence n <-> total angular momentum
  !
  COMPLEX(DP), ALLOCATABLE :: & ! variables for spin-orbit/noncolinear case:
       dvan_so(:,:,:,:)         ! D_{nm}
  !
  REAL(DP), ALLOCATABLE :: &
       beta(:,:,:)           ! beta functions for CP (without struct.factor)
  !
  !
END MODULE intw_uspp

