!
MODULE intw_eqv
  USE kinds, ONLY :  DP
  !
  ! ... The wavefunctions at point k+q
  !
  SAVE
!haritz
  ! variables
  public :: evq, dvpsi, dpsi, drhoscfs, dmuxc, vlocq, eprec
  !
  private
!haritz
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  !
  ! ... The variable describing the linear response problem
  !
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:), drhoscfs (:,:,:)
  ! the product of dV psi
  ! the change of the wavefunctions
  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin),
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the derivative of the xc potential
  ! the local potential at q+G
  REAL (DP), ALLOCATABLE :: eprec(:,:) ! needed for preconditioning
  !
END MODULE intw_eqv
!
!
!
MODULE intw_phus
  USE kinds, ONLY :  DP
  USE intw_becmod, ONLY : bec_type
  !
  ! ... These are additional variables needed for the linear response
  ! ... program with the US pseudopotentials
  !
  SAVE
!haritz
  !variables
  public :: alphasum, dpqq, int1, int2, int3, int3_paw, int4, int5, &
            int1_nc, int2_so, int3_nc, int4_nc, int5_so, becsum_nc, &
            becsumort, alphasum_nc, dpqq_so,becp1, alphap

  private
!haritz
  !
  REAL (DP), ALLOCATABLE :: &
       alphasum(:,:,:,:),   &! nhm*(nhm+1)/2,3,nat,nspin)
                             ! used to compute modes
       dpqq(:,:,:,:)         ! (nhm, nhm, 3, ntyp)
  ! alphasum contains \sum_i <psi_i| d/du (|\beta_n><beta_m|) | psi_i> + (m-n)
  ! dipole moment of each Q
  COMPLEX (DP), ALLOCATABLE :: &
       int1(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2(:,:,:,:,:),     &! nhm, nhm, 3,nat, nat),&
       int3(:,:,:,:,:),     &! nhm, nhm, npert, nat, nspin),&
       int3_paw(:,:,:,:,:), &! nhm, nhm, npert, nat, nspin),&
       int4(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nspin),&
       int5(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nat),&
       int1_nc(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2_so(:,:,:,:,:,:),   &! nhm, nhm, 3, nat,nat,nspin),&
       int3_nc(:,:,:,:,:),     &! nhm, nhm, npert, nat, nspin),&
       int4_nc(:,:,:,:,:,:),   &! nhm, nhm, 3, 3, nat, nspin),&
       int5_so(:,:,:,:,:,:,:), &! nhm*(nhm+1)/2, 3, 3, nat, nat, nspin),&
!
!  These variables contains the five integrals defined in PRB 64, 35118 (2001)
!  int1 -> \int V_eff d/du (Q) d^3r
!  int2 -> \int d/du (V_loc) Q d^3r
!  int3 -> \int d\du (V_Hxc) Q d^3r
!  int4 -> \int V_eff d^2/dudu (Q) d^3r
!  int5 -> \int d/du (V_loc) d/du (Q) d^3r
!
!  int3_paw contains d/du (D^1-\tilde D^1)
!
!
       becsum_nc(:,:,:,:),     &! nhm*(nhm+1)/2,nat,npol,npol)
       becsumort(:,:,:,:),     &! nhm*(nhm+1)/2,nat,nspin,3*nat)
       alphasum_nc(:,:,:,:,:), &! nhm*(nhm+1)/2,3,nat,npol,npol)
       dpqq_so(:,:,:,:,:)       ! nhm, nhm, nspin, 3, ntyp
!
!  becsum contains \sum_i <\psi_i | \beta_n><\beta_m| \psi_i > + (m-n)
!  besumort contains alphasum+\sum_i <\psi_i | \beta_n><\beta_m| \delta \psi_i >
!  dpqq_so dipole moment of each Q multiplied by the fcoef factors
!
  type (bec_type),  ALLOCATABLE, TARGET :: &
       becp1(:)              ! (nksq); (nkbtot, nbnd)
  !
  ! becp1 contains < beta_n | \psi_i >
  !
  type (bec_type),  ALLOCATABLE, TARGET :: &
       alphap(:,:)           ! nkbtot, nbnd, 3, nksq)
  !
  ! alphap contains < d\du (\beta_n) | psi_i>
  !
END MODULE intw_phus
!
!

