module intw_allwfcs
  ! 1. Modulo honetan uhin funtzioak irakurri transformatu eta berriro inv-transformatuko ditugu igk, igkq lista
  ! berdinak erabiliz. Era horretan sinplifikazio nabarmenta lortuko dugu.
  ! 2. Uhin funtzio guztiak (irreduzibleak) memorian gordeko ditugu, ordena magnitude batzuk azkarrago dena.
  !    Irreduzibleetatik lortuko digu denak baina azken hauek memorian gordeta izanik.
  ! 3. dv potentzial irreduzibleak ere RAM memorian dauzkagu.
  use kinds, only: dp

  implicit none

  !
  ! variables
  public :: wfc_k_irr_all, list_iG_all, QE_eig_irr_all,ngk_all
  !
  ! subroutines
  public :: allocate_and_get_all_irreducible_wfc, get_psi_general_k_all_wfc
  !
  private

  complex(dp), allocatable :: wfc_k_irr_all(:,:,:,:)
  integer, allocatable :: list_iG_all(:,:), ngk_all(:)
  real(dp), allocatable :: QE_eig_irr_all(:,:)

contains
  !
  subroutine allocate_and_get_all_irreducible_wfc()

    use intw_reading, only: nG_max, nkpoints_QE, get_K_folder_data_with_nG, nspin, nbands
    use intw_useful_constants, only: cmplx_0

    implicit none

    integer, allocatable :: list_iG(:)
    real(dp), allocatable :: QE_eig(:)
    complex(dp), allocatable :: wfc_g(:,:,:)
    integer :: i_folder, ipol, iG, ibnd,  ngk


    ! allocate what is useful
    !

    if (allocated(list_iG)) deallocate (list_iG)
    allocate   (list_iG(nG_max))
    !
    if (allocated(QE_eig)) deallocate (QE_eig)
    allocate   (QE_eig(nbands))
    !
    if (allocated(wfc_g)) deallocate (wfc_g)
    allocate   (wfc_g(nG_max,nbands,nspin))
    !
    if (allocated(wfc_k_irr_all)) deallocate (wfc_k_irr_all)
    allocate (wfc_k_irr_all(nkpoints_QE,nG_max,nbands,nspin))
    !
    if (allocated(QE_eig_irr_all)) deallocate (QE_eig_irr_all)
    allocate (QE_eig_irr_all(nkpoints_QE,nbands))
    !
    if (allocated(list_iG_all)) deallocate(list_iG_all)
    allocate (list_iG_all(nkpoints_QE,nG_max))

    if (allocated(ngk_all)) deallocate(ngk_all)
    allocate (ngk_all(nkpoints_QE))
    !
    wfc_k_irr_all=cmplx_0
    list_iG_all=0

    do i_folder=1,nkpoints_QE
      !
      call get_K_folder_data_with_nG(i_folder,list_iG,wfc_g,QE_eig, ngk_all(i_folder))
      !
      QE_eig_irr_all(i_folder,1:nbands)=QE_eig(1:nbands)
      !
      list_iG_all(i_folder,1:ngk_all(i_folder)  )=list_iG(1:ngk_all(i_folder))
      !
      do iG=1,ngk_all(i_folder)
        do ipol=1,nspin
          do ibnd=1,nbands
            !
            wfc_k_irr_all(i_folder,iG,ibnd,ipol)=wfc_g(iG,ibnd,ipol)
            !
          enddo ! ibnd
        enddo !ipol
      enddo !iG
      !
    enddo ! i_folder
    !
    return

  end subroutine allocate_and_get_all_irreducible_wfc

  subroutine get_psi_general_k_all_wfc(add_G,kpoint,list_iG,wfc_k,QE_eig,G_plus) !, ngk, igk, gk )

    use intw_reading, only: s, ftau, nG_max, nspin
    use intw_input_parameters, only: nk1, nk2, nk3
    use intw_symmetries, only: sym_G, full_mesh, inverse_indices, symlink, &
                               QE_folder_sym, apply_TR_to_wfc, rotate_wfc_test
    use intw_utility, only: find_k_1BZ_and_G, switch_indices
    use intw_fft, only: wfc_by_expigr
    use w90_parameters, only: num_bands

    implicit none

    !I/O variables

    logical, intent(in) :: add_G
    real(dp), intent(in) :: kpoint(3)
    integer, intent(out) :: list_iG(nG_max)
    real(dp), intent(out) :: QE_eig(num_bands)
    complex(dp), intent(out) :: wfc_k(nG_max,num_bands,nspin)
    integer, intent(out) :: G_plus(3)


    !local variables

    integer :: ikpt, i_folder
    integer :: i_sym, TR
    integer :: i_1bz, j_1bz, k_1bz
    integer :: sym(3,3), G_sym(3)
    real(dp) :: ftau_sym(3)
    real(dp) :: kpoint_1BZ(3)
    integer :: list_iG_irr(nG_max)
    complex(dp) :: wfc_k_irr(nG_max,num_bands,nspin)


    call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i_1bz,j_1bz,k_1bz,kpoint_1bz,G_plus)
    !
    call switch_indices(nk1,nk2,nk3,ikpt,i_1bz,j_1bz,k_1bz,+1)
    !
    if (full_mesh) then
      !
      ! Use the full BZ, no symmetry!
      !
      wfc_k(:,:,:)=wfc_k_irr_all(ikpt,:,:,:)
      list_iG(:)=list_iG_all(ikpt,:)
      QE_eig(:)=QE_eig_irr_all(ikpt,:)
      !
    else
      !
      ! Use the IBZ and symmetry!
      !
      ! identify the right folder
      !
      i_folder=QE_folder_sym(ikpt)
      !
      ! The symmetry which takes k_irr to k
      !
      i_sym    = symlink(ikpt,1)
      TR       = symlink(ikpt,2)
      G_sym    = sym_G(:,ikpt) !
      ftau_sym = ftau(:,i_sym)
      sym      = s(:,:,inverse_indices(i_sym))
      !
      ! debugging test
      !        call rotate_wfc (wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
      !			inverse_indices(i_sym), sym, ftau_sym, G_sym)
      !
      wfc_k_irr(:,:,:)=wfc_k_irr_all(i_folder,:,:,:)
      list_iG_irr(:)=list_iG_all(i_folder,:)
      QE_eig(:)=QE_eig_irr_all(i_folder,:)
      !
      call rotate_wfc_test (wfc_k_irr,list_iG_irr,wfc_k,list_iG,i_sym, &
                                                      sym,ftau_sym,G_sym)
      !
      ! If time-reversal is present, the wavefunction currently stored
      ! in wfc_k is actually for (-k). Complex conjugation must now
      ! be applied to recover wfc_k.
      !
      if (TR==1) then
        !
        call apply_TR_to_wfc(wfc_k,list_iG)
        !
      endif
      !
      list_iG_irr=list_iG
      !
      if (add_G) then
        !
        call wfc_by_expigr(kpoint,num_bands,nspin,ng_max,list_iG_irr,list_iG,wfc_k,G_plus)
        !
      endif
      !
    endif
    !
    return

  end subroutine  get_psi_general_k_all_wfc

!---------------------------------------------
end module intw_allwfcs
!---------------------------------------------
