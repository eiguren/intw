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

    use intw_reading, only: nG_max, nkpoints_QE, get_K_folder_data, nspin, num_bands_intw
    use intw_useful_constants, only: cmplx_0

    implicit none

    integer, allocatable :: list_iG(:)
    real(dp), allocatable :: QE_eig(:)
    complex(dp), allocatable :: wfc_g(:,:,:)
    integer :: i_folder, ispin, iG, ibnd


    ! allocate what is useful
    !
    if (allocated(list_iG)) deallocate (list_iG)
    allocate   (list_iG(nG_max))
    !
    if (allocated(QE_eig)) deallocate (QE_eig)
    !allocate   (QE_eig(nbands))
    allocate   (QE_eig(num_bands_intw))
    !
    if (allocated(wfc_g)) deallocate (wfc_g)
    !allocate   (wfc_g(nG_max,nbands,nspin))
    allocate   (wfc_g(nG_max,num_bands_intw,nspin))
    !
    if (allocated(wfc_k_irr_all)) deallocate (wfc_k_irr_all)
    !allocate (wfc_k_irr_all(nkpoints_QE,nG_max,nbands,nspin))
    allocate (wfc_k_irr_all(nkpoints_QE,nG_max,num_bands_intw,nspin))
    !
    if (allocated(QE_eig_irr_all)) deallocate (QE_eig_irr_all)
    !allocate (QE_eig_irr_all(nkpoints_QE,nbands))
    allocate (QE_eig_irr_all(nkpoints_QE,num_bands_intw))
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
      call get_K_folder_data(i_folder,list_iG,wfc_g,QE_eig, ngk_all(i_folder))
      !
      !QE_eig_irr_all(i_folder,1:nbands)=QE_eig(1:nbands)
      QE_eig_irr_all(i_folder,:)=QE_eig(:)
      !
      list_iG_all(i_folder,1:ngk_all(i_folder)  )=list_iG(1:ngk_all(i_folder))
      !
      do iG=1,ngk_all(i_folder)
        do ispin=1,nspin
          !do ibnd=1,nbands
          do ibnd=1, num_bands_intw
            !
            wfc_k_irr_all(i_folder,iG,ibnd,ispin)=wfc_g(iG,ibnd,ispin)
            !
          enddo ! ibnd
        enddo !ispin
      enddo !iG
      !
    enddo ! i_folder
    !
    return

  end subroutine allocate_and_get_all_irreducible_wfc


  subroutine get_psi_general_k_all_wfc(kpoint, ngk, list_iG, wfc_k, QE_eig)

    use intw_reading, only: s, ftau, nG_max, nspin, kpoints_QE, num_bands_intw
    use intw_input_parameters, only: nk1, nk2, nk3
    use intw_symmetries, only: full_mesh, symlink, QE_folder_sym, QE_folder_nosym, &
                               apply_TR_to_wfc, rotate_wfc_test
    use intw_utility, only: find_k_1BZ_and_G, triple_to_joint_index_g
    use intw_fft, only: wfc_by_expigr
    use intw_useful_constants, only : eps_5, ZERO

    implicit none

    !I/O variables
    real(dp), intent(in) :: kpoint(3)
    integer, intent(out) :: ngk, list_iG(nG_max)
    complex(dp), intent(out) :: wfc_k(nG_max,num_bands_intw,nspin)
    real(dp), intent(out), optional :: QE_eig(num_bands_intw)

    !local variables
    integer :: ikpt, i_folder
    integer :: i_sym, TR
    integer :: i_1bz, j_1bz, k_1bz
    integer :: sym(3,3), G_sym(3), G_plus(3)
    real(dp) :: ftau_sym(3)
    real(dp) :: kpoint_1BZ(3), k_QE(3), ktest(3)
    integer :: list_iG_irr(nG_max)
    complex(dp) :: wfc_k_irr(nG_max,num_bands_intw,nspin)


    call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i_1bz,j_1bz,k_1bz,kpoint_1bz,G_plus)
    !
    call triple_to_joint_index_g(nk1,nk2,nk3,ikpt,i_1bz,j_1bz,k_1bz)
    !
    if (full_mesh) then
      !
      ! Use the full BZ, no symmetry!
      !
      i_folder = QE_folder_nosym(ikpt)
      k_QE = kpoints_QE(:,i_folder)
      !
      ngk = ngk_all(ikpt)
      list_iG(:) = list_iG_all(ikpt,:)
      wfc_k(:,:,:) = wfc_k_irr_all(ikpt,:,:,:)
      if (present(QE_eig)) QE_eig(:) = QE_eig_irr_all(ikpt,:)
      !
      G_sym = kpoint - k_QE
      !
      call wfc_by_expigr(num_bands_intw, nspin, G_sym, list_iG, wfc_k)
      !
    else
      !
      ! Use the IBZ and symmetry!
      !
      ! This is the irreducible point
      ! connected to aimed kpoint
      i_folder = QE_folder_sym(ikpt)
      k_QE = kpoints_QE(:,i_folder)
      !
      ! The symmetry which takes kpoints_QE(:,i_folder) into aimed kpoint.
      ! sym * kpoints_QE =  kpoint
      i_sym = symlink(ikpt,1)
      TR = symlink(ikpt,2)
      ftau_sym = ftau(:,(i_sym))
      sym = s(:,:,(i_sym))
      !
      ! Load the corresponding irreducible wfcs in kpoints_QE
      ngk = ngk_all(i_folder)
      list_iG_irr(:) = list_iG_all(i_folder,:)
      wfc_k_irr(:,:,:) = wfc_k_irr_all(i_folder,:,:,:)
      if (present(QE_eig)) QE_eig(:) = QE_eig_irr_all(i_folder,:)
      !
      ktest = matmul(sym ,k_QE)
      !
      if (TR==1) then
        ! If TR needed
        ! G_sym  = aimed_kpoint -TR[S*QE_kpoint] = aimed_kpoint + S*QE_kpoint, such that
        ! aimed_kpoint = -S*QE_kpoint + G_ sym
        G_sym = nint(kpoint + ktest)
      else
        ! G_sym  = aimed_kpoint -   S*QE_kpoint, such that
        ! aimed_kpoint =  S*QE_kpoint + G_ sym
        G_sym = nint(kpoint - ktest)
      end if
      !
      call rotate_wfc_test(wfc_k_irr,list_iG_irr,wfc_k,list_iG,i_sym,sym,ftau_sym,(/0,0,0/))
      !
      ! If time-reversal is present, the wavefunction currently stored
      ! in wfc_k is actually for (-k). Complex conjugation must now
      ! be applied to recover wfc_k.
      !
      if (TR==1) then
        ktest = -ktest
        !
        call apply_TR_to_wfc(wfc_k,list_iG)
        !
      endif
      !
      if ( sum( abs( ktest + dble(G_sym) - kpoint ) )>eps_5 ) then
          write(*,*) "ERROR in get_psi_general_k_all_wfc:"
          write(*,"(a,100f12.6)") "Aimed kpoint is     :", kpoint
          write(*,"(a,100f12.6)") "Symmetry induced is :", ktest + G_sym
      end if
      !
      call wfc_by_expigr(num_bands_intw, nspin, G_sym, list_iG, wfc_k)
      !
    endif

  end subroutine get_psi_general_k_all_wfc

!---------------------------------------------
end module intw_allwfcs
!---------------------------------------------
