!
! Copyright (C) 2024 INTW group
!
! This file is part of INTW.
!
! INTW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! INTW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
module intw_allwfcs

  !---------------------------------------------------------------------------------!
  ! This module contains variables and subroutines to obtain wave functions.        !
  ! Wave functions for irreducible k points are read and stored, and wave functions !
  ! for general k points are obtained using symmetry.                               !
  !---------------------------------------------------------------------------------!

  use kinds, only: dp

  implicit none

  ! variables
  public :: eig_all_irr, wfc_k_all_irr, list_iG_all_irr, ngk_all_irr

  ! subroutines
  public :: allocate_and_get_all_irreducible_wfc, get_psi_general_k_all_wfc

  private

  real(dp), allocatable :: eig_all_irr(:,:)
  ! Eigenvalues for each band and irreducible k point
  complex(dp), allocatable :: wfc_k_all_irr(:,:,:,:)
  ! Wave function components (also spin component) for each band and irreducible k point
  integer, allocatable :: list_iG_all_irr(:,:)
  ! G vector indices of the wave function components for each irreducible k point
  integer, allocatable :: ngk_all_irr(:)
  ! Number of components of the wave functions for each irreducible k point

contains

  subroutine allocate_and_get_all_irreducible_wfc()

    use intw_reading, only: nG_max, nkpoints_QE, get_K_folder_data, nspin, num_bands_intw
    use intw_useful_constants, only: cmplx_0, ZERO

    implicit none

    integer :: list_iG(nG_max)
    real(dp) :: QE_eig(num_bands_intw)
    complex(dp) :: wfc_g(nG_max,num_bands_intw,nspin)
    integer :: i_folder, ispin, iG, ibnd


    !
    ! Allocate eigenvalue and wfc related arrays
    !
    if (allocated(eig_all_irr)) deallocate (eig_all_irr)
    allocate(eig_all_irr(nkpoints_QE,num_bands_intw))
    eig_all_irr = ZERO
    !
    if (allocated(wfc_k_all_irr)) deallocate (wfc_k_all_irr)
    allocate(wfc_k_all_irr(nkpoints_QE,nG_max,num_bands_intw,nspin))
    wfc_k_all_irr = cmplx_0
    !
    if (allocated(list_iG_all_irr)) deallocate(list_iG_all_irr)
    allocate(list_iG_all_irr(nkpoints_QE,nG_max))
    list_iG_all_irr = 0
    !
    if (allocated(ngk_all_irr)) deallocate(ngk_all_irr)
    allocate(ngk_all_irr(nkpoints_QE))
    ngk_all_irr = 0

    !
    ! Read all irreducible k points
    !
    do i_folder = 1, nkpoints_QE
      !
      call get_K_folder_data(i_folder, list_iG, wfc_g, QE_eig, ngk_all_irr(i_folder))
      !
      eig_all_irr(i_folder,:) = QE_eig(:)
      !
      list_iG_all_irr(i_folder, 1:ngk_all_irr(i_folder)) = list_iG(1:ngk_all_irr(i_folder))
      !
      do iG = 1, ngk_all_irr(i_folder)
        do ispin = 1, nspin
          do ibnd = 1, num_bands_intw
            !
            wfc_k_all_irr(i_folder,iG,ibnd,ispin) = wfc_g(iG,ibnd,ispin)
            !
          enddo ! ibnd
        enddo ! ispin
      enddo ! iG
      !
    enddo ! i_folder

  end subroutine allocate_and_get_all_irreducible_wfc


  subroutine get_psi_general_k_all_wfc(kpoint, ngk, list_iG_k, wfc_k, eig_k)

    use intw_reading, only: s, ftau, nG_max, nspin, kpoints_QE, num_bands_intw
    use intw_input_parameters, only: nk1, nk2, nk3
    use intw_symmetries, only: full_mesh, symlink, QE_folder_sym, QE_folder_nosym, &
                               apply_TR_to_wfc, rotate_wfc_test
    use intw_utility, only: find_k_1BZ_and_G, triple_to_joint_index_g
    use intw_fft, only: wfc_by_expigr
    use intw_useful_constants, only : eps_5, ZERO

    implicit none

    !I/O variables
    real(dp), intent(in) :: kpoint(3) ! General k point in crystal coordinates for which wave functions are desired
    integer, intent(out) :: ngk ! Number of components of the wave functions
    integer, intent(out) :: list_iG_k(nG_max) ! G vector indices of the wave function components
    complex(dp), intent(out) :: wfc_k(nG_max,num_bands_intw,nspin) ! Wave function components (also spin component) for all bands
    real(dp), intent(out), optional :: eig_k(num_bands_intw) ! Eigenvalues for all band

    !local variables
    integer :: ikpt, i_folder
    integer :: i_sym, TR
    integer :: i_1bz, j_1bz, k_1bz
    integer :: sym(3,3), G_sym(3), G_plus(3)
    real(dp) :: ftau_sym(3)
    real(dp) :: kpoint_1BZ(3), k_irr(3), ktest(3)
    integer :: list_iG_irr(nG_max)
    complex(dp) :: wfc_k_irr(nG_max,num_bands_intw,nspin)


    call find_k_1BZ_and_G(kpoint, nk1, nk2, nk3, i_1bz, j_1bz, k_1bz, kpoint_1bz, G_plus)
    !
    call triple_to_joint_index_g(nk1, nk2, nk3, ikpt, i_1bz, j_1bz, k_1bz)
    !
    if (full_mesh) then
      !
      ! Use the full BZ, no symmetry!
      !
      i_folder = QE_folder_nosym(ikpt)
      k_irr = kpoints_QE(:,i_folder)
      !
      ngk = ngk_all_irr(ikpt)
      list_iG_k(:) = list_iG_all_irr(ikpt,:)
      wfc_k(:,:,:) = wfc_k_all_irr(ikpt,:,:,:)
      if (present(eig_k)) eig_k(:) = eig_all_irr(ikpt,:)
      !
      G_sym = kpoint - k_irr
      !
      call wfc_by_expigr(num_bands_intw, nspin, G_sym, list_iG_k, wfc_k)
      !
    else
      !
      ! Use the IBZ and symmetry!
      !
      ! This is the irreducible point
      ! connected to aimed kpoint
      i_folder = QE_folder_sym(ikpt)
      k_irr = kpoints_QE(:,i_folder)
      !
      ! The symmetry which takes kpoints_QE(:,i_folder) into aimed kpoint.
      ! sym * kpoints_QE = kpoint
      i_sym = symlink(ikpt,1)
      TR = symlink(ikpt,2)
      ftau_sym = ftau(:,(i_sym))
      sym = s(:,:,(i_sym))
      !
      ! Load the corresponding irreducible wfcs in kpoints_QE
      ngk = ngk_all_irr(i_folder)
      list_iG_irr(:) = list_iG_all_irr(i_folder,:)
      wfc_k_irr(:,:,:) = wfc_k_all_irr(i_folder,:,:,:)
      if (present(eig_k)) eig_k(:) = eig_all_irr(i_folder,:)
      !
      ktest = matmul(sym ,k_irr)
      !
      if (TR==1) then
        ! If TR needed
        ! G_sym = aimed_kpoint -TR[S*QE_kpoint] = aimed_kpoint + S*QE_kpoint, such that
        ! aimed_kpoint = -S*QE_kpoint + G_ sym
        G_sym = nint(kpoint + ktest)
      else
        ! G_sym = aimed_kpoint - S*QE_kpoint, such that
        ! aimed_kpoint = S*QE_kpoint + G_ sym
        G_sym = nint(kpoint - ktest)
      end if
      !
      call rotate_wfc_test(wfc_k_irr, list_iG_irr, wfc_k, list_iG_k, i_sym, sym, ftau_sym, (/0, 0, 0/))
      !
      ! If time-reversal is present, the wavefunction currently stored
      ! in wfc_k is actually for (-k). Complex conjugation must now
      ! be applied to recover wfc_k.
      !
      if (TR==1) then
        ktest = -ktest
        !
        call apply_TR_to_wfc(wfc_k, list_iG_k)
        !
      endif
      !
      if ( sum( abs( ktest + dble(G_sym) - kpoint ) ) > eps_5 ) then
          write(*,*)            "ERROR in get_psi_general_k_all_wfc:"
          write(*,"(a,3f12.6)") "Aimed kpoint is     :", kpoint
          write(*,"(a,3f12.6)") "Symmetry induced is :", ktest + G_sym
      end if
      !
      call wfc_by_expigr(num_bands_intw, nspin, G_sym, list_iG_k, wfc_k)
      !
    endif

  end subroutine get_psi_general_k_all_wfc

end module intw_allwfcs
