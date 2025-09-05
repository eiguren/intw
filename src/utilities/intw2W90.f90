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
program intw2W90

  ! The purpose of this utility is to perform the same tasks as the
  ! program "pw2wannier90" which is part of the QE distribution, but
  ! utilizing a minimum set of (QE generated) Bloch functions, using
  ! symmetry.

  use kinds, only: dp
  use intw_version, only: print_intw_version
  use intw_intw2wannier, only: nnkp_exclude_bands, read_nnkp_file, output_nnkp_file, nnkp_n_proj, &
                               intw2w90_check_mesh, generate_mmn_using_allwfc, &
                               generate_amn_using_allwfc, deallocate_nnkp
  use intw_symmetries, only: full_mesh, IBZ, QE_folder_nosym, QE_folder_sym, &
                             nosym_G, sym_G, symlink, find_size_of_irreducible_k_set, &
                             find_the_irreducible_k_set, allocate_symmetry_related_k, &
                             find_inverse_symmetry_matrices_indices, &
                             allocate_and_build_spin_symmetry_matrices, &
                             set_symmetry_relations, deallocate_symmetry_related_k, &
                             deallocate_spin_symmetry_matrices
  use intw_fft, only: allocate_fft, deallocate_fft, generate_nl
  use intw_input_parameters, only: outdir, prefix, intw2W_fullzone, intw2W_method, &
                                   nk1, nk2, nk3, compute_mmn, compute_amn, &
                                   read_input
  use intw_reading, only: kpoints_QE, nspin, lspin, nkpoints_QE, nsym, &
                          read_parameters_data_file, get_gvec, &
                          read_kpoints_data_file, deallocate_reading_variables, &
                          num_wann_intw, num_exclude_bands_intw, &
                          set_num_bands
  USE intw_allwfcs, only: allocate_and_get_all_irreducible_wfc
  use intw_utility, only: get_timing, print_date_time, generate_kmesh

!================================================================================
! Declare the variables
!================================================================================
  implicit none

  integer :: nk_irr , nkmesh

  logical :: read_status, have_nnkp

  logical :: k_points_consistent

  character(256) :: nnkp_file

  real(dp), allocatable :: kmesh(:,:)

  real(dp) :: time1, time2


  20 format(A)
  30 format(A,F8.2,6X,A)
  !
  !
  !================================================================================
  ! Beginning
  !================================================================================
  !
  call get_timing(time1)
  !
  write(*,20) '====================================================='
  write(*,20) '|                  program intw2W90                 |'
  write(*,20) '|         ---------------------------------         |'
  call print_intw_version()
  call print_date_time("Start of execution")
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Read the necessary information from standard input
  !================================================================================
  !
  call read_input(read_status)
  !
  if (read_status ) stop
  !
  if (trim(intw2W_method) == 'CONVOLUTION') then
    write(*,20) '| - intw2W_method   = CONVOLUTION                   |'
  else if (trim(intw2W_method) == 'FFT') then
    write(*,20) '| - intw2W_method   = FFT                           |'
  else
    write(*,20) '***********************************************'
    write(*,20) '* UNKNOWN COMPUTATION METHOD:'
    write(*,20) '* Only "CONVOLUTION" and "FFT" available'
    write(*,20) '***********************************************'
    stop
  end if
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Check that $prefix.nnkp is present
  !================================================================================
  !
  nnkp_file = trim(outdir)//trim(prefix)//".nnkp"
  !
  inquire(file=nnkp_file,exist=have_nnkp)
  !
  if(.not. have_nnkp) then
    write(*,20) '**********************************************************'
    write(*,20) '* Could not find the file '//trim(nnkp_file)
    write(*,20) '* Did you run W90 -pp $seed to get the parameter file?   '
    write(*,20) '**********************************************************'
    stop
  end if
  write(*,20) '| - .nnkp file found                                |'
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Read the parameters from the SCF calculation
  !================================================================================
  !
  call read_parameters_data_file()
  !
  !
  !================================================================================
  ! Set the number of wave functions
  !================================================================================
  !
  call set_num_bands()
  !
  !
  !================================================================================
  ! Print spin information
  !================================================================================
  !
  if (lspin) then
    write(*,20) '| - Spin-polarized calculation nspin = 2            |'
  else
    write(*,20) '| - Paramagnetic calculation nspin = 1              |'
  endif
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Set symmetry arrays
  !================================================================================
  !
  write(*,20) '| - Setting symmetry arrays...                      |'
  !
  ! Compute the indices of the inverse rotation matrices
  call find_inverse_symmetry_matrices_indices()
  !
  ! Set up spin_symmetry_matrices, needed to rotate wave functions and indueced potential for non-colinear calculations
  call allocate_and_build_spin_symmetry_matrices(nsym)
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Read the .nnkp file
  !================================================================================
  !
  write(*,20) '| - Reading nnkp file...                            |'
  !
  call read_nnkp_file(nnkp_file)
  !
  if (nnkp_exclude_bands /= num_exclude_bands_intw) then
    write(*,20) '**********************************************************'
    write(*,20) '* Mismatch in number of excluded bands                    '
    write(*,20) '**********************************************************'
    stop
  else if (nnkp_n_proj /= num_wann_intw) then
    write(*,20) '**********************************************************'
    write(*,20) '* Mismatch in number of projections/wannier functions     '
    write(*,20) '**********************************************************'
    stop
  end if
  !
  ! just as a test; can be removed later
  call output_nnkp_file()
  !
  !
  !================================================================================
  ! Set up the gvec array and all FFT variables
  !================================================================================
  !
  write(*,20) '| - Reading G vectors...                            |'
  !
  call get_gvec()
  !
  ! Allocate useful variables
  call allocate_fft()
  !
  ! Generate some important indices for FFT
  call generate_nl()
  !
  !
  !================================================================================
  ! Read the kpoints from the calculation
  !================================================================================
  !
  write(*,20) '| - Reading k-points...                             |'
  !
  allocate(kpoints_QE(3,nkpoints_QE))
  !
  call read_kpoints_data_file(kpoints_QE)
  !
  !
  !================================================================================
  ! Build the wave function's k-mesh
  !================================================================================
  !
  write(*,20) '| - Building k-mesh...                              |'
  !
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh, nk1, nk2, nk3)
  !
  ! Find the size of the irreducible set of k-points (IBZ)
  call find_size_of_irreducible_k_set(nk1, nk2, nk3, nk_irr)
  !
  !
  !================================================================================
  ! Check that kmesh and nnkp_kpoints are consistent
  !================================================================================
  !
  call intw2W90_check_mesh(nkmesh, kmesh)
  !
  write(*,20) '| - The mesh in the Wannier90 input.win file        |'
  write(*,20) '|   and the intw mesh are equal.                    |'
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Set symmetry relations between irreducible k-points and full k-mesh
  !================================================================================
  !
  ! Allocate arrays
  call allocate_symmetry_related_k(nk1, nk2, nk3)
  !
  ! Fill the symmetry arrays
  call set_symmetry_relations(nk1, nk2, nk3, nkpoints_QE, kpoints_QE, &
                              QE_folder_nosym, nosym_G, QE_folder_sym, sym_G, &
                              symlink, full_mesh, IBZ )
  !
  !
  !================================================================================
  ! Check that the number of kpoints corresponds to either a full mesh or the IBZ
  !================================================================================
  !
  if (full_mesh .and. IBZ) then
    write(*,20) '| - The kpoints present in the QE folders           |'
    write(*,20) '|   are consistent with a full 1BZ and a            |'
    write(*,20) '|   IBZ has also been found.                        |'
    write(*,20) '|         ---------------------------------         |'
  else if(IBZ) then
    write(*,20) '| - The kpoints present in the QE folders           |'
    write(*,20) '|   are consistent with an IBZ.                     |'
    write(*,20) '|         ---------------------------------         |'
  else
    write(*,*) '**********************************************************'
    write(*,*) '* The kpoints present in the QE folders are not consistent'
    write(*,*) '* with the parameters of the input file!                 '
    write(*,*) '**********************************************************'
    write(*,*) '* debug information:                                *'
    write(*,*) '*        nkpoints_QE = ',nkpoints_QE
    write(*,*) '*        nkmesh      = ',nkmesh
    write(*,*) '*        nk_irr      = ',nk_irr
    stop
  end if
  !
  !
  !================================================================================
  ! Check that the requested calculation is possible
  !================================================================================
  !
  if (intw2W_fullzone) then
    write(*,20) '| - intw2W_fullzone = .true.                        |'
    if (full_mesh) then
      write(*,20) '|   all k-points are explicitely calculated         |'
      write(*,20) '|   no symmetry is assumed.                         |'
      write(*,20) '|   (This is mostly for testing)                    |'
      write(*,20) '|         ---------------------------------         |'
    else
      write(*,20) '**********************************************************'
      write(*,20) '* A full mesh is not present in the QE folders!          '
      write(*,20) '* The requested calculation is impossible.               '
      write(*,20) '*                   program stops.                       '
      write(*,20) '**********************************************************'
      stop
    end if
  else
    write(*,20) '| - intw2W_fullzone = .false.                       |'
    write(*,20) '|   Symmetries will be utilized.                    |'
    write(*,20) '|         ---------------------------------         |'
  end if
  !
  !
  !================================================================================
  ! Read all wave functions
  !================================================================================
  !
  write(*,20) '| - Reading wave functions...                       |'
  !
  call allocate_and_get_all_irreducible_wfc()
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Compute the mmn file
  !================================================================================
  !
  if (compute_mmn) then
    write(*,20) '| - Computing the file prefix.mmn and               |'
    write(*,20) '|   writing the file prefix.eig...                  |'
    write(*,20) '|  (this is labor intensive and may take some time) |'
    write(*,20) '|         ---------------------------------         |'

    call generate_mmn_using_allwfc(intw2W_fullzone,intw2W_method)
  end if
  !
  !
  !================================================================================
  ! Compute the amn file
  !================================================================================
  !
  if (compute_amn) then
    write(*,20) '| - Computing the file prefix.amn...                |'
    write(*,20) '|  (this is labor intensive and may take some time) |'
    write(*,20) '|         ---------------------------------         |'

    call generate_amn_using_allwfc(intw2W_fullzone,intw2W_method)
  end if
  !
  !
  !================================================================================
  ! Clean up
  !================================================================================
  !
  call deallocate_symmetry_related_k()
  call deallocate_nnkp()
  call deallocate_fft()
  call deallocate_reading_variables()
  call deallocate_spin_symmetry_matrices()
  deallocate(kpoints_QE)
  deallocate(kmesh)
  !
  !
  !================================================================================
  ! Finish
  !================================================================================
  !
  call get_timing(time2)
  !
  write(*,20) '|                      ALL DONE                     |'
  write(*,30) '|     Total time: ',time2-time1,' seconds            |'
  call print_date_time('End of execution  ')
  write(*,20) '====================================================='

end program intw2W90
