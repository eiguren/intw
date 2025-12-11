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
program interpolate_phonons

  !! display: none
  !!
  !! Interpolate phonon dispersion and DOS using Fourier interpolation.
  !!
  !! ### Details
  !!
  !! Uses atom-pair-adapted WS vectors to interpolate dynamical matrices.
  !!
  !! MBR June 2024
  !!
  !! #### Input parameters
  !!
  !! ```{.txt}
  !! &input
  !!     outdir                = 'directory'
  !!     prefix                = 'prefix'
  !!     TR_symmetry           = T or F
  !!     use_exclude_bands     = 'none', 'wannier' or 'custom'
  !!     include_bands_initial = integer
  !!     include_bands_final   = integer
  !! /
  !! &ph
  !!     qlist           = 'file'
  !!     read_for_dynmat = 'fc' or 'dynq'
  !!     fc_mat          = 'file'
  !!     nq1             = integer
  !!     nq2             = integer
  !!     nq3             = integer
  !!     nqirr           = integer
  !!     apply_asr       = T or F
  !! /
  !! &DOS_ph
  !!     nq1_dosph = integer
  !!     nq2_dosph = integer
  !!     nq3_dosph = integer
  !!     nomega    = integer
  !!     omega_ini = real
  !!     omega_fin = real
  !!     osmear_q  = real
  !! /
  !! Q_PATH
  !!     nqpath nqspecial
  !!     label(1) qspecial_x(1) qspecial_y(1) qspecial_z(1)
  !!     label(2) qspecial_x(2) qspecial_y(2) qspecial_z(2)
  !!     ...
  !!     label(nqspecial) qspecial_x(nqspecial) qspecial_y(nqspecial) qspecial_z(nqspecial)
  !! ```
  !!
  !! See [[intw_input_parameters]] module for the description of each parameter.
  !!

#ifdef _OPENMP
  use omp_lib, only: omp_set_max_active_levels
#endif

  use kinds, only: dp

  use intw_version, only: print_intw_version

  use intw_useful_constants, only: Ha_to_eV, tpi

  use intw_utility, only: get_timing, print_threads, print_date_time, find_free_unit, &
                          generate_kmesh, cryst_to_cart, generate_and_allocate_kpath


  use intw_input_parameters, only: read_input, read_cards, &
                                   outdir, prefix, &
                                   nq1, nq2, nq3, nqirr, fc_mat, &
                                   exist_qpath, nqpath, nqspecial, qspecial, &
                                   nq1_dosph, nq2_dosph, nq3_dosph, &
                                   nomega, omega_ini, omega_fin, osmear_q, read_for_dynmat

  use intw_reading, only: read_parameters_data_file, &
                          nat, bg, at, tpiba, nsym, tau

  use intw_ph, only: nqmesh, qmesh, QE_folder_nosym_q, QE_folder_sym_q, &
                     symlink_q, q_irr_cryst, &
                     read_ph_information

  use intw_ph_interpolate, only: dyn_q, w2_q, u_q, dyn_diagonalize_1q, &
                                 dyn_q_to_dyn_r, dyn_interp_1q, &
                                 allocate_and_build_ws_irvec_qtau, &
                                 allocate_and_build_dyn_qmesh, &
                                 allocate_and_build_dyn_qmesh_from_fc, &
                                 interpolated_phonon_DOS

  use intw_symmetries, only: rtau, rtau_cryst, rtau_index, rot_atoms, find_size_of_irreducible_k_set, &
                             set_symmetry_relations, find_inverse_symmetry_matrices_indices

  implicit none

  character(256) :: phband_file_name
  logical :: read_status
  logical :: full_mesh_q, IBZ_q
  integer :: iq, iomega, imode
  integer :: qmesh_nqirr
  integer :: ph_unit
  integer, allocatable :: qspecial_indices(:)
  real(dp) :: omega_step, omega
  real(dp) :: qpoint(3)
  real(dp), allocatable :: qpath(:,:), dqpath(:), dosph(:,:)
  real(dp), allocatable :: w2_qint(:), w_qint(:)
  complex(dp), allocatable :: dyn_qint(:,:), u_qint(:,:)

  ! timing
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
  write(*,20) '|            program interpolate_phonons            |'
  write(*,20) '|         ---------------------------------         |'
  call print_intw_version()
#ifdef _OPENMP
  call omp_set_max_active_levels(1) ! This utility usea a single active parallel level
#endif
  call print_threads()
  call print_date_time("Start of execution")
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Read the necessary information from standard input file
  !================================================================================
  !
  call read_input(read_status)
  !
  if (read_status) stop
  !
  ! generate q-list for phonon bands plot with special points in Q_PATH
  call read_cards()
  if (.not. exist_qpath) then
    write(*,*) 'Q_PATH not found. Phonon bands/DOS cannot be interpolated. Stopping.'
    stop
  end if
  !
  !
  !================================================================================
  ! Read the parameters from the DFT calculation
  !================================================================================
  !
  write(*,20) '| - Reading calculation parameters...               |'
  !
  call read_parameters_data_file()
  !
  !
  !================================================================================
  ! Read phonon information
  !================================================================================
  !
  write(*,20) '| - Reading phonon info...                          |'
  !
  ! Read q-points and irreducible patterns
  call read_ph_information()
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Build the phonon q-mesh
  !================================================================================
  !
  write(*,20) '| - Building q-mesh...                              |'
  !
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  !
  call generate_kmesh(qmesh, nq1, nq2, nq3)
  !
  !
  !================================================================================
  ! Set symmetry arrays
  !================================================================================
  !
  write(*,20) '| - Setting symmetry arrays...                      |'
  !
  ! Set the rotation table for each atom and symmetry
  allocate(rtau_index(nat,nsym))
  allocate(rtau(3,nsym,nat))
  allocate(rtau_cryst(3,nsym,nat))
  !
  call rot_atoms(nat, nsym, tau)
  !
  ! Compute the indices of the inverse rotation matrices
  call find_inverse_symmetry_matrices_indices()
  !
  !
  !================================================================================
  ! Symmetry relations between irreducible q-points and full q-mesh
  !================================================================================
  !
  allocate(QE_folder_nosym_q(nqmesh))
  allocate(QE_folder_sym_q(nqmesh))
  allocate(symlink_q(nqmesh,2))
  !
  ! Find the size of the irreducible set of q-points (IBZ)
  call find_size_of_irreducible_k_set(nq1, nq2, nq3, qmesh_nqirr)
  !
  call set_symmetry_relations(nq1, nq2, nq3, nqirr, q_irr_cryst, &
                              QE_folder_nosym_q, QE_folder_sym_q, symlink_q, &
                              full_mesh_q, IBZ_q)
  !
  !
  !================================================================================
  ! Check that the number of kpoints corresponds to either a full mesh or the IBZ
  !================================================================================
  !
  if (full_mesh_q .and. IBZ_q) then
    write(*,20) '| - The qpoints present in the QE folders           |'
    write(*,20) '|   are consistent with a full 1BZ and a            |'
    write(*,20) '|   IBZ has also been found.                        |'
  else if(IBZ_q) then
    write(*,20) '| - The qpoints present in the QE folders           |'
    write(*,20) '|   are consistent with an IBZ.                     |'
  else
    write(*,*) '**********************************************************'
    write(*,*) '* The qpoints present in the QE folders are not consistent'
    write(*,*) '* with the parameters of the input file!                 '
    write(*,*) '**********************************************************'
    write(*,*) '* debug information:                                *'
    write(*,*) '*        nqpoints_QE = ', nqirr
    write(*,*) '*        nqmesh      = ', nqmesh
    write(*,*) '*        qmesh_nqirr = ', qmesh_nqirr
    stop
  end if
  !
  write(*,20) '|         ---------------------------------         |'
    !
  !
  !================================================================================
  ! Get dynamical matrices. Two options to get dyn_q:
  !  - Read dyn files.
  !  - Read the force constants and transfor to q space
  !================================================================================
  !
  write(*,20) '| - Reading dynamical matrices...                   |'
  !
  if (read_for_dynmat == 'fc') then ! read force constants
    call allocate_and_build_dyn_qmesh_from_fc(fc_mat)
  else if (read_for_dynmat == 'dynq') then ! read dyn files
    call allocate_and_build_dyn_qmesh()
  end if
  !
  ! diagonalize
  do iq=1,nqmesh
    qpoint = qmesh(:,iq)
    call dyn_diagonalize_1q(3*nat, dyn_q(:,:,iq), u_q(:,:,iq), w2_q(:,iq))
  end do
  !
  !
  !================================================================================
  ! Wigner-Seitz cells
  !================================================================================
  !
  write(*,20) '| - Building WS cells...                            |'
  !
  call allocate_and_build_ws_irvec_qtau()
  !
  !
  !================================================================================
  ! Transform dynamical matrices to real space (force constants)
  !================================================================================
  !
  write(*,20) '| - Computing force constants...                    |'
  !
  call dyn_q_to_dyn_r()
  !
  ! test decay of dyn_r elements with distance
  ! do ir=1,nrpts_q
  !   rcart = real(irvec_q(:,ir),dp)
  !   call cryst_to_cart(1, rcart, at, 1)
  !   rcart = rcart * alat ! bohr units
  !   write(1000,'(i5,f16.6,8e16.4)') ir, sqrt ( sum(rcart*rcart) ), &
  !           abs(dyn_r(1,1,ir)), abs(dyn_r(1,2,ir)), abs(dyn_r(1,4,ir)), abs(dyn_r(1,5,ir))
  ! end do
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Build qpoint path to plot bands.
  ! The nqpath number of points from the input might fluctuate.
  ! Use qspecial_indices option to print out the special q-points
  ! along the path (useful for plotting).
  !================================================================================
  !
  write(*,20) '| - Building q-path...                              |'
  !
  call generate_and_allocate_kpath(at, bg, tpiba, nqpath, nqspecial, qspecial, &
                                   qpath, dqpath, qspecial_indices)
  !
  !
  !================================================================================
  ! Interpolate on qpath for bands plot
  !================================================================================
  !
  write(*,20) '| - Computing phonon dispersion...                  |'
  !
  allocate(dyn_qint(3*nat,3*nat), u_qint(3*nat,3*nat), w2_qint(3*nat), w_qint(3*nat))
  !
  phband_file_name = trim(outdir)//trim(prefix)//".qbnd_int"
  ph_unit = find_free_unit()
  open(ph_unit, file=phband_file_name, status='unknown')
  write(ph_unit,'(A)') '# q-point   omega(imode=1)[meV]  omega(2)[meV]   omega(3)[meV] ...'
  do iq=1,nqpath
    qpoint = qpath(:,iq)
    call dyn_interp_1q(qpoint, dyn_qint)
    call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint) ! freqs are given in a.u
    w_qint = sign(sqrt(abs(w2_qint)), w2_qint) * Ha_to_eV*1000.0_dp
    write(ph_unit,'(100e14.6)') dqpath(iq), (w_qint(imode), imode=1,3*nat) ! meV
    ! write(ph_unit,'(100e14.6)') dqpath(iq)/tpiba, (w_qint(imode)*8.065610_dp, imode=1,3*nat) ! Matdyn (cm^-1)
    ! write(ph_unit,'(100e14.6)') dqpath(iq)/tpi, (w_qint(imode)/4.135665538536_dp, imode=1,3*nat) ! Phonopy (tHz)
  end do
  !
  deallocate(dyn_qint, u_qint, w2_qint, w_qint)
  !
  ! Print special q-points information in the phonon bands file
  write(ph_unit,*) '#'
  write(ph_unit,*) '#Special q-points in the .qbnd_int file are:'
  do iq=1,nqspecial
    write(ph_unit,'(a,3f10.4,a,i4,e14.6)') '#', qspecial(:,iq), ' --> ', qspecial_indices(iq), dqpath(qspecial_indices(iq))
  end do
  write(ph_unit,*) '#'
  !
  close(ph_unit)
  !
  write(*,20) '|   Phonon dispersion computed and written to:      |'
  write(*,20) "|   "//phband_file_name(1:max(47,len(trim(phband_file_name))))//" |"
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Interpolate on fine q-grid for bands plot
  ! Parameters of DOS plot from namelist /DOS_ph/
  !================================================================================
  !
  write(*,20) '| - Computing phonon DOS...                         |'
  !
  allocate(dosph(3*nat,nomega))
  !
  call interpolated_phonon_DOS(nq1_dosph, nq2_dosph, nq3_dosph, omega_ini, omega_fin, osmear_q, nomega, dosph)
  !
  ! Write DOS to file
  phband_file_name = trim(outdir)//trim(prefix)//".qdos_int"
  ph_unit = find_free_unit()
  open(ph_unit, file=phband_file_name, status='unknown')
  !
  write(ph_unit,'(A)') '# omega[Ry]  phonon-DOS(total)  PDOS(imode=1)  PDOS(imode=2)  PDOS(imode=3) ...'
  !
  omega_step = (omega_fin-omega_ini)/real(nomega-1,dp)
  do iomega=1,nomega
    omega = omega_ini + omega_step*real(iomega-1,dp)
    write(ph_unit,'(100e14.6)') omega, sum(dosph(:,iomega)), (dosph(imode,iomega), imode=1,3*nat)
  end do
  !
  close(ph_unit)
  !
  write(*,20) '| - DOS sum test:                                   |'
  write(*,'(A19,I10,23X,A1)') '|   Number of modes = ', 3*nat, '|'
  write(*,'(A19,F10.6,23X,A1)') '|   DOS integral = ', omega_step*sum(dosph(:,:)), '|'
  !
  write(*,20) '|   Phonon DOS computed and written to:             |'
  write(*,20) '|   '//phband_file_name(1:max(47,len(trim(phband_file_name))))//' |'
  !
  write(*,20) '====================================================='
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

end program interpolate_phonons
