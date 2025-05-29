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
program interpolatephonons

  ! MBR June 2024

  ! Interpolate dynamical matrices using atom-pair-adapted WS vectors

  use kinds, only: dp
  use intw_useful_constants, only: Ha_to_eV, Ha_to_Ry, tpi
  use intw_utility, only: find_free_unit, &
                          generate_kmesh, cryst_to_cart, smeared_delta, &
                          generate_and_allocate_kpath
  use intw_matrix_vector, only: ainv
  use intw_input_parameters, only: read_input, read_cards, &
                                   mesh_dir, prefix, &
                                   nq1, nq2, nq3, nqirr, fc_mat, &
                                   exist_qpath, nqpath, nqspecial, qspecial, &
                                   nq1_dosph, nq2_dosph, nq3_dosph, &
                                   nomega, omega_ini, omega_fin, osmear_q, read_for_dynmat
  use intw_reading, only: read_parameters_data_file_xml, &
                          nat, bg, at, tpiba, nsym, tau
  use intw_ph, only: nqmesh, qmesh, QE_folder_nosym_q, QE_folder_sym_q, &
                     nosym_G_q, sym_G_q, symlink_q, q_irr, q_irr_cryst, &
                     read_ph_information_xml

  use intw_ph_interpolate, only: dyn_q, w2_q, u_q, dyn_diagonalize_1q, &
                                 dyn_q_to_dyn_r, dyn_interp_1q, &
                                 allocate_and_build_ws_irvec_qtau, &
                                 allocate_and_build_dyn_qmesh, &
                                 allocate_and_build_dyn_qmesh_from_fc
  use intw_symmetries, only: rtau, rtau_cryst, rtau_index, rot_atoms, find_size_of_irreducible_k_set, &
                             set_symmetry_relations, find_inverse_symmetry_matrices_indices

  implicit none

  character(256) :: phband_file_name
  logical :: read_status
  logical :: q_points_consistent, full_mesh_q, IBZ_q
  integer :: iq, iq1, iq2, iq3, iomega, imode
  integer :: qmesh_nqirr
  integer :: ph_unit
  integer , allocatable :: qspecial_indices(:)
  real(dp) :: omega_step, omega, rfacq
  real(dp) :: qpoint(3)
  real(dp), allocatable :: qpath(:,:), dqpath(:), dosph(:,:)
  real(dp), allocatable :: w2_qint(:), w_qint(:)
  complex(dp), allocatable :: dyn_qint(:,:), u_qint(:,:)


  20 format(A)
  !
  !
  !================================================================================
  ! Begining
  !================================================================================
  !
  write(*,'(A)') '====================================================='
  write(*,'(A)') '|         program phonons interpolate               |'
  write(*,'(A)') '|        ---------------------------------          |'
  write(*,'(A)') '====================================================='
  write(*,'(A)') '|    waiting for input file...                      |'
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
  if ( .not. exist_qpath ) then
    write(*,*) ' Q_PATH not found. Phonon bands/DOS cannot be interpolated. Stopping.'
    stop
  end if
  !
  !
  !================================================================================
  ! Read the parameters from the DFT calculation
  !================================================================================
  !
  write(*,20) '|       - Reading calculation parameters...         |'
  !
  call read_parameters_data_file_xml()
  !
  !
  !================================================================================
  ! Read phonon information
  !================================================================================
  !
  write(*,20) '|       - Reading q-points...                      |'
  !
  ! Read q-points and irreducible patterns
  call read_ph_information_xml()
  !
  allocate(q_irr_cryst(3,nqirr))
  do iq=1,nqirr
    q_irr_cryst(:,iq) = matmul(ainv(bg),q_irr(:,iq))
  enddo
  !
  !
  !================================================================================
  ! Build the phonon q-mesh
  !================================================================================
  !
  write(*,20) '|       - Building q-mesh...                        |'
  !
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  !
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !
  !
  !================================================================================
  ! Set symmetry arrays
  !================================================================================
  !
  write(*,20) '|       - Setting symmetry arrays...                |'
  !
  ! Set the rotation table for each atom and symmetry
  allocate(rtau_index(nat,nsym))
  allocate(rtau(3,nsym,nat))
  allocate(rtau_cryst(3,nsym,nat))
  !
  call rot_atoms(nat,nsym,tau)
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
  allocate(nosym_G_q(3,nqmesh))
  allocate(sym_G_q(3,nqmesh))
  allocate(symlink_q(nqmesh,2))
  !
  ! Find the size of the irreducible set of q-points (IBZ)
  call find_size_of_irreducible_k_set(nq1,nq2,nq3,qmesh,qmesh_nqirr)
  !
  call set_symmetry_relations(nq1, nq2, nq3, nqirr, q_irr_cryst, qmesh, q_points_consistent, &
                              QE_folder_nosym_q, QE_folder_sym_q, &
                              nosym_G_q, sym_G_q, symlink_q, full_mesh_q, IBZ_q)
  !
  !
  !================================================================================
  ! Check that the number of kpoints corresponds to either a full mesh or the IBZ
  !================================================================================
  !
  if (full_mesh_q .and. IBZ_q) then
    write(*,20) '|       - the qpoints present in the QE folders     |'
    write(*,20) '|         are consistent with a full 1BZ and a      |'
    write(*,20) '|         IBZ has also been found.                  |'
    write(*,20) '|           ---------------------------------       |'
  else if(IBZ_q) then
    write(*,20) '|       - the qpoints present in the QE folders     |'
    write(*,20) '|         are consistent with an IBZ.               |'
    write(*,20) '|           ---------------------------------       |'
  else
    write(*,20) '**********************************************************'
    write(*,20) '* The qpoints present in the QE folders are not consistent'
    write(*,20) '* with the parameters of the input file!                 '
    write(*,20) '**********************************************************'
    write(*,20) '* debug information:                                *'
    write(*,*) '*        nqpoints_QE = ', nqirr
    write(*,*) '*        nqmesh      = ', nqmesh
    write(*,*) '*        qmesh_nqirr = ', qmesh_nqirr
    stop
  end if
  !
  !
  !================================================================================
  !   Build qpoint path to plot bands.
  !   The nqpath number of points from the input might fluctuate.
  !   Use qspecial_indices option to print out the special q-points
  !   along the path (useful for plotting).
  !================================================================================
  !
  write(*,20) '|       - Building q-path...                        |'
  !
  call generate_and_allocate_kpath(at, bg, tpiba, nqpath, nqspecial, qspecial, qpath, dqpath, &
          qspecial_indices)
  !
  !
  !================================================================================
  !   Wigner-Seitz cells
  !================================================================================
  !
  write(*,20) '|       - Building WS cells...                      |'
  !
  call allocate_and_build_ws_irvec_qtau()
  !
  !
  !================================================================================
  ! Get dynamical matrices. Two options to get dyn_q:
  !  - Read dyn files.
  !  - Read the force constants and transfor to q space
  !================================================================================
  !
  write(*,20) '|       - Reading dynamical matrices...             |'
  !
  if (read_for_dynmat == 'fc' ) then ! read force constants
    call allocate_and_build_dyn_qmesh_from_fc(fc_mat)
  else if (read_for_dynmat == 'dynq' ) then ! read dyn files
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
  ! Transform dynamical matrices to real space (force constants)
  !================================================================================
  !
  write(*,20) '|       - Compute force constants...                |'
  !
  call dyn_q_to_dyn_r()
  !
  ! test decay of dyn_r elements with distance
  ! do ir=1,nrpts_q
  !   rcart = real(irvec_q(:,ir),dp)
  !   call cryst_to_cart (1, rcart, at, 1)
  !   rcart = rcart * alat  ! bohr units
  !   write(520,'(i5,f16.6,8e16.4)') ir,  sqrt ( sum(rcart*rcart) ), &
  !           abs(dyn_r(1,1,ir)), abs(dyn_r(1,2,ir)), abs(dyn_r(1,4,ir)), abs(dyn_r(1,5,ir))
  ! end do
  !
  !
  !================================================================================
  ! Interpolate on qpath for bands plot
  !================================================================================
  !
  write(*,20) '|       - Interpolate dynamical matrices...         |'
  !
  allocate( dyn_qint(3*nat,3*nat), u_qint(3*nat,3*nat), w2_qint(3*nat), w_qint(3*nat) )
  !
  phband_file_name = trim(mesh_dir)//trim(prefix)//".qbnd_int"
  ph_unit = find_free_unit()
  open(ph_unit,file=phband_file_name,status='unknown')
  write(ph_unit,'(A)') '# q-point   omega(imode=1)[meV]  omega(2)[meV]   omega(3)[meV] ...'
  do iq=1,nqpath
    qpoint = qpath(:,iq)
    call dyn_interp_1q(qpoint, dyn_qint)
    call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint) ! freqs are given in a.u
    w_qint = sign(sqrt(abs(w2_qint)),w2_qint) * Ha_to_eV*1000.0_dp
    write(ph_unit,'(20e14.6)') dqpath(iq), w_qint ! meV
    ! write(ph_unit,'(20e14.6)') dqpath(iq)/tpiba, w_qint*8.065610_dp ! Matdyn (cm^-1)
    ! write(ph_unit,'(20e14.6)') dqpath(iq)/tpi, w_qint/4.135665538536_dp ! Phonopy (tHz)
  end do
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
  write(*,*)' Phonon bands interpolation finished and written to file', phband_file_name
  !
  !
  !================================================================================
  ! Interpolate on fine q-grid for bands plot
  ! Parameters of DOS plot from namelist /DOS_ph/
  !================================================================================
  !
  allocate (dosph(3*nat,nomega))
  dosph = 0.0_dp
  !
  omega_step = (omega_fin-omega_ini)/real(nomega-1,dp)
  !
  ! Fine q-grid
  do iq1 = 1, nq1_dosph
    qpoint(1) = real(iq1-1,dp) / real(nq1_dosph,dp)
    do iq2 = 1, nq2_dosph
      qpoint(2) = real(iq2-1,dp) / real(nq2_dosph,dp)
      do iq3 = 1, nq3_dosph
        qpoint(3) = real(iq3-1,dp) / real(nq3_dosph,dp)
        !
        ! Interpolate frequency in qpoint
        call dyn_interp_1q(qpoint, dyn_qint)
        !
        call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint)
        w_qint=sign(sqrt(abs(w2_qint)),w2_qint)
        !
        ! Phonon frequency in a.u.: pass to Ry
        w_qint = w_qint*Ha_to_Ry
        !
        ! Smear omega(q) for DOS (gaussian)
        do imode = 1,3*nat
          do iomega=1, nomega ! Frequencies in Ry
            omega = omega_ini + omega_step*real(iomega-1,dp)
            rfacq = smeared_delta(omega-w_qint(imode), osmear_q)
            dosph(imode,iomega) = dosph(imode,iomega) + rfacq
          end do
        end do
        !
      end do
    end do
  end do
  !
  dosph = dosph / real(nq1_dosph*nq2_dosph*nq3_dosph,dp) ! Normalize for Nq points
  !
  ! Write DOS to file
  phband_file_name = trim(mesh_dir)//trim(prefix)//".qdos_int"
  ph_unit = find_free_unit()
  open(ph_unit,file=phband_file_name,status='unknown')
  !
  write(ph_unit,'(A)') '# omega[Ry]  phonon-DOS(total)  PDOS(imode=1)  PDOS(imode=2)  PDOS(imode=3) ...'
  !
  do iomega=1,nomega
    omega = omega_ini + omega_step*real(iomega-1,dp)
    write(ph_unit,'(40e14.6)') omega, sum(dosph(:,iomega)), dosph(:,iomega)
  end do
  !
  close(ph_unit)
  !
  write(*,'(A)') '|  DOS sum test:                                    |'
  write(*,*)'       DOS integral (trapeze) = ', omega_step*sum(dosph(:,:))
  write(*,*)'       Number of modes =', 3*nat
  write(*,*)' Phonon DOS interpolation finished and written to file', phband_file_name
  !
  !
  !================================================================================
  !       clean up and finish
  !================================================================================
  !
  write(*,'(A)') '====================================================='
  write(*,'(A)') '|               end program phonons                 |'
  write(*,'(A)') '|        ---------------------------------          |'
  write(*,'(A)') '====================================================='

end program interpolatephonons
