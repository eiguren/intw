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
program ep_melements

  use kinds, only: dp
  use intw_version, only: print_intw_version
  use intw_input_parameters, only: nk1, nk2, nk3, &
                                   nq1, nq2, nq3, nqirr, outdir, &
                                   ep_mat_file, &
                                   read_input, &
                                   ep_bands, ep_bands_initial, ep_bands_final
  use intw_reading, only: nkpoints_QE, kpoints_QE, nspin, lspin, nsym, &
                          s, nGk_max, nat, tau, &
                          nr1, nr2, nr3, num_bands_intw, &
                          read_parameters_data_file, &
                          get_gvec, &
                          read_kpoints_data_file, &
                          set_num_bands
  use intw_pseudo, only: read_all_pseudo
  use intw_pseudo_local, only: calculate_local_part_dv, dvqpsi_local
  use intw_pseudo_non_local, only: init_KB_PP, &
                                   multiply_psi_by_dvKB
  use intw_utility, only: get_timing, print_threads, print_date_time, &
                          find_free_unit, &
                          generate_kmesh, &
                          conmesurate_and_coarser
  use intw_matrix_vector, only: ainv
  use intw_useful_constants, only: cmplx_0, cmplx_1
  use intw_symmetries, only: full_mesh, IBZ, QE_folder_nosym, QE_folder_sym, symlink, &
                             symtable, rtau_index, rtau, rtau_cryst, rot_atoms, &
                             find_size_of_irreducible_k_set, &
                             allocate_symmetry_related_k, &
                             find_inverse_symmetry_matrices_indices, &
                             allocate_and_build_spin_symmetry_matrices, &
                             set_symmetry_relations, multable
  use intw_fft, only: generate_nl, &
                      allocate_fft
  use intw_ph, only: nqmesh, qmesh, QE_folder_nosym_q, QE_folder_sym_q, &
                     symlink_q, q_irr_cryst, &
                     read_ph_information, &
                     read_allq_dvr, &
                     get_dv
  use intw_allwfcs, only: allocate_and_get_all_irreducible_wfc, &
                          get_psi_general_k_all_wfc

  !================================================================================
  ! Declare the variables
  !================================================================================
  implicit none

  !k point related variables
  integer                  :: ik, kmesh_nkirr, nkmesh
  real(dp), allocatable    :: kmesh(:,:)
  real(dp)                 :: kpoint(3)

  !q point related variables
  real(dp)                 :: qpoint(3)
  integer                  :: iq, qmesh_nqirr
  logical                  :: full_mesh_q, IBZ_q
  character(len=4)         :: iq_loc

  !wave function realted variables
  integer                  :: nGk, nGkq
  integer, allocatable     :: list_igk(:)
  integer, allocatable     :: list_igkq(:)
  complex(dp), allocatable :: wfc_k(:,:,:)
  complex(dp), allocatable :: wfc_kq(:,:,:)

  !phonon related variables
  complex(dp), allocatable :: dvq_local(:,:,:,:)

  !ep related variables
  integer                  :: num_bands_ep
  complex(dp), allocatable :: dvpsi(:,:,:,:,:)
  complex(dp), allocatable :: ep_mat_el(:,:,:,:,:,:)
  integer                  :: ep_unit, record_length, ierr

  !local/aux variables
  real(dp)                 :: time1, time2
  integer                  :: ibnd, jbnd, ispin, jspin, imode
  logical                  :: read_status

  complex(dp), external :: zdotc


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
  write(*,20) '|                program ep_melements               |'
  write(*,20) '|         ---------------------------------         |'
  call print_intw_version()
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
  !
  !================================================================================
  ! Check if the k-mesh and q-mesh are conmesurate
  !================================================================================
  !
  if (.not.conmesurate_and_coarser(nk1,nk2,nk3,nq1,nq2,nq3)) then
    write(*,20) '**********************************************************'
    write(*,20) '*ERROR                                                  '
    write(*,20) '*   the electron k and phonon q are not                 '
    write(*,20) '*   conmesurate and the k grid does not contain         '
    write(*,20) '*   the phonon q grid                                   '
    write(*,20) '**********************************************************'
    stop
  endif
  !
  !
  !================================================================================
  ! Read the parameters from the SCF calculation
  !================================================================================
  !
  write(*,20) '| - Reading calculation parameters...               |'
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
  ! Select bands for which ep elements will be calculated
  !
  if (trim(ep_bands) .eq. 'custom') then
    num_bands_ep = ep_bands_final-ep_bands_initial+1
    write(*,*) ' ep_bands == custom chosen.'
    write(*,*) ' ep elements to be calculated for bands ', ep_bands_initial, &
               ' to ', ep_bands_final, ' of the num_bands_intw list'
  else ! all available bands
    num_bands_ep = num_bands_intw
  end if
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
  ! Calculate the multiplication talble for symmetry operations
  allocate(symtable(nsym,nsym))
  call multable(nsym, s, symtable)
  !
  ! Set up spin_symmetry_matrices, needed to rotate wave functions and indueced potential for non-colinear calculations
  call allocate_and_build_spin_symmetry_matrices(nsym)
  !
  write(*,20) '|         ---------------------------------         |'
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
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Read PPs
  !================================================================================
  !
  write(*,20) '| - Reading pseudopotentials...                     |'
  !
  call read_all_pseudo()
  !
  write(*,20) '|                    PPs are OK                     |'
  !
  ! Allocate and set PP variables
  call init_KB_PP()
  !
  write(*,20) '|         ---------------------------------         |'
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
  call find_size_of_irreducible_k_set(nk1, nk2, nk3, kmesh_nkirr)
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
                              QE_folder_nosym, QE_folder_sym, symlink, &
                              full_mesh, IBZ)
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
    write(*,20) '**********************************************************'
    write(*,20) '* The kpoints present in the QE folders are not consistent'
    write(*,20) '* with the parameters of the input file!                 '
    write(*,20) '**********************************************************'
    write(*,20) '* debug information:                                *'
    write(*,*) '*        nkpoints_QE = ', nkpoints_QE
    write(*,*) '*        nkmesh      = ', nkmesh
    write(*,*) '*        kmesh_nkirr = ', kmesh_nkirr
    stop
  end if
  !
  !
  !================================================================================
  ! Read phonon information
  !================================================================================
  !
  write(*,20) '| - Reading q-points...                             |'
  !
  ! Read q-points and irreducible patterns
  call read_ph_information()
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
  ! Find the size of the irreducible set of q-points (IBZ)
  call find_size_of_irreducible_k_set(nq1, nq2, nq3, qmesh_nqirr)
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
  call set_symmetry_relations(nq1,nq2,nq3, nqirr, q_irr_cryst, &
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
    write(*,20) '|         ---------------------------------         |'
  else if(IBZ_q) then
    write(*,20) '| - The qpoints present in the QE folders           |'
    write(*,20) '|   are consistent with an IBZ.                     |'
    write(*,20) '|         ---------------------------------         |'
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
  ! Read the induced potentials
  !================================================================================
  !
  write(*,20) '| - Reading induced potentials...                   |'
  !
  call read_allq_dvr()
  !
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Compute matrix elements
  !================================================================================
  !
  write(*,20) '| - Computing matrix elements...                    |'
  !
  ! Allocate wfc related variables
  allocate(list_igk(nGk_max))
  allocate(list_igkq(nGk_max))
  allocate(wfc_k(nGk_max,num_bands_intw,nspin))
  allocate(wfc_kq(nGk_max,num_bands_intw,nspin))
  !
  ! We will calculate num_bands_ep bands if ep_bands=custom
  ! indexed 1:num_bands_ep (ep_bands_initial to ep_bands_final)
  !
  ! Allocate induced potential related variables
  allocate(dvq_local(nr1*nr2*nr3,3*nat,nspin,nspin))
  allocate(dvpsi(nGk_max,num_bands_ep,nspin,nspin,3*nat))
  !
  ! Allocate matrix elements variable
  allocate(ep_mat_el(nkmesh,num_bands_ep,num_bands_ep,nspin,nspin,3*nat))
  !
  do iq=1,nqmesh
    !
    ep_mat_el = cmplx_0
    !
    qpoint = qmesh(:,iq)
    write(*,"(a,i4,a,a)") "|                    qpoint ", iq, "                    |"
    write(*,"(a,3f10.5,a)") "|         q =", qpoint, "         |"
    !
    if (                iq <   10) write(iq_loc,"(i1)") iq
    if ( 10 <= iq .and. iq <  100) write(iq_loc,"(i2)") iq
    if (100 <= iq .and. iq < 1000) write(iq_loc,"(i3)") iq
    !
    ep_unit = find_free_unit()
    inquire(iolength=record_length) ep_mat_el
    open(unit=ep_unit, iostat=ierr, &
         file=trim(outdir)//trim(ep_mat_file)//trim('_')//adjustl(iq_loc), &
         form='unformatted', status='unknown', access='direct', recl=record_length)
    if (ierr /= 0 ) stop 'Error opening ep_mat_file'
    !
    ! Get induced potential for q
    dvq_local = cmplx_0
    call get_dv(qpoint, dvq_local)
    !
    ! Add local part of the KB-PP to the induced potential
    call calculate_local_part_dv(qpoint, dvq_local)
    !
    do ik=1,nkmesh
      !
      kpoint = kmesh(:,ik)
      !
      write(*,'(a,i4,a,3f10.5,a)') "|   kpoint ", ik, ': k =', kpoint, "  |"
      !
      ! Get wave functions for k and k+q
      call get_psi_general_k_all_wfc(       kpoint,  nGk,  list_iGk,  wfc_k)
      call get_psi_general_k_all_wfc(kpoint+qpoint, nGkq, list_iGkq, wfc_kq)
      !
      !
      ! Multiply induced potential + local part of KB-PP with wave function:
      ! dvpsi: dv_q^local x | psi_k > (G)
      !
      if ( trim(ep_bands) .eq. 'intw') then
        call dvqpsi_local(num_bands_intw, list_iGk, list_iGkq, wfc_k, dvq_local, dvpsi)
      else if ( trim(ep_bands) .eq. 'custom') then
        call dvqpsi_local(num_bands_ep, list_iGk, list_iGkq, &
                          wfc_k(:,ep_bands_initial:ep_bands_final,:), dvq_local, dvpsi)
      end if
      !
      ! Add non-local part of KB-PP multiplied with the wave function:
      ! dvpsi: --> dvpsi + dq^mode [ KB projectors ] x | psi_k > (G)
      !            local + non-local
      !
      if ( trim(ep_bands) .eq. 'intw') then
        call multiply_psi_by_dvKB(kpoint, qpoint, list_iGk, list_iGkq, num_bands_intw, wfc_k, dvpsi)
      else if ( trim(ep_bands) .eq. 'custom') then
        call multiply_psi_by_dvKB(kpoint, qpoint, list_iGk, list_iGkq, &
                                  num_bands_ep,  wfc_k(:,ep_bands_initial:ep_bands_final,:), dvpsi)
      end if
      !
      do imode=1,3*nat ! This are Cartesian modes, not real phonon modes
        !
        ! Compute matrix elements:
        ! ep_mat_el: < psi_{k+q} | x dvpsi
        do jspin=1,nspin
          do ispin=1,nspin
            do jbnd=1,num_bands_ep
              do ibnd=1,num_bands_ep
                !
                ep_mat_el(ik,ibnd,jbnd,ispin,jspin,imode) = zdotc( nGk_max, wfc_kq(:,ibnd,ispin), 1, dvpsi(:,jbnd,ispin,jspin,imode), 1 )
                !
              enddo !ibnd (intw or custom)
            enddo !jbnd (intw or custom)
          enddo !is
        enddo !js
        !
      enddo !imode
      !
    enddo !ik
    !
    write(unit=ep_unit, rec = 1, iostat = ierr) ep_mat_el(:,:,:,:,:,:)
    !
    close(unit=ep_unit)
    !

#ifdef DEBUG
      ! Write the matrix elements to a formatted file to compare with QE matrix elements
      write(123400+iq,"(3f10.6)") qpoint
      write(123400+iq,"(i4)") nkmesh
      do ik=1,nkmesh
        kpoint = kmesh(:,ik)
        write(123400+iq,"(i4,3f10.6)") ik, kpoint
        write(123400+iq,"(i4,3f10.6)") 2*(ik-1)+2, kpoint+qpoint
      enddo
      write(123400+iq,"(4i4,2f20.15)")((((ik, ibnd, jbnd, imode, sum(ep_mat_el(ik,ibnd,jbnd,:,:,imode)), ibnd=1,num_bands_intw), jbnd=1,num_bands_intw), ik=1,nkmesh), imode=1,3*nat)
#endif

    !
  enddo !iq

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

end program ep_melements
