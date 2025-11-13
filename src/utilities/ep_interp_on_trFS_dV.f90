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
program ep_on_trFS_dV

  ! MBR 24/04/2024
  !
  ! This utility calculates electron-phonon matrix elements interpolated on a Fermi surface triangulation
  ! following the method of:
  !
  ! A. Eiguren, C. Ambrosch-Draxl, Phys. Rev. B 78, 045124 (2008)
  !
  ! For that, the induced potentials calculated with QE on the q-grid are read in
  ! and interpolated to the non-uniform q-list corresponding to the triangulation (q=k'-k)
  ! Needed wavefunctions on the triangulated k-points are calculated on the fly by calling QE and stored.
  ! If such calculations preexist, those wavefunctions are read.
  !
  ! 1st part:
  !
  !    This utility uses < intw.in as others.
  !    It uses outdir, prefix, and info on which FS sheets are needed.
  !    Unlike with Wannier g-element interpolation, here
  !    the location of pw.x and pw2intw.x executables is needed,
  !    which is also given in the intw.in file.
  !
  !    1.- From files prefix.N_FS_tri.off, where N is the Fermi surface (FS) sheet,
  !        it reads the lists of k-points at the vertices, and the faces, in order to calculate
  !        areas.
  !
  !    2.- It writes a prefix-nscf.nscf.in input for pw.x by reading prefix.scf.in
  !        and replacing the k-points by the FS triangularization information.
  !
  !    3.- It writes a prefix-nscf.pw2intw.in file, which will convert wfc data from
  !        prefix-nscf.save directory.
  !
  !    4.- It sets a directory prefix-nscf.save and initializes it with the scf charge density.
  !        It runs pw.x and applies pw2intw.x to the obtained wfcs. In case the prefix-nscf.save/wfc
  !        files exist, it skips this and warns the user.
  !
  !    (In the case of using SIESTA step 2 and 3 are done with siesta2intw.x)
  !
  ! 2nd part:
  !
  !    Read derivative of local potential (dvq_local) on the real space unit cell for
  !    q-vectors of full BZ and Fourier transform to R (Wigner-Seitz)
  !    as the first step of the interpolation.
  !    The non-local part will be added on the fly when interpolating over the triangulated k-points.
  !
  ! 3rd part:
  !
  !    Loop over FS (one over irreducible, one over full).
  !    Interpolate the potential.
  !    Calculate ep elements (local + non-local) as done in ep_melements.f90 utility.
  !    Write to file

  use kinds, only: dp

  use intw_version, only: print_intw_version

  use intw_useful_constants, only: cmplx_1, cmplx_0, cmplx_i, Ha_to_eV, tpi, eps_8

  use intw_utility, only: get_timing, print_threads, print_date_time, find_free_unit, &
                          cryst_to_cart, generate_kmesh, &
                          joint_to_triple_index_r

  use intw_matrix_vector, only: area_vec, ainv

  use intw_input_parameters, only: outdir, prefix, read_input, &
                                   intw2W_method, intw2W_fullzone, nk1, nk2, nk3, chemical_potential, &
                                   nq1, nq2, nq3, nqirr, ph_dir, &
                                   ep_interp_method, ep_interp_bands, nfs_sheets_initial, nfs_sheets_final, &
                                   nscf_code

  use intw_reading, only: num_bands_intw, set_num_bands, read_parameters_data_file, &
                          nGk_max, get_gvec, &
                          nspin, lspin, &
                          at, nr1, nr2, nr3, &
                          nat, tau, nsym, s, &
                          get_K_folder_data

  use intw_pseudo, only: read_all_pseudo

  use intw_pseudo_local, only: calculate_local_part_dv, dvqpsi_local

  use intw_pseudo_non_local, only: init_KB_PP, &
                                   multiply_psi_by_dvKB

  use intw_symmetries, only: set_symmetry_relations, multable, rot_atoms, &
                             symtable, rtau_index, rtau, rtau_cryst, &
                             allocate_and_build_spin_symmetry_matrices, &
                             find_inverse_symmetry_matrices_indices, &
                             find_size_of_irreducible_k_set

  use intw_fft, only: generate_nl, allocate_fft

  use intw_ph, only: nqmesh, qmesh, read_ph_information, &
                     q_irr_cryst, read_allq_dvr, get_dv, &
                     QE_folder_nosym_q, QE_folder_sym_q, symlink_q

  use intw_ph_interpolate, only: irvec_q, nrpts_q, ndegen_q, allocate_and_build_ws_irvec_q

  implicit none

  ! for part I
  logical :: read_status
  character(5) :: is_loc, comenta
  character(100) :: file_off
  integer :: unit_off
  integer :: nkpt_tr_tot, nkpt_tr_ibz_tot
  integer :: is, js, ik, ik1, i, j, k, iface, iks, ish, iksp, ishp, ibp
  integer :: nfs_sheets_tot ! number of sheets considered
  integer, allocatable :: nfs_sheet(:), & ! band indices of the sheets (num_bands_intw set)
                          nkpt_tr(:), & ! number of kpoints in each FS sheet
                          nkpt_tr_ibz(:), & ! number of kpoints in each FS sheet irreducible BZ wedge
                          nface_tr(:) ! number of faces in each FS sheet
  real(dp) :: k1(3), k2(3), k3(3), kwei
  real(dp), allocatable :: kpts_tr(:,:), kpts_tr_area(:)

  ! for part II
  logical :: full_mesh_q, IBZ_q
  integer :: qmesh_nqirr
  integer :: iq, ir, jr, ir1, ir2, ir3
  real(dp) :: qpoint(3), rvec(3)
  complex(dp) :: facq
  complex(dp), allocatable :: dvq_local_R(:,:,:,:,:)

  ! for part III
  logical :: have_ep
  character(256) :: altprefix, file_ep
  integer :: ib, iat, ikp, unit_ep
  integer :: nG, nG_p
  integer, allocatable :: list_iG(:), list_iG_p(:)
  real(dp) :: kpoint(3), kpoint_p(3)
  real(dp), allocatable :: QE_eig(:), QE_eig_p(:)
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:)
  complex(dp), allocatable :: wfc(:,:,:), wfc_p(:,:,:)
  complex(dp), allocatable :: dvq_local(:,:,:,:)
  complex(dp), allocatable :: dvpsi(:,:,:,:,:)

  ! timing
  real(dp) :: time1, time2

  complex(dp), external :: zdotc


  20 format(A)
  30 format(A,F8.2,6X,A)

  !================================================================================
  ! Beginning
  !================================================================================

  call get_timing(time1)

  write(*,20) '====================================================='
  write(*,20) '|               program ep_on_triFS_dV              |'
  write(*,20) '|         e-p calculation on triangulated FS        |'
  write(*,20) '|  using interpolation of dV on the triangulated q  | '
  write(*,20) '|         ---------------------------------         |'
  call print_intw_version()
  call print_threads()
  call print_date_time("Start of execution")
  write(*,20) '====================================================='


  !================================================================================
  ! Read the input file
  !================================================================================

  call read_input(read_status)

  if (read_status) stop

  ! check right method is chosen in intw.in
  if ( ep_interp_method /= 'dV_interpolate' ) then
    write(*,*) 'ep_interp_method /= dV_interpolate in input. Stopping.'
    stop
  end if


  !================================================================================
  ! Read the parameters from the SCF calculation
  !================================================================================

  write(*,20) '| - Reading calculation parameters...               |'

  call read_parameters_data_file()


  !================================================================================
  ! Set the number of bands
  !================================================================================

  ! this will read num_wann_intw and num_bands_intw dimensions
  call set_num_bands()


  !================================================================================
  ! Print spin information
  !================================================================================

  if (lspin) then
    write(*,20) '| - Spin-polarized calculation nspin = 2            |'
  else
    write(*,20) '| - Paramagnetic calculation nspin = 1              |'
  endif

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Set symmetry arrays
  !================================================================================

  write(*,20) '| - Setting symmetry arrays...                      |'

  ! Set the rotation table for each atom and symmetry
  allocate(rtau_index(nat,nsym))
  allocate(rtau(3,nsym,nat))
  allocate(rtau_cryst(3,nsym,nat))

  call rot_atoms(nat, nsym, tau)

  ! Compute the indices of the inverse rotation matrices
  call find_inverse_symmetry_matrices_indices()

  ! Calculate the multiplication talble for symmetry operations
  allocate(symtable(nsym, nsym))
  call multable(nsym, s, symtable)

  ! Set up spin_symmetry_matrices, needed to rotate wave functions and indueced potential for non-colinear calculations
  call allocate_and_build_spin_symmetry_matrices(nsym)

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Set up the gvec array and all FFT variables
  !================================================================================

  write(*,20) '| - Reading G vectors...                            |'

  call get_gvec()

  ! Allocate useful variables
  call allocate_fft()

  ! Generate some important indices for FFT
  call generate_nl()

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read PPs
  !================================================================================

  write(*,20) '| - Reading pseudopotentials...                     |'

  call read_all_pseudo()

  write(*,20) '|                    PPs are OK                     |'

  ! Allocate and set PP variables
  call init_KB_PP()


  !================================================================================
  ! Choose Fermi surface sheets according to ep_interp_bands
  !================================================================================

  if ( ep_interp_bands == 'intw_bands' ) then
    nfs_sheets_tot = num_bands_intw
    allocate(nfs_sheet(nfs_sheets_tot))
    do ib = 1, num_bands_intw
      nfs_sheet(ib) = ib
    end do
  else if ( ep_interp_bands == 'ef_crossing' ) then
    nfs_sheets_tot = nfs_sheets_final - nfs_sheets_initial + 1
    allocate(nfs_sheet(nfs_sheets_tot))
    do ib = 1, nfs_sheets_tot
      nfs_sheet(ib) = nfs_sheets_initial + ib-1
    end do
  end if


  !******************************** Part I ************************************

  !================================================================================
  ! Read .off files
  !================================================================================

  write(*,20) '====================================================='
  write(*,20) '| - Reading .off files...                           |'

  allocate(nkpt_tr(nfs_sheets_tot), nface_tr(nfs_sheets_tot))
  allocate(nkpt_tr_ibz(nfs_sheets_tot))

  ! open all sheet files just to see dimensions of kpoint lists
  do is=1,nfs_sheets_tot

    if (                is <   10) write(is_loc,"(i1)") nfs_sheet(is)
    if ( 10 <= is .and. is <  100) write(is_loc,"(i2)") nfs_sheet(is)

    file_off = trim(outdir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_tri.off')
    write(*,'(A)') '|     '//file_off(1:max(45,len(trim(file_off))))//' |'

    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')
    read(unit_off,*) comenta
    read(unit_off,*) nkpt_tr(is), nface_tr(is), k ! number of vertices and faces (ignore edges)
    close(unit_off)

    ! open the IBZ off file and search for dimension nkpt_tr_ibz(is).
    ! Its vertices coincide with the first nkpt_tr_ibz(is) vertices of the full off vertex list.

    file_off = trim(outdir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_IBZ_FS_tri.off')

    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')
    read(unit_off,*) comenta
    read(unit_off,*) nkpt_tr_ibz(is), j, k ! number of vertices (ignore faces and edges)
    close(unit_off)

  end do

  ! total number of k-points to be calculated
  nkpt_tr_tot = sum(nkpt_tr)
  nkpt_tr_ibz_tot = sum(nkpt_tr_ibz)

  write(*,20) '|   Number of k-points (total vertices):            |'
  write(*,'(A1,3X,I6,42X,A1)') '|', nkpt_tr_tot, '|'
  write(*,20) '|   Number of k-points in IBZ (total vertices):     |'
  write(*,'(A1,3X,I6,42X,A1)') '|', nkpt_tr_ibz_tot, '|'
  write(*,20) '|   Number of faces (total triangles):              |'
  write(*,'(A1,3X,I6,42X,A1)') '|', sum(nface_tr), '|'


  allocate(kpts_tr(3,nkpt_tr_tot), kpts_tr_area(nkpt_tr_tot))


  ! open .off files again to read k-points
  ik1 = 0
  do is=1,nfs_sheets_tot

    if (                is <   10) write(is_loc,"(i1)") nfs_sheet(is)
    if ( 10 <= is .and. is <  100) write(is_loc,"(i2)") nfs_sheet(is)

    ! .off file for this sheet

    file_off = trim(outdir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_tri.off')
    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')

    read(unit_off,*) comenta
    read(unit_off,*) i, j, k ! number vertices, faces and edges (I will ignore edges)
    ! read(unit_off,'(/)') ! DUDA... This will depend on how the line break is written in the file, I think...
    if ( (i /= nkpt_tr(is)) .or. (j /= nface_tr(is)) ) then
      write(*,*) 'Error reading ', file_off, '. Stopping.'
      stop
    end if

    ! read vertices
    do ik=1,nkpt_tr(is)
      read(unit_off,*) kpts_tr(:,ik1+ik) ! units in the trFS.off file are cartesian 2pi/alat ("tpiba" for QE)
    end do

    ! Read (triangular) faces on this sheet.
    ! Each face contributes with 1/3 of its area to the effective area of each of its vertices.
    ! Calculate the are on the go and add the contribution to each vertex, storing for global indices (i.e. ik1+ik).
    do iface = 1,nface_tr(is)
      read(unit_off,*) i, ir1, ir2, ir3 ! indices ik of the vertices of the face, indexed from 0
      ir1 = ir1 + 1
      ir2 = ir2 + 1
      ir3 = ir3 + 1 ! now, ik of the vertices of the face, indexed from 1
      if ( i /= 3 ) then
        write(*,*) 'Error reading ', file_off, 'Only triangles allowed. Stopping.'
        stop
      end if
      ! triangle vertex vectors (cartesian 2pi/alat)
      k1 = kpts_tr(:,ik1+ir1)
      k2 = kpts_tr(:,ik1+ir2)
      k3 = kpts_tr(:,ik1+ir3)
      ! get spanned area and add contribution to each vertex
      ! function copied from FSH/modules/geometry.f90
      kwei = area_vec(k2-k1,k3-k1)/3.0_dp
      kpts_tr_area(ik1+ir1) = kpts_tr_area(ik1+ir1) + kwei
      kpts_tr_area(ik1+ir2) = kpts_tr_area(ik1+ir2) + kwei
      kpts_tr_area(ik1+ir3) = kpts_tr_area(ik1+ir3) + kwei
    end do

    close(unit_off)

    ! accumulate ik global index for the reading of next sheet
    ik1 = ik1 + nkpt_tr(is)

  end do

  write(*,20) '|   .... reading done                               |'

  write(*,20) '|   Total FS area:                                  |'
  write(*,'(A1,3X,F12.6,A19,17X,A1)') '|', sum(kpts_tr_area), ' (2 x pi / alat)^2 ', '|'

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Run nscf calculation for the FS k-points
  !================================================================================

  if ( nscf_code == "QE" ) then
    call QE_nscf()
  else if ( nscf_code == "SIESTA" ) then
    call SIESTA_nscf()
  else
    stop "ERROR: Invalid value for nscf_code variable."
  endif

  write(*,20) '|                                                   |'
  write(*,20) '| ---------------- Part I completed --------------- |'
  write(*,20) '|                                                   |'
  write(*,20) '====================================================='


  !******************************** Part II ************************************

  !================================================================================
  ! Read phonon information
  !================================================================================

  write(*,20) '| - Reading phonon info...                          |'

  ! Read irreducible q-points and irreducible patterns
  call read_ph_information()

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Build the phonon q-mesh
  !================================================================================

  write(*,20) '| - Building q-mesh...                              |'

  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))

  call generate_kmesh(qmesh, nq1, nq2, nq3)

  ! Find the size of the irreducible set of q-points (IBZ)
  call find_size_of_irreducible_k_set(nq1, nq2, nq3, qmesh_nqirr)


  !================================================================================
  ! Set symmetry relations between irreducible q-points and full q-mesh
  !================================================================================

  allocate(QE_folder_nosym_q(nqmesh))
  allocate(QE_folder_sym_q(nqmesh))
  allocate(symlink_q(nqmesh,2))

  call set_symmetry_relations(nq1, nq2, nq3, nqirr, q_irr_cryst, &
                              QE_folder_nosym_q, QE_folder_sym_q, symlink_q, &
                              full_mesh_q, IBZ_q)


  !================================================================================
  ! Check that the number of q-points corresponds to either a full mesh or the IBZ
  !================================================================================

  if (full_mesh_q .and. IBZ_q) then
    write(*,20) '| - The qpoints present in the QE folders           |'
    write(*,20) '|   are consistent with a full 1BZ and a            |'
    write(*,20) '|   IBZ has also been found.                        |'
  else if(IBZ_q) then
    write(*,20) '| - The qpoints present in the QE folders           |'
    write(*,20) '|   are consistent with an IBZ.                     |'
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


  !================================================================================
  ! Obtain dV(q) (as in the ep_elements utility written by Haritz).
  !================================================================================

  write(*,20) '| - Reading induced potentials...                   |'

  call read_allq_dvr()

  write(*,20) '====================================================='


  !================================================================================
  ! Fourier transform a la Wannier of q->R of dvq_local using the Wigner-Seitz supercell grid.
  ! Loop over whole mesh: obtain dvq_local for each q-point and add it to the Fourier transform
  !================================================================================

  ! Wigner-Seitz R points
  write(*,20) '| - Building WS mesh...                             |'

  call allocate_and_build_ws_irvec_q()

  write(*,20) '| - Fourier transform induced potential...          |'

  ! Potential in the supercell using Wigner-Seitz R points
  allocate(dvq_local_R(nrpts_q,nr1*nr2*nr3,3*nat,nspin,nspin))

  allocate(dvq_local(nr1*nr2*nr3,3*nat,nspin,nspin))
  dvq_local_R = cmplx_0

  do iq = 1, nqmesh

    qpoint = qmesh(:,iq)
    write(*,'(A12,I5,3f11.5,A3)') '|     qpoint', iq, qpoint, "  |"

    ! Compute induced potential using symmetries
    dvq_local = cmplx_0
    call get_dv(qpoint, dvq_local)

    ! Compute local part of the derivative of the KB PP
    call calculate_local_part_dv(qpoint, dvq_local)

    ! transform with phase: iq*(r-R)

    do ir = 1,nr1*nr2*nr3 ! unit cell coordinates(1:nr1)
      ! ir to (ir1,ir2,ir3), 3rd index fastest
      call joint_to_triple_index_r(nr1, nr2, nr3, ir, ir1, ir2, ir3)
      rvec(1) = real(ir1-1,dp)/real(nr1,dp) ! r-vector in fractional coord.
      rvec(2) = real(ir2-1,dp)/real(nr2,dp)
      rvec(3) = real(ir3-1,dp)/real(nr3,dp)
      do jr = 1,nrpts_q
        facq = exp(cmplx_i*tpi*dot_product(qpoint, rvec(:)-irvec_q(:,jr)))
        dvq_local_R(jr,ir,:,:,:) = dvq_local_R(jr,ir,:,:,:) + facq * dvq_local(ir,:,:,:)
      end do ! R in WS
    end do ! e in unit cell

  end do
  dvq_local_R = dvq_local_R / real(nqmesh,dp) ! normalize Fourier transform

  write(*,20) '|                                                   |'
  write(*,20) '|  --------------- Part II completed -------------- |'
  write(*,20) '|                                                   |'
  write(*,20) '====================================================='

  ! DUDA dvq_local_R in the supercell should be real?


  !******************************* Part III *************************************

  write(*,20) '| - Interpolating e-p elements...                   |'

  ! For the a2F integral we will need a double loop on kpoints of the FS triangulation.
  ! In this part the final ep element is calculated on the go for q=k'-k after
  ! backtransform of dV local:
  ! dvq_local(k+q,k) = dvq_local(k',k)

  ! allocate dvpsi for only one band, as band index is included in the k indices and thus
  ! run one-by-one in the loops below
  allocate(dvpsi(nGk_max,1,nspin,nspin,3*nat))

  allocate(aep_mat_el(nkpt_tr_tot,nkpt_tr_ibz_tot,nspin,nspin,3*nat))

  allocate(list_iG(nGk_max), list_iG_p(nGk_max))
  allocate(QE_eig(num_bands_intw), QE_eig_p(num_bands_intw))
  allocate(wfc(nGk_max,num_bands_intw,nspin), wfc_p(nGk_max,num_bands_intw,nspin))


  altprefix = trim(prefix)//'-nscf'
  file_ep = trim(outdir)//trim(prefix)//trim('_ep_interp.dat')

  inquire(file=file_ep,exist=have_ep)

  if (.not.have_ep) then ! calculate interpolated ep elements and write to file_ep

    unit_ep = find_free_unit()
    open(unit_ep, file=file_ep, status='unknown')
    write(unit_ep,*)'# ik(irr)   jk(full)    is js   g(canonical modes)'

    ! ik, ikp indices implicitly contain the FS sheet index, i.e. the band indices ib, ib'
    ! to be selected, so instead of iterating over nkpt_tr_tot, I separate over sheets
    ! (aep_mat_el elements are stored only for the needed pair of sheets for a given kk')

    ik = 0
    ik1 = 0
    do ish = 1, nfs_sheets_tot

      ib = nfs_sheet(ish) ! band index for k

      do iks = 1, nkpt_tr_ibz(ish)

        ! ik is the k-index over nkpt_tr_tot in the Irreducible BZ
        ! ik1 is the corresponding k-index in the full BZ kpts_tr list (at the end of the loop I add the rest of nkpt_tr(ish))
        ik = ik + 1
        ik1 = ik1 + 1

        write(*, '(A14,I4,A1,I4,A6,I5,A19)') '|     ik_IBZ: ', ik, "/", nkpt_tr_ibz_tot, ' (ik: ', ik1, ")                 |"

        kpoint = kpts_tr(:,ik1) ! this is cartesians x 2pi/alat. Transform to cryst.
        call cryst_to_cart(1, kpoint, at, -1)

        ! Read wavefunction k
        call get_K_folder_data(ik1, list_iG, wfc, QE_eig, nG, altprefix)

        do ishp = 1, nfs_sheets_tot

          ibp = nfs_sheet(ishp) ! band index for k'

          !$omp parallel do &
          !$omp firstprivate(iks, ik, kpoint, list_iG, wfc, QE_eig, nG) &
          !$omp private(iksp, ikp, kpoint_p, list_iG_p, wfc_p, QE_eig_p, nG_p) &
          !$omp private(ir, jr, ir1, ir2, ir3, rvec, facq, dvq_local) &
          !$omp private(qpoint, dvpsi) &
          !$omp private(iat, is, js)
          do iksp = 1, nkpt_tr(ishp)

            ikp = iksp + sum(nkpt_tr(:ishp-1)) ! k'-index over nkpt_tr_tot

            kpoint_p = kpts_tr(:,ikp) ! this is cartesians x 2pi/alat. Transform to cryst.
            call cryst_to_cart(1, kpoint_p, at, -1)

            ! Read wavefunction k'
            !$omp critical
            call get_K_folder_data(ikp, list_iG_p, wfc_p, QE_eig_p, nG_p, altprefix)
            !$omp end critical

            qpoint = kpoint_p-kpoint

            ! Fourier backtransform a la Wannier: interpolation of dV local at qpoint
            ! with e^{-iq.(r-R)}, where R is a lattice vector of the WS supercell, as
            ! we would do in Wannier.

            dvq_local = cmplx_0 ! watch out: I am reusing this array in this step of k,k' double loop
            do ir = 1,nr1*nr2*nr3 ! unit cell coordinates(1:nr1)
              ! ir to (ir1,ir2,ir3), 3rd index fastest
              call joint_to_triple_index_r(nr1, nr2, nr3, ir, ir1, ir2, ir3)
              rvec(1) = real(ir1-1,dp)/real(nr1,dp) ! r-vector in fractional coord.
              rvec(2) = real(ir2-1,dp)/real(nr2,dp)
              rvec(3) = real(ir3-1,dp)/real(nr3,dp)
              do jr = 1,nrpts_q
                facq = exp(-cmplx_i*tpi*dot_product(qpoint, rvec(:)-irvec_q(:,jr)))/real(ndegen_q(jr),dp)
                dvq_local(ir,:,:,:) = dvq_local(ir,:,:,:) + facq * dvq_local_R(jr,ir,:,:,:)
              end do ! R in WS
            end do ! r in unit cell

            ! Multiply psi_k with induced potential + local part of the KB PP:
            !   dvpsi_{k+q} = dv_local_{q} x |psi_{k}>
            call dvqpsi_local(1, list_iG, list_iG_p, wfc(:,ib:ib,:), dvq_local, dvpsi)

            ! Add non-local contribution: Multiply psi_k with non-local part of the KB PP:
            !   dvpsi_{k+q} --> dvpsi_{k+q} + d_{q} [ KB ] |psi_{k}>
            !                   (local)     + (non-local)
            call multiply_psi_by_dvKB(kpoint, qpoint, list_iG, list_iG_p, 1, wfc(:,ib:ib,:), dvpsi)

            ! Compute matrix elements: Multiply dvpsi with psi_{k+q}
            do iat=1,3*nat
              do js=1,nspin
                do is=1,nspin
                  aep_mat_el(ikp,ik, js,is, iat) = &
                            zdotc( nGk_max, wfc_p(:,ibp,js), 1, dvpsi(:,1,js,is,iat), 1 )
                enddo ! is
              enddo ! js
            enddo ! iat

          end do ! k'
          !$omp end parallel do

          ! Save interpolated matrix elements for kpoint
          do iksp = 1, nkpt_tr(ishp)
            !
            ikp = iksp + sum(nkpt_tr(:ishp-1)) ! k'-index over nkpt_tr_tot
            !
            do js=1,nspin
              do is=1,nspin
                write(unit_ep,fmt="(6i6,100e16.6)") ibp, iksp, ikp, ib, iks, ik,  &
                    (aep_mat_el(ikp,ik, js,is, iat), iat=1,3*nat)
              end do
            end do
            !
          end do ! k'

        end do ! sheet'

      end do ! k

      ik1 = ik1 + nkpt_tr(ish) - nkpt_tr_ibz(ish)

    end do ! sheet

    close(unit_ep)

    write(*,20) '| - e-p elements interpolated and written to file:  |'
    write(*,20) "|   "//file_ep(1:max(47,len(trim(file_ep))))//" |"

  else ! have_ep--> read interpolated matrix elements

    write(*,20) '| - interpolated e-p file already exists. Check!    |'
    write(*,20) "|   "//file_ep(1:max(47,len(trim(file_ep))))//" |"

  end if

  write(*,20) '|                                                   |'
  write(*,20) '|  -------------- Part III completed -------------- |'
  write(*,20) '|                                                   |'
  write(*,20) '====================================================='


  !================================================================================
  ! Finish
  !================================================================================

  call get_timing(time2)

  write(*,20) '|                      ALL DONE                     |'
  write(*,30) '|     Total time: ',time2-time1,' seconds            |'
  call print_date_time('End of execution  ')
  write(*,20) '====================================================='


contains

  subroutine QE_nscf()

    implicit none

    character(100) :: file_nscf_in, file_nscf_out
    character(100) :: file_pw2intw_in, file_pw2intw_out
    character(100) :: dir_nscf, dir_nscf_intw
    logical :: have_nscf_out, have_nscf_wfcN, have_nscf_wfc1
    character(5) :: ik_loc

    !================================================================================
    ! Read prefix.scf.in file line by line and dump information to prefix-nscf.nscf.in.
    ! For that, the lines containing calculation and prefix are modified.
    ! When arriving at KPOINTS, reading stops and the triangle list is dumped to the
    ! prefix-nscf.nscf.in file.
    ! Note: this is done only in case the nscf files do not exist (otherwise,
    ! it is assumed that the wfcs have already been calculated and transformed to intw format)
    !================================================================================

    !================================================================================
    ! Check pw2intw
    ! Check if first and last wfc are present. If so, it probably means that the nscf
    ! calculation was already made, so skip this step
    !================================================================================

    file_pw2intw_in = trim(outdir)//trim(prefix)//'-nscf.pw2intw.in'
    file_pw2intw_out = trim(outdir)//trim(prefix)//'-nscf.pw2intw.out'
    dir_nscf_intw = trim(outdir)//trim(prefix)//"-nscf.save.intw/"

    inquire(file=file_pw2intw_out, exist=have_nscf_out)
    write(ik_loc,"(i5.5)") 1
    inquire(file=trim(dir_nscf_intw)//'wfc'//trim(adjustl(ik_loc))//'.dat', exist=have_nscf_wfc1)
    write(ik_loc,"(i5.5)") nkpt_tr_tot
    inquire(file=trim(dir_nscf_intw)//'wfc'//trim(adjustl(ik_loc))//'.dat', exist=have_nscf_wfcN)

    if (have_nscf_out .and. have_nscf_wfc1 .and. have_nscf_wfcN) then
      write(*,'(A)') '| - INTW nscf wfcs seem to be already available.    |'
      write(*,'(A)') '|   Skipping nscf calculation!                      |'
      return
    endif

    !================================================================================
    ! Check QE nscf
    ! Check if first and last wfc are present. If so, it probably means that the nscf
    ! calculation was already made, so skip this step
    !================================================================================

    file_nscf_in = trim(outdir)//trim(prefix)//'-nscf.nscf.in'
    file_nscf_out = trim(outdir)//trim(prefix)//'-nscf.nscf.out'
    dir_nscf = trim(outdir)//trim(prefix)//"-nscf.save/"

    inquire(file=file_nscf_out, exist=have_nscf_out)
    write(ik_loc,"(i5)") 1
    inquire(file=trim(dir_nscf)//'wfc'//trim(adjustl(ik_loc))//'.dat', exist=have_nscf_wfc1)
    write(ik_loc,"(i5)") nkpt_tr_tot
    inquire(file=trim(dir_nscf)//'wfc'//trim(adjustl(ik_loc))//'.dat', exist=have_nscf_wfcN)

    if (have_nscf_out .and. have_nscf_wfc1 .and. have_nscf_wfcN) then
      write(*,'(A)') '| - QE nscf wfcs seem to be already available.      |'
      write(*,'(A)') '|   Skipping QE nscf calculation!                   |'
      call run_pw2intw(file_pw2intw_in, file_pw2intw_out)
    else
      write(*,'(A)') '| - Calculating nscf wfcs with QE...                |'
      call run_pw(file_nscf_in, file_nscf_out, dir_nscf)
      call run_pw2intw(file_pw2intw_in, file_pw2intw_out)
    end if

  end subroutine QE_nscf

  subroutine run_pw(file_nscf_in, file_nscf_out, dir_nscf)
    ! Create the input file and run the nscf calculation with pw.x

    use intw_input_parameters, only: command_pw, file_pw

    implicit none

    character(100), intent(in) :: file_nscf_in, file_nscf_out
    character(100), intent(in) :: dir_nscf

    logical :: exists_charge_density, exists_schema, exists_charge_density_qe, exists_schema_qe
    character(70) :: line, lleft
    character(100) :: dir_scf, dir_scf_intw
    character(200) :: comando, datafile
    integer :: unit_pw, unit_nscf, ios
    integer :: ikpt


    ! Create input file for nscf calculation

    unit_pw = find_free_unit()
    open(unit_pw, file=trim(file_pw), status='old', iostat=ios)
    if (ios /= 0) stop "ERROR opening file_pw! Stopping."

    unit_nscf = find_free_unit()
    open(unit_nscf, file=trim(file_nscf_in), status='unknown', iostat=ios)
    if (ios /= 0) stop "ERROR opening file_nscf_in! Stopping."

    do
      read(unit_pw,'(a)', iostat=ios) line
      if (ios /= 0) stop "ERROR reading file_pw! Stopping."
      lleft = trim(adjustl(line))
      if ( lleft(1:11) == 'calculation' ) then
        write(unit_nscf,*) " calculation = 'nscf'"
      else if ( lleft(1:6) == 'prefix' ) then
        write(unit_nscf,*) " prefix = '"//trim(prefix)//"-nscf'"
      else if ( lleft(1:8) == 'K_POINTS' ) then
        exit
      else
        write(unit_nscf,*) line
      end if
    end do

    write(unit_nscf,*) " K_POINTS {tpiba} "
    write(unit_nscf,*) nkpt_tr_tot
    do ikpt = 1, nkpt_tr_tot
      write(unit_nscf,'(3f16.9,2x,f4.1)') kpts_tr(:,ikpt), 1.0_dp
    end do

    close(unit_pw)
    close(unit_nscf)

    ! Create directory for nscf calculation prefix-nscf.save

    comando = 'mkdir -p '//trim(dir_nscf) ! it will not create it if it already exists
    call execute_command_line(comando)

    ! copy scf input files (force copy) to prefix-nscf.save

    ! If charge-density.dat and data-file-schema.xml have already been stored in
    ! prefix.save.intw directory, first try to get them from there.
    ! Otherwise, attempt to get them from prefix.save (the scf QE calculation directory,
    ! if it still exists)

    ! Try INTW directory (.save.intw)
    dir_scf_intw = trim(outdir)//trim(prefix)//".save.intw/"
    datafile = trim(dir_scf_intw)//'charge-density.dat'
    inquire(file=datafile, exist=exists_charge_density)
    datafile = trim(dir_scf_intw)//'data-file-schema.xml'
    inquire(file=datafile, exist=exists_schema)

    if ( exists_schema .and. exists_charge_density ) then
      ! Continuation files found in INTW directory
      dir_scf = dir_scf_intw
    else
      ! Try searching them in QE scf directory (prefix.save)
      write(*,'(A)') '| - Continuation files for nscf not present in:     |'
      write(*,'(A)') '|   '//dir_scf_intw(1:max(47,len(trim(dir_scf_intw))))//" |"

      dir_scf = trim(outdir)//trim(prefix)//".save/"
      write(*,'(A)') '| - Searching in:                                   |'
      write(*,'(A)') '|   '//dir_scf(1:max(47,len(trim(dir_scf))))//" |"

      datafile = trim(dir_scf)//'charge-density.dat'
      inquire(file=datafile, exist=exists_charge_density_qe)
      datafile = trim(dir_scf)//'data-file-schema.xml'
      inquire(file=datafile, exist=exists_schema_qe)
      if ( (.not. exists_schema_qe) .or. (.not. exists_charge_density_qe) ) then
        write(*,*) 'Continuation files for nscf neither present in: ', trim(dir_scf)
        write(*,*) 'Error. Stopping'
        stop
      end if
    end if

    write(*,'(A)') '| - Using continuation files for nscf from:         |'
    write(*,'(A)') '|   '//dir_scf(1:max(47,len(trim(dir_scf))))//" |"

    ! copy scf input files (force copy)

    comando = 'cp -f '//trim(dir_scf)//'charge-density.dat '//trim(dir_nscf)
    call execute_command_line(comando)
    comando = 'cp -f '//trim(dir_scf)//'data-file-schema.xml '//trim(dir_nscf)
    call execute_command_line(comando)

    ! Invoke QE pw.x to do nscf.
    ! command_pw may contain running options (e.g. mpirun)

    if ( command_pw == "unassigned" ) then
      stop "ERROR: Unassigned command_pw variable."
    else
      comando = trim(command_pw)//" < "//trim(file_nscf_in)//" > "//trim(file_nscf_out)
    end if
    write(*,'(A)') "| - Running QE with command:                        |"
    write(*,'(A)') "|   "//comando(1:max(47,len(trim(comando))))//" |"
    call execute_command_line(comando)

  end subroutine run_pw

  subroutine run_pw2intw(file_pw2intw_in, file_pw2intw_out)
    ! Create input file and run pw2intw.x

    use intw_input_parameters, only: command_pw2intw

    implicit none

    character(100), intent(in) :: file_pw2intw_in, file_pw2intw_out

    character(200) :: comando
    integer :: unit_pw2intw, ios


    ! Create input file for pw2intw

    unit_pw2intw = find_free_unit()
    open(unit_pw2intw, file=trim(file_pw2intw_in), status='unknown', iostat=ios)
    if (ios /= 0) stop "ERROR opening file_pw2intw_in! Stopping."

    write(unit_pw2intw, *) trim(prefix)
    write(unit_pw2intw, *) "&inputpp"
    write(unit_pw2intw, *) "  outdir = './'"
    write(unit_pw2intw, *) "  prefix = '"//trim(prefix)//"-nscf'"
    write(unit_pw2intw, *) "  phonons = .false."
    write(unit_pw2intw, *) "/"

    close(unit_pw2intw)

    ! Invoke pw2intw.x to transform nscf data to INTW format
    ! command_pw2intw may contain running options (e.g. mpirun)

    if ( command_pw2intw == "unassigned" ) then
      stop "ERROR: Unassigned command_pw2intw variable."
    else
      comando = trim(command_pw2intw)//" < "//trim(file_pw2intw_in)//" > "//trim(file_pw2intw_out)
    end if
    write(*,'(A)') "| - Running pw2intw with command:                   |"
    write(*,'(A)') "|   "//comando(1:max(47,len(trim(comando))))//" |"

    call execute_command_line(comando)

  end subroutine run_pw2intw


  subroutine SIESTA_nscf()

    implicit none

    character(100) :: file_s2intw_nscf_in, file_s2intw_nscf_out
    character(100) :: dir_nscf_intw
    logical :: have_nscf_out, have_nscf_wfc1, have_nscf_wfcN
    character(5) :: ik_loc

    !================================================================================
    ! Read s2intw.in file line by line and dump information to s2intw.in_nscf.
    ! nscf data will be written in prefix-nscf.save.intw directory.
    ! For that, the prefix is modified.
    ! The triangle list is added to the s2intw.in_nscf file.
    ! Note: this is done only in case the siesta2intw-nscf.out does not exist (otherwise,
    ! it is assumed that the wfcs have already been calculated and transformed to intw format)
    !================================================================================

    !================================================================================
    ! Check siesta2intw
    ! Check if first and last wfc are present. If so, it probably means that the nscf
    ! calculation was already made, so skip this step
    !================================================================================

    file_s2intw_nscf_in = trim(outdir)//trim(prefix)//'-nscf.siesta2intw.in'
    file_s2intw_nscf_out = trim(outdir)//trim(prefix)//'-nscf.siesta2intw.out'
    dir_nscf_intw = trim(outdir)//trim(prefix)//"-nscf.save.intw/"

    inquire(file=file_s2intw_nscf_out, exist=have_nscf_out)
    write(ik_loc,"(i5.5)") 1
    inquire(file=trim(dir_nscf_intw)//'wfc'//trim(adjustl(ik_loc))//'.dat', exist=have_nscf_wfc1)
    write(ik_loc,"(i5.5)") nkpt_tr_tot
    inquire(file=trim(dir_nscf_intw)//'wfc'//trim(adjustl(ik_loc))//'.dat', exist=have_nscf_wfcN)

    if (have_nscf_out .and. have_nscf_wfc1 .and. have_nscf_wfcN) then
      write(*,'(A)') '| - INTW nscf wfcs seem to be already available.    |'
      write(*,'(A)') '|   Skipping nscf calculation!                      |'
    else
      write(*,'(A)') '| - Calculating nscf wfcs with siesta2intw...       |'
      call run_siest2intw(file_s2intw_nscf_in, file_s2intw_nscf_out)
    endif

  end subroutine SIESTA_nscf

  subroutine run_siest2intw(file_s2intw_nscf_in, file_s2intw_nscf_out)
    ! Create input file and run siesta2intw.x

    use intw_input_parameters, only: command_siesta2intw, file_siesta2intw

    implicit none

    character(100), intent(in) :: file_s2intw_nscf_in, file_s2intw_nscf_out

    logical :: exists_density_matrix
    character(70) :: line, lleft
    logical :: write_nk
    character(200) :: comando, datafile
    integer :: unit_s2intw, unit_s2intw_nscf, ios
    integer :: ikpt
    real(dp), dimension(3) :: kpts_tr_cryst


    ! Create input file for siesta2intw

    unit_s2intw = find_free_unit()
    open(unit_s2intw, file=trim(file_siesta2intw), status='old', iostat=ios)
    if (ios /= 0) stop "ERROR opening file_siesta2intw! Stopping."

    unit_s2intw_nscf = find_free_unit()
    open(unit_s2intw_nscf, file=trim(file_s2intw_nscf_in), status='unknown', iostat=ios)
    if (ios /= 0) stop "ERROR opening file_s2intw_nscf_in! Stopping."

    write_nk = .true.
    do
      read(unit_s2intw,'(a)', iostat=ios) line
      if (ios /= 0) stop "ERROR reading file_siesta2intw! Stopping."
      lleft = trim(adjustl(line))
      if ( lleft(1:6) == 'prefix' ) then
        write(unit_s2intw_nscf,*) " prefix = '"//trim(prefix)//"-nscf'"
      else if ( lleft(1:7) == 'phonons' ) then
        write(unit_s2intw_nscf,*) " phonons = .false."
      else if ( lleft(1:3) == 'nk1' .or. lleft(1:3) == 'nk2' .or. lleft(1:3) == 'nk3' ) then
        if (write_nk) write(unit_s2intw_nscf,*) " nk1=0, nk2=0, nk3=0"
        write_nk = .false.
      else if ( lleft(1:1) == '/' ) then
        exit
      else
        write(unit_s2intw_nscf,*) trim(line)
      end if
    end do

    close(unit_s2intw)

    write(unit_s2intw_nscf,*) "/"
    write(unit_s2intw_nscf,*) "KPOINTS"
    write(unit_s2intw_nscf,*) nkpt_tr_tot
    do ikpt = 1, nkpt_tr_tot
      kpts_tr_cryst = kpts_tr(:,ikpt) ! this is cartesians x 2pi/alat. Transform to cryst.
      call cryst_to_cart(1, kpts_tr_cryst, at, -1)
      write(unit_s2intw_nscf,'(3f16.9,2x,f4.1)') kpts_tr_cryst
    end do

    close(unit_s2intw_nscf)

    ! Check if density matrix file is present

    datafile = trim(outdir)//trim(prefix)//".DM"
    inquire(file=datafile, exist=exists_density_matrix)

    if ( .not. exists_density_matrix ) then
      write(*,*) 'Density matrix file needed for nscf not present: ', trim(datafile)
      write(*,*) 'Error. Stopping'
      stop
    end if

    write(*,'(A)') '| - Using density matrix file:                      |'
    write(*,'(A)') '|   '//datafile(1:max(47,len(trim(datafile))))//" |"

    ! Invoke siesta2intw.x to do nscf and transform nscf data to INTW format
    ! command_siesta2intw may contain running options (e.g. mpirun)

    if ( command_siesta2intw == "unassigned" ) then
      stop "ERROR: Unassigned command_siesta2intw variable."
    else
      comando = trim(command_siesta2intw)//" < "//trim(file_s2intw_nscf_in)//" > "//trim(file_s2intw_nscf_out)
    end if
    write(*,'(A)') "| - Running siesta2intw with command:               |"
    write(*,'(A)') "|   "//comando(1:max(47,len(trim(comando))))//" |"
    call execute_command_line(comando)

  end subroutine run_siest2intw

end program ep_on_trFS_dV
