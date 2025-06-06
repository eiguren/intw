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
program a2F_on_trFS

  ! MBR 26/04/2024

  ! Calculate Eliashberg function from electron-phonon coupling
  ! matrix elements interpolated on a triangulated Fermi surface.

  ! a2F integral = double loop on kpoints of the FS triangulation.
  ! The ep element calculated for q=k'-k is read from a file.
  ! The needed phonons and dynamical matrices calculated by QE are
  ! read in and then interpolated as in the method of:
  ! F. Giustino et al, Phys. Rev. B 76, 165108 (2007)

  use kinds, only: dp

  use intw_useful_constants, only: cmplx_1, cmplx_0, cmplx_i, Ha_to_Ry, tpi, eps_8, eps_6

  use intw_utility, only: get_timing, find_free_unit, cryst_to_cart, &
                          smeared_delta, smeared_lorentz, &
                          generate_kmesh

  use intw_matrix_vector, only: ainv, area_vec

  use intw_input_parameters, only: mesh_dir, prefix, read_input, &
                                   nk1, nk2, nk3, nq1, nq2, nq3, nqirr, ph_dir, fc_mat, &
                                   nomega, omega_ini, omega_fin, osmear_q, &
                                   read_for_dynmat, &
                                   ep_interp_bands, nfs_sheets_initial, nfs_sheets_final

  use intw_reading, only: read_parameters_data_file_xml, &
                          nspin, noncolin, at, bg, volume0, alat, nat, ityp, amass, &
                          num_bands_intw, nsym, tau

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

  logical :: read_status
  logical :: q_points_consistent, full_mesh_q, IBZ_q
  character(5) :: is_loc, js_loc, comenta
  character(100) :: file_off, file_a2f
  integer :: unit_off, unit_a2f
  integer :: nkpt_tr_tot, nkpt_tr_ibz_tot
  integer :: qmesh_nqirr
  integer :: is, js, ik, ik1, iq, i, j, k, iface, ir1, ir2, ir3
  integer :: nfs_sheets_tot ! number of sheets considered
  integer, allocatable :: nfs_sheet(:), & ! band indices of the sheets (num_bands_intw set)
                          nkpt_tr(:), & ! number of kpoints in each FS sheet
                          nkpt_tr_ibz(:), & ! number of kpoints in each FS sheet irreducible BZ wedge
                          nface_tr(:), & ! number of faces in each FS sheet
                          nface_tr_ibz(:) ! number of faces in each IFS
  real(dp) :: k1(3), k2(3), k3(3), kwei, vol1bz
  real(dp), allocatable :: kpts_tr(:,:), kpts_tr_area(:), vk_tr(:,:), vabsk_tr(:), &
                           kpts_tr_ibz(:,:), kpts_tr_ibz_area(:), vk_tr_ibz(:,:), vabsk_tr_ibz(:)
  real(dp), allocatable :: area_ibz(:), area_fbz(:), factor_area_ibz(:)

  logical :: have_ep
  character(256) :: file_ep
  integer :: ib, iat, imode, ikp, unit_ep

  integer :: iomega, iks, ish, iksp, ishp, ibp
  real(dp) :: omega, omega_step, rfacq, dosef, dsk_vk_2
  real(dp) :: qpoint(3), kpoint(3), kpoint_p(3)
  real(dp), allocatable :: alpha2F(:,:,:,:), w2_qint(:), w_qint(:), lambda(:)
  complex(dp), allocatable :: gep_int(:)
  complex(dp), allocatable :: dyn_qint(:,:), u_qint(:,:)
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:)

  ! From mat_inv_four_t, see IGG's comments
  ! TODO add to useful_constants
  real(dp), parameter :: pmass = 1822.88848426_dp, aumev = 27211.396d0

  ! timing
  real(dp) :: time1, time2


  20 format(A)
  30 format(A,F8.2,6X,A)


  !================================================================================
  ! Beginning
  !================================================================================

  call get_timing(time1)

  write(*,20) '====================================================='
  write(*,20) '|                program a2F_on_triFS               |'
  write(*,20) '|         Eliashberg function calculation           |'
  write(*,20) '|    from interpolated ep elements on the FS        |'
  write(*,20) '|         ---------------------------------         |'
  write(*,20) '====================================================='


  !================================================================================
  ! Read the input file
  !================================================================================

  call read_input(read_status)

  if (read_status) stop


  !================================================================================
  ! Read the parameters from the SCF calculation
  !================================================================================

  write(*,20) '| - Reading calculation parameters...               |'

  call read_parameters_data_file_xml()


  !================================================================================
  ! Print spin information
  !================================================================================

  if (nspin==1) then
    write(*,20) '| - Paramagnetic calculation nspin=1                |'
  elseif (nspin==2) then
    write(*,20) '| - Spin-polarized calculation nspin = 2            |'
    if (noncolin) write(*,20) '| - Non-collinear spin calculation                  |'
  else
    write(*,20) '*****************************************************'
    write(*,20) '* ERROR: Allowed values for nspin are 1 or 2        *'
    write(*,20) '*            program stops.                         *'
    write(*,20) '*****************************************************'
    stop
  endif

  write(*,20) '|         ---------------------------------         |'


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

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read prefix.off and velocity file(s)
  !================================================================================

  write(*,20) '====================================================='
  write(*,20) '| - Reading .off and v_k files...                   |'

  allocate(nkpt_tr(nfs_sheets_tot), nface_tr(nfs_sheets_tot))
  allocate(nkpt_tr_ibz(nfs_sheets_tot), nface_tr_ibz(nfs_sheets_tot))

  ! open all sheet files just to see dimensions of kpoint lists

  do is=1,nfs_sheets_tot

    if (                is <  10) write(is_loc,"(i1)") nfs_sheet(is)
    if ( 10 <= is .and. is < 100) write(is_loc,"(i2)") nfs_sheet(is)

    file_off = trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_tri.off')
    write(*,'(A)') '|     '//file_off(1:max(45,len(trim(file_off))))//' |'

    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')
    read(unit_off,*) comenta
    read(unit_off,*) nkpt_tr(is), nface_tr(is), k ! number of vertices and faces (ignore edges)
    close(unit_off)

    ! open the IBZ off file and search for dimension nkpt_tr_ibz(is).
    ! Its vertices coincide with the first nkpt_tr_ibz(is) vertices of the full off vertex list.

    file_off = trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_IBZ_FS_tri.off')
    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')
    read(unit_off,*) comenta
    read(unit_off,*) nkpt_tr_ibz(is), nface_tr_ibz(is), k ! number of vertices and faces (ignore rest)
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
  write(*,20) '|   Number of faces in IBZ (total triangles):       |'
  write(*,'(A1,3X,I6,42X,A1)') '|', sum(nface_tr_ibz), '|'


  allocate(kpts_tr_ibz(3,nkpt_tr_ibz_tot), kpts_tr_ibz_area(nkpt_tr_ibz_tot))
  allocate(kpts_tr(3,nkpt_tr_tot), kpts_tr_area(nkpt_tr_tot))
  allocate(vk_tr_ibz(3,nkpt_tr_ibz_tot), vabsk_tr_ibz(nkpt_tr_ibz_tot))
  allocate(vk_tr(3,nkpt_tr_tot), vabsk_tr(nkpt_tr_tot))


  kpts_tr_ibz_area = 0.0_dp
  kpts_tr_area = 0.0_dp

  ! open .off files again to read k-points and read also velocity files
  ik1 = 0
  do is=1,nfs_sheets_tot

    if (                is <  10) write(is_loc,"(i1)") nfs_sheet(is)
    if ( 10 <= is .and. is < 100) write(is_loc,"(i2)") nfs_sheet(is)

    ! .off file for this sheet

    file_off = trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_tri.off')
    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')

    read(unit_off,*) comenta
    read(unit_off,*) i, j, k ! number vertices, faces and edges (I will ignore edges)
    ! read(unit_off,'(/)') ! DUDA... esto dependdera de como está escrito el salto de línea en el fichero, creo...
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
    do iface = 1, nface_tr(is)
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

    ! velocity for this sheet (use same unit)

    file_off = trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_tri_v_k.dat')
    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')

    read(unit_off,*) i
    if ( i /= nkpt_tr(is) ) then
      write(*,*) 'Error reading ', file_off, '. Stopping.'
      stop
    end if

    do ik=1,nkpt_tr(is)
      read(unit_off,*) i, vk_tr(:,ik1+ik), vabsk_tr(ik1+ik) ! velocity xyz and its modulus. DUDA 2pi/alat???
    end do

    close(unit_off)

    ! accumulate ik global index for the reading of next sheet
    ik1 = ik1 + nkpt_tr(is)

  end do


  ! MBR correction 100724
  ! now do the same with IBZ .off files to obtain the triangle areas there
  ik1 = 0
  do is=1,nfs_sheets_tot

    if (                is <  10) write(is_loc,"(i1)") nfs_sheet(is)
    if ( 10 <= is .and. is < 100) write(is_loc,"(i2)") nfs_sheet(is)

    ! .off file for this sheet

    file_off = trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_IBZ_FS_tri.off')
    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')

    read(unit_off,*) comenta
    read(unit_off,*) i, j, k ! number vertices, faces and edges (I will ignore edges)
    ! read(unit_off,'(/)') ! DUDA... esto dependdera de como está escrito el salto de línea en el fichero, creo...
    if ( (i /= nkpt_tr_ibz(is)) .or. (j /= nface_tr_ibz(is)) ) then
      write(*,*) 'Error reading ', file_off, '. Stopping.'
      stop
    end if

    ! read vertices (actually, these are repeated from the first nkpt_tr_ibz(is) lines of the full BZ file)
    do ik=1,nkpt_tr_ibz(is)
      read(unit_off,*) kpts_tr_ibz(:,ik1+ik) ! units in the trFS.off file are cartesian 2pi/alat ("tpiba" for QE)
    end do

    ! Read (triangular) faces on this sheet.
    ! Each face contributes with 1/3 of its area to the effective area of each of its vertices.
    ! Calculate the are on the go and add the contribution to each vertex, storing for global indices (i.e. ik1+ik).
    do iface = 1, nface_tr_ibz(is)
      read(unit_off,*) i, ir1, ir2, ir3 ! indices ik of the vertices of the face, indexed from 0
      ir1 = ir1 + 1
      ir2 = ir2 + 1
      ir3 = ir3 + 1 ! now, ik of the vertices of the face, indexed from 1
      if ( i /= 3 ) then
        write(*,*) 'Error reading ', file_off, 'Only triangles allowed. Stopping.'
        stop
      end if
      ! triangle vertex vectors (cartesian 2pi/alat)
      k1 = kpts_tr_ibz(:,ik1+ir1)
      k2 = kpts_tr_ibz(:,ik1+ir2)
      k3 = kpts_tr_ibz(:,ik1+ir3)
      ! get spanned area and add contribution to each vertex
      ! function copied from FSH/modules/geometry.f90
      kwei = area_vec(k2-k1,k3-k1)/3.0_dp
      kpts_tr_ibz_area(ik1+ir1) = kpts_tr_ibz_area(ik1+ir1) + kwei
      kpts_tr_ibz_area(ik1+ir2) = kpts_tr_ibz_area(ik1+ir2) + kwei
      kpts_tr_ibz_area(ik1+ir3) = kpts_tr_ibz_area(ik1+ir3) + kwei
    end do

    close(unit_off)

    ! velocity for this sheet (use same unit)

    file_off = trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_IBZ_FS_tri_v_k.dat')
    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')

    read(unit_off,*) i
    if ( i /= nkpt_tr_ibz(is) ) then
      write(*,*) 'Error reading ', file_off, '. Stopping.'
      stop
    end if

    do ik=1,nkpt_tr_ibz(is)
      read(unit_off,*) i, vk_tr_ibz(:,ik1+ik), vabsk_tr_ibz(ik1+ik) ! velocity xyz and its modulus (in units 2pi/alat, Hartree a.u.)
    end do

    close(unit_off)

    ! accumulate ik global index for the reading of next sheet
    ik1 = ik1 + nkpt_tr_ibz(is)

  end do

  write(*,20) '|   ...reading done                                 |'

  ! Calculate triangle areas for calculating the a2F integral later
  allocate(area_fbz(nfs_sheets_tot), area_ibz(nfs_sheets_tot), factor_area_ibz(nfs_sheets_tot))
  area_ibz = 0.0_dp
  ik = 0
  do ish = 1, nfs_sheets_tot
    do iks = 1, nkpt_tr_ibz(ish)
      ik = ik + 1
      area_ibz(ish) = area_ibz(ish) + kpts_tr_ibz_area(ik)
    end do
  end do

  area_fbz = 0.0_dp
  ik = 0
  do ish = 1, nfs_sheets_tot
    do iks = 1, nkpt_tr(ish)
      ik = ik + 1
      area_fbz(ish) = area_fbz(ish) + kpts_tr_area(ik)
    end do
  end do

  do ish=1,nfs_sheets_tot
    factor_area_ibz(ish) = area_fbz(ish)/area_ibz(ish)
  end do

  write(*,'(A36,F12.6,4X,A1)') '|   FS area in full BZ (tpiba^2) =   ', sum(area_fbz(:)), '|'
  do ish=1,nfs_sheets_tot
    write(*,'(A10,I3,F12.6,27X,A1)') '|      is ', ish, area_fbz(ish), '|'
  enddo

  write(*,'(A36,F12.6,4X,A1)') '|   FS area in irr  BZ (tpiba^2) =   ', sum(area_ibz(:)), '|'
  do ish=1,nfs_sheets_tot
    write(*,'(A10,I3,F12.6,27X,A1)') '|      is ', ish, area_ibz(ish), '|'
  enddo

  write(*,'(A36,F12.6,4X,A1)') '|   FS area ratio (full / irr BZ) = ', sum(factor_area_ibz(:)), '|'
  do ish=1,nfs_sheets_tot
    write(*,'(A10,I3,F12.6,27X,A1)') '|      is ', ish, factor_area_ibz(ish), '|'
  enddo

  ! N(EF) from sum over vertices, all sheets, using full triangulated mesh

  dosef = 0.0_dp
  do ik = 1, nkpt_tr_tot
    dosef = dosef + kpts_tr_area(ik) * (tpi/alat)**2 / ( vabsk_tr(ik) * Ha_to_Ry ) ! velocities in Hartree a.u., pass to Ry
  end do
  vol1bz = tpi**3 / volume0
  dosef = dosef / vol1bz ! in Ry^-1
  write(*,'(A29,F12.6,11X,A1)') '|   Volume of BZ (bohr^-3) = ', vol1bz, '|'
  write(*,'(A29,F12.6,11X,A1)') '|   DOS at FS (Ry^-1) =      ', dosef, '|'

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read all the information about phonons and prepare for interpolation of phonons:
  ! Calculate dynamical matrices in qmesh and transform to real space
  !================================================================================

  write(*,20) '| - Reading phonon info...                          |'

  ! Read q-points and irreducible patterns
  call read_ph_information_xml()

  allocate(q_irr_cryst(3,nqirr))
  do iq=1,nqirr
    q_irr_cryst(:,iq) = matmul(ainv(bg), q_irr(:,iq))
  enddo


  !================================================================================
  ! Build the phonon q-mesh
  !================================================================================

  write(*,20) '| - Building q-mesh...                              |'

  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  call generate_kmesh(qmesh, nq1, nq2, nq3)

  ! Find the size of the irreducible set of q-points (IBZ)
  call find_size_of_irreducible_k_set(nq1, nq2, nq3, qmesh, qmesh_nqirr)


  !================================================================================
  ! Set symmetry relations between irreducible q-points and full q-mesh
  !================================================================================

  allocate(QE_folder_nosym_q(nqmesh))
  allocate(QE_folder_sym_q(nqmesh))
  allocate(nosym_G_q(3,nqmesh))
  allocate(sym_G_q(3,nqmesh))
  allocate(symlink_q(nqmesh,2))

  call set_symmetry_relations(nq1, nq2, nq3, nqirr, q_irr_cryst, qmesh, q_points_consistent, &
                              QE_folder_nosym_q, QE_folder_sym_q, &
                              nosym_G_q, sym_G_q, symlink_q, full_mesh_q, IBZ_q)


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
    write(*,*) '*        nq_irr      = ', qmesh_nqirr
    stop
  end if


  !================================================================================
  ! Get dynamical matrices. Two options to get dyn_q:
  !  - Read dyn files.
  !  - Read the force constants and transfor to q space
  !================================================================================

  write(*,20) '| - Reading dynamical matrices...                   |'

  if (read_for_dynmat == 'fc') then ! read force constants
    call allocate_and_build_dyn_qmesh_from_fc(fc_mat)
  else if (read_for_dynmat == 'dynq') then ! read dyn files
    call allocate_and_build_dyn_qmesh()
  end if

  ! diagonalize
  do iq=1,nqmesh
    qpoint = qmesh(:,iq)
    call dyn_diagonalize_1q(3*nat, dyn_q(:,:,iq), u_q(:,:,iq), w2_q(:,iq))
  end do


  !================================================================================
  ! Wigner-Seitz cells
  !================================================================================

  write(*,20) '| - Building WS cells...                            |'

  ! This Wigner-Seitz mesh will be needed to interpolate phonons below
  ! (q->R space will be done on that grid of nrpts_q lattice points).
  call allocate_and_build_ws_irvec_qtau()


  !================================================================================
  ! Transform dynamical matrices to real space (force constants)
  !================================================================================

  write(*,20) '| - Computing force constants...                    |'

  call dyn_q_to_dyn_r()

  ! test decay of dyn_r elements with distance
  ! do ir=1,nrpts_q
  !   rcart = real(irvec_q(:,ir),dp)
  !   call cryst_to_cart(1, rcart, at, 1)
  !   rcart = rcart * alat ! bohr units
  !   write(520,'(i5,f16.6,8e16.4)') ir, sqrt ( sum(rcart*rcart) ), &
  !           abs(dyn_r(1,1,ir)), abs(dyn_r(1,2,ir)), abs(dyn_r(1,4,ir)), abs(dyn_r(1,5,ir))
  ! end do

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read ep elements file
  !================================================================================

  write(*,20) '| - Reading ep_mat files...                         |'

  file_ep = trim(mesh_dir)//trim(prefix)//trim('_ep_interp.dat')
  inquire(file=file_ep, exist=have_ep)

  if (.not.have_ep) then

    write(*,*) '!      e-p elements interpolated and written to file:         !'
    write(*,*) trim(file_ep)
    write(*,*) '!      e-p elements file does not exist:                      !'
    write(*,*) '!      Stopping                                               !'
    stop

  else

    ! Note, indices ikp,ik include all the calculated sheets
    allocate(aep_mat_el(nkpt_tr_tot,nkpt_tr_ibz_tot,nspin,nspin,3*nat))

    unit_ep = find_free_unit()
    open(unit_ep, file=file_ep, status='old')
    read(unit_ep,*) comenta

    do i = 1, nkpt_tr_tot
      do j = 1, nkpt_tr_ibz_tot
        do js=1,nspin
          do is=1,nspin
            read(unit_ep, fmt="(6i6,100e16.6)") ibp, iksp, ikp, ib, iks, ik,&
                  (aep_mat_el(ikp,ik, js,is,iat), iat=1,3*nat)
          end do
        end do
      end do
    end do

    close(unit_ep)

  end if

  write(*,20) '|   ...reading done                                 |'
  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Calculate a2F integral
  !================================================================================

  ! phonon arrays
  allocate(dyn_qint(3*nat,3*nat), u_qint(3*nat,3*nat), w2_qint(3*nat), w_qint(3*nat))
  ! ep elements
  allocate(gep_int(3*nat))
  ! integrals
  allocate(alpha2F(nspin,nspin,3*nat,nomega), lambda(3*nat))


  alpha2F = 0.0_dp
  lambda = 0.0_dp

  ! Energy range for Eliashberg is read from intw.in
  ! with Ry units throughout.
  ! Step in energy for integrals and phonon DOS:
  omega_step = (omega_fin-omega_ini)/real(nomega-1,dp)

  write(*,20) '| - Calculating a2F...                              |'

  ! ik, ikp indices implicitly contain the FS sheet index, i.e. the band indices ib, ib'
  ! to be selected, so instead of iterating over nkpt_tr_tot, I separate over sheets

  ik = 0
  ik1 = 0
  do ish = 1, nfs_sheets_tot

    ib = nfs_sheet(ish) ! band index for k

    do iks = 1, nkpt_tr_ibz(ish)

      ! ik is the k-index over nkpt_tr_tot in the Irreducible BZ
      ik = ik + 1

      kpoint = kpts_tr_ibz(:,ik) ! this is cartesians x 2pi/alat. Transform to cryst.
      call cryst_to_cart(1, kpoint, at, -1)

      ikp = 0
      do ishp = 1, nfs_sheets_tot

        ibp = nfs_sheet(ishp) ! band index for k'

        do iksp = 1, nkpt_tr(ishp)

          ikp = ikp + 1 ! k'-index over nkpt_tr_tot

          kpoint_p = kpts_tr(:,ikp) ! this is cartesians x 2pi/alat. Transform to cryst.
          call cryst_to_cart(1, kpoint_p, at, -1)

          qpoint = kpoint_p-kpoint

          ! interpolate phonon:
          call dyn_interp_1q(qpoint, dyn_qint)
          call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint)
          w_qint = sign(sqrt(abs(w2_qint)), w2_qint)
          ! u_qint contains the polarization vector (each column is a mode)
          ! phonon frequency in a.u.: pass to Ry
          w_qint = w_qint * Ha_to_Ry

          ! Generate ep elements with the added mass factor and polarization vector
          ! component for each mode and sum over modes.
          ! In the eigenvector matrix u_qint, each column is a mode. Rows are individual atom displacements.
          ! The matrix elements g are given for atom displacements and I want them as a
          ! function of modes to compare with QE a2F.dos* files

          ! band indices to be used are the FS sheet indices given at the beginning of the loop.
          ! I select only those ones for the interpolated elements

          do is=1,nspin
            do js=1,nspin
              !
              gep_int = cmplx_0
              !
              do imode = 1, 3*nat
                !
                if (w_qint(imode)>eps_8) then ! stable mode frequency
                  do iat=1,3*nat
                    k = (iat-1)/3 + 1 ! atom index
                    rfacq = 2.0_dp * w_qint(imode) * amass(ityp(k)) * (pmass / Ha_to_Ry) ! pmass is in Ha a.u., so pass to Ry
                    gep_int(imode) = gep_int(imode) + &
                        aep_mat_el(ikp,ik,is,js,iat) * u_qint(iat,imode) / sqrt(rfacq)
                  end do
                end if
                !
                ! for testing:
                write(999, fmt="(5i6,2e16.6)") ikp, ik, ibp, ib, imode, w_qint(imode), abs(gep_int(imode))
                !
              end do

              ! Sum over modes and weight with areas and velocities v_k * v_k' for contribution of k',k pair to a2F.
              ! Mind, the areas were in (tpi/alat)**2 units.
              ! Add weight of this irreducible wedge (IBZ) to the integral using factor_area_ibz.

              dsk_vk_2 = kpts_tr_ibz_area(ik) * factor_area_ibz(ish) * kpts_tr_area(ikp) * &
                        (tpi/alat)**4 / ( vabsk_tr_ibz(ik) * vabsk_tr(ikp) * Ha_to_Ry * Ha_to_Ry ) ! velocities in Hartree a.u., pass to Ry

              do iomega = 1, nomega ! frequencies
                omega = omega_ini + omega_step*real(iomega-1,dp)
                do imode = 1, 3*nat

                  ! smear with a gaussian
                  rfacq = smeared_delta(omega-w_qint(imode), osmear_q) ! smear omega(q)
                  rfacq = rfacq - smeared_delta(omega+w_qint(imode), osmear_q) ! add antisymmetric term at -wq to remove divergence at w=0

                  ! or with a lorentzian
                  ! rfacq = smeared_lorentz(omega-w_qint(imode), osmear_q) ! smear omega(q)
                  ! rfacq = rfacq - smeared_lorentz(omega+w_qint(imode), osmear_q) ! add antisymmetric term at -wq to remove divergence at w=0

                  alpha2F(is,js,imode,iomega) = alpha2F(is,js,imode,iomega) + &
                          (abs(gep_int(imode)))**2 * rfacq * dsk_vk_2

                end do ! imode
              end do ! iomega

            end do ! spin js
          end do ! spin is

        end do ! k' (kpoint)
      end do ! k' (sheet)

    end do ! k (kpoint)
  end do ! k (sheet)

  ! divide by the N(E_F) factor and by 1BZ-volume^2 (because of k and k' integrals)
  alpha2F = alpha2F/dosef/vol1bz**2

  ! print out by spin blocks (one file per block if nspin=2).
  do is=1,nspin
    do js=1,nspin

      if (nspin == 1) then
        file_a2f = trim(mesh_dir)//trim(prefix)//trim('_a2F_interp.dat')
        unit_a2f = find_free_unit()
        open(unit_a2f, file=file_a2f, status='unknown')
        write(unit_a2f,'(A)') '#omega(Ry)  alpha2F(total)    alpha2F(1:nmode)'
      else
        write(is_loc,"(i1)") is
        write(js_loc,"(i1)") js
        file_a2f = trim(mesh_dir)//trim(prefix)//trim('_s')//trim(adjustl(is_loc))//trim(adjustl(js_loc))//trim('_a2F_interp.dat')
        unit_a2f = find_free_unit()
        open(unit_a2f, file=file_a2f, status='unknown')
        write(unit_a2f,'(A)') '#omega(Ry)  alpha2F(total)    alpha2F(1:nmode)'
        write(unit_a2f,*) '#spins = ', is,js
      end if

      ! print a2F(omega)
      omega = omega_ini
      do iomega=1,nomega ! frequencies
        write(unit_a2f,'(e16.4,2x,15e16.4)') omega, sum(alpha2F(is,js,:,iomega)), &
                                                    alpha2F(is,js,:,iomega)
        omega = omega + omega_step
      end do

      ! calculate lambda for spin diagonal

      if (is == js) then
        lambda = 0.0_dp
        do imode = 1, 3*nat
          omega = omega_ini + omega_step
          do iomega=2,nomega ! skip omega=0 value
            lambda(imode) = lambda(imode) + alpha2F(is,js,imode,iomega) / omega
            omega = omega + omega_step
          end do
          lambda(imode) = lambda(imode) * 2.0_dp * omega_step
        end do
        write(unit_a2f,*) '# lambda = 2 \int dw a2F(w)/w for each mode'
        write(unit_a2f,'(a,i3,a,20e14.4)') '#imode = ', imode, ', lambda = ', lambda(:)
        write(unit_a2f,*) '#total lambda = ', sum(lambda)
      end if

      close(unit_a2f)

      write(*,20) '| - a2F calculated and written to file:             |'
      write(*,20) "|   "//file_a2f(1:max(47,len(trim(file_a2f))))//" |"

    end do ! js
  end do ! is

  write(*,20) '====================================================='


  !================================================================================
  ! Finish
  !================================================================================

  call get_timing(time2)

  write(*,20) '|                      ALL DONE                     |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='

end program a2F_on_trFS
