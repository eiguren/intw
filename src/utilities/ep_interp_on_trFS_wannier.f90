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
program ep_on_trFS_wannier

  ! MBR 24/04/2024

  ! This utility reads the electron-phonon matrix elements calculated by the utility
  ! ep_melements.f90 of INTW and interpolates them on a triangulated Fermi surface,
  ! following the method of:
  ! F. Giustino et al, Phys. Rev. B 76, 165108 (2007)
  ! Finally, elements are printed to file.

#ifdef _OPENMP
  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
#endif

  use kinds, only: dp

  use intw_version, only: print_intw_version

  use intw_useful_constants, only: cmplx_0, cmplx_i, Ha_to_eV, tpi

  use intw_utility, only: get_timing, print_threads, print_date_time, find_free_unit, &
                          find_k_1BZ_and_G, triple_to_joint_index_g, &
                          generate_kmesh, cryst_to_cart

  use intw_matrix_vector, only: area_vec

  use intw_w90_setup, only: nrpts, &
                            allocate_and_read_ham_r, allocate_and_read_u_mesh, &
                            wann_rotate_matrix, wann_IFT_1index, wann_FT_1index_1k, &
                            interpolate_1k

  use intw_input_parameters, only: nk1, nk2, nk3, outdir, prefix, read_input, &
                                   nq1, nq2, nq3, &
                                   ep_bands, ep_bands_initial, ep_bands_final, &
                                   use_exclude_bands, &
                                   ep_interp_method, ep_interp_bands, nfs_sheets_initial, nfs_sheets_final, &
                                   ep_mat_file

  use intw_reading, only: read_parameters_data_file, set_num_bands, &
                          nat, nspin, lspin, num_bands_intw, num_wann_intw, at

  use intw_intw2wannier, only: nnkp_kpoints

  use intw_ph, only: nqmesh, qmesh, read_ph_information

  use intw_ph_interpolate, only: allocate_and_build_ws_irvec_q, &
                                 irvec_q, nrpts_q, ndegen_q, &
                                 wann_IFT_1index_q

  implicit none

  ! FS triangulation
  integer :: nfs_sheets_tot ! number of FS sheets considered
  integer :: nkpt_tr_tot, nkpt_tr_ibz_tot ! total number of kpoints in the FS and in the irreducible BZ wedge
  integer, allocatable :: nfs_sheet(:), & ! band indices of the FS sheets (from num_bands_intw set)
                          nkpt_tr(:), & ! number of kpoints in each FS sheet
                          nkpt_tr_ibz(:), & ! number of kpoints in each FS sheet irreducible BZ wedge
                          nface_tr(:) ! number of faces in each FS sheet
  real(dp), allocatable :: kpts_tr(:,:), & ! list of all kpoints
                           kpts_tr_area(:) ! area of each kpoint in the FS
  integer, allocatable :: ikibz_2_ik(:), & ! index of the ikibz kpoint inside the irreducible BZ wedge in the list of all kpoints
                          ik_2_ish(:), & ! FS sheet index of kpts_tr(:,ik)
                          ik_2_iks(:) ! index of kpts_tr(:,ik) in its FS sheet kpoints list

  ! Part I
  logical :: read_status
  character(5) :: ish_loc, comment
  character(100) :: file_off
  integer :: unit_off
  real(dp) :: k1(3), k2(3), k3(3), kwei

  ! Part II
  character(4) :: iq_loc
  integer :: record_lengh, ierr, ep_unit
  integer :: ikq, irq
  integer :: nkmesh, Gkq_1bz(3)
  integer, allocatable :: kqmap(:,:)
  real(dp) :: kpoint(3), qpoint(3), kqpoint(3), kq_1bz(3)
  real(dp), allocatable :: kmesh(:,:), kqmesh(:,:,:)
  real(dp), allocatable :: eig_kint(:), eig_kqint(:)
  complex(dp) :: facq
  complex(dp), allocatable :: u_kint(:,:), u_kqint(:,:), u_kint_all(:,:,:)
  complex(dp), allocatable :: gmatkqk_wann(:,:,:,:), gmatL_wann(:,:,:,:,:,:,:), gmat_aux(:,:,:,:), &
                              gmat_int(:,:), gmat_int_rot(:,:), gmat_aux1(:,:,:), gep_int(:,:)
  complex(dp), allocatable :: ep_mat_el_coarse(:,:,:,:,:,:,:)
  integer :: i_start, i_end
#ifdef _OPENMP
  integer :: n, n_remaining
  integer :: thread_id, thread_num
#endif

  ! Part III
  logical :: have_ep
  character(256) :: file_ep
  integer :: unit_ep
  real(dp) :: kpoint_p(3)
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:)

  ! loop variables and indices
  integer :: ik, ikp, ikibz, ikibz_global, ikibz_do, iks, iksp
  integer :: ish, ishp, ib, ibp
  integer :: ik1, ik2, ik3
  integer :: iface, iedge
  integer :: is, js
  integer :: ir
  integer :: iat
  integer :: iq

  ! timing
  real(dp) :: time1, time2


  20 format(A)
  30 format(A,F8.2,6X,A)


  !================================================================================
  ! Beginning
  !================================================================================

  call get_timing(time1)

  write(*,20) '====================================================='
  write(*,20) '|            program ep_on_triFS_wannier            |'
  write(*,20) '|       Wannier interpolation of e-p elements       |'
  write(*,20) '|               on the triangulated FS              |'
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

  if ( ep_interp_method /= 'wannier' ) then
    write(*,*) 'ep_interp_method /= wannier in input. Stopping.'
    stop
  end if

  if ( use_exclude_bands /= "wannier" ) then
    write(*,*) 'use_exclude_bands /= wannier in input. Stopping.'
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
  ! Control: for the Wannier interpolation we need all coarse ep elements in the
  ! ep_bands == 'intw' set (otherwise we cannot rotate with the U matrices).
  ! Therefore, check if we will have all the elements, even though custom_bands
  ! was set.
  !================================================================================

  if ( ep_bands == 'custom' .and. &
      ( (ep_bands_initial /= 1) .or. (ep_bands_final /= num_bands_intw) ) ) then
    write(*,*) 'Coarse ep elements where calculated in the', &
                ep_bands_initial, ' to', ep_bands_final, ' band subset only.'
    write(*,*) 'This means that we cannot Wannier interpolate. Stopping.'
    stop
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
  do ish=1,nfs_sheets_tot

    if (                 ish <  10) write(ish_loc,"(i1)") nfs_sheet(ish)
    if ( 10 <= ish .and. ish < 100) write(ish_loc,"(i2)") nfs_sheet(ish)

    file_off = trim(outdir)//trim(prefix)//trim('.')//trim(adjustl(ish_loc))//trim('_FS_tri.off')
    write(*,'(A)') '|     '//file_off(1:max(45,len(trim(file_off))))//' |'

    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')
    read(unit_off,*) comment
    read(unit_off,*) nkpt_tr(ish), nface_tr(ish), iedge ! number of vertices and faces (ignore edges)
    close(unit_off)

    ! open the IBZ off file and search for dimension nkpt_tr_ibz(ish).
    ! Its vertices coincide with the first nkpt_tr_ibz(ish) vertices of the full off vertex list.

    file_off = trim(outdir)//trim(prefix)//trim('.')//trim(adjustl(ish_loc))//trim('_IBZ_FS_tri.off')

    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')
    read(unit_off,*) comment
    read(unit_off,*) nkpt_tr_ibz(ish), iface, iedge ! number of vertices (ignore faces and edges)
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
  allocate(ikibz_2_ik(nkpt_tr_ibz_tot))
  allocate(ik_2_ish(nkpt_tr_tot))
  allocate(ik_2_iks(nkpt_tr_tot))
  ikibz_2_ik = -10
  ik_2_ish = -10
  ik_2_iks = -10


  ! open .off files again to read k-points
  do ish=1,nfs_sheets_tot

    if (                 ish <  10) write(ish_loc,"(i1)") nfs_sheet(ish)
    if ( 10 <= ish .and. ish < 100) write(ish_loc,"(i2)") nfs_sheet(ish)

    ! .off file for this sheet

    file_off = trim(outdir)//trim(prefix)//trim('.')//trim(adjustl(ish_loc))//trim('_FS_tri.off')
    unit_off = find_free_unit()
    open(unit_off, file=file_off, status='old')

    read(unit_off,*) comment
    read(unit_off,*) ik, iface, iedge ! number vertices, faces and edges (I will ignore edges)
    ! read(unit_off,'(/)') ! DUDA... This will depend on how the line break is written in the file, I think...
    if ( (ik /= nkpt_tr(ish)) .or. (iface /= nface_tr(ish)) ) then
      write(*,*) 'Error reading ', file_off, '. Stopping.'
      stop
    end if

    ! read vertices
    do iks=1,nkpt_tr(ish)
      ik = iks + sum(nkpt_tr(:ish-1))
      ik_2_ish(ik) = ish
      ik_2_iks(ik) = iks
      if (iks<=nkpt_tr_ibz(ish)) then
        ! kpoint is in the irreducible BZ wedge
        ikibz = iks + sum(nkpt_tr_ibz(:ish-1))
        ikibz_2_ik(ikibz) = ik
      endif
      read(unit_off,*) kpts_tr(:,ik) ! units in the trFS.off file are cartesian 2pi/alat ("tpiba" for QE)
    end do

    ! Read (triangular) faces on this sheet.
    ! Each face contributes with 1/3 of its area to the effective area of each of its vertices.
    ! Calculate the are on the go and add the contribution to each vertex, storing for global indices (i.e. ik).
    do iface = 1, nface_tr(ish)
      read(unit_off,*) ik, ik1, ik2, ik3 ! indices ik of the vertices of the face, indexed from 0
      if ( ik /= 3 ) then
        write(*,*) 'Error reading ', file_off, 'Only triangles allowed. Stopping.'
        stop
      end if
      ik1 = ik1 + 1
      ik2 = ik2 + 1
      ik3 = ik3 + 1 ! now, ik of the vertices of the face, indexed from 1
      ik1 = ik1 + sum(nkpt_tr(:ish-1))
      ik2 = ik2 + sum(nkpt_tr(:ish-1))
      ik3 = ik3 + sum(nkpt_tr(:ish-1)) ! now, ik in the global ik list
      ! triangle vertex vectors (cartesian 2pi/alat)
      k1 = kpts_tr(:,ik1)
      k2 = kpts_tr(:,ik2)
      k3 = kpts_tr(:,ik3)
      ! get spanned area and add contribution to each vertex
      ! function copied from FSH/modules/geometry.f90
      kwei = area_vec(k2-k1,k3-k1)/3.0_dp
      kpts_tr_area(ik1) = kpts_tr_area(ik1) + kwei
      kpts_tr_area(ik2) = kpts_tr_area(ik2) + kwei
      kpts_tr_area(ik3) = kpts_tr_area(ik3) + kwei
    end do

    close(unit_off)

  end do

  if (any(ikibz_2_ik == -10)) stop "ERROR ikibz_2_ik"
  if (any(ik_2_ish == -10)) stop "ERROR ik_2_ib"
  if (any(ik_2_iks == -10)) stop "ERROR ik_2_iks"

  write(*,20) '|   .... reading done                               |'

  write(*,20) '|   Total FS area:                                  |'
  write(*,'(A1,3X,F12.6,A19,17X,A1)') '|', sum(kpts_tr_area), ' (2 x pi / alat)^2 ', '|'

  write(*,20) '|                                                   |'
  write(*,20) '| ---------------- Part I completed --------------- |'
  write(*,20) '|                                                   |'
  write(*,20) '====================================================='


  !******************************** Part II ************************************

  !================================================================================
  ! Read u_mesh file from w902intw
  !================================================================================

  write(*,20) '| - Reading Wannier U matrix...                     |'

  call allocate_and_read_u_mesh()


  !================================================================================
  ! Read ham_r file from w902intw, to be used in band interpolation
  ! (theis also allocates and reads the irvec list)
  !================================================================================

  write(*,20) '| - Reading Wannier H(R)...                         |'

  call allocate_and_read_ham_r()

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read phonon information
  !================================================================================

  write(*,20) '| - Reading phonon info...                          |'

  ! Read irreducible q-points and irreducible patterns
  call read_ph_information()


  !================================================================================
  ! Generate coarse meshes
  !================================================================================

  write(*,20) '| - Building coarse k-mesh...                       |'

  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh, nk1, nk2, nk3)

  write(*,20) '| - Building coarse q-mesh...                       |'

  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  call generate_kmesh(qmesh, nq1, nq2, nq3)


  ! DUDA kmesh is equal to nnkp_kpoints, but I do not know if this is general
  ! if not, we have to use kmesh explicitly in the (I)FTs
  ! do ik=1,nkmesh
  !   write(*,'(i3,3f7.3)') ik, kmesh(:,ik)-nnkp_kpoints(:,ik)
  ! end do

  ! k+q mesh and index correspondence
  allocate(kqmesh(3,nqmesh,nkmesh))
  allocate(kqmap(nqmesh,nkmesh))
  do iq=1,nqmesh
    qpoint = qmesh(:,iq)
    do ik=1,nkmesh
      kpoint = kmesh(:,ik) ! nnkp_kpoints(:,ik) is equivalent
      kqpoint = kpoint+qpoint
      kqmesh(:,iq,ik) = kqpoint
      ! elements k+q,k are: ep_mat_el_coarse(iq,ik,:,:,is,js,imode)
      ! locate (k+q)-point index ikq in the kmesh
      call find_k_1BZ_and_G(kqpoint, nk1, nk2, nk3, ik1, ik2, ik3, kq_1bz, Gkq_1bz)
      call triple_to_joint_index_g(nk1, nk2, nk3, ikq, ik1, ik2, ik3)
      kqmap(iq,ik) = ikq
    end do
  end do

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read ep mat elements
  !================================================================================

  write(*,20) '| - Reading ep_mat files...                         |'

  allocate(ep_mat_el_coarse(nqmesh,nkmesh,num_bands_intw,num_bands_intw,nspin,nspin,3*nat))
  inquire(iolength=record_lengh) ep_mat_el_coarse(1,:,:,:,:,:,:)

  do iq = 1, nqmesh
    ep_unit = find_free_unit()
    if (                iq <   10) write(iq_loc,"(i1)")iq
    if ( 10 <= iq .and. iq <  100) write(iq_loc,"(i2)")iq
    if (100 <= iq .and. iq < 1000) write(iq_loc,"(i3)")iq
    open(unit=ep_unit, file=trim(ep_mat_file)//trim('_')//adjustl(iq_loc), &
         status='old', form='unformatted', access='direct', recl=record_lengh)
    read(unit=ep_unit, rec=1, iostat=ierr) ep_mat_el_coarse(iq,:,:,:,:,:,:)
    close(unit=ep_unit)
  end do
  write(*,20) '|                       ...done                     |'
  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Get the WS vectors for the coarse q-mesh Wannier interpolation
  !================================================================================

  write(*,20) '| - Building WS mesh...                             |'

  call allocate_and_build_ws_irvec_q()


  !================================================================================
  ! ep elements in real space
  ! For each nspin,nspin,3*nat:
  !    1. rotate matrix g(k+q,k,m,n) to wannier gauge with u_mesh
  !    2. For k+q and k indices, inverse Fourier index by index.
  !       Note, k and q are fouriered on different grids.
  !================================================================================

  write(*,20) '| - Fourier transform ep mat elements...            |'

  write(*,'(A42,I4,A7)') '|   Number of WS vectors for k transform: ', nrpts, '      |'
  write(*,'(A42,I4,A7)') '|   Number of WS vectors for q transform: ', nrpts_q, '      |'

  allocate(gmatkqk_wann(num_wann_intw,num_wann_intw,nqmesh,nkmesh))
  allocate(gmat_aux(num_wann_intw,num_wann_intw,nqmesh,nrpts))
  allocate(gmatL_wann(3*nat,num_wann_intw,num_wann_intw,nspin,nspin,nrpts_q,nrpts))

  gmatL_wann = cmplx_0

  do iat=1,3*nat

    do is=1,nspin
      do js=1,nspin

        ! 1.
        ! rotate U^dagger(k+q) * gmat * U(k)
        gmatkqk_wann = cmplx_0
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(is, js, iat) &
        !$omp shared(nqmesh, nkmesh, kqmap) &
        !$omp shared(ep_mat_el_coarse, gmatkqk_wann) &
        !$omp private(ik, ikq)
        do iq=1,nqmesh
          !$omp parallel do &
          !$omp default(none) &
          !$omp shared(is, js, iat) &
          !$omp shared(iq, nkmesh, kqmap) &
          !$omp shared(ep_mat_el_coarse, gmatkqk_wann) &
          !$omp private(ikq)
          do ik=1,nkmesh
            ikq = kqmap(iq,ik)
            call wann_rotate_matrix(ikq, ik, ep_mat_el_coarse(iq,ik,:,:,is,js,iat), gmatkqk_wann(:,:,iq,ik))
          end do
          !$omp end parallel do
        end do
        !$omp end parallel do

        ! 2.
        ! Inverse Fourier over k-index, using kmesh and irvec grids
        gmat_aux = cmplx_0
        !$omp parallel &
        !$omp default(none) &
        !$omp shared(nqmesh, nkmesh, kmesh) &
        !$omp shared(gmatkqk_wann, gmat_aux) &
        !$omp shared(thread_num, n, n_remaining) &
        !$omp private(thread_id, i_start, i_end)
        !
        ! Calculate the range of iterations for this thread.
        ! If nqmesh is a multiple of thread_num, each thread
        ! will run n iterations.
        ! Otherwise, the first n_remaining threads will run
        ! an extra iteration.
        !
#ifdef _OPENMP
        !$omp single
        thread_num = omp_get_num_threads()
        n = int(nqmesh/thread_num) ! Number of iterations for each thread
        n_remaining = mod(nqmesh, thread_num) ! Remainig iterations that need to be distributed
        !$omp end single
        !
        thread_id = omp_get_thread_num()
        i_start = n * thread_id + min(thread_id, n_remaining) + 1
        i_end = n * (thread_id + 1) + min(thread_id + 1, n_remaining)
#else
        i_start = 1
        i_end = nqmesh
#endif
        !
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(nqmesh, nkmesh, kmesh) &
        !$omp shared(gmatkqk_wann, gmat_aux) &
        !$omp shared(i_start, i_end)
        do iq=i_start,i_end
          call wann_IFT_1index(nkmesh, kmesh, gmatkqk_wann(:,:,iq,:), gmat_aux(:,:,iq,:))
        end do
        !$omp end parallel do
        !$omp end parallel

        ! 3.
        ! Inverse Fourier over q-index, using qmesh and irvec_q grids
        !$omp parallel &
        !$omp default(none) &
        !$omp shared(iat, is, js) &
        !$omp shared(nrpts, nqmesh, qmesh) &
        !$omp shared(gmat_aux, gmatL_wann) &
        !$omp shared(thread_num, n, n_remaining) &
        !$omp private(thread_id, i_start, i_end)
        !
        ! Calculate the range of iterations for this thread.
        ! If nrpts is a multiple of thread_num, each thread
        ! will run n iterations.
        ! Otherwise, the first n_remaining threads will run
        ! an extra iteration.
        !
#ifdef _OPENMP
        !$omp single
        thread_num = omp_get_num_threads()
        n = int(nrpts/thread_num) ! Number of iterations for each thread
        n_remaining = mod(nrpts, thread_num) ! Remainig iterations that need to be distributed
        !$omp end single
        !
        thread_id = omp_get_thread_num()
        i_start = n * thread_id + min(thread_id, n_remaining) + 1
        i_end = n * (thread_id + 1) + min(thread_id + 1, n_remaining)
#else
        i_start = 1
        i_end = nrpts
#endif
        !
        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(iat, is, js) &
        !$omp shared(nrpts, nqmesh, qmesh) &
        !$omp shared(gmat_aux, gmatL_wann) &
        !$omp shared(i_start, i_end)
        do ir=i_start,i_end
          call wann_IFT_1index_q(nqmesh, qmesh, gmat_aux(:,:,:,ir), gmatL_wann(iat,:,:,is,js,:,ir))
        end do
        !$omp end parallel do
        !$omp end parallel

      end do ! spin
    end do ! spin

    write(*,'(A38,I4,10X,A1)') '|     IFT on WS done for displacement ', iat, '|'

  end do ! iat

  deallocate(ep_mat_el_coarse)
  deallocate(gmat_aux, gmatkqk_wann)

  write(*,20) '|                                                   |'
  write(*,20) '| --------------- Part II completed --------------- |'
  write(*,20) '|                                                   |'
  write(*,20) '====================================================='


  !******************************* Part III *************************************

  !================================================================================
  ! Interpolate bands
  !================================================================================

  write(*,20) '| - Interpolating Wannier U matrix...               |'

  ! interpolated bands in nkpt_tr_tot vertices
  ! u_kint_all, contains the rotation matrices (row is bloch and column is wannier,
  ! i.e. the "orbital weights" for each band is in the columns)

  ! band arrays
  allocate(eig_kint(num_wann_intw), eig_kqint(num_wann_intw))
  allocate(u_kint(num_wann_intw,num_wann_intw), u_kqint(num_wann_intw,num_wann_intw))
  allocate(u_kint_all(num_wann_intw,num_wann_intw,nkpt_tr_tot))

  do ik=1,nkpt_tr_tot
    ! this is cartesians x 2pi/alat. Transform to crystal before interpolation
    kpoint = kpts_tr(:,ik)
    call cryst_to_cart(1, kpoint, at, -1)
    call interpolate_1k(kpoint, eig_kint, u_kint)
    ! eig_kint = eig_kint * 2.0_dp / Ha_to_eV
    u_kint_all(:,:,ik) = u_kint
  end do

  deallocate(eig_kint)


  !================================================================================
  ! Interpolate elements
  ! Loop over k'=k+q, k electrons. Get interpolated ep elements and write out
  !================================================================================

  write(*,20) '| - Interpolating e-p elements...                   |'

  ! the interpolation grid is formed by the triangulated FS points
  ! In principle, in the trianglulated k-point list up to nkpt_tr_tot
  ! the index carries the information about the band implicitly.
  ! We will retrieve U rotation matrices and eigenvalues on the
  ! full num_wann_intw list for those k-points, as the w90_setup
  ! routines give that info, but then we will intepolate the ep element only
  ! for the bands index at Fermi given by nfs_sheet().

  ! ep arrays
  allocate(gmat_aux1(num_wann_intw,num_wann_intw,nrpts))
  allocate(gmat_int(num_wann_intw,num_wann_intw))
  allocate(gmat_int_rot(num_wann_intw,num_wann_intw))
  allocate(gep_int(num_wann_intw,num_wann_intw))
  allocate(aep_mat_el(nkpt_tr_tot,nkpt_tr_ibz_tot,nspin,nspin,3*nat))

  file_ep = trim(outdir)//trim(prefix)//trim('_ep_interp.dat')

  inquire(file=file_ep, exist=have_ep)

  if (.not.have_ep) then ! calculate interpolated ep elements and write to file_ep

    ikibz_global = 0
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(at, nat, nspin) &
    !$omp shared(kpts_tr, nkpt_tr_tot, nkpt_tr_ibz_tot, nfs_sheet) &
    !$omp shared(ikibz_2_ik, ik_2_iks, ik_2_ish) &
    !$omp shared(ndegen_q, irvec_q, nrpts_q, gmatL_wann) &
    !$omp shared(u_kint_all) &
    !$omp shared(aep_mat_el) &
    !$omp shared(ikibz_global) &
    !$omp private(ikibz, ik, iks, ish, ib, kpoint) &
    !$omp private(ikp, iksp, ishp, ibp, kpoint_p, qpoint) &
    !$omp private(gmat_aux1, gmat_int, gmat_int_rot) &
    !$omp private(u_kint, u_kqint) &
    !$omp private(irq, facq) &
    !$omp private(iat, is, js)
    do ikibz_do= 1, nkpt_tr_ibz_tot

      !$omp critical
      ikibz_global = ikibz_global + 1
      ikibz = ikibz_global ! ikibz is the k-index over kpoints in the irreducible BZ wedge
      ik = ikibz_2_ik(ikibz) ! ik is the corresponding k-index over kpoints in the full BZ kpts_tr list
      iks = ik_2_iks(ik) ! iks is the corresponding k-index over the kpoints in the FS sheet
      !
      ish = ik_2_ish(ik) ! FS sheet index for k
      ib = nfs_sheet(ish) ! band index for k

      write(*, '(A14,I4,A1,I4,A6,I5,A19)') '|     ik_IBZ: ', ikibz, "/", nkpt_tr_ibz_tot, ' (ik: ', ik, ")                 |"

      kpoint = kpts_tr(:,ik) ! this is cartesians x 2pi/alat. Transform to cryst.
      call cryst_to_cart(1, kpoint, at, -1)

      ! interpolated bands in k and k' calculated above
      ! u_kint_all contains the rotation matrices (row is band and column is wannier,
      ! i.e. the "orbital weights" for each band is in the columns)
      u_kint = transpose(conjg(u_kint_all(:,:,ik))) ! we will need "dagger" below for U(k+q) g U^dagger(k)
      !$omp end critical


      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(at, nat, nspin) &
      !$omp shared(ik_2_iks, ik_2_ish, nfs_sheet) &
      !$omp shared(kpts_tr, nkpt_tr_tot) &
      !$omp shared(ndegen_q, irvec_q, nrpts_q, gmatL_wann) &
      !$omp shared(u_kint_all, u_kint) &
      !$omp shared(aep_mat_el) &
      !$omp shared(ikibz, ib, kpoint) &
      !$omp private(iksp, ishp, ibp, kpoint_p, qpoint) &
      !$omp private(gmat_aux1, gmat_int, gmat_int_rot) &
      !$omp private(u_kqint) &
      !$omp private(irq, facq) &
      !$omp private(iat, is, js)
      do ikp = 1, nkpt_tr_tot

        iksp = ik_2_iks(ikp) ! iksp is the corresponding k-index over the kpoints in the FS sheet
        !
        ishp = ik_2_ish(ikp) ! FS sheet index for k
        ibp = nfs_sheet(ishp) ! band index for k'

        kpoint_p = kpts_tr(:,ikp) ! this is cartesians x 2pi/alat. Transform to cryst.
        call cryst_to_cart(1, kpoint_p, at, -1)

        ! interpolated bands in k and k' calculated above
        ! u_kint_all contains the rotation matrices (row is band and column is wannier,
        ! i.e. the "orbital weights" for each band is in the columns)
        u_kqint = u_kint_all(:,:,ikp)

        qpoint = kpoint_p - kpoint

        ! In the loop over atoms and directions below, I will
        ! interpolate matrix elements for each displacement.
        ! 1. interpolate irvec_q on this qpoint and store in gmat_aux1(:,:,nrpts)
        ! 2. interpolate irvec on kpoint and store in gmat_int.
        ! 3. Rotate Wannier as gmat_int_rot = U(k+q) * gmat_int * U(k)^dagger
        ! 4. The ibp,ib elements are the ep element needed (in canonical atom displacement coordinates)

        ! TODO some of these can be done out of the 3nat and spin loops to go faster
        do js=1,nspin
          do is=1,nspin

            do iat=1,3*nat
              gmat_aux1 = cmplx_0
              gmat_int = cmplx_0
              gmat_int_rot = cmplx_0
              do irq=1,nrpts_q
                facq = exp(cmplx_i*tpi*dot_product(qpoint(:), irvec_q(:,irq)))/real(ndegen_q(irq),dp)
                ! TODO DUDA the spin order in gmatL_wann is correct???
                gmat_aux1(:,:,:) = gmat_aux1(:,:,:) + facq*gmatL_wann(iat,:,:,js,is,irq,:)
              end do
              call wann_FT_1index_1k(kpoint, gmat_aux1(:,:,:), gmat_int(:,:))
              gmat_int_rot(:,:) = matmul(u_kqint, matmul(gmat_int, u_kint))
              aep_mat_el(ikp,ikibz,js,is,iat) = gmat_int_rot(ibp,ib) ! TODO DUDA exclude bands ?? (ver junto al comentario JLB en w90_setup)
            end do

          end do
        end do ! spins

      end do ! k'
      !$omp end parallel do

    end do ! k

    deallocate(gmat_aux1, gmat_int, gmat_int_rot, gep_int)

    ! Save interpolated matrix elements

    unit_ep = find_free_unit()
    open(unit_ep, file=file_ep, status='unknown')
    write(unit_ep,*)'# ik(irr)   jk(full)    is js   g(canonical modes)'
    !
    do ikibz = 1, nkpt_tr_ibz_tot
      !
      ik = ikibz_2_ik(ikibz) ! ik is the corresponding k-index over kpoints in the full BZ kpts_tr list
      iks = ik_2_iks(ik) ! iks is the corresponding k-index over the kpoints in the FS sheet
      ish = ik_2_ish(ik) ! FS sheet index for k
      ib = nfs_sheet(ish) ! band index for k
      !
      do ikp = 1, nkpt_tr_tot
        !
        iksp = ik_2_iks(ikp) ! iksp is the corresponding k-index over the kpoints in the FS sheet
        ishp = ik_2_ish(ikp) ! FS sheet index for k
        ibp = nfs_sheet(ishp) ! band index for k'
        !
        do js=1,nspin
          do is=1,nspin
            write(unit_ep,fmt="(6i6,100e16.6)") ibp, iksp, ikp, ib, iks, ikibz, &
                                                (aep_mat_el(ikp,ikibz, js,is, iat), iat=1,3*nat)
          end do
        end do
        !
      end do ! k'
      !
    end do ! k
    !
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

end program ep_on_trFS_wannier
