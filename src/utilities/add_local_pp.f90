program ep_melements

  use kinds, only: dp
  use intw_input_parameters, only: intw2W_method, nk1, nk2, nk3, &
                                   nq1, nq2, nq3, nqirr, mesh_dir, ph_dir, fc_mat, &
                                   calc_epmat, ep_mat_file, TR_symmetry, &
                                   read_input
  use intw_reading, only: nkpoints_QE, kpoints_QE, nspin, noncolin, gvec, ngm, nsym, &
                          spinorb_mag, can_use_TR, s, nbands, nG_max, nat, alat, &
                          ntyp, amass, nspin, tau, bg, nr1, nr2, nr3, &
                          read_parameters_data_file_xml, get_gvec, &
                          read_kpoints_data_file_xml, &
                          num_bands_intw, set_num_bands
  use intw_pseudo, only: read_all_pseudo, init_pp, phq_init
  use intw_utility, only: get_timing, find_free_unit, switch_indices, &
                          generate_kmesh, conmesurate_and_coarser, ainv, &
                          find_k_1BZ_and_G
  use intw_useful_constants, only: cmplx_0, cmplx_1
  use intw_symmetries, only: full_mesh, IBZ, QE_folder_nosym, QE_folder_sym, nosym_G, sym_G, symlink, &
                             symtable, rtau_index, rtau, rtau_cryst, &
                             rot_atoms, &
                             find_size_of_irreducible_k_set, &
                             find_the_irreducible_k_set, &
                             allocate_symmetry_related_k, &
                             find_inverse_symmetry_matrices_indices, &
                             allocate_and_build_spin_symmetry_matrices, &
                             set_symmetry_relations, multable
  use intw_fft, only: generate_nl, allocate_fft, nl
  use intw_ph, only: nqmesh, qmesh, dvq_local, dvpsi, QE_folder_nosym_q, QE_folder_sym_q, &
                     nosym_G_q, sym_G_q, symlink_q, q_irr, q_irr_cryst, frc, &
                     read_ph_information_xml, readfc, read_allq_dvr, get_dv, &
                     mat_inv_four_t, calculate_local_part_dv
  use intw_allwfcs, only: allocate_and_get_all_irreducible_wfc

  !================================================================================
  !       Declare the variables
  !================================================================================
  implicit none

  !k point related variables
  logical                  :: k_points_consistent
  integer                  :: nk_irr, ik, nkmesh
  real(dp), allocatable     :: kmesh(:, :)
  integer                  :: ikpt_k, ikpt_kq
  real(dp), allocatable     :: kpoints_irr  (:, :)
  real(dp)                 :: kpoint(3), kpoint_cart(3), kqpoint_in_1bz(3), kqpoint_cart(3)
  real(dp)                 :: kpoint_in_1bz(3)

  !q point related variables
  logical                  :: q_points_consistent
  real(dp)                 :: qpoint(3)
  integer                  :: iq, nq1_, nq2_, nq3_
  logical                  :: full_mesh_q, IBZ_q

  !phonon related variables
  character(len=256)       :: fc_file_name
  real(dp)                 :: at_frc(3, 3)
  integer                  :: imode
  real(dp), allocatable    :: qstar(:, :)

  !symmetry variables
  integer, allocatable     :: nsym_sgk(:), sindex_sgk(:, :)
  complex(dp), allocatable :: unit_sym_sgk( :, :, :, :), umat(:, :)

  !time related variables
  real(dp)                 :: time1, time2

  !ep interaction related variables
  integer                  :: ep_unit
  complex(dp), allocatable :: aep_mat_el(:, :, :, :, :, :, :), ep_mat_el(:, :, :, :, :, :)
  complex(dp), allocatable :: fr(:, :, :), fg(:, :, :)

  !wave function realted variables information
  real(dp), allocatable    :: QE_eig_k(:)
  real(dp), allocatable    :: QE_eig_kq(:)
  integer, allocatable     :: list_igk(:)
  integer, allocatable     :: list_igkq(:)
  integer, allocatable     :: list_igk_aux (:)
  integer, allocatable     :: list_igk_orig (:)
  complex(dp), allocatable :: wfc_k (:, :, :) ! nG_max is defined in reading
  complex(dp), allocatable :: wfc_kq (:, :, :)
  complex(dp), allocatable :: wfc_k_aux (:, :, :)
  complex(dp), allocatable :: wfc_k_orig (:, :, :)


  !fft related
  integer                  :: nr(3)
  integer                  :: GKQ_bz(3), G_plusk(3), G_pluskq(3)

  !ep fitxa
  integer                  :: record_lengh, ierr

  !local/aux variables
  integer                  :: nbands_loc
  integer                  :: npw, npwq

  integer                  :: i, j, k
  integer                  :: ig, ibnd, jbnd, ipol, jpol
  logical                  :: read_status
  character(256)           :: method
  character(len=4)         :: iq_loc

  complex(dp), allocatable :: wfc_k_r(:)
  integer :: nG

  complex(dp), external :: zdotc
  !
  real(dp) :: r(3)
  integer :: ir

  complex(dp), allocatable    :: dvq_local_pp    (:, :, :, :)


  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  call get_timing(time1)
  write(*, 20) '====================================================='
  write(*, 20) '|                  program me                       |'
  write(*, 20) '|        ---------------------------------          |'
  write(*, 20) '====================================================='
  write(*, 20) '|    waiting for input file...                      |'
  !
  !================================================================================
  !       read the input file
  !       Read in the necessary information from standard input
  !================================================================================
  !
  call read_input(read_status)
  !
  method=trim(intw2W_method)
  !
  if (read_status) then
     !
     stop
     !
  endif
  !
  !================================================================================
  !       read the parameters from the SCF QE calculation, in
  !       particular, read in the symmetry matrices!!!
  !================================================================================
  !
  call read_parameters_data_file_xml()
  !
  !-call test_symmetry_axis_angle()
  !
  !================================================================================
  !      set up the gvec array, which will contain the global
  !      G-vectors, as well as the FFT code, which is necessary to
  !      generate g_fft_map, which in turn is necessary in the
  !      wavefunction rotation code!
  !================================================================================
  !
  call get_gvec()
  !
  !allocate useful variables
  !
  call allocate_fft()
  !
  !generate some important indices for FFT
  !
  call generate_nl()
  !
  if (method=='CONVOLUTION') then
     !
     write(*, 20) '|       - intw2W_method   = CONVOLUTION             |'
     !
  elseif (method=='FFT') then
     !
     write(*, 20) '|       - intw2W_method   = FFT                     |'
     !
  else
     !
     write(*, 20) '***********************************************'
     write(*, 20) '* UNKNOWN COMPUTATION METHOD:'
     write(*, 20) '* Only "CONVOLUTION" and "FFT" available'
     write(*, 20) '***********************************************'
     !
     stop
     !
  endif
  !
  write(*, 20) '|           ---------------------------------       |'
  !
  !================================================================================
  !       read in the kpoints from the QE folders
  !================================================================================
  !
  ! allocate(kpoints_QE(3, nkpoints_QE))
  !
  ! call read_kpoints_data_file_xml(kpoints_QE)
  !
  !================================================================================
  ! Build the kmesh corresponding to the parameters in the input file
  !================================================================================
  !
  ! nkmesh = nk1*nk2*nk3
  ! allocate(kmesh(3, nkmesh))
  ! call generate_kmesh(kmesh, nk1, nk2, nk3)
  !
  !
  if (nspin==1) then
     !
     write(*, 20) '|       - The calculation is paramagnetic nspin=1   |'
     write(*, 20) '|                                                   |'
     write(*, 20) '|           ---------------------------------       |'
     !
  elseif (nspin==2) then
     !
     write(*, 20) '|       - Spin calculation nspin = 2                |'
     !
     if (noncolin) then
        !
        write(*, 20) '|         Non-collinear Spin calculation            |'
        !
     endif
     !
     write(*, 20) '|           ---------------------------------       |'
     !
  else
     !
     write(*, 20) '*****************************************************'
     write(*, 20) '* ERROR: Allowed values for nspin are 1 or 2        *'
     write(*, 20) '*            program stops.                         *'
     write(*, 20) '*****************************************************'
     !
     stop
     !
  endif
  !
  !================================================================================
  !    Find the size of the irreducible set of k-points (IBZ)
  !    and check that the number of kpoints corresponds to either
  !    a full mesh or the IBZ.
  !================================================================================
  !
  ! call find_size_of_irreducible_k_set(nk1, nk2, nk3, kmesh, nk_irr)
  !
  !This is only for testing: The result for nk_irr is different in both.
  !
  ! allocate(kpoints_irr(3, nk1*nk2*nk3))
  ! call find_the_irreducible_k_set (nk1, nk2, nk3, kmesh, kpoints_irr, nk_irr)
  !
  !
  !================================================================================
  !      allocate the symmetry arrays
  !================================================================================
  !
  ! call allocate_symmetry_related_k(nk1, nk2, nk3, nsym)
  !
  !================================================================================
  !     Compute the indices of the inverse rotation matrices.
  !     This will be useful later in the exectution.
  !================================================================================
  !
  ! call find_inverse_symmetry_matrices_indices()
  !
  !================================================================================
  !     Set up the array spin_symmetry_matrices, which contain
  !     the matrices which must be used to rotate spin
  !================================================================================
  !
  call allocate_and_build_spin_symmetry_matrices(nsym)
  !
  !================================================================================
  !      Fill the symmetry arrays appropriately
  !================================================================================
  !
!Peio
  ! if (spinorb_mag) then
  !    can_use_TR=.true.
  ! endif
!Peio
  !
  ! call set_symmetry_relations(nk1, nk2, nk3, nkpoints_QE, kpoints_QE, kmesh, k_points_consistent, &
  !                             QE_folder_nosym, QE_folder_sym, &
  !                             nosym_G, sym_G, symlink, full_mesh, IBZ)
  !
  !================================================================================
  !       Tell the user what is in the QE folders
  !================================================================================
  !
  ! if (full_mesh .and. IBZ) then
  !   write(*, 20) '|       - the kpoints present in the QE folders     |'
  !   write(*, 20) '|         are consistent with a full 1BZ and a      |'
  !   write(*, 20) '|         IBZ has also been found.                  |'
  !   write(*, 20) '|           ---------------------------------       |'
  ! else if(IBZ) then
  !   write(*, 20) '|       - the kpoints present in the QE folders     |'
  !   write(*, 20) '|         are consistent with an IBZ.               |'
  !   write(*, 20) '|           ---------------------------------       |'
  ! else
  !   write(*, 20) '**********************************************************'
  !   write(*, 20) '* The kpoints present in the QE folders are not consistent'
  !   write(*, 20) '* with the parameters of the input file!                 '
  !   write(*, 20) '**********************************************************'
  !   write(*, 20) '* debug information:                                *'
  !   write(*, *) '*        nkpoints_QE = ', nkpoints_QE
  !   write(*, *) '*        nkmesh      = ', nkmesh
  !   write(*, *) '*        nk_irr      = ', nk_irr, nsym, tr_symmetry
  !   write(*, *) '*        IBZ         = ', IBZ
  !   stop
  ! end if
  !
  !================================================================================
  !      Calculate the multiplication talble for symmetry operations
  !================================================================================
  !
  ! call multable(nsym, s, symtable)
  !
  !================================================================================
  !       Tell the user what is in the QE folders
  !================================================================================
  !
  ! if (full_mesh .and. IBZ) then
  !   !
  !   write(*, 20) '|       - the kpoints present in the QE folders     |'
  !   write(*, 20) '|         are consistent with a full 1BZ and a      |'
  !   write(*, 20) '|         IBZ has also been found.                  |'
  !   write(*, 20) '|           ---------------------------------       |'
  !   !
  ! elseif(IBZ) then
  !   !
  !   write(*, 20) '|       - the kpoints present in the QE folders     |'
  !   write(*, 20) '|         are consistent with an IBZ.               |'
  !   write(*, 20) '|           ---------------------------------       |'
  !   !
  ! else
  !   !
  !   write(*, 20) '**********************************************************'
  !   write(*, 20) '* The kpoints present in the QE folders are not consistent'
  !   write(*, 20) '* with the parameters of the input file!                 '
  !   write(*, 20) '**********************************************************'
  !   write(*, 20) '* debug information:                                *'
  !   write(*, *) '*        nkpoints_QE = ', nkpoints_QE
  !   write(*, *) '*        nkmesh      = ', nkmesh
  !   write(*, *) '*        nk_irr      = ', nk_irr, nsym, tr_symmetry
  !   write(*, *) '*        IBZ         = ', IBZ
  !   !
  !   stop
  !   !
  ! endif
  !
  !================================================================================
  !       Check that the requested calculation is possible
  !================================================================================
  !
  write(*, 20) '|       - reading pseudopotentials from UPF files   |'
  write(*, 20) '|            (defined in .save data files)          |'
  write(*, 20) '|                                                   |'
  !
  call read_all_pseudo ()
  !if (.not.lspinorb) call average_pp(ntyp)
  !
  write(*, 20) '|                    PPs are OK                     |'
  write(*, 20) '|           ---------------------------------       |'
  !
  ! if (.not.conmesurate_and_coarser(nk1, nk2, nk3, nq1, nq2, nq3)) then
  !    !
  !    write(*, 20) '**********************************************************'
  !    write(*, 20) '*ERROR                                                  '
  !    write(*, 20) '*   the electron k and phonon q are not                 '
  !    write(*, 20) '*   conmesurate and the k grid does not contain         '
  !    write(*, 20) '*   the phonon q grid                                   '
  !    write(*, 20) '**********************************************************'
  !    !
  !    stop
  !    !
  ! endif
  ! !
  !
  !================================================================================
  !       Allocate wfc related files
  !================================================================================
  !
  !!!JLB: To be discussed
  !num_bands = nbands
  !nbands_loc = num_bands
  !
  ! call set_num_bands()
  ! nbands_loc = num_bands_intw ! JLB: Why do we need this?
  !
  allocate(list_igk(nG_max))
  allocate(list_igkq(nG_max))
  allocate(list_igk_aux(nG_max))
  allocate(list_igk_orig(nG_max))
  !
  ! allocate(wfc_k(nG_max, num_bands_intw, nspin))
  ! allocate(wfc_kq(nG_max, num_bands_intw, nspin))
  ! allocate(wfc_k_aux(nG_max, num_bands_intw, nspin))
  ! allocate(wfc_k_orig(nG_max, num_bands_intw, nspin))
  ! !
  ! allocate(wfc_k_r(nr1*nr2*nr3))
  ! !
  ! allocate(QE_eig_k(num_bands_intw))
  ! allocate(QE_eig_kq(num_bands_intw))
  ! !
  ! allocate(fr(nr1*nr2*nr3, num_bands_intw, nspin))
  ! allocate(fg(nG_max, num_bands_intw, nspin))
  !
  !================================================================================
  !       Read all the information about phonons and pseudopotentials
  !================================================================================
  !
  call read_ph_information_xml()
  !
  !================================================================================
  !       Read the force constant matrix from the QE directory
  !================================================================================
  !
  ! fc_file_name = trim(trim(mesh_dir)//trim(ph_dir)//trim(fc_mat))
  !
!  if (.not.lspinorb) then
!     call readfc(fc_file_name, nq1_, nq2_, nq3_, nat, alat, at_frc, ntyp, amass)
!  else
!     call read_fc_from_XML ()
!  endif
  !
  ! call readfc(fc_file_name, nq1_, nq2_, nq3_, nat, alat, at_frc, ntyp, amass)
  !
  !================================================================================
  ! Build the phonon qmesh corresponding to the parameters in the input file
  !================================================================================
  !
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3, nqmesh))
  !
  !hauek behekoak ph_module.mod-n definituta daude eta aldagai globalak dira kontuz.
  !
  allocate(dvq_local(nr1*nr2*nr3, 3*nat, nspin, nspin))
  allocate(dvq_local_pp(nr1*nr2*nr3, 3*nat, nspin, nspin))
  ! allocate(dvpsi(nG_max, nbands_loc, nspin, nspin, 3*nat))
  !
  call generate_kmesh(qmesh, nq1, nq2, nq3)
  !
  ! allocate(ep_mat_el(nk1*nk2*nk3, num_bands_intw, num_bands_intw, nspin, nspin, 3*nat))
  ! allocate(aep_mat_el(nqmesh, nk1*nk2*nk3, num_bands_intw, num_bands_intw, nspin, nspin, 3*nat))
  !
  ! aep_mat_el(:, :, :, :, :, :, :) = cmplx_0
  !
  !================================================================================
  ! Below, the table of rotations for each atom and symmetry.
  !================================================================================
  !
  allocate(rtau_index(nat, nsym))
  allocate(rtau(3, nsym, nat))
  allocate(rtau_cryst(3, nsym, nat))
  !
  nr = (/nr1, nr2, nr3/)
  !
  call rot_atoms(nat, nsym, tau)
  !
  !
  !================================================================================
  ! Symmetry realtions between irreducible q directory file number and the full mesh.
  !================================================================================
  !
  allocate(QE_folder_nosym_q(nqmesh))
  allocate(QE_folder_sym_q(nqmesh))
  allocate(nosym_G_q(3, nqmesh))
  allocate(sym_G_q(3, nqmesh))
  allocate(symlink_q(nqmesh, 2))
  allocate(q_irr_cryst(3, nqirr))
  allocate(qstar(3, 48))
  !
  do iq=1, nqirr
     q_irr_cryst(:, iq) = matmul(ainv(bg), q_irr(:, iq))
  enddo
  !
  !================================================================================
  ! Read de dv induced potentials from QE directories.
  !================================================================================
  !
  call read_allq_dvr(nqirr, 3*nat)
  !
  !================================================================================
  ! We need the symmtry relations between irreducible/full in order to obtain
  ! the induced potential for any q point.
  !================================================================================
  !
  call set_symmetry_relations(nq1, nq2, nq3, nqirr, q_irr_cryst, qmesh, q_points_consistent, &
                              QE_folder_nosym_q, QE_folder_sym_q, &
                              nosym_G_q, sym_G_q, symlink_q, full_mesh_q, IBZ_q)
  !
  !================================================================================
  ! testing ...dv and mat_ep
  !================================================================================
  !
  ! call allocate_and_get_all_irreducible_wfc()
  !
  !================================================================================
  !      Calculate the small group of symm. operations for the irreducible k points
  !================================================================================
  !
  ! allocate(nsym_sgk(nk_irr))
  ! allocate(sindex_sgk(nk_irr, nsym))
  ! allocate(unit_sym_sgk(nkmesh, nsym, num_bands_intw, num_bands_intw))
  ! allocate(umat(num_bands_intw, num_bands_intw))
  ! !
  ! unit_sym_sgk = cmplx_0
  ! !
  ! do ik=1, nkmesh
  !    do ibnd=1, num_bands_intw
  !       !
  !       unit_sym_sgk(ik, 1:nsym, ibnd, ibnd) = cmplx_1
  !       !
  !    enddo !ibnd
  ! enddo !ik
  !
  call allocate_nlpot
  !
  call allocate_phq
  !
  call init_pp

  ep_mat_el = cmplx_0
  !
  !-------------------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------------------
  !
    !
    do iq=1, nqmesh
      !
      if (                iq <   10) write(iq_loc, "(i1)")iq
      if ( 10 <= iq .and. iq <  100) write(iq_loc, "(i2)")iq
      if (100 <= iq .and. iq < 1000) write(iq_loc, "(i3)")iq
      !
      ep_unit=find_free_unit()
      inquire(iolength=record_lengh) ep_mat_el
      open(unit=ep_unit, file=trim(trim(mesh_dir)//trim(ep_mat_file)//trim('_')//adjustl(iq_loc)), iostat=ierr, &
      form='unformatted', status='unknown', access='direct', recl=record_lengh)
      !
      qpoint(:)=qmesh(:, iq)
      write(*, "(a, i4, 100f12.6)") "qpoint", iq, qpoint
      !
      !-matrize dinamikoa kalkulatu q puntu jakin batentzat.
      !
      call set_nqxq(qpoint)
      !
      !-potentzialaren alde induzitua kalkulatu simetria erabiliz (errotazioz beharrezkoa izanez).
      !
      dvq_local = cmplx_0
      call get_dv(iq, qpoint, 3*nat, nspin, dvq_local)
      !
      !-alde induzituari (goian), KB pseudopotentzialaren(pp) deribatuaren ALDE LOKALA gehitu.
      !
      ! qpoint(1) = 0.25_dp
      ! qpoint(3) = -0.5_dp
      ! call phq_init(matmul(bg, qpoint))
      ! dvq_local_pp = cmplx_0
      ! call calculate_local_part_dv(qpoint, nat, nspin, dvq_local_pp)
      !
      ! do k=1, nr3
      !   if (k/=1) cycle
      !   do j=1, nr2
      !     if (j/=1) cycle
      !     do i=1, nr1
      !       r(1) = real(i-1, dp)/nr1
      !       r(2) = real(j-1, dp)/nr2
      !       r(3) = real(k-1, dp)/nr3
      !       ir = i + nr1*(j-1) + nr1*nr2*(k-1)
      !       write(4000+iq, "(f10.6, 2f12.8)") sqrt(sum(r**2)), dvq_local(ir, 1, 1, 1)
      !     enddo
      !   enddo
      ! enddo
      do imode=1, 3*nat
        if (nspin==2) then
          write(4000+iq, "(i10, 4f12.8)") ( nr1*nr2*nr3*(imode-1) + ir, abs(dvq_local(ir, imode, 1, 1)), abs(dvq_local(ir, imode, 1, 2)), abs(dvq_local(ir, imode, 2, 1)), abs(dvq_local(ir, imode, 2, 2)), ir=1, nr1*nr2*nr3)
        else
          write(4000+iq, "(i10, 1f12.8)") ( nr1*nr2*nr3*(imode-1) + ir, abs(dvq_local(ir, imode, 1, 1)), ir=1, nr1*nr2*nr3)
        endif
        ! write(5000+iq, "(i6, 1f12.8)") ( nr1*nr2*nr3*(imode-1) + ir, abs(dvq_local(ir, imode, 1, 1) + dvq_local_pp(ir, imode, 1, 1)), ir=1, nr1*nr2*nr3)
        ! write(6000+iq, "(i6, 1f12.8)") ( nr1*nr2*nr3*(imode-1) + ir, abs(dvq_local_pp(ir, imode, 1, 1)), ir=1, nr1*nr2*nr3)
      enddo
      ! dvq_local_pp(:, 1, 1, 1) = cmplx_0
      ! call calculate_local_part_v( dvq_local_pp(:, 1, 1, 1))
      ! write(6000+iq, "(i6, 1f12.8)") ( nr1*nr2*nr3*(imode-1) + ir, abs(dvq_local_pp(ir, 1, 1, 1)), ir=1, nr1*nr2*nr3)
      ! dvq_local_pp(:, 1, 1, 1) = cmplx_0
      ! call calculate_local_part_ni( dvq_local_pp(:, 1, 1, 1))
      ! write(7000+iq, "(i6, 1f12.8)") ( nr1*nr2*nr3*(imode-1) + ir, abs(dvq_local_pp(ir, 1, 1, 1)), ir=1, nr1*nr2*nr3)

      ! dvq_local_pp = cmplx_0
      ! call calculate_local_part_dv(qpoint, nat, nspin, dvq_local_pp)
      ! do imode=1, 3*nat
      !   write(8000+iq, "(i8, 1f18.8)") ( nr1*nr2*nr3*(imode-1) + ir, aimag(dvq_local_pp(ir, imode, 1, 1)), ir=1, nr1*nr2*nr3)
      ! enddo
      ! dvq_local_pp = cmplx_0
      ! call calculate_local_part_dv_ni(qpoint, nat, nspin, dvq_local_pp)
      ! do imode=1, 3*nat
      !   write(9000+iq, "(i8, 1f18.8)") ( nr1*nr2*nr3*(imode-1) + ir, aimag(dvq_local_pp(ir, imode, 1, 1)), ir=1, nr1*nr2*nr3)
      ! enddo


      !
    enddo !iq


  ! deallocate (ep_mat_el)
  ! deallocate (aep_mat_el)

  ! deallocate(nsym_sgk)
  ! deallocate(sindex_sgk)
  ! deallocate(kpoints_irr )

  !call deallocate_upfeak ()
  !================================================================================
  !       Finish
  !================================================================================
  call get_timing(time2)
  write(*, 20) '|                     ALL DONE                       |'
  write(*, 30) '|     total time: ', time2-time1, ' seconds            |'
  write(*, 20) '====================================================='


20 format(A)
30 format(A, F8.2, 6X, A)

contains

  subroutine calculate_local_part_ni( v_local )
    ! We have dV_scf as input and we add to it the derivative of the PP   !

    use intw_useful_constants, only: tpi, cmplx_i
    use intw_reading, only: nr1, nr2, nr3, ngm, ityp, alat
    use intw_fft, only: nl, gvec_cart
    use intw_pseudo, only: vlocq

    implicit none

    external :: cfftnd

    !I/O variables

    complex(dp), intent(inout) :: v_local(nr1*nr2*nr3) ! spin idependentea da baina koherentzia mantenduko dugu.

    !local variables

    integer :: ia, ig, nt
    complex(dp) :: aux(nr1*nr2*nr3), gtau

    do ia=1,nat
      !
      nt = ityp(ia)
      aux = cmplx_0
      !
      do ig=1,ngm
        gtau = exp(-cmplx_i*sum(gvec_cart(:,ig)*tau(:,ia))*tpi)
        aux(nl(ig)) = aux(nl(ig)) + vlocq(ig, nt)*gtau
      enddo !ig
      !
      !
      call cfftnd(3,(/nr1,nr2,nr3/),1,aux)
      !
      v_local = v_local + aux
      !
    enddo !ia

  end subroutine calculate_local_part_ni


  subroutine calculate_local_part_ni_2( v_local )
    ! We have dV_scf as input and we add to it the derivative of the PP   !

    use intw_useful_constants, only: tpi, cmplx_i
    use intw_reading, only: nr1, nr2, nr3, ngm, ityp, alat
    use intw_fft, only: nl, gvec_cart
    use intw_pseudo, only: vlocq

    implicit none

    external :: cfftnd

    !I/O variables

    complex(dp), intent(inout) :: v_local(nr1*nr2*nr3) ! spin idependentea da baina koherentzia mantenduko dugu.

    !local variables

    integer :: ia, ig, nt
    complex(dp) :: aux(nr1*nr2*nr3), gtau

    do ia=1,nat
      !
      nt = ityp(ia)
      aux = cmplx_0
      !
      do ig=1,ngm
        gtau = exp(-cmplx_i*sum(gvec_cart(:,ig)*tau(:,ia))*tpi)
        aux(nl(ig)) = aux(nl(ig)) + vlocq(ig, nt)*gtau
      enddo !ig
      !
      !
      call cfftnd(3,(/nr1,nr2,nr3/),1,aux)
      !
      v_local = v_local + aux
      !
    enddo !ia

  end subroutine calculate_local_part_ni_2


  subroutine calculate_local_part_dv_ni( qpoint, nat, nspin, dv_local )
    ! We have dV_scf as input and we add to it the derivative of the PP   !

    use intw_useful_constants, only: tpi, cmplx_i
    use intw_reading, only: nr1, nr2, nr3, ngm, ityp, alat, bg, tpiba
    use intw_fft, only: nl, gvec_cart, eigts1, eigts2, eigts3, mill
    use intw_pseudo, only: vlocq
    use intw_ph, only: eigqts

    implicit none

    external :: cfftnd

    !I/O variables
    real(kind=dp),intent(in) :: qpoint(1:3) ! crystal coord.
    integer, intent(in) :: nat, nspin
    complex(dp), intent(inout) :: dv_local(nr1*nr2*nr3, 3*nat, nspin, nspin) ! spin idependentea da baina koherentzia mantenduko dugu.

    !local variables
    integer :: imode, ia, ig, nt, is, na, ipol
    real(dp) :: qcart(3) ! qpoint in cart.
    complex(dp) :: aux(nr1*nr2*nr3), qtau, gtau
    ! real(kind=dp), parameter :: dx = 0.01_dp
    ! real(kind=dp) :: tau_p(3,nat), tau_m(3,nat)
    ! complex(dp) :: aux1(nr1*nr2*nr3), gtau1
    ! complex(dp) :: aux2(nr1*nr2*nr3), gtau2

    qcart = matmul(bg, qpoint)

    do imode = 1, 3*nat
      !
      na = (imode-1)/3+1
      ipol = modulo(imode-1, 3)+1
      nt = ityp(na)
      !
      ! tau_p = tau
      ! tau_p(ipol, na) = tau(ipol, na) + dx/alat
      ! tau_m = tau
      ! tau_m(ipol, na) = tau(ipol, na) - dx/alat
      !
      aux = cmplx_0
      ! aux1 = cmplx_0
      ! aux2 = cmplx_0
      !
      qtau = exp(-cmplx_i*sum(qcart(:)*tau(:,na))*tpi)
      do ig=1,ngm
        gtau = exp(-cmplx_i*sum(gvec_cart(:,ig)*tau(:,na))*tpi)
        ! gtau = eigts1(mill(1,ig),na)*eigts2(mill(2,ig),na)*eigts3(mill(3,ig),na)
        ! print*, gtau, eigts1(mill(1,ig),na)*eigts2(mill(2,ig),na)*eigts3(mill(3,ig),na)
        aux(nl(ig)) = aux(nl(ig)) -tpiba*cmplx_i*(qcart(ipol)+gvec_cart(ipol, ig))*qtau*gtau*vlocq(ig, nt)
        ! gtau1 = exp(-cmplx_i*sum(gvec_cart(:,ig)*tau_p(:,na))*tpi)
        ! aux1(nl(ig)) = aux1(nl(ig)) + vlocq(ig, nt)*gtau1
        ! gtau2 = exp(-cmplx_i*sum(gvec_cart(:,ig)*tau_m(:,na))*tpi)
        ! aux2(nl(ig)) = aux2(nl(ig)) + vlocq(ig, nt)*gtau2
      enddo !ig
      !
      !
      call cfftnd(3, (/nr1, nr2, nr3/), 1, aux)
      ! call cfftnd(3,(/nr1,nr2,nr3/),1,aux1)
      ! call cfftnd(3,(/nr1,nr2,nr3/),1,aux2)
      !
      do is=1,nspin
        dv_local(:, imode, is, is) = dv_local(:, imode, is, is) + aux(:)
      enddo
      ! dv_local(:,imode) = dv_local(:,imode) + (aux1 - aux2)/(2.0_dp*dx)
      !
    enddo !imode

  end subroutine calculate_local_part_dv_ni


  subroutine calculate_local_part_v( v_local )
    ! We have dV_scf as input and we add to it the derivative of the PP   !

    use intw_useful_constants, only: tpi, cmplx_i
    use intw_reading, only: nr1, nr2, nr3, ngm, ityp, alat
    use intw_fft, only: eigts1, eigts2, eigts3, nl, mill, gvec_cart
    use intw_pseudo, only: vlocq

    implicit none

    external :: cfftnd

    !I/O variables

    complex(dp),intent(inout) :: v_local(nr1*nr2*nr3) ! spin idependentea da baina koherentzia mantenduko dugu.

    !local variables

    integer :: na, ig, nt
    complex(dp) :: aux(nr1*nr2*nr3), gtau

    !
    do na=1,nat
       !
       nt = ityp(na)
       !
       aux = cmplx_0
       !
       do ig=1,ngm
          !
          gtau = eigts1(mill(1,ig),na)*eigts2(mill(2,ig),na)*eigts3(mill(3,ig),na)
          aux(nl(ig)) = aux(nl(ig)) + vlocq(ig,nt)*gtau
          !
       enddo !ig
       !

       call cfftnd(3, (/nr1,nr2,nr3/), 1, aux)
       !
       v_local = v_local + aux
       !
    enddo !na
    !
  end subroutine calculate_local_part_v


end program ep_melements
