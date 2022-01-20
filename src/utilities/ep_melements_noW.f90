program ep_melements

  use kinds, only: dp
  use intw_input_parameters, only: intw2W_method, intw2W_fullzone, nk1, nk2, nk3, &
                                   nq1, nq2, nq3, nqirr, mesh_dir, ph_dir, fc_mat, &
                                   calc_epmat, ep_mat_file, TR_symmetry, &
                                   read_input
  use intw_reading, only: nkpoints_QE, nspin, noncolin, gvec, ngm, nsym, &
                          spinorb_mag, can_use_TR, s, nbands, nG_max, nat, alat, &
                          ntyp, amass, npol, tau, bg, nr1, nr2, nr3, &
                          read_parameters_data_file_xml, get_ngm, get_gvec, &
                          read_kpoints_data_file_xml
  use intw_pseudo, only: vkb, vkqb, read_all_pseudo
  use intw_utility, only: get_timing, find_free_unit, switch_indices, &
                          generate_kmesh, conmesurate_and_coarser, ainv, &
                          find_k_1BZ_and_G
  use intw_useful_constants, only: cmplx_0, cmplx_1
  use intw_symmetries, only: full_mesh, IBZ, QE_folder_sym, sym_G, symlink, &
                             symtable, rtau_index, rtau, tau_cryst, rtau_cryst, &
                             rot_atoms, &
                             find_size_of_irreducible_k_set, &
                             find_the_irreducible_k_set, &
                             allocate_symmetry_related_k, &
                             find_inverse_symmetry_matrices_indices, &
                             allocate_and_build_spin_symmetry_matrices, &
                             set_symmetry_relations, multable
  use intw_fft, only: generate_nl, allocate_fft, nl
  use intw_ph, only: nqmesh, qmesh, dvq_local, dvpsi, QE_folder_sym_q, sym_G_q, &
                     symlink_q, q_irr, q_irr_cryst, frc, &
                     read_ph_information_xml, readfc, read_allq_dvr, get_dv, &
                     mat_inv_four_t, calculate_local_part_dv
  ! use intw_w90
  use w90_io, only: io_error,io_file_unit,io_time,io_stopwatch
  use w90_parameters, only: num_bands
  use intw_allwfcs, only: allocate_and_get_all_irreducible_wfc, get_psi_general_k_all_wfc
  use intw_matrix_elements, only: get_elec_phon_matrix_element_convolution
  use intw_allwfcs, only: get_psi_general_k_all_wfc

  !================================================================================
  !       Declare the variables
  !================================================================================
  implicit none

  !k point related variables
  logical                  :: k_points_consistent
  integer                  :: nk_irr , ik, nkmesh
  real(dp),allocatable     :: kmesh(:,:), kpoints_QE(:,:)
  integer                  :: ikpt_k, ikpt_kq
  real(dp),allocatable     :: kpoints_irr  (:,:)
  real(dp)                 :: kpoint(3), kpoint_cart(3), kqpoint_in_1bz(3), kqpoint_cart(3)
  real(dp)                 :: kpoint_in_1bz(3)

  !q point related variables
  real(dp)                 :: qpoint(3)

  !phonon related variables
  character(len=256)       :: fc_file_name
  real(dp)                 :: at_frc(3,3)
  integer                  :: imode
  complex(dp), allocatable :: dyn_q(:,:)
  real(dp), allocatable    :: w2(:)
  real(dp), allocatable    :: qstar(:,:)

  !symmetry variables
  integer, allocatable     :: nsym_sgk(:), sindex_sgk(:,:)
  complex(dp), allocatable :: unit_sym_sgk( :,:,:,:), umat(:,:)

  !time related variables
  real(dp)                 :: time1, time2

  !q point related variables
  integer                  :: iq, nq1_, nq2_, nq3_

  !ep interaction related variables
  integer                  :: ep_unit
  logical                  :: q_points_consistent
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:,:,:), ep_mat_el(:,:,:,:,:,:)
  complex(dp), allocatable :: fr(:,:,:), fg(:,:,:)

  !wave function realted variables information
  integer, allocatable     :: list_igk (:)
  integer, allocatable     :: list_igkq(:)
  integer, allocatable     :: list_igk_aux (:)
  integer, allocatable     :: list_igk_orig (:)

  complex(dp), allocatable :: wfc_k (:,:,:) ! nG_max is defined in reading
  complex(dp), allocatable :: wfc_kq (:,:,:)
  complex(dp), allocatable :: wfc_k_aux (:,:,:)
  complex(dp), allocatable :: wfc_k_orig (:,:,:)

  real(dp), allocatable    :: QE_eig_k(:)
  real(dp), allocatable    :: QE_eig_kq(:)

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

  complex(dp),allocatable :: wfc_k_r(:)
  integer :: nG

  complex(dp), external :: zdotc

  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  call get_timing(time1)
  write(*,20) '====================================================='
  write(*,20) '|                  program me                       |'
  write(*,20) '|        ---------------------------------          |'
  write(*,20) '====================================================='
  write(*,20) '|    waiting for input file...                      |'
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
!haritz: ngm irakurri da read_parameters_data_file_xml()-en, eztao get_ngm() erabili beharrik
  call get_ngm()
  allocate(gvec(3,ngm)) ! it is deallocated in the bottom of this main code.
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
     write(*,20) '|       - intw2W_method   = CONVOLUTION             |'
     !
  elseif (method=='FFT') then
     !
     write(*,20) '|       - intw2W_method   = FFT                     |'
     !
  else
     !
     write(*,20) '***********************************************'
     write(*,20) '* UNKNOWN COMPUTATION METHOD:'
     write(*,20) '* Only "CONVOLUTION" and "FFT" available'
     write(*,20) '***********************************************'
     !
     stop
     !
  endif
  !
  write(*,20) '|           ---------------------------------       |'
  !
  !================================================================================
  !       read in the kpoints from the QE folders
  !================================================================================
  !
  allocate(kpoints_QE(3,nkpoints_QE))
  !
  call read_kpoints_data_file_xml(kpoints_QE)
  !
  !================================================================================
  ! Build the kmesh corresponding to the parameters in the input file
  !================================================================================
  !
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)
  !
  !
  if (nspin==1) then
     !
     write(*,20) '|       - The calculation is paramagnetic nspin=1   |'
     write(*,20) '|                                                   |'
     write(*,20) '|           ---------------------------------       |'
     !
  elseif (nspin==2) then
     !
     write(*,20) '|       - Spin calculation nspin = 2                |'
     !
     if (noncolin) then
        !
        write(*,20) '|         Non-collinear Spin calculation            |'
        !
     endif
     !
     write(*,20) '|           ---------------------------------       |'
     !
  else
     !
     write(*,20) '*****************************************************'
     write(*,20) '* ERROR: Allowed values for nspin are 1 or 2        *'
     write(*,20) '*            program stops.                         *'
     write(*,20) '*****************************************************'
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
  call find_size_of_irreducible_k_set(nk1,nk2,nk3,kmesh,nk_irr)
  !
  !This is only for testing: The result for nk_irr is different in both.
  !
  allocate(kpoints_irr(3,nk1*nk2*nk3))
  call find_the_irreducible_k_set (nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
  !
  if (nkpoints_QE/=nkmesh .and. nkpoints_QE/=nk_irr) then
     !
     ! the points in the folder are not consistent with the
     ! input file. Stop the program!
     !
     write(*,20) '*****************************************************'
     write(*,20) '*      The number of kpoints present in the QE      *'
     write(*,20) '*      folders are not consistent with a full       *'
     write(*,20) '*      Brillouin Zone or an irreducible Brillouin   *'
     write(*,20) '*      zone! Review your input...                   *'
     write(*,20) '*                   Program stops.                  *'
     write(*,20) '*****************************************************'
     write(*,20) '* debug information:                                *'
     write(*,*) '*        nkpoints_QE = ',nkpoints_QE
     write(*,*) '*        nkmesh      = ',nkmesh
     write(*,*) '*        nk_irr      = ',nk_irr, nsym, tr_symmetry
     !
     stop
     !
  elseif (nkpoints_QE==nk_irr) then
     !
     ! The points in the QE folders *could* be a valid choice for the IBZ;
     ! this must be checked!
     !
     full_mesh = .false.
     IBZ       = .true.
     !
  elseif (nkpoints_QE==nkmesh) then
     !
     ! The points in the QE folders *could* be consistent with a full mesh;
     ! this must be checked!
     !
     full_mesh = .true.
     IBZ       = .false.
     !
  endif
  !
  !================================================================================
  !      allocate the symmetry arrays
  !      CAREFUL! the subroutine needs to know the global value of "full_mesh",
  !      so it is crucial that this allocation occurs AFTER setting full_mesh
  !================================================================================
  !
  call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)
  !
  !================================================================================
  !     Compute the indices of the inverse rotation matrices.
  !     This will be useful later in the exectution.
  !================================================================================
  !
  call find_inverse_symmetry_matrices_indices()
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
  if (spinorb_mag) then
     can_use_TR=.true.
  endif
!Peio
  !
  call set_symmetry_relations(nk1,nk2,nk3,nkpoints_QE,kpoints_QE,kmesh, &
                         k_points_consistent,QE_folder_sym,sym_G,symlink)
  !
  !================================================================================
  !      Calculate the multiplication talble for symmetry operations
  !================================================================================
  !
  call multable(nsym,s,symtable)
  !
  !================================================================================
  !       Tell the user what is in the QE folders
  !================================================================================
  !
  if (full_mesh .and. IBZ) then
     !
     write(*,20) '|       - the kpoints present in the QE folders     |'
     write(*,20) '|         are consistent with a full 1BZ and a      |'
     write(*,20) '|         IBZ has also been found.                  |'
     write(*,20) '|           ---------------------------------       |'
     !
  elseif(IBZ) then
     !
     write(*,20) '|       - the kpoints present in the QE folders     |'
     write(*,20) '|         are consistent with an IBZ.               |'
     write(*,20) '|           ---------------------------------       |'
     !
  else
     !
     write(*,20)      '**********************************************************'
     write(*,20)      '* The kpoints present in the QE folders are not consistent'
     write(*,20)      '* with the parameters of the input file!                 '
     write(*,20)      '**********************************************************'
     write(*,20) '* debug information:                                *'
     write(*,*) '*        nkpoints_QE = ',nkpoints_QE
     write(*,*) '*        nkmesh      = ',nkmesh
     write(*,*) '*        nk_irr      = ',nk_irr, nsym, tr_symmetry
     write(*,*) '*        IBZ         = ',IBZ
     !
     stop
     !
  endif
  !
  !================================================================================
  !       Check that the requested calculation is possible
  !================================================================================
  !
  if (intw2W_fullzone) then
     !
     write(*,20) '|       - intw2W_fullzone = .true.                  |'
     !
     if (full_mesh) then
        !
        write(*,20) '|         all k-points are explicitely calculated   |'
        write(*,20) '|         no symmetry is assumed.                   |'
        write(*,20) '|             (This is mostly for testing)          |'
        write(*,20) '|           ---------------------------------       |'
        !
     else
        !
        write(*,20)      '**********************************************************'
        write(*,20)      '* A full mesh is not present in the QE folders!          '
        write(*,20)      '* The requested calculation is impossible.               '
        write(*,20)      '*                   program stops.                       '
        write(*,20)      '**********************************************************'
        !
        stop
        !
     endif
     !
  else
     !
     write(*,20) '|       - intw2W_fullzone = .false.                 |'
     write(*,20) '|         Symmetries will be utilized.              |'
     write(*,20) '|           ---------------------------------       |'
     !
  endif
  !
  write(*,20) '|       - reading pseudopotentials from UPF files   |'
  write(*,20) '|            (defined in .save data files)          |'
  write(*,20) '|                                                   |'
  !
  call read_all_pseudo ()
  !if (.not.lspinorb) call average_pp(ntyp)
  !
  write(*,20) '|                    PPs are OK                     |'
  write(*,20) '|           ---------------------------------       |'
  !
  if (.not.conmesurate_and_coarser(nk1,nk2,nk3,nq1,nq2,nq3)) then
     !
     write(*,20)      '**********************************************************'
     write(*,20)       '*ERROR                                                  '
     write(*,20)       '*   the electron k and phonon q are not                 '
     write(*,20)       '*   conmesurate and the k grid does not contain         '
     write(*,20)       '*   the phonon q grid                                   '
     write(*,20)      '**********************************************************'
     !
     stop
     !
  endif
  !
  !
  !================================================================================
  !       Allocate wfc related files
  !================================================================================
  !
  num_bands = nbands
  nbands_loc=num_bands
  !
  allocate (list_igk(nG_max))
  allocate (list_igkq(nG_max))
  allocate (list_igk_aux(nG_max))
  allocate (list_igk_orig(nG_max))
  !
  allocate (wfc_k(nG_max,num_bands,nspin))
  allocate (wfc_kq(nG_max,num_bands,nspin))
  allocate (wfc_k_aux(nG_max,num_bands,nspin))
  allocate (wfc_k_orig(nG_max,num_bands,nspin))
  !
  allocate (wfc_k_r(nr1*nr2*nr3))
  !
  allocate (QE_eig_k(num_bands))
  allocate (QE_eig_kq(num_bands))
  !
  allocate (fr(nr1*nr2*nr3,num_bands,nspin))
  allocate (fg(nG_max,num_bands,nspin))
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
  fc_file_name=trim(trim(mesh_dir)//trim(ph_dir)//trim(fc_mat))
  !
!  if (.not.lspinorb) then
!     call readfc(fc_file_name,nq1_,nq2_,nq3_,nat,alat,at_frc,ntyp,amass)
!  else
!     call read_fc_from_XML ()
!  endif
  !
  call readfc(fc_file_name,nq1_,nq2_,nq3_,nat,alat,at_frc,ntyp,amass)
  !
  !================================================================================
  ! Build the phonon qmesh corresponding to the parameters in the input file
  !================================================================================
  !
  nqmesh=nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  !
  !hauek behekoak ph_module.mod-n definituta daude eta aldagai globalak dira kontuz.
  !
  allocate(dvq_local(nr1*nr2*nr3,3*nat,npol,npol))
  allocate(dvpsi(nG_max,nbands_loc,npol,npol,3*nat))
  !
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !
  allocate(ep_mat_el(nk1*nk2*nk3,num_bands,num_bands,nspin,nspin,3*nat))
  allocate(aep_mat_el(nqmesh,nk1*nk2*nk3,num_bands,num_bands,nspin,nspin,3*nat))
  !
  aep_mat_el(:,:,:,:,:,:,:)=cmplx_0
  !
  !================================================================================
  ! Below, the table of rotations for each atom and symmetry.
  !================================================================================
  !
  allocate(rtau_index(nat,nsym))
  allocate(rtau(3,nsym,nat))
  allocate(tau_cryst(3,nat))
  allocate(rtau_cryst(3,nsym,nat))
  !
  nr=(/nr1,nr2,nr3/)
  !
  call rot_atoms(nat,nsym,tau)
  !
  !================================================================================
  ! Allocate dynamical matrix and phonon energy array
  !================================================================================
  !
  allocate(dyn_q(3*nat,3*nat),w2(3*nat))
  !
  !================================================================================
  ! Symmetry realtions between irreducible q directory file number and the full mesh.
  !================================================================================
  !
  allocate(QE_folder_sym_q(nqmesh))
  allocate(sym_G_q(3,nqmesh))
  allocate(symlink_q(nqmesh,2))
  allocate(q_irr_cryst(3,nqirr))
  allocate(qstar(3,48))
  !
  do iq=1,nqirr
     q_irr_cryst(:,iq)=matmul(ainv(bg),q_irr(:,iq))
  enddo
  !
  !================================================================================
  ! Read de dv induced potentials from QE directories.
  !================================================================================
  !
  call read_allq_dvr(nqirr,3*nat)
  !
  !================================================================================
  ! We need the symmtry relations between irreducible/full in order to obtain
  ! the induced potential for any q point.
  !================================================================================
  !
  call set_symmetry_relations(nq1,nq2,nq3,nqirr,q_irr_cryst,qmesh, &
              q_points_consistent,QE_folder_sym_q,sym_G_q,symlink_q)
  !
  !================================================================================
  ! testing ...dv and mat_ep
  !================================================================================
  !
  call allocate_and_get_all_irreducible_wfc()
  !
  !================================================================================
  !      Calculate the small group of symm. operations for the irreducible k points
  !================================================================================
  !
  allocate(nsym_sgk(nk_irr))
  allocate(sindex_sgk(nk_irr,nsym))
  allocate(unit_sym_sgk(nkmesh,nsym,num_bands,num_bands))
  allocate(umat(num_bands,num_bands))
  !
  unit_sym_sgk=cmplx_0
  !
  do ik=1,nkmesh
     do ibnd=1,num_bands
        !
        unit_sym_sgk(ik,1:nsym,ibnd,ibnd)=cmplx_1
        !
     enddo !ibnd
  enddo !ik
  !
  call allocate_nlpot1
  !
  call allocate_phq
  !
  call init_pp
  
  ep_mat_el = cmplx_0
  !
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
  !
  if (calc_epmat) then
     !
     do iq=1,nqirr
        !
        if (                iq <   10) write(iq_loc,"(i1)")iq
        if ( 10 <= iq .and. iq <  100) write(iq_loc,"(i2)")iq
        if (100 <= iq .and. iq < 1000) write(iq_loc,"(i3)")iq
        !
        ep_unit=find_free_unit()
        inquire(iolength=record_lengh) ep_mat_el
        open(unit=ep_unit,file=trim(trim(mesh_dir)//trim(ep_mat_file)//trim('_')//adjustl(iq_loc)),iostat=ierr,&
        form='unformatted',status='unknown',access='direct',recl=record_lengh)
        !
        qpoint(:)=qmesh(:,iq)
        write(*,"(a,i4,100f12.6)")"qpoint",iq, qpoint
        !
        !-matrize dinamikoa kalkulatu q puntu jakin batentzat.
        !
        call mat_inv_four_t(qpoint,nq1,nq2,nq3,3*nat,frc,dyn_q)
        !
        !diagonalizatu goian kalkulatutako matrize dinamikoa: polarizazio bektore eta energiak lortu ditugu ia.
        !
        call diagonalize_cmat (3*nat, dyn_q,w2)
        !
        !-q puntuaren menpeko aldagai bat allocatu (hau kenduko dugu laster).
        !
        call allocate_nlpot2 (qpoint)
        !
        !-potentzialaren alde induzitua kalkulatu simetria erabiliz (errotazioz beharrezkoa izanez).
        !
        dvq_local=cmplx_0
        call get_dv(iq,qpoint,3*nat,dvq_local)
        !
        !-alde induzituari (goian), KB pseudopotentzialaren(pp) deribatuaren ALDE LOKALA gehitu.
        !
        call phq_init(matmul(bg,qpoint))
        call calculate_local_part_dv(qpoint,nat,npol,dvq_local)
        !
        !-bi subroutina hauek (goikoak), biak batera joan behar dira beti).
        !
        do ik=1,nkmesh
           !
           write(*,'(a,i4,a,3(f15.8))') "ik= ", ik,' k= ',matmul(bg,kmesh(:,ik))
           kpoint=kmesh(:,ik)
           !
           kpoint_cart=matmul(bg,kpoint)
           kqpoint_cart=matmul(bg,kpoint+qpoint)
           !
           !-uhina lortu RAM memorian dauden uhin irreduzibleak errotatuta. (k+q)
           !
           call get_psi_general_k_all_wfc(.true.,kpoint       ,list_iGk ,wfc_k ,QE_eig_k ,G_plusk)
           call get_psi_general_k_all_wfc(.true.,kpoint+qpoint,list_iGkq,wfc_kq,QE_eig_kq,G_pluskq)
           !
           !-k eta k+q puntuen indizeak topatu. ikpt_1 ik-aldagaiaren berdina izan
           ! behar da, baina hurrengo eran segurtasuna irabazten dugu.
           ! k-ren indizea: ikpt_k
           !
           call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
           call switch_indices(nk1,nk2,nk3,ikpt_k,i,j,k,+1)
           !
           ! kq-ren indizea: ikpt_kq
           !
           call find_k_1BZ_and_G(kpoint+qpoint,nk1,nk2,nk3,i,j,k,kqpoint_in_1bz,GKQ_bz)
           call switch_indices(nk1,nk2,nk3,ikpt_kq,i,j,k,+1)
           !
           ! Ordenan jartzen ditugu G bektoreak, k=0-n nola dauden ordenatuta arabera
           !
           nG=0
           do iG=1,nG_max
              if (list_iGk(iG)==0) exit
              nG=nG+1
           enddo
           npw=nG
           nG=0
           do iG=1,nG_max
              if (list_iGkq(iG)==0) exit
              nG=nG+1
           enddo
           npwq=nG
           !
           !-hemen KB potentzial ez lokaleko |beta> funtzioak kalkulatzen dira (k eta k+q puntuetarako hurrenez hurren).
           !
           vkb =cmplx_0
           vkqb=cmplx_0
           !
           call init_KB_projectors(npw ,nG_max,list_iGk ,kpoint       ,vkb )
           call init_KB_projectors(npwq,nG_max,list_iGkq,kpoint+qpoint,vkqb)
           !
           !-psi_k uhinak, potentzial induzitua + KB pp-ren ALDE LOKALAREN
           ! batuketarekin biderkatu (output-a:dvpsi): dv_local x |psi_k> (G)
           !
           call dvqpsi_local (3*nat, nG_max, nbands_loc, npol, list_iGk, list_iGkq, wfc_k, dvq_local, &
                                                    dvpsi(1:nG_max,1:nbands_loc,1:npol,1:npol,1:3*nat))
           !
           !-psi_k uhinak KB potentzialaren alde ez lokalarekin biderkatu eta emaitza dvpsi aldagaiari gehitu:
           !                    dvpsi^q_k --> dvpsi^q_k + D^q_mode [ KB ] |psi_k> (G)
           !                                  (lokala) + (ez lokala)
           call multiply_psi_by_KB( kpoint, qpoint, npol, nbands, nG_max, list_iGk, list_iGkq, wfc_k, dvpsi)
           !
           !-QE-ren subroutina goikoaren berdina egiteko.
           !
           do imode=1,3*nat ! Osagai kanonikoak, ez dira moduak, kontuz.
              !matrize elementuak kalkulatu
              !
              do jpol=1,nspin
                 do ipol=1,nspin
                    do jbnd=1,num_bands
                       do ibnd=1,num_bands
                          !
                          aep_mat_el(iq,ikpt_k,ibnd,jbnd,ipol,jpol,imode) = zdotc( nG_max, wfc_kq(:,ibnd,ipol), 1, dvpsi(:,jbnd,ipol,jpol,imode), 1 )
                          !
                       enddo !ibnd
                    enddo !jbnd
                 enddo !is
              enddo !js
              !
           enddo !imode
           !
        enddo !ik
        !
        write(unit=ep_unit,  rec = 1, iostat = ierr)  aep_mat_el      (iq,:,:,:,:,:,:)
        !
        close(unit=ep_unit)
        !
     enddo !iq

     else

        do iq=1,1
           qpoint(:)=qmesh(:,iq)
           !
           if (                iq <   10) write(iq_loc,"(i1)")iq
           if ( 10 <= iq .and. iq <  100) write(iq_loc,"(i2)")iq
           if (100 <= iq .and. iq < 1000) write(iq_loc,"(i3)")iq
           !
           ep_unit=find_free_unit()
           inquire(iolength=record_lengh) ep_mat_el
           open(unit=ep_unit,&
           file=trim(trim(mesh_dir)//trim('./')//trim(ep_mat_file)//trim('_')//adjustl(iq_loc)),&
           iostat=ierr,form='unformatted',status='unknown',access='direct',recl=record_lengh)
           !
           write(*,*)'Reading Reciprocal Bloch matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh,ierr
           read(unit=ep_unit,rec=1,iostat=ierr) aep_mat_el(iq,:,:,:,:,:,:)
           !
           close(unit=ep_unit)
           !
        enddo !iq

  end if ! calc_epmat

  deallocate (ep_mat_el)
  deallocate (aep_mat_el)

  deallocate(nsym_sgk)
  deallocate(sindex_sgk)
  deallocate(kpoints_irr )

  !call deallocate_upfeak ()
  !================================================================================
  !       Finish
  !================================================================================
  call get_timing(time2)
  write(*,20) '|                     ALL DONE                       |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='


20 format(A)
30 format(A,F8.2,6X,A)

contains

  subroutine dvqpsi_local(nmode,nG_max,nbands,npol,list_iGk,list_iGkq,wfc_k,dvq_local,dvpsi_local)

    implicit none

    !I/O variables

    integer,intent(in) :: nmode,nG_max,nbands,npol,list_iGk(nG_max),list_iGkq(nG_max)
    complex(kind=dp),intent(in) :: dvq_local(nr1*nr2*nr3,nmode,npol,npol),wfc_k(nG_max,nbands,npol)
    complex(kind=dp),intent(inout) :: dvpsi_local(nG_max,nbands,npol,npol,nmode)

    !local variables

    integer :: ibnd,ipol,ig,imode,ir
    complex(kind=dp) :: wfc_r(nr1*nr2*nr3,npol,npol),wfc_r1(nr1*nr2*nr3,npol)

    dvpsi_local(:,:,:,:,:)=cmplx_0
    !
    do imode =1,nmode !3*nat
      do ibnd=1,nbands
          !
          wfc_r1(:,:)=cmplx_0
          wfc_r(:,:,:)=cmplx_0
          !
          do ipol=1,npol !elektroi uhin funtzioak espazio errelarela pasatu. SO kasuan ipol osagaiz osagai.
             do ig=1,nG_max
                !
                if (list_iGk(ig)==0) exit
                wfc_r1(nl(list_iGk(ig)),ipol)=wfc_k(ig,ibnd,ipol)
                !
             enddo !ig
             !
             call cfftnd(3,(/nr1,nr2,nr3/),1,wfc_r1(:,ipol))
             !
          enddo !ipol

          !
!          if ((npol==2).and.spinorb_mag) then
          if (npol==2) then
             !
             do ir=1,nr1*nr2*nr3
                !
                do ipol=1,npol
                   do jpol=1,npol
                   !
                   wfc_r(ir,ipol,jpol)=dvq_local(ir,imode,ipol,jpol)*wfc_r1(ir,jpol)
                   !
                   enddo !jpol
                enddo !ipol
                !
             enddo !ir
             !
          else !npol
             !
             do ir=1,nr1*nr2*nr3
                !
                wfc_r(ir,1,1)=dvq_local(ir,imode,1,1)*wfc_r1(ir,1)
                !
             enddo !ir
             !
          endif !npol
          !
          do ipol=1,npol
             do jpol=1,npol
                !
                call cfftnd(3,(/nr1,nr2,nr3/),-1,wfc_r(:,ipol,jpol))
                !
                do ig=1,nG_max
                   !
                   if (list_iGkq(ig)==0) exit
                   !
!                   dvpsi_local(ig,ibnd,ipol,jpol,imode)=wfc_r(nl(ig),ipol,jpol)
                   dvpsi_local(ig,ibnd,ipol,jpol,imode)=wfc_r(nl(list_iGkq(ig)),ipol,jpol)
                   !
                enddo !ig
             enddo !jpol
          enddo !ipol
          !
       enddo !ibnd
    enddo !imode
    !
    return

  end subroutine dvqpsi_local

  !----------------------------------------

  subroutine diagonalize_cmat (n,a,w)

    integer, intent(in)  :: n
    complex(dp),intent(inout) :: a(n,n)
    real(dp),intent(out) :: w(n)

    complex(dp) :: a_pack(n*(n+1)/2)

    integer :: i,j, nfound

    complex(dp) :: cwork(2*n)
    real   (dp) :: rwork(7*n)
    integer     :: iwork(5*n), ifail(n), info


    do j=1,n
       do i=1,j
          a_pack(i+((j-1)*j)/2)=a(i,j)
       enddo
    enddo

    ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
    call ZHPEVX('V', 'A', 'U', n, a_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp,  &
         nfound, w, a, n,  cwork,  rwork,       iwork,          ifail, info)

  end subroutine diagonalize_cmat


end program ep_melements
