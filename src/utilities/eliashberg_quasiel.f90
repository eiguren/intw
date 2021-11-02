program eliashberg_quasiel

  use intw_input_parameters
  use intw_reading
  use intw_pseudo
  use intw_useful_constants
  use intw_intw2wannier
  use intw_symmetries
  use intw_fft
  use intw_band_crossing
  use intw_ph
  use intw_w90
  use w90_io, only: io_error,stdout,io_file_unit,seedname,io_time,io_stopwatch
  use w90_parameters, only: num_wann,bands_num_points,real_metric,recip_metric, &
                                            bands_num_spec_points,timing_level, &
                               bands_spec_points,bands_label,bands_plot_format, &
                   bands_plot_mode,num_bands_project,bands_plot_project,num_bands
  use w90_hamiltonian, only: irvec,nrpts,ndegen,ham_r
  use intw_allwfcs
  use intw_uspp, only: nkb, vkb, vkqb, nhtol, nhtolm, indv

  !================================================================================
  !       Declare the variables 
  !================================================================================
  implicit none

  !k point related variables
  logical                  :: k_points_consistent
  integer                  :: nk_irr , ik, nkmesh
  real(dp),allocatable     :: kmesh(:,:), kpoints_QE(:,:)
  integer                  :: ikpt, iqpt, ikqpt, ikpt_1, ikpt_2, ikpt_k, ikpt_kq
  real(dp),allocatable     :: kpoints_irr  (:,:)
  real(dp)                 :: kpoint(3), kpoint_cart(3), kqpoint(3), kqpoint_in_1bz(3)
  real(dp)                 :: kpoint_in_1bz(3), kpoint_rot(3)

  !q point related variables
  real(dp)                 :: qpoint(3), qpoint_1bz(3), qpoint_rot(3)
  integer                  :: ikpt_q

  !phonon related variables
  character(len=256)       :: fc_file_name
  real(dp)                 :: at_frc(3,3)
  integer                  :: mode, q_index_irr, q_index, imode, jmode
  complex(dp), allocatable :: dyn_q(:,:)
  real(dp), allocatable    :: w2(:)
  real(dp), allocatable    :: qstar(:,:)

  !symmetry variables
  integer                  :: isym, jsym, s_index, imq, i_sym, j_sym
  integer, allocatable     :: nsym_sgk(:), sindex_sgk(:,:)
  complex(dp), allocatable :: unit_sym_sgk( :,:,:,:), umat(:,:)

  !time related variables   
  real(dp)                 :: time, time1, time2

  !q point related variables
  integer                  :: iq, iqq, nq1_, nq2_, nq3_

  !path for bands plot
  integer                  :: number_of_special_points, nk_vec_path
  real(dp),allocatable     :: list_k_path(:,:)
  real(dp),allocatable     :: list_x_path(:)
  real(dp),allocatable     :: list_xticks(:)

  !ep interaction related variables
  integer                  :: ep_unit, Wep_unit
  logical                  :: q_points_consistent
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:,:,:), ep_mat_el(:,:,:,:,:,:)
  complex(dp), allocatable :: aW_rot_ep_mat_el(:,:,:,:,:,:,:), W_rot_ep_mat_el(:,:,:,:,:,:)
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
  integer                  :: ngk, ngkq, nr(3)
  integer                  :: G(1:3), GKQ_bz(3), G_plusk(3), G_pluskq(3)

  !ep fitxa
  integer                  :: record_lengh, ierr
  integer                  :: i_folder

  !Wigner-Seitz zones (for el (k) and ph (q))
  integer :: nrr_k, nrr_q
  integer,allocatable :: irr_k(:,:),irr_q(:,:)
  integer,allocatable :: ndegen_k(:), ndegen_q(:)

  !kmesh + qmesh finak
  integer :: nkf1,nkf2,nkf3,nkfmesh,nkf_irr,ik_irr
  integer,allocatable :: list_ikf_irr(:)
  integer :: nqf1,nqf2,nqf3,nqfmesh
  real(dp),allocatable :: kfmesh(:,:),kfpoints_irr(:,:),weightf_irr(:)
  real(dp),allocatable :: qfmesh(:,:),w_ik(:)

  !matrix-elements
  integer :: g_realw_unit
  real(dp) :: g2
  complex(dp),allocatable :: g_ReRp(:,:,:,:,:,:,:),ag_ReRp(:,:,:,:,:,:)
  complex(dp),allocatable :: g_kq(:,:,:,:,:,:,:),g_Req(:,:,:,:,:,:,:)
  complex(dp),allocatable :: g_Rp(:,:,:,:,:,:),g_wannier(:,:,:,:,:)
  complex(dp),allocatable :: g_elph(:,:,:,:,:)

  !Reciprocal Bloch to Real Wannier
  logical :: kqb_2_ReRpw_me,kqb_2_kqw_me,h_kb_2_Rew,dm_q_2_Rp

  !selfenergy + eliashberg
  integer :: npts,wpts,ipts
  real(dp) :: temp ! temperature of the system
  real(dp) :: trshld,w_min,w_max,d_w,dkvec(3)
  real(dp) :: ekq,ekk,delta_energ
  real(dp),allocatable :: eliash_f_so(:,:,:),eliash_f(:)

  !Hamiltonian (Electrons)
  logical :: write_el
  integer :: h_realw_unit,el_fine_unit,eigen_fine_unit
  real(dp),allocatable :: eigen(:,:),eigen_w(:)
  complex(dp),allocatable :: ham_rw(:,:,:),a_h_rw(:,:),U_el(:,:,:),U_el_w(:,:)

  !Dynamical Matrix (Phonons)
  logical :: write_ph,iwant
  integer :: ph_real_unit,ph_fine_unit,omega_fine_unit
  real(dp),allocatable :: omega(:,:),a_omega_f(:),omega_phon(:)
  complex(dp),allocatable :: dyn_r(:,:,:),a_ph_r(:,:),a_ph_q_f(:,:),U_ph(:,:,:),U_phon(:,:)

  !Filer of q-ak
  integer :: iq_cont,iq_2,nqfmesh_2
  integer,allocatable :: list_iq(:,:)

  !local/aux variables
  real(dp)                 :: sigma_w,term,sigma_e
  integer                  :: nbands_loc,irk,irq
  integer                  :: i,j,k
  integer                  :: ii,jj,kk
  integer                  :: jr,kr
  integer                  :: ig,ibnd,jbnd,ipol,jpol,iibnd,jjbnd,ifs,nkk,nfs
  integer                  :: switch
  logical                  :: read_status,have_nnkp,confirm
  character(256)           :: nnkp_file,method
  character(len=4)         :: iq_loc,ifs_loc
  character(1) :: ik_loc1
  character(2) :: ik_loc2
  character(3) :: ik_loc3
  character(4) :: ik_loc4


  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(*,*) '====================================================='
  write(*,*) '|                  program me                       |'
  write(*,*) '|        ---------------------------------          |'
  write(*,*) '====================================================='
  write(*,*) '|    waiting for input file...                      |'
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
  call get_ngm()
  allocate (gvec(3,ngm)) ! it is deallocated in the bottom of this main code.
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
  if (method == 'CONVOLUTION') then
     !
     write(*,*) '|       - intw2W_method   = CONVOLUTION             |'
     !
  elseif (method == 'FFT') then
     !
     write(*,*) '|       - intw2W_method   = FFT                     |'
     !
  else
     !
     write(*,*) '***********************************************'
     write(*,*) '* UNKNOWN COMPUTATION METHOD:'
     write(*,*) '* Only "CONVOLUTION" and "FFT" available'
     write(*,*) '***********************************************'
     !
     stop
     !
  endif
  !
  write(*,*) '|           ---------------------------------       |'
  !
  !================================================================================
  !       Check that $prefix.nnkp is present
  !================================================================================
  !
  nnkp_file = trim(mesh_dir)//trim(prefix)//".nnkp"
  !
  inquire(file=nnkp_file,exist=have_nnkp)
  !
  if (.not.have_nnkp) then
     !
     write(*,*)      '**********************************************************'
     write(*,*)      '* Could not find the file '//trim(nnkp_file)
     write(*,*)      '* Did you run W90 -pp $seed to get the parameter file?   '
     write(*,*)      '**********************************************************'
     !
     stop
     !
  endif
  !
  write(*,*) '|       - .nnkp file found                          |'
  write(*,*) '|           ---------------------------------       |'
  !
  !================================================================================
  !       read the parameters in the .nnkp file 
  !================================================================================
  !
  call read_nnkp_file(nnkp_file)
  !
  ! just as a test; can be removed later
  !
  call output_nnkp_file()
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
  !================================================================================
  ! check that kmesh and nnkp_kpoints are consistent with one another 
  !       This insures that the Wannier data is consistent with intw data.
  !================================================================================
  !
  call intw2W90_check_mesh(nkmesh,kmesh)
  !
  write(*,*) '|       - The mesh in the Wannier90 input.win file  |'
  write(*,*) '|         and the intw mesh are equal.              |'
  write(*,*) '|           ---------------------------------       |'
  !
  if (nspin==1) then
     !
     write(*,*) '|       - The calculation is paramagnetic nspin=1   |'
     write(*,*) '|                                                   |'
     write(*,*) '|           ---------------------------------       |'
     !
  elseif (nspin==2) then
     !
     write(*,*) '|       - Spin calculation nspin = 2                |'
     !
     if (noncolin) then
        !
        write(*,*) '|         Non-collinear Spin calculation            |'
        !
     endif
     !
     write(*,*) '|           ---------------------------------       |'
     !
  else
     !
     write(*,*) '*****************************************************'
     write(*,*) '* ERROR: Allowed values for nspin are 1 or 2        *'
     write(*,*) '*            program stops.                         *'
     write(*,*) '*****************************************************'
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
  call find_the_irreducible_k_set(nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
  !
  if (nkpoints_QE/=nkmesh .and. nkpoints_QE/=nk_irr) then
     !
     ! the points in the folder are not consistent with the
     ! input file. Stop the program!
     !
     write(*,*) '*****************************************************'
     write(*,*) '*      The number of kpoints present in the QE      *'
     write(*,*) '*      folders are not consistent with a full       *'
     write(*,*) '*      Brillouin Zone or an irreducible Brillouin   *'
     write(*,*) '*      zone! Review your input...                   *'
     write(*,*) '*                   Program stops.                  *'
     write(*,*) '*****************************************************'
     write(*,*) '* debug information:                                *'
     write(*,*) '*        nkpoints_QE = ',nkpoints_QE
     write(*,*) '*        nkmesh      = ',nkmesh
     write(*,*) '*        nk_irr      = ',nk_irr, nsym, tr_symmetry
     !
     stop
     !
  elseif (nkpoints_QE == nk_irr) then
     !
     ! The points in the QE folders *could* be a valid choice for the IBZ;
     ! this must be checked!
     !
     full_mesh = .false.
     IBZ       = .true.
     !
  else if (nkpoints_QE == nkmesh ) then
     !
     ! The points in the QE folders *could* be consistent with a full mesh;
     ! this must be checked!
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
  call allocate_and_build_spin_symmetry_matrices(nsym)
  !
  !================================================================================
  !      Fill the symmetry arrays appropriately
  !================================================================================
  call set_symmetry_relations(nk1, nk2, nk3, nkpoints_QE, kpoints_QE,kmesh, k_points_consistent &
       , QE_folder_sym, sym_G, symlink)
  !
  !================================================================================
  !      Calculate the multiplication talble for symmetry operations 
  !================================================================================
  call multable(nsym, s, symtable)
  !
  !================================================================================
  !       Tell the user what is in the QE folders
  !================================================================================
  if (full_mesh .and. IBZ) then
     write(*,*) '|       - the kpoints present in the QE folders     |'
     write(*,*) '|         are consistent with a full 1BZ and a      |'
     write(*,*) '|         IBZ has also been found.                  |'
     write(*,*) '|           ---------------------------------       |'
  else if(IBZ) then
     write(*,*) '|       - the kpoints present in the QE folders     |'
     write(*,*) '|         are consistent with an IBZ.               |'
     write(*,*) '|           ---------------------------------       |'
  else
     write(*,*)      '**********************************************************'
     write(*,*)      '* The kpoints present in the QE folders are not consistent'
     write(*,*)      '* with the parameters of the input file!                 '
     write(*,*)      '**********************************************************'
     write(*,*) '* debug information:                                *'
     write(*,*) '*        nkpoints_QE = ',nkpoints_QE
     write(*,*) '*        nkmesh      = ',nkmesh
     write(*,*) '*        nk_irr      = ',nk_irr, nsym, tr_symmetry
     write(*,*) '*        IBZ         = ',IBZ
     stop
  end if
  !
  !================================================================================
  !       Check that the requested calculation is possible
  !================================================================================
  if (intw2W_fullzone ) then
     write(*,*) '|       - intw2W_fullzone = .true.                  |'
     if ( full_mesh ) then
        write(*,*) '|         all k-points are explicitely calculated   |'
        write(*,*) '|         no symmetry is assumed.                   |'
        write(*,*) '|             (This is mostly for testing)          |'
        write(*,*) '|           ---------------------------------       |'
     else
        write(*,*)      '**********************************************************'
        write(*,*)      '* A full mesh is not present in the QE folders!          '
        write(*,*)      '* The requested calculation is impossible.               '
        write(*,*)      '*                   program stops.                       '
        write(*,*)      '**********************************************************'
        stop
     end if
  else
     write(*,*) '|       - intw2W_fullzone = .false.                 |'
     write(*,*) '|         Symmetries will be utilized.              |'
     write(*,*) '|           ---------------------------------       |'
  end if

  write(*,*) '|       - reading pseudopotentials from UPF files   |'
  write(*,*) '|            (defined in .save data files)          |'
  write(*,*) '|                                                   |'

  call read_all_pseudo ()
  if (.not.lspinorb) call average_pp(ntyp)

  write(*,*) '|                    PPs are OK                     |'
  write(*,*) '|           ---------------------------------       |'

  if (.not.conmesurate_and_coarser (nk1, nk2, nk3, nq1, nq2, nq3)) then
     write(*,*)      '**********************************************************'
     write(*,*)       '*ERROR                                                  '
     write(*,*)       '*   the electron k and phonon q are not                 '
     write(*,*)       '*   conmesurate and the k grid does not contain         '
     write(*,*)       '*   the phonon q grid                                   '
     write(*,*)      '**********************************************************'
     stop
  end if
  !
  !================================================================================
  !       set up W90 data and generate the u_mesh array 
  !================================================================================
  !       CAREFUL!! It is CRUCIAL to allocate and read the W90 stuff BEFORE
  !       referencing it. num_wann, for example, exists but is not properly
  !       DEFINED before.
  call allocate_and_read_W90()
  !
  !       NOW you can allocate
  !
  ! extract the nband x num_wann Wannier projection+rotation matrices.
  call produce_u_mesh()
  !
  !================================================================================
  !       Get the special points from the file
  !================================================================================
  if (bands_points /= 0 ) then
     call get_interpolation_special_points( number_of_special_points)
  end if
  !
  !================================================================================
  !       Allocate wfc related files
  !================================================================================
  nbands_loc=num_bands
  !
  allocate (list_igk (nG_max))
  allocate (list_igkq(nG_max))
  allocate (list_igk_aux(nG_max))
  allocate (list_igk_orig(nG_max))
  !
  allocate (wfc_k(nG_max,nbands_loc,nspin))
  allocate (wfc_kq(nG_max,nbands_loc,nspin))
  allocate (wfc_k_aux(nG_max,nbands_loc,nspin))
  allocate (wfc_k_orig(nG_max,nbands_loc,nspin))
  !
  allocate (QE_eig_k (nbands_loc))
  allocate (QE_eig_kq(nbands_loc))
  !
  allocate (fr (nr1*nr2*nr3,nbands_loc,nspin))
  allocate (fg (nG_max,nbands_loc,nspin))
  !
  !================================================================================
  !       Read all the information about phonons and pseudopotentials
  !================================================================================
  call read_ph_information_xml()
  !
  !================================================================================
  !       Read the force constant matrix from the QE directory
  !================================================================================
  fc_file_name=trim(trim(mesh_dir)//"/"//trim(ph_dir)//"/"//trim(fc_mat))
  !
  call readfc(fc_file_name,nq1_,nq2_,nq3_,nat,alat,at_frc,ntyp,amass)
  !
  !================================================================================
  ! Build the phonon qmesh corresponding to the parameters in the input file 
  !================================================================================
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  !
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !
  allocate (aep_mat_el(nqmesh,nk1*nk2*nk3,nbands_loc,nbands_loc,nspin,nspin,3*nat))  ! q bakoitzerako guztia
  allocate (aW_rot_ep_mat_el(nqmesh,nk1*nk2*nk3,num_wann,num_wann,nspin,nspin,3*nat))
  allocate (W_rot_ep_mat_el(nk1*nk2*nk3,num_wann,num_wann,nspin,nspin,3*nat))
  allocate (ep_mat_el(nk1*nk2*nk3,nbands_loc,nbands_loc,nspin,nspin,3*nat))
  !
  !================================================================================
  ! Below, the table of rotations for each atom and symmetry.
  !================================================================================
  allocate ( rtau_index(nat, nsym))
  allocate ( rtau(3,nsym,nat) )
  allocate ( tau_cryst(3,nat) )
  allocate ( rtau_cryst(3,nsym,nat) )
  !
  nr=(/nr1, nr2, nr3/)
  !
  call rot_atoms(nat,nsym, tau)
  !
  !================================================================================
  ! Allocate dynamical matrix and phonon energy array
  !================================================================================
  allocate(dyn_q(3*nat,3*nat), w2(3*nat))
  !
  !================================================================================
  ! Symmetry realtions between irreducible q directory file number and the full mesh.
  !================================================================================
  allocate(QE_folder_sym_q(nqmesh))
  allocate(sym_G_q      (3,nqmesh))
  allocate(symlink_q    (nqmesh,2))
  allocate(q_irr_cryst   (3,nqirr))
  allocate(qstar(3,48)            )
  !
  do iq=1,nqirr
     q_irr_cryst(:,iq) = matmul(ainv(bg), q_irr(:,iq))
  enddo
  !
  !================================================================================
  ! We need the symmtry relations between irreducible/full in order to obtain
  ! the induced potential for any q point.
  !================================================================================
  call set_symmetry_relations(nq1, nq2, nq3, nqirr, q_irr_cryst ,qmesh,       &
       q_points_consistent, QE_folder_sym_q, sym_G_q, symlink_q)
  !
  !================================================================================
  ! testing ...dv and mat_ep
  !================================================================================
  call allocate_and_get_all_irreducible_wfc()
  !
  !================================================================================
  !      Calculate the small group of symm. operations for the irreducible k points
  !================================================================================
  !
  call allocate_nlpot1
  !
  call allocate_phq
  !
  call init_us_1
  !
  !================================================================================
  ! Starting the el-ph calculation
  !================================================================================
  kqb_2_ReRpw_me=.false.
  kqb_2_kqw_me=.false.
  h_kb_2_Rew=.false.
  dm_q_2_Rp=.false.
  !
  !================================================================================
  ! Determination of the Wigner-Seitz points 
  !================================================================================
  nrr_k=0
  nrr_q=0
  !
  call count_wigner_seitz_points(nk1,nk2,nk3,nrr_k)
  call count_wigner_seitz_points(nq1,nq2,nq3,nrr_q)
  !
  allocate(ndegen_k(nrr_k),ndegen_q(nrr_q))
  allocate(irr_k(3,nrr_k),irr_q(3,nrr_q))
  !
  ndegen_k(:)=0.0d0
  ndegen_q(:)=0.0d0
  !
  call wigner_seitz_points(nk1,nk2,nk3,nrr_k,irr_k,ndegen_k)
  call wigner_seitz_points(nq1,nq2,nq3,nrr_q,irr_q,ndegen_q)
  !
  allocate(g_ReRp(nrr_q,nrr_k,num_wann,num_wann,nspin,nspin,3*nat))
  g_ReRp=cmplx_0
  !
  allocate(ag_ReRp(nrr_k,num_wann,num_wann,nspin,nspin,3*nat)) ! For writing
  !
  !===================================================================================================
  ! Reciprocal Bloch coarse to Real Wannier (Wannier Rotation + Fourier Transform): Matrix Elements
  !===================================================================================================
  if (kqb_2_ReRpw_me) then
     !
     ! If we have not yet calculated the matrix elements in the Reciprocal Wannier Representation, we do it
     !
     !=============================================================================================
     ! We read the Reciprocal Wannier matrix elements in the coarse mesh
     ! If not calculated, we transform (rotate) the Reciprocal Bloch matrix elements with
     ! Wannier Rotation 
     !=============================================================================================
     allocate(g_kq(nqmesh,nkmesh,num_wann,num_wann,nspin,nspin,3*nat))
     g_kq=cmplx_0
     !
     if (kqb_2_kqw_me) then
        !
        ! We read the matrix elements in Reciprocal Bloch Representation
        ! (and canonical-cartesian representation for phonons)
        !
        do iq=1,nqmesh
           qpoint(:)=qmesh(:,iq)        
           !
           if (iq<10)        write(iq_loc,"(i1)")iq
           if (10<=iq<100)   write(iq_loc,"(i2)")iq
           if (100<=iq<1000) write(iq_loc,"(i3)")iq        
           !
           ep_unit=find_free_unit()
           inquire(iolength=record_lengh) ep_mat_el
           open(unit=ep_unit,&
           file=trim(trim(mesh_dir)//trim('EP_DATA/')//trim(ep_mat_file)//trim('_')//adjustl(iq_loc)),&
           iostat=ierr,form='unformatted',status='unknown',access='direct',recl=record_lengh)
           !
           write(*,*)'Reading Reciprocal Bloch matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh,ierr
           read(unit=ep_unit,rec=1,iostat=ierr) aep_mat_el(iq,:,:,:,:,:,:)
           !
           close(unit=ep_unit)
           !
           Wep_unit=find_free_unit()
           inquire(iolength=record_lengh) W_rot_ep_mat_el
           open (unit=Wep_unit,&
           file=trim(trim(mesh_dir)//trim(ep_mat_file)//trim("W")//trim('_')//adjustl(iq_loc)),&
           iostat=ierr,form='unformatted',status='unknown',access='direct',recl=record_lengh)
           !
           write(*,*)'Calculating Reciprocal Wannier matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh,ierr
           !
           do ik=1,nkmesh
              kpoint(:)=kmesh(:,ik)
              !
              call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
              call switch_indices(nk1,nk2,nk3,ikpt_k,i,j,k,+1)
              call find_k_1BZ_and_G(kpoint+qpoint,nk1,nk2,nk3,i,j,k,kqpoint_in_1bz,GKQ_bz)
              call switch_indices(nk1,nk2,nk3,ikpt_kq,i,j,k,+1)
              !
              ! Now, we rotate (smooth) the matrix elements using the U Wannier unitary matrix
              !
              do imode=1,3*nat
                 call Wrotate_plane_wave_matrix_elements(ikpt_kq,ikpt_k,u_mesh,&
                 aep_mat_el(iq,ikpt_k,1:nbands_loc,1:nbands_loc,1:nspin,1:nspin,imode),&
                 aW_rot_ep_mat_el(iq,ikpt_k,1:num_wann,1:num_wann,1:nspin,1:nspin,imode))
              enddo !imode
           enddo !ik
           !
           write(*,*)'Writing Reciprocal Wannier matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh
           write(unit=Wep_unit,rec=1,iostat=ierr) aW_rot_ep_mat_el(iq,:,:,:,:,:,:)
           close(unit=Wep_unit)
           !
           write(*,*)' '
        enddo !iq
        !
        write(*,*)'****************************************************************'
!Peio
!do iq=30,30
!do ik=240,240
!do imode=1,3*nat
!   write(7000,"(a,i4,a,i4)") ' ik= ',ik,' imode= ',imode
!   write(7500,"(a,i4,a,i4)") ' ik= ',ik,' imode= ',imode
!   term=cmplx_0
!   do ipol=1,npol
!      do jpol=1,npol
!         write(7000,"(a,i4,a,i4)")' ipol= ',ipol,' jpol= ',jpol
!         write(7500,"(a,i4,a,i4)")' ipol= ',ipol,' jpol= ',jpol
!         do ibnd=1,nbands_loc
!            write(7000,"(100f12.6)") (abs(aep_mat_el(iq,ikpt_k,ibnd,jbnd,ipol,jpol,imode)),jbnd=1,nbands_loc)
!            write(7500,"(100(a,f12.6,a,f12.6,a))") (' ',real(aep_mat_el(iq,ikpt_k,ibnd,jbnd,ipol,jpol,imode)),&
!                        " + ",aimag(aep_mat_el(iq,ikpt_k,ibnd,jbnd,ipol,jpol,imode)),"*I,", jbnd=1,nbands_loc)
!            write(7000,*)"------------------------------------------------"
!            write(7500,*)"------------------------------------------------"
!            do jbnd=1,nbands_loc
!               term=term+aep_mat_el(iq,ikpt_k,ibnd,jbnd,ipol,jpol,imode)
!            enddo !jbnd
!         enddo !ibnd
!      enddo !jpol
!   enddo !ipol
!   write(333,*) abs(term)
!enddo !imode
!enddo
!enddo
!stop
!Peio
     else
        !
        ! In this case, the matrix elements are already smooth, we directly read them
        ! The matrix elements are already in Reciprocal Wannier Representation
        !
        do iq=1,nqmesh
           if (iq<10)        write(iq_loc,"(i1)")iq
           if (10<=iq<100)   write(iq_loc,"(i2)")iq
           if (100<=iq<1000) write(iq_loc,"(i3)")iq
           !
           Wep_unit=find_free_unit()
           inquire(iolength=record_lengh) W_rot_ep_mat_el
           open(unit=Wep_unit,&
           file=trim(trim(mesh_dir)//trim('EPW_DATA/')//trim(ep_mat_file)//trim("W")//trim('_')//adjustl(iq_loc)),&
           iostat=ierr,form='unformatted',status='unknown',access='direct',recl=record_lengh)
           !
           write(*,*)'Reading Reciprocal Wannier matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh
           read(unit=Wep_unit,rec=1,iostat=ierr)  aW_rot_ep_mat_el(iq,:,:,:,:,:,:)
           close(unit=Wep_unit)
           !
           write(*,*)' '
        enddo
        !
        write(*,*)'****************************************************************'
     endif
     !
     ! We change of name variable, for comodity
     !
     g_kq(:,:,:,:,:,:,:)=aW_rot_ep_mat_el(:,:,:,:,:,:,:)
     !
     deallocate(aep_mat_el,ep_mat_el,aW_rot_ep_mat_el,W_rot_ep_mat_el)
     !     
     !==============================================================================================
     ! Reciprocal coarse Wannier Representation to Real Wannier Representation: g(k,q) -> g(Re,Rp)
     !==============================================================================================
     !
     ! We perform this transformation in two steps.
     ! 1) g(k,q) -> g(Re,q)
     ! 2) g(Re,q) -> g(Re,Rp)
     !
     allocate(g_Req(nqmesh,nrr_k,num_wann,num_wann,nspin,nspin,3*nat))
     g_Req=cmplx_0
     !
     write(*,*)'------------------------------------------------------------------'
     write(*,*)'------------------------------------------------------------------'
     write(*,*)'Transformation from Reciprocal Wannier on the coarse mesh to Real Wannier'
     write(*,*)'First step: g(k,q) -> g(Re,q)'
     !
     call G_kqw2Reqw(nrr_k,irr_k,kmesh,g_kq,g_Req)
     !
     deallocate(g_kq)
     !
     write(*,*)'------------------------------------------------------------------'
     write(*,*)' '
     write(*,*)'Second step: g(Re,q) -> g(Re,Rp)'
     !
     call G_Reqw2ReRpw(nrr_q,nrr_k,irr_q,irr_k,qmesh,g_Req,g_ReRp)
     !
     deallocate(g_Req)
     !
     write(*,*)'------------------------------------------------------------------'
     write(*,*)' '
     !
     !================================================================================
     ! Writing the matrix elements in Real Wannier Representation
     !================================================================================
     g_realw_unit=find_free_unit()
     inquire(iolength=record_lengh) ag_ReRp
     open(unit=g_realw_unit,file=trim(trim(mesh_dir)//trim(g_realw_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     do irq=1,nrr_q
        write(*,*)'Writing Real Wannier matrix elements irq = ',irq,'/',nrr_q
        write(unit=g_realw_unit,rec=irq,iostat=ierr) g_ReRp(irq,:,:,:,:,:,:)
     enddo
     !
     close(unit=g_realw_unit)
  else
  !
  ! If we have already calculated matrix elements in Real Wannier Representation, we read them  
  !
  !================================================================================
  ! Reading the matrix elements in Real Wannier Representation
  !================================================================================
     g_realw_unit=find_free_unit()
     inquire(iolength=record_lengh) ag_ReRp
     open(unit=g_realw_unit,file=trim(trim(mesh_dir)//trim(g_realw_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     write(*,*)'------------------------------------------------------------------'
     write(*,*)'------------------------------------------------------------------'
     !
     do irq=1,nrr_q
        write(*,*)'Reading Real Wannier matrix elements irq = ',irq,'/',nrr_q
        read(unit=g_realw_unit,rec=irq,iostat=ierr) g_ReRp(irq,:,:,:,:,:,:)
     enddo
     !
     close(unit=g_realw_unit)
     !
     write(*,*)'------------------------------------------------------------------'
  endif
  !
  !====================================================================================
  ! Reciprocal Bloch coarse to Real Wannier Transform: 1) electrons (Wannier + Fourier)
  ! Reciprocal coarse to Real Transform: 2) phonons (Fourier)
  !====================================================================================
  write(*,*)' '
  write(*,*)'------------------------------------------------------------------'
  write(*,*)'------------------------------------------------------------------'
  write(*,*)'Reciprocal coarse to Real transform for 1) Electrons and 2) Phonons'
  write(*,*)' '
  !
  ! 1) Electrons
  !
  write(*,*)'******************************************************************'
  write(*,*)'Reciprocal Bloch to Real Wannier: 1) Electrons'
  write(*,*)' '
  !
  allocate(ham_rw(num_wann,num_wann,nrr_k))
  allocate(a_h_rw(num_wann,num_wann))
  ham_rw(:,:,:)=cmplx_0
  !
  if (h_kb_2_Rew) then
     write(*,*)'We must calculate the Hamiltonians in the Real Wannier Representation'
     !
     call H_kb2Rew(nrr_k,irr_k,kmesh,ham_rw)
     !
     write(*,*)' '
     !
     !===============================================================================
     ! Writing the Hamiltoninas in the Real Wannier Representation
     !===============================================================================    
     h_realw_unit=find_free_unit()
     inquire(iolength=record_lengh) a_h_rw
     open(unit=h_realw_unit,file=trim(trim(mesh_dir)//trim(h_realw_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     do irk=1,nrr_k
        write(*,*)'Writing Real Wannier Hamiltonians irk = ',irk,'/',nrr_k
        write(unit=h_realw_unit,rec=irk,iostat=ierr) ham_rw(:,:,irk)
     enddo
     !
     close(unit=h_realw_unit)
  else
  !
  ! If we have already calculated the Hamiltonians in Real Wannier Representation, we read them
  !
  !=================================================================================
  ! Reading the Hamiltonians in Real Wannier Representation
  !=================================================================================
     h_realw_unit=find_free_unit()
     inquire(iolength=record_lengh) a_h_rw
     open(unit=h_realw_unit,file=trim(trim(mesh_dir)//trim(h_realw_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     do irk=1,nrr_k
        write(*,*)'Reading Real Wannier Hamiltonians irk = ',irk,'/',nrr_k
        read(unit=h_realw_unit,rec=irk,iostat=ierr) ham_rw(:,:,irk)
     enddo
     !
     close(unit=h_realw_unit)
  endif
  deallocate(a_h_rw)
  !
  write(*,*)' '
  !
  ! 2) Phonons
  !
  write(*,*)'******************************************************************'
  write(*,*)'2) Phonons'
  write(*,*)' '
  !
  allocate(dyn_r(3*nat,3*nat,nrr_q))
  allocate(a_ph_r(3*nat,3*nat))
  dyn_r(:,:,:)=cmplx_0
  !
  if (dm_q_2_Rp) then
     write(*,*)'We must calculate the Real Dynamical Matrices'
     !
     call DM_q2Rp(nrr_q,irr_q,dyn_r)
     !
     write(*,*)' '
     !
     !=============================================================================
     ! Writing the Dynamical Matrices in the Real representation
     !=============================================================================
     ph_real_unit=find_free_unit()
     inquire(iolength=record_lengh) a_ph_r
     open(unit=ph_real_unit,file=trim(trim(mesh_dir)//trim(ph_real_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     do irq=1,nrr_q
        write(*,*)'Writing Real Dynamical Matrices irq = ',irq,'/',nrr_q
        write(unit=ph_real_unit,rec=irq,iostat=ierr) dyn_r(:,:,irq)
     enddo
     !
     close(unit=ph_real_unit)
  else
  !
  ! If we have already calculated the Dynamical Matrices in Real representation, we read them
  !
  !================================================================================
  ! Reading the Dynamical Matrices in Real representation
  !================================================================================
     ph_real_unit=find_free_unit()
     inquire(iolength=record_lengh) a_ph_r
     open(unit=ph_real_unit,file=trim(trim(mesh_dir)//trim(ph_real_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     do irq=1,nrr_q
        write(*,*)'Reading Real Dynamical Matrices irq = ',irq,'/',nrr_q
        read(unit=ph_real_unit,rec=irq,iostat=ierr) dyn_r(:,:,irq)
     enddo
     !
     close(unit=ph_real_unit)
  endif
  deallocate(a_ph_r)
  !
  write(*,*)' '
  !
  !=========================================================================================
  ! Definition of the electrons and phonons fine grids (Hamiltonians + energies for each k,
  ! Dynamical Matrix + phonon modes for each q,...)
  !
  ! Note: these variables (fine el mesh [nkf1,nkf2,nkf3] + ph mesh [nfq1,nqf2,nqf3]) 
  ! would have to be called from the input files, by the moment we defne by hand here
  !=========================================================================================
  !
  ! We define the electron fine grid
  !
  write(*,*)'Definition of the electron fine grid...'
  !
  nkf1=1000
  nkf2=1000
  nkf3=1
  nkfmesh=nkf1*nkf2*nkf3
  allocate(kfmesh(3,nkfmesh))   
  call generate_kmesh(kfmesh,nkf1,nkf2,nkf3)
  !
  ! We calculate and diagonalize the Hamiltonian of all the k points of the fine grid
  ! We stock them in the memory, in that manner we won't waste time calculating + diagonalizing
  ! them once and again every time we meet a new k point
  ! This operation takes some time, but we will thank it aftrwards
  !
  allocate(U_el_w(num_wann,num_wann))
  allocate(eigen_w(num_wann))
  allocate(U_el(num_wann,num_wann,nkfmesh))
  allocate(eigen(num_wann,nkfmesh))
  U_el(:,:,:)=cmplx_0
  eigen(:,:)=0.0d0
  write_el=.false.
  !
  if (write_el) then
     !
     do ik=1,nkfmesh
        !
        if ((ik/10000)*10000==ik) then
           write(*,*) 'ik = ',ik,'/',nkfmesh
        endif
        !
        kpoint(:)=kfmesh(:,ik)
        call H_Rew2kb(nrr_k,irr_k,ndegen_k,kpoint,ham_rw,U_el(:,:,ik),eigen(:,ik))
        !
     enddo !ik
     !
     ! We write it, for future calculations
     !
     el_fine_unit=find_free_unit()
     inquire(iolength=record_lengh) U_el_w
     open(unit=el_fine_unit,file=trim(trim(mesh_dir)//trim(el_fine_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     eigen_fine_unit=find_free_unit()
     inquire(iolength=record_lengh) eigen_w
     open(unit=eigen_fine_unit,file=trim(trim(mesh_dir)//trim(eigen_fine_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     do ik=1,nkfmesh
        !
        if ((ik/10000)*10000==ik) then
           write(*,*)'Writing eigenvectors & eigenvalues ik = ',ik,'/',nkfmesh
        endif
        !
        write(unit=el_fine_unit,rec=ik,iostat=ierr) U_el(:,:,ik)
        write(unit=eigen_fine_unit,rec=ik,iostat=ierr) eigen(:,ik)
        !
     enddo !ik
     !
     close(unit=el_fine_unit)
     close(unit=eigen_fine_unit)
     !
  else
     !
     ! We have already calculated, we read them
     !
     el_fine_unit=find_free_unit()
     inquire(iolength=record_lengh) U_el_w
     open(unit=el_fine_unit,file=trim(trim(mesh_dir)//trim(el_fine_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     eigen_fine_unit=find_free_unit()
     inquire(iolength=record_lengh) eigen_w
     open(unit=eigen_fine_unit,file=trim(trim(mesh_dir)//trim(eigen_fine_file)),iostat=ierr,&
     form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     do ik=1,nkfmesh
        !
        if ((ik/10000)*10000==ik) then
           write(*,*)'Reading eigenvectors & eigenvalues ik = ',ik,'/',nkfmesh
        endif
        !
        read(unit=el_fine_unit,rec=ik,iostat=ierr) U_el(:,:,ik)
        read(unit=eigen_fine_unit,rec=ik,iostat=ierr) eigen(:,ik)
        !
     enddo !ik
     !
     close(unit=el_fine_unit)
     close(unit=eigen_fine_unit)
     !
  endif
  !
!THIS IS FOR AVERAGE PROGRAM
!  ! We calculate the number of irreducible k points in the k fine grid and its respective weights and
!  ! associate in the all k points fine grid list
!  ! Useful for caltulating if summing in the k variable. Instead of summing all the k points, just the
!  ! irreducibles -> less time of computing! Can't be done for summing in the q variable.
!  !
!  call find_size_of_irreducible_k_set(nkf1,nkf2,nkf3,kfmesh,nkf_irr)
!  allocate(kfpoints_irr(3,nkf_irr))
!  allocate(weightf_irr(nkf_irr))
!  allocate(list_ikf_irr(nkf_irr))
!  call find_the_irreducible_k_set_irr(nkf1,nkf2,nkf3,kfmesh,nkf_irr,list_ikf_irr,kfpoints_irr,weightf_irr)
!THIS IS FOR AVERAGE PROGRAM
  !
  ! We define the phonon fine grid
  !
  write(*,*)'Definition of the phonon fine grid...'
  !
  nqf1=1000
  nqf2=1000
  nqf3=1
  nqfmesh=nqf1*nqf2*nqf3
  allocate(qfmesh(3,nqfmesh))
  call generate_kmesh(qfmesh,nqf1,nqf2,nqf3)
  !
  ! We can calculate and diagonalize now the DM in order to save time for the next. That is a good
  ! stratedy, specially regarding average calculations (very long) (iwant=.true.)
  !
  ! If yes, we calculate and diagonalize the Dynamical Matrix of all the q points of the fine grid
  ! Two options: we calculate it now (write_ph=.true.) or we calculated it in a previous calculation 
  ! and we read them (write_ph=.false.) -> adecuate in the input file
  ! In that manner, we stock in the memory the DM and the omega and we won't calculate them, fact that
  ! requires a lot of time
  !
  allocate(U_phon(3*nat,3*nat))
  allocate(omega_phon(3*nat))
  iwant=.true.
  write_ph=.false.
  !
  if (iwant) then
     !
     allocate(U_ph(3*nat,3*nat,nqfmesh))
     allocate(omega(3*nat,nqfmesh))
     U_ph(:,:,:)=0.0d0
     omega(:,:)=0.0d0
     !
     allocate(a_ph_q_f(3*nat,3*nat))
     allocate(a_omega_f(3*nat))
     !
     if (write_ph) then
        !
        ! We calculate in the whole BZ
        ! 
        do iq=1,nqfmesh
           !
           if ((iq/1000)*1000==iq) then
              write(*,*)'iq = ',iq,' / ',nqfmesh
           endif
           !
           qpoint(:)=qfmesh(:,iq)
           call DM_Rp2q(nrr_q,irr_q,ndegen_q,qpoint,dyn_r,U_ph(:,:,iq),omega(:,iq))
           !
        enddo
        !
        ! We write them, for future calculations
        !
        ph_fine_unit=find_free_unit()
        inquire(iolength=record_lengh) a_ph_q_f
        open(unit=ph_fine_unit,file=trim(trim(mesh_dir)//trim(ph_fine_file)),iostat=ierr,&
        form='unformatted',status='unknown',access='direct',recl=record_lengh)
        !
        omega_fine_unit=find_free_unit()
        inquire(iolength=record_lengh) a_omega_f
        open(unit=omega_fine_unit,file=trim(trim(mesh_dir)//trim(omega_fine_file)),iostat=ierr,&
        form='unformatted',status='unknown',access='direct',recl=record_lengh)
        !
        do iq=1,nqfmesh
           !
           if ((iq/1000)*1000==iq) then
              write(*,*)'Writing Dynamical Matrices & phonons iq = ',iq,'/',nqfmesh
           endif
           !
           write(unit=ph_fine_unit,rec=iq,iostat=ierr) U_ph(:,:,iq)
           write(unit=omega_fine_unit,rec=iq,iostat=ierr) omega(:,iq)
           !
        enddo !iq
        !
        close(unit=ph_fine_unit)
        close(unit=omega_fine_unit)
        !
     else
        !
        ! We have already calculated, so we read it
        !
        ph_fine_unit=find_free_unit()
        inquire(iolength=record_lengh) a_ph_q_f
        open(unit=ph_fine_unit,file=trim(trim(mesh_dir)//trim(ph_fine_file)),iostat=ierr,&
        form='unformatted',status='unknown',access='direct',recl=record_lengh)
        !
        omega_fine_unit=find_free_unit()
        inquire(iolength=record_lengh) a_omega_f
        open(unit=omega_fine_unit,file=trim(trim(mesh_dir)//trim(omega_fine_file)),iostat=ierr,&
        form='unformatted',status='unknown',access='direct',recl=record_lengh)
        !
        do iq=1,nqfmesh
           !
           if ((iq/1000)*1000==iq) then
              write(*,*)'Reading Dynamical Matrices & phonons iq = ',iq,'/',nqfmesh
           endif
           !
           read(unit=ph_fine_unit,rec=iq,iostat=ierr) U_ph(:,:,iq)
           read(unit=omega_fine_unit,rec=iq,iostat=ierr) omega(:,iq)
           !
        enddo !iq
        !
        close(unit=ph_fine_unit)
        close(unit=omega_fine_unit)
        !
     endif !write_ph
     !
     deallocate(a_ph_q_f)
     deallocate(a_omega_f)
     !
  else
     !
     allocate(U_ph(3*nat,3*nat,1))
     allocate(omega(3*nat,1))
     !
  endif !iwant
  !
  !=========================================================================================
  ! Allocation of needed variables in oder to make the transformation:
  !
  ! g(Re,Rp) Real Wannier -> g(k_f,q_f) Reciprocal Bloch
  !=========================================================================================
  allocate(g_Rp(nrr_q,num_wann,num_wann,nspin,nspin,3*nat))
  allocate(g_wannier(num_wann,num_wann,nspin,nspin,3*nat))
  allocate(g_elph(num_wann,num_wann,nspin,nspin,3*nat))
  !
  !==========================================================================
  ! Set-up the environment of our el-ph calculation
  !==========================================================================
  !
  temp=0.0d0 ! Temperature of the system
  !
  if (lspinorb) then
     trshld=-2.88073485d0
  else
     trshld=-2.60126050d0
  endif
  !
  w_min=0.0d0
  w_max=0.070d0
  wpts=3501 ! how many points for our w variable in eliashberg function a^2F(w)
  !
  if (lspinorb) allocate(eliash_f_so(wpts,npol,npol))
  !
  allocate(eliash_f(wpts))
  !
  d_w=abs(w_max-w_min)/(wpts-1)
  sigma_w=10.d0*d_w
  sigma_e=0.0005
  !
  ! We choose in which band and direction we will make the calculation and how
  ! many points will be in our calculation:
  !
  jbnd=2
  nkk=200
  !
  dkvec(1)=0.001d0
  dkvec(2)=0.001d0
  dkvec(3)=0.d0
  !
!  dkvec(1)=0.001d0
!  dkvec(2)=-0.002d0
!  dkvec(3)=0.d0
  !
  !==============================================================================================
  ! Calculation of the state dependent Eliashberg Quasi Elastic Function in fonction the
  ! (k,n) states through a given direction in the BZ (normaly a special direction)
  ! It depends on the state we want to research
  !
  ! Here we perform the transformation from Real Wannier to Reciprocal Bloch respresentation
  ! of the matrix elements, so that: g(Re,Rp) -> g(k',q')
  !==============================================================================================
  !
  ! We know:
  !
  ! a^2F(k,n,si,sj,w)=Sum_q_m_nu[|g(q,k,m,n,si,sj,nu)|^2*Delta[w-w(q,nu)]*Delta[E(k,n)-E(k+q,m)]]
  !
  ! So, we can see that in this calculation a lot of k+q points will contribute negligibly, because 
  ! states (k+q,m) had eigenvalues very far away from the Fermi Level (that is E(k,n)) (Quasi Elastic
  ! Approximation)
  ! We simulate the deltas with gaussian functions, but gaussian tails "die" very fast from the center
  ! of the peak:
  !
  ! Delta[x-xo]=Exp[-(x-xo)^2/sigma^2]/(sigma*Sqrt[Pi])
  !
  ! When we move 4*sigma away from the center we have a contribution of 10^-7 normalized to the value
  ! of the peak. When we move 5*sigma we have a contribution of 10^-11. So for states k+q points with
  ! eigenvalues as: |E(k,n)-E(k+q,m)|>4*sigma, we don't take them into account, and then we ease up
  ! very much the calculation
  !
  allocate(list_iq(nqfmesh*num_wann,2))
  !
  do ik=1,nkk
     !
     write(*,*) 'ik = ',ik,'/',nkk
     !
     ! We choose the kpoint in the direction we chosed before and we identify it within the
     ! kfmesh we defined before.
     !
     kpoint(:)=0.0d0+dkvec(:)*(ik-1)
     !
!     if (ik.lt.101) then
!        jbnd=2
!        dkvec(1)=0.001d0
!        dkvec(2)=0.001d0
!        dkvec(3)=0.d0
!        kpoint(1)=0.234d0+dkvec(1)*(ik-1)
!        kpoint(2)=0.234d0+dkvec(2)*(ik-1)
!        kpoint(3)=0.0d0
!     else
!        jbnd=1
!        dkvec(1)=0.001d0
!        dkvec(2)=-0.002d0
!        dkvec(3)=0.d0
!        kpoint(1)=0.333d0+dkvec(1)*(ik-100)
!        kpoint(2)=0.333d0+dkvec(2)*(ik-100)
!        kpoint(3)=0.0d0
!     endif
     !
     call find_k_1BZ_and_G(kpoint,nkf1,nkf2,nkf3,i,j,k,kpoint_in_1bz,GKQ_bz)
     call switch_indices(nkf1,nkf2,nkf3,ikpt_k,i,j,k,+1)
     kpoint(:)=kfmesh(:,ikpt_k)
     kpoint_cart=matmul(bg,kpoint)
     !
     ! We define the eigenvalues of the state (k,n)
     !
     ekk=eigen(jbnd,ikpt_k)
     !
     write(1314999,'(i10,3f8.5,f12.6,100f12.4)') ik,kpoint(:), &
           sqrt(kpoint_cart(1)**2+kpoint_cart(2)**2+kpoint_cart(3)**2),eigen(:,ikpt_k)+trshld
     !
     ! We initialize the eliashberg function for this state
     !
     if (lspinorb) eliash_f_so(:,:,:)=0.0d0
     !
     eliash_f(:)=0.0d0
     !
     ! We make the transformation: g(Re,Rp) -> g(k,Rp)
     !
     g_Rp(:,:,:,:,:,:)=cmplx_0
     call G_ReRpw2kRpw(nrr_k,nrr_q,irr_k,ndegen_k,kpoint,g_ReRp,g_Rp)
     !
     ! We start the summ with respect to q (qmesh) and to m (ibnd). We limit this addition,
     ! we just take into account the states (q,m) we don't misprise
     !
     iq_cont=0
     list_iq(:,:)=0
     !
     do iq=1,nqfmesh
        !
        if ((iq/10000)*10000==iq) then
           write(*,*) 'iq = ',iq,'/',nqfmesh,' '
        endif 
        !
        qpoint(:)=qfmesh(:,iq)
        call find_k_1BZ_and_G(kpoint+qpoint,nkf1,nkf2,nkf3,i,j,k,kqpoint_in_1bz,GKQ_bz)
        call switch_indices(nkf1,nkf2,nkf3,ikpt_kq,i,j,k,+1)
        !
        do ibnd=1,num_wann
           !
           ! We define the eigenvalues of (q,m)
           !
           ekq=eigen(ibnd,ikpt_kq)
           !
           ! We create the list of terms we are interested in
           !
           if (abs(ekk-ekq).gt.4.d0*sigma_e) then
              !
           else
              !
              iq_cont=iq_cont+1
              list_iq(iq_cont,1)=iq
              list_iq(iq_cont,2)=ibnd
              !
           endif !abs(ekk-ekq)
           !
        enddo !ibnd
     enddo !iq
     !
     nqfmesh_2=iq_cont !how many states (q,m) we take into account in the addition
     !
     do iq_2=1,nqfmesh_2
        !
        if ((iq_2/500)*500==iq_2) then
           write(*,*) 'iq = ',iq_2,'/',nqfmesh_2,' ik = ',ik,'/',nkk
        endif
        !
        ! For each (q,m) term of the addition, we have (k+q,m) state. We calculate the
        ! energy of this state and his delta value: Delta[E(k,n)-E(k+q,m)]
        !
        iq=list_iq(iq_2,1)
        ibnd=list_iq(iq_2,2)
        !
        qpoint(:)=qfmesh(:,iq)
        call find_k_1BZ_and_G(kpoint+qpoint,nkf1,nkf2,nkf3,i,j,k,kqpoint_in_1bz,GKQ_bz)
        call switch_indices(nkf1,nkf2,nkf3,ikpt_kq,i,j,k,+1)
        ekq=eigen(ibnd,ikpt_kq)
        !
        if (iwant) then
           !
           U_phon(:,:)=U_ph(:,:,iq)
           omega_phon(:)=omega(:,iq)
           !
        else
           !
           call DM_Rp2q(nrr_q,irr_q,ndegen_q,qpoint,dyn_r,U_ph(:,:,1),omega(:,1))
           U_phon(:,:)=U_ph(:,:,1)
           omega_phon(:)=omega(:,1)
           !
        endif !iwant
        !
        ! We make the transformation: g(k,Rp) -> g(k,q)
        !
        g_wannier(:,:,:,:,:)=cmplx_0
        g_elph(:,:,:,:,:)=cmplx_0
        call G_kRpw2kqw(nrr_q,irr_q,ndegen_q,qpoint,g_Rp,g_wannier)
        call G_kqw2kqb(U_el(:,:,ikpt_kq),U_el(:,:,ikpt_k),U_phon(:,:),g_wannier,g_elph)
        !
        delta_energ=exp(-((ekk-ekq)**2.d0)/((sigma_e)**2.d0))/(sigma_e*sqrt(pi))
        !
        if (lspinorb) call add_q_m_nu_term_qe_so(w_min,w_max,wpts,sigma_w,delta_energ,omega_phon,g_elph(ibnd,jbnd,:,:,:),eliash_f_so)
        !
        call add_q_m_nu_term_qe(w_min,w_max,wpts,sigma_w,delta_energ,omega_phon,g_elph(ibnd,jbnd,:,:,:),eliash_f)
        !
     enddo !iq_2
     !
     ! We normalized the Eliashberg Function
     !
     if (lspinorb) eliash_f_so(:,:,:)=eliash_f_so(:,:,:)/nqfmesh
     !
     eliash_f(:)=eliash_f(:)/nqfmesh
     !
     ! We write our results, for next calculations (in .dat but also in .bin)
     !
     if (ik<10) write(ik_loc1,'(i1)')ik
     if (ik<100) write(ik_loc2,'(i2)')ik
     if (ik<1000) write(ik_loc3,'(i3)')ik
     if (ik<10000) write(ik_loc4,'(i4)')ik
     !
     nfs=nspin*nspin
     !
     if (lspinorb) then
        !
        do ifs=1,nfs
           !
           write(ifs_loc,'(i1)')ifs
           !
           if (ik<10) then
              open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc1)//trim('_')//trim(ifs_loc)//trim('.dat')))
              open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc1)//trim('_')//trim(ifs_loc)//trim('.bin')),form='unformatted')
           elseif (ik<100) then
              open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc2)//trim('_')//trim(ifs_loc)//trim('.dat')))
              open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc2)//trim('_')//trim(ifs_loc)//trim('.bin')),form='unformatted')
           elseif (ik<1000) then
              open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc3)//trim('_')//trim(ifs_loc)//trim('.dat')))
              open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc3)//trim('_')//trim(ifs_loc)//trim('.bin')),form='unformatted')
           elseif (ik<10000) then
              open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc4)//trim('_')//trim(ifs_loc)//trim('.dat')))
              open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_so_')//trim(ik_loc4)//trim('_')//trim(ifs_loc)//trim('.bin')),form='unformatted')
           endif
           !
           do ipts=1,wpts
              !
              if (ifs==1) then
                 !
                 write(10*ik,'(100f16.8)') w_min+abs(w_max-w_min)*(ipts-1)/(wpts-1),eliash_f_so(ipts,1,1)
                 write(10*ik+100000) eliash_f_so(ipts,1,1)
                 !
              elseif (ifs==2) then
                 !
                 write(10*ik,'(100f16.8)') w_min+abs(w_max-w_min)*(ipts-1)/(wpts-1),eliash_f_so(ipts,1,2)
                 write(10*ik+100000) eliash_f_so(ipts,1,2)
                 !
              elseif (ifs==3) then
                 !
                 write(10*ik,'(100f16.8)') w_min+abs(w_max-w_min)*(ipts-1)/(wpts-1),eliash_f_so(ipts,2,1)
                 write(10*ik+100000) eliash_f_so(ipts,2,1)
                 !
              elseif (ifs==4) then
                 !
                 write(10*ik,'(100f16.8)') w_min+abs(w_max-w_min)*(ipts-1)/(wpts-1),eliash_f_so(ipts,2,2)
                 write(10*ik+100000) eliash_f_so(ipts,2,2)
                 !
              endif !ifs
              !
           enddo !ipts
           !
           close(10*ik)
           close(10*ik+100000)
           !
        enddo !ifs
        !
     endif !lspinorb
     !
     if (ik<10) then
        open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc1)//trim('.dat')))
        open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc1)//trim('.bin')),form='unformatted')
     elseif (ik<100) then
        open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc2)//trim('.dat')))
        open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc2)//trim('.bin')),form='unformatted')
     elseif (ik<1000) then
        open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc3)//trim('.dat')))
        open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc3)//trim('.bin')),form='unformatted')
     elseif (ik<10000) then
        open(unit=10*ik,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc4)//trim('.dat')))
        open(unit=10*ik+100000,file=trim(trim(mesh_dir)//trim('eliash_')//trim(ik_loc4)//trim('.bin')),form='unformatted')
     endif
     !
     do ipts=1,wpts
        !
        write(10*ik,'(100f16.8)') w_min+abs(w_max-w_min)*(ipts-1)/(wpts-1),eliash_f(ipts)
        write(10*ik+100000) eliash_f(ipts)
        !
     enddo !ipts
     !
     close(10*ik)
     close(10*ik+100000)
     !
  enddo !ik
  deallocate(list_iq)

end program eliashberg_quasiel
