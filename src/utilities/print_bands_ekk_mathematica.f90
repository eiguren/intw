program eliashberg_full

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
  real(dp)                 :: qpoint(3), qpoint_1bz(3), qpoint_rot(3), qpoint_cart(3)
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
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:,:,:),ep_mat_el(:,:,:,:,:,:)
  complex(dp), allocatable :: aW_rot_ep_mat_el(:,:,:,:,:),aep_mat_el_2(:,:,:,:,:),W_rot_ep_mat_el(:,:,:,:)
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
  complex(dp),allocatable :: g_ReRp(:,:,:,:,:),ag_ReRp(:,:,:,:)
  complex(dp),allocatable :: g_kq(:,:,:,:,:),g_Req(:,:,:,:,:)
  complex(dp),allocatable :: g_Rp(:,:,:,:),g_wannier(:,:,:)
  complex(dp),allocatable :: g_elph(:,:,:),g_me(:)

  !Reciprocal Bloch to Real Wannier
  logical :: kqb_2_ReRpw_me,kqb_2_kqw_me,h_kb_2_Rew,dm_q_2_Rp

  !selfenergy + eliashberg
  integer :: npts,wpts,ipts
  real(dp) :: temp ! temperature of the system
  real(dp) :: trshld,w_min,w_max,d_w,dkvec(3)
  real(dp) :: ekq,ekk,delta_eplus,delta_eminus
  real(dp),allocatable :: eliash_fplus(:),eliash_fminus(:)

  !Hamiltonian (Electrons)
  logical :: write_el
  integer :: h_realw_unit,el_fine_unit,eigen_fine_unit
  real(dp),allocatable :: eigen(:,:),eigen_w(:),eigen_bands(:,:,:)
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
  real(dp)                 :: sigma_w,term,sigma_e,modk,modkq,sintita,costita,k1(3),k2(3),k3(3),k4(3)
  real(dp) :: k1mod,k2mod,k3mod,k4mod
  integer                  :: nbands_loc,irk,irq,iq_w_min_K,iq_w_min_M,iq_trshld_K,iq_trshld_M,iqx,iqy
  integer                  :: i,j,k,ikx,iky
  integer                  :: ii,jj,kk
  integer                  :: jr,kr
  integer                  :: ig,ibnd,jbnd,is,js,iibnd,jjbnd,ifs,nkk,nfs
  integer                  :: switch
  logical                  :: read_status,have_nnkp,confirm,deg_k,deg_kq
  logical,allocatable :: deg_iq(:)
  character(256)           :: nnkp_file,method
  character(len=4)         :: iq_loc,ifs_loc
  character(1) :: ibnd_loc


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
  else if (method == 'FFT') then
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
  elseif (nkpoints_QE == nkmesh ) then
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
  h_kb_2_Rew=.false.
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
  write_el=.true.
  !
  if (write_el) then
     !
     !$omp parallel default(none) &
     !$omp shared(nkfmesh,kfmesh,nrr_k,irr_k,ndegen_k,ham_rw,U_el,eigen) &
     !$omp private(ik,kpoint) 
     !
     !$omp do
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
     !$omp end parallel
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
  if (lspinorb) then
     trshld=-2.57509105d0
  else
     trshld=-2.47036625d0
  endif
  !
  allocate(eigen_bands(-nkf1:nkf1,-nkf2:nkf2,num_wann))
  !
  do ik=1,nkfmesh
     !
     kpoint(:)=kfmesh(:,ik)
     call find_k_1BZ_and_G(kpoint,nkf1,nkf2,nkf3,i,j,k,kpoint_in_1bz,GKQ_bz)
     !
     eigen_bands((i-1)-nkf1,(j-1)-nkf2,:)=eigen(:,ik)+trshld
     eigen_bands((i-1)-nkf1,j-1,:)=eigen(:,ik)+trshld
     eigen_bands(i-1,(j-1)-nkf2,:)=eigen(:,ik)+trshld
     eigen_bands(i-1,j-1,:)=eigen(:,ik)+trshld
     !
  enddo !ik
  !
  do ibnd=2,2
     !
     write(ibnd_loc,'(i1)')ibnd
     !
     open(unit=11,file=trim(trim(mesh_dir)//trim('band_')//(trim(ibnd_loc))//trim('.dat')))
     !
     do ikx=-nkf1,nkf1
        do iky=-nkf2,nkf2
           !
           kpoint(1)=1.d0*ikx/nkf1
           kpoint(2)=1.d0*iky/nkf2
           !
           write(11,'(100f16.8)') kpoint(1),kpoint(2),eigen_bands(ikx,iky,ibnd)
           !
        enddo !iky
     enddo !ikx
     !
     close(11)
     !
  enddo !ibnd

end program eliashberg_full
