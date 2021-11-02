program imaginary_selfen

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
  use intw_utility

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
  real(dp)                 :: kpoint_in_1bz(3), kpoint_rot(3), qpoint_in_1bz(3)

  !q point related variables
  real(dp)                 :: qpoint(3), qpoint_1bz(3), qpoint_rot(3), qpoint_cart(3)
  integer                  :: ikpt_q

  !phonon related variables
  character(len=256)       :: fc_file_name
  real(dp)                 :: at_frc(3,3)
  integer                  :: mode, q_index_irr, q_index, imode, jmode
  complex(dp), allocatable :: dyn_q(:,:),dyn_q_coarse(:,:,:)
  real(dp), allocatable    :: w2(:),w22(:)
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
  integer                  :: ep_unit, Wep_unit, Wep_unit_ext, ep_unit_ext
  logical                  :: q_points_consistent
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:,:,:),aep_mat_el_ext(:,:,:,:,:,:,:),ep_mat_el(:,:,:,:,:,:)
  complex(dp), allocatable :: aep_mat_el_2(:,:,:,:,:),aep_mat_el_2_ext(:,:,:,:,:)
  complex(dp), allocatable :: aW_rot_ep_mat_el(:,:,:,:,:),aW_rot_ep_mat_el_ext(:,:,:,:,:),W_rot_ep_mat_el(:,:,:,:)
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
  integer :: nkf1,nkf2,nkf3,nkfmesh
  integer :: nqf1,nqf2,nqf3,nqfmesh
  real(dp),allocatable :: kfmesh(:,:)
  real(dp),allocatable :: qfmesh(:,:)

  !matrix-elements KS
  integer :: g_realw_unit
  real(dp) :: g2
  complex(dp),allocatable :: g_ReRp(:,:,:,:,:),ag_ReRp(:,:,:,:)
  complex(dp),allocatable :: g_kq(:,:,:,:,:),g_Req(:,:,:,:,:)
  complex(dp),allocatable :: g_Re(:,:,:,:),g_wannier(:,:,:)
  complex(dp),allocatable :: g_elph(:,:,:),g_me(:)

  !matrix-elements EXT
  real(dp) :: g2_ext
  complex(dp),allocatable :: g_ReRp_ext(:,:,:,:,:)
  complex(dp),allocatable :: g_kq_ext(:,:,:,:,:),g_Req_ext(:,:,:,:,:)
  complex(dp),allocatable :: g_Re_ext(:,:,:,:),g_wannier_ext(:,:,:)
  complex(dp),allocatable :: g_elph_ext(:,:,:),g_me_ext(:)

  !Reciprocal Bloch to Real Wannier
  logical :: kqb_2_ReRpw_me,kqb_2_kqw_me,kqb_2_ReRpw_me_ext,kqb_2_kqw_me_ext,h_kb_2_Rew,dm_q_2_Rp

  !selfenergy
  real(dp) :: temp ! temperature of the system
  real(dp) :: trshld,trshld1,trshld2,dqvec(3),wqv
  real(dp) :: ekq,ekk
  real(dp) :: ef
  real(dp),allocatable :: phonren(:,:)

  !Hamiltonian (Electrons)
  logical :: write_el
  integer :: h_realw_unit,el_fine_unit,eigen_fine_unit
  real(dp),allocatable :: eigen(:,:),eigen_w(:)
  complex(dp),allocatable :: ham_rw(:,:,:),a_h_rw(:,:),U_el(:,:,:),U_el_w(:,:)

  !Dynamical Matrix (Phonons)
  logical :: write_ph,iwant
  integer :: ph_real_unit,ph_fine_unit,omega_fine_unit
  real(dp),allocatable :: omega(:,:),a_omega_f(:),omega_phon(:)
  complex(dp),allocatable :: dyn_r(:,:,:),a_ph_r(:,:),a_ph_q_f(:,:),U_ph(:,:,:),U_phon(:,:),c_nodn(:,:),pol1(:),pol2(:)

  !local/aux variables
  real(dp)                 :: term,sigma_e,arg
  complex(dp) :: fac
  integer                  :: nbands_loc,irk,irq,na,naimode,najmode
  integer                  :: i,j,k,is,js
  integer                  :: ii,jj,kk
  integer                  :: jr,kr
  integer                  :: ig,ibnd,jbnd,ipol,jpol,iibnd,jjbnd,ifs,nkk,nqq
  integer                  :: switch
  logical                  :: read_status,have_nnkp,confirm
  character(256)           :: nnkp_file,method
  character(len=4)         :: iq_loc,ifs_loc
  character(1) :: ibnd_loc,jbnd_loc,imode_loc
  real(dp),parameter :: rytomev=13605.698066d0 , amutory=911.44424213227238468637d0

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
  allocate (aep_mat_el(nqmesh,nk1*nk2*nk3,nbands_loc,nbands_loc,nspin,nspin,3*nat))  ! q bakoitzerako guztia
  allocate (aep_mat_el_2(nqmesh,nk1*nk2*nk3,nbands_loc,nbands_loc,3*nat))
  allocate (aep_mat_el_ext(nqmesh,nk1*nk2*nk3,nbands_loc,nbands_loc,nspin,nspin,3*nat))  ! q bakoitzerako guztia
  allocate (aep_mat_el_2_ext(nqmesh,nk1*nk2*nk3,nbands_loc,nbands_loc,3*nat))
  allocate (ep_mat_el(nk1*nk2*nk3,nbands_loc,nbands_loc,nspin,nspin,3*nat))
  !
  aep_mat_el(:,:,:,:,:,:,:)=cmplx_0
  aep_mat_el_2(:,:,:,:,:)=cmplx_0
  aep_mat_el_ext(:,:,:,:,:,:,:)=cmplx_0
  aep_mat_el_2_ext(:,:,:,:,:)=cmplx_0
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
  allocate(dyn_q(3*nat,3*nat), w2(3*nat),w22(3*nat))
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
  ! Determination of the Wigner-Seitz points 
  !================================================================================
  nrr_q=0
  !
  call count_wigner_seitz_points(nq1,nq2,nq3,nrr_q)
  !
  allocate(ndegen_q(nrr_q))
  allocate(irr_q(3,nrr_q))
  !
  ndegen_q(:)=0.0d0
  !
  call wigner_seitz_points(nq1,nq2,nq3,nrr_q,irr_q,ndegen_q)
  !
  !=============================================================================================
  ! We read the Reciprocal Wannier matrix elements in the coarse mesh
  !=============================================================================================
  allocate(g_kq(nqmesh,nkmesh,num_wann,num_wann,3*nat))
  allocate(g_kq_ext(nqmesh,nkmesh,num_wann,num_wann,3*nat))
  g_kq=cmplx_0
  g_kq_ext=cmplx_0
  !
  allocate(g_me(3*nat),g_me_ext(3*nat))
  g_me=cmplx_0
  g_me_ext=cmplx_0
  !
  allocate(dyn_q_coarse(3*nat,3*nat,nqmesh))
  dyn_q_coarse=cmplx_0
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
     write(*,*)'Calculating Reciprocal Bloch NO SPIN matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh,ierr
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
           !
           do ibnd=1,num_bands
              do jbnd=1,num_bands
                 !
                 do is=1,nspin
                    do js=1,nspin
                       !
                       aep_mat_el_2(iq,ikpt_k,ibnd,jbnd,imode)=aep_mat_el_2(iq,ikpt_k,ibnd,jbnd,imode)+ &
                                                               aep_mat_el(iq,ikpt_k,ibnd,jbnd,is,js,imode)
                       !
                    enddo !js
                 enddo !is
                 !
              enddo !jbnd
           enddo !ibnd
           !
           if (iq.eq.1.and.imode.lt.4) aep_mat_el_2(iq,ikpt_k,:,:,imode)=cmplx_0
           !
        enddo !imode
     enddo !ik
     !
     ep_unit=find_free_unit()
     inquire(iolength=record_lengh) ep_mat_el
     open(unit=ep_unit,&
     file=trim(trim(mesh_dir)//trim('EP_DATA_EXT/')//trim(ep_mat_file)//trim('_')//adjustl(iq_loc)),&
     iostat=ierr,form='unformatted',status='unknown',access='direct',recl=record_lengh)
     !
     write(*,*)'Reading Reciprocal Bloch matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh,ierr
     read(unit=ep_unit,rec=1,iostat=ierr) aep_mat_el_ext(iq,:,:,:,:,:,:)
     !
     close(unit=ep_unit)
     !
     write(*,*)'Calculating Reciprocal Bloch NO SPIN matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh,ierr
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
           !
           do ibnd=1,num_bands
              do jbnd=1,num_bands
                 !
                 do is=1,nspin
                    do js=1,nspin
                       !
                       aep_mat_el_2_ext(iq,ikpt_k,ibnd,jbnd,imode)=aep_mat_el_2_ext(iq,ikpt_k,ibnd,jbnd,imode)+ &
                                                                   aep_mat_el_ext(iq,ikpt_k,ibnd,jbnd,is,js,imode)
                       !
                    enddo !js
                 enddo !is
                 !
              enddo !jbnd
           enddo !ibnd
           !
           if (iq.eq.1.and.imode.lt.4) aep_mat_el_2_ext(iq,ikpt_k,:,:,imode)=cmplx_0
           !
        enddo !imode
     enddo !ik
     !
     write(*,*)' '
     write(*,*)'Calculating Dynamical Matrix in the coarse mesh for iq = ',iq,'/',nqmesh
     !
     call mat_inv_four_t(qpoint,nq1,nq2,nq3,3*nat,frc,dyn_q)
     !
     dyn_q_coarse(:,:,iq)=dyn_q(:,:)
     !
  enddo !iq
  !
  write(*,*)'****************************************************************'
  !
  ! We change of name variable, for comodity
  !
  g_kq(:,:,:,:,:)=aep_mat_el_2(:,:,:,:,:)
  g_kq_ext(:,:,:,:,:)=aep_mat_el_2_ext(:,:,:,:,:)
  !
  deallocate(aep_mat_el,aep_mat_el_2)
  deallocate(aep_mat_el_ext,aep_mat_el_2_ext)
  !
  if (lspinorb) then
     trshld1=-2.57509105d0 !TlSi
     trshld2=-2.39934425d0
!     trshld=1.9594 !C graph 
  else
     trshld1=-2.57838555d0 !TlSi
     trshld2=-2.47036625d0
!     trshld=1.9594 !C graph
  endif
  !
  ef=0.350d0
  !
  ! We transform in to ARU the atom masses
  !
  amass(:)=amass(:)*amutory
  !
  ! Elementu guztiak kargatuta, coarseko matrix elements eta coarseko dynamical matrix
  ! Goazen orain dynamical matrix-en aldatuta kalkulatzera
  !
  write(*,*)''
  write(*,*)'Calculating changed Dynamical Matrix in the coarse mesh'
  !
  do iq=1,nqmesh
     !
     write(*,*)'Calculating changed Dynamical Matrix in the coarse mesh',iq,'/',nqmesh
     !
     qpoint(:)=qmesh(:,iq)
     !
     do ik=1,nkmesh
        !
        kpoint(:)=kmesh(:,ik)
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
        do ibnd=1,num_wann
           !
           if (ibnd.lt.2) ekq=QE_eig_kq(ibnd)+trshld1
           if (ibnd.gt.1) ekq=QE_eig_kq(ibnd)+trshld2
           ! 
           do jbnd=1,num_wann
              !
              if (jbnd.lt.2) ekk=QE_eig_k(jbnd)+trshld1
              if (jbnd.gt.1) ekk=QE_eig_k(jbnd)+trshld2
              !
              if (ekk.lt.ef.and.ekq.gt.ef) then
                 !
                 g_me(:)=g_kq(iq,ikpt_k,ibnd,jbnd,:)
                 g_me_ext(:)=g_kq_ext(iq,ikpt_k,ibnd,jbnd,:)
                 !
                 do imode=1,3*nat
                    !
                    naimode=((imode-1)/3)+1
                    !
                    do jmode=1,3*nat
                       !
                       najmode=((jmode-1)/3)+1
                       !
                       g2=2*real(g_me_ext(jmode)*conjg(g_me(imode)))*rytomev**2.d0
!                       g2=2*real(g_me(jmode)*conjg(g_me(imode)))
                       !
                       dyn_q_coarse(imode,jmode,iq)=dyn_q_coarse(imode,jmode,iq)+&
                                                    (1.d0/dble(nkmesh))*(1.d0/sqrt(amass(ityp(naimode))*amass(ityp(najmode))))*&
                                                    (g2/(ekk-ekq))
                       !
                    enddo !jmode
                    !
                 enddo !imode
                 !
              endif ! ekk<ef, ekq>ef
              !
              if (ekk.gt.ef.and.ekq.lt.ef) then
                 !
                 g_me(:)=g_kq(iq,ikpt_k,ibnd,jbnd,:)
                 g_me_ext(:)=g_kq_ext(iq,ikpt_k,ibnd,jbnd,:)
                 !
                 do imode=1,3*nat
                    !
                    naimode=((imode-1)/3)+1
                    !
                    do jmode=1,3*nat
                       !
                       najmode=((jmode-1)/3)+1
                       !
                       g2=2*real(g_me_ext(jmode)*conjg(g_me(imode)))*rytomev**2.d0
!                       g2=2*real(g_me(jmode)*conjg(g_me(imode)))
                       !
                       dyn_q_coarse(imode,jmode,iq)=dyn_q_coarse(imode,jmode,iq)-&
                                                    (1.d0/dble(nkmesh))*(1.d0/sqrt(amass(ityp(naimode))*amass(ityp(najmode))))*&
                                                    (g2/(ekk-ekq))
                       !
                    enddo !jmode
                    !
                 enddo !imode
                 !
              endif ! ekk>ef, ekq<ef
              !
           enddo !jbnd
           !
        enddo !ibnd
        !
        ekq=QE_eig_kq(2)
        ekk=QE_eig_k(1)
        !
        do imode=1,3*nat
           !
           naimode=((imode-1)/3)+1
           !
           do jmode=1,3*nat
              !
              najmode=((jmode-1)/3)+1
              !
              g2=2*real(g_kq_ext(iq,ikpt_k,2,1,jmode)*conjg(g_kq(iq,ikpt_k,2,1,imode)))*rytomev**2.d0
!              g2=2*real(g_kq(iq,ikpt_k,2,1,jmode)*conjg(g_kq(iq,ikpt_k,2,1,imode)))
!              g2=0.0d0
              !
              dyn_q_coarse(imode,jmode,iq)=dyn_q_coarse(imode,jmode,iq)-&
                                           (1.d0/dble(nkmesh))*(1.d0/sqrt(amass(ityp(naimode))*amass(ityp(najmode))))*&
                                           (g2/(ekk-ekq))
              !
           enddo !jmode
           !
        enddo !imode
        !
        ekq=QE_eig_kq(1)
        ekk=QE_eig_k(2)
        !
        do imode=1,3*nat
           !
           naimode=((imode-1)/3)+1
           !
           do jmode=1,3*nat
              !
              najmode=((jmode-1)/3)+1
              !
              g2=2*real(g_kq_ext(iq,ikpt_k,1,2,jmode)*conjg(g_kq(iq,ikpt_k,1,2,imode)))*rytomev**2.d0
!              g2=2*real(g_kq(iq,ikpt_k,1,2,jmode)*conjg(g_kq(iq,ikpt_k,1,2,imode)))
!              g2=0.0d0
              !
              dyn_q_coarse(imode,jmode,iq)=dyn_q_coarse(imode,jmode,iq)+&
                                           (1.d0/dble(nkmesh))*(1.d0/sqrt(amass(ityp(naimode))*amass(ityp(najmode))))*&
                                           (g2/(ekk-ekq))
              !
           enddo !jmode
           !
        enddo !imode
        !
     enddo !ik
     !
     dyn_q_coarse(:,:,iq)=0.5d0*(dyn_q_coarse(:,:,iq)+conjg(transpose(dyn_q_coarse(:,:,iq))))
     !
  enddo !iq
  !
  write(*,*)' '
  write(*,*)'------------------------------------------------------------------'
  write(*,*)'------------------------------------------------------------------'
  write(*,*)'Reciprocal coarse to Real transform for Phonons'
  write(*,*)' '
  !  
  allocate(dyn_r(3*nat,3*nat,nrr_q))
  allocate(a_ph_r(3*nat,3*nat))
  dyn_r(:,:,:)=cmplx_0
  !
  write(*,*)'Fourier transform: Reciprocal Dynamical Matrix -> Real Dynamical Matrix'
  !
  write(*,*)' '
  !
  do irq=1,nrr_q
     write(*,*)'irq = ',irq,'/',nrr_q
     !
     do iq=1,nqmesh
        !
        arg=tpi*dot_product(qmesh(:,iq),dble(irr_q(:,irq)))
        fac=exp(-cmplx_i*arg)/dble(nqmesh)
        dyn_r(:,:,irq)=dyn_r(:,:,irq)+fac*dyn_q_coarse(:,:,iq)
        !
     enddo !iq
     !
  enddo !irq
  !
  ! We interpolate the dynamical matrix to the path we want to plot
  !
  nqq=501
  !
  allocate(U_ph(3*nat,3*nat,nqq))
  allocate(omega(3*nat,nqq))
  U_ph(:,:,:)=cmplx_0
  omega(:,:)=0.d0  
  !
  ! We choose in which band and direction we will make the calculation and how
  ! many points will be in our calculation:
  !
  open(unit=1314,file=trim(trim(mesh_dir)//trim('path_phon_ren_DM.dat')))
  !
  do iq=1,nqq
     !
     ! We choose the qpoint in the direction we chosed before and we identify it within the
     ! qfmesh we defined before.
     !
     write(*,*) 'iq = ',iq,'/',nqq
     !
     if (iq.lt.168) then
        !
        dqvec(1)=0.002d0
        dqvec(2)=0.002d0
        dqvec(3)=0.000d0
        !
        qpoint(1)=0.000d0+dqvec(1)*(iq-1)
        qpoint(2)=0.000d0+dqvec(2)*(iq-1)
        qpoint(3)=0.000d0+dqvec(3)*(iq-1)
        !
     elseif (iq.gt.167.and.iq.lt.251)then
        !
        dqvec(1)=0.002d0
        dqvec(2)=0.004d0
        dqvec(3)=0.000d0
        !
        qpoint(1)=0.334d0+dqvec(1)*(iq-168)
        qpoint(2)=0.332d0-dqvec(2)*(iq-168)
        qpoint(3)=0.000d0+dqvec(3)*(iq-168)
        !
     elseif (iq.gt.250) then
        !
        dqvec(1)=0.002d0
        dqvec(2)=0.000d0
        dqvec(3)=0.000d0
        !
        qpoint(1)=0.500d0-dqvec(1)*(iq-251)
        qpoint(2)=0.000d0+dqvec(2)*(iq-251)
        qpoint(3)=0.000d0+dqvec(3)*(iq-251)
        !
     endif
     !
!     call find_k_1BZ_and_G(qpoint,nqf1,nqf2,nqf3,i,j,k,qpoint_in_1bz,GKQ_bz)
!     call switch_indices(nqf1,nqf2,nqf3,iqpt,i,j,k,+1)
!     qpoint(:)=qfmesh(:,iqpt)
     qpoint_cart=matmul(bg,qpoint)
     !
     call DM_Rp2q(nrr_q,irr_q,ndegen_q,qpoint,dyn_r,U_ph(:,:,iq),omega(:,iq))
     !
!     write(1314,'(i10,3f8.5,f12.6,100f15.8)') iq,qpoint(:), &
!           sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2),eigen(:,iqpt)
     !
     if (iq.lt.168) then
        !
        write(1314,'(i10,3f8.5,f12.6,100f15.8)') iq,qpoint(:), &
           sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2),omega(:,iq)!abs(omega(:,iqpt))
        !
     elseif (iq.gt.167.and.iq.lt.251) then
        !
        write(1314,'(i10,3f8.5,f12.6,100f15.8)') iq,qpoint(:), &
           1.d0-sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2-1.d0/3.d0),omega(:,iq)!abs(omega(:,iqpt))
        !
     elseif (ik.gt.250) then
        !
        write(1314,'(i10,3f8.5,f12.6,100f15.8)') iq,qpoint(:), &
           1.d0+sqrt(1.d0/3.d0)-sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2),omega(:,iq)!abs(omega(:,iqpt))
        !
     endif
     !
  enddo !iq
  !
  close(1314)

end program imaginary_selfen
