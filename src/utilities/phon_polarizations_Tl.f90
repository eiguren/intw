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
  real(dp)                 :: kpoint(3), kpoint_cart(3), kqpoint(3), kqpoint_cart(3), kqpoint_in_1bz(3)
  real(dp)                 :: kpoint_in_1bz(3), kpoint_rot(3), qpoint_in_1bz(3)

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
  complex(dp), allocatable :: aep_mat_el(:,:,:,:), ep_mat_el(:,:,:)
  complex(dp), allocatable :: aW_rot_ep_mat_el(:,:,:,:), W_rot_ep_mat_el(:,:,:)
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
  integer :: nqf1,nqf2,nqf3,nqfmesh
  real(dp),allocatable :: kfmesh(:,:)
  real(dp),allocatable :: qfmesh(:,:)

  !matrix-elements
  integer :: g_realw_unit
  complex(dp) :: g1
  real(dp) :: g2
  complex(dp),allocatable :: g_ReRp(:,:,:,:),ag_ReRp(:,:,:)
  complex(dp),allocatable :: g_kq(:,:,:,:),g_Req(:,:,:,:)
  complex(dp),allocatable :: g_Re(:,:,:),g_wannier(:,:)
  complex(dp),allocatable :: g_elph(:,:)
  real(dp),allocatable :: sumg1(:),sumg2(:),sumg1_tex(:,:),sumg2_tex(:,:)

  !Reciprocal Bloch to Real Wannier
  logical :: kqb_2_ReRpw_me,kqb_2_kqw_me,h_kb_2_Rew,dm_q_2_Rp

  !selfenergy
  real(dp) :: temp ! temperature of the system
  real(dp) :: trshld,dkvec(3),wqv,w_max
  real(dp) :: ekq,ekk,d_w,sigma_w
  real(dp) :: ef
  real(dp),allocatable :: imagse(:)
  integer :: iq_cont_plus,iq_cont_minus,iq_2,nqfmesh_plus,nqfmesh_minus,ikx,iky
  integer,allocatable :: list_iq_plus(:,:),list_iq_minus(:,:)
  real(dp) :: termplus,termminus

  !Tl polarizations
  complex(dp),allocatable :: tl_vibration(:,:,:),tl_pol(:,:,:,:)

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

  !local/aux variables
  real(dp)                 :: term,sigma_e,arg
  complex(dp) :: fac
  integer                  :: nbands_loc,irk,irq,iGkq,iGk,iG_1,iG_2,ispin
  integer                  :: i,j,k,is,js
  integer                  :: ii,jj,kk
  integer                  :: jr,kr
  integer                  :: ig,ibnd,jbnd,ipol,jpol,iibnd,jjbnd,ifs,nkk,nfs
  integer                  :: switch
  logical                  :: read_status,have_nnkp,confirm
  logical                  :: deg_k_12,deg_k_34,deg_k_56
  logical                  :: deg_kq_12,deg_kq_34,deg_kq_56
  character(256)           :: nnkp_file,method
  character(len=4)         :: iq_loc,ifs_loc,imode_loc
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
  allocate (aep_mat_el(nqmesh,nk1*nk2*nk3,nbands_loc,nbands_loc))  ! q bakoitzerako guztia
  allocate (aW_rot_ep_mat_el(nqmesh,nk1*nk2*nk3,num_wann,num_wann))
  allocate (W_rot_ep_mat_el(nk1*nk2*nk3,num_wann,num_wann))
  allocate (ep_mat_el(nk1*nk2*nk3,nbands_loc,nbands_loc))
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
  nrr_q=0
  !
  call count_wigner_seitz_points(nq1,nq2,nq3,nrr_q)
  allocate(ndegen_q(nrr_q))
  allocate(irr_q(3,nrr_q))
  ndegen_q(:)=0.0d0
  call wigner_seitz_points(nq1,nq2,nq3,nrr_q,irr_q,ndegen_q)
  !
  dm_q_2_Rp=.true.
  !
  write(*,*)'******************************************************************'
  write(*,*)'Phonons'
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
  write(*,*)'Definition of the phonon fine grid...'
  !
  nqf1=24
  nqf2=24
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
  write_ph=.true.
  !
  if (iwant) then
     !
     allocate(U_ph(3*nat,3*nat,nqfmesh))
     allocate(omega(3*nat,nqfmesh))
     U_ph(:,:,:)=cmplx_0
     omega(:,:)=0.0d0
     !
     allocate(a_ph_q_f(3*nat,3*nat))
     allocate(a_omega_f(3*nat))
     !
     if (write_ph) then
        !
        ! We calculate in the whole BZ
        !
        !$omp parallel default(none) &
        !$omp shared(nqfmesh,qfmesh,nrr_q,irr_q,ndegen_q,dyn_r,U_ph,omega) &
        !$omp private(iq,qpoint) 
        !
        !$omp do
        ! 
        do iq=1,nqfmesh
           !
           if ((iq/10000)*10000==iq) then
              write(*,*)'iq = ',iq,' / ',nqfmesh
           endif
           !
           qpoint(:)=qfmesh(:,iq)
           call DM_Rp2q(nrr_q,irr_q,ndegen_q,qpoint,dyn_r,U_ph(:,:,iq),omega(:,iq))
           !
        enddo
        !
        !$omp end parallel
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
           if ((iq/10000)*10000==iq) then
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
           if ((iq/10000)*10000==iq) then
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
  allocate(tl_vibration(nqfmesh,3,3))
  !
  do iq=1,nqfmesh
     !
     qpoint(:)=qfmesh(:,iq)
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
write(1789,*) iq,qpoint,omega_phon(:)
write(1789,*) ''
write(1789,*) '------'
     !
     do imode=1,3
        !
        do ipol=1,3
           !
           tl_vibration(iq,ipol,imode)=U_phon(ipol,imode)
           !
        enddo !ipol
        !
     enddo !imode
     !
  enddo !iq
  !
  allocate(tl_pol(-nqf1:nqf1,-nqf2:nqf2,3,3))
  tl_pol(:,:,:,:)=cmplx_0
  !
  do ipol=1,3
     do imode=1,3
        !
        do iq=1,nqfmesh
           !
           qpoint(:)=qfmesh(:,iq)
           call find_k_1BZ_and_G(qpoint,nqf1,nqf2,nqf3,i,j,k,kpoint_in_1bz,GKQ_bz)
           !
           tl_pol((i-1)-nqf1,(j-1)-nqf2,ipol,imode)=tl_vibration(iq,ipol,imode)
           tl_pol((i-1)-nqf1,j-1,ipol,imode)=tl_vibration(iq,ipol,imode)
           tl_pol(i-1,(j-1)-nqf2,ipol,imode)=tl_vibration(iq,ipol,imode)
           tl_pol(i-1,j-1,ipol,imode)=tl_vibration(iq,ipol,imode)   
           !
        enddo !iq
        !
        do iq=-nqf1,nqf1-1
           !
           tl_pol(iq,nqf2,ipol,imode)=tl_pol(iq,0,ipol,imode)
           !
        enddo !ik
        !
        do iq=-nqf2,nqf2-1
           !
           tl_pol(nqf1,iq,ipol,imode)=tl_pol(0,iq,ipol,imode)
           !
        enddo !ik
        !
        tl_pol(nqf1,nqf2,ipol,imode)=tl_pol(0,0,ipol,imode)
        !
     enddo !imode
  enddo !ipol
  !
  do imode=1,3
     !
     write(imode_loc,'(i1)')imode
     !
     open(unit=11,file=trim(trim(mesh_dir)//trim('tl_vibration_')//(trim(imode_loc))//trim('.dat')))
     !
     do ikx=-nqf1,nqf1
        do iky=-nqf2,nqf2
           !
           qpoint(1)=1.d0*ikx/nqf1
           qpoint(2)=1.d0*iky/nqf2
           !
           write(11,'(100f16.8)') qpoint(1),qpoint(2),(tl_pol(ikx,iky,ipol,imode),ipol=1,3)
           !
        enddo !iky
     enddo !ikx
     !
     close(11)
     !
  enddo !imode
  !

end program imaginary_selfen
