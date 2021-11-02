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
  complex(dp), allocatable :: trans(:,:,:),trans_scal(:,:,:)
  complex(dp), allocatable :: fr(:,:,:), fg(:,:,:)

  !wave function realted variables information
  integer, allocatable     :: list_igk (:)
  integer, allocatable     :: list_igkq(:)
  integer, allocatable     :: list_igk_aux (:)
  integer, allocatable     :: list_igk_orig (:)

  complex(dp), allocatable :: wfc_k (:,:,:),wfc_r_k(:,:),wfc_r_kq(:,:) ! nG_max is defined in reading  
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
  real(dp)                 :: term,sigma_e,arg,modkq,modkk,moduloa
  complex(dp) :: fac,akq,akk,bkq,bkk,phikq,phikk,chikq(2),chikk(2),a,b
  complex(dP),allocatable :: a_r(:,:),b_r(:,:),a_kq(:,:),b_kq(:,:),a_kk(:,:),b_kk(:,:)
  integer                  :: nbands_loc,irk,irq,iGkq,iGk,iG_1,iG_2,ispin
  integer                  :: i,j,k,is,js
  integer                  :: ii,jj,kk
  integer                  :: jr,kr
  integer                  :: ig,ibnd,jbnd,ipol,jpol,iibnd,jjbnd,ifs,nkk,nfs
  integer                  :: switch
  logical                  :: read_status,have_nnkp,confirm,calculate
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
  allocate(trans(nkmesh,num_bands,num_bands))
  trans(:,:,:)=cmplx_0
  allocate(wfc_r_kq(nr1*nr2*nr3,nspin),wfc_r_k(nr1*nr2*nr3,nspin))
  allocate(a_r(nr1*nr2*nr3,num_bands),b_r(nr1*nr2*nr3,num_bands))
  allocate(a_kq(nG_max,num_bands),b_kq(nG_max,num_bands))
  allocate(a_kk(nG_max,num_bands),b_kk(nG_max,num_bands))
  allocate(trans_scal(nkmesh,4,4))
  trans_scal(:,:,:)=cmplx_0
  !
  qpoint(1)=0.3333333d0
  qpoint(2)=0.3333333d0
  qpoint(3)=0.0d0
  !
  ep_unit=find_free_unit()
  inquire(iolength=record_lengh) trans
  open(unit=ep_unit,file=trim(trim(mesh_dir)//trim('trans_sin_spinor.dat_scal')),&
  iostat=ierr,form='unformatted',status='unknown',access='direct',recl=record_lengh)
  read(unit=ep_unit,rec=1,iostat=ierr) trans_scal(:,:,:)
  close(unit=ep_unit)
  !
!  calculate=.false.
!  if (calculate) then
     !
     do ik=1,nkmesh
        !
        print*, "ik= ", ik,' k= ',matmul(bg,kmesh(:,ik))
        kpoint=kmesh(:,ik)
        !
        !-uhina lortu RAM memorian dauden uhin irreduzibleak errotatuta. (k+q)
        !
        call get_psi_general_k_all_wfc(.true.,kpoint       ,list_iGk ,wfc_k ,QE_eig_k ,G_plusk)
        call get_psi_general_k_all_wfc(.true.,kpoint+qpoint,list_iGkq,wfc_kq,QE_eig_kq,G_pluskq)
        !
!        !$omp parallel default(none) &
!        !$omp shared(nkmesh,kmesh,kpoint,ik,list_iGk,list_iGkq,wfc_k,wfc_kq,qpoint,num_bands,nG_max,trans,cmplx_1,nspin,trans_scal,cmplx_0) &
!        !$omp private(ibnd,jbnd,iGkq,iGk,iG_1,iG_2,a_kq,b_kq,a_kk,b_kk,moduloa,a,b,chikq,chikk,iibnd,jjbnd)
!        !$omp do
        !
        do ibnd=1,num_bands
           !
           a_kq(:,:)=cmplx_0
           b_kq(:,:)=cmplx_0
           !
           do iGkq=1,nG_max
              !
              iG_1=list_iGkq(iGkq)
              !
              if (iG_1==0) exit
              !
              a=cmplx_1
              b=wfc_kq(iGkq,ibnd,2)/wfc_kq(iGkq,ibnd,1)
              moduloa=sqrt(conjg(a)*a+conjg(b)*b)
              !
              a_kq(iGkq,ibnd)=a/moduloa
              b_kq(iGkq,ibnd)=b/moduloa
              !
           enddo !iGkq
           !
           chikq(1)=sum(a_kq(:,ibnd))/nG_max
           chikq(2)=sum(b_kq(:,ibnd))/nG_max
           !
           do jbnd=1,num_bands
              !
              a_kk(:,:)=cmplx_0
              b_kk(:,:)=cmplx_0
              !
              do iGk=1,nG_max
                 !
                 iG_2=list_iGk(iGk)
                 !
                 if (iG_2==0) exit
                 !
                 a=cmplx_1
                 b=wfc_k(iGk,jbnd,2)/wfc_k(iGk,jbnd,1)
                 moduloa=sqrt(conjg(a)*a+conjg(b)*b)
                 !
                 a_kk(iGk,jbnd)=a/moduloa
                 b_kk(iGk,jbnd)=b/moduloa
                 !
              enddo !iGk
              !
              chikk(1)=sum(a_kk(:,jbnd))/nG_max
              chikk(2)=sum(b_kk(:,jbnd))/nG_max
              !
              if (ibnd.lt.3) iibnd=1
              if (ibnd.gt.2.and.ibnd.lt.5) iibnd=2
              if (ibnd.gt.4) iibnd=3
              !
              if (jbnd.lt.3) jjbnd=1
              if (jbnd.gt.2.and.jbnd.lt.5) jjbnd=2
              if (jbnd.gt.4) jjbnd=3              
              !
              trans(ik,ibnd,jbnd)=trans_scal(ik,iibnd,jjbnd)*(conjg(chikq(1))*chikk(1)+conjg(chikq(2))*chikk(2))
              !
!              do iGkq=1,nG_max
!                 !
!                 iG_1=list_iGkq(iGkq)
!                 !
!                 if (iG_1==0) exit
!                 !
!                 do iGk=1,nG_max
!                    !
!                    iG_2=list_iGk(iGk)
!                    !
!                    if (iG_1==iG_2) then
!                       !
!                       do is=1,nspin
!                          !
!                          trans(ik,ibnd,jbnd)=trans(ik,ibnd,jbnd)+conjg(wfc_kq(iGkq,ibnd,is))*wfc_k(iGk,jbnd,is)
!                          !
!                       enddo !is
!                       !
!                       akq=cmplx_1
!                       bkq=wfc_kq(iGkq,ibnd,2)/wfc_kq(iGkq,ibnd,1)
!                       modkq=sqrt(conjg(akq)*akq+conjg(bkq)*bkq)
!                       phikq=wfc_kq(iGkq,ibnd,1)*modkq
!                       !
!                       akk=cmplx_1
!                       bkk=wfc_k(iGk,jbnd,2)/wfc_k(iGk,jbnd,1)
!                       modkk=sqrt(conjg(akk)*akk+conjg(bkk)*bkk)
!                       phikk=wfc_k(iGk,jbnd,1)*modkk
!                       !
!                       trans(ik,ibnd,jbnd)=trans(ik,ibnd,jbnd)+conjg(phikq)*phikk*(conjg(akq)*akk+conjg(bkq)*bkk)/(modkk*modkq)
!                       trans(ik,ibnd,jbnd)=trans(ik,ibnd,jbnd)+conjg(phikq)*phikk
!                       !
!                       do is=1,nspin
!                          !
!                          trans(ik,ibnd,jbnd)=trans(ik,ibnd,jbnd)+sum(conjg(wfc_r_kq(:,is))*wfc_r_k(:,is))/(nr1*nr2*nr3)
!                          !
!                       enddo !is
!                       !
!                    endif
!                    !
!                 enddo !iGkk
!              enddo !iGkq
!              !
!              do is=1,nspin
!                 !
!                 trans(ik,ibnd,jbnd)=trans(ik,ibnd,jbnd)+conjg(chikq(is))*chikk(is)
!                 !
!              enddo
!              trans(ik,ibnd,jbnd)=cmplx_1
              !
           enddo !jbnd
        enddo !ibnd
        !
!        !$omp end parallel
        !
     enddo !ik
     !
!  endif
!stop
  !
  ep_unit=find_free_unit()
  inquire(iolength=record_lengh) trans
  open(unit=ep_unit,file=trim(trim(mesh_dir)//trim('trans_sin_spinor.dat')),&
  iostat=ierr,form='unformatted',status='unknown',access='direct',recl=record_lengh)
  write(unit=ep_unit,rec=1,iostat=ierr) trans(:,:,:)
!  read(unit=ep_unit,rec=1,iostat=ierr) trans(:,:,:)
  close(unit=ep_unit)
  !
  temp=0.0d0 ! Temperature of the system
  !
  if (lspinorb) then
     trshld=-2.88073485d0
!     trshld=1.6087
  else
     trshld=-2.60126050d0
!     trshld=1.6087
  endif
  !
  allocate(sumg1(nkmesh),sumg2(nkmesh))
  !
  ef=-trshld
  !
  sumg1(:)=0.d0
  sumg2(:)=0.d0
  !
  do ik=1,nkmesh
     !
     call get_psi_general_k_all_wfc(.true.,kpoint       ,list_iGk ,wfc_k ,QE_eig_k ,G_plusk)
     call get_psi_general_k_all_wfc(.true.,kpoint+qpoint,list_iGkq,wfc_kq,QE_eig_kq,G_pluskq)
     kpoint(:)=kmesh(:,ik)
     call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
     call switch_indices(nk1,nk2,nk3,ikpt_k,i,j,k,+1)
     kpoint(:)=kmesh(:,ikpt_k)
     call find_k_1BZ_and_G(kpoint+qpoint,nk1,nk2,nk3,i,j,k,kqpoint_in_1bz,GKQ_bz)
     call switch_indices(nk1,nk2,nk3,ikpt_kq,i,j,k,+1)
     !
     do ibnd=1,num_wann
        !
        ekq=QE_eig_kq(ibnd)
        !
        do jbnd=1,num_wann
           !
           ekk=QE_eig_k(jbnd)
           !
           if (ekk.lt.ef.and.ekq.gt.ef) then
              !
              g1=trans(ik,ibnd,jbnd)
              g2=(abs(g1)**2.d0)
              !
              sumg1(ik)=sumg1(ik)+g2
              sumg2(ik)=sumg2(ik)+g2/(ekk-ekq)
              !
           elseif (ekk.gt.ef.and.ekq.lt.ef) then
              !
              g1=trans(ik,ibnd,jbnd)
              g2=(abs(g1)**2.d0)
              !
              sumg1(ik)=sumg1(ik)+g2
              sumg2(ik)=sumg2(ik)+g2/(ekq-ekk)              
              !
           endif !ekk->ekq + ekq->ekk
           !
        enddo !jbnd
     enddo !ibnd
     !
  enddo !ik
  !
  allocate(sumg1_tex(-nk1:nk1,-nk2:nk2),sumg2_tex(-nk1:nk1,-nk2:nk2))
  sumg1_tex(:,:)=0.0d0
  sumg2_tex(:,:)=0.0d0
  !
  do ik=1,nkmesh
     !
     kpoint(:)=kmesh(:,ik)
     call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
     !
     sumg1_tex((i-1)-nk1,(j-1)-nk2)=sumg1(ik)
     sumg1_tex((i-1)-nk1,j-1)=sumg1(ik)
     sumg1_tex(i-1,(j-1)-nk2)=sumg1(ik)
     sumg1_tex(i-1,j-1)=sumg1(ik)
     !
     sumg2_tex((i-1)-nk1,(j-1)-nk2)=sumg2(ik)
     sumg2_tex((i-1)-nk1,j-1)=sumg2(ik)
     sumg2_tex(i-1,(j-1)-nk2)=sumg2(ik)
     sumg2_tex(i-1,j-1)=sumg2(ik)
     !
  enddo !ik
  !
  do ik=-nk1,nk1-1
     !
     sumg1_tex(ik,nk2)=sumg1_tex(ik,0)
     sumg2_tex(ik,nk2)=sumg2_tex(ik,0)
     !
  enddo !ik
  !
  do ik=-nk2,nk2-1
     !
     sumg1_tex(nk1,ik)=sumg1_tex(0,ik)
     sumg2_tex(nk1,ik)=sumg2_tex(0,ik)
     !
  enddo !ik
  !
  sumg1_tex(nk1,nk2)=sumg1_tex(0,0)
  sumg2_tex(nk1,nk2)=sumg2_tex(0,0)
  !
  open(unit=11,file=trim(trim(mesh_dir)//trim('trans_tex_sin_spinor.dat')))
  open(unit=22,file=trim(trim(mesh_dir)//trim('trans_zati_w_tex_sin_spinor.dat')))
  !
  do ikx=-nk1,nk1
     do iky=-nk2,nk2
        !
        kpoint(1)=1.d0*ikx/nk1
        kpoint(2)=1.d0*iky/nk2
        !
        write(11,'(2f16.8,f20.15)') kpoint(1),kpoint(2),sumg1_tex(ikx,iky)
        write(22,'(2f16.8,f20.15)') kpoint(1),kpoint(2),sumg2_tex(ikx,iky)
        !
     enddo !iky
  enddo !ikx
  !
  close(11)
  close(22)

end program imaginary_selfen
