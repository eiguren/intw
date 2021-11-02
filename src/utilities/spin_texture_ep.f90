program ep_melements

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
  use w90_parameters, only: num_wann,bands_num_points,real_metric,recip_metric,&
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
  real(dp)                 :: kpoint(3), kpoint_cart(3), kqpoint(3), kqpoint_in_1bz(3), kqpoint_cart(3)
  real(dp)                 :: kpoint_in_1bz(3), kpoint_rot(3) 

  !q point related variables
  real(dp)                 :: qpoint(3), qpoint_1bz(3), qpoint_rot(3)

  !phonon related variables
  character(len=256)       :: fc_file_name
  real(dp)                 :: at_frc(3,3)
  integer                  :: mode, q_index_irr, q_index, imode
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
  integer                  :: ep_unit
  logical                  :: q_points_consistent 

  !wave function realted variables information
  integer, allocatable     :: list_igk (:)
  integer, allocatable     :: list_igkq(:)
  integer, allocatable     :: list_igk_aux (:)
  integer, allocatable     :: list_igk_orig (:)

  complex(dp), allocatable :: wfc_k (:,:,:) ! nG_max is defined in reading  
  complex(dp), allocatable :: wfc_kq (:,:,:)

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
  real(dp),allocatable :: kfmesh(:,:)

  !spin components
  integer :: spin_kb_unit,spin_Rw_unit,ikx,iky
  complex(dp) :: a,b
  real(dp),allocatable :: spin_m(:,:,:,:)
  complex(dp),allocatable :: spin_kb(:,:,:,:),a_spin_kb(:,:)
  complex(dp),allocatable :: spin_Rw(:,:,:,:),a_spin_Rw(:,:,:)
  complex(dp),allocatable :: dvpsi_spin(:,:,:,:)
  real(dp),allocatable :: suma(:,:,:)

  !local/aux variables
  integer                  :: nstar, ntot, na
  integer                  :: suma_i
  complex(dp)              :: suma_c,term
  real(dp)                 :: suma_r 
  integer                  :: LEN
  integer                  :: nbands_loc
  integer                  :: npw, npwq
  
  integer                  :: i, j, k,is,iG11,iG22,jG
  integer                  :: ii, jj, kk
  integer                  :: jr, kr
  integer                  :: ig, ibnd, jbnd, ipol, jpol
  integer                  :: switch
  logical                  :: read_status, have_nnkp
  character(256)           :: nnkp_file, method
  character(len=4)         :: iq_loc,ibnd_loc,imode_loc

  complex(dp),allocatable :: wfc_k_r(:)
  integer :: igg1, igg2, gg1(3), gg2(3), nG, ikb
  logical :: found
  complex(dp),allocatable :: auxil(:,:,:)
  integer,allocatable :: list_aux(:)

  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
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
  !       Check that $prefix.nnkp is present
  !================================================================================
  !
  nnkp_file = trim(mesh_dir)//trim(prefix)//".nnkp"
  !
  inquire(file=nnkp_file,exist=have_nnkp)
  !
  if (.not.have_nnkp) then
     !
     write(*,20)      '**********************************************************'
     write(*,20)      '* Could not find the file '//trim(nnkp_file)
     write(*,20)      '* Did you run W90 -pp $seed to get the parameter file?   ' 
     write(*,20)      '**********************************************************'
     !
     stop
     !
  endif
  !
  write(*,20) '|       - .nnkp file found                          |'
  write(*,20) '|           ---------------------------------       |'
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
  write(*,20) '|       - The mesh in the Wannier90 input.win file  |'
  write(*,20) '|         and the intw mesh are equal.              |'
  write(*,20) '|           ---------------------------------       |'
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
  if (.not.lspinorb) call average_pp(ntyp)
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
  !================================================================================
  !       set up W90 data and generate the u_mesh array 
  !================================================================================
  !
  ! CAREFUL!! It is CRUCIAL to allocate and read the W90 stuff BEFORE
  ! referencing it. num_wann, for example, exists but is not properly
  ! DEFINED before.
  !
  call allocate_and_read_W90()
  !
  ! extract the nband x num_wann Wannier projection+rotation matrices.
  !
  call produce_u_mesh()
  !
  !================================================================================
  !       Get the special points from the file
  !================================================================================
  !
  if (bands_points/=0) then
     !
     call get_interpolation_special_points(number_of_special_points)
     !
  endif
  !
  !================================================================================
  !       Allocate wfc related files
  !================================================================================
  !
  nbands_loc=num_bands
  !
  allocate (list_igk(nG_max))
  allocate (list_igkq(nG_max))
  !
  allocate (wfc_k(nG_max,num_bands,nspin))
  allocate (wfc_kq(nG_max,num_bands,nspin))
  !
  allocate (wfc_k_r(nr1*nr2*nr3))
  !
  allocate (QE_eig_k(num_bands))
  allocate (QE_eig_kq(num_bands))
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
  fc_file_name=trim(trim(mesh_dir)//"/"//trim(ph_dir)//"/"//trim(fc_mat))
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
  allocate(dvpsi_spin(nG_max,nbands_loc,npol,3*nat))
  !
  call generate_kmesh(qmesh,nq1,nq2,nq3)
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
  call allocate_nlpot1
  !
  call allocate_phq
  !
  call init_us_1
  !
  iq=1
  !
  qpoint(1)=0.3333333d0
  qpoint(2)=0.3333333d0
  qpoint(3)=0.0000000d0
  !
  call allocate_nlpot2 (qpoint)
  !
  dvq_local=cmplx_0
  call get_dv(iq,qpoint,3*nat,dvq_local)
  !
  call phq_init(matmul(bg,qpoint))
  call calculate_local_part_dv(qpoint,nat,npol,dvq_local)
  !
  allocate(spin_kb(nkmesh,nbands_loc,3*nat,3))
  allocate(suma(nbands_loc,nspin,3*nat))
  spin_kb(:,:,:,:)=cmplx_0
  !
  ep_unit=find_free_unit()
  inquire(iolength=record_lengh) spin_kb
  open(unit=ep_unit,file=trim(trim(mesh_dir)//trim('spin_kb.dat')),iostat=ierr,&
  form='unformatted',status='unknown',access='direct',recl=record_lengh)
  !
  !!$omp parallel default(none) &
  !!$omp shared(kmesh,nkmesh,bg,qpoint,nk1,nk2,nk3,nG_max,nat,nbands_loc,npol,dvq_local,spin_kb,cmplx_0,num_bands,nspin) &
  !!$omp private(ik,kpoint,kpoint_cart,kqpoint_cart,list_iGk,wfc_k,QE_eig_k,G_plusk,list_iGkq,wfc_kq,QE_eig_kq,G_pluskq,i,j,k,kpoint_in_1bz,GKQ_bz,ikpt_k) &
  !!$omp private(kqpoint_in_1bz,ikpt_kq,nG,iG,vkb,vkqb,npw,npwq,dvpsi,ibnd,imode,is,dvpsi_spin)
  !!
  !!$omp do
  !
!  do ik=1,nkmesh
!     !
!     dvpsi(:,:,:,:,:)=cmplx_0
!     dvpsi_spin(:,:,:,:)=cmplx_0
!     !
!     print*, "ik= ", ik,' k= ',matmul(bg,kmesh(:,ik))
!     kpoint=kmesh(:,ik)
!     !
!     kpoint_cart=matmul(bg,kpoint)
!     kqpoint_cart=matmul(bg,kpoint+qpoint)
!     !
!     print*, "ikq= ", ik,' kq= ',kqpoint_cart
!     !
!     call get_psi_general_k_all_wfc(.true.,kpoint       ,list_iGk ,wfc_k ,QE_eig_k ,G_plusk)
!     call get_psi_general_k_all_wfc(.true.,kpoint+qpoint,list_iGkq,wfc_kq,QE_eig_kq,G_pluskq)
!     !
!     call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
!     call switch_indices(nk1,nk2,nk3,ikpt_k,i,j,k,+1)
!     !
!     call find_k_1BZ_and_G(kpoint+qpoint,nk1,nk2,nk3,i,j,k,kqpoint_in_1bz,GKQ_bz)
!     call switch_indices(nk1,nk2,nk3,ikpt_kq,i,j,k,+1)
!     !
!     nG=0
!     do iG=1,nG_max
!        if (list_iGk(iG)==0) exit
!        nG=nG+1
!     enddo
!     npw=nG
!     nG=0
!     do iG=1,nG_max
!        if (list_iGkq(iG)==0) exit
!        nG=nG+1
!     enddo
!     npwq=nG
!     !
!     vkb=cmplx_0
!     vkqb=cmplx_0
!     call init_us_2(npw ,nG_max,list_iGk ,kpoint       ,vkb)
!     call init_us_2(npwq,nG_max,list_iGkq,kpoint+qpoint,vkqb)
!     !
!     call dvqpsi_local (3*nat, nG_max, nbands_loc, npol, list_iGk, wfc_k, dvq_local, dvpsi(1:nG_max,1:nbands_loc,1:npol,1:npol,1:3*nat))
!     !
!     call dvqpsi_us_only (matmul(bg, kpoint), matmul(bg, qpoint), nat, nG_max, nbands_loc, npol, list_iGk, list_iGkq, wfc_k, dvpsi)
!     !
!     do iG=1,nG_max
!        do ibnd=1,num_bands
!           do imode=1,3*nat
!              do is=1,nspin
!                 !
!                 dvpsi_spin(iG,ibnd,is,imode)=dvpsi(iG,ibnd,is,is,imode)
!                 !
!              enddo !is
!           enddo !imode
!        enddo !ibnd
!     enddo !iG
!     !
!     suma(:,:,:)=0.d0
!     !
!     do ibnd=1,num_bands
!        do is=1,nspin
!           do imode=1,3*nat
!              !
!              do iG=1,nG_max
!                 !
!                 suma(ibnd,is,imode)=suma(ibnd,is,imode)+0.5d0*conjg(dvpsi_spin(iG,ibnd,is,imode))* &
!                                                                         dvpsi_spin(iG,ibnd,is,imode)
!                 !
!              enddo
!              !
!           enddo
!        enddo
!     enddo
!     !
!     do imode=1,3*nat
!        do ibnd=1,num_bands
!           do jbnd=1,num_bands
!              !
!              do iG=1,nG_max
!                 !
!                 iG11=list_iGkq(iG)
!                 !
!                 if (iG11==0) exit
!                 !
!                 do jG=1,nG_max
!                    !
!                    iG22=list_iGk(jG)
!                    !
!                    if (iG22==0) exit
!                    !
!                    if (iG11==iG22) then
!                       !
!                       do is=1,nspin
!                          !
!                          suma(ibnd,jbnd,imode)=suma(ibnd,jbnd,imode)+conjg(wfc_kq(iG,ibnd,is))*dvpsi_spin(jG,jbnd,is,imode)
!                          !
!                       enddo !is
!                       !
!                    endif !iG=jG
!                    !
!                 enddo !jG
!              enddo !iG
!              !
!           enddo !jbnd
!        enddo !ibnd
!     enddo !imode
!     !
!     do imode=1,3*nat
!        !
!        write(1340,'(a,i3)') 'phonon imode = ',imode
!        write(1340,*) ''
!        write(1340,*) ''
!        !
!        do ibnd=1,num_bands
!           !
!           write(1340,'(a,i2)') 'electron ibnd = ',ibnd
!           write(1340,*) ''
!           !
!           do is=1,nspin
!              !
!              write(1340,'(G12.6,a,G12.6)') real(suma(ibnd,is,imode))
!              !
!           enddo
!           !
!           write(1340,*) ''
!           write(1340,*) '------------------------------------------------------------------------------' 
!           !
!        enddo
!        !
!        write(1340,*) '          jbnd = 1    jbnd = 2    jbnd = 3    jbnd = 4    jbnd = 5    jbnd = 6'
!        write(1340,*) ''
!        !
!        do ibnd=1,num_bands
!           !
!           write(1340,'(a,i1,a,6(x,G10.3,x))') 'ibnd = ',ibnd,' ',(abs(suma(ibnd,jbnd,imode)),jbnd=1,num_bands)
!           write(1340,*) ''
!           !
!        enddo
!        !
!        write(1340,*) ''
!        write(1340,*) '******************************************************************************'
!        !
!     enddo
!stop
!     !
!     do imode=1,3*nat
!        !
!        call get_spin_component_convolution(list_iGk,dvpsi_spin(1:nG_max,1:num_bands,1:nspin,imode),spin_kb(ik,1:nbands_loc,imode,1:3))
!        !
!     enddo 
!     !
!  enddo !ik
  !
  !!$omp end parallel
  !
!  write(unit=ep_unit,rec=1,iostat=ierr) spin_kb(:,:,:,:)
  read(unit=ep_unit,rec=1,iostat=ierr) spin_kb(:,:,:,:)
  close(unit=ep_unit)
  !
  allocate(spin_m(-nk1:nk1,-nk2:nk2,num_bands,3))
  !
  do imode=1,3*nat
     !
     spin_m=0.0d0
     do ipol=1,3
        do ibnd=1,num_bands
           !
           do ik=1,nkmesh
              !
              kpoint(:)=kmesh(:,ik)
              call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
              !
              spin_m((i-1)-nk1,(j-1)-nk2,ibnd,ipol)=real(spin_kb(ik,ibnd,imode,ipol))
              spin_m((i-1)-nk1,j-1,ibnd,ipol)=real(spin_kb(ik,ibnd,imode,ipol))
              spin_m(i-1,(j-1)-nk2,ibnd,ipol)=real(spin_kb(ik,ibnd,imode,ipol))
              spin_m(i-1,j-1,ibnd,ipol)=real(spin_kb(ik,ibnd,imode,ipol))
              !
           enddo !ik
           !
           do ik=-nk1,nk1-1
              !
              spin_m(ik,nk2,ibnd,ipol)=spin_m(ik,0,ibnd,ipol)
              !
           enddo !ik
           !
           do ik=-nk2,nk2-1
              !
              spin_m(nk1,ik,ibnd,ipol)=spin_m(0,ik,ibnd,ipol)
              !
           enddo !ik
           !
           spin_m(nk1,nk2,ibnd,ipol)=spin_m(0,0,ibnd,ipol)
           !
        enddo !ibnd
     enddo !ipol
     !
     do ibnd=1,num_bands
        !
        write(ibnd_loc,'(i1)')ibnd
        if (imode.lt.10) write(imode_loc,'(i1)')imode
        if (imode.ge.10) write(imode_loc,'(i2)')imode 
        !
        open(unit=11,file=trim(trim(mesh_dir)//trim('spin_fine_ep')//(trim(ibnd_loc))//(trim('_'))//(trim(imode_loc))//trim('.dat')))
        !
        do ikx=-nk1,nk1
           do iky=-nk2,nk2
              !
              kpoint(1)=1.d0*ikx/nk1
              kpoint(2)=1.d0*iky/nk2
              !
              write(11,'(100f16.8)') kpoint(1),kpoint(2),(spin_m(ikx,iky,ibnd,ipol),ipol=1,3)
              !
           enddo !iky
        enddo !ikx
        !
        close(11)
        !
     enddo !ibnd
     !
  enddo !imode
  !
  stop

  !================================================================================
  !       Define the interpolation path     
  !================================================================================

  nk_vec_path = bands_points

  allocate(list_k_path(3,nk_vec_path))
  allocate(list_x_path(nk_vec_path))
  allocate(list_xticks(number_of_special_points))

  call build_list_kpath(number_of_special_points,nk_vec_path,    &
       list_k_path,list_x_path, list_xticks)


  !io_unit=find_free_unit()
  !open(unit=io_unit, file=trim(ph_bands_file),status="unknown")

  !do i=1,nk_vec_path
  ! qpoint(1:3) =  list_k_path(1:3,i)
  ! call mat_inv_four_t ( qpoint, nq1, nq2, nq3,  3*nat , frc , dyn_q )
  ! call diagonalize_cmat (3*nat, dyn_q,w2)
  ! write(io_unit, "(i5,10000f18.12)")  i, sqrt(abs(w2))
  !enddo

  !close(unit=io_unit)

  deallocate(dyn_q, w2)

  call deallocate_upfeak ()
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
   
!********************************************************************************************************
!-----------------------------------------------------------------------------------------------------
  subroutine dvqpsi_local(nmode,nG_max,nbands,npol,list_iGk,wfc_k,dvq_local,dvpsi_local)
!-----------------------------------------------------------------------------------------------------

    implicit none

    !I/O variables
     
    integer,intent(in) :: nmode,nG_max,nbands,npol,list_iGk(nG_max)
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
                   dvpsi_local(ig,ibnd,ipol,jpol,imode)=wfc_r(nl(ig),ipol,jpol)
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
!********************************************************************************
  !----------------------------------------
  subroutine diagonalize_cmat (n,a,w)
    !----------------------------------------
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

  !--------------------------------------
  subroutine wfc_bands_to_wann(wfc_k,U_k,wfc_k_W)
   !----------------------------------------

   integer :: is, nG, i, j
   complex(dp),intent(in) :: wfc_k(nG_max,nbands_loc,nspin)
   complex(dp),intent(in) :: U_k(nbands_loc,num_wann)
   complex(dp),intent(out) :: wfc_k_W(nG_max,nbands_loc,nspin)
   
   wfc_k_W=(0.d0,0.d0)
   
   do nG=1,nG_max
   do is=1,nspin
   do i=1,num_wann
      do j=1,nbands_loc
         wfc_k_W(nG,i,is)=wfc_k_W(nG,i,is)+U_k(j,i)*wfc_k(nG,j,is)
      enddo
   enddo
   enddo
   enddo   

  end subroutine wfc_bands_to_wann

end program ep_melements

