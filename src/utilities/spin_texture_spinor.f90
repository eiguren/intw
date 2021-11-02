program spin_texture_spinor

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
  real(dp),allocatable :: kfmesh(:,:)

  !spin components
  integer :: spin_kb_unit,spin_Rw_unit,ikx,iky
  complex(dp) :: a,b,chi_kk_spinstate(2)
  real(dp),allocatable :: spin_m(:,:,:,:),spin_m_tot(:,:,:),grad_spin(:,:,:),spin_kb(:,:,:)
  complex(dp),allocatable :: a_spin_kb(:,:)
  complex(dp),allocatable :: spin_kb_xy(:,:,:,:)
  complex(dp),allocatable :: spin_Rw(:,:,:,:),a_spin_Rw(:,:,:)
  complex(dp),allocatable :: chi_kk(:,:,:)

  !Reciprocal Bloch to Real Wannier
  logical :: kb_2_Rw_spin,kb_spin,h_kb_2_Rew

  !Hamiltonian (Electrons)
  logical :: write_el
  integer :: h_realw_unit,el_fine_unit,eigen_fine_unit
  real(dp),allocatable :: eigen(:,:),eigen_w(:)
  complex(dp),allocatable :: ham_rw(:,:,:),a_h_rw(:,:),U_el(:,:,:),U_el_w(:,:)

  !local/aux variables
  real(dp)                 :: sigma_w,term,sigma_e,tita,alpha,moduloa
  integer                  :: nbands_loc,irk,irq,ikq,nGnorm,iGk,iG_1
  integer                  :: i,j,k
  integer                  :: ii,jj,kk
  integer                  :: jr,kr
  integer                  :: ig,ibnd,jbnd,ipol,jpol,iibnd,jjbnd,ifs,nkk,nfs
  integer                  :: switch,ik1,ik2
  logical                  :: read_status,have_nnkp,confirm
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
     call get_interpolation_special_points(number_of_special_points)
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
  ! variables related with the wavefunction // spin-state
  !
  allocate(chi_kk(nG_max,nbands_loc,2))
  !
  allocate(spin_kb(nkmesh,nbands_loc,4))
  allocate(spin_kb_xy(nk1,nk2,nbands_loc,4))
  allocate(a_spin_kb(nbands_loc,4))
  spin_kb(:,:,:)=0.d0
  !
  do ik=1,nkmesh
     !
     ! first thing we want, for each k momentum try to ger just the spinor state from the
     ! whole spinor wave function.
     !
     print*, "ik= ", ik,' k= ',matmul(bg,kmesh(:,ik))
     kpoint=kmesh(:,ik)
     !
     call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
     call switch_indices(nk1,nk2,nk3,ikpt_k,i,j,k,+1) 
     !
     !-uhina lortu RAM memorian dauden uhin irreduzibleak errotatuta. (k+q)
     !
     call get_psi_general_k_all_wfc(.true.,kpoint       ,list_iGk ,wfc_k ,QE_eig_k ,G_plusk)
     !
     do ibnd=1,num_bands
        !
        chi_kk(:,:,:)=cmplx_0
        nGnorm=0
        chi_kk_spinstate(:)=cmplx_0
        !
        do iGk=1,nG_max
           !
           iG_1=list_iGk(iGk)
           !
           if (iG_1==0) exit
           !
!           if (abs(wfc_k(iGk,ibnd,2)).lt.abs(wfc_k(iGk,ibnd,1))) then
              !
              a=cmplx_1
              b=wfc_k(iGk,ibnd,2)/wfc_k(iGk,ibnd,1)
              moduloa=sqrt(conjg(a)*a+conjg(b)*b)
              !
!           elseif (abs(wfc_k(iGk,ibnd,1)).lt.abs(wfc_k(iGk,ibnd,2))) then
!              !
!              a=wfc_k(iGk,ibnd,1)/wfc_k(iGk,ibnd,2)
!              b=cmplx_1
!              moduloa=sqrt(conjg(a)*a+conjg(b)*b)
!              !
!           endif
           !
           chi_kk(iGk,ibnd,1)=a/moduloa
           chi_kk(iGk,ibnd,2)=b/moduloa
           !
           nGnorm=nGnorm+1
           !
        enddo !iGk
        !
        chi_kk_spinstate(1)=sum(chi_kk(:,ibnd,1))/nGnorm
        chi_kk_spinstate(2)=sum(chi_kk(:,ibnd,2))/nGnorm
        !
        spin_kb(ik,ibnd,1)=2.d0*real(conjg(chi_kk_spinstate(1))*chi_kk_spinstate(2))
        spin_kb(ik,ibnd,2)=2.d0*aimag(conjg(chi_kk_spinstate(1))*chi_kk_spinstate(2))
        spin_kb(ik,ibnd,3)=abs(chi_kk_spinstate(1))**2.d0-abs(chi_kk_spinstate(2))**2.d0
        spin_kb(ik,ibnd,4)=abs(chi_kk_spinstate(1))**2.d0+abs(chi_kk_spinstate(2))**2.d0
        !
!        do iGk=1,nG_max
!           !
!           spin_kb(ik,ibnd,1)=spin_kb(ik,ibnd,1)+2.d0*real(conjg(chi_kk(iGk,ibnd,1))*chi_kk(iGk,ibnd,2))
!           spin_kb(ik,ibnd,2)=spin_kb(ik,ibnd,2)+2.d0*aimag(conjg(chi_kk(iGk,ibnd,1))*chi_kk(iGk,ibnd,2))
!           spin_kb(ik,ibnd,3)=spin_kb(ik,ibnd,3)+abs(chi_kk(iGk,ibnd,1))**2.d0-abs(chi_kk(iGk,ibnd,2))**2.d0
!           spin_kb(ik,ibnd,4)=spin_kb(ik,ibnd,4)+abs(chi_kk(iGk,ibnd,1))**2.d0+abs(chi_kk(iGk,ibnd,2))**2.d0
!           !
!        enddo
!        !
!        spin_kb(ik,ibnd,:)=spin_kb(ik,ibnd,:)/nGnorm
        !
     enddo !ibnd 
     !
     spin_kb_xy(i,j,:,:)=spin_kb(ik,:,:)
     !
  enddo !ik
  !
  allocate(spin_m(-nk1:nk1,-nk2:nk2,num_bands,4))
  allocate(spin_m_tot(-nk1:nk1,-nk2:nk2,num_bands))
  allocate(grad_spin(-nk1:nk1,-nk2:nk2,num_bands))
  spin_m(:,:,:,:)=0.0d0
  spin_m_tot(:,:,:)=0.0d0
  do ipol=1,4
     do ibnd=1,num_bands
        !
        do ik=1,nkmesh
           !
           kpoint(:)=kmesh(:,ik)
           call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
           !
           spin_m((i-1)-nk1,(j-1)-nk2,ibnd,ipol)=spin_kb(ik,ibnd,ipol)
           spin_m((i-1)-nk1,j-1,ibnd,ipol)=spin_kb(ik,ibnd,ipol)
           spin_m(i-1,(j-1)-nk2,ibnd,ipol)=spin_kb(ik,ibnd,ipol)
           spin_m(i-1,j-1,ibnd,ipol)=spin_kb(ik,ibnd,ipol)
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
     !
     open(unit=11,file=trim(trim(mesh_dir)//trim('spin_spinor_unique_fine_')//(trim(ibnd_loc))//trim('.dat')))
     !
     do ikx=-nk1,nk1
        do iky=-nk2,nk2
           !
           kpoint(1)=1.d0*ikx/nk1
           kpoint(2)=1.d0*iky/nk2
           !
           write(11,'(100f16.8)') kpoint(1),kpoint(2),(spin_m(ikx,iky,ibnd,ipol),ipol=1,4)
           !
        enddo !iky
     enddo !ikx
     !
     close(11)
     !
  enddo !ibnd
  
end program spin_texture_spinor
