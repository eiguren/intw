!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       program intw2W90 
!       -----------------------------------
!
!       This is a "utility" program which is part of the intw project.
!
!       The purpose of this utility is to perform the same tasks as the
!       program "pw2wannier90" which is part of the QE distribution, but
!       utilizing a minimum set of (QE generated) Bloch functions, using
!       symmetry.
!
!       The code is heavily inspired from pw2wannier90. 
!       Mad props to the QE people; their code is a bit above my paygrade,
!       however, so I've taken some bits and modified them so that a monkey
!       could understand.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program intw2W90 

  use intw_useful_constants
  use intw_intw2wannier
  use intw_tests
  use intw_symmetries
  use intw_fft
  use intw_input_parameters
  use w90_parameters, only: num_bands,num_exclude_bands,exclude_bands

!================================================================================
!       Declare the variables 
!================================================================================
  implicit none

  integer :: ikpt, i_folder, i, j, k, ik, ibnd, switch
  integer :: nk_irr, nkmesh 
  logical :: read_status, have_nnkp, k_points_consistent
  character(256) :: nnkp_file, method
  real(dp),allocatable :: kpoints_QE(:,:), kmesh(:,:)
  real(dp) :: time1, time2
  real(dp),allocatable :: kpoints_irr(:,:) 
  

  call get_timing(time1)
  
!================================================================================
!       Talk to the user
!================================================================================
  
  write(*,20) '====================================================='
  write(*,20) '|                  program intw2W90                 |'
  write(*,20) '|        ---------------------------------          |'
  write(*,20) '====================================================='
  write(*,20) '|    waiting for input file...                      |'
  
!================================================================================
!       read the input file 
!       Read in the necessary information from standard input
!================================================================================
  
  call read_input(read_status)

  !
  method = trim(intw2W_method)
  !
  if (read_status) then
     !
     stop 
     !
  endif        
  
!================================================================================
!       read the parameters from the SCF QE calculation, in  
!       particular, read in the symmetry matrices!!! 
!================================================================================
 
  call read_parameters_data_file_xml()

!-call test_symmetry_axis_angle()

!================================================================================
!      set up the gvec array, which will contain the global
!      G-vectors, as well as the FFT code, which is necessary to
!      generate g_fft_map, which in turn is necessary in the 
!      wavefunction rotation code!
!================================================================================

  call get_ngm()

  allocate(gvec(3,ngm)) ! it is deallocated in the bottom of this main code.

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

!================================================================================
!       Check that $prefix.nnkp is present
!================================================================================

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

!================================================================================
!       read the parameters in the .nnkp file 
!================================================================================

  call read_nnkp_file(nnkp_file)
  !
  num_bands=nbands-nnkp_exclude_bands
  !
  ! just as a test; can be removed later
  !
  call output_nnkp_file()

!================================================================================
!       read in the kpoints from the QE folders
!================================================================================

  allocate(kpoints_QE(3,nkpoints_QE))
  call read_kpoints_data_file_xml(kpoints_QE) 

!================================================================================
! Build the kmesh corresponding to the parameters in the input file 
!================================================================================

  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)

!================================================================================
! check that kmesh and nnkp_kpoints are consistent with one another 
!       This insures that the Wannier data is consistent with intw data.
!================================================================================

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

!================================================================================
!    Find the size of the irreducible set of k-points (IBZ)
!    and check that the number of kpoints corresponds to either
!    a full mesh or the IBZ.
!================================================================================

  call find_size_of_irreducible_k_set(nk1,nk2,nk3,kmesh,nk_irr)
  !
  !This is only for testing: The result for nk_irr is different in both.
  !
  allocate(kpoints_irr(3,nk1*nk2*nk3))
  call find_the_irreducible_k_set(nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
  deallocate(kpoints_irr )
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

!================================================================================
!      allocate the symmetry arrays 
!      CAREFUL! the subroutine needs to know the global value of "full_mesh",
!      so it is crucial that this allocation occurs AFTER setting full_mesh
!================================================================================

  call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)

!================================================================================
!     Compute the indices of the inverse rotation matrices.
!     This will be useful later in the exectution.
!================================================================================

  call find_inverse_symmetry_matrices_indices()

!================================================================================
!     Set up the array spin_symmetry_matrices, which contain
!     the matrices which must be used to rotate spin
!================================================================================

  call allocate_and_build_spin_symmetry_matrices(nsym)

!================================================================================
!      Fill the symmetry arrays appropriately
!================================================================================

!Peio
  if (spinorb_mag) then
     can_use_TR=.true.
  endif
!Peio

  call set_symmetry_relations(nk1,nk2,nk3,nkpoints_QE,kpoints_QE,kmesh, &
                         k_points_consistent,QE_folder_sym,sym_G,symlink)


!================================================================================
!       Tell the user what is in the QE folders
!================================================================================

  if (full_mesh .and. IBZ) then
     !
     write(*,20) '|       - the kpoints present in the QE folders     |'
     write(*,20) '|         are consistent with a full 1BZ and a      |'
     write(*,20) '|         IBZ has also been found.                  |'
     write(*,20) '|           ---------------------------------       |'
     !
  elseif (IBZ) then
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

!================================================================================
!       Check that the requested calculation is possible
!================================================================================

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

!================================================================================
!       Compute the mmn file      
!================================================================================

  if (compute_mmn) then
     !
     write(*,20) '|       - writing the file prefix.eig and           |'
     write(*,20) '|         computing the file prefix.mmn...          |'
     write(*,20) '| (this is labor intensive and may  take some time) |'
     write(*,20) '|           ---------------------------------       |'
     !
     call generate_mmn(intw2W_fullzone,method)
     !
  endif

!================================================================================
!       Compute the amn file      
!================================================================================

  if (compute_amn) then
     !
     write(*,20) '|       - computing the file prefix.amn...          |'
     write(*,20) '| (this is labor intensive and may  take some time) |'
     write(*,20) '|           ---------------------------------       |'
     !
     call generate_amn(intw2W_fullzone,method)
     !
  endif

!================================================================================
!       clean up 
!================================================================================

  call deallocate_symmetry_related_k()
  call deallocate_nnkp()
  call deallocate_fft()
  call deallocate_reading_variables()
  call deallocate_spin_symmetry_matrices()
  !
  deallocate(kpoints_QE)
  deallocate (gvec)
  deallocate(kmesh)

!================================================================================
!       Finish 
!================================================================================

  call get_timing(time2)
  !
  write(*,20) '|                     ALL DONE                      |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='
  !
  20 format(A)
  30 format(A,F8.2,6X,A)

end program intw2W90 
