program ep_melements

  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use intw_intw2wannier
  use intw_symmetries
  use intw_w90
  use w90_io, only: io_error,stdout,io_file_unit,seedname,io_time,io_stopwatch
  use w90_parameters, only: num_wann,bands_num_points,real_metric,recip_metric,&
                                            bands_num_spec_points,timing_level, &
                               bands_spec_points,bands_label,bands_plot_format, &
                   bands_plot_mode,num_bands_project,bands_plot_project,num_bands
  use w90_hamiltonian, only: irvec,nrpts,ndegen,ham_r
  use intw_allwfcs

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

  integer                  :: nbands_loc

  complex(dp), allocatable :: wfc_k (:,:,:) ! nG_max is defined in reading
  complex(dp), allocatable :: wfc_kq (:,:,:)
  complex(dp), allocatable :: wfc_k_aux (:,:,:)
  complex(dp), allocatable :: wfc_k_orig (:,:,:)

  real(dp), allocatable    :: QE_eig_k(:)
  real(dp), allocatable    :: QE_eig_kq(:)

  !fft related
  integer                  :: ngk, ngkq, nr(3)
  integer                  :: G(1:3), GKQ_bz(3), G_plusk(3), G_pluskq(3)

  integer                  :: i, j, k
  integer                  :: ii, jj, kk
  integer                  :: jr, kr
  integer                  :: ig, ibnd, jbnd, ipol, jpol
  integer                  :: switch
  logical                  :: read_status, have_nnkp
  character(256)           :: nnkp_file, method
  character(len=4)         :: iq_loc

  complex(dp),allocatable :: wfc_k_r(:)
  integer :: igg1, igg2, gg1(3), gg2(3), nG, ikb
  logical :: found
  complex(dp),allocatable :: auxil(:,:,:)
  integer,allocatable :: list_aux(:)


  complex(dp),allocatable :: u_mesh_irr(:,:,:), uu(:,:)

  integer :: nd
  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(*,20) '====================================================='
  write(*,20) '|                  program read_wan                 |'
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
  
  call allocate_and_read_W90()
  !
  ! extract the nband x num_wann Wannier projection+rotation matrices.
  !
  call produce_u_mesh()
  !
  !================================================================================
  !       Allocate wfc related files
  !================================================================================
  !
  nbands_loc=num_bands

  
  allocate(u_mesh_irr(num_bands,num_wann,num_kpts))

  allocate(uu(num_wann,num_wann)) 
  u_mesh_irr=(0.0d0,0.0d0)
  nd=0

  print*,"hemen nago", num_kpts, nk_irr

  do ik=1,num_kpts
    if (QE_folder_sym(ik)/=18) then 
        cycle
    end if

    print*, ik
    nd=nd+1
    u_mesh_irr(:,:,18) = u_mesh_irr(:,:,18) + u_mesh(:,:,ik)

    uu=matmul(transpose(conjg(u_mesh(:,:,ik))),u_mesh(:,:,ik))
    do i=1, num_bands
      write(unit=44,fmt="( 100(2f12.6,x) )") ( u_mesh(i,j,ik), j=1, num_wann )
    end do
    write(unit=44,fmt="(a)") "----------"

    do i=1, num_wann
      write(unit=45,fmt="( 100(2f12.6,x) )") ( uu(i,j), j=1, num_wann )
    end do
    write(unit=45,fmt="(a)") "----------"

  end do!ik

    write(unit=44,fmt="(a)") " mean----------"

    do i=1, num_bands
      write(unit=44,fmt="( 100(2f12.6,x) )") ( u_mesh_irr(i,j,18)/nd, j=1, num_wann )
    end do
    write(unit=44,fmt="(a)") "----------"

    uu=matmul(transpose(conjg(u_mesh_irr(:,:,18))),u_mesh_irr(:,:,18))
 
   write(unit=44,fmt="(a)") " uu ---------"

    do i=1, num_wann
      write(unit=44,fmt="( 100(2f12.6,x) )") ( uu(i,j)/nd**2, j=1, num_wann )
    end do
    write(unit=44,fmt="(a)") "----------"
 
    

   
stop
  
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
  !
  !================================================================================
  call get_timing(time2)
  write(*,20) '|                     ALL DONE                       |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='


20 format(A)
30 format(A,F8.2,6X,A)

contains

  subroutine calculate_rotation_unitary_matrices()

  complex(dp)  :: unit_mat(nbands_loc, nbands_loc)

  unit_mat = cmplx_0

  do ibnd=1,nbands_loc
    unit_mat(ibnd, ibnd) = cmplx_1
  enddo


    do ik=1,nkmesh
!    do ik=1,1

       kpoint(:) = kmesh(:,ik)

       call get_psi_general_k_all_wfc(.true., kpoint, list_iGk_orig , wfc_k_orig ,  QE_eig_k,  G_plusk)

       do isym=1,nsym

          kpoint_rot  =matmul(s(:,:, isym ), kpoint(:) )

          call find_k_1BZ_and_G(kpoint_rot,nk1, nk2, nk3, i ,j, k, kpoint_in_1bz, GKQ_bz)
          call switch_indices (nk1, nk2, nk3, ikpt_k, i, j, k, +1)

          call get_psi_general_k_all_wfc(.true., kpoint_in_1bz, list_iGk , wfc_k ,  QE_eig_k,  G_plusk)

          list_iGk_aux = list_iGk_orig
          wfc_k_aux    = wfc_k_orig

          G_pluskq=nint(kmesh(:,ikpt_k)-kpoint_rot)

          call rotate_wfc_test (wfc_k_aux,list_iGk_aux, wfc_kq, list_iGkq,         &
               isym, s(:,:,isym),  ftau(:,isym) , (/0,0,0/))

          call get_plane_wave_matrix_element_convolution  ((G_pluskq-G_plusk)*0, list_iGk, list_iGkq, wfc_k, wfc_kq , &
               unit_sym_sgk( ik, isym, 1:nbands_loc,1:nbands_loc)  )

          umat=  unit_sym_sgk( ikpt_k, isym, 1:nbands_loc,1:nbands_loc)

          write(12334,*)ik,isym
          do i=1,nbands_loc
             write(12334,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,nbands_loc)
          enddo

          !Test. Determinant of Det[u^{-1}.u] must be = 1.
          umat=matmul(umat,conjg(transpose(umat)))

          write(10000,*)ik,isym
          do i=1,nbands_loc
             write(10000,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,nbands_loc)
          enddo


          !call diagonalize_cmat (nbands_loc, umat, QE_eig_k )
          !suma_c=cmplx_1
          !do i=1,nbands_loc
          !   suma_c=suma_c*QE_eig_k(i)
          !enddo
          !if (abs(suma_c-cmplx_1)>1E-3) then
             !write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
             !stop
          !end if

          if (sum(abs(unit_mat - umat))>0.001) then
             write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
             write(*,*) ik, isym
             write(*,*) sum(abs(unit_mat - umat))
             write(10000,*) '*****ERROREA*****'
             write(4325,*) ik, isym
!             stop
          end if

       enddo ! isym

    end do !ik

  end subroutine calculate_rotation_unitary_matrices

  subroutine calculate_rotation_unitary_matrices_wannier()

  complex(dp)  :: unit_mat(nbands_loc, nbands_loc)
  complex(dp)  :: wfc_k_W(nG_max,nbands_loc,nspin)
  complex(dp)  :: wfc_kq_W(nG_max,nbands_loc,nspin)
  complex(dp)  :: wfc_k_orig_W(nG_max,nbands_loc,nspin)
  complex(dp)  :: U_k(nbands_loc,num_wann)

  unit_mat = cmplx_0

  do ibnd=1,num_wann
    unit_mat(ibnd, ibnd) = cmplx_1
  enddo


    do ik=1,nkmesh

       kpoint(:) = kmesh(:,ik)

       call get_psi_general_k_all_wfc(.true., kpoint, list_iGk_orig , wfc_k_orig ,  QE_eig_k,  G_plusk)

       call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3, i ,j, k, kpoint_in_1bz, GKQ_bz)
       call switch_indices (nk1, nk2, nk3, ikpt_k, i, j, k, +1)

       U_k(:,:)=u_mesh(:,:,ikpt_k)
       call wfc_bands_to_wann(wfc_k_orig,U_k,wfc_k_orig_W)

       do isym=1,nsym

          kpoint_rot  =matmul(s(:,:, isym ), kpoint(:) )

          call find_k_1BZ_and_G(kpoint_rot,nk1, nk2, nk3, i ,j, k, kpoint_in_1bz, GKQ_bz)
          call switch_indices (nk1, nk2, nk3, ikpt_k, i, j, k, +1)

          call get_psi_general_k_all_wfc(.true., kpoint_in_1bz, list_iGk , wfc_k ,  QE_eig_k,  G_plusk)

          U_k(:,:)=u_mesh(:,:,ikpt_k)
          call wfc_bands_to_wann(wfc_k,U_k,wfc_k_W)

          list_iGk_aux = list_iGk_orig
          wfc_k_aux    = wfc_k_orig_W

          G_pluskq=nint(kmesh(:,ikpt_k)-kpoint_rot)

          call rotate_wfc_test (wfc_k_aux,list_iGk_aux, wfc_kq_W, list_iGkq,         &
               isym, s(:,:,isym),  ftau(:,isym) , (/0,0,0/))

          do i=num_wann+1,nbands_loc
             wfc_kq_W(:,i,:)=(0.d0,0.d0)
             wfc_k_W(:,i,:)=(0.d0,0.d0)
          enddo

          call get_plane_wave_matrix_element_convolution  ((G_pluskq-G_plusk)*0, list_iGk, list_iGkq, wfc_k_W, wfc_kq_W , &
               unit_sym_sgk( ik, isym, 1:num_wann,1:num_wann)  )

!          umat=  unit_sym_sgk( ikpt_k, isym, 1:nbands_loc,1:nbands_loc)
          umat=  unit_sym_sgk( ik, isym, 1:num_wann,1:num_wann)

          write(12334,*)ik,isym
          do i=1,num_wann
             write(12334,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,num_wann)
          enddo

          !Test. Determinant of Det[u^{-1}.u] must be = 1.
          umat=matmul(umat,conjg(transpose(umat)))

          write(10000,*)ik,isym
          do i=1,num_wann
             write(10000,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,num_wann)
          enddo


          !call diagonalize_cmat (nbands_loc, umat, QE_eig_k )
          !suma_c=cmplx_1
          !do i=1,nbands_loc
          !   suma_c=suma_c*QE_eig_k(i)
          !enddo
          !if (abs(suma_c-cmplx_1)>1E-3) then
             !write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
             !stop
          !end if

          if (sum(abs(unit_mat - umat))>0.001) then
             write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
             write(*,*) ik, isym
             write(*,*) sum(abs(unit_mat - umat))
             write(10000,*) '*****ERROREA*****'
             write(4326,*) ik, isym
!             stop
          end if

       enddo ! isym

    end do !ik

  end subroutine calculate_rotation_unitary_matrices_wannier

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
