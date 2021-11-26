!----------------------------------------------------------------------------!
!  intw project.
!
!----------------------------------------------------------------------------!
!
module intw_W90
!
!----------------------------------------------------------------------------!
!
!       The subroutines in this module call subroutines from the W90
!  distribution in order to READ the prefix.chk file, which contains
!  the U matrices.
!
!----------------------------------------------------------------------------!
  !
  use kinds, only: dp
  !
  !
  implicit none
  !
  save
  !
  ! variables
  public :: u_mesh, excluded_bands
  !
  !subroutines
  public :: allocate_and_read_W90, &
            deallocate_W90, &
            produce_u_mesh, &
            interpolate_u_matrices_and_eigenvalues, &
            interpolate_u_matrices_and_eigenvalues_proj, &
            Wrotate_plane_wave_matrix_elements, &
            Wrotate_matrix_elements_no_spin, &
            convert_ham_r_to_H_R, &
            interpolate_u_matrices_and_eigenvalues2, &
            compute_derivative_of_Hk, &
            compute_iR_H_R, &
            remove_useless_W90, &
            deallocate_hr_W90
  !
  private
  !
  complex(dp), allocatable :: u_mesh(:,:,:)

  logical(dp), allocatable :: excluded_bands(:)
  !
  !
contains
!
!--------------------------------------
  subroutine allocate_and_read_W90()
    !----------------------------------------------------------------------------!
    ! Call *obscure* routines from the W90 distribution, which will read
    ! all the info in $SEED.chk, including the U matrices.
    !
    ! Let seedname contain the path of the wannier directory (which is
    ! also the mesh directory). This way, no need to locally copy
    ! large $SEED.chk, $SEED.eig etc files to computation folder!
    !
    ! NOTE: I (Bruno) have modified Wannier90 so that the string "seedname" is
    ! large enough to accomodate a path. This is because I don't understand
    ! the W90 routines, and they require seedname to lead to the file.
    !----------------------------------------------------------------------------!
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit
    use intw_reading, only: nbands
    use w90_parameters, only: param_read, param_read_chkpt, num_exclude_bands, &
                              exclude_bands, num_kpts, num_bands, eigval, num_wann
    use w90_io, only: seedname
    use w90_hamiltonian, only: hamiltonian_setup, hamiltonian_get_hr
    use w90_kmesh, only: kmesh_get

    implicit none

    integer :: io_unit, k, n, i, j

    seedname = trim(mesh_dir)//trim(prefix)
    !
    ! read the parameters from wannier90
    !
    call param_read()
    call kmesh_get()
    call param_read_chkpt()
    !
    ! define the excluded bands from the .chk file
    !
    allocate(excluded_bands(nbands))
    excluded_bands(:)=.false.
    !
    do i=1,num_exclude_bands
      !
      n=exclude_bands(i)
      excluded_bands(n)=.true.
      !
    enddo
    !
    ! Explicitely read the eigenvalues from $SEED.eig (which must be present!)
    ! this code snippet was taken from the W90 distribution and slightly
    ! modified.
    !
    io_unit=find_free_unit()
    !
    open(unit=io_unit,file=trim(seedname)//'.eig',form='formatted',status='old')
    !
    do k=1,num_kpts
      do n=1,num_bands
        !
        read(io_unit,*) i, j, eigval(n,k)
        !
      enddo
    enddo
    !
    ! set up the hamiltonian in R space by calling W90 routines
    !
    call hamiltonian_setup()
    call hamiltonian_get_hr()
    !
    ! remove allocated variables which are of no use to intw
    !
    call remove_useless_W90()
    !
    allocate(u_mesh(num_bands,num_wann,num_kpts))
    !
    return

  end subroutine allocate_and_read_W90
!-------------------------------------------
!*******************************************
!-------------------------------------------
  subroutine deallocate_W90()
    !----------------------------------------------------------------------------!
    ! This subroutine uses deallocation functions from the W90 distribution
    ! to clean up the W90 environment
    !----------------------------------------------------------------------------!
    use w90_parameters, only: param_dealloc
    use w90_kmesh, only: kmesh_dealloc
    use w90_overlap, only: overlap_dealloc

    implicit none

    call kmesh_dealloc()
    call param_dealloc()
    call overlap_dealloc()
    !
    deallocate(u_mesh)
    !
    return

  end subroutine deallocate_W90
!---------------------------------
!*********************************
!---------------------------------
  subroutine produce_u_mesh()
    !----------------------------------------------------------------------------!
    ! This subroutine computes the "rotation" matrices which must be applied
    ! to the matrix elements to take them from band space to Wannier space.
    !
    ! If disentanglement is used, u_mesh    = u_matrix_opt x u_matrix
    ! If no disentanglement is used, u_mesh = u_matrix
    !
    ! note that the number of k-points in the coarse mesh, in Wannier90,
    ! is defined in the variable num_kpts. Consistency demands that this number
    ! be equal to nkmesh, the same quantity in the intw program.
    !----------------------------------------------------------------------------!
    use intw_reading, only: nbands
    use w90_parameters, only: u_matrix_opt, u_matrix, num_kpts, num_wann

    implicit none

    integer :: i, n1, n2, nb1, ikpt


    u_mesh = (0.0_dp,0.0_dp)
    !
    if (allocated(u_matrix_opt)) then
      !
      do ikpt=1,num_kpts
        !
        n1 = 0
        !
        do nb1=1,nbands
          !
          if (excluded_bands(nb1)) cycle
          !
          n1 = n1 + 1
          !
          do n2=1,num_wann
            do i=1,num_wann
              !
              u_mesh(n1,n2,ikpt) = u_mesh(n1,n2,ikpt) + &
                                   u_matrix_opt(n1,i,ikpt)*u_matrix(i,n2,ikpt)
                   !
            enddo !i
          enddo !n2
        enddo !nb1
        !
      enddo !ikpt
      !
    else
      !
      u_mesh = u_matrix
      !
    endif
    !
    ! deallocate these arrays from W90, which are no longer useful
    !
    if (allocated(u_matrix_opt)) deallocate(u_matrix_opt)
    if (allocated(u_matrix)) deallocate(u_matrix)
    !
    return

  end subroutine produce_u_mesh
!-----------------------------------------------------------------------------------
!***********************************************************************************
!-----------------------------------------------------------------------------------
  subroutine interpolate_u_matrices_and_eigenvalues(u_int, eig_int, kvec_int, H_int, nk_max, nk)
  !----------------------------------------------------------------------!
  !                                                                      !
  !   This  subroutine is strongly based on "plot_interpolate_bands"     !
  !   from the W90 distribution.                                         !
  !   This subroutine has been modified to suit the needs of the intw    !
  !   project. The goal of this subroutine is to compute the interpolated!
  !   eigenvalues and u_matrices for k points not on the coarse mesh.    !
  !                                                                      !
  !   CAREFUL!                                                           !
  !   It appears that the eigenvectors of H are stored as the columns of !
  !   the "u" array; this is not the order expected in the Wannier       !
  !   formalism! This subroutine thus transpose the u matrix.            !
  !----------------------------------------------------------------------!
  use intw_useful_constants, only: cmplx_0, cmplx_i, twopi
  use w90_io, only: io_error, io_file_unit, io_time, io_stopwatch
  use w90_parameters, only: num_wann
  use w90_hamiltonian, only: irvec, nrpts, ndegen, ham_r

  implicit none

  external :: ZHPEVX

  integer :: nk_max ! the size of the kpoint dimension in arrays
  integer :: nk ! the number of meaningful kpoints in the arrays
  integer :: info, i, j
  integer :: loop_rpt, nfound, loop_kpt
  integer :: ierr

  integer, allocatable :: iwork(:), ifail(:)

  complex(kind=dp), allocatable :: ham_pack(:)
  complex(kind=dp), allocatable :: ham_kprm(:,:)
  complex(kind=dp), allocatable :: cwork(:)

  complex(kind=dp) :: fac
  complex(kind=dp) :: u_dagger(num_wann,num_wann)
  complex(kind=dp) :: u_int(num_wann,num_wann,nk_max)
  complex(kind=dp) :: H_int(num_wann,num_wann,nk_max)


  real(kind=dp), allocatable :: rwork(:)

  real(kind=dp) :: eig_int(num_wann,nk_max)
  real(kind=dp) :: kvec_int(3,nk_max)
  real(kind=dp) :: rdotk

  !
  !
  !
  allocate(ham_pack((num_wann*(num_wann+1))/2),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating ham_pack in plot_interpolate_bands')
  allocate(ham_kprm(num_wann,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating ham_kprm in plot_interpolate_bands')
  allocate(cwork(2*num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating cwork in plot_interpolate_bands')
  allocate(rwork(7*num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating rwork in plot_interpolate_bands')
  allocate(iwork(5*num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating iwork in plot_interpolate_bands')
  allocate(ifail(num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating ifail in plot_interpolate_bands')


  eig_int=0.0_dp
  !
  !
  ! Interpolate the Hamiltonian at each meaningful kpoint
  !

  do loop_kpt=1,nk

     ham_kprm=cmplx_0

     do loop_rpt=1,nrpts
       rdotk    = twopi*dot_product(kvec_int(:,loop_kpt),irvec(:,loop_rpt))
       fac      = exp(cmplx_i*rdotk)/real(ndegen(loop_rpt),dp)
       ham_kprm = ham_kprm+fac*ham_r(:,:,loop_rpt)
     enddo

     H_int(:,:,loop_kpt) = ham_kprm

     ! Diagonalise H_k (->basis of eigenstates)
     do j=1,num_wann
       do i=1,j
         ham_pack(i+((j-1)*j)/2)=ham_kprm(i,j)
       enddo
     enddo

     ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
     call ZHPEVX( 'V', 'A', 'U', num_wann, ham_pack, &
                  0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                  nfound, eig_int(:,loop_kpt), u_dagger, &
                  num_wann, cwork, rwork, iwork, ifail, info )

        u_int(:,:,loop_kpt) = conjg(transpose( u_dagger))

     if(info < 0) then
       write(*,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
       call io_error('Error in plot_interpolate_bands')
     endif

     if(info > 0) then
       write(*,'(i3,a)') info,' EIGENVECTORS FAILED TO CONVERGE'
       call io_error('Error in "interpolate_u_matrices_and_eigenvalues"')
     endif
  !
  enddo


  end subroutine interpolate_u_matrices_and_eigenvalues

  subroutine interpolate_u_matrices_and_eigenvalues_proj(numz, u_int, eig_int, kvec_int, H_int, nk_max, nk)
  !----------------------------------------------------------------------!
  !                                                                      !
  !   This  subroutine is strongly based on "plot_interpolate_bands"     !
  !   from the W90 distribution.                                         !
  !   This subroutine has been modified to suit the needs of the intw    !
  !   project. The goal of this subroutine is to compute the interpolated!
  !   eigenvalues and u_matrices for k points not on the coarse mesh.    !
  !                                                                      !
  !   CAREFUL!                                                           !
  !   It appears that the eigenvectors of H are stored as the columns of !
  !   the "u" array; this is not the order expected in the Wannier       !
  !   formalism! This subroutine thus transpose the u matrix.            !
  !----------------------------------------------------------------------!
  use intw_useful_constants, only: cmplx_i, cmplx_0, twopi
  use w90_io, only: io_error, io_file_unit, io_time, io_stopwatch
  use w90_parameters, only: num_wann
  use w90_hamiltonian, only: irvec, nrpts, ndegen, ham_r

  implicit none

  external :: ZHPEVX

  integer :: nk_max, numz ! the size of the kpoint dimension in arrays
  integer :: nk ! the number of meaningful kpoints in the arrays
  integer :: info, i, j
  integer :: loop_rpt, nfound, loop_kpt
  integer :: ierr

  integer, allocatable :: iwork(:), ifail(:)

  complex(kind=dp), allocatable :: ham_pack(:)
  complex(kind=dp), allocatable :: ham_kprm(:,:)
  complex(kind=dp), allocatable :: cwork(:)

  complex(kind=dp) :: fac
  complex(kind=dp) :: u_dagger(num_wann,num_wann)
  complex(kind=dp) :: u_int(num_wann,num_wann,nk_max)
  complex(kind=dp) :: H_int(num_wann,num_wann,nk_max)


  real(kind=dp), allocatable :: rwork(:)

  real(kind=dp) :: eig_int(num_wann,nk_max), eig_int_all(num_wann*numz,nk_max)
  real(kind=dp) :: kvec_int(3,nk_max)
  real(kind=dp) :: rdotk

  integer :: ikz, ibnd
  !
  !
  !
  allocate(ham_pack((num_wann*(num_wann+1))/2),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating ham_pack in plot_interpolate_bands')
  allocate(ham_kprm(num_wann,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating ham_kprm in plot_interpolate_bands')
  allocate(cwork(2*num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating cwork in plot_interpolate_bands')
  allocate(rwork(7*num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating rwork in plot_interpolate_bands')
  allocate(iwork(5*num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating iwork in plot_interpolate_bands')
  allocate(ifail(num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating ifail in plot_interpolate_bands')


  eig_int=0.0_dp
  !
  !
  ! Interpolate the Hamiltonian at each meaningful kpoint
  !
  do ikz=1,numz

  do loop_kpt=1,nk

     ham_kprm=cmplx_0

     do loop_rpt=1,nrpts
       rdotk    = twopi*dot_product(kvec_int(:,loop_kpt),irvec(:,loop_rpt))
       fac      = exp(cmplx_i*rdotk)/real(ndegen(loop_rpt),dp)
       ham_kprm = ham_kprm+fac*ham_r(:,:,loop_rpt)
     enddo

     H_int(:,:,loop_kpt) = ham_kprm

     ! Diagonalise H_k (->basis of eigenstates)
     do j=1,num_wann
       do i=1,j
         ham_pack(i+((j-1)*j)/2)=ham_kprm(i,j)
       enddo
     enddo

     ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
     call ZHPEVX( 'V', 'A', 'U', num_wann, ham_pack, &
                  0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                  nfound, eig_int(:,loop_kpt), u_dagger, &
                  num_wann, cwork, rwork, iwork, ifail, info )

        do ibnd=1, num_wann
          eig_int_all(ibnd+(ikz-1)*numz,loop_kpt) = eig_int_all(ibnd,loop_kpt)
        enddo

        u_int(:,:,loop_kpt) = conjg(transpose( u_dagger))

     if(info < 0) then
       write(*,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
       call io_error('Error in plot_interpolate_bands')
     endif

     if(info > 0) then
       write(*,'(i3,a)') info,' EIGENVECTORS FAILED TO CONVERGE'
       call io_error('Error in "interpolate_u_matrices_and_eigenvalues"')
     endif
  !
  enddo

  enddo !ikz

  do loop_kpt=1,nk
     write(112,'(10000f12.6)') (eig_int_all(ibnd,loop_kpt), ibnd=1,num_wann*numz)
  enddo


  end subroutine interpolate_u_matrices_and_eigenvalues_proj
!*******************************************************************************
!-------------------------------------------------------------------------------
  subroutine Wrotate_plane_wave_matrix_elements(ikpt_1, ikpt_2, u_mesh, pw_matrix_elements, W_rot_mat_el)
    !----------------------------------------------------------------------!
    !   This  subroutine W-rotate the matrix elements with respect to both
    !   k1 and k2 indices, using the appropriate u_matrices.
    !
    !   input:  pw_matrix_elements(nbands,nbands)
    !             given by  < n1 k1 | e^{-i (q+G) r} | n2 k2>
    !             with ikpt_1 the singlet index of k1,
    !             and ikpt_2 the singlet index of k2.
    !
    !           u_mesh(num_bands,num_wann,num_kpts)
    !              num_wann <= nbands: if disentanglement is present,
    !              num_wann < nbands.  In this case, the matrices are
    !              rectangular.
    !
    !   output: W_rot_mat_el(num_wann,num_wann)
    !              the Wannier rotated (and projected) matrix elements.
    !----------------------------------------------------------------------!
    use intw_useful_constants, only: cmplx_0
    use intw_reading, only: npol
    use w90_parameters, only: num_bands, num_wann, num_kpts

    implicit none

    integer, intent(in) :: ikpt_1, ikpt_2 ! the singlet index of k1 and k2
    complex(dp), intent(in) :: pw_matrix_elements(num_bands,num_bands,npol,npol)
    complex(dp), intent(in) :: u_mesh(num_bands,num_wann,num_kpts)
    complex(dp), intent(out) :: W_rot_mat_el(num_wann,num_wann,npol,npol)

    integer :: ipol, jpol ! spin indices.
    integer :: nb1, nb2 ! bands loop indices
    integer :: nbW1, nbW2 ! wannier projection loop indices
    complex(dp) :: u1(num_bands,num_wann), u2(num_bands,num_wann)

    ! fish out the matrices
    u1 = u_mesh(:,:,ikpt_1)
    u2 = u_mesh(:,:,ikpt_2)
    !
    ! initialize
    W_rot_mat_el = cmplx_0
    !
    !       PROJECT + ROTATE on bands
    !       Let's be careful about the order of the loops to avoid
    !    "cache misses". Remember: In fortran, the first index
    !    rolls faster.
    !
    do jpol=1,npol
      do ipol=1,npol
        do nbW2 = 1,num_wann
          do nbW1 = 1,num_wann
            !
            do nb2=1,num_bands
              do nb1=1,num_bands
                !
                W_rot_mat_el(nbW1,nbW2,ipol,jpol)=W_rot_mat_el(nbW1,nbW2,ipol,jpol)+ &
                                    conjg(u1(nb1,nbW1))*pw_matrix_elements(nb1,nb2,ipol,jpol)*u2(nb2,nbW2)
                !
              enddo !nb1
            enddo !nb2
            !
          enddo ! nbW1
        enddo  !nbW2
      enddo !ipol
    enddo !jpol

  end subroutine Wrotate_plane_wave_matrix_elements
!----------------------------------------------------------------------------------------------------------
!**********************************************************************************************************
!----------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine Wrotate_matrix_elements_no_spin(ikpt_1, ikpt_2, u_mesh, pw_matrix_elements, W_rot_mat_el)
!-------------------------------------------------------------------------------
!
!----------------------------------------------------------------------!
!   This  subroutine W-rotate the matrix elements with respect to both
!   k1 and k2 indices, using the appropriate u_matrices.
!
!   input:  pw_matrix_elements(nbands,nbands)
!             given by  < n1 k1 | e^{-i (q+G) r} | n2 k2>
!             with ikpt_1 the singlet index of k1,
!             and ikpt_2 the singlet index of k2.
!
!           u_mesh(num_bands,num_wann,num_kpts)
!              num_wann <= nbands: if disentanglement is present,
!              num_wann < nbands.  In this case, the matrices are
!              rectangular.
!
!   output: W_rot_mat_el(num_wann,num_wann)
!              the Wannier rotated (and projected) matrix elements.
!----------------------------------------------------------------------!
  use intw_useful_constants, only: cmplx_0
  use w90_parameters, only: num_bands, num_wann, num_kpts

  implicit none

  integer, intent(in) :: ikpt_1, ikpt_2 ! the singlet index of k1 and k2
  complex(dp), intent(in) :: pw_matrix_elements(num_bands,num_bands)
  complex(dp), intent(in) :: u_mesh(num_bands,num_wann,num_kpts)
  complex(dp), intent(out) :: W_rot_mat_el(num_wann,num_wann)

  integer :: nb1, nb2 ! bands loop indices
  integer :: nbW1, nbW2 ! wannier projection loop indices
  complex(dp) :: u1(num_bands,num_wann), u2(num_bands,num_wann)

  ! fish out the matrices
  u1 = u_mesh(:,:,ikpt_1)
  u2 = u_mesh(:,:,ikpt_2)
  !
  ! initialize
  W_rot_mat_el = cmplx_0
  !
  !       PROJECT + ROTATE on bands
  !       Let's be careful about the order of the loops to avoid
  !    "cache misses". Remember: In fortran, the first index
  !    rolls faster.
  !
  do nbW2 = 1,num_wann
     do nbW1 = 1,num_wann
        !
        do nb2=1,num_bands
           do nb1=1,num_bands
              !
              W_rot_mat_el(nbW1,nbW2)=W_rot_mat_el(nbW1,nbW2)+ &
                            conjg(u1(nb1,nbW1))*pw_matrix_elements(nb1,nb2)*u2(nb2,nbW2)
              !
           enddo !nb1
        enddo !nb2
        !
     enddo ! nbW1
  enddo  !nbW2

  end subroutine Wrotate_matrix_elements_no_spin
!-----------------------------------------------------------------------------------
!***********************************************************************************
!-----------------------------------------------------------------------------------

  subroutine convert_ham_r_to_H_R(num_pack, H_R)
  !----------------------------------------------------------------------!
  !                                                                      !
  !   This  subroutine converts the real-space Hamiltonian (force        !
  !   constants) from the original Wannier90 format to a more palatable  !
  !   (and useful) format for FFT interpolation. Since the Hamiltonian   !
  !   is Hermitian, there is absolutely no reason to not store it as an  !
  !   upper triangular matrix (this is what "pack" will mean).           !
  !----------------------------------------------------------------------!
  use intw_input_parameters, only: nk1, nk2, nk3
  use intw_useful_constants, only: cmplx_0
  use w90_hamiltonian, only: irvec, nrpts, ham_r
  use w90_parameters, only: num_wann

  implicit none

  !input variables
  integer :: num_pack
  complex(dp) :: H_R(nk1,nk2,nk3,num_pack)

  ! local variables
  logical :: R_is_found(nk1,nk2,nk3)
  integer :: ir_loop, ir1, ir2, ir3
  integer :: m_bnd, n_bnd, i_pack

  ! initialize
  H_R(:,:,:,:)      = cmplx_0
  R_is_found(:,:,:) = .false.

  ! loop on all R vectors
  do ir_loop = 1, nrpts ! this is a Wannier90 variable, not necessarily = nk1*nk2*nk3

     ! The indices of the R vector in Wannier90 are centered at Gamma
     ! return to off-centered indices
     ir1 = modulo(irvec(1,ir_loop),nk1)+1
     ir2 = modulo(irvec(2,ir_loop),nk2)+1
     ir3 = modulo(irvec(3,ir_loop),nk3)+1

  ! if this point hasn't already been found, process it
     if ( .not.  R_is_found(ir1,ir2,ir3) ) then

    R_is_found(ir1,ir2,ir3) = .true.
    ! upper triangular packing

          do n_bnd =1,num_wann
              do m_bnd = 1, n_bnd

                i_pack = m_bnd+((n_bnd-1)*n_bnd)/2

                H_R(ir1,ir2,ir3,i_pack) = ham_r(m_bnd,n_bnd,ir_loop)

              enddo
           enddo

  endif

  enddo

  end subroutine convert_ham_r_to_H_R

  subroutine interpolate_u_matrices_and_eigenvalues2(nk_vec_max, nk_vec, num_pack, H, u, eig)
  !----------------------------------------------------------------------!
  !                                                                      !
  !   This  subroutine, given Hamiltonians in upper triangular form,     !
  !   computes the eigenvalues and the eigenvectors (stored in u         !
  !   matrices).
  !                                                                      !
  !   CAREFUL!                                                           !
  !   It appears that the eigenvectors of H are stored as the columns of !
  !   the "u" array; this is not the order expected in the Wannier       !
  !   formalism! This subroutine thus transposes the u matrix.           !
  !----------------------------------------------------------------------!
  use intw_useful_constants, only: cmplx_0, ZERO
  use w90_parameters, only: num_wann

  implicit none

  external :: ZHPEV

  ! input variables
  integer :: nk_vec_max ! the size of the kpoint dimension in arrays
  integer :: nk_vec ! the number of meaningful kpoints in the  arrays
  integer :: num_pack ! array size, given upper triangular storage

  complex(dp) :: H   (num_pack,nk_vec_max)

  ! output variables
  complex(dp) :: u   (num_wann,num_wann,nk_vec_max)
  real(dp) :: eig (num_wann,nk_vec_max)


  ! local variables
  integer :: info
  integer :: i_kpt

  complex(dp) :: u_dagger(num_wann,num_wann)

  ! work arrays for the diagonalization routine

  real(dp) :: rwork(3*num_wann-2)
  complex(dp) :: cwork(2*num_wann-1)

  ! initialize
  eig  = ZERO
  u    = cmplx_0

  do i_kpt=1,nk_vec

     ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
     !call ZHPEVX(    'V',    'A',    'U',    num_wann,       ham_pack,       &
     !             0.0_dp, 0.0_dp,      0,           0,        -1.0_dp,       &
     !             nfound,    eig_int(:,loop_kpt),                            &
     !                        u_dagger           ,                            &
     !           num_wann,  cwork,  rwork,       iwork,          ifail, info)

     ! Diagonalizing routine from mkl. See manual.
     call zhpev('V', 'U', num_wann, H(:,i_kpt), &
                eig(:,i_kpt), u_dagger, num_wann, &
                cwork, rwork, info )


        u(:,:,i_kpt) = conjg(transpose(u_dagger))

     if(info < 0) then
       write(*,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEV HAD AN ILLEGAL VALUE'
       write(*,*) 'Error in "interpolate_u_matrices_and_eigenvalues2"'
     endif

     if(info > 0) then
       write(*,'(i3,a)') info,' EIGENVECTORS FAILED TO CONVERGE'
       write(*,*) 'Error in "interpolate_u_matrices_and_eigenvalues2"'
     endif
  enddo ! i_kpt

  end subroutine interpolate_u_matrices_and_eigenvalues2


  subroutine compute_derivative_of_Hk(kvec, dHk)
  !----------------------------------------------------------------------------!
  !     This subroutine computes  d/dk_alpha H(k) for a given k
  !     using DFT. This routine will be for testing.
  !----------------------------------------------------------------------------!
  use intw_useful_constants, only: cmplx_0, cmplx_i, twopi
  use w90_parameters, only: num_wann
  use w90_hamiltonian, only: irvec, nrpts, ham_r, ndegen

  implicit none

  ! input variables
  real(dp) :: kvec(3) ! in crystal coordinates


  ! output variable
  complex(dp) :: dHk(num_wann,num_wann,3)

  ! local variables
  integer :: loop_rpt
  real(dp) :: rdotk
  complex(dp) :: fac

  dHk(:,:,:) = cmplx_0

      write(45,'(A,3F20.12)') ' kvec = ', kvec(:)
  do loop_rpt=1,nrpts

      rdotk    = twopi*dot_product(kvec(:),irvec(:,loop_rpt))
      fac      = exp(cmplx_i*rdotk)/real(ndegen(loop_rpt),dp)

      dHk(:,:,1) = dHk(:,:,1) + cmplx_i*irvec(1,loop_rpt)*fac*ham_r(:,:,loop_rpt)
      dHk(:,:,2) = dHk(:,:,2) + cmplx_i*irvec(2,loop_rpt)*fac*ham_r(:,:,loop_rpt)
      dHk(:,:,3) = dHk(:,:,3) + cmplx_i*irvec(3,loop_rpt)*fac*ham_r(:,:,loop_rpt)

      write(45,'(I8,8X,3I8)') ndegen(loop_rpt), irvec(:,loop_rpt)
  enddo



  end subroutine compute_derivative_of_Hk


  subroutine compute_iR_H_R(num_pack, iR_H_R)
    !----------------------------------------------------------------------!
    !                                                                      !
    !   This  subroutine computes the real space matrices iR_alpha H(R),   !
    !   which will be useful to compute the band derivatives by FFT.       !
    !   This in turn will be used to implement an adaptive broadening      !
    !   scheme.                                                            !
    !----------------------------------------------------------------------!
    use intw_input_parameters, only: nk1, nk2, nk3
    use intw_useful_constants, only: cmplx_0, cmplx_i
    use w90_hamiltonian, only: irvec, nrpts, ham_r
    use w90_parameters, only: num_wann

    implicit none

    !input variables
    integer :: num_pack


    !output variables
    complex(dp) :: iR_H_R(nk1,nk2,nk3,3,num_pack)

    ! local variables
    logical :: R_is_found(nk1,nk2,nk3)
    integer :: ir_loop, ir1, ir2, ir3
    integer :: m_bnd, n_bnd, i_pack
    integer :: alpha

    ! initialize
    iR_H_R(:,:,:,:,:) = cmplx_0
    R_is_found(:,:,:) = .false.

    ! loop on all R vectors
    do ir_loop = 1, nrpts ! this is a Wannier90 variable, not necessarily = nk1*nk2*nk3

      ! The indices of the R vector in Wannier90 are centered at Gamma
      ! return to off-centered indices
      ir1 = modulo(irvec(1,ir_loop),nk1)+1
      ir2 = modulo(irvec(2,ir_loop),nk2)+1
      ir3 = modulo(irvec(3,ir_loop),nk3)+1

      ! if this point hasn't already been found, process it
       if ( .not.  R_is_found(ir1,ir2,ir3) ) then

        R_is_found(ir1,ir2,ir3) = .true.
        ! upper triangular packing

        do n_bnd =1,num_wann
          do m_bnd = 1, n_bnd

            i_pack = m_bnd+((n_bnd-1)*n_bnd)/2

            do alpha = 1, 3
              iR_H_R(ir1,ir2,ir3,alpha,i_pack) = cmplx_i*irvec(alpha,ir_loop)*ham_r(m_bnd,n_bnd,ir_loop)
            enddo

          enddo
        enddo

      endif

    enddo


  end subroutine compute_iR_H_R

  subroutine remove_useless_W90()
  !----------------------------------------------------------------------------!
  ! This subroutine deallocates variables that were allocated by the W90
  ! subroutines, but are not useful to intw.
  !----------------------------------------------------------------------------!
  use w90_parameters, only: proj_zona, proj_site, proj_radial, proj_x, proj_z, &
                            proj_l, proj_m, atoms_symbol, atoms_label, atoms_species_num, &
                            atoms_pos_cart, atoms_pos_frac, bands_spec_points, &
                            bands_label, bands_plot_project, nncell, nnlist, &
                            kpt_latt, kpt_cart, wb, neigh, lwindow, ndimwin, &
                            bka, bk, shell_list, wannier_spreads, wannier_centres, &
                            m_matrix, m_matrix_orig, a_matrix, eigval, wannier_plot_list

  implicit none

!  if (allocated(exclude_bands)) then
!      deallocate(exclude_bands)
!  endif
  if (allocated(wannier_plot_list)) then
      deallocate(wannier_plot_list)
  endif
  if (allocated(bands_plot_project)) then
      deallocate(bands_plot_project)
  endif
  if (allocated(shell_list)) then
      deallocate(shell_list)
  endif
  if (allocated(kpt_latt)) then
      deallocate(kpt_latt)
  endif
  if (allocated(atoms_pos_frac)) then
      deallocate(atoms_pos_frac)
  endif
  if (allocated(atoms_pos_cart)) then
      deallocate(atoms_pos_cart)
  endif
  if (allocated(atoms_species_num)) then
      deallocate(atoms_species_num)
  endif
  if (allocated(atoms_label)) then
      deallocate(atoms_label)
  endif
  if (allocated(atoms_symbol)) then
      deallocate(atoms_symbol)
  endif
  if (allocated(proj_site)) then
      deallocate(proj_site)
  endif
  if (allocated(proj_l)) then
      deallocate(proj_l)
  endif
  if (allocated(proj_m)) then
      deallocate(proj_m)
  endif
  if (allocated(proj_z)) then
      deallocate(proj_z)
  endif
  if (allocated(proj_x)) then
      deallocate(proj_x)
  endif
  if (allocated(proj_radial)) then
      deallocate(proj_radial)
  endif
  if (allocated(proj_zona)) then
      deallocate(proj_zona)
  endif
  if (allocated(bands_label)) then
      deallocate(bands_label)
  endif
  if (allocated(bands_spec_points)) then
      deallocate(bands_spec_points)
  endif
  if (allocated(kpt_cart)) then
      deallocate(kpt_cart)
  endif
  if (allocated(nnlist)) then
      deallocate(nnlist)
  endif
  if (allocated(neigh)) then
      deallocate(neigh)
  endif
  if (allocated(nncell)) then
      deallocate(nncell)
  endif
  if (allocated(wb)) then
      deallocate(wb)
  endif
  if (allocated(bk)) then
      deallocate(bk)
  endif
  if (allocated(bka)) then
      deallocate(bka)
  endif
  if (allocated(ndimwin)) then
      deallocate(ndimwin)
  endif
  if (allocated(lwindow)) then
      deallocate(lwindow)
  endif
  if (allocated(a_matrix)) then
      deallocate(a_matrix)
  endif
  if (allocated(m_matrix_orig)) then
      deallocate(m_matrix_orig)
  endif
  if (allocated(eigval)) then
      deallocate(eigval)
  endif
  if (allocated(m_matrix)) then
      deallocate(m_matrix)
  endif
  if (allocated(wannier_centres)) then
      deallocate(wannier_centres)
  endif
  if (allocated(wannier_spreads)) then
      deallocate(wannier_spreads)
  endif

  end subroutine remove_useless_W90

  subroutine deallocate_hr_W90()
  !----------------------------------------------------------------------------!
  ! This subroutine deallocates the W90 h(R) hamiltonian, which is
  ! no longer useful past a certain point.
  !----------------------------------------------------------------------------!
  use w90_parameters
  use w90_io
  use w90_hamiltonian
  use w90_kmesh
  use w90_disentangle
  use w90_overlap
  use w90_wannierise
  use w90_plot
  use w90_transport
  implicit none

  if (allocated(ham_r)) then
      deallocate(ham_r)
  endif

  end subroutine deallocate_hr_W90

!----------------------------------------------------------------------------!
!
end module intw_W90
!
!----------------------------------------------------------------------------!
