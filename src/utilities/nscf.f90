program nscf

  use kinds, only: dp
  use intw_input_parameters, only: intw2W_method, read_input
  use intw_reading, only: nspin, noncolin, ngm, nsym, &
       spinorb_mag, can_use_TR, s, nbands, nG_max, nat, alat, &
       ntyp, amass, nspin, tau, bg, nr1, nr2, nr3, nbands, nkpoints_QE, &
       get_gvec, &
       read_parameters_data_file_xml, &
       read_kpoints_data_file_xml, &
       get_K_folder_data
  use intw_pseudo, only: read_all_pseudo
  use intw_utility, only: get_timing, find_free_unit, switch_indices

  use intw_useful_constants, only: cmplx_0, cmplx_1

  use intw_fft, only: generate_nl, allocate_fft, nl

  use intw_allwfcs, only: allocate_and_get_all_irreducible_wfc
  use intw_allwfcs, only: wfc_k_irr_all, list_iG_all, QE_eig_irr_all


  !================================================================================
  !       Declare the variables
  !================================================================================
  implicit none

  integer, allocatable     :: list_igk (:,:), ngk(:)

  complex(dp), allocatable :: wfc_k (:,:,:,:) ! nG_max is defined in reading
  real(dp), allocatable    :: QE_eig_k(:,:)
  !fft related
  integer                  :: nr(3)

  !local/aux variables
  integer                  :: nbands_loc
  integer                  :: npw

  integer                  :: i, j, k, ik
  integer                  :: ig, ibnd, jbnd, ipol, jpol
  logical                  :: read_status
  character(256)           :: method

  complex(dp),allocatable :: wfc_k_r(:)
  integer :: nG
  complex(dp), external :: zdotc
  real(kind=dp) :: time1, time2

  real(kind=dp), allocatable            :: kpoints_QE(:,:)

  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  call get_timing(time1)
  write(*,20) '====================================================='
  write(*,20) '|                  program nscf                     |'
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
  write(*,20) '|       - reading pseudopotentials from KBPP files  |'
  write(*,20) '|            (defined in .save data files)          |'
  write(*,20) '|                                                   |'
  !
  call read_all_pseudo ()
  !if (.not.lspinorb) call average_pp(ntyp)
  !
  write(*,20) '|                    PPs are OK                     |'
  write(*,20) '|           ---------------------------------       |'
  !
  allocate (list_igk(nG_max,nkpoints_QE))
  allocate (ngk(nkpoints_QE))
  !
  allocate (wfc_k(nG_max,nbands,nspin,nkpoints_QE))
  !
  allocate (wfc_k_r(nr1*nr2*nr3))
  !
  allocate (QE_eig_k(nbands,nkpoints_QE))
  !
  !================================================================================
  !       read in the kpoints from the QE folders
  !================================================================================
  !
  allocate(kpoints_QE(3,nkpoints_QE))
  !
  call read_kpoints_data_file_xml(kpoints_QE)
  write(*,20) '|      k point list (crystal)                       |'
  do ik=1, nkpoints_QE
     write(*,"(a,i4,3f12.6,a)")'|',ik,kpoints_QE(:,i),'           |'
     call get_K_folder_data(ik,list_iGk(:,ik),wfc_k(:,:,:,ik),QE_eig_k(:,ik),ngk(ik))
     write(*,20) '|      Energies are:                                |'
     !write(*,"(5f10.2)")QE_eig_k
     write(*,*)maxval(list_iGk(:,ik)), "ngk", ngk(ik)
  end do

  write(*,*)maxval(list_iGk(:,:))

  !call allocate_and_get_all_irreducible_wfc()
  !
  !call deallocate_upfeak ()
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
    !--------------------------------------------------------------------------------------------------
  subroutine convolution_of_two_functions (list_iG_1,ngk1,list_iG_2,ngk2,wfc_1,wfc_2, product_wfc)
    !--------------------------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    !	Given two wavefunctions wfc_1  and wfc_2, this subroutine computes
    !              the convolution of them as:
    !
    !	          pw_mat_el_ij (G)    \sum_GP  wfc_1_i(G-GP) * wfc_2_j(GP)
    !
    !     The G-vectors are referenced by their indices in list_iG1, list_iG2,
    !     which refer to the global list gvec(3,ngm), which should be defined
    !     BEFORE using this subroutine....ç
    !
    !
    !--------------------------------------------------------------------------------
    use intw_useful_constants, only: zero, one, cmplx_0
    use intw_reading, only: nG_max, gvec, nspin, nbands, ngm

    implicit none

    !I/O variables in (F fortran style). For example,
    !Input list_iG_1(:), instead of list_iG_1(nG_max).

    integer,intent(in)      :: list_iG_1(:),ngk1,list_iG_2(:),ngk2
    complex(dp),intent(in)  :: wfc_1(:,:,:), wfc_2(:,:,:) !wfc_1(nG_max,num_bands,nspin),wfc_2(nG_max,num_bands,nspin)

    ! In output, we have nbndxnbnd functions in (G), but G in the full ngm list

    complex(dp),intent(out) :: product_wfc (:,:,:,:,:) !(num_bands,num_bands,nspin,nspin,ngm)

    !local variables
    integer :: nbnd_l
    integer :: nG_max_l
    integer :: nspin_l

    integer :: G2pG1(3)
    integer :: i,j,ibnd,jbnd,iG_1,iG_2, iG
    integer :: nG_max_non_zero
    !complex(dp) :: pw_mat_el_local(num_bands,num_bands,nspin,nspin)
    logical :: found

    nG_max_l  = size(wfc_1,1) !
    nbnd_l    = size(wfc_1,2)
    nspin_l    = size(wfc_1,3)

    product_wfc =  (0.d0,0.0d0)

    do iG_1=1,nGk1    !This is G-Gp
       do iG_2=1,nGk2 !This is Gp

          G2pG1 = gvec(:,list_iG_1(iG_1)) + gvec(:,list_iG_2(iG_2)) ! This is G
          call find_iG(G2pG1,iG) !  This is the index for G

          do ibnd=1, nbnd_l
             do jbnd=1, nbnd_l
                do ipol=1,nspin_l
                   do jpol=1,nspin_l

                      product_wfc (ibnd,jbnd,ipol,jpol,iG) = product_wfc (ibnd,jbnd,ipol,jpol,iG) + &
                           conjg(wfc_2 (iG_2,jbnd,jpol)) * wfc_1 (iG_1,ibnd,ipol)

                   end do !jpol
                end do !ipol
             end do !jbnd
          end do !ibnd

       end do !iG_2
    end do !iG_1

  end subroutine convolution_of_two_functions

  subroutine diagonalize_cmat (n,a,w)

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


end program nscf
