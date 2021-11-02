!----------------------------------------------------------------------------!
!       intw project.
!
!----------------------------------------------------------------------------!
module intw_plane_wave_setup
!----------------------------------------------------------------------------!
!
!       This module contains subroutines which will set up the interpolation
!       of plane wave matrix elements, of the type
!
!               < n1 k | e^{-i (q+G)*r} | n2 k+q >.
!
!       For a given q, and for fixed n1,n2,G, the above can be thought of
!       as a function of (k,k+q), where k lives in the 1BZ. Its domain is
!       thus a 3D subspace of the more general 6D space (k1,k2). In order       
!       to use hexacubic or hexalinear interpolation, it will be necessary
!       to identify and properly tabulate the relevant "supercubes".
!
!----------------------------------------------------------------------------!

use intw_useful_constants

use iotk_module
use intw_utility
use intw_input_parameters
use intw_new_interpolate
use intw_W90
use intw_matrix_elements
use intw_symmetries

!Peio
!Number of Bloch Original states for the wannierization
use w90_parameters, only: num_bands
!Peio

  !
  implicit none
  !
  save
  !

contains

  subroutine build_sub_blocks(ikpt,qpt,sub_blocks,occupied_blocks)
  !----------------------------------------------------------------------------!
  !    Consider all 6D vectors of the form (k1,k2), where both k1 and k2
  !    belong to the coarse mesh. Each vector forms the identifying corner
  !    of a supercube in 6D space. 
  !
  !    For a given k point in the coarse mesh, identified by ikpt, a sub-set
  !    of the fine mesh {ks} can be defined as all the ks points which fit inside
  !    the 3D cube defined by the k point.
  !
  !    The image of this subset upon addition of the qpoint, qpt, {ks +qpt}
  !    will in general overlap up to 8 coarse mesh cubes. The direct product
  !    of  {ks} and {ks+qpt} thus forms an object which lives in 6 dimensions;
  !    in general, this object will live inside up to 8 6D supercubes.
  ! 
  !    This subroutine computes the 6D coordinates of these touched supercubes,
  !    and tabulates which ones are actually occupied. These coordinates will
  !    be of the form
  !                i such that ki = (i-1)/nki for each dimension
  !    However, ki will not be restrained to be in the 1BZ.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: ikpt
  real(dp)       :: qpt(3)
  logical        :: occupied_blocks(8)
  integer        :: sub_blocks(8,6)

  ! local variables
  real(dp)       :: kq(3)

  integer        :: i1, i2, i3, i4, i5, i6
  integer        :: i1s, i2s, i3s
  integer        :: icube

  logical        :: already_found
  integer        :: ikpts
  integer        :: norm
  integer        :: iloop, jloop

  integer        :: switch1, switch2 
  integer        :: list_corners_ikpts(8)

  switch1 = -1   ! singlet to triplet
  switch2 =  1   ! triplet to singlet


  ! initialize occupied_blocks
  occupied_blocks  =  .false.
  
  ! find the coordinates of ikpt 
  call switch_indices(nk1,nk2,nk3,ikpt,i1,i2,i3,switch1)

  ! find the 8 corners in the smooth mesh inside the coarse cube ikpt
  call find_submesh_corners(ikpt,list_corners_ikpts)

  icube = 0
  ! loop on the corners in the smooth mesh; see if they land in 
  ! distinct cubes
  do iloop = 1, 8 
     ikpts = list_corners_ikpts(iloop)

     ! find k+q
     call switch_indices(nk1s,nk2s,nk3s,ikpts,i1s,i2s,i3s,switch1)

     kq(1) = dble(i1s-1)/dble(nk1s)+qpt(1)
     kq(2) = dble(i2s-1)/dble(nk2s)+qpt(2)
     kq(3) = dble(i3s-1)/dble(nk3s)+qpt(3)

     ! find the corner in the coarse mesh associated with this k-point
     i4    = floor(dble(nk1)*kq(1)+1.0_dp)
     i5    = floor(dble(nk2)*kq(2)+1.0_dp)
     i6    = floor(dble(nk3)*kq(3)+1.0_dp)

     ! see if it is already accounted for
     already_found = .false.
     do jloop = 1,icube    
                norm =   (i4-sub_blocks(jloop,4))**2+ &
                         (i5-sub_blocks(jloop,5))**2+ &
                         (i6-sub_blocks(jloop,6))**2  

                if (norm == 0) already_found = .true.

     end do

     if ( .not. already_found) then
        icube = icube + 1 
        occupied_blocks(icube) = .true.

        sub_blocks(icube,:)   = (/ i1,i2,i3,i4,i5,i6 /)
     end if
  end do

  end subroutine build_sub_blocks


  subroutine find_is_min_is_max(nk,nks,i,is_min,is_max)
  !-----------------------------------------------------
  ! convenience routine
  !-----------------------------------------------------
  integer        :: iwork,  nk,  nks,  i,  is_min, is_max

  iwork = nks*(i-1)+nk
  if ( modulo(iwork,nk) == 0 ) then
        is_min = iwork/nk
  else
        is_min = iwork/nk+1
  end if

  iwork = iwork+nks
  if ( modulo(iwork,nk) == 0 ) then
        is_max = iwork/nk-1
  else
        is_max = iwork/nk
  end if

  
  end subroutine find_is_min_is_max


  subroutine find_submesh_corners(ikpt,list_corners_ikpts)
  !----------------------------------------------------------------------------!
  ! Given a kpoint in the coarse mesh, identified by ikpt, which defines a
  ! cube in the 1BZ lattice, returns the list of indices of the kpoints in 
  ! the fine mesh which form the corners of the sub-fine mesh which fits 
  ! in that cube.
  !
  ! The algorithm only manipulates integers to avoid error. The basic idea
  ! is:
  !             ks = (is-1)/nks,    k = (i-1)/nk
  !
  !               (i-1)/nk <=  (is-1)/nks <  i/nk
  !
  !                     nks*(i-1) + nk <= nk*is
  !                     nks*  i   + nk  > nk*is
  !
  !  Find is_min and is_max by being careful with the integers.
  !----------------------------------------------------------------------------!
  implicit none

 
  integer        :: switch 
  integer        :: i1, i2, i3, i1s, i2s, i3s
  integer        :: i1s_min, i2s_min, i3s_min
  integer        :: i1s_max, i2s_max, i3s_max

  integer        :: ikpt, ikpts

  integer        :: list_corners_ikpts(8)

  integer        :: list_i1s(2), list_i2s(2), list_i3s(2)

  integer        :: iloop, index_i1s, index_i2s, index_i3s


  ! identify the 8 corner points in the fine mesh which define the
  ! the cube which lands inside the cube defined by ikpt.

  switch = -1   ! singlet to triplet
  call switch_indices(nk1,nk2,nk3,ikpt,i1,i2,i3,switch)

  call find_is_min_is_max(nk1,nk1s,i1,i1s_min,i1s_max)
  call find_is_min_is_max(nk2,nk2s,i2,i2s_min,i2s_max)
  call find_is_min_is_max(nk3,nk3s,i3,i3s_min,i3s_max)

  list_i1s = (/i1s_min,i1s_max/)
  list_i2s = (/i2s_min,i2s_max/)
  list_i3s = (/i3s_min,i3s_max/)

  switch = 1   ! triplet to singlet

  iloop = 0
  do index_i3s = 1,2
     do index_i2s = 1,2
        do index_i1s = 1,2

           iloop = iloop+1
           i1s   = list_i1s(index_i1s)
           i2s   = list_i2s(index_i2s)
           i3s   = list_i3s(index_i3s)

           ! find the points in the smooth mesh
           call switch_indices(nk1s,nk2s,nk3s,ikpts,i1s,i2s,i3s,switch)

           list_corners_ikpts(iloop) = ikpts

        end do
     end do
  end do

  end subroutine find_submesh_corners

  subroutine compute_matrix_elements_at_l_scalar       &
             (N_scheme,D_scheme,sze_scheme,use_IBZ,    &
              l_scalar,X_l,W_rotated_matrix_elements,  &
              pw_matrix_elements,wfc_1,wfc_2,QE_eig,   &
              list_iG_1,list_iG_2)

  !------------------------------------------------------------------
  ! 
  ! This subroutine computes all the plane wave matrix elements
  ! corresponding to a given scalar index l_scalar.
  !
  !------------------------------------------------------------------

  implicit none

  integer        ::      N_scheme,D_scheme,sze_scheme
  integer        ::      X_l(sze_scheme,D_scheme)
  logical        ::      use_IBZ
  integer        ::      l_scalar

  integer        ::      i1, i2, i3, i4, i5, i6

  integer        ::      G1(3), G2(3), G(3)

  integer        ::      iG
  integer        ::      ikpt1, ikpt2

  integer        ::      switch

  complex(dp)    :: W_rotated_matrix_elements(nG_shell_max,num_wann,num_wann,nspin,nspin)
  complex(dp)    :: pw_matrix_elements(num_bands,num_bands,nspin,nspin)
!  complex(dp)    :: pw_matrix_elements(nbands,nbands,nspin,nspin)

  complex(dp)    :: wfc_1(nG_max,num_bands,nspin), wfc_2(nG_max,num_bands,nspin)
!  complex(dp)    :: wfc_1(nG_max,nbands,nspin), wfc_2(nG_max,nbands,nspin)

  real(dp)       :: QE_eig(num_bands)
!  real(dp)       :: QE_eig(nbands)

  integer        :: list_iG_1(nG_max), list_iG_2(nG_max) 

  switch  =  1 ! triplet to singlet
  W_rotated_matrix_elements = cmplx_0

  ! find k1 and k2, as well as the appropriate G vectors
  i1 = modulo(X_l(l_scalar,1)-1,nk1)+1
  i2 = modulo(X_l(l_scalar,2)-1,nk2)+1
  i3 = modulo(X_l(l_scalar,3)-1,nk3)+1

  G1(1) = (X_l(l_scalar,1)-i1)/nk1
  G1(2) = (X_l(l_scalar,2)-i2)/nk2
  G1(3) = (X_l(l_scalar,3)-i3)/nk3

  call switch_indices(nk1,nk2,nk3,ikpt1,i1,i2,i3,switch)

  i4 = modulo(X_l(l_scalar,4)-1,nk1)+1
  i5 = modulo(X_l(l_scalar,5)-1,nk2)+1
  i6 = modulo(X_l(l_scalar,6)-1,nk3)+1

  G2(1) = (X_l(l_scalar,4)-i4)/nk1
  G2(2) = (X_l(l_scalar,5)-i5)/nk2
  G2(3) = (X_l(l_scalar,6)-i6)/nk3

  call switch_indices(nk1,nk2,nk3,ikpt2,i4,i5,i6,switch)

  ! get the wavefunctions
  call get_psi_k(ikpt1,use_IBZ,list_iG_1,wfc_1,QE_eig)
  call get_psi_k(ikpt2,use_IBZ,list_iG_2,wfc_2,QE_eig)

  ! loop on G vectors
  do iG = 1,nG_shell_max

        ! identify the G vector that goes in the matrix element
        G = gvec(:,iG) +G2(:)-G1(:)
        ! compute the plane wave matrix elements

!	ORIGINAL CONVOLUTION ROUTINE: 
!		apparently faster than the FFT routines (even FFTW)
!		this is because my tests only involve 1 G vector:
!		FFT may get the advantage for many G shells.
!        call get_plane_wave_matrix_element_convolution      &
!                (G,list_iG_1,list_iG_2, wfc_1,wfc_2,        &
!                                 pw_matrix_elements)

!	IMPROVED CONVOLUTION ROUTINE: 
!	this is similar in spirit to the original convolution routine,
!	but utilizes mkl matrix-matrix products. Is this faster?
!	PRELIMINARY profiling suggests that this is significantly
!	faster. MAKE SURE THERE IS NO BUGS!
        call get_plane_wave_matrix_element_convolution_alt  &
                (G,list_iG_1,list_iG_2, wfc_1,wfc_2,        &
                                 pw_matrix_elements)

!	FFT ROUTINE: 
!	does not currently produce functional code, probably
!	due to bad linking to the FFT libraries (gives seg faults)
!       call get_plane_wave_matrix_element_FFT              &
!                (G,list_iG_1,list_iG_2, wfc_1,wfc_2,       &
!                                 pw_matrix_elements)

        ! W rotate
        ! Asier && Idoia 15 07 2014: We add nspin variable as input below. 
        call Wrotate_plane_wave_matrix_elements(ikpt1,ikpt2, &
                                  u_mesh,pw_matrix_elements, &
                                W_rotated_matrix_elements(iG,:,:,:,:))

  end do

  end subroutine compute_matrix_elements_at_l_scalar

  subroutine read_write_checkpoint_xml(N_scheme,D_scheme,sze_scheme,X_l, what)
  !------------------------------------------------------------------
  ! 
  ! This subroutine writes a checkpoint file, in xml format, 
  ! which allows the program to keep track of "where it is".
  !
  ! This file will be intimately linked to a direct-access binary file
  ! which will contain the actual data. This xml file will allow the
  ! the program to check for consistency.
  !
  ! The file should contain the (N+1)^D D dimensional integer vectors
  ! which identify the points in the 6D coarse mesh for which 
  ! the matrix elements were computed. This information will be
  ! very useful in order to recycle data from one sub block to the next.
  ! It should also contain a minimum set of data to allow consistency check.
  ! 
  !------------------------------------------------------------------
  implicit none

  integer        :: N_scheme, D_scheme, sze_scheme

  integer        :: N_scheme_read, D_scheme_read, sze_scheme_read
  integer        :: nk1_read, nk2_read, nk3_read, nG_shell_max_read

  integer        :: X_l(sze_scheme,D_scheme), X_l_read(D_scheme,sze_scheme)


  character(*),parameter :: filename='checkpoint.xml'
  character(*)   :: what

  logical        :: consistency_test

  integer        :: io_unit

  io_unit = find_free_unit()

  if ( what .eq. 'write' ) then

     call iotk_open_write (io_unit,filename,binary=.false.)
         call iotk_write_comment(io_unit,                      &
                'This file contains information about the '    &
                        //'data stored in the binary file.')

         call iotk_write_dat  (io_unit,'D_scheme',D_scheme)
         call iotk_write_dat  (io_unit,'N_scheme',N_scheme)
         call iotk_write_dat  (io_unit,'nk1',nk1)
         call iotk_write_dat  (io_unit,'nk2',nk2)
         call iotk_write_dat  (io_unit,'nk3',nk3)
         call iotk_write_dat  (io_unit,'nG_shell_max',nG_shell_max)
!         call iotk_write_dat  (io_unit,'X_l',X_l,columns=sze_scheme)

         call iotk_write_dat  (io_unit,'X_l',transpose(X_l),columns=D_scheme)

       call iotk_close_write (io_unit)

   else if ( what .eq. 'read' ) then

     call iotk_open_read (io_unit,filename,binary=.false.)

        call iotk_scan_dat (io_unit,'D_scheme', D_scheme_read)
        call iotk_scan_dat (io_unit,'N_scheme', N_scheme_read)
        call iotk_scan_dat (io_unit,'nk1', nk1_read)
        call iotk_scan_dat (io_unit,'nk2', nk2_read)
        call iotk_scan_dat (io_unit,'nk3', nk3_read)
        call iotk_scan_dat (io_unit,'nG_shell_max',nG_shell_max_read)

        consistency_test  = (nk1_read == nk1)          .and.    &     
                            (nk2_read == nk2)          .and.    &     
                            (nk3_read == nk3)          .and.    &     
                       (D_scheme_read == D_scheme )    .and.    &     
                       (N_scheme_read == N_scheme )    .and.    &     
                   (nG_shell_max_read == nG_shell_max)

        if ( .not. consistency_test ) then
                write(*,*) '******************************************'        

                write(*,*) '* The data in the checkpoint.xml file    *'        
                write(*,*) '* is not consistent with internal state; *'        
                write(*,*) '* the computation cannot continue.       *'        
                write(*,*) '******************************************'        
                stop
        end if                


        call iotk_scan_dat  (io_unit,'X_l',X_l_read)
        X_l = transpose(X_l_read)

     call iotk_close_read(io_unit)

   end if

  end subroutine read_write_checkpoint_xml


  subroutine get_record_index(l_scalar,iG,record_index)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes the record index, given that the 
  ! data is specified by two labels: l_scalar and iG. The same
  ! old scheme will be used to go from one to the other.
  !------------------------------------------------------------------
  implicit none
  integer        :: l_scalar, iG, record_index
  
  ! This scheme puts all the G shells together for a given l.      
  ! This choice doesn't matter since the file is direct access.

  record_index = (l_scalar-1)*nG_shell_max+iG  
  
  end subroutine get_record_index


  subroutine create_rotated_matrix_elements_file &
                (N_scheme,D_scheme,sze_scheme,X_l,use_IBZ)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes all the plane wave matrix elements
  ! necessary for interpolation for a specified sub-block of the coarse
  ! mesh, and saves the information in the binary file.
  !
  !------------------------------------------------------------------

  implicit none

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio

  integer        ::      N_scheme,D_scheme,sze_scheme
  integer        ::      X_l(sze_scheme,D_scheme)
  logical        ::      use_IBZ

  integer        ::      l_scalar, iG

  ! variables to get record index
  integer        ::      record_index
  integer        ::      record_length
  integer        ::      io_unit_bin

  complex(dp),allocatable  :: W_rotated_matrix_elements(:,:,:,:,:)
  complex(dp),allocatable  :: pw_matrix_elements(:,:,:,:)
  complex(dp),allocatable  :: wfc_1(:,:,:), wfc_2(:,:,:)

  real(dp),allocatable     :: QE_eig(:)

  integer,allocatable      :: list_iG_1(:), list_iG_2(:) 

  character(*),parameter  :: filename = 'W_matrix_elements.bin'

  nbands_loc=num_bands

  io_unit_bin = find_free_unit()
  ! open a file
  record_length = num_wann*num_wann* nspin**2 *direct_io_factor_cmplx 

  open(unit = io_unit_bin, file = filename, form = 'unformatted', &
     status = 'unknown', access = 'direct', recl = record_length, &
     action = 'readwrite')

  allocate(W_rotated_matrix_elements(nG_shell_max,num_wann,num_wann,nspin,nspin))
  allocate(pw_matrix_elements(nbands_loc,nbands_loc,nspin,nspin))
!  allocate(pw_matrix_elements(nbands,nbands,nspin,nspin))
  allocate(wfc_1(nG_max,nbands_loc,nspin))
!  allocate(wfc_1(nG_max,nbands,nspin))
  allocate(wfc_2(nG_max,nbands_loc,nspin))
!  allocate(wfc_2(nG_max,nbands,nspin))
  allocate(list_iG_1(nG_max))
  allocate(list_iG_2(nG_max))
  allocate(QE_eig(nbands_loc))
!  allocate(QE_eig(nbands))
  
  ! TEST
!  write(36,*) '# rotated matrix elements'
!  write(36,*) '# l  Xl(6)  W_pw(real) W_pw(im)'

  do l_scalar = 1, sze_scheme

     ! compute the matrix elements
     call compute_matrix_elements_at_l_scalar          &
             (N_scheme,D_scheme,sze_scheme,use_IBZ,    &
              l_scalar,X_l,W_rotated_matrix_elements,  &
              pw_matrix_elements,wfc_1,wfc_2,QE_eig,   &
              list_iG_1,list_iG_2)


     do iG=1,nG_shell_max
        ! get the record index
        call get_record_index(l_scalar,iG,record_index)
        ! write to file!
        write(io_unit_bin,rec=record_index) W_rotated_matrix_elements(iG,:,:,:,:)
!        print *, 'W = ', sqrt(sum(abs(W_rotated_matrix_elements(iG,:,:))**2))
     end do

!     write(36,'(I4,6I6,2F12.8)')  l_scalar, X_l(l_scalar,:), W_rotated_matrix_elements(1,1,1)
  end do

  

  deallocate(pw_matrix_elements)
  deallocate(W_rotated_matrix_elements)
  deallocate(wfc_1)
  deallocate(wfc_2)
  deallocate(list_iG_1)
  deallocate(list_iG_2)
  deallocate(QE_eig)

  close(io_unit_bin)

  end subroutine create_rotated_matrix_elements_file


  subroutine update_rotated_matrix_elements_file &
                (N_scheme,D_scheme,sze_scheme,X_l,use_IBZ,n_found, &
                 reading_time, writing_time, computing_time)

  !------------------------------------------------------------------
  ! 
  ! This subroutine produces the binary file necessary for interpolation,
  ! assuming such a file already exists for a previous block. 
  ! The subroutine recycles as much data as possible and only 
  ! computes matrix elements for the new points. Once the new 
  ! binary file is obtained, it is moved to the original file.
  !
  ! Clearly, it would be better to modify the file "in place";
  ! this however will take time to implement.
  !------------------------------------------------------------------

  implicit none

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio

  integer        ::      N_scheme,D_scheme,sze_scheme
  integer        ::      X_l(sze_scheme,D_scheme)
  logical        ::      use_IBZ

  integer        ::      X_l_file(sze_scheme,D_scheme)

  logical        ::      found

  integer        ::      n_found

  integer        ::      l_scalar, iG, l_scalar_file

  integer        ::      io_unit_old_bin, io_unit_new_bin
  integer        ::      record_index, record_length

  complex(dp),allocatable  :: W_rotated_matrix_elements(:,:,:,:,:)
  complex(dp),allocatable  :: pw_matrix_elements(:,:,:,:)
  complex(dp),allocatable  :: wfc_1(:,:,:), wfc_2(:,:,:)

  real(dp),allocatable     :: QE_eig(:)

  integer,allocatable      :: list_iG_1(:), list_iG_2(:) 

  character(256)           :: oldfilename, newfilename
  real(dp)                 :: time1, time2
  real(dp)                 :: reading_time, writing_time, computing_time

  nbands_loc=num_bands

  reading_time   = 0.0_dp
  writing_time   = 0.0_dp
  computing_time = 0.0_dp

  oldfilename = trim('W_matrix_elements.bin')
  newfilename = trim('W_matrix_elements_new.bin')

  allocate(W_rotated_matrix_elements(nG_shell_max,num_wann,num_wann,nspin,nspin))
  allocate(pw_matrix_elements(nbands_loc,nbands_loc,nspin,nspin))
!  allocate(pw_matrix_elements(nbands,nbands,nspin,nspin))
  allocate(wfc_1(nG_max,nbands_loc,nspin))
!  allocate(wfc_1(nG_max,nbands,nspin))
  allocate(wfc_2(nG_max,nbands_loc,nspin))
!  allocate(wfc_2(nG_max,nbands,nspin))
  allocate(list_iG_1(nG_max))
  allocate(list_iG_2(nG_max))
  allocate(QE_eig(nbands_loc))
!  allocate(QE_eig(nbands))
  
  ! read X_l stored in the file
  call read_write_checkpoint_xml(N_scheme,D_scheme,sze_scheme,X_l_file, 'read')

  ! open files
  record_length = num_wann*num_wann*nspin**2*direct_io_factor_cmplx

  io_unit_old_bin = find_free_unit()

  open(unit = io_unit_old_bin, file = oldfilename, form = 'unformatted', &
     status = 'unknown', access = 'direct', recl = record_length, &
     action = 'read')

  io_unit_new_bin = find_free_unit()

  open(unit = io_unit_new_bin, file = newfilename, form = 'unformatted', &
     status = 'unknown', access = 'direct', recl = record_length, &
     action = 'write')

  n_found = 0
  ! for every l_scalar, check if data is already available; if not, compute

  do l_scalar = 1,sze_scheme
        call find_data(N_scheme,D_scheme,sze_scheme, found, &
                        X_l,X_l_file,l_scalar,l_scalar_file)

        if (found ) then         
           n_found = n_found + 1 
           ! read the data and write it in the new file
           call get_timing(time1)
           do iG = 1, nG_shell_max
              ! read 
              call get_record_index(l_scalar_file,iG,record_index)

              read(io_unit_old_bin,rec=record_index)    &
                                        W_rotated_matrix_elements(iG,:,:,:,:)
           end do
           call get_timing(time2)
           reading_time = reading_time+time2-time1

        else
        ! compute the data since it is not available
           call get_timing(time1)
           call compute_matrix_elements_at_l_scalar          &
                   (N_scheme,D_scheme,sze_scheme,use_IBZ,    &
                    l_scalar,X_l,W_rotated_matrix_elements,  &
                    pw_matrix_elements,wfc_1,wfc_2,QE_eig,   &
                    list_iG_1,list_iG_2)
           call get_timing(time2)
           computing_time = computing_time+time2-time1
        end if

        ! write data to new file
        call get_timing(time1)
        do iG = 1, nG_shell_max
           ! write 
           call get_record_index(l_scalar,iG,record_index)
           write(io_unit_new_bin,rec=record_index)    &
                                  W_rotated_matrix_elements(iG,:,:,:,:)
 !         print *, 'W = ', sqrt(sum(abs(W_rotated_matrix_elements(iG,:,:))**2))
        end do
        call get_timing(time2)
        writing_time = writing_time + time2-time1
  end do

  ! finish up
 
  close(io_unit_old_bin)
  close(io_unit_new_bin)

  ! replace the old file with the new file
  call system('mv '//trim(newfilename)//' '//trim(oldfilename))

  ! update the xml file
  call read_write_checkpoint_xml(N_scheme,D_scheme,sze_scheme,X_l, 'write')

  deallocate(pw_matrix_elements)
  deallocate(W_rotated_matrix_elements)
  deallocate(wfc_1)
  deallocate(wfc_2)
  deallocate(list_iG_1)
  deallocate(list_iG_2)
  deallocate(QE_eig)


  end subroutine update_rotated_matrix_elements_file

 
  subroutine find_data(N_scheme,D_scheme,sze_scheme, found, &
                        X_l,X_l_old,l_scalar,l_scalar_old)
  !------------------------------------------------------------------
  ! 
  ! This subroutine finds the common points in the mesh arrays X_l
  ! and X_l_old
  !------------------------------------------------------------------
  integer        :: N_scheme, D_scheme, sze_scheme

  integer        :: X_l(sze_scheme,D_scheme), X_l_old(sze_scheme,D_scheme)

  integer        :: l_scalar, l_scalar_old, test

  logical        :: found

  found = .false.

  do l_scalar_old = 1, sze_scheme
       test = sum((X_l(l_scalar,:)-X_l_old(l_scalar_old,:))**2) 

       if ( test == 0 ) then
                found = .true.
                exit
       end if

  end do
  if ( .not. found ) l_scalar_old  = 0

  end subroutine find_data

  subroutine find_boundaries(nk,nks,i,is_min,is_max)
  !----------------------------------------------------------------------------!
  ! convenience subroutine
  !----------------------------------------------------------------------------!
  implicit none

  integer       :: nk, nks, i, is_min, is_max

  if ( modulo(nks*(i-1),nk) == 0) then
        ! careful here! You are dividing INTEGERS
        is_min = (nks*(i-1))/nk+1
  else
        is_min = ceiling(dble(nks)/dble(nk)*(dble(i)-1.0_dp)+1.0_dp)
  end if

  if ( modulo(nks*i,nk) == 0) then
        ! careful here! You are dividing INTEGERS
        is_max = (nks*i)/nk
  else
        is_max = floor(dble(nks)/dble(nk)*dble(i)+1.0_dp)
  end if 

  end subroutine find_boundaries

  subroutine find_ks_and_ksq_in_sub_block(X0,qpt,list_ks,list_ksq,   &
                                          list_t,nk_vec_max,nk_vec)
  !----------------------------------------------------------------------------!
  !    Given a sub-block identified by its corner X0, this subroutine 
  !    computes the vectors ks and ks+q which are inside the sub-block, where
  !    ks is in the smooth mesh. It also computes the list of t vectors, where
  !    t[i] = x[i]/Delta x[i] - [x[i]/Delta x[i]] is the relative coordinate.
  !    
  !    nk_vec is the number of such k vectors; the arrays will be larger,
  !    of size nk_vec_max, for convenience.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: X0(6)
  integer        :: nk_vec_max
  real(dp)       :: qpt(3)

  
  ! output variables
  integer        :: nk_vec
  real(dp)       :: list_ks(3,nk_vec_max), list_ksq(3,nk_vec_max)
  real(dp)       :: list_t(6,nk_vec_max)


  ! local variables
  real(dp)       :: ks(3), ksq(3)

  integer        :: i1, i2, i3, i4, i5, i6
  integer        :: i4_test, i5_test, i6_test

  integer        :: i1s, i2s, i3s
  integer        :: i1s_min, i1s_max
  integer        :: i2s_min, i2s_max
  integer        :: i3s_min, i3s_max

  logical        :: right_corner

  ! initialize 
  list_ks  = 0.0_dp
  list_ksq = 0.0_dp
  list_t   = 0.0_dp
  nk_vec   = 0

  ! Identify the corner of the cube in the 1BZ in which ks is located
  i1    = X0(1)
  i2    = X0(2)
  i3    = X0(3)

  ! Identify the corner of the cube in which ks+q  is located
  i4    = X0(4)
  i5    = X0(5)
  i6    = X0(6)

  ! identify the boundaries of the triplet indices of ks so that it
  ! is located inside the cube identified by X0(1:3)
  call find_boundaries(nk1,nk1s,i1,i1s_min,i1s_max)
  call find_boundaries(nk2,nk2s,i2,i2s_min,i2s_max)
  call find_boundaries(nk3,nk3s,i3,i3s_min,i3s_max)

  ! loop through the ks thus defined, and find if ks+q is in the cube X0(4:6)
  do i1s=i1s_min, i1s_max
     do i2s=i2s_min, i2s_max
        do i3s=i3s_min, i3s_max
           ! find ks
           ks(1)  = dble(i1s-1)/dble(nk1s)
           ks(2)  = dble(i2s-1)/dble(nk2s)
           ks(3)  = dble(i3s-1)/dble(nk3s)

           ! find ks+q
           ksq(1) = ks(1)+qpt(1)
           ksq(2) = ks(2)+qpt(2)
           ksq(3) = ks(3)+qpt(3)

           ! find the corner in the coarse mesh associated with ksq
           i4_test = floor(dble(nk1)*ksq(1)+1.0_dp)
           i5_test = floor(dble(nk2)*ksq(2)+1.0_dp)
           i6_test = floor(dble(nk3)*ksq(3)+1.0_dp)

           ! check if this is the right corner
           right_corner =       (i4_test == i4 ) &
                          .and. (i5_test == i5 ) &
                          .and. (i6_test == i6 )

           ! if it is the right corner, add these points to the list
           if (right_corner) then
                nk_vec = nk_vec+1

                list_ks (:,nk_vec)  = ks(:)
                list_ksq(:,nk_vec)  = ksq(:)

                list_t(1,nk_vec)    = dble(nk1)*ks(1)-floor(dble(nk1)*ks(1))
                list_t(2,nk_vec)    = dble(nk2)*ks(2)-floor(dble(nk2)*ks(2))
                list_t(3,nk_vec)    = dble(nk3)*ks(3)-floor(dble(nk3)*ks(3))

                list_t(4,nk_vec)    = dble(nk1)*ksq(1)-floor(dble(nk1)*ksq(1))
                list_t(5,nk_vec)    = dble(nk2)*ksq(2)-floor(dble(nk2)*ksq(2))
                list_t(6,nk_vec)    = dble(nk3)*ksq(3)-floor(dble(nk3)*ksq(3))

           end if

       end do
     end do
  end do

  end subroutine find_ks_and_ksq_in_sub_block

  subroutine create_coefficients_file(N_scheme,D_scheme,sze_scheme)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes the C coefficients necessary for the hexa
  ! interpolation scheme for a given sub-block, and stores the data
  ! in a binary file. This is because the amount of data can get quite
  ! large and cannot be kept in RAM. 
  !
  ! It is unclear at this point what is the best mkl BLAS routine
  ! to use to improve performance. 
  !------------------------------------------------------------------

  implicit none


  integer        ::      N_scheme, D_scheme,sze_scheme


  integer        ::      l_scalar, iG

  integer        ::      nbnd1, nbnd2 

  ! variables to get record index
  integer        ::      record_index
  integer        ::      record_length
  integer        ::      io_unit_bin, io_unit_coeff

  ! CAREFUL! the array where the rotated matrix elements are 
  ! stored now contains every l_scalar value. Do not be confused:
  ! when building the bin file, this array had nG_shell_max as 
  ! a first dimension: here the purpose is different.
  complex(dp)    ::      W_rotated_matrix_elements(sze_scheme,num_wann,num_wann,nspin,nspin)
  complex(dp)    ::      f(sze_scheme), c(sze_scheme)
  real   (dp)    ::      re_matf(sze_scheme,num_wann), re_matc(sze_scheme,num_wann)
  real   (dp)    ::      im_matf(sze_scheme,num_wann), im_matc(sze_scheme,num_wann)
  complex(dp)    ::      matf(sze_scheme,num_wann), matc(sze_scheme,num_wann)

  complex(dp)    ::      c_coefficients(sze_scheme,num_wann,num_wann,nspin,nspin)

  real(dp)       ::      time1, time2
  real(dp)       ::      computing_time, writing_time, reading_time

  character(*),parameter  :: filename_bin   = 'W_matrix_elements.bin'
  character(*),parameter  :: filename_coeff = 'coefficients.bin'

  character      ::  matdescra(6)
  integer        :: ipol, jpol

  record_length = num_wann*num_wann*nspin**2*direct_io_factor_cmplx

  ! open the bin file
  io_unit_bin = find_free_unit()
  open(unit = io_unit_bin, file = filename_bin, form = 'unformatted', &
     status = 'unknown', access = 'direct', recl = record_length, &
     action = 'read')

  ! open the coefficients file
  io_unit_coeff = find_free_unit()
  open(unit = io_unit_coeff, file = filename_coeff, form = 'unformatted', &
     status = 'unknown', access = 'direct', recl = record_length, &
     action = 'write')

  reading_time   = ZERO
  writing_time   = ZERO
  computing_time = ZERO

  ! loop on G vectors
  do iG=1,nG_shell_max

     ! extract a "G-slice" of the rotated matrix elements,
     ! ie the elements for every corner for a given G.
     call get_timing(time1)
     do l_scalar = 1, sze_scheme
        ! get the record index
        call get_record_index(l_scalar,iG,record_index)

        ! read the matrix elements from the bin file
        read(io_unit_bin,rec=record_index) W_rotated_matrix_elements(l_scalar,:,:,:,:)

!        print *, 'l = ',l_scalar,' |W| =', sqrt(sum(abs(W_rotated_matrix_elements(l_scalar,:,:))**2))

     end do
     call get_timing(time2)
     reading_time = reading_time+time2-time1

     ! compute the coefficients using the interpolation scheme, for every
     ! bands
!     do nbnd1 = 1,num_wann
!        do nbnd2 = 1,num_wann
       

!	TRY SPARSE BLAS!
!            f(:)  = W_rotated_matrix_elements(:,nbnd1,nbnd2)

!           call get_timing(time1)
!           call mkl_zcoogemv('n', sze_scheme,  MA_sparse_array,          &
!                              MA_row_indices, MA_col_indices,            &
!                              sparse_sze,f,c)
!
!           call zgemv('n', sze_scheme, sze_scheme, cmplx_1,     &
!                       cmplx_MA, sze_scheme, f, 1, cmplx_0, c, 1)

!           c = matmul(cmplx_MA,f)

!           call get_timing(time2)

!           c_coefficients(:,nbnd1,nbnd2) = c(:)

!           computing_time = computing_time+time2-time1
!        end do
!     end do

     matdescra(1) = 'G'
     matdescra(4) = 'F'
     do nbnd2 = 1,num_wann

      do ipol=1,nspin
       do jpol=1,nspin

          matf = W_rotated_matrix_elements(:,:,nbnd2,ipol,jpol)

          re_matf = real(matf)
          im_matf = aimag(matf)

          call get_timing(time1)

          call mkl_dcoomm('n',sze_scheme, num_wann, sze_scheme,        &
                     ONE, matdescra,                                   &
                     MA_sparse_array, MA_row_indices, MA_col_indices,  &
                     sparse_sze,                                       &
                     re_matf, sze_scheme,ZERO,re_matc, sze_scheme)

          call mkl_dcoomm('n',sze_scheme, num_wann, sze_scheme,        &
                     ONE, matdescra,                                   &
                     MA_sparse_array, MA_row_indices, MA_col_indices,  &
                     sparse_sze,                                       &
                     im_matf, sze_scheme,ZERO,im_matc, sze_scheme)

!           matc = cmplx_1*re_matc+cmplx_i*im_matc

           call get_timing(time2)
           computing_time = computing_time+time2-time1
           
!          c_coefficients(:,:,nbnd2) = matc
           c_coefficients(:,:,nbnd2,ipol,jpol) = cmplx_1*re_matc+cmplx_i*im_matc
         end do !jpol
       end do !ipol
     end do


     ! write the coefficients to file
     call get_timing(time1)
     do l_scalar = 1, sze_scheme
        call get_record_index(l_scalar,iG,record_index)
        write(io_unit_coeff,rec=record_index) c_coefficients(l_scalar,:,:,:,:)
!        print *, 'l = ',l_scalar,' |c| =', sqrt(sum(abs(c_coefficients(l_scalar,:,:))**2))
     end do
     call get_timing(time2)
     writing_time = writing_time+time2-time1

  end do


  close(io_unit_bin)
  close(io_unit_coeff)

  end subroutine create_coefficients_file

  subroutine fetch_c_coefficients(sze_scheme,iG,c_coefficients)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes the C coefficients necessary for the hexa
  ! interpolation scheme for a given sub-block, and stores the data
  ! in a binary file. This is because the amount of data can get quite
  ! large and cannot be kept in RAM. 
  !
  !------------------------------------------------------------------

  implicit none


  integer        ::      sze_scheme


  integer        ::      l_scalar, iG

  ! variables to get record index
  integer        ::      record_index
  integer        ::      record_length
  integer        ::      io_unit_coeff

  complex(dp)    ::      c_coefficients(sze_scheme,num_wann,num_wann,nspin,nspin)

  character(*),parameter  :: filename_coeff = 'coefficients.bin'

  record_length = num_wann*num_wann*nspin**2*direct_io_factor_cmplx

  ! open the coefficients file
  io_unit_coeff = find_free_unit()

  open(unit = io_unit_coeff, file = filename_coeff, form = 'unformatted', &
     status = 'unknown', access = 'direct', recl = record_length, &
     action = 'read')

  ! extract a "G-slice" of the rotated matrix elements,
  ! ie the elements for every corner for a given G.
  do l_scalar = 1, sze_scheme
     ! get the record index
     call get_record_index(l_scalar,iG,record_index)

     ! read the matrix elements from the bin file
     read(io_unit_coeff,rec=record_index) c_coefficients(l_scalar,:,:,:,:)
  end do

  close(io_unit_coeff)

  end subroutine fetch_c_coefficients

  subroutine computed_interpolated_matrix_elements(nk_vec_max,nk_vec, &
                sze_scheme,c_coefficients,list_mu,u_int_ks,u_int_ksq, &
                pw_matrix_elements)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes interpolated matrix elements for all
  ! ks belonging to a given sub-block, for a given G. 
  !
  !------------------------------------------------------------------
  implicit none

  ! input variables
  integer       ::  sze_scheme
  integer       ::  nk_vec_max, nk_vec
  complex(dp)   ::  c_coefficients(sze_scheme,num_wann,num_wann,nspin,nspin)
  complex(dp)   ::  u_int_ks (num_wann,num_wann,nk_vec_max)
  complex(dp)   ::  u_int_ksq(num_wann,num_wann,nk_vec_max)
  real(dp)      ::  list_mu(sze_scheme,nk_vec_max)

  ! output variable
  complex(dp)   ::  pw_matrix_elements(num_wann,num_wann,nspin,nspin,nk_vec_max)


  ! local variables
  complex(dp)   ::  pw_work(num_wann,num_wann,nspin,nspin)

  integer       ::  nbnd1, nbnd2, nbnd3, nbnd4 
  integer       ::  nk_loop, l_scalar, ipol, jpol,  ipolp, jpolp 

  ! initialize
  pw_matrix_elements  =  cmplx_0

  do nk_loop = 1,nk_vec
     ! interpolate the ROTATED matrix elements for every band,
     ! for very ks in the sub-block
     do nbnd1 =1, num_wann
        do nbnd2 =1, num_wann
          do ipol=1,nspin
           do jpol=1,nspin
            do l_scalar = 1, sze_scheme
              pw_matrix_elements(nbnd1,nbnd2,ipol,jpol,nk_loop)  =    &         
              pw_matrix_elements(nbnd1,nbnd2,ipol,jpol,nk_loop)  +    &
                c_coefficients(l_scalar,nbnd1,nbnd2,ipol,jpol)*list_mu(l_scalar,nk_loop)
            end do
           end do
          end do

        end do
     end do
     ! ANTI ROTATE
     pw_work  =  cmplx_0

     do nbnd1 =1, num_wann
        do nbnd2 =1, num_wann
           do nbnd3 =1, num_wann
              do nbnd4 =1, num_wann
                 
                do ipol =1, nspin 
                  do jpol =1, nspin 

                         pw_work(nbnd1,nbnd2,ipol,jpol) = pw_work(nbnd1,nbnd2,ipol,jpol) +      &
                                       u_int_ks(nbnd1,nbnd3,nk_loop)*                           &
                                       pw_matrix_elements(nbnd3,nbnd4,ipol,jpol,nk_loop)*       &
                                     conjg(u_int_ksq(nbnd2,nbnd4,nk_loop))

                   enddo
                  enddo

              end do
           end do
        end do
     end do
    ! store back the matrix elements in the original array
    pw_matrix_elements(:,:,:,:,nk_loop)  =   pw_work(:,:,:,:)
  end do  ! nk_loop

  end subroutine computed_interpolated_matrix_elements

  subroutine computed_interpolated_matrix_elements_no_back_rotate     &
                (nk_vec_max,nk_vec,sze_scheme,c_coefficients,         &
                 list_mu,pw_matrix_elements)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes interpolated matrix elements for all
  ! ks belonging to a given sub-block, for a given G, but 
  ! doesn't anti W_rotate the result.
  !------------------------------------------------------------------
  implicit none

  ! input variables
  integer       ::  sze_scheme
  integer       ::  nk_vec_max, nk_vec
  complex(dp)   ::  c_coefficients(sze_scheme,num_wann,num_wann,nspin,nspin)
  real(dp)      ::  list_mu(sze_scheme,nk_vec_max)

  ! output variable
  complex(dp)   ::  pw_matrix_elements(num_wann,num_wann,nspin,nspin,nk_vec_max)


  ! local variables
  integer       ::  nk_loop, nbnd1, nbnd2,ipol,jpol

  ! initialize
  pw_matrix_elements  =  cmplx_0

  do nk_loop = 1,nk_vec
     ! interpolate the ROTATED matrix elements for every band,
     ! for very ks in the sub-block
     do nbnd1 =1, num_wann
        do nbnd2 =1, num_wann
          do ipol=1,nspin
           do jpol=1,nspin
              pw_matrix_elements(nbnd1,nbnd2,ipol,jpol,nk_loop)  =    &         
                sum(c_coefficients(:,nbnd1,nbnd2,ipol,jpol)*list_mu(:,nk_loop))
           enddo
          enddo

        end do
     end do
  end do

  end subroutine computed_interpolated_matrix_elements_no_back_rotate 

!--------------------------------------------------------------------------------
!
end module intw_plane_wave_setup
!
!--------------------------------------------------------------------------------

