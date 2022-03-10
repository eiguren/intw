!llwfc----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
!
module intw_intw2wannier
!
!----------------------------------------------------------------------------!
!       This module contains subroutines which will perform the same tasks
!       as the Quantum Espresso utility program "pw2wannier", but utilizing
!       symmetries.
!
!       In particular, this module will produce the file $seed.mmn
!       using only the wavefunctions in the IBZ, rotating them appropriately.
!----------------------------------------------------------------------------!
  use kinds, only: dp
  use intw_reading, only: nbands, nG_max, ngm, nspin
  use intw_useful_constants, only: bohr, pi, tpi, fpi, eps_8, ZERO, cmplx_0, cmplx_i

  implicit none
  !
  ! variables
  public :: nnkp_exclude_bands, nnkp_real_lattice, nnkp_recip_lattice, &
            nnkp_num_kpoints, nnkp_nnkpts, nnkp_n_proj, nnkp_kpoints, &
            nnkp_Wcenters, nnkp_proj_x, nnkp_proj_z, nnkp_proj_zona, &
            nnkp_proj_n, nnkp_proj_l, nnkp_proj_m, nnkp_list_ikpt_nn, &
            nnkp_list_G, nnkp_excluded_bands
  !
  ! subroutines
  public :: deallocate_nnkp, scan_file_to, read_nnkp_file, output_nnkp_file, &
            intw2W90_check_mesh, generate_header,  &
            generate_mmn_using_allwfc, generate_amn_using_allwfc, &
            generate_guiding_function, &
            get_guiding_function_overlap_FFT, &
            get_guiding_function_overlap_convolution, get_radial_part, &
            get_angular_part, ylm_wannier
  !
  private
  !
  save
  ! declare some global variables; the prefix "nnkp" will indicate
  ! variables read from the $prefix.nnkp file, which must be checked with
  ! the data read from the input file for consistency.

  external :: s, pz_func, px, py, dz2, dxz, dyz, dx2my2, dxy, fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2

  integer     :: nnkp_exclude_bands     !how many bands are excluded?

  real(dp)    :: nnkp_real_lattice(3,3) ! the lattice vectors
                                        ! These are in angstrom coordinates in the file
                                        ! but will be transformed to alat cartesian
                                        ! coordinates

  real(dp)    :: nnkp_recip_lattice(3,3)! the reciprocal lattice vectors
                                        ! These are in angstrom^-1 coordinates in the file
                                        ! but will be transformed to 2pi/alat cartesian
                                        ! coordinates

  integer     :: nnkp_num_kpoints       !the number of k-points in the 1BZ

  integer     :: nnkp_nnkpts            !the number near neighbor k-points

  integer     :: nnkp_n_proj            !the number of Wannier functions used.

  real(dp),allocatable :: nnkp_kpoints(:,:)   ! the k vectors in the 1BZ
  real(dp),allocatable :: nnkp_Wcenters(:,:)  ! the Wannier centers, for trial projection
  real(dp),allocatable :: nnkp_proj_x(:,:)    ! x projection, for trial
  real(dp),allocatable :: nnkp_proj_z(:,:)    ! z projection, for trial
  real(dp),allocatable :: nnkp_proj_zona(:)   ! Z on a: gauge broadness of the hydrogenoid trial wfc.

  ! The hydrogenoid quantum numbers, psi_{nlm}
  integer,allocatable  :: nnkp_proj_n(:)     ! the radial-pojection, for trial
  integer,allocatable  :: nnkp_proj_l(:)     ! the l-pojection, for trial
  integer,allocatable  :: nnkp_proj_m(:)     ! the m-pojection, for trial

  integer,allocatable  :: nnkp_list_ikpt_nn(:,:) ! the neighbors of ikpt1
  integer,allocatable  :: nnkp_list_G(:,:,:)   ! the G vectors linking one point to another

  logical,allocatable  :: nnkp_excluded_bands(:) ! what bands are excluded

contains
  !
!---------------------------------
  subroutine deallocate_nnkp()
!---------------------------------

  implicit none

  deallocate(nnkp_kpoints)
  deallocate(nnkp_Wcenters)
  deallocate(nnkp_proj_l)
  deallocate(nnkp_proj_m)
  deallocate(nnkp_proj_n)
  deallocate(nnkp_proj_x)
  deallocate(nnkp_proj_z)
  deallocate(nnkp_proj_zona)
  deallocate(nnkp_list_ikpt_nn)
  deallocate(nnkp_list_G)
  deallocate(nnkp_excluded_bands)

  end subroutine deallocate_nnkp
!----------------------------------------------
!**********************************************
!----------------------------------------------
  subroutine scan_file_to (nnkp_unit,keyword)
!----------------------------------------------
!
!----------------------------------------------------------------------!
! This subroutine reads a file all the way to the line
!            begin $keyword
! This is useful when extracting parameters from the ascii file $seed.nnkp.
! The subroutine is heavily inspired by the subroutine
!     QE/PP/pw2pwannier.f90{scan_file_to}
!-----------------------------------------------------------------------!

  implicit none

  character(len=*)   :: keyword
  character(len=256) :: word1, word2
  logical            :: found, test
  integer            :: ios,   nnkp_unit

  found=.false.
  !
  do
     !
     read(nnkp_unit,*,iostat=ios) word1, word2
     !
     test=(trim(word1).eq.'begin').and.(trim(word2).eq.keyword)
     !
     if (test) exit
     !
     if (ios.ne.0) then
        !
        write (*,*) keyword," data-block missing "
        stop
        !
     endif
     !
  enddo
  !
  return

  end subroutine scan_file_to
!------------------------------------------
!******************************************
!------------------------------------------
  subroutine read_nnkp_file(nnkp_file)
!------------------------------------------
!
!----------------------------------------------------------------------------!
! This subroutine reads the .nnkp input file produced by W90.
! It assumes the file has been checked for existence and is already
! opened.
!----------------------------------------------------------------------------!

  use intw_reading, only: alat
  use intw_utility, only: find_free_unit
  use w90_parameters, only: num_exclude_bands, exclude_bands

  implicit none

  integer :: nnkp_unit
  integer :: i, j, nn, dummy
  character(*) :: nnkp_file

  ! read in the various quantities stored in the .nnkp parameters file.

  nnkp_unit=find_free_unit()
  !
  !===============================
  ! open the file, which exists!
  !===============================
  !
  open(unit=nnkp_unit,file=nnkp_file,status='old')
  !
  !=======================
  ! real lattice vectors
  !=======================
  !
  call scan_file_to (nnkp_unit,'real_lattice')
  !
  do j=1,3
     !
     read(nnkp_unit,*) (nnkp_real_lattice(i,j),i=1,3)
     !
  enddo
  !
  ! convert to alat coordinates
  !
  nnkp_real_lattice=nnkp_real_lattice/(alat*bohr)
  !
  !==============================
  ! reciprocal lattice vectors
  !==============================
  !
  call scan_file_to (nnkp_unit,'recip_lattice')
  !
  do j=1,3
     !
     read(nnkp_unit,*) (nnkp_recip_lattice(i,j),i=1,3)
     !
  enddo
  !
  ! convert to 2pi/alat coordinates
  !
  nnkp_recip_lattice=nnkp_recip_lattice*(alat*bohr)/tpi
  !
  !======================================
  ! kpoints in the zone (crystal units)
  !======================================
  !
  call scan_file_to (nnkp_unit,'kpoints')
  read(nnkp_unit,*) nnkp_num_kpoints
  !
  allocate(nnkp_kpoints(3,nnkp_num_kpoints))
  !
  do j=1,nnkp_num_kpoints
     !
     read(nnkp_unit,*) (nnkp_kpoints(i,j),i=1,3)
     !
  enddo
  !
  !==========================
  ! projection information
  !==========================
  !
  call scan_file_to (nnkp_unit,'projections')
  !
  read(nnkp_unit,*) nnkp_n_proj
  !
  allocate(nnkp_Wcenters(3,nnkp_n_proj))
  allocate(nnkp_proj_l(nnkp_n_proj))
  allocate(nnkp_proj_m(nnkp_n_proj))
  allocate(nnkp_proj_n(nnkp_n_proj))
  allocate(nnkp_proj_x(3,nnkp_n_proj))
  allocate(nnkp_proj_z(3,nnkp_n_proj))
  allocate(nnkp_proj_zona(nnkp_n_proj))
  !
  do j=1,nnkp_n_proj
     !
     read(nnkp_unit,*) (nnkp_Wcenters(i,j),i=1,3), &
                        nnkp_proj_l(j), &
                        nnkp_proj_m(j), &
                        nnkp_proj_n(j)
     !
     read(nnkp_unit,*) (nnkp_proj_z(i,j),i=1,3), &
                       (nnkp_proj_x(i,j),i=1,3), &
                        nnkp_proj_zona(j)
     !
  enddo
  !
  !==========================================
  ! connection information (near neighbors)
  !==========================================
  !
  call scan_file_to (nnkp_unit,'nnkpts')
  !
  read(nnkp_unit,*) nnkp_nnkpts
  !
  allocate(nnkp_list_ikpt_nn(nnkp_nnkpts,nnkp_num_kpoints))
  allocate(nnkp_list_G(3,nnkp_nnkpts,nnkp_num_kpoints))
  !
  do j=1,nnkp_num_kpoints
     do nn=1,nnkp_nnkpts
        !
        read(nnkp_unit,*) dummy,nnkp_list_ikpt_nn(nn,j),(nnkp_list_G(i,nn,j),i=1,3)
        !
     enddo
  enddo
  !
  !=================
  ! band exclusion
  !=================
  !
  call scan_file_to (nnkp_unit,'exclude_bands ')
  !
  read(nnkp_unit,*) nnkp_exclude_bands
  !
  num_exclude_bands=nnkp_exclude_bands
  !
  allocate(exclude_bands(num_exclude_bands))
  allocate(nnkp_excluded_bands(nbands))
  !
  nnkp_excluded_bands(:)=.false.
  !
  do j=1,nnkp_exclude_bands
     !
     read(nnkp_unit,*) nn
     nnkp_excluded_bands(nn)=.true.
     exclude_bands(j)=nn
     !
  enddo
  !
  !=============
  ! close file
  !=============
  !
  close(unit=nnkp_unit)
  !
  return

  end subroutine read_nnkp_file
!--------------------------------------
!**************************************
!--------------------------------------
  subroutine output_nnkp_file()
!--------------------------------------
!
!----------------------------------------------------------------------------!
! This subroutine simply outputs what was read from the nnkp file, for
! testing.
!----------------------------------------------------------------------------!
  use intw_utility, only: find_free_unit

  implicit none

  integer  :: io_unit
  integer  :: i, j, nn

  io_unit = find_free_unit()

  open(unit=io_unit,file='nnkp.test',status='unknown')

  write(io_unit,*) '====================================================='
  write(io_unit,*) '|            content of prefix.nnkp                 |'
  write(io_unit,*) '|        ---------------------------------          |'
  write(io_unit,*) '|        This file is generated for testing         |'
  write(io_unit,*) '|        purposes; it should be compared to         |'
  write(io_unit,*) '|        the actual .nnkp file to check             |'
  write(io_unit,*) '|        consistency.                               |'
  write(io_unit,*) '====================================================='


  write(io_unit,*) 'nnkp_real_lattice'
  write(io_unit,*) '-----------------'
  do j=1,3
        write(io_unit,'(3F8.4)') (nnkp_real_lattice(i,j),i=1,3)
  end do
  write(io_unit,*) ''

  write(io_unit,*) 'nnkp_recip_lattice'
  write(io_unit,*) '------------------'
  do j=1,3
        write(io_unit,'(3F8.4)') (nnkp_recip_lattice(i,j),i=1,3)
  end do
  write(io_unit,*) ''

  write(io_unit,'(a,I6)') 'nnkp_num_kpoints = ',nnkp_num_kpoints
  write(io_unit,*) ''

  write(io_unit,'(a,I6)') 'k-points'
  write(io_unit,'(a,I6)') '--------'
  do j=1,nnkp_num_kpoints
        write(io_unit,'(3F8.4)') (nnkp_kpoints(i,j),i=1,3)
  end do


  write(io_unit,'(a,I6)') 'nnkp_n_proj= ',nnkp_n_proj
  write(io_unit,*) ''
  write(io_unit,*) 'projection information'
  write(io_unit,*) '----------------------'
  write(io_unit,*) '  n          center           l   m   n          proj_z                  proj_x            Z/a'
  write(io_unit,*) '---------------------------------------------------------------------------------------------------'
  do j=1,nnkp_n_proj

     write(io_unit,100) j,(nnkp_Wcenters(i,j),i=1,3),   &
                           nnkp_proj_l(j),              &
                           nnkp_proj_m(j),              &
                           nnkp_proj_n(j),              &
                          (nnkp_proj_z  (i,j),i=1,3),   &
                          (nnkp_proj_x  (i,j),i=1,3),   &
                           nnkp_proj_zona( j)

  end do

  write(io_unit,*) ''
  write(io_unit,*) 'connection information'
  write(io_unit,*) '----------------------'
  write(io_unit,*) 'ikpt1     ikpt2       G'
  write(io_unit,*) '----------------------------'

  do j=1,nnkp_num_kpoints
     do nn = 1,nnkp_nnkpts

       write(io_unit,200) j,nnkp_list_ikpt_nn(nn,j),(nnkp_list_G(i,nn,j),i=1,3)
     end do
  end do

  close(unit=io_unit)

  100 format(I4,3F8.4,3I4,7F8.4)
  200 format(I4,6x,I4,2x,3I4)

  end subroutine output_nnkp_file
!-----------------------------------------------
!***********************************************
!-----------------------------------------------
  subroutine intw2W90_check_mesh(nkmesh,kmesh)
!-----------------------------------------------
!
!--------------------------------------------------------------------!
! This subroutine checks that the irreducible kpoints and the mesh
! read from the nnkp file are related.
!--------------------------------------------------------------------!

  implicit none

  !I/O variables

  integer,intent(in) :: nkmesh
  real(dp),intent(in) :: kmesh(3,nkmesh)

  !local variables

  integer   :: ikpt
  real(dp)  :: kpt(3),norm

  if (nkmesh.ne.nnkp_num_kpoints) then
     !
     write(*,*)   '**********************************************************'
     write(*,*)   '                INCONSISTENCY ERROR                       '
     write(*,*)   '       The number of points in the MP mesh corresponding  '
     write(*,*)   '       to the input DOES NOT MATCH the number of k-points '
     write(*,*)   '       in the Wannier90 .nnkp file.                       '
     write(*,*)   '  nkmesh                                                  '
     write(*,*)   '         ',nkmesh
     write(*,*)   '  nnkp_num_kpoints                                                    '
     write(*,*)   '         ',nnkp_num_kpoints
     write(*,*)   '**********************************************************'
     !
     stop
     !
  endif
  !
  do ikpt=1,nkmesh
     !
     kpt=kmesh(:,ikpt)-nnkp_kpoints(:,ikpt)
     !
     norm=sqrt(kpt(1)**2 + kpt(2)**2 + kpt(3)**2)
     !
     if (norm>eps_8) then
        !
        write(*,*)    '**********************************************************'
        write(*,*)    '                INCONSISTENCY ERROR                   '
        write(*,*)    '  The MP mesh generated from input does not match     '
        write(*,*)    '  the MP mesh expected by Wannier90. Are the points   '
        write(*,*)    '  properly ordered in the .win file?                  '
        write(*,*)    '**********************************************************'
        !
     endif
     !
  enddo !ikpt
  !
  return

  end subroutine intw2W90_check_mesh
!--------------------------------------------
!********************************************
!--------------------------------------------
  subroutine generate_header(method,header)
!--------------------------------------------
!
!-----------------------------------------------------------
! This is a utility routine which creates a date and time
! header to provide a time stamp in our files
!-----------------------------------------------------------

  character(256) :: header
  character(8)   :: cdate
  character(10)  :: ctime
  character(*)   :: method

  call date_and_time(DATE=cdate,TIME=ctime)
  !
  header = 'intw2W90::'//trim(method)//' time:: '         &
        //trim(ctime(1:2))//':'//trim(ctime(3:4))//':'//trim(ctime(5:6)) &
        //' day:: '//trim(cdate(7:8))//'/'//trim(cdate(5:6))//'/'//trim(cdate(1:4))
  !
  return

  end subroutine generate_header
!----------------------------------------------------
!------------------------------------------------------------------
!******************************************************************
!------------------------------------------------------------------
  subroutine generate_mmn_using_allwfc (intw2W_fullzone,method)
  !----------------------------------------------------------------------------!
  ! This subroutine computes the plane wave matrix elements  needed by Wannier90
  ! by using symmetry. It fetches the wfc in the IBZ, rotates them, and computes
  ! the needed matrix elements.
  !
  !----------------------------------------------------------------------------!
  USE intw_allwfcs, only: get_psi_general_k_all_wfc, ngk_all
  use intw_utility, only: find_free_unit
  use intw_matrix_elements, only: get_plane_wave_matrix_element_FFT, get_plane_wave_matrix_element_convolution
  use intw_matrix_elements, only: get_plane_wave_matrix_element_convolution_map
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3
  use w90_parameters, only: num_bands
  use intw_symmetries, only:  QE_folder_sym
  

  implicit none

  logical        :: intw2W_fullzone
  character(*)   :: method


  integer        :: io_unit_mmn, io_unit_eig
  integer        :: nkmesh

  integer        :: nn

  character(256) :: filename

  integer        :: ikpt_1, ikpt_2
  integer        :: nb1, nb2, nb
  integer        :: G(3), G_plus(3)

  integer        :: ngk1, ngk2

  integer        :: list_iG_1(nG_max), list_iG_2(nG_max)
  complex(dp)    :: wfc_1(nG_max,num_bands,nspin), wfc_2(nG_max,num_bands,nspin)
  real(dp)       :: QE_eig(num_bands)

  complex(dp)    :: pw_mat_el(num_bands,num_bands,nspin,nspin)

  character(256) :: header

  integer        :: nbands_loc

  nbands_loc=num_bands
  !-----------------------------------
  ! Open all the needed files
  !-----------------------------------
  io_unit_mmn = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('.mmn')
  open(unit=io_unit_mmn,file=filename,status='unknown')

  io_unit_eig = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('.eig')
  open(unit=io_unit_eig,file=filename,status='unknown')

  !-----------------------------------
  ! define a few useful variables
  !-----------------------------------
  nkmesh = nk1*nk2*nk3

  if (intw2W_fullzone) then
          call generate_header(trim(method)//trim('-fullzone'),header)
  else
          call generate_header(trim(method)//trim('-IBZ'),header)
  end if

  write(io_unit_mmn,*) trim(header)
  write(io_unit_mmn,'(3i12)') nbands-nnkp_exclude_bands, nnkp_num_kpoints , nnkp_nnkpts

  !loop on all points
  do ikpt_1 = 1, nkmesh
    ! fetch the data
    call get_psi_general_k_all_wfc( .true., nnkp_kpoints(:,ikpt_1),list_iG_1, wfc_1, QE_eig, G_plus)

    ngk1= ngk_all(QE_folder_sym(ikpt_1))

    ! print out the eigenvalues
    do nb =1, nbands_loc
       write(io_unit_eig,'(2I5,F18.12)') nb, ikpt_1, QE_eig(nb)
    end do

    ! loop on neighbors
    do nn = 1, nnkp_nnkpts
      G      = nnkp_list_G(:,nn,ikpt_1)
      ikpt_2 = nnkp_list_ikpt_nn(nn,ikpt_1)

      write(io_unit_mmn,'(5I7)')  ikpt_1,ikpt_2,G

      ! fetch data
      call get_psi_general_k_all_wfc(.true.,  nnkp_kpoints(:, ikpt_2) + G  ,list_iG_2, wfc_2, QE_eig, G_plus)

      ngk2= ngk_all(QE_folder_sym(ikpt_2))
      ! Compute the matrix elements
      if ( trim(method) == 'CONVOLUTION' ) then
          call get_plane_wave_matrix_element_convolution_map      &
                        (G,list_iG_1,ngk1,list_iG_2,ngk2, wfc_1,wfc_2,pw_mat_el)
          !call get_plane_wave_matrix_element_convolution      &
          !(G,list_iG_1,list_iG_2,wfc_1,wfc_2,pw_mat_el)


      else if ( trim(method) == 'FFT' ) then
          call get_plane_wave_matrix_element_FFT              &
                        (G,list_iG_1,list_iG_2, wfc_1,wfc_2,pw_mat_el)
      else
          write(*,*) 'ERROR in generate_mmn'
          stop
      end if

      if (nspin==2) then
         do nb1 = 1,nbands_loc
!         do nb1 = 1,nbands
           do nb2 = 1,nbands_loc
!           do nb2 = 1,nbands
              write(io_unit_mmn,'(2F18.12)')   pw_mat_el(nb2,nb1,1,1) + pw_mat_el(nb2,nb1,2,2)
           end do
         end do
      else if  (nspin==1) then
         do nb1 = 1,nbands_loc
!         do nb1 = 1,nbands
           do nb2 = 1,nbands_loc
!           do nb2 = 1,nbands
              write(io_unit_mmn,'(2F18.12)')   pw_mat_el(nb2,nb1,1,1)
           end do
         end do
      endif

    end do

  end do

  close(io_unit_mmn)
  close(io_unit_eig)

  end subroutine generate_mmn_using_allwfc

  subroutine generate_amn_using_allwfc (intw2W_fullzone,method)
  !----------------------------------------------------------------------------!
  ! This subroutine computes the overlap with trial functions, thus producing
  ! the amn file needed by Wannier90.
  !----------------------------------------------------------------------------!
   USE intw_allwfcs, only: get_psi_general_k_all_wfc
   use intw_utility, only: find_free_unit
   use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3
   use w90_parameters, only: num_bands

  implicit none

  logical        :: intw2W_fullzone
  character(256) :: method

  integer        :: io_unit_amn
  integer        :: nkmesh


  character(256) :: filename

  integer        :: ikpt
  integer        :: nb, n_proj

  integer        :: list_iG(nG_max)

  complex(dp)    :: wfc(nG_max,num_bands,nspin)

  complex(dp)    :: guiding_function(ngm)

  real(dp)       :: QE_eig(num_bands)

  complex(dp)    :: amn(num_bands)

  character(256) :: header

  integer :: G_plus(3)

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio

  nbands_loc=num_bands
  nkmesh = nk1*nk2*nk3

  io_unit_amn = find_free_unit()
  filename    = trim(mesh_dir)//trim(prefix)//trim('.amn')
  open(unit=io_unit_amn,file=filename,status='unknown')

  call generate_header(method,header)
  write(io_unit_amn,*) trim(header)


  write(io_unit_amn,'(3I12)') nbands-nnkp_exclude_bands, nnkp_num_kpoints , nnkp_n_proj


  !loop on all k-points
  do ikpt = 1, nnkp_num_kpoints

    ! fetch the data
    !call get_psi_k(ikpt,.not.intw2W_fullzone,list_iG,wfc,QE_eig)
!     call get_psi_general_k_all_wfc(.true.,  nnkp_kpoints(:, ikpt)   , nspin,list_iG, wfc, QE_eig, G_plus)
    call get_psi_general_k_all_wfc(.true.,  nnkp_kpoints(:, ikpt)   , list_iG, wfc, QE_eig, G_plus)

    !loop on all bands and all trial functions
    do n_proj = 1,nnkp_n_proj

    ! Generate the fourier transform of the trial function (called guiding
    ! function, just like in pw2wannier). In order
    ! to insure that it is normalized, the coefficients must be computed
    ! on the whole gvec mesh, not just the G vectors where the wavefunction
    ! does not vanish. It will be assumed that gvec is large enough to
    ! insure proper normalization.


      call generate_guiding_function(ikpt,n_proj,guiding_function)

      if (trim(method) == 'CONVOLUTION') then
              call get_guiding_function_overlap_convolution(list_iG,wfc,   &
                                                        guiding_function,amn)
      else if (trim(method) == 'FFT') then
              call get_guiding_function_overlap_FFT(list_iG,wfc,   &
                                                        guiding_function,amn)
      else
          write(*,*) 'ERROR in generate_amn'
          stop
      end if

      !Write result to file $prefix.amn.

      do nb = 1,nbands_loc
!      do nb = 1,nbands
        write(io_unit_amn,'(3I7,2F18.12)') nb, n_proj, ikpt, amn(nb)
      end do !nb

    end do ! n_proj


  end do !ikpt

  close(io_unit_amn)


  end subroutine generate_amn_using_allwfc
!--------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
  subroutine generate_guiding_function(ikpt,n_proj,guiding_function)
!----------------------------------------------------------------------
!
!---------------------------------------------------------------------------
! This subroutine computes the normalized guiding function in reciprocal
! space, over the global list of G vectors gvec.
!
! This subroutine is heavily inspired by a similar routine in pw2wannier.
!---------------------------------------------------------------------------

  use intw_reading, only: gvec, ngm, alat

  implicit none

  !I/O variables

  integer,intent(in) :: ikpt, n_proj
  complex(dp),intent(out) :: guiding_function(ngm)

  !local variables

  integer :: i, mu, iG
  integer :: proj_nr, proj_l, proj_m
  real(dp) :: zona, zaxis(3), xaxis(3), ylm(ngm)
  real(dp) :: k_cryst(3), tau_cryst(3), tau_cart(3)
  real(dp) :: k_plus_G_cart(3,ngm)
  real(dp) :: norm2, x, y
  complex(dp) :: four_pi_i_l

  ! get all the relevant vectors in reciprocal space
  ! express the vectors in bohrs^-1 cartesian coordinates
  !
  k_cryst(:)=nnkp_kpoints(:,ikpt)
  k_plus_G_cart=ZERO
  !
  do iG=1,ngm
     do mu=1,3
        do i=1,3
           !
           k_plus_G_cart(mu,iG)=k_plus_G_cart(mu,iG) &
                     +nnkp_recip_lattice(mu,i)*(k_cryst(i)+dble(gvec(i,iG)))
           !
        enddo !i
     enddo !mu
  enddo !iG
  !
  k_plus_G_cart=k_plus_G_cart*tpi/alat ! bohr^-1
  !
  ! TEST
  !  write(111,'(a,3F8.4)') 'k_crsyt = ',k_cryst
  !  do iG = 1, ngm
  !          write(111,'(a,3I4,a,3F8.4)') 'G = ',gvec(:,iG),' {k+G}_cart =',k_plus_G_cart(:,iG)
  !  end do
  !  stop
  !
  ! compute the guiding function
  !
  guiding_function=cmplx_0
  !
  ! get the part from the radial integration
  !
  proj_nr=nnkp_proj_n(n_proj) ! the radial projection parameter
  zona=nnkp_proj_zona(n_proj) ! Z/a, the diffusive parameter
  !
  call get_radial_part(proj_nr,zona,k_plus_G_cart,guiding_function)
  !
  ! multiply by 4 pi (-i)^l e^{-i q * Wcenter}
  !
  proj_l=nnkp_proj_l(n_proj)
  !
  tau_cryst=nnkp_Wcenters(:,n_proj) ! crystal coordinates
  tau_cart=ZERO
  !
  do mu=1,3
     do i=1,3
        !
        tau_cart(mu)=tau_cart(mu)+alat*nnkp_real_lattice(mu,i)*tau_cryst(i)
        !
     enddo !i
  enddo !mu
  !
  ! TEST
  !  write(111,'(a,3F8.4)') 'Wcenter = ',nnkp_Wcenters(:,n_proj)
  !  write(111,'(a,3F8.4)') 'tau_cart = ',tau_cart(:)
  !
  ! compute the guiding function
  !
  four_pi_i_l=fpi*(-cmplx_i)**proj_l
  !
  do iG=1,ngm
     guiding_function(iG)=guiding_function(iG)*four_pi_i_l*exp(-cmplx_i* &
                                        (k_plus_G_cart(1,iG)*tau_cart(1) &
                                       + k_plus_G_cart(2,iG)*tau_cart(2) &
                                       + k_plus_G_cart(3,iG)*tau_cart(3)))
  end do !iG
  !
  ! get angular part
  !
  proj_m=nnkp_proj_m(n_proj)
  zaxis(:)=nnkp_proj_z(:,n_proj)
  xaxis(:)=nnkp_proj_x(:,n_proj)
  !
  call get_angular_part(proj_l,proj_m,xaxis,zaxis,k_plus_G_cart,ylm)
  !
  do iG=1,ngm
     !
     guiding_function(iG)=guiding_function(iG)*ylm(iG)
     !
  enddo !iG
  !
  ! normalize
  !
  norm2=ZERO
  !
  do iG=1,ngm
     !
     x=real(guiding_function(iG))
     y=aimag(guiding_function(iG))
     !
     norm2=norm2+x**2+y**2
     !
  enddo !iG
  !
  guiding_function=guiding_function/sqrt(norm2)
  !
  return

  end subroutine generate_guiding_function
!---------------------------------------------------------------------------------
!*********************************************************************************
!---------------------------------------------------------------------------------
  subroutine get_guiding_function_overlap_FFT(list_iG,wfc,guiding_function,amn)
!---------------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!  This subroutine computes the overlap between a given wavefunction
!  wfc and a given guiding_function (assumed normalized)
!
!                    amn(band) =  < wfc(band) |  guiding_function > .
!
!  The computation is done over all bands using FFT.
!------------------------------------------------------------------------

  use intw_reading, only: ngm, nG_max, nr1, nr2, nr3
  use intw_fft, only: nl, find_iG, func_from_g_to_r, func_from_r_to_g, &
                      wfc_from_g_to_r
  use w90_parameters, only: num_bands

  implicit none

  !I/O variables

  integer,intent(in) :: list_iG(nG_max)
  complex(dp),intent(in) :: wfc(nG_max,num_bands,nspin)
  complex(dp),intent(in) :: guiding_function(ngm)
  complex(dp),intent(out) :: amn(num_bands)

  !local variables

  integer :: ibnd, ir, is, iG0, iG0_fft, G0(3)
  complex(dp) :: wfc_r(nr1*nr2*nr3), fr(nr1*nr2*nr3)
  complex(dp) :: fg(ngm)

  amn=cmplx_0
  !
  G0=0
  !
  call find_iG(G0,iG0)
  !
  ! find its scalar FFT index
  !
  iG0_fft=nl(iG0)
  !
  do ibnd=1,num_bands
     !
     fr(:)=cmplx_0
     !
     do is=1,nspin
        !
        call wfc_from_g_to_r (list_iG,wfc(:,ibnd,is),wfc_r)
        call func_from_g_to_r (1,ngm,guiding_function,fr)
        !
        do ir=1,nr1*nr2*nr3
           !
           fr(ir)=fr(ir) + conjg(wfc_r(ir))*fr(ir)
           !
        enddo !ir
     enddo !is
     !
     call func_from_r_to_g (1,ngm,fr,fg)
     !
     amn(ibnd)=fg(iG0_fft)
     !
  enddo !ibnd
  !
  return

  end subroutine get_guiding_function_overlap_FFT
!----------------------------------------------------------------------------------------
!****************************************************************************************
!----------------------------------------------------------------------------------------
  subroutine get_guiding_function_overlap_convolution(list_iG,wfc,guiding_function,amn)
!----------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
!     This subroutine computes the overlap between a given wavefunction
!     wfc and a given guiding_function (assumed normalized)
!
!	                amn(band) =  < wfc(band) |  guiding_function > .
!
!     The computation is done over all bands.
!
!     The G-vectors are referenced by their indices in list_iG
!     which refer to the global list gvec(3,ngm).
!
!--------------------------------------------------------------------------

  use intw_reading, only: ngm, nG_max
  use w90_parameters, only: num_bands

  implicit none

  !I/O variables

  integer,intent(in) :: list_iG(nG_max)
  complex(dp),intent(in) :: wfc(nG_max,num_bands,nspin)
  complex(dp),intent(in) :: guiding_function(ngm)
  complex(dp),intent(out) :: amn(num_bands)

  !local variables

  integer :: i, ibnd, is, iG, nG_max_non_zero
  integer        :: list_iG1(nG_max)
  complex(dp) :: amn_local(num_bands)
!  real(dp)       :: wrong_norm2

  amn = cmplx_0
  !
  nG_max_non_zero = 0
  !
  do i=1,nG_max
     !
     iG=list_iG(i)
     !
     ! list_iG is zero-padded at the end.
     ! When this point is reached, the computation is over.
     !
     if (iG==0) exit
     !
     nG_max_non_zero=nG_max_non_zero+1
     list_iG1(nG_max_non_zero)=iG
     !
  enddo !i
  !
  !$omp parallel default(none)                          &
  !$omp shared(nG_max_non_zero,nG_max,list_iG,list_iG1,cmplx_0)  &
  !$omp shared(nbands,num_bands,nspin,wfc,amn,guiding_function)   &
  !$omp private(i,iG,ibnd,is,amn_local)
  !
  amn_local = cmplx_0
  !
  ! First, build the pw_mat_el_local arrays, on each thread.
  !$omp do
  !
  do i=1,nG_max_non_zero
     !
     iG=list_iG1(i)
     !
     ! note that the guiding function is defined for all G in gvec;
     ! it is thus properly indexed by iG, not i.
     !
     do ibnd=1,num_bands
        do is=1,nspin
           !
           amn_local(ibnd)=amn_local(ibnd)+CONJG(wfc(i,ibnd,is))*guiding_function(iG)
           !
        enddo !is
     enddo !ibnd
     !
  enddo !i loop
  !
  !$omp end do
  !
  do ibnd=1,num_bands
     !
     !$omp atomic
     !
     amn(ibnd)=amn(ibnd)+amn_local(ibnd)
     !
  enddo
  !
  !$omp end parallel
  !
  return

  end subroutine get_guiding_function_overlap_convolution
!--------------------------------------------------------------------------
!**************************************************************************
!--------------------------------------------------------------------------
  subroutine get_radial_part(proj_nr,zona,k_plus_G_cart,guiding_function)
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
!
!     This subroutine computes the overlap between a given wavefunction
!     wfc and a given guiding_function (assumed normalized)
!
!	                amn(band) =  < wfc(band) |  guiding_function > .
!
!     The computation is done over all bands.
!
!     The G-vectors are referenced by their indices in list_iG
!     which refer to the global list gvec(3,ngm).
!
!--------------------------------------------------------------------------------

  use intw_reading, only: ngm

  implicit none

  !I/O variables

  integer,intent(in) :: proj_nr
  real(dp),intent(in) :: zona, k_plus_G_cart(3,ngm)
  complex(dp),intent(inout) :: guiding_function(ngm)

  !local variables

  integer :: iG
  real(dp) :: z, z2, z4, z6, z52, sqrt_z, pref
  real(dp) :: p(ngm)

  z=zona
  z2=z*z
  z4=z2*z2
  z6=z4*z2
  sqrt_z=sqrt(z)
  z52=z2*sqrt_z
  !
  !find the norms
  !
  do iG=1,ngm
     !
     p(iG)=sqrt( k_plus_G_cart(1,iG)**2 &
               + k_plus_G_cart(2,iG)**2 &
               + k_plus_G_cart(3,iG)**2 )
     !
  enddo !iG
  !
  ! the radial functions
  ! Their functional forms are obtained from analytical integrals
  ! done using MATHEMATICA (this is not a guarantee, however! always
  ! check for bugs...).
  !
  if (proj_nr==1) then
     !
     pref=4.0_dp*z52
     !
     guiding_function=pref/(p**2+z2)**2
     !
     ! CHEAT
     ! there appears to be a BUG in pw2wannier, on line 2106.
     ! There, the numerical integration of r R_{nl}(r) j_l(kr) is done,
     ! BUT IT IS THE INTEGRAL OF r^2 R_{nl}(r) j_l(kr) WHICH IS NEEDED.
     ! Below is the analytical expression for this WRONG integral;
     ! it yields better agreement between Amn computed by this code and
     ! pw2wannier.
     ! guiding_function=2.0_dp*z2*sqrt_z/(p**2+z2)
     !
  elseif (proj_nr==2) then
     !
     pref=16.0_dp*sqrt(2.0_dp)*z52
     !
     guiding_function=pref*(4.0_dp*p**2-z2)/(4.0_dp*p**2+z2)**3
     !
  elseif (proj_nr==3) then
     !
     pref=36.0_dp*sqrt(3.0_dp)*z52
     !
     guiding_function=pref*(81.0_dp*p**4-30.0_dp*p**2*z2+z4)  &
                                           /(9.0_dp*p**2+z2)**4
     !
  elseif (proj_nr == 4) then
     !
     pref=128.0_dp*z52
     !
     guiding_function=pref*(4096.0_dp*p**6-1792.0_dp*p**4*z2  &
                                     + 112.0_dp*p**2*z4 - z6) &
                                          /(16.0_dp*p**2+z2)**5
     !
  else
     !
     write(*,*) 'ERROR in intw2W90: this radial projection is not implemented'
     !
  endif
  !
  return

  end subroutine get_radial_part
!---------------------------------------------------------------------------
!***************************************************************************
!---------------------------------------------------------------------------
  subroutine get_angular_part(proj_l,proj_m,xaxis,zaxis,k_plus_G_cart,ylm)
!
!--------------------------------------------------------------------------
! This subroutine computes appropriate spherical harmonic corresponding
! to
!
!                   amn(band) =  < wfc(band) |  guiding_function > .
!
! The computation is done over all bands.
!
! The G-vectors are referenced by their indices in list_iG
! which refer to the global list gvec(3,ngm).
!--------------------------------------------------------------------------

  use intw_reading, only: ngm

  implicit none

  !I/O variables

  integer,intent(inout) :: proj_l,proj_m
  real(dp),intent(inout) :: xaxis(3), zaxis(3)
  real(dp),intent(inout) :: k_plus_G_cart(3,ngm)
  real(dp),intent(inout) :: ylm(ngm)

  !local variables

  integer :: iG
  real(dp) :: yaxis(3),q(3,ngm),norm_y


  ! produce the yaxis using z cross x. THEY REALLY SHOULD BE ORTHOGONAL.
  !
  yaxis(1)=  zaxis(2)*xaxis(3) - zaxis(3)*xaxis(2)
  yaxis(2)= -zaxis(1)*xaxis(3) + zaxis(3)*xaxis(1)
  yaxis(3)=  zaxis(1)*xaxis(2) - zaxis(2)*xaxis(1)
  !
  norm_y=sqrt(yaxis(1)**2+yaxis(2)**2+yaxis(3)**2)
  !
  yaxis=yaxis/norm_y
  !
  ! project k_plus_G_cart onto these axes.
  !
  do iG=1,ngm
     !
     q(1,iG)=k_plus_G_cart(1,iG)*xaxis(1) &
           + k_plus_G_cart(2,iG)*xaxis(2) &
           + k_plus_G_cart(3,iG)*xaxis(3)
     !
     q(2,iG)=k_plus_G_cart(1,iG)*yaxis(1) &
           + k_plus_G_cart(2,iG)*yaxis(2) &
           + k_plus_G_cart(3,iG)*yaxis(3)
     !
     q(3,iG)=k_plus_G_cart(1,iG)*zaxis(1) &
           + k_plus_G_cart(2,iG)*zaxis(2) &
           + k_plus_G_cart(3,iG)*zaxis(3)
     !
  enddo !iG
  !
  call ylm_wannier(ylm,proj_l,proj_m,q,ngm)
  !
  return

  end subroutine get_angular_part
!----------------------------------------
!****************************************
!----------------------------------------
  subroutine ylm_wannier(ylm,l,mr,r,nr)
!----------------------------------------
!
!--------------------------------------------------------------------------
!     This subroutine has been taken from pw2wannier and modified
!     to suit my needs.
!     The original comments follow.
!--------------------------------------------------------------------------
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
! of the spherical harmonic identified  by indices (l,mr)
! in table 3.1 of the wannierf90 specification.
!
! No reference to the particular ylm ordering internal to quantum-espresso
! is assumed.
!
! If ordering in wannier90 code is changed or extended this should be the
! only place to be modified accordingly
!--------------------------------------------------------------------------

   implicit none

   !I/O variables

   integer,intent(inout) :: l, mr, nr
   real(dp),intent(inout) :: ylm(nr), r(3,nr)

   !local variables

   real(dp), external :: s, pz_func, px, py, dz2, dxz, dyz, dx2my2, dxy, &
                         fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(dp) :: rr, cost, phi
   integer :: ir
   real(dp) :: bs2, bs3, bs6, bs12

   bs2  = 1.d0/sqrt(2.d0)
   bs3  = 1.d0/sqrt(3.d0)
   bs6  = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)
   !
   if ( l>3 .OR. l<-5 ) then
      !
      write(*,*)' ylm_wannier: l out of range! '
      stop
      !
   endif
   !
   if (l>=0) then
      !
      if ( mr<1 .OR. mr>2*l+1 ) then
         !
         write(*,*)' ylm_wannier: m out of range! '
         stop
         !
      endif
      !
   else
      !
      if ( mr<1 .OR. mr>abs(l)+1 ) then
         !
         write(*,*)' ylm_wannier: m out of range! '
         stop
         !
      endif
      !
   endif
   !
   do ir=1,nr
      !
      rr=sqrt( r(1,ir)*r(1,ir) + r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      !
      cost=r(3,ir)/rr
      !
      ! beware the arc tan, it is defined modulo pi
      !
      if (r(1,ir) > eps_8) then
         !
         phi=atan( r(2,ir)/r(1,ir) )
         !
      elseif (r(1,ir) < -eps_8 ) then
         !
         phi=atan( r(2,ir)/r(1,ir) ) + pi
         !
      else
         !
         phi = sign( pi/2.d0,r(2,ir) )
         !
      endif
      !
      ! if the norm of r is very small, just arbitrarily pick
      ! the angle to be theta = 0 , phi = 0
      !
      if (rr < eps_8) then
         !
         cost=1.0_dp
         phi=0.0_dp
         !
      endif
      !
      if (l==0) then   ! s orbital
         !
         ylm(ir)=s(cost,phi)
         !
      endif
      !
      if (l==1) then   ! p orbitals
         !
         if (mr==1) ylm(ir)=pz_func(cost,phi)
         if (mr==2) ylm(ir)=px(cost,phi)
         if (mr==3) ylm(ir)=py(cost,phi)
         !
      endif
      !
      if (l==2) then   ! d orbitals
         !
         if (mr==1) ylm(ir)=dz2(cost,phi)
         if (mr==2) ylm(ir)=dxz(cost,phi)
         if (mr==3) ylm(ir)=dyz(cost,phi)
         if (mr==4) ylm(ir)=dx2my2(cost,phi)
         if (mr==5) ylm(ir)=dxy(cost,phi)
         !
      endif
      !
      if (l==3) then   ! f orbitals
         !
         if (mr==1) ylm(ir)=fz3(cost,phi)
         if (mr==2) ylm(ir)=fxz2(cost,phi)
         if (mr==3) ylm(ir)=fyz2(cost,phi)
         if (mr==4) ylm(ir)=fzx2my2(cost,phi)
         if (mr==5) ylm(ir)=fxyz(cost,phi)
         if (mr==6) ylm(ir)=fxx2m3y2(cost,phi)
         if (mr==7) ylm(ir)=fy3x2my2(cost,phi)
         !
      endif
      !
      if (l==-1) then  !  sp hybrids
         !
         if (mr==1) ylm(ir)=bs2 * ( s(cost,phi) + px(cost,phi) )
         if (mr==2) ylm(ir)=bs2 * ( s(cost,phi) - px(cost,phi) )
         !
      endif
      !
      if (l==-2) then  !  sp2 hybrids
         !
         if (mr==1) ylm(ir)= bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         if (mr==2) ylm(ir)= bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         if (mr==3) ylm(ir)= bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
         !
      endif
      !
      if (l==-3) then  !  sp3 hybrids
         !
         if (mr==1) ylm(ir)=0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+pz_func(cost,phi))
         if (mr==2) ylm(ir)=0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-pz_func(cost,phi))
         if (mr==3) ylm(ir)=0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-pz_func(cost,phi))
         if (mr==4) ylm(ir)=0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+pz_func(cost,phi))
         !
      endif
      !
      if (l==-4) then  !  sp3d hybrids
         !
         if (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         if (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         if (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
         if (mr==4) ylm(ir) = bs2*pz_func(cost,phi)+bs2*dz2(cost,phi)
         if (mr==5) ylm(ir) =-bs2*pz_func(cost,phi)+bs2*dz2(cost,phi)
         !
      endif
      !
      if (l==-5) then  ! sp3d2 hybrids
         !
         if (mr==1) ylm(ir) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
         if (mr==2) ylm(ir) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
         if (mr==3) ylm(ir) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
         if (mr==4) ylm(ir) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
         if (mr==5) ylm(ir) = bs6*s(cost,phi)-bs2*pz_func(cost,phi)+bs3*dz2(cost,phi)
         if (mr==6) ylm(ir) = bs6*s(cost,phi)+bs2*pz_func(cost,phi)+bs3*dz2(cost,phi)
         !
      endif
      !
   enddo !ir
   !
   return

   end subroutine ylm_wannier
!----------------------------------------------------------------------------!
!
!
end module intw_intw2wannier
!
!
!----------------------------------------------------------------------------!


!======== l = 0 =====================================================================
function s(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) :: s, cost, phi

   s = 1.d0/ sqrt(fpi)
   return
end function s

!======== l = 1 =====================================================================
function pz_func(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::pz_func, cost,phi
   pz_func =  sqrt(3.d0/fpi) * cost
   return
end function pz_func

function px(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::px, cost, phi, sint

   sint = sqrt(abs(1.d0 - cost*cost))
   px   =  sqrt(3.d0/fpi) * sint * cos(phi)
   return
end function px

function py(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py   =  sqrt(3.d0/fpi) * sint * sin(phi)
   return
end function py
!======== l = 2 =====================================================================
function dz2(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::dz2, cost, phi
   dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)
   return
end function dz2

function dxz(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)
   return
end function dxz

function dyz(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz  =  sqrt(15.d0/fpi) * sint*cost * sin(phi)
   return
end function dyz

function dx2my2(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::dx2my2, cost, phi, sint
   sint   = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)
   return
end function dx2my2

function dxy(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy  =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)
   return
end function dxy
!======== l = 3 =====================================================================
function fz3(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::fz3, cost, phi
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
   return
end function fz3

function fxz2(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
   return
end function fxz2

function fyz2(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
   return
end function fyz2

function fzx2my2(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
   return
end function fzx2my2

function fxyz(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
   return
end function fxyz

function fxx2m3y2(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
   return
end function fxx2m3y2

function fy3x2my2(cost,phi)
  use kinds, only: dp
  use intw_useful_constants, only: pi,tpi,fpi
   implicit none
   real(dp) ::fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
   return
end function fy3x2my2
