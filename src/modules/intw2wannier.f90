!
! Copyright (C) 2024 INTW group
!
! This file is part of INTW.
!
! INTW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! INTW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
module intw_intw2wannier

  !---------------------------------------------------------------------------!
  ! This module contains subroutines which will perform the same tasks        !
  ! as the Quantum Espresso utility program "pw2wannier", but utilizing       !
  ! symmetries.                                                               !
  !                                                                           !
  ! In particular, this module will produce the $seed.mmn and $seed.amn files !
  ! using only the wavefunctions in the IBZ, rotating them appropriately.     !
  !---------------------------------------------------------------------------!

  use kinds, only: dp
  use intw_reading, only: nbands, nGk_max, ngm, nspin
  use intw_useful_constants, only: bohr, pi, tpi, fpi, eps_8, ZERO, cmplx_0, cmplx_i

  implicit none
  !
  ! variables
  public :: nnkp_exclude_bands, nnkp_real_lattice, nnkp_recip_lattice, &
            nnkp_num_kpoints, nnkp_nnkpts, nnkp_n_proj, nnkp_kpoints, &
            nnkp_Wcenters, nnkp_proj_x, nnkp_proj_z, nnkp_proj_zona, &
            nnkp_proj_n, nnkp_proj_l, nnkp_proj_m, &
            nnkp_proj_s, nnkp_proj_spin_axis, &
            nnkp_list_ikpt_nn, nnkp_list_G, nnkp_excluded_bands
  !
  ! subroutines
  public :: deallocate_nnkp, read_nnkp_file, output_nnkp_file, &
            intw2W90_check_mesh, generate_header,  &
            generate_mmn_using_allwfc, generate_amn_using_allwfc, &
            generate_guiding_function, &
            get_guiding_function_overlap_FFT, &
            get_guiding_function_overlap_convolution, get_radial_part, &
            get_angular_part, ylm_wannier  !, get_bvec_list
  !
  private
  !
  save
  ! declare some global variables; the prefix "nnkp" will indicate
  ! variables read from the $prefix.nnkp file, which must be checked with
  ! the data read from the input file for consistency.

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

  integer     :: nnkp_n_proj            !the number of projections specified
                                        !(should be the same as number of Wannier functions)

  real(dp),allocatable :: nnkp_kpoints(:,:)   ! the k vectors in the 1BZ
  real(dp),allocatable :: nnkp_Wcenters(:,:)  ! the Wannier centers, for trial projection
  real(dp),allocatable :: nnkp_proj_x(:,:)    ! x projection, for trial
  real(dp),allocatable :: nnkp_proj_z(:,:)    ! z projection, for trial
  real(dp),allocatable :: nnkp_proj_zona(:)   ! Z on a: gauge broadness of the hydrogenoid trial wfc.

  ! The hydrogenoid quantum numbers, psi_{nlm}
  integer,allocatable  :: nnkp_proj_n(:)     ! the radial-pojection, for trial
  integer,allocatable  :: nnkp_proj_l(:)     ! the l-pojection, for trial
  integer,allocatable  :: nnkp_proj_m(:)     ! the m-pojection, for trial

  ! JLB: Extra variables for spinor projections
  real(dp),allocatable :: nnkp_proj_spin_axis(:,:) ! spin quantization axis for spinor projection
  integer,allocatable  :: nnkp_proj_s(:) ! up or down spin

  integer,allocatable  :: nnkp_list_ikpt_nn(:,:) ! the neighbors of ikpt1
  integer,allocatable  :: nnkp_list_G(:,:,:)   ! the G vectors linking one point to another

  logical,allocatable  :: nnkp_excluded_bands(:) ! what bands are excluded


contains
  !
!---------------------------------
  subroutine deallocate_nnkp()
!---------------------------------

  use intw_reading, only: noncolin

  implicit none

  deallocate(nnkp_kpoints)
  deallocate(nnkp_Wcenters)
  deallocate(nnkp_proj_l)
  deallocate(nnkp_proj_m)
  deallocate(nnkp_proj_n)
  deallocate(nnkp_proj_x)
  deallocate(nnkp_proj_z)
  deallocate(nnkp_proj_zona)
  if (noncolin) deallocate(nnkp_proj_spin_axis)
  if (noncolin) deallocate(nnkp_proj_s)
  deallocate(nnkp_list_ikpt_nn)
  deallocate(nnkp_list_G)
  deallocate(nnkp_excluded_bands)

  end subroutine deallocate_nnkp
!----------------------------------------------
!**********************************************
!------------------------------------------
  subroutine read_nnkp_file(nnkp_file)
!------------------------------------------
!
!----------------------------------------------------------------------------!
! This subroutine reads the .nnkp input file produced by W90.
! It assumes the file has been checked for existence and is already
! opened.
!----------------------------------------------------------------------------!

  use intw_reading, only: alat, noncolin, scan_file_to
  use intw_utility, only: find_free_unit

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
  if (noncolin) then
    call scan_file_to (nnkp_unit,'spinor_projections')
    read(nnkp_unit,*) nnkp_n_proj
    allocate(nnkp_proj_s(nnkp_n_proj))
    allocate(nnkp_proj_spin_axis(3,nnkp_n_proj))
  else
    call scan_file_to (nnkp_unit,'projections')
    read(nnkp_unit,*) nnkp_n_proj
  end if
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
     if (noncolin) then
       read(nnkp_unit,*) nnkp_proj_s(j), nnkp_proj_spin_axis(:,j)
       ! Check
       if (abs(nnkp_proj_s(j)) /= 1) then
         write(*,*) "Error in spinor projection! Should be +-1"
         write(*,*) "Stopping..."
         stop
       else if (abs(nnkp_proj_spin_axis(1,j))>eps_8 &
                .or. abs(nnkp_proj_spin_axis(2,j))>eps_8 &
                .or. abs(nnkp_proj_spin_axis(3,j)-1)>eps_8) then
         write(*,*) "Currently, only (/0, 0, 1/) spin-axis orientation implemented"
         write(*,*) "Stopping..."
         stop
       end if
       !
     end if
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
  ! JLB: Now this is done in reading.f90 -> set_num_bands
  !      Here only poulating nnkp_* variables,
  !      then check for consistency with intw_* in main program / utility
  !num_exclude_bands=nnkp_exclude_bands
  !allocate(exclude_bands(num_exclude_bands))
  !
  allocate(nnkp_excluded_bands(nbands))
  !
  nnkp_excluded_bands(:)=.false.
  !
  do j=1,nnkp_exclude_bands
     !
     read(nnkp_unit,*) nn
     nnkp_excluded_bands(nn)=.true.
     !exclude_bands(j)=nn
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
  use intw_reading, only: noncolin
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

      if (noncolin) then
         write(io_unit,'(I6, 3F12.6)') nnkp_proj_s(j), &
                                       (nnkp_proj_spin_axis(i,j), i=1,3)
      end if

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


  subroutine generate_mmn_using_allwfc (intw2W_fullzone,method)
    !----------------------------------------------------------------------------!
    ! This subroutine computes the plane wave matrix elements  needed by Wannier90
    ! by using symmetry. It fetches the wfc in the IBZ, rotates them, and computes
    ! the needed matrix elements.
    !
    !----------------------------------------------------------------------------!
    use intw_allwfcs, only: get_psi_general_k_all_wfc
    use intw_utility, only: find_free_unit
    use intw_matrix_elements, only: get_plane_wave_matrix_element_FFT, get_plane_wave_matrix_element_convolution_map
    use intw_input_parameters, only: outdir, prefix, nk1, nk2, nk3
    use intw_reading, only: num_bands_intw

    implicit none

    logical        :: intw2W_fullzone
    character(*)   :: method

    integer        :: io_unit_mmn, io_unit_eig
    integer        :: nkmesh

    integer        :: nn

    character(256) :: filename

    integer        :: ikpt_1, ikpt_2
    integer        :: nb1, nb2, nb
    integer        :: G(3)

    integer        :: ngk1, ngk2

    integer        :: list_iG_1(nGk_max), list_iG_2(nGk_max)
    complex(dp)    :: wfc_1(nGk_max,num_bands_intw,nspin), wfc_2(nGk_max,num_bands_intw,nspin)
    real(dp)       :: QE_eig(num_bands_intw)

    complex(dp)    :: pw_mat_el(num_bands_intw,num_bands_intw,nspin,nspin)

    character(256) :: header


    !-----------------------------------
    ! Open all the needed files
    !-----------------------------------
    io_unit_mmn = find_free_unit()
    filename = trim(outdir)//trim(prefix)//trim('.mmn')
    open(unit=io_unit_mmn,file=filename,status='unknown')

    io_unit_eig = find_free_unit()
    filename = trim(outdir)//trim(prefix)//trim('.eig')
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
      call get_psi_general_k_all_wfc(nnkp_kpoints(:,ikpt_1), ngk1, list_iG_1, wfc_1, QE_eig)

      ! print out the eigenvalues
      do nb =1, num_bands_intw
        write(io_unit_eig,'(2(I6,x),F18.12)') nb, ikpt_1, QE_eig(nb)
      end do

      ! loop on neighbors
      do nn = 1, nnkp_nnkpts
        G      = nnkp_list_G(:,nn,ikpt_1)
        ikpt_2 = nnkp_list_ikpt_nn(nn,ikpt_1)

        write(io_unit_mmn,'(5I7)')  ikpt_1,ikpt_2,G

        ! fetch data
        call get_psi_general_k_all_wfc(nnkp_kpoints(:, ikpt_2) + G, ngk2, list_iG_2, wfc_2, QE_eig)

        ! Compute the matrix elements
        if ( trim(method) == 'CONVOLUTION' ) then
          call get_plane_wave_matrix_element_convolution_map      &
                        ((/0, 0, 0/),list_iG_1,ngk1,list_iG_2,ngk2, wfc_1,wfc_2,pw_mat_el)


        else if ( trim(method) == 'FFT' ) then
          call get_plane_wave_matrix_element_FFT              &
                        ((/0, 0, 0/),list_iG_1,list_iG_2, wfc_1,wfc_2,pw_mat_el)
        else
          write(*,*) 'ERROR in generate_mmn'
          stop
        end if

        if (nspin==2) then
          do nb1 = 1,num_bands_intw
            do nb2 = 1,num_bands_intw
              write(io_unit_mmn,'(2F18.12)')   pw_mat_el(nb2,nb1,1,1) + pw_mat_el(nb2,nb1,2,2)
            end do
          end do
        else if  (nspin==1) then
          do nb1 = 1,num_bands_intw
            do nb2 = 1,num_bands_intw
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

    use intw_allwfcs, only: get_psi_general_k_all_wfc
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: cmplx_0
    use intw_reading, only : noncolin, num_bands_intw
    use intw_input_parameters, only: outdir, prefix, nk1, nk2, nk3

    implicit none

    logical        :: intw2W_fullzone
    character(256) :: method

    integer        :: io_unit_amn
    integer        :: nkmesh

    character(256) :: filename

    integer        :: ikpt
    integer        :: nb, n_proj

    integer        :: ngk, list_iG(nGk_max)

    complex(dp)    :: wfc(nGk_max,num_bands_intw,nspin)

    complex(dp)    :: guiding_function(nGk_max,nspin)

    real(dp)       :: QE_eig(num_bands_intw)

    complex(dp)    :: amn(num_bands_intw)

    character(256) :: header


    nkmesh = nk1*nk2*nk3

    io_unit_amn = find_free_unit()
    filename    = trim(outdir)//trim(prefix)//trim('.amn')
    open(unit=io_unit_amn,file=filename,status='unknown')

    call generate_header(method,header)
    write(io_unit_amn,*) trim(header)
    write(io_unit_amn,'(3I12)') nbands-nnkp_exclude_bands, nnkp_num_kpoints , nnkp_n_proj

    !loop on all k-points
    do ikpt = 1, nnkp_num_kpoints

      ! fetch the data
      call get_psi_general_k_all_wfc(nnkp_kpoints(:, ikpt), ngk, list_iG, wfc, QE_eig)

      !loop on all bands and all trial functions
      do n_proj = 1,nnkp_n_proj

        ! Generate the fourier transform of the trial function (called guiding
        ! function, just like in pw2wannier).
        call generate_guiding_function(ikpt, ngk, list_iG, n_proj, guiding_function(:,1))

        !JLB spinor projection. Should be generalized to quantization axis /= z
        if (noncolin) then
          if (nnkp_proj_s(n_proj) < 0) then
            guiding_function(:,2) = guiding_function(:,1)
            guiding_function(:,1) = cmplx_0
          else
            guiding_function(:,2) = cmplx_0
          end if
        end if

        if (trim(method) == 'CONVOLUTION') then
          call get_guiding_function_overlap_convolution(ngk, wfc, guiding_function, amn)
        else if (trim(method) == 'FFT') then
          call get_guiding_function_overlap_FFT(list_iG, wfc, guiding_function, amn)
        else
          write(*,*) 'ERROR in generate_amn'
          stop
        end if

        !Write result to file $prefix.amn.
        do nb = 1,num_bands_intw
          write(io_unit_amn,'(3I7,2F18.12)') nb, n_proj, ikpt, amn(nb)
        end do !nb

      end do ! n_proj

    end do !ikpt

    close(io_unit_amn)

  end subroutine generate_amn_using_allwfc


  subroutine generate_guiding_function(ikpt, ngk, list_iG, n_proj, guiding_function)
    !---------------------------------------------------------------------------
    ! This subroutine computes the normalized guiding function in reciprocal
    ! space.
    !
    ! This subroutine is heavily inspired by a similar routine in pw2wannier.
    !---------------------------------------------------------------------------

    use intw_reading, only: gvec, alat

    implicit none

    !I/O variables

    integer, intent(in) :: ikpt, ngk, n_proj, list_iG(nGk_max)
    complex(dp), intent(out) :: guiding_function(nGk_max)

    !local variables

    integer, parameter :: lmax=3, lmmax=(lmax+1)**2 ! max Ylm implemented in wannier90
    integer :: i, mu, iG
    integer :: proj_nr, proj_l, proj_m
    integer :: l, m, lm
    real(dp) :: zona, zaxis(3), xaxis(3), ylm(nGk_max)
    real(dp) :: k_cryst(3), tau_cryst(3), tau_cart(3)
    real(dp) :: k_plus_G_cart(3,ngk)
    real(dp) :: norm2
    real(dp) :: radial_l(nGk_max, 0:lmax), coef(lmmax)
    complex(dp) :: four_pi_i_l

    !
    ! get all the relevant vectors in reciprocal space
    ! express the vectors in bohrs^-1 cartesian coordinates
    !
    k_cryst(:)=nnkp_kpoints(:,ikpt)
    k_plus_G_cart=ZERO
    !
    do iG=1,ngk
      do mu=1,3
        do i=1,3
          !
          k_plus_G_cart(mu,iG) = k_plus_G_cart(mu,iG) &
                                 + nnkp_recip_lattice(mu,i)*(k_cryst(i)+dble(gvec(i,list_ig(iG))))
          !
        enddo !i
      enddo !mu
    enddo !iG
    !
    k_plus_G_cart = k_plus_G_cart*tpi/alat ! bohr^-1
    !
    ! compute the guiding function
    !
    guiding_function = cmplx_0
    !
    ! Set projection l and m
    proj_l = nnkp_proj_l(n_proj)
    proj_m = nnkp_proj_m(n_proj)
    ! JLB note: Not sure whether this works currently, needs to be tested
    zaxis(:) = nnkp_proj_z(:,n_proj)
    xaxis(:) = nnkp_proj_x(:,n_proj)
    !
    ! JLB: Expansion coefficients of this projection in lm-s (needed for hybrids)
    call projection_expansion(proj_l, proj_m, coef)
    !
    ! get the part from the radial integration
    proj_nr = nnkp_proj_n(n_proj) ! the radial projection parameter
    zona = nnkp_proj_zona(n_proj) ! Z/a, the diffusive parameter
    !
    ! MBR, JLB: radial integrals based on pw2wannier
    !call get_radial_part(proj_nr,zona,k_plus_G_cart,guiding_function)
    call get_radial_part_numerical(lmax, coef, proj_nr, zona, ngk, k_plus_G_cart, radial_l)

    !
    do l=0, lmax
      do m=1, 2*l+1
        !
        ! Check which lm are needed in this projection
        lm = l**2 + m
        if (abs(coef(lm)) < 1.0d-8) cycle
        !
        call get_angular_part(l, m, xaxis, zaxis, ngk, k_plus_G_cart, ylm)
        !
        four_pi_i_l = fpi*(-cmplx_i)**l
        !
        do iG=1,ngk
          !
          guiding_function(iG) = guiding_function(iG) + coef(lm)*radial_l(iG, l)*ylm(iG)*four_pi_i_l
          !
        enddo !iG
        !
      end do !m
    end do !l
    !
    ! Add extra phase coming from Wannier center positions
    tau_cryst = nnkp_Wcenters(:,n_proj) ! crystal coordinates
    tau_cart = ZERO
    do mu=1,3
      do i=1,3
        !
        tau_cart(mu) = tau_cart(mu) + alat*nnkp_real_lattice(mu,i)*tau_cryst(i)
        !
      enddo !i
    enddo !mu
    !
    do iG=1,ngk
      guiding_function(ig) = guiding_function(ig) &
                             * exp(-cmplx_i*(  k_plus_G_cart(1,ig)*tau_cart(1) &
                                             + k_plus_G_cart(2,ig)*tau_cart(2) &
                                             + k_plus_G_cart(3,ig)*tau_cart(3) ))
    end do !iG
    !
    ! Normalize
    norm2 = sum(abs(guiding_function)**2)
    guiding_function = guiding_function/sqrt(norm2)
    !
    return

  end subroutine generate_guiding_function


  subroutine get_guiding_function_overlap_FFT(list_iG, wfc, guiding_function, amn)
    !------------------------------------------------------------------------
    !  This subroutine computes the overlap between a given wavefunction
    !  wfc and a given guiding_function (assumed normalized)
    !
    !                    amn(band) =  < wfc(band) |  guiding_function > .
    !
    !  The computation is done over all bands using FFT.
    !------------------------------------------------------------------------

    use intw_reading, only: nGk_max, nr1, nr2, nr3, num_bands_intw
    use intw_fft, only: wfc_from_g_to_r, wfc_from_r_to_g

    implicit none

    !I/O variables

    integer,intent(in) :: list_iG(nGk_max)
    complex(dp),intent(in) :: wfc(nGk_max,num_bands_intw,nspin)
    complex(dp),intent(in) :: guiding_function(nGk_max,nspin)
    complex(dp),intent(out) :: amn(num_bands_intw)

    !local variables

    integer :: ibnd, ir, is
    complex(dp) :: wfc_r(nr1*nr2*nr3), fr(nr1*nr2*nr3)
    complex(dp) :: fg(nGk_max)


    amn = cmplx_0
    !
    do ibnd=1,num_bands_intw
      !
      do is=1,nspin
        !
        call wfc_from_g_to_r(list_iG, guiding_function(:,is), fr)
        call wfc_from_g_to_r(list_iG, wfc(:,ibnd,is), wfc_r)
        !
        do ir=1,nr1*nr2*nr3
          !
          fr(ir) = conjg(wfc_r(ir))*fr(ir)
          !
        enddo !ir
        !
        call wfc_from_r_to_g(list_iG, fr, fg)
        !
        amn(ibnd) = amn(ibnd) + fg(1)
        !
      enddo !is
      !
    enddo !ibnd
    !
    return

  end subroutine get_guiding_function_overlap_FFT


  subroutine get_guiding_function_overlap_convolution(ngk, wfc, guiding_function, amn)
    !--------------------------------------------------------------------------
    !     This subroutine computes the overlap between a given wavefunction
    !     wfc and a given guiding_function (assumed normalized)
    !
    !	                amn(band) =  < wfc(band) |  guiding_function > .
    !
    !     The computation is done over all bands.
    !--------------------------------------------------------------------------

    use intw_reading, only: nGk_max, num_bands_intw

    implicit none

    !I/O variables

    integer, intent(in) :: ngk
    complex(dp), intent(in) :: wfc(nGk_max,num_bands_intw,nspin)
    complex(dp), intent(in) :: guiding_function(nGk_max,nspin)
    complex(dp), intent(out) :: amn(num_bands_intw)

    !local variables

    integer :: ibnd, is, iG
    complex(dp) :: amn_local(num_bands_intw)


    amn = cmplx_0
    !
    !
    !$omp parallel default(none) &
    !$omp shared(num_bands_intw,nspin,wfc,amn,guiding_function,ngk) &
    !$omp private(iG,ibnd,is,amn_local)
    !
    amn_local = cmplx_0
    !
    ! First, build the pw_mat_el_local arrays, on each thread.
    !$omp do
    !
    do iG=1,ngk
      !
      do ibnd=1,num_bands_intw
          do is=1,nspin
            !
            amn_local(ibnd) = amn_local(ibnd) + CONJG(wfc(iG,ibnd,is))*guiding_function(iG,is)
            !
          enddo !is
      enddo !ibnd
      !
    enddo !i loop
    !
    !$omp end do
    !
    do ibnd=1,num_bands_intw
      !
      !$omp atomic
      !
      amn(ibnd) = amn(ibnd) + amn_local(ibnd)
      !
    enddo
    !
    !$omp end parallel
    !
    return

  end subroutine get_guiding_function_overlap_convolution


  subroutine get_radial_part_numerical(lmax, coef, proj_nr, zona, ngk, k_plus_G_cart, radial_l)
    ! MBR
    ! Numerical integration
    ! JLB: Extended to multiple l, needed for hybrid projections
    !
    ! This subroutine is based on radialpart subroutine distributed as part of the Quantum Espresso
    ! code and has been adapted to INTW:
    !   Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use intw_reading, only: nGk_max, volume0
    use intw_utility, ONLY : simpson, sphb
    use intw_useful_constants, only: fpi, ZERO

    implicit none

    !I/O variables

    integer, intent(in) :: lmax, proj_nr, ngk
    real(dp), intent(in) :: coef((lmax+1)**2), zona, k_plus_G_cart(3,ngk)
    real(dp), intent(out) :: radial_l(nGk_max, 0:lmax)

    !local variables

    integer :: iG, l
    real(dp) :: z, z2, z4, z6, z52, sqrt_z
    real(dp) :: p(ngk)

    ! from pw2wannier

    integer :: mesh_r, ir
    real(DP), PARAMETER :: xmin=-6.d0, dx=0.025d0, rmax=10.d0
    real(DP) :: x, rad_int
    real(DP), ALLOCATABLE :: bes(:), func_r(:), r(:), rij(:), aux(:)


    z = zona
    z2 = z*z
    z4 = z2*z2
    z6 = z4*z2
    sqrt_z = sqrt(z)
    z52 = z2*sqrt_z
    !
    !find the norms
    !
    do iG=1,ngk
     !
     p(iG) = sqrt( k_plus_G_cart(1,iG)**2 &
                 + k_plus_G_cart(2,iG)**2 &
                 + k_plus_G_cart(3,iG)**2 )
     !
    enddo !iG

    !
    ! from pw2intw:
    !
    mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
    allocate(bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r))
    allocate(aux(mesh_r))
    !
    !    compute the radial mesh
    !
    do ir = 1, mesh_r
     x = xmin + dble(ir - 1) * dx
     r(ir) = exp(x) / zona
     rij(ir) = dx * r(ir)
    end do
    !
    !
    if (proj_nr==1) then
      func_r(:) = 2.d0 * zona**(3.d0/2.d0) * exp(-zona*r(:))
    else if (proj_nr==2) then
      func_r(:) = 1.d0/sqrt(8.d0) * zona**(3.d0/2.d0) * &
                  (2.0d0 - zona*r(:)) * exp(-zona*r(:)*0.5d0)
    else if (proj_nr==3) then
      func_r(:) = sqrt(4.d0/27.d0) * zona**(3.0d0/2.0d0) * &
                (1.d0 - 2.0d0/3.0d0*zona*r(:) + 2.d0*(zona*r(:))**2/27.d0) * exp(-zona*r(:)/3.0d0)
    else
      write(*,*) 'ERROR in intw2W90: this radial projection is not implemented'
    endif

    radial_l = ZERO
    do l=0, lmax
      ! JLB: Check which l-s are used in this projection
      if ( any(coef(l**2+1:l**2+1+2*l) > 1.0d-8) ) then
        !
        do iG=1,ngk
          aux = r*r*sphb(l,p(iG)*r)*func_r
          call simpson(mesh_r, aux, rij, rad_int)
          radial_l(iG, l) = rad_int * fpi/sqrt(volume0)
        end do
        !
      end if
    end do

    deallocate(bes,func_r,r,rij,aux)

    return

  end subroutine get_radial_part_numerical
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
  subroutine get_angular_part(proj_l, proj_m, xaxis, zaxis, ngk, k_plus_G_cart, ylm)
    !--------------------------------------------------------------------------
    ! This subroutine computes appropriate spherical harmonic corresponding
    ! to
    !
    !                   amn(band) =  < wfc(band) |  guiding_function > .
    !
    ! The computation is done over all bands.
    !--------------------------------------------------------------------------

    use intw_reading, only: nGk_max

    implicit none

    !I/O variables

    integer, intent(in) :: ngk, proj_l, proj_m
    real(dp), intent(in) :: xaxis(3), zaxis(3)
    real(dp), intent(in) :: k_plus_G_cart(3,ngk)
    real(dp), intent(out) :: ylm(nGk_max)

    !local variables

    integer :: iG
    real(dp) :: yaxis(3), q(3,ngk), norm_y


    ! produce the yaxis using z cross x. THEY REALLY SHOULD BE ORTHOGONAL.
    !
    yaxis(1) =  zaxis(2)*xaxis(3) - zaxis(3)*xaxis(2)
    yaxis(2) = -zaxis(1)*xaxis(3) + zaxis(3)*xaxis(1)
    yaxis(3) =  zaxis(1)*xaxis(2) - zaxis(2)*xaxis(1)
    !
    norm_y = sqrt(yaxis(1)**2+yaxis(2)**2+yaxis(3)**2)
    !
    yaxis=yaxis/norm_y
    !
    ! project k_plus_G_cart onto these axes.
    !
    do iG=1,ngk
      !
      q(1,iG) = k_plus_G_cart(1,iG)*xaxis(1) &
              + k_plus_G_cart(2,iG)*xaxis(2) &
              + k_plus_G_cart(3,iG)*xaxis(3)
      !
      q(2,iG) = k_plus_G_cart(1,iG)*yaxis(1) &
              + k_plus_G_cart(2,iG)*yaxis(2) &
              + k_plus_G_cart(3,iG)*yaxis(3)
      !
      q(3,iG) = k_plus_G_cart(1,iG)*zaxis(1) &
              + k_plus_G_cart(2,iG)*zaxis(2) &
              + k_plus_G_cart(3,iG)*zaxis(3)
      !
    enddo !iG
    !
    call ylm_wannier(ylm,proj_l,proj_m,q,ngk)
    !
    return

  end subroutine get_angular_part


  subroutine ylm_wannier(ylm,l,mr,r,nr)
    !
    ! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
    ! of the spherical harmonic identified  by indices (l,mr)
    ! in table 3.1 of the wannierf90 specification.
    !
    ! No reference to the particular ylm ordering internal to quantum-espresso
    ! is assumed.
    !
    ! If ordering in wannier90 code is changed or extended this should be the
    ! only place to be modified accordingly
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    implicit none

    !I/O variables

    integer,intent(in) :: l, mr, nr
    real(dp), intent(in) :: r(3,nr)
    real(dp),intent(out) :: ylm(nr)

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
      rr = sqrt( r(1,ir)*r(1,ir) + r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      !
      cost = r(3,ir)/(rr+eps_8)
      !
      ! beware the arc tan, it is defined modulo pi
      !
      if (r(1,ir) > eps_8) then
        !
        phi = atan( r(2,ir)/r(1,ir) )
        !
      elseif (r(1,ir) < -eps_8) then
        !
        phi = atan( r(2,ir)/r(1,ir) ) + pi
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
        cost = 1.0_dp
        phi = 0.0_dp
        !
      endif
      !
      if (l==0) then ! s orbital
        !
        ylm(ir) = s(cost,phi)
        !
      endif
      !
      if (l==1) then ! p orbitals
        !
        if (mr==1) ylm(ir) = pz_func(cost,phi)
        if (mr==2) ylm(ir) = px(cost,phi)
        if (mr==3) ylm(ir) = py(cost,phi)
        !
      endif
      !
      if (l==2) then ! d orbitals
        !
        if (mr==1) ylm(ir) = dz2(cost,phi)
        if (mr==2) ylm(ir) = dxz(cost,phi)
        if (mr==3) ylm(ir) = dyz(cost,phi)
        if (mr==4) ylm(ir) = dx2my2(cost,phi)
        if (mr==5) ylm(ir) = dxy(cost,phi)
        !
      endif
      !
      if (l==3) then ! f orbitals
        !
        if (mr==1) ylm(ir) = fz3(cost,phi)
        if (mr==2) ylm(ir) = fxz2(cost,phi)
        if (mr==3) ylm(ir) = fyz2(cost,phi)
        if (mr==4) ylm(ir) = fzx2my2(cost,phi)
        if (mr==5) ylm(ir) = fxyz(cost,phi)
        if (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
        if (mr==7) ylm(ir) = fy3x2my2(cost,phi)
        !
      endif
      !
      if (l==-1) then ! sp hybrids
        !
        if (mr==1) ylm(ir) = bs2 * ( s(cost,phi) + px(cost,phi) )
        if (mr==2) ylm(ir) = bs2 * ( s(cost,phi) - px(cost,phi) )
        !
      endif
      !
      if (l==-2) then ! sp2 hybrids
        !
        if (mr==1) ylm(ir) = bs3*s(cost,phi) - bs6*px(cost,phi) + bs2*py(cost,phi)
        if (mr==2) ylm(ir) = bs3*s(cost,phi) - bs6*px(cost,phi) - bs2*py(cost,phi)
        if (mr==3) ylm(ir) = bs3*s(cost,phi) + 2.d0*bs6*px(cost,phi)
        !
      endif
      !
      if (l==-3) then ! sp3 hybrids
        !
        if (mr==1) ylm(ir) = 0.5d0*( s(cost,phi) + px(cost,phi) + py(cost,phi) + pz_func(cost,phi) )
        if (mr==2) ylm(ir) = 0.5d0*( s(cost,phi) + px(cost,phi) - py(cost,phi) - pz_func(cost,phi) )
        if (mr==3) ylm(ir) = 0.5d0*( s(cost,phi) - px(cost,phi) + py(cost,phi) - pz_func(cost,phi) )
        if (mr==4) ylm(ir) = 0.5d0*( s(cost,phi) - px(cost,phi) - py(cost,phi) + pz_func(cost,phi) )
        !
      endif
      !
      if (l==-4) then ! sp3d hybrids
        !
        if (mr==1) ylm(ir) = bs3*s(cost,phi) - bs6*px(cost,phi) + bs2*py(cost,phi)
        if (mr==2) ylm(ir) = bs3*s(cost,phi) - bs6*px(cost,phi) - bs2*py(cost,phi)
        if (mr==3) ylm(ir) = bs3*s(cost,phi) + 2.d0*bs6*px(cost,phi)
        if (mr==4) ylm(ir) = bs2*pz_func(cost,phi) + bs2*dz2(cost,phi)
        if (mr==5) ylm(ir) =-bs2*pz_func(cost,phi) + bs2*dz2(cost,phi)
        !
      endif
      !
      if (l==-5) then ! sp3d2 hybrids
        !
        if (mr==1) ylm(ir) = bs6*s(cost,phi) - bs2*px(cost,phi) - bs12*dz2(cost,phi) + 0.5d0*dx2my2(cost,phi)
        if (mr==2) ylm(ir) = bs6*s(cost,phi) + bs2*px(cost,phi) - bs12*dz2(cost,phi) + 0.5d0*dx2my2(cost,phi)
        if (mr==3) ylm(ir) = bs6*s(cost,phi) - bs2*py(cost,phi) - bs12*dz2(cost,phi) - 0.5d0*dx2my2(cost,phi)
        if (mr==4) ylm(ir) = bs6*s(cost,phi) + bs2*py(cost,phi) - bs12*dz2(cost,phi) - 0.5d0*dx2my2(cost,phi)
        if (mr==5) ylm(ir) = bs6*s(cost,phi) - bs2*pz_func(cost,phi) + bs3*dz2(cost,phi)
        if (mr==6) ylm(ir) = bs6*s(cost,phi) + bs2*pz_func(cost,phi) + bs3*dz2(cost,phi)
        !
      endif
      !
    enddo !ir
    !
    return

  end subroutine ylm_wannier


  subroutine projection_expansion(l, mr, coef)
!----------------------------------------
!
!--------------------------------------------------------------------------
!     Outputs expansion coefficients for hybrid projections,
!     following wannier90 user guide table 3.2
!--------------------------------------------------------------------------

   use intw_useful_constants, only: ZERO

   implicit none

   !I/O variables
   integer, intent(in)  :: l, mr
   real(dp), intent(out) :: coef(:)
   !
   !local variables
   integer :: lm
   real(dp) :: fac1, fac2, fac3, fac4, fac5

   coef = ZERO
   !
   ! Check if l and m are within what's implemented in wannier90
   if (l>3 .or. l<-5) then
      !
      write(*,*)' projection_expansion: l out of range! '
      stop
      !
   end if
   !
   ! Compute coefficients
   !
   if (l>-1) then   ! single orbitals
      !
      ! Double-check if l and mr make sense
      if ( mr<1 .OR. mr>2*l+1 ) then
         !
         write(*,*)' ylm_wannier: m out of range! '
         stop
         !
      endif
      !
      lm = l**2 + mr
      coef(lm) = 1.d0
      !
   else ! hybrid orbitals
      !
      ! Double check if l and mr make sense
      if (mr < 1 .or. mr > abs(l)+1 ) then
         !
         write(*,*)' ylm_wannier: m out of range! '
         stop
         !
      end if
      !
      if (l==-1) then  !  sp hybrids
         !
         fac1 = 1.d0/sqrt(2.d0)
         !
         if (mr==1) then ! sp-1
            coef(1) = fac1
            coef(3) = fac1
         else if (mr==2) then ! sp-2
            coef(1) =  fac1
            coef(3) = -fac1
         end if
         !
      else if (l==-2) then  !  sp2 hybrids
         !
         fac1 = 1.d0/sqrt(3.d0)
         fac2 = 1.d0/sqrt(6.d0)
         fac3 = 1.d0/sqrt(2.d0)
         fac4 = 2.d0/sqrt(6.d0)
         !
         if (mr==1) then ! sp2-1
            coef(1) =  fac1
            coef(3) = -fac2
            coef(4) =  fac3
         else if (mr==2) then ! sp2-2
            coef(1) =  fac1
            coef(3) = -fac2
            coef(4) = -fac3
         else if (mr==3) then ! sp2-2
            coef(1) =  fac1
            coef(3) =  fac4
         end if
         !
      else if (l==-3) then  !  sp3 hybrids
         !
         fac1 = 1.d0/2.d0
         !
         if (mr==1) then ! sp3-1
            coef(1) =  fac1
            coef(2) =  fac1
            coef(3) =  fac1
            coef(4) =  fac1
         else if (mr==2) then ! sp3-2
            coef(1) =  fac1
            coef(2) = -fac1
            coef(3) =  fac1
            coef(4) = -fac1
         else if (mr==3) then ! sp3-3
            coef(1) =  fac1
            coef(2) = -fac1
            coef(3) = -fac1
            coef(4) =  fac1
         else if (mr==4) then ! sp3-4
            coef(1) =  fac1
            coef(2) =  fac1
            coef(3) = -fac1
            coef(4) = -fac1
         end if
         !
      else if (l==-4) then  !  sp3d hybrids
         !
         fac1 = 1.d0/sqrt(3.d0)
         fac2 = 1.d0/sqrt(6.d0)
         fac3 = 1.d0/sqrt(2.d0)
         fac4 = 2.d0/sqrt(6.d0)
         !
         if (mr==1) then ! sp3d-1
            coef(1) =  fac1
            coef(3) = -fac2
            coef(4) =  fac3
         else if (mr==2) then ! sp3d-2
            coef(1) =  fac1
            coef(3) = -fac2
            coef(4) = -fac3
         else if (mr==3) then ! sp3d-3
            coef(1) =  fac1
            coef(3) =  fac4
         else if (mr==4) then ! sp3d-4
            coef(2) =  fac3
            coef(5) =  fac3
         else if (mr==5) then ! sp3d-5
            coef(2) = -fac3
            coef(5) =  fac3
         end if
         !
      else if (l==-5) then  !  sp3d2 hybrids
         !
         fac1 = 1.d0/sqrt(6.d0)
         fac2 = 1.d0/sqrt(2.d0)
         fac3 = 1.d0/sqrt(12.d0)
         fac4 = 1.d0/2.d0
         fac5 = 1.d0/sqrt(3.d0)
         !
         if (mr==1) then ! sp3d2-1
            coef(1) =  fac1
            coef(3) = -fac2
            coef(5) = -fac3
            coef(8) =  fac4
         else if (mr==2) then ! sp3d2-2
            coef(1) =  fac1
            coef(3) =  fac2
            coef(5) = -fac3
            coef(8) =  fac4
         else if (mr==3) then ! sp3d2-3
            coef(1) =  fac1
            coef(3) = -fac2
            coef(5) = -fac3
            coef(8) = -fac4
         else if (mr==4) then ! sp3d2-4
            coef(1) =  fac1
            coef(3) =  fac2
            coef(5) = -fac3
            coef(8) = -fac4
         else if (mr==5) then ! sp3d2-5
            coef(1) =  fac1
            coef(2) = -fac2
            coef(5) =  fac5
         else if (mr==6) then ! sp3d2-6
            coef(1) =  fac1
            coef(2) =  fac2
            coef(5) =  fac5
         end if
         !
      end if
      !
   end if

   end subroutine



!===================================================
            ! MBR 18/04/24: currently unused!!!
    subroutine get_bvec_list (bvec)
! Reads the first k-vector and its neighbours shells from nnkp, and
! returns the list of the nnkp_nnkpts b-vectors used by W90
! in fractional coordinates

  implicit none
  real(dp) , intent(out) :: bvec(3,nnkp_nnkpts)
  integer   :: nn, ik, G(3)
  real(dp)  :: kpoint(3), kpoint_plus_b(3)

  kpoint = nnkp_kpoints(:,1)
  do nn = 1, nnkp_nnkpts
      G = nnkp_list_G(:,nn,1)
      ik = nnkp_list_ikpt_nn(nn,1)
      kpoint_plus_b =  nnkp_kpoints(:,ik) + real(G,dp)
      bvec(:,nn) = kpoint_plus_b - kpoint
  end do

  return
  end subroutine get_bvec_list
!===================================================



!----------------------------------------------------------------------------!
!
!
end module intw_intw2wannier
!
!
!----------------------------------------------------------------------------!
!
! This functions are originally distributed as part of the Quantum Espresso code and has
! been adapted to INTW:
!   Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
!   Distributed under the terms of the GNU General Public License.
!   See the LICENSE file in the original Quantum Espresso source for license details.
!   For the original source visit: https://www.quantum-espresso.org/
!
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
