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
program nscf

  use kinds, only: dp
  use intw_version, only: print_intw_version
  use intw_input_parameters, only: intw2W_method, read_input
  use intw_reading, only: nspin, &
                          nbands, nGk_max, &
                          nr1, nr2, nr3, nbands, nkpoints_QE, &
                          get_gvec, &
                          read_parameters_data_file_xml, &
                          read_kpoints_data_file_xml, &
                          get_K_folder_data
  use intw_pseudo, only: read_all_pseudo
  use intw_utility, only: get_timing, print_date_time

  use intw_useful_constants, only: cmplx_0

  use intw_fft, only: generate_nl, allocate_fft, find_iG

  use intw_allwfcs, only: allocate_and_get_all_irreducible_wfc


  !================================================================================
  !       Declare the variables
  !================================================================================
  implicit none

  integer, allocatable     :: list_igk(:,:), ngk(:)

  complex(dp), allocatable :: wfc_k(:,:,:,:) ! nGk_max is defined in reading
  real(dp), allocatable    :: QE_eig_k(:,:)

  !local/aux variables
  integer                  :: ik
  integer                  :: ipol, jpol
  logical                  :: read_status
  character(256)           :: method

  complex(dp),allocatable :: wfc_k_r(:)
  real(kind=dp) :: time1, time2
  real(kind=dp), allocatable :: kpoints_QE(:,:)

  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  call get_timing(time1)
  write(*,20) '====================================================='
  write(*,20) '|                  program nscf                     |'
  write(*,20) '|         ---------------------------------         |'
  call print_intw_version()
  call print_date_time("Start of execution")
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
  !
  write(*,20) '|                    PPs are OK                     |'
  write(*,20) '|           ---------------------------------       |'
  !
  allocate (list_igk(nGk_max,nkpoints_QE))
  allocate (ngk(nkpoints_QE))
  !
  allocate (wfc_k(nGk_max,nbands,nspin,nkpoints_QE))
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
     write(*,"(a,i4,3f12.6,a)")'|',ik,kpoints_QE(:,ik),'           |'
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
  ! Finish
  !================================================================================

  call get_timing(time2)

  write(*,20) '|                      ALL DONE                     |'
  write(*,30) '|     Total time: ',time2-time1,' seconds            |'
  call print_date_time('End of execution  ')
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
    use intw_reading, only: gvec

    implicit none

    !I/O variables in (F fortran style). For example,
    !Input list_iG_1(:), instead of list_iG_1(nGk_max).

    integer,intent(in)      :: list_iG_1(:),ngk1,list_iG_2(:),ngk2
    complex(dp),intent(in)  :: wfc_1(:,:,:), wfc_2(:,:,:)

    ! In output, we have nbndxnbnd functions in (G), but G in the full ngm list

    complex(dp),intent(out) :: product_wfc (:,:,:,:,:)

    !local variables
    integer :: nbnd_l
    integer :: nGk_max_l
    integer :: nspin_l

    integer :: G2pG1(3)
    integer :: ibnd, jbnd, iG_1, iG_2, iG


    nGk_max_l  = size(wfc_1,1) !
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

end program nscf
