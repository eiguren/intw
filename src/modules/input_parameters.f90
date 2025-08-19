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
module intw_input_parameters

  !-----------------------------------------------------------------------!
  ! This module contains the definitions of all input parameters for intw !
  ! as well as subroutines to read and test the input.                    !
  !-----------------------------------------------------------------------!

  use kinds, only: dp
  !
  implicit none
  !
  save
  !
  ! variables
  ! &input
  public :: input, outdir, prefix, nk1, nk2, nk3, TR_symmetry, chemical_potential, &
            use_exclude_bands, include_bands_initial, include_bands_final
  ! &intw2W
  public :: intw2W, intw2W_fullzone, intw2W_method, compute_mmn, compute_amn
  ! &ph
  public :: ph, ph_dir, qlist, read_for_dynmat, fc_mat, &
            nq1, nq2, nq3, nqirr, apply_asr
  ! &DOS
  public :: DOS, ne_dos, eini_dos, efin_dos, esmear_dos, ktsmear, nk1_dos, nk2_dos, nk3_dos
  ! &DOS_ph
  public :: DOS_ph, nq1_dosph, nq2_dosph, nq3_dosph, nomega, omega_ini, omega_fin, osmear_q, omega_cut
  ! &elphon
  public :: elphon, ep_mat_file, ep_bands, ep_bands_initial, ep_bands_final, ep_interp_method, &
            ep_interp_bands, nfs_sheets_initial, nfs_sheets_final, nscf_code, &
            command_pw, command_pw2intw, file_pw, &
            command_siesta2intw, file_siesta2intw
  ! K_PATH
  public :: exist_kpath, nkpath, nkspecial, kspecial
  ! Q_PATH
  public :: exist_qpath, nqpath, nqspecial, qspecial
  !
  ! subroutines
  public :: read_input, read_cards
  !
  private

  !----------------------------------------------------------------------------!
  ! &input namelist variables
  !----------------------------------------------------------------------------!

  character(len=256) :: outdir = 'unassigned'
  ! The directory where the DFT calculations are stored

  character(len=256) :: prefix = 'unassigned'
  ! The prefix of the DFT calculation

  integer :: nk1 = 0, nk2 = 0, nk3 = 0
  ! Monkhorst-Pack mesh indices for the coarse k-mesh

  real(dp) :: chemical_potential = 0.0_dp
  ! The value which determines the occupation factors, in eV

  logical :: TR_symmetry
  ! If TR symmetry is present TR_symmetry = .true.

  character(256) :: use_exclude_bands="unassigned"
  ! Three options to exclude bands from INTW
  ! This flags is to control the bands used in various utilities:
  !     use_exclude_bands="none": we don't exclude any band and therefore we use all bands from the DFT calculation (nbands)
  !          (Error if nnkp file is present)
  !          TODO: Haritz 17/07/2025: Should we raise this error?
  !     use_exclude_bands="wannier": we exclude the bands indicated by Wannier90 in nnkp file
  !          (Error if nnkp file is NOT present)
  !     use_exclude_bands="custom": will read a range of bands out of the whole nbands list and exclude the rest of bands
  !          (Error if include_bands_initial, include_bands_final not present)
  integer ::  include_bands_initial = 0, include_bands_final = 0

  !----------------------------------------------------------------------------!
  ! &intw2W namelist variables
  !----------------------------------------------------------------------------!

  logical :: intw2W_fullzone = .False.
  ! If True, the code wil assume that a full zone DFT calculation
  ! has been performed and that wavefunctions for every k-point
  ! are available. This is mostly for testing and directly comparing
  ! the results of intw2W90 and pw2wannier

  character(256) :: intw2W_method = 'CONVOLUTION'
  ! What method should be used to compute matrix elements:
  ! CONVOLUTION or FFT?

  logical :: compute_mmn = .true.
  ! If True, the code produces the $prefix.mmn and $prefix.eig files

  logical :: compute_amn = .true.
  ! If True, the code produces the $prefix.amn file

  !----------------------------------------------------------------------------!
  ! &ph namelist variables
  !----------------------------------------------------------------------------!

  character(256) :: ph_dir = './'
  ! Haritz 22/05/2024: At this moment this is only used as the directory where qlist and fc_mat are present
  ! I think that could be removed

  character(256) :: qlist = 'qlist.txt'
  ! Name of the file containing the irreducible q-points list

  character(256) :: fc_mat = '--.fc'
  ! Name of the force constants matrix

  character(256) :: read_for_dynmat = 'dynq'
  ! read_for_dynmat = 'fc' --> this will read force constants from fc_mat file
  ! read_for_dynmat = 'dynq' --> this will read .dyn files from the phonon calculations

  logical :: apply_asr = .true.
  ! Whether apply Acoustic Sum Rule or not

  integer :: nq1 = -1, nq2 = -1, nq3 = -1, nqirr = -1
  ! Monkhorst-Pack mesh indices for the coarse q-mesh and number of irreducible q-points

  !----------------------------------------------------------------------------!
  ! &DOS namelist variables for Wannier interpolatied DOS plot
  !----------------------------------------------------------------------------!

  integer :: nk1_dos = 20, nk2_dos = 20, nk3_dos = 20
  ! Interpolation k-grid for DOS plot

  integer :: ne_dos = 100
  ! Number of energy points in DOS plot

  real(dp) :: eini_dos = -10.0_dp, efin_dos = 10.0_dp, esmear_dos = 0.05_dp, ktsmear = 0.01_dp
  ! Energy range and smearing width, Fermi level smearing (kT). In eV!!!

  !----------------------------------------------------------------------------!
  ! &DOS_ph namelist variables for intepolated phonon DOS plot
  ! and also used for the Eliashberg function plot
  !----------------------------------------------------------------------------!

  integer :: nq1_dosph = 20, nq2_dosph = 20, nq3_dosph = 20
  ! Interpolation q-grid for phonon DOS plot

  integer :: nomega = 100
  ! Number of energy points in DOS plot

  real(dp) :: omega_ini = 0.0_dp, omega_fin = 0.005_dp, osmear_q = 0.000075
  ! Energy range and smearing width for phonon DOS plot (Units: Ry)

  real(dp) :: omega_cut = 0.0_dp
  ! Phonon energy cut-off for removing w -> 0 peak in a2F (Units: Ry)

  !----------------------------------------------------------------------------!
  ! &elphon namelist variables for ep elements calculation
  !----------------------------------------------------------------------------!

  character(256) :: ep_mat_file = "ep_mat.dat"
  ! Name of the electron-phonon matrix elements file

  character(256) :: ep_bands = 'intw'
  ! Will calculate ep elements for all bands (num_bands_intw, as by set_num_bands) or a
  ! selected substet:
  !   ep_bands = 'intw'   : all bands
  !   ep_bands = 'custom' : subset selected with ep_bands_initial and ep_bands_final

  integer ::  ep_bands_initial = 0, ep_bands_final = 0
  ! If ep_bands = 'custom', use the range ep_bands_initial:ep_bands_final subset
  ! from the num_bands_intw set

  ! For ep elements interpolation utilities:
  character(256) :: ep_interp_method = 'unassigned'
  !    ep_interp_method = 'wannier'
  !    ep_interp_method = 'dV_interpolate'

  character(256) :: ep_interp_bands = 'intw_bands'
  ! Will interpolate ep elements for all bands (num_bands_intw, as by set_num_bands) or a
  ! selected substet formed by those that cross the Fermi level (note, these indices
  ! must be provided by the user)
  !   ep_interp_bands = 'intw_bands'  : all bands
  !   ep_interp_bands = 'ef_crossing' : subset selected with nfs_sheets_initial and nfs_sheets_final

  integer :: nfs_sheets_initial=0, nfs_sheets_final=0
  ! If ef_crossing, use the range nfs_sheets_initial:nfs_sheets_final
  ! subset from the num_bands_intw set

  character(len=256) :: nscf_code = 'unassigned'
  ! The DFT code used for running the nscf calculation (only needed if
  ! ep_interp_method = 'dV_interpolate' is chosen).
  !   nscf_code = 'QE'
  !   nscf_code = 'SIESTA'

  character(len=256) :: command_pw = 'unassigned'
  character(len=256) :: command_pw2intw = 'unassigned'
  ! pw.x and pw2intw.x executables, with optional running options
  ! that go ahead of the executable like 'mpirun -np N', 'nice', etc.

  character(len=256) :: file_pw = 'unassigned'
  ! The input file for pw.x.

  character(len=256) :: command_siesta2intw = 'unassigned'
  ! siesta2intw.x executable, with optional running options
  ! that go ahead of the executable like 'mpirun -np N', 'nice', etc.

  character(len=256) :: file_siesta2intw = 'unassigned'
  ! The input file for siesta2intw.x.

  !----------------------------------------------------------------------------!
  ! K_PATH card variables
  !----------------------------------------------------------------------------!

  logical :: exist_kpath
  ! Whether K_PATH card for plotting electron bands exists or not

  integer :: nkpath, nkspecial
  ! Number of total k-points and number of special k-points in the path

  real(dp) , allocatable :: kspecial(:,:)
  ! Special k-points in crystal coordinates

  !----------------------------------------------------------------------------!
  ! Q_PATH card variables
  !----------------------------------------------------------------------------!

  logical :: exist_qpath
  ! Whether Q_PATH card for plotting phonons bands exists or not

  integer :: nqpath, nqspecial
  ! Number of total q-points and number of special q-points in the path

  real(dp) , allocatable :: qspecial(:,:)
  ! Special q-points in crystal coordinates

  !----------------------------------------------------------------------------!
  ! Define namelists
  !----------------------------------------------------------------------------!

  NAMELIST / input / outdir, prefix, nk1, nk2, nk3, &
                     TR_symmetry, chemical_potential, &
                     use_exclude_bands, include_bands_initial, include_bands_final

  NAMELIST / intw2W / intw2W_fullzone, intw2W_method, &
                      compute_mmn, compute_amn

  NAMELIST / ph / ph_dir, qlist, read_for_dynmat, fc_mat, ep_mat_file, &
                  nq1, nq2, nq3, nqirr, apply_asr

  NAMELIST / DOS / nk1_dos, nk2_dos, nk3_dos, ne_dos, &
                   eini_dos, efin_dos, esmear_dos, ktsmear

  NAMELIST / DOS_ph / nq1_dosph, nq2_dosph, nq3_dosph, nomega, &
                      omega_ini, omega_fin, osmear_q, omega_cut

  NAMELIST / elphon / ep_bands, ep_bands_initial, ep_bands_final, &
                      ep_interp_method, ep_interp_bands, &
                      nfs_sheets_initial, nfs_sheets_final, nscf_code, &
                      command_pw, command_pw2intw, file_pw, &
                      command_siesta2intw, file_siesta2intw

  ! ----------------------------------------------------------------------
  !  END namelist
  ! ----------------------------------------------------------------------

contains

  subroutine read_input(read_status)
    !------------------------------------------------------------------
    ! This subroutine reads outdir and prefix from standard input
    !------------------------------------------------------------------
    use intw_useful_constants, only: eps_10

    implicit none

    integer :: ios, strlen
    logical :: read_status

    read_status = .false.


    write(*,20) '| - Reading standard input file...                  |'
    write(*,20) '|         namelist             ios                  |'
    write(*,20) '|         --------             ----                 |'


    read( 5, input, iostat = ios)
    write(*,22) '|           &input             ',ios,'                   |'

    read( 5, intw2W, iostat = ios)
    write(*,22) '|           &intw2W            ',ios,'                   |'

    read( 5, ph, iostat = ios)
    write(*,22) '|           &ph                ',ios,'                   |'

    read( 5, DOS, iostat = ios)
    write(*,22) '|           &DOS               ',ios,'                   |'

    read( 5, DOS_ph, iostat = ios)
    write(*,22) '|           &DOS_ph            ',ios,'                   |'

    read( 5, elphon, iostat = ios)
    write(*,22) '|           &elphon            ',ios,'                   |'

    write(*,20) '====================================================='


    !  Test the various read parameters

    if ( outdir == 'unassigned' ) then
      read_status = .true.
      write(*,*) 'MISSING outdir!'
    end if

    if ( prefix == 'unassigned' ) then
      read_status = .true.
      write(*,*) 'MISSING prefix!'
    end if

    if ( nk1 == 0 ) then
      read_status = .true.
      write(*,*) 'MISSING nk1!'
    end if

    if ( nk2 == 0 ) then
      read_status = .true.
      write(*,*) 'MISSING nk2!'
    end if

    if ( nk3 == 0 ) then
      read_status = .true.
      write(*,*) 'MISSING nk3!'
    end if

    if ( trim(use_exclude_bands) == 'unassigned' ) then
      read_status = .true.
      write(*,*) 'MISSING use_exclude_bands'
    else if (      trim(use_exclude_bands) /= 'none' &
             .and. trim(use_exclude_bands) /= 'custom' &
             .and. trim(use_exclude_bands) /= 'wannier' ) then
      read_status = .true.
      write(*,*) 'Error: Wrong value for use_exclude_bands!'
    end if

    if ( trim(use_exclude_bands) == 'custom' .and. &
        (include_bands_final*include_bands_initial .le. 0 .or. &
          include_bands_final .lt. include_bands_initial ) ) then
      read_status = .true.
      write(*,*) 'Error: Invalid use_exclude_bands custom selection!'
    end if

    if (      trim(ep_bands) /= 'intw' &
        .and. trim(ep_bands) /= 'custom' ) then
      read_status = .true.
      write(*,*) 'Error: Wrong value for ep_bands!'
    end if

    if ( trim(ep_bands) == 'custom' &
         .and. &
         (      ep_bands_initial*ep_bands_final <= 0 &
           .or. ep_bands_final < ep_bands_initial ) ) then
      read_status = .true.
      write(*,*) 'Error: Invalid ep_bands custom selection!'
    end if

    if (      trim(ep_interp_bands) /= 'intw_bands' &
        .and. trim(ep_interp_bands) /= 'ef_crossing' ) then
      read_status = .true.
      write(*,*) 'Error: Wrong value for ep_interp_bands!'
    end if

    if ( trim(ep_interp_bands) == 'ef_crossing' &
         .and. &
         ( nfs_sheets_initial*nfs_sheets_final < 0 &
           .or. nfs_sheets_final < nfs_sheets_initial ) ) then
      read_status = .true.
      write(*,*) 'Error: Invalid ep_interp_bands custom selection!'
    end if

    if (trim(ep_interp_method) /= 'unassigned') then
      if (      trim(ep_interp_method) /= 'wannier' &
              .and. trim(ep_interp_method) /= 'dV_interpolate' ) then
        read_status = .true.
        write(*,*) 'Error: Wrong value for ep_interp_method!'
      endif
    end if

    if ( trim(ep_interp_method) == 'dV_interpolate' ) then
      !
      if ( trim(nscf_code) == 'unassigned' ) then
        read_status = .true.
        write(*,*) 'Error: dV_interpolate method chosen, but nscf_code not specified!'
      else if ( trim(nscf_code) == 'QE' ) then
        if ( trim(command_pw) == 'unassigned' ) then
          read_status = .true.
          write(*,*) 'Error: dV_interpolate method chosen with QE, but command_pw not specified!'
        endif
        if ( trim(command_pw2intw) == 'unassigned' ) then
          read_status = .true.
          write(*,*) 'Error: dV_interpolate method chosen with QE, but command_pw2intw not specified!'
        endif
        if ( trim(file_pw) == 'unassigned' ) then
          read_status = .true.
          write(*,*) 'Error: dV_interpolate method chosen with QE, but file_pw not specified!'
        endif
      else if ( trim(nscf_code) == 'SIESTA' ) then
        if ( trim(command_siesta2intw) == 'unassigned' ) then
          read_status = .true.
          write(*,*) 'Error: dV_interpolate method chosen with SIESTA, but command_siesta2intw not specified!'
        endif
        if ( trim(file_siesta2intw) == 'unassigned' ) then
          read_status = .true.
          write(*,*) 'Error: dV_interpolate method chosen with SIESTA, but file_siesta2intw not specified!'
        endif
      else
        read_status = .true.
        write(*,*) 'Error: dV_interpolate method chosen with unknown nscf_code!'
      endif
      !
    end if


    if (      trim(read_for_dynmat) /= 'dynq' &
        .and. trim(read_for_dynmat) /= 'fc' ) then
      read_status = .true.
      write(*,*) 'Error: Wrong value for read_for_dynmat!'
    end if


    if ( read_status ) then
      write(*,*) "The input should be of the form:"
      write(*,*) "&input"
      write(*,*) "             outdir                = 'directory'"
      write(*,*) "             prefix                = 'prefix'"
      write(*,*) "             nk1                   = integer"
      write(*,*) "             nk2                   = integer"
      write(*,*) "             nk3                   = integer"
      write(*,*) "             TR_symmetry           = T or F"
      write(*,*) "             chemical_potential    = real"
      write(*,*) "             use_exclude_bands     = 'none', 'wannier' or 'custom'"
      write(*,*) "             include_bands_initial = integer"
      write(*,*) "             include_bands_final   = integer"
      write(*,*) "/"
      write(*,*) "&intw2W"
      write(*,*) "             intw2W_fullzone = T or F"
      write(*,*) "             intw2W_method   = CONVOLUTION or FFT"
      write(*,*) "             compute_amn     = T or F"
      write(*,*) "             compute_mmn     = T or F"
      write(*,*) "/"
      write(*,*) "&ph"
      write(*,*) "             ph_dir          = 'directory'"
      write(*,*) "             qlist           = 'file'"
      write(*,*) "             read_for_dynmat = 'fc' or 'dynq' (D)"
      write(*,*) "             fc_mat          = 'file'"
      write(*,*) "             nq1             = integer"
      write(*,*) "             nq2             = integer"
      write(*,*) "             nq3             = integer"
      write(*,*) "             nqirr           = integer"
      write(*,*) "             apply_asr       = T(D) or F"
      write(*,*) "/"
      write(*,*) "&DOS"
      write(*,*) "             nk1_dos    = integer"
      write(*,*) "             nk2_dos    = integer"
      write(*,*) "             nk3_dos    = integer"
      write(*,*) "             ne_dos     = integer"
      write(*,*) "             eini_dos   = real"
      write(*,*) "             efin_dos   = real"
      write(*,*) "             esmear_dos = real"
      write(*,*) "             kTsmear    = real"
      write(*,*) "/"
      write(*,*) "&DOS_ph"
      write(*,*) "             nq1_dosph = integer"
      write(*,*) "             nq2_dosph = integer"
      write(*,*) "             nq3_dosph = integer"
      write(*,*) "             nomega    = integer"
      write(*,*) "             omega_ini = real"
      write(*,*) "             omega_fin = real"
      write(*,*) "             osmear_q  = real"
      write(*,*) "/"
      write(*,*) "&elphon"
      write(*,*) "             ep_mat_file         = 'file'"
      write(*,*) "             ep_bands            = 'intw' or 'custom'"
      write(*,*) "             ep_bands_initial    = integer"
      write(*,*) "             ep_bands_final      = integer"
      write(*,*) "             ep_interp_method    = 'wannier' or 'dV_interpolate'"
      write(*,*) "             ep_interp_bands     = 'intw_bands' or 'ef_crossing' "
      write(*,*) "             nfs_sheets_initial  = integer"
      write(*,*) "             nfs_sheets_final    = integer"
      write(*,*) "             nscf_code           = 'QE' or 'SIESTA'"
      write(*,*) "             command_pw          = 'string'"
      write(*,*) "             command_pw2intw     = 'string'"
      write(*,*) "             file_pw             = 'file'"
      write(*,*) "             command_siesta2intw = 'string'"
      write(*,*) "             file_siesta2intw    = 'file'"
      write(*,*) "/"
    end if

    strlen = len_trim(outdir)
    if ( outdir(strlen:strlen+1) .ne. "/" ) outdir(strlen+1:strlen+2) = "/"
    strlen = len_trim(ph_dir)
    if ( ph_dir(strlen:strlen+1) .ne. "/" ) ph_dir(strlen+1:strlen+2) = "/"

    return

20 format(A)
22 format(A,I2,A)

  end subroutine read_input


  subroutine read_cards ()
    ! MBR 03/05/2024

    ! Namelists have already been read from unit 5.
    ! Continue parsing unit until finding cards.

    ! IMPORTANT: some cards may be missing, but the present
    ! ones must be ordered as described below.

    ! Format:
    !   K_PATH
    !      nkpath  nkspecial
    !      k1 k2 k3
    !      k1 k2 k3
    !      ....
    !   Q_PATH
    !      nqpath  nqspecial
    !      q1 q2 q3
    !      q1 q2 q3
    !      ....
    ! where nkpath are the # points sought and nkspecial are
    ! the number ofintermediate special points, given below in
    ! FRACTIONAL coordinates.

    implicit none

    character(80) :: cardname
    character(1) :: klabel, qlabel
    integer :: i


    exist_kpath = .false.
    exist_qpath = .false.

    write(*,20) '| - Reading standard input file...                  |'
    write(*,20) '|   card:                                           |'

    ! parse until K_PATH, Q_PATH
    do
      !
      read(5,*,end=100) cardname
      !
      if ( trim(cardname) == 'K_PATH') then
        !
        write(*,20) '|             K_PATH                                |'
        exist_kpath = .true.
        read(5,*) nkpath, nkspecial
        !
        if (nkspecial .lt. 2) then
          write(*,20) 'More than 2 k-points are needed to set up the path. Stopping.'
          stop
        end if
        !
        allocate (kspecial(3,nkspecial))
        do i=1,nkspecial
          read(5,*) klabel, kspecial(:,i)
        end do
        !
      else if ( trim(cardname) == 'Q_PATH') then
        !
        write(*,20) '|             Q_PATH                                |'
        exist_qpath = .true.
        read(5,*) nqpath, nqspecial
        !
        if (nqspecial .lt. 2) then
          write(*,20) 'More than 2 q-points are needed to set up the path. Stopping.'
          stop
        end if
        !
        allocate (qspecial(3,nqspecial))
        do i=1,nqspecial
          read(5,*) qlabel, qspecial(:,i)
        end do
        !
      end if
      !
    end do

100 write(*,20) '|   EOF reached                                     |'
    if ( .not. exist_kpath .and. .not. exist_qpath) then
      write(*,20) '|   no cards found                                  |'
    end if

    write(*,20) '====================================================='

    20 format(A)

   end subroutine read_cards

END MODULE intw_input_parameters
