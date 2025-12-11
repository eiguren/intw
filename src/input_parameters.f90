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

  !! display: public
  !!
  !! This module contains the definitions of all input parameters for INTW
  !! as well as subroutines for reading and checking the input.
  !!
  !! ### Details
  !!
  !! #### Input parameters
  !!
  !! ```{.txt}
  !! &input
  !!     outdir                = 'directory'
  !!     prefix                = 'prefix'
  !!     nk1                   = integer
  !!     nk2                   = integer
  !!     nk3                   = integer
  !!     tr_symmetry           = T or F (Default: T)
  !!     chemical_potential    = real (Default: 0.0 eV)
  !!     use_exclude_bands     = 'none', 'wannier' or 'custom'
  !!     include_bands_initial = integer
  !!     include_bands_final   = integer
  !! /
  !! &intw2W
  !!     intw2w_fullzone = T or F (Default: F)
  !!     intw2w_method   = 'CONVOLUTION' or 'FFT' (Default: 'CONVOLUTION')
  !!     compute_amn     = T or F (Default: T)
  !!     compute_mmn     = T or F (Default: T)
  !! /
  !! &ph
  !!     qlist           = 'file' (Default: 'qlist.txt')
  !!     read_for_dynmat = 'fc' or 'dynq' (Default: 'dynq')
  !!     fc_mat          = 'file'
  !!     nq1             = integer
  !!     nq2             = integer
  !!     nq3             = integer
  !!     nqirr           = integer
  !!     apply_asr       = T or F (Default: T)
  !! /
  !! &DOS
  !!     nk1_dos    = integer (Default: 20)
  !!     nk2_dos    = integer (Default: 20)
  !!     nk3_dos    = integer (Default: 20)
  !!     ne_dos     = integer (Default: 100)
  !!     eini_dos   = real (Default: -10.0 eV)
  !!     efin_dos   = real (Default: 10.0 eV)
  !!     esmear_dos = real (Default: 0.05 eV)
  !!     ktsmear    = real (Default: 0.01 eV)
  !! /
  !! &DOS_ph
  !!     nq1_dosph = integer (Default: 20)
  !!     nq2_dosph = integer (Default: 20)
  !!     nq3_dosph = integer (Default: 20)
  !!     nomega    = integer (Default: 100)
  !!     omega_ini = real (Default: 0.0 Ry)
  !!     omega_fin = real (Default: 0.005 Ry)
  !!     osmear_q  = real (Default: 0.000075 Ry)
  !!     omega_cut = real (Default: 0.0 Ry)
  !! /
  !! &elphon
  !!     ep_mat_file         = 'file' (Default: 'ep_mat.dat')
  !!     ep_bands            = 'intw' or 'custom' (Default: 'intw')
  !!     ep_bands_initial    = integer
  !!     ep_bands_final      = integer
  !!     ep_interp_method    = 'wannier' or 'dV_interpolate'
  !!     ep_interp_bands     = 'intw_bands' or 'ef_crossing' (Default: 'intw_bands')
  !!     nfs_sheets_initial  = integer
  !!     nfs_sheets_final    = integer
  !!     nscf_code           = 'QE' or 'SIESTA'
  !!     command_pw          = 'command'
  !!     command_pw2intw     = 'command'
  !!     file_pw             = 'file'
  !!     command_siesta2intw = 'command'
  !!     file_siesta2intw    = 'file'
  !! /
  !! K_PATH
  !!     nkpath nkspecial
  !!     label(1) kspecial_x(1) kspecial_y(1) kspecial_z(1)
  !!     label(2) kspecial_x(2) kspecial_y(2) kspecial_z(2)
  !!     ...
  !!     label(nkspecial) kspecial_x(nkspecial) kspecial_y(nkspecial) kspecial_z(nkspecial)
  !! Q_PATH
  !!     nqpath nqspecial
  !!     label(1) qspecial_x(1) qspecial_y(1) qspecial_z(1)
  !!     label(2) qspecial_x(2) qspecial_y(2) qspecial_z(2)
  !!     ...
  !!     label(nqspecial) qspecial_x(nqspecial) qspecial_y(nqspecial) qspecial_z(nqspecial)
  !! ```
  !!
  !! `K_PATH` and `Q_PATH` cards are used to generate the path along the Brilloin zone
  !! to plot electron band structure and phonon frequency dispersion, respectively.
  !!
  !! Variables without default values must be explicitly set in the input
  !! file; otherwise, an error will be raised.
  !!

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
  public :: ph, qlist, read_for_dynmat, fc_mat, &
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
  !! The directory where the DFT calculations are stored

  character(len=256) :: prefix = 'unassigned'
  !! The prefix of the DFT calculation

  integer :: nk1 = -1, nk2 = -1, nk3 = -1
  !! Monkhorst-Pack mesh indices for the coarse k-mesh

  real(dp) :: chemical_potential = 0.0_dp
  !! The value which determines the occupation factors (Units: eV)

  logical :: TR_symmetry = .true.
  !! If TR symmetry is present TR_symmetry = .true.

  character(256) :: use_exclude_bands = "unassigned"
  !! Three options to select the bands used in various utilities:
  !!
  !! - `use_exclude_bands='none'`: we don't exclude any band and therefore we use all bands from the DFT calculation
  !!   (Error if nnkp file is present)
  !    NOTE: Haritz 17/07/2025: Should we raise this error?
  !!
  !! - `use_exclude_bands='wannier'`: we exclude the bands indicated by Wannier90 in nnkp file
  !!   (Error if nnkp file is NOT present)
  !!
  !! - `use_exclude_bands='custom'`: use the subset of bands (`include_bands_initial:include_bands_final`)
  !!   from the bands considered by the DFT calculation and exclude the rest of bands
  !!   (Error if include_bands_initial, include_bands_final not present)

  integer :: include_bands_initial = 0, include_bands_final = 0
  !! The initial and final band indices for `use_exclude_bands='custom'`

  !----------------------------------------------------------------------------!
  ! &intw2W namelist variables
  !----------------------------------------------------------------------------!

  logical :: intw2W_fullzone = .false.
  !! If True, the code wil assume that a full zone DFT calculation
  !! has been performed and that wavefunctions for every k-point
  !! are available. This is mostly for testing and directly comparing
  !! the results of intw2W90 and pw2wannier

  character(256) :: intw2W_method = 'CONVOLUTION'
  !! Two options to select which method use to compute matrix elements:
  !!
  !! - `intw2W_method = 'CONVOLUTION'`
  !! - `intw2W_method = 'FFT'`

  logical :: compute_mmn = .true.
  !! If True, the code produces the $prefix.mmn and $prefix.eig files

  logical :: compute_amn = .true.
  !! If True, the code produces the $prefix.amn file

  !----------------------------------------------------------------------------!
  ! &ph namelist variables
  !----------------------------------------------------------------------------!

  character(256) :: qlist = 'qlist.txt'
  !! Name of the file containing the irreducible q-points list (relative to outdir)

  character(256) :: fc_mat = '--.fc'
  !! Name of the force constants matrix (relative to outdir)

  character(256) :: read_for_dynmat = 'dynq'
  !! Two options to choose:
  !!
  !! - `read_for_dynmat = 'fc'`: to read force constants from fc_mat file
  !! - `read_for_dynmat = 'dynq'`: to read .dyn files from the phonon calculations

  logical :: apply_asr = .true.
  !! Whether apply Acoustic Sum Rule or not

  integer :: nq1 = -1, nq2 = -1, nq3 = -1
  !! Monkhorst-Pack mesh indices for the coarse q-mesh

  integer :: nqirr = -1
  !! Number of irreducible q-points

  !----------------------------------------------------------------------------!
  ! &DOS namelist variables for Wannier interpolatied DOS plot
  !----------------------------------------------------------------------------!

  integer :: nk1_dos = 20, nk2_dos = 20, nk3_dos = 20
  !! Interpolation k-grid for DOS plot

  integer :: ne_dos = 100
  !! Number of energy points in DOS plot

  real(dp) :: eini_dos = -10.0_dp, efin_dos = 10.0_dp
  !! Energy range for electron DOS plot (Units: eV)

  real(dp) :: esmear_dos = 0.05_dp
  !! Smearing width for electron DOS plot (Units: eV)

  real(dp) :: ktsmear = 0.01_dp
  !! Smearing (kT) for Fermi level electron DOS calculation (Units: eV)

  !----------------------------------------------------------------------------!
  ! &DOS_ph namelist variables for intepolated phonon DOS plot
  ! and also used for the Eliashberg function plot
  !----------------------------------------------------------------------------!

  integer :: nq1_dosph = 20, nq2_dosph = 20, nq3_dosph = 20
  !! Interpolation q-grid for phonon DOS plot

  integer :: nomega = 100
  !! Number of energy points in DOS plot

  real(dp) :: omega_ini = 0.0_dp, omega_fin = 0.005_dp
  !! Energy range for phonon DOS plot (Units: Ry)

  real(dp) :: osmear_q = 0.000075
  !! Smearing width for phonon DOS plot (Units: Ry)

  real(dp) :: omega_cut = -1.0_dp
  !! Phonon frequency cut-off for removing w -> 0 peak in a2F (Units: Ry)

  !----------------------------------------------------------------------------!
  ! &elphon namelist variables for ep elements calculation
  !----------------------------------------------------------------------------!

  character(256) :: ep_mat_file = "ep_mat.dat"
  !! Name of the electron-phonon matrix elements file

  character(256) :: ep_bands = 'intw'
  !! Two options to select the subset of bands for computing the electron-phonon
  !! matrix elements:
  !!
  !! - `ep_bands = 'intw'`: compute matrix elements for all bands considered by INTW (See [[use_exclude_bands]] variable)
  !! - `ep_bands = 'custom'`: compute matrix elements for the custom subset of bands (ep_bands_initial:ep_bands_final)
  !!   from the bands considered by INTW (See [[use_exclude_bands]] variable)

  integer :: ep_bands_initial = 0
  !! The initial and final band indices for `ep_bands = 'custom'`
  integer :: ep_bands_final = 0
  !! The initial and final band indices for `ep_bands = 'custom'`

  ! For ep elements interpolation utilities:
  character(256) :: ep_interp_method = 'unassigned'
  !! Two options to choose the interpolation method:
  !!
  !! - `ep_interp_method = 'wannier'`
  !! - `ep_interp_method = 'dV_interpolate'`

  character(256) :: ep_interp_bands = 'intw_bands'
  !! Two options to select the subset of bands to interpolate electron-phonon matrix elements:
  !!
  !! - `ep_interp_bands = 'intw_bands'`: interpolate matrix elements for all bands considered by INTW (See use_exclude_bands variable)
  !! - `ep_interp_bands = 'ef_crossing'`: interpolate matrix elements for the substet of bands that cross the Fermi level (nfs_sheets_initial:nfs_sheets_final)
  !!   from the bands considered by INTW (See use_exclude_bands variable)

  integer :: nfs_sheets_initial=0
  !! The initial band index for `ep_interp_bands = 'ef_crossing'`
  integer :: nfs_sheets_final=0
  !! The final band index for `ep_interp_bands = 'ef_crossing'`

  character(len=256) :: nscf_code = 'unassigned'
  !! The DFT code used for running the nscf calculation (only needed if
  !! `ep_interp_method = 'dV_interpolate'` is chosen):
  !!
  !! - `nscf_code = 'QE'`
  !! - `nscf_code = 'SIESTA'`

  character(len=256) :: command_pw = 'unassigned'
  !! pw.x executable, with optional running options
  !! that go ahead of the executable like 'mpirun -np N', 'nice', etc.

  character(len=256) :: command_pw2intw = 'unassigned'
  !! pw2intw.x executable, with optional running options
  !! that go ahead of the executable like 'mpirun -np N', 'nice', etc.

  character(len=256) :: file_pw = 'unassigned'
  !! The input file for pw.x.

  character(len=256) :: command_siesta2intw = 'unassigned'
  !! siesta2intw.x executable, with optional running options
  !! that go ahead of the executable like 'mpirun -np N', 'nice', etc.

  character(len=256) :: file_siesta2intw = 'unassigned'
  !! The input file for siesta2intw.x.

  !----------------------------------------------------------------------------!
  ! K_PATH card variables
  !----------------------------------------------------------------------------!

  logical :: exist_kpath
  ! Whether K_PATH card for plotting electron bands exists or not

  integer :: nkpath
  !! Number of total k-points in the path

  integer :: nkspecial
  !! Number of special k-points in the path

  real(dp) , allocatable :: kspecial(:,:)
  !! List of special k-points in crystal coordinates

  !----------------------------------------------------------------------------!
  ! Q_PATH card variables
  !----------------------------------------------------------------------------!

  logical :: exist_qpath
  ! Whether Q_PATH card for plotting phonons bands exists or not

  integer :: nqpath
  !! Number of total q-points in the path

  integer :: nqspecial
  !! Number of special q-points in the path

  real(dp) , allocatable :: qspecial(:,:)
  !! List of special q-points in crystal coordinates

  !----------------------------------------------------------------------------!
  ! Define namelists
  !----------------------------------------------------------------------------!

  NAMELIST / input / outdir, prefix, nk1, nk2, nk3, &
                     TR_symmetry, chemical_potential, &
                     use_exclude_bands, include_bands_initial, include_bands_final

  NAMELIST / intw2W / intw2W_fullzone, intw2W_method, &
                      compute_mmn, compute_amn

  NAMELIST / ph / qlist, read_for_dynmat, fc_mat, ep_mat_file, &
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

    !! This subroutine reads all namelists from standard input,
    !! and checks their validity and consistency.

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

    return

20 format(A)
22 format(A,I2,A)

  end subroutine read_input


  subroutine read_cards ()

    ! Namelists have already been read from unit 5.
    ! Continue parsing unit until finding cards.

    !! This subroutine reads K_PATH and Q_PATH cards from standard input if present.
    !!
    !! IMPORTANT: some cards may be missing, but the present
    !! ones must be ordered as described below.
    !!
    !! ```{.txt}
    !! K_PATH
    !!     label(1) k_x(1) k_y(1) k_z(1)
    !!     label(2) k_x(2) k_y(2) k_z(2)
    !!     ...
    !!     label(nkspecial) k_x(nkspecial) k_y(nkspecial) k_z(nkspecial)
    !! Q_PATH
    !!     nqpath nqspecial
    !!     label(1) q_x(1) q_y(1) q_z(1)
    !!     label(2) q_x(2) q_y(2) q_z(2)
    !!     ...
    !!     label(nqspecial) q_x(nqspecial) q_y(nqspecial) q_z(nqspecial)
    !! ```
    !!
    !! where `nkpath` indicate the desired total number of k-points along the path,
    !! and `nkspecial` specifies the number of intermediate special points,
    !! which are given below in crystal coordinates.
    !!
    !! MBR 03/05/2024

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
