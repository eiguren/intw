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
module intw_reading

  ! TODO: Add a description.

  use kinds, only: dp

  implicit none

  save

  ! variables
  public :: alat, tpiba, tpiba2, at, volume0, bg, &
            ntyp, atom_labels, amass, atom_pfile, &
            nat, tau, tau_cryst, ityp, &
            nsym, s, ftau, TR, &
            nkpoints_QE, kpoints_QE, &
            nbands, num_bands_intw, num_wann_intw, &
            num_exclude_bands_intw, band_excluded_intw, &
            nspin, lspin, lspinorb, lmag, &
            nr1, nr2, nr3, ecutwfc, ecutrho, &
            nG, gvec, nGk_max, gamma_only


  ! subroutines
  public :: read_parameters_data_file, &
            read_kpoints_data_file, &
            get_gvec, &
            get_K_folder_data, &
            deallocate_reading_variables, &
            set_num_bands, &
            scan_file_to

  private

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! unit cell variables
  ! system / configuration variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp) :: alat
  ! The lattice parameter TODOD: Units

  real(dp):: tpiba
  ! tpi/alat

  real(dp) :: tpiba2
  ! tpiba**2

  real(dp) :: at(3, 3)
  ! The a1 a2 and a3 lattice vectors column wise TODO: Units

  real(dp) :: volume0
  ! The volume of the unit cell TODO: units

  real(dp) :: bg(3, 3)
  ! The reciprocal lattice vectors column wise TODO: units

  integer :: ntyp
  ! The number types of atoms in the cell

  character(len=3), allocatable :: atom_labels(:)
  ! atom_label( j ) = name of the j-th atomic type (or species)

  real(dp), allocatable :: amass(:)
  ! amass(1:ntyp) = atomic masses TODOD: units

  character(len=80), allocatable :: atom_pfile(:)
  ! atom_pfile( j ) = name of pseudopotential file for
  ! the j-th atomic type (or species)

  integer :: nat
  ! The number of atoms in the cell

  real(dp), allocatable :: tau(:, :)
  ! tau( 1:3, i ) = position of the i-th atom (in cartesian, alat units)

  real(dp), allocatable :: tau_cryst(:, :)
  ! tau( 1:3, i ) = position of the i-th atom (in Crystal units)

  integer, allocatable :: ityp(:)
  ! ityp( i ) = type of the i-th atom

  !!!!!!!!!!!!!!!!!!!!!!
  ! Symmetry variables
  !!!!!!!!!!!!!!!!!!!!!!

  integer :: nsym
  ! Number of crystal symmetries

  integer, allocatable :: s(:, :, :)
  ! Symmetry matrices, in crystal coordinates,
  ! CAREFUL!!! since the symmetry matrices are expressed
  ! in the crystal coordinates, they may not *look* like
  ! rotation matrices: convert to cartesian coordinates
  ! and they'll have a sensible form.

  real(dp), allocatable :: ftau(:, :)
  ! Fractional translations, in crystal coordinates.

  logical, allocatable :: TR(:)
  ! This array indicates if TR can (should) be used
  ! with a given point group operation. TR is optional
  ! if the ground state is not magnetically polarized, but
  ! can become mandatory; for a magnetic ground state,
  ! certain rotations can only be permitted in conjunction
  ! with time reversal symmetry (ie: the rotation flips B, and
  ! TR flips it back). This array will take note of
  ! which point group symmetries require TR.

  !!!!!!!!!!!!!!!!!!!
  ! k-point variables
  !!!!!!!!!!!!!!!!!!!

  integer :: nkpoints_QE
  ! Number of k points used in the DFT calculation

  real(dp), allocatable :: kpoints_QE(:, :)
  ! The kpoints used in the DFT calculation TODO: crystal or Cartesian?

  !!!!!!!!!!!!!!!!!!!
  ! Bands variables
  !!!!!!!!!!!!!!!!!!!

  integer :: nbands
  ! The number of bands computed in the DFT calculation

  !! JLB: These are not read from the DFT calculation, but I think it's cleanest to put them here.
  !!      They are set in the subroutine "set_num_bands" in this module.
  !!      To be discussed.
  integer :: num_bands_intw
  ! Number of bands forming the subspace from which Wannier functions are to be obtained
  ! Might be different from nbands if exclude_bands is used

  integer :: num_wann_intw
  ! Number of Wannier functions, or equivalently bands/KS states to be interpolated
  ! Might be different from num_bands_intw if disentanglement is used

  integer :: num_exclude_bands_intw
  ! Number of bands to be excluded

  logical, allocatable :: band_excluded_intw(:)
  ! Array determining the bands to be excluded

  !!!!!!!!!!!!!!!!!!!
  ! Spin variables
  !!!!!!!!!!!!!!!!!!!

  integer :: nspin
  ! nspin = 1 for non-polarized calculations, nspin = 2 for spin-polarized ones.
  ! NOTE: Colinear spin-polarized calculations are transformed to a non-colinear
  !       spinor format by pw2intw or siesta2intw.

  logical :: lspin
  ! If a spin-polarized calculation (colinear or non-collinear) lspin = T

  logical :: lspinorb
  ! If spin-orbit non-collinear calculation lspinorb = T

  logical :: lmag
  ! If the calculation is magnetic (time-reversal broken) lmag = T
  ! NOTE: In QE, if the system has TR symmetry (no magnetization),
  !       only one spin component of the induced potential is calculated
  !       and saved to the file, even if the calculation is non-collinear.

  !!!!!!!!!!!!!!!!!!!
  ! Plane-waves variables
  !!!!!!!!!!!!!!!!!!!

  integer :: nr1, nr2, nr3
  ! The FFT grid

  real(dp) :: ecutwfc
  ! Cutoff energy for wave functions (Rydberg)

  real(dp) :: ecutrho
  ! Cutoff energy for charge density (Rydberg)

  integer :: nG
  ! The number of g vectors

  integer, allocatable :: gvec(:, :)
  ! The global g vectors (to which indices refer). TODO: crystal or Cartesian?

  integer :: nGk_max
  ! The maximum number of G vectors for any k point

  logical :: gamma_only
  ! If Gamma is the only k-point and Gamma only tricks are used (half of the G-vectors are used)
  ! gamma_only = .true.

contains

  subroutine read_parameters_data_file()
    !------------------------------------------------------------------
    ! This subroutine reads in atomic positions and composition, symmetries
    ! and all data from the prefix.save.intw/crystal.dat file.
    !------------------------------------------------------------------
    use intw_input_parameters, only: outdir, prefix, TR_symmetry
    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: ainv
    use intw_useful_constants, only: tpi

    implicit none

    integer :: io_unit, ierr ! input/output variables
    integer :: i, j, ii ! loop variables
    character(len=256) :: dummy

    character(256) :: datafile ! full path of the data-file.xml file in the .xml file
    ! note: outdir and prefix are defined in the input_parameter module.

    ! variables for the magnetic case.
    integer :: max_nsym ! maximum possible value of nsym. Relevant
    integer :: i_sym
    integer, allocatable, dimension(:, :, :) :: ss
    real(dp), allocatable, dimension(:, :) :: fftau
    integer, allocatable, dimension(:) :: trev


    ! Read the data related to the crystal structure
    datafile = trim(trim(outdir)//trim(prefix)//".save.intw/"//"crystal.dat")
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="formatted", iostat=ierr)
    if (ierr/=0) stop "ERROR: read_parameters_data_file: Error opening crystal.dat file!"

    !ALAT
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) alat

    !AT
    read(unit=io_unit, fmt=*) dummy
    do i = 1, 3
      read(unit=io_unit, fmt=*) ( at(i, j), j = 1, 3 )
    enddo

    !BG
    read(unit=io_unit, fmt=*) dummy
    do i = 1, 3
      read(unit=io_unit, fmt=*) ( bg(i, j), j = 1, 3 )
    enddo

    !NTYP
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) ntyp

    allocate(atom_labels(ntyp))
    allocate(atom_pfile(ntyp))
    allocate(amass(ntyp))

    !ATOM LABELS and PP files
    read(unit=io_unit, fmt=*) dummy
    do i = 1, ntyp
      read(unit=io_unit, fmt=*) atom_labels(i), amass(i), atom_pfile(i)
    end do

    !NAT
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nat

    allocate(ityp(nat))
    allocate(tau(3, nat))
    allocate(tau_cryst(3, nat))

    !POSITIONS
    read(unit=io_unit, fmt=*) dummy
    do i = 1, nat
      read(unit=io_unit, fmt=*) dummy, ityp(i), tau(:, i)
      tau_cryst(:, i) = matmul(ainv(at), tau(:, i))
    end do

    !NSYM
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nsym

    allocate(ss(3, 3, nsym))
    allocate(fftau(3, nsym))
    allocate(trev(nsym))

    !SYM
    do ii = 1, nsym
      read(unit=io_unit, fmt=*) dummy
      read(unit=io_unit, fmt=*) i_sym
      do i = 1, 3
        read(unit=io_unit, fmt=*) ( ss(i, j, ii), j = 1, 3 )
      enddo
      read(unit=io_unit, fmt=*) ( fftau(j, ii), j = 1, 3 )
      read(unit=io_unit, fmt=*) trev(ii)
    end do !ii

    !LSPIN
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) lspin

    !LSPINORB
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) lspinorb

    !LMAG
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) lmag

    if (lspin) then
       nspin = 2
    else
      nspin = 1
    end if

    !NKS
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nkpoints_QE

    !NBAND
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nbands

    !FFT GRID NR1 NR2 NR3
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nr1, nr2, nr3

    !ECUTWFC
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) ecutwfc

    !ECUTRHO
    read(unit=io_unit, fmt=*)dummy
    read(unit=io_unit, fmt=*) ecutrho

    !NG
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nG

    !NGK_MAX
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nGk_max

    close(unit=io_unit)

    ! Set some variables

    tpiba = tpi/alat
    tpiba2 = tpiba*tpiba

    volume0 = alat**3*abs( at(1, 1)*(at(2, 2)*at(3, 3)-at(2, 3)*at(3, 2)) + &
                           at(1, 2)*(at(2, 3)*at(3, 1)-at(2, 1)*at(3, 3)) + &
                           at(1, 3)*(at(2, 1)*at(3, 2)-at(2, 2)*at(3, 1)) )

    ! Remove symmetry opreations that are not wanted

    if ( lmag ) then

      ! Count valid symmetries
      max_nsym = nsym
      nsym = 0
      do ii = 1, max_nsym
        if ( trev(ii) == 0 .or. TR_symmetry ) nsym = nsym + 1
      end do

      allocate(s(3, 3, nsym))
      allocate(ftau(3, nsym))
      allocate(TR(nsym))

      ! Select only valid symmetries
      i_sym = 0
      do ii = 1, max_nsym
        if ( trev(ii) == 0 .or. TR_symmetry ) then
          i_sym = i_sym + 1
          s(:, :, i_sym) = ss(:, :, ii)
          ftau(:, i_sym) = fftau(:, ii)
          TR(i_sym) = trev(ii) == 1
        endif
      end do

    else

      ! This is the easy case. Simply read in the symmetries

      allocate(s(3, 3, nsym))
      allocate(ftau(3, nsym))
      allocate(TR(nsym))

      s = ss
      ftau = fftau
      TR = TR_symmetry ! If .not. lmag trev does not mean anything, it is allways 0

    end if

    deallocate(ss)
    deallocate(fftau)
    deallocate(trev)

  end subroutine read_parameters_data_file


  subroutine read_kpoints_data_file(kpoints_cryst)
    !------------------------------------------------------------------------
    ! This subroutine reads the prefix.save.intw/kpoints.dat file.
    ! The k-points are read in cartesian 2pi/a units: they will be converted
    ! to the more convenient crystal units
    !------------------------------------------------------------------------
    use intw_input_parameters, only: outdir, prefix
    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: ainv

    implicit none

    integer :: io_unit
    character(256) :: datafile

    real(dp) :: k(3), kpoints_cryst(3, nkpoints_QE)

    integer :: i


    datafile = trim(trim(outdir)//trim(prefix)//".save.intw/"//"kpoints.dat")
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="formatted")

    do i = 1, nkpoints_QE
      read(unit=io_unit, fmt=*) k(1:3)
      kpoints_cryst(:, i) = matmul(ainv(bg), k)
    end do

    close(unit=io_unit)

  end subroutine read_kpoints_data_file


  subroutine get_gvec()
    !------------------------------------------------------------------------
    ! TODO: Add description
    !------------------------------------------------------------------------
    use intw_input_parameters, only: outdir, prefix
    use intw_utility, only: find_free_unit

    implicit none

    integer :: io_unit
    character(256) :: datafile
    integer :: ig, nG_read


    datafile = trim(trim(outdir)//trim(prefix)//".save.intw/"//"gvectors.dat")
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="unformatted")

    read(unit=io_unit) nG_read

    if (nG/=nG_read) stop "ERROR: get_gvec: Wrong number of G vectors!"

    allocate(gvec(3, nG))

    do ig = 1, nG
      read(unit=io_unit) gvec(1:3, ig)
    end do

    close(unit=io_unit)

  end subroutine get_gvec


  subroutine get_K_folder_data(ik, list_iG, wfc, eig, nGk, altprefix)
    !------------------------------------------------------------------------
    ! For the kpoint labeled by ik, this subroutine reads all the
    ! wave functions for bands 1, .., nbands and stores them in the array
    ! wfc(nGk_max, nbands). It reads the G vectors index array list_iG,
    ! which refers to the global list of G vectors gvecs. It also reads
    ! the eigenvalues.
    !
    ! Do not be confused! G means "reciprocal lattice vector".
    !                     K means "point in the 1BZ"
    !
    ! NOTE: In fortran, "the first index varies fastest". I take this to mean
    ! that arrays are stored column-wise, namely, for a matrix
    !
    ! M = [ m_{11}   m_{12}   m_{13} ]
    !     [ m_{21}   m_{22}   m_{23} ]
    !     [ m_{31}   m_{32}   m_{33} ]
    !
    ! m_{21} is closer in memory to m_{11} than m_{12}, and in fact
    ! M is stored as
    ! M ~ [ m_{11} m_{21} m_{31} m_{12} m_{22} m_{32} m_{13} m_{23} m_{33}]
    ! Thus, it makes GOOD SENSE to put the G index first, as this is
    ! the index that will be used to perform inner products.
    !
    ! JLB: Excluded bands are taken care for at this point.
    !      Only relevant "num_bands_intw" wave functions are stored on output.
    !      This subroutine is only called in symmetries.f90,
    !      where only the wfc-s within "num_bands_intw" are rotated.
    !------------------------------------------------------------------------
    use intw_input_parameters, only: outdir, prefix
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: ZERO, cmplx_0

    implicit none

    !I/O variables
    integer, intent(in) :: ik
    integer, intent(out) :: list_iG(nGk_max)
    real(dp), intent(out) :: eig(num_bands_intw)
    complex(dp), intent(out) :: wfc(nGk_max, num_bands_intw, nspin)
    integer, intent(out) :: nGk
    character(256), optional, intent(in) :: altprefix

    !logical variables

    character(256) :: wfc_file, datafile
    integer :: io_unit, ibnd, is, n_yes
    real(dp) :: eig_all(nbands)

    !
    ! Initialize the arrays to zero (zero will be broadcasted)
    list_iG(:) = 0
    wfc(:, :, :) = cmplx_0
    eig(:) = ZERO
    !
    write(wfc_file, 100) ik
    100 format('wfc'I5.5'.dat')
    !
    if (present(altprefix)) then
      datafile = trim(trim(outdir)//trim(altprefix)//".save.intw/"//trim(wfc_file))
    else
      datafile = trim(trim(outdir)//trim(prefix)//".save.intw/"//trim(wfc_file))
    end if
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="unformatted")
    !
    ! Read data
    read(unit=io_unit) nGk
    !
    read(unit=io_unit) list_iG(1:nGk)
    !
    read(unit=io_unit) eig_all(1:nbands)
    !
    n_yes = 0
    do ibnd = 1, nbands
      !
      if (band_excluded_intw(ibnd)) then
        read(unit=io_unit)
      else
        n_yes = n_yes + 1
        read(unit=io_unit) ( wfc(1:nGk, n_yes, is), is = 1, nspin )
        eig(n_yes) = eig_all(ibnd)
      endif
        !
    enddo
    !
    close(io_unit)

  end subroutine get_K_folder_data


  subroutine deallocate_reading_variables()
    !------------------------------------------------------------------------
    ! TODO: Add description
    !------------------------------------------------------------------------

    deallocate(atom_labels)
    deallocate(atom_pfile)
    deallocate(amass)
    deallocate(ityp)
    deallocate(tau)
    deallocate(s)
    deallocate(ftau)

    ! JLB: Not sure whether this should go somewhere else
    if (allocated(band_excluded_intw)) deallocate(band_excluded_intw)

  end subroutine deallocate_reading_variables


  subroutine scan_file_to(nnkp_unit, keyword)
    !----------------------------------------------------------------------
    ! This subroutine reads a file all the way to the line
    !            begin $keyword
    ! This is useful when extracting parameters from the ascii file $seed.nnkp.
    !
    ! JLB 07/2023: Moved this subroutine here from intw2wannier.f90
    !----------------------------------------------------------------------

    implicit none

    integer, intent(in) :: nnkp_unit
    character(len=*), intent(in) :: keyword

    character(len=256) :: word1, word2
    logical            :: found, test
    integer            :: ios


    found = .false.
    !
    do
      !
      read(nnkp_unit, *, iostat=ios) word1, word2
      !
      test = (trim(word1).eq.'begin') .and. (trim(word2).eq.keyword)
      !
      if (test) exit
      !
      if (ios/=0) then
        !
        write(*, *) keyword, " data-block missing"
        stop
        !
      endif
      !
    enddo

  end subroutine scan_file_to


  subroutine set_num_bands()
    !----------------------------------------------------------------------
    ! This subroutine sets the sizes of the different band subspaces.
    ! Very important for all subsequent calculations,
    ! should be called at the beginning of all utilities.
    !
    ! For the moment it detects whether a .nnkp file exists,
    ! and in that case it reads them from there.
    ! If there's no such file, the number of bands is set to nbands
    ! (the number of bands in the DFT calculation).
    ! JLB: I think it would be useful to add the option to set them on input.
    !
    ! WARNING: This subroutine must always be called after "read_parameters_data_file"
    !          so that nbands is set.
    !
    ! MBR 20/05/2024: three options for exclude_bands
    !------------------------------------------------------------------------
    use intw_utility, only: find_free_unit
    use intw_input_parameters, only: prefix, use_exclude_bands, &
                                     include_bands_initial, include_bands_final

    implicit none

    character(len=256) :: nnkp_filename
    logical :: have_nnkp
    integer :: nnkp_unit, i, nn


    ! Check if .nnkp file exists
    nnkp_filename = trim(prefix)//trim('.nnkp')
    inquire(file=nnkp_filename, exist=have_nnkp)

    ! control consistency with use_exclude_bands flag
    if ( have_nnkp .and. trim(use_exclude_bands) .eq. 'none' ) then
            write(*, '(A)') '| - use_exclude_bands = none in input but           |'
            write(*, '(A)') '| - .nnkp file found. Inconsistency!!               |'
            write(*, '(A)') ' Stopping '
            stop
    else if ( have_nnkp .and. trim(use_exclude_bands) .eq. 'custom' ) then
            write(*, '(A)') '| - use_exclude_bands = custom in input but         |'
            write(*, '(A)') '| - .nnkp file found. Inconsistency!!               |'
            write(*, '(A)') ' Stopping '
            stop
    else if ( .not. have_nnkp .and. trim(use_exclude_bands) .eq. 'wannier' ) then
            write(*, '(A)') '| - use_exclude_bands = wannier in input but        |'
            write(*, '(A)') '| - .nnkp file not found. Inconsistency!!           |'
            write(*, '(A)') ' Stopping '
            stop
    end if

    ! Set number of bands

    ! from .nnkp
    if ( trim(use_exclude_bands) .eq. 'wannier' ) then
      !
      write(*, '(A)') '| - .nnkp file found                                |'
      write(*, '(A)') '| - Setting number of bands from .nnkp              |'
      write(*, '(A)') '|           ---------------------------------       |'
      !
      nnkp_unit = find_free_unit()
      open(unit=nnkp_unit, file=nnkp_filename, status='old')
      !
      ! Number of wannier functions (after disentanglement)
      ! must be the same as number of projections
      if (lspin) then
        call scan_file_to(nnkp_unit, 'spinor_projections')
        read(nnkp_unit, *) num_wann_intw
      else
        call scan_file_to(nnkp_unit, 'projections')
        read(nnkp_unit, *) num_wann_intw
      end if
      !
      ! Excluded bands, if any
      call scan_file_to(nnkp_unit, 'exclude_bands ')
      read(nnkp_unit, *) num_exclude_bands_intw
      allocate(band_excluded_intw(nbands))
      band_excluded_intw(:) = .false.
      do i = 1, num_exclude_bands_intw
         read(nnkp_unit, *) nn
         band_excluded_intw(nn) = .true.
      enddo
      !
      close(nnkp_unit)
      !
      ! Number of bands (before disentanglement)
      num_bands_intw = nbands - num_exclude_bands_intw
      !
    ! all bands from DFT
    else if ( trim(use_exclude_bands) .eq. 'none' ) then
      write(*, '(A)') '| - Setting number of bands from calculation        |'
      write(*, '(A)') '|           ---------------------------------       |'
      !
      allocate(band_excluded_intw(nbands))
      band_excluded_intw(:) = .false.
      num_bands_intw = nbands
      num_wann_intw = nbands
      !
    ! custom bands from input
    else if ( trim(use_exclude_bands) .eq. 'custom' ) then
      write(*, '(A)') '| - Setting custom number of bands:                 |'
      write(*, '(A9,I4,A4,I4,31X,A1)') '|   From ', include_bands_initial, ' to ', include_bands_final, '|'
      write(*, '(A)') '|           ---------------------------------       |'
      allocate(band_excluded_intw(nbands))
      num_bands_intw = include_bands_final - include_bands_initial + 1
      num_wann_intw = num_bands_intw
      num_exclude_bands_intw = nbands - num_bands_intw
      band_excluded_intw(:) = .true.
      do i = include_bands_initial, include_bands_final
          band_excluded_intw(i) = .false.
      end do
      !
    end if

  end subroutine

end module intw_reading
