!----------------------------------------------------------------------------!
! intw project.
!----------------------------------------------------------------------------!
!
module intw_reading
  !
  use kinds, only: dp
  !
  implicit none
  !
  save
  !
  !
  ! variables
  public :: nbands, num_bands_intw, num_wann_intw, &
            num_exclude_bands_intw, band_excluded_intw, &
            s, can_use_TR, ftau, nsym, &
            atom_pfile, atom_labels, &
            at, bg, alat, tpiba, tpiba2, volume0, &
            nat, ntyp, ityp, tau, tau_cryst, amass, &
            nkpoints_QE, kpoints_QE, &
            nr1, nr2, nr3, ecutwfc, ecutrho, &
            ngm, gvec, nG_max, gamma_only, &
            nspin, lsda, noncolin, lspinorb, spinorb_mag

  !
  ! subroutines
  public :: read_parameters_data_file_xml, &
            read_kpoints_data_file_xml, &
            get_gvec, &
            get_K_folder_data, &
            deallocate_reading_variables, &
            set_num_bands, &
            scan_file_to
  !
  private
  !

  !----------------------------------------------------------------------------!
  ! Variables strongly inspired by QE variables
  ! (but not necessarily exactly equivalent!)
  !----------------------------------------------------------------------------!
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

  !!!!!!!!!!!!!!!!!!!!!!
    ! Symmetry variables
  !!!!!!!!!!!!!!!!!!!!!!

  integer, allocatable :: s(:, :, :)
  ! Symmetry matrices, in crystal coordinates,
  ! CAREFUL!!! since the symmetry matrices are expressed
  ! in the crystal coordinates, they may not *look* like
  ! rotation matrices: convert to cartesian coordinates
  ! and they'll have a sensible form.

  logical, allocatable :: can_use_TR(:)
  ! This array indicates if TR can (should) be used
  ! with a given point group operation. TR is optional
  ! if the ground state is not magnetically polarized, but
  ! can become mandatory; for a magnetic ground state,
  ! certain rotations can only be permitted in conjunction
  ! with time reversal symmetry (ie: the rotation flips B, and
  ! TR flips it back). This array will take note of
  ! which point group symmetries require TR.

  real(dp), allocatable :: ftau(:, :)
  ! Fractional translations, in crystal coordinates.

  integer :: nsym
  ! Number of crystal symmetries

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! system / configuration variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  character(len=80), allocatable :: atom_pfile(:)
  ! atom_pfile( j ) = name of pseudopotential file for
  ! the j-th atomic type (or species)

  character(len=3), allocatable :: atom_labels(:)
  ! atom_label( j ) = name of the j-th atomic type (or species)

  real(dp) :: at(3, 3)
  ! The a1 a2 and a3 lattice vectors column wise

  real(dp) :: bg(3, 3)
  ! The reciprocal lattice vectors column wise

  real(dp) :: alat
  ! The lattice parameter

  real(dp):: tpiba
  ! tpi/alat

  real(dp) :: tpiba2
  ! tpiba**2

  real(dp) :: volume0
  ! The volume of the unit cell

  integer :: nat
  ! The number of atoms in the cell

  integer :: ntyp
  ! The number types of atoms in the cell

  integer, allocatable :: ityp(:)
  ! ityp( i ) = type of the i-th atom

  real(dp), allocatable :: tau(:, :)
  ! tau( 1:3, i ) = position of the i-th atom (in cartesian, alat units)

  real(dp), allocatable :: tau_cryst(:, :)
  ! tau( 1:3, i ) = position of the i-th atom (in Crystal units)

  real(dp), allocatable :: amass(:)
  ! amass(1:ntyp) = atomic masses

  integer :: nkpoints_QE
  ! Number of k points used in the DFT calculation

  real(dp), allocatable :: kpoints_QE(:, :)
  ! The kpoints used in the DFT calculation

  integer :: nr1, nr2, nr3
  ! The FFT grid

  real(dp) :: ecutwfc
  ! Cutoff energy for wave functions (Hartree)

  real(dp) :: ecutrho
  ! Cutoff energy for charge density (Hartree)

  integer :: ngm
  ! The number of g vectors

  integer, allocatable :: gvec(:, :)
  ! The global g vectors (to which indices refer).

  integer :: nG_max
  ! The maximum number of G vectors for any k point

  logical :: gamma_only
  ! If Gamma is the only k-point and Gamma only tricks are used (half of the G-vectors are used)
  ! gamma_only = .true.

  integer :: nspin
  ! nspin = 1 for non-polarized calculations, nspin = 2 for spin-polarized ones.
  ! NOTE: Colinear spin-polarized calculations are transformed to a non-colinear
  !       spinor format by pw2intw or siesta2intw.

  logical :: lsda
  ! If a colinear calculation lsda = T

  logical :: noncolin
  ! If a noncolinear calculation noncolin = T

  logical :: lspinorb
  ! If spin-orbit non-collinear calculation lspinorb = T

  logical :: spinorb_mag ! TODO: I think that we should change the name of this variable, as at this moments only indicates if the induced potential is a 2x2 matrix or not.
  ! If spinorb_mag = T It is a situation with TR broken and spinor calculation.
  ! If spinorb_mag = F TR sym. is present!

contains

  subroutine read_parameters_data_file_xml()
    !------------------------------------------------------------------
    ! This subroutine reads in atomic positions and composition,
    ! as well as symmetry data from the prefix.save.intw folder.
    !------------------------------------------------------------------
    use intw_input_parameters, only: mesh_dir, prefix, TR_symmetry
    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: ainv
    use intw_useful_constants, only: tpi

    implicit none

    integer :: io_unit ! input/output file unit
    integer :: ii ! a looping variable

    character(256) :: datafile ! full path of the data-file.xml file in the .xml file
    ! note: mesh_dir and prefix are defined in the input_parameter module.

    integer :: max_nsym ! maximum possible value of nsym. Relevant
    ! variable for the magnetic case.

    integer :: i_sym

    integer, dimension(48) :: trev

    logical :: datafile_exists ! for checking if files in .save.intw folder exist

    integer, allocatable, dimension(:, :, :) :: ss
    real(dp), allocatable, dimension(:, :) :: fftau
    character(len=256) :: dummy

    integer :: i, j


    ! First, read the value of the parameter ngm
    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"gvectors.dat")

    inquire(file=datafile, exist=datafile_exists)

    if ( .not. datafile_exists ) then
      write(*, "(a)") '*****************************************************'
      write(*, "(a)") '*   ERROR                                           *'
      write(*, "(a)") '*   The file gvectors.dat cannot be found,          *'
      write(*, "(a)") '*   check if the folder .save.intw                  *'
      write(*, "(a)") '*   exists or if the prefix in intw input file      *'
      write(*, "(a)") '*   has been correctly chosen.                      *'
      write(*, "(a)") '*   STOPPING                                        *'
      write(*, "(a)") '*****************************************************'
      stop
    end if

    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="unformatted")
    read(unit=io_unit)ngm
    close(unit=io_unit)

    !read the data related to the crystal structure
    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"crystal.dat")
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="formatted")
    !ALAT
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) alat
    !AT
    read(unit=io_unit, fmt=*) dummy
    do i = 1, 3
      read(unit=io_unit, fmt=*) ( at(i, j), j = 1, 3 )
    enddo
    volume0 = alat**3*abs( at(1, 1)*(at(2, 2)*at(3, 3)-at(2, 3)*at(3, 2)) + &
                           at(1, 2)*(at(2, 3)*at(3, 1)-at(2, 1)*at(3, 3)) + &
                           at(1, 3)*(at(2, 1)*at(3, 2)-at(2, 2)*at(3, 1)) )

    !BG
    read(unit=io_unit, fmt=*) dummy
    do i = 1, 3
      read(unit=io_unit, fmt=*) ( bg(i, j), j = 1, 3 )
    enddo
    !FFT grid nr1 nr2 nr3
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nr1, nr2, nr3

    !ECUT_WFCS
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) ecutwfc
    !ECUT_POT
    read(unit=io_unit, fmt=*)dummy
    read(unit=io_unit, fmt=*) ecutrho

    tpiba = tpi/alat
    tpiba2 = tpiba*tpiba

    !LSDA
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) lsda
    !NONCOLIN
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) noncolin
    !LSO
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) lspinorb
    !DOMAG
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) spinorb_mag
    if (noncolin) then
       nspin = 2
    else
      nspin = 1
    end if

    !NAT
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nat

    !NTYP
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) ntyp
    if (ntyp<1) stop "ERROR: read_parameters_data_file_xml: ntyp < 1"

    allocate(atom_labels(ntyp))
    allocate(atom_pfile(ntyp))
    allocate(ityp(nat))
    allocate(tau(3, nat))
    allocate(tau_cryst(3, nat))
    allocate(amass(ntyp))
    !ATOM LABELS and PP files
    read(unit=io_unit, fmt=*) dummy
    do i = 1, ntyp
      read(unit=io_unit, fmt=*) atom_labels(i), amass(i), atom_pfile(i)
    end do

    !POSITIONS
    read(unit=io_unit, fmt=*) dummy
    do i = 1, nat
      read(unit=io_unit, fmt=*) dummy, ityp(i), tau(:, i)
      tau_cryst(:, i) = matmul(ainv(at), tau(:, i))
    end do

    !NSYM
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nsym
    allocate(s(3, 3, nsym))
    allocate(ftau(3, nsym))
    allocate(can_use_TR(nsym))
    !
    do ii = 1, nsym
      read(unit=io_unit, fmt=*) i_sym
      do i = 1, 3
        read(unit=io_unit, fmt=*) ( s(i, j, ii), j = 1, 3 )
      enddo
      read(unit=io_unit, fmt=*) ( ftau(j, ii), j = 1, 3 )
      read(unit=io_unit, fmt=*) trev(ii)
    end do !ii


    if ( (nspin == 1 .or. nspin == 2) .and. TR_symmetry ) then
      ! This is the easy case. Simply read in the symmetries

      do ii = 1, nsym
        if ( nspin == 1 ) then
          can_use_TR(ii) = TR_symmetry
        else if ( nspin == 2 ) then

          if ( trev(ii) == 1 ) then
            can_use_TR(ii) = .true.
          else
            can_use_TR(ii) = .false.
            if ( .not. spinorb_mag ) can_use_TR(ii) = .true.
          end if
        end if !nspin 1 or 2
      end do
    else if ( nspin == 2 .and. (.not. TR_symmetry) ) then
      ! This case is slightly more complicated. QE by default seems to
      ! assume time reversal symmetry is sometimes ok even in the presence
      ! of a B field. This is only true if B is colinear! I don't think it is
      ! true in the non-colinear case. Thus, rotations combined with TR must
      ! be suppressed!

      max_nsym = nsym

      nsym = 0

      do ii = 1, max_nsym
        if ( trev(ii) == 0 ) nsym = nsym + 1
      end do

      ! next, populate the symmetry array
      allocate(ss(3, 3, max_nsym))
      allocate(fftau(3, max_nsym))

      can_use_TR(:) = .false.
      ss(:, :, 1:max_nsym) = s(:, :, 1:max_nsym)
      fftau(1:3, 1:max_nsym) = ftau(1:3, 1:max_nsym)

      s = 0
      ftau = 0.0_dp

      i_sym = 0
      do ii = 1, max_nsym
        if ( trev(ii) == 0 ) then
          i_sym = i_sym + 1
          s(:, :, i_sym) = ss(:, :, ii)
          ftau(:, i_sym) = fftau(:, ii)
        end if
      end do

    end if !nspin

    !-KONTUZ
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nkpoints_QE

    !GAMMA_ONLY
    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) gamma_only
    if (gamma_only) stop "ERROR: intw does not support gamma_only calculations yet"

    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nG_max

    read(unit=io_unit, fmt=*) dummy
    read(unit=io_unit, fmt=*) nbands

    close(unit=io_unit)

  end subroutine read_parameters_data_file_xml


  subroutine read_kpoints_data_file_xml(kpoints_cryst)
    !------------------------------------------------------------------------
    ! This subroutine reads the prefix.save.intw/kpoints.dat file.
    ! The k-points are read in cartesian 2pi/a units: they will be converted
    ! to the more convenient crystal units
    !------------------------------------------------------------------------
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: ainv

    implicit none

    integer :: io_unit
    character(256) :: datafile

    real(dp) :: k(3), kpoints_cryst(3, nkpoints_QE)

    integer :: i


    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"kpoints.dat")
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="formatted")

    do i = 1, nkpoints_QE
      read(unit=io_unit, fmt=*) k(1:3)
      kpoints_cryst(:, i) = matmul(ainv(bg), k)
    end do

    close(unit=io_unit)

  end subroutine read_kpoints_data_file_xml


  subroutine get_gvec()
    !------------------------------------------------------------------------
    ! TODO: Add description
    !------------------------------------------------------------------------
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit

    implicit none

    integer :: io_unit
    character(256) :: datafile
    integer :: ig


    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"gvectors.dat")
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="unformatted")
    read(unit=io_unit) ngm

    allocate(gvec(3, ngm))

    do ig = 1, ngm
      read(unit=io_unit) gvec(1:3, ig)
    end do

    close(unit=io_unit)

  end subroutine get_gvec


  subroutine get_K_folder_data(ik, list_iG, wfc, eig, nG, altprefix)
    !------------------------------------------------------------------------
    ! For the kpoint labeled by ik, this subroutine reads all the
    ! wave functions for bands 1, .., nbands and stores them in the array
    ! wfc(nG_max_k, nbands). It reads the G vectors index array list_iG,
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
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: ZERO, cmplx_0

    implicit none

    !I/O variables
    integer, intent(in) :: ik
    integer, intent(out) :: list_iG(nG_max)
    real(dp), intent(out) :: eig(num_bands_intw)
    complex(dp), intent(out) :: wfc(nG_max, num_bands_intw, nspin)
    integer, intent(out) :: nG
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
      datafile = trim(trim(mesh_dir)//trim(altprefix)//".save.intw/"//trim(wfc_file))
    else
      datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//trim(wfc_file))
    end if
    io_unit = find_free_unit()
    open(unit=io_unit, file=datafile, status="unknown", action="read", form="unformatted")
    !
    ! Read data
    read(unit=io_unit) nG
    !
    read(unit=io_unit) list_iG(1:nG)
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
        read(unit=io_unit) ( wfc(1:nG, n_yes, is), is = 1, nspin )
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
    ! WARNING: This subroutine must always be called after "read_parameters_data_file_xml"
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
    if ( have_nnkp .and. trim(use_exclude_bands) .eq. 'all' ) then
            write(*, '(A)') '|  -  use_exclude_bands = all in input but          |'
            write(*, '(A)') '|       - .nnkp file found. Inconsistency!!         |'
            write(*, '(A)') ' Stopping '
            stop
    else if ( have_nnkp .and. trim(use_exclude_bands) .eq. 'custom' ) then
            write(*, '(A)') '|  -  use_exclude_bands = custom in input but       |'
            write(*, '(A)') '|       - .nnkp file found. Inconsistency!!         |'
            write(*, '(A)') ' Stopping '
            stop
    else if ( .not. have_nnkp .and. trim(use_exclude_bands) .eq. 'wannier' ) then
            write(*, '(A)') '|  -  use_exclude_bands = wannier in input but      |'
            write(*, '(A)') '|       - .nnkp file not found. Inconsistency!!     |'
            write(*, '(A)') ' Stopping '
            stop
    end if

    ! Set number of bands

    ! from .nnkp
    if ( trim(use_exclude_bands) .eq. 'wannier' ) then
      !
      write(*, '(A)') '|       - .nnkp file found                          |'
      write(*, '(A)') '|       - Setting number of bands from .nnkp        |'
      write(*, '(A)') '|           ---------------------------------       |'
      !
      nnkp_unit = find_free_unit()
      open(unit=nnkp_unit, file=nnkp_filename, status='old')
      !
      ! Number of wannier functions (after disentanglement)
      ! must be the same as number of projections
      if (noncolin) then
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
      ! Number of bands (before disentanglement)
      num_bands_intw = nbands - num_exclude_bands_intw
      !
    ! all bands from DFT
    else if ( trim(use_exclude_bands) .eq. 'all' ) then
      write(*, '(A)') '|       - Setting number of bands from calculation  |'
      write(*, '(A)') '|           ---------------------------------       |'
      !
      allocate(band_excluded_intw(nbands))
      band_excluded_intw(:) = .false.
      num_bands_intw = nbands
      num_wann_intw = nbands
      !
    ! custom bands from input
    else if ( trim(use_exclude_bands) .eq. 'custom' ) then
      write(*, '(A)') '|       - Setting custom number of bands from:      |'
      write(*, '(A)') '|', include_bands_initial, ' to ', include_bands_final
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
