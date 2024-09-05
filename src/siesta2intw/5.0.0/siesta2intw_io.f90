module siesta2intw_io

  use precision, only: dp

  implicit none

  public :: stdout, siesta_output_file, intwdir
  public :: outdir, prefix, phdir, dvscfdir, &
            nk1, nk2, nk3, cutoff, nbnd_initial, nbnd_final, &
            phonons, nqirr, dynxml, use_sym

  public :: read_input, find_free_unit

  private

  integer :: stdout
  character(len=22), parameter :: siesta_output_file = "siesta2intw.siesta.out"

  character(len=256) :: intwdir

  ! input variables
  character(len=256) :: outdir, prefix, phdir, dvscfdir
  integer :: nk1, nk2, nk3
  real(kind=dp) :: cutoff
  integer :: nbnd_initial, nbnd_final
  logical:: phonons
  integer :: nqirr
  logical :: dynxml
  logical :: use_sym

  namelist / inputpp / outdir, prefix,  phdir, dvscfdir, &
                       nk1, nk2, nk3, cutoff, nbnd_initial, nbnd_final, &
                       phonons, nqirr, dynxml, use_sym

  ! Input variables not in pw2intw: nk1, nk2, nk3, symfile, nsym, cutoff, nbnd_initial, nbnd_final

contains

  subroutine read_input()
    !
    ! Read input file s2intw.in
    !

    implicit none

    integer :: ios, iunit
    integer :: strlen

    ! Set default values
    outdir = "./"
    prefix = "siesta"
    phdir ="ph/"
    dvscfdir = "dvscf/"
    nk1 = 1
    nk2 = 1
    nk3 = 1
    cutoff = 300.0_dp
    nbnd_initial = -1
    nbnd_final = -1
    phonons = .false.
    nqirr = 1
    dynxml = .false.
    use_sym = .true.

    iunit = find_free_unit()
    open(unit=iunit, file="s2intw.in", status="unknown", iostat=ios)
    if (ios/=0) stop "ERROR opening s2intw.in file"

    read(iunit, inputpp, iostat=ios)
    if (ios/=0) stop "ERROR reading s2intw.in file"

    close(unit=iunit, iostat=ios)
    if (ios/=0) stop "ERROR closing s2intw.in file"

    ! Some checks
    strlen = len_trim(outdir)
    if ( outdir(strlen:strlen+1) .ne. "/" ) outdir(strlen+1:strlen+2) = "/"
    strlen = len_trim(dvscfdir)
    if ( dvscfdir(strlen:strlen+1) .ne. "/" ) dvscfdir(strlen+1:strlen+2) = "/"
    strlen = len_trim(phdir)
    if ( phdir(strlen:strlen+1) .ne. "/" ) phdir(strlen+1:strlen+2) = "/"

  end subroutine read_input


  function find_free_unit()
    !
    !
    !

    use precision, only: dp

    implicit none

    integer :: find_free_unit, ios
    real(kind=dp) :: r
    logical :: opend

    call random_seed()
    do
       call random_number(r)
       find_free_unit = nint(r*1000)
       inquire ( unit = find_free_unit, opened = opend, iostat = ios)
       if ((.not.opend).and.(ios==0) ) exit
    enddo

  end function find_free_unit


end module siesta2intw_io