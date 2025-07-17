module pw2intw_io

  implicit none

  ! functions and subroutines
  public :: read_pw2intw_input

  ! variables
  public :: prefix, outdir, phonons, &
            phdir, nqirr, dynxml, fildyn, fildvscf, files4nscf

  ! Input variables
  character(len=256) :: prefix = " "
  character(len=256) :: outdir = "./"
  logical :: phonons = .false.
  character(len=256) :: phdir = "./"
  integer :: nqirr = 0
  logical :: dynxml = .false.
  character(len=256) :: fildyn = "matdyn"
  character(len=256) :: fildvscf = "dvscf"
  logical :: files4nscf = .false.

  namelist / inputpp / prefix, outdir, phonons, &
                       phdir, nqirr, dynxml, fildyn, fildvscf, files4nscf

contains

  subroutine read_pw2intw_input()
    !
    ! Read input file
    !
    use, intrinsic :: iso_fortran_env, only: stdin=>input_unit

    implicit none

    integer :: ios, strlen

    external :: errore


    read(stdin, inputpp, iostat=ios)
    if (ios /= 0) call errore( "pw2intw", "ERROR: pw2intw: error reading inputpp", ios )

    ! Some checks
    strlen = len_trim(outdir)
    if ( outdir(strlen:strlen+1) .ne. "/" ) outdir(strlen+1:strlen+2) = "/"
    strlen = len_trim(phdir)
    if ( phdir(strlen:strlen+1) .ne. "/" ) phdir(strlen+1:strlen+2) = "/"

  end subroutine read_pw2intw_input

end module pw2intw_io