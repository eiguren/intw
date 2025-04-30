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
module siesta2ph_io

  use precision, only: dp

  implicit none

  public :: stdin, stdout
  public :: outdir, prefix, phdir, v0dir, dvscfdir, kpath_file, &
            nr1, nr2, nr3, dx, lpm, disp_along_cart, dv_precision, &
            irreducible_q, xsf, full_xsf, verbose, use_sym

  public :: read_input, write_fdf, find_free_unit

  private

  ! i/o variables
  integer, parameter :: stdin = 5
  integer, parameter :: stdout = 6

  ! input variables
  character(len=256) :: outdir, prefix, phdir, v0dir, dvscfdir, kpath_file
  integer :: nr1, nr2, nr3
  real(kind=dp) :: dx
  logical :: lpm
  logical :: disp_along_cart
  logical :: use_sym
  logical :: irreducible_q
  logical :: xsf, full_xsf
  character(len=2) :: dv_precision
  logical :: verbose
  !
  namelist / input / outdir, prefix, phdir, v0dir, dvscfdir, kpath_file, &
                     nr1, nr2, nr3, dx, lpm, disp_along_cart, dv_precision, &
                     use_sym, irreducible_q, xsf, full_xsf, verbose


contains

  subroutine read_input()
    !
    ! Read the input namelist from the stdin file.
    !
    implicit none
    !
    integer :: ios


    ! Set default values for variables in namelist
    outdir = "./"
    prefix = ""
    phdir = "ph/"
    v0dir = "v0/"
    dvscfdir = "dvscf/"
    kpath_file = ""
    nr1 = 1
    nr2 = 1
    nr3 = 1
    dx = 0.01_dp
    lpm = .true.
    disp_along_cart = .false.
    use_sym = .true.
    irreducible_q = .true.
    xsf = .false.
    full_xsf = .false.
    dv_precision = "dp"
    verbose = .false.
    !
    ! Reading the namelist inputpp
    !
    ios = 0
    read(stdin, input, iostat=ios)
    !
    ! Check of namelist variables
    !
    if (ios /= 0) stop "ERROR reading input namelist"
    !
    if (dv_precision /= "dp" .and. dv_precision /= "sp") stop "ERROR: Invalid value for precision"
    !
    if (outdir(len(trim(outdir)):len(trim(outdir))+1) /= "/") outdir(len(trim(outdir))+1:len(trim(outdir))+2) = "/"
    if (phdir(len(trim(phdir)):len(trim(phdir))+1) /= "/") phdir(len(trim(phdir))+1:len(trim(phdir))+2) = "/"
    if (v0dir(len(trim(v0dir)):len(trim(v0dir))+1) /= "/") v0dir(len(trim(v0dir))+1:len(trim(v0dir))+2) = "/"
    if (dvscfdir(len(trim(dvscfdir)):len(trim(dvscfdir))+1) /= "/") dvscfdir(len(trim(dvscfdir))+1:len(trim(dvscfdir))+2) = "/"

  end subroutine read_input


  subroutine write_fdf(nat, tau, ityp, cell, reference_fdffile, fdffilename)
    !
    ! Write a fdf file taking an existing fdf as a reference and replacing
    ! the unit cell by nat, tau, ityp and cell.
    !
    use, intrinsic :: iso_fortran_env, only: iostat_end ! end of file
    !
    use siesta2ph_utils, only: word_in_string
    !
    implicit none
    !
    integer, intent(in) :: nat
    real(kind=dp), dimension(:,:), intent(in) :: tau
    integer, dimension(:), intent(in) :: ityp
    real(kind=dp), dimension(3,3), intent(in) :: cell
    character(len=*), intent(in) :: reference_fdffile, fdffilename
    !
    character(len=256) :: line
    integer :: iounit_reference, iounit
    integer :: ia, id
    integer :: ios


    ! Check input array dimensions
    if (size(tau, dim=1) /= 3) stop "write_fdf: tau must have dimensions (3,nat)"
    if (size(tau, dim=2) /= nat) stop "write_fdf: tau must have dimensions (3,nat)"
    if (size(ityp, dim=1) /= nat) stop "write_fdf: ityp must have dimensions (nat)"
    if (size(cell, dim=1) /= 3 .or. size(cell, dim=2) /= 3) stop "write_fdf: cell must have dimensions (3,3)"
    !
    iounit_reference = find_free_unit()
    open(unit=iounit_reference, file=trim(outdir)//trim(reference_fdffile), status="old", action="read", iostat=ios)
    if ( ios /= 0 ) stop "write_fdf: Error opening file outdir/reference_fdffile"
    !
    iounit = find_free_unit()
    open(unit=iounit, file=trim(outdir)//trim(fdffilename), status="replace", action="write", iostat=ios)
    if ( ios /= 0 ) stop "write_fdf: Error opening file outdir/filename"
    !
    do
      !
      ! Read each line of the input file
      read(unit=iounit_reference, fmt="(a)", iostat=ios) line
      !
      ! Exit the do loop if it is the last line
      if ( ios == iostat_end ) exit
      !
      ! Check what is specified in the line
      if (word_in_string("NumberOfAtoms",line)) then
        !
        write(iounit,"(a14,i4)") "NumberOfAtoms ", nat
        !
      elseif (word_in_string("LatticeVectors",line) .and. word_in_string("block",line)) then
        !
        do
          read(unit=iounit_reference, fmt="(a)", iostat=ios) line
          if (word_in_string("LatticeVectors",line) .and. word_in_string("endblock",line)) exit
        enddo
        !
        write(iounit,"(a)") "%block LatticeVectors"
        do id=1,3
          write(iounit,"(3f15.9,i4)") cell(:,id)
        enddo
        write(iounit,"(a)") "%endblock LatticeVectors"
        !
      elseif (word_in_string("AtomicCoordinatesFormat",line)) then
        !
        ! Specify the AtomicCoordinatesFormat as ScaledCartesian (atomic coordinates will be written in this format)
        write(iounit,"(a)") "AtomicCoordinatesFormat    ScaledCartesian"
        !
      elseif (word_in_string("AtomicCoordinatesAndAtomicSpecies",line) .and. word_in_string("block",line)) then
        !
        do
          read(unit=iounit_reference, fmt="(a)", iostat=ios) line
          if (word_in_string("AtomicCoordinatesAndAtomicSpecies",line) .and. word_in_string("endblock",line)) exit
        enddo
        !
        ! Write the coordinates of the atoms
        write(iounit,"(a)") "%block AtomicCoordinatesAndAtomicSpecies"
        do ia=1,nat
          write(iounit,"(3f15.9,i4)") tau(:,ia), ityp(ia)
        enddo
        write(iounit,"(a)") "%endblock AtomicCoordinatesAndAtomicSpecies"
        !
      else
        !
        write(iounit,"(a)") trim(line)
        !
      endif
      !
    end do
    !
    close(iounit)
    close(iounit_reference)

  end subroutine write_fdf


  function find_free_unit()
    !
    ! This function finds a free input/output unit.
    !
    implicit none
    !
    integer :: find_free_unit
    integer :: io_unit
    logical :: opnd


    do io_unit = 100, 1000
      !
      inquire( unit = io_unit, opened = opnd )
      if ( .not. opnd ) then
        find_free_unit = io_unit
        return
      end if
      !
    end do
    !
    stop "ERROR: find_free_unit: NO free units available!"

  end function find_free_unit

end module siesta2ph_io
