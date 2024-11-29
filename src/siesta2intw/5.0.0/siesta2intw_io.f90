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
module siesta2intw_io

  use precision, only: dp

  implicit none

  public :: stdout, siesta_output_file, intwdir
  public :: outdir, prefix, phdir, dvscfdir, &
            nk1, nk2, nk3, cutoff, nbnd_initial, nbnd_final, &
            phonons, nqirr, dynxml, use_sym, &
            kmesh, nkpoints, kpoints

  public :: read_input, find_free_unit

  private

  integer :: stdout
  character(len=22), parameter :: siesta_output_file = "siesta2intw.siesta.out"

  character(len=256) :: intwdir

  ! input variables
  character(len=256) :: outdir, prefix, phdir, dvscfdir
  integer :: nk1, nk2, nk3 ! k mesh, if they are negative, a list of k points given in the input is used instead of a regular k mesh.
  real(kind=dp) :: cutoff
  integer :: nbnd_initial, nbnd_final
  logical:: phonons
  integer :: nqirr
  logical :: dynxml
  logical :: use_sym
  !
  logical :: kmesh ! .true. if nk1, nk2, nk3 are specified
  integer :: nkpoints
  real(kind=dp), allocatable, dimension(:,:) :: kpoints

  namelist / inputpp / outdir, prefix, phdir, dvscfdir, &
                       nk1, nk2, nk3, cutoff, nbnd_initial, nbnd_final, &
                       phonons, nqirr, dynxml, use_sym

  ! Input variables not in pw2intw: nk1, nk2, nk3, symfile, nsym, cutoff, nbnd_initial, nbnd_final

contains

  subroutine read_input()
    !
    ! Read input file s2intw.in
    !

    use, intrinsic :: iso_fortran_env, only: iostat_end ! end of file

    implicit none

    integer :: ios, iunit
    integer :: strlen
    character(80) :: cardname
    integer :: i
    logical :: exist_kpoints


    ! Set default values
    outdir = "./"
    prefix = "siesta"
    phdir ="ph/"
    dvscfdir = "dvscf/"
    nk1 = 0
    nk2 = 0
    nk3 = 0
    cutoff = 300.0_dp
    nbnd_initial = -1
    nbnd_final = -1
    phonons = .false.
    nqirr = 1
    dynxml = .false.
    use_sym = .true.
    !
    nkpoints = 0
    kpoints = 0.0_dp

    iunit = find_free_unit()
    open(unit=iunit, file="s2intw.in", status="unknown", iostat=ios)
    if (ios/=0) stop "ERROR opening s2intw.in file"

    read(iunit, inputpp, iostat=ios)
    if (ios/=0) stop "ERROR reading s2intw.in file"

    ! Some checks
    strlen = len_trim(outdir)
    if ( outdir(strlen:strlen+1) .ne. "/" ) outdir(strlen+1:strlen+2) = "/"
    strlen = len_trim(dvscfdir)
    if ( dvscfdir(strlen:strlen+1) .ne. "/" ) dvscfdir(strlen+1:strlen+2) = "/"
    strlen = len_trim(phdir)
    if ( phdir(strlen:strlen+1) .ne. "/" ) phdir(strlen+1:strlen+2) = "/"

    ! Read KPOINTS
    ! Format:
    !   KPOINTS
    !      nkpoints
    !      k1 k2 k3
    !      k1 k2 k3
    !      ....
    ! where k points are given in crystal coordinates.

    exist_kpoints = .false.
    ! Read until KPOINTS is found
    do
      !
      read(iunit, *, iostat=ios) cardname
      if ( ios == iostat_end ) exit
      !
      if ( trim(cardname) == 'KPOINTS') then
        !
        exist_kpoints = .true.
        read(iunit,*) nkpoints
        !
        allocate (kpoints(3,nkpoints))
        do i=1,nkpoints
          read(iunit,*) kpoints(:,i)
        end do
        !
      end if
      !
    end do
    !
    ! Some checks
    if ( nk1 <= 0 .and. nk2 <= 0 .and. nk3 <= 0 ) then
      kmesh = .false.
      if ( .not. exist_kpoints ) then
        stop "ERROR: You must specify (nk1, nk2, nk3) or KPOINTS."
      end if
    else
      kmesh = .true.
      if ( exist_kpoints ) then
        write(stdout, *) "WARNINRG: Both (nk1, nk2, nk3) and KPOINTS are specified, (nk1, nk2, nk3) will be used."
      endif
    endif
    !
    close(unit=iunit, iostat=ios)
    if (ios/=0) stop "ERROR closing s2intw.in file"

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