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
module siesta2ph_utils

  use precision, only: dp

  implicit none

  public :: pi, tpi, cmplx_0, cmplx_i, ZERO, &
            pmass, &
            eV2meV, meV2eV, Hartree2Ry, Ry2Hartree, &
            Ry2eV, eV2Ry, Ry2meV, meV2Ry, &
            Hartree2eV, eV2Hartree, Hartree2meV, meV2Hartree, &
            Bohr2Ang, Ang2Bohr, &
            au2THz, au2cm1

  public :: xsf_crystal_structure, xsf_datagrid_2d, xsf_datagrid_3d, &
            eqvect, word_in_string

  private

  ! useful constants
  real(kind=dp), parameter :: pi = dacos(-1.0_dp)
  real(kind=dp), parameter :: tpi = 2*dacos(-1.0_dp)
  complex(kind=dp), parameter :: cmplx_0 = (0.0_dp, 0.0_dp)
  complex(kind=dp), parameter :: cmplx_i = (0.0_dp, 1.0_dp)
  real(kind=dp), parameter :: ZERO = 0.0_dp
  !
  ! unit conversion constants
  ! mass
  real(kind=dp), parameter :: pmass = 1822.88848426_dp ! 1 Dalton in m_e
  ! energy
  real(kind=dp), parameter :: eV2meV = 1000.0_dp ! 1 eV in meV
  real(kind=dp), parameter :: meV2eV = 0.001_dp ! 1 meV in eV
  real(kind=dp), parameter :: Ry2eV = 13.605698066_dp ! 1 Ry in eV
  real(kind=dp), parameter :: eV2Ry = 1.0_dp/Ry2eV ! 1 eV in Ry
  real(kind=dp), parameter :: Ry2meV = 13605.698066_dp ! 1 Ry in meV
  real(kind=dp), parameter :: meV2Ry = 1.0_dp/Ry2meV ! 1 meV in Ry
  real(kind=dp), parameter :: Hartree2eV = 27.211396132_dp ! 1 Hartree in eV
  real(kind=dp), parameter :: eV2Hartree = 1.0_dp/Hartree2eV ! 1 eV in Hartree
  real(kind=dp), parameter :: Hartree2meV = 27211.396132_dp ! 1 Hartree in meV
  real(kind=dp), parameter :: meV2Hartree = 1.0_dp/Hartree2meV ! 1 meV in Hartree
  real(kind=dp), parameter :: Hartree2Ry = 2.0_dp ! 1 Hartree in Ry
  real(kind=dp), parameter :: Ry2Hartree = 0.5_dp ! 1 Ry in Hartree
  ! lenght
  real(kind=dp), parameter :: Bohr2Ang = 0.5291772_dp ! 1 Ang in Bohr
  real(kind=dp), parameter :: Ang2Bohr = 1.0_dp/Bohr2Ang ! 1 Ang in Bohr
  ! freq
  real(kind=dp), parameter :: au2THz = 6579.6839205_dp ! 1 freq. a.u. in THz
  real(kind=dp), parameter :: au2cm1 = Hartree2meV*8.065610_dp ! 1 freq. a.u. in cm^-1



contains

  subroutine xsf_crystal_structure(alat, at, nat, tau, atom_labels, ityp, iounit, filename)
    !
    ! Write the crystal structure in xsf format
    !
    implicit none
    !
    integer, intent(in) :: nat, ityp(nat)
    character(len=3), intent(in) :: atom_labels(*)
    real(kind=dp), intent(in) :: alat, tau(3, nat), at(3, 3)
    integer, intent(in) :: iounit
    character(len=*), intent(in) :: filename
    !
    integer :: in, ios


    open(unit=iounit, file=trim(filename), iostat=ios, status="replace", action="write")
    if ( ios /= 0 ) stop "ERROR: xsf_crystal_structure: couldn't open file"

    write(iounit, *) "CRYSTAL"
    ! lattice vectors in Angstroms
    write(iounit, *) "PRIMVEC"
    write(iounit, "(2(3F15.9/), 3f15.9)") at*alat*Bohr2Ang
    ! positions in Angstroms
    write(iounit, *) "PRIMCOORD"
    write(iounit, *) nat, 1
    do in=1, nat
      write(iounit, "(a3, 3x, 3f15.9)") atom_labels(ityp(in)), tau(:, in)*alat*Bohr2Ang
    enddo

    close(iounit)

  end subroutine xsf_crystal_structure


  subroutine xsf_datagrid_2d(rho, x0, e1, e2, nr1, nr2, at, alat, nat, tau, atom_labels, ityp, iounit, filename)
    !
    ! Write a 2D datagrid in xsf format
    !
    implicit none
    !
    integer, intent(in) :: nr1, nr2, nat, ityp(nat)
    character(len=3), intent(in) :: atom_labels(*)
    real(kind=dp), intent(in) :: x0(3), e1(3), e2(3), alat, tau(3, nat), at(3, 3), rho(nr1, nr2)
    integer, intent(in) :: iounit
    character(len=*), intent(in) :: filename
    !
    integer :: ix, iy, count, i, ind_x(10), ind_y(10)
    integer :: ios


    open(unit=iounit, file=trim(filename), iostat=ios, status="replace", action="write")
    if ( ios /= 0 ) stop "ERROR: xsf_datagrid_2d: couldn't open file"

    write(iounit, *) "CRYSTAL"
    ! lattice vectors in Angstroms
    write(iounit, *) "PRIMVEC"
    write(iounit, "(2(3F15.9/), 3f15.9)") at*alat*Bohr2Ang
    ! positions in Angstroms
    write(iounit, *) "PRIMCOORD"
    write(iounit, *) nat, 1
    do i=1, nat
      write(iounit, "(a3, 3x, 3f15.9)") atom_labels(ityp(i)), tau(:, i)*alat*Bohr2Ang
    enddo


    write(iounit,"(a)") "BEGIN_BLOCK_DATAGRID_2D"
    write(iounit,"(a)") "2D_PWSCF"
    write(iounit,"(a)") "DATAGRID_2D_UNKNOWN"

    ! real space mesh
    write(iounit, *) nr1, nr2
    ! origin in Angstoroms
    write(iounit, "(3f10.6)") x0*alat*Bohr2Ang
    ! spanning vectors in Angstoroms
    write(iounit, "(3f10.6)") e1*alat*Bohr2Ang
    write(iounit, "(3f10.6)") e2*alat*Bohr2Ang
    ! rho
    count=0
    do iy=1,nr2
      do ix=1,nr1
        !
        count = count + 1
        !
        ind_x(count) = ix
        ind_y(count) = iy
        !
        if (count == 6) then
          write(iounit,"(6e13.5)") (rho(ind_x(i), ind_y(i)), i=1, 6)
          count = 0
        endif
      enddo
    enddo

    write(iounit,"(6e13.5:)") (rho(ind_x(i), ind_y(i)), i=1, count)
    write(iounit,"(a)") "END_DATAGRID_2D"
    write(iounit,"(a)") "END_BLOCK_DATAGRID_2D"

    close(iounit)

  END subroutine xsf_datagrid_2d


  subroutine xsf_datagrid_3d(rho, nr1, nr2, nr3, at, alat, nat, tau, atom_labels, ityp, iounit, filename)
    !
    ! Write a 3D datagrid in xsf format
    !
    implicit none
    !
    integer, intent(in) :: nr1, nr2, nr3, nat, ityp(nat)
    character(len=3), intent(in) :: atom_labels(*)
    real(kind=dp), intent(in) :: alat, tau(3, nat), at(3, 3), rho(nr1, nr2, nr3)
    integer, intent(in) :: iounit
    character(len=*), intent(in) :: filename
    !
    integer :: i1, i2, i3, ix, iy, iz, count, i, ind_x(10), ind_y(10), ind_z(10)
    integer :: ios


    open(unit=iounit, file=trim(filename), iostat=ios, status="replace", action="write")
    if ( ios /= 0 ) stop "ERROR: xsf_datagrid_3d: couldn't open file"

    write(iounit, *) "CRYSTAL"
    ! lattice vectors in Angstroms
    write(iounit, *) "PRIMVEC"
    write(iounit, "(2(3F15.9/), 3f15.9)") at*alat*Bohr2Ang
    ! positions in Angstroms
    write(iounit, *) "PRIMCOORD"
    write(iounit, *) nat, 1
    do i=1, nat
      write(iounit, "(a3, 3x, 3f15.9)") atom_labels(ityp(i)), tau(:, i)*alat*Bohr2Ang
    enddo


    write(iounit, "(a)") "BEGIN_BLOCK_DATAGRID_3D"
    write(iounit, "(a)") "3D_PWSCF"
    write(iounit, "(a)") "DATAGRID_3D_UNKNOWN"
    ! real space mesh
    write(iounit, *) nr1+1, nr2+1, nr3+1
    ! origin in Angstoroms
    write(iounit, "(3f10.6)") 0.0d0, 0.0d0, 0.0d0
    ! lattice vectors in Angstoroms
    write(iounit, "(3f10.6)") at(:, 1)*alat*Bohr2Ang
    write(iounit, "(3f10.6)") at(:, 2)*alat*Bohr2Ang
    write(iounit, "(3f10.6)") at(:, 3)*alat*Bohr2Ang
    ! rho
    count=0
    do i3=0, nr3
      !
      iz = mod(i3, nr3) + 1
      !
      do i2=0, nr2
        !
        iy = mod(i2, nr2) + 1
        !
        do i1=0, nr1
          !
          ix = mod(i1, nr1) + 1
          !
          count = count + 1
          !
          ind_x(count) = ix
          ind_y(count) = iy
          ind_z(count) = iz
          !
          if (count==6) then
            write(iounit, "(6e13.5)") (rho(ind_x(i), ind_y(i), ind_z(i)), i=1, 6)
            count=0
          endif
          !
        enddo
        !
      enddo
      !
    enddo
    write(iounit, "(6e13.5:)") (rho(ind_x(i), ind_y(i), ind_z(i)), i=1, count)


    write(iounit, "(a)") "END_DATAGRID_3D"
    write(iounit, "(a)") "END_BLOCK_DATAGRID_3D"

    close(iounit)

  end subroutine xsf_datagrid_3d


  logical function eqvect(x, y)
    !
    ! TODO: Add description
    !
    use precision, only: dp
    !
    implicit none
    !
    real(kind=dp), dimension(3), intent(in) :: x, y ! input: input vectors
    real(kind=dp), parameter :: accep = 1.0d-4 ! acceptance parameter


    eqvect = abs( x(1)-y(1) - nint(x(1)-y(1)) ) < accep .and. &
             abs( x(2)-y(2) - nint(x(2)-y(2)) ) < accep .and. &
             abs( x(3)-y(3) - nint(x(3)-y(3)) ) < accep

  end function eqvect


  function word_in_string(word,string)
    !
    ! Check if a word is present in a string.
    !
    implicit none
    !
    character(*), intent(in) :: word, string
    logical :: word_in_string
    !
    character(len=len(word)) :: tmp_word
    character(len=len(string)) :: tmp_string
    integer :: len_word, len_string


    tmp_word = string_tolower(adjustl(word))
    len_word = len(trim(tmp_word))
    tmp_string = string_tolower(adjustl(string))
    len_string = len(trim(tmp_string))
    !
    if ( index(tmp_string,tmp_word) == 0 ) then
      word_in_string = .false.
    else
      word_in_string = .true.
    endif

  end function word_in_string


  function string_tolower(string)
    !
    ! Return the string in lowercase.
    !
    implicit none
    !
    character(*) :: string
    character(len=len(string)) :: string_tolower
    !
    integer :: i, k, length


    length = len(string)
    string_tolower = string
    do i = 1,len(string)
      k = iachar(string(i:i))
      if ( k >= iachar("A") .and. k <= iachar("Z") ) then
        k = k + iachar("a") - iachar("A")
        string_tolower(i:i) = achar(k)
      endif
    enddo

  end function string_tolower

end module siesta2ph_utils
