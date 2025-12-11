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
module intw_version

  !! display: none
  !!
  !! This module handles the version information of INTW.
  !!

  ! subroutines
  public :: print_intw_version

#ifdef VERSION_STRING
  character(25), parameter, private :: version = VERSION_STRING
  ! VERSION_STRING is passed as a pre-processing variable.
  ! This variable will be set by CMake from the version number specified in
  ! cmake project() (for example, v1.0.0), or, if possible, using 'git describe'
  ! (for example, v1.0.0-502-gff5fbda: commit ff5fbda, 502 commits ahead from parent tag v0.0.0).
#endif

contains

subroutine print_intw_version()

  implicit none

#ifdef VERSION_STRING
  write(*,'("|     INTW version: ",A25,"       |")') version
#else
    write(*,'("|     INTW version: ",A25,"       |")') "Unknown"
#endif

end subroutine print_intw_version

end module intw_version
