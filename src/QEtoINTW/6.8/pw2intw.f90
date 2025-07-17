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
program pw2intw
  !-----------------------------------------------------------------------
  !
  use kinds
  use mp_global, only: mp_startup
  use io_files, only: prefix_QE => prefix
  use environment, only: environment_start, environment_end
  use basis, only: starting_wfc
  use control_flags, only: gamma_only

  use pw2intw_io
  use pw2intw_write

  implicit none

  external :: errore, read_file, wfcinit

  !
  ! Initialize MPI and environment
#if defined(__MPI)
  call mp_startup( )
#endif

  call environment_start( "pw2intw" )

  !
  ! Read input file
  call read_pw2intw_input()

  !
  ! Create intwdir
  intwdir = trim(outdir)//trim(prefix)//".save.intw/"
  call execute_command_line("mkdir -p "//intwdir)

  !
  ! Read QE data
  prefix_QE = prefix
  starting_wfc = "file"
  call read_file()
  call wfcinit()

  !
  ! Some checks
  if (gamma_only) call errore( "pw2intw", "gamma_only calculations are not implemented in INTW", 1 )

  !
  ! Write all data in intw format
  call write_pp_intw()

  call write_crystal_info()

  call scf_v_and_rho()

  call write_FFT_information ()

  call write_wfc()

  if (phonons) call write_phonon_info()

  if (files4nscf) call copy_nscf_files()

  !
  ! End job
  call environment_end( "pw2intw" )

end program pw2intw
