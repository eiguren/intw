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
program siesta2ph
  !
  ! This program reads an input fdf file and calculates the irreducible
  ! displacements needed to calculate the full force constant matrix.
  ! Finally, it writtes the new input fdf files to run the calculations.
  !
  use precision, only: dp
  use m_timestamp, only: timestamp
#ifdef MPI
  use mpi_siesta, only: MPI_Comm_World
#endif
  use parallel, only: Node, Nodes
  !
  use siesta2ph_io, only: outdir, v0dir, phdir, prefix, nr1, nr2, nr3, lpm, verbose, stdout
  use siesta2ph_system, only: nat, at, alat, tau, ityp
  use siesta2ph_symmetry, only: irred_atm, irred_disp
  !
  use siesta2ph_io, only: read_input, write_fdf
  use siesta2ph_system, only: read_unit_cell_data, print_unit_cell_data
  use siesta2ph_symmetry, only: compute_symmetry, print_symmetry_data, &
                                find_irreducible_atoms, find_irreducible_displacements, &
                                find_site_symmetry
  !
  implicit none
  !
#ifdef MPI
  ! MPI variables
  logical :: initialized
  integer :: MPIerror
#endif
  !
  ! supercell variables
  real(kind=dp), allocatable, dimension(:,:) :: tau_sc
  real(kind=dp), dimension(3,3) :: at_sc
  integer, allocatable, dimension(:) :: ityp_sc
  integer :: nat_sc
  !
  ! Monkhorst-Pack kgrid for the unit cell and supercell
  integer :: kgrid_uc(3,3), kgrid_sc(3,3)
  real(kind=dp) :: kgrid_disp(3)
  !
  ! Real space mesh for the unit cell and supercell
  integer, dimension(3) :: rmesh_uc, rmesh_sc
  !
  !
  ! Initialise environment
  !
#ifdef MPI
  call MPI_Initialized( initialized, MPIerror )
  if (.not.initialized) then
    call MPI_Init( MPIerror )
  endif ! (.not.initialized)
  call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )

  if (Nodes>1) then
    call die("siesta2ph: ERROR: This code is not prepared to run in multiple nodes")
  endif
#endif
  !
  Node = 1 ! to turn off output from siesta subroutines
  !
  call timestamp("Start of run")
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Read input file
  call read_input()
  !
  ! Read the unit cell data from the fdf file
  call read_unit_cell_data()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute symmetries
  call compute_symmetry()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Find the irreducible atoms and atom mapping
  call find_irreducible_atoms()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Find site symmetries for each irreducible atom
  call find_site_symmetry()
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Print all data
  if (verbose) then
    call print_unit_cell_data()
    call print_symmetry_data()
  endif
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Find the irreducible displacements
  call find_irreducible_displacements()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Create the input file of the supercell
  call create_supercell()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Create the input file for each irreducible displacements
  call create_displacements()
  !
  !
#ifdef MPI
  call MPI_Finalize( MPIerror )
#endif

contains

  subroutine create_supercell()
    !
    ! Creates the lattice vectors and atomic positions of the supercell
    ! and writes the fdf file
    !
    implicit none
    !
    real(kind=dp), dimension(3) :: R
    integer :: ia, ja, ir1, ir2, ir3


    write(stdout,*) "- Creating the supercell file..."
    !
    nat_sc = nat*nr1*nr2*nr3
    !
    allocate(tau_sc(3,nat_sc), ityp_sc(nat_sc))
    ja = 0
    do ir3=1,nr3
      do ir2=1,nr2
        do ir1=1,nr1
          !
          R = (ir1-1)*at(:,1) + (ir2-1)*at(:,2) + (ir3-1)*at(:,3)
          !
          do ia=1,nat
            ja = ja + 1
            tau_sc(:, ja) = tau(:, ia) + R
            ityp_sc(ja) = ityp(ia)
          enddo
          !
        enddo
      enddo
    enddo
    !
    at_sc(:,1) = at(:,1)*nr1
    at_sc(:,2) = at(:,2)*nr2
    at_sc(:,3) = at(:,3)*nr3
    !
    call execute_command_line("mkdir -p "//trim(outdir)//trim(v0dir))
    call write_fdf(nat_sc, tau_sc, ityp_sc, at_sc, prefix, trim(v0dir)//"supercell-"//trim(prefix))
    !
    ! Check kgrid
    call check_kgrid(kgrid_uc, kgrid_disp)
    kgrid_sc(:,1) = kgrid_uc(:,1)/nr1
    kgrid_sc(:,2) = kgrid_uc(:,2)/nr2
    kgrid_sc(:,3) = kgrid_uc(:,3)/nr3
    !
    if (any(kgrid_sc /= 0)) then
      call modify_kgrid(kgrid_sc, kgrid_disp, trim(v0dir)//"supercell-"//trim(prefix))
    endif
    !
    ! Check mesh sizes
    call check_mesh_sizes(rmesh_uc)
    rmesh_sc = rmesh_uc*(/nr1, nr2, nr3/)
    !
    if (any(rmesh_sc > 0)) then
      call modify_MeshSizes(rmesh_sc, trim(v0dir)//"supercell-"//trim(prefix))
    endif

  end subroutine create_supercell


  subroutine create_displacements()
    !
    ! Displaces the irreducible atoms of the supercell along the irreducible
    ! directions and writes the fdf file for each displacement
    !
    use siesta2ph_symmetry, only: set_displacements_direction
    !
    implicit none
    !
    character(len=256) :: dispp_folder
    character(len=256) :: dispn_folder
    integer :: ia, id, iirred
    real(kind=dp), dimension(3,3) :: disp


    write(stdout,*) "- Creating the ireducible displacement files..."
    !
    call set_displacements_direction(disp)
    !
    iirred = 0
    do ia=1,nat
      !
      if ( .not. irred_atm(ia) ) cycle
      !
      do id=1,3
        !
        if ( .not. irred_disp(id,ia) ) cycle
        !
        iirred = iirred + 1
        write(stdout,"(a,i3,2(a5,i3))") "irreducible displacement", iirred, ": ia=", ia, ", id=", id
        !
        ! Positive displacement
        if (lpm) then
          write(dispp_folder,"(a5,i4.4,a1)") "disp-", 2*iirred-1, "/"
        else
          write(dispp_folder,"(a5,i4.4,a1)") "disp-", iirred, "/"
        endif
        !
        tau_sc(:,ia) = tau(:,ia) + disp(:,id)/alat ! add the positive displacement to the atom
        !
        call execute_command_line("mkdir -p "//trim(outdir)//trim(phdir)//trim(dispp_folder))
        call write_fdf(nat_sc, tau_sc, ityp_sc, at_sc, prefix, trim(phdir)//trim(dispp_folder)//"supercell-"//trim(prefix))
        !
        if (any(kgrid_sc /= 0)) then
          call modify_kgrid(kgrid_sc, kgrid_disp, trim(phdir)//trim(dispp_folder)//"supercell-"//trim(prefix))
        endif
        !
        if (any(rmesh_sc > 0)) then
          call modify_MeshSizes(rmesh_sc, trim(phdir)//trim(dispp_folder)//"supercell-"//trim(prefix))
        endif
        !
        if (lpm) then
          !
          ! Negative displacement
          write(dispn_folder,"(a5,i4.4,a1)") "disp-", 2*iirred, "/"
          !
          tau_sc(:,ia) = tau(:,ia) - disp(:,id)/alat ! add the negative displacement to the atom
          !
          call execute_command_line("mkdir -p "//trim(outdir)//trim(phdir)//trim(dispn_folder))
          call write_fdf(nat_sc, tau_sc, ityp_sc, at_sc, prefix, trim(phdir)//trim(dispn_folder)//"supercell-"//trim(prefix))
          !
          if (any(kgrid_sc /= 0)) then
            call modify_kgrid(kgrid_sc, kgrid_disp, trim(phdir)//trim(dispn_folder)//"supercell-"//trim(prefix))
          endif
          !
          if (any(rmesh_sc > 0)) then
            call modify_MeshSizes(rmesh_sc, trim(phdir)//trim(dispn_folder)//"supercell-"//trim(prefix))
          endif
          !
        endif
        !
        tau_sc(id,ia) = tau(id,ia) ! restore equilibrium position
        !
      enddo
    enddo

  end subroutine create_displacements


  subroutine check_kgrid(kscell, displ)
    !
    ! Check if the Monkhorst-Pack kgrid is defined. If it is defined returns it
    ! in kscell and displ, if it is not defined, there is no need to modify it
    ! for the supercell and returns empty arrays kscell and displ.
    !
    use fdf, only: block_fdf, parsed_line, fdf_block, fdf_bline, fdf_bclose, &
                   fdf_bintegers, fdf_bnvalues, fdf_bvalues
    !
    implicit none
    !
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline
    !
    integer, intent(out) :: kscell(3,3)
    real(kind=dp), intent(out) :: displ(3)
    !
    integer :: i


    ! Default values
    kscell = 0
    displ = 0.0_dp

    ! Read kgridMonkhorstPack block
    if ( fdf_block("kgridMonkhorstPack",bfdf) ) then
      !
      do i= 1, 3
        !
        if ( fdf_bline(bfdf,pline) ) then
          kscell(1,i) = fdf_bintegers(pline,1)
          kscell(2,i) = fdf_bintegers(pline,2)
          kscell(3,i) = fdf_bintegers(pline,3)
          if ( fdf_bnvalues(pline) > 3 ) then
            displ(i) = mod(fdf_bvalues(pline,4), 1._dp)
          else
            displ(i) = 0._dp
          end if
        else
          stop "kgridinit: ERROR in kgridMonkhorstPack block"
        endif
        !
      enddo
      call fdf_bclose(bfdf)
    endif

  end subroutine check_kgrid


  subroutine modify_kgrid(kscell, displ, fdffilename, backup_in)
    !
    ! Modifys the fdf fdffilename to set kgridMonkhorstPack to kscell and displ
    !
    use, intrinsic :: iso_fortran_env, only: iostat_end ! end of file
    use fdf, only: fdf_defined
    !
    use siesta2ph_io, only: find_free_unit
    use siesta2ph_utils, only: word_in_string
    !
    implicit none
    !
    integer, intent(in) :: kscell(3,3)
    real(kind=dp), intent(in) :: displ(3)
    character(len=*), intent(in) :: fdffilename
    logical, optional, intent(in) :: backup_in
    !
    character(len=256) :: line
    integer :: iounit_reference, iounit, ios
    logical :: backup, defined, found
    integer :: id


    ! Default value for optional variable
    if (present(backup_in)) then
      backup = backup_in
    else
      backup = .false.
    endif
    !
    ! Check if kgridMonkhorstPack is defined in the fdf
    if (fdf_defined("kgridMonkhorstPack")) then
      defined = .true.
    else
      defined = .false.
    endif
    !
    ! Backup the original fdf
    call execute_command_line("cp "//trim(outdir)//trim(fdffilename)//" "//trim(outdir)//trim(fdffilename)//"_backup")
    !
    ! Modify kgrid
    iounit_reference = find_free_unit()
    open(unit=iounit_reference, file=trim(outdir)//trim(fdffilename)//"_backup", status="old", action="read", iostat=ios)
    if ( ios /= 0 ) stop "modify_kgrid: Error opening backup fdf file"
    !
    iounit = find_free_unit()
    open(unit=iounit, file=trim(outdir)//trim(fdffilename), status="replace", action="write", iostat=ios)
    if ( ios /= 0 ) stop "modify_kgrid: Error opening new fdf file"
    !
    found = .false.
    do
      !
      ! Read each line of the input file
      read(unit=iounit_reference, fmt="(a)", iostat=ios) line
      !
      ! Exit the do loop if it is the last line
      if ( ios == iostat_end ) exit
      !
      if (defined .and. word_in_string("kgridMonkhorstPack",packlabel(line)) .and. word_in_string("%block",line)) then
        !
        found = .true.
        write(iounit,"(a)") trim(line)
        do id=1,3
          write(iounit,"(3i4,f6.1)") kscell(:,id), displ(id)
        enddo
        !
        do
          read(unit=iounit_reference, fmt="(a)", iostat=ios) line
          if (word_in_string("kgridMonkhorstPack",packlabel(line)) .and. word_in_string("%endblock",line)) exit
        enddo
        !
        write(iounit,"(a)") trim(line)
        !
      else
        !
        write(iounit,"(a)") trim(line)
        !
      endif
      !
    end do
    !
    if ( defined .and. (.not. found) ) stop "ERROR: modify_kgrid: kgridMonkhorstPack defined but not found"
    !
    close(iounit)
    close(iounit_reference)
    !
    if (backup) then
      write(stdout,*) "- Backup "//trim(outdir)//trim(fdffilename)//"..."
    else
      call execute_command_line("rm "//trim(outdir)//trim(fdffilename)//"_backup")
    endif

  end subroutine modify_kgrid


  subroutine check_mesh_sizes(suggested_mesh_sizes)
    !
    ! Check if the unit cell and the supercell have a compatible
    ! Mesh Sizes. If the Mesh Sizes are compatible returns an
    ! empty suggested_mesh_sizes array, if they are not compatible,
    ! suggested_mesh_sizes returns a compatible Mesh Sizes for the
    ! unit cell.
    !
    use fdf, only: block_fdf, parsed_line, fdf_defined, fdf_islist, fdf_list, &
                   fdf_block, fdf_bline, fdf_bclose, fdf_bintegers, fdf_get
    !
    use siesta2ph_io, only: stdout, nr1, nr2, nr3
    !
    use siesta2ph_linalg, only: ainv
    !
    implicit none
    !
    integer, dimension(3), intent(out) :: suggested_mesh_sizes
    !
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline
    !
    real(kind=dp) :: suggested_g2cut
    integer, dimension(3) :: mesh_sizes_uc, mesh_sizes_sc, mesh_sizes
    integer :: size
    real(kind=dp) :: g2cut
    logical :: compatible


    ! Default output value
    suggested_mesh_sizes = (/0, 0, 0/)
    !
    if ( fdf_defined("MeshSizes") ) then
      !
      ! Mesh sizes is given in the fdf
      !
      if (fdf_islist("MeshSizes") ) then
        !
        size = -1
        call fdf_list("MeshSizes", size, suggested_mesh_sizes)
        if ( size /= 3 ) stop "check_mesh_sizes: ERROR in MeshSizes list"
        ! Do the actual read of the mesh sizes
        call fdf_list("MeshSizes", size, suggested_mesh_sizes)
        !
      else if ( fdf_block("MeshSizes", bfdf) ) then
        !
        if ( fdf_bline(bfdf, pline) ) then
          suggested_mesh_sizes(1) = fdf_bintegers(pline,1)
          suggested_mesh_sizes(2) = fdf_bintegers(pline,2)
          suggested_mesh_sizes(3) = fdf_bintegers(pline,3)
        else
            stop "check_mesh_sizes: ERROR in MeshSizes block"
        endif
        call fdf_bclose(bfdf)
        !
      endif
      !
      return
      !
    else
      !
      ! Mesh sizes is not given in the fdf
      !
      ! Get mesh cutoff
      g2cut = fdf_get("MeshCutoff", 300.0_dp, "Ry")
      !
      ! Compute mesh sizes for the unit cell and the supercell
      call get_mesh_sizes(at, g2cut, mesh_sizes_uc)
      call get_mesh_sizes(at_sc, g2cut, mesh_sizes_sc)
      mesh_sizes = int(dble(mesh_sizes_sc)/(/nr1, nr2, nr3/))
      !
      ! Check if mesh sizes is compatible with the supercell
      compatible = ( any(mesh_sizes == mesh_sizes_uc) .and. &
                     any(abs(dble(mesh_sizes_sc)/(/nr1, nr2, nr3/) - dble(mesh_sizes)) < 0.0001_dp) )
      !
      if (compatible) return
      !
      ! Print a warning message
      write(stdout, "(a)") "WARNING: mesh sizes are not compatible:"
      write(stdout, "(a15,3i5)") "  mesh_sizes_uc:", mesh_sizes_uc
      write(stdout, "(a15,3i5)") "  mesh_sizes_sc:", mesh_sizes_sc
      !
      ! Try to find a MeshCutoff that has compatible mesh sizes
      suggested_g2cut = g2cut
      do
        !
        suggested_g2cut = suggested_g2cut + 1.0_dp
        call get_mesh_sizes(at, suggested_g2cut, mesh_sizes_uc)
        call get_mesh_sizes(at_sc, suggested_g2cut, mesh_sizes_sc)
        mesh_sizes = int(dble(mesh_sizes_sc)/(/nr1, nr2, nr3/))
        !
        compatible = ( any(mesh_sizes == mesh_sizes_uc) .and. &
                       any(abs(dble(mesh_sizes_sc)/(/nr1, nr2, nr3/) - dble(mesh_sizes)) < 0.0001_dp) )
        !
        if (compatible) exit
        if (suggested_g2cut > 2*g2cut) exit
        !
      enddo
      !
      if (compatible) then
        write(stdout, "(a36,f6.1,a3,f6.1)") "WARNING: Increaseing MeshCutoff from", g2cut, " to", suggested_g2cut
        call modify_MeshCutoff(suggested_g2cut, trim(outdir)//trim(prefix), .true.)
        return
      endif
      !
      ! If MeshCutoff was not found try to find mesh sizes manually
      suggested_mesh_sizes = nint(dble(mesh_sizes_uc)/(2*3*5))*2*3*5
      !
      write(stdout, "(a27,3i5,a1)") "WARNING: Using MeshSizes [", suggested_mesh_sizes, "]"
      call modify_MeshSizes(suggested_mesh_sizes, trim(outdir)//trim(prefix), .true.)
      !
    endif

  end subroutine check_mesh_sizes


  subroutine modify_MeshCutoff(g2cut, fdffilename, backup_in)
    !
    ! Modifys the fdf fdffilename to set the mesh cutoff to g2cut
    !
    use, intrinsic :: iso_fortran_env, only: iostat_end ! end of file
    use fdf, only: fdf_defined
    !
    use siesta2ph_io, only: find_free_unit
    use siesta2ph_utils, only: word_in_string
    !
    implicit none
    !
    real(kind=dp), intent(in) :: g2cut
    character(len=*), intent(in) :: fdffilename
    logical, optional, intent(in) :: backup_in
    !
    character(len=256) :: line
    integer :: iounit_reference, iounit, ios
    logical :: backup, defined, found


    ! Default value for optional variable
    if (present(backup_in)) then
      backup = backup_in
    else
      backup = .false.
    endif
    !
    ! Check if MeshCutoff is defined in the fdf
    if (fdf_defined("MeshCutoff")) then
      defined = .true.
    else
      defined = .false.
    endif
    !
    ! Backup the original fdf
    call execute_command_line("cp "//trim(outdir)//trim(fdffilename)//" "//trim(outdir)//trim(fdffilename)//"_backup")
    !
    ! Modify MeshCutoff
    iounit_reference = find_free_unit()
    open(unit=iounit_reference, file=trim(outdir)//trim(fdffilename)//"_backup", status="old", action="read", iostat=ios)
    if ( ios /= 0 ) stop "modify_MeshCutoff: Error opening backup fdf file"
    !
    iounit = find_free_unit()
    open(unit=iounit, file=trim(outdir)//trim(fdffilename), status="replace", action="write", iostat=ios)
    if ( ios /= 0 ) stop "modify_MeshCutoff: Error opening new fdf file"
    !
    found = .false.
    do
      !
      ! Read each line of the input file
      read(unit=iounit_reference, fmt="(a)", iostat=ios) line
      !
      ! Exit the do loop if it is the last line
      if ( ios == iostat_end ) exit
      !
      if (defined .and. word_in_string("MeshCutoff",packlabel(line))) then
        !
        found = .true.
        write(iounit,"(a10,f8.1,a3)") "MeshCutoff", g2cut, " Ry"
        !
      else
        !
        write(iounit,"(a)") trim(line)
        !
      endif
      !
    end do
    !
    if ( defined .and. (.not. found) ) stop "ERROR: modify_MeshCutoff: MeshCutoff defined but not found"
    !
    if (.not. found) write(iounit,"(a,f8.1,a3)") "MeshCutoff", g2cut, " Ry"
    !
    close(iounit)
    close(iounit_reference)
    !
    if (backup) then
      write(stdout,*) "- Backup "//trim(outdir)//trim(fdffilename)//"..."
    else
      call execute_command_line("rm "//trim(outdir)//trim(fdffilename)//"_backup")
    endif

  end subroutine modify_MeshCutoff


  subroutine modify_MeshSizes(mesh_sizes, fdffilename, backup_in)
    !
    ! Modifys the fdf fdffilename to set the mesh sizes to mesh_sizes
    !
    use, intrinsic :: iso_fortran_env, only: iostat_end ! end of file
    use fdf, only: fdf_islist, fdf_isblock
    !
    use siesta2ph_io, only: find_free_unit
    use siesta2ph_utils, only: word_in_string
    !
    implicit none
    !
    integer, dimension(3), intent(in) :: mesh_sizes
    character(len=*), intent(in) :: fdffilename
    logical, optional, intent(in) :: backup_in
    !
    character(len=256) :: line
    integer :: iounit_reference, iounit, ios
    logical :: backup, defined_list, defined_block, found


    ! Default value for optional variable
    if (present(backup_in)) then
      backup = backup_in
    else
      backup = .false.
    endif
    !
    ! Check if Mesh.Sizes is defined in the fdf
    defined_list = fdf_islist("MeshSizes")
    defined_block = fdf_isblock("MeshSizes")
    !
    ! Backup the original fdf
    call execute_command_line("cp "//trim(outdir)//trim(fdffilename)//" "//trim(outdir)//trim(fdffilename)//"_backup")
    !
    ! Modify MeshCutoff
    iounit_reference = find_free_unit()
    open(unit=iounit_reference, file=trim(outdir)//trim(fdffilename)//"_backup", status="old", action="read", iostat=ios)
    if ( ios /= 0 ) stop "modify_MeshSizes: Error opening backup fdf file"
    !
    iounit = find_free_unit()
    open(unit=iounit, file=trim(outdir)//trim(fdffilename), status="replace", action="write", iostat=ios)
    if ( ios /= 0 ) stop "modify_MeshSizes: Error opening new fdf file"
    !
    found = .false.
    do
      !
      ! Read each line of the input file
      read(unit=iounit_reference, fmt="(a)", iostat=ios) line
      !
      ! Exit the do loop if it is the last line
      if ( ios == iostat_end ) exit
      !
      if (word_in_string("MeshSizes",packlabel(line))) then
        !
        if (defined_block .and. word_in_string("%block", line)) then
          !
          found = .true.
          !
          write(iounit,"(a)") trim(line)
          write(iounit,"(3i5)") mesh_sizes
          do
            read(unit=iounit_reference, fmt="(a)", iostat=ios) line
            if (word_in_string("%endblock", line)) exit
          end do
          write(iounit,"(a)") trim(line)
          !
        else if (defined_list .and. (.not.word_in_string("%block", line))) then
          !
          found = .true.
          !
          write(iounit,"(a,3i5,a)") "MeshSizes [", mesh_sizes, "]"
          !
        else
          write(iounit,"(a)") trim(line)
        endif
        !
      else
        !
        write(iounit,"(a)") trim(line)
        !
      endif
      !
    end do
    !
    if ( (defined_list.or.defined_block) .and. (.not. found) ) stop "ERROR: modify_MeshSizes: MeshSizes defined but not found"
    !
    if (.not. found ) write(iounit,"(a,3i5,a)") "MeshSizes [", mesh_sizes, "]"
    !
    close(iounit)
    close(iounit_reference)
    !
    if (backup) then
      write(stdout,*) "- Backup "//trim(outdir)//trim(fdffilename)//"..."
    else
      call execute_command_line("rm "//trim(outdir)//trim(fdffilename)//"_backup")
    endif

  end subroutine modify_MeshSizes


  subroutine get_mesh_sizes(cell, cutoff, ntm)
    !
    ! Computes the mesh dimensions ntm for a given cell
    ! and mesh cutoff
    !
    use precision, only: dp
    use m_chkgmx, only: chkgmx ! Checks planewave cutoff of a mesh
    use fft1d, only: nfft ! Finds allowed value for 1-D FFT
    use fdf, only: fdf_integer
    !
    use siesta2ph_io, only: stdout
    use siesta2ph_utils, only: tpi
    !
    use siesta2ph_linalg, only: ainv
    !
    implicit none
    !
    real(kind=dp), intent(in) :: cell(3,3)
    real(kind=dp), intent(in) :: cutoff
    integer, intent(out) :: ntm(3)
    !
    real(kind=dp) :: RealCutoff
    integer :: nsm
    real(kind=dp) :: rcell(3,3), vecmod
    integer :: i
    real(kind=dp), parameter :: k0(3) = (/ 0.0, 0.0, 0.0 /)


    ! Get mesh sub-divisions
    nsm = fdf_integer("MeshSubDivisions", 2)
    nsm = max(nsm, 1)

    ! Find number of mesh intervals for each cell vector.
    ! The reciprocal vectors of the mesh unit cell (cell/ntm)
    ! are rcell*ntm, and must be larger than 2*cutoff
    rcell = transpose(ainv(cell))*tpi/alat
    do i = 1,3
      vecmod = sqrt(dot_product(rcell(:,i),rcell(:,i)))
      ntm(i) = 2 * sqrt(cutoff) / vecmod + 1
    enddo

    ! Impose that mesh cut-off is large enough
    impose_cutoff: do

      ! Require that ntm is suitable for FFTs and a multiple of nsm
      do i = 1,3
        impose_fft: do
          ! nfft increases ntm to the next integer suitable for FFTs
          call nfft( ntm(i) )
          if ( mod( ntm(i), nsm )==0 ) exit impose_fft
          ntm(i) = ntm(i) + 1
        end do impose_fft
      enddo ! i

      ! Check that effective cut-off is large enough as for non-right
      ! angled unit cells this is not guaranteed to be the case.
      ! If cut-off needs to be larger, increase ntm and try again.
      ! chkgmx decreases G2mesh to the right cutoff value of the mesh
      RealCutoff = huge(1.0_dp)
      call chkgmx( k0, rcell, ntm, RealCutoff )
      if (RealCutoff >= cutoff) then
        exit impose_cutoff
      else
        ntm(1:3) = ntm(1:3) + 1
      end if

    end do impose_cutoff

  end subroutine get_mesh_sizes


  function packlabel(string)
    !
    ! Interface function to SIESTA's packable subroutine
    !
    use fdf_utils, only: packlabel_siesta => packlabel
    !
    implicit none
    !
    character(*), intent(in) :: string
    character(len(string)) :: packlabel

    call packlabel_siesta(string, packlabel)

  end function packlabel

end program siesta2ph
