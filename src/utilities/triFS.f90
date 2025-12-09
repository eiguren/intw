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
module triFS_input_parameters

  !! display: public
  !!
  !! This module contains the definitions of input parameters for [[triFS]] utility.
  !!
  !! ### Details
  !!
  !! #### Input parameters
  !!
  !! ```{.txt}
  !! &tri_FS
  !!     n1           = integer
  !!     n2           = integer
  !!     n3           = integer
  !!     volume_nodes = T or F (Default: T)
  !!     volnodfac    = real (Default: 1.0)
  !!     hr_file      = 'file'
  !!     ef           = real (Default: 0.0 eV)
  !!     verbose      = T or F (Default: F)
  !!     plot_BZ      = T or F (Default: T)
  !!     dos          = T or F (Default: T)
  !!     eps_dupv     = real (Default: 1.0E-06)
  !! /
  !! &FS_opt
  !!     collapse          = T or F (Default: T)
  !!     collapse_criteria = real (Default: 0.2)
  !!     relax             = T or F (Default: T)
  !!     relax_iter        = integer (Default: 1000)
  !!     newton_raphson    = 0: not applied, 1: only in the end, 2: beginning and end (Default: 2)
  !!     newton_iter       = integer (Default: 10)
  !!     relax_vinface     = T or F (Default: F)
  !!     eps_vinface       = real (Default: 1.0E-5)
  !! /
  !! ```
  !!
  !! Variables without default values must be explicitly set in the input
  !! file; otherwise, an error will be raised.
  !!

  use kinds, only: dp

  implicit none

  public :: tri_FS, TR_sym, n1, n2, n3, volume_nodes, volnodfac, hr_file, ef, verbose, plot_BZ, dos, eps_dupv
  public :: FS_opt, collapse, collapse_criteria, relax, relax_iter, newton_raphson, newton_iter, relax_vinface, eps_vinface

  private

  ! tri_FS

  logical :: TR_sym = .true.
  !! TR symmetry operation
  integer :: n1 = -1, n2 = -1, n3 = -1
  !! BZ sampling for creating tetrahedra
  logical :: volume_nodes = .true.
  !! .true. if nodes are to be added inside the IBZ volume as Steiner points
  real(kind=dp) :: volnodfac = 1.0_dp
  !! Factor with which multiply n1, n2, n3 to add points within IBZ volume
  character(len=256) :: hr_file = 'unassigned'
  !! Wannier90 $seedname.hr file
  real(kind=dp) :: ef = 0.0_dp
  !! Fermi level (Units: eV)
  logical :: verbose = .false.
  !! Verbose output
  logical :: plot_BZ = .true.
  !! Plot BZ edges
  logical :: dos = .true.
  !! Compute Fermi level DOS
  real(kind=dp) :: eps_dupv = 1.0E-06_dp
  !! Threshold parameter to detect duplicated vertices

  ! FS_opt

  logical :: collapse = .true.
  !! Collapse edges of triangulated mesh
  real(kind=dp) :: collapse_criteria = 0.2_dp
  !! Shortest edge to longest edge ratio to proceed with collapse
  logical :: relax = .true.
  !! Tangentially relax edges of triangulated mesh
  integer :: relax_iter = 1000
  !! Number of tangential relax iterations
  integer :: newton_raphson = 2
  !! Newthon-Raphson relax edges of triangulated mesh:
  !!
  !! - `0`: not applied
  !! - `1`: applied only at the end
  !! - `2`: applied at beginning and end
  integer :: newton_iter = 10
  !! Number of Newthon-Raphson relax iterations
  logical :: relax_vinface = .false.
  !! Relax vertices on faces. May give erros in some examples.
  real(kind=dp) :: eps_vinface = 1.0E-5_dp
  !! Threshold to detect vertices on BZ faces

  ! Input namelists
  NAMELIST / tri_FS / TR_sym, n1, n2, n3, volume_nodes, volnodfac, hr_file, ef, verbose, plot_BZ, dos, eps_dupv
  NAMELIST / FS_opt / collapse, collapse_criteria, relax, relax_iter, newton_raphson, newton_iter, relax_vinface, eps_vinface

end module triFS_input_parameters
!
program triFS

  !! display: none
  !!
  !! Triangulate Fermi surface.
  !!
  !! ### Details
  !!
  !! #### Input parameters
  !!
  !! ```{.txt}
  !! &input
  !!     outdir = 'directory'
  !!     prefix = 'prefix'
  !! /
  !! &tri_FS
  !!     TR_sym = T or F
  !!     n1           = integer
  !!     n2           = integer
  !!     n3           = integer
  !!     volume_nodes = T or F
  !!     volnodfac    = real
  !!     hr_file      = 'file'
  !!     ef           = real
  !!     verbose      = T or F
  !!     plot_BZ      = T or F
  !!     dos          = T or F
  !!     eps_dupv     = real
  !! /
  !! &FS_opt
  !!     collapse          = T or F
  !!     collapse_criteria = real
  !!     relax             = T or F
  !!     relax_iter        = integer
  !!     newton_raphson    = 0: not applied, 1: only in the end, 2: beginning and end
  !!     newton_iter       = integer
  !!     relax_vinface     = T or F
  !!     eps_vinface       = real
  !! /
  !! ```
  !!
  !! `triFS.x` utility uses special `tri_FS` and `FS_opt` namelists,
  !! see [[triFS_input_parameters]] module for the description of each parameter in those namelists.
  !! See [[intw_input_parameters]] module for the description of each parameter in the `input` namelist.
  !!

  use kinds, only: dp
  use triFS_input_parameters
  use triFS_geometry, only: wigner_seitz_cell, plot_poly, polyhedra_off, tetrasym, plot_tetra_off, compact_tetra, &
                      tetraIBZ_2_vert_faces_edges, irr_faces, triangulate_faces, add_nodes_IBZ_volume
  use triFS_isosurface, only: vert_veloc_rot, vert_veloc, vert_index_rot, vert_coord_rot, ntri_rot, nvert_rot, &
                        vert_index, vert_coord, nvert, ntri, &
                        read_tetrahedra, create_isosurface_IBZ, write_full_isosurface, DOS_isosurface, &
                        write_IBZ_isosurface
  use triFS_mesh_opt, only: mesh_optimization
  use intw_version, only: print_intw_version
  use intw_utility, only: find_free_unit, get_timing, print_threads, print_date_time
  use intw_input_parameters, only: input, outdir, prefix
  use intw_reading, only: read_parameters_data_file, alat, at, bg, volume0, nsym, s

  implicit none

  ! Tetrahedra variables
  integer, parameter :: ntetmax = 400
  integer :: n_BZ_tetra_irr, n_BZ_tetra_all, tetra_equiv(ntetmax), tetra_symlink(ntetmax,1:2)
  real(kind=dp), dimension(1:3, 1:4, ntetmax) :: BZ_tetra_irr, BZ_tetra_all
  integer :: nvert_IBZ, nfaces_IBZ, nedges_IBZ, nfaces_irr
  real(kind=dp), allocatable :: vert_IBZ(:,:)
  integer, allocatable :: faces_IBZ_as_vert(:,:), edges_IBZ(:,:)
  integer, allocatable :: faces_Gsymlink(:,:,:,:,:), faces_indx(:), faces_inv_indx(:)

  ! Wannier and isosurface variables
  integer :: hr_unit, num_wann, irpt, nrpts, iwann, jwann
  integer, dimension(:), allocatable :: ndegen
  integer, dimension(:,:), allocatable :: irvec
  complex(kind=dp), dimension(:,:,:), allocatable :: ham_r

  integer :: i, j
  character(len=25) :: tag_in
  integer :: ios1, ios2, ios3
  real(dp) :: time1, time2


  20 format(A)
  30 format(A,F8.2,6X,A)
  !
  !
  !================================================================================
  ! Begining
  !================================================================================
  !
  call get_timing(time1)
  !
  write(*,20) '====================================================='
  write(*,20) '|                   program triFS                   |'
  write(*,20) '|         ---------------------------------         |'
  call print_intw_version()
  call print_threads()
  call print_date_time("Start of execution")
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Read the necessary information from standard input file
  !================================================================================
  !
  22 format(A,I2,A)
  !
  write(*,20) "| - Reading standard input file...                  |"
  write(*,20) "|         namelist             ios                  |"
  write(*,20) "|         --------             ----                 |"
  !
  ! Read input file
  READ(5, nml=input, iostat=ios1)
  write(*,22) "|           &input             ", ios1, "                   |"
  !
  READ(5, nml=tri_FS, iostat=ios2)
  write(*,22) "|           &triFS             ", ios2, "                   |"
  !
  READ(5, nml=FS_opt, iostat=ios3)
  write(*,22) "|           &FS_opt            ", ios3, "                   |"
  !
  ! Check input
  if (ios1 /= 0) then
    write(*,*) "PLEASE CHECK INPUT NAMELIST  &input  as it is not OK!"
    stop
  else if (ios2 /=0) then
    write(*,*) "PLEASE CHECK INPUT NAMELIST  &triFS  as it is not OK!"
    stop
  else if (ios3 /= 0) then
    write(*,*) "PLEASE CHECK INPUT NAMELIST  &FS_opt  as it is not OK!"
    stop
  end if
  !
  ! Check input parameters
  if ( outdir == 'unassigned' ) stop 'MISSING outdir!'
  !
  if ( prefix == 'unassigned' ) stop 'MISSING prefix!'
  !
  if ( hr_file == 'unassigned' ) stop 'MISSING hr_file!'
  !
  if ( n1 == -1 .or. n2 == -1 .or. n3 == -1) stop 'MISSING n1, n2, n3!'
  !
  if ( newton_raphson /= 0 .and. newton_raphson /= 1 .and. newton_raphson /= 2 ) stop 'INVALID newton_raphson!'
  !
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Read the parameters from the SCF calculation
  !================================================================================
  !
  write(*,20) "| - Reading calculation parameters...               |"
  !
  call read_parameters_data_file()
  !
  !
  !================================================================================
  ! Read _hr file
  !================================================================================
  !
  write(*,20) "| - Reading Wannier H(R)...                         |"
  !
  hr_unit = find_free_unit()
  open(unit=hr_unit, file=hr_file, status="unknown", action="read")
  read(unit=hr_unit, fmt=*)
  read(unit=hr_unit, fmt=*) num_wann
  read(unit=hr_unit, fmt=*) nrpts
  allocate(ndegen(nrpts), irvec(1:3,nrpts), ham_r(num_wann,num_wann,nrpts))
  read(unit=hr_unit, fmt='(15I5)') (ndegen(i), i=1,nrpts)
  do irpt=1,nrpts
     do i=1,num_wann
        do j=1,num_wann
           ! JL: This is the accuracy on output hr file of wannier90, it seems too low...
           read(unit=hr_unit, fmt='(5I5,2F12.6)') irvec(:,irpt), jwann, iwann, ham_r(j,i,irpt)
           !read(unit=hr_unit, fmt='(5I5,2ES18.10)') irvec(:,irpt), jwann, iwann, ham_r(j,i,irpt)
        end do
     end do
  end do
  close(unit=hr_unit)
  !
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Create BZ, tetrahedralize with symmetric tetra, and find IBZ
  !================================================================================
  !
  ! Create BZ
  !
  write(*,20) "| - Creating WS BZ...                               |"
  !
  call wigner_seitz_cell(transpose(bg), verbose)
  !
  ! Write BZ border lines in plotting format if needed
  if(plot_BZ) call plot_poly()
  !
  ! Write BZ in OFF format (reads polyhedra.dat) ! JL this should be improved to write directly in OFF format...
  call polyhedra_off()
  !
  ! Find tetrahedra forming the irreducible BZ volume
  !
  write(*,20) "| - Finding IBZ...                                  |"
  !
  call tetrasym(bg, nsym, s, TR_sym, ntetmax, n_BZ_tetra_all, BZ_tetra_all, n_BZ_tetra_irr, BZ_tetra_irr, tetra_equiv, tetra_symlink)
  !
  if (verbose) then
    ! Write tetrahedized full BZ
    tag_in = "tetra_BZ.off"
    call plot_tetra_off(tag_in, n_BZ_tetra_all, BZ_tetra_all, plot_BZ)
  end if
  !
  ! Compact tetrahedra on IBZ
  call compact_tetra(bg, nsym, s, TR_sym, ntetmax, n_BZ_tetra_irr, BZ_tetra_irr, n_BZ_tetra_all, BZ_tetra_all, tetra_equiv, tetra_symlink)
  !
  ! Write compact tetrahedralized irreducible BZ
  tag_in = "IBZ.off"
  call plot_tetra_off(tag_in, n_BZ_tetra_irr, BZ_tetra_irr, .false.)
  !
  write(*,20) "|   ...done                                         |"
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Detect irreducible faces within IBZ with S+G symmetries
  !================================================================================
  !
  write(*,20) "| - Detecting equivalent faces of IBZ...            |"
  !
  allocate(vert_IBZ(1:3,4*n_BZ_tetra_irr), faces_IBZ_as_vert(1:3,4*n_BZ_tetra_irr), edges_IBZ(1:2,6*n_BZ_tetra_irr))
  call tetraIBZ_2_vert_faces_edges(n_BZ_tetra_irr, BZ_tetra_irr, verbose, nvert_IBZ, nfaces_IBZ, nedges_IBZ, vert_IBZ, faces_IBZ_as_vert, edges_IBZ)
  !
  allocate(faces_Gsymlink(4*n_BZ_tetra_irr,2,-1:1,-1:1,-1:1), faces_indx(4*n_BZ_tetra_irr), faces_inv_indx(4*n_BZ_tetra_irr))
  call irr_faces(n_BZ_tetra_irr, nsym, s, TR_sym, bg, verbose, vert_IBZ, nfaces_IBZ, faces_IBZ_as_vert, nfaces_irr, faces_Gsymlink, faces_indx, faces_inv_indx)
  !
  write(*,20) "|   ...done                                         |"
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Split edges and perform triangulation of irr faces
  !================================================================================
  !
  write(*,20) "| - Triangulating faces of IBZ...                   |"
  !
  ! This uses Triangle executable
  call triangulate_faces(n_BZ_tetra_irr, nfaces_IBZ, faces_Gsymlink, nsym, s, faces_indx, faces_inv_indx, n1, n2, n3, bg, &
                         faces_IBZ_as_vert, vert_IBZ, verbose)
  !
  ! Remove unnecessary files
  CALL EXECUTE_COMMAND_LINE("rm face_split_edges.*")
  !
  write(*,20) "|   ...done                                         |"
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Add extra n1, n2, n3 nodes within IBZ volume
  !================================================================================
  !
  if (volume_nodes) then
    !
    write(*,20) "| - Adding n1, n2, n3 nodes to IBZ volume...        |"
    !
    call add_nodes_IBZ_volume(n1, n2, n3, volnodfac, eps_vinface, bg, n_BZ_tetra_irr, BZ_tetra_irr, verbose)
    !
  write(*,20) "|   ...done                                         |"
  write(*,20) "====================================================="
    !
  end if
  !
  !
  !================================================================================
  ! Pass triangulated IBZ border and extra nodes to TetGen and tetrahedralize
  !================================================================================
  !
  write(*,20) "| - Creating tetrahedra in IBZ volume...            |"
  !
  ! Call TetGen
  ! -Y  doesn't let adding points on surface
  ! -i uses adds extra nodes from .node file
  ! -v verbose output
  if (verbose) write(*,20) "|   Running tetgen with command:                    |"
  if (volume_nodes) then
    if (verbose) write(*,20) "|  tetgen -Ykv -i Triangulated_IBZ.off > tetgen.out |"
    CALL EXECUTE_COMMAND_LINE("tetgen -Ykv -i Triangulated_IBZ.off > tetgen.out")
  else
    if (verbose) write(*,20) "|   tetgen -Ykv Triangulated_IBZ.off > tetgen.out   |"
    CALL EXECUTE_COMMAND_LINE("tetgen -Ykv Triangulated_IBZ.off")
  end if
  !
  ! Re-name files
  CALL EXECUTE_COMMAND_LINE("rm Triangulated_IBZ.1.smesh")
  CALL EXECUTE_COMMAND_LINE("mv Triangulated_IBZ.1.vtk Tetrahedralized_IBZ.vtk")
  CALL EXECUTE_COMMAND_LINE("mv Triangulated_IBZ.1.node Tetrahedralized_IBZ.node")
  CALL EXECUTE_COMMAND_LINE("mv Triangulated_IBZ.1.face Tetrahedralized_IBZ.face")
  CALL EXECUTE_COMMAND_LINE("mv Triangulated_IBZ.1.ele Tetrahedralized_IBZ.ele")
  CALL EXECUTE_COMMAND_LINE("mv Triangulated_IBZ.1.edge Tetrahedralized_IBZ.edge")
  !
  write(*,20) "|   ...done                                         |"
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Read tetrahedra output and obtain isosurface
  !================================================================================
  !
  write(*,20) '| - Creating FS on IBZ...                           |'
  write(*,20) "|         ---------------------------------         |"
  !
  ! Read small tetrahedra from output files. Output variables defined on isosurface module
  call read_tetrahedra()
  !
  write(*,20) "|         ---------------------------------         |"
  !
  ! Create isosurface on IBZ
  call create_isosurface_IBZ(ef, num_wann, nrpts, ndegen, irvec, ham_r, alat, at, bg, nsym, s, TR_sym, verbose, eps_dupv)
  !
  ! Write triangulated mesh in OFF format
  write(*,20) "|         ---------------------------------         |"
  tag_in = "initial_IBZ_FS_tri"
  call write_IBZ_isosurface(tag_in, num_wann, .false.)
  !
  write(*,20) "|         ---------------------------------         |"
  write(*,20) "|   ...done                                         |"
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Mesh optimization
  !================================================================================
  !
  write(*,20) "| - Optimizing FS...                                |"
  !
  if(collapse .or. relax .or. (newton_raphson>0)) then
    !
    call mesh_optimization(collapse, relax, newton_raphson, collapse_criteria, relax_iter, newton_iter, relax_vinface, &
                           eps_vinface, eps_dupv, verbose, ef, nrpts, irvec, ndegen, ham_r, alat, at, bg, nsym, s, TR_sym, &
                           num_wann, nfaces_IBZ, faces_IBZ_as_vert, vert_IBZ, ntri, nvert, vert_coord, vert_index, vert_veloc)
    !
    write(*,20) "|         ---------------------------------         |"
    write(*,20) '|   ...done                                         |'
    !
  else
    !
    write(*,20) "|   ...nothing to do                                |"
    !
  endif
  !
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Rotate isosurface to full BZ and write to file
  !================================================================================
  !
  write(*,20) "| - Writing final triangulated FS...                |"
  !
  write(*,20) "|         ---------------------------------         |"
  !
  ! write isosurface on IBZ
  tag_in = "IBZ_FS_tri"
  call write_IBZ_isosurface(tag_in, num_wann, .true., prefix)
  !
  write(*,20) "|         ---------------------------------         |"
  !
  ! write isosurface on full BZ
  tag_in = "FS_tri"
  call write_full_isosurface(bg, nsym, s, TR_sym, num_wann, verbose, eps_dupv, tag_in, prefix)
  !
  write(*,20) "|         ---------------------------------         |"
  write(*,20) "|   ...done                                         |"
  write(*,20) "====================================================="
  !
  !
  !================================================================================
  ! Compute DOS at isosurface
  !================================================================================
  !
  if (dos) then
    !
    write(*,20) "| - Computing DOS at FS...                          |"
    !
    call DOS_isosurface(alat, alat**3/volume0, num_wann, nvert, ntri, vert_coord, vert_index, vert_veloc, nvert_rot, ntri_rot, &
                        vert_coord_rot, vert_index_rot, vert_veloc_rot)
    !
    write(*,20) "|   ...done                                         |"
    write(*,20) "====================================================="
  end if
  !
  !
  !
  !================================================================================
  ! Finish
  !================================================================================
  !
  call get_timing(time2)
  !
  write(*,20) '|                      ALL DONE                     |'
  write(*,30) '|     Total time: ',time2-time1,' seconds            |'
  call print_date_time('End of execution  ')
  write(*,20) '====================================================='

end program triFS
