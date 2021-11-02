!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
module intw_setup
!----------------------------------------------------------------------------!
!
!       This module contains the subroutine "setup", which basically sets the
!       stage for further computations. 
!       
!       It also contains subroutines which perform tests for the
!       setup subprogram, insuring the proper behavior of tested subroutines.
!----------------------------------------------------------------------------!

use intw_input_parameters
use intw_utility
use intw_useful_constants
use intw_reading
use intw_symmetries
use intw_fft
use intw_W90
  !
  implicit none
  !
  save
  !
  real(dp), allocatable :: kpoints_QE(:,:)
  real(dp), allocatable :: kmesh(:,:)

  
contains

  subroutine setup_QE(time)
  !--------------------------------------------------------------------------------
  !
  !       The purpose of this subroutine is to set up the necessary information
  !       structures for computation, only using Quantum Espresso data. 
  !       It will read in the input file and generate the data stuctures necessary 
  !	  for the reading and rotation of wave functions. It will not
  !       use any data structures from Wannier90. This can be useful 
  !	  when performing tests when a W90 calculation is not necessary.
  !--------------------------------------------------------------------------------
  
  !================================================================================
  !       Declare the variables 
  !================================================================================
  implicit none
  ! time variables
  real(dp)      :: time, time1, time2
  
  ! computation variables
  logical       :: read_status
  integer       :: iloop
  
  ! mesh variables
  integer               :: nkmesh
  
  ! consistency variable
  logical               :: k_points_consistent 
  
  !--------------------------------------------------------------------------------
  
  ! start global clock
  call get_timing(time1)
  
  !================================================================================
  !        read the input file 
  !================================================================================
  !  Read in the necessary information from standard input




  call read_input(read_status)
  if (read_status ) then
         stop 

  end if        

  !================================================================================
  !       read the global parameters from the xml files in the QE folders
  !================================================================================
  call read_parameters_data_file_xml()

  
  !================================================================================
  !      Allocate and fetch the global list of G vectors 
  !      Careful! This must be done BEFORE initializing the FFT 
  !      code because it relies on the existence of this global list.
  !================================================================================
  !Here we obtain 'ngm', the number of G vectors in the global list

  call get_ngm() 

  call get_nGmax ()   
  ! allocate and populate the global list of all G vectors
  allocate (gvec(3,ngm)) 
  call get_gvec ()  
stop  
  !================================================================================
  !      set up FFT variables 
  !================================================================================
  !allocate useful variables
stop
  call allocate_fft()

  !generate some important indices for FFT
  call generate_nl()

  
  !================================================================================
  !	read in the k vectors stored in the QE folder.
  !       CAREFUL! This could be the entire MP coarse mesh, or only
  !       the IBZ. We don't know at this point. Thus, distinguish the variables
  !       (nkpoints_QE, kpoints_QE) and (nkmesh, kmesh)
  !================================================================================
  ! allocate the kpoints array, which contains the actual k-points 
  ! present in the QE folder; this could be the IBZ! 
  allocate(kpoints_QE(3,nkpoints_QE))
  
  call read_kpoints_data_file_xml(kpoints_QE)
  
  !================================================================================
  ! Build the kmesh corresponding to the parameters in the input file 
  !================================================================================
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)
 

  !================================================================================
  !      allocate the symmetry arrays 
  !      CAREFUL! the subroutine needs to know the global value of "full_mesh",
  !      so it is crucial that this allocation occurs AFTER setting full_mesh
  !================================================================================

  call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)
  !================================================================================
  !     Compute the indices of the inverse rotation matrices.
  !     This will be useful later in the exectution.
  !================================================================================
  call find_inverse_symmetry_matrices_indices()

  !================================================================================
  !     Set up the array spin_symmetry_matrices, which contain
  !     the matrices which must be used to rotate spin
  !================================================================================
  call allocate_and_build_spin_symmetry_matrices(nsym)


  !================================================================================
  !      Fill the symmetry arrays appropriately
  !     the subroutine "set_symmetry_relations" also performs tests
  !     on the k-points to make sure they are consistent with the
  !     parameters in the input file
  !================================================================================
  
  call set_symmetry_relations(nk1, nk2, nk3, nkpoints_QE, kpoints_QE,kmesh, k_points_consistent, &
             QE_folder_sym, sym_G, symlink )
  
  ! We need to be careful here. There are many possibilities:
  !  1) A full zone calculation is stored in the QE 
  !     and kpoints_QE = kmesh
  !  2) Only the IBZ data is available, and thus kpoints_QE != kmesh.
  !     However, the IBZ points are consistent with kmesh. 
  !  3) none of the above; throw an error, what the hell is going on!?!

  if(.not. k_points_consistent .or. (.not. full_mesh .and. .not. IBZ)) then
  if (mpi_rank == 0) then
        write(*,20) '*****************************************************'
        write(*,20) '*      The kpoints present in the QE folders are    *'
        write(*,20) '*      not consistent with a full Brillouin Zone    *'
        write(*,20) '*      or an irreducible Brillouin zone!            *'
        write(*,20) '*      Review your input...                         *'
        write(*,20) '*          program continues at your own risks!     *'
        write(*,20) '*****************************************************'
        write(*,20) '  Debug information:                                 '
        write(*,*) 'k_points_consistent = ',k_points_consistent 
        write(*,*) 'full_mesh           = ',full_mesh
        write(*,*) 'IBZ                 = ',IBZ
        do iloop = 1, nkpoints_QE
                write(*,'(3F8.4)') kpoints_QE(:,iloop)
        end do
    end if
  else if (full_mesh .and. IBZ) then
       if (mpi_rank == 0) then
          write(*,20) '|      The kpoints present in the QE folder are     |'
          write(*,20) '|      consistent with a full 1BZ and contain       |'
          write(*,20) '|      an irreducible first Brillouin zone.         |'
       end if
  
  else if(.not. full_mesh .and. IBZ) then
       if (mpi_rank == 0) then
          write(*,20) '|      The kpoints present in the QE folder are     |'
          write(*,20) '|      consistent with the IBZ.                     |'
       end if
  else
       if (mpi_rank == 0) then
          write(*,*) 'nonsense: check code'
       end if
          stop
  end if
  if (mpi_rank == 0) then
    write(*,20) '====================================================='
  end if
  
  
  !================================================================================
  ! Build the list of G shells
  !================================================================================
  call get_G_shells(number_G_shells,nG_shell_max)
  
  
  ! stop global clock
  call get_timing(time2)
 
  time = time2-time1 

  20 format(A) 
  end subroutine setup_QE

  subroutine setup(time)
  !--------------------------------------------------------------------------------
  !
  !       The purpose of this subroutine is to set up the necessary information
  !       structures for computation. In particular, it will read in the 
  !       input file and generate the data stuctures necessary for the 
  !       reading and rotation of wave functions.
  !
  !--------------------------------------------------------------------------------
  
  !================================================================================
  !       Declare the variables 
  !================================================================================
  implicit none
  ! time variables
  real(dp)      :: time, time1, time2
  
  ! computation variables
  logical       :: read_status
  integer       :: iloop
  
  ! mesh variables
  integer               :: nkmesh
  
  ! consistency variable
  logical               :: k_points_consistent 

    !--------------------------------------------------------------------------------
  

  ! start global clock
  call get_timing(time1)
  
  !================================================================================
  !        read the input file 
  !================================================================================
  !  Read in the necessary information from standard input, on the master node
  if (mpi_rank == 0) then
    call read_input(read_status)
    if (read_status ) then
         stop 
    end if        
  end if

  if ( mpi_nproc > 1) then 
    call broadcast_input()
  end if
  !================================================================================
  !       read the global parameters from the xml files in the QE folders
  !================================================================================
  call read_parameters_data_file_xml()
 
 
  !================================================================================
  !      Allocate and fetch the global list of G vectors 
  !      Careful! This must be done BEFORE initializing the FFT 
  !      code because it relies on the existence of this global list.
  !================================================================================
  !Here we obtain 'ngm', the number of G vectors in the global list
  call get_ngm() 
  
  ! allocate and populate the global list of all G vectors
  allocate (gvec(3,ngm)) 
  call get_gvec ()  
  
  !================================================================================
  !      set up FFT variables 
  !================================================================================
  !allocate useful variables
  call allocate_fft()
  !generate some important indices for FFT
  call generate_nl()



  
  !================================================================================
  !	read in the k vectors stored in the QE folder.
  !       CAREFUL! This could be the entire MP coarse mesh, or only
  !       the IBZ. We don't know at this point. Thus, distinguish the variables
  !       (nkpoints_QE, kpoints_QE) and (nkmesh, kmesh)
  !================================================================================
  ! allocate the kpoints array, which contains the actual k-points 
  ! present in the QE folder; this could be the IBZ! 
  allocate(kpoints_QE(3,nkpoints_QE))
  
  call read_kpoints_data_file_xml(kpoints_QE)
  
  !================================================================================
  ! Build the kmesh corresponding to the parameters in the input file 
  !================================================================================
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)
 

  !================================================================================
  !      allocate the symmetry arrays 
  !      CAREFUL! the subroutine needs to know the global value of "full_mesh",
  !      so it is crucial that this allocation occurs AFTER setting full_mesh
  !================================================================================
  call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)

  !================================================================================
  !     Compute the indices of the inverse rotation matrices.
  !     This will be useful later in the exectution.
  !================================================================================
  call find_inverse_symmetry_matrices_indices()

  !================================================================================
  !     Set up the array spin_symmetry_matrices, which contain
  !     the matrices which must be used to rotate spin
  !================================================================================
  call allocate_and_build_spin_symmetry_matrices(nsym)

  !================================================================================
  !      Fill the symmetry arrays appropriately
  !     the subroutine "set_symmetry_relations" also performs tests
  !     on the k-points to make sure they are consistent with the
  !     parameters in the input file
  !================================================================================
  
  call set_symmetry_relations(nk1, nk2, nk3, nkpoints_QE, kpoints_QE,kmesh, k_points_consistent, &
          QE_folder_sym, sym_G, symlink )
  
  ! We need to be careful here. There are many possibilities:
  !  1) A full zone calculation is stored in the QE 
  !     and kpoints_QE = kmesh
  !  2) Only the IBZ data is available, and thus kpoints_QE != kmesh.
  !     However, the IBZ points are consistent with kmesh. 
  !  3) none of the above; throw an error, what the hell is going on!?!

  if(.not. k_points_consistent .or. (.not. full_mesh .and. .not. IBZ)) then
      if (mpi_rank == 0) then
        write(*,20) '*****************************************************'
        write(*,20) '*      The kpoints present in the QE folders are    *'
        write(*,20) '*      not consistent with a full Brillouin Zone    *'
        write(*,20) '*      or an irreducible Brillouin zone!            *'
        write(*,20) '*      Review your input...                         *'
        write(*,20) '*                   Program stops.                  *'
        write(*,20) '*****************************************************'
        write(*,20) '  Debug information:                                 '
        write(*,*) 'k_points_consistent = ',k_points_consistent 
        write(*,*) 'full_mesh           = ',full_mesh
        write(*,*) 'IBZ                 = ',IBZ
        do iloop = 1, nkpoints_QE
                write(*,'(3F8.4)') kpoints_QE(:,iloop)
        end do
      end if
        stop
  else if (full_mesh .and. IBZ) then
       if (mpi_rank == 0) then
          write(*,20) '|      The kpoints present in the QE folder are     |'
          write(*,20) '|      consistent with a full 1BZ and contain       |'
          write(*,20) '|      an irreducible first Brillouin zone.         |'
       end if
  
  else if(.not. full_mesh .and. IBZ) then
       if (mpi_rank == 0) then
          write(*,20) '|      The kpoints present in the QE folder are     |'
          write(*,20) '|      consistent with the IBZ.                     |'
       end if
  else
       if (mpi_rank == 0) then
          write(*,*) 'nonsense: check code'
       end if
          stop
  end if
  if (mpi_rank == 0) then
    write(*,20) '====================================================='
  end if
  
  
  !================================================================================
  ! Build the list of G shells
  !================================================================================
  call get_G_shells(number_G_shells,nG_shell_max)
  
  !================================================================================
  !       set up W90 data and generate the u_mesh array 
  !================================================================================
  !       CAREFUL!! It is CRUCIAL to allocate and read the W90 stuff BEFORE
  !       referencing it. num_wann, for example, exists but is not properly
  !       DEFINED before.
  
  call allocate_and_read_W90()
  !
  !       NOW you can allocate
  !
  ! extract the nband x num_wann Wannier projection+rotation matrices.
  call produce_u_mesh()
  
  
  ! stop global clock
  call get_timing(time2)
 
  time = time2-time1 

  20 format(A) 
  end subroutine setup


  subroutine cleanup() 
  !--------------------------------------------------------------------------------
  !        This subroutine deallocates and cleans up.
  !--------------------------------------------------------------------------------
  deallocate (gvec)
  deallocate (kpoints_QE)
  deallocate (kmesh)

  call deallocate_fft()
  call deallocate_reading_variables()
  call deallocate_symmetry_related_k()

  end subroutine cleanup



!------------------------------------------------------------------------
  
  subroutine write_input(log_file)
  !------------------------------------------------------------------
  ! This subroutine simply echoes the input parameters, 
  ! as read by the subroutine read_input in intw_input_parameters.
  ! 
  !------------------------------------------------------------------
 
  implicit none
 
  character(*) ::    log_file
  integer      ::    io_unit
  logical      ::    file_exists

  io_unit = find_free_unit()

  call create_or_append_to_file(io_unit,log_file)

  write(io_unit,10) "========================================="
  write(io_unit,10) "= echo of data read from input file     ="
  write(io_unit,10) "========================================="
  write(io_unit,10) "&input" 
  write(io_unit,20) "  mesh_dir    = ",trim(mesh_dir)
  write(io_unit,20) "  prefix      = ",trim(prefix)
  write(io_unit,25) "  TR_symmetry = ",TR_symmetry
  write(io_unit,30) "  nk1         = ",nk1
  write(io_unit,30) "  nk2         = ",nk2
  write(io_unit,30) "  nk3         = ",nk3
  write(io_unit,30) "  nk1s        = ",nk1s
  write(io_unit,30) "  nk2s        = ",nk2s
  write(io_unit,30) "  nk3s        = ",nk3s
  write(io_unit,30) "  G shells    = ",number_G_shells 
  write(io_unit,10) "&intw2W"
  write(io_unit,25) "  fullzone    = ",intw2W_fullzone 
  write(io_unit,20) "  method      = ",trim(intw2W_method)
  write(io_unit,10) "========================================="
  
  close(io_unit)     

  10 format(A) 
  20 format(A,A) 
  25 format(A,L5) 
  30 format(A,I4) 

  end subroutine write_input


  subroutine write_parameters_data_file_xml(log_file)
  !------------------------------------------------------------------
  ! This subroutine simply echoes the data read from the xml file,
  ! using the subroutine read_parameters_data_file_xml in intw_reading.
  !------------------------------------------------------------------
  implicit none

  character(*) ::    log_file
  integer      ::    io_unit
  logical      ::    file_exists

  integer      ::  ilp

  io_unit = find_free_unit()

  call create_or_append_to_file(io_unit,log_file)

  write(io_unit,10) "========================================="
  write(io_unit,10) "= echo of data read from data-file.xml  ="
  write(io_unit,10) "========================================="
  
  write(io_unit,10) "-----------------------------------------"
  write(io_unit,10) "---         atomic information        ---"
  write(io_unit,10) "-----------------------------------------"
  write(io_unit,30) "nat      = ",nat
  write(io_unit,30) "ntyp     = ",ntyp
  write(io_unit,10) ""
 
  write(io_unit,10) "-----------------------------------------"
  write(io_unit,10) "---               types               ---"
  write(io_unit,10) "-----------------------------------------"
  do ilp=1,ntyp
      write(io_unit,30) "type ",ilp
      write(io_unit,35)"atom_labels(",ilp,") = ",atom_labels(ilp)
      write(io_unit,35)"atom_pfile (",ilp,") = ",atom_pfile (ilp)
  end do

  write(io_unit,10) "-----------------------------------------"
  write(io_unit,10) "---               atoms               ---"
  write(io_unit,10) "-----------------------------------------"
  do ilp=1,nat
      write(io_unit,30) "atom ",ilp
      write(io_unit,40) "ityp(",ilp,") = ",ityp(ilp)
      write(io_unit,45)  "tau (",ilp,") = (",tau(:,ilp),")"
  end do
  write(io_unit,10) ""

  write(io_unit,10) "-----------------------------------------"
  write(io_unit,10) "---        Bands information          ---"
  write(io_unit,10) "-----------------------------------------"

  write(io_unit,30) "nkpoints_QE = ",nkpoints_QE
  write(io_unit,30) "ngm         = ",ngm
  write(io_unit,30) "nG_max      = ",nG_max
  write(io_unit,30) "nbands      = ",nbands
  write(io_unit,10) ""

  write(io_unit,10) "-----------------------------------------"
  write(io_unit,10) "---      symmetry information         ---"
  write(io_unit,10) "-----------------------------------------"
  do ilp=1,nsym
     write(io_unit,50) "s(",ilp,")  = (",   &
           s(1,1,ilp),s(1,2,ilp),s(1,3,ilp),")  f(",ilp," ) = (",ftau(1,ilp),")"

     write(io_unit,55) "(",                 &
           s(2,1,ilp),s(2,2,ilp),s(2,3,ilp),")           (",ftau(2,ilp),")"

     write(io_unit,55) "(",                  &
           s(3,1,ilp),s(3,2,ilp),s(3,3,ilp),")           (",ftau(3,ilp),")"
  
     write(io_unit,10) ""

  end do
  write(io_unit,10) "========================================="

  close(io_unit)


  10 format(A) 
  20 format(A,A) 
  25 format(A,L5) 
  30 format(4X,A,I4) 
  35 format(8X,A,I2,A,A) 
  40 format(8X,A,I2,A,I2) 
  45 format(8X,A,I2,A,3F8.4,A) 
  50 format(A,I2,A,3I5,A,I2,A,F8.4,A)
  55 format(9X,A,3I5,A,F8.4,A)
  end subroutine write_parameters_data_file_xml


!----------------------------------------------------------------------------!
!
!
end module intw_setup
!
!
!----------------------------------------------------------------------------!
