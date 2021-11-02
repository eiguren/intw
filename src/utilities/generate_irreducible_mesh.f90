!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       program generate_irreducible_mesh 
!       -----------------------------------
!
!       This is a "utility" program which is part of the intw project.
!
!       The purpose of this utility is to generate a set of irreducible
!       k points in crystal coordinates which can then be used to peform
!       a Quantum Espresso non-self consistent calculation. 
!
!       This program will read the same input file as the main program
!       intw.x, and will deduce the irreducible k-points. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program generate_irreducible_mesh 

 
  use intw_useful_constants
  use intw_input_parameters
  use intw_tests
  use intw_symmetries
  use intw_tests

!================================================================================
!       Declare the variables 
!================================================================================
implicit none

integer       ::      ikpt
integer       ::      i, j, k 

integer       ::      io_unit

logical       ::      read_status

integer       ::      nk_irr , nkmesh

real(dp),allocatable:: kpoints_irr(:,:), kmesh(:,:)



!--------------------------------------------------------------------------------

!================================================================================
!       read the input file 
!       Read in the necessary information from standard input
!================================================================================
call read_input(read_status)

if (read_status ) then
       stop 
end if        

!================================================================================
!       read the parameters from the SCF QE calculation, in  
!       particular, read in the symmetry matrices!!! 
!================================================================================
call read_parameters_data_file_xml()

!================================================================================
!       Talk to the user
!================================================================================

write(*,*) '====================================================='
write(*,*) '|        program generate_irreducible_mesh          |'
write(*,*) '|        ---------------------------------          |'
write(*,*) '====================================================='
write(*,*) '|        This program is part of the intw           |'
write(*,*) '|        project. It computes the irreducible       |'
write(*,*) '|        k-points (the k-points in the wedge)       |'
write(*,*) '|        at which the wavefunctions are needed.     |'
write(*,*) '|                                                   |'
write(*,*) '|        The ideal usage of this program is to      |'
write(*,*) '|        1) Perform a scf QE calculation            |'
write(*,*) '|        2) Produce a intw input file               |'
write(*,*) '|        3) use generate_irreducible_mesh.x         |'
write(*,*) '|           to produce the appropriate k-points     |'
write(*,*) '|           for the nscf QE calculation.            |'
write(*,*) '|                                                   |'
write(*,*) '|        The k-points are stored in the file        |'
write(*,*) '|                   irr_mesh.dat                    |'
write(*,*) '|        The symmetry matrices used are found in    |'
write(*,*) '|   	test_point_group_symmetries.test        |'
write(*,*) '====================================================='
          
 

!================================================================================
!	generate the k vectors for the coarse mesh	
!       according to the parameters read in
!================================================================================

call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)

nkmesh = nk1*nk2*nk3
allocate(kmesh(3,nkmesh))
allocate(kpoints_irr(3,nkmesh))
call generate_kmesh(kmesh,nk1,nk2,nk3)

call find_the_irreducible_k_set(nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)


!================================================================================
!       print out the irreducible k-points, in crystal units 
!       in a format suitable for Quantum Espresso
!================================================================================
io_unit = find_free_unit()
open(unit=io_unit,file='irr_mesh.dat',status='unknown')

write(io_unit,*) 'K_POINTS crystal'
write(io_unit,50) nk_irr 

do ikpt = 1, nk_irr
        write(io_unit,100) kpoints_irr(:,ikpt),1.0
end do



call test_point_group_symmetries()
!================================================================================
!       clean up 
!================================================================================
close(io_unit)
call deallocate_symmetry_related_k()
deallocate(kmesh)
deallocate(kpoints_irr)

50  format(I6)
100 format(4x,3F16.12,2x,F4.2)

end program generate_irreducible_mesh 
