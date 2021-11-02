!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       program generate_Wannier_spin_overlaps
!       -----------------------------------
!
!       This is a "sub program" which is part of the intw project.
!
!       The purpose of this program is to compute the overlap of Wannier
!	functions with the purpose of projecting approximately the wavefunctions
!	onto spin components. 
!
!	This will be useful when performing a non-collinear calculation 
!	for a system which can be solved using the collinear approach.
!	The use of spinors is somewhat artificial in this case, and spin 
!	should be a good quantum number.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program generate_Wannier_spin_overlaps
  
  use intw_useful_constants
  use intw_product_basis
  use intw_reading
  use intw_setup
  use intw_symmetries
  use intw_fft
  use intw_W90
  use intw_tests


!================================================================================
!       Declare the variables 
!================================================================================
implicit none


! time variables
real(dp) :: time1, time2, time, total_time1, total_time2
real(dp) :: loop_time1, loop_time2, loop_time
real(dp) :: read_wannier_time, compute_overlap_time


logical  :: use_IBZ
logical  :: W_file_exists

! netcdf file units
integer  :: W_nc_id


!--------------------------------
! the supercell
!--------------------------------
! dimensions of the supercell
integer    :: snr1, snr2, snr3, snr
integer    :: nr


integer  ::  ipol_up, ipol_down
integer  ::  nb_W1, nb_W2, L_index
integer  ::  bare_basis_size



! block variables
integer              :: max_bands_per_block
integer              :: i_loop
integer              :: number_of_blocks, i_block, rest_block, nb_per_block
integer              :: nb_min, nb_max
integer,allocatable  :: block_bounds(:,:)


integer              :: nb1_min, nb1_max
integer              :: nb2_min, nb2_max

integer              :: i_block1, i_block2


complex(dp),allocatable    :: wannier_function_read(:,:)

complex(dp),allocatable    :: wannier_function1(:,:,:)
complex(dp),allocatable    :: wannier_function2(:,:,:)
complex(dp),allocatable    :: wannier_spin_overlaps(:,:,:)

real(dp)  :: volume
!--------------------------------------------------------------------------------


! start global clock
call get_timing(total_time1)

!================================================================================
!        Talk to user
!================================================================================
write(*,20) '====================================================='
write(*,20) '|    program  generate_Wannier_spin_overlaps        |'
write(*,20) '|    ---------------------------------------        |'
write(*,20) '====================================================='
write(*,20) '|                                                   |'
write(*,20) '|   This program computes the spin-resolved inner   |'
write(*,20) '|   product of Wannier functions and stores the     |'
write(*,20) '|   result to file. This should be useful to        |'
write(*,20) '|   identify the spin direction which can be        |'
write(*,20) '|   associated to every band in the case that spin  |'
write(*,20) '|   remains a good quantum number.                  |'
write(*,20) '====================================================='

!================================================================================
!       call the setup subroutine
!================================================================================
call setup(time)

! ham_r, from W90, is no longer needed.
call deallocate_hr_W90()

write(*,20) '====================================================='
write(*,10) '|           setup     - time::',time,' seconds    |'
write(*,20) '====================================================='

use_IBZ = .not. intw2W_fullzone

!================================================================================
!        Tell the user if the system has spin
!================================================================================

if     (nspin==1) then
	write(*,20) '|       - The calculation is paramagnetic nspin=1   |'
        npol = 1
	if (magnon) magnon = .false.
	write(*,20) '*****************************************************'
	write(*,20) '* ERROR: Using this program is pointless in this    *'
	write(*,20) '*       case. This program only applies to magnetic *'
	write(*,20) '*       systems.                                    *'
	write(*,20) '*           program stops.                          *'
	write(*,20) '*****************************************************'
	stop
elseif (nspin==2) then
	write(*,20) '|       - Non-collinear Spin calculation  nspin=2   |'
        npol = 2
else
	write(*,20) '*****************************************************'
	write(*,20) '* ERROR: Allowed values for nspin are 1 or 2        *'
	write(*,20) '*****************************************************'
	stop
endif

if (magnon) then
	write(*,20) '|       - This is a magnon calculation              |'
	write(*,20) '|         ==> nspin = 2 and npol = 1                |'
        npol = 1
endif

write(*,20) '|           ---------------------------------       |'


!================================================================================
!        Open the netcdf file containing the Wannier functions, 
!	 and extract the relevant parameters.
!================================================================================

! check that the file exists

inquire(file=wannier_functions_filename, exist=W_file_exists) 

if (.not. W_file_exists) then
  write(*,20) '*******************************************************'
  write(*,20) '*   ERROR: The file which should contain the Wannier  *'
  write(*,20) '*          functions does not appear to exist.        *'
  write(*,20) '*          Review input.                              *'
  write(*,20) '*               program stops.                        *'
  write(*,20) '*******************************************************'
  write(*,*)  ' wannier_functions_filename = ',wannier_functions_filename 
  stop
end if

! read in the relevant parameters
call open_wannier_functions_netcdf(W_nc_id,snr1,snr2,snr3,snr)


write(*,20) '====================================================='
write(*,10) '|    The file Wannier function file was found.      |'
write(*,10) '|    The parameters are:                            |'
write(*,22) '|               nr1     = ',nr1,'                      |'
write(*,22) '|               nr2     = ',nr2,'                      |'
write(*,22) '|               nr3     = ',nr3,'                      |'
write(*,22) '|              snr1     = ',snr1,'                      |'
write(*,22) '|              snr2     = ',snr2,'                      |'
write(*,22) '|              snr3     = ',snr3,'                      |'
write(*,22) '|              num_wann = ',num_wann,'                      |'
write(*,20) '====================================================='

nr = nr1*nr2*nr3


!================================================================================
! determine the size of the Wannier functions sub-blocks which will maximize
! RAM usage
!================================================================================

max_bands_per_block = int(MAX_RAM/dble(2*nspin*double_complex*snr))

i_loop = 0
do
        i_loop = i_loop + 1

        number_of_blocks = i_loop

        rest_block = modulo(num_wann,number_of_blocks)

        nb_per_block = (num_wann-rest_block)/number_of_blocks

        if (nb_per_block <= max_bands_per_block )  exit
end do


! build the blocks

allocate(block_bounds(number_of_blocks,2))

block_bounds(:,:) = 0

nb_min = 1
do i_block = 1, number_of_blocks

   block_bounds(i_block,1) = nb_min

   if ( rest_block == 0) then
        nb_max     = nb_min+nb_per_block-1
   else
        nb_max     = nb_min+nb_per_block
        rest_block = rest_block -1
   end if

   if ( nb_max < num_wann) then
        block_bounds(i_block,2) = nb_max
   else if ( nb_max >= num_wann) then
        block_bounds(i_block,2) = num_wann
        exit
   end if
   nb_min = nb_max+1
end do

! tell user about RAM management
write(*,20) '====================================================='
write(*,20) '|     Typically, the Wannier functions do not       |'
write(*,20) '|     fit in RAM, and it is time consuming to       |'
write(*,20) '|     read them over and over.                      |'
write(*,20) '|                                                   |'
write(*,20) '|     They will thus be read in blocks,             |'
write(*,20) '|     according to:                                 |'
write(*,20) '|                                                   |'
write(*,18) '|         max_ram            = ',max_ram,' bytes       |'
write(*,11) '|       number_of_blocks     = ',number_of_blocks,'             |'
write(*,11) '|         nb_per_block(+/-1) = ',nb_per_block,'             |'
write(*,20) '|    block bounds:                                  |'
write(*,20) '|    i block    nb_min             nb_max           |'
do i_block = 1, number_of_blocks
     write(*,21) i_block, block_bounds(i_block,1), block_bounds(i_block,2) 
end do
write(*,20) '|                                                   |'
write(*,20) '====================================================='


!================================================================================
!        Build the spin resolved overlaps 
!================================================================================

!-------------------------------------------------------------------------------
!
!	This program essentially computes
!
!	Wannier_overlap_{m1m2}^{sigma} = int dr  W_{m1}^{sigma,*}(r) W_{m2}^{sigma}(r) 
!
!	This integral allows us to compute the overlap of the Bloch states,
!
!			< Psi_{nk} | P | Psi_{nk} > = [ | psi_nk^up|^2 ]
!			                              [ | psi_nk^dn|^2 ]
!
!	where Psi_{nk} = [  psi_{nk}^up ]
!	                 [  psi_{nk}^dn ]
!
!		and P = [ 1 0 ]        or P = [ 0 0 ]
!		        [ 0 0 ]               [ 0 1 ]
!
!-------------------------------------------------------------------------------
L_index   = 0

ipol_up   = 1
ipol_down = 2


allocate(wannier_spin_overlaps(num_wann,num_wann,nspin))


volume = alat**3*abs(                                  &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))


write(*,20) '====================================================='
write(*,20) '|    Computing the spin resolved overlaps :         |'
write(*,20) '====================================================='

bare_basis_size = num_wann**2

read_wannier_time    = ZERO
compute_overlap_time = ZERO

call get_timing(loop_time1)

allocate(wannier_function_read(snr,nspin))

! loop over the blocks

do i_block1 = 1,number_of_blocks

     ! read in the first block
     nb1_min = block_bounds(i_block1,1)
     nb1_max = block_bounds(i_block1,2)

     allocate(wannier_function1(snr,nb1_max-nb1_min+1,nspin))

     call get_timing(time1)
     do nb_W1 = nb1_min,nb1_max
        call read_wannier_functions_netcdf  &
                (W_nc_id,nb_W1,snr,wannier_function_read)

         wannier_function1(:,nb_W1-nb1_min+1,:) = wannier_function_read(:,:)
	
     end do 

     call get_timing(time2)
     read_wannier_time = read_wannier_time + time2-time1


     do i_block2 = 1,number_of_blocks

        L_index = L_index + 1
        write(*,16) '| block     ',L_index,' of ',number_of_blocks**2,'                        |'


     	! read in the second block
     	nb2_min = block_bounds(i_block2,1)
     	nb2_max = block_bounds(i_block2,2)

     	allocate(wannier_function2(snr,nb2_max-nb2_min+1,nspin))

     	call get_timing(time1)
     	do nb_W2 = nb2_min,nb2_max
        	call read_wannier_functions_netcdf  &
                	(W_nc_id,nb_W2,snr,wannier_function_read)

         	wannier_function2(:,nb_W2-nb2_min+1,:) = wannier_function_read(:,:)
	
     	end do 

     	call get_timing(time2)
     	read_wannier_time = read_wannier_time + time2-time1

        ! compute product of Wannier functions

        call get_timing(time1)
     	do nb_W1 = nb1_min,nb1_max
     	   do nb_W2 = nb2_min,nb2_max

		wannier_spin_overlaps(nb_W1,nb_W2,ipol_up) =                 &
     		  sum(conjg(wannier_function1(:,nb_W1-nb1_min+1,ipol_up))*   &
     	                    wannier_function2(:,nb_W2-nb2_min+1,ipol_up))

		wannier_spin_overlaps(nb_W1,nb_W2,ipol_down) =               &
     		  sum(conjg(wannier_function1(:,nb_W1-nb1_min+1,ipol_down))* &
     	                    wannier_function2(:,nb_W2-nb2_min+1,ipol_down))
           end do ! nb_W2 
        end do ! nb_W1

        call get_timing(time2)

        compute_overlap_time = compute_overlap_time + time2-time1

     	deallocate(wannier_function2)
     end do ! i_block2

     deallocate(wannier_function1)
end do !i_block1

call get_timing(loop_time2)
loop_time = loop_time2-loop_time1

!close file
call close_W_functions_netcdf(W_nc_id)

! normalize at the end
wannier_spin_overlaps(:,:,:) = wannier_spin_overlaps(:,:,:)*volume/nr 

write(*,20) '====================================================='
write(*,18) '|  - reading Wannier       : ',read_wannier_time, ' sec.          |'
write(*,18) '|  - computing overlap     : ',compute_overlap_time, ' sec.          |'
write(*,18) '|  - total loop time       : ',loop_time,' sec.          |'
write(*,20) '====================================================='



write(*,20) '====================================================='
write(*,20) '|    write overlaps to netcdf file...               |'
call get_timing(time1)
call output_spin_resolved_Wannier_overlaps_netcdf &
                                (num_wann,wannier_spin_overlaps)
call get_timing(time2)

write(*,10) '|                     - time::',time,' seconds    |'
write(*,20) '====================================================='





!================================================================================
!        clean up and finish
!================================================================================

call cleanup()

! good bye!
call get_timing(total_time2)

write(*,30) '|        DONE - total time::',     &
                                total_time2-total_time1,' seconds      |'
write(*,20) '====================================================='

10 format(A,F10.1,A)
11 format(A,I8,A)
16 format(A,I6,A,I6,A)
18 format(A,ES8.2,A)
20 format(A)
21 format(5X,I6,I12,10X,I12)
22 format(A,I4,A)
23 format(A,I6,A)
30 format(A,F10.3,A)

end program generate_Wannier_spin_overlaps
