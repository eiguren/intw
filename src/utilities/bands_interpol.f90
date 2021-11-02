!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       program bands_interpol
!       -----------------------------------
!
!       This is a "sub program" which is part of the intw project.
!
!       The purpose of this program is to interpolate the bands using Wannier  
!       interpolation. The Wannier90 code already has this functionality,  
!       but this sub-program will also implement routines which can resolve
!       band crossing/anti-crossings. This will be important when integrating
!       using the "improved" tetrahedron method.                                  
!
!       Furthermore, this sub-program will interpolate the W_rotated matrix
!       elements on a path: this will be useful for testing.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program bands_interpol 


  use intw_useful_constants
  use intw_setup

  use intw_input_parameters
  use intw_reading
  use intw_symmetries
  use intw_W90
  use intw_tests
  use intw_plane_wave_setup
  use intw_new_interpolate 
  use intw_band_crossing
!  use intw_product_basis
!  use intw_chi_0


!================================================================================
!       Declare the variables 
!================================================================================
implicit none

! keeping track of time!
real(dp)       ::      time, time1, time2
real(dp)       ::      total_time1, total_time2


integer        ::      io_unit  ,io_unit1, io_unit2  

logical        ::      use_IBZ

integer        ::      nk_vec_path

integer        ::      iloop,  ik_loop  

integer        ::      nb


integer        ::      number_of_special_points


real(dp),allocatable     :: list_k_path(:,:)
real(dp),allocatable     :: list_x_path(:)
real(dp),allocatable     :: list_xticks(:)

real(dp),allocatable     :: eig_int_k(:,:)
complex(dp),allocatable  :: u_int_k(:,:,:)
complex(dp),allocatable  :: H_int_k(:,:,:)

complex(dp),allocatable  :: wannier_spin_overlaps(:,:,:)


real(dp),allocatable     :: psi2_k(:,:,:)
complex(dp),allocatable  :: U(:,:), Ud(:,:)
complex(dp),allocatable  :: Mup(:,:), Mdn(:,:)
integer                  :: i,ipol_up, ipol_down

call get_timing(total_time1)
!================================================================================
!        Talk to user
!================================================================================


use_IBZ = .true.

write(*,20) '====================================================='
write(*,20) '|         program bands_interpol                    |'
write(*,20) '|  ---------------------------------                |'
write(*,20) '====================================================='
write(*,20) '|                                                   |'
write(*,20) '|  This sub-program computes the interpolated bands |'
write(*,20) '|  using Wannier interpolation. Its main purpose is |'
write(*,20) '|  to test the quality of the MLWFs.                |'
write(*,20) '====================================================='

!================================================================================
!        call "setup" routine, which sets up the necessary environment.
!================================================================================

print*, "kaix1"
call setup(time)
print*, "kaix2"



write(*,20) '====================================================='
write(*,10) '|            setup    - time::',time,' seconds    |'
write(*,20) '====================================================='

!================================================================================
!       Get the special points from the file
!================================================================================

if (bands_points /= 0 ) then
        call get_interpolation_special_points( number_of_special_points)
end if


!================================================================================
!       Define the interpolation path     
!================================================================================

nk_vec_path = bands_points

allocate(list_k_path(3,nk_vec_path))
allocate(list_x_path(nk_vec_path))
allocate(list_xticks(number_of_special_points))

call build_list_kpath(number_of_special_points,nk_vec_path,    &
                      list_k_path,list_x_path, list_xticks)

!================================================================================
!      Establish if the matrix elements must be interpolated.
!      Allocate necessary arrays
!================================================================================

write(*,20) '|      Method for disentanglement of bands:         |'
write(*,21) '|                  ',crossing_method,'              |'
write(*,20) '====================================================='

!================================================================================
!       First, Compute the interpolated eigenvalues and U matrices on path       
!================================================================================
allocate(u_int_k    (num_wann,num_wann,nk_vec_path))
allocate(H_int_k    (num_wann,num_wann,nk_vec_path))
allocate(eig_int_k  (num_wann,nk_vec_path))


call get_timing(time1)

call interpolate_u_matrices_and_eigenvalues(u_int_k,eig_int_k,list_k_path,    &
                                            H_int_k,nk_vec_path,nk_vec_path)

call get_timing(time2)
write(*,15)    '           Wannier time : ',time2-time1,' seconds'

! ham_r, from W90, is no longer needed.
call deallocate_hr_W90()

!================================================================================
!       Output bands 
!================================================================================

io_unit = find_free_unit()
open(unit=io_unit,file=trim(bands_file),status='unknown')

write(io_unit,20) '########################################################'
write(io_unit,20) '# Bands interpolated using Wannier interpolation        '
write(io_unit,20) '# ----------------------------------------------        '
write(io_unit,20) '# The first column corresponds to | k |                 '
write(io_unit,20) '# subsequent columns are the band eigenvalues, in eV.   '
write(io_unit,20) '# '
write(io_unit,20) '# Ticks:'

do iloop = 1, number_of_special_points
        write(io_unit,'(A,A,5X,F20.16)') '#! ',interpolation_labels(iloop),list_xticks(iloop)
end do
write(io_unit,20) '########################################################'



if (crossing_method == 'simple') then
   ! order the bands

  call resolve_band_ordering(nk_vec_path,num_wann,    &
                        list_k_path,H_int_k,eig_int_k,u_int_k)

end if


! output the bands
do ik_loop = 1,nk_vec_path
   write(io_unit,50) list_x_path(ik_loop), eig_int_k(:,ik_loop)
end do 


close(unit=io_unit)

!================================================================================
! output spin projections, if need be
!================================================================================


if (spin_projection) then

	io_unit1 = find_free_unit()
	open(unit=io_unit1,file=trim('projections_up.dat'),status='unknown')

	io_unit2 = find_free_unit()
	open(unit=io_unit2,file=trim('projections_down.dat'),status='unknown')

	write(io_unit1,20) '########################################################'
	write(io_unit1,20) '# Approximate projection of bands onto spin, using      '
	write(io_unit1,20) '# Wannier functions.                                    '
	write(io_unit1,20) '# ----------------------------------------------        '
	write(io_unit1,20) '# The first column corresponds to | k |                 '
	write(io_unit1,20) '# subsequent columns are the |psi_{nk}^up|^2, in the    '
	write(io_unit1,20) '# same order as the bands.                              '
	write(io_unit1,20) '# '
	write(io_unit1,20) '########################################################'

	write(io_unit2,20) '########################################################'
	write(io_unit2,20) '# Approximate projection of bands onto spin, using      '
	write(io_unit2,20) '# Wannier functions.                                    '
	write(io_unit2,20) '# ----------------------------------------------        '
	write(io_unit2,20) '# The first column corresponds to | k |                 '
	write(io_unit2,20) '# subsequent columns are the |psi_{nk}^dn|^2, in the    '
	write(io_unit2,20) '# same order as the bands.                              '
	write(io_unit2,20) '# '
	write(io_unit2,20) '########################################################'


	! output the projections
	do ik_loop = 1,nk_vec_path
	   write(io_unit1,50) list_x_path(ik_loop), psi2_k(:,ipol_up,ik_loop)
	   write(io_unit2,50) list_x_path(ik_loop), psi2_k(:,ipol_down,ik_loop)
	end do 

	close(unit=io_unit1)
	close(unit=io_unit2)
end if
!================================================================================
!       finish up
!================================================================================

deallocate(list_k_path)
deallocate(list_x_path)
deallocate(list_xticks)

call cleanup()

call get_timing(total_time2)
time = total_time2-total_time1

write(*,20) '====================================================='
write(*,10) '|       DONE     -total time::',time,' seconds    |'
write(*,20) '====================================================='

10 format(A,F10.1,A)
15 format(A,es12.4,A)
16 format(A,I6)
17 format(15X,A,ES12.3)
20 format(A)
21 format(A,A,A)
24 format(A,F8.4,A,F8.4,A,F8.4,A)
25 format(A,F8.4,A,F8.4,A,F8.4,A,3X,A,F8.4,A,F8.4,A,F8.4,A)
26 format(I4,4X,F8.4,4X,F8.4,4X,F8.4,4X,F8.4,4X,F8.4,4X,F8.4,4X,F8.4)
27 format(I8,4X,100I8)
30 format(A,I4,A)
40 format(A,I4,A,I4)
50 format(1000F20.12)
55 format(1000F8.4)
56 format(1000I4)
60 format(1000F8.4)
61 format(1000(2F8.4,4x))

end program bands_interpol 
