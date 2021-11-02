!       program generate_optimal_basis
!       -----------------------------------
!
!       This is a "sub program" which is part of the intw project.
!
!       The purpose of this program is to generate the "optimal basis", namely
!	the combination of bare product basis functions {|B>} which generate
!	an orthonormal basis {|b>} which :
!		- spans the functional space {1, mz, Delta, [other functions]},
!		- is generated using the QR decomposition, which allows 
!		  easy convergence testing with basis size.
!
!	This program is mainly useful for magnon computations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program generate_optimal_basis
  
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

character(10000) :: header_string, line_string,G_string

logical  :: use_IBZ
logical  :: mz_file_exists, Delta_file_exists

integer  :: ir, nr
integer  :: ir1, ir2, ir3
integer  :: iG
integer  :: switch_singlet_to_triplet
integer  :: number_of_symmetry_G_shells
integer  :: ishell
integer  :: bare_basis_size, ortho_basis_size, special_functions_basis_size
integer  :: L_index, I_index
integer  :: info, lwork,lrwork
integer  :: pw_dim
integer  :: B_nc_id
integer  :: nb_W1, nb_W2
integer  :: io_unit

integer,allocatable  :: list_m1(:), list_m2(:)
integer,allocatable  :: shells_table(:,:)


real(dp) :: time1, time2, time, total_time1, total_time2
real(dp) :: qvec(3), rvec(3)
real(dp) :: max_Lambda
real(dp) :: real_norm
real(dp) :: q_direction(3)

real(dp),allocatable    :: eta_array(:)
real(dp),allocatable    :: Lambda(:)

complex(dp)  :: eta_normalization
complex(dp)  :: norm

complex(dp),allocatable :: overlaps(:,:)
complex(dp),allocatable :: overlaps_b(:,:)
complex(dp),allocatable :: Vo(:,:)
complex(dp),allocatable :: Pc(:,:)
complex(dp),allocatable :: PT(:,:)
complex(dp),allocatable :: work(:)
complex(dp),allocatable :: T_functional_projectors(:,:)
complex(dp),allocatable :: field_r(:)
complex(dp),allocatable :: bare_basis_kernel(:,:)
complex(dp),allocatable :: C_basis_kernel(:,:)
complex(dp),allocatable :: b_basis_kernel(:,:)
complex(dp),allocatable :: one_r(:)
complex(dp),allocatable :: mz_r(:)
complex(dp),allocatable :: Delta_r(:)
complex(dp),allocatable :: plane_waves(:,:)
complex(dp),allocatable :: special_functions_B_overlaps(:,:)
complex(dp),allocatable :: A_QR_decomposition(:,:)
complex(dp),allocatable :: tau_QR_decomposition(:)
complex(dp),allocatable :: special_functions_C_coefficients(:,:)
complex(dp),allocatable :: special_functions_b_coefficients(:,:)
complex(dp),allocatable :: b_basis_C_coefficients(:,:)
complex(dp),allocatable :: B_function(:)
complex(dp),allocatable :: g_functions_B_coefficients(:,:)


!--------------------------------------------------------------------------------


! start global clock
call get_timing(total_time1)

!================================================================================
!        Talk to user
!================================================================================
write(*,20) '====================================================='
write(*,20) '|          program generate_optimal_basis           |'
write(*,20) '|        ---------------------------------          |'
write(*,20) '====================================================='
write(*,20) '|                                                   |'
write(*,20) '|        This program computes the optimal basis    |'
write(*,20) '|        for magnon calculations.                   |'
write(*,20) '====================================================='

!================================================================================
!       call the setup subroutine
!================================================================================
call setup(time)
switch_singlet_to_triplet = -1
! ham_r, from W90, is no longer needed.
call deallocate_hr_W90()

write(*,20) '====================================================='
write(*,10) '|           setup     - time::',time,' seconds    |'
write(*,20) '====================================================='

use_IBZ = .not. intw2W_fullzone

!================================================================================
!        Tell the user if the system has spin
!================================================================================

if     (nspin==1 .or. .not. magnon) then
	write(*,20) '*****************************************************'
	write(*,20) '*   The calculation is paramagnetic nspin=1         *'
	write(*,20) '* ERROR: This program is for magnons only!          *'
	write(*,20) '*****************************************************'
	stop
elseif (nspin==2) then
	write(*,20) '|       - Non-collinear Spin calculation  nspin=2   |'
        npol = 2
else
	write(*,20) '*****************************************************'
	write(*,20) '* ERROR: This program is for magnons only!          *'
	write(*,20) '*****************************************************'
	stop
endif

if (magnon) then
	write(*,20) '|       - This is a magnon calculation              |'
	write(*,20) '|         ==> nspin = 2 and npol = 1                |'
        npol = 1
endif

write(*,20) '====================================================='


!-------------------------------------------------------------------------------
! Some explanations are in order, or else all these basis functions get
! confusing.
!	Define:
!
!                   | B_L >  : the bare "product basis"; functions are linearly
!			       dependent. There are bare_basis_size = num_wann**2 
!				of them.
!
! O_{L1L2} = < B_L1 | B_L2 > :  overlap of bare product basis functions
!
!				O = Vo * Lambda * Vo^dagger,   Vo unitary
!
!				Pc_{LI} = Vo_{LI}/sqrt{Lambda_{II}}, Lambda_{II} > tolerance
!
!		    | C_I >  : temporary orthogonal basis functions, built 
!				by orthogonalizing the overlaps O and throwing
!				out combinations with eigenvalues smaller than
!				a given tolerance.
!
! 				| C_I > = sum_L | B_L > Pc_{LI} 
!
!		     | g_i > : special functions necessary to impose the acoustic
!			       condition.
!
!				i       function
!                               -------------------
!				1        one
!				2        mz(r)
!				3        Delta(r)
!				...	plane waves-like functions
!
!				| g'_i > = sum_I | C_I > g_{Ii}
!
!				g_{Ii}  = <C_I | g_i > 
!				        = sum_{L} Pc^dagger_{IL} <B_L | G_i > 
!
!				The difference between | g'_i > and |g_i > 
!				indicates the shortcomings of the basis in describing
!				the special functions.
!
!		     | b_i > : orthonormal basis functions obtained by QR 
!				decomposition of the g_{Ii} matrix.
!
!				This basis is orthogonal and spans a subspace
!				containing the special functions as much as possible.
!
!-------------------------------------------------------------------------------



! We are neglecting the dependence in R of the basis function;
! that could be improved in the future.
bare_basis_size = num_wann**2


!================================================================================
!        Read the bare basis XC kernel and the overlap matrix.
!================================================================================

allocate(overlaps(bare_basis_size,bare_basis_size ))
allocate(bare_basis_kernel(bare_basis_size,bare_basis_size ))

call read_overlaps_netcdf(bare_basis_size,overlaps,bare_basis_kernel,qvec)


!================================================================================
!        Diagonalize the overlap and throw away eigenvectors corresponding
!	 to small eigenvalues in order to extract a temporary
!	 orthogonal basis set.
!================================================================================

allocate(Vo(bare_basis_size,bare_basis_size))
allocate(Lambda(bare_basis_size))

call diagonalize_matrix(bare_basis_size,overlaps,Vo,Lambda)



max_Lambda = maxval(Lambda(:))

! Find the size of the temporary basis
ortho_basis_size = 0

do L_index = 1, bare_basis_size
	if ( Lambda(L_index)/max_Lambda  >= tolerance ) then
		ortho_basis_size = ortho_basis_size + 1
	end if
end do

! Build the coefficients of the C basis 

allocate(Pc(bare_basis_size,ortho_basis_size))

I_index = 0
do L_index = 1,bare_basis_size
	if ( Lambda(L_index)/max_Lambda  >= tolerance ) then
		I_index = I_index + 1
		Pc(:,I_index) = &
			Vo(:,L_index)/sqrt(Lambda(L_index))
	end if
end do

! get rid of arrays no longer needed
deallocate(Vo)
deallocate(Lambda)


write(*,20) '|           For the given tolerance, the            |'
write(*,22) '|           orthonormal basis has dimension',ortho_basis_size ,'     |'
write(*,20) '====================================================='

!================================================================================
!        Read the ascii files containing the ground state densities. 
!================================================================================

! check that the file exists

inquire(file=mz_r_filename,    exist=mz_file_exists) 
inquire(file=Delta_r_filename, exist=Delta_file_exists) 

if (.not. mz_file_exists .or. .not. Delta_file_exists ) then
  write(*,20) '*******************************************************'
  write(*,20) '*   ERROR: The file which should contain the ground   *'
  write(*,20) '*          state magnetization and/or the difference  *'
  write(*,20) '*          of XC potentials do/does not exist!        *'
  write(*,20) '*               program stops.                        *'
  write(*,20) '*******************************************************'
  stop
else
	write(*,20) '|       - The files containing the mz(r), Delta(r)  |'
	write(*,20) '|         have been found.                          |'
	write(*,20) '|           ---------------------------------       |'
end if


nr = nr1*nr2*nr3

allocate(one_r(nr))
allocate(mz_r(nr))
allocate(Delta_r(nr))

allocate(field_r(nr))


one_r(:) = cmplx_1

! read in the data from the files
! CAREFUL!!!
!
!	The data in the ascii files was generated using a hack from
!	Quantum Espresso; it is not very clear, but the most likely
!	order for the r-mesh is Fortran style (ie, the first dimension loops
!	fastest, "column-major").

!	HOWEVER, in this code the real space functions are stored using C style,
!	where the last dimension loops fastest "row-major". There is no 
!	good reason for this (I just started doing it this way, I can't remember why);
!	it is simpler to account for it here than to modify everything else in the code.

call read_field_ascii_file(mz_r_filename,nr,field_r)
call change_order(nr,field_r,mz_r)

call read_field_ascii_file(Delta_r_filename,nr,field_r)
call change_order(nr,field_r,Delta_r)


deallocate(field_r)

!================================================================================
!	Get the integrals of the relevant functions with the 
!	bare B basis functions.
!================================================================================

if (number_G_shells >1) then
   call get_G_shells(number_G_shells,nG_shell_max)


   if (trim(plane_wave_combination) == 'symmetric'  .or.    &
       trim(plane_wave_combination) == 'mz_times_symmetric') then

      allocate(shells_table(nG_shell_max,nsym))

      q_direction(:) = (/qpt_dir_1,qpt_dir_2,qpt_dir_3/)
      call sort_G_vectors_in_symmetry_shells(nG_shell_max,q_direction,shells_table)

      ! output the results to an ascii file, for checking!
      call output_sorted_G_vectors_in_symmetry_shells(nG_shell_max,q_direction,shells_table)

      number_of_symmetry_G_shells = 0
      do I_index = 1,nG_shell_max
	 if (shells_table(I_index,1) == 0) exit
	
         number_of_symmetry_G_shells = number_of_symmetry_G_shells +1
      end do

      pw_dim = number_of_symmetry_G_shells -1

      write(*,20) '|       - The input file specifies that the basis   |'
      write(*,20) '|         should contain symmetric combinations of  |'
      write(*,20) '|         plane waves. The number of such functions |'
      write(*,20) '|         excluding "1" for the specified parameter |'
      write(*,20) '|         number_G_shells is                        |'
      write(*,22) '|            ',pw_dim,'.                                  |'
      write(*,20) '|           ---------------------------------       |'

      optimal_basis_size = number_of_symmetry_G_shells+2
      write(*,20) '|       - The parameter optimal_basis_size has been |'
      write(*,20) '|         set to                                    |'
      write(*,22) '|            ',optimal_basis_size,'.                                  |'
      write(*,20) '|         WARNING: it is now obsolete to specify    |'
      write(*,20) '|                  this number in the input file.   |'
      write(*,20) '|           ---------------------------------       |'


      if (optimal_basis_size > ortho_basis_size ) then
	   ! the basis size needs to be no larger than the number of available
           ! orthogonal basis vectors!
  	   write(*,20) '*******************************************************'
  	   write(*,20) '*   ERROR: The optimal_basis_size parameter is larger *'
           write(*,20) '*          than the number of available orthonormal   *'
           write(*,20) '*          basis vectors! This cannot be accomodated. *'
           write(*,20) '*               program stops.                        *'
           write(*,20) '*******************************************************'
           stop
      end if

      write(*,20) '|       - Generating plane waves on real space      |'
      write(*,20) '|         grid ...                                  |'

      call get_timing(time1)
      allocate(plane_waves(nr,pw_dim))

      plane_waves(:,:) = cmplx_0

      do ir = 1, nr
	   call switch_indices(nr1,nr2,nr3,ir,ir1,ir2,ir3,switch_singlet_to_triplet)

	   rvec(1) = dble(ir1-1)/dble(nr1)
	   rvec(2) = dble(ir2-1)/dble(nr2)
	   rvec(3) = dble(ir3-1)/dble(nr3)

           ! skip the first shell, which is the function 1.
           do ishell = 2, number_of_symmetry_G_shells 

              do I_index =1, nsym
	         iG = shells_table(ishell,I_index)

		 if (iG == 0) exit

		 plane_waves(ir,ishell-1) =                   &
		 		plane_waves(ir,ishell-1) +    &
				exp(cmplx_i*tpi*dot_product(rvec(:),gvec(:,iG)))

	      end do ! I_index


           end do ! ishell 
	
      end do ! ir

      ! multiply the plane waves by mz(r) if requested by user
      if (trim(plane_wave_combination) == 'mz_times_symmetric') then
        do ishell = 2, number_of_symmetry_G_shells 
	   plane_waves(:,ishell-1) = mz_r(:)*plane_waves(:,ishell-1)
        end do
      end if

      ! normalize
      do ishell = 1, pw_dim
	 norm = sqrt(abs(sum(conjg(plane_waves(:,ishell))*plane_waves(:,ishell))/nr))
         plane_waves(:,ishell) = plane_waves(:,ishell)/norm
      end do


      call get_timing(time2)
      write(*,30) '|                  time    :',time2-time1,' seconds      |'
      write(*,20) '|           ---------------------------------       |'


   else if (trim(plane_wave_combination) == 'all') then


        pw_dim = nG_shell_max-1

	write(*,20) '|       - The input file specifies that the basis   |'
	write(*,20) '|         should contain plane waves. For the       |'
	write(*,20) '|         specified parameter number_G_shells,      |'
	write(*,20) '|         the number of plane waves not equal to 1  |'
	write(*,22) '|         is ',pw_dim,'.                                  |'
	write(*,20) '|           ---------------------------------       |'


        optimal_basis_size = pw_dim+3
        write(*,20) '|       - The parameter optimal_basis_size has been |'
        write(*,20) '|         set to                                    |'
        write(*,22) '|            ',optimal_basis_size,'.                                  |'
        write(*,20) '|         WARNING: it is now obsolete to specify    |'
        write(*,20) '|                  this number in the input file.   |'
        write(*,20) '|           ---------------------------------       |'

        if (optimal_basis_size > ortho_basis_size ) then
	   ! the basis size needs to be no larger than the number of available
           ! orthogonal basis vectors!
  	   write(*,20) '*******************************************************'
  	   write(*,20) '*   ERROR: The optimal_basis_size parameter is larger *'
           write(*,20) '*          than the number of available orthonormal   *'
           write(*,20) '*          basis vectors! This cannot be accomodated. *'
           write(*,20) '*               program stops.                        *'
           write(*,20) '*******************************************************'
           stop
        end if


	write(*,20) '|       - Generating plane waves on real space      |'
	write(*,20) '|         grid ...                                  |'

	call get_timing(time1)
	allocate(plane_waves(nr,pw_dim))

	do ir = 1, nr
	   call switch_indices(nr1,nr2,nr3,ir,ir1,ir2,ir3,switch_singlet_to_triplet)

	   rvec(1) = dble(ir1-1)/dble(nr1)
	   rvec(2) = dble(ir2-1)/dble(nr2)
	   rvec(3) = dble(ir3-1)/dble(nr3)

	   do iG = 2, nG_shell_max
		plane_waves(ir,iG-1) = exp(cmplx_i*tpi*dot_product(rvec(:),gvec(:,iG)))
	   end do !iG
	
	end do ! ir
	call get_timing(time2)
	write(*,30) '|                  time    :',time2-time1,' seconds      |'
	write(*,20) '|           ---------------------------------       |'

   end if

end if 

!================================================================================
!	Compute the coefficients of the special functions in the C basis,
!	g_q^{Ii}....
!================================================================================


write(*,20) '|       - computing overlaps with B functions...    |'
call get_timing(time1)

! order the special functions as {1, mz, Delta, [plane waves]}.
special_functions_basis_size = 3+pw_dim 

allocate(special_functions_B_overlaps(bare_basis_size,special_functions_basis_size ))

allocate(B_function(nr))

! open to read (but not to write!)
call open_B_functions_netcdf(B_nc_id)

do L_index = 1, bare_basis_size

        call read_B_function_netcdf(B_nc_id,L_index,B_function,nr)


	special_functions_B_overlaps(L_index,1) =    &
				sum(conjg(B_function(:))*one_r(:))/dble(nr)

	special_functions_B_overlaps(L_index,2) =    &
				sum(conjg(B_function(:))*mz_r(:))/dble(nr)

	special_functions_B_overlaps(L_index,3) =    &
				sum(conjg(B_function(:))*Delta_r(:))/dble(nr)

	if (nG_shell_max > 1) then
	   do iG = 1,pw_dim
		special_functions_B_overlaps(L_index,3+iG) =    &
				sum(conjg(B_function(:))*plane_waves(:,iG))/dble(nr)
           end do
	end if

end do ! L_index
call close_B_functions_netcdf(B_nc_id)
deallocate(B_function)


if (nG_shell_max > 1) deallocate(plane_waves)


call get_timing(time2)
write(*,30) '|                  time    :',time2-time1,' seconds      |'
write(*,20) '|           ---------------------------------       |'




! Express the special functions in the orthogonal | C_I >  basis

allocate(special_functions_C_coefficients(ortho_basis_size,special_functions_basis_size))

special_functions_C_coefficients(:,:) =             &
		matmul(transpose(conjg(Pc(:,:))),special_functions_B_overlaps(:,:))


!================================================================================
!        Extract the coefficients of the optimal basis using QR decomposition.
!================================================================================

write(*,20) '|       - Extracting optimal basis from QR          |'
write(*,20) '|         decomposition...                          |'
call get_timing(time1)


lwork = special_functions_basis_size
allocate(work(lwork))

allocate(A_QR_decomposition(ortho_basis_size,special_functions_basis_size))
allocate(tau_QR_decomposition(ortho_basis_size))

A_QR_decomposition(:,:) = special_functions_C_coefficients(:,:)


! A = QR decomposition
! The A matrix is the matrix special_functions_C_coefficients
! zgeqrf(m, n, a, lda, tau, work, lwork, info)
call zgeqrf(        ortho_basis_size,        & ! number of rows of the matrix
	special_functions_basis_size,        & ! number of columns of the matrix
                  A_QR_decomposition,        & ! the A matrix
  		    ortho_basis_size,        & ! the first dimension of A; modified on output
                tau_QR_decomposition,        & ! contains information about Q matrix
			work, lwork, info)     ! extra info

if     (info /=0 ) then
	write(*,20) '*****************************************************'
	write(*,20) '* ERROR: QR decomposition failed! Something must be *'
	write(*,20) '*        very wrong... program stops!               *'
	write(*,20) '*****************************************************'
	stop
end if

! generate the Q matrix explicitely
!call zungqr(m, p, p, a, lda, tau, work, lwork, info)

call zungqr(        ortho_basis_size,        & ! number of rows of the Q matrix
	special_functions_basis_size,        & ! number of required columns of the Q matrix
	special_functions_basis_size,        & ! number of required columns of the Q matrix
                  A_QR_decomposition,        & ! the A matrix
  		    ortho_basis_size,        & ! the first dimension of A; modified on output
                tau_QR_decomposition,        & ! contains information about Q matrix
	 	    work, lwork, info)

if     (info /=0 ) then
	write(*,20) '*****************************************************'
	write(*,20) '* ERROR: Constructing Q matrix failed!Something must*'
	write(*,20) '*        be very wrong... program stops!            *'
	write(*,20) '*****************************************************'
	stop
end if


! The array A_QR_decomposition now contains the first 
! special_functions_basis_size columns of the Q matrix.
! These columns correspond to the coefficients of the optimal
! basis in the C basis

allocate(b_basis_C_coefficients(ortho_basis_size,special_functions_basis_size))

b_basis_C_coefficients(:,:) = A_QR_decomposition(:,1:special_functions_basis_size)

! clean out the trash
deallocate(work)
deallocate(A_QR_decomposition)
deallocate(tau_QR_decomposition)


! make sure the b_basis is orthonormal

allocate(overlaps_b(special_functions_basis_size,special_functions_basis_size))

overlaps_b(:,:) = matmul(transpose(conjg(b_basis_C_coefficients)),b_basis_C_coefficients)

do I_index = 1, special_functions_basis_size
	overlaps_b(I_index,I_index) = overlaps_b(I_index,I_index)-cmplx_1 
end do

real_norm = sqrt(real(sum(conjg(overlaps_b(:,:))*overlaps_b(:,:))))

write(*,31) '|       | <b_i| b_j> - delta_ij | = ',real_norm,' |'

call get_timing(time2)
write(*,30) '|                  time    :',time2-time1,' seconds      |'
write(*,20) '|           ---------------------------------       |'


!================================================================================
!	Compute the kernel in the b basis
!================================================================================

allocate(C_basis_kernel(ortho_basis_size,ortho_basis_size))
allocate(b_basis_kernel(special_functions_basis_size,special_functions_basis_size))


! build the kernel in the | b > sub-basis 
!
! Note that  | C_I > = sum_{L} | B_L > Pc_{LI}
! Thus,  KC_{IJ} = < C_I | fxc | C_J > = sum_{L1L2} Pc^dagger_{IL1} KB_{L1L2} Pc_{L2J}, 
!
!
! Also,      | b_i > = sum_{I} | C_I > b_{Ii} 
!

! Thus,     < b_i | fxc | b_j > = sum_{I1I2} b^dagger_{iI1} KC_{I1I2} b_{I2j}, 

C_basis_kernel(:,:) = matmul(transpose(conjg(Pc)),   &
				matmul(bare_basis_kernel,Pc))

b_basis_kernel(:,:) = matmul(transpose(conjg(b_basis_C_coefficients)),   &
				matmul(C_basis_kernel,b_basis_C_coefficients))

! output the b_basis_kernel, for testing and checking purposes.
call output_full_projected_kernel_netcdf(qvec,special_functions_basis_size,b_basis_kernel)

!================================================================================
!	Create the functional projectors onto the sub-basis containing
!	optimal_basis_size functions.
!================================================================================

! The functional projectors are given by
!		T_{iL} = < b_i | B_L > 
!		       = sum_{I} b^dagger_{iI} < C_I | B_L > 
!		       = sum_{I} b^dagger_{iI} sum_{L1} Pc^dagger_{IL1} < B_L1 | B_L > 
! thus
!		T_{iL} = sum_{I,L1} b^dagger_{iI} Pc^dagger_{IL1} O_{L1L}


allocate(PT(bare_basis_size,optimal_basis_size))

! truncate the t basis to the specified size
PT(:,:) = matmul(Pc,b_basis_C_coefficients)

!        Finally, compute the functional overlaps, which will be used
!	to compute chi_KS
allocate(T_functional_projectors(optimal_basis_size,bare_basis_size))

T_functional_projectors(:,:) = matmul(transpose(conjg(PT)),overlaps)


!================================================================================
!        Output the projectors and other necessary info to netcdf file.
!================================================================================

! create the arrays which define the Wanner band indices
L_index   = 0

allocate(list_m1(bare_basis_size))
allocate(list_m2(bare_basis_size))

do nb_W1 = 1,num_wann
   do nb_W2 = 1,num_wann
	L_index  = L_index + 1

	list_m1(L_index) = nb_W1
	list_m2(L_index) = nb_W2

   end do
end do


! project the special functions onto the | b > basis
! The functions are ordered as { 1, mz, Delta, [plane waves]}. 
!
!	We have | g_j > = sum_i | b_i> g_{ij}, such that
!	         g_{ij} =  < b_i | g_j >. 
!
!	But this is just the R matrix of the QR decomposition.
!	Let's compute it explicitely, as a check.
!

! compute the coefficients in the b basis!

allocate(special_functions_b_coefficients(optimal_basis_size,special_functions_basis_size))

special_functions_b_coefficients = matmul(transpose(conjg(b_basis_C_coefficients)), &
						special_functions_C_coefficients)

!================================================================================
!        Output functional projectors to netcdf file
!================================================================================


write(*,20) '|       - Outputing the functional projectors to    |'
write(*,20) '|         netcdf file...                            |'

call output_functional_projectors_netcdf(                           &
		qvec,bare_basis_size,special_functions_basis_size,  &
		optimal_basis_size,T_functional_projectors,         &
		list_m1,list_m2,special_functions_b_coefficients,   &
		PT,b_basis_kernel)


!================================================================================
!        Output some useful information to an ascii file
!================================================================================

io_unit = find_free_unit()
open(io_unit,file='optimal_basis.dat')

write(io_unit,20) '#======================================================='
write(io_unit,20) '# This file contains useful information related to the  ' 
write(io_unit,20) '# Computation of the optimal basis.                     ' 
write(io_unit,20) '#======================================================='
write(io_unit,20) ''
write(io_unit,22) ' The size of the optimal basis is    ',optimal_basis_size   ,' '
write(io_unit,20) ''
write(io_unit,22) ' The size of the orthogonal basis is ',ortho_basis_size   ,' '
write(io_unit,20) ''
write(io_unit,20) ''
write(io_unit,20) ' The special functions are represented as '
write(io_unit,20) '     | g_j > = sum_{i}  |b_i> g_{ij}.     '
write(io_unit,20) ''
write(io_unit,20) ' The coefficients are given by:           '
write(io_unit,20) ''
write(io_unit,20) '# orthonormal basis  |                    special functions                   '

header_string = '#    | b_i >         |        one(r)        |         mz(r)        |       Delta(r)        |' 
line_string   = '#==========================================================================================='


if (number_G_shells >1) then

   if (trim(plane_wave_combination) == 'symmetric'  .or.    &
       trim(plane_wave_combination) == 'mz_times_symmetric') then


     do iG = 1, pw_dim

	write(G_string,22) ' shell #',iG+1,'          |'

	header_string = trim(header_string)//trim(G_string)
	line_string   = trim(line_string)//trim('=======================')
     end do ! iG 

   else if (trim(plane_wave_combination) == 'all') then
     do iG = 2, nG_shell_max

	write(G_string,34) ' G = (',gvec(1,iG),',',     &
                                    gvec(2,iG),',',     &
                                    gvec(3,iG),')    |'

	header_string = trim(header_string)//trim(G_string)
	line_string   = trim(line_string)//trim('=======================')
     end do ! iG 

   end if
end if

write(io_unit,20) trim(header_string)
write(io_unit,20) trim(line_string)


do I_index = 1,optimal_basis_size
	write(io_unit,32) I_index,special_functions_b_coefficients(I_index,:)
end do


! build the g(r) functions explicitely and compare with the G(r) functions 


! first, extract their coefficients in the | B_L > basis
allocate(g_functions_B_coefficients(bare_basis_size,special_functions_basis_size))

g_functions_B_coefficients(:,:) = matmul(PT, special_functions_b_coefficients)


! Next, compute the overlaps

allocate(eta_array(special_functions_basis_size))




eta_normalization = sum(conjg(one_r(:))*one_r(:))/dble(nr)
eta_array(1)      = sqrt(abs(                                       &
			sum(conjg(g_functions_B_coefficients(:,1))* &
		                special_functions_B_overlaps(:,1))  &
				                      /eta_normalization ))

eta_normalization = sum(conjg(mz_r(:))*mz_r(:))/dble(nr)
eta_array(2)      = sqrt(abs(                                       &
			sum(conjg(g_functions_B_coefficients(:,2))* &
		                special_functions_B_overlaps(:,2))  &
				                      /eta_normalization ))

eta_normalization = sum(conjg(Delta_r(:))*Delta_r(:))/dble(nr)
eta_array(3)      = sqrt(abs(                                       &
			sum(conjg(g_functions_B_coefficients(:,3))* &
		                special_functions_B_overlaps(:,3))  &
				                      /eta_normalization ))

if (number_G_shells >1) then
   do iG = 1, pw_dim
	eta_normalization = ONE
	eta_array(3+iG) = sqrt(abs(                                      &
			sum(conjg(g_functions_B_coefficients(:,3+iG))*   &
		                special_functions_B_overlaps(:,3+iG))    &
				                      /eta_normalization ))

   end do ! iG 
end if

write(io_unit,20) ''
write(io_unit,20) ' The quality of the basis can be evaluated using'
write(io_unit,20) '        eta_i = sqrt{| < G_i | g_i > / <G_i | G_i > | } '
write(io_unit,20) '  where G_i is the original function and g_i its projection. '
write(io_unit,20) ''
write(io_unit,20) '       function                  eta_i '
write(io_unit,20) '#============================================================================='

write(io_unit,33) '          1          ',eta_array(1)
write(io_unit,33) '          mz         ',eta_array(2)
write(io_unit,33) '         Delta       ',eta_array(3)

if (number_G_shells >1) then


   if (trim(plane_wave_combination) == 'symmetric'  .or.    &
       trim(plane_wave_combination) == 'mz_times_symmetric') then

     do iG = 1, pw_dim
	write(io_unit,35) '         shell',iG+1,eta_array(3+iG)
     end do ! iG

   else if (trim(plane_wave_combination) == 'all') then
     do iG = 2, nG_shell_max
	write(io_unit,34) ' G = (',gvec(1,iG),',',     &
                                   gvec(2,iG),',',     &
                                   gvec(3,iG),')',     &
				   eta_array(3+iG-1)
     end do ! iG 
   end if

end if


close(io_unit)

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
16 format(A,I6,A,I6,A)
17 format(F12.6)
18 format(A,ES8.2,A)
19 format(ES18.6)
20 format(A)
22 format(A,I4,A)
30 format(A,F10.3,A)
31 format(A,ES15.3,A)
32 format(6X,I6,10X,1000(F8.4,2X,F8.4,5x))
33 format(A,5x,F12.8)
34 format(A,I3,A,I3,A,I3,A,8X,F12.8)
35 format(A,I4,8X,F12.8)

end program generate_optimal_basis
