!       program generate_optimal_basis_plasmons
!       -------------------------------------------
!
!       This is a "sub program" which is part of the intw project.
!
!       The purpose of this program is to generate the "optimal basis", namely
!	the combination of bare product basis functions {|B>} which generate
!	an orthonormal basis {|b>} which :
!		- spans at least the functional space {1},
!		- diagonalizes the Hartree-exchange-correlation kernel as much as possible
!		- is ordered to make the Hartree-exchange-correlation kernel as block
!		  diagonal as possible.
!
!	This program is mainly useful for plasmons computations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program generate_optimal_basis_plasmons
  
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


real(dp) :: qvec(3)
real(dp) :: q_plus_G(3), q_plus_G_cart(3)

logical  :: use_IBZ
logical  :: mz_file_exists, Delta_file_exists

integer  :: ir, nr
integer  :: ir1, ir2, ir3
integer  :: iG
integer  :: switch_singlet_to_triplet
real(dp) :: rvec(3)


integer  :: number_of_symmetry_G_shells

integer  :: ishell
integer,allocatable  :: shells_table(:,:)

integer  :: n
integer  :: bare_basis_size, ortho_basis_size, special_functions_basis_size
integer  :: specific_basis_size,basis_size_m3
integer  :: F_basis_size
integer  :: L_index, I_index, L2_index

integer  :: I1_index, I2_index

complex(dp)             :: graham_schmidt_overlap

complex(dp),allocatable :: overlaps(:,:)


complex(dp),allocatable :: overlaps_d(:,:)


complex(dp),allocatable :: Vc(:,:)
real(dp),allocatable    :: Sc(:)
complex(dp),allocatable :: Pc(:,:)

complex(dp),allocatable :: VF(:,:)
real(dp),allocatable    :: SF(:)
complex(dp),allocatable :: PF(:,:)

complex(dp),allocatable :: Pd(:,:), Pd_inverse(:,:)
complex(dp),allocatable :: Pb(:,:)
complex(dp),allocatable :: Pt(:,:)

complex(dp),allocatable :: norm_checking_array(:,:)
complex(dp),allocatable :: A(:,:), pseudo_inverse_A(:,:)
real(dp),allocatable    :: Sigma_A(:)
complex(dp),allocatable :: matrix_Sigma_A(:,:)
complex(dp),allocatable :: inverse_Sigma_A(:,:)
complex(dp),allocatable :: work(:)
complex(dp),allocatable :: U_A(:,:), V_dagger_A(:,:)
real(dp),allocatable    :: rwork(:)

integer                 :: info, lwork,lrwork


complex(dp),allocatable :: P_tmp(:,:)

complex(dp),allocatable :: T_functional_projectors(:,:)

complex(dp),allocatable :: field_r(:)

! function in r, integrals with bare basis, coefficient in optimal basis

complex(dp),allocatable :: bare_basis_kernel(:,:)
complex(dp),allocatable :: bare_basis_kernel_xc(:,:)
complex(dp),allocatable :: bare_basis_kernel_Ha(:,:)

complex(dp),allocatable :: C_basis_kernel(:,:)
complex(dp),allocatable :: C_basis_kernel_xc(:,:)
complex(dp),allocatable :: C_basis_kernel_Ha(:,:)



complex(dp),allocatable :: d_basis_kernel(:,:)
complex(dp),allocatable :: F_basis_kernel(:,:)
complex(dp),allocatable :: h_basis_kernel(:,:)

complex(dp),allocatable :: t_basis_kernel(:,:)
complex(dp),allocatable :: t_basis_kernel_xc(:,:)
complex(dp),allocatable :: t_basis_kernel_Ha(:,:)

complex(dp),allocatable :: optimal_basis_kernel(:,:)
complex(dp),allocatable :: optimal_basis_kernel_xc(:,:)
complex(dp),allocatable :: optimal_basis_kernel_Ha(:,:)

real(dp),   allocatable :: kernel_to_be_sorted(:)
real(dp),   allocatable :: wing_column(:)
integer ,   allocatable :: P_sort(:)

complex(dp),allocatable :: permutation(:,:)


complex(dp),allocatable :: one_r(:)
complex(dp),allocatable :: mz_r(:)
complex(dp),allocatable :: Delta_r(:)

integer                 :: pw_dim
integer                 :: dim1, dim2
complex(dp),allocatable :: plane_waves(:,:)

complex(dp),allocatable :: special_functions_B_overlaps(:,:)

complex(dp),allocatable :: special_functions_C_coefficients(:,:)

complex(dp),allocatable :: special_functions_t_coefficients(:,:)

complex(dp),allocatable :: d_basis_C_coefficients(:,:)
complex(dp),allocatable :: F_basis_C_coefficients(:,:)
complex(dp),allocatable :: h_basis_C_coefficients(:,:)
complex(dp),allocatable :: b_basis_C_coefficients(:,:)
complex(dp),allocatable :: t_basis_C_coefficients(:,:)


complex(dp),allocatable :: overlaps_g(:,:)
complex(dp),allocatable :: Vg(:,:)
real(dp),allocatable    :: Sg(:)

complex(dp),allocatable :: Vk(:,:)
real(dp),allocatable    :: Sk(:)

complex(dp)             :: eta_normalization
complex(dp)             :: norm
real(dp)                :: real_norm

real(dp)                :: q_direction(3)
integer  :: B_nc_id

complex(dp),allocatable :: B_function(:)
complex(dp),allocatable :: d_basis_functions(:,:)
complex(dp),allocatable :: b_basis_functions(:,:)

complex(dp),allocatable :: g_functions_B_coefficients(:,:)
real(dp),allocatable    :: eta_array(:)


integer              :: nb_W1, nb_W2
integer,allocatable  :: list_m1(:), list_m2(:)

integer              :: io_unit

character(10000)     :: header_string, line_string,G_string
!--------------------------------------------------------------------------------


! start global clock
call get_timing(total_time1)

!================================================================================
!        Talk to user
!================================================================================
write(*,20) '====================================================='
write(*,20) '|   program generate_optimal_basis_plasmons         |'
write(*,20) '|   --------------------------------------          |'
write(*,20) '====================================================='
write(*,20) '|                                                   |'
write(*,20) '|        This program computes the optimal basis    |'
write(*,20) '|        for plasmons calculations.                 |'
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
!        Tell the user if the system has spin; if so, stop computation!
!================================================================================

if (nspin==2 .or. magnon) then
   write(*,20) '*****************************************************'
   write(*,20) '* ERROR: This program is for plasmons only!         *'
   write(*,20) '* nspin must 1, and magnon = .false.                *'
   write(*,20) '*  program stops.                                   *'
   write(*,20) '*****************************************************'
   stop
endif

write(*,20) '====================================================='
write(*,20) '|   The calculation is paramagnetic nspin=1         |'
write(*,20) '|   Computing optimal basis for plasmon calculation |'
write(*,20) '|                                                   |'
write(*,20) '|                 ----------------                  |'
write(*,20) '|   This program is still in development. It assumes|'
write(*,20) '|   the basis is entirely made of plane waves. In   |'
write(*,20) '|   order to use it, the parameters must be         |'
write(*,20) '|                                                   |'
write(*,20) '|    optimal_basis_size     = number of plane waves |'
write(*,20) '|    plane_wave_combination = "all"                 |'
write(*,20) '|                                                   |'
write(*,20) '|     Also, the vector q cannot be Gamma.           |'
write(*,20) '====================================================='

if (trim(plane_wave_combination) /= 'all' ) then 
   write(*,20) '*****************************************************'
   write(*,20) '* ERROR: Only plane_wave_combination = all          *'
   write(*,20) '*        is implemented at this time!               *'
   write(*,20) '*        Program stops.                             *'
   write(*,20) '*****************************************************'
   stop
end if

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
!				O = Vc * Sc * Vc^dagger,   Vc unitary
!
!				Pc_{LI} = Vc_{LI}/sqrt{Sc_{II}}, Sc_{II} > tolerance
!
!		    | C_I >  : temporary orthogonal basis functions, built 
!				by orthogonalizing the overlaps O and throwing
!				out combinations with eigenvalues smaller than
!				a given tolerance.
!
! 				| C_I > = sum_L | B_L > Pc_{LI} 
!
!		     | G_i > : special functions. They are different 
!			       depending if its a magnon or plasmon
!			       calculation.
!
!		        --->   Magnon calculation:
!                              The | G_i > are the special functions 
!			       necessary to impose the acoustic
!			       condition:
!
!				i       function
!                               -------------------
!				1        one
!				2        mz(r)
!				3        Delta(r)
!
!		         --->   Plasmon calculation:
!                               The | G_i > are chose to be:
!				i       function
!                               -------------------
!				i        exp(iGr)
!				
!				The following applies for both 
!				plasmon and magnon cases:
!
!				| g_i > = sum_I | C_I > g_{Ii}
!
!				g_{Ii}  = <C_I | G_i > 
!				        = sum_{L} Pc^dagger_{IL} <B_L | G_i > 
!
!				The difference between | G_i > and |g_i > 
!				indicates the shortcomings of the basis in describing
!				the special functions.
!
!				
!
!  Og_{ij} = < g_i | g_j >   : overlap of special functions projections onto
!				the orthonormal basis.
!
!				Og = Vg * Sg * Vg^dagger
!				Pd_{ij} = Vg_{ij}/sqrt{Sg_{jj}}
!
!		     | d_i > : orthonormal basis functions, obtained by diagonalizing
!				Og;
!				| d_i > = sum_{j} |g_j > Pd_{ji}
!
!				This basis is orthogonal and spans a subspace
!				containing the special functions as much as possible.
!				However, there remains some freedom. Let's use this
!				freedom to diagonalize the kernel!
!
!  K_{ij} = < d_i | fxc | d_j > : exchange-correlation kernel.
!				 K = Vk * Sk * Vk^dagger
!
!		     | b_i > : orthonormal basis functions which diagonalize the kernel,
!				the optimal basis!
!				| b_i > = sum_{j} |d_j > Vk_{ji}
!
!	Thus, we have
!		| b_i > = sum_L |B_L > Pb_{Li}
!
!		Pb_{Li} = sum_{I,l,j} Pc_{LI} g_{Il} Pd_{lj} Vk_{ji}
!
!
!	In order to go beyond a three function basis, let's build orthogonal
!	functions in the remaining subspace.
!
!		    | F_I >  : temporary orthogonal basis functions orthogonal
!			       to the {|b_i > } set, and built from the | C_I >
!			       basis functions.
!
!	Define the projector:
!
!		P_basis = sum_I | C_I >< C_I |  - \sum_{i=1}^{3} | b_i > < b_i | 
!
!		Clearly, 
!				P_basis | b_i > = 0
!				P_basis | F_I > = | F_I > 
!
!		The projector P_basis should have three eigenvectors with
!		eigenvalues zero (the | b_i > !), and ortho_basis_size-3
!		basis functions with eigenvalue 1 ( the | F_I > ).
!
!		Once again, there is a lot of arbitariness in this basis of |F_I>;
!		Diagonalize the kernel in this subspace!
!
!  KK_{IJ} = < F_I | fxc | F_J > : exchange-correlation kernel.
!				 KK = Vkk * Skk * Vkk^dagger
!
!		| f_I > : orthonormal basis functions which diagonalize the kernel,
!			  in this subspace
!			| f_I > = sum_{J} |F_J > Vkk_{JI}
!
!
!  We thus obtain a larger basis, { | b_i > , | f_I > }, where the { | b_i > }
!  span the sub-space of the special functions and { | f_I > } its orthogonal
!  subspace; the kernel is diagonalized in both subspaces.
!
!  This basis set is very large, however (of dimension ortho_basis_size!), and 
!  the user will certainly want to truncate it. It would be nice to find a really
!  clever, physically motivated way of cherry picking basis functions, but
!  until such a criterion is obtained (if ever), we'll have to content ourselves
!  with sorting the basis functions and hoping for the best.
!
!  Define the temporary basis { | h_I > } = { | b_i > , | f_I > };
!
!  the kernel in this basis is given by
!
!	Kh_{IJ} = < h_I | fxc | h_J > 
!
!	Kh      = [       |                ]
!	          [ diag  |     Wing       ]
!	          [------------------------]
!	          [       |                ]
!	          [       |                ]
!	          [ Wing  |   Diagonal     ]
!	          [       |                ]
!	          [       |                ]
!
!	Finally define the basis as 
!
!		| t_I > = sum_{J} | h_J > Permutation_{JI}
!
!	where the permutation insures that the Wing has is weight as close to the
!	diagonal as possible.
!	The { | t_I > } basis should be orthogonal, should have { | b_i > } has
!	its first three basis functions and the rest of the functions should
!	be somewhat optimized. Again, a clever idea on how to select basis 
!	functions is still lacking.
!-------------------------------------------------------------------------------


! test that the chosen optimal basis size isn't larger than the bare basis size!

! We are neglecting the dependence in R of the basis function;
! that could be improved in the future.
bare_basis_size = num_wann**2


!================================================================================
!        Read the bare basis XC kernel and the overlap matrix to get Pc.
!================================================================================
allocate(overlaps(bare_basis_size,bare_basis_size ))

allocate(bare_basis_kernel(bare_basis_size,bare_basis_size ))
allocate(bare_basis_kernel_xc(bare_basis_size,bare_basis_size ))
allocate(bare_basis_kernel_Ha(bare_basis_size,bare_basis_size ))

allocate(Vc(bare_basis_size,bare_basis_size))
allocate(Sc(bare_basis_size))


call read_overlaps_plasmons_netcdf(bare_basis_size,       &
				   overlaps,              &
		                   bare_basis_kernel,     &
				   bare_basis_kernel_xc,  &
				   bare_basis_kernel_Ha,qvec)


call diagonalize_matrix(bare_basis_size,overlaps,Vc,Sc)


ortho_basis_size = 0

do L_index = 1, bare_basis_size
	if ( Sc(L_index)/maxval(Sc(:))  >= tolerance ) then
		ortho_basis_size = ortho_basis_size + 1
	end if
end do

allocate(Pc(bare_basis_size,ortho_basis_size))

I_index = 0
do L_index = 1,bare_basis_size
	if ( Sc(L_index)/maxval(Sc(:))  >= tolerance ) then
		I_index = I_index + 1
		Pc(:,I_index) = &
			Vc(:,L_index)/sqrt(Sc(L_index))
	end if
end do

deallocate(Vc)
deallocate(Sc)


write(*,20) '|           For the given tolerance, the            |'
write(*,22) '|           orthonormal basis has dimension',ortho_basis_size ,'     |'
write(*,20) '====================================================='

if (optimal_basis_size > ortho_basis_size ) then
  write(*,20) '*******************************************************'
  write(*,20) '*   ERROR: The optimal_basis_size specified on input  *'
  write(*,20) '*          is larger than the ortho_basis_size; this  *'
  write(*,20) '*          is impossible to accomodate. Pick a smaller*'
  write(*,20) '*          optimal_basis_size!                        *'
  write(*,20) '*               program stops.                        *'
  write(*,20) '*******************************************************'
  stop
end if


nr = nr1*nr2*nr3

! In density calculation, the optimal
! basis includes the function 1 and
! a finite number of plane waves.

allocate(one_r(nr))
one_r(:) = cmplx_1


!================================================================================
!	Get the integrals of the relevant functions with the 
!	bare B basis functions.
!================================================================================

if (number_G_shells == 1) then
	pw_dim = 0

else if (number_G_shells >1) then
   call get_G_shells(number_G_shells,nG_shell_max)

   if (optimal_basis_size /= nG_shell_max) then
      write(*,20) '*****************************************************'
      write(*,20) '* ERROR: the optimal basis size must be set exactly *'
      write(*,20) '*        to the number of plane waves!              *'
      write(*,20) '*        Input value for optimal_basis_size:        *'
      write(*,22) '*            ',optimal_basis_size,'                                   *'
      write(*,20) '*        Number of plane waves:                     *'
      write(*,22) '*            ',nG_shell_max,'                                   *'
      write(*,20) '*        Program stops.                             *'
      write(*,20) '*****************************************************'
      stop
   end if


   pw_dim = nG_shell_max-1


   write(*,20) '|       - The input file specifies that the basis   |'
   write(*,20) '|         should contain plane waves. For the       |'
   write(*,20) '|         specified parameter number_G_shells,      |'
   write(*,20) '|         the number of plane waves not equal to 1  |'
   write(*,22) '|         is ',pw_dim,'.                                  |'
   write(*,20) '|           ---------------------------------       |'
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

write(*,20) '|       - computing overlaps with B functions...    |'
call get_timing(time1)


! order the special functions as {1, [plane waves]}.
special_functions_basis_size = 1+pw_dim 

! Prepare for computing overlaps between the special 
! functions and the bare basis functions.
allocate(special_functions_B_overlaps(bare_basis_size,special_functions_basis_size ))

allocate(B_function(nr))

! open to read (but not to write!)
call open_B_functions_netcdf(B_nc_id)


do L_index = 1, bare_basis_size

   call read_B_function_netcdf(B_nc_id,L_index,B_function,nr)


   special_functions_B_overlaps(L_index,1) =    &
				sum(conjg(B_function(:))*one_r(:))/dble(nr)

   if (nG_shell_max > 1) then
      do iG = 1,pw_dim
         special_functions_B_overlaps(L_index,1+iG) =    &
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
!        Create the "d" basis, which spans the space of the special functions,
!	expressed in the orthogonal |C_I> basis.
!================================================================================
! This variable is defined only to make it very explicit 
! the size of the special functions WITHOUT plane waves

! it represents the size of the basis {1}
specific_basis_size = 1
write(*,20) '|       - Generating basis functions spanning the   |'
write(*,20) '|         subspace {1} ...                          |'



! Make sure the first basis function span {1}
allocate(overlaps_g(specific_basis_size ,specific_basis_size ))
allocate(Vg(specific_basis_size,specific_basis_size))
allocate(Sg(specific_basis_size))
allocate(Pd(specific_basis_size,specific_basis_size))

allocate(d_basis_C_coefficients(ortho_basis_size,special_functions_basis_size))

call get_timing(time1)

overlaps_g(:,:) = matmul(                                                &
		transpose(conjg(                                         &
		special_functions_C_coefficients(:,1:specific_basis_size))),    &
		special_functions_C_coefficients(:,1:specific_basis_size))

! Diagonalize the overlap matrix using LAPACK within MKL

call diagonalize_matrix(specific_basis_size,overlaps_g,Vg,Sg)


! create the projectors to this basis

do I_index = 1, specific_basis_size
	Pd(:,I_index) = Vg(:,I_index)/sqrt(Sg(I_index))
end do


! create the first basis function, which span {1}
d_basis_C_coefficients(:,1:specific_basis_size) =   &
	 matmul(special_functions_C_coefficients(:,1:specific_basis_size),Pd)

deallocate(overlaps_g)
deallocate(Vg)
deallocate(Pd)

call get_timing(time2)
write(*,30) '|                  time    :',time2-time1,' seconds      |'
write(*,20) '|           ---------------------------------       |'

! define the other basis functions, if need be, using Graham-Schmidt 
! orthogonalization. 

if (number_G_shells >1) then
  call get_timing(time1)

   write(*,20) '|       - Generating basis functions spanning the   |'
   write(*,20) '|         subspace of the plane waves ...           |'


   do iG = 1,pw_dim

      I1_index = specific_basis_size+iG


      ! initialize 
      d_basis_C_coefficients(:,I1_index) =                  &
		special_functions_C_coefficients(:,I1_index)

      ! subtract the components from basis functions which are already in 
      ! the basis
      do I2_index = 1,I1_index-1

	graham_schmidt_overlap =                               &
	        sum(conjg(d_basis_C_coefficients(:,I2_index))* &
		special_functions_C_coefficients(:,I1_index))


	d_basis_C_coefficients(:,I1_index) =                 &
		d_basis_C_coefficients(:,I1_index)           &
		-graham_schmidt_overlap*d_basis_C_coefficients(:,I2_index)

      end do ! I2_index

      ! normalize
      norm = sqrt(sum(                                       &
		conjg(d_basis_C_coefficients(:,I1_index))*   &
	              d_basis_C_coefficients(:,I1_index)))

      d_basis_C_coefficients(:,I1_index) = d_basis_C_coefficients(:,I1_index)/norm

   end do ! iG

   ! make sure the d_basis is orthonormal

   allocate(overlaps_d(special_functions_basis_size,special_functions_basis_size))

   overlaps_d(:,:) = matmul(transpose(conjg(d_basis_C_coefficients)),d_basis_C_coefficients)

   do I_index = 1, special_functions_basis_size
	overlaps_d(I_index,I_index) = overlaps_d(I_index,I_index)-cmplx_1 
   end do

   real_norm = sqrt(real(sum(conjg(overlaps_d(:,:))*overlaps_d(:,:))))

   write(*,31) '|       | <d_i| d_j> - delta_ij | = ',real_norm,' |'

   call get_timing(time2)
   write(*,30) '|                  time    :',time2-time1,' seconds      |'
   write(*,20) '|           ---------------------------------       |'

end if

!================================================================================
!	Compute the kernel in the d basis, and extract the
!	optimal basis.
!================================================================================

allocate(C_basis_kernel(ortho_basis_size,ortho_basis_size))
allocate(C_basis_kernel_xc(ortho_basis_size,ortho_basis_size))
allocate(C_basis_kernel_Ha(ortho_basis_size,ortho_basis_size))



allocate(d_basis_kernel(special_functions_basis_size,special_functions_basis_size))

allocate(b_basis_C_coefficients(ortho_basis_size,special_functions_basis_size))
   
   
! build the kernel in the | d > sub-basis 
!
! Note that  | C_I > = sum_{L} | B_L > Pc_{LI}
! Thus,  KC_{IJ} = < C_I | vc+fxc | C_J > = sum_{L1L2} Pc^dagger_{IL1} KB_{L1L2} Pc_{L2J}, 
!
!
! Also,      | d_i > = sum_{I} | C_I > d_{Ii} 
!
! Thus,     < d_i | vc+fxc | d_j > = sum_{I1I2} d^dagger_{iI1} KC_{I1I2} d_{I2j}, 
   
C_basis_kernel(:,:)    = matmul(transpose(conjg(Pc)),   &
   				matmul(bare_basis_kernel,Pc))

C_basis_kernel_xc(:,:) = matmul(transpose(conjg(Pc)),   &
   				matmul(bare_basis_kernel_xc,Pc))

C_basis_kernel_Ha(:,:) = matmul(transpose(conjg(Pc)),   &
   				matmul(bare_basis_kernel_Ha,Pc))


   
d_basis_kernel(:,:) = matmul(transpose(conjg(d_basis_C_coefficients)),   &
   				matmul(C_basis_kernel,d_basis_C_coefficients))
   
   
!        Diagonalize the kernel in the two subspaces of {1} and
!        the plane waves
   
allocate(Vk(specific_basis_size,specific_basis_size))
allocate(Sk(specific_basis_size))
   
call diagonalize_matrix(specific_basis_size,                           &
       d_basis_kernel(1:specific_basis_size,1:specific_basis_size),Vk,Sk)
   
! compute the coefficients of the optimal basis in the | C_I > basis
   
b_basis_C_coefficients(:,1:specific_basis_size) = &
   		matmul(d_basis_C_coefficients(:,1:specific_basis_size),Vk)
   

deallocate(Vk)
deallocate(Sk)
   
if (number_G_shells >1) then
   
   allocate(Vk(pw_dim,pw_dim))
   allocate(Sk(pw_dim))
   
   dim1 = specific_basis_size+1
   dim2 = specific_basis_size+pw_dim
   call diagonalize_matrix(pw_dim, d_basis_kernel(dim1:dim2,dim1:dim2), Vk,Sk)
   
   b_basis_C_coefficients(:,dim1:dim2) = matmul(d_basis_C_coefficients(:,dim1:dim2),Vk)
   
   deallocate(Vk)
   deallocate(Sk)
end if
   
   
! make sure the b basis is still orthonormal
allocate(norm_checking_array(special_functions_basis_size,special_functions_basis_size))
   
norm_checking_array(:,:) = matmul(transpose(conjg(b_basis_C_coefficients)),b_basis_C_coefficients)
   
   
do I_index = 1, special_functions_basis_size
   norm_checking_array(I_index,I_index) = norm_checking_array(I_index,I_index)-cmplx_1 
end do
   
real_norm = sqrt(real(sum(conjg(norm_checking_array(:,:))*norm_checking_array(:,:))))

write(*,31) '|       | <b_i| b_j> - delta_ij | = ',real_norm,' |'
write(*,20) '|           ---------------------------------       |'


deallocate(norm_checking_array)

deallocate(d_basis_kernel)
deallocate(d_basis_C_coefficients)

!================================================================================
!	Build the sub-space orthogonal to the optimal basis.
!================================================================================

! Build the projector operator, which is orthogonal to { | b_i > }
allocate(P_tmp(ortho_basis_size,ortho_basis_size))
    
allocate(F_basis_C_coefficients(ortho_basis_size,       &
			ortho_basis_size-special_functions_basis_size))

allocate(VF(ortho_basis_size,ortho_basis_size))
allocate(SF(ortho_basis_size))

P_tmp(:,:) = cmplx_0
do I_index = 1, ortho_basis_size
   P_tmp(I_index,I_index) =  cmplx_1
end do

do I1_index = 1, ortho_basis_size
   do I2_index = 1, ortho_basis_size
	P_tmp(I1_index,I2_index) =  P_tmp(I1_index,I2_index)    & 
				    -sum(                       &
                       b_basis_C_coefficients(I1_index,:)       &
                *conjg(b_basis_C_coefficients(I2_index,:)))
   end do
end do

! extract the | F_I > basis

call diagonalize_matrix(ortho_basis_size,P_tmp,VF,SF)


F_basis_size = 0
do I_index = 1, ortho_basis_size
   ! if the eigenvalue is small, skip the rest!
   if ( SF(I_index) > 0.1) then
   	F_basis_size = F_basis_size +1
   end if
end do


if (F_basis_size /=  ortho_basis_size-special_functions_basis_size) then
   write(*,*) '*************************************************'
   write(*,*) '* ERROR: the size of the | F > basis is not     *'
   write(*,*) '*        what it is supposed to be...           *'
   write(*,*) '*        This can happen if the bare basis      *'
   write(*,*) '*        nc file is corrupt...                  *'
   write(*,*) '*************************************************'
   write(*,*) 'F_basis_size = ' ,F_basis_size 
   stop
end if


I1_index = 0
do I_index = 1, ortho_basis_size
   ! if the eigenvalue is small, skip the rest!
   if ( SF(I_index) > 0.1) then
      I1_index = I1_index + 1
      F_basis_C_coefficients(:,I1_index) = VF(:,I_index)
   end if
end do
    
    
    
deallocate(P_tmp)
deallocate(VF)
deallocate(SF)

!================================================================================
!	Diagonalize kernel in | F > basis, and directly create the 
!       h_basis
!================================================================================
allocate(F_basis_kernel(F_basis_size ,F_basis_size ))
    
allocate(h_basis_C_coefficients(ortho_basis_size, ortho_basis_size))
    
allocate(VF(F_basis_size,F_basis_size))
allocate(SF(F_basis_size))
    
F_basis_kernel(:,:) = matmul(transpose(conjg(F_basis_C_coefficients)),   &
			matmul(C_basis_kernel,F_basis_C_coefficients))
    

h_basis_C_coefficients(:, 1:special_functions_basis_size) = b_basis_C_coefficients(:,:)
    
call diagonalize_matrix(F_basis_size,F_basis_kernel,VF,SF)
    
    
h_basis_C_coefficients(:, special_functions_basis_size+1:ortho_basis_size) =    &
				matmul(F_basis_C_coefficients,VF)
    
deallocate(F_basis_kernel)
deallocate(b_basis_C_coefficients)
deallocate(VF)
deallocate(SF)

!================================================================================
!	Compute the kernel in the | h_I > basis, and sort! 
!================================================================================
allocate(h_basis_kernel(ortho_basis_size,ortho_basis_size))
    
h_basis_kernel(:,:) = matmul(transpose(conjg(h_basis_C_coefficients)),   &
    				matmul(C_basis_kernel,h_basis_C_coefficients))
    
    
! Build a permutation matrix which will order the basis appropriately
allocate(permutation(ortho_basis_size,ortho_basis_size))
permutation = cmplx_0
    
! maintain the order of the first special functions
do I_index = 1, special_functions_basis_size
   permutation(I_index,I_index) = cmplx_1
end do
    
    
allocate(kernel_to_be_sorted(F_basis_size))
    
! Decide what should be sorted; the wings, or the diagonal?
if (trim(basis_order) == 'wings') then
    
    	allocate(wing_column(special_functions_basis_size))
    	allocate(P_sort(special_functions_basis_size))
    
    	! first, find the largest contribution in each column of the wing
    	do I_index = 1, F_basis_size
    		wing_column(:) = abs(h_basis_kernel(1:special_functions_basis_size,  &
    						special_functions_basis_size+I_index))
    
    		call HPSORT_real(special_functions_basis_size,wing_column,P_sort)
    
    		! the largest value is the last one
    		kernel_to_be_sorted(I_index) = wing_column(special_functions_basis_size)
    
    
    	end do 
    
    	deallocate(wing_column)
    	deallocate(P_sort)
    	
    
else if (trim(basis_order) == 'diagonal') then
    
    	! Sort according to magnitude on the diagonal!
    	do I_index = 1, F_basis_size
    		I1_index = special_functions_basis_size+I_index
    
    		kernel_to_be_sorted(I_index) =         &
    				abs(h_basis_kernel(I1_index,I1_index))
    	end do 
    
end if
    
! next, sort the array kernel_to_be_sorted; what matters here is the order!
allocate(P_sort(F_basis_size))
call HPSORT_real(F_basis_size,kernel_to_be_sorted,P_sort)
    
do I_index = 1,F_basis_size
   ! The sort is in ascending order!
   I1_index = P_sort(F_basis_size+1-I_index)+ special_functions_basis_size
   I2_index =                       I_index + special_functions_basis_size

   permutation(I1_index,I2_index) = cmplx_1
end do
    
deallocate(P_sort)

!================================================================================
!	Finally, permute the h basis to obtain the definitive basis,
!	which we'll call the t basis.
!================================================================================


allocate(t_basis_C_coefficients(ortho_basis_size,ortho_basis_size))

t_basis_C_coefficients(:,:) = matmul(h_basis_C_coefficients, permutation)

! make sure the t basis is still orthonormal
allocate(norm_checking_array(ortho_basis_size,ortho_basis_size))

norm_checking_array(:,:) = matmul(transpose(conjg(t_basis_C_coefficients)),t_basis_C_coefficients)


do I_index = 1, ortho_basis_size
   norm_checking_array(I_index,I_index) = norm_checking_array(I_index,I_index)-cmplx_1 
end do

real_norm = sqrt(real(sum(conjg(norm_checking_array(:,:))*norm_checking_array(:,:))))

write(*,31) '|       | <t_i| t_j> - delta_ij | = ',real_norm,' |'
write(*,20) '|           ---------------------------------       |'


deallocate(norm_checking_array)


! compute the kernel explicitely in this basis; this will be useful 
! afterwards to make sure the kernel behaves as expected. There is a lot
! of linear algebra in the code above; it will be important to make sure
! no mistakes were made!


allocate(t_basis_kernel(ortho_basis_size,ortho_basis_size))
allocate(t_basis_kernel_xc(ortho_basis_size,ortho_basis_size))
allocate(t_basis_kernel_Ha(ortho_basis_size,ortho_basis_size))

t_basis_kernel(:,:) = matmul(transpose(conjg(t_basis_C_coefficients)),   &
				matmul(C_basis_kernel,t_basis_C_coefficients))

t_basis_kernel_xc(:,:) = matmul(transpose(conjg(t_basis_C_coefficients)),   &
				matmul(C_basis_kernel_xc,t_basis_C_coefficients))

t_basis_kernel_Ha(:,:) = matmul(transpose(conjg(t_basis_C_coefficients)),   &
				matmul(C_basis_kernel_Ha,t_basis_C_coefficients))



! output the full t_basis_kernel, for testing and checking purposes.

call output_full_projected_kernel_netcdf(qvec,ortho_basis_size,t_basis_kernel)

!================================================================================
!	Create the functional projectors onto the sub-basis containing
!	optimal_basis_size functions.
!================================================================================

! The functional projectors are given by
!		T_{IL} = < t_I | B_L > 
!		       = sum_{I1} t^dagger_{II1} < C_I1 | B_L > 
!		       = sum_{I1} t^dagger_{II1} sum_{L1} Pc^dagger_{I1L1} < B_L1 | B_L > 
! thus
!		T_{IL} = sum_{I1,L1} t^dagger_{II1} Pc^dagger_{I1L1} O_{L1L}


allocate(Pt(bare_basis_size,optimal_basis_size))

! truncate the t basis to the specified size
Pt(:,:) = matmul(Pc,t_basis_C_coefficients(:,1:optimal_basis_size))


!allocate(norm_checking_array(optimal_basis_size,optimal_basis_size))
!norm_checking_array(:,:) = matmul(transpose(conjg(Pt)),matmul(overlaps,Pt))
!do I_index = 1, optimal_basis_size
!	norm_checking_array(I_index,I_index) = norm_checking_array(I_index,I_index) -cmplx_1
!end do
!real_norm = sqrt(real(sum(conjg(norm_checking_array(:,:))*norm_checking_array(:,:))))
!write(*,31) '|        | t^daqgger*O*t- identity |  = ',real_norm,'     |'
!deallocate(norm_checking_array)


!        Finally, compute the functional overlaps, which will be used
!	to compute chi_KS
allocate(T_functional_projectors(optimal_basis_size,bare_basis_size))

T_functional_projectors(:,:) = matmul(transpose(conjg(Pt)),overlaps)


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


! project the special functions onto the | t > basis
! The functions are ordered as { 1, mz, Delta, [plane waves]}. 
!
!	We have | g_i > = sum_I | t_I> g_{Ii}, such that
!	         g_{Ii} =  < t_I | g_i >. 
!
!	We have already computed < B_L | g_i > , which is given by
!
!		< B_L | g_i > = sum_I < B_L | t_I> g_{Ii}
!		< B_L | g_i > = sum_I T^dagger_{LI} g_{Ii}
!
!	Using the singular value decomposition, 
!
!		T^dagger = U * Sigma * V^dagger
!
!		pseudo_interse_T^dagger = V * Sigma^{-1} * U^dagger




! first, build the singular value decomposition of T^dagger
write(*,20) '|    Computing the pseudo-inverse of A= T^dagger... |'

allocate(A(bare_basis_size,optimal_basis_size))

A = transpose(conjg(T_functional_projectors))


!  singular value decomposition: 
!			A = U * S * V^dagger 

lwork = 2*optimal_basis_size+bare_basis_size
lrwork = 5*optimal_basis_size

allocate(Sigma_A(optimal_basis_size))
allocate(U_A(bare_basis_size,bare_basis_size))
allocate(V_dagger_A(optimal_basis_size,optimal_basis_size))
allocate(work(lwork))
allocate(rwork(lrwork))

call zgesvd(                       & ! Driver routine to compute SVD
		'A',               & ! compute all of the matrix U
		'A',               & ! compute all of the matrix V^dagger
    bare_basis_size,               & ! number of rows of A
 optimal_basis_size,               & ! number of columns of A
                  A,               & ! the A matrix on input, destroyed on output
    bare_basis_size,               & ! first dimension of A
            Sigma_A,               & ! singular values of A on output
                U_A,               & ! U matrix on output
    bare_basis_size,               & ! leading dimension of U
         V_dagger_A,               & ! V^dagger matrix on output
 optimal_basis_size,               & ! leading dimension of V^dagger
	 work,lwork,               & ! work array and dimensions
	      rwork,               & ! real work array
               info)


! Remember! zgesvd destroys the A matrix; define it again for convenience ...
A = transpose(conjg(T_functional_projectors))


write(*,22) '|          SVD info :  ',info,'                         |'

allocate(norm_checking_array(bare_basis_size,optimal_basis_size))
allocate(matrix_Sigma_A(bare_basis_size,optimal_basis_size))


matrix_Sigma_A(:,:) = cmplx_0

do I_index =1, optimal_basis_size
	matrix_Sigma_A(I_index,I_index) = cmplx_1*Sigma_A(I_index)
end do 

norm_checking_array(:,:) = matmul(U_A,matmul(matrix_Sigma_A,V_dagger_A))- A(:,:)

real_norm = sqrt(real(sum(conjg(norm_checking_array(:,:))*norm_checking_array(:,:))))



write(*,31) '|        | A- U*S*V^dagger |  = ',real_norm,'     |'
deallocate(norm_checking_array)
deallocate(matrix_Sigma_A)



! compute the pseudo inverse

allocate(pseudo_inverse_A(optimal_basis_size,bare_basis_size))
allocate(inverse_Sigma_A(optimal_basis_size,bare_basis_size))

inverse_Sigma_A(:,:) = cmplx_0

do I_index = 1, optimal_basis_size
	inverse_Sigma_A(I_index,I_index) = cmplx_1/Sigma_A(I_index)
end do


pseudo_inverse_A = matmul(transpose(conjg(V_dagger_A)),      &
		   matmul(inverse_Sigma_A,                   &
		          transpose(conjg(U_A))))


allocate(norm_checking_array(optimal_basis_size,optimal_basis_size))

norm_checking_array(:,:)= matmul(pseudo_inverse_A ,A)

do I_index = 1, optimal_basis_size
	norm_checking_array(I_index,I_index) = norm_checking_array(I_index,I_index)-cmplx_1 
end do

real_norm = sqrt(real(sum(conjg(norm_checking_array(:,:))*norm_checking_array(:,:))))

write(*,31) '|       | A^-1*A - identity | = ',real_norm,'     |'
write(*,20) '|           ---------------------------------       |'

deallocate(A)
deallocate(Sigma_A)
deallocate(inverse_Sigma_A)
deallocate(U_A)
deallocate(V_dagger_A)
deallocate(norm_checking_array)
deallocate(work)
deallocate(rwork)



! compute the coefficients in the t basis!

allocate(special_functions_t_coefficients(optimal_basis_size,special_functions_basis_size))

!special_functions_t_coefficients = matmul(pseudo_inverse_A,special_functions_B_overlaps)

special_functions_t_coefficients = matmul(transpose(conjg(Pt)),special_functions_B_overlaps)

deallocate(pseudo_inverse_A )


!================================================================================
!        Output functional projectors to netcdf file
!================================================================================

! truncate the kernel to the right basis size

allocate(optimal_basis_kernel(optimal_basis_size,optimal_basis_size))
allocate(optimal_basis_kernel_xc(optimal_basis_size,optimal_basis_size))
allocate(optimal_basis_kernel_Ha(optimal_basis_size,optimal_basis_size))



optimal_basis_kernel(:,:)    = t_basis_kernel(1:optimal_basis_size,1:optimal_basis_size)
optimal_basis_kernel_xc(:,:) = t_basis_kernel_xc(1:optimal_basis_size,1:optimal_basis_size)
optimal_basis_kernel_Ha(:,:) = t_basis_kernel_Ha(1:optimal_basis_size,1:optimal_basis_size)

write(*,20) '|       - Outputing the functional projectors to    |'
write(*,20) '|         netcdf file...                            |'


call output_functional_projectors_plasmons_netcdf(          &
	qvec,bare_basis_size,special_functions_basis_size,  &
	optimal_basis_size,T_functional_projectors,         &
	list_m1,list_m2,special_functions_t_coefficients,   &
	Pt,optimal_basis_kernel,optimal_basis_kernel_xc,optimal_basis_kernel_Ha)


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
write(io_unit,20) ' The eigenvalues of Og = < g | g > are given by '
write(io_unit,20) ''
do I_index = 1, specific_basis_size
	write(io_unit,19) Sg(I_index)
end do

write(io_unit,20) ''
write(io_unit,20) ' The special functions are represented as '
write(io_unit,20) '     | g_i > = sum_{I}  |t_I> g_{Ii}.     '
write(io_unit,20) ''
write(io_unit,20) ' The coefficients are given by:           '
write(io_unit,20) ''
write(io_unit,20) '# orthonormal basis  |                    special functions                   '

header_string = '#    | t_I >         |        one(r)        |' 
line_string   = '#============================================'


if (number_G_shells >1) then

     do iG = 2, nG_shell_max

	write(G_string,34) ' G = (',gvec(1,iG),',',     &
                                    gvec(2,iG),',',     &
                                    gvec(3,iG),')    |'

	header_string = trim(header_string)//trim(G_string)
	line_string   = trim(line_string)//trim('=======================')
     end do ! iG 

end if

write(io_unit,20) trim(header_string)
write(io_unit,20) trim(line_string)


do I_index = 1,optimal_basis_size
	write(io_unit,32) I_index,special_functions_t_coefficients(I_index,:)
end do


! build the g(r) functions explicitely and compare with the G(r) functions 

! first, extract their coefficients in the | B_L > basis
allocate(g_functions_B_coefficients(bare_basis_size,special_functions_basis_size))

g_functions_B_coefficients(:,:) = matmul(Pt, special_functions_t_coefficients)


! Next, compute the overlaps

allocate(eta_array(special_functions_basis_size))




eta_normalization = sum(conjg(one_r(:))*one_r(:))/dble(nr)
eta_array(1)      = sqrt(abs(                                       &
			sum(conjg(g_functions_B_coefficients(:,1))* &
		                special_functions_B_overlaps(:,1))  &
				                      /eta_normalization ))

if (number_G_shells >1) then
   do iG = 1, pw_dim
	eta_normalization = ONE
	eta_array(specific_basis_size +iG) = sqrt(abs(                                      &
			sum(conjg(g_functions_B_coefficients(:,specific_basis_size+iG))*   &
		                special_functions_B_overlaps(:,specific_basis_size+iG))    &
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

if (number_G_shells >1) then

     do iG = 2, nG_shell_max
	write(io_unit,34) ' G = (',gvec(1,iG),',',     &
                                   gvec(2,iG),',',     &
                                   gvec(3,iG),')',     &
				   eta_array(specific_basis_size+iG-1)
     end do ! iG 
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

end program generate_optimal_basis_plasmons
