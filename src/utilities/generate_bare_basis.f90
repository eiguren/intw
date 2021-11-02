!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       program generate_bare_basis
!       -----------------------------------
!
!       This is a "sub program" which is part of the intw project.
!
!       The purpose of this program is to compute and store the "bare basis"
!	functions which are products of Wannier functions.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program generate_bare_basis
  
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
real(dp) :: read_wannier_time, compute_B_time, write_B_time
real(dp) :: read_B_time, compute_overlap_time
real(dp) :: FFT_B_time


logical  :: use_IBZ
logical  :: W_file_exists
logical  :: fxc_file_exists

logical  :: compute_Coulomb_kernel

! netcdf file units
integer  :: B_nc_id, W_nc_id

! file units
integer  :: B_FFT_io_unit

integer  :: record_length 

!--------------------------------
! the supercell
!--------------------------------
! dimensions of the supercell
integer    :: snr1, snr2, snr3, snr

integer    :: sir1, sir2, sir3, sir
integer    ::  ir1,  ir2, ir3, ir
integer    ::  nr
integer    ::  iG

integer    ::  switch_singlet_to_triplet
integer    ::  switch_triplet_to_singlet

! r-vectors in the supercell
real(dp),allocatable :: rmesh_SC(:,:) ! supercell, centered at origin
real(dp),allocatable :: rmesh_WS(:,:) ! single cell, Fourier like, not centered
real(dp)             :: list_r_WS(3,8)
real(dp)             :: r_SC(3), r_WS(3), r_UC(3) 
integer              :: Rlat(3)


! the q vector for the basis functions

real(dp) :: qvec(3)

real(dp) :: q_plus_G(3), q_plus_G_cart(3)


integer  :: sparse_array_size, sparse_array_index

complex(dp),allocatable :: e_iqr(:)

integer,allocatable     :: sparse_row_indices(:),sparse_col_indices(:)


integer  ::  ipol_up, ipol_down
integer  ::  nb_W1, nb_W2, L_index
integer  ::  L1, L2
integer  ::  bare_basis_size

complex(dp),allocatable    :: wannier_function1(:,:)
complex(dp),allocatable    :: wannier_function2(:,:)
complex(dp),allocatable    :: wannier_product(:)
complex(dp),allocatable    :: B_function(:)
complex(dp),allocatable    :: B_function1(:),B_function2(:)
complex(dp),allocatable    :: B_G1(:),B_G2(:)
complex(dp),allocatable    :: Coulomb_kernel(:)


complex(dp),allocatable :: overlaps(:,:)
complex(dp),allocatable :: kernel_xc(:,:), kernel_Ha(:,:)


! function in r, integrals with bare basis, coefficient in optimal basis
complex(dp),allocatable :: fxc_r(:), field_r(:)




real(dp)  :: volume
!--------------------------------------------------------------------------------


! start global clock
call get_timing(total_time1)

!================================================================================
!        Talk to user
!================================================================================
write(*,20) '====================================================='
write(*,20) '|          program generate_bare_basis              |'
write(*,20) '|        ---------------------------------          |'
write(*,20) '====================================================='
write(*,20) '|                                                   |'
write(*,20) '|        This program computes the product basis    |'
write(*,20) '|        using Maximally localized Wannier          |'
write(*,20) '|        functions.                                 |'
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
compute_Coulomb_kernel = .false.

if     (nspin==1) then
	write(*,20) '|       - The calculation is paramagnetic nspin=1   |'
        npol = 1
	compute_Coulomb_kernel = .true.

	if (magnon) then
	  magnon = .false.
	  write(*,20)'*         WARNING                                   *' 
	  write(*,20)'*         The input file contained magnon=.true.    *' 
	  write(*,20)'*         This is a paramagnetic system so the      *' 
	  write(*,20)'*         program sets magnon=.false.               *' 
	end if


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


!================================================================================
!        prepare the netcdf file which will contain the 
!	 bare product basis.
!================================================================================
switch_singlet_to_triplet  = -1
switch_triplet_to_singlet  =  1

qvec(1) = qpt_basis_1
qvec(2) = qpt_basis_2
qvec(3) = qpt_basis_3

nr      = nr1*nr2*nr3
allocate(rmesh_WS(3,nr))

do ir =1, nr

	call switch_indices(nr1,nr2,nr3,ir,ir1,ir2,ir3, &
					switch_singlet_to_triplet)

	   rmesh_WS(1,ir) =  dble(ir1-1)/dble(nr1)
	   rmesh_WS(2,ir) =  dble(ir2-1)/dble(nr2)
	   rmesh_WS(3,ir) =  dble(ir3-1)/dble(nr3)
end do


call prepare_B_functions_netcdf(qvec,rmesh_WS,nr,B_nc_id)


deallocate(rmesh_WS)
!================================================================================
!        Build the e_iqR array, which contains the phase factors
!	 relevant for the given q vector.
!================================================================================
write(*,20) '====================================================='
write(*,20) '|    Preparing real space arrays and exponential    |'
write(*,20) '|    factors ...                                    |'

call get_timing(time1)


! allocate and read the supercell real space mesh
allocate(rmesh_SC(3,snr))

call read_supercell_rmesh_from_wannier_functions_netcdf  &
             				(W_nc_id,snr,rmesh_SC)


!-------------------------------------------------------------------------------
! CAREFUL! There is great potential for confusion here.
!
!	Define the B functions as
!
!	B_{I=m1m2}(r) = Omega/N sum_{R} e^{-i q*(r+R)} W^*_m1(r+R) W_m2(r+R)
!
!	The Wannier functions are expressed on a SuperCell (SC) centered at the
!	origin; by approximation, they are taken to be identically zero
!	outside this supercell.
!
!	On the other hand, it is easy to see that B_I(r) is a periodic function.
!	Thus, it can safely be stored only in a single cell, which in crystal
!	coordinates can be defined as the domain r in [0,1]x[0,1]x[0,1] = Unit Cell (UC).
!
!
!	Hence, let's define three coordinate grids:
!
!		r_SC :  the Gamma-centered supercell coordinate, 
!			relevant to describe the Wannier functions
!
!		r_UC :  a coordinate in the off-centered (sometimes referred to as
!			FFT grid) unit cell coordinate
!
!	We have:
!		r_SC  = r_UC+ R,                R some lattice vector
!
!	Clearly, many different r_SC will contribute to the same r_UC (one for each R!).
!
!	Thus, define a Matrix,
!		e_iqr[ r_UC, r_SC] = e^{-iq*r_SC} , r_SC = r_UC+R for some R.
!
!	Then,
!		B[r_UC] = sum_{r_SC} e_iqR[r_UC,r_SC] [W^*_m1(r_SC)W_m2(r_SC)]
!
!	Clearly, the matrix e_iqr will be sparse; for a given r_UC, only
!	the components with r_SC such that  r_SC = r_UC+R will be non-zero.
!
!	The only task left is to determine the relationship bewteen r_SC and r_UC.
!
!-------------------------------------------------------------------------------

! build the sparse matrix

sparse_array_size = snr ! (this variable is only defined for clarity)

! the actual size of the arrays is unknown, but will be no larger than the maximum
allocate(             e_iqr(sparse_array_size))
allocate(sparse_row_indices(sparse_array_size))
allocate(sparse_col_indices(sparse_array_size))

sparse_array_index = 0

! simple test
!write(121,*) '#---------------------------------------------------------'
!write(121,*) '# r supercell      r Wigner Seitz     r unit cell         '
!write(121,*) '#---------------------------------------------------------'
do sir = 1,snr

   ! increment the sparse_array index
   sparse_array_index = sparse_array_index + 1

   ! compute the exponential contribution
   e_iqr(sparse_array_index) = exp(-cmplx_i*tpi*sum(r_SC(:)*qvec(:)))

   ! identify the column index
   sparse_col_indices(sparse_array_index) = sir

   ! getting the row index is more tricky!

   ! vector in supercell
   r_SC(:) = rmesh_SC(:,sir)

   ! r_UC = r_WS+Rlat
   call find_k_1BZ_and_G(r_SC,nr1,nr2,nr3,ir1,ir2,ir3,r_UC,Rlat)

   ! find singlet index for WS cell index
   call switch_indices(nr1,nr2,nr3,ir,ir1,ir2,ir3,switch_triplet_to_singlet)

   sparse_row_indices(sparse_array_index) =  ir


   ! simple test
!   write(121,'(3(3F8.4,4x))')  r_SC(:),r_WS(:),r_UC(:)


end do ! sir


call get_timing(time2)
time = time2-time1

write(*,10) '|                     - time::',time,' seconds    |'
write(*,20) '====================================================='

!================================================================================
!        Build the bare product functions
!================================================================================
L_index   = 0


! This is pertinent only in the case where spin is present
ipol_up   = 1
ipol_down = 2

allocate(wannier_function1(snr,nspin))
allocate(wannier_function2(snr,nspin))
allocate(  wannier_product(snr))
allocate(       B_function(nr))


volume = alat**3*abs(                                  &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))


write(*,20) '====================================================='
write(*,20) '|    Computing the bare basis functions, B_L(r):    |'
write(*,20) '====================================================='

bare_basis_size = num_wann**2

read_wannier_time  = ZERO
compute_B_time     = ZERO
write_B_time       = ZERO

call get_timing(loop_time1)
do nb_W1 = 1,num_wann

   call get_timing(time1)
   call read_wannier_functions_netcdf  &
                (W_nc_id,nb_W1,snr,wannier_function1)
   call get_timing(time2)
   read_wannier_time = read_wannier_time + time2-time1

   do nb_W2 = 1,num_wann

	L_index  = L_index + 1

        write(*,16) '| L_index = ',L_index,' of ',bare_basis_size,'                        |'

        call get_timing(time1)
	! Read the Wannier functions from the netcdf file
        call read_wannier_functions_netcdf  &
                (W_nc_id,nb_W2,snr,wannier_function2)
        call get_timing(time2)
        read_wannier_time = read_wannier_time + time2-time1

        call get_timing(time1)
        ! compute product of Wannier functions, for the Magnon case

	if (nspin.eq.1) then

	  wannier_product(:) = conjg(wannier_function1(:,1))*    &
				     wannier_function2(:,1)

	else if (magnon) then

	  wannier_product(:) = conjg(wannier_function1(:,ipol_up))*    &
				     wannier_function2(:,ipol_down)

	else
	
          write(*,20) '*****************************************************'
          write(*,20) '* ERROR: General spin case not yet implemented      *'
          write(*,20) '* Program stops!                                    *'
          write(*,20) '*****************************************************'

	  stop

	end if

        ! Compute the product y = A*x, where A is a sparse matrix
	!call mkl_zcoogemv(transa, m, val, rowind, colind, nnz, x, y)

	call mkl_zcoogemv('N',        & ! Normal order
	           	   nr,        & ! Number of rows in the A matrix
			e_iqr,        & ! The sparse matrix A
           sparse_row_indices,        & ! row indices where A!=0
           sparse_col_indices,        & ! col indices where A!=0
	    sparse_array_size,        & ! number of non-vanishing elements in A
              wannier_product,	      & ! the x array
		   B_function)          ! the y array

        ! normalize by the volume at the end
	B_function(:) = B_function(:)*volume
        call get_timing(time2)
        compute_B_time = compute_B_time + time2-time1

        call get_timing(time1)
	call write_B_functions_netcdf(B_nc_id,L_index,B_function,nr)
        call get_timing(time2)
        write_B_time = write_B_time + time2-time1




   end do ! nb_W2 
end do ! nb_W1
call get_timing(loop_time2)
loop_time = loop_time2-loop_time1


write(*,20) '====================================================='
write(*,18) '|  - reading Wannier       : ',read_wannier_time, ' sec.          |'
write(*,18) '|  - computing basis       : ',compute_B_time, ' sec.          |'
write(*,18) '|  - writing basis to file : ',write_B_time,   ' sec.          |'
write(*,20) '|---------------------------------------------------|'
write(*,18) '|  - total B function time : ',loop_time,' sec.          |'
write(*,20) '====================================================='


deallocate(wannier_function1)
deallocate(wannier_function2)
deallocate(  wannier_product)
deallocate(       B_function)


!close to flush!
call close_B_functions_netcdf(B_nc_id)
call close_W_functions_netcdf(W_nc_id)
!================================================================================
!        compute the overlaps
!================================================================================


write(*,20) '====================================================='
write(*,20) '|    Computing the overlap matrix,                  |'
write(*,20) '|          O_{L1L2} = < B_L1 | B_L2 >               |'
write(*,20) '|    and the kernel matrix,                         |'
if (compute_Coulomb_kernel) then
	write(*,20) '|     kernel_xc_{L1L2} = < B_L1 | fxc | B_L2 >      |'
	write(*,20) '|     kernel_Ha_{L1L2} = < B_L1 |  vc | B_L2 >      |'
	write(*,20) '|            kernel = kernel_xc + kernel_Ha         |'
else
	write(*,20) '|     kernel_{L1L2} = < B_L1 | fxc | B_L2 >         |'
end if
write(*,20) '====================================================='


allocate(fxc_r(nr))

! In magnon case, the exchange-correlation kernel
! is read in real space from existing file.

if (magnon) then
   allocate(field_r(nr))

   ! check that the file exists

   inquire(file=fxc_r_filename, exist=fxc_file_exists)

   if (.not.  fxc_file_exists) then
      write(*,20) '*******************************************************'
      write(*,20) '*   ERROR: The file which should contain the XC       *'
      write(*,20) '*          kernel does not exist!                     *'
      write(*,20) '*               program stops.                        *'
      write(*,20) '*******************************************************'
      stop
   else
      write(*,20) '|       - The files containing fxc(r) has been      |'
      write(*,20) '|         found.                                    |'
      write(*,20) '|           ---------------------------------       |'
   end if

   ! read in the data from the files
   ! CAREFUL!!!
   !
   !       The data in the ascii files was generated using a hack from
   !       Quantum Espresso; it is not very clear, but the most likely
   !       order for the r-mesh is Fortran style (ie, the first dimension loops
   !       fastest, "column-major").
   !       HOWEVER, in this code the real space functions are stored using C style,
   !       where the last dimension loops fastest "row-major". There is no 
   !       good reason for this (I just started doing it this way, I can't remember why);
   !       it is simpler to account for it here than to modify everything else in the code.


   call read_field_ascii_file(fxc_r_filename,nr,field_r)
   call change_order(nr,field_r,fxc_r)

   deallocate(field_r)

   ! In density calculation, the Coulomb part of 
   ! the xc kernel will be computed in k space.

else if ( nspin == 1 ) then

   ! It might be useful to introduce an RPA switch, which 
   ! allows the user to use RPA instead of LDA
   if (RPA) then 
      fxc_r(:) = cmplx_0 
   else
      inquire(file=fxc_r_filename, exist=fxc_file_exists)

      if (.not.  fxc_file_exists) then
         write(*,20) '*******************************************************'
         write(*,20) '*   ERROR: The file which should contain the XC       *'
         write(*,20) '*          kernel does not exist!                     *'
         write(*,20) '*               program stops.                        *'
         write(*,20) '*******************************************************'
         stop
      else
         write(*,20) '|       - The files containing fxc(r) has been      |'
         write(*,20) '|         found.                                    |'
         write(*,20) '|           ---------------------------------       |'
      end if

      allocate(field_r(nr))

      call read_field_ascii_file(fxc_r_filename,nr,field_r)
      call change_order(nr,field_r,fxc_r)
      deallocate(field_r)
   end if

end if


! open to read (but not to write!)
call open_B_functions_netcdf(B_nc_id)

allocate(B_function1(nr))
allocate(B_function2(nr))
allocate(overlaps(bare_basis_size,bare_basis_size))

allocate(kernel_xc(bare_basis_size,bare_basis_size))

kernel_xc(:,:) = cmplx_0

if (compute_Coulomb_kernel) then
   write(*,20) '====================================================='
   write(*,20) '|    Computing the Coulomb kernel and FFT of        |'
   write(*,20) '|    bare basis functions...                        |'
   !----------------------------------------------------------------
   ! Computing the Coulomb kernel is more efficient in G space. 
   ! Indeed:
   !		v_q(r,r') = \sum_{G} e^{iG(r-r')} 4pi e^2/|q+G|^2
   !
   !	< Bi_q| v_q | Bj_q > = sum_{G} Bi_q(G)^* Bj_q(G) 4pi e^2/|q+G|^2
   !
   ! The functions must thus be Fourier transformed to G space
   !
   !----------------------------------------------------------------

   allocate(B_G1(ngm))
   allocate(B_G2(ngm))
   allocate(Coulomb_kernel(ngm))
   allocate(kernel_Ha(bare_basis_size,bare_basis_size))

   kernel_Ha(:,:) = cmplx_0

   ! compute the Coulomb kernel
   do iG = 1, ngm 
      q_plus_G(:) = qvec(:)+gvec(:,iG)

      q_plus_G_cart(:) = 2.0_dp*pi/alat*(           &
                         q_plus_G(1)*bg(:,1)+       &
                         q_plus_G(2)*bg(:,2)+       &
                         q_plus_G(3)*bg(:,3))

      Coulomb_kernel(iG) = 4.0_dp*pi/sum(q_plus_G_cart(:)**2)*Ha_in_eV


   end do !iG
  

   read_B_time  = ZERO
   FFT_B_time   = ZERO
   write_B_time = ZERO

   ! find a free unit
  B_FFT_io_unit = find_free_unit()


  ! open a temporary file which will hold the Fourier transformed B functions
  record_length = ngm*direct_io_factor_cmplx

  open(unit = B_FFT_io_unit, form = 'unformatted',                  &
     status = 'scratch', access = 'direct', recl = record_length,   &
     action = 'readwrite')

   do L1 = 1, bare_basis_size

      ! compute all the Fourier transforms, and store to scratch file
      call get_timing(time1)
      call read_B_function_netcdf(B_nc_id,L1,B_function1,nr)
      call get_timing(time2)
      read_B_time = read_B_time + time2-time1

      call get_timing(time1)
      call func_from_r_to_g (1, ngm, B_function1, B_G1)
      call get_timing(time2)
      FFT_B_time = FFT_B_time + time2-time1

      call get_timing(time1)
      write(B_FFT_io_unit,rec=L1)  B_G1(:)
      call get_timing(time2)
      write_B_time = write_B_time + time2-time1

   end do

   write(*,18) '|  - read B functions (nc) : ',read_B_time, ' sec.          |'
   write(*,18) '|  - FFT B functions       : ',FFT_B_time,' sec.          |'
   write(*,18) '|  - write B functions(bin): ',write_B_time,' sec.          |'
   write(*,20) '====================================================='


end if

do L1 = 1, bare_basis_size
  write(*,20) '|---------------------------------------------------|'
  call get_timing(loop_time1)
  write(*,23) '| L1  = ',L1,'                                      |'

  read_B_time           = ZERO
  FFT_B_time            = ZERO
  compute_overlap_time  = ZERO


  call get_timing(time1)
  call read_B_function_netcdf(B_nc_id,L1,B_function1,nr)
  call get_timing(time2)
  read_B_time = read_B_time + time2-time1

  if (compute_Coulomb_kernel) then
     call get_timing(time1)
     read(B_FFT_io_unit,rec=L1)  B_G1(:)
     call get_timing(time2)
     FFT_B_time = FFT_B_time + time2-time1
  end if

   do L2 = 1, bare_basis_size

        call get_timing(time1)
        call read_B_function_netcdf(B_nc_id,L2,B_function2,nr)
        call get_timing(time2)
        read_B_time = read_B_time + time2-time1
  
        if (compute_Coulomb_kernel)  then
           call get_timing(time1)
           read(B_FFT_io_unit,rec=L2)  B_G2(:)
           call get_timing(time2)
           FFT_B_time = FFT_B_time + time2-time1
        end if

        call get_timing(time1)
	overlaps(L1,L2) = sum(conjg(B_function1(:))*B_function2(:))/dble(nr)


        ! First, the contribution of the exchange-correlation kernel
	kernel_xc(L1,L2) = sum(conjg(B_function1(:))*fxc_r(:)*B_function2(:))/dble(nr)
        if (compute_Coulomb_kernel) then 
	   kernel_Ha(L1,L2)= sum(conjg(B_G1(:))*Coulomb_kernel(:)*B_G2(:))
        end if


        call get_timing(time2)
        compute_overlap_time = compute_overlap_time + time2-time1


   end do ! L2
   call get_timing(loop_time2)
   loop_time = loop_time2-loop_time1

   write(*,18) '|  - reading B function    : ',read_B_time, ' sec.          |'
   if (compute_Coulomb_kernel)  then
       write(*,18) '|  - reading FFT B function: ',FFT_B_time, ' sec.          |'
   end if
   write(*,18) '|  - computing overlaps    : ',compute_overlap_time, ' sec.          |'
   write(*,18) '|  - total loop time       : ',loop_time,' sec.          |'
end do ! L1


if (compute_Coulomb_kernel) then
   deallocate(B_G1)
   deallocate(B_G2)
   deallocate(Coulomb_kernel)
  close(unit = B_FFT_io_unit)
end if



write(*,20) '====================================================='
write(*,20) '|    write overlaps to netcdf file...               |'
call get_timing(time1)
if (compute_Coulomb_kernel) then
   call output_overlaps_plasmons_netcdf(overlaps,kernel_xc,kernel_Ha,qvec,bare_basis_size)
else
   call output_overlaps_netcdf(overlaps,kernel_xc,qvec,bare_basis_size)
end if
call get_timing(time2)

write(*,10) '|                     - time::',time,' seconds    |'


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
18 format(A,ES8.2,A)
20 format(A)
22 format(A,I4,A)
23 format(A,I6,A)
30 format(A,F10.3,A)

end program generate_bare_basis
