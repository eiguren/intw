!       program generate_wannier_functions
!       -----------------------------------
!
!       This is a utility program which is part of the intw project.
!
!       The purpose of this program is to compute the Wannier functions 
!	on a real space mesh:
!
!		W^sigma_{n}(r) = sqrt{N/V} 1/N sum_{mk} U_{mn}(k) e^{ikr} u^sigma_{mk}(r)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program generate_wannier_functions

  use intw_useful_constants
  use intw_setup

  use intw_input_parameters
  use intw_reading
  use intw_symmetries
  use intw_dft_fftw
  use intw_W90
  use intw_plane_wave_setup
  use intw_matrix_elements
  use intw_product_basis

!Peio
!Number of Bloch Original states for the wannierization
  use w90_parameters, only: num_bands
!Peio

!================================================================================
!       Declare the variables 
!================================================================================
implicit none


integer       :: ikpt
integer       :: ir1, ir2, ir3
integer       :: sir1, sir2, sir3, sir
integer       :: ik1, ik2, ik3
integer       :: switch_triplet_to_singlet,switch_singlet_to_triplet

integer       :: snr1, snr2, snr3, snr
integer       :: ir, nr

integer       :: W_io_unit
integer       :: record_length , record_index

integer       :: nb, nb_W, ipol

real(dp)      :: volume
complex(dp)   :: normalization

complex(dp),allocatable    :: fr(:,:,:), fg(:,:,:)

real(dp) ::      rvec(3), rvec_in_WS(3)
integer  ::      Rlat(3)

integer, allocatable    :: list_ir_of_sir(:)

! keeping time
real(dp) ::      setup_time

real(dp) ::      total_time1, total_time2, total_time

real(dp) ::      loop_time1, loop_time2, loop_time
real(dp) ::      exp_time1, exp_time2, exp_time
real(dp) ::      psi_time1, psi_time2, psi_time
real(dp) ::      d_time1, d_time2, d_time
real(dp) ::      contract_time1, contract_time2, contract_time
real(dp) ::      fft_time1, fft_time2, fft_time
real(dp) ::      ud_time1, ud_time2, ud_time
real(dp) ::      sum_time


logical        ::      use_IBZ



real(dp),allocatable       ::  rmesh(:,:)
complex(dp),allocatable    ::  e_ikr(:)

real(dp)  ::  kpt(3)

! wavefunction variables
complex(dp),allocatable    :: wfc_k(:,:,:)
complex(dp),allocatable    :: wfc_r(:,:,:)

integer,allocatable        :: list_iG(:)
real(dp),allocatable       :: QE_eig(:)
complex(dp),allocatable    :: wfc_k_3D(:,:,:,:,:)

complex(dp),allocatable    :: Uk_wfc_r(:,:,:)

! Wannier variables
complex(dp),allocatable    :: Wannier_functions(:,:,:)
complex(dp),allocatable    :: Wf(:)
complex(dp),allocatable    :: U_k(:,:)


! block variables
integer              :: number_of_blocks, i_block, rest_block, nr_per_block
integer,allocatable  :: block_bounds(:,:)
integer              :: sir_min, sir_max, nr_block_size
integer              :: sir_loop


character    ::  matdescra(6) 

integer, allocatable     :: sparse_rows(:), sparse_columns(:)
integer                  :: sparse_size

!Peio
!The variable we use instead of num_bands
integer :: nbands_loc
!Peio

nbands_loc=num_bands

call get_timing(total_time1)


!================================================================================
!        Talk to user
!================================================================================

  write(*,20) '====================================================='
  write(*,20) '|         program generate_wannier_functions        |'
  write(*,20) '|  ---------------------------------------------    |'
  write(*,20) '====================================================='
  write(*,20) '|                                                   |'
  write(*,20) '|  This sub-program produces the Wannier functions  |'
  write(*,20) '|  on the real space mesh.                          |'
  write(*,20) '|                                                   |'
  write(*,20) '====================================================='



switch_triplet_to_singlet =  1
switch_singlet_to_triplet = -1

matdescra(1) = 'G' ! general matrix
matdescra(4) = 'F' ! fortran indexing


!================================================================================
!        call "setup" routine, which sets up the necessary environment.
!================================================================================

call setup(setup_time)
! ham_r, from W90, is no longer needed.
call deallocate_hr_W90()

volume = alat**3*abs(                                  &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))

normalization = cmplx_1/sqrt(volume)/(nk1*nk2*nk3)


write(*,20) '====================================================='
write(*,10) '|            setup    - time::',setup_time,' seconds    |'
write(*,20) '====================================================='



use_IBZ = .not. intw2W_fullzone

!================================================================================
!        Tell the user if the system has spin
!================================================================================

if     (nspin==1) then
	write(*,20) '|       - The calculation is paramagnetic nspin=1   |'
        npol = 1
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
!   prepare the arrays to build the Wannier functions
!================================================================================

call setup_fftw_3D(nr1,nr2,nr3)

nr   = nr1*nr2*nr3

snr1 = super_factor_1*nr1+1
snr2 = super_factor_2*nr2+1
snr3 = super_factor_3*nr3+1

snr  = snr1*snr2*snr3

! determine the size of the sub-blocks of the dense r-mesh
call split_rmesh_in_subblocks_intw(number_of_blocks, nr, snr ,rest_block,nr_per_block)

! build the blocks

allocate(block_bounds(number_of_blocks,2))


block_bounds(:,:) = 0

sir_min = 1

do i_block = 1, number_of_blocks

   block_bounds(i_block,1) = sir_min

   if ( rest_block == 0) then
        sir_max     = sir_min+nr_per_block-1
   else
        sir_max     = sir_min+nr_per_block
        rest_block  = rest_block -1
   end if

   if ( sir_max < snr) then
        block_bounds(i_block,2) = sir_max
   else if ( sir_max >= snr) then
        block_bounds(i_block,2) = snr
        exit
   end if
   sir_min = sir_max+1
end do


! tell user about RAM management

  write(*,20) '|                                                   |'
  write(*,20) '|  Given that the calculation uses a lot of RAM,    |'
  write(*,20) '|  the real-space mesh is split into blocks.        |'
  write(*,20) '|  Summary of parameters:                           |'
  write(*,20) '|                                                   |'
  write(*,11) '|             nr1 = ',nr1,'                        |'
  write(*,11) '|             nr2 = ',nr2,'                        |'
  write(*,11) '|             nr3 = ',nr3,'                        |'
  write(*,11) '|    ===>     nr  = ',nr,'                        |'
  write(*,20) '|                                                   |'
  write(*,11) '|            snr1 = ',snr1,'                        |'
  write(*,11) '|            snr2 = ',snr2,'                        |'
  write(*,11) '|            snr3 = ',snr3,'                        |'
  write(*,11) '|    ===>    snr  = ',snr,'                        |'
  write(*,20) '|                                                   |'
  write(*,18) '|         max_ram            = ',max_ram,' bytes       |'
  write(*,11) '|       number_of_blocks     = ',number_of_blocks,'             |'
  write(*,11) '|         nr_per_block(+/-1) = ',nr_per_block,'             |'
  write(*,20) '|    block bounds:                                  |'
  write(*,20) '|    i block    sir_min            sir_max          |'
  do i_block = 1, number_of_blocks
     write(*,21) i_block, block_bounds(i_block,1), block_bounds(i_block,2) 
  end do
  write(*,20) '|                                                   |'
  write(*,20) '====================================================='



! allocate some global arrays
allocate(wfc_k(nG_max,nbands_loc,nspin))
allocate(wfc_r(nr,nbands_loc,nspin))
allocate(wfc_k_3D(nr1,nr2,nr3,nbands_loc,nspin))
!allocate(wfc_k(nG_max,nbands,nspin))
!allocate(wfc_r(nr,nbands,nspin))
!allocate(wfc_k_3D(nr1,nr2,nr3,nbands,nspin))
allocate(Uk_wfc_r(nr,num_wann,nspin))
allocate(QE_eig(nbands_loc))
!allocate(QE_eig(nbands))
allocate(list_iG(nG_max))
allocate(U_k(nbands_loc,num_wann))
!allocate(U_k(nbands,num_wann))
allocate(fr(nr1,nr2,nr3))
allocate(fg(nr1,nr2,nr3))



! find a free unit
W_io_unit = find_free_unit()

! open a temporary file which will hold the Wannier functions

! BE VERY CAREFUL HERE. IS THE RECORD UNIT OF LENGTH A BYTE
! OR A LONGWORD, WHICH IS 4 BYTES???
! THIS IS SET BY COMPILER OPTIONS


record_length = snr*direct_io_factor_cmplx

!record_length = snr*double_complex

! open a scratch file which will contain the Wannier data;
! this doesn't fit in RAM!

open(unit = W_io_unit ,        form = 'unformatted', &
   status = 'scratch',       access = 'direct',      &
     recl = record_length,   action = 'readwrite')

! initialize
allocate(Wf(snr))

Wf(:) = cmplx_0

record_index = 0
do ipol =1, nspin
  do nb_W = 1,num_wann

	record_index = record_index + 1 
	write(W_io_unit,rec=record_index) Wf

  end do
end do


deallocate(Wf)


!================================================================================
!   Loop on the blocks of the super-cell r mesh
!================================================================================


do i_block = 1, number_of_blocks
        write(*,20) '|===================================================|'
        write(*,16) '| block ',i_block,' of ',number_of_blocks,'                            |'
        write(*,20) '|===================================================|'

   	sir_min = block_bounds(i_block,1)
   	sir_max = block_bounds(i_block,2)

	nr_block_size = sir_max-sir_min+1


        ! allocate the needed arrays
	allocate(list_ir_of_sir(nr_block_size))
	allocate(Wannier_functions(nr_block_size,num_wann,nspin))
	allocate(e_ikr(nr_block_size))
	allocate(rmesh(3,nr_block_size))

        Wannier_functions(:,:,:) = cmplx_0

        ! generate r-mesh, which is a super cell of the Wigner seitz 
        ! cell, and which is centered at (0,0,0).
        ! Also build a list tracking indices.

	do sir_loop =1, nr_block_size

	   sir = sir_min + sir_loop-1

           call switch_indices(snr1,snr2,snr3,sir,sir1,sir2,sir3, &
					switch_singlet_to_triplet)

   	   rvec(1) = dble(sir1-snr1/2-1)/dble(nr1)
           rvec(2) = dble(sir2-snr2/2-1)/dble(nr2)
           rvec(3) = dble(sir3-snr3/2-1)/dble(nr3)

           rmesh(:,sir_loop) = rvec(:)


   	   call find_k_1BZ_and_G(rvec,nr1,nr2,nr3,ir1,ir2,ir3,rvec_in_WS,Rlat)


   	   call switch_indices(nr1,nr2,nr3,ir,ir1,ir2,ir3,switch_triplet_to_singlet)


   	   list_ir_of_sir(sir_loop) = ir

        end do

        ! loop on k points in the coarse mesh
        do ikpt = 1, nk1*nk2*nk3
           call get_timing(loop_time1)

           write(*,16) '| ik = ',ikpt,' of ',nk1*nk2*nk3,'                             |'

           call switch_indices(nk1,nk2,nk3,ikpt,ik1,ik2,ik3,switch_singlet_to_triplet)

           ! k-point in crystal coordinates
           kpt(1) = dble(ik1-1)/dble(nk1)
           kpt(2) = dble(ik2-1)/dble(nk2)
           kpt(3) = dble(ik3-1)/dble(nk3)

           ! compute the exponential
           call get_timing(exp_time1)
           e_ikr(:) = exp(  twopi*cmplx_i*(	        &
	       		     kpt(1)*rmesh(1,:) + 	&
			     kpt(2)*rmesh(2,:) + 	&
			     kpt(3)*rmesh(3,:) ))
           call get_timing(exp_time2)
           exp_time  = exp_time2 -exp_time1

           ! extract the U matrix
           U_k(:,:) = u_mesh(:,:,ikpt)


           ! Get wavefunction from the QE folders
           call get_timing(psi_time1)
           call get_psi_k(ikpt,use_IBZ,list_iG,wfc_k,QE_eig)
           call get_timing(psi_time2)
           psi_time = psi_time2 -psi_time1


           !---------------------------------------------------------------------
           ! Fourier transform wavefunction to real space
           !---------------------------------------------------------------------

           !--- Code below uses FFTW, which is fast, but requires going to a 3D mesh ----

           ! put wavefunction on a 3D G-space mesh
           call get_timing(d_time1)
           call wfc_G_from_1D_to_3D (list_iG,wfc_k,wfc_k_3D)
           call get_timing(d_time2)
           d_time   = d_time2 -d_time1

           call get_timing(fft_time1)
           do ipol = 1, nspin
              do nb = 1, nbands_loc
!              do nb = 1, nbands

                  fg(:,:,:) = wfc_k_3D(:,:,:,nb,ipol)

                  call func_from_g_to_r_3D_FFTW(nr1,nr2,nr3,fr,fg)

                  do ir = 1, nr

          		call switch_indices(nr1,nr2,nr3,ir,ir1,ir2,ir3,&
						switch_singlet_to_triplet)

                        wfc_r(ir,nb,ipol) = fr(ir1,ir2,ir3)
                  end do ! ir

              end do ! nb
           end do ! ipol
           call get_timing(fft_time2)
           fft_time = fft_time2-fft_time1


           !---------------------------------------------------------------------
           ! contract the wave functions with the Wannier rotation matrices
           !---------------------------------------------------------------------

           call get_timing(contract_time1)
           Uk_wfc_r(:,:,:) = cmplx_0

           do ipol = 1,nspin
              ! matrix-matrix product
              ! C = alpha*A*B + beta*C
              ! with
              !  A = wfc_r(:,:,ipol)    of dimensions nr     x nbands
              !  B = Uk(:,:,ipol)       of dimensions nbands x num_wann
              !  C = Uk_wfc_r(:,:,ipol) of dimensions nr x num_wann

	      call  zgemm('N',    & ! A is in the normal order 
                          'N',    & ! B is in the normal order 
                           nr,    & ! first dimension of A
                     num_wann,    & ! second dimension of B and C
                   nbands_loc,    & ! second dimension of A
!                       nbands,    & ! second dimension of A
                      cmplx_1,    & ! alpha
              wfc_r(:,:,ipol),    & ! A matrix
                           nr,    & ! first dimension of A
                     U_k(:,:),    & ! B matrix
                   nbands_loc,    & ! first dimension of B
!                       nbands,    & ! first dimension of B
                      cmplx_0,    & ! beta
           Uk_wfc_r(:,:,ipol),    & ! C matrix
  		          nr)       ! first dimension of C


           end do 

           call get_timing(contract_time2)
           contract_time = contract_time2- contract_time1

           !---------------------------------------------------------------------
           ! build the Wannier functions IN THIS SUB-BLOCK
           !---------------------------------------------------------------------

           ! prepare the Sparse matrix for updating the Wannier function

           call get_timing(ud_time1)

           sparse_size = nr_block_size
	   allocate(sparse_rows(sparse_size))
	   allocate(sparse_columns(sparse_size))

	   do sir_loop =1, nr_block_size
                ir = list_ir_of_sir(sir_loop)

		sparse_rows(sir_loop)    = sir_loop
		sparse_columns(sir_loop) = ir
           end do


           do ipol =1, nspin
              ! sparse matrix-matrix product
              ! C = alpha*A*B + beta*C
              ! with
              !  A = the sparse matrix which contains e^ikr and which 
              !      sorts out the relationship between sir and ir.
              !      A has dimensions nr_block_size x nr
              !  B = Uk_wfc_r(:,:,ipol) of dimensions nr x num_wann
              !  C = Wannier_functions(:,:,ipol) of dimensions nr_block_size x num_wann
              

	     call mkl_zcoomm(        'N',        & ! normal order
                           nr_block_size,        & ! first dimension of A matrix
                                num_wann,        & ! second dimension of C matrix
                                      nr,        & ! second dimension of A matrix
                                 cmplx_1,        & ! alpha = 1
                               matdescra,        & ! info about matrix
                                   e_ikr,        & ! values of sparse matrix.
                             sparse_rows,        & ! row indices, COO format
                          sparse_columns,        & ! column indices, COO format
                             sparse_size,        & ! nnz variable, number of non-zero elements
		      Uk_wfc_r(:,:,ipol),        & ! B matrix
                                      nr,        & ! first dimension of B matrix
                                 cmplx_1,        & ! beta = 1; we are updating C!!
             Wannier_functions(:,:,ipol),        & ! C matrix
                           nr_block_size)          ! first dimension of C matrix


           end do

	   deallocate(sparse_rows)
	   deallocate(sparse_columns)

	   call get_timing(ud_time2)
           ud_time = ud_time2-ud_time1


           call get_timing(loop_time2)

           loop_time = loop_time2-loop_time1
           sum_time  = exp_time+psi_time+contract_time+d_time+fft_time+ud_time

           write(*,18) '|  - computing exponential : ',exp_time, ' sec.          |'
           write(*,18) '|  - getting wfc           : ',psi_time, ' sec.          |'
           write(*,18) '|  - going from 1d to 3d   : ',d_time,   ' sec.          |'
           write(*,18) '|  - performing FFTW       : ',fft_time, ' sec.          |'
           write(*,18) '|  - contracting Uk and wfc: ',contract_time,' sec.          |'
           write(*,18) '|  - building Wannier      : ',ud_time,  ' sec.          |'
           write(*,18) '|  - sum of the above      : ',sum_time, ' sec.          |'
           write(*,20) '|  ---------------------------------------------    |'
           write(*,18) '|  - total measured time   : ',loop_time,' sec.          |'
           write(*,20) '|---------------------------------------------------|'

        end do ! ikpt

        ! deallocate temperorary arrays, which might change dimensions at next
        ! iteration


        write(*,20) '| Updating the Wannier functions in scratch file... |'

        call get_timing(ud_time1)
	deallocate(list_ir_of_sir)
	deallocate(e_ikr)
	deallocate(rmesh)

        ! write Wannier function to file!
        allocate(Wf(snr))

     
        record_index = 0
        do ipol =1, nspin
          do nb_W = 1,num_wann

	     record_index = record_index + 1 
	     read(W_io_unit,rec=record_index) Wf

	     do sir_loop =1, nr_block_size
	        sir = sir_min + sir_loop-1

                Wf(sir) = Wf(sir) + Wannier_functions(sir_loop,nb_W,ipol)

	     end do

	     write(W_io_unit,rec=record_index) Wf
          end do
       end do


        ! finally, deallocate the Wannier functions
	deallocate(Wf)
	deallocate(Wannier_functions)

        call get_timing(ud_time2)
	ud_time = ud_time2-ud_time1

        write(*,18) '|  -  updating time       : ',ud_time,' sec.           |'

end do! i_block 

!--------------------------------------
! deallocate all that is no longer
! needed
!--------------------------------------

deallocate(wfc_k)
deallocate(wfc_k_3D)
deallocate(wfc_r)
deallocate(QE_eig)
deallocate(list_iG)
deallocate(U_k)


call deallocate_spin_symmetry_matrices()
call cleanup()
call cleanup_fftw_3D()


! normalize the Wannier functions 

write(*,20) '| Normalizing the Wannier functions ...             |'
call get_timing(ud_time1)
allocate(Wf(snr))

record_index = 0
do ipol =1, nspin
  do nb_W = 1,num_wann
	record_index = record_index + 1 
	read(W_io_unit,rec=record_index) Wf

	Wf(:) = Wf(:)*normalization

        write(W_io_unit,rec=record_index) Wf

  end do
end do
deallocate(Wf)
call get_timing(ud_time2)
ud_time = ud_time2-ud_time1

write(*,18) '|  -  normalizing time    : ',ud_time,' sec.           |'


allocate(rmesh(3,snr))

do sir = 1, snr

        call switch_indices(snr1,snr2,snr3,sir,sir1,sir2,sir3,switch_singlet_to_triplet)

   	rvec(1) = dble(sir1-snr1/2-1)/dble(nr1)
        rvec(2) = dble(sir2-snr2/2-1)/dble(nr2)
        rvec(3) = dble(sir3-snr3/2-1)/dble(nr3)

	rmesh(:,sir) = rvec(:)

end do


!================================================================================
!       output to netcdf file
!================================================================================
write(*,20) '| Writing to netcdf file...                         |'

call get_timing(ud_time1)
call output_wannier_functions_netcdf(W_io_unit, snr1,snr2,snr3,snr, rmesh)
call get_timing(ud_time2)
ud_time = ud_time2-ud_time1
write(*,18) '|  -  writing to file     : ',ud_time,' sec.           |'

!================================================================================
!       finish up
!================================================================================


call get_timing(total_time2)

total_time = total_time2-total_time1

write(*,20) '====================================================='
write(*,10) '|       DONE     -total time::',total_time,' seconds    |'
write(*,20) '====================================================='




10 format(A,F10.1,A)
11 format(A,I8,A)
15 format(A,I6)
16 format(A,I6,A,I6,A)
17 format(15X,A,ES12.3)
18 format(A,ES8.2,A)
19 format(A,F8.3,A)
20 format(A)
21 format(5X,I6,I12,10X,I12)
30 format(A,I4,A)
40 format(A,I4,A,I4)



end program generate_wannier_functions

