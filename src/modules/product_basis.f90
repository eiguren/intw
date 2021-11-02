!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
module intw_product_basis
!----------------------------------------------------------------------------!
!
!       This module contains subroutines which compute the product basis
!	using Maximally localized Wannier functions.
!
!----------------------------------------------------------------------------!

use intw_useful_constants

contains


  subroutine output_wannier_functions_netcdf(W_io_unit,snr1,snr2,snr3,snr,rmesh)
  !------------------------------------------------------------------------
  ! This subroutine prints out the wannier functions in nectdf format, 
  ! along with relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer        :: W_io_unit 
  integer        :: snr1, snr2, snr3, snr
  real(dp)       :: rmesh(3,snr)



  character(256) :: filename

  ! computation variables
  complex(dp)    :: Wannier_function(snr)


  integer        :: record_length , record_index

  integer        :: switch_singlet_to_triplet
  integer        :: nb_W, ipol
  integer        :: sir

  integer,allocatable :: start(:)
  integer,allocatable :: count(:)
  

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
 
  integer :: ierr  
  ! Wannier functions
  integer :: var_re_W_id, var_im_W_id
  character (*), parameter :: re_W_name     = "real_Wannier_functions"
  character (*), parameter :: im_W_name     = "imaginary_Wannier_functions"
  character (*), parameter :: W_units       = "1/a_0^{3/2}"


  ! r mesh
  integer :: var_r_mesh_id
  character (*), parameter :: rmesh_name  = "real_space_mesh"
  character (*), parameter :: rmesh_units = "crystal_basis"


  ! grid dimensions
  integer :: var_grid_id
  character (*), parameter :: grid_name  = "real_space_grid_triplet"
  character (*), parameter :: grid_units = "crystal_basis"

  ! grid dimensions
  integer :: var_supergrid_id
  character (*), parameter :: supergrid_name  = "real_space_supergrid_triplet"
  character (*), parameter :: supergrid_units = "crystal_basis"


  ! spin
  integer :: dim_sigma_id, dim_sigma
  character (*), parameter :: sigma_name    = "sigma"

  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "real_space_coordinates"


  ! crystal directions
  integer :: dim_A1A2A3_id, dim_A1A2A3
  character (*), parameter :: A1A2A3_name       = "a1_x_a2_x_a3"

  ! Wannier basis 
  integer :: dim_wann_id, dim_wann
  character (*), parameter :: wann_name       = "Wannier_index"


  ! Define the units
  character (*), parameter :: UNITS = "units"



  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! prepare the netcdf file, writing meta-data
  ! create file
  filename = trim(wannier_functions_filename)

  call check_netcdf(nf90_create(filename, NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Wannier functions on a regular real space mesh") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"The real space mesh is defined on a supercell and is centered at the  "//&
	" origin."))


  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation"))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )


  ! create the A1, A2, A3 dimensions
  dim_A1A2A3 = snr
  call check_netcdf( nf90_def_dim(nc_id,A1A2A3_name, dim_A1A2A3, dim_A1A2A3_id) )

  if ( nspin /= 1) then
     ! create the "spin" dimensions
	dim_sigma = 2
  	call check_netcdf( nf90_def_dim(nc_id,sigma_name, dim_sigma, dim_sigma_id) )
  end if

  ! create the wannier basis dimension
  dim_wann = num_wann
  call check_netcdf( nf90_def_dim(nc_id,wann_name, dim_wann, dim_wann_id) )

  !------------------------------
  ! define the variables
  !------------------------------

  ! create the grid and supergrid triplet arrays
  call check_netcdf( nf90_def_var(nc_id,grid_name, NF90_INT, &
					dim_space_id, var_grid_id) )

  call check_netcdf( nf90_def_var(nc_id,supergrid_name, NF90_INT, &
					dim_space_id, var_supergrid_id) )

  ! create the "r-mesh" array
  call check_netcdf( nf90_def_var(nc_id,rmesh_name, NF90_DOUBLE, &
			(/dim_space_id,dim_A1A2A3_id/), var_r_mesh_id) )


  ! create the real and imaginary part of the response function 

  if (nspin ==1 ) then

     call check_netcdf( nf90_def_var(nc_id,re_W_name, NF90_DOUBLE, &
     		(/dim_A1A2A3_id,dim_wann_id/), var_re_W_id) )

     call check_netcdf( nf90_def_var(nc_id,im_W_name, NF90_DOUBLE, &
     		(/dim_A1A2A3_id,dim_wann_id/), var_im_W_id) )


  else
     call check_netcdf( nf90_def_var(nc_id,re_W_name, NF90_DOUBLE, &
     (/dim_A1A2A3_id,dim_wann_id,dim_sigma_id/),  var_re_W_id) )

     call check_netcdf( nf90_def_var(nc_id,im_W_name, NF90_DOUBLE, &
     (/dim_A1A2A3_id,dim_wann_id,dim_sigma_id/),  var_im_W_id) )

  end if

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------


  ! grids
  call check_netcdf( nf90_put_att(nc_id, var_grid_id, UNITS, grid_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_supergrid_id, UNITS, supergrid_UNITS) )


  ! r-mesh
  call check_netcdf( nf90_put_att(nc_id, var_r_mesh_id, UNITS, rmesh_UNITS) )

  ! Wannier functions
  call check_netcdf( nf90_put_att(nc_id, var_im_W_id, UNITS, W_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_W_id, UNITS, W_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------

  ! write the grid parameters
  call check_netcdf( nf90_put_var(nc_id, var_grid_id, (/nr1,nr2,nr3/)))
  call check_netcdf( nf90_put_var(nc_id, var_supergrid_id, (/snr1,snr2,snr3/)))

  ! write the r-mesh
  call check_netcdf( nf90_put_var(nc_id, var_r_mesh_id, rmesh ))

  ! write the real and imaginary Wannier functions
  switch_singlet_to_triplet = -1

  if (nspin == 1) then
    allocate(start(2))
    allocate(count(2))
    ! dump the 1D array in a dummy 3D array. This is necessary
    ! to avoid overloading the RAM!

    count(:) = (/dim_A1A2A3,1/)
    
    record_index = 0
    do nb_W = 1, num_wann

       start(:) = (/1,nb_W/)

       record_index = record_index + 1 
       read(W_io_unit,rec=record_index) Wannier_function

       call check_netcdf( nf90_put_var(nc_id, var_re_W_id,  &
		real (Wannier_function(:)),start,count))

       call check_netcdf( nf90_put_var(nc_id, var_im_W_id,  &
                aimag (Wannier_function(:)),start,count))

    end do

  else


    allocate(start(3))
    allocate(count(3))

    count(:) = (/dim_A1A2A3,1,1/)

    record_index = 0
    do ipol = 1, nspin
      do nb_W = 1, num_wann

         start(:) = (/1,nb_W,ipol/)

	 record_index = record_index + 1 
	 read(W_io_unit,rec=record_index) Wannier_function


         call check_netcdf( nf90_put_var(nc_id, var_re_W_id, &
			real(Wannier_function(:)),start,count))

         call check_netcdf( nf90_put_var(nc_id, var_im_W_id, &
			aimag(Wannier_function(:)),start,count))


      end do 
    end do 

   
  end if                    

  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )
  close(unit = W_io_unit )

  end subroutine output_wannier_functions_netcdf


  subroutine split_rmesh_in_subblocks_intw(number_of_blocks, &
		nr, snr ,rest_block,nr_per_block)
  !----------------------------------------------------------------------------!
  ! It is impossible to hold all the Wannier functions on a super-cell r mesh
  ! in RAM. On the other hand, it is very slow to read-update-write the Wannier
  ! functions from a scratch file at every k-point.
  ! 
  ! It now seems much more efficient to use the RAM as much as possible by
  ! splitting the super-cell r mesh into blocks, to keep the Wannier functions
  ! on that block in RAM and to sum on all k in that block. Only write at the
  ! end of the block!
  !
  ! This subroutine determines how many such blocks there should be, and
  ! what their boundaries are.
  !
  ! Basically, we have:
  !			MAX_RAM ~ A*size_of_block+ B
  ! where A accounts for the arrays that depend on the size of the blocks
  ! and B represents the arrays which are independent of the size of the 
  ! blocks
  !----------------------------------------------------------------------------!
  use w90_parameters
  use intw_useful_constants
  use intw_reading
  implicit none

!Peio
!The variable we use instead of num_bands
    integer :: nbands_loc
!Peio

  ! input parameters

  integer  :: nr, snr

  ! output parameters
  integer :: number_of_blocks, rest_block, nr_per_block

  ! computation parameters
  real(dp) :: A_ram, B_ram

  integer  :: max_nr_per_block

  integer  :: i_loop


  nbands_loc=num_bands

  ! Rough estimate of the total memory requirement which is independent
  ! of the dense r-mesh.
  B_ram  = dble(double_complex)*dble(2*nr)                + &  ! the fr and fg arrays, for FFT
!           dble(double_complex)*dble(nG_max*nbands*nspin) + &  ! the wfc_k array
           dble(double_complex)*dble(nG_max*nbands_loc*nspin) + &  ! the wfc_k array
!           dble(double_complex)*dble(nr*nbands*nspin)     + &  ! the wfc_r array
           dble(double_complex)*dble(nr*nbands_loc*nspin)     + &  ! the wfc_r array
!           dble(double_complex)*dble(nr*nbands*nspin)     + &  ! the wfc_k_3D array
           dble(double_complex)*dble(nr*nbands_loc*nspin)     + &  ! the wfc_k_3D array
           dble(4)*dble(nG_max)                           + &   ! the list_iG array of integers
           dble(double_complex)*dble(nr*num_wann*nspin)     ! the Uk_wfc_r array
           	

  A_ram  = dble(4)                      + &   ! the list_ir_of_sir array of integers
	   dble(double_complex)         + &   ! the e_ikr array
	   dble(double_real   )*dble(3) + &   ! the r mesh
	   dble(double_complex)*dble(num_wann*nspin) ! the Wannier_functions array


  max_nr_per_block = int ( (max_ram-B_ram)/A_ram )

  ! we want the number of blocks to be such that each block has the
  ! same number of r-points (roughly). We also want 
  !  	nr_per_block <= max_nr_per_block .
  !
  ! The following identifies the optimal parameters

  i_loop = 0 
  do 
	i_loop = i_loop + 1

	number_of_blocks = i_loop

	rest_block = modulo(snr,number_of_blocks)

	nr_per_block = (snr-rest_block)/number_of_blocks

	if (nr_per_block <= max_nr_per_block )  exit
  end do


  end subroutine split_rmesh_in_subblocks_intw

  subroutine open_wannier_functions_netcdf(W_nc_id,snr1,snr2,snr3,snr)
  !------------------------------------------------------------------------
  ! This subroutine opens the netcdf file which contains the Wannier
  ! functions and reads the relevant parameters from this file.
  ! By leaving the file opened, it will be straightforward to read
  ! in the Wannier functions afterwards.
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading

  implicit none

  ! output variables
  integer        :: W_nc_id

  integer        :: snr1, snr2, snr3, snr


  ! local variables

  ! a dummy string
  character (256) :: name

  ! crystal directions
  integer :: dim_A1A2A3_id, dim_A1A2A3
  character (*), parameter :: A1A2A3_name       = "a1_x_a2_x_a3"

  ! grid dimensions
  integer :: var_supergrid_id
  character (*), parameter :: supergrid_name  = "real_space_supergrid_triplet"
  character (*), parameter :: supergrid_units = "crystal_basis"
  integer                  :: supergrid_triplets(3)


  ! Open the file in read-only mode, which is assumed to exist
  call check_netcdf(nf90_open(trim(wannier_functions_filename),  &
					NF90_NOWRITE, W_nc_id))

  !----------------------------------------------------------------
  ! Read in the dimensions of the supercell on which the 
  ! Wannier functions are defined.
  !----------------------------------------------------------------

  call check_netcdf(nf90_inq_dimid(W_nc_id, A1A2A3_name, dim_A1A2A3_id))
  call check_netcdf(nf90_inquire_dimension(W_nc_id, dim_A1A2A3_id, name, dim_A1A2A3))
  snr = dim_A1A2A3

  call check_netcdf(nf90_inq_varid(W_nc_id,supergrid_name,var_supergrid_id))
  call check_netcdf(nf90_get_var(W_nc_id, var_supergrid_id,supergrid_triplets))

  snr1 = supergrid_triplets(1)
  snr2 = supergrid_triplets(2)
  snr3 = supergrid_triplets(3)


  end subroutine open_wannier_functions_netcdf


  subroutine read_supercell_rmesh_from_wannier_functions_netcdf  &
	     (W_nc_id,snr,rmesh_SC)
  !------------------------------------------------------------------------
  ! This subroutine reads the supercell real space mesh from the 
  ! netcdf file containing the Wannier functions, which is assumed to 
  ! already be opened.
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading

  implicit none

  ! input variables
  integer        :: W_nc_id
  integer        :: snr

  ! output variables
  real(dp)       ::  rmesh_SC(3,snr)

  ! local variables

  ! r mesh
  integer :: var_r_mesh_id
  character (*), parameter :: rmesh_name  = "real_space_mesh"

  !----------------------------------------------------------------
  ! Read in the supercell real space mesh.
  !----------------------------------------------------------------

  call check_netcdf(nf90_inq_varid(W_nc_id,rmesh_name,var_r_mesh_id))
  call check_netcdf(nf90_get_var(W_nc_id, var_r_mesh_id,rmesh_SC ))


  end subroutine read_supercell_rmesh_from_wannier_functions_netcdf  

  subroutine prepare_B_functions_netcdf(qvec,rmesh_WS,nr,B_nc_id)
  !------------------------------------------------------------------------
  ! This subroutine prepares the netcdf file which will contain the 
  ! bare product basis, namely basis functions composed of the 
  ! products of Wannier functions.
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer        :: nr
  real(dp)       :: qvec(3)
  real(dp)       :: rmesh_WS(3,nr)


  ! output variables
  integer        :: B_nc_id

  integer :: ierr  

  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "crystal_basis"

  ! product basis dimension
  integer :: dim_L_id, dim_L

  ! crystal dimension
  integer :: dim_A1A2A3_id, dim_A1A2A3
  character (*), parameter :: A1A2A3_name       = "a1_x_a2_x_a3"

  ! spin
  integer :: dim_sigma
  integer :: dim_sigma1_id, dim_sigma2_id
  character (*), parameter :: sigma1_name    = "sigma1"
  character (*), parameter :: sigma2_name    = "sigma2"

  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"

  ! r mesh
  integer :: var_r_mesh_id
  character (*), parameter :: rmesh_name  = "r_mesh"
  character (*), parameter :: rmesh_units = "crystal_basis"


  character (*), parameter :: L_name       = "L"


  ! product basis functions
  integer :: var_re_B_id, var_im_B_id
  character (*), parameter :: re_B_name     = "real_B"
  character (*), parameter :: im_B_name     = "imaginary_B"
  character (*), parameter :: B_units       = "(none)"


  ! Define the units
  character (*), parameter :: UNITS = "units"

  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim(product_basis_filename), NF90_64BIT_OFFSET, B_nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(B_nc_id, NF90_GLOBAL, "title",  &
        		"Product Basis on a regular real space mesh") )

  call check_netcdf( nf90_put_att(B_nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(B_nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(B_nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))


  if (magnon) then
    call check_netcdf(nf90_put_att(B_nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation"))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(B_nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  else
    call check_netcdf(nf90_put_att(B_nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------
  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(B_nc_id,space_name, dim_space, dim_space_id) )

  ! create the A1xA2xA3 dimension, which is the flattened 3D grid
  dim_A1A2A3 = nr
  call check_netcdf( nf90_def_dim(B_nc_id,A1A2A3_name, dim_A1A2A3, dim_A1A2A3_id) )

  ! create the product basis dimension
  dim_L = num_wann**2
  call check_netcdf( nf90_def_dim(B_nc_id,L_name, dim_L, dim_L_id) )


  ! create a spin variable, if necessary
  if ( nspin /= 1 .and. .not. magnon) then
     ! create the "spin" dimensions
	dim_sigma = 2
  	call check_netcdf( nf90_def_dim(B_nc_id,sigma1_name, dim_sigma, dim_sigma1_id) )
  	call check_netcdf( nf90_def_dim(B_nc_id,sigma2_name, dim_sigma, dim_sigma2_id) )
  end if

  !------------------------------
  ! define the variables
  !------------------------------

  ! create the "q_vector" array
  call check_netcdf( nf90_def_var(B_nc_id,qvec_name, NF90_DOUBLE, &
			dim_space_id, var_qvec_id) )

  ! create the "r-mesh" array
  call check_netcdf( nf90_def_var(B_nc_id,rmesh_name, NF90_DOUBLE, &
			(/dim_space_id,dim_A1A2A3_id/), var_r_mesh_id) )


  ! create the real and imaginary part of the product basis

  if ( nspin /= 1 .and. .not. magnon) then
     call check_netcdf( nf90_def_var(B_nc_id,re_B_name, NF90_DOUBLE, &
     (/dim_A1A2A3_id,dim_L_id,dim_sigma1_id,dim_sigma2_id/),  var_re_B_id) )

     call check_netcdf( nf90_def_var(B_nc_id,im_B_name, NF90_DOUBLE, &
     (/dim_A1A2A3_id,dim_L_id,dim_sigma1_id,dim_sigma2_id/),  var_im_B_id) )

  else
     call check_netcdf( nf90_def_var(B_nc_id,re_B_name, NF90_DOUBLE, &
     (/dim_A1A2A3_id,dim_L_id/), var_re_B_id) )

     call check_netcdf( nf90_def_var(B_nc_id,im_B_name, NF90_DOUBLE, &
     (/dim_A1A2A3_id,dim_L_id/),  var_im_B_id) )

  end if

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! r-mesh
  call check_netcdf( nf90_put_att(B_nc_id, var_r_mesh_id, UNITS, rmesh_units) )

  ! Wannier functions
  call check_netcdf( nf90_put_att(B_nc_id, var_im_B_id, UNITS, B_UNITS) )
  call check_netcdf( nf90_put_att(B_nc_id, var_re_B_id, UNITS, B_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(B_nc_id) )

  !--------------------------------------------------------------------
  ! Write some data to the file
  !--------------------------------------------------------------------
  ! write the q-vector
  call check_netcdf( nf90_put_var(B_nc_id, var_qvec_id, qvec))

  ! write the r-mesh
  call check_netcdf( nf90_put_var(B_nc_id, var_r_mesh_id, rmesh_WS))


  end subroutine prepare_B_functions_netcdf

  subroutine close_B_functions_netcdf(B_nc_id)
  !-------------------------------------------------------------------------!
  ! This subroutine properly closes the product basis file so that 
  ! all buffers are flushed.
  !-------------------------------------------------------------------------!
  use netcdf

  use intw_utility

  implicit none

  integer  :: B_nc_id

  call check_netcdf( nf90_close(B_nc_id) )

  end subroutine close_B_functions_netcdf

  subroutine open_B_functions_netcdf(B_nc_id)
  !-------------------------------------------------------------------------!
  ! This subroutine properly closes the product basis file so that 
  ! all buffers are flushed.
  !-------------------------------------------------------------------------!
  use netcdf

  use intw_utility
  use intw_input_parameters

  implicit none

  integer  :: B_nc_id

  call check_netcdf(nf90_open(trim(product_basis_filename), &
					NF90_NOWRITE, B_nc_id))

  end subroutine open_B_functions_netcdf


  subroutine close_W_functions_netcdf(W_nc_id)
  !-------------------------------------------------------------------------!
  ! This subroutine properly closes the Wannier functions file .
  !-------------------------------------------------------------------------!
  use netcdf

  use intw_utility

  implicit none

  integer  :: W_nc_id

  call check_netcdf( nf90_close(W_nc_id) )

  end subroutine close_W_functions_netcdf

  subroutine read_wannier_functions_netcdf  &
	     	(W_nc_id,nb_W,snr,wannier_function)
  !------------------------------------------------------------------------
  ! This subroutine reads the supercell real space mesh from the 
  ! netcdf file containing the Wannier functions, which is assumed to 
  ! already be opened.
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading

  implicit none

  ! input variables
  integer        :: W_nc_id
  integer        :: nb_W
  integer        :: snr

  ! output variables
  complex(dp)    ::  wannier_function(snr,nspin)

  real(dp)       ::  re_WF(snr)
  integer        ::  ipol
  ! local variables

  integer,allocatable :: start(:)
  integer,allocatable :: count(:)

  ! r mesh
  integer :: ierr  
  ! Wannier functions
  integer :: var_re_W_id, var_im_W_id
  character (*), parameter :: re_W_name     = "real_Wannier_functions"
  character (*), parameter :: im_W_name     = "imaginary_Wannier_functions"
  character (*), parameter :: W_units       = "1/a_0^{3/2}"

  !----------------------------------------------------------------
  ! Read in the supercell real space mesh.
  !----------------------------------------------------------------
  call check_netcdf(nf90_inq_varid(W_nc_id,re_W_name,var_re_W_id))
  call check_netcdf(nf90_inq_varid(W_nc_id,im_W_name,var_im_W_id))

  if (nspin == 1) then
    allocate(start(2))
    allocate(count(2))

    start(:) = (/1,nb_W/)
    count(:) = (/snr,1/)

    ipol = 1
    call check_netcdf(nf90_get_var(W_nc_id, var_re_W_id,re_WF, start,count ))
    wannier_function(:,ipol) = cmplx_1*re_WF(:)

    call check_netcdf(nf90_get_var(W_nc_id, var_im_W_id,re_WF, start,count ))
    wannier_function(:,ipol) = wannier_function(:,ipol)+ cmplx_i*re_WF(:)


  else

    allocate(start(3))
    allocate(count(3))


    do ipol =1, nspin

    	start(:) = (/1,nb_W,ipol/)
    	count(:) = (/snr,1,1/)

    	call check_netcdf(nf90_get_var(W_nc_id, var_re_W_id,re_WF, start,count ))
    	wannier_function(:,ipol) = cmplx_1*re_WF(:)

    	call check_netcdf(nf90_get_var(W_nc_id, var_im_W_id,re_WF, start,count ))
    	wannier_function(:,ipol) = wannier_function(:,ipol)+ cmplx_i*re_WF(:)


    end do

  end if



  deallocate(start)
  deallocate(count)





  end subroutine read_wannier_functions_netcdf  

  subroutine write_B_functions_netcdf  &
	     (B_nc_id,L_index,B_function,nr)
  !------------------------------------------------------------------------
  ! This subroutine writes the product basis function to file
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading

  implicit none

  ! input variables
  integer        :: B_nc_id
  integer        :: L_index, nr
  complex(dp)    :: B_function(nr)



  ! local variables
  integer,allocatable :: start(:)
  integer,allocatable :: count(:)

  integer :: ierr  

  ! product basis functions
  integer :: var_re_B_id, var_im_B_id
  character (*), parameter :: re_B_name     = "real_B"
  character (*), parameter :: im_B_name     = "imaginary_B"



  if ( nspin /= 1 .and. .not. magnon) then
	write(*,*) 'NOT IMPLEMENTED! REVIEW CODE!'
	stop
  else
    allocate(start(2))
    allocate(count(2))

    start(:) = (/1,L_index/)
    count(:) = (/nr,1/)

  end if


  call check_netcdf(nf90_inq_varid(B_nc_id,re_B_name,var_re_B_id))
  call check_netcdf(nf90_inq_varid(B_nc_id,im_B_name,var_im_B_id))

  call check_netcdf(nf90_put_var(B_nc_id, var_re_B_id, &
				real(B_function), start,count ))

  call check_netcdf(nf90_put_var(B_nc_id, var_im_B_id, &
				aimag(B_function), start,count ))


  deallocate(start)
  deallocate(count)


  end subroutine write_B_functions_netcdf  

  subroutine read_B_function_netcdf  &
	     (B_nc_id,L_index,B_function,nr)
  !------------------------------------------------------------------------
  ! This subroutine writes the product basis function to file
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading

  implicit none

  ! input variables
  integer        :: B_nc_id
  integer        :: L_index, nr
  complex(dp)    :: B_function(nr)




  ! local variables

  real(dp)       :: re_B_function(nr)


  integer :: start(2)
  integer :: count(2)

  integer :: ierr  

  ! product basis functions
  integer :: var_re_B_id, var_im_B_id
  character (*), parameter :: re_B_name     = "real_B"
  character (*), parameter :: im_B_name     = "imaginary_B"



  if ( nspin /= 1 .and. .not. magnon) then
	write(*,*) 'NOT IMPLEMENTED! REVIEW CODE!'
	stop
  end if

  start(:) = (/1,L_index/)
  count(:) = (/nr,1/)

  call check_netcdf(nf90_inq_varid(B_nc_id,re_B_name,var_re_B_id))
  call check_netcdf(nf90_inq_varid(B_nc_id,im_B_name,var_im_B_id))

  call check_netcdf(nf90_get_var(B_nc_id, var_re_B_id, &
			re_B_function, start,count ))

  B_function(:) = cmplx_1*re_B_function(:)

  call check_netcdf(nf90_get_var(B_nc_id, var_im_B_id, &
				re_B_function, start,count ))

  B_function(:) = B_function(:)+cmplx_i*re_B_function(:)


  end subroutine read_B_function_netcdf  


  subroutine output_overlaps_netcdf(overlaps,kernel,qvec,L_size)
  !------------------------------------------------------------------------
  ! This subroutine prints out the wannier functions in nectdf format, 
  ! along with relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer           :: L_size
  real(dp)          :: qvec(3)
  complex(dp)       :: overlaps(L_size,L_size)
  complex(dp)       :: kernel(L_size,L_size)


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
 
  integer :: ierr  

  ! overlaps
  integer :: var_re_O_id, var_im_O_id
  character (*), parameter :: re_O_name     = "real_O"
  character (*), parameter :: im_O_name     = "imaginary_O"
  character (*), parameter :: O_units       = "none"

  ! kernel
  integer :: var_re_K_id, var_im_K_id
  character (*), parameter :: re_K_name     = "real_kernel"
  character (*), parameter :: im_K_name     = "imaginary_kernel"
  character (*), parameter :: K_units       = "eV_a0^3"

  ! product basis dimensions
  integer :: dim_L1_id, dim_L1
  integer :: dim_L2_id, dim_L2
  character (*), parameter :: L1_name       = "L1"
  character (*), parameter :: L2_name       = "L2"


  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "crystal_basis"


  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"
  character (*), parameter :: qvec_units = "crystal_basis"

  ! Define the units
  character (*), parameter :: UNITS = "units"



  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim(overlaps_filename), NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Bare basis functions overlaps and exchange-correlation kernel") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))


  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation"))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the "L1" dimension 
  dim_L1  = L_size
  call check_netcdf( nf90_def_dim(nc_id,L1_name, dim_L1, dim_L1_id) )

  ! create the "L2" dimension 
  dim_L2  = L_size
  call check_netcdf( nf90_def_dim(nc_id,L2_name, dim_L2, dim_L2_id) )

  !------------------------------
  ! define the variables
  !------------------------------

  ! create the "q_vector" array
  call check_netcdf( nf90_def_var(nc_id,qvec_name, NF90_DOUBLE, &
			dim_space_id, var_qvec_id) )

  ! create the "overlaps" array

  call check_netcdf( nf90_def_var(nc_id,re_O_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_re_O_id) )

  call check_netcdf( nf90_def_var(nc_id,im_O_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_im_O_id) )

  ! create the "kernel" array

  call check_netcdf( nf90_def_var(nc_id,re_K_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_re_K_id) )

  call check_netcdf( nf90_def_var(nc_id,im_K_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_im_K_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! overlaps
  call check_netcdf( nf90_put_att(nc_id, var_im_O_id, UNITS, O_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_O_id, UNITS, O_UNITS) )

  ! kernel
  call check_netcdf( nf90_put_att(nc_id, var_im_K_id, UNITS, K_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_K_id, UNITS, K_UNITS) )

  ! q vector
  call check_netcdf( nf90_put_att(nc_id, var_qvec_id, UNITS, qvec_units) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------


  ! write the q vector
  call check_netcdf( nf90_put_var(nc_id, var_qvec_id, qvec))

  ! write the real and imaginary overlaps
  call check_netcdf( nf90_put_var(nc_id, var_re_O_id,real (overlaps)))
  call check_netcdf( nf90_put_var(nc_id, var_im_O_id,aimag(overlaps)))

  ! write the real and imaginary kernel
  call check_netcdf( nf90_put_var(nc_id, var_re_K_id,real (kernel)))
  call check_netcdf( nf90_put_var(nc_id, var_im_K_id,aimag(kernel)))


  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_overlaps_netcdf




  subroutine read_field_ascii_file(filename,nr,field)
  !----------------------------------------------------------------------------!
  ! This subroutine reads a complex field written to an ascii file, 
  ! assuming it contains 3 columns, namely ir, real[field], imag[field].
  ! It assumes the comments begin with "#".
  !----------------------------------------------------------------------------!
  use w90_parameters
  use intw_useful_constants
  use intw_reading
  implicit none

  ! input parameters
  character(*)  :: filename
  integer       :: nr

  ! output parameters
  complex(dp)   :: field(nr)

  ! computation parameters
  integer        :: io_unit
  integer        :: iline, skip_line


  character(256) :: header_string

  integer        :: ir
  real(dp)       :: x, y

  io_unit = find_free_unit()

  open(io_unit,file=filename)

  ! read the header

  skip_line = 0
  do 

     read(io_unit,*) header_string
     if (header_string(1:1) /= '#') exit
     skip_line = skip_line + 1
  end do
  rewind(io_unit)

  do iline = 1, skip_line
     read(io_unit,*) header_string
  end do

  do iline = 1, nr
     read(io_unit,35) ir,x,y
     field(ir) = cmplx_1*x+cmplx_i*y
  end do

  close(io_unit)

35 format(I6,2E20.12)

  end subroutine read_field_ascii_file

  subroutine change_order(nr,field_wrong_order,field_right_order)
  !----------------------------------------------------------------------------!
  ! This subroutine goes from C ordering to F ordering.
  !----------------------------------------------------------------------------!
  use intw_useful_constants
  use intw_utility
  use intw_reading
  implicit none

  ! input parameters
  integer        :: nr
  complex(dp)    :: field_wrong_order(nr)

  ! output parameters
  complex(dp)    :: field_right_order(nr)

  ! local parameters
  integer        :: ir_wrong, ir_right
  integer        :: switch_singlet_to_triplet
  integer        :: switch_triplet_to_singlet

  integer        :: ir1, ir2, ir3

  switch_singlet_to_triplet = -1
  switch_triplet_to_singlet =  1

  do ir_wrong = 1 , nr

  	call switch_indices(nr3,nr2,nr1,ir_wrong,ir3,ir2,ir1,switch_singlet_to_triplet)

  	call switch_indices(nr1,nr2,nr3,ir_right,ir1,ir2,ir3,switch_triplet_to_singlet)

	field_right_order(ir_right) = field_wrong_order(ir_wrong)
  end do

  end subroutine change_order


  subroutine read_overlaps_netcdf(L_size,overlaps,kernel,qvec)
  !------------------------------------------------------------------------
  ! This subroutine prints out the wannier functions in nectdf format, 
  ! along with relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer           :: L_size

  !output variables
  complex(dp)       :: overlaps(L_size,L_size)
  complex(dp)       :: kernel(L_size,L_size)
  real(dp)          :: qvec(3)


  ! local variables
  real(dp)       :: real_field(L_size,L_size)

  integer :: nc_id        
  ! netcdf file id
 
  integer :: ierr  

  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"


  ! overlaps
  integer :: var_re_O_id, var_im_O_id
  character (*), parameter :: re_O_name     = "real_O"
  character (*), parameter :: im_O_name     = "imaginary_O"

  ! kernel
  integer :: var_re_K_id, var_im_K_id
  character (*), parameter :: re_K_name     = "real_kernel"
  character (*), parameter :: im_K_name     = "imaginary_kernel"


  !------------------------------
  ! read overlaps
  !------------------------------
  call check_netcdf(nf90_open(trim(overlaps_filename), NF90_NOWRITE, nc_id))


  call check_netcdf(nf90_inq_varid(nc_id,qvec_name,var_qvec_id))
  call check_netcdf(nf90_get_var(nc_id, var_qvec_id,qvec))

  call check_netcdf(nf90_inq_varid(nc_id,re_O_name,var_re_O_id))
  call check_netcdf(nf90_inq_varid(nc_id,im_O_name,var_im_O_id))


  call check_netcdf(nf90_inq_varid(nc_id,re_K_name,var_re_K_id))
  call check_netcdf(nf90_inq_varid(nc_id,im_K_name,var_im_K_id))


  call check_netcdf(nf90_get_var(nc_id, var_re_O_id,real_field))
  overlaps(:,:) = cmplx_1*real_field(:,:)

  call check_netcdf(nf90_get_var(nc_id, var_im_O_id,real_field))
  overlaps(:,:) = overlaps(:,:) + cmplx_i*real_field(:,:)


  call check_netcdf(nf90_get_var(nc_id, var_re_K_id,real_field))
  kernel(:,:) = cmplx_1*real_field(:,:)

  call check_netcdf(nf90_get_var(nc_id, var_im_K_id,real_field))
  kernel(:,:) = kernel(:,:) + cmplx_i*real_field(:,:)

  call check_netcdf( nf90_close(nc_id) )

  end subroutine read_overlaps_netcdf

  subroutine diagonalize_matrix(dim,matrix,V_matrix,Lambda)
  !------------------------------------------------------------------------
  ! This subroutine is a wrapper around the LAPACK routine zheevd,
  ! which diagonalizes an hermitian matrix, returning its eigenvalues and 
  ! eigenvectors.
  !	The eigenvectors of matrix are stored in the columns of V_matrix,
  !	and the eivenvalues in Lambda.
  !------------------------------------------------------------------------


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer           :: dim
  complex(dp)       :: matrix(dim,dim)

  !output variables
  complex(dp)       :: V_matrix(dim,dim)
  real(dp)          :: Lambda(dim)

  ! local variables
  integer   :: lwork,  liwork, lrwork
  integer   :: info, n

  complex(dp),allocatable :: work(:)
  real(dp),   allocatable :: rwork(:)
  integer ,   allocatable :: iwork(:)



  V_matrix(:,:) = matrix(:,:)

!call zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)

  n = dim
  lwork  =   n**2+ 2*n
  lrwork = 2*n**2+ 5*n+ 1
  liwork =         5*n+ 3


  allocate(work (lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))

  call zheevd('V',             & ! Eigenvalues and Eigenvectors are computed
              'U',             & ! doesn't matter
                n,             & ! size of the matrix 
         V_matrix,             & ! Matrix to be diagonalized on input, eigenvectors on outpu
                n,             & ! size of the matrix 
           Lambda,             & ! sorted eigenvalues
             work, lwork, rwork, lrwork, iwork, liwork, info) ! junk

  if ( info /= 0 ) then
        write(*,*) 'info = ',info
        write(*,*) '*******************************************'
        write(*,*) '* ERROR: diagonalization failed!          *'
        write(*,*) '*******************************************'
        stop
  end if

  deallocate(work )
  deallocate(rwork)
  deallocate(iwork)



  end subroutine diagonalize_matrix

  subroutine output_functional_projectors_netcdf(qvec,bare_basis_size, &
		special_functions_basis_size, optimal_basis_size_local,&
		T_functional_projectors,list_m1,list_m2, &
		special_functions_t_coefficients,Pt,kernel)
  !------------------------------------------------------------------------
  ! This subroutine outputs to netcdf all the relevant information to
  ! perform chi_KS calculation in the optimal basis.
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_input_parameters
  use intw_reading
  use intw_utility

  implicit none

  ! input variables
  real(dp)       :: qvec(3)
  integer        :: bare_basis_size, optimal_basis_size_local
                    ! define optimal_basis_size_local to not conflict with the
                    ! global variable!

  integer        :: special_functions_basis_size

  integer        :: list_m1(bare_basis_size), list_m2(bare_basis_size) 

  complex(dp)    :: T_functional_projectors(optimal_basis_size_local,bare_basis_size)
  complex(dp)    :: special_functions_t_coefficients(optimal_basis_size_local,  &
						special_functions_basis_size)
  complex(dp)    :: kernel(optimal_basis_size_local,optimal_basis_size_local)

  ! | b_i > = sum_L | B_L > Pt_{Li}
  complex(dp)    :: Pt(bare_basis_size,optimal_basis_size_local)

  ! local variables
  integer        :: nc_id

  integer :: ierr  

  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "crystal_basis"

  ! product basis dimension
  integer :: dim_L_id, dim_L
  character (*), parameter :: L_name       = "L"

  ! optimal basis dimension
  integer :: dim_I1_id, dim_I1
  character (*), parameter :: I1_name       = "I1_projected"

  integer :: dim_I2_id, dim_I2
  character (*), parameter :: I2_name       = "I2_projected"

  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"
  character (*), parameter :: q_UNITS = "crystal_basis"


  ! product basis functions
  integer :: var_list_m1_id, var_list_m2_id
  character (*), parameter :: list_m1_name  = "list_m1"
  character (*), parameter :: list_m2_name  = "list_m2"

  ! product basis functions
  integer :: var_re_projector_id, var_im_projector_id
  character (*), parameter :: re_projector_name  = "real_projector"
  character (*), parameter :: im_projector_name  = "imaginary_projector"
  character (*), parameter :: projector_UNITS = "none"

  ! kernel 
  integer :: var_re_kernel_id, var_im_kernel_id
  character (*), parameter :: re_kernel_name  = "real_projected_kernel"
  character (*), parameter :: im_kernel_name  = "imaginary_projected_kernel"
  character (*), parameter :: kernel_UNITS    = "eV*a_0^3"

  !  < 1 | b_i > 
  integer :: var_re_b_integral_id, var_im_b_integral_id
  character (*), parameter :: re_b_integral_name  = "real_b_integral"
  character (*), parameter :: im_b_integral_name  = "imaginary_b_integral"
  character (*), parameter :: integral_UNITS    = "none"


  !  < b_i | 1 > 
  integer :: var_re_1_coeff_id, var_im_1_coeff_id
  character (*), parameter :: re_1_coeff_name  = "real_1_b_coefficients"
  character (*), parameter :: im_1_coeff_name  = "imaginary_1_b_coefficients"
  character (*), parameter :: one_UNITS    = "none"

  !  < b_i | mz > 
  integer :: var_re_mz_coeff_id, var_im_mz_coeff_id
  character (*), parameter :: re_mz_coeff_name  = "real_mz_b_coefficients"
  character (*), parameter :: im_mz_coeff_name  = "imaginary_mz_b_coefficients"
  character (*), parameter :: mz_UNITS    = "1/a_0^3"

  !  < b_i | Delta > 
  integer :: var_re_Delta_coeff_id, var_im_Delta_coeff_id
  character (*), parameter :: re_Delta_coeff_name  = "real_Delta_b_coefficients"
  character (*), parameter :: im_Delta_coeff_name  = "imaginary_Delta_b_coefficients"
  character (*), parameter :: Delta_UNITS    = "eV"

  !  | b_i > = sum_L |B_L > b_{Li}
  integer :: var_re_basis_B_coeff_id, var_im_basis_B_coeff_id
  character (*), parameter :: re_basis_B_coeff_name  = "real_basis_B_coefficients"
  character (*), parameter :: im_basis_B_coeff_name  = "imaginary_basis_B_coefficients"

  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! Define the units
  character (*), parameter :: UNITS = "units"

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim(projectors_filename), NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Functional projectors and optimal basis coefficients") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))


  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation"))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------
  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the bare product basis dimension
  dim_L = bare_basis_size
  call check_netcdf( nf90_def_dim(nc_id,L_name, dim_L, dim_L_id) )

  ! create the optimal basis dimensions
  dim_I1 = optimal_basis_size_local
  call check_netcdf( nf90_def_dim(nc_id,I1_name, dim_I1, dim_I1_id) )

  dim_I2 = optimal_basis_size_local
  call check_netcdf( nf90_def_dim(nc_id,I2_name, dim_I2, dim_I2_id) )

  !------------------------------
  ! define the variables
  !------------------------------

  ! create the "q_vector" array
  call check_netcdf( nf90_def_var(nc_id,qvec_name, NF90_DOUBLE, &
			dim_space_id, var_qvec_id) )


  ! create the list_m1, list_m2 array
  call check_netcdf( nf90_def_var(nc_id,list_m1_name, NF90_INT, &
			dim_L_id, var_list_m1_id) )

  call check_netcdf( nf90_def_var(nc_id,list_m2_name, NF90_INT, &
			dim_L_id, var_list_m2_id) )



  ! create the projectors array
  call check_netcdf( nf90_def_var(nc_id,re_projector_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_re_projector_id) )

  call check_netcdf( nf90_def_var(nc_id,im_projector_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_im_projector_id) )


  ! create the kernel array
  call check_netcdf( nf90_def_var(nc_id,re_kernel_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_re_kernel_id) )

  call check_netcdf( nf90_def_var(nc_id,im_kernel_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_im_kernel_id) )

  ! create the integrals array
  call check_netcdf( nf90_def_var(nc_id,re_b_integral_name, NF90_DOUBLE, &
			dim_I1_id, var_re_b_integral_id) )

  call check_netcdf( nf90_def_var(nc_id,im_b_integral_name, NF90_DOUBLE, &
			dim_I1_id, var_im_b_integral_id) )


  ! create the various coefficients
  call check_netcdf( nf90_def_var(nc_id,re_1_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_re_1_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,im_1_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_im_1_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,re_mz_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_re_mz_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,im_mz_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_im_mz_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,re_Delta_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_re_Delta_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,im_Delta_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_im_Delta_coeff_id) )


  ! The basis coefficients in the original B basis
  call check_netcdf( nf90_def_var(nc_id,re_basis_B_coeff_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_re_basis_B_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,im_basis_B_coeff_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_im_basis_B_coeff_id) )


  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, var_qvec_id, UNITS, q_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_projector_id, UNITS, projector_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_projector_id, UNITS, projector_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_kernel_id, UNITS, kernel_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_kernel_id, UNITS, kernel_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_b_integral_id, UNITS, integral_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_b_integral_id, UNITS, integral_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_1_coeff_id, UNITS, one_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_1_coeff_id, UNITS, one_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_mz_coeff_id, UNITS, mz_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_mz_coeff_id, UNITS, mz_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_Delta_coeff_id, UNITS, Delta_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_Delta_coeff_id, UNITS, Delta_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write some data to the file
  !--------------------------------------------------------------------
  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_qvec_id, qvec))


  call check_netcdf( nf90_put_var(nc_id,var_list_m1_id, list_m1))
  call check_netcdf( nf90_put_var(nc_id,var_list_m2_id, list_m2))

  call check_netcdf( nf90_put_var(nc_id,var_re_projector_id,  &
				real(transpose(T_functional_projectors))))

  call check_netcdf( nf90_put_var(nc_id,var_im_projector_id,  &
				aimag(transpose(T_functional_projectors))))



  call check_netcdf( nf90_put_var(nc_id,var_re_kernel_id, real(kernel)))
  call check_netcdf( nf90_put_var(nc_id,var_im_kernel_id, aimag(kernel)))

  call check_netcdf( nf90_put_var(nc_id,var_re_b_integral_id, &
			real(conjg(special_functions_t_coefficients(:,1)))))
  call check_netcdf( nf90_put_var(nc_id,var_im_b_integral_id, &
			aimag(conjg(special_functions_t_coefficients(:,1)))))

  call check_netcdf( nf90_put_var(nc_id,var_re_1_coeff_id, &
					real(special_functions_t_coefficients(:,1))))
  call check_netcdf( nf90_put_var(nc_id,var_im_1_coeff_id, &
					aimag(special_functions_t_coefficients(:,1))))

  call check_netcdf( nf90_put_var(nc_id,var_re_mz_coeff_id, &
					real(special_functions_t_coefficients(:,2))))
  call check_netcdf( nf90_put_var(nc_id,var_im_mz_coeff_id, &
					aimag(special_functions_t_coefficients(:,2))))

  call check_netcdf( nf90_put_var(nc_id,var_re_Delta_coeff_id, &
					real(special_functions_t_coefficients(:,3))))
  call check_netcdf( nf90_put_var(nc_id,var_im_Delta_coeff_id, &
					aimag(special_functions_t_coefficients(:,3))))

  call check_netcdf( nf90_put_var(nc_id,var_re_basis_B_coeff_id,real(Pt(:,:))))
  call check_netcdf( nf90_put_var(nc_id,var_im_basis_B_coeff_id,aimag(Pt(:,:))))

  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_functional_projectors_netcdf

  subroutine output_optimal_basis_netcdf(qvec, nr,  &
		optimal_basis_size_local, b_basis_functions)
  !------------------------------------------------------------------------
  ! This subroutine outputs to netcdf all the relevant information to
  ! perform chi_KS calculation in the optimal basis.
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_input_parameters
  use intw_reading
  use intw_utility

  implicit none

  ! input variables
  real(dp)       :: qvec(3)
  integer        :: nr
  integer        :: optimal_basis_size_local
  complex(dp)    :: b_basis_functions(nr,optimal_basis_size_local)


  ! local variables
  integer        :: nc_id

  integer :: ierr  

  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "crystal_basis"

  ! optimal basis dimension
  integer :: dim_I1_id, dim_I1
  character (*), parameter :: I1_name       = "I1_projected"

  ! crystal directions
  integer :: dim_A1A2A3_id, dim_A1A2A3
  character (*), parameter :: A1A2A3_name       = "a1_x_a2_x_a3"

  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"
  character (*), parameter :: q_UNITS = "crystal_basis"

  ! optimal basis functions
  integer :: var_re_optimal_basis_id, var_im_optimal_basis_id
  character (*), parameter :: re_optimal_basis_name  = "real_optimal_basis"
  character (*), parameter :: im_optimal_basis_name  = "imaginary_optimal_basis"


  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! Define the units
  character (*), parameter :: UNITS = "units"

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim('optimal_basis_functions.nc'), &
						NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
	"Optimal basis functions, explicitely as function of r.") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))


  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation"))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------
  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the optimal basis dimensions
  dim_I1 = optimal_basis_size
  call check_netcdf( nf90_def_dim(nc_id,I1_name, dim_I1, dim_I1_id) )

  ! create the optimal basis dimensions
  dim_A1A2A3 = nr
  call check_netcdf( nf90_def_dim(nc_id,A1A2A3_name, dim_A1A2A3, dim_A1A2A3_id) )


  !------------------------------
  ! define the variables
  !------------------------------

  ! create the "q_vector" array
  call check_netcdf( nf90_def_var(nc_id,qvec_name, NF90_DOUBLE, &
					dim_space_id, var_qvec_id) )



  ! create the projectors array
  call check_netcdf( nf90_def_var(nc_id,re_optimal_basis_name, NF90_DOUBLE, &
			(/dim_A1A2A3_id,dim_I1_id/), var_re_optimal_basis_id) )

  call check_netcdf( nf90_def_var(nc_id,im_optimal_basis_name, NF90_DOUBLE, &
			(/dim_A1A2A3_id,dim_I1_id/), var_im_optimal_basis_id) )


  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, var_qvec_id, UNITS, q_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write some data to the file
  !--------------------------------------------------------------------
  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_qvec_id, qvec))


  call check_netcdf( nf90_put_var(nc_id,var_re_optimal_basis_id,  &
						real(b_basis_functions)))

  call check_netcdf( nf90_put_var(nc_id,var_im_optimal_basis_id,  &
						aimag(b_basis_functions)))


  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_optimal_basis_netcdf

  subroutine output_full_projected_kernel_netcdf(qvec,ortho_basis_size,     &
								b_basis_kernel)
  !------------------------------------------------------------------------
  ! This subroutine outputs to netcdf the kernel in the full orthonormalized
  ! and "optimized" basis.
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_input_parameters
  use intw_reading
  use intw_utility

  implicit none

  ! input variables
  real(dp)       :: qvec(3)
  integer        :: ortho_basis_size
  complex(dp)    :: b_basis_kernel(ortho_basis_size,ortho_basis_size)


  ! local variables
  integer        :: nc_id

  integer :: ierr  

  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "crystal_basis"

  ! optimal basis dimension
  integer :: dim_I1_id, dim_I1
  character (*), parameter :: I1_name       = "I1"

  integer :: dim_I2_id, dim_I2
  character (*), parameter :: I2_name       = "I2"


  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"
  character (*), parameter :: q_UNITS = "crystal_basis"

  ! kernel 
  integer :: var_re_kernel_id, var_im_kernel_id
  character (*), parameter :: re_kernel_name  = "real_full_projected_kernel"
  character (*), parameter :: im_kernel_name  = "imaginary_full_projected_kernel"
  character (*), parameter :: kernel_UNITS    = "eV*a_0^3"

  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! Define the units
  character (*), parameter :: UNITS = "units"

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim('full_projected_kernel.nc'),   &
						NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------


  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Coefficients of the kernel in the final orthogonal basis") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  if (magnon) then
     call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"The first three basis functions span the space of {1, mz, Delta}.  The"//&
	" other functions are picked to diagonalize the kernel and sort the wings."))
  else if( nspin == 1) then
     call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"The basis functions are assumed to be plane waves, and the Kernel "//&
	"is given by the matrix elements of the Coulomb potential between plane waves."))
  end if


  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation"))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------
  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the optimal basis dimensions
  dim_I1 = ortho_basis_size
  call check_netcdf( nf90_def_dim(nc_id,I1_name, dim_I1, dim_I1_id) )

  dim_I2 = ortho_basis_size
  call check_netcdf( nf90_def_dim(nc_id,I2_name, dim_I2, dim_I2_id) )

  !------------------------------
  ! define the variables
  !------------------------------

  ! create the "q_vector" array
  call check_netcdf( nf90_def_var(nc_id,qvec_name, NF90_DOUBLE, &
			dim_space_id, var_qvec_id) )


  ! create the kernel array
  call check_netcdf( nf90_def_var(nc_id,re_kernel_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_re_kernel_id) )

  call check_netcdf( nf90_def_var(nc_id,im_kernel_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_im_kernel_id) )


  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, var_qvec_id, UNITS, q_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_kernel_id, UNITS, kernel_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_kernel_id, UNITS, kernel_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write some data to the file
  !--------------------------------------------------------------------
  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_qvec_id, qvec))


  call check_netcdf( nf90_put_var(nc_id,var_re_kernel_id, real(b_basis_kernel)))
  call check_netcdf( nf90_put_var(nc_id,var_im_kernel_id, aimag(b_basis_kernel)))


  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_full_projected_kernel_netcdf

  subroutine output_spin_resolved_Wannier_overlaps_netcdf &
				(num_wann2,wannier_spin_overlaps)
  !------------------------------------------------------------------------
  ! This subroutine prints out the spin resolved overlap of wannier functions 
  ! in nectdf format, along with relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer           :: num_wann2
  complex(dp)       :: wannier_spin_overlaps(num_wann2,num_wann2,nspin)


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
 
  integer :: ierr  


  character (*), parameter :: UNITS = "units"

  ! overlaps
  integer :: var_re_O_id, var_im_O_id
  character (*), parameter :: re_O_name     = "real_Wannier_overlap"
  character (*), parameter :: im_O_name     = "imaginary_Wannier_overlap"
  character (*), parameter :: O_units       = "none"

  ! Wannier band indices dimensions
  integer :: dim_nb1_id, dim_nb1
  integer :: dim_nb2_id, dim_nb2

  character (*), parameter :: nb1_name       = "m1"
  character (*), parameter :: nb2_name       = "m2"

  ! spin
  integer :: dim_sigma_id
  integer :: dim_sigma
  character (*), parameter :: sigma_name     = "sigma"


  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim('Wannier_overlaps.nc'), NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Wannier functions overlaps") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))


  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation"))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------

  dim_nb1  = num_wann2
  call check_netcdf( nf90_def_dim(nc_id,nb1_name, dim_nb1, dim_nb1_id) )

  dim_nb2  = num_wann2
  call check_netcdf( nf90_def_dim(nc_id,nb2_name, dim_nb2, dim_nb2_id) )


  dim_sigma  = nspin
  call check_netcdf( nf90_def_dim(nc_id,sigma_name, dim_sigma, dim_sigma_id) )

  ! create the "overlaps" array

  call check_netcdf( nf90_def_var(nc_id,re_O_name, NF90_DOUBLE, &
     		(/dim_nb1_id,dim_nb2_id,dim_sigma_id/), var_re_O_id) )

  call check_netcdf( nf90_def_var(nc_id,im_O_name, NF90_DOUBLE, &
     		(/dim_nb1_id,dim_nb2_id,dim_sigma_id/), var_im_O_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! overlaps
  call check_netcdf( nf90_put_att(nc_id, var_im_O_id, UNITS, O_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_O_id, UNITS, O_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------

  ! write the real and imaginary overlaps
  call check_netcdf( nf90_put_var(nc_id, var_re_O_id,real (wannier_spin_overlaps)))
  call check_netcdf( nf90_put_var(nc_id, var_im_O_id,aimag(wannier_spin_overlaps)))


  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_spin_resolved_Wannier_overlaps_netcdf 



  subroutine read_spin_resolved_Wannier_overlaps_netcdf &
				(num_wann2,wannier_spin_overlaps)
  !------------------------------------------------------------------------
  ! This subroutine prints out the wannier functions in nectdf format, 
  ! along with relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading

  implicit none


  ! input variables
  integer           :: num_wann2

  ! output variables
  complex(dp)       :: wannier_spin_overlaps(num_wann2,num_wann2,nspin)


  ! local variables
  real(dp)       :: real_field(num_wann2,num_wann2,nspin)


  character(256) :: wannier_overlap_filename 

  integer :: nc_id        
  ! netcdf file id
 
  integer :: ierr  


  ! overlaps
  integer :: var_re_O_id, var_im_O_id
  character (*), parameter :: re_O_name     = "real_Wannier_overlap"
  character (*), parameter :: im_O_name     = "imaginary_Wannier_overlap"

  !------------------------------
  ! read overlaps
  !------------------------------

  wannier_overlap_filename = 'Wannier_overlaps.nc'
  call check_netcdf(nf90_open(trim(wannier_overlap_filename), NF90_NOWRITE, nc_id))

  call check_netcdf(nf90_inq_varid(nc_id,re_O_name,var_re_O_id))
  call check_netcdf(nf90_inq_varid(nc_id,im_O_name,var_im_O_id))


  call check_netcdf(nf90_get_var(nc_id, var_re_O_id,real_field))
  wannier_spin_overlaps(:,:,:) = cmplx_1*real_field(:,:,:)

  call check_netcdf(nf90_get_var(nc_id, var_im_O_id,real_field))

  wannier_spin_overlaps(:,:,:) = wannier_spin_overlaps(:,:,:) +  &
				          cmplx_i*real_field(:,:,:)

  call check_netcdf( nf90_close(nc_id) )

  end subroutine read_spin_resolved_Wannier_overlaps_netcdf 

  subroutine output_functional_projectors_plasmons_netcdf(qvec,bare_basis_size, &
		special_functions_basis_size, optimal_basis_size_local,&
		T_functional_projectors,list_m1,list_m2, &
		special_functions_t_coefficients,Pt,kernel,kernel_xc,kernel_Ha)
  !------------------------------------------------------------------------
  ! This subroutine outputs to netcdf all the relevant information to
  ! perform chi_KS calculation in the optimal basis, using the plasmon
  ! formalism.
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_input_parameters
  use intw_reading
  use intw_utility

  implicit none

  ! input variables
  real(dp)       :: qvec(3)
  integer        :: bare_basis_size, optimal_basis_size_local
                    ! define optimal_basis_size_local to not conflict with the
                    ! global variable!

  integer        :: special_functions_basis_size

  integer        :: list_m1(bare_basis_size), list_m2(bare_basis_size) 

  complex(dp)    :: T_functional_projectors(optimal_basis_size_local,bare_basis_size)
  complex(dp)    :: special_functions_t_coefficients(optimal_basis_size_local,  &
						special_functions_basis_size)
  complex(dp)    :: kernel(optimal_basis_size_local,optimal_basis_size_local)
  complex(dp)    :: kernel_xc(optimal_basis_size_local,optimal_basis_size_local)
  complex(dp)    :: kernel_Ha(optimal_basis_size_local,optimal_basis_size_local)

  ! | b_i > = sum_L | B_L > Pt_{Li}
  complex(dp)    :: Pt(bare_basis_size,optimal_basis_size_local)

  ! local variables
  integer        :: nc_id

  integer :: ierr  

  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "crystal_basis"

  ! product basis dimension
  integer :: dim_L_id, dim_L
  character (*), parameter :: L_name       = "L"

  ! optimal basis dimension
  integer :: dim_I1_id, dim_I1
  character (*), parameter :: I1_name       = "I1_projected"

  integer :: dim_I2_id, dim_I2
  character (*), parameter :: I2_name       = "I2_projected"

  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"
  character (*), parameter :: q_UNITS = "crystal_basis"


  ! product basis functions
  integer :: var_list_m1_id, var_list_m2_id
  character (*), parameter :: list_m1_name  = "list_m1"
  character (*), parameter :: list_m2_name  = "list_m2"

  ! product basis functions
  integer :: var_re_projector_id, var_im_projector_id
  character (*), parameter :: re_projector_name  = "real_projector"
  character (*), parameter :: im_projector_name  = "imaginary_projector"
  character (*), parameter :: projector_UNITS = "none"

  ! kernels
  integer :: var_re_kernel_id, var_im_kernel_id
  character (*), parameter :: re_kernel_name  = "real_projected_kernel"
  character (*), parameter :: im_kernel_name  = "imaginary_projected_kernel"
  character (*), parameter :: kernel_UNITS    = "eV*a_0^3"

  integer :: var_re_kernel_xc_id, var_im_kernel_xc_id
  character (*), parameter :: re_kernel_xc_name  = "real_projected_kernel_xc"
  character (*), parameter :: im_kernel_xc_name  = "imaginary_projected_kernel_xc"

  integer :: var_re_kernel_Ha_id, var_im_kernel_Ha_id
  character (*), parameter :: re_kernel_Ha_name  = "real_projected_kernel_Hartree"
  character (*), parameter :: im_kernel_Ha_name  = "imaginary_projected_kernel_Hartree"




  !  < 1 | b_i > 
  integer :: var_re_b_integral_id, var_im_b_integral_id
  character (*), parameter :: re_b_integral_name  = "real_b_integral"
  character (*), parameter :: im_b_integral_name  = "imaginary_b_integral"
  character (*), parameter :: integral_UNITS    = "none"


  !  < b_i | 1 > 
  integer :: var_re_1_coeff_id, var_im_1_coeff_id
  character (*), parameter :: re_1_coeff_name  = "real_1_b_coefficients"
  character (*), parameter :: im_1_coeff_name  = "imaginary_1_b_coefficients"
  character (*), parameter :: one_UNITS    = "none"

  !  | b_i > = sum_L |B_L > b_{Li}
  integer :: var_re_basis_B_coeff_id, var_im_basis_B_coeff_id
  character (*), parameter :: re_basis_B_coeff_name  = "real_basis_B_coefficients"
  character (*), parameter :: im_basis_B_coeff_name  = "imaginary_basis_B_coefficients"

  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! Define the units
  character (*), parameter :: UNITS = "units"

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim(projectors_filename), NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Functional projectors and optimal basis coefficients") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))


  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
					"This is a density calculation" ))

  !------------------------------
  ! define the dimensions 
  !------------------------------
  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the bare product basis dimension
  dim_L = bare_basis_size
  call check_netcdf( nf90_def_dim(nc_id,L_name, dim_L, dim_L_id) )

  ! create the optimal basis dimensions
  dim_I1 = optimal_basis_size_local
  call check_netcdf( nf90_def_dim(nc_id,I1_name, dim_I1, dim_I1_id) )

  dim_I2 = optimal_basis_size_local
  call check_netcdf( nf90_def_dim(nc_id,I2_name, dim_I2, dim_I2_id) )

  !------------------------------
  ! define the variables
  !------------------------------

  ! create the "q_vector" array
  call check_netcdf( nf90_def_var(nc_id,qvec_name, NF90_DOUBLE, &
			dim_space_id, var_qvec_id) )


  ! create the list_m1, list_m2 array
  call check_netcdf( nf90_def_var(nc_id,list_m1_name, NF90_INT, &
			dim_L_id, var_list_m1_id) )

  call check_netcdf( nf90_def_var(nc_id,list_m2_name, NF90_INT, &
			dim_L_id, var_list_m2_id) )


  ! create the projectors array
  call check_netcdf( nf90_def_var(nc_id,re_projector_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_re_projector_id) )

  call check_netcdf( nf90_def_var(nc_id,im_projector_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_im_projector_id) )


  ! create the kernels array
  call check_netcdf( nf90_def_var(nc_id,re_kernel_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_re_kernel_id) )

  call check_netcdf( nf90_def_var(nc_id,im_kernel_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_im_kernel_id) )

  call check_netcdf( nf90_def_var(nc_id,re_kernel_xc_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_re_kernel_xc_id) )

  call check_netcdf( nf90_def_var(nc_id,im_kernel_xc_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_im_kernel_xc_id) )

  call check_netcdf( nf90_def_var(nc_id,re_kernel_Ha_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_re_kernel_Ha_id) )

  call check_netcdf( nf90_def_var(nc_id,im_kernel_Ha_name, NF90_DOUBLE, &
			(/dim_I1_id,dim_I2_id/), var_im_kernel_Ha_id) )


  ! create the integrals array
  call check_netcdf( nf90_def_var(nc_id,re_b_integral_name, NF90_DOUBLE, &
			dim_I1_id, var_re_b_integral_id) )

  call check_netcdf( nf90_def_var(nc_id,im_b_integral_name, NF90_DOUBLE, &
			dim_I1_id, var_im_b_integral_id) )


  ! create the various coefficients
  call check_netcdf( nf90_def_var(nc_id,re_1_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_re_1_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,im_1_coeff_name, NF90_DOUBLE, &
			dim_I1_id, var_im_1_coeff_id) )

  ! The basis coefficients in the original B basis
  call check_netcdf( nf90_def_var(nc_id,re_basis_B_coeff_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_re_basis_B_coeff_id) )

  call check_netcdf( nf90_def_var(nc_id,im_basis_B_coeff_name, NF90_DOUBLE, &
			(/dim_L_id,dim_I1_id/), var_im_basis_B_coeff_id) )


  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, var_qvec_id, UNITS, q_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_projector_id, UNITS, projector_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_projector_id, UNITS, projector_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_kernel_id, UNITS, kernel_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_kernel_id, UNITS, kernel_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_kernel_xc_id, UNITS, kernel_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_kernel_xc_id, UNITS, kernel_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_kernel_Ha_id, UNITS, kernel_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_kernel_Ha_id, UNITS, kernel_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_b_integral_id, UNITS, integral_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_b_integral_id, UNITS, integral_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_1_coeff_id, UNITS, one_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_1_coeff_id, UNITS, one_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write some data to the file
  !--------------------------------------------------------------------
  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_qvec_id, qvec))


  call check_netcdf( nf90_put_var(nc_id,var_list_m1_id, list_m1))
  call check_netcdf( nf90_put_var(nc_id,var_list_m2_id, list_m2))

  call check_netcdf( nf90_put_var(nc_id,var_re_projector_id,  &
				real(transpose(T_functional_projectors))))

  call check_netcdf( nf90_put_var(nc_id,var_im_projector_id,  &
				aimag(transpose(T_functional_projectors))))



  call check_netcdf( nf90_put_var(nc_id,var_re_kernel_id, real(kernel)))
  call check_netcdf( nf90_put_var(nc_id,var_im_kernel_id, aimag(kernel)))

  call check_netcdf( nf90_put_var(nc_id,var_re_kernel_xc_id, real(kernel_xc)))
  call check_netcdf( nf90_put_var(nc_id,var_im_kernel_xc_id, aimag(kernel_xc)))

  call check_netcdf( nf90_put_var(nc_id,var_re_kernel_Ha_id, real(kernel_Ha)))
  call check_netcdf( nf90_put_var(nc_id,var_im_kernel_Ha_id, aimag(kernel_Ha)))


  call check_netcdf( nf90_put_var(nc_id,var_re_b_integral_id, &
			real(conjg(special_functions_t_coefficients(:,1)))))
  call check_netcdf( nf90_put_var(nc_id,var_im_b_integral_id, &
			aimag(conjg(special_functions_t_coefficients(:,1)))))

  call check_netcdf( nf90_put_var(nc_id,var_re_1_coeff_id, &
					real(special_functions_t_coefficients(:,1))))
  call check_netcdf( nf90_put_var(nc_id,var_im_1_coeff_id, &
					aimag(special_functions_t_coefficients(:,1))))

  call check_netcdf( nf90_put_var(nc_id,var_re_basis_B_coeff_id,real(Pt(:,:))))
  call check_netcdf( nf90_put_var(nc_id,var_im_basis_B_coeff_id,aimag(Pt(:,:))))

  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_functional_projectors_plasmons_netcdf

  subroutine output_overlaps_plasmons_netcdf(overlaps,kernel_xc,kernel_Ha,qvec,L_size)
  !------------------------------------------------------------------------
  ! This subroutine prints out the wannier functions in nectdf format, 
  ! along with relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer           :: L_size
  real(dp)          :: qvec(3)
  complex(dp)       :: overlaps(L_size,L_size)
  complex(dp)       :: kernel_xc(L_size,L_size)
  complex(dp)       :: kernel_Ha(L_size,L_size)

  ! local variables
  complex(dp)       :: kernel(L_size,L_size)


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
 
  integer :: ierr  

  ! overlaps
  integer :: var_re_O_id, var_im_O_id
  character (*), parameter :: re_O_name     = "real_O"
  character (*), parameter :: im_O_name     = "imaginary_O"
  character (*), parameter :: O_units       = "none"

  ! kernels
  integer :: var_re_K_id, var_im_K_id
  character (*), parameter :: re_K_name     = "real_kernel"
  character (*), parameter :: im_K_name     = "imaginary_kernel"
  character (*), parameter :: K_units       = "eV_a0^3"

  integer :: var_re_K_xc_id, var_im_K_xc_id
  character (*), parameter :: re_K_xc_name     = "real_kernel_xc"
  character (*), parameter :: im_K_xc_name     = "imaginary_kernel_xc"

  integer :: var_re_K_Ha_id, var_im_K_Ha_id
  character (*), parameter :: re_K_Ha_name     = "real_kernel_Hartree"
  character (*), parameter :: im_K_Ha_name     = "imaginary_kernel_Hartree"


  ! product basis dimensions
  integer :: dim_L1_id, dim_L1
  integer :: dim_L2_id, dim_L2
  character (*), parameter :: L1_name       = "L1"
  character (*), parameter :: L2_name       = "L2"


  ! 3D space
  integer :: dim_space, dim_space_id
  character (*), parameter :: space_name    = "crystal_basis"


  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"
  character (*), parameter :: qvec_units = "crystal_basis"

  ! Define the units
  character (*), parameter :: UNITS = "units"



  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! prepare the netcdf file, writing meta-data
  ! create file

  call check_netcdf(nf90_create(trim(overlaps_filename), NF90_64BIT_OFFSET, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Bare basis functions overlaps and Hartree-exchange-correlation kernel") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))


  if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the "L1" dimension 
  dim_L1  = L_size
  call check_netcdf( nf90_def_dim(nc_id,L1_name, dim_L1, dim_L1_id) )

  ! create the "L2" dimension 
  dim_L2  = L_size
  call check_netcdf( nf90_def_dim(nc_id,L2_name, dim_L2, dim_L2_id) )

  !------------------------------
  ! define the variables
  !------------------------------

  ! create the "q_vector" array
  call check_netcdf( nf90_def_var(nc_id,qvec_name, NF90_DOUBLE, &
			dim_space_id, var_qvec_id) )

  ! create the "overlaps" array

  call check_netcdf( nf90_def_var(nc_id,re_O_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_re_O_id) )

  call check_netcdf( nf90_def_var(nc_id,im_O_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_im_O_id) )

  ! create the "kernel" array

  call check_netcdf( nf90_def_var(nc_id,re_K_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_re_K_id) )

  call check_netcdf( nf90_def_var(nc_id,im_K_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_im_K_id) )

  !----
  call check_netcdf( nf90_def_var(nc_id,re_K_xc_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_re_K_xc_id) )

  call check_netcdf( nf90_def_var(nc_id,im_K_xc_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_im_K_xc_id) )
  !----
  call check_netcdf( nf90_def_var(nc_id,re_K_Ha_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_re_K_Ha_id) )

  call check_netcdf( nf90_def_var(nc_id,im_K_Ha_name, NF90_DOUBLE, &
     		(/dim_L1_id,dim_L2_id/), var_im_K_Ha_id) )


  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! overlaps
  call check_netcdf( nf90_put_att(nc_id, var_im_O_id, UNITS, O_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_O_id, UNITS, O_UNITS) )

  ! kernel
  call check_netcdf( nf90_put_att(nc_id, var_im_K_id, UNITS, K_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_K_id, UNITS, K_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_im_K_xc_id, UNITS, K_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_K_xc_id, UNITS, K_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_im_K_Ha_id, UNITS, K_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_K_Ha_id, UNITS, K_UNITS) )

  ! q vector
  call check_netcdf( nf90_put_att(nc_id, var_qvec_id, UNITS, qvec_units) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------


  ! write the q vector
  call check_netcdf( nf90_put_var(nc_id, var_qvec_id, qvec))

  ! write the real and imaginary overlaps
  call check_netcdf( nf90_put_var(nc_id, var_re_O_id,real (overlaps)))
  call check_netcdf( nf90_put_var(nc_id, var_im_O_id,aimag(overlaps)))

  ! write the real and imaginary kernel

  kernel(:,:) = kernel_Ha(:,:) + kernel_xc(:,:)

  call check_netcdf( nf90_put_var(nc_id, var_re_K_id,real (kernel)))
  call check_netcdf( nf90_put_var(nc_id, var_im_K_id,aimag(kernel)))

  call check_netcdf( nf90_put_var(nc_id, var_re_K_xc_id,real (kernel_xc)))
  call check_netcdf( nf90_put_var(nc_id, var_im_K_xc_id,aimag(kernel_xc)))

  call check_netcdf( nf90_put_var(nc_id, var_re_K_Ha_id,real (kernel_Ha)))
  call check_netcdf( nf90_put_var(nc_id, var_im_K_Ha_id,aimag(kernel_Ha)))


  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_overlaps_plasmons_netcdf

  subroutine read_overlaps_plasmons_netcdf(L_size,overlaps,kernel,\
							kernel_xc,kernel_Ha,qvec)
  !------------------------------------------------------------------------
  ! This subroutine prints out the wannier functions in nectdf format, 
  ! along with relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf


  use intw_useful_constants
  use intw_utility
  use intw_input_parameters
  use intw_reading
  use intw_W90

  implicit none

  ! input variables
  integer           :: L_size

  !output variables
  complex(dp)       :: overlaps(L_size,L_size)
  complex(dp)       :: kernel(L_size,L_size)
  complex(dp)       :: kernel_xc(L_size,L_size)
  complex(dp)       :: kernel_Ha(L_size,L_size)
  real(dp)          :: qvec(3)


  ! local variables
  real(dp)       :: real_field(L_size,L_size)

  integer :: nc_id        
  ! netcdf file id
 
  integer :: ierr  

  ! q vector
  integer :: var_qvec_id
  character (*), parameter :: qvec_name  = "q_vector"


  ! overlaps
  integer :: var_re_O_id, var_im_O_id
  character (*), parameter :: re_O_name     = "real_O"
  character (*), parameter :: im_O_name     = "imaginary_O"

  ! kernel
  integer :: var_re_K_id, var_im_K_id
  character (*), parameter :: re_K_name     = "real_kernel"
  character (*), parameter :: im_K_name     = "imaginary_kernel"

  integer :: var_re_K_xc_id, var_im_K_xc_id
  character (*), parameter :: re_K_xc_name  = "real_kernel_xc"
  character (*), parameter :: im_K_xc_name  = "imaginary_kernel_xc"

  integer :: var_re_K_Ha_id, var_im_K_Ha_id
  character (*), parameter :: re_K_Ha_name  = "real_kernel_Hartree"
  character (*), parameter :: im_K_Ha_name  = "imaginary_kernel_Hartree"


  !------------------------------
  ! read overlaps
  !------------------------------
  call check_netcdf(nf90_open(trim(overlaps_filename), NF90_NOWRITE, nc_id))


  call check_netcdf(nf90_inq_varid(nc_id,qvec_name,var_qvec_id))
  call check_netcdf(nf90_get_var(nc_id, var_qvec_id,qvec))

  call check_netcdf(nf90_inq_varid(nc_id,re_O_name,var_re_O_id))
  call check_netcdf(nf90_inq_varid(nc_id,im_O_name,var_im_O_id))


  call check_netcdf(nf90_inq_varid(nc_id,re_K_name,var_re_K_id))
  call check_netcdf(nf90_inq_varid(nc_id,im_K_name,var_im_K_id))

  call check_netcdf(nf90_inq_varid(nc_id,re_K_xc_name,var_re_K_xc_id))
  call check_netcdf(nf90_inq_varid(nc_id,im_K_xc_name,var_im_K_xc_id))

  call check_netcdf(nf90_inq_varid(nc_id,re_K_Ha_name,var_re_K_Ha_id))
  call check_netcdf(nf90_inq_varid(nc_id,im_K_Ha_name,var_im_K_Ha_id))


  call check_netcdf(nf90_get_var(nc_id, var_re_O_id,real_field))
  overlaps(:,:) = cmplx_1*real_field(:,:)

  call check_netcdf(nf90_get_var(nc_id, var_im_O_id,real_field))
  overlaps(:,:) = overlaps(:,:) + cmplx_i*real_field(:,:)


  call check_netcdf(nf90_get_var(nc_id, var_re_K_id,real_field))
  kernel(:,:) = cmplx_1*real_field(:,:)

  call check_netcdf(nf90_get_var(nc_id, var_im_K_id,real_field))
  kernel(:,:) = kernel(:,:) + cmplx_i*real_field(:,:)


  call check_netcdf(nf90_get_var(nc_id, var_re_K_xc_id,real_field))
  kernel_xc(:,:) = cmplx_1*real_field(:,:)

  call check_netcdf(nf90_get_var(nc_id, var_im_K_xc_id,real_field))
  kernel_xc(:,:) = kernel_xc(:,:) + cmplx_i*real_field(:,:)


  call check_netcdf(nf90_get_var(nc_id, var_re_K_Ha_id,real_field))
  kernel_Ha(:,:) = cmplx_1*real_field(:,:)

  call check_netcdf(nf90_get_var(nc_id, var_im_K_Ha_id,real_field))
  kernel_Ha(:,:) = kernel_Ha(:,:) + cmplx_i*real_field(:,:)




  call check_netcdf( nf90_close(nc_id) )

  end subroutine read_overlaps_plasmons_netcdf


!----------------------------------------------------------------------------!
!
!
end module intw_product_basis
!
!
!----------------------------------------------------------------------------!

