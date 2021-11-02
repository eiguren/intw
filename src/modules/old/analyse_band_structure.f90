!----------------------------------------------------------------------------!
!       intw project.
!
!----------------------------------------------------------------------------!
module intw_analyse_band_structure
!----------------------------------------------------------------------------!
!
!       This module contains various useful subroutines to analyse and 
!	test various quantities related to the band structure.
!	This module is mostly for post-processing / analysis.
!
!----------------------------------------------------------------------------!

use intw_useful_constants
use intw_input_parameters


  !
  implicit none

contains



  subroutine output_band_structure(band_structure_ks_up,band_structure_ks_down)
  !------------------------------------------------------------------------
  ! This subroutine prints out the band structure in nectdf format, along with
  ! relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_setup

  implicit none

   real(dp)  :: band_structure_ks_up(num_wann,nk1s,nk2s,nk3s)
   real(dp)  :: band_structure_ks_down(num_wann,nk1s,nk2s,nk3s)

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
   
  integer :: var_eps_up_id,var_eps_down_id
  ! variables id: band structure


  integer :: dim_bands_id, dim_k1_id, dim_k2_id, dim_k3_id
  ! dimension id: number of bands, kx,ky,kz

  integer :: dim_bands, dim_k1, dim_k2, dim_k3

  character (*), parameter :: bands_name    = "bands"
  character (*), parameter :: k1_name       = "k1"
  character (*), parameter :: k2_name       = "k2"
  character (*), parameter :: k3_name       = "k3"

  character (*), parameter :: eps_up_name    = "band_structure_spin_up"
  character (*), parameter :: eps_down_name  = "band_structure_spin_down"


  ! Define the units
  character (*), parameter :: UNITS      = "units"
  character (*), parameter :: eps_units  = "eV"
  character (*), parameter :: k_UNITS    = "crystal_basis"


  character(256) :: time_stamp, filename
  integer        :: time_values(8)


  ! create file
  filename = trim('band_structure.nc')

  call check_netcdf(nf90_create(filename, NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Band structure in the 1BZ") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "Coordinates",  &
	"The first brillouin zone is represented on a 3D grid in crystal coordinates"))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note_2",  &
	"Spin is determined approximately using Wannier functions; bands are    "//&
	"initialized to (minus big number) so that they cannot be confused with "//&
	"the Fermi energy"))

  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnetic system; spin is well defined. "))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------
  dim_bands = num_wann
  call check_netcdf( nf90_def_dim(nc_id,bands_name, dim_bands, dim_bands_id) )

  dim_k1 = nk1s
  call check_netcdf( nf90_def_dim(nc_id,k1_name, dim_k1, dim_k1_id) )

  dim_k2 = nk2s
  call check_netcdf( nf90_def_dim(nc_id,k2_name, dim_k2, dim_k2_id) )

  dim_k3 = nk3s
  call check_netcdf( nf90_def_dim(nc_id,k3_name, dim_k3, dim_k3_id) )

  !------------------------------
  ! define the variables
  !------------------------------
  ! create the band structure variables

  call check_netcdf( nf90_def_var(nc_id,eps_up_name, NF90_DOUBLE, &
 	      (/dim_bands_id,dim_k1_id,dim_k2_id,dim_k3_id /), &
						var_eps_up_id) )

  call check_netcdf( nf90_def_var(nc_id,eps_down_name, NF90_DOUBLE, &
 	      (/dim_bands_id,dim_k1_id,dim_k2_id,dim_k3_id /), &
						var_eps_down_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, var_eps_up_id, UNITS, eps_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_eps_down_id, UNITS, eps_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------


  call check_netcdf( nf90_put_var(nc_id, var_eps_up_id,  band_structure_ks_up))
  call check_netcdf( nf90_put_var(nc_id, var_eps_down_id,  band_structure_ks_down))


  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------

 
  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_band_structure

  subroutine determine_spin(wannier_spin_overlaps,U,spin_up,spin_down,psi_up,psi_down)
  !------------------------------------------------------------------------
  ! This subroutine prints out the band structure in nectdf format, along with
  ! relevant information to interpret the data
  !------------------------------------------------------------------------

  use intw_useful_constants
  use intw_setup

  implicit none

   ! input variables
   complex(dp) :: wannier_spin_overlaps(num_wann,num_wann,nspin)
   complex(dp) :: U(num_wann,num_wann)

   ! output variables
   logical     :: spin_up(num_wann)
   logical     :: spin_down(num_wann)
   real(dp)    :: psi_up(num_wann), psi_down(num_wann)

   ! local variables
   integer     :: ipol_up, ipol_down
   integer     :: nb
   real(dp)    :: tol
   complex(dp) :: Ud(num_wann,num_wann)
   complex(dp) :: Mup(num_wann,num_wann)
   complex(dp) :: Mdown(num_wann,num_wann)


   real(dp)    :: psi2_up, psi2_down


   ipol_up   = 1
   ipol_down = 2

   Ud(:,:)  = transpose(conjg(U))

   Mup(:,:)   = matmul(U,matmul(wannier_spin_overlaps(:,:,ipol_up),Ud))
   Mdown(:,:) = matmul(U,matmul(wannier_spin_overlaps(:,:,ipol_down),Ud))


   spin_up(:)     = .false.
   spin_down(:)   = .false.

   do nb = 1 ,num_wann

        psi2_up   = abs(Mup(nb,nb))
        psi2_down = abs(Mdown(nb,nb))


        psi_up(nb)    = sqrt(psi2_up)
        psi_down(nb)  = sqrt(psi2_down)


 	if (psi2_up > psi2_down) then
		spin_up(nb) = .true.
	else
		spin_down(nb) = .true.
	end if


   end do
 
  end subroutine determine_spin


  subroutine output_band_structure2(band_structure_ks,projection_ks_up,projection_ks_down) 
  !------------------------------------------------------------------------
  ! This subroutine prints out the band structure in nectdf format, along with
  ! relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_setup

  implicit none

   real(dp)  :: band_structure_ks(num_wann,nk1s,nk2s,nk3s)
   real(dp)  :: projection_ks_up(num_wann,nk1s,nk2s,nk3s)
   real(dp)  :: projection_ks_down(num_wann,nk1s,nk2s,nk3s)

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
   
  integer :: var_eps_id,var_proj_up_id, var_proj_down_id
  ! variables id: band structure


  integer :: dim_bands_id, dim_k1_id, dim_k2_id, dim_k3_id
  ! dimension id: number of bands, kx,ky,kz

  integer :: dim_bands, dim_k1, dim_k2, dim_k3

  character (*), parameter :: bands_name    = "bands"
  character (*), parameter :: k1_name       = "k1"
  character (*), parameter :: k2_name       = "k2"
  character (*), parameter :: k3_name       = "k3"

  character (*), parameter :: eps_name      = "band_structure"
  character (*), parameter :: proj_up_name  = "projection_spin_up"
  character (*), parameter :: proj_down_name= "projection_spin_down"


  ! Define the units
  character (*), parameter :: UNITS      = "units"
  character (*), parameter :: eps_units  = "eV"
  character (*), parameter :: proj_units = "(unitless)"
  character (*), parameter :: k_UNITS    = "crystal_basis"


  character(256) :: time_stamp, filename
  integer        :: time_values(8)


  ! create file
  filename = trim('band_structure2.nc')

  call check_netcdf(nf90_create(filename, NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Band structure in the 1BZ") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "Coordinates",  &
	"The first brillouin zone is represented on a 3D grid in crystal coordinates"))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note_2",  &
	"Spin is determined approximately using Wannier functions."))

  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnetic system; spin is well defined. "))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------
  dim_bands = num_wann
  call check_netcdf( nf90_def_dim(nc_id,bands_name, dim_bands, dim_bands_id) )

  dim_k1 = nk1s
  call check_netcdf( nf90_def_dim(nc_id,k1_name, dim_k1, dim_k1_id) )

  dim_k2 = nk2s
  call check_netcdf( nf90_def_dim(nc_id,k2_name, dim_k2, dim_k2_id) )

  dim_k3 = nk3s
  call check_netcdf( nf90_def_dim(nc_id,k3_name, dim_k3, dim_k3_id) )

  !------------------------------
  ! define the variables
  !------------------------------
  ! create the band structure variables

  call check_netcdf( nf90_def_var(nc_id,eps_name, NF90_DOUBLE, &
 	      (/dim_bands_id,dim_k1_id,dim_k2_id,dim_k3_id /), &
						var_eps_id) )

  call check_netcdf( nf90_def_var(nc_id,proj_up_name, NF90_DOUBLE, &
 	      (/dim_bands_id,dim_k1_id,dim_k2_id,dim_k3_id /), &
						var_proj_up_id) )



  call check_netcdf( nf90_def_var(nc_id,proj_down_name, NF90_DOUBLE, &
 	      (/dim_bands_id,dim_k1_id,dim_k2_id,dim_k3_id /), &
						var_proj_down_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, var_eps_id, UNITS, eps_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_proj_up_id  , UNITS, proj_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_proj_down_id, UNITS, proj_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------


  call check_netcdf( nf90_put_var(nc_id, var_eps_id      , band_structure_ks))
  call check_netcdf( nf90_put_var(nc_id, var_proj_up_id  , projection_ks_up))
  call check_netcdf( nf90_put_var(nc_id, var_proj_down_id, projection_ks_down))


  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------

 
  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_band_structure2


  subroutine output_nesting_function(nesting_Iq,nesting_integrant_up_ks,    &
						nesting_integrant_down_ksq)
  !------------------------------------------------------------------------
  ! This subroutine outputs the nesting function, as well as its 
  ! integrant for analysis purposes.
  !------------------------------------------------------------------------
  use netcdf

  use intw_useful_constants
  use intw_setup


  implicit none

  real(dp)       :: nesting_Iq

  real(dp)       :: nesting_integrant_up_ks    (nk1s,nk2s,nk3s)
  real(dp)       :: nesting_integrant_down_ksq (nk1s,nk2s,nk3s)

  integer :: nc_id        

  integer :: var_Iq_id, var_q_id, var_integrant_up_id, var_integrant_down_id
  ! variables id: nesting, q vector, integrant up and down

  integer ::  dim_Iq_id, dim_k1_id, dim_k2_id, dim_k3_id, dim_space_id
  ! dimension id: Iq, kx,ky,kz, 3D space

  integer :: dim_Iq, dim_k1, dim_k2, dim_k3, dim_space

  character (*), parameter :: q_name        = "q_vector"
  character (*), parameter :: k1_name       = "k1"
  character (*), parameter :: k2_name       = "k2"
  character (*), parameter :: k3_name       = "k3"
  character (*), parameter :: space_name    = "space_coordinates"
  character (*), parameter :: Iq_name       = "Nesting_function"
  character (*), parameter :: int_up_name   = "integrant_k_up"
  character (*), parameter :: int_dn_name   = "integrant_kq_down"

  ! Define the units
  character (*), parameter :: UNITS      = "units"
  character (*), parameter :: Iq_units   = "1 / eV^2 a_0^3"
  character (*), parameter :: int_units  = "1 / eV^2"
  character (*), parameter :: k_UNITS    = "crystal_basis"
  character (*), parameter :: space_UNITS = "crystal_basis"


  character(256) :: time_stamp, filename
  integer        :: time_values(8)

  ! netcdf file id
  filename = trim('intw_nesting_function')//'.nc'

  call check_netcdf(nf90_create(filename, NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Band structure in the 1BZ") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "Coordinates",  &
	"The first brillouin zone is represented on a 3D grid in crystal coordinates"))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note_2",  &
	"The integrant allows to determine the 'hot spots' on the Fermi surface."))

  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnetic system; spin is well defined. "))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension, 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  dim_Iq = 1
  call check_netcdf( nf90_def_dim(nc_id, Iq_name, dim_Iq, dim_Iq_id) )

  dim_k1 = nk1s
  call check_netcdf( nf90_def_dim(nc_id,k1_name, dim_k1, dim_k1_id) )

  dim_k2 = nk2s
  call check_netcdf( nf90_def_dim(nc_id,k2_name, dim_k2, dim_k2_id) )

  dim_k3 = nk3s
  call check_netcdf( nf90_def_dim(nc_id,k3_name, dim_k3, dim_k3_id) )

  !------------------------------
  ! define the variables
  !------------------------------
  ! create the band structure variables

  ! create the "q-vector" array
  call check_netcdf( nf90_def_var(nc_id,q_name, NF90_DOUBLE, dim_space_id, var_q_id) )

  ! create the Iq array
  call check_netcdf( nf90_def_var(nc_id,Iq_name, NF90_DOUBLE, dim_Iq_id, var_Iq_id) )


  call check_netcdf( nf90_def_var(nc_id,int_up_name, NF90_DOUBLE, &
 	      (/dim_k1_id,dim_k2_id,dim_k3_id /),var_integrant_up_id) )

  call check_netcdf( nf90_def_var(nc_id, int_dn_name, NF90_DOUBLE, &
 	      (/dim_k1_id,dim_k2_id,dim_k3_id /), var_integrant_down_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! q vector
  call check_netcdf( nf90_put_att(nc_id, var_q_id, UNITS, space_UNITS) )

  ! Iq
  call check_netcdf( nf90_put_att(nc_id, var_Iq_id, UNITS, Iq_UNITS) )

  ! integrant
  call check_netcdf( nf90_put_att(nc_id, var_integrant_up_id, UNITS, int_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_integrant_down_id, UNITS, int_UNITS) )


  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_q_id,(/qpt1,qpt2,qpt3/)) )

  call check_netcdf( nf90_put_var(nc_id, var_Iq_id, nesting_Iq))

  call check_netcdf( nf90_put_var(nc_id, var_integrant_up_id, nesting_integrant_up_ks))
  call check_netcdf( nf90_put_var(nc_id, var_integrant_down_id, nesting_integrant_down_ksq ))

  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )


  5   format(A)
  10  format(A,3G18.8E3)
  20  format(F18.8)
 
  end subroutine output_nesting_function


!--------------------------------------------------------------------------------
!
end module intw_analyse_band_structure
!
!--------------------------------------------------------------------------------
