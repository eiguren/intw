!----------------------------------------------------------------------------!
!	intw project.
!
!	This module is strongly inspired from the module "input_parameters"
!	found in the QE distribution.
!
!----------------------------------------------------------------------------!
!
module intw_input_parameters
!
!----------------------------------------------------------------------------!
!
!  this module contains the definitions of all input parameters for intw
!  as well as subroutines to read and test the input.
!
!----------------------------------------------------------------------------!
  !
!haritz
  use kinds, only: dp
!  use intw_useful_constants
!haritz
  !
  implicit none
  !
  save
  !
!haritz
  ! variables
  public :: mesh_dir, prefix, projectors_filename , nk1, nk2, nk3, &
            number_G_shells, nG_shell_max, TR_symmetry, use_symmetry, &
            magnon, RPA, max_ram, nk1s, nk2s, nk3s, compute_mmn, &
            compute_amn, intw2W_fullzone, intw2W_method, qpt1, qpt2, &
            qpt3, chemical_potential, hw_min, hw_max, n_hw, &
            delta_width, adaptive_width_coeff, chi0_file_prefix, &
            broadening_method, broadening_function, &
            matrix_elements_method, real_part_chi0, use_KK, &
            skip_KK, el_test, el_number_G_shells, &
            interpolation_scheme, bands_points, crossing_method, &
            bands_file, output_format, spin_projection, &
            super_factor_1, super_factor_2, super_factor_3, &
            qpt_basis_1, qpt_basis_2, qpt_basis_3, qpt_dir_1, &
            qpt_dir_2, qpt_dir_3, wannier_functions_filename, &
            product_basis_filename, overlaps_filename, &
            mz_r_filename, Delta_r_filename, fxc_r_filename, &
            basis_order, plane_wave_combination, tolerance, &
            optimal_basis_size, elph, calc_epmat, ph_dir, &
            dvscf_dir, dvscf_name, data_dir, data_name, qlist, &
            fc_mat, ph_bands_file, ep_mat_file, g_realw_file, &
            g_realw_ext_file, h_realw_file, ph_real_file, &
            ph_fine_file, omega_fine_file, el_fine_file, &
            eigen_fine_file, spin_kb_file, spin_Rw_file, &
            nq1, nq2, nq3, nqirr, fsh_symmetry, fsh_sparse, &
            nfsh, num_modes, newton_iter1, newton_iter2, Emin, Emax

  ! subroutines
  public :: read_input, broadcast_input
  !
  private
!haritz

!----------------------------------------------------------------------------!
!  &input namelist nnput parameters
!----------------------------------------------------------------------------!
!
        character(len=256) :: mesh_dir= 'unassigned'
          ! the directory where the coarse mesh QE calculations are stored. 

        character(len=256) :: prefix = 'unassigned'
          ! the QE prefix for the computations 
        
        character(len=256) :: projectors_filename= 'functional_projectors.nc'
          ! name of the netcdf file which contains the projectors
	  ! which describe the optimal product basis.


        integer :: nk1 = 0, nk2 = 0, nk3 = 0
          ! Monkhorst-Pack mesh indices for the coarse mesh 

        integer :: number_G_shells  = 0
          ! number of complete shells of G vectors should be used
          ! in the computation of the response function. 
          ! These shells will be built from the global array gvec.

        integer :: nG_shell_max  = 0
          ! Number of G vectors in the specified number of shells. 
          ! Also, index of the last G vector in gvec belonging to the
          ! shells. This is true because gvec is sorted by QE. 
          ! nG_shell_max will not actually be read from input, but
          ! it is logically related to number_G_shells, and it will be
          ! useful to have it as a global variable. 

        logical :: TR_symmetry 
          !If TR symmetry is present TR_symmetry=.true.

        logical :: use_symmetry = .true.
          ! should the IBZ be used if possible?

        logical :: magnon = .false.
          ! If true, only the spectral weight part of the uddu and duud components
          ! will be computed, assuming that nspin /= 1.

        logical :: RPA = .false.
          ! If true, the exchange-correlation kernel will be set to zero
          ! and only the Coulomb contribution to the kernel will be used.
          ! This is pertinent for a plasmon calculation.


	real(dp) :: max_ram = 1e9
          ! a rough measure (say up to a factor of 2) of how much ram
          ! the program pw_turbo should use. 

        integer :: nk1s = 0, nk2s = 0, nk3s = 0
          ! Monkhorst-Pack mesh indices for the fine, interpolation mesh 

	logical :: compute_mmn= .true.
          ! If  True, the code produces the $prefix.mmn  and $prefix.eig
          ! files.

	logical :: compute_amn= .true.
          ! If  True, the code produces the $prefix.amn file.

	logical :: intw2W_fullzone = .False.
          ! If  True, the code wil assume that a full zone QE calculation
          ! has been performed and that wavefunctions for every k-point
          ! are available. This is mostly for testing and directly comparing 
          ! the results of intw2W90 and pw2wannier.

        character(256) :: intw2W_method = 'CONVOLUTION' 
          ! What method should be used to compute matrix elements;
          ! CONVOLUTION or FFT? 

        real(dp):: qpt1 = 0.0_dp
        real(dp):: qpt2 = 0.0_dp
        real(dp):: qpt3 = 0.0_dp
          ! coordinates of the q point, in crystal coordinates, for the
          ! response function calculation.

        real(dp)      :: chemical_potential = 0.0_dp
        ! The  value which determines the occupation factors, in eV.

        real(dp)      :: hw_min =0.0_dp
        real(dp)      :: hw_max =0.0_dp
        ! minimum and maximum frequencies for response calculation

        integer       :: n_hw  = 0
        ! how many frequencies to use

        real(dp)      :: delta_width = 0.0_dp
        ! The  "infinitesimal" broadening going in the response function,
        ! in eV.

        real(dp)      :: adaptive_width_coeff = 1.0_dp
        ! The  adaptive width is defined as delta = a_coeff * |de/dk|*dk
        ! where a_coeff = adaptive_width_coeff. This parameter is
        ! unitless

        character(256):: chi0_file_prefix = 'chi0'
          ! the prefix for the name of the file in which the chi0
          ! data is stored.

        character(256) :: broadening_method = 'adaptive' 
          ! What method should be used to integrate chi0;
          ! a fixed width, or an adaptive scheme? 

        character(256) :: broadening_function= 'gaussian' 
          ! What function should be used to broaden the integral?
          ! A lorentzian or a gaussian? The function used for 
          ! the real part in the case "gaussian" will actually be
          ! an hermite polynomial times the gaussian. See appropriate
          ! part of the code for details

        character(256) :: matrix_elements_method = 'CONVOLUTION' 
          ! What method should be used to compute matrix elements
          ! in the interpolation scheme


        character(256) :: real_part_chi0= 'kramers-kronig' 
        logical        :: use_KK        = .true.
        logical        :: skip_KK       = .false.
          ! Should the "real part" (not necessarily real) be computed
          ! directly, using the Kramers-Kronig relations, or not at all?
          ! (not at all is useful for testing purposes)


        logical ::  el_test    = .false. 
          ! empty_lattice_test:
          ! should the code compute the response function for an empty 
          ! lattice of free electrons? This is for testing purposes.
          ! Note that the empty lattice matrix elements CANNOT BE
          ! INTERPOLATED. They are not smooth; they go from 1 to 0
          ! discontinuously at the 1BZ boundary.

        integer ::  el_number_G_shells = 1 
          ! number of complete shells of G vectors that should be used 
          ! to compute the empty lattice response function. Each G
          ! corresponds to a band.
       
        character(256) :: interpolation_scheme = 'linear' 
          ! What method should be used to interpolate matrix elements;
          ! hexaclinear, hexacubic, or none? If none, simply set
          ! all matrix elements to 1.


        integer :: bands_points = 0
          ! For how many points will the band structure be interpolated

        character(256):: crossing_method = 'simple'
          ! For how many points will the band structure be interpolated

        character(256):: bands_file = 'interpolated_bands.dat'
          ! For how many points will the band structure be interpolated

        character(256):: output_format = 'ascii'
          ! How should the response function be output to file?
          ! ascii or netcdf.

        logical ::   spin_projection = .false. 
          ! Indicates if we want to calculate the spin projection by using
          ! Wannier function overlaps



        integer       :: super_factor_1 = 2
        integer       :: super_factor_2 = 2
        integer       :: super_factor_3 = 2
        ! Supercell factors for the generation of the Wannier functions

        real(dp):: qpt_basis_1 = 0.0_dp
        real(dp):: qpt_basis_2 = 0.0_dp
        real(dp):: qpt_basis_3 = 0.0_dp
          ! coordinates of the q point, in crystal coordinates, for the
          ! computation of the Basis function. In theory, this should
          ! be the same as [qpt1,qpt2,qpt3], but it is probably a 
          ! good approximation to neglect the q dependence and set
          ! q=Gamma. It is important to be able to change this and
          ! check the validity of the approximation, hence these new
          ! parameters.

        real(dp):: qpt_dir_1 = 0.0_dp
        real(dp):: qpt_dir_2 = 0.0_dp
        real(dp):: qpt_dir_3 = 0.0_dp
          ! qpt which sets the direction which will be used to 
          ! build the symmetrized plane wave basis.


        character(256):: wannier_functions_filename = 'Wannier_functions.nc'
          ! name of the file where the Wannier functions will be saved/read.

        character(256):: product_basis_filename = 'B_functions.nc'
          ! name of the file where the product basis functions will be saved/read.

        character(256):: overlaps_filename = 'overlaps.nc'
          ! name of the file where the overlaps of product basis functions 
          ! will be saved/read.

        character(256):: mz_r_filename = 'QE_mz_r.out'
          ! name of the ascii file where the ground state magnetization is written
        character(256):: Delta_r_filename = 'QE_Delta_r.out'
          ! name of the ascii file where the ground state xc potential difference
          ! is written
        character(256):: fxc_r_filename = 'QE_fxc_r.out'
          ! name of the ascii file where the ground state xc kernel is written

        character(256):: basis_order = 'diagonal'
          ! How should the basis be ordered? To maximize the wing elements, or the 
	  ! diagonal elements?

        character(256):: plane_wave_combination = 'mz_times_symmetric'
          ! If plane waves are used, should they be combined in symmetric functions,
	  ! or should they be used independently?

        real(dp)     :: tolerance = 1.0e-8
          ! eigenvectors of the overlap matrix with
          ! eigenvalues smaller than tolerance should be discarded 

	integer      :: optimal_basis_size = 3
	  ! size of the optimal basis which will be used to compute magnons. 
 	  ! The size should be at least 3!

       logical ::   elph=.false. 
          !Indicates if we want to calculate something related to ph.

       logical ::   calc_epmat=.false. 
          !Indicates if we want to calculate ep matrix elements 

       character (256) :: ph_dir  ='./', dvscf_dir='./', dvscf_name='dvscf_q'
       character (256) :: data_dir='data_dir'  , data_name='data-file.1'
       character (256) :: qlist='qlist.txt'
       character (256) :: fc_mat='--.fc', ph_bands_file='ph_bands.dat'
            
       character (256) :: ep_mat_file="ep_mat.dat" 
       character (256) :: g_realw_file="g_real.dat"
       character (256) :: g_realw_ext_file="g_real_ext.dat"
       character (256) :: h_realw_file="h_real.dat"
       character (256) :: ph_real_file="ph_real.dat"
       character (256) :: ph_fine_file="ph_fine.dat"
       character (256) :: omega_fine_file="omega_fine.dat"
       character (256) :: el_fine_file="el_fine.dat"
       character (256) :: eigen_fine_file="eigen_fine.dat"
       character (256) :: spin_kb_file="spin_coarse.dat" 
       character (256) :: spin_Rw_file="spin_real.dat"

       integer :: nq1=-1, nq2=-1, nq3=-1, nqirr=-1

       logical :: fsh_symmetry=.true., fsh_sparse=.true.
       integer :: nfsh=401
       integer :: num_modes=101
       integer :: newton_iter1 = 2
       integer :: newton_iter2 = 2
       real(dp) :: Emin =-0.01_dp, Emax= 1000.0_dp

        NAMELIST / input / mesh_dir, prefix, nk1, nk2, nk3,             &
			   number_G_shells, TR_symmetry, use_symmetry,  &
                           magnon, RPA,nk1s, nk2s, nk3s,                &
                           qpt1,qpt2,qpt3, max_ram,projectors_filename, &
                           chemical_potential

        NAMELIST / intw2W / intw2W_fullzone,intw2W_method, compute_mmn, &
			    compute_amn

        NAMELIST / response / interpolation_scheme,  &
                              hw_min, hw_max, n_hw, delta_width,        &
                              chi0_file_prefix, broadening_method,      &
			      broadening_function,real_part_chi0,       &
			      adaptive_width_coeff,el_test,             &
			      output_format,el_number_G_shells,         &
			      matrix_elements_method 
                              
        NAMELIST / fsh_input / nfsh, fsh_symmetry, fsh_sparse, newton_iter1,newton_iter2 , emin, emax,num_modes

        NAMELIST / ph / nq1, nq2, nq3, nqirr,  elph, fc_mat,  ph_dir, data_dir, &
             data_name, dvscf_dir, dvscf_name, ph_bands_file, qlist, ep_mat_file, calc_epmat
                        
        NAMELIST / bands /  bands_points, crossing_method ,bands_file,  &
		            spin_projection 

        NAMELIST / product_basis /                                      &
			super_factor_1,super_factor_2, super_factor_3,  &
			qpt_basis_1,qpt_basis_2,qpt_basis_3,            &
			qpt_dir_1,qpt_dir_2,qpt_dir_3,                  &
			product_basis_filename,overlaps_filename,       &
         		mz_r_filename, Delta_r_filename,fxc_r_filename, &
			wannier_functions_filename, optimal_basis_size, &
			basis_order, plane_wave_combination, tolerance
                            
        
! ----------------------------------------------------------------------
!  END namelist 
! ----------------------------------------------------------------------

contains

  subroutine read_input(read_status)
  !------------------------------------------------------------------
  ! This subroutine reads outdir and prefix from standard input 
  !------------------------------------------------------------------
!haritz
  use intw_useful_constants, only: eps_10
!haritz

  implicit none

  integer         :: ios
  logical         :: read_status

  read_status = .false. 


  write(*,20) '|  Reading namelist from input files is sometimes   |'
  write(*,20) '|  buggy; make sure all namelists were read         |'
  write(*,20) '|  correctly!                                       |'
  write(*,20) '|                                                   |'
  write(*,20) '|  Reading input file from standard input ...       |'
  write(*,20) '|                                                   |'
  write(*,20) '|     namelist             ios                      |'
  write(*,20) '|     --------             ----                     |'


  read( 5, input,        iostat = ios)
  write(*,22) '|       &input             ',ios,'                       |'

  read( 5, intw2W,       iostat = ios)
  write(*,22) '|       &intw2W            ',ios,'                       |'

  read( 5, response,     iostat = ios)
  write(*,22) '|       &response          ',ios,'                       |'

  read( 5, bands,        iostat = ios)
  write(*,22) '|       &bands             ',ios,'                       |'

  read( 5, fsh_input,        iostat = ios)
!haritz
!  write(*,22) '|       &fsh                ',ios,'                       |'
  write(*,22) '|       &fsh_input         ',ios,'                       |'
!haritz

  read( 5, ph,        iostat = ios)
  write(*,22) '|       &ph                ',ios,'                       |'

  read( 5, product_basis,iostat = ios)
  write(*,22) '|       &product_basis     ',ios,'                       |'
  write(*,20) '====================================================='


  !	Test the various read parameters

  if      (     mesh_dir .eq. 'unassigned'  ) then
        read_status = .true. 
        write(*,*) 'MISSING mesh_dir'
  end if

  if      (    prefix .eq. 'unassigned'  ) then
        read_status = .true. 
        write(*,*) 'MISSING prefix'
  end if

  if      (    nk1 .eq. 0 ) then
        read_status = .true. 
        write(*,*) 'MISSING nk1'
  end if

  if      (    nk2 .eq. 0 ) then
        read_status = .true. 
        write(*,*) 'MISSING nk2'
  end if

  if      (    nk3 .eq. 0 ) then
        read_status = .true. 
        write(*,*) 'MISSING nk3'
  end if
  
  if      (    nk1s .eq. 0 ) then
        read_status = .true. 
        write(*,*) 'MISSING nk1s'
  end if
  
  if      (    nk2s .eq. 0 ) then
        read_status = .true. 
        write(*,*) 'MISSING nk2s'
  end if
  
  if      (    nk3s .eq. 0 ) then
        read_status = .true. 
        write(*,*) 'MISSING nk3s'
  end if

  if      (   number_G_shells .eq. 0 ) then
        read_status = .true. 
        write(*,*) 'MISSING number_G_shells'
  end if


  if      (   trim(interpolation_scheme) /= 'none'   .and.   &
              trim(interpolation_scheme) /= 'linear' .and.   &
              trim(interpolation_scheme) /= 'cubic' ) then
        read_status = .true. 
        write(*,*) 'UNKNOWN interpolation scheme'
  end if

  if      (   broadening_method  /= 'adaptive' .and.        &
              broadening_method /= 'fixed' )  then
              read_status = .true. 
        write(*,*) 'UNKNOWN broadening method'
  end if

  if      (   output_format /= 'ascii' .and.        &
              output_format /= 'netcdf' )  then
              read_status = .true. 
        write(*,*) 'UNKNOWN output format'
  end if


  if      (   broadening_method  == 'direct' ) then
        read_status = .true. 
        write(*,*) '*************************************************'
        write(*,*) '* OBSOLESCENCE ERROR:                           *'
        write(*,*) '*************************************************'
        write(*,*) '* The broadening method requested, direct,      *'
        write(*,*) '* is implemented, but has not been updated with *'
        write(*,*) '* respect to many changes that have occured in  *'
        write(*,*) '* the development of the code. The relevant     *'
        write(*,*) '* subroutines are thus most probably inconsis-  *'
        write(*,*) '* tent with the rest of the program. If you     *'
        write(*,*) '* really really want to use this broadening     *'
        write(*,*) '* method, update the code.                      *'
        write(*,*) '*              program stops.                   *'
        write(*,*) '*************************************************'
  end if

  if      (   broadening_function /= 'gaussian'   .and.        &
              broadening_function /= 'lorentzian' .and.        &
              broadening_function /= 'cold_smearing')  then
              read_status = .true. 
        write(*,*) 'UNKNOWN broadening function'
  end if


  if      (   broadening_function == 'lorentzian' )  then
              read_status = .true. 
        write(*,*) '*************************************************'
        write(*,*) '* OBSOLESCENCE ERROR:                           *'
        write(*,*) '*************************************************'
        write(*,*) '* The broadening function requested, lorentizan,*'
        write(*,*) '* is implemented, but has not been updated with *'
        write(*,*) '* respect to many changes that have occured in  *'
        write(*,*) '* the development of the code. The relevant     *'
        write(*,*) '* subroutines are thus most probably inconsis-  *'
        write(*,*) '* tent with the rest of the program. If you     *'
        write(*,*) '* really really want to use this broadening     *'
        write(*,*) '* function, update the code.                    *'
        write(*,*) '*              program stops.                   *'
        write(*,*) '*************************************************'

  end if

        
  if      (   matrix_elements_method /= 'CONVOLUTION'  .and.  &
              matrix_elements_method /= 'FFT' ) then
              read_status = .true. 
        write(*,*) 'UNKNOWN matrix_element_method'
  end if

  if      (   real_part_chi0 /= 'kramers-kronig' .and.        &
              real_part_chi0 /= 'direct'         .and.        &
              real_part_chi0 /= 'none' )  then
              read_status = .true. 
        write(*,*) 'UNKNOWN real_part_chi0'
  end if

  !---------------------------------------------------------------------------
  ! If use_KK = true,  compute the imaginary part 
  ! If use_KK = false, compute the imaginary part  and real part directly
  !
  ! If skip_KK = true,  do not perform the kramers-kronig integral
  ! If use_KK = false,  perform the kramers-kronig integral
  !---------------------------------------------------------------------------
  if      (real_part_chi0 == 'direct') then
	 use_KK  = .false.
	 skip_KK = .true.

  else if (real_part_chi0 == 'kramers-kronig' ) then
	 use_KK  = .true.
	 skip_KK = .false.

  else if (real_part_chi0 == 'none' ) then
	 use_KK  = .true.
	 skip_KK = .true.
  end if
  if (.not. skip_KK) then
     if ( modulo(n_hw,2) /= 1) then
        write(*,*) '*****************************************'
        write(*,*) '* ERROR:                                *'
        write(*,*) '* In order to use the Kramers-Kronig    *'
        write(*,*) '* relations, n_hw must be equal to      *'
        write(*,*) '* an odd number (ie n_hw = 2 N + 1 ).   *'
        write(*,*) '*****************************************'
        read_status = .true. 
     end if
     if ( abs(hw_max+hw_min) > eps_10 .and. abs(hw_min) > eps_10) then
        write(*,*) '*****************************************'
        write(*,*) '* ERROR:                                *'
        write(*,*) '* In order to use the Kramers-Kronig    *'
        write(*,*) '* relations, set                        *'
        write(*,*) '*     hw_min = 0 (program assumes       *'
        write(*,*) '*                 antisymmetry about    *'
        write(*,*) '*                 hw=0)                 *'
        write(*,*) '*             OR                        *'
        write(*,*) '*     hw_min = -hw_max                  *'
        write(*,*) '*****************************************'
        read_status = .true. 
     end if
  end if

  if (	modulo(super_factor_1,2) /= 0 .or.             &
      	modulo(super_factor_2,2) /= 0 .or.             &
      	modulo(super_factor_3,2) /= 0 ) then
        write(*,*) '*****************************************'
        write(*,*) '* ERROR:                                *'
        write(*,*) '* The supercell factors used to build   *'
        write(*,*) '* the Wannier functions should be even. *'
        write(*,*) '*****************************************'

        read_status = .true. 
  end if


  if ( trim(basis_order) /= 'wings' .and. trim(basis_order) /= 'diagonal' ) then
        write(*,*) '*****************************************'
        write(*,*) '* ERROR:                                *'
        write(*,*) '* The basis order can only be "wings"   *'
        write(*,*) '* or "diagonal". Review input!          *'
        write(*,*) '*****************************************'
        read_status = .true. 
  end if

  if (optimal_basis_size < 3 .and. magnon) then
        write(*,*) '*****************************************'
        write(*,*) '* ERROR:                                *'
        write(*,*) '* The optimal basis size should be at   *'
        write(*,*) '* least 3 to accomodate the acoustic    *'
        write(*,*) '* condition.                            *'
        write(*,*) '*****************************************'
        read_status = .true. 
  end if


  if ( trim(plane_wave_combination) /= 'mz_times_symmetric' .and.    &
       trim(plane_wave_combination) /= 'symmetric'          .and.    &
       trim(plane_wave_combination) /= 'all'              ) then
        write(*,*) '*****************************************'
        write(*,*) '* ERROR:                                *'
        write(*,*) '* The keyword "plane_wave_combination"  *'
        write(*,*) '* can only be "symmetric" ,  "all"      *'
        write(*,*) '* or "mz_times_symmetric".              *'
        write(*,*) '*    Review your input.                 *'
        write(*,*) '*****************************************'
        read_status = .true. 
  end if


  if ( read_status  ) then
        write(*,*) "PROBLEM!: the input should be of the form:"
        write(*,*) "&input" 
        write(*,*) "             mesh_dir        = 'directory'"
        write(*,*) "             prefix          = 'prefix'"
        write(*,*) "             nk1             = integer "
        write(*,*) "             nk2             = integer "
        write(*,*) "             nk3             = integer "
        write(*,*) "             TR_symmetry     = T or F  " 
        write(*,*) "             magnon          = T or F  " 
        write(*,*) "             RPA             = T or F  " 
        write(*,*) "             nk1s            = integer"
        write(*,*) "             nk2s            = integer"
        write(*,*) "             nk3s            = integer"
        write(*,*) "             number_G_shells = integer"
        write(*,*) "             qpt1            = real"
        write(*,*) "             qpt2            = real"
        write(*,*) "             qpt3            = real"
        write(*,*) "             max_ram         = real"
        write(*,*) "             projectors_filename = (string)"
        write(*,*) "/" 
        write(*,*) "&intw2W" 
        write(*,*) "             intw2W_fullzone = T or F"
        write(*,*) "             intw2W_method   = CONVOLUTION or FFT"
        write(*,*) "             compute_amn     = T or F"
        write(*,*) "             compute_mmn     = T or F"
        write(*,*) "/" 
        write(*,*) "&response" 
        write(*,*) "          interpolation_scheme  = (none, linear, cubic)"
        write(*,*) "             chemical_potential = (in eV)"
        write(*,*) "                      hw_min    = (in eV)"
        write(*,*) "                      hw_max    = (in eV)"
        write(*,*) "                      n_hw      = integer"
        write(*,*) "                   delta_width  = (in eV)"
        write(*,*) "          matrix_element_method = (CONVOLUTION or FFT)"
        write(*,*) "          real_part_chi0        = (direct, kramers-kronig, none)"
        write(*,*) "          broadening_function   = (lorentzian, gaussian, cold_smearing)"
        write(*,*) "          broadening_method     = (fixed, adaptive)"
        write(*,*) "          delta_width           = (if fixed: in eV)"
        write(*,*) "          adaptive_width_coeff  = (if adaptive: unitless, of order 1)"
        write(*,*) "          chi0_file_prefix      = (string)"
        write(*,*) "          el_number_of_G_shells = (integer number of G shells)"
        write(*,*) "          output_format         = (ascii or netcdf)"
        write(*,*) "/" 
        write(*,*) "&bands" 
        write(*,*) "                   bands_points = (integer)"
        write(*,*) "/"
        write(*,*) "/"  
        write(*,*) "&product_basis" 
        write(*,*) "                 super_factor_1 = (even integer)"
        write(*,*) "                 super_factor_2 = (even integer)"
        write(*,*) "                 super_factor_3 = (even integer)"
        write(*,*) "                    qpt_basis_1 = real"
        write(*,*) "                    qpt_basis_2 = real"
        write(*,*) "                    qpt_basis_3 = real"
        write(*,*) "                    qpt_dir_1   = real"
        write(*,*) "                    qpt_dir_2   = real"
        write(*,*) "                    qpt_dir_3   = real"
        write(*,*) "     product_basis_filename     = (string)"
        write(*,*) "     wannier_functions_filename = (string)"
        write(*,*) "     overlaps_filename          = (string)"
        write(*,*) "     mz_r_filename              = (string)"
        write(*,*) "     Delta_r_filename           = (string)"
        write(*,*) "     fxc_r_filename             = (string)"
        write(*,*) "     basis_order                = wings or diagonal "
        write(*,*) "     optimal_basis_size         = integer >= 3"
        write(*,*) "     plane_wave_combination     = symmetric,all or mz_times_symmetric"
        write(*,*) "     tolerance                  = real"
        write(*,*) "/" 
  end if

  return

20 format(A)
22 format(A,I2,A)

  end subroutine read_input

  subroutine broadcast_input()
  !---------------------------------------------------------------------------
  ! This subroutine broadcasts all the read input to the other processors
  !---------------------------------------------------------------------------
  implicit none
!haritz
!  #ifdef _MPI
!  include 'mpif.h'
!  #endif
#ifdef _MPI
  include 'mpif.h'
!haritz

  integer   :: mpi_master_node


  mpi_master_node = 0

  !------------------- character strings ------------------
  call MPI_BCAST(mesh_dir               ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(prefix                 ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(interpolation_scheme   ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(chi0_file_prefix       ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(broadening_method      ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(broadening_function    ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(real_part_chi0         ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(matrix_elements_method ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(output_format          ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(crossing_method        ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(bands_file             ,256,MPI_CHARACTER,mpi_master_node, MPI_COMM_WORLD,mpi_err)


  !------------------- integers ------------------
  call MPI_BCAST(nk1                ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(nk2                ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(nk3                ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(number_G_shells    ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(nk1s               ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(nk2s               ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(nk3s               ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(n_hw               ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(el_number_G_shells ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(bands_points       ,1,MPI_INTEGER,mpi_master_node, MPI_COMM_WORLD,mpi_err)

  !------------------- logicals ------------------
  call MPI_BCAST(TR_symmetry           ,1,MPI_BYTE,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(use_symmetry          ,1,MPI_BYTE,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(magnon                ,1,MPI_BYTE,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(intw2W_fullzone       ,1,MPI_BYTE,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(intw2W_method         ,1,MPI_BYTE,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(el_test               ,1,MPI_BYTE,mpi_master_node, MPI_COMM_WORLD,mpi_err)


  !------------------- doubles ------------------
  call MPI_BCAST(qpt1                 ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(qpt2                 ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(qpt3                 ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(max_ram              ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(chemical_potential   ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(hw_min               ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(hw_max               ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(delta_width          ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)
  call MPI_BCAST(adaptive_width_coeff ,1,MPI_DOUBLE_PRECISION,mpi_master_node, MPI_COMM_WORLD,mpi_err)

!haritz
#endif
!haritz

  end subroutine broadcast_input


!----------------------------------------------------------------------------!
!
END MODULE intw_input_parameters
!
!----------------------------------------------------------------------------!
