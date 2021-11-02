!----------------------------------------------------------------------------!
!       intw project.
!
!----------------------------------------------------------------------------!
module intw_density_of_states
!----------------------------------------------------------------------------!
!
!       This module contains all necessary components to compute the 
!       density of states, DOS.
!
!----------------------------------------------------------------------------!

use intw_useful_constants
use intw_input_parameters
use intw_W90
use intw_band_crossing
use intw_kramers_kronig
use intw_chi_0

  !
  implicit none
  !
  save
  !
        real(dp),allocatable :: DOS(:)
        ! The density of states


contains

  subroutine create_DOS()
  !------------------------------------------------------------------------
  ! This subroutine initializes the chi0_intw array
  !------------------------------------------------------------------------

  implicit none

  complex(dp)   :: i_delta

  integer       :: iloop


  ! allocate
  allocate(list_hbar_omega(n_hw))

  allocate(DOS(n_hw))

  
  average_width   = ZERO
  number_of_width = 0

  average_width_deg   = ZERO
  number_of_width_deg = 0

  DOS = ZERO
  ! generate the frequency list
  ! generate this global variable
  dhw  = (hw_max-hw_min)/(dble(n_hw-1))

  do iloop = 1,n_hw
     list_hbar_omega(iloop) = hw_min+dble(iloop-1)*dhw
  end do

  end subroutine create_DOS


  subroutine output_DOS()
  !------------------------------------------------------------------------
  ! This subroutine prints out the DOS
  !------------------------------------------------------------------------
  implicit none


  real(dp)       :: spin_sum_factor

  integer        :: io_unit

  integer        :: iw_loop

  character(256) :: filename


  ! normalize at the very end, to avoid truncation error

  ! In the paramagnetic case, all spin states are degenerate
  ! and there is a factor of 2 coming from a sum on spin. 
  ! This factor should not be present in the case of a non-paramagnetic 
  ! system.
  if (nspin == 1) then
	spin_sum_factor = 2.0_dp
  else 
	spin_sum_factor = 1.0_dp
  end if


  DOS  = spin_sum_factor*DOS/dble(nk1s*nk2s*nk3s)


  io_unit = find_free_unit()

  filename = trim('intw_DOS.dat')


  open(unit = io_unit, file = filename)

  write(io_unit,5) '#=========================================================='
  write(io_unit,5) '#   This file contains the computed DOS       '
  write(io_unit,5) '#=========================================================='
  write(io_unit,5) '#    hw (eV)            DOS  (1/eV)'                
  write(io_unit,5) '#=========================================================='

  do iw_loop = 1,n_hw

     ! write to file
     write(io_unit,20) list_hbar_omega(iw_loop), DOS(iw_loop)

  end do 


  5   format(A)
  20  format(2F18.8)
 

  end subroutine output_DOS


  subroutine compute_adaptive_widths_DOS( nk_vec_max,nk_vec,num_pack,   &
                        		eig_ks, u_ks, dH_ks, adaptive_widths)
  !------------------------------------------------------------------------
  ! This subroutine computes the widths delta_{n1n2}(k) that should be
  ! used at a given k point in the computation of chi0. 
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max, num_pack
  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  
  complex(dp)  :: u_ks (num_wann,num_wann,nk_vec_max)
  complex(dp)  :: dH_ks  (num_pack,3,nk_vec_max)

  ! output variables
  real(dp)     :: adaptive_widths (num_wann,nk_vec_max)

  ! computation variables
  integer      :: degenerate_sets_k (num_wann,num_wann)
  integer      :: number_of_sets_k

  integer      :: f_set_k(num_wann)
  real(dp)     :: gradients_k (num_wann,num_wann,3)

  integer      :: i_set_k, i_subset_k 
  integer      :: n_subset_max_k

  real(dp)     :: gr(3)
  real(dp)     :: dk (3)
  real(dp)     :: delta

  integer      :: k_vec_loop
  integer      :: i_pack

  integer      :: m_bnd, n_bnd


  dk(1)       = twopi/dble(nk1s)
  dk(2)       = twopi/dble(nk2s)
  dk(3)       = twopi/dble(nk3s)

  ! loop on all k-points in the block
  do k_vec_loop = 1,nk_vec

      ! resolve degeneracies 
      call find_degenerate_eigenvalue_sets(       &
                eig_ks (:,k_vec_loop),            &
                degenerate_sets_k,number_of_sets_k)



      call get_gradients_and_occupations(num_pack,     &
                eig_ks  (:,k_vec_loop),                &
                u_ks    (:,:,k_vec_loop),              &
                dH_ks   (:,:,k_vec_loop),              &
                degenerate_sets_k,number_of_sets_k,f_set_k,gradients_k)


      ! loop on sets and build widths

      do i_set_k=1, number_of_sets_k
         call find_n_subset_max(degenerate_sets_k,i_set_k,n_subset_max_k)

         delta = ZERO

         do i_subset_k = 1, n_subset_max_k

             gr(:) = gradients_k(i_set_k,i_subset_k,:)
        
             delta = delta  +           &
             	     (gr(1)*dk(1))**2 + &
		     (gr(2)*dk(2))**2 + &
	             (gr(3)*dk(3))**2

                        
         end do !i_subset_k

         delta = adaptive_width_coeff*sqrt(delta/dble(3*n_subset_max_k))

         do i_subset_k = 1, n_subset_max_k
            m_bnd = degenerate_sets_k(i_set_k,i_subset_k)

            adaptive_widths(m_bnd,k_vec_loop) = delta

            average_width   = average_width + delta
            number_of_width = number_of_width + 1 

            if ( n_subset_max_k /= 1) then
                average_width_deg   = average_width_deg + delta
                number_of_width_deg = number_of_width_deg + 1 
            end if

                       
         end do !i_subset_k

      end do !i_set_k

  end do !k_vec_loop 

  end subroutine compute_adaptive_widths_DOS



  subroutine compute_dDOS(d_DOS, cycle_out,i_hw_min, i_hw_max,xi,delta)
  !------------------------------------------------------------------------
  ! This subroutine dispatches the computation of dDOS to 
  ! the appropriate specialized subroutine corresponding to the 
  ! calculation requested by the user: gaussian smearing, lorentzian
  ! smearing, etc....
  !
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  real(dp)     :: delta  ! width of gaussian
  real(dp)     :: xi     ! e_k-eF

  ! output variables
  real(dp)     :: d_DOS(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max



  if ( broadening_function == 'gaussian' ) then
		call compute_dDOS_gaussian &
			(d_DOS, cycle_out,i_hw_min, i_hw_max,xi,delta)
  else if(broadening_function == 'cold_smearing' ) then
  		call compute_dDOS_cold_smearing &
			(d_DOS, cycle_out,i_hw_min, i_hw_max,xi,delta)
  else 
     write(*,*) '**************************************'
     write(*,*) '*  Broadening method not implemented.*'
     write(*,*) '*     program stops.                 *'
     write(*,*) '**************************************'
     stop
  end if
	


 end subroutine compute_dDOS


  subroutine compute_dDOS_gaussian &
	(d_DOS, cycle_out,i_hw_min, i_hw_max,xi,delta)
  !------------------------------------------------------------------------
  ! Compute the contribution to DOS using gaussian width.
  !------------------------------------------------------------------------
  implicit none

  ! log_small = sqrt( - ln(10^{-16}) )
  ! thus, if |x| > log_small, e^{-x^2} < 10^{-16}, and can be neglected!
  real(dp),parameter  :: log_small = 6.0697085175405858

  ! input variables
  real(dp)     :: delta  ! width of gaussian
  real(dp)     :: xi      ! e_k-eF

  ! output variables
  real(dp)     :: d_DOS(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max


  ! computation variables
  real(dp)     :: e_y2, y, y2
  real(dp)     :: sqrt_pi_delta 
  integer      :: i_hw
  real(dp)     :: const_hw1 , const_hw2 

  sqrt_pi_delta = one/(sqrt_pi*delta)
  ! only compute the contribution to the delta-function
  ! where the function doesn't virtually vanish!
  const_hw1 = xi-hw_min
  const_hw2 = delta*log_small



  i_hw_min  = nint((const_hw1-const_hw2)/dhw)+1
  i_hw_max  = nint((const_hw1+const_hw2)/dhw)+1

  ! test the bounds: if the numbers are outside
  ! the range of validity, cycle out (as long as we are using KK)!

  cycle_out = .false.
  if ( i_hw_min  < 1) then
	i_hw_min = 1
  else if ( i_hw_min > n_hw ) then
 	cycle_out = .true.
	i_hw_min  = 0
  end if

  if ( i_hw_max  > n_hw) then
	i_hw_max = n_hw
  else if ( i_hw_max < 1 ) then
 	cycle_out = .true.
        i_hw_max  = 0
  end if

  if ( cycle_out ) return


  ! only compute the contributions to d_DOS which
  ! are larger than machine-precision
  do i_hw = i_hw_min, i_hw_max

     y    = (xi-list_hbar_omega(i_hw))/delta
     y2   = y**2
     e_y2 = dexp(-y2)

     d_DOS(i_hw) = e_y2*sqrt_pi_delta 
  end do

  end subroutine compute_dDOS_gaussian

  subroutine compute_dDOS_cold_smearing &
	(d_DOS, cycle_out,i_hw_min, i_hw_max,xi,delta)
  !------------------------------------------------------------------------
  ! Compute the contribution to DOS using gaussian width.
  !------------------------------------------------------------------------
  implicit none

  ! delta(x) < 10^-16*delta_max if x < x_minus or x > x_plus
  ! the parameters were obtained numerically.

  real(dp),parameter  :: x_minus = -5.524654195658513
  real(dp),parameter  :: x_plus  =  6.920282881677849


  ! input variables
  real(dp)     :: delta  ! width of gaussian
  real(dp)     :: xi      ! e_k-eF

  ! output variables
  real(dp)     :: d_DOS(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max


  ! computation variables
  real(dp)     :: e_y2, y, y2
  real(dp)     :: sqrt_pi_delta 
  real(dp)     :: sqrt2, one_sqrt2
  integer      :: i_hw
  real(dp)     :: const_hw1 , const_hw2 ,const_hw3


  ! only compute the contribution to the delta-function
  ! where the function doesn't virtually vanish!
  const_hw1 = xi-hw_min
  const_hw2 = delta*x_minus 
  const_hw3 = delta*x_plus

  sqrt2     = sqrt(2.0_dp)
  one_sqrt2 = ONE/sqrt2

  sqrt_pi_delta = one/(sqrt_pi*delta)

  i_hw_min  = nint((const_hw1+const_hw2)/dhw)+1
  i_hw_max  = nint((const_hw1+const_hw3)/dhw)+1

  ! test the bounds: if the numbers are outside
  ! the range of validity, cycle out (as long as we are using KK)!

  cycle_out = .false.
  if ( i_hw_min  < 1) then
	i_hw_min = 1
  else if ( i_hw_min > n_hw ) then
 	cycle_out = .true.
	i_hw_min  = 0
  end if

  if ( i_hw_max  > n_hw) then
	i_hw_max = n_hw
  else if ( i_hw_max < 1 ) then
 	cycle_out = .true.
        i_hw_max  = 0
  end if

  if ( cycle_out ) return


  ! only compute the contributions to d_DOS which
  ! are larger than machine-precision
  do i_hw = i_hw_min, i_hw_max

     y    = (xi-list_hbar_omega(i_hw))/delta
     y2   = (y-one_sqrt2)**2
     e_y2 = dexp(-y2)

     d_DOS(i_hw) = e_y2*(2.0_dp-sqrt2*y)*sqrt_pi_delta 
  end do

  end subroutine compute_dDOS_cold_smearing 


  subroutine update_DOS_omp(nk_vec_max,nk_vec,eig_ks, widths)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to chi0_intw, the
  ! array of coefficients in the MLWF basis.
  !
  ! The details of exactly how the contribution to chi0_intw is calculated
  ! (eg, only imaginary part vs. real and imaginary; gaussian vs. lorentzian
  !  broadening,etc...) is left to another subroutine.
  !
  ! This subroutine will be explicit (ie, not many parts will be hidden
  ! in other utility subroutines) so that we can parallelize on the k
  ! loop using openmp.
  !
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max
  real(dp)     :: eig_ks (num_wann,nk_vec_max)

  real(dp)     :: widths(num_wann,nk_vec_max)

  ! computation variables
  real(dp)     :: delta, xi 

  integer      :: k_vec_loop

  integer      :: kn
  integer      :: nbnd
  real(dp)     :: ek
  real(dp)     :: d_DOS(n_hw)

  integer      :: i_hw, i_hw_min, i_hw_max

  logical      :: cycle_out


  real(dp)     :: time1, time2

  real(dp)  :: DOS_local(n_hw)

  ! combine the loops in k, n1 and n2

  !$omp parallel default(none)                                         &
  !$omp shared(nk_vec,num_wann, eig_ks,chemical_potential)             &
  !$omp shared(widths,nk_vec_max,n_hw,DOS)                             &
  !$omp private(DOS_local,kn,nbnd, ek,delta,d_DOS,cycle_out)           &
  !$omp private(i_hw_min, i_hw_max, k_vec_loop,i_hw,xi)


  DOS_local(:) = 0.0_dp

  !$omp do 
  
  do kn = 1,nk_vec*num_wann

     ! extract k_vec_loop, nbnd1, nbnd2

     nbnd        = modulo(kn-1,num_wann)+1
     k_vec_loop  = (kn-nbnd)/num_wann+1

     ek    =  eig_ks (nbnd,k_vec_loop)
     xi    =  ek-chemical_potential

     delta = widths(nbnd,k_vec_loop)

     ! compute the contribution to DOS. Let the subroutine
     ! decide what method to use!

     call compute_dDOS(d_DOS, cycle_out,i_hw_min, i_hw_max,xi,delta)
     if (cycle_out) cycle

     do i_hw = i_hw_min, i_hw_max
     	   DOS_local(i_hw)  =  DOS_local(i_hw)+d_DOS(i_hw)
     end do !i_hw

  end do  ! kn
  !$omp end do

  !$omp barrier
  ! Next, dump DOS_local into DOS;
  ! each thread should dump all its components here, so don't
  ! parallelize over the loop variables! We must still be in the
  ! "parallel" environment, however, or else DOS_local becomes
  ! undefined
  do i_hw = 1, n_hw
     !$omp atomic
     DOS(i_hw)  = DOS(i_hw) + DOS_local(i_hw)
  end do ! i_hw
  !$omp end parallel

  end subroutine update_DOS_omp




!--------------------------------------------------------------------------------
!
end module intw_density_of_states
!
!--------------------------------------------------------------------------------
