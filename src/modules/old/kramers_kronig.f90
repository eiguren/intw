!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
module intw_kramers_kronig
!----------------------------------------------------------------------------!
!
!       This module contains subroutines which implement the kramers-kronig
!       relations, using FFT to perform the Fourier integrals.
!       
! Using convolutions:
!		F_r(w_i) = \sum_j K[i-j] F_i(w_j)
!	which implies  tilde F_r = tilde K * tilde F_i
!	for "tilde" meaning the fourier transform. 
!	As is well known, the convolution of two functions is the
!	Fourier transform of the product of their Fourier transforms!
!	
!
! Using Fourier integrals:
!	**********************************************************
!	* The method below, which was first implemented, is not  *
!	* optimal or as clever as I thought. It is a big pain    *
!	* to implement, and since there are two numerical        *
!	* integrals, the error is actually larger than with      *
!	* a simpler scheme!                                      *
!	**********************************************************
!	Performing Fourier-type integrals is very well described in        
!	Numerical Recipes, chapter 13. The scheme presented there must
!	however be modified to account for the fact that the function to       
!	be transformed is not smooth. Indeed,
!
!		F_r(w) = int dw' e^{-iw't} L(w-w') * F_i(w')
!
!		===> F_r(t) = L(t) * F_i (t)  (Convolution theorem)
!
!		===> F_r(w) = int dt e^{iwt} F_r(t)
!
!		Where L(t) = sign(t)
!
!	Thus sign(t) * F_i(t) is potentially not smooth at t=0!
!       
!----------------------------------------------------------------------------!

!haritz
use kinds, only: dp
!haritz
use intw_useful_constants


!
implicit none
!
include 'fftw3.f'  

save
!
  ! parameters defining dimensions of the w-space and t-space arrays
  integer  ::  Mfft, Nfft

  ! parameter which separates "small" and "large" values of theta 
  real(dp) ::  th_limit = 0.1

  ! Fourier conjugate variables, t and w (not necessarily physically meaningful)
  real(dp),allocatable :: list_t_FFT(:)
  real(dp),allocatable :: list_w_FFT(:)

  ! the following arrays are made global to be easily accessed and probed
  complex(dp),allocatable :: list_Fi_t(:)
  complex(dp),allocatable :: list_Fr_t(:)

  ! useful arrays
  real(dp),allocatable :: list_theta(:)
  real(dp),allocatable :: list_gamma(:)
  real(dp),allocatable :: list_W_of_theta(:)
  real(dp),allocatable :: list_W_of_gamma(:)

  complex(dp),allocatable :: list_alpha_jmin_of_theta(:)
  complex(dp),allocatable :: list_alpha_jmax_of_theta(:)
  complex(dp),allocatable :: list_alpha_jmin_of_gamma(:)
  complex(dp),allocatable :: list_alpha_jmax_of_gamma(:)

  complex(dp),allocatable :: list_bm2_of_gamma(:)
  complex(dp),allocatable :: list_bm1_of_gamma(:)
  complex(dp),allocatable :: list_b0_of_gamma(:)
  complex(dp),allocatable :: list_b1_of_gamma(:)
  complex(dp),allocatable :: list_b2_of_gamma(:)

  ! parameters for the FFT scheme
  real(dp)   :: Delta_w, tmin, tmax, Delta_t 


  ! parameters for FFTW
  
  !*****************************************************************
  ! For the Fourier integral algorithm, which is obsolete
  !*****************************************************************
  integer*8   :: fftw_t_to_w_plan,  fftw_w_to_t_plan

  complex(dp),allocatable :: fftw_KK_t_to_w_in(:), fftw_KK_t_to_w_out(:)
  complex(dp),allocatable :: fftw_KK_w_to_t_in(:), fftw_KK_w_to_t_out(:)


  !*****************************************************************
  ! For the convolution algorithm, inspired by Thomas Frederiksen
  !*****************************************************************

  integer     :: convolution_size

  integer*8   :: fftw_convolution_plan,  fftw_deconvolution_plan
  complex(dp),allocatable :: fftw_convolution_in(:), fftw_convolution_out(:)
  complex(dp),allocatable :: fftw_deconvolution_in(:), fftw_deconvolution_out(:)

  complex(dp),allocatable :: convolution_kernel(:)

  ! for testing
  complex(dp),allocatable :: kernel_direct(:)



contains




   subroutine setup_FFTW_kramers_kronig_plans()
   !------------------------------------------------------------------ 
   ! This subroutine sets up the environment needed by FFTW.
   !
   ! According to the FFTW documentation :
   !		BACKWARD ==> exp{+iwt} --> FROM t TO w
   !		FORWARD  ==> exp{-iwt} --> FROM w TO t
   !------------------------------------------------------------------ 

   implicit none

   integer   ::  size_of_fft

   size_of_fft = Nfft+1

   allocate(fftw_KK_t_to_w_in (size_of_fft))
   allocate(fftw_KK_t_to_w_out(size_of_fft))

   allocate(fftw_KK_w_to_t_in (size_of_fft))
   allocate(fftw_KK_w_to_t_out(size_of_fft))


   CALL  DFFTW_PLAN_DFT_1D(fftw_t_to_w_plan,size_of_fft,      &
                        fftw_KK_t_to_w_in,fftw_KK_t_to_w_out, &
			FFTW_BACKWARD, FFTW_MEASURE)

   CALL  DFFTW_PLAN_DFT_1D(fftw_w_to_t_plan,size_of_fft,      &
                        fftw_KK_w_to_t_in,fftw_KK_w_to_t_out, &
			FFTW_FORWARD, FFTW_MEASURE)


   end subroutine setup_FFTW_kramers_kronig_plans


   subroutine setup_kramers_kronig_arrays(wmin,wmax,Mfft_input)
   !------------------------------------------------------------------ 
   ! This subroutine sets up the necessary arrays for the Kramers-
   ! Kronig procedure.                               
   !------------------------------------------------------------------ 
   implicit none

   ! inputs
   real(dp)   :: wmin, wmax
   integer    :: Mfft_input

   ! working parameters
   integer    :: factor_M_to_N
   real(dp)   :: w, t

   integer    :: iloop

   factor_M_to_N = 8 ! use this for now; is it optimal?

   Mfft = Mfft_input
   Nfft = factor_M_to_N*Mfft


   Delta_w = (wmax-wmin)/dble(Mfft)

   Delta_t = twopi/(Delta_w*dble(Nfft+1))

   ! The following way of computing tmax avoids overflowing,
   ! which is what happens if you compute Mfft*Nfft (a large integer!) 
   tmax    = 0.5_dp*pi/wmax
   tmax    = tmax*dble(Mfft)
   tmax    = tmax*dble(Nfft)/dble(Nfft+1)

   tmin    = -tmax ! assume symmetric domain!


   allocate(list_w_FFT(Mfft+1))
   allocate(list_gamma(Mfft+1))
   allocate(list_W_of_gamma(Mfft+1))
   allocate(list_bm2_of_gamma(Mfft+1))
   allocate(list_bm1_of_gamma(Mfft+1))
   allocate(list_b0_of_gamma (Mfft+1))
   allocate(list_b1_of_gamma (Mfft+1))
   allocate(list_b2_of_gamma (Mfft+1))
   allocate(list_alpha_jmin_of_gamma(Mfft+1))
   allocate(list_alpha_jmax_of_gamma(Mfft+1))


   allocate(list_t_FFT     (Nfft+1))
   allocate(list_theta     (Nfft+1))
   allocate(list_W_of_theta(Nfft+1))

   allocate(list_Fi_t(Nfft+1))
   allocate(list_Fr_t(Nfft+1))

   allocate(list_alpha_jmin_of_theta(Nfft+1))
   allocate(list_alpha_jmax_of_theta(Nfft+1))

   ! set up the arrays

   t = tmin
   do iloop = 1,Nfft+1
        list_t_FFT(iloop) = t
        t                 = t+Delta_t
   end do

   w = wmin
   do iloop = 1,Mfft+1
        list_w_FFT(iloop) = w
        w                 = w+Delta_w
   end do

   list_theta(:) = list_t_FFT(:)*Delta_w
   list_gamma(:) = list_w_FFT(:)*Delta_t

   end subroutine setup_kramers_kronig_arrays


   subroutine deallocate_kramers_kronig_arrays()
   !------------------------------------------------------------------ 
   ! This is a cleanup subroutine .
   !------------------------------------------------------------------ 
   implicit none

   deallocate(list_w_FFT)
   deallocate(list_gamma)
   deallocate(list_W_of_gamma)
   deallocate(list_alpha_jmin_of_gamma)
   deallocate(list_alpha_jmax_of_gamma)

   deallocate(list_bm2_of_gamma)
   deallocate(list_bm1_of_gamma)
   deallocate(list_b0_of_gamma )
   deallocate(list_b1_of_gamma )
   deallocate(list_b2_of_gamma )

   deallocate(list_t_FFT)
   deallocate(list_theta)
   deallocate(list_W_of_theta)
   deallocate(list_alpha_jmin_of_theta)
   deallocate(list_alpha_jmax_of_theta)


   end subroutine deallocate_kramers_kronig_arrays

   subroutine get_W_trapezoid()
   !------------------------------------------------------------------ 
   ! This subroutine computes W(theta) and W(gamma) (careful, these
   ! arrays have different dimensions...). These arrays are useful
   ! for the FFT based Kramers-Kronig algorithm.
   !------------------------------------------------------------------ 
   implicit none

   integer :: iloop

   real(dp):: th, gma, W

   do iloop = 1, Nfft+1
      th = list_theta(iloop)
      if ( abs(th) < th_limit) then
	call W_th_small_trap(th,W)
      else
	call W_th_large_trap(th,W)
      end if
      list_W_of_theta(iloop) = W
   end do

   do iloop = 1,Mfft+1
      gma = list_gamma(iloop)
      if (abs(gma) < th_limit) then
	call W_th_small_trap(gma,W)
      else
	call W_th_large_trap(gma,W)
      end if
      list_W_of_gamma(iloop) = W
   end do

   end subroutine get_W_trapezoid


   subroutine W_th_small_trap(th,W)
   !------------------------------------------------------------------ 
   ! This function computes the value of W(theta) using the 
   ! trapezoidal scheme when theta is "small".
   !------------------------------------------------------------------ 
   implicit none 

   real(dp),intent(in)   :: th
   real(dp),intent(out)  :: W


   real(dp)              :: th2
   real(dp)              :: th4,th6

   ! parameters
 
   real(dp),parameter    :: c0  =   1.0_dp
   real(dp),parameter    :: c2  =  -1.0_dp/  12.0_dp
   real(dp),parameter    :: c4  =   1.0_dp/  360.0_dp
   real(dp),parameter    :: c6  =  -1.0_dp/20160.0_dp

   th2 = th*th
   !th4 = th2*th2
   !th6 = th4*th2

   W = c0+th2*(c2+th2*(c4+th2*c6))
   !W = c0+th2*c2+th4*c4+th6*c6

   end subroutine W_th_small_trap

   subroutine W_th_large_trap(th,W)
   !------------------------------------------------------------------ 
   ! This function computes the value of W(theta) using the 
   ! trapezoidal scheme when theta is "small".
   !------------------------------------------------------------------ 
   implicit none 

   real(dp),intent(in)   :: th
   real(dp),intent(out)  :: W

   real(dp)              :: th2

   th2 = th*th

   W = 2.0_dp*(1.0_dp-dcos(th))/th2

   end subroutine W_th_large_trap

   subroutine get_alpha_trapezoid()
   !------------------------------------------------------------------ 
   ! This subroutine computes the end-point corrections alpha(-theta) 
   ! and alpha(gamma)  in the trapezoidal scheme.
   !------------------------------------------------------------------ 
   implicit none

   integer :: iloop

   real(dp)    :: th, gma
   complex(dp) :: phase, a, alpha0


   do iloop = 1, Nfft+1
      th = list_theta(iloop)
      if ( abs(th) < th_limit) then
	call alpha0_th_small_trap(th,alpha0)
      else
	call alpha0_th_large_trap(th,alpha0)
      end if

      phase = exp(-cmplx_i*dble(Mfft/2)*th)
      a     = alpha0*phase
      list_alpha_jmax_of_theta(iloop) = a
      list_alpha_jmin_of_theta(iloop) = conjg(a)
   end do

   do iloop = 1,Mfft+1
      gma = list_gamma(iloop)
      if (abs(gma) < th_limit) then
	call alpha0_th_small_trap(gma,alpha0)
      else
	call alpha0_th_large_trap(gma,alpha0)
      end if

      phase = exp(cmplx_i*dble(Nfft/2)*gma)
      a     = alpha0*phase
      list_alpha_jmax_of_gamma(iloop) = a
      list_alpha_jmin_of_gamma(iloop) = conjg(a)
   end do

   end subroutine get_alpha_trapezoid

   subroutine alpha0_th_small_trap(th,alpha0)
   !------------------------------------------------------------------ 
   ! This function computes the value of alpha0(theta) using the 
   ! trapezoidal scheme when theta is "large". This value can then
   ! be used to obtain the appropriate end point corrections for the
   ! Fourier integration scheme.
   !------------------------------------------------------------------ 
   implicit none 
   
   real(dp),    intent(in)  :: th
   complex(dp),intent(out)  :: alpha0 
   
   
   real(dp),parameter :: c0  =-1.0_dp/     2.0_dp
   real(dp),parameter :: c2  = 1.0_dp/    24.0_dp
   real(dp),parameter :: c4  =-1.0_dp/   720.0_dp
   real(dp),parameter :: c6  = 1.0_dp/ 40320.0_dp
   
   real(dp),parameter :: c1  = 1.0_dp/     6.0_dp
   real(dp),parameter :: c3  =-1.0_dp/   120.0_dp
   real(dp),parameter :: c5  = 1.0_dp/  5040.0_dp
   real(dp),parameter :: c7  =-1.0_dp/362880.0_dp


   complex(dp) :: th2
   complex(dp) :: ith
   
   th2 = cmplx_1*th*th 
   ith = cmplx_i*th
   
   alpha0 = c0+th2*(c2+th2*(c4+th2*c6))
   alpha0 = alpha0+ith*(c1+th2*(c3+th2*(c5+th2*c7)))
   
   
   end subroutine alpha0_th_small_trap

   subroutine alpha0_th_large_trap(th,alpha0)
   !------------------------------------------------------------------ 
   ! This function computes the value of alpha0(theta) using the 
   ! trapezoidal scheme when theta is "small". This value can then
   ! be used to obtain the appropriate end point corrections for the
   ! Fourier integration scheme.
   ! This function is defined exactly as in Numerical Recipes.
   !------------------------------------------------------------------ 
   implicit none 

   real(dp),    intent(in)  :: th
   complex(dp),intent(out)  :: alpha0 

   real(dp)                 :: th2


   th2 = th*th

   alpha0  =  (cmplx_i*(th-sin(th))-cmplx_1*(cmplx_1-cos(th)))/th2


   end subroutine alpha0_th_large_trap

   subroutine beta_gma_small_trap(gma,beta_m2,beta_m1,beta_0,beta_1,beta_2)
   !------------------------------------------------------------------ 
   ! This function computes the correction coefficients beta
   ! which implement a sided scheme near t = 0; this is necessary
   ! because the function to be Fourier integrated has a discontinuity
   ! at t = 0!
   !------------------------------------------------------------------ 
   implicit none 

   real(dp),intent(in)      :: gma
   complex(dp),intent(out)  :: beta_m2,beta_m1,beta_0,beta_1,beta_2


   complex(dp) :: gma2

   ! c0_2 = -1/2
   ! c1_2 = -1/6   j
   ! c2_2 =  1/24
   ! c3_2 =  1/120 j
   ! c4_2 = -1/720
   ! c5_2 = -1/5040j
 
   complex(dp),parameter :: c0_2  =(-0.5_dp ,                   0.0_dp )
   complex(dp),parameter :: c1_2  =( 0.0_dp ,-0.1666666666666666667_dp )
   complex(dp),parameter :: c2_2  =( 0.0416666666666666667_dp , 0.0_dp )
   complex(dp),parameter :: c3_2  =( 0.0_dp , 0.0083333333333333333_dp )
   complex(dp),parameter :: c4_2  =(-0.0013888888888888889_dp,  0.0_dp ) 
   complex(dp),parameter :: c5_2  =( 0.0_dp ,-0.0001984126984126984_dp)

   ! c0_0 = -1
   ! c2_0 =  1/12
   ! c4_0 = -1/360
   ! c6_0 =  1/20160

   complex(dp),parameter :: c0_0   =(-1.0_dp                  ,0.0_dp)
   complex(dp),parameter :: c2_0   =( 0.0833333333333333333_dp,0.0_dp)
   complex(dp),parameter :: c4_0   =(-0.0027777777777777778_dp,0.0_dp)
   complex(dp),parameter :: c6_0   =( 0.0000496031746031746_dp,0.0_dp)


   gma2 = cmplx_1*gma*gma


   beta_2  = c0_2+gma*(c1_2+gma*(c2_2+gma*(c3_2+gma*(c4_2+gma*c5_2))))
   beta_1  = -2.0_dp*beta_2

   beta_0  = c0_0+gma2*(c2_0+gma2*(c4_0+gma2*c6_0))

   beta_m1 = conjg(beta_1)
   beta_m2 = conjg(beta_2)


   end subroutine beta_gma_small_trap

   subroutine beta_gma_large_trap(gma,beta_m2,beta_m1,beta_0,beta_1,beta_2)
   !------------------------------------------------------------------ 
   ! This function computes the correction coefficients beta
   ! which implement a sided scheme near t = 0; this is necessary
   ! because the function to be Fourier integrated has a discontinuity
   ! at t = 0!
   !------------------------------------------------------------------ 
   implicit none 

   real(dp),intent(in)      :: gma
   complex(dp),intent(out)  :: beta_m2,beta_m1,beta_0,beta_1,beta_2


   real(dp)    :: gma2
   complex(dp) :: egamma

   gma2    = gma*gma
   egamma  = exp(cmplx_i*gma)

   beta_2  = (egamma-cmplx_1-cmplx_i*gma)/gma2
   beta_1  = -2.0_dp*beta_2
   beta_0  = -2.0_dp*(cmplx_1-cos(gma))/gma2
   beta_m1 = conjg(beta_1)
   beta_m2 = conjg(beta_2)


   end subroutine beta_gma_large_trap

   subroutine get_beta_trapezoid()
   !------------------------------------------------------------------ 
   ! This subroutine computes the beta(gamma) arrays, for the 
   ! trapezoidal scheme.
   !------------------------------------------------------------------ 
   implicit none

   integer :: iloop

   real(dp)    :: gma
   complex(dp) ::  bm2,bm1,b0,b1,b2

   do iloop = 1,Mfft+1
      gma = list_gamma(iloop)
      if (abs(gma) < th_limit) then
        call beta_gma_small_trap(gma,bm2,bm1,b0,b1,b2)
      else
        call beta_gma_large_trap(gma,bm2,bm1,b0,b1,b2)
      end if
      list_bm2_of_gamma(iloop) = bm2
      list_bm1_of_gamma(iloop) = bm1
      list_b0_of_gamma(iloop)  = b0
      list_b1_of_gamma(iloop)  = b1
      list_b2_of_gamma(iloop)  = b2
   end do

   end subroutine get_beta_trapezoid



   subroutine apply_kramers_kronig_FFT(list_Fi_w,list_Fr_w)
   !----------------------------------------------------------------------
   ! This subroutine applies the FFT scheme to compute the Kramers-Kronig
   ! integral equation to Fi, in order to obtain Fr. Fi is not 
   ! necessarily purely imaginary, nor is Fr purely real. This subroutine
   ! requires the setup of all necessary arrays to have been performed.
   !----------------------------------------------------------------------
   implicit none

   ! inputs
   complex(dp)   :: list_Fi_w(Mfft+1), list_Fr_w(Mfft+1)


   complex(dp)   :: Fi_w_jmin, Fi_w_jmax
   complex(dp)   :: Fr_t_jmin, Fr_t_jmax
   complex(dp)   :: Fr_t_m2, Fr_t_m1, Fr_t_0, Fr_t_1, Fr_t_2  

   ! computation variables
   integer       :: iloop
   integer       :: ishift

   integer       :: test_corrections

   ! First, set the in array to zero. This will account for zero padding
   fftw_KK_w_to_t_in(:) = cmplx_0

   ! Second, properly shift the data inside the array, to make sure that
   ! it is in the appropriate order for FFT.
   
   ! put positive frequencies at the begining of the array
   ishift = Mfft/2
   do iloop = 1,Mfft/2+1
      fftw_KK_w_to_t_in(iloop) = list_Fi_w(ishift+iloop)
   end do

   ! put negative frequencies at the end of the array
   ishift = Nfft+1-Mfft/2
   do iloop = 1,Mfft/2
      fftw_KK_w_to_t_in(ishift+iloop) = list_Fi_w(iloop)
   end do

   ! perform the FFT!
   CALL DFFTW_EXECUTE(fftw_w_to_t_plan)

   ! compute F_i(t) by unshifting the FFT data
   list_Fi_t(:) = cmplx_0

   ishift = Nfft/2+1
   do iloop = 1,Nfft/2
      list_Fi_t(iloop) = fftw_KK_w_to_t_out(ishift+iloop) 
   end do

   ishift = Nfft/2
   do iloop = 1,Nfft/2+1
      list_Fi_t(ishift+iloop) = fftw_KK_w_to_t_out(iloop) 
   end do

   ! normalize the array to get the proper Fourier integral
   Fi_w_jmin = list_Fi_w(1)
   Fi_w_jmax = list_Fi_w(Mfft+1)


   test_corrections = 0

   if (test_corrections == 2) then
      ! no corrections at all  (FOR TESTING)
      list_Fi_t(:) = Delta_w/twopi*(list_Fi_t(:))
   else if (test_corrections == 1) then
      ! no end point corrections (FOR TESTING)
      list_Fi_t(:) = Delta_w/twopi*(                 &
		  list_W_of_theta(:)*list_Fi_t(:) )
   else
      ! full correction
      list_Fi_t(:) = Delta_w/twopi*(                              &
			list_W_of_theta(:)*list_Fi_t(:)        +  &
			list_alpha_jmin_of_theta(:)*Fi_w_jmin  +  &
			list_alpha_jmax_of_theta(:)*Fi_w_jmax)
   end if


   ! Set up Fr_t, which is just Fi_t times the sign of t

   list_Fr_t(Nfft/2+1) = cmplx_0
   do iloop = 1,Nfft/2
      list_Fr_t(iloop) =  -list_Fi_t(iloop)
   end do
   do iloop = Nfft/2+2,Nfft+1
      list_Fr_t(iloop) =   list_Fi_t(iloop)
   end do

   ! prepare the array for FFT
   fftw_KK_t_to_w_in(:) = cmplx_0

   ! properly shift the data inside the array, to make sure that
   ! it is in the appropriate order for FFT.

   ! put positive times at the begining of the FFT array
   ishift = Nfft/2
   do iloop = 1,Nfft/2+1
      fftw_KK_t_to_w_in(iloop) = list_Fr_t(ishift+iloop)
   end do

   ! put negative times at the end of the FFT array
   ishift = Nfft/2+1
   do iloop = 1,Nfft/2
      fftw_KK_t_to_w_in(ishift+iloop) = list_Fr_t(iloop)
   end do

   ! perform the FFT!
   CALL DFFTW_EXECUTE(fftw_t_to_w_plan)

   ! compute F_r(w) by unshifting the FFT data, cropping the data
   ! appropriately
   list_Fr_w(:) = cmplx_0

   ishift = Nfft-Mfft/2+1
   do iloop = 1,Mfft/2
      list_Fr_w(iloop) = fftw_KK_t_to_w_out(ishift+iloop) 
   end do

   ishift = Mfft/2
   do iloop = 1,Mfft/2+1
      list_Fr_w(ishift+iloop) = fftw_KK_t_to_w_out(iloop) 
   end do

   ! normalize the array to get the proper Fourier integral
   
   ! end point corrections
   Fr_t_jmin = list_Fr_t(1)
   Fr_t_jmax = list_Fr_t(Nfft+1)

   Fr_t_m2   = list_Fr_t(Nfft/2-1)
   Fr_t_m1   = list_Fr_t(Nfft/2  )
   Fr_t_0    = list_Fr_t(Nfft/2+1)
   Fr_t_1    = list_Fr_t(Nfft/2+2)
   Fr_t_2    = list_Fr_t(Nfft/2+3)

   ! central to sided scheme corrections

   test_corrections = 0

   if (test_corrections == 3) then
      ! no corrections at all  (FOR TESTING)
      list_Fr_w(:) = Delta_t*list_Fr_w(:)
   else if (test_corrections == 2) then
      ! no end point corrections, no central scheme corrections (FOR TESTING)
      list_Fr_w(:) = Delta_t*(                      &
		  list_W_of_gamma(:)*list_Fr_w(:) )
   else if (test_corrections == 1) then
      ! no central scheme corrections (FOR TESTING)
      list_Fr_w(:) = Delta_t*(                              &
		  list_W_of_gamma(:)*list_Fr_w(:)        +  &
	          list_alpha_jmin_of_gamma(:)*Fr_t_jmin  +  &
		  list_alpha_jmax_of_gamma(:)*Fr_t_jmax )
   else 
      ! full corrections
      list_Fr_w(:) = Delta_t*(                              &
		  list_W_of_gamma(:)*list_Fr_w(:)        +  &
	          list_alpha_jmin_of_gamma(:)*Fr_t_jmin  +  &
		  list_alpha_jmax_of_gamma(:)*Fr_t_jmax  +  &
    		  Fr_t_m2*list_bm2_of_gamma(:)           +  &
                  Fr_t_m1*list_bm1_of_gamma(:)           +  &
                  Fr_t_0*list_b0_of_gamma(:)             +  &
                  Fr_t_1*list_b1_of_gamma(:)             +  &
                  Fr_t_2*list_b2_of_gamma(:))
   end if

   end subroutine apply_kramers_kronig_FFT

   subroutine setup_FFTW_kramers_kronig_plans_convolution()
   !------------------------------------------------------------------ 
   ! This subroutine sets up the environment needed by FFTW.
   !
   ! According to the FFTW documentation :
   !		BACKWARD ==> exp{+iwt} --> convolution
   !		FORWARD  ==> exp{-iwt} --> deconvolution
   !------------------------------------------------------------------ 

   implicit none


   allocate(fftw_convolution_in (convolution_size))
   allocate(fftw_convolution_out(convolution_size))
   allocate(fftw_deconvolution_in (convolution_size))
   allocate(fftw_deconvolution_out(convolution_size))

   CALL  DFFTW_PLAN_DFT_1D(fftw_convolution_plan,convolution_size, &
                        fftw_convolution_in,fftw_convolution_out,  &
			FFTW_BACKWARD, FFTW_MEASURE)

   CALL  DFFTW_PLAN_DFT_1D(fftw_deconvolution_plan,convolution_size,   &
                        fftw_deconvolution_in,fftw_deconvolution_out,  &
			FFTW_FORWARD, FFTW_MEASURE)


   end subroutine setup_FFTW_kramers_kronig_plans_convolution


   subroutine build_convolution_kernel()
   !----------------------------------------------------------------------
   ! This subroutine computes the array "convolution_kernel" which contains the
   ! appropriate information to perform the Kramers-Kronig integral
   ! through a convolution.
   !
   ! The discretized Kramers-Kronig integral can be expressed as
   !     R_i = sum_{j} K_{ij} * I_j
   !
   ! where the kernel K is actually of the form K_{ij} = k[i-j].
   ! Thus, the sum can be preformed using a convolution, which can be 
   ! peformed efficiently with FFT!
   !
   !----------------------------------------------------------------------
   implicit none

   complex(dp)  :: one_on_ipi


   integer      :: l, lp1
   complex(dp)  :: kernel
   integer      :: n_hw_local



   n_hw_local = convolution_size/2

   one_on_ipi = -cmplx_i/pi
   ! CAUTION! convolution_size should already have been defined!
   allocate(convolution_kernel(convolution_size))
   convolution_kernel(:)     = cmplx_0

  ! Build the linear kernel in direct space
  do lp1 = 1,n_hw_local
	l = lp1-1
        if (l == 0 ) then
                kernel = cmplx_0
        else if (l == 1 ) then
                kernel = -cmplx_1*2.0_dp*log(2.0_dp)
        else
                kernel =  cmplx_1*                         &
                         (2.0_dp*l*log(abs(ONE*l))         &
                        - (ONE+l)*log(abs(ONE+l))          &
                        + (ONE-l)*log(abs(ONE-l)))
        end if

        convolution_kernel(lp1)     = one_on_ipi * kernel
  end do
  do lp1 = n_hw_local+2,convolution_size

        convolution_kernel(lp1) = -convolution_kernel(convolution_size-lp1+2)

  end do

  allocate (kernel_direct(convolution_size))
  kernel_direct(:) = convolution_kernel(:)

  ! perform the convolution of the kernel, going to Fourier space
   
  ! Properly shift the data inside the array, to make sure that
  ! it is in the appropriate order for FFT.


   ! perform the FFT!
   fftw_convolution_in(:) = convolution_kernel(:)
   CALL DFFTW_EXECUTE(fftw_convolution_plan)
   convolution_kernel(:)  = fftw_convolution_out(:) 


   end subroutine build_convolution_kernel


   subroutine apply_kramers_kronig_FFT_convolution(list_Fi_w,list_Fr_w)
   !----------------------------------------------------------------------
   ! This subroutine performs the convolution using FFTW.
   !----------------------------------------------------------------------
   implicit none

   ! input variables
   complex(dp) :: list_Fi_w(convolution_size)

   ! output variables
   complex(dp) :: list_Fr_w(convolution_size)

   ! local variables
   complex(dp) :: KI(convolution_size)


   ! First, convolve list_Fi_w
   fftw_convolution_in(:) = list_Fi_w(:)
   CALL DFFTW_EXECUTE(fftw_convolution_plan)
   KI(:) = fftw_convolution_out(:)


   ! compute the product of the kernel and the "imaginary part" array
   KI(:) = KI(:)*convolution_kernel(:)

   ! deconvolve the array; remember to normalize!
   fftw_deconvolution_in(:) = KI(:)
   CALL DFFTW_EXECUTE(fftw_deconvolution_plan)
   list_Fr_w(:) = fftw_deconvolution_out(:)/dble(convolution_size)

   end subroutine apply_kramers_kronig_FFT_convolution

   subroutine deallocate_kramers_kronig_convolution()
   !------------------------------------------------------------------ 
   ! This subroutine deallocates arrays which are no longer needed.
   !------------------------------------------------------------------ 

   implicit none

   deallocate(fftw_convolution_in )
   deallocate(fftw_convolution_out)
   deallocate(fftw_deconvolution_in)
   deallocate(fftw_deconvolution_out)


   end subroutine deallocate_kramers_kronig_convolution

!----------------------------------------------------------------------------!
!
!
end module intw_kramers_kronig
!
!
!----------------------------------------------------------------------------!
