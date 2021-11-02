!----------------------------------------------------------------------------!
!       intw project.
!
!----------------------------------------------------------------------------!
module intw_plane_wave_interpolate
!----------------------------------------------------------------------------!
!
!       This module contains subroutines which will set up the interpolation
!       of plane wave matrix elements, of the type
!
!         O(k;q)  ~  < n1 k | e^{-i (q+G)*r} | n2 k+q >.
!
!        For a given fixed q, the above is periodic with respect to
!        k and can thus be submitted to Fourier interpolation:
!
!                O(k;q)  ~  sum_R e^{ikR} O(R;q).
!
!        For a fixed k, O(k;q) is not periodic with q. However, it should
!        still be a smooth function of k'=k+q, and thus can be interpolated
!        using a 3D neighbor-based scheme.
!
!        Utilizing these ideas should yield a method which is much faster
!        than a 6D neighbor-based interpolation scheme.
!
!
!----------------------------------------------------------------------------!

use intw_useful_constants

use iotk_module
use intw_utility
use intw_reading
use intw_input_parameters
use intw_new_interpolate
use intw_W90
use intw_matrix_elements
use intw_symmetries
use intw_plane_wave_setup

!Peio
!Number of Bloch Original states for the wannierization
use w90_parameters, only: num_bands
!Peio

  !
  implicit none
  !
  save
  !

contains

  subroutine find_X0(ikpt,qpt,X0)
  !----------------------------------------------------------------------------!
  !    Given the index of a k-point in the coarse mesh and a q-point,
  !    this subroutine computes the 3D coordinates of the corner of the 
  !    cube in the coarse grid where kpt+qpt lands.  These coordinates will
  !    be of the form 
  !          X0 = (i1,i2,i3) such that ki = (i-1)/nki for each dimension.
  !    However, ki will not be restrained to be in the 1BZ.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: ikpt
  real(dp)       :: qpt(3)

  ! output variables
  integer        :: X0(3)

  ! local variables
  real(dp)       :: kq(3)

  integer        :: i1, i2, i3

  integer        :: switch


  switch = -1   ! singlet to triplet

  ! find the coordinates of ikpt 
  call switch_indices(nk1,nk2,nk3,ikpt,i1,i2,i3,switch)

  ! find k+q
  kq(1) = dble(i1-1)/dble(nk1)+qpt(1)
  kq(2) = dble(i2-1)/dble(nk2)+qpt(2)
  kq(3) = dble(i3-1)/dble(nk3)+qpt(3)

  ! find the corner in the coarse mesh associated with this k-point
  X0(1) = floor(dble(nk1)*kq(1)+1.0_dp)
  X0(2) = floor(dble(nk2)*kq(2)+1.0_dp)
  X0(3) = floor(dble(nk3)*kq(3)+1.0_dp)

 
  end subroutine find_X0

  subroutine find_t(ikpt,qpt,X0,t)
  !----------------------------------------------------------------------------!
  !    Given the index of a k-point in the coarse mesh and a q-point
  !    as well as the coordinate of the corner where kpt+qpt lands (stored in X0),
  !    this subroutine computes the relative position "t" appropriate for
  !    interpolation.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: ikpt
  real(dp)       :: qpt(3)
  integer        :: X0(3)

  ! output variables
  real(dp)       :: t(3)

  ! local variables
  real(dp)       :: kq(3)

  integer        :: i1, i2, i3
  integer        :: switch


  switch = -1   ! singlet to triplet

  ! find the coordinates of ikpt 
  call switch_indices(nk1,nk2,nk3,ikpt,i1,i2,i3,switch)

  ! find k+q
  kq(1) = dble(i1-1)/dble(nk1)+qpt(1)
  kq(2) = dble(i2-1)/dble(nk2)+qpt(2)
  kq(3) = dble(i3-1)/dble(nk3)+qpt(3)

  ! find t
  t(1)  = dble(nk1)*kq(1)+dble(1-X0(1))
  t(2)  = dble(nk2)*kq(2)+dble(1-X0(2))
  t(3)  = dble(nk3)*kq(3)+dble(1-X0(3))

  end subroutine find_t



  subroutine compute_matrix_elements_at_l_scalar_3D         &
             (N_scheme,D_scheme,sze_scheme,use_IBZ,         &
              ikpt,l_scalar,X_l,W_rotated_matrix_elements,  &
              wfc_1,wfc_2,QE_eig, list_iG_1,list_iG_2)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes all the plane wave matrix elements
  ! corresponding to a given scalar index l_scalar. It is assumed
  ! that this is applied to a 3D scheme (D_scheme = 3) and that
  ! X_l represents the 2nd k-point.
  !------------------------------------------------------------------

  implicit none

  integer        :: ikpt  ! the actual k-point in the coarse mesh

  integer        :: N_scheme,D_scheme,sze_scheme
  integer        :: X_l(sze_scheme,D_scheme)
  logical        :: use_IBZ
  integer        :: l_scalar

  integer        :: i1, i2, i3, ipol, jpol

  integer        :: G2(3), G(3)
  integer        :: list_G(3,nG_shell_max)

  integer        :: iG
  integer        :: ikpt_prime ! neighbor of k+q in the coarse mesh

  integer        :: switch

  complex(dp)    :: W_rotated_matrix_elements(nG_shell_max,num_wann,num_wann,npol,npol)

  complex(dp)    :: wfc_1(nG_max,num_bands,nspin), wfc_2(nG_max,num_bands,nspin)
!  complex(dp)    :: wfc_1(nG_max,nbands,nspin), wfc_2(nG_max,nbands,nspin)

  real(dp)       :: QE_eig(num_bands)
!  real(dp)       :: QE_eig(nbands)

  integer        :: list_iG_1(nG_max), list_iG_2(nG_max) 

  complex(dp)    :: pw_matrix_elements(nG_shell_max,num_bands,num_bands,npol,npol)
!  complex(dp)    :: pw_matrix_elements(nG_shell_max,nbands,nbands,npol,npol)

  switch  =  1 ! triplet to singlet
  W_rotated_matrix_elements = cmplx_0

  ! find ikpt_prime, corresponding to the 2nd k-point 

  i1 = modulo(X_l(l_scalar,1)-1,nk1)+1
  i2 = modulo(X_l(l_scalar,2)-1,nk2)+1
  i3 = modulo(X_l(l_scalar,3)-1,nk3)+1

  G2(1) = (X_l(l_scalar,1)-i1)/nk1
  G2(2) = (X_l(l_scalar,2)-i2)/nk2
  G2(3) = (X_l(l_scalar,3)-i3)/nk3

  call switch_indices(nk1,nk2,nk3,ikpt_prime,i1,i2,i3,switch)


  ! get the wavefunctions
  call get_psi_k(ikpt      ,use_IBZ,list_iG_1,wfc_1,QE_eig)
  call get_psi_k(ikpt_prime,use_IBZ,list_iG_2,wfc_2,QE_eig)

  ! create an array of the G vectors of interest
  do iG = 1,nG_shell_max
        G(:) = gvec(:,iG) +G2(:)
	list_G(:,iG) = G(:)
  end do



  if ( matrix_elements_method == 'FFT' ) then
   call get_plane_wave_matrix_element_FFTW_improved               &
                        (list_G,list_iG_1,list_iG_2, wfc_1,wfc_2, &
			 pw_matrix_elements)

  else if ( matrix_elements_method == 'CONVOLUTION' ) then
  ! loop on G vectors
    do iG = 1,nG_shell_max

        ! identify the G vector that goes in the matrix element
        G(:) = list_G(:,iG)
        ! compute the plane wave matrix elements

!        ORIGINAL CONVOLUTION ROUTINE: 
!                apparently faster than the FFT routines (even FFTW)
!                this is because my tests only involve 1 G vector:
!                FFT may get the advantage for many G shells.
!		omp parallelization

!        call get_plane_wave_matrix_element_convolution      &
!                (G,list_iG_1,list_iG_2, wfc_1,wfc_2,        &
!                            pw_matrix_elements(iG,:,:,:,:) )

!        IMPROVED CONVOLUTION ROUTINE: 
!        this is similar in spirit to the original convolution routine,
!        but utilizes mkl matrix-matrix products. Is this faster?
        call get_plane_wave_matrix_element_convolution_alt  &
                (G,list_iG_1,list_iG_2, wfc_1,wfc_2,        &
                            pw_matrix_elements(iG,:,:,:,:) )

    end do
  end if

  do iG = 1,nG_shell_max
        ! W rotate 
        ! Asier && Idoia 15 07 2014: We add nspin variable as input below.
        call Wrotate_plane_wave_matrix_elements(ikpt,ikpt_prime, u_mesh, &
                                pw_matrix_elements(iG,:,:,:,:),          &
                                W_rotated_matrix_elements(iG,:,:,:,:))

  end do

  end subroutine compute_matrix_elements_at_l_scalar_3D


  subroutine get_matrix_elements_force_constants(nkmesh,W_int_k,W_int_R)
  !----------------------------------------------------------------------!
  !                                                                      !
  !   This  subroutine computes the matrix elements force constants,     !
  !   defined as                                                         !
  !             O(k) = < m k | e^{-i (q+G)r | n k+q >                    !
  !                   == >                                               !
  !             O(R) = 1/N_k sum_{k} e^{i k R} O(k).                     !
  !                                                                      !
  !   These force constants can then be used to Fourier interpolate      !
  !   the matrix elements everywhere!                                    !
  !                                                                      !
  !   NOTE : this can probably be done better with FFT                   !
  !----------------------------------------------------------------------!
  
  
  implicit none

  ! input variables
  integer       ::  nkmesh
  complex(dp)   ::  W_int_k(nkmesh,num_wann,num_wann,nspin,nspin)

  ! output variables
  complex(dp)   ::  W_int_R(nkmesh,num_wann,num_wann,nspin,nspin)

  ! local variables
  integer  :: switch

  integer  :: ikpt, ik1, ik2, ik3
  integer  :: irpt, ir1, ir2, ir3

  real(dp)   :: rdotk
  complex(dp)::  fac

  switch = -1 ! singlet to triplet

  W_int_R(:,:,:,:,:) = cmplx_0

  do ikpt = 1,nkmesh
        call switch_indices(nk1,nk2,nk3,ikpt,ik1,ik2,ik3,switch)
        do irpt = 1,nkmesh
                ! The real space vectors are just a simple mesh of
                ! similar shape as the reciprocal space mesh.

                call switch_indices(nk1,nk2,nk3,irpt,ir1,ir2,ir3,switch)
                rdotk  = twopi*(                                    &
                                 dble((ik1-1)*(ir1-1))/dble(nk1) +  &
                                 dble((ik2-1)*(ir2-1))/dble(nk2) +  &
                                 dble((ik3-1)*(ir3-1))/dble(nk3))
 
                 fac    = exp(-cmplx_i*rdotk)
  
                 W_int_R(irpt,:,:,:,:) = W_int_R(irpt,:,:,:,:) +  &
                                     fac*W_int_k(ikpt,:,:,:,:)

        enddo
  enddo

  W_int_R(:,:,:,:,:) = W_int_R(:,:,:,:,:)/dble(nkmesh)


  end subroutine get_matrix_elements_force_constants


  subroutine get_interpolated_W_matrix_elements(nkmesh,kvec,W_int_R,Wk)
  !----------------------------------------------------------------------!
  !                                                                      !
  !   This  subroutine computes the Fourier interpolated matrix elements !
  !   using the force constants.                                         !
  !----------------------------------------------------------------------!
  
  
  implicit none

  ! input variables
  integer       ::  nkmesh
  complex(dp)   ::  W_int_R(nkmesh,num_wann,num_wann)
  real   (dp)   ::  kvec(3)

  ! output
  complex(dp)   ::  Wk(num_wann,num_wann)

  ! local variables
  integer  :: switch

  integer  :: ikpt, ik1, ik2, ik3
  integer  :: irpt, ir1, ir2, ir3
  integer  :: irc1, irc2, irc3

  real(dp)   :: rdotk
  complex(dp)::  fac

  switch  = -1 ! singlet to triplet

  Wk(:,:) = cmplx_0


  do irpt = 1,nkmesh
     ! The real space vectors are just a simple mesh of
     ! similar shape as the reciprocal space mesh.

     call switch_indices(nk1,nk2,nk3,irpt,ir1,ir2,ir3,switch)

     ! coordinates centered at the origin
     irc1 = modulo(ir1-1 + nk1/2,nk1)-nk1/2
     irc2 = modulo(ir2-1 + nk2/2,nk2)-nk2/2
     irc3 = modulo(ir3-1 + nk3/2,nk3)-nk3/2
 
     rdotk = twopi*(dble(irc1)*kvec(1)+  &
                    dble(irc2)*kvec(2)+  &
                    dble(irc3)*kvec(3))

     fac   = exp(cmplx_i*rdotk)
  
     Wk(:,:) = Wk(:,:)+ fac*W_int_R(irpt,:,:)

  enddo


  end subroutine get_interpolated_W_matrix_elements


  subroutine get_matrix_elements_force_constants_FFT(nkmesh,W_k,W_R)
  !----------------------------------------------------------------------!
  !                                                                      !
  !   This  subroutine computes the matrix elements force constants,     !
  !   defined as                                                         !
  !             O(k) = < m k | e^{-i (q+G)r | n k+q >                    !
  !                   == >                                               !
  !             O(R) = 1/N_k sum_{k} e^{i k R} O(k).                     !
  !                                                                      !
  !   The matrix elements actually being interpolated, labelled "W"      !
  !   are Wannier rotated, and thus "smooth" with respect to k.          !
  !   The force constants can then be used to Fourier interpolate        !
  !   the matrix elements everywhere!                                    !
  !                                                                      !
  !   The computation is done using FFTW, which is available through     !
  !   the MKL library.                                                   !
  !                                                                      !
  !       NOTE ON FFTW:                                                  !
  !                                                                      !
  ! According to the mkl manual, the Fourier transform is defined as     !
  !                                                                      !
  ! Z(k1,...,kd) = S \sum_{j1=0}^{n1-1}... W(j1,..,jd)                   !
  !                                            Exp[i delta 2pi jl*kl/nl]    !
  !                                                                      !
  !       with        delta = -1  FORWARD  transform                     !
  !                   delta =  1  BACKWARD transform                     !
  !               S = 1 by default on both sides                         !
  !                                                                      !
  ! Looking at the W90 code suggest the following convention is used:    !
  !                                                                      !
  !       H(R) = 1/Nk \sum_k Exp[-i k*R] H(k)                            !
  !                                                                      !
  !       THUS, THE APPROPRIATE CHOICE FOR INTERPOLATION TO W(R) IS      !
  !       FORWARD TRANSFORMATION AND  S=1/Nk.                            !
  !----------------------------------------------------------------------!

  implicit none
  include 'fftw3.f'
  ! input variables
  integer       ::  nkmesh
  complex(dp)   ::  W_k(num_wann,num_wann,npol,npol,nkmesh)
  complex(dp)   ::  W_R(num_wann,num_wann,npol,npol,nkmesh)

  complex(dp)   ::  W_k_work(nk1,nk2,nk3,num_wann,num_wann,npol,npol)
  complex(dp)   ::  W_R_work(nk1,nk2,nk3,num_wann,num_wann,npol,npol)

  ! local variables
  integer       ::  n_bnd, m_bnd, ipol, jpol
  integer       :: switch
  integer       :: ikpt, ik1, ik2, ik3

  ! FFTW variables
  integer*8     ::  fftw_forward_plan
  complex(dp)   ::  fftw_in(nk1,nk2,nk3), fftw_out(nk1,nk2,nk3)


  switch = -1 ! singlet to triplet

  ! First, convert to the 3D labelling scheme
  do ikpt = 1,nkmesh
        call switch_indices(nk1,nk2,nk3,ikpt,ik1,ik2,ik3,switch)
        W_k_work(ik1,ik2,ik3,:,:,:,:) = W_k(:,:,:,:,ikpt)
  end do

  ! give me a plan!
  call DFFTW_PLAN_DFT_3D(fftw_forward_plan,nk1,nk2,nk3,fftw_in,fftw_out, &
                                      FFTW_FORWARD, FFTW_MEASURE)

  ! do your beautiful magic, mister FFT
  do jpol =1, npol
   do ipol =1, npol
    do n_bnd = 1,num_wann
     do m_bnd = 1,num_wann

        fftw_in(:,:,:) = W_k_work(:,:,:,m_bnd,n_bnd,ipol,jpol)

          call DFFTW_EXECUTE(fftw_forward_plan)

        W_R_work(:,:,:,m_bnd,n_bnd,ipol,jpol) = fftw_out(:,:,:)/dble(nkmesh)

     end do !m_bnd
    end do !n_bnd
   end do !ipol
  end do !jpol

  ! Return to the singlet labelling scheme
  do ikpt = 1,nkmesh
        call switch_indices(nk1,nk2,nk3,ikpt,ik1,ik2,ik3,switch)
        W_R(:,:,:,:,ikpt) = W_R_work(ik1,ik2,ik3,:,:,:,:)
  end do


  end subroutine get_matrix_elements_force_constants_FFT



  subroutine Fourier_interpolate_matrix_elements                    &
                (fftw_in,fftw_out,fftw_backward_plan,               &
                 nk_vec_max,ikpts_min, ikpts_max,u_ks, u_ksq,       &
                 pw_matrix_elements)
  !----------------------------------------------------------------------------!
  !     This subroutine uses FFT to compute the interpolated matrix elements
  !     on the fine mesh, and then extracts only the relevant terms.
  !----------------------------------------------------------------------------!
  implicit none

  ! input variables
  integer     :: nk_vec_max
  integer     :: ikpts_min, ikpts_max

  complex(dp) :: u_ks (num_wann,num_wann,nk_vec_max)
  complex(dp) :: u_ksq(num_wann,num_wann,nk_vec_max)


  integer*8   :: fftw_backward_plan
  complex(dp) :: fftw_in (nk1s,nk2s,nk3s)
  complex(dp) :: fftw_out(nk1s,nk2s,nk3s)
        
  ! output variables 
  complex(dp) :: pw_matrix_elements(nG_shell_max,num_wann,num_wann,npol,npol,nk_vec_max)

  ! local variables
  complex(dp) :: W_R_work (nk1, nk2, nk3 )
  complex(dp) :: W_Rs_work(nk1s,nk2s,nk3s)
  complex(dp) :: W_ks_work(nk1s,nk2s,nk3s)
  complex(dp) :: W_R(num_wann,num_wann,npol,npol,nk1*nk2*nk3)

  complex(dp) :: tmp(num_wann,num_wann)

  integer     :: iG_1
  integer     :: m_bnd, n_bnd, ipol, jpol
  integer     :: irpt, ir1, ir2, ir3

  integer     :: ikpts, block_ikpts 
  integer     :: ik1s,  ik2s,  ik3s

  integer     :: switch


  switch = -1 ! singlet to triplet

  ! read in the coarse mesh force constants for this G vector

  do iG_1 = 1, nG_shell_max
     call read_WR_from_netcdf_database(iG_1,W_R)
 
     W_Rs_work = cmplx_0
     W_R_work  = cmplx_0
     ! compute interoplated matrix elements using FFT
     do jpol =1, npol
      do ipol =1, npol
       do n_bnd = 1, num_wann
        do m_bnd = 1, num_wann

           do irpt = 1, nk1*nk2*nk3 

              call switch_indices(nk1,nk2,nk3,irpt,ir1,ir2,ir3,switch)
              W_R_work(ir1,ir2,ir3) = W_R(m_bnd,n_bnd,ipol,jpol,irpt)

           end do !irpt

           ! put the force constants on the fine mesh 
           call coarse_to_smooth(nk1, nk2, nk3,  W_R_work,     &
                              nk1s,nk2s,nk3s, W_Rs_work)

           ! Fourier interpolate        
           fftw_in   = W_Rs_work
           CALL DFFTW_EXECUTE(fftw_backward_plan)
           W_ks_work = fftw_out

           ! fish out only the elements that are needed, anti-rotating
           ! along the way
           do ikpts = ikpts_min, ikpts_max
              ! find the triplet index of ks
              call switch_indices(nk1s,nk2s,nk3s,ikpts,ik1s,ik2s,ik3s,switch)

              ! go to an index that starts at 1
              block_ikpts = ikpts - ikpts_min+1

              ! pw_matrix elements now contains the W-ROTATED matrix elements 
              pw_matrix_elements &
                (iG_1,m_bnd,n_bnd,ipol,jpol,block_ikpts) = W_ks_work(ik1s,ik2s,ik3s)
           end do !ikpts 

        end do ! m_bnd
       end do ! n_bnd 
      end do ! ipol
     end do !jpol

     ! ANTI-ROTATE the matrix elements

     do ikpts = ikpts_min, ikpts_max
        block_ikpts = ikpts - ikpts_min+1

        do jpol =1, npol
         do ipol =1, npol
           tmp(:,:) = pw_matrix_elements(iG_1,:,:,ipol,jpol,block_ikpts)

           pw_matrix_elements(iG_1,:,:,ipol,jpol,block_ikpts)  =    &
                      matmul( u_ks(:,:,block_ikpts),                &
                      matmul(tmp, transpose(conjg(u_ksq(:,:,block_ikpts)))))
         end do !ipol
        end do !jpol

     end do !ikpts 


  end do ! iG_1 

  end subroutine Fourier_interpolate_matrix_elements


  subroutine compute_rotated_matrix_elements_3D                    &
                (ikpt,N_scheme,D_scheme,sze_scheme,X_l,use_IBZ,    &
  				      W_rotated_matrix_elements)
  !------------------------------------------------------------------
  ! 
  ! This subroutine performs essentially the same task as 
  ! subroutine create_rotated_matrix_elements_file_3D; the difference
  ! is that the matrix elements are kept in RAM instead of being
  ! written to a temporary .bin file.
  !
  ! This subroutine computes all the plane wave matrix elements
  ! necessary for interpolation for a specified block of the coarse
  ! mesh. It assumes a 3D polynomial interpolation scheme.
  !
  ! If the logical variable magnon = .true., the matrix elements
  ! computed are < down | ... | up >.
  !------------------------------------------------------------------

  implicit none


  ! input variables
  integer  :: ikpt
  integer  :: N_scheme,D_scheme,sze_scheme
  integer  :: X_l(sze_scheme,D_scheme)
  logical  :: use_IBZ

  ! looping variables
  integer  :: l_scalar, iG

  ! variables to read QE folders
  complex(dp)    :: wfc_1(nG_max,num_bands,nspin), wfc_2(nG_max,num_bands,nspin)
!  complex(dp)    :: wfc_1(nG_max,nbands,nspin), wfc_2(nG_max,nbands,nspin)
  real   (dp)    :: QE_eig(num_bands)
!  real   (dp)    :: QE_eig(nbands)
  integer        :: list_iG_1(nG_max), list_iG_2(nG_max) 

  ! arrays containing computed results
  complex(dp)    :: W_rotated_matrix_elements(sze_scheme,nG_shell_max, &
						num_wann,num_wann,npol,npol)

  complex(dp)    :: Wm(nG_shell_max,num_wann,num_wann,npol,npol)



  do l_scalar = 1, sze_scheme


        call compute_matrix_elements_at_l_scalar_3D         &
             (N_scheme,D_scheme,sze_scheme,use_IBZ,         &
              ikpt,l_scalar,X_l,Wm, wfc_1,wfc_2,            &
	      QE_eig,list_iG_1,list_iG_2)

  	W_rotated_matrix_elements(l_scalar,:,:,:,:,:) = Wm(:,:,:,:,:)
        
  end do


  end subroutine compute_rotated_matrix_elements_3D 

  subroutine compute_interpolation_coefficients(N_scheme,D_scheme,sze_scheme, &
                		W_rotated_matrix_elements, iG,c_coefficients)
  !------------------------------------------------------------------
  ! 
  ! This subroutine computes the C coefficients necessary for the 
  ! 3D neighbor interpolation scheme. It must fetch data from 
  ! a file, as not all data can be kept in RAM for a large system. 
  !
  ! It is unclear at this point what is the best mkl BLAS routine
  ! to use to improve performance. 
  !------------------------------------------------------------------
  implicit none


  ! input variables 
  integer :: N_scheme, D_scheme,sze_scheme
  integer :: iG

  complex(dp) :: W_rotated_matrix_elements(sze_scheme,nG_shell_max, &
					num_wann,num_wann,npol,npol)

  ! output variable
  complex(dp) :: c_coefficients(sze_scheme,num_wann,num_wann,npol,npol)
  complex(dp) :: Wm(sze_scheme,num_wann,num_wann,npol,npol)

  ! looping variables 
  integer :: l_scalar
  integer :: nbnd1, nbnd2 
  integer :: ipol, jpol


  real   (dp) :: re_matf(sze_scheme,num_wann), re_matc(sze_scheme,num_wann)
  real   (dp) :: im_matf(sze_scheme,num_wann), im_matc(sze_scheme,num_wann)
  complex(dp) :: matf(sze_scheme,num_wann), matc(sze_scheme,num_wann)


  character      ::  matdescra(6)

  ! extract a "G-slice" of the rotated matrix elements,
  ! ie the elements for every corner for a given G.
  Wm(:,:,:,:,:) =  W_rotated_matrix_elements(:,iG,:,:,:,:)

  ! compute the interpolation coefficients 

  matdescra(1) = 'G'
  matdescra(4) = 'F'

  do jpol= 1,npol
    do ipol= 1,npol
      do nbnd2 = 1,num_wann
          matf = Wm(:,:,nbnd2,ipol,jpol)
            
          ! the matrix is split into real and imaginary because
          ! the somewhat exotic zcoomm is not always available 
          ! (ex. on the DIPC cluster)

          re_matf = real(matf)
          im_matf = aimag(matf)

          call mkl_dcoomm('n',sze_scheme, num_wann, sze_scheme,        &
                     ONE, matdescra,                                   &
                     MA_sparse_array, MA_row_indices, MA_col_indices,  &
                     sparse_sze,                                       &
                     re_matf, sze_scheme,ZERO,re_matc, sze_scheme)

          call mkl_dcoomm('n',sze_scheme, num_wann, sze_scheme,        &
                     ONE, matdescra,                                   &
                     MA_sparse_array, MA_row_indices, MA_col_indices,  &
                     sparse_sze,                                       &
                     im_matf, sze_scheme,ZERO,im_matc, sze_scheme)

!           matc = cmplx_1*re_matc+cmplx_i*im_matc
!          c_coefficients(:,:,nbnd2) = matc

          c_coefficients(:,:,nbnd2,ipol,jpol) = cmplx_1*re_matc+cmplx_i*im_matc
      end do !nbnd2
    end do !ipol
  end do !jpol


  end subroutine compute_interpolation_coefficients


  subroutine create_netcdf_database()
  !------------------------------------------------------------------------
  ! This subroutine will create a nectdf file which will contain
  ! the matrix element interpolation in R and k space, as well as 
  ! check information.
  !------------------------------------------------------------------------
  use netcdf
  implicit none



  character (*), parameter :: filename = 'intw_db.nc'

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id

  !-------------------------------------------------------
  !   variables:     
  !	q_vec,  the q vector at which the calculation is
  !		conducted.
  !
  !	Wk ,  the interpolated matrix elements in k space
  !
  !	WR ,  the force constants in R space associated to Wk
  !
  !	HR ,  the force constants in R space associated to 
  !           the interpolated Hamiltonian.
  !
  !  list_R,  the R vectors, in crystal coordinates and 
  !	      centered at the origin, to which HR and WR 
  !	      correspond
  !
  !  list_G,  the G vectors, in crystal coordinates, at 
  !	      which the matrix elements are evaluated.
  !-------------------------------------------------------

  ! ------------------ q --------------------
  character (*), parameter :: q_name        = "q_vector"

  integer :: var_q_id

  ! ------------------ HR --------------------
  character (*), parameter :: re_HR_name    = "real_HR"
  character (*), parameter :: im_HR_name    = "imaginary_HR"
  integer :: var_re_HR_id, var_im_HR_id

  ! ------------------ Wk --------------------
  character (*), parameter :: re_Wk_name    = "real_Wk"
  character (*), parameter :: im_Wk_name    = "imaginary_Wk"

  integer :: var_re_Wk_id, var_im_Wk_id

  ! ------------------ WR --------------------
  character (*), parameter :: re_WR_name    = "real_WR"
  character (*), parameter :: im_WR_name    = "imaginary_WR"

  integer :: var_re_WR_id, var_im_WR_id

  ! ------------------ list_R --------------------
  character (*), parameter :: Rvecs_name    = "R_vectors"
  integer :: var_Rvecs_id

  ! ------------------ list_G --------------------
  character (*), parameter :: Gvecs_name    = "G_vectors"
  integer :: var_Gvecs_id

  !-------------------------------------------------------
  !   dimensions
  !	space_dimension,  dimension of 3D space (3!)
  !
  !	coarse mesh, the number of coarse k points, which
  !	is also the number of R vectors.
  !
  !	bands,  the number of wannier bands
  !
  !	spins,  the number of spin polarizations considered in the
  !             problem.
  !
  !  list_R,  the R vectors, in crystal coordinates and 
  !	      centered at the origin, to which HR and WR 
  !	      correspond
  !
  !  list_G,  the G vectors, in crystal coordinates, at 
  !	      which the matrix elements are evaluated.
  !-------------------------------------------------------

  !----------------  space ----------------------
  character (*), parameter :: space_name    = "space_coordinates"
  integer :: dim_space_id
  integer :: dim_space

  !----------------  coarse mesh  ---------------
  character (*), parameter :: coarse_name   = "k_mesh_and_R_mesh"
  integer :: dim_coarse_mesh_id
  integer :: dim_coarse_mesh


  !----------------  bands ----------------------
  character (*), parameter :: bands1_name   = "bands_1"
  character (*), parameter :: bands2_name   = "bands_2"

  integer :: dim_bands1_id, dim_bands2_id
  integer :: dim_bands


  !----------------  spins ----------------------
  character (*), parameter :: sigma1_name   = "sigma_1"
  character (*), parameter :: sigma2_name   = "sigma_2"

  integer :: dim_s1_id, dim_s2_id


  !----------------  G vectors ----------------------
  character (*), parameter :: G_name        = "iG"
  integer :: dim_G_id
  integer :: dim_G



  ! Define the units
  character (*), parameter :: UNITS = "units"
  character (*), parameter :: SPIN_COM1   = "first_spin_component"
  character (*), parameter :: SPIN_COM2   = "second_spin_component"

  character (*), parameter :: real_space_UNITS  = "crystal_basis_(R_space)"
  character (*), parameter :: reci_space_UNITS  = "crystal_basis_(K_space)"

  character (*), parameter :: no_UNITS  = "(unitless)"

  character (*), parameter :: energy_UNITS  = "eV"



  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! prepare the netcdf file, writing meta-data

  ! create file

  call check_netcdf(nf90_create(filename, NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"INTW database") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This file contains arrays which are necessary to the computation "//&
	"of the response function, chi0. You should already know what system "//&
        "this data applies to."))


  if (magnon) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a magnon calculation: spin combination up-down "))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation: no spin" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed." ))
  end if

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "file_status", &
			"potentially corrupted or incomplete: do not use to restart!" ))

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension, 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create one "G-vector" dimension, the number of G vectors
  dim_G = nG_shell_max
  call check_netcdf( nf90_def_dim(nc_id,G_name, dim_G, dim_G_id) )

  call check_netcdf( nf90_def_dim(nc_id,sigma1_name, npol, dim_s1_id) )
  call check_netcdf( nf90_def_dim(nc_id,sigma2_name, npol, dim_s2_id) )

  ! create two bands dimensions, the number of wannier bands
  dim_bands = num_wann
  call check_netcdf( nf90_def_dim(nc_id,bands1_name, dim_bands, dim_bands1_id) )
  call check_netcdf( nf90_def_dim(nc_id,bands2_name, dim_bands, dim_bands2_id) )

  ! create the coarse mesh dimension
  dim_coarse_mesh = nk1*nk2*nk3
  call check_netcdf( nf90_def_dim(nc_id,coarse_name, dim_coarse_mesh, dim_coarse_mesh_id))

  !------------------------------
  ! define the variables
  !------------------------------
  ! create the "q-vector" array
  call check_netcdf( nf90_def_var(nc_id,q_name, NF90_DOUBLE, dim_space_id, var_q_id) )

  ! create the R vectors array
  call check_netcdf( nf90_def_var(nc_id,Rvecs_name, NF90_INT, &
			(/dim_space_id,dim_coarse_mesh_id/), var_Rvecs_id) )

  ! create the G vectors array
  call check_netcdf( nf90_def_var(nc_id,Gvecs_name, NF90_INT, &
			(/dim_space_id,dim_G_id/), var_Gvecs_id) )


  ! create the data arrays
  call check_netcdf( nf90_def_var(nc_id,re_Wk_name, NF90_DOUBLE,   &
  (/ dim_bands1_id, dim_bands2_id, dim_s1_id, dim_s2_id,           &
                                   dim_G_id, dim_coarse_mesh_id/), &
				                     var_re_Wk_id) )

  call check_netcdf( nf90_def_var(nc_id,im_Wk_name, NF90_DOUBLE,   &
  (/ dim_bands1_id, dim_bands2_id, dim_s1_id, dim_s2_id,           &
                                   dim_G_id, dim_coarse_mesh_id/), &
						     var_im_Wk_id) )

  call check_netcdf( nf90_def_var(nc_id,re_WR_name, NF90_DOUBLE,   &
  (/ dim_bands1_id, dim_bands2_id, dim_s1_id, dim_s2_id,           &
                                   dim_G_id, dim_coarse_mesh_id/), &
						     var_re_WR_id) )

  call check_netcdf( nf90_def_var(nc_id,im_WR_name, NF90_DOUBLE,   &
  (/ dim_bands1_id, dim_bands2_id, dim_s1_id, dim_s2_id,           &
                                   dim_G_id, dim_coarse_mesh_id/), &
						     var_im_WR_id) )


  call check_netcdf( nf90_def_var(nc_id,re_HR_name, NF90_DOUBLE,   &
  (/ dim_bands1_id, dim_bands2_id, dim_coarse_mesh_id/),           &
				                   var_re_HR_id) )

  call check_netcdf( nf90_def_var(nc_id,im_HR_name, NF90_DOUBLE,   &
  (/ dim_bands1_id, dim_bands2_id, dim_coarse_mesh_id/),           &
				                   var_im_HR_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! q vector
  call check_netcdf( nf90_put_att(nc_id, var_q_id, UNITS, reci_space_UNITS) )

  ! G vectors
  call check_netcdf( nf90_put_att(nc_id, var_Gvecs_id, UNITS, reci_space_UNITS) )

  ! R vectors
  call check_netcdf( nf90_put_att(nc_id, var_Rvecs_id, UNITS, real_space_UNITS) )

  ! Hamiltonian
  call check_netcdf( nf90_put_att(nc_id, var_re_HR_id, UNITS, energy_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_HR_id, UNITS, energy_UNITS) )

  ! matrix elements
  call check_netcdf( nf90_put_att(nc_id, var_re_WR_id, UNITS, no_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_WR_id, UNITS, no_UNITS) )

  call check_netcdf( nf90_put_att(nc_id, var_re_Wk_id, UNITS, no_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_im_Wk_id, UNITS, no_UNITS) )



  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write some basic data to the file
  !--------------------------------------------------------------------

  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_q_id,(/qpt1,qpt2,qpt3/)) )

  ! write the G-vectors
  call check_netcdf( nf90_put_var(nc_id, var_Gvecs_id, gvec(:,1:nG_shell_max) ))


  call check_netcdf( nf90_put_var(nc_id, var_Gvecs_id, gvec(:,1:nG_shell_max) ))


  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )

  end subroutine create_netcdf_database

  subroutine write_HR_to_netcdf_database(num_pack,H_R)
  !------------------------------------------------------------------------
  ! This subroutine will create a nectdf file which will contain
  ! the matrix element interpolation in R and k space, as well as 
  ! check information.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  ! input variables
  integer        ::	num_pack
  complex(dp)    ::     H_R(nk1,nk2,nk3,num_pack)

  ! computation variables
  character (*), parameter :: filename = 'intw_db.nc'

  integer	::	switch, irpt
  integer	::	ir1, ir2, ir3
  integer	::	irc1, irc2, irc3
  integer	::	n_bnd, m_bnd, i_pack

  integer       :: list_R(3,nk1*nk2*nk3)
  complex       :: list_HR(num_wann,num_wann,nk1*nk2*nk3)

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id

  ! ------------------ HR --------------------
  character (*), parameter :: re_HR_name    = "real_HR"
  character (*), parameter :: im_HR_name    = "imaginary_HR"
  integer :: var_re_HR_id, var_im_HR_id

  ! ------------------ list_R --------------------
  character (*), parameter :: Rvecs_name    = "R_vectors"
  integer :: var_Rvecs_id

  ! open the file, which we assume exists 
  call check_netcdf( nf90_open(filename, NF90_WRITE, nc_id) )


  ! Get the id of the relevant variables, based on names.
  call check_netcdf( nf90_inq_varid(nc_id,re_HR_name , var_re_HR_id) )
  call check_netcdf( nf90_inq_varid(nc_id,im_HR_name , var_im_HR_id) )
  call check_netcdf( nf90_inq_varid(nc_id,Rvecs_name , var_Rvecs_id) )


  ! generate and write the R vectors as well as the HR matrices
  switch = 1 ! triplet to singlet
  do ir1 = 1,nk1
     irc1 = modulo(ir1-1 + nk1/2,nk1)-nk1/2
     do ir2 = 1,nk2
        irc2 = modulo(ir2-1 + nk2/2,nk2)-nk2/2
        do ir3 = 1,nk3
           irc3 = modulo(ir3-1 + nk3/2,nk3)-nk3/2

           call switch_indices(nk1,nk2,nk3,irpt,ir1,ir2,ir3,switch)

	   list_R(1,irpt) = irc1
	   list_R(2,irpt) = irc2
	   list_R(3,irpt) = irc3

	   do n_bnd =1,num_wann
              do m_bnd = 1, n_bnd

                 i_pack = m_bnd+((n_bnd-1)*n_bnd)/2


		 list_HR(m_bnd,n_bnd,irpt) =  H_R(ir1,ir2,ir3,i_pack)
		 list_HR(n_bnd,m_bnd,irpt) =  conjg(H_R(ir1,ir2,ir3,i_pack))

               end do
           end do

        end do ! irc3
     end do ! irc2
  end do ! irc1


  call check_netcdf( nf90_put_var(nc_id, var_Rvecs_id, list_R) )
  call check_netcdf( nf90_put_var(nc_id, var_re_HR_id, real(list_HR)) )
  call check_netcdf( nf90_put_var(nc_id, var_im_HR_id, aimag(list_HR)) )

  ! close the file
  call check_netcdf( nf90_close(nc_id) )

  end subroutine write_HR_to_netcdf_database


  subroutine write_Wk_to_netcdf_database(ikpt,W_int_me)
  !------------------------------------------------------------------------
  ! This subroutine will create a nectdf file which will contain
  ! the matrix element interpolation in R and k space, as well as 
  ! check information.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  ! input variables
  integer        ::     ikpt
  complex(dp)    ::     W_int_me(num_wann,num_wann,npol,npol,nG_shell_max)

  ! computation variables
  integer        ::     start(6), count(6)

  character (*), parameter :: filename = 'intw_db.nc'


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id

  ! ------------------ Wk --------------------
  character (*), parameter :: re_Wk_name    = "real_Wk"
  character (*), parameter :: im_Wk_name    = "imaginary_Wk"
  integer :: var_re_Wk_id, var_im_Wk_id



  ! open the file, which we assume exists 
  call check_netcdf( nf90_open(filename, NF90_WRITE, nc_id) )


  ! Get the id of the relevant variables, based on names.
  call check_netcdf( nf90_inq_varid(nc_id,re_Wk_name , var_re_Wk_id) )
  call check_netcdf( nf90_inq_varid(nc_id,im_Wk_name , var_im_Wk_id) )


  ! write the variable
  start(:) = (/1,1,1,1,1,ikpt/)
  count(:) = (/num_wann,num_wann,npol,npol,nG_shell_max,1/)

  call check_netcdf( nf90_put_var(nc_id, var_re_Wk_id, &
				real(W_int_me), start,count))

  call check_netcdf( nf90_put_var(nc_id, var_im_Wk_id, &
				aimag(W_int_me), start,count))


  ! close the file
  call check_netcdf( nf90_close(nc_id) )

  end subroutine write_Wk_to_netcdf_database

  subroutine read_Wk_from_netcdf_database(iG,W_k)
  !------------------------------------------------------------------------
  ! This subroutine reads a G slice of the Wk data in the nectdf database 
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  ! input variables
  integer        ::     iG

  ! output variables
  complex(dp)    ::     W_k(num_wann,num_wann,npol,npol,nk1*nk2*nk3)

  ! computation variables
  integer        ::     start(6), count(6) ! control the data slice in netcdf

  real(dp)       ::     re_W_k(num_wann,num_wann,npol,npol,nk1*nk2*nk3)
  real(dp)       ::     im_W_k(num_wann,num_wann,npol,npol,nk1*nk2*nk3)


  character (*), parameter :: filename = 'intw_db.nc'

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id

  ! ------------------ Wk --------------------
  character (*), parameter :: re_Wk_name    = "real_Wk"
  character (*), parameter :: im_Wk_name    = "imaginary_Wk"
  integer :: var_re_Wk_id, var_im_Wk_id


  ! open the file, which we assume exists 
  call check_netcdf( nf90_open(filename, NF90_NOWRITE, nc_id) )


  ! Get the id of the relevant variables, based on names.
  call check_netcdf( nf90_inq_varid(nc_id,re_Wk_name , var_re_Wk_id) )
  call check_netcdf( nf90_inq_varid(nc_id,im_Wk_name , var_im_Wk_id) )


  ! read the variable

  start(:) = (/1,1,1,1,iG,1/)
  count(:) = (/num_wann,num_wann,npol,npol,1,nk1*nk2*nk3/)


  call check_netcdf( nf90_get_var(nc_id, var_re_Wk_id, &
					re_W_k, start,count))

  call check_netcdf( nf90_get_var(nc_id, var_im_Wk_id, &
					im_W_k, start,count))


  ! close the file
  call check_netcdf( nf90_close(nc_id) )

  W_k = cmplx_1*re_W_k+cmplx_i*im_W_k

  end subroutine read_Wk_from_netcdf_database

  subroutine write_WR_to_netcdf_database(iG,W_R)
  !------------------------------------------------------------------------
  ! This subroutine writes WR to the netcdf database file.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  ! input variables
  integer        ::     iG
  complex(dp)    ::     W_R(num_wann,num_wann,npol,npol,nk1*nk2*nk3)

  ! computation variables
  integer        ::     start(6), count(6)

  character (*), parameter :: filename = 'intw_db.nc'


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id

  ! ------------------ WR --------------------
  character (*), parameter :: re_WR_name    = "real_WR"
  character (*), parameter :: im_WR_name    = "imaginary_WR"
  integer :: var_re_WR_id, var_im_WR_id



  ! open the file, which we assume exists 
  call check_netcdf( nf90_open(filename, NF90_WRITE, nc_id) )


  ! Get the id of the relevant variables, based on names.
  call check_netcdf( nf90_inq_varid(nc_id,re_WR_name , var_re_WR_id) )
  call check_netcdf( nf90_inq_varid(nc_id,im_WR_name , var_im_WR_id) )


  ! write the variable
  start(:) = (/1,1,1,1,iG,1/)
  count(:) = (/num_wann,num_wann,npol,npol,1,nk1*nk2*nk3/)

  call check_netcdf( nf90_put_var(nc_id, var_re_WR_id, &
				real(W_R), start,count))

  call check_netcdf( nf90_put_var(nc_id, var_im_WR_id, &
				aimag(W_R), start,count))


  ! close the file
  call check_netcdf( nf90_close(nc_id) )

  end subroutine write_WR_to_netcdf_database

  subroutine read_WR_from_netcdf_database(iG,W_R)
  !------------------------------------------------------------------------
  ! This subroutine reads WR from the netcdf database file.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  ! input variables
  integer        ::     iG
  complex(dp)    ::     W_R(num_wann,num_wann,npol,npol,nk1*nk2*nk3)

  ! computation variables
  integer        ::     start(6), count(6)


  character (*), parameter :: filename = 'intw_db.nc'
  real(dp)       ::   re_W_R(num_wann,num_wann,npol,npol,nk1*nk2*nk3)
  real(dp)       ::   im_W_R(num_wann,num_wann,npol,npol,nk1*nk2*nk3)


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id

  ! ------------------ WR --------------------
  character (*), parameter :: re_WR_name    = "real_WR"
  character (*), parameter :: im_WR_name    = "imaginary_WR"
  integer :: var_re_WR_id, var_im_WR_id



  ! open the file, which we assume exists 
  call check_netcdf( nf90_open(filename, NF90_NOWRITE, nc_id) )


  ! Get the id of the relevant variables, based on names.
  call check_netcdf( nf90_inq_varid(nc_id,re_WR_name , var_re_WR_id) )
  call check_netcdf( nf90_inq_varid(nc_id,im_WR_name , var_im_WR_id) )


  start(:) = (/1,1,1,1,iG,1/)
  count(:) = (/num_wann,num_wann,npol,npol,1,nk1*nk2*nk3/)

  call check_netcdf( nf90_get_var(nc_id, var_re_WR_id, &
				re_W_R, start,count))

  call check_netcdf( nf90_get_var(nc_id, var_im_WR_id, &
				im_W_R, start,count))

  W_R = cmplx_1*re_W_R+cmplx_i*im_W_R

  ! close the file
  call check_netcdf( nf90_close(nc_id) )
  end subroutine read_WR_from_netcdf_database

  subroutine change_status_of_netcdf_database()
  !------------------------------------------------------------------------
  ! This subroutine changes the status of the netcdf file to confirm
  ! that all operations have been conducted successfully and the file
  ! is (hopefully) safe to use.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  character (*), parameter :: filename = 'intw_db.nc'


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id

  ! open the file, which we assume exists 
  call check_netcdf( nf90_open(filename, NF90_WRITE, nc_id) )

  ! enter into redefine mode, to change the attribute
  call check_netcdf( nf90_redef(nc_id) )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "file_status",  &
		"ok" ))


  ! close the file
  call check_netcdf( nf90_close(nc_id) )

  end subroutine change_status_of_netcdf_database

  subroutine check_netcdf_database(netcdf_file_exists,netcdf_file_ok)
  !------------------------------------------------------------------------
  ! This subroutine checks if the netcdf database file is present,
  ! and if so it checks if it contains the right q point and if it is
  ! corrupted.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  character (*), parameter :: filename = 'intw_db.nc'
  character (*), parameter :: q_name = "q_vector"

  logical      :: netcdf_file_exists, netcdf_file_ok
  character (2):: status_line

  integer      :: var_q_id

  real(dp)     :: file_qpt(3), dq


  integer      :: nc_id        

  ! check if the file is present
  inquire(file=filename,exist=netcdf_file_exists)

  if  (netcdf_file_exists) then
     ! open the file, which we assume exists 
     call check_netcdf( nf90_open(filename, NF90_NOWRITE, nc_id) )

     ! extract the q point
     call check_netcdf( nf90_inq_varid(nc_id,q_name, var_q_id) )
     call check_netcdf( nf90_get_var(nc_id, var_q_id, file_qpt))

     ! get status of file
     call check_netcdf( nf90_get_att(nc_id, NF90_GLOBAL, "file_status",  &
						status_line))
     ! close the file
     call check_netcdf( nf90_close(nc_id) )
  end if


  netcdf_file_ok = .true.
  if ( status_line /= 'ok' ) netcdf_file_ok = .false.


  dq = sqrt((file_qpt(1)-qpt1)**2 + &
            (file_qpt(2)-qpt2)**2 + & 
            (file_qpt(3)-qpt3)**2 )

  if ( dq > eps_8 ) netcdf_file_ok = .false.

  end subroutine check_netcdf_database

!--------------------------------------------------------------------------------
!
end module intw_plane_wave_interpolate
!
!--------------------------------------------------------------------------------
