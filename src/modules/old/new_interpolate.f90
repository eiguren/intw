!----------------------------------------------------------------------------!
!       intw project.
!
!----------------------------------------------------------------------------!
module intw_new_interpolate
!----------------------------------------------------------------------------!
!
!       This module contains subroutines which will compute Mm1 for the 
!       hexacubic scheme. 
!
!       Hard coding the array is NOT A GOOD IDEA. It is very large,
!       and compilation takes for ever. 
!
!       The code in this module is heavily inspired from the PYTHON code
!       which I created to perform the same task.
!
!       This module will eventually replace "interpolate". It is still 
!       in a design phase.
!----------------------------------------------------------------------------!

!haritz
use kinds, only: dp
!haritz
use intw_useful_constants

!
  implicit none
  !
  save
  !
  real(dp),allocatable    :: MA(:,:)
  real(dp),allocatable    :: MA_sparse_array(:)
  complex(dp),allocatable :: cmplx_MA(:,:)
  integer ,allocatable    :: MA_row_indices(:)
  integer ,allocatable    :: MA_col_indices(:)

  integer :: sparse_sze



contains

  subroutine allocate_MA(sze_scheme)
  !----------------------------------------------------------------------------!
  !     This subroutine simply allocate MA, which is now a global 
  !     array. I do this to avoid strange seg faults when passing MA around
  !     to subroutines
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: sze_scheme

  allocate(MA(sze_scheme,sze_scheme))
  allocate(cmplx_MA(sze_scheme,sze_scheme))


  end subroutine allocate_MA

  subroutine deallocate_MA()
  !----------------------------------------------------------------------------!
  !	function is obvious.
  !----------------------------------------------------------------------------!
  implicit none
   deallocate(MA)

   deallocate(cmplx_MA)

   if (allocated(MA_sparse_array)) then
      deallocate(MA_sparse_array)
      deallocate(MA_row_indices)
      deallocate(MA_col_indices)
   end if

  end subroutine deallocate_MA
  

  subroutine convert_MA_to_sparse_coo_format(sze_scheme)
  !----------------------------------------------------------------------!
  ! This subroutine converts the MA matrix, which is sparse, to 
  ! the sparse coordinates format, appropriate for sparse BLAS 
  ! routines.
  !----------------------------------------------------------------------!
  implicit none

  integer    :: sze_scheme

  integer    :: i, j, cnt

  ! first, find how many non-zero elements there are 
  sparse_sze = 0
  do j=1, sze_scheme
      do i=1, sze_scheme
        if (abs(MA(i,j)) > eps_8) then
           sparse_sze = sparse_sze + 1
        end if
      end do
  end do

  ! allocate
  allocate(MA_sparse_array(sparse_sze))
  allocate(MA_row_indices(sparse_sze))
  allocate(MA_col_indices(sparse_sze))

  cnt = 0
  ! fill the arrays
  do j=1, sze_scheme
      do i=1, sze_scheme
        if (abs(MA(i,j)) > eps_8) then
           cnt = cnt + 1
           MA_sparse_array(cnt) = MA(i,j)
           MA_row_indices (cnt) = i
           MA_col_indices (cnt) = j
        end if
      end do
  end do

  end subroutine convert_MA_to_sparse_coo_format



  subroutine scalar_to_vector_index(D,vector_index,scalar_index,i_min,i_max,switch)
  !----------------------------------------------------------------------------!
  !     This subroutine goes from scalar to vector indices and vice-versa. 
  !     It is similar in spirit to "switch_indices", the utility subroutine. 
  !     
  !     This subroutine is useful because it insures CONSISTENCY in ordering
  !     vectors.
  !
  !     VARIABLES:
  !             D       : the dimension of the vector index 
  !        vector_index : an integer array of dimension D, the vector index
  !        scalar_index : an integer; the scalar index
  !     i_min / i_max   : integers; the minimum/ maximum value taken by a 
  !                       component of the vector index. This is assumed to be
  !                       the same for every dimension.
  !             switch  : integer which indicates whether we want to go from
  !                       scalar to vector, or vice-versa.
  !
  !     switch = +1   :  vector - to - scalar 
  !
  !     switch = -1   :   scalar - to - vector 
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: D
  integer        :: scalar_index, vector_index(D)
  integer        :: i_min,       i_max
  integer        :: switch 

  ! local variables

  integer        :: N, N_power

  integer        :: dd

  integer        :: left_over  , p



  N = i_max-i_min+1 ! the span of the indices of the vector index


  if (switch .eq. 1) then

        scalar_index = 1
        N_power      = 1       ! N^0

        do dd=1, D
                scalar_index = scalar_index + (vector_index(dd)-i_min)*N_power

                N_power = N_power*N
        end do 

  else if (switch .eq. -1) then
        left_over  = scalar_index - 1
  
        N_power = N**(D-1)

        do dd = D,1,-1
                p = left_over/N_power
                left_over = left_over-p*N_power
                N_power   = N_power/N

                vector_index(dd) = p+i_min
        end do

  end if

  end subroutine scalar_to_vector_index


subroutine  build_little_m(N,der_max,little_m)
  !----------------------------------------------------------------------!
  !      This subroutine computes
  !         little_m[k, der ,p ] =  d^der/d t^der  t^p |_t = k 
  !
  !     for der = 0,...,der_max;
  !         k = 0,1;
  !         p = 0,...,N;
  !
  !     the tensor little_m will thus contain all the information 
  !     necessary to compute the effect of the differentiation scheme
  !     on a single variable.
  !----------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: N         ! the order of the scheme
  integer        :: der_max   ! the maximum order of derivation of a given variable
  integer        :: little_m(0:1,0:der_max,0:N)

  ! local variables
  integer        ::  der, k, p 

  integer        ::  power, factor
  integer        ::  factorial

  do p = 0, N
        do der = 0, der_max 

           if ( der > p ) then
                little_m(0,der,p) = 0
                little_m(1,der,p) = 0
           else
                power  = p-der
                factor = factorial(p)/factorial(p-der)
                do k = 0, 1
                   little_m(k,der,p) = factor*k**power
                end do

           end if
        end do
  end do

end subroutine  build_little_m


subroutine build_derivatives(N,D,derivatives)
  !----------------------------------------------------------------------!
  ! simply returns the arrays representing the derivatives pertaining to
  ! the D-N interpolation scheme, where D is the dimension of space, 
  ! and N the order of the scheme. Only "linear" and "cubic" are 
  ! implemented.
  !----------------------------------------------------------------------!
  implicit none

  integer    :: N, D

  integer    :: i_scalar, i_min, i_max, switch

  integer    :: derivatives(((N+1)/2)**D,D) 

  switch =  -1 ! scalar to vector
  i_min  =   0
  i_max  =   1
 
  do i_scalar = 1, ((N+1)/2)**D
     call scalar_to_vector_index(D,derivatives(i_scalar,:),     &
                               i_scalar,i_min,i_max,switch)
  end do
        
end subroutine build_derivatives

subroutine  build_corners(D, corners)
  !----------------------------------------------------------------------!
  !      This function builds the 2-dimensional array corners, where
  !      each row is a vector of length D representing the position of 
  !      a corner. There are 2^D such rows (the total number of corners).
  !----------------------------------------------------------------------!
  implicit none

  integer    :: D
  integer    :: corners(2**D,D) 

  integer    :: k_scalar, i_min, i_max, switch
  

  switch =  -1 ! scalar to vector
  i_min  =   0
  i_max  =   1

  do k_scalar = 1, 2**D
     call scalar_to_vector_index(D,corners(k_scalar,:),     &
                               k_scalar,i_min,i_max,switch)
  end do

end subroutine  build_corners

subroutine  build_powers(N,D,powers)
  !----------------------------------------------------------------------!
  !      This function builds a 2-dimensional array, where every row is
  !      the vector (p_1,..,p_D), the powers of the corresponding monomer. 
  !      The vectors are ordered according to a scalar index. 
  !----------------------------------------------------------------------!
  implicit none

  integer    :: N,D
  integer    :: powers((N+1)**D,D) 

  integer    :: p_scalar
  integer    :: i_min, i_max, switch

  switch = -1 ! scalar to vector
  i_min  =  0
  i_max  =  N

  do p_scalar = 1, (N+1)**D
        call scalar_to_vector_index(D,powers(p_scalar,:),p_scalar,i_min,i_max,switch)
  end do


end subroutine  build_powers

subroutine build_subM(N,D,corners,powers,der_max,little_m,derivatives,subM)
  !----------------------------------------------------------------------!
  ! Given an array of length D representing the order of the derivative 
  ! for each variable r_1,..,r_D, as well as the ordered corner 
  ! coordinates, This function builds the sub-block of the M matrix 
  ! corresponding to the 2^D equations corresponding to evaluating the 
  ! scheme at the 2^D corners.
  !----------------------------------------------------------------------!
  implicit none

  integer    :: N,D, der_max

  integer    :: corners(2**D,D) 
  integer    :: powers((N+1)**D,D) 
  integer    :: little_m(0:1,0:der_max,0:N)
  integer    :: derivatives(D) 

  integer    :: subM(2**D,(N+1)**D)


  integer    :: k_scalar, p_scalar

  integer    :: dm 

  subM = 1

  do p_scalar = 1, (N+1)**D
   

     do k_scalar = 1, 2**D

        do dm =1, D
                subM(k_scalar,p_scalar)  = subM(k_scalar,p_scalar)  * &
                                 little_m( corners(k_scalar,dm),      &
                                                derivatives(dm),      &
                                            powers(p_scalar,dm))
        end do 

     end do
  end do

end subroutine build_subM

subroutine build_M(N,D,corners,powers,list_derivatives,M)
  !----------------------------------------------------------------------!
  ! Given a list of derivatives array, this subroutine builds the M
  ! matrix. It assumes that the maximum order of derivation is 1,
  ! which is appropriate for the D-cubic schemes.
  !----------------------------------------------------------------------!
  implicit none

  integer,parameter :: der_max = 1

  integer    :: N,D
  integer    :: corners(2**D,D) 
  integer    :: powers((N+1)**D,D) 
  integer    :: list_derivatives(((N+1)/2)**D,D) 
  integer    :: derivatives(D) 

  integer    :: subM(2**D,(N+1)**D)

  integer    :: M((N+1)**D,(N+1)**D)

  integer    :: little_m(0:1,0:der_max,0:N)

  integer    :: i_corner, i_scalar, block
  integer    :: i,  j, offset


  call build_little_m(N,der_max,little_m)

  do block = 1, ((N+1)/2)**D
        derivatives = list_derivatives(block,:) 

        offset = (block-1)*2**D
        call build_subM(N,D,corners,powers,der_max,     &
                        little_m,derivatives,M(offset+1:offset+2**D,:))


  end do

end subroutine build_M


subroutine invert_M(sze,M,Mm1,info_lu,info_inv)
  implicit none

  integer    :: sze
  integer    :: info_lu, info_inv
  integer    :: lwork, ipiv_lapack(sze)

  real(dp)   :: M(sze,sze), Mm1(sze,sze), fake_work(sze)

  real(dp),allocatable  :: work(:)

 
  Mm1 = M

  ! The matrix M is symmetric. Invert it using the right routine

  ! LU decomposition
  call DGETRF(sze,sze,Mm1,sze,ipiv_lapack, info_lu)

  ! find out the best value for lwork
  call DGETRI(sze,Mm1,sze, ipiv_lapack, fake_work, -1,info_lu)

  lwork = int(fake_work(1))

  allocate(work(lwork))
  call DGETRI(sze,Mm1,sze, ipiv_lapack, work, lwork,info_inv)
  deallocate(work)


end subroutine invert_M

  subroutine read_and_write_MA(N,D,sze,what)
  !------------------------------------------------------------------
  ! 
  ! This subroutine reads or writes the matrix MA to xml format. 
  !
  !------------------------------------------------------------------
  use iotk_module
  use intw_utility,    only: find_free_unit
  implicit none

  integer        :: io_unit
  integer        :: record_length, record_index
  integer        :: N,D, sze

  integer        :: read_N,read_D

  character(256) :: filename
  character(49 ) :: attr
  character(1)   :: what


  filename = trim('MA'//trim('.xml'))
  io_unit = find_free_unit()
  
  if ( what .eq. 'w' ) then
          call iotk_open_write (io_unit,filename,binary=.true.)
             call iotk_write_comment(io_unit,  &
                'This file contains the matrix MA relevant to the interpolation scheme.')
             call iotk_write_attr (attr,'D',D)
             call iotk_write_attr (attr,'N',N)
             call iotk_write_empty (io_unit,"dimension_and_order_of_scheme",attr)

             call iotk_write_dat  (io_unit,'MA',MA,columns=sze)
          call iotk_close_write (io_unit)

   else if (what .eq. 'r' ) then
          call iotk_open_read (io_unit,filename,binary=.true.)

             call iotk_scan_empty (io_unit,"dimension_and_order_of_scheme",attr)
             call iotk_scan_attr(attr,"D", read_D)
             call iotk_scan_attr(attr,"N", read_N)

             if (read_D /= D .or. read_N /= N ) then
                write(*,*) '*********************************'
                write(*,*) '**ERROR: read MA matrix has    **'
                write(*,*) '**       inconsistent size     **'
                write(*,*) '**       with input.           **'
                write(*,*) '*********************************'
                stop
             end if
             
             call iotk_scan_dat  (io_unit,'MA',MA)

          call iotk_close_read(io_unit)
   else 
                write(*,*) '*********************************'
                write(*,*) '**ERROR: unknown action in     **'
                write(*,*) '**       read_and_write_MA     **'
                write(*,*) '*********************************'
   end if

  end subroutine read_and_write_MA



subroutine build_delta_offset(D,theta_order,derivatives,delta_offset)
  !----------------------------------------------------------------------!
  ! Given a derivatives array of dimension D with each value either 0 or 1
  ! and with sum of component equal to theta_order(this is not checked) ,
  ! this subroutine returns the offsets corresponding to the centered 
  ! finite difference scheme.
  !----------------------------------------------------------------------!

  implicit none


  integer    :: D
  integer    :: theta_order
  integer    :: derivatives(D) 
  integer    :: delta_offset(2**theta_order,D) 


  integer    :: dm, s
  integer    :: switch, i_min, i_max
  integer    :: J_scalar

  integer    :: p_vector(theta_order)

  ! build the offset arrays.
  i_min  = 0;  i_max = 1
  switch = -1 ! scalar to vector

  do J_scalar = 1, 2**theta_order
     call scalar_to_vector_index(theta_order,p_vector,J_scalar,i_min,i_max,switch)

     s = 1 ! s accounds for the fact that the p_vector does not have dimension D.
     do dm = 1,D
        if (derivatives(dm) == 0) then
            delta_offset(J_scalar,dm) = 0    
        else
            delta_offset(J_scalar,dm) = (-1)**p_vector(s)
            s = s+1
        end if

     end do
  end do

end subroutine build_delta_offset

subroutine build_subA(N,D,corners,theta_order,delta_offset,subA)
  !----------------------------------------------------------------------!
  ! This subroutine builds the subA matrix, which corresponds to the
  ! equation stored in delta_offset.
  !
  !----------------------------------------------------------------------!
  implicit none

  integer    :: N,D

  integer    :: theta_order
  integer    :: corners(2**D,D) 
  integer    :: delta_offset(2**theta_order,D) 
  integer    :: offset(D) 

  real(dp)   :: subA(2**D,(N+1)**D)

  real(dp)   :: a, frac, phase


  integer    :: k_scalar, l_scalar, j_scalar
  integer    :: k_vector(D), l_vector(D)

  integer    :: i_min, i_max, switch
  integer    :: dm

  subA = ZERO

  if     (N == 3 ) then 
          i_min  = -1
          i_max  = 2
  else if (N == 1) then
          i_min  = 0 
          i_max  = 1
  end if

  switch = 1 ! vector to scalar

  frac = 1.0_dp/(2**theta_order)

  do j_scalar = 1, 2**theta_order

     offset  =  delta_offset(j_scalar,:)

     phase = 1
     do dm = 1, D 
        if( offset(dm) /= 0) phase = phase*offset(dm)
     end do 
     a = frac*phase

     do k_scalar = 1, 2**D

        k_vector = corners(k_scalar,:)
     
        l_vector = k_vector+offset
        call scalar_to_vector_index(D,l_vector,l_scalar,i_min,i_max,switch)
        subA(k_scalar,l_scalar) = subA(k_scalar,l_scalar) + a 
     end do

  end do

end subroutine build_subA

subroutine build_A(N,D,corners,list_derivatives,A)
  !----------------------------------------------------------------------!
  ! This subroutine builds the A matrix from its subcomponents subA 
  !----------------------------------------------------------------------!
  implicit none

  integer    :: N,D
  integer    :: corners(2**D,D) 
  integer    :: list_derivatives(((N+1)/2)**D,D) 
  integer    :: derivatives(D) 

  integer    :: delta_offset(2**D,D) 

  integer    :: theta_order

  real(dp)   :: A((N+1)**D,(N+1)**D)

  integer    :: block, offset


  do block = 1, ((N+1)/2)**D
        derivatives = list_derivatives(block,:) 
        theta_order = sum(derivatives)
        
        call build_delta_offset(D,theta_order,derivatives, &
                            delta_offset(1:2**theta_order,:))

        offset = (block-1)*2**D
        
        call build_subA(N,D,corners,theta_order,   &
                delta_offset(1:2**theta_order,:),A(offset+1:offset+2**D,:))
  end do

end subroutine build_A

subroutine build_mu_power_array(N,D,powers,t,mu)
  !----------------------------------------------------------------------!
  ! Given a a vector of variables t(D) = (t1,...,tD), with 0 <= td <= 1, 
  ! This routine returns the array mu(p_scalar), with 
  !  mu(p_scalar) = t1^p1 * t2^p2 *...* tD^pD,   where (p1,..,pD) 
  !                                              is the vector array 
  !                                             corresponding to p_scalar.
  !----------------------------------------------------------------------!
  implicit none

  ! input variables
  integer    :: N,D
  integer    :: powers((N+1)**D,D) 

  real(dp)   :: t(D)
  real(dp)   :: mu((N+1)**D)

  ! local variables
  integer    :: p_scalar
  integer    :: dm, p
  real(dp)   :: t_powers(0:N,D)

  t_powers(0,:) = ONE
  ! compute the powers
  do p = 1,N
     do dm = 1,D
        t_powers(p,dm) = t_powers(p-1,dm)*t(dm)
     end do
  end do
 
  ! compute the mu array 
  mu = ONE

  do p_scalar = 1,(N+1)**D

     do dm = 1, D

        p  = powers(p_scalar,dm) 
        mu(p_scalar) = mu(p_scalar) * t_powers(p,dm)

     end do 

  end do

end subroutine build_mu_power_array

subroutine find_mesh_offset(N,D,dX)
  !----------------------------------------------------------------------!
  ! Given the dimension of space (D) and the order of the scheme (N),
  ! This routine returns the integer array
  !             dX((N+1)**D,D) 
  ! where each row corresponds to the offset from the X0 corner corresponding
  ! to l_scalar .
  !----------------------------------------------------------------------!
  implicit none
  integer  :: N, D
  integer  :: dX((N+1)**D,D)

  integer  :: i_min, i_max, switch
  integer  :: l_scalar
  integer  :: l_vector(D)

  if     (N == 3 ) then
          i_min  = -1
          i_max  = 2
  else if (N == 1) then
          i_min  = 0
          i_max  = 1
  end if

  switch = -1 ! scalar to vector 
  do l_scalar =1, (N+1)**D
      call scalar_to_vector_index(D,l_vector,l_scalar,i_min,i_max,switch)
      dX(l_scalar,:) = l_vector
  end do

end subroutine find_mesh_offset


subroutine find_determinant(sze,M,determinant,info_lu)
  implicit none

  integer    :: sze
  integer    :: info_lu, i
  integer    :: ipiv_lapack(sze)

  real(dp)   :: M(sze,sze),local_M(sze,sze)
  real(dp)   :: determinant
  
  local_M = M
  ! LU decomposition
  call DGETRF(sze,sze,local_M,sze,ipiv_lapack, info_lu)
 
  determinant = ONE
  do i = 1, sze
     determinant = determinant * local_M(i,i)
  end do 

end subroutine find_determinant


subroutine build_MA(N,D)
  !----------------------------------------------------------------------!
  ! This is a driver subroutine which performs all the steps 
  ! necessary for the construction of the MA matrix, necessary
  ! for neighbor-based interpolation.
  !----------------------------------------------------------------------!
  implicit none 

  integer  ::  N, D
  integer  ::  sze, sub_sze
  integer  ::  info_lu, info_inv

  integer,allocatable  ::      corners(:,:), powers(:,:), derivatives(:,:)

  integer,allocatable  ::      M(:,:)
  real(dp),allocatable ::      A(:,:), Mreal(:,:), Mm1(:,:)


  !---------------------------

  sze     = (N+1)**D
  sub_sze = ((N+1)/2)**D
  
  ! allocate various arrays. Allocating inside this subroutine
  ! allows passing N_scheme and D_scheme, the more fundamental
  ! variables, instead of sze.

  allocate(corners(2**D,D))
  allocate(powers(sze,D))
  allocate(derivatives(sub_sze,D))
  allocate(M(sze,sze))
  allocate(A(sze,sze))
  allocate(Mreal(sze,sze))
  allocate(Mm1(sze,sze))


  call build_corners(D, corners)
  call build_powers(N,D,powers)
  call build_derivatives(N,D,derivatives)


  call build_M(N,D,corners,powers,derivatives,M)

  Mreal(:,:)  = dble(M)

  call invert_M(sze,Mreal,Mm1,info_lu,info_inv)

  if ( info_lu /= 0 .or. info_inv /= 0 ) then
     write(*,*) '*******************************************'
     write(*,*) '*  PROBLEM GENERATING MATRIX MA:          *'
     write(*,*) '*  REVIEW CODE, SUBROUTINE build_MA       *'
     write(*,*) '*            program stops.               *'
     write(*,*) '*******************************************'
     stop
  end if

  call build_A(N,D,corners,derivatives,A)
  
  ! get the product MA = M^-1 * A

  call DGEMM ( 'N', 'N', sze, sze, sze, ONE, Mm1,      &
                        sze, A, sze, ZERO, MA, sze)

  deallocate(corners)
  deallocate(powers)
  deallocate(derivatives)
  deallocate(M)
  deallocate(A)
  deallocate(Mreal)
  deallocate(Mm1)

end subroutine build_MA

!================================================================================
!  Older subroutines, still useful for testing
!================================================================================

subroutine old_build_subM(N,D,corners,powers,derivatives,subM)
  !----------------------------------------------------------------------!
  ! Given an array of length D representing the order of the derivative 
  ! for each variable r_1,..,r_D, as well as the ordered corner 
  ! coordinates, This function builds the sub-block of the M matrix 
  ! corresponding to the 2^D equations corresponding to evaluating the 
  ! scheme at the 2^D corners.
  !----------------------------------------------------------------------!
  implicit none

  integer    :: N,D

  integer    :: corners(2**D,D) 
  integer    :: derivatives(D) 
  integer    :: powers((N+1)**D,D) 

  integer    :: subM(2**D,(N+1)**D)

  integer    :: i_corner, i_scalar

  integer    :: der, derivative 
  integer    :: dm, coeff, power, value, K
  integer    :: differentiate

  do i_corner =1, 2**D
     do i_scalar = 1, (N+1)**D

        value = 1
        
        do dm = 1,D         

                derivative   = derivatives(dm)
                power        = powers(i_scalar,dm) 
                K            = corners(i_corner,dm) 

                value = value*differentiate(K,power,derivative)

        end do

        subM(i_corner,i_scalar)  =  value

     end do 
  end do 

end subroutine old_build_subM


!--------------------------------------------------------------------------------
!
end module intw_new_interpolate
!
!--------------------------------------------------------------------------------

integer function differentiate(x,power,derivative_order)
        implicit none
        INTEGER, INTENT(IN) :: x,power,derivative_order
        INTEGER             :: coeff,term
      
        coeff = 1 

        if ( derivative_order > power) then
                 differentiate = 0
        else
                do term = power-derivative_order+1, power
                        coeff = coeff*term
                end do
                differentiate = coeff*x**(power-derivative_order)
        end if

end function differentiate

integer function factorial(n)
        implicit none
        integer, intent(in) :: n
        integer             :: i
    
        factorial = 1 
        do i = 1,n
                factorial = factorial*i
        end do

end function factorial
