!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       program interpolator
!       -----------------------------------
!
!       This is a "utility" program which is part of the intw project.
!
!       The purpose of this utility is to test the interpolation scheme. 
!
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program interpolator


  use intw_utility
  use intw_useful_constants
  use intw_new_interpolate


!================================================================================
!       Declare the variables 
!================================================================================
implicit none


integer              ::      N, D, sze
integer              ::      i, j, offset, l_scalar

integer              ::      i_c, i_mu, i_cmu

real(dp),allocatable ::      t(:), x(:), mu(:), f(:), c(:)
real(dp),allocatable ::      x_corner_1(:),x_corner_2(:)
real(dp),allocatable ::      list_x(:,:), dx_l(:,:)
integer,allocatable  ::      powers(:,:)


real(dp)            ::      exact_f, int_f 
real(dp)            ::      test_function
real(dp)            ::      time1, time2
real(dp)            ::      time_c1, time_c2, total_time_c
real(dp)            ::      time_mu1, time_mu2, total_time_mu
real(dp)            ::      time_cmu1, time_cmu2, total_time_cmu
real(dp)            ::      DDOT

integer              ::      io_unit
!================================================================================
!       ask the user for input 
!================================================================================
write(*,*) 'what is the dimension of space?'
read(*,*) D
write(*,*) 'what is the order of the scheme?'
read(*,*) N

!================================================================================
!       allocate the relevant arrays 
!================================================================================
sze = (N+1)**D

call allocate_MA(sze)
allocate(powers(sze,D))
allocate(dx_l(sze,D))
allocate(t(D))
allocate(x(D))
allocate(x_corner_1(D))
allocate(x_corner_2(D))
allocate(mu(sze))
allocate(f(sze))
allocate(c(sze))


call build_powers(N,D,powers)

call find_dx_l(N,D,dx_l)


!================================================================================
!       read in the matrix 
!================================================================================
write(*,*) 'reading the MA matrix...'
call get_timing(time1)
call read_and_write_MA(N,D,sze,'r')
call get_timing(time2)
write(*,'(a,F8.4,a)') 'it took ',time2-time1,' seconds'

! TEST
!i_mu = 0
!do i=1, sze
!   do j=1, sze
!        if (abs(MA(i,j)) > eps_8) i_mu = i_mu + 1
!   end do
!end do
!
!write(*,*) 'there is a total of ',sze*sze,' elements in MA'
!write(*,*) 'there are ',i_mu,' non-zero elements in MA'
!write(*,*) 'Non-zero fraction: ',dble(i_mu)/(sze*sze)
!stop
!================================================================================
!       interpolate 
!================================================================================

io_unit =find_free_unit()
open(unit=io_unit,file='output.dat',status='unknown')
write(io_unit,'(a)') '#   x(1)              exact f           interpolated f'

x_corner_1 = ZERO

i_c  = 0
i_mu = 0
i_cmu= 0
total_time_c  = ZERO
total_time_mu = ZERO
total_time_cmu= ZERO

write(*,*) 'interpolating...'
call get_timing(time1)
! find fl
do l_scalar =1 ,sze
   f(l_scalar) = test_function(D,x_corner_1+dx_l(l_scalar,:))
end do
! find c
call DGEMV ( 'N', sze, sze, ONE, MA, sze, f, 1, ZERO, c, 1)

do i=1,1000

        x = dble(i)/dble(100)

       
        call find_corner(D,x,x_corner_2,t)

        if ( sum(x_corner_2-x_corner_1)**2 > 1e-8 ) then
                x_corner_1 = x_corner_2

                ! find fl
                do l_scalar =1 ,sze
                        f(l_scalar) = test_function(D,x_corner_1+dx_l(l_scalar,:))
                end do
                ! find c
                call get_timing(time_c1)
                call DGEMV ( 'N', sze, sze, ONE, MA, sze, f, 1, ZERO, c, 1)
                call get_timing(time_c2)
                i_c = i_c + 1
                total_time_c = total_time_c +time_c2-time_c1 
        end if
        ! find mu

        call get_timing(time_mu1)
        call build_mu_power_array(N,D,powers,t,mu)
        call get_timing(time_mu2)
        i_mu = i_mu + 1
        total_time_mu = total_time_mu +time_mu2-time_mu1 

        call get_timing(time_cmu1)
        int_f   = DDOT(sze,c,1,mu,1)
        call get_timing(time_cmu2)
        i_cmu = i_cmu + 1
        total_time_cmu = total_time_cmu +time_cmu2-time_cmu1 

        exact_f = test_function(D,x)

        write(io_unit,'(3F18.12)') x(1), exact_f, int_f 
end do

call get_timing(time2)
close(unit=io_unit)

write(*,'(a,F8.4,a)') 'it took ',time2-time1,' seconds'

write(*,'(a)') '------------   STATS  --------------'
write(*,'(a,F8.4,a,e12.4,a)') 'Computing the c array :: ',total_time_c,' seconds (',total_time_c/i_c,' seconds per call)'
write(*,'(a,F8.4,a,e12.4,a)') 'Computing the mu array:: ',total_time_mu,' seconds (',total_time_mu/i_mu,' seconds per call)'
write(*,'(a,F8.4,a,e12.4,a)') 'Computing c*mu        :: ',total_time_cmu,' seconds (',total_time_cmu/i_cmu,' seconds per call)'

!do offset = 0,5


!   t    = (dble(offset)-1.0)
!   f(1) = test_function(D,t)
!
!   t    = (dble(offset))
!   f(2) = test_function(D,t)
!
!   t    = (dble(offset)+1.0)
!   f(3) = test_function(D,t)
!
!   t    = (dble(offset)+2.0)
!   f(4) = test_function(D,t)
!
!   ! find c
!   call DGEMV ( 'N', sze, sze, ONE, MA, sze, f, 1, ZERO, c, 1)
!
!        do i = 0, 100
!
!           
!           t = (dble(offset)+dble(i)/dble(100))
!           exact_f = test_function(D,t)
!
!           t = (dble(i)/dble(100))
!
!           call build_mu_power_array(N,D,powers,t,mu)
!           int_f = DDOT(sze,c,1,mu,1)
!
!           write(111,'(3F18.12)') dble(offset)+t(1), exact_f, int_f 
!        end do
!end do

!================================================================================
!       clean up 
!================================================================================
call deallocate_MA()
deallocate(powers)
deallocate(t)
deallocate(mu)
deallocate(x)
deallocate(dx_l)
deallocate(x_corner_1)
deallocate(x_corner_2)
deallocate(f)
deallocate(c)

end program interpolator

!================================================================================
!       some driving subroutines
!================================================================================

subroutine find_corner(D,x,x_corner,t)
  !----------------------------------------------------------------------!
  ! Given a a vector of variables x(D) = (x1,...,xD)
  ! This routine returns the arrays
  !      x_corner(D) = (X1,..,XD) 
  !      t(D)        = (t1,..,tD) 
  ! where Xi is the nearest integer below xi and  ti the remainder
  !----------------------------------------------------------------------!
  use intw_useful_constants
  implicit none
  integer  :: D 
  real(dp) :: x(D),t(D),x_corner(D)

  integer  :: dm

  x_corner = floor(x)
  t        = x - x_corner

end subroutine find_corner

subroutine find_dx_l(N,D,dx_l)
  !----------------------------------------------------------------------!
  ! Given the dimension of space (D) and the order of the scheme (N),
  ! This routine returns the array
  !             dx_l((N+1)**D,D) 
  ! where each row corresponds to the offset from the X0 corner corresponding
  ! to l_scalar .
  !----------------------------------------------------------------------!
  use intw_useful_constants
  use intw_new_interpolate
  implicit none
  integer  :: N, D 
  real(dp) :: dx_l((N+1)**D,D)

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
      dx_l(l_scalar,:) = dble(l_vector)
  end do

end subroutine find_dx_l

!================================================================================
!       some test interpolation function 
!================================================================================
real(dp) function test_function(D,t)
        use intw_useful_constants
        implicit none
        integer  :: D , dm
        real(dp) :: t(D) 

        test_function  = ONE
        do dm = 1,D
           test_function = test_function + t(dm)*cos(t(dm)) 
        end do

end function test_function


