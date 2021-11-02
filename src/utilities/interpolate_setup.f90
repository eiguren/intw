!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       program interpolate_setup 
!       -----------------------------------
!
!       This is a "utility" program which is part of the intw project.
!
!       The purpose of this utility is to generate the matrix relevant 
!       to a given interpolation scheme.
!
!       The interpolation formula has the following form:
!
!               f(x) = sum_i C[...]  x1^i1 ... xD^iD
!       where D is the dimension of space, and 0 <= i <= N, where N
!       is the order of the scheme.
!
!       The coefficients C[...] are related to the exact values of f
!       on a mesh of points according to 
!
!                             M * C = F = A * f
!       where the exact definitions are provided elsewhere. 
!
!       This program computes the matrix M^{-1}*A, which allows
!       the direct computation of C through:
!                       C = [ M^{-1}*A ] f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program interpolate_setup 


  use intw_utility
  use intw_new_interpolate


!================================================================================
!       Declare the variables 
!================================================================================
implicit none

integer              ::      N, D 
integer              ::      i, j, non_zero
integer              ::      sze , sub_sze
integer              ::      info_lu, info_inv

integer,allocatable  ::      corners(:,:), powers(:,:)
integer,allocatable  ::      M(:,:)



real(dp),allocatable::      Mreal(:,:), Mm1(:,:), A(:,:),  prod(:,:)
real(dp),allocatable::      MA_matmul(:,:)
real(dp),allocatable::      MA_test(:,:)

real(dp)            ::      dMA

integer,allocatable  ::      derivatives(:,:)

real(dp)            ::      total_time1, total_time2, time1, time2


call get_timing(total_time1)
!================================================================================
!       Talk to the user 
!================================================================================

write(*,*) '====================================================='
write(*,*) '|              program interpolate_setup            |'
write(*,*) '|        ---------------------------------          |'
write(*,*) '====================================================='
write(*,*) '|                                                   |'
write(*,*) '|        This program computes the matrix           |'
write(*,*) '|                     M^{-1} * A                    |'
write(*,*) '|        relevant to the interpolation scheme.      |'
write(*,*) '|        The matrix is then saved in binary         |'
write(*,*) '|        format for easy use by the main program.   |'
write(*,*) '====================================================='


!================================================================================
!       Ask the user what order and dimension he needs
!================================================================================

write(*,*) '|        what is the dimension of space?            |'
read(*,*) D
write(*,*) '|        what is the order of the scheme?           |'
read(*,*) N

write(*,'(a,I3,a,I3,a)') ' |        the user has chosen: D =',D,', N =',N,'        |'

if (N /= 1 .and. N /=3 ) then

                write(*,*) '************************************'
                write(*,*) '**ERROR: only N = 1 (linear )    ***'
                write(*,*) '**       or   N = 3 (cubic )     ***'
                write(*,*) '**       schemes are implemented ***'
                write(*,*) '************************************'
        
                stop
end if

!================================================================================
!       allocate the relevant arrays 
!================================================================================
sze     = (N+1)**D
sub_sze = ((N+1)/2)**D

allocate(corners(2**D,D))
allocate(powers(sze,D))
allocate(M(sze,sze))
allocate(A(sze,sze))
call allocate_MA(sze)
allocate(MA_matmul(sze,sze))
allocate(MA_test(sze,sze))
allocate(Mm1(sze,sze))
allocate(Mreal(sze,sze))
allocate(prod(sze,sze))
allocate(derivatives(sub_sze,D))

!================================================================================
!       build MA using the subroutine
!================================================================================

call build_MA(N,D)
MA_test = MA
MA      = ZERO

!================================================================================
!       build the necessary arrays 
!================================================================================


call build_corners(D, corners)
call build_powers(N,D,powers)
call build_derivatives(N,D,derivatives)

!================================================================================
!       Perform computation of the M matrix 
!================================================================================

write(*,*) '|        - finding the M matrix ...                 |'
call get_timing(time1)
call build_M(N,D,corners,powers,derivatives,M) 
call get_timing(time2)

write(*,'(a,F6.2,a)') ' |                  time :: ',time2-time1,  &
                                  ' seconds           |'



write(*,*) '|        - finding the M^{-1} matrix ...            |'

Mreal(:,:)  = dble(M)
call get_timing(time1)
call invert_M(sze,Mreal,Mm1,info_lu,info_inv)
call get_timing(time2)

write(*,'(a,F6.2,a)') ' |                  time :: ',time2-time1,  &
                                  ' seconds           |'


write(*,*) '|        - finding the A matrix...                  |'
call get_timing(time1)
call build_A(N,D,corners,derivatives,A)
call get_timing(time2)

write(*,'(a,F6.2,a)') ' |                  time :: ',time2-time1,  &
                                  ' seconds           |'


write(*,*) '|        - finding the MA = M^{-1}*A matrix...      |'
call get_timing(time1)
call DGEMM ( 'N', 'N', sze, sze, sze, ONE, Mm1,      &
                        sze, A, sze, ZERO, MA, sze)
call get_timing(time2)
write(*,'(a,F6.2,a)') ' |          DGEMM   time :: ',time2-time1,  &
                                  ' seconds           |'

call get_timing(time1)
MA_matmul = matmul(Mm1,A)
call get_timing(time2)
write(*,'(a,F6.2,a)') ' |          matmul  time :: ',time2-time1,  &
                                  ' seconds           |'


write(*,*) '====================================================='
write(*,*) '==========             checks                ========'
write(*,*) '====================================================='
write(*,'(a,I4,a)') ' |       - info(LU decomposition) =  ',info_lu,'            |'
write(*,'(a,I4,a)') ' |       - info( inversion      ) =  ',info_inv,'            |'

dMA = sqrt(sum(MA(:,:)-MA_matmul(:,:))**2)
write(*,'(a,E12.4,a)') ' |       - |MA - MA_matmul |      =  ',dMA,'    |'

dMA = sqrt(sum(MA(:,:)-MA_test(:,:))**2)
write(*,'(a,E12.4,a)') ' |       - |MA - MA_test   |      =  ',dMA,'    |'


! test inverse is full of integers

do i = 1, sze
   do j = 1, sze
        if ( abs(Mm1(i,j) - nint(Mm1(i,j))) > 1e-8)  then

                write(*,*) '***********************************'
                write(*,*) '**ERROR: the inverse matrix is  ***'
                write(*,*) '**       not full of integers!  ***'
                write(*,*) '***********************************'
                stop
        end if
   end do
end do

write(*,*) '|       - The inverse matrix is composed of integers|' 


write(*,*) '|       - finding the product M * M^{-1} ...        |'
call get_timing(time1)
call DGEMM ( 'N', 'N', sze, sze, sze, ONE, Mreal,       &
                        sze, Mm1, sze, ZERO, prod, sze)
call get_timing(time2)
write(*,'(a,F6.2,a)') ' |                  time :: ',time2-time1,  &
                                  ' seconds           |'

! check that the product of M and M^{-1} is the identity
do i = 1, sze
   do j = 1, sze

        if ( (i == j .and. abs(prod(i,j) - ONE) > 1e-8 )        & 
                           .or.                                 &                                                        
             (i /= j .and.     abs(prod(i,j))   > 1e-8 )        &
                                                        )  then

                write(*,*) '***********************************'
                write(*,*) '**ERROR: the product M * M^{-1} ***'
                write(*,*) '**       is not the identity !  ***'
                write(*,*) '***********************************'
                stop
        end if
   end do
end do


write(*,*) '|       -  M*M^{-1} = Identity                      |' 
write(*,*) '====================================================='
write(*,*) '==========             writing...            ========'
write(*,*) '====================================================='

call read_and_write_MA(N,D,sze,'w')

!================================================================================
!       clean up 
!================================================================================

deallocate(corners)
deallocate(powers)
deallocate(M)
call deallocate_MA()
deallocate(A)
deallocate(Mm1)
deallocate(prod)
deallocate(derivatives)

call get_timing(total_time2)
write(*,'(a,F10.3,a)') ' |     DONE - total time::',total_time2-total_time1,' seconds         |' 
write(*,*) '====================================================='

end program interpolate_setup 
