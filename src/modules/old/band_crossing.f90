!       intw project.
!
!----------------------------------------------------------------------------!
module intw_band_crossing
!----------------------------------------------------------------------------!
!
!       This module contains subroutines which help to disentangle band 
!       crossing and anti-crossing information, using  Wannier interpolation
!       data.
!
!----------------------------------------------------------------------------!

use intw_useful_constants

use intw_utility
use intw_input_parameters
use intw_W90

  !
  implicit none
  !
  save
 
  character(1),          allocatable :: interpolation_labels(:)
  real(kind=dp),         allocatable :: interpolation_special_points(:,:)


contains

  subroutine find_degenerate_eigenvalue_sets(eig,degenerate_sets,number_of_sets)
  !----------------------------------------------------------------------------!
  !   Given a list of eigenvalues, this subroutine finds the degenerate subsets
  !   and the band indices corresponding to each subset.
  !   The eigenvalues are NOT assumed to be ordered.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  real(dp)       :: eig(num_wann)

  ! output variables
  integer        :: degenerate_sets(num_wann,num_wann)
  integer        :: number_of_sets

  ! local variables
  integer        :: n_bnd, m_bnd

  integer        :: i_set, n_subset
  
  logical        :: found(num_wann)

  real(dp)       :: tol

  tol = eps_8
!  tol = ONE

  ! initialize
  degenerate_sets = 0

  n_subset    = 0
  i_set       = 0


  n_bnd       = 1
  found       = .false.

  ! tabulate the degeneracies
  do
     ! loop exit condition
     if (n_bnd == num_wann+1) exit

     ! n_bnd is found
     found(n_bnd) = .true.

     ! initialize a new set!
     i_set     = i_set+1 
     n_subset  = 1
  
     degenerate_sets(i_set,n_subset) = n_bnd

     ! fish out all remaining terms which belong to this set
     do m_bnd = n_bnd+1,num_wann
        if (.not. found(m_bnd) .and. abs(eig(n_bnd)-eig(m_bnd)) < tol ) then
             found(m_bnd) = .true.
             n_subset     = n_subset + 1
             degenerate_sets(i_set,n_subset) = m_bnd
        end if
     end do

     do while ( n_bnd <= num_wann)
        if  ( found(n_bnd) ) then
                n_bnd = n_bnd + 1
        else
                exit
        end if
     end do

!     CAREFUL! n_bnd can be out of bounds below! Bug found by Ion!
!     do while ( found(n_bnd) .and. n_bnd <= num_wann)
!        n_bnd = n_bnd + 1
!     end do

  end do 

  number_of_sets = i_set

  end subroutine find_degenerate_eigenvalue_sets


  subroutine find_n_subset_max(degenerate_sets_k,i_set,n_subset_max)
  !------------------------------------------------------------------
  ! This subroutine simply finds the number of elements in a set
  !------------------------------------------------------------------

  integer       :: degenerate_sets_k(num_wann,num_wann)
  integer       :: n_subset_max
  integer       :: i_set

  n_subset_max = 1
  do
     if (degenerate_sets_k(i_set,n_subset_max) == 0 &
                     .or. n_subset_max == num_wann) exit

      n_subset_max = n_subset_max + 1
   end do
   n_subset_max = n_subset_max - 1

  end subroutine find_n_subset_max

  subroutine find_number_of_sets(degenerate_sets_k,number_of_sets)
  !------------------------------------------------------------------
  ! This subroutine simply finds the number sets stored in 
  ! degenerate_sets_k.
  !------------------------------------------------------------------

  integer       :: degenerate_sets_k(num_wann,num_wann)
  integer       :: number_of_sets

  integer       :: i_set

  number_of_sets = num_wann
  do
     if (degenerate_sets_k(number_of_sets,1) /= 0 ) exit
     number_of_sets  = number_of_sets - 1
   end do

  end subroutine find_number_of_sets



  subroutine find_eigenvalue_degeneracies(eig_array,degeneracy_array, &
                                                number_of_degeneracies)
  !----------------------------------------------------------------------------!
  !   Given a list of eigenvalues, this subroutine establishes if there
  !   are degeneracies. It assumes that the eigenvalues are sorted.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  real(dp)       :: eig_array(num_wann)
  integer        :: degeneracy_array(num_wann)
  integer        :: number_of_degeneracies

  ! local variables
  integer        :: n_bnd
  integer        :: degenerate_index
  real(dp)       :: eigenvalue

  number_of_degeneracies = 0
  
  degeneracy_array       = 0
  
  degenerate_index       = 1

  eigenvalue             = eig_array(1)
  degeneracy_array(1)    = degenerate_index

  do n_bnd = 2,num_wann
        if (abs(eigenvalue-eig_array(n_bnd)) < eps_8) then
                number_of_degeneracies  = number_of_degeneracies + 1
        else
                degenerate_index        = degenerate_index  +  1
                eigenvalue              = eig_array(n_bnd)
        end if
        degeneracy_array(n_bnd) = degenerate_index
  end do

  end subroutine find_eigenvalue_degeneracies



  subroutine create_permutation_matrix(P_k,Permutations)
  !----------------------------------------------------------------------------!
  ! This subroutine creates a permutation matrix for the eigenvalues, 
  ! based on the absolute values of the elements in the matrix of square 
  ! coefficients.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  real(dp)       :: Permutations(num_wann,num_wann),P_k(num_wann,num_wann)

  ! local variables
  integer        :: n_bnd, m_bnd, max_bnd
  integer        :: n_bnd_max, m_bnd_max
  real(dp)       :: p, threshold_A
  logical        :: allowed_column(num_wann), allowed_row(num_wann)

  integer        :: left_over

  allowed_column = .true.
  allowed_row    = .true.
  Permutations   = 0.0_dp

  ! step A: a first pass to find large elements
  threshold_A = 0.7

  left_over   = num_wann

  do m_bnd = 1, num_wann
     do n_bnd = 1, num_wann
        if ( P_k(m_bnd,n_bnd) > threshold_A .and. allowed_column(n_bnd) ) then
                Permutations(m_bnd,n_bnd)   = 1.0_dp
                allowed_column(n_bnd)       = .false.
                allowed_row   (m_bnd)       = .false.
                left_over = left_over - 1
                exit
        end if
     end do
  end do
 
  
  ! step B: iterate on what is left over

  do 
     if (left_over == 0) exit
        write(*,*) 'step B'

     ! find the first acceptable point
     do m_bnd = 1, num_wann
        if( allowed_row(m_bnd)) then
                m_bnd_max = m_bnd
                exit
        end if
     end do

     do n_bnd = 1, num_wann
        if( allowed_column(n_bnd)) then
                n_bnd_max = n_bnd
                exit
        end if
     end do

     p  =  P_k(m_bnd_max,n_bnd_max)

     do m_bnd = 1, num_wann
        if( allowed_row(m_bnd)) then
            do n_bnd = 1, num_wann
               if( allowed_column(n_bnd) ) then
                        if ( P_k(m_bnd,n_bnd) >= p) then
                             n_bnd_max = n_bnd
                             m_bnd_max = m_bnd
                             p         = P_k(m_bnd,n_bnd)
                        end if
               end if
            end do

        end if 
     end do

     Permutations(m_bnd_max,n_bnd_max) = 1.0_dp
     allowed_row   (m_bnd_max)  =  .false.
     allowed_column(n_bnd_max)  =  .false.
 
     left_over = left_over - 1 
  end do

  end subroutine create_permutation_matrix

  subroutine check_permutation_matrix_old(Permutations,test_passed)
  !----------------------------------------------------------------------------!
  ! This subroutine checks that the matrix provided contains exactly one 1
  ! in every row and every column.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  real(dp)       :: Permutations(num_wann,num_wann)
  logical        :: test_passed
  ! local variables

  integer        :: columns(num_wann), rows(num_wann) 
  integer        :: n_bnd, m_bnd
  real(dp)       :: p


  test_passed = .true.

  columns = 0
  rows    = 0

  do m_bnd = 1, num_wann
     do n_bnd = 1, num_wann
        if (abs(Permutations(m_bnd,n_bnd)-1.0_dp) < eps_8) then 
               columns(n_bnd)  = columns(n_bnd) +1
               rows   (m_bnd)  = rows   (m_bnd) +1
        else if (abs(Permutations(m_bnd,n_bnd)) > eps_8) then 
               write(*,*) ' Elements not 1 or 0, m =',m_bnd,' n = ',n_bnd
               test_passed = .false.
        end if
     end do
  end do
 
  do m_bnd = 1, num_wann
       if ( columns(m_bnd) /= 1 ) then
               write(*,*) ' column #',m_bnd,' doesnt have a single 1'
               test_passed = .false.
       else if ( rows(m_bnd) /= 1 ) then
               test_passed = .false.
               write(*,*) ' row    #',m_bnd,' doesnt have a single 1'
       end if
  end do 

  end subroutine check_permutation_matrix_old


  subroutine get_interpolation_special_points(number_of_special_points)
  !----------------------------------------------------------------------------!
  !   This subroutine reads in the input file  and reads the special k points
  !----------------------------------------------------------------------------!

    implicit none

    integer           :: stdin_unit
    integer           :: ierr

    integer           :: number_of_special_points, ipoint

    logical           :: found_card
    character(len=maxlen) :: test_line, line


    stdin_unit  =   5

    rewind(stdin_unit)  

    found_card = .false.

    do
        read(5,'(A)',iostat=ierr) line

        test_line = adjustl(line)
        if (test_line(1:11) == 'Kpoint_Path') then
                found_card = .true.                     
                exit 
        end if

        if(  ierr /= 0) exit 
    end do      


    if(found_card) then

        read(stdin_unit,*,iostat=ierr) number_of_special_points 

        allocate(interpolation_labels(number_of_special_points))
        allocate(interpolation_special_points(3,number_of_special_points))

        do ipoint = 1,number_of_special_points
                read(stdin_unit,*,iostat=ierr)                            &
                                        interpolation_labels(ipoint),                &
                                        interpolation_special_points(1,ipoint),      &
                                        interpolation_special_points(2,ipoint),      &
                                        interpolation_special_points(3,ipoint)
                                        
        end do
    end if

  end subroutine get_interpolation_special_points

  subroutine build_list_kpath(number_of_special_points,nk_vec_path,    &
                                list_k_path, list_x_path, list_xticks)
  !----------------------------------------------------------------------------!
  !   This subroutine creates a list of k vectors along the path specified
  !   in the input file. The array list_k_path will contain the k-vectors
  !   in CRYSTAL coordinates.
  !
  !   The subroutine also computes a distance "x" from one k-point to another.
  !   To compute this quantity, CARTESIAN coordinates will be used, since
  !   this coordinate system preserves distances.
  !----------------------------------------------------------------------------!

    implicit none

    integer           :: number_of_special_points,nk_vec_path
    integer           :: iloop, jloop

    integer           :: branch_points, total_branch_points
    real(dp)          :: total_path_length
    

    real(dp)          :: list_xticks(number_of_special_points)
    real(dp)          :: list_x_path(nk_vec_path)
    real(dp)          :: list_k_path(3,nk_vec_path)

    real(dp)          :: k_cart(3), dk_cart(3), dk(3)
    real(dp)          :: norm_dk, norm_dk_cart
    real(dp)          :: x,  dx
    integer           :: ix

    total_path_length   = 0.0_dp
    total_branch_points = 0

    do iloop = 2, number_of_special_points

        dk = interpolation_special_points(:,iloop)   &
            -interpolation_special_points(:,iloop-1)

        norm_dk = sqrt(dk(1)**2+dk(2)**2+dk(3)**2)

        total_path_length  = total_path_length  + norm_dk  
                                
    end do

    ix = 1
    list_x_path(ix)   = 0.0_dp    

    list_k_path(:,1)  = interpolation_special_points(:,1)


    list_xticks(1)  = 0.0_dp    

    do iloop = 2, number_of_special_points-1

        dk(:)   = interpolation_special_points(:,iloop)   &
                 -interpolation_special_points(:,iloop-1)


        norm_dk    = sqrt(dot_product(dk,dk))

        branch_points =  int(dble(nk_vec_path-1)*      &
                         norm_dk/dble(total_path_length))

        total_branch_points = total_branch_points + branch_points

        dk = dk/dble(branch_points)
        norm_dk    = sqrt(dot_product(dk,dk))

        dk_cart(:)   = matmul(bg,dk)
        norm_dk_cart = sqrt(dot_product(dk_cart,dk_cart))

        do jloop=1,branch_points
                ix = ix+1
                list_k_path(:,ix) = list_k_path(:,ix-1)+dk
                list_x_path(ix)   = list_x_path(ix-1)  +norm_dk_cart
        end do
        list_xticks(iloop)  = list_x_path(ix)
    end do

    branch_points  = nk_vec_path-1-total_branch_points 

    dk   = interpolation_special_points(:,number_of_special_points)   &
           -interpolation_special_points(:,number_of_special_points-1)

    dk      = dk/dble(branch_points)
    norm_dk = sqrt(dot_product(dk,dk))

    dk_cart(:)   = matmul(bg,dk)
    norm_dk_cart = sqrt(dot_product(dk_cart,dk_cart))

    do jloop=1,branch_points
       ix = ix+1
       list_k_path(:,ix) = list_k_path(:,ix-1)+dk
       list_x_path(ix)   = list_x_path(ix-1)  +norm_dk_cart
    end do
    list_xticks(number_of_special_points)  = list_x_path(ix)


  end subroutine build_list_kpath


  subroutine diagonalize_perturbation(u_dagger,dh,d_eig,degenerate_sets,number_of_sets)
  !----------------------------------------------------------------------------!
  !   Given a perturbation Hamiltonian and the eigenvectors stored as columns, 
  !   this subroutine
  !             -  modifies the u_dagger array so that degenerate subspaces 
  !                diagonalize the perturbation
  !		-  computes the change in eigenvalues to 1st order.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  complex(dp)    :: dh(num_wann,num_wann)
  integer        :: degenerate_sets(num_wann,num_wann)
  integer        :: number_of_sets

  ! input/output variables
  complex(dp)    :: u_dagger(num_wann,num_wann)
  
  ! output variables
  real(dp)       :: d_eig(num_wann)


  ! local variables
  integer        :: i_set, n_subset, n_subset_max
  integer        :: n_bnd, n_bnd_current


  ! diagonalize the perturbation, for all degeneracies
 
  n_bnd_current = 0 

  do i_set = 1, number_of_sets

        call find_n_subset_max(degenerate_sets,i_set,n_subset_max)

        if (n_subset_max == 1)  then
           ! Subspace is not degenerate. No need to change the 
           ! eigenvectors, the change in eigenvalue is trivial.

            n_bnd_current = n_bnd_current + 1

            n_bnd  = degenerate_sets(i_set,1)

            ! remember! dot_product conjugates its first argument...
            d_eig(n_bnd_current) =  dot_product(u_dagger(:,n_bnd),  &
                                    matmul(dh,u_dagger(:,n_bnd)))

        else
            ! the subspace is degenerate! Apply diagonalization algorithm

            call rotate_sub_u_2(n_subset_max,                           &
                                degenerate_sets(i_set,1:n_subset_max),  &
                                u_dagger,dh,                            &
                     d_eig(n_bnd_current + 1:n_bnd_current+n_subset_max))
            
            n_bnd_current = n_bnd_current+n_subset_max


        end if

  end do

  end subroutine diagonalize_perturbation


  subroutine diagonalize_degeneracies(eig,u_dagger,dh,degenerate)
  !----------------------------------------------------------------------------!
  !   Given a list of eigenvalues, a perturbation hamiltonian, and 
  !   the eigenvectors stored as columns, this subroutine:
  !		- modifies the u_dagger array so that degenerate subspaces 
  !               diagonalize the perturbation.
  !		- computes the eigenvalue perturbation
  !		- computes the states perturbation
  !
  !   It cannot be assumed that the eigenvalues are in order!!!
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  real   (dp)    :: eig(num_wann)
  complex(dp)    :: u_dagger(num_wann,num_wann)
  complex(dp)    :: dh(num_wann,num_wann)
  logical        :: degenerate(num_wann)


  ! local variables
  integer        :: n_bnd,   m_bnd
  integer        :: i_index, i_block
  integer        :: degenerate_block(num_wann,num_wann)
  integer        :: degenerate_size(num_wann)
  logical        :: first_step

  logical        :: found(num_wann) 
  integer        :: iloop

  ! initialize
  found            = .false.
  degenerate       = .false.
  degenerate_block = 0
  degenerate_size  = 0
  i_block          = 0
  i_index          = 0


  n_bnd = 1
  ! tabulate the degeneracies
  do
     if (n_bnd == num_wann) exit

     do while ( found(n_bnd) .and. n_bnd < num_wann)
        n_bnd = n_bnd + 1
     end do

     found(n_bnd) = .true.
     first_step   = .true. 

     do m_bnd = n_bnd+1,num_wann

        if (.not. found(m_bnd) .and. abs(eig(n_bnd)-eig(m_bnd)) < eps_6 ) then

             found(m_bnd) = .true.
             
             if (first_step) then

                   first_step = .false.
                   i_block    = i_block + 1
                   degenerate_size(i_block)    = 2
                   degenerate_block(i_block,1) = n_bnd
                   degenerate_block(i_block,2) = m_bnd
                   i_index    = 2

                   degenerate(n_bnd)   = .true.
                   degenerate(m_bnd)   = .true.
              else

                   i_index = i_index + 1 
                   degenerate_size(i_block)          = degenerate_size(i_block)+1
                   degenerate_block(i_block,i_index) = m_bnd
                   degenerate(m_bnd) = .true.
          
              end if

        end if

     end do


  end do

  ! diagonalize the perturbation, for all degeneracies

  i_block = 0
  do i_block = 1, num_wann
        if (degenerate_size(i_block) /= 0)  then

                write(*,*) 'deg size = ',degenerate_size(i_block)

                call rotate_sub_u( degenerate_size(i_block),      &
                           degenerate_block(i_block,:),           &
                           u_dagger,dh)
        end if
  end do

  end subroutine diagonalize_degeneracies



  subroutine exponentiate_perturbation(u_k_dagger,du_k_dagger,u_kdk_dagger)
  !----------------------------------------------------------------------------!
  !   The change in the states obtained from perturbation theory
  !   yield an anti-hermitian matrix. It is thus possible to express
  !   
  !  	U^dagger(k+dk) ~  U^dagger(k) [ 1 -i {i U(k)*dU^dagger(k) } ]
  !  	               ~  U^dagger(k) exp[ U(k)*dU^dagger(k)] (This is now unitary!)
  !                              
  !	Define    M = i U(k)*dU^dagger(k)
  !			M is hermitian 
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  complex(dp)    :: u_k_dagger(num_wann,num_wann)
  complex(dp)    :: du_k_dagger(num_wann,num_wann)

  ! ouptut variable
  complex(dp)    :: u_kdk_dagger(num_wann,num_wann)

  ! local variables
  complex(dp)    :: Mat(num_wann,num_wann)

  integer        :: n_bnd

  ! diagonalization variables
  integer        :: nfound, info
  real(dp)       :: eig(num_wann)

  complex(dp)    :: T(num_wann,num_wann)
  complex(dp)    :: T_dagger(num_wann,num_wann)
  complex(dp)    :: tmp(num_wann,num_wann)

  integer        :: lwork, lrwork, liwork

  complex(dp)    :: work(2*num_wann)
  real   (dp)    :: rwork(24*num_wann)
  integer        :: iwork(10*num_wann)
  integer        :: isuppz(2*num_wann)


  ! size of work arrays (already hard coded above! this is just for clarity)
  lwork  = 2*num_wann
  liwork = 10*num_wann
  lrwork = 24*num_wann


  ! build the Hermitian matrix iM 

  Mat(:,:) = cmplx_i*matmul(transpose(conjg( u_k_dagger(:,:))),du_k_dagger(:,:))

  ! Diagonalizing routine from lapack

  call zheevr('V', 'A', 'U', num_wann, Mat, num_wann,         &
               0.0_dp, 0.0_dp, 0, 0, 0.0_dp, nfound,          &
               eig, T, num_wann,                              &
               isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)

  T_dagger(:,:) = transpose(conjg(T(:,:)))

  ! initialize
  u_kdk_dagger(:,:)  = cmplx_0

  ! fill with exponentiated eigenvalues
  do n_bnd= 1,num_wann
     u_kdk_dagger(n_bnd,n_bnd)  = exp(-cmplx_i*eig(n_bnd))
  end do

  ! return to original basis
  tmp(:,:)          = matmul(matmul(T(:,:), u_kdk_dagger(:,:)),T_dagger(:,:))
  u_kdk_dagger(:,:) = matmul(u_k_dagger(:,:),tmp(:,:))

  end subroutine exponentiate_perturbation


  subroutine rotate_sub_u_2(deg_size,deg_indices,u_dagger,dh,d_eig)
  !----------------------------------------------------------------------------!
  !   Given the indices of degenerate eigenvectors, stored as columns in u_dagger,
  !   this subroutine modifies the eigenvectors so that they diagonalize 
  !   the perturbation.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: deg_size
  integer        :: deg_indices(deg_size)
  complex(dp)    :: u_dagger(num_wann,num_wann)
  complex(dp)    :: dh(num_wann,num_wann)

  complex(dp)    :: new_columns(num_wann,deg_size)

  ! local variables
  complex(dp)    :: u_dh_u(deg_size,deg_size)
  integer        :: iloop, jloop, lloop, kloop
  integer        :: n_bnd1, n_bnd2

  ! diagonalization variables
  integer        :: nfound, info
  real(dp)       :: d_eig(deg_size)

  complex(dp)    :: diag_coefficients(deg_size,deg_size)

  integer        :: lwork, lrwork, liwork

  complex(dp)    :: work(2*deg_size)
  real   (dp)    :: rwork(24*deg_size)
  integer        :: iwork(10*deg_size)
  integer        :: isuppz(2*deg_size)


  ! size of work arrays (already hard coded above! this is just for clarity)
  lwork  = 2*deg_size
  liwork = 10*deg_size
  lrwork = 24*deg_size


  ! build the perturbation Hamiltonian in the degenerate sub-space

  do iloop = 1,deg_size
     n_bnd1 = deg_indices(iloop)

     do jloop = 1,deg_size
        n_bnd2 = deg_indices(jloop)

                u_dh_u(iloop,jloop)  =                         &
                        dot_product(u_dagger(:,n_bnd1),        &
                          matmul(dh,u_dagger(:,n_bnd2)))

     end do
  end do


  ! Diagonalizing routine from lapack

  call zheevr('V', 'A', 'U', deg_size, u_dh_u, deg_size,      &
               0.0_dp, 0.0_dp, 0, 0, 0.0_dp, nfound,          &
               d_eig, diag_coefficients, deg_size,            &
               isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)

  ! modify the u_dagger matrix according to the diagonalized states
  new_columns = cmplx_0

  do iloop = 1,deg_size
     do jloop = 1,deg_size

          n_bnd1 = deg_indices(jloop)
          new_columns(:,iloop) = new_columns(:,iloop) + &
                diag_coefficients(jloop,iloop)*u_dagger(:,n_bnd1)

     end do
  end do

  do iloop = 1,deg_size
     n_bnd1             = deg_indices(iloop)
     u_dagger(:,n_bnd1) = new_columns(:,iloop)
  end do


  end subroutine rotate_sub_u_2


  subroutine rotate_sub_u(deg_size,deg_indices,u_dagger,dh)
  !----------------------------------------------------------------------------!
  !   Given the indices of degenerate eigenvectors, stored as columns in u_dagger,
  !   this subroutine modifies the eigenvectors so that they diagonalize 
  !   the perturbation.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  integer        :: deg_size
  integer        :: deg_indices(deg_size)
  complex(dp)    :: u_dagger(num_wann,num_wann)
  complex(dp)    :: dh(num_wann,num_wann)

  complex(dp)    :: new_columns(num_wann,deg_size)

  ! local variables
  complex(dp)    :: u_dh_u(deg_size,deg_size)

  integer        :: iloop, jloop, lloop, kloop
  integer        :: n_bnd1, n_bnd2

  ! diagonalization variables
  integer        :: nfound, info
  complex(dp)    :: diag_coefficients(deg_size,deg_size)

  complex(dp),allocatable    :: u_dh_u_pack(:)

  real(dp),   allocatable    :: delta_eig(:) 
  real(dp)   ,allocatable    :: rwork(:)
  complex(dp),allocatable    :: cwork(:)

  integer,    allocatable    :: iwork(:), ifail(:)

  ! allocate

  allocate(u_dh_u_pack(deg_size*(deg_size+1)/2))
  allocate(rwork(7*deg_size))
  allocate(cwork(2*deg_size))
  allocate(iwork(5*deg_size))
  allocate(ifail(deg_size))
  allocate(delta_eig(deg_size))

  ! build the Hamiltonian in the degenerate sub-space

  do iloop = 1,deg_size
     n_bnd1 = deg_indices(iloop)

     do jloop = 1,deg_size
        n_bnd2 = deg_indices(jloop)

                u_dh_u(iloop,jloop)  =                         &
                        dot_product(u_dagger(:,n_bnd1),        &
                          matmul(dh,u_dagger(:,n_bnd2)))

     end do
  end do

  ! pack the perturbation, for the diagonalization routine
     do jloop=1,deg_size
       do iloop=1,jloop
         u_dh_u_pack(iloop+((jloop-1)*jloop)/2) =  u_dh_u(iloop,jloop)
       enddo
     enddo


  ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
  call ZHPEVX(    'V',    'A',    'U',    deg_size,    u_dh_u_pack,       &
               0.0_dp, 0.0_dp,      0,           0,        -1.0_dp,       &
               nfound,    delta_eig,                                      &
                          diag_coefficients  ,                            &
             deg_size,  cwork,  rwork,       iwork,          ifail, info)



  ! modify the u_dagger matrix according to the diagonalized states
  new_columns = cmplx_0

  do iloop = 1,deg_size
     do jloop = 1,deg_size

          n_bnd1 = deg_indices(jloop)
          new_columns(:,iloop) = new_columns(:,iloop) + &
                diag_coefficients(jloop,iloop)*u_dagger(:,n_bnd1)

     end do
  end do

  do iloop = 1,deg_size
     n_bnd1             = deg_indices(iloop)
     u_dagger(:,n_bnd1) = new_columns(:,iloop)
  end do

  deallocate(u_dh_u_pack)
  deallocate(rwork)
  deallocate(cwork)
  deallocate(iwork)
  deallocate(ifail)
  deallocate(delta_eig)

  end subroutine rotate_sub_u


  subroutine check_permutation_matrix(permutation_matrix,sze,test)
  !------------------------------------------------------------------
  ! This subroutine checks if a given matrix is a permutation matrix
  !------------------------------------------------------------------

  integer       :: sze, p
  integer       :: permutation_matrix(sze,sze)
  logical       :: test

  integer       :: i,j 
  integer       :: test_matrix(sze,sze)
  test = .true.

  do i = 1, sze
     do j = 1, sze
        p = permutation_matrix(i,j)
        if (p /= 0 .and. p /= 1) then
           test = .false.
           exit
        end if
     end do
     if ( .not. test) exit
  end do

  if ( test) then
     test_matrix(:,:) = matmul(transpose(permutation_matrix(:,:)),   &
                                         permutation_matrix(:,:))
     
     do i=1,sze
        test_matrix(i,i) = test_matrix(i,i) - 1 
     end do


        if ( sum(test_matrix(:,:)**2) /= 0) then
           test = .false.
        end if
  end if

  end subroutine check_permutation_matrix




  subroutine resolve_band_ordering(nk_vec_path,num_wann,    &
			list_k_path,H_int_k,eig_int_k,u_int_k)
  !------------------------------------------------------------------
  ! This subroutine orders the bands
  !------------------------------------------------------------------

  implicit none

  ! input variables
  integer        :: nk_vec_path, num_wann
  real(dp)       :: list_k_path(3,nk_vec_path)
  complex(dp)    :: u_int_k(num_wann,num_wann,nk_vec_path)
  complex(dp)    :: H_int_k(num_wann,num_wann,nk_vec_path)


  ! input/output variables
  real(dp)       :: eig_int_k(num_wann,nk_vec_path)

  ! local variables

  real(dp)       :: Ek1(num_wann), Ek2(num_wann) 

  complex(dp)    :: u_k1 (num_wann,num_wann), u_k2(num_wann,num_wann)
  complex(dp)    :: u_k1_dagger (num_wann,num_wann), u_k2_dagger(num_wann,num_wann)

  complex(dp)    :: H_k1 (num_wann,num_wann), H_k2 (num_wann,num_wann)
  complex(dp)    :: d_Hk (num_wann,num_wann)
  complex(dp)    :: A_k1(num_wann,num_wann), E_Ak1(num_wann,num_wann)

  complex(dp)    :: V_dagger(num_wann,num_wann)

  integer        :: permutation(num_wann,num_wann)
  logical        :: secondary_step

  real(dp)       :: k1(3), k2(3)
  real(dp)       :: dk(3)

  integer        :: io_unit_resolve

  integer        :: ik_loop
  integer        :: nb1, nb2

  complex(dp)    :: tmp_matrix(num_wann,num_wann)
  real(dp)       :: error 


  integer        :: number_of_sets_k1, number_of_sets_k2
  integer        :: degenerate_sets_k1(num_wann,num_wann)
  integer        :: degenerate_sets_k2(num_wann,num_wann)


  integer        :: n_subset_max
  integer        :: i_set, n_subset




   io_unit_resolve = find_free_unit()
   open(unit=io_unit_resolve,file=trim('band_crossing_resolution.test'),status='unknown')

   write(io_unit_resolve,20) '########################################################'
   write(io_unit_resolve,20) '# Bands Crossing resolution '
   write(io_unit_resolve,20) '# ----------------------------------------------        '
   write(io_unit_resolve,20) '# This file provides information regarding the band     ' 
   write(io_unit_resolve,20) '# connection algorithm. It will be useful for debugging.'
   write(io_unit_resolve,20) '# '
   write(io_unit_resolve,20) '########################################################'
   write(io_unit_resolve,20) '#'
   write(io_unit_resolve,20) '#'


   ! check for band crossings.
   do ik_loop = 2,nk_vec_path

      !-------------- extract some useful variables ---------------------
      Ek1(:)  =  eig_int_k(:,ik_loop-1)
      Ek2(:)  =  eig_int_k(:,ik_loop)

      k1(:)   =  list_k_path(:,ik_loop-1)
      k2(:)   =  list_k_path(:,ik_loop)
      dk(:)   =  k2(:)-k1(:)

      u_k1(:,:) =  u_int_k(:,:,ik_loop-1)
      u_k2(:,:) =  u_int_k(:,:,ik_loop)


      u_k1_dagger(:,:) = transpose(conjg(u_k1))
      u_k2_dagger(:,:) = transpose(conjg(u_k2))

      H_k1(:,:) =  H_int_k(:,:,ik_loop-1)
      H_k2(:,:) =  H_int_k(:,:,ik_loop)


      !-------------- talk to user ---------------------

      write(io_unit_resolve,20) '#-------------------------------------------------------'
      write(io_unit_resolve,24) 'ik_loop = ', ik_loop 
      write(io_unit_resolve,25) '  k = (',k1(1),',',k1(2),',',k1(3),')', &
                               'k+dk = (',k2(1),',',k2(2),',',k2(3),')'
      write(io_unit_resolve,20) ''


      !-------------- make sure matrices are aligned as expected ---------------------
      ! Make sure that u*H*u^dagger = E
      tmp_matrix(:,:) = matmul(u_k1,matmul(H_k1,u_k1_dagger))

      do nb1 = 1,num_wann
	tmp_matrix(nb1,nb1) = tmp_matrix(nb1,nb1) - Ek1(nb1) 
      end do

      error = sqrt(sum(matmul(transpose(conjg(tmp_matrix)),tmp_matrix)))

      write(io_unit_resolve,21) ' | u*H*u^d - E|  : ', error


      !-------------- build perturbation ---------------------
      d_Hk (:,:) =  H_k2(:,:) - H_k1(:,:)


      !-------------- find degenerate eigenvalue sets --------
      call find_degenerate_eigenvalue_sets(       &
               Ek1,degenerate_sets_k1,number_of_sets_k1)

      call find_degenerate_eigenvalue_sets(       &
               Ek2,degenerate_sets_k2,number_of_sets_k2)



      !-------------- output the degenerate sets --------

      call find_number_of_sets(degenerate_sets_k1,number_of_sets_k1)

      if (number_of_sets_k1 == num_wann) then
      	write(io_unit_resolve,20) ' no eigenvalue degeneracy'
      else 
      	write(io_unit_resolve,24) ' DEGENERATE: number of sets = ',number_of_sets_k1
      end if

      
      	write(io_unit_resolve,20) 'set index ---- eigenvalue (eV) -----  bands indices'
        
	do i_set = 1,number_of_sets_k1 
      	   call find_n_subset_max(degenerate_sets_k1,i_set,n_subset_max)
	
           nb1 = degenerate_sets_k1(i_set,1)
      	   write(io_unit_resolve,30) i_set, Ek1(nb1),degenerate_sets_k1(i_set,1:n_subset_max)
	end do

        !-------------- build parallel transport connection --------
        call build_parallel_transport_connection(u_k1_dagger,Ek1, d_Hk,  &
				degenerate_sets_k1, number_of_sets_k1,A_k1)


      	write(io_unit_resolve,20) 'Connection A_{mn} :'
        do nb1 = 1,num_wann
      	   write(io_unit_resolve,22) A_k1(nb1,:)
        end do

        !-------------- Exponentiate the connection --------
	call exponentiate_connection(A_k1,E_Ak1)

        !------------------------------------------------------
	! We expect, approximately, that 
        !  | PSI_{n k2} > = \sum_m | psi_{m k2} > V_{mn}^dagger
        ! where the PSI are parallel transported from k1.
        ! Thus,
        !  < psi_{m k1} | PSI_{n k2} > = 
	!	\sum_{l} <psi_{m k1}| psi_{l k2} > V_{l n}^dagger
        ! or
        !  e^{A} ~ u_k1 * u_k2^dagger * V^dagger
        !------------------------------------------------------


        V_dagger(:,:) = matmul(u_k2,matmul(u_k1_dagger,E_Ak1))

      	write(io_unit_resolve,20) 'unitary transform V_dagger :'
        do nb1 = 1,num_wann
      	   write(io_unit_resolve,22) V_dagger(nb1,:)
        end do

        !------- find closest permutation to V_dagger --------
  	call find_permutation(V_dagger,permutation,secondary_step)

        if (secondary_step) then
      		write(io_unit_resolve,20) 'permutation matrix (secondary step was necessary!)'
        else
      		write(io_unit_resolve,20) 'permutation matrix :'
        end if
        do nb1 = 1,num_wann
      	   write(io_unit_resolve,26) permutation(nb1,:)
        end do

        !------- Apply the permutation to the states and the eigenvalues--------


        eig_int_k(:,ik_loop) = matmul(Ek2(:),permutation(:,:))

        u_int_k(:,:,ik_loop) = matmul(transpose(permutation),u_k2)
   end do


  20 format(A)
  21 format(A,E16.8)
  22 format(100(2F8.4,5X))
  23 format(1000F8.4)
  24 format(A,I4)
  25 format(A,F8.4,A,F8.4,A,F8.4,A,3X,A,F8.4,A,F8.4,A,F8.4,A)
  26 format(100I2)
  30 format(I4,13X,F8.4,10X,100I4)

  end subroutine resolve_band_ordering

  subroutine build_parallel_transport_connection(u_k_dagger,Ek, d_Hk,  &
				degenerate_sets_k, number_of_sets_k,A_k)
  !----------------------------------------------------------------------------!
  ! This subroutine builds the parallel transport connection for the
  ! change in the wavefunctions. Care must be taken with respect to 
  ! states with degenerate eigenvalues.
  !
  !	The connection is given by 
  !		A_{mn}(k) = < psi_{m k} | d psi_{n k} >.
  !
  !	In the parallel transport gauge, we have
  !
  !    for non degenerate eigenvalues:
  !	A_{mn}(k) = < psi_{mk} | d H_k | psi_{nk} >
  !                -------------------------------
  !                  epsilon_{nk} - epsilon_{mk}
  !  
  !    for degenerate eigenvalues, A_{mn}(k) = 0, but the unperturbed states
  !    must be chosen to diagonalize the perturbation.
  !
  !	A build this way should be anti-hermitian.
  !----------------------------------------------------------------------------!
  implicit none

  ! input variables
  complex(dp)    :: d_Hk (num_wann,num_wann)
  real(dp)       :: Ek(num_wann)
  integer        :: number_of_sets_k
  integer        :: degenerate_sets_k(num_wann,num_wann)

  ! input/output variables
  complex(dp)    :: A_k(num_wann,num_wann)
  complex(dp)    :: u_k_dagger (num_wann,num_wann)

  ! local variables
  complex(dp)    :: u_k(num_wann,num_wann)
  real(dp)       :: dE (num_wann)
  integer        :: nb1, nb2
  real(dp)       :: e1, e2 

  complex(dp)    :: M_k(num_wann,num_wann)

  ! first, diagonalize the degenerate states  
  call diagonalize_perturbation(u_k_dagger,d_Hk,dE,degenerate_sets_k,number_of_sets_k)

  u_k(:,:) = conjg(transpose(u_k_dagger))
  ! then, compute the perturbation overlaps

  M_k(:,:) = matmul( u_k, matmul(d_Hk,u_k_dagger))

  ! next, build the connection
  A_k(:,:) = cmplx_0

  do  nb2 = 1, num_wann
    e2 = Ek(nb2)
    do  nb1 = 1, num_wann
        e1 = Ek(nb1)

	if (abs(e2-e1) < eps_8) cycle
		
	A_k(nb1,nb2) = M_k(nb1,nb2)/(e1-e2)

    end do
  end do

  end subroutine build_parallel_transport_connection

  subroutine exponentiate_connection(A_k,E_Ak)
  !----------------------------------------------------------------------------!
  !   
  !     Computes E_Ak = exp(A_k), which is a unitary matrix.
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  complex(dp)    :: A_k(num_wann,num_wann)

  ! ouptut variable
  complex(dp)    :: E_Ak(num_wann,num_wann)


  ! local variables
  complex(dp)    :: B_k(num_wann,num_wann)
  complex(dp)    :: B_pack(num_wann*(num_wann+1)/2)

  complex(dp)    :: V_dagger(num_wann,num_wann), V(num_wann,num_wann) 

  real(dp)       :: lambdas(num_wann)

  integer        :: i, j

  ! diagonalization variables
  integer        :: nfound, ifail, info


  complex(dp)    :: cwork(2*num_wann)
  real(dp)       :: rwork(7*num_wann)
  integer        :: iwork(5*num_wann)


  complex(dp)    :: tmp(num_wann,num_wann)

  ! define B = -iA, which is hermitan

  B_k(:,:) = -cmplx_i*A_k(:,:)

  ! size of work arrays (already hard coded above! this is just for clarity)



  ! Diagonalise B_k 
  do j=1,num_wann
     do i=1,j
        B_pack(i+((j-1)*j)/2)=B_k(i,j)
     enddo
  enddo

  ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
  call ZHPEVX(    'V',    'A',    'U',    num_wann,         B_pack,       &
               0.0_dp, 0.0_dp,      0,           0,        -1.0_dp,       &
               nfound,    lambdas,                                        &
                          V_dagger           ,                            &
             num_wann,  cwork,  rwork,       iwork,          ifail, info)


  V = conjg(transpose(V_dagger))
  ! initialize to zero
  tmp(:,:) = cmplx_0
 
  ! fill with diagonal eigenvalues
  do i= 1,num_wann
     tmp(i,i)  = exp(cmplx_i*lambdas(i))
  end do

  ! rotate back to original basis
  E_Ak(:,:)    = matmul(V_dagger , matmul(tmp,V))


  end subroutine exponentiate_connection

  subroutine find_permutation(V_dagger,permutation,secondary_step)
  !----------------------------------------------------------------------------!
  !   
  !     Find the best permutation approximation to V_dagger
  !----------------------------------------------------------------------------!
  implicit none
 
  ! input variables
  complex(dp)    :: V_dagger(num_wann,num_wann)

  ! ouptut variable
  integer        :: permutation(num_wann,num_wann)
  logical        :: secondary_step


  ! local variables
  integer        :: m_bnd, n_bnd

  logical        :: col_is_found(num_wann)
  logical        :: row_is_found(num_wann)
  logical        :: test

  real(dp)       :: tolerance
  real(dp)       :: x, x_max
  real(dp)       :: m_bnd_max, n_bnd_max

  ! a number larger than tolerance must be the largest number in 
  ! a column
  tolerance        = 1./sqrt(2.)

  col_is_found(:)  = .false.
  row_is_found(:)  = .false.

  ! initialize permutation
  permutation(:,:) = 0

  ! do a first pass in all the columns to see if
  ! there is always an element larger than tolerance




  do m_bnd = 1,num_wann
    do n_bnd = 1, num_wann
        if ( col_is_found(n_bnd) ) cycle

	x = abs(V_dagger(m_bnd,n_bnd))

	if (x > tolerance) then
    		permutation(m_bnd,n_bnd) = 1
		col_is_found(n_bnd) = .true.
		row_is_found(m_bnd) = .true.
        end if
    end do
  end do

  ! if all the columns are found, we are done! If not, let's zero in
  ! on what is left


  test = all(col_is_found(:))
  if (test) then
     secondary_step = .false.
  else
     secondary_step = .true.
  end if

  do 

    ! if all columns are found, we are done!
    test = all(col_is_found(:))
    if ( test ) exit
 
    x_max = ZERO
    ! find the largest element left over 
    do m_bnd = 1,num_wann
      if (row_is_found(m_bnd) ) cycle
      do n_bnd = 1, num_wann
         if (col_is_found(n_bnd) ) cycle

	 x = abs(V_dagger(m_bnd,n_bnd))

         if (x >= x_max) then
	    x_max = x 
            m_bnd_max = m_bnd
            n_bnd_max = n_bnd
         end if

      end do
    end do

!haritz
!    permutation(m_bnd_max,n_bnd_max) = 1
!    col_is_found(n_bnd_max) = .true.
!    row_is_found(m_bnd_max) = .true.
    permutation(int(m_bnd_max),int(n_bnd_max)) = 1
    col_is_found(int(n_bnd_max)) = .true.
    row_is_found(int(m_bnd_max)) = .true.
!haritz

  end do
  

  end subroutine find_permutation

!--------------------------------------------------------------------------------
!
end module intw_band_crossing
!
!--------------------------------------------------------------------------------
