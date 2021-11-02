!----------------------------------------------------------------------------!
!       intw project.
!
!----------------------------------------------------------------------------!
module intw_chi_0
!----------------------------------------------------------------------------!
!
!       This module contains all necessary components to compute the 
!       response function, chi0, and the dielectric function.
!
!----------------------------------------------------------------------------!

use intw_useful_constants
use intw_input_parameters
use intw_W90
use intw_band_crossing
use intw_kramers_kronig

  !
  implicit none
  !
  save
  !
        integer                 :: D_scheme, N_scheme, sze_scheme
          ! The dimension of the interpolation scheme, and the
          ! the order. 
        complex(dp),allocatable :: chi0(:,:,:,:,:,:,:)
          ! The response function chi0

        real(dp),allocatable :: re_chi0_tmp(:,:,:,:,:,:,:)
        real(dp),allocatable :: im_chi0_tmp(:,:,:,:,:,:,:)
          ! buffers for computations


        complex(dp),allocatable :: chi0_intw(:,:,:)
          ! The response function coefficients in the MLWF product basis


        real(dp),allocatable    :: list_hbar_omega(:)

          ! the frequency interval, in eV
        real(dp)                :: dhw

        real(dp) :: average_width, average_width_deg 
          ! The average width used by the adaptive scheme, in eV.
        integer  :: number_of_width, number_of_width_deg

 	! rectangular matrix containing information about the basis 
	complex(dp), allocatable  :: projectors_intw(:,:)
  	integer ,allocatable      :: list_m1_intw(:), list_m2_intw(:)

	integer                   :: basis_size_intw

contains

  subroutine create_chi0()
  !------------------------------------------------------------------------
  ! This subroutine initializes the chi0 array
  !------------------------------------------------------------------------

  implicit none

  complex(dp)   :: i_delta

  integer       :: iloop


  ! allocate
  allocate(list_hbar_omega(n_hw))

  allocate(chi0(n_hw,nG_shell_max,nG_shell_max,npol,npol,npol,npol))

  
  average_width   = ZERO
  number_of_width = 0

  average_width_deg   = ZERO
  number_of_width_deg = 0

  chi0 = cmplx_0
  ! generate the frequency list
  ! generate this global variable
  dhw  = (hw_max-hw_min)/(dble(n_hw-1))

  do iloop = 1,n_hw
     list_hbar_omega(iloop) = hw_min+dble(iloop-1)*dhw
  end do

  end subroutine create_chi0


  subroutine create_chi0_intw()
  !------------------------------------------------------------------------
  ! This subroutine initializes the chi0_intw array
  !------------------------------------------------------------------------

  implicit none

  complex(dp)   :: i_delta

  integer       :: iloop


  ! allocate
  allocate(list_hbar_omega(n_hw))

  allocate(chi0_intw(n_hw,basis_size_intw,basis_size_intw))

  
  average_width   = ZERO
  number_of_width = 0

  average_width_deg   = ZERO
  number_of_width_deg = 0

  chi0_intw = cmplx_0
  ! generate the frequency list
  ! generate this global variable
  dhw  = (hw_max-hw_min)/(dble(n_hw-1))

  do iloop = 1,n_hw
     list_hbar_omega(iloop) = hw_min+dble(iloop-1)*dhw
  end do

  end subroutine create_chi0_intw


  subroutine output_chi0()
  !------------------------------------------------------------------------
  ! This subroutine prints out chi0
  !------------------------------------------------------------------------
  implicit none


  real(dp)            :: volume,  rho_0, prefactor, spin_sum_factor

  real(dp),allocatable :: output_array(:)

  integer        :: io_unit
  integer        :: iG_1,  iG_2
  integer        :: ipol1, ipol2, ipol3, ipol4
  integer        :: ipol_loop

  integer        :: iw_loop

  character(256) :: filename, iG_1_str, iG_2_str

  ! header lines for spin response function. The exact length
  ! is 18(characters) x (1(hw) + 16(re) +16(im) ) =  594 
  character(594) :: header_line1, header_line2
  character(36)  :: header_bit

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

  allocate(output_array(2*nspin**4))

  volume = alat**3*abs(                                &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))

  rho_0     = 1.0_dp/volume

  

  chi0      = chi0/dble(nk1s*nk2s*nk3s)*rho_0*spin_sum_factor 


  do iG_1 = 1,nG_shell_max
     do iG_2 = 1,nG_shell_max

        io_unit = find_free_unit()

        write(iG_1_str,100) iG_1
        write(iG_2_str,100) iG_2

        filename = trim(chi0_file_prefix)//'_'//       &
                        trim(adjustl(iG_1_str))//'_'// &
                        trim(adjustl(iG_2_str))//'.dat'


        open(unit = io_unit, file = filename)

        write(io_unit,5) '#=========================================================='
        write(io_unit,5) '#   This file contains the computed response  '
        write(io_unit,5) '#   function chi0 for the following parameters:'
        write(io_unit,5) '#                                             '
        write(io_unit,5) '#   The q vector, in crystal coordinates:     '
        write(io_unit,5) '#                                             '
        write(io_unit,10)'#         q(1)  = ', qpt1
        write(io_unit,10)'#         q(2)  = ', qpt2
        write(io_unit,10)'#         q(3)  = ', qpt3
        write(io_unit,5) '#                                             '
        write(io_unit,5) '#   The 1st G-vector, in crystal coordinates:'
        write(io_unit,5) '#                                             '
        write(io_unit,15)'#         G1(1) = ', gvec(1,iG_1)
        write(io_unit,15)'#         G1(2) = ', gvec(2,iG_1)
        write(io_unit,15)'#         G1(3) = ', gvec(3,iG_1)
        write(io_unit,5) '#                                             '
        write(io_unit,5) '#   The 2nd G-vector, in crystal coordinates:'
        write(io_unit,5) '#                                             '
        write(io_unit,15)'#         G2(1) = ', gvec(1,iG_2)
        write(io_unit,15)'#         G2(2) = ', gvec(2,iG_2)
        write(io_unit,15)'#         G2(3) = ', gvec(3,iG_2)
        write(io_unit,5) '#                                             '


	if ( nspin ==1 ) then
           write(io_unit,5) '# nspin = 1: paramagnetic case                '
           write(io_unit,5) '#    The sum on degenerate spins yields a     '
           write(io_unit,5) '#    factor of 2, which is included in chi0   '
           write(io_unit,5) '#                                             '
           write(io_unit,5) '#  1st   column: hw (eV)                      '
           write(io_unit,5) '#  2st   column: Re[chi0](1/a0^3eV)           '
           write(io_unit,5) '#  3rd   column: Im[chi0](1/a0^3eV)           '                     
           write(io_unit,5) '#                                             '                     
           write(io_unit,5) '#=========================================================='
           write(io_unit,5) '#    hw (eV)            Re[chi0]          Im[chi0]'
           write(io_unit,5) '#=========================================================='

	else if ( nspin ==2 .and. .not. magnon) then
           write(io_unit,5) '# nspin = 2: non-paramagnetic case            '
           write(io_unit,5) '#    Spin is no longer degenerate; chi0 has   '
           write(io_unit,5) '#    four spin indices, thus 16 complex       '
           write(io_unit,5) '#    components:                              '
           write(io_unit,5) '#    chi_0_{s1s2,s3s4}                        '
           write(io_unit,5) '#    where si = 1,2 for up, down              '
           write(io_unit,5) '#                                             '
           write(io_unit,5) '#  1st   column: hw (eV)                      '
           write(io_unit,5) '#  even  column: Re[chi0_{s1s2s3s4}](1/a0^3eV)'
           write(io_unit,5) '#  odd   column: Im[chi0_{s1s2s3s4}](1/a0^3eV)'
           write(io_unit,5) '#                                             '
           write(io_unit,5) '#  from column to column, s4 loop fastest     '
           write(io_unit,5) '#  and s1 slowest.                            '
           write(io_unit,5) '#=========================================================='

	   ! build an informative string which will tell the user what she's looking at
	   header_line1 = '#    hw (eV)       |'
	   header_line2 = '#                  |'

           do ipol1=1,nspin
            do ipol2=1,nspin 
             do ipol3=1,nspin
              do ipol4=1,nspin 

	        write(header_bit,25) ipol1,ipol2,ipol3,ipol4
                header_line1 = trim(header_line1)//header_bit
                header_line2 = trim(header_line2)//'     Re[chi0]          Im[chi0]    |'

              end do ! ipol4
             end do !ipol3
            end do !ipol2
           end do !ipol1

           write(io_unit,5) header_line1 
           write(io_unit,5) header_line2 
           write(io_unit,5) '#=========================================================='

	else if ( nspin ==2 .and.  magnon) then
           write(io_unit,5) '# nspin = 2: non-paramagnetic case            '
           write(io_unit,5) '#    Spin is no longer degenerate; chi0 has   '
           write(io_unit,5) '#    four spin indices, thus 16 complex       '
           write(io_unit,5) '#    components:                              '
           write(io_unit,5) '#    chi_0_{s1s2,s3s4}                        '
           write(io_unit,5) '#    where si = 1,2 for up, down              '
           write(io_unit,5) '# magnon = True:                              '
           write(io_unit,5) '#    only chi_0_{12,21} and  chi_0_{21,12}    '
           write(io_unit,5) '#    will be written to file.                 '
           write(io_unit,5) '#                                             '
           write(io_unit,5) '#    units of chi0:  1/a0^3eV                 '
           write(io_unit,5) '#                                             '
           write(io_unit,5) '#======================================================='
	   ! build an informative string which will tell the user what she's looking at
	   header_line1 = '#    hw (eV)       |'
	   header_line2 = '#                  |'

           ipol1=1; ipol2=2; ipol3=2; ipol4=1
	   write(header_bit,25) ipol1,ipol2,ipol3,ipol4
           header_line1 = trim(header_line1)//header_bit
           header_line2 = trim(header_line2)//'     Re[chi0]          Im[chi0]    |'

           write(io_unit,5) trim(header_line1)
           write(io_unit,5) trim(header_line2)
           write(io_unit,5) '#======================================================='

        end if

	if ( .not.  magnon) then
          do iw_loop = 1,n_hw

          ! organize the response function for output
            ipol_loop = 0
            do ipol1=1,nspin
             do ipol2=1,nspin 
              do ipol3=1,nspin
               do ipol4=1,nspin 

                ipol_loop = ipol_loop + 1 
  
                output_array(ipol_loop ) =   &
                   real(chi0(iw_loop,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4))

                ipol_loop = ipol_loop + 1 
                output_array(ipol_loop ) =   &
                   aimag(chi0(iw_loop,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4))

               end do !ipol4
              end do !ipol3
             end do !ipol2
            end do !ipol1


            ! write to file
            write(io_unit,20) list_hbar_omega(iw_loop), output_array(:)

          end do ! iw_loop

	else if ( magnon) then

          do iw_loop = 1,n_hw

            ipol_loop = 0

            ! Only account for the magnon response functions
            ipol_loop = ipol_loop + 1 
            output_array(ipol_loop ) =   &
                   real(chi0(iw_loop,iG_1,iG_2,1,1,1,1))

            ipol_loop = ipol_loop + 1 
            output_array(ipol_loop ) =   &
                   aimag(chi0(iw_loop,iG_1,iG_2,1,1,1,1))

            ! write to file
            write(io_unit,20) list_hbar_omega(iw_loop), output_array(1:ipol_loop)

          end do

	end if 

        close(io_unit)
     end do !iG_2
  end do !iG_1


  5   format(A)
  10  format(A,3G18.8E3)
  15  format(A,I6)
  20  format(100G18.8E3)
 
  25  format('             chi0_{',4I1,'}           |')

  100 format(I4)

  end subroutine output_chi0



  subroutine get_gradients_and_occupations(num_pack,eig, u, dH_pack,   &
             degenerate_sets,number_of_sets,f,gradients)
  !------------------------------------------------------------------------
  ! This convenience subroutine performs various tasks pertaining
  ! to computing the gradient of the Hamiltonian and determining
  ! what set of states are occupied. 
  !
  ! The gradients array will have the same shape as degenerate_sets,
  ! with an extra dimension for alpha = x,y,z.
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: num_pack

  real(dp)     :: eig (num_wann)

  complex(dp)  :: u (num_wann,num_wann)

  complex(dp)  :: dH_pack (num_pack,3)

  integer      :: degenerate_sets (num_wann,num_wann)

  integer      :: number_of_sets

  ! output variables
  integer      :: f(num_wann)
  real(dp)     :: gradients(num_wann,num_wann,3)


  ! computation variables
  complex(dp)  :: dH(num_wann,num_wann,3)

  real(dp)     :: e

  integer      :: i_set, alpha
  integer      :: i_pack
  integer      :: m_bnd, n_bnd
  integer      :: i_subset, j_subset

  integer      :: n_subset_max

  integer      :: num_pack_sub_dH, i_pack_sub_dH

  complex(dp),allocatable  :: sub_dH_pack(:)


  ! diagonalization variables
  complex(dp),allocatable  :: work(:)
  real   (dp),allocatable  :: rwork(:)
  integer    ,allocatable  :: iwork(:) 
  integer                  :: lwork, lrwork,  liwork
  integer                  :: info


  complex(dp)  :: z_garbage(1,1)

  ! initialize
  f(:)      =   0
  gradients =   ZERO

  ! unpack gradient of hamiltonian

  do n_bnd =1,num_wann
     do m_bnd = 1, n_bnd
           i_pack = m_bnd+((n_bnd-1)*n_bnd)/2

           dH(m_bnd, n_bnd,:) = dH_pack(i_pack,:)
           dH(n_bnd, m_bnd,:) = conjg(dH(m_bnd, n_bnd,:))
     end do

     dH(n_bnd, n_bnd,:) = real(dH(n_bnd, n_bnd,:))
  end do

  ! iterate on sets
  do i_set = 1, number_of_sets

        call find_n_subset_max(degenerate_sets,i_set,n_subset_max)

        ! get the eigenvalue for this subset
        m_bnd = degenerate_sets(i_set,1)
        e     = eig (m_bnd)

        ! is this set occupied?
        if (e  <= chemical_potential ) then
                f(i_set) = 1
        end if

        ! get the gradients
        if ( n_subset_max == 1) then
            ! non-degenerate case, no need to diagonalize

            do alpha =1,3 
                gradients (i_set,1,alpha)  =   &
                ! remember, dot_product conjugates the first argument!
                dot_product(conjg(u(m_bnd,:)),     &
                            matmul(dH(:,:,alpha),conjg(u(m_bnd,:))))
            end do 

        else
            ! degeneracies are present: the perturbation must be diagonalized

            num_pack_sub_dH = n_subset_max*(n_subset_max+1)/2

            lwork  = n_subset_max
            lrwork = n_subset_max
            liwork = 1

            allocate(sub_dH_pack(num_pack_sub_dH))
            allocate(work (lwork))
            allocate(rwork(lrwork))
            allocate(iwork(liwork))


            do alpha =1,3 

               sub_dH_pack(:) = cmplx_0

               do j_subset =1,n_subset_max
                     n_bnd = degenerate_sets(i_set,j_subset)

                  do i_subset =1,j_subset
                        m_bnd = degenerate_sets(i_set,i_subset)

                        i_pack_sub_dH = i_subset+((j_subset-1)*j_subset)/2

                        sub_dH_pack( i_pack_sub_dH ) =              &
                            dot_product(conjg(u(m_bnd,:)),          &
                            matmul(dH(:,:,alpha),conjg(u(n_bnd,:))))

                  end do ! i_subset
               end do ! j_subset

               ! diagonalize ! Routine only calculates eigenvalues.
               ! see MKL manual for meaning of various variables
               call zhpevd('N', 'U', n_subset_max, sub_dH_pack,   &
                        gradients (i_set,1:n_subset_max,alpha),   &
                        z_garbage, 1, work, lwork, rwork, lrwork, &
                        iwork, liwork, info)


                if ( info /= 0) then
                        write (*,*) '***********************************'
                        write (*,*) '* info /= 0 in zhpevd!            *'
                        write (*,*) '* problem in subroutine           *'
                        write (*,*) '* get_gradients_and_occupations   *'
                        write (*,*) '***********************************'
                end if
            end do ! alpha

            deallocate(sub_dH_pack)
            deallocate(work )
            deallocate(rwork)
            deallocate(iwork)

        end if
           


     end do  !i_set

  end subroutine get_gradients_and_occupations



  subroutine compute_adaptive_widths_and_occupations(            &
                        nk_vec_max,nk_vec,num_pack,              &
                        eig_ks, eig_ksq, u_ks, u_ksq,            &
                        dH_ks, dH_ksq, adaptive_widths,fk_m_fkq)
  !------------------------------------------------------------------------
  ! This subroutine computes the widths delta_{n1n2}(k) that should be
  ! used at a given k point in the computation of chi0. 
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max, num_pack

  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)

  complex(dp)  :: u_ks (num_wann,num_wann,nk_vec_max)
  complex(dp)  :: u_ksq(num_wann,num_wann,nk_vec_max)

  complex(dp)  :: dH_ks  (num_pack,3,nk_vec_max)
  complex(dp)  :: dH_ksq (num_pack,3,nk_vec_max)

  ! output variables
  real(dp)     :: adaptive_widths (num_wann,num_wann,nk_vec_max)
  integer      :: fk_m_fkq        (num_wann,num_wann,nk_vec_max)

  ! computation variables
  integer      :: degenerate_sets_k (num_wann,num_wann)
  integer      :: degenerate_sets_kq(num_wann,num_wann)

  integer      :: number_of_sets_k, number_of_sets_kq


  real(dp)     :: gradients_k (num_wann,num_wann,3)
  real(dp)     :: gradients_kq(num_wann,num_wann,3)

  real(dp)     :: gr(3)

  real(dp)     :: dk (3)
  real(dp)     :: delta

  integer      :: k_vec_loop

  integer      :: i_set_k, i_set_kq
  integer      :: i_subset_k, i_subset_kq
  integer      :: i_pack

  integer      :: n_subset_max_k, n_subset_max_kq 
  integer      :: m_bnd, n_bnd

  integer      :: f_set_k(num_wann), f_set_kq(num_wann) 
  integer      :: f_m_f

  ! initialize 
  adaptive_widths  = ZERO
  fk_m_fkq         =  0 


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
                eig_ks(:,k_vec_loop),                  &
                u_ks    (:,:,k_vec_loop),              &
                dH_ks   (:,:,k_vec_loop),              &
                degenerate_sets_k,number_of_sets_k,f_set_k,gradients_k)


      call find_degenerate_eigenvalue_sets(       &
                     eig_ksq(:,k_vec_loop),            &
                degenerate_sets_kq,number_of_sets_kq)

      call get_gradients_and_occupations(num_pack,      &
                eig_ksq    (:,k_vec_loop),              &
                u_ksq    (:,:,k_vec_loop),              &
                dH_ksq   (:,:,k_vec_loop),              &
                degenerate_sets_kq,number_of_sets_kq,f_set_kq,gradients_kq)

      ! loop on sets and build widths

      do i_set_k=1, number_of_sets_k
         call find_n_subset_max(degenerate_sets_k,i_set_k,n_subset_max_k)

         do i_set_kq=1, number_of_sets_kq
            call find_n_subset_max(degenerate_sets_kq,i_set_kq,n_subset_max_kq)

            f_m_f = f_set_k(i_set_k)- f_set_kq(i_set_kq)

            if ( f_m_f /= 0 ) then 

                    delta = ZERO

                    do i_subset_k = 1, n_subset_max_k
                          do i_subset_kq = 1, n_subset_max_kq

                        	gr(:) = gradients_k(i_set_k,i_subset_k,:)     &
                                      - gradients_kq(i_set_kq,i_subset_kq,:)
        
                          	delta = delta  +           &
                        		(gr(1)*dk(1))**2 + &
		                        (gr(2)*dk(2))**2 + &
	                                (gr(3)*dk(3))**2

                        
                          end do ! i_subset_kq
                end do !i_subset_k

                delta = adaptive_width_coeff*            &
                        sqrt(delta/dble(3*n_subset_max_k*n_subset_max_kq))


                do i_subset_k = 1, n_subset_max_k
                   m_bnd = degenerate_sets_k(i_set_k,i_subset_k)
                   do i_subset_kq = 1, n_subset_max_kq
                      n_bnd = degenerate_sets_kq(i_set_kq,i_subset_kq)

                        fk_m_fkq(m_bnd,n_bnd,k_vec_loop)        = f_m_f
                        adaptive_widths(m_bnd,n_bnd,k_vec_loop) = delta

                        average_width   = average_width + delta
                        number_of_width = number_of_width + 1 

                      if (n_subset_max_kq /= 1 .or. n_subset_max_k /= 1) then
                           average_width_deg   = average_width_deg + delta
                           number_of_width_deg = number_of_width_deg + 1 
                      end if

                       
                   end do ! i_subset_kq
                end do !i_subset_k

            end if

         end do !i_set_kq
      end do !i_set_k

  end do !k_vec_loop 

  end subroutine compute_adaptive_widths_and_occupations



  subroutine compute_fixed_widths_and_occupations(nk_vec_max,nk_vec, &
                        eig_ks, eig_ksq, widths,fk_m_fkq,delta)
  !------------------------------------------------------------------------
  ! This subroutine computes the widths delta_{n1n2}(k) that should be
  ! used at a given k point in the computation of chi0. 
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max, num_pack

  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)

  ! output variables
  real(dp)     :: widths    (num_wann,num_wann,nk_vec_max)
  integer      :: fk_m_fkq  (num_wann,num_wann,nk_vec_max)

  ! computation variables

  real(dp)     :: delta

  integer      :: k_vec_loop

  integer      :: m_bnd, n_bnd

  integer      :: f_set_k(num_wann), f_set_kq(num_wann) 
  integer      :: f_n, f_m, f_m_f

  ! initialize 
  widths    = ZERO
  fk_m_fkq  =  0 

  ! loop on all k-points in the block
  do k_vec_loop = 1,nk_vec

     do m_bnd = 1,num_wann

        f_m = 0
        if (eig_ks (m_bnd ,k_vec_loop) <= chemical_potential) then
           f_m = 1
        end if

        do n_bnd = 1,num_wann
           f_n = 0
           if (eig_ksq(n_bnd ,k_vec_loop) <= chemical_potential) then
              f_n = 1
           end if

           f_m_f = f_m-f_n
           if (f_m_f /= 0) then
               fk_m_fkq  (m_bnd,n_bnd,k_vec_loop) =  f_m_f
               widths    (m_bnd,n_bnd,k_vec_loop) =  delta
           end if
              

        end do ! n_bnd
     end do ! m_bnd
  end do ! k_vec_loop

  end subroutine compute_fixed_widths_and_occupations

  subroutine apply_kramers_kronig_to_chi0()
  !------------------------------------------------------------------------
  ! This subroutine computes the "real part" of chi0 using the
  ! "imaginary part", using the Kramers-Kronig (KK) FFT algorithm. 
  ! Note that since the KK equations are linear, they can be applied to
  ! complex functions, namely the "real part" need not be real, nor does
  ! the "imaginary part" need be imaginary.
  !------------------------------------------------------------------------
  implicit none


  real(dp)       :: wmin, wmax

  complex(dp),allocatable  :: list_Fi_w(:)
  complex(dp),allocatable  :: list_Fr_w(:)

  logical  :: double_data
  integer  :: iw_loop
  integer  :: Mfft_input
  integer  :: iG_1, iG_2
  integer  :: ipol1, jpol1, ipol2, jpol2


  ! set up the Kramers-Kronig scheme


  wmax = hw_max

  if (hw_min == -hw_max) then
     ! use the data as is
     wmin =  hw_min
     double_data = .False.
     Mfft_input  = n_hw-1
  else if( abs(hw_min) < eps_10) then
     ! symmetrize below w = 0
     wmin = -hw_max
     double_data = .True.
     Mfft_input  = 2*(n_hw-1)
  else
     ! something is wrong. This should be tested for 
     ! at the begining of the execution of the program...
     write(*,*) '**************************************'
     write(*,*) '*  Problems with frequencies:        *'
     write(*,*) '*  cannot perform Kramers-Kronig     *'
     write(*,*) '*     program stops.                 *'
     write(*,*) '**************************************'
     stop
  end if

 
  ! set up the KK scheme 
  call setup_kramers_kronig_arrays(wmin,wmax,Mfft_input)
  call setup_FFTW_kramers_kronig_plans()
  call get_W_trapezoid()
  call get_alpha_trapezoid()
  call get_beta_trapezoid()


  ! allocate the working arrays for the KK scheme
  allocate(list_Fi_w(Mfft+1))
  allocate(list_Fr_w(Mfft+1))

  ! loop on all components of the response function, which
  ! is assumed to only contain the "imaginary part" 


    do jpol2=1,npol
     do ipol2=1,npol
      do jpol1=1,npol
       do ipol1=1,npol
        do iG_2 = 1,nG_shell_max
         do iG_1 = 1,nG_shell_max

               ! dump the imaginary part of chi0 in the working array

               if ( double_data == .False. ) then
		  list_Fi_w(:) =  chi0(:,iG_1,iG_2,ipol1,jpol1,ipol2,jpol2)
	       else 
		  list_Fi_w(Mfft/2+1:Mfft+1) =  &
				chi0(:,iG_1,iG_2,ipol1,jpol1,ipol2,jpol2)

                  if ( magnon == .True.) then
		     ! if this is a magnon calculation, zero pad
                     ! on the negative part of the axis.
                  	do iw_loop =1, Mfft/2
		     		list_Fi_w(iw_loop)  =  cmplx_0
                  	end do
                  else
                  	do iw_loop =1, Mfft/2
		     		list_Fi_w(iw_loop)  =  -list_Fi_w(Mfft+2-iw_loop)
                  	end do
                  end if


               end if

               ! apply the KK-FFT scheme
               call apply_kramers_kronig_FFT(list_Fi_w,list_Fr_w)


               ! add the "real part" of chi0 to the "imaginary part"
               if ( double_data == .False. ) then
		  chi0(:,iG_1,iG_2,ipol1,jpol1,ipol2,jpol2) =      &
		      chi0(:,iG_1,iG_2,ipol1,jpol1,ipol2,jpol2) +  &
					          	list_Fr_w(:)
               else
		  chi0(:,iG_1,iG_2,ipol1,jpol1,ipol2,jpol2) =      &
		      chi0(:,iG_1,iG_2,ipol1,jpol1,ipol2,jpol2) +  &
			                  list_Fr_w(Mfft/2+1:Mfft+1)
               end if

         end do !iG_1
        end do !iG_2
       end do !ipol1
      end do !jpol1
     end do !ipol2
    end do !jpol2


  ! clean up
  deallocate(list_Fi_w)
  deallocate(list_Fr_w)

  call deallocate_kramers_kronig_arrays()

  end subroutine apply_kramers_kronig_to_chi0


  subroutine output_chi0_netcdf()
  !------------------------------------------------------------------------
  ! This subroutine prints out chi0 in nectdf format, along with
  ! relevant information to interpret the data
  !------------------------------------------------------------------------
  use netcdf
  implicit none


  real(dp)       :: volume,  rho_0, prefactor, spin_sum_factor


  integer        :: io_unit
  integer        :: iG_1,  iG_2
  integer        :: ipol1, ipol2, ipol3, ipol4
  integer        :: ipol_loop

  integer        :: iw_loop

  character(256) :: filename, iG_1_str, iG_2_str


  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
   
  integer :: var_hw_id, var_re_chi0_id,var_im_chi0_id, var_q_id, var_Gvect_id
  ! variables id: frequency, response function, q vector, G vectors 


  integer :: dim_hw_id, dim_chi0_id
  integer :: dim_G1_id, dim_G2_id, dim_space_id
  integer :: dim_s1_id, dim_s2_id,dim_s3_id, dim_s4_id
  ! dimension id: frequency, response function, 
  !               2 G vectors dimensions, 3D space dimension
  !		  spin dimensions

  integer :: dim_hw, dim_chi0, dim_G, dim_space
  ! dimensions: frequency, response function, q vector, G vectors
  !	        dimension of 3D space (3!)

  character (*), parameter :: hw_name       = "frequency"
  character (*), parameter :: G1_name       = "iG1"
  character (*), parameter :: G2_name       = "iG2"
  character (*), parameter :: Gvect_name    = "G_vectors"
  character (*), parameter :: q_name        = "q_vector"
  character (*), parameter :: re_chi0_name  = "real_chi_KS"
  character (*), parameter :: im_chi0_name  = "imaginary_chi_KS"
  character (*), parameter :: space_name    = "space_coordinates"
  character (*), parameter :: sigma1_name   = "sigma1"
  character (*), parameter :: sigma2_name   = "sigma2"
  character (*), parameter :: sigma3_name   = "sigma3"
  character (*), parameter :: sigma4_name   = "sigma4"

  ! Define the units
  character (*), parameter :: UNITS = "units"
  character (*), parameter :: SPIN_COM1   = "first_spin_component"
  character (*), parameter :: SPIN_COM2   = "second_spin_component"

  character (*), parameter :: hw_UNITS    = "eV"
  character (*), parameter :: chi0_UNITS  = "1/(eV*a_0^3)"
  character (*), parameter :: space_UNITS = "crystal_basis"

  character (*), parameter :: uddu  = "up-dn-dn-up"
  character (*), parameter :: duud  = "dn-up-up-dn"



  character(256) :: time_stamp, hostname
 integer         :: time_values(8)

  ! header lines for spin response function. The exact length
  ! is 18(characters) x (1(hw) + 16(re) +16(im) ) =  594 
  character(594) :: header_line1, header_line2
  character(36)  :: header_bit

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

  volume = alat**3*abs(                                &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))

  rho_0     = 1.0_dp/volume


  chi0      = chi0/dble(nk1s*nk2s*nk3s)*rho_0*spin_sum_factor 


  ! prepare the netcdf file, writing meta-data

  ! create file
  filename = trim(chi0_file_prefix)//'.nc'

  call check_netcdf(nf90_create(filename, NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Kohn-Sham susceptibility function") )

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
			"This is a magnon calculation: chi_{up-dn-dn-up} "))
  else if( nspin == 1) then
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"This is a density calculation: chi_{nn}" ))
  else
    call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "calculation", &
			"All spin contributions are computed: chi_{s1s2s3s4}" ))
  end if

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension, 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create one "G-vector" dimension, the number of G vectors
  dim_G = nG_shell_max
  call check_netcdf( nf90_def_dim(nc_id,G1_name, dim_G, dim_G1_id) )
  call check_netcdf( nf90_def_dim(nc_id,G2_name, dim_G, dim_G2_id) )

  ! create the "frequency" dimension, as well as the frequency variable
  dim_hw = n_hw
  call check_netcdf( nf90_def_dim(nc_id,hw_name, dim_hw, dim_hw_id) )

  if (.not. magnon .and. nspin /= 1) then
     ! create the "spin" dimensions
  	call check_netcdf( nf90_def_dim(nc_id,sigma1_name, 2, dim_s1_id) )
  	call check_netcdf( nf90_def_dim(nc_id,sigma2_name, 2, dim_s2_id) )
  	call check_netcdf( nf90_def_dim(nc_id,sigma3_name, 2, dim_s3_id) )
  	call check_netcdf( nf90_def_dim(nc_id,sigma4_name, 2, dim_s4_id) )
  end if

  !------------------------------
  ! define the variables
  !------------------------------
  ! create the "q-vector" array
  call check_netcdf( nf90_def_var(nc_id,q_name, NF90_DOUBLE, dim_space_id, var_q_id) )

  ! create the "G-vectors" array
  call check_netcdf( nf90_def_var(nc_id,Gvect_name, NF90_INT, &
			(/dim_space_id,dim_G1_id/), var_Gvect_id) )

  ! create the "frequency" array
  call check_netcdf( nf90_def_var(nc_id,hw_name, NF90_DOUBLE, dim_hw_id, var_hw_id) )

  ! create the real and imaginary part of the response function 

  if (nspin ==1 .or. magnon) then
     call check_netcdf( nf90_def_var(nc_id,re_chi0_name, NF90_DOUBLE, &
  			   (/dim_hw_id,dim_G1_id,dim_G2_id /),        &
						var_re_chi0_id) )

     call check_netcdf( nf90_def_var(nc_id,im_chi0_name, NF90_DOUBLE, &
			   (/dim_hw_id,dim_G1_id,dim_G2_id /),        &
						var_im_chi0_id) )

  else
     call check_netcdf( nf90_def_var(nc_id,re_chi0_name, NF90_DOUBLE,            &
     (/dim_hw_id,dim_G1_id,dim_G2_id,dim_s1_id,dim_s2_id,dim_s3_id,dim_s4_id /), &
						var_re_chi0_id) )

     call check_netcdf( nf90_def_var(nc_id,im_chi0_name, NF90_DOUBLE,            &
     (/dim_hw_id,dim_G1_id,dim_G2_id,dim_s1_id,dim_s2_id,dim_s3_id,dim_s4_id /), &
						var_im_chi0_id) )

  end if

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! frequency
  call check_netcdf( nf90_put_att(nc_id, var_hw_id, UNITS, hw_UNITS) )

  ! G vectors
  call check_netcdf( nf90_put_att(nc_id, var_Gvect_id, UNITS, space_UNITS) )

  ! q vector
  call check_netcdf( nf90_put_att(nc_id, var_q_id, UNITS, space_UNITS) )

  ! response function
  call check_netcdf( nf90_put_att(nc_id, var_im_chi0_id, UNITS, chi0_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_chi0_id, UNITS, chi0_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------

  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_q_id,(/qpt1,qpt2,qpt3/)) )

  ! write the G-vectors
  call check_netcdf( nf90_put_var(nc_id, var_Gvect_id, gvec(:,1:nG_shell_max) ))

  ! write the frequencies
  call check_netcdf( nf90_put_var(nc_id, var_hw_id, list_hbar_omega ))

  ! write the real and imaginary response functions

  if (magnon .or. nspin == 1) then

    call check_netcdf( nf90_put_var(nc_id, var_re_chi0_id,  &
			real (chi0(:,:,:,1,1,1,1))))

    call check_netcdf( nf90_put_var(nc_id, var_im_chi0_id,  &
	                aimag (chi0(:,:,:,1,1,1,1))))

  else


    call check_netcdf( nf90_put_var(nc_id, var_re_chi0_id, real(chi0)))

    call check_netcdf( nf90_put_var(nc_id, var_im_chi0_id,aimag(chi0)))

   
  end if                    

  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )

  end subroutine output_chi0_netcdf



  subroutine update_chi0_openmp\
	(nk_vec_max,nk_vec, eig_ks,eig_ksq, fk_m_fkq, widths,pw_matrix_elements)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to chi0.
  ! The details of exactly how the contribution to chi0 is calculated
  ! (eg, only imaginary part vs. real and imaginary; gaussian vs. lorentzian
  !  broadening,etc...) is left to another subroutine.
  !
  ! NOTE ABOUT TIMING:
  ! This is a resource intensive sub-routine: a large fraction of the computation
  ! time is spent here. It is tempting to time everything here with cpu_time
  ! to know exactly what takes time.  However, empirical tests show that 
  ! calling cpu_time all over the place actually MAKES THE ROUTINE MUCH SLOWER.
  ! It is thus probably much better not to do such refined timing in order to 
  ! maintain speed of execution.                                                 
  !------------------------------------------------------------------------
  use omp_lib
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max
  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)

  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)
  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)

  complex(dp)  :: pw_matrix_elements(nG_shell_max,num_wann,num_wann, &
					   npol,npol,nk_vec_max)

  complex(dp)  :: m1(npol,npol), m2(npol,npol)

  real(dp)     :: delta,  de

  ! computation variables
  integer      :: iG_1, iG_2

  integer      :: ipol1, ipol2, ipol3, ipol4

  integer      :: i_hw, i_hw_min, i_hw_max
  integer      :: k_vec_loop, iw_loop
  integer      :: nbnd1, nbnd2
  real(dp)     :: ek, ekq, f_m_f, hw
  complex(dp)  :: numerator
  complex(dp)  :: d_chi0(n_hw)

  integer      :: I_band

  logical      :: cycle_out

  complex(dp)  ::chi0_local(n_hw,nG_shell_max,nG_shell_max,npol,npol,npol,npol)


  !$omp parallel default(none)               &
  !$omp shared(nk_vec, num_wann, fk_m_fkq, eig_ks, eig_ksq, widths)   &
  !$omp shared(hw_min, dhw, list_hbar_omega, nG_shell_max)            &
  !$omp shared(pw_matrix_elements, npol,chi0,cmplx_0,n_hw)            &
  !$omp private( k_vec_loop, I_band , nbnd1, nbnd2, ek )              &
  !$omp private( ekq, de, delta, i_hw_min, i_hw_max, f_m_f )          &
  !$omp private( i_hw , d_chi0, iG_2, m2, iG_1, m1,cycle_out)         &
  !$omp private( ipol1, ipol2, ipol3, ipol4, numerator, chi0_local) 

  chi0_local(:,:,:,:,:,:,:) = cmplx_0

  ! First, build the chi0_local arrays, on each thread.
  !$omp do 
  do k_vec_loop = 1,nk_vec
     do I_band = 1,num_wann*num_wann

          nbnd1  = modulo(I_band-1,num_wann)+1
          nbnd2  = (I_band-nbnd1)/num_wann+1

          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle


          ek    =  eig_ks (nbnd1,k_vec_loop)
          ekq   =  eig_ksq(nbnd2,k_vec_loop)
          de    =  ek-ekq

          f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 

          delta = widths(nbnd1,nbnd2,k_vec_loop)


          ! compute the contribution to chi0. Let the subroutine
          ! decide what method to use!
  	  call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)

	  if (cycle_out) cycle


          do iG_2 = 1, nG_shell_max
             ! to avoid constantly multiplying by f_m_f, do it here
	       ! take transpose to account for comment above
               m2(:,:)=transpose(conjg(           &
			pw_matrix_elements(iG_2,nbnd1,nbnd2,:,:,k_vec_loop)))*f_m_f



             do iG_1 = 1, nG_shell_max
               !---------------------------------------------------- 
               ! Careful! 
               ! There is potential confusion here. The matrix 
               ! elements are defined as:
               ! 	pw_matrix_elements(nbnd1,nbnd2,ipol1,ipol2,k,G)
               !                  =
               !    < psi_{nbnd1,k}^{ipol1} | e^{-i(q+G)*r} | psi_{nbnd2,k+q}^{ipol2} >
               !
               !    and the response function 
               !
               ! 	chi_{ipol1,ipol2,ipol3,ipol4} = ...
               !    < psi_{nbnd1,k}^{ipol1}   | e^{-i(q+G)*r}  | psi_{nbnd2,k+q}^{ipol2} >
               !    < psi_{nbnd2,k+q}^{ipol3} | e^{ i(q+G')*r} | psi_{nbnd1,k}^{ipol4} >
               !                                      = ...
               !    < psi_{nbnd1,k}^{ipol1} | e^{-i(q+G)*r}  | psi_{nbnd2,k+q}^{ipol2} >
               !   [< psi_{nbnd1,k}^{ipol4} | e^{-i(q+G')*r} | psi_{nbnd2,k+q}^{ipol3} >]^*
               !
               !   So, careful when taking matrix elements!
               !   Note that npol = 1 for the magnon case and the matrix
               !   elements already account for < up | ... | down >.
               !---------------------------------------------------- 

               m1(:,:)=pw_matrix_elements(iG_1,nbnd1,nbnd2,:,:,k_vec_loop)
             
               do ipol4 = 1,npol
                do ipol3 = 1,npol
                 do ipol2 = 1,npol
                  do ipol1 = 1,npol


                     numerator = m1(ipol1,ipol2)*m2(ipol3,ipol4)


                     do i_hw = i_hw_min, i_hw_max
			   chi0_local(i_hw,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4)=     &
                     	   chi0_local(i_hw,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4)    + &
                                                d_chi0(i_hw)*numerator
		     end do ! i_hw

                  enddo
                 enddo
                enddo
               enddo

             end do 
          end do 

     end do  ! I_band
  end do  ! kvec_loop
  !$omp end do

  !$omp barrier


  ! Next, dump chi0_local into chi0;
  ! each thread should dump all its components here, so don't
  ! parallelize over the loop variables! We must still be in the
  ! "parallel" environment, however, or else chi0_local becomes
  ! undefined
  do ipol4 = 1,npol
    do ipol3 = 1,npol
      do ipol2 = 1,npol
        do ipol1 = 1,npol
          do iG_2 = 1, nG_shell_max
            do iG_1 = 1, nG_shell_max
              do i_hw = 1, n_hw
                !$omp atomic
     		chi0(i_hw,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4) =          &
     			chi0(i_hw,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4) +  &
        		chi0_local(i_hw,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4)
  	      end do ! i_hw
  	    end do ! iG_1
  	  end do ! iG_2
        end do ! ipol1
      end do ! ipol2
    end do ! ipol3
  end do ! ipol4

  !$omp end parallel

  end subroutine update_chi0_openmp


  subroutine update_chi0(nk_vec_max,nk_vec, eig_ks,eig_ksq, &
			fk_m_fkq, widths,pw_matrix_elements)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to chi0.
  ! The details of exactly how the contribution to chi0 is calculated
  ! (eg, only imaginary part vs. real and imaginary; gaussian vs. lorentzian
  !  broadening,etc...) is left to another subroutine.
  !
  ! NOTE ABOUT TIMING:
  ! This is a resource intensive sub-routine: a large fraction of the computation
  ! time is spent here. It is tempting to time everything here with cpu_time
  ! to know exactly what takes time.  However, empirical tests show that 
  ! calling cpu_time all over the place actually MAKES THE ROUTINE MUCH SLOWER.
  ! It is thus probably much better not to do such refined timing in order to 
  ! maintain speed of execution.                                                 
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max
  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)

  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)
  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)

  complex(dp)  :: pw_matrix_elements(nG_shell_max,num_wann,num_wann, &
					   npol,npol,nk_vec_max)

  complex(dp)  :: m1(npol,npol), m2(npol,npol)

  real(dp)     :: delta,  de

  ! computation variables
  integer      :: iG_1, iG_2

  integer      :: ipol1, ipol2, ipol3, ipol4

  integer      :: i_hw, i_hw_min, i_hw_max
  integer      :: k_vec_loop, iw_loop
  integer      :: nbnd1, nbnd2
  real(dp)     :: ek, ekq, f_m_f, hw
  complex(dp)  :: numerator
  complex(dp)  :: d_chi0(n_hw)

  integer      :: I_band

  logical      :: cycle_out

  ! netcdf variables
  !character(256) :: chi0_tmp_filename
  !character (*), parameter :: re_chi0_name  = "real_chi_KS"
  !character (*), parameter :: im_chi0_name  = "imaginary_chi_KS"
  !integer      :: nc_id
  !integer      :: var_re_chi0_id, var_im_chi0_id

  !chi0_tmp_filename = 'chi0_tmp.nc'
  !call check_netcdf( nf90_open(trim(chi0_tmp_filename), NF90_WRITE, nc_id) )
  !call check_netcdf( nf90_inq_varid(nc_id,re_chi0_name , var_re_chi0_id) )
  !call check_netcdf( nf90_inq_varid(nc_id,im_chi0_name , var_im_chi0_id) )

  do k_vec_loop = 1,nk_vec

     do I_band = 1,num_wann*num_wann

          nbnd1  = modulo(I_band-1,num_wann)+1
          nbnd2  = (I_band-nbnd1)/num_wann+1

          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle


          ek    =  eig_ks (nbnd1,k_vec_loop)
          ekq   =  eig_ksq(nbnd2,k_vec_loop)
          de    =  ek-ekq

          f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 

          delta = widths(nbnd1,nbnd2,k_vec_loop)


          ! compute the contribution to chi0. Let the subroutine
          ! decide what method to use!
  	  call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)

	  if (cycle_out) cycle

  	 call add_d_chi0 (d_chi0,f_m_f,k_vec_loop,nk_vec_max, &
		i_hw_min, i_hw_max,nbnd1, nbnd2, pw_matrix_elements)


        ! The subroutine below is much too slow! chi0 must be kept
        ! in RAM... 

  	! call add_d_chi0_netcdf( d_chi0,f_m_f ,k_vec_loop,nk_vec_max, &
	!		        i_hw_min, i_hw_max, nbnd1, nbnd2,      &
	!			pw_matrix_elements,                    &
	!			nc_id,var_re_chi0_id,var_im_chi0_id)


     end do  ! I_band
  end do  ! kvec_loop

  ! close netcdf file
  !call check_netcdf( nf90_close(nc_id) )

  end subroutine update_chi0





  subroutine compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
  !------------------------------------------------------------------------
  ! This subroutine dispatches the computation of d_chi0 to 
  ! the appropriate specialized subroutine corresponding to the 
  ! calculation requested by the user: gaussian smearing, lorentzian
  ! smearing, kramers-kronig or direct, etc....
  !
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  real(dp)     :: delta  ! width of gaussian
  real(dp)     :: de     ! e_k-e_k+q

  ! output variables
  complex(dp)  :: d_chi0(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max



  if ( broadening_function == 'gaussian' ) then
  	if( use_KK ) then
		call compute_dchi0_gaussian_KK &
			(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
	else
		call compute_dchi0_gaussian_direct &
			(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
	end if

  else if ( broadening_function == 'lorentzian' ) then
  	if( use_KK ) then
		call compute_dchi0_lorentzian_KK &
			(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
	else
		call compute_dchi0_lorentzian_direct &
			(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
	end if
  else if ( broadening_function == 'cold_smearing' ) then
  	if( use_KK ) then
  		call compute_dchi0_cold_smearing_KK &
			(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)

	else
		write(*,*) '*************************************'
		write(*,*) '* ERROR:                            *'
		write(*,*) '* COLD SMEARING IS ONLY IMPLEMENTED *'
		write(*,*) '* WITH KRAMERS-KRONIG.              *'
		write(*,*) '*      program stops.               *'
		write(*,*) '*************************************'
		stop
	end if



  end if

 end subroutine compute_dchi0


  subroutine compute_dchi0_gaussian_KK &
	(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
  !------------------------------------------------------------------------
  ! Compute the contribution to chi0 using gaussian width;
  ! only the imaginary part is calculated and the real part will be 
  ! obtained from Kramers-Kronig.
  !------------------------------------------------------------------------
  implicit none

  ! log_small = sqrt( - ln(10^{-16}) )
  ! thus, if |x| > log_small, e^{-x^2} < 10^{-16}, and can be neglected!
  real(dp),parameter  :: log_small = 6.0697085175405858

  ! input variables
  real(dp)     :: delta  ! width of gaussian
  real(dp)     :: de     ! e_k-e_k+q

  ! output variables
  complex(dp)  :: d_chi0(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max


  ! computation variables
  real(dp)     :: e_y2, y, y2
  complex(dp)  :: i_sqrt_pi_on_delta
  integer      :: i_hw
  real(dp)     :: const_hw1 , const_hw2 

  ! only compute the contribution to the delta-function
  ! where the function doesn't virtually vanish!
  const_hw1 = -de-hw_min
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


  i_sqrt_pi_on_delta = -cmplx_i*sqrt_pi/delta

  ! only compute the contributions to d_chi0 which
  ! are larger than machine-precision
  do i_hw = i_hw_min, i_hw_max

     y    = (list_hbar_omega(i_hw) + de)/delta
     y2   = y**2
     e_y2 = dexp(-y2)

     d_chi0(i_hw) = e_y2*i_sqrt_pi_on_delta
  end do

  end subroutine compute_dchi0_gaussian_KK 


  subroutine compute_dchi0_lorentzian_KK &
	(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
  !------------------------------------------------------------------------
  ! Compute the contribution to chi0 using lorentzian width;
  ! only the imaginary part is calculated and the real part will be 
  ! obtained from Kramers-Kronig.
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  real(dp)     :: delta  ! width of lorentzian
  real(dp)     :: de     ! e_k-e_k+q

  ! output variables
  complex(dp)  :: d_chi0(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max


  ! computation variables
  real(dp)     :: y, y2
  complex(dp)  :: i_sqrt_pi_delta
  real(dp)     :: delta2
  integer      :: i_hw


  ! don't test the bounds. It's a pain in the neck, and you'll never
  ! use lorentzians for a serious calculation anyways.
  cycle_out = .false.
  i_hw_min  = 1
  i_hw_max  = n_hw 


  i_sqrt_pi_delta = -cmplx_i*delta

  delta2 = delta**2

  ! only compute the contributions to d_chi0 which
  ! are larger than machine-precision
  do i_hw = i_hw_min, i_hw_max

     y    = list_hbar_omega(i_hw) + de
     y2   = y**2

     d_chi0(i_hw) = i_sqrt_pi_delta/(y2+delta2)
  end do

  end subroutine compute_dchi0_lorentzian_KK 



  subroutine compute_dchi0_gaussian_direct &
	(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
  !------------------------------------------------------------------------
  ! Compute the contribution to chi0 using gaussian width;
  ! both the real and imaginary parts are calculated directly.
  !------------------------------------------------------------------------
  implicit none

  ! log_small = sqrt( - ln(10^{-16}) )
  ! thus, if |x| > log_small, e^{-x^2} < 10^{-16}, and can be neglected!
  real(dp),parameter  :: log_small = 6.0697085175405858

  ! dawson function
  real(dp)     :: daw

  ! input variables
  real(dp)     :: delta  ! width of gaussian
  real(dp)     :: de     ! e_k-e_k+q


  ! output variables
  complex(dp)  :: d_chi0(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max


  ! computation variables
  real(dp)     :: e_y2, y, y2
  complex(dp)  :: i_sqrt_pi_on_delta, two_on_delta
  integer      :: i_hw
  real(dp)     :: const_hw1 , const_hw2 

  ! only compute the contribution to the delta-function
  ! where the function doesn't virtually vanish!
  const_hw1 = -de-hw_min
  const_hw2 = delta*log_small

  i_hw_min  = nint((const_hw1-const_hw2)/dhw)+1
  i_hw_max  = nint((const_hw1+const_hw2)/dhw)+1

  ! test the bounds: if the numbers are outside
  ! the range of validity, cycle out (as long as we are using KK)!

  d_chi0(:) = cmplx_0

  ! cannot cycle out because of real part!
  cycle_out = .false.

  if ( i_hw_min  < 1)    i_hw_min = 1
  if ( i_hw_max  > n_hw) i_hw_max = n_hw

  i_sqrt_pi_on_delta = -cmplx_i*sqrt_pi/delta

  ! only compute the imaginary contributions to d_chi0 which
  ! are larger than machine-precision
  do i_hw = i_hw_min, i_hw_max

     y    = (list_hbar_omega(i_hw) + de)/delta
     y2   = y**2
     e_y2 = dexp(-y2)

     d_chi0(i_hw) = e_y2*i_sqrt_pi_on_delta
  end do

  i_hw_min = 1
  i_hw_max = n_hw

  ! compute the proper analytical real part, which
  ! is related to the so-called "Dawson integral"
  ! only perform these very expensive computations
  ! if instructed NOT to use Kramers-Kronig
  two_on_delta = 2.0_dp/delta    
  do i_hw = 1, n_hw
      y  = (list_hbar_omega(i_hw) + de)/delta

      d_chi0(i_hw) = d_chi0(i_hw)+two_on_delta*daw(y)
  end do

  end subroutine compute_dchi0_gaussian_direct

  subroutine compute_dchi0_lorentzian_direct &
	(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
  !------------------------------------------------------------------------
  ! Compute the contribution to chi0 using lorentzian width;
  ! only the imaginary part is calculated and the real part will be 
  ! obtained from Kramers-Kronig.
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  real(dp)     :: delta  ! width of lorentzian
  real(dp)     :: de     ! e_k-e_k+q

  ! output variables
  complex(dp)  :: d_chi0(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max


  ! computation variables
  real(dp)     :: y
  complex(dp)  :: i_sqrt_pi_delta
  integer      :: i_hw


  ! don't test the bounds. It's a pain in the neck, and you'll never
  ! use lorentzians for a serious calculation anyways.
  cycle_out = .false.
  i_hw_min  = 1
  i_hw_max  = n_hw 

  i_sqrt_pi_delta = cmplx_i*delta

  ! only compute the contributions to d_chi0 which
  ! are larger than machine-precision
  do i_hw = i_hw_min, i_hw_max

     y    = list_hbar_omega(i_hw) + de

     d_chi0(i_hw) = cmplx_1/(y+i_sqrt_pi_delta)
  end do

  end subroutine compute_dchi0_lorentzian_direct


  subroutine create_chi0_tmp_netcdf(filename)
  !------------------------------------------------------------------------
  ! This subroutine creates a file where to temporarily store 
  ! chi0 in nectdf format, along with relevant information to 
  ! interpret the data.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  integer        :: iw_loop

  character(256) :: filename 

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
   
  integer :: var_re_chi0_id, var_im_chi0_id
  ! variables id: response function

  integer :: dim_hw_id, dim_chi0_id
  integer :: dim_G1_id, dim_G2_id
  integer :: dim_s1_id, dim_s2_id,dim_s3_id, dim_s4_id
  ! dimension id: frequency, response function, 
  !               2 G vectors dimensions, 3D space dimension
  !		  spin dimensions

  integer :: dim_hw, dim_chi0, dim_G
  ! dimensions: frequency, response function, q vector, G vectors
  !	        dimension of 3D space (3!)

  character (*), parameter :: hw_name       = "frequency"
  character (*), parameter :: G1_name       = "iG1"
  character (*), parameter :: G2_name       = "iG2"
  character (*), parameter :: re_chi0_name  = "real_chi_KS"
  character (*), parameter :: im_chi0_name  = "imaginary_chi_KS"
  character (*), parameter :: sigma1_name   = "sigma1"
  character (*), parameter :: sigma2_name   = "sigma2"
  character (*), parameter :: sigma3_name   = "sigma3"
  character (*), parameter :: sigma4_name   = "sigma4"


  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  ! header lines for spin response function. The exact length
  ! is 18(characters) x (1(hw) + 16(re) +16(im) ) =  594 
  character(594) :: header_line1, header_line2
  character(36)  :: header_bit

  integer        :: start(7), count(7)

  ! normalize at the very end, to avoid truncation error

  ! prepare the netcdf file, writing meta-data

  call check_netcdf(nf90_create(trim(filename), NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"TEMPORARY UNORMALIZED Kohn-Sham susceptibility function") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This is a storage file for a running calculation: DO NOT USE THIS"))

  !------------------------------
  ! define the dimensions 
  !------------------------------
  ! create one "G-vector" dimension, the number of G vectors
  dim_G = nG_shell_max
  call check_netcdf( nf90_def_dim(nc_id,G1_name, dim_G, dim_G1_id) )
  call check_netcdf( nf90_def_dim(nc_id,G2_name, dim_G, dim_G2_id) )

  ! create the "frequency" dimension, as well as the frequency variable
  dim_hw = n_hw
  call check_netcdf( nf90_def_dim(nc_id,hw_name, dim_hw, dim_hw_id) )

  ! create the "spin" dimensions
  call check_netcdf( nf90_def_dim(nc_id,sigma1_name, npol, dim_s1_id) )
  call check_netcdf( nf90_def_dim(nc_id,sigma2_name, npol, dim_s2_id) )
  call check_netcdf( nf90_def_dim(nc_id,sigma3_name, npol, dim_s3_id) )
  call check_netcdf( nf90_def_dim(nc_id,sigma4_name, npol, dim_s4_id) )

  !------------------------------
  ! define the variables
  !------------------------------

  call check_netcdf( nf90_def_var(nc_id,re_chi0_name, NF90_DOUBLE,            &
  (/dim_hw_id,dim_G1_id,dim_G2_id,dim_s1_id,dim_s2_id,dim_s3_id,dim_s4_id /), &
						var_re_chi0_id) )

  call check_netcdf( nf90_def_var(nc_id,im_chi0_name, NF90_DOUBLE,            &
  (/dim_hw_id,dim_G1_id,dim_G2_id,dim_s1_id,dim_s2_id,dim_s3_id,dim_s4_id /), &
						var_im_chi0_id) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Initialize the data to zero
  !--------------------------------------------------------------------
  
  !count(:) = (/n_hw,nG_shell_max,nG_shell_max,npol,npol,npol,npol/)
  !start(:) = (/1,1,1,1,1,1,1/)

  call check_netcdf( nf90_put_var(nc_id, var_re_chi0_id,  &
			real (chi0(:,:,:,1,1,1,1))))

  call check_netcdf( nf90_put_var(nc_id, var_im_chi0_id,  &
	                aimag (chi0(:,:,:,1,1,1,1))))

  !--------------------------------------------------------------------
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_close(nc_id) )

  end subroutine create_chi0_tmp_netcdf



  subroutine add_d_chi0 ( d_chi0,f_m_f ,k_vec_loop,nk_vec_max, &
			i_hw_min, i_hw_max,nbnd1, nbnd2, pw_matrix_elements)
  !------------------------------------------------------------------------
  ! This subroutine adds the contribution of d_chi0 to chi0.
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  complex(dp)  :: d_chi0(n_hw)
  integer      :: k_vec_loop, nk_vec_max
  integer      :: nbnd1, nbnd2
  integer      :: i_hw_min, i_hw_max
  real(dp)     :: f_m_f
  complex(dp)  :: pw_matrix_elements(nG_shell_max,num_wann,num_wann, &
					   npol,npol,nk_vec_max)

  ! local variables
  complex(dp)  :: m1(npol,npol), m2(npol,npol)
  integer      :: iG_1, iG_2
  integer      :: ipol1, ipol2, ipol3, ipol4
  integer      :: i_hw
  complex(dp)  :: numerator



  do iG_2 = 1, nG_shell_max
     ! to avoid constantly multiplying by f_m_f, do it here
     ! take transpose to account for comment above
     m2(:,:)=transpose(conjg(           &
	     pw_matrix_elements(iG_2,nbnd1,nbnd2,:,:,k_vec_loop)))*f_m_f

     do iG_1 = 1, nG_shell_max
        !---------------------------------------------------- 
        ! Careful! 
        ! There is potential confusion here. The matrix 
        ! elements are defined as:
        ! 	pw_matrix_elements(nbnd1,nbnd2,ipol1,ipol2,k,G)
        !                  =
        !    < psi_{nbnd1,k}^{ipol1} | e^{-i(q+G)*r} | psi_{nbnd2,k+q}^{ipol2} >
        !
        !    and the response function 
        !
        ! 	chi_{ipol1,ipol2,ipol3,ipol4} = ...
        !    < psi_{nbnd1,k}^{ipol1}   | e^{-i(q+G)*r}  | psi_{nbnd2,k+q}^{ipol2} >
        !    < psi_{nbnd2,k+q}^{ipol3} | e^{ i(q+G')*r} | psi_{nbnd1,k}^{ipol4} >
        !                                      = ...
        !    < psi_{nbnd1,k}^{ipol1} | e^{-i(q+G)*r}  | psi_{nbnd2,k+q}^{ipol2} >
        !   [< psi_{nbnd1,k}^{ipol4} | e^{-i(q+G')*r} | psi_{nbnd2,k+q}^{ipol3} >]^*
        !
        !   So, careful when taking matrix elements!
        !   Note that npol = 1 for the magnon case and the matrix
        !   elements already account for < up | ... | down >.
        !---------------------------------------------------- 
        m1(:,:)=pw_matrix_elements(iG_1,nbnd1,nbnd2,:,:,k_vec_loop)
             
        do ipol4 = 1,npol
          do ipol3 = 1,npol
            do ipol2 = 1,npol
              do ipol1 = 1,npol

                 numerator = m1(ipol1,ipol2)*m2(ipol3,ipol4)

                 do i_hw = i_hw_min, i_hw_max
		    chi0(i_hw,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4)=     &
                    chi0(i_hw,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4)    + &
                                            d_chi0(i_hw)*numerator
		 end do ! i_hw

              enddo
            enddo
          enddo
        enddo

     end do 
   end do 

  end subroutine add_d_chi0


  subroutine add_d_chi0_netcdf( d_chi0,f_m_f ,k_vec_loop,nk_vec_max, &
			        i_hw_min, i_hw_max, nbnd1, nbnd2,    &
				pw_matrix_elements,                  &
				nc_id,var_re_chi0_id,var_im_chi0_id)
  !------------------------------------------------------------------------
  ! This subroutine adds the contribution of d_chi0 to chi0, stored
  ! in a netcdf file.
  !------------------------------------------------------------------------
  use netcdf
  implicit none
  ! input variables
  complex(dp)  :: d_chi0(n_hw)
  integer      :: k_vec_loop, nk_vec_max
  integer      :: nbnd1, nbnd2
  integer      :: i_hw_min, i_hw_max
  real(dp)     :: f_m_f
  complex(dp)  :: pw_matrix_elements(nG_shell_max,num_wann,num_wann, &
					   npol,npol,nk_vec_max)
  integer      :: nc_id
  integer      :: var_re_chi0_id, var_im_chi0_id


  ! local variables
  complex(dp)  :: m1(npol,npol), m2(npol,npol)
  integer      :: iG_1, iG_2
  integer      :: ipol1, ipol2, ipol3, ipol4
  integer      :: i_hw
  complex(dp)  :: numerator

  real(dp)     :: re_data(i_hw_max-i_hw_min+1)
  real(dp)     :: im_data(i_hw_max-i_hw_min+1)

  ! netcdf variables

  integer      :: start(7),count(7)


  do iG_2 = 1, nG_shell_max
     ! to avoid constantly multiplying by f_m_f, do it here
     ! take transpose to account for comment above
     m2(:,:)=transpose(conjg(           &
	     pw_matrix_elements(iG_2,nbnd1,nbnd2,:,:,k_vec_loop)))*f_m_f

     do iG_1 = 1, nG_shell_max
        !---------------------------------------------------- 
        ! Careful! 
        ! There is potential confusion here. The matrix 
        ! elements are defined as:
        ! 	pw_matrix_elements(nbnd1,nbnd2,ipol1,ipol2,k,G)
        !                  =
        !    < psi_{nbnd1,k}^{ipol1} | e^{-i(q+G)*r} | psi_{nbnd2,k+q}^{ipol2} >
        
        !    and the response function 
        !
        ! 	chi_{ipol1,ipol2,ipol3,ipol4} = ...
        !    < psi_{nbnd1,k}^{ipol1}   | e^{-i(q+G)*r}  | psi_{nbnd2,k+q}^{ipol2} >
        !    < psi_{nbnd2,k+q}^{ipol3} | e^{ i(q+G')*r} | psi_{nbnd1,k}^{ipol4} >
        !                                      = ...
        !    < psi_{nbnd1,k}^{ipol1} | e^{-i(q+G)*r}  | psi_{nbnd2,k+q}^{ipol2} >
        !   [< psi_{nbnd1,k}^{ipol4} | e^{-i(q+G')*r} | psi_{nbnd2,k+q}^{ipol3} >]^*
        !
        !   So, careful when taking matrix elements!
        !   Note that npol = 1 for the magnon case and the matrix
        !   elements already account for < up | ... | down >.
        !---------------------------------------------------- 
        m1(:,:)=pw_matrix_elements(iG_1,nbnd1,nbnd2,:,:,k_vec_loop)
             
        do ipol4 = 1,npol
          do ipol3 = 1,npol
            do ipol2 = 1,npol
              do ipol1 = 1,npol

                 numerator = m1(ipol1,ipol2)*m2(ipol3,ipol4)

		 ! read the data from netcdf file
     		 start(:) = (/i_hw_min,iG_1,iG_2,ipol1,ipol2,ipol3,ipol4/)
  		 count(:) = (/i_hw_max-i_hw_min+1,1,1,1,1,1,1/)

  		 call check_netcdf( nf90_get_var(nc_id, var_re_chi0_id, &
						re_data, start,count))
  		 call check_netcdf( nf90_get_var(nc_id, var_im_chi0_id, &
						im_data, start,count))

		 ! compute the change
                 do i_hw = i_hw_min, i_hw_max

		    re_data(i_hw-i_hw_min+1) =               &
		    		re_data(i_hw-i_hw_min+1) +   &
				real(d_chi0(i_hw)*numerator)

		    im_data(i_hw-i_hw_min+1) =               &
		    		im_data(i_hw-i_hw_min+1) +   &
				aimag(d_chi0(i_hw)*numerator)

		 end do ! i_hw

		 ! write to file
  		  call check_netcdf( nf90_put_var(nc_id, var_re_chi0_id, &
						re_data, start,count))
  		  call check_netcdf( nf90_put_var(nc_id, var_im_chi0_id, &
						im_data, start,count))

              enddo
            enddo
          enddo
        enddo

     end do 
   end do 

  end subroutine add_d_chi0_netcdf

  subroutine read_chi0_tmp_netcdf(chi0_tmp_filename,re_chi0_tmp,im_chi0_tmp)
  !------------------------------------------------------------------------
  ! This subroutine creates a file where to temporarily store 
  ! chi0 in nectdf format, along with relevant information to 
  ! interpret the data.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  integer        :: iw_loop

  character(256) :: chi0_tmp_filename 

  real(dp)       :: re_chi0_tmp(n_hw,nG_shell_max,nG_shell_max,npol,npol,npol,npol)
  real(dp)       :: im_chi0_tmp(n_hw,nG_shell_max,nG_shell_max,npol,npol,npol,npol)

  ! netcdf variables

  integer :: nc_id        
  ! netcdf file id
   
  integer :: var_re_chi0_id, var_im_chi0_id
  ! variables id: response function

  character (*), parameter :: re_chi0_name  = "real_chi_KS"
  character (*), parameter :: im_chi0_name  = "imaginary_chi_KS"


  call check_netcdf( nf90_open(trim(chi0_tmp_filename), NF90_NOWRITE, nc_id) )

  call check_netcdf( nf90_inq_varid(nc_id,re_chi0_name , var_re_chi0_id) )
  call check_netcdf( nf90_inq_varid(nc_id,im_chi0_name , var_im_chi0_id) )



   call check_netcdf( nf90_get_var(nc_id, var_re_chi0_id, re_chi0_tmp))
   call check_netcdf( nf90_get_var(nc_id, var_im_chi0_id, im_chi0_tmp))


  call check_netcdf( nf90_close(nc_id) )

  end subroutine read_chi0_tmp_netcdf

!
!  subroutine update_chi0_intw_netcdf(nc_id, var_re_chi0_id, var_im_chi0_id , &
!		nk_vec_max,nk_vec, eig_ks,eig_ksq, u_ks,u_ksq,fk_m_fkq, widths)
!  !------------------------------------------------------------------------
!  ! This subroutine adds a sub-block contribution to chi0.
!  ! The details of exactly how the contribution to chi0 is calculated
!  ! (eg, only imaginary part vs. real and imaginary; gaussian vs. lorentzian
!  !  broadening,etc...) is left to another subroutine.
!  !------------------------------------------------------------------------
!  use omp_lib
!  use netcdf
!  implicit none
!
!  ! input variables
!  integer      :: nc_id, var_re_chi0_id, var_im_chi0_id 
!  integer      :: nk_vec, nk_vec_max
!
!  real(dp)     :: eig_ks (num_wann,nk_vec_max)
!  real(dp)     :: eig_ksq(num_wann,nk_vec_max)
!
!  complex(dp)  :: u_ks (num_wann,num_wann,nk_vec_max)
!  complex(dp)  :: u_ksq(num_wann,num_wann,nk_vec_max)
!
!  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)
!  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)
!
!
!  ! computation variables
!  integer      :: start(3), count(3)
!  integer      :: Iproj1, Iproj2
!  integer      :: m1, m2, m3, m4
!  real(dp)     :: delta,  de
!
!
!  integer      :: i_hw, i_hw_min, i_hw_max
!  integer      :: k_vec_loop, iw_loop
!  integer      :: nbnd1, nbnd2
!  real(dp)     :: ek, ekq, f_m_f, hw
!  complex(dp)  :: numerator
!  complex(dp)  :: d_chi0(n_hw), partial_chi0(n_hw)
!  real(dp)     :: real_partial_chi0(n_hw)
!
!
!
!  complex(dp),allocatable  :: matrix_data(:,:,:)
!!  complex(dp)  :: matrix_data(n_hw)
!
!  complex(dp),allocatable  :: b_projectors(:,:,:)
!
!  complex(dp)  :: b_projector, b_I1_projector, b_I2_projector
!
!  complex(dp)  :: work_u1(num_wann), work_u2(num_wann), work_V(num_wann*num_wann)
!
!  integer      :: I_band, I_band_loop, I_band_max, I_band_local_max
!  integer      :: I_m1m2
!
!  integer,allocatable  :: list_nbnd1(:,:), list_nbnd2(:,:) 
!  integer,allocatable  :: list_i_hw_min(:,:), list_i_hw_max(:,:) 
!
!  logical      :: cycle_out
!
!
!  real(dp)     :: time1, time2
!
!  !--------------------------------------------------------------------------
!  ! First, build an array that will contain
!  !		[f_{n1 k}-f_{n2 k+q}]* delta(hw+ E_{n1 k} -E_{n2 k+q})
!  ! for k, (n1n2), hw. In order to maximize RAM usage, pack the 
!  ! (n1n2) indices such that only terms where [f_{n1 k}-f_{n2 k+q}] /= 0
!  ! appear.
!  !--------------------------------------------------------------------------
!
!  ! First, probe the fk_m_fkq array to see how large the (n1n2) dimension 
!  ! must be
!
!  I_band_max = 0
!
!  write(142,5) '#==========================================================='
!  write(142,5) '#                          new block                        '
!  write(142,5) '#-----------------------------------------------------------'
!  
!
!  write(142,5) '#                                                           '
!  write(142,5) '#  Tabulating bands and frequency indices ...               '
!  call get_timing(time1)
!  do k_vec_loop = 1,nk_vec
!
!     ! find the number of components (n1n2) which are not zero, for a given kpoint
!     I_band_local_max = 0
!
!     do I_band = 1,num_wann*num_wann
!
!          nbnd1  = modulo(I_band-1,num_wann)+1
!          nbnd2  = (I_band-nbnd1)/num_wann+1
!
!          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle
!          I_band_local_max = I_band_local_max + 1
!
!     end do
!     
!     ! find maximum for all kpoints
!     if (I_band_local_max > I_band_max) I_band_max = I_band_local_max 
!  end do
!
!
!  ! next, allocate an array that will hold the data, and tabulate
!  ! the band indices
!  allocate(list_nbnd1(I_band_max ,nk_vec))
!  allocate(list_nbnd2(I_band_max ,nk_vec))
!
!  allocate(list_i_hw_min(I_band_max ,nk_vec))
!  allocate(list_i_hw_max(I_band_max ,nk_vec))
!  call get_timing(time2)
!  write(142,10) '#         time :',time2-time1,' seconds '
!  write(142, 5) '#-----------------------------------------------------------'
!
!
!  ! This can be very large! Make sure it fits in RAM.
!  !allocate(matrix_data(n_hw,I_band_max ,nk_vec))
!
!  ! initialize arrays
!
!  !matrix_data = cmplx_0
!
!  list_nbnd1    = 0
!  list_nbnd2    = 0
!  list_i_hw_min = 0
!  list_i_hw_max = 0
!
!  write(142,5) '#                                                           '
!  write(142,5) '#  Building frequency dependent data array ...              '
!  call get_timing(time1)
!  do k_vec_loop = 1,nk_vec
!     !write(142,*) '# k_vec_loop = ',k_vec_loop,' of ',nk_vec
!
!     I_band_loop = 0
!
!     do I_band = 1,num_wann*num_wann
!
!          nbnd1  = modulo(I_band-1,num_wann)+1
!          nbnd2  = (I_band-nbnd1)/num_wann+1
!
!          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle
!
!          I_band_loop  = I_band_loop + 1
!        
!  	  list_nbnd1(I_band_loop ,k_vec_loop) =  nbnd1
!          list_nbnd2(I_band_loop ,k_vec_loop) =  nbnd2
!
!
!          ek    =  eig_ks (nbnd1,k_vec_loop)
!          ekq   =  eig_ksq(nbnd2,k_vec_loop)
!          de    =  ek-ekq
!
!          f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 
!
!          delta = widths(nbnd1,nbnd2,k_vec_loop)
!
!  	  call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
!
!
!
!          ! if cycle_out, there is no fequency contribution. Do not
!          ! update indices, they are already set to zero 
!	  if (cycle_out) cycle
!
!  	  list_i_hw_min(I_band_loop ,k_vec_loop) = i_hw_min
!  	  list_i_hw_max(I_band_loop ,k_vec_loop) = i_hw_max
!
!          !write(*,'(A,I4,A,I4,A,2I4)') &
!	  !	      ' k_vec_loop  = ',k_vec_loop,      &
!          !            ' I_band_loop = ', I_band_loop,  &  
!	  !	      ' imin,imax   = ',i_hw_min,i_hw_max 
!
!          matrix_data(i_hw_min:i_hw_max,I_band_loop ,k_vec_loop)  =  &
!				f_m_f*d_chi0(i_hw_min:i_hw_max)
!
!     end do
!  end do
!  call get_timing(time2)
!  write(142,10) '#         time :',time2-time1,' seconds '
!  write(142, 5) '#-----------------------------------------------------------'
!
!  !--------------------------------------------------------------------------
!  ! Next, compute the projector functions to project onto the 
!  ! projected Wannier basis.
!  !--------------------------------------------------------------------------
!  write(142,5) '#                                                           '
!  write(142,5) '#  Building projector arrays ...                            '
!
!  call get_timing(time1)
!  allocate(b_projectors(basis_size_intw,I_band_max ,nk_vec))
!
!  b_projectors(:,:,:) = cmplx_0
!
!  do k_vec_loop = 1,nk_vec
!     do I_band_loop = 1,I_band_max
!
!  	nbnd1 = list_nbnd1(I_band_loop ,k_vec_loop)
!        nbnd2 = list_nbnd2(I_band_loop ,k_vec_loop)
!
!        i_hw_min = list_i_hw_min(I_band_loop ,k_vec_loop)  
!	!i_hw_max = list_i_hw_max(I_band_loop ,k_vec_loop) 
!
!        ! don't bother to calculate if there is no contribution!
!	if (nbnd1 == 0 .or. i_hw_min == 0 ) exit
!
!        do Iproj1 = 1, basis_size_intw
!
!           b_projector = cmplx_0
!
!           work_V (:) = conjg(V_intw_bar(:,Iproj1))
!	   work_u1(:) = u_ksq(nbnd2,:,k_vec_loop)
!	   work_u2(:) = conjg(u_ks (nbnd1,:,k_vec_loop))
!
!
!           do I_m1m2 = 1, num_wann**2   
!
!              m1 = list_m1_intw(I_m1m2)
!              m2 = list_m2_intw(I_m1m2)
!
!	      b_projector = b_projector   +        & 
!			work_V(I_m1m2)*work_u1(m1)*work_u2(m2)
!
!!	      b_projector = b_projector   +        & 
!!			conjg(V_intw_bar(I_m1m2,Iproj1))     &
!!		      *(      u_ksq(nbnd2,m1,k_vec_loop)     &
!!                       *conjg(u_ks (nbnd1,m2,k_vec_loop)))
!
!           end do ! I_m1m2
!        
!           b_projectors(Iproj1,I_band_loop,k_vec_loop) = b_projector 
!
!       end do ! Iproj1
!
!     end do ! I_band_loop
!  end do ! k_vec_loop
!  call get_timing(time2)
!  write(142,10) '#         time :',time2-time1,' seconds '
!  write(142, 5) '#-----------------------------------------------------------'
!
!
!  !--------------------------------------------------------------------------
!  ! Next, compute the contribution  to each combination {Iproj1,Iproj2}
!  ! of the projected Wannier basis.
!  !--------------------------------------------------------------------------
!
!  count(:) = (/n_hw,1,1/)
!
!  write(142,5) '#                                                           '
!  write(142,5) '#  read, update, write chi0 array to file...                '
!
!  call get_timing(time1)
!  do Iproj2 = 1,basis_size_intw
!        do Iproj1 = 1,basis_size_intw
!  
!        write(142,15) '#  Iproj1 = ',Iproj1,' , Iproj2 = ',Iproj2
!
!           partial_chi0(:) = cmplx_0
!
!           do k_vec_loop = 1,nk_vec
!              do I_band_loop = 1,I_band_max
!
!
!  	         nbnd1 = list_nbnd1(I_band_loop ,k_vec_loop)
!                 nbnd2 = list_nbnd2(I_band_loop ,k_vec_loop)
!
!  	         i_hw_min = list_i_hw_min(I_band_loop ,k_vec_loop)
!  	         i_hw_max = list_i_hw_max(I_band_loop ,k_vec_loop)
!
!                 ! if there is no contribution, move on to next iteration
!		 if (nbnd1 == 0 ) exit
!		 if (i_hw_min == 0) cycle
!              
!                 ! get projectors to go to the projected basis
!
!		 b_I1_projector =  b_projectors(Iproj1,I_band_loop,k_vec_loop)
!		 b_I2_projector =  b_projectors(Iproj2,I_band_loop,k_vec_loop)
!
!
!
!                 partial_chi0(i_hw_min:i_hw_max) =                        &
!                   partial_chi0(i_hw_min:i_hw_max) +                      &
!		   matrix_data(i_hw_min:i_hw_max,I_band_loop ,k_vec_loop) &
!		   *conjg(b_I1_projector)*b_I2_projector 
!
!
!              end do !I_band_loop
!           end do !k_vec_loop
!
!
!           ! Write to netcdf file!
!
!
!           start(:) = (/1,Iproj1,Iproj2/)
!	   call check_netcdf( nf90_get_var(nc_id,var_re_chi0_id, &
!					real_partial_chi0(:), start, count) ) 
!
!           real_partial_chi0(:) = real_partial_chi0(:) + real(partial_chi0(:))
!
!	   call check_netcdf( nf90_put_var(nc_id,var_re_chi0_id, &
!					real_partial_chi0(:), start, count) ) 
!
!	   call check_netcdf( nf90_get_var(nc_id,var_im_chi0_id, &
!					real_partial_chi0(:), start, count) ) 
!
!           real_partial_chi0(:) = real_partial_chi0(:) + aimag(partial_chi0(:))
!
!	   call check_netcdf( nf90_put_var(nc_id,var_im_chi0_id, &
!					real_partial_chi0(:), start, count) ) 
!
!           !call check_netcdf( nf90_sync(nc_id) ) 
!           call get_timing(time2)
!
!        end do ! Iproj1 
!      end do !Iproj2
!      call get_timing(time2)
!      write(142,10) '#         time :',time2-time1,' seconds '
!      write(142, 5) '#==========================================================='
!
!  5   format(A)
!  10  format(A,F8.4,A)
!  15  format(A,I4,A,I4)
!
!  end subroutine update_chi0_intw_netcdf


  subroutine normalize_chi0_intw_netcdf(nc_id,var_re_chi0_id,var_im_chi0_id)
  !------------------------------------------------------------------------
  ! This subroutine normalizes the response function in the netcdf file 
  ! and closes the file.
  !------------------------------------------------------------------------
  use netcdf
  implicit none


  integer        :: nc_id, var_re_chi0_id, var_im_chi0_id

  real(dp)       :: volume,  rho_0

  integer        :: Iproj1, Iproj2

  integer        :: start(3), count(3)

  real(dp)       :: partial_chi0(n_hw)

  !------------------------------------------------------------------

  volume = alat**3*abs(                                &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))

  rho_0     = 1.0_dp/volume
  

  !--------------------------------------------------------------------
  ! read data, normalize and write back
  !--------------------------------------------------------------------

  count(:) = (/n_hw,1,1/)
  do Iproj1 = 1,basis_size_intw
     do Iproj2 = 1,basis_size_intw

           start(:) = (/1,Iproj1,Iproj2/)

           ! normalize the real part
	   call check_netcdf( nf90_get_var(nc_id,var_re_chi0_id, &
					partial_chi0(:), start, count) ) 

           partial_chi0(:) =  partial_chi0(:)*rho_0

	   call check_netcdf( nf90_put_var(nc_id,var_re_chi0_id, &
					partial_chi0(:), start, count) ) 

           ! normalize the imaginary part
	   call check_netcdf( nf90_get_var(nc_id,var_im_chi0_id, &
					partial_chi0(:), start, count) ) 

           partial_chi0(:) =  partial_chi0(:)*rho_0

	   call check_netcdf( nf90_put_var(nc_id,var_im_chi0_id, &
					partial_chi0(:), start, count) ) 

    end do 
  end do 


  end subroutine normalize_chi0_intw_netcdf


  subroutine open_chi0_intw_netcdf(nc_id,var_re_chi0_id,var_im_chi0_id)
  !------------------------------------------------------------------------
  ! This subroutine simply opens the netcdf file, assuming it exists,
  ! and returns the necessary IDs.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  character(256) :: filename

  character (*), parameter :: re_chi0_name  = "real_chi_KS"
  character (*), parameter :: im_chi0_name  = "imaginary_chi_KS"

  ! netcdf variables
  integer :: nc_id        
  ! netcdf file id

  integer :: var_re_chi0_id,var_im_chi0_id

  filename = trim(chi0_file_prefix)//'.nc'

  call check_netcdf(nf90_open(filename, NF90_WRITE, nc_id))

  call check_netcdf(nf90_inq_varid(nc_id, re_chi0_name ,var_re_chi0_id))
  call check_netcdf(nf90_inq_varid(nc_id, im_chi0_name ,var_im_chi0_id))



  end subroutine open_chi0_intw_netcdf


  subroutine create_chi0_intw_bin(chi0_io_unit)
  !------------------------------------------------------------------------
  ! This subroutine creates and initializes a fortran .bin file 
  ! which will contain the chi0 array.  The purpose of this file
  ! is to store the very large dataset without having it stay in RAM!
  ! The use of a fortran native format instead of netcdf is to improve
  ! io performance.
  !------------------------------------------------------------------------
  implicit none

  ! input-output
  integer :: chi0_io_unit


  ! computation variable
  character (*), parameter :: filename = "chi_KS_tmp.bin"

  integer   :: Iproj1, Iproj2
  integer   :: record_index, record_length 

  complex(dp) :: chi0_tmp(n_hw)

  ! find a free unit
  chi0_io_unit = find_free_unit()


  ! open a temporary file which will hold the total chi0 array
  record_length = n_hw*direct_io_factor_cmplx

  open(unit = chi0_io_unit , file = filename, form = 'unformatted', &
     status = 'scratch', access = 'direct', recl = record_length,   &
     action = 'readwrite')

  ! initialize data to zero. Make sure the writing is as sequential
  ! as possible !

  chi0_tmp(:)   = cmplx_0
  record_index = 0

  do Iproj2 = 1,basis_size_intw
     do Iproj1 = 1,basis_size_intw
	   record_index = record_index + 1

	   write(chi0_io_unit ,rec=record_index) chi0_tmp(:)

    end do 
  end do 

  ! leave the file open so that other subroutines can write to it! 
  !call check_netcdf( nf90_close(nc_id) )

  end subroutine create_chi0_intw_bin
!
!  subroutine update_chi0_intw_bin(chi0_io_unit,  &
!		nk_vec_max,nk_vec, eig_ks,eig_ksq, u_ks,u_ksq,fk_m_fkq, widths)
!  !------------------------------------------------------------------------
!  ! This subroutine adds a sub-block contribution to chi0.
!  ! The details of exactly how the contribution to chi0 is calculated
!  ! (eg, only imaginary part vs. real and imaginary; gaussian vs. lorentzian
!  !  broadening,etc...) is left to another subroutine.
!  ! The data is written to a binary file which is assumed to be already
!  ! opened.
!  !------------------------------------------------------------------------
!  use omp_lib
!  implicit none
!
!  ! input variables
!  integer      :: chi0_io_unit
!  integer      :: nk_vec, nk_vec_max
!
!  real(dp)     :: eig_ks (num_wann,nk_vec_max)
!  real(dp)     :: eig_ksq(num_wann,nk_vec_max)
!
!  complex(dp)  :: u_ks (num_wann,num_wann,nk_vec_max)
!  complex(dp)  :: u_ksq(num_wann,num_wann,nk_vec_max)
!
!  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)
!  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)
!
!
!  ! computation variables
!
!  integer      :: record_index_chi0
!
!  integer      :: Iproj1, Iproj2
!  integer      :: m1, m2, m3, m4
!  real(dp)     :: delta,  de
!
!
!  integer      :: i_hw, i_hw_min, i_hw_max
!  integer      :: k_vec_loop, iw_loop
!  integer      :: nbnd1, nbnd2
!  real(dp)     :: ek, ekq, f_m_f, hw
!
!  real(dp)     :: d_i_hw_avg
!
!  complex(dp)  :: chi0_tmp(n_hw)
!  complex(dp)  :: d_chi0(n_hw), partial_chi0(n_hw)
!
!
!  ! read-write the array "matrix_data" from a bin file
!  integer ::  matrix_data_io_unit, record_length_M, record_index_M
!
!!  complex(dp),allocatable  :: matrix_data(:,:,:)
!  complex(dp)  :: matrix_data(n_hw)
!
!  complex(dp),allocatable  :: b_projectors(:,:,:)
!
!  complex(dp)  :: b_projector, b_I1_projector, b_I2_projector
!
!  complex(dp)  :: work_u1(num_wann), work_u2(num_wann), work_V(num_wann*num_wann)
!
!  integer      :: I_band, I_band_loop, I_band_max, I_band_local_max
!  integer      :: I_m1m2
!
!  integer,allocatable  :: list_nbnd1(:,:), list_nbnd2(:,:) 
!  integer,allocatable  :: list_i_hw_min(:,:), list_i_hw_max(:,:) 
!
!  logical      :: cycle_out
!
!
!  real(dp)     :: time1, time2
!  real(dp)     :: read_time1, read_time2, read_time
!  real(dp)     :: c_time1, c_time2, c_time
!
!  !--------------------------------------------------------------------------
!  ! First, build an array that will contain
!  !		[f_{n1 k}-f_{n2 k+q}]* delta(hw+ E_{n1 k} -E_{n2 k+q})
!  ! for k, (n1n2), hw. In order to maximize RAM usage, pack the 
!  ! (n1n2) indices such that only terms where [f_{n1 k}-f_{n2 k+q}] /= 0
!  ! appear.
!  !--------------------------------------------------------------------------
!
!  ! First, probe the fk_m_fkq array to see how large the (n1n2) dimension 
!  ! must be
!
!  I_band_max = 0
!
!  write(142,5) '#==========================================================='
!  write(142,5) '#                          new block                        '
!  write(142,5) '#-----------------------------------------------------------'
!  
!
!  write(142,5) '#                                                           '
!  write(142,5) '#  Tabulating bands and frequency indices ...               '
!  call get_timing(time1)
!  do k_vec_loop = 1,nk_vec
!
!     ! find the number of components (n1n2) which are not zero, for a given kpoint
!     I_band_local_max = 0
!
!     do I_band = 1,num_wann*num_wann
!
!          nbnd1  = modulo(I_band-1,num_wann)+1
!          nbnd2  = (I_band-nbnd1)/num_wann+1
!
!          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle
!          I_band_local_max = I_band_local_max + 1
!
!     end do
!     
!     ! find maximum for all kpoints
!     if (I_band_local_max > I_band_max) I_band_max = I_band_local_max 
!  end do
!
!
!  ! next, allocate an array that will hold the data, and tabulate
!  ! the band indices
!  allocate(list_nbnd1(I_band_max ,nk_vec))
!  allocate(list_nbnd2(I_band_max ,nk_vec))
!
!  allocate(list_i_hw_min(I_band_max ,nk_vec))
!  allocate(list_i_hw_max(I_band_max ,nk_vec))
!  call get_timing(time2)
!  write(142,10) '#         time :',time2-time1,' seconds '
!  write(142, 5) '#-----------------------------------------------------------'
!
!
!  ! This can be very large! Make sure it fits in RAM.
!  !allocate(matrix_data(n_hw,I_band_max ,nk_vec))
!
!  ! initialize arrays
!
!  matrix_data = cmplx_0
!
!  list_nbnd1    = 0
!  list_nbnd2    = 0
!  list_i_hw_min = 0
!  list_i_hw_max = 0
!
!  write(142,5) '#                                                           '
!  write(142,5) '#  Building frequency dependent data array ...              '
!
!
!  ! this can be a very sparse array! Only store what you need!
!  matrix_data_io_unit = find_free_unit()
!  record_length_M     = direct_io_factor_cmplx
!
!  open(unit = matrix_data_io_unit ,  form = 'unformatted',            &
!     status = 'scratch', access = 'direct', recl = record_length_M,   &
!     action = 'readwrite')
!
!  record_index_M = 0
!  call get_timing(time1)
!  do k_vec_loop = 1,nk_vec
!     !write(142,*) '# k_vec_loop = ',k_vec_loop,' of ',nk_vec
!
!     I_band_loop = 0
!
!     do I_band = 1,num_wann*num_wann
!
!          nbnd1  = modulo(I_band-1,num_wann)+1
!          nbnd2  = (I_band-nbnd1)/num_wann+1
!
!          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle
!
!          I_band_loop  = I_band_loop + 1
!        
!  	  list_nbnd1(I_band_loop ,k_vec_loop) =  nbnd1
!          list_nbnd2(I_band_loop ,k_vec_loop) =  nbnd2
!
!
!          ek    =  eig_ks (nbnd1,k_vec_loop)
!          ekq   =  eig_ksq(nbnd2,k_vec_loop)
!          de    =  ek-ekq
!
!          f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 
!
!          delta = widths(nbnd1,nbnd2,k_vec_loop)
!
!  	  call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
!
!
!
!          ! if cycle_out, there is no fequency contribution. Do not
!          ! update indices, they are already set to zero 
!	  if (cycle_out) cycle
!
!  	  list_i_hw_min(I_band_loop ,k_vec_loop) = i_hw_min
!  	  list_i_hw_max(I_band_loop ,k_vec_loop) = i_hw_max
!
!          !write(*,'(A,I4,A,I4,A,2I4)') &
!	  !	      ' k_vec_loop  = ',k_vec_loop,      &
!          !            ' I_band_loop = ', I_band_loop,  &  
!	  !	      ' imin,imax   = ',i_hw_min,i_hw_max 
!
!
!  	  matrix_data  = cmplx_0
!
!          matrix_data(i_hw_min:i_hw_max)  =  &
!	  			f_m_f*d_chi0(i_hw_min:i_hw_max)
!
!          do i_hw = i_hw_min, i_hw_max 
!             	record_index_M = record_index_M + 1 
!	  	write(matrix_data_io_unit ,rec=record_index_M) matrix_data(i_hw)
!          end do
!
!          !matrix_data(i_hw_min:i_hw_max,I_band_loop ,k_vec_loop)  =  &
!	  !			f_m_f*d_chi0(i_hw_min:i_hw_max)
!
!     end do
!  end do
!  call get_timing(time2)
!  write(142,10) '#         time :',time2-time1,' seconds '
!  write(142, 5) '#-----------------------------------------------------------'
!
!  !--------------------------------------------------------------------------
!  ! Next, compute the projector functions to project onto the 
!  ! projected Wannier basis.
!  !--------------------------------------------------------------------------
!  write(142,5) '#                                                           '
!  write(142,5) '#  Building projector arrays ...                            '
!
!  call get_timing(time1)
!  !allocate(b_projectors(basis_size_intw,I_band_max ,nk_vec))
!  allocate(b_projectors(I_band_max,nk_vec, basis_size_intw))
!
!  b_projectors(:,:,:) = cmplx_0
!
!  do k_vec_loop = 1,nk_vec
!     do I_band_loop = 1,I_band_max
!
!  	nbnd1 = list_nbnd1(I_band_loop ,k_vec_loop)
!        nbnd2 = list_nbnd2(I_band_loop ,k_vec_loop)
!
!        i_hw_min = list_i_hw_min(I_band_loop ,k_vec_loop)  
!	!i_hw_max = list_i_hw_max(I_band_loop ,k_vec_loop) 
!
!        ! don't bother to calculate if there is no contribution!
!	if (nbnd1 == 0 .or. i_hw_min == 0 ) exit
!
!        do Iproj1 = 1, basis_size_intw
!
!           b_projector = cmplx_0
!
!           work_V (:) = conjg(V_intw_bar(:,Iproj1))
!	   work_u1(:) = u_ksq(nbnd2,:,k_vec_loop)
!	   work_u2(:) = conjg(u_ks (nbnd1,:,k_vec_loop))
!
!
!           do I_m1m2 = 1, num_wann**2   
!
!              m1 = list_m1_intw(I_m1m2)
!              m2 = list_m2_intw(I_m1m2)
!
!	      b_projector = b_projector   +        & 
!			work_V(I_m1m2)*work_u1(m1)*work_u2(m2)
!
!!	      b_projector = b_projector   +        & 
!!			conjg(V_intw_bar(I_m1m2,Iproj1))     &
!!		      *(      u_ksq(nbnd2,m1,k_vec_loop)     &
!!                       *conjg(u_ks (nbnd1,m2,k_vec_loop)))
!
!           end do ! I_m1m2
!        
!           b_projectors(I_band_loop,k_vec_loop,Iproj1) = b_projector 
!
!       end do ! Iproj1
!
!     end do ! I_band_loop
!  end do ! k_vec_loop
!  call get_timing(time2)
!  write(142,10) '#         time :',time2-time1,' seconds '
!  write(142, 5) '#-----------------------------------------------------------'
!
!
!  !--------------------------------------------------------------------------
!  ! Next, compute the contribution  to each combination {Iproj1,Iproj2}
!  ! of the projected Wannier basis.
!  !--------------------------------------------------------------------------
!
!  write(142,5) '#                                                           '
!  write(142,5) '#  read, update, write chi0 array to file...                '
!  write(142,7) '#      - nk_vec     = ',nk_vec
!  write(142,7) '#      - I_band_max = ',I_band_max
!
!
!  record_index_chi0 = 0
!  call get_timing(time1)
!  do Iproj2 = 1,basis_size_intw
!        do Iproj1 = 1,basis_size_intw
!  
!
!           call get_timing(c_time1)
!
!           partial_chi0(:) = cmplx_0
!
!	   d_i_hw_avg = ZERO
!
!           record_index_M = 0
!           do k_vec_loop = 1,nk_vec
!              do I_band_loop = 1,I_band_max
!
!
!  	         nbnd1 = list_nbnd1(I_band_loop ,k_vec_loop)
!                 nbnd2 = list_nbnd2(I_band_loop ,k_vec_loop)
!
!  	         i_hw_min = list_i_hw_min(I_band_loop ,k_vec_loop)
!  	         i_hw_max = list_i_hw_max(I_band_loop ,k_vec_loop)
!
!                 ! if there is no contribution, move on to next iteration
!		 if (nbnd1 == 0 ) exit
!		 if (i_hw_min == 0) cycle
!
!		 d_i_hw_avg = d_i_hw_avg + &
!			dble(i_hw_max-i_hw_min)/dble(nk_vec*I_band_max)
!              
!                 ! get projectors to go to the projected basis
!
!		 !b_I1_projector =  b_projectors(Iproj1,I_band_loop,k_vec_loop)
!		 !b_I2_projector =  b_projectors(Iproj2,I_band_loop,k_vec_loop)
!
!		 b_I1_projector =  b_projectors(I_band_loop,k_vec_loop,Iproj1)
!		 b_I2_projector =  b_projectors(I_band_loop,k_vec_loop,Iproj2)
!
!          	 do i_hw = i_hw_min, i_hw_max 
!             	    record_index_M = record_index_M + 1 
!	  	    read(matrix_data_io_unit ,rec=record_index_M) matrix_data(i_hw)
!                 end do
!
!                 partial_chi0(i_hw_min:i_hw_max) =         &
!                   partial_chi0(i_hw_min:i_hw_max) +       &
!		   	matrix_data(i_hw_min:i_hw_max)     &
!		   	*conjg(b_I1_projector)*b_I2_projector 
!
!                ! partial_chi0(i_hw_min:i_hw_max) =                        &
!                !   partial_chi0(i_hw_min:i_hw_max) +                      &
!		!   matrix_data(i_hw_min:i_hw_max,I_band_loop ,k_vec_loop) &
!		!   *conjg(b_I1_projector)*b_I2_projector 
!
!
!              end do !I_band_loop
!           end do !k_vec_loop
!           call get_timing(c_time2)
!	   c_time = c_time2-c_time1
!
!
!           ! read-update-write to bin file!
!
!           call get_timing(read_time1)
!           record_index_chi0 = record_index_chi0 + 1 
!	   read(chi0_io_unit ,rec=record_index_chi0) chi0_tmp(:)
!
!           chi0_tmp(:) = chi0_tmp(:) + partial_chi0(:)
!
!	   write(chi0_io_unit ,rec=record_index_chi0) chi0_tmp(:)
!           call get_timing(read_time2)
!           read_time = read_time2-read_time1
!
!           write(142,15) '# I1,I2 = ',Iproj1,Iproj2,                  &
!                         ' <i_hw_max-i_hw_min> = ',nint(d_i_hw_avg),  &
!			 ': c time = ',c_time,' sec, read_time = ',read_time,' sec'
!        end do ! Iproj1 
!      end do !Iproj2
!      call get_timing(time2)
!
!      write(142,10) '#         time :',time2-time1,' seconds '
!      write(142, 5) '#==========================================================='
!
!      close(matrix_data_io_unit)
!
!
!  5   format(A)
!  7   format(A,I4)
!  10  format(A,F8.4,A)
!  15  format(A,2I4,A,I4,A,E10.4,A,E10.4,A)
!
!  end subroutine update_chi0_intw_bin


  subroutine normalize_chi0_intw_bin(chi0_io_unit)
  !------------------------------------------------------------------------
  ! This subroutine normalizes the response function in the netcdf file 
  ! and closes the file.
  !------------------------------------------------------------------------
  implicit none


  ! input
  integer        :: chi0_io_unit


  ! computation variables
  real(dp)       :: volume,  rho_0

  integer        :: Iproj1, Iproj2
  integer        :: record_index, i_hw
  integer        :: test_io_unit

  complex(dp)    :: chi0_tmp(n_hw)
  character(20)  :: filename

  !------------------------------------------------------------------

  volume = alat**3*abs(                                &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))

  rho_0     = 1.0_dp/volume
  
  !--------------------------------------------------------------------
  ! read data, normalize and write back
  !--------------------------------------------------------------------
  record_index = 0
  do Iproj2 = 1,basis_size_intw
        do Iproj1 = 1,basis_size_intw

	      record_index = record_index + 1
	      read(chi0_io_unit ,rec=record_index) chi0_tmp(:)

	      ! cheap test: is chi_KS really complex?
	      !write(filename,100) Iproj1,Iproj2
	      !test_io_unit = find_free_unit()
              !open(unit = test_io_unit, file = filename)
	      !do i_hw = 1, n_hw
	      !	 write(test_io_unit,10) chi0_tmp(i_hw)
              !end do
              !close(test_io_unit)


              ! normalize with respect to sum on k-points!
	      chi0_tmp(:) = chi0_tmp(:)*rho_0/dble(nk1s*nk2s*nk3s)

	      write(chi0_io_unit ,rec=record_index) chi0_tmp(:)

    end do 
  end do 


  10   format(2F12.6)
  100  format('test_chi0_',I1,'_',I1,'.dat')

  end subroutine normalize_chi0_intw_bin


  subroutine apply_kramers_kronig_to_chi0_bin(chi0_io_unit)
  !------------------------------------------------------------------------
  ! This subroutine computes the "real part" of chi0 using the
  ! "imaginary part", using the Kramers-Kronig (KK) FFT algorithm. 
  ! Note that since the KK equations are linear, they can be applied to
  ! complex functions, namely the "real part" need not be real, nor does
  ! the "imaginary part" need be imaginary.
  !------------------------------------------------------------------------
  implicit none

  integer :: chi0_io_unit
  real(dp)       :: wmin, wmax

  complex(dp),allocatable  :: list_Fi_w(:)
  complex(dp),allocatable  :: list_Fr_w(:)

  logical  :: double_data
  integer  :: iw_loop
  integer  :: Mfft_input
  integer  :: Iproj1, Iproj2
  integer  :: record_index, i_hw

  complex(dp)  :: partial_chi0(n_hw)
  complex(dp)  :: chi0_tmp(n_hw)



  ! set up the Kramers-Kronig scheme

  wmax = hw_max

  if (hw_min == -hw_max) then
     ! use the data as is
     wmin =  hw_min
     double_data = .False.
     Mfft_input  = n_hw-1
  else if( abs(hw_min) < eps_10) then
     ! symmetrize below w = 0
     wmin = -hw_max
     double_data = .True.
     Mfft_input  = 2*(n_hw-1)
  else
     ! something is wrong. This should be tested for 
     ! at the begining of the execution of the program...
     write(*,*) '**************************************'
     write(*,*) '*  Problems with frequencies:        *'
     write(*,*) '*  cannot perform Kramers-Kronig     *'
     write(*,*) '*     program stops.                 *'
     write(*,*) '**************************************'
     stop
  end if

 
  ! set up the KK scheme 
  call setup_kramers_kronig_arrays(wmin,wmax,Mfft_input)
  call setup_FFTW_kramers_kronig_plans()
  call get_W_trapezoid()
  call get_alpha_trapezoid()
  call get_beta_trapezoid()


  ! allocate the working arrays for the KK scheme
  allocate(list_Fi_w(Mfft+1))
  allocate(list_Fr_w(Mfft+1))

  ! loop on all components of the response function, which
  ! is assumed to only contain the "imaginary part" 


  record_index = 0
  do Iproj2 = 1,basis_size_intw
        do Iproj1 = 1,basis_size_intw

           record_index = record_index  + 1

           ! read from file; we assume only the "imaginary part"
           ! is present.
	   read(chi0_io_unit ,rec=record_index) partial_chi0(:)

           ! dump the "imaginary part" of chi0 in the working array

           if ( double_data == .False. ) then
		  list_Fi_w(:) =  partial_chi0(:) 
	   else 
                  !**************************************
                  ! I'm not sure what is below is valid.
                  ! It comes from an earlier version of
                  ! the algorithm. For now, it is best
                  ! not to use double_data...
                  !**************************************
		  list_Fi_w(Mfft/2+1:Mfft+1) =  partial_chi0(:) 
				
                  if ( magnon == .True.) then
		     ! if this is a magnon calculation, zero pad
                     ! on the negative part of the axis.
                  	do iw_loop =1, Mfft/2
		     		list_Fi_w(iw_loop)  =  cmplx_0
                  	end do
                  else
                  	do iw_loop =1, Mfft/2
		     		list_Fi_w(iw_loop)  =  -list_Fi_w(Mfft+2-iw_loop)
                  	end do
                  end if


           end if

           ! apply the KK-FFT scheme
           call apply_kramers_kronig_FFT(list_Fi_w,list_Fr_w)


           ! add the "real part" of chi0 to the "imaginary part"
           if ( double_data == .False. ) then

               partial_chi0(:) = partial_chi0(:) + list_Fr_w(:)

           else
               partial_chi0(:) = partial_chi0(:) + list_Fr_w(Mfft/2+1:Mfft+1)
           end if


           ! write back to file
	   write(chi0_io_unit ,rec=record_index) partial_chi0(:)


    end do !Iproj2
  end do !Iproj1



  ! clean up
  deallocate(list_Fi_w)
  deallocate(list_Fr_w)

  call deallocate_kramers_kronig_arrays()

  end subroutine apply_kramers_kronig_to_chi0_bin


  subroutine read_projectors()
  !------------------------------------------------------------------------
  ! This subroutine reads in the V_bar matrix from a netcdf file.
  !------------------------------------------------------------------------
  use netcdf
  implicit none


  ! netcdf variables
  integer :: nc_id        
  ! netcdf file id

  real(dp),allocatable :: re_projectors(:,:)

  integer   :: re_projectors_id, im_projectors_id, L_id, basis_size_id
  integer   :: list_m1_id, list_m2_id
  integer   :: dim_L

  character(256) :: name

  !------------------------------------------------------------------
  ! Read data from the projectors file, which gives relevant information
  ! about the projected basis.
  !------------------------------------------------------------------
  call check_netcdf(nf90_open(projectors_filename, NF90_NOWRITE, nc_id))

  call check_netcdf(nf90_inq_dimid(nc_id, "L", L_id))
  call check_netcdf(nf90_inq_dimid(nc_id, "I1_projected", basis_size_id))

  call check_netcdf(nf90_inquire_dimension(nc_id, L_id, name, dim_L))
  call check_netcdf(nf90_inquire_dimension(nc_id, basis_size_id, name, basis_size_intw))

  call check_netcdf(nf90_inq_varid(nc_id, "real_projector", re_projectors_id))
  call check_netcdf(nf90_inq_varid(nc_id, "imaginary_projector", im_projectors_id))
  call check_netcdf(nf90_inq_varid(nc_id, "list_m1", list_m1_id))
  call check_netcdf(nf90_inq_varid(nc_id, "list_m2", list_m2_id))


  allocate(list_m1_intw(dim_L))
  allocate(list_m2_intw(dim_L))
  call check_netcdf(nf90_get_var(nc_id, list_m1_id, list_m1_intw))
  call check_netcdf(nf90_get_var(nc_id, list_m2_id, list_m2_intw))


  allocate(projectors_intw(dim_L,basis_size_intw))
  allocate(re_projectors(dim_L,basis_size_intw))

  call check_netcdf(nf90_get_var(nc_id, re_projectors_id, re_projectors))

  projectors_intw(:,:) = cmplx_1*re_projectors(:,:)

  call check_netcdf(nf90_get_var(nc_id, im_projectors_id, re_projectors))
  projectors_intw(:,:) = projectors_intw(:,:) +  &
 			 	cmplx_i*re_projectors(:,:)

  deallocate(re_projectors)

  call check_netcdf(nf90_close(nc_id))

  end subroutine read_projectors


  subroutine output_chi0_intw_to_netcdf(chi0_io_unit)
  !------------------------------------------------------------------------
  ! This subroutine creates and initializes the netcdf file which will
  ! contain chi0 in nectdf format, in the Wannier basis,
  ! along with relevant information to interpret the data.
  !------------------------------------------------------------------------
  use netcdf
  implicit none

  ! input variable
  integer :: chi0_io_unit


  ! computation variables
  character(256) :: filename

  ! netcdf variables
  integer :: nc_id        
  ! netcdf file id

  integer :: start(3), count(3)
  integer :: Iproj1, Iproj2
   
  integer :: var_hw_id, var_re_chi0_id,var_im_chi0_id, var_q_id
  ! variables id: frequency, response function, q vector

  integer :: dim_hw_id, dim_chi0_id
  integer :: dim_space_id
  integer :: dim_Iproj1_id, dim_Iproj2_id
  ! dimension id: frequency, response function, 
  !               3D space dimension, Wannier basis indices dimensions.

  integer :: dim_hw, dim_chi0,  dim_space, dim_Iproj
  ! dimensions: frequency, response function, 
  !	        dimension of 3D space (3!), wannier basis dimensions

  character (*), parameter :: hw_name       = "frequency"
  character (*), parameter :: q_name        = "q_vector"
  character (*), parameter :: re_chi0_name  = "real_chi_KS"
  character (*), parameter :: im_chi0_name  = "imaginary_chi_KS"
  character (*), parameter :: space_name    = "space_coordinates"
  character (*), parameter :: Iproj1_name   = "I_projected_1"
  character (*), parameter :: Iproj2_name   = "I_projected_2"

  ! Define the units
  character (*), parameter :: UNITS       = "units"

  character (*), parameter :: hw_UNITS    = "eV"
  character (*), parameter :: chi0_UNITS  = "1/(eV*a_0^3)"
  character (*), parameter :: space_UNITS = "crystal_basis"

  character(256) :: name


  character(256) :: time_stamp, hostname
  integer        :: time_values(8)

  complex(dp)    :: chi0_tmp(n_hw)

  integer        :: record_index

  !------------------------------------------------------------------
  ! open and prepare a netcdf file for the response function.
  !------------------------------------------------------------------

  ! prepare the netcdf file, writing meta-data

  ! create file
  filename = trim(chi0_file_prefix)//'.nc'

  call check_netcdf(nf90_create(filename, NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Kohn-Sham susceptibility function in Wannier basis") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension, 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the "frequency" dimension, as well as the frequency variable
  dim_hw = n_hw
  call check_netcdf( nf90_def_dim(nc_id,hw_name, dim_hw, dim_hw_id) )

  dim_Iproj = basis_size_intw
  ! create the "wannier bands" dimensions
  call check_netcdf( nf90_def_dim(nc_id,Iproj1_name, dim_Iproj, dim_Iproj1_id) )
  call check_netcdf( nf90_def_dim(nc_id,Iproj2_name, dim_Iproj, dim_Iproj2_id) )

  !------------------------------
  ! define the variables
  !------------------------------
  ! create the "q-vector" array
  call check_netcdf( nf90_def_var(nc_id,q_name, NF90_DOUBLE, dim_space_id, var_q_id) )

  ! create the "frequency" array
  call check_netcdf( nf90_def_var(nc_id,hw_name, NF90_DOUBLE, dim_hw_id, var_hw_id) )

  ! create the real and imaginary part of the response function 

  call check_netcdf( nf90_def_var(nc_id,re_chi0_name, NF90_DOUBLE,            &
     (/dim_hw_id,dim_Iproj1_id,dim_Iproj2_id /), var_re_chi0_id) )

  call check_netcdf( nf90_def_var(nc_id,im_chi0_name, NF90_DOUBLE,            &
     (/dim_hw_id,dim_Iproj1_id,dim_Iproj2_id /), var_im_chi0_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! frequency
  call check_netcdf( nf90_put_att(nc_id, var_hw_id, UNITS, hw_UNITS) )

  ! q vector
  call check_netcdf( nf90_put_att(nc_id, var_q_id, UNITS, space_UNITS) )

  ! response function
  call check_netcdf( nf90_put_att(nc_id, var_im_chi0_id, UNITS, chi0_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_chi0_id, UNITS, chi0_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------

  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_q_id,(/qpt1,qpt2,qpt3/)) )

  ! write the frequencies
  call check_netcdf( nf90_put_var(nc_id, var_hw_id, list_hbar_omega ))

  ! write the response function to file

  count(:) = (/n_hw,1,1/)

  record_index = 0
  do Iproj2 = 1,basis_size_intw
     do Iproj1 = 1,basis_size_intw
	   record_index = record_index + 1

           ! read the data from the bin file
	   read(chi0_io_unit ,rec=record_index) chi0_tmp(:)

           start(:) = (/1,Iproj1,Iproj2/)

	   call check_netcdf( nf90_put_var(nc_id,var_re_chi0_id, &
					real(chi0_tmp(:)), start, count) ) 

	   call check_netcdf( nf90_put_var(nc_id,var_im_chi0_id, &
					aimag(chi0_tmp(:)), start, count) ) 
    end do 
  end do 

  ! close netcdf file
  call check_netcdf( nf90_close(nc_id) )

  ! close scratch bin file!
  close(chi0_io_unit )

  end subroutine output_chi0_intw_to_netcdf

  subroutine update_chi0_intw_sparse_matrix_algorithm(chi0_io_unit,  &
		nk_vec_max,nk_vec, eig_ks,eig_ksq, u_ks,u_ksq,fk_m_fkq, widths)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to chi0.
  !
  ! The details of exactly how the contribution to chi0 is calculated
  ! (eg, only imaginary part vs. real and imaginary; gaussian vs. lorentzian
  !  broadening,etc...) is left to another subroutine.
  !
  ! The data is written to a binary file which is assumed to be already
  ! opened.
  !
  ! This subroutine attempts to be as efficient as possible by
  ! using the fact that the result that is seeked is the product
  ! of sparse matrices!
  !------------------------------------------------------------------------
  use omp_lib
  implicit none

  ! input variables
  integer      :: chi0_io_unit
  integer      :: nk_vec, nk_vec_max

  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)

  complex(dp)  :: u_ks (num_wann,num_wann,nk_vec_max)
  complex(dp)  :: u_ksq(num_wann,num_wann,nk_vec_max)

  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)
  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)


  ! computation variables

  integer      :: record_index_chi0

  integer      :: Iproj1, Iproj2
  integer      :: m1, m2, m3, m4
  real(dp)     :: delta,  de


  integer      :: i_hw, i_hw_min, i_hw_max
  integer      :: k_vec_loop, iw_loop
  integer      :: nbnd1, nbnd2

  integer      :: k_nn, k_nn_max

  real(dp)     :: ek, ekq, f_m_f, hw

  real(dp)     :: d_i_hw_avg

  complex(dp)  :: chi0_tmp(n_hw)
  complex(dp)  :: d_chi0(n_hw), partial_chi0(n_hw)


  ! read-write the array "matrix_data" from a bin file
  integer ::  matrix_data_io_unit, record_length_M, record_index_M

!  complex(dp),allocatable  :: matrix_data(:,:,:)
  integer :: matrix_data_nnz, loop_nnz
  integer,allocatable     :: matrix_data_rows(:),matrix_data_columns(:)
  complex(dp),allocatable :: matrix_data_coo_format_values(:)
  complex(dp),allocatable :: vector(:)

  complex(dp)  :: matrix_data(n_hw)

  complex(dp),allocatable  :: b_projectors(:,:)
!  complex(dp),allocatable  :: b_projector_product(:,:)

  complex(dp)  :: b_projector, b_I1_projector, b_I2_projector

  complex(dp)  :: work_u1(num_wann), work_u2(num_wann), work_V(num_wann*num_wann)

  integer      :: I_band, I_band_loop, I_band_max, I_band_local_max
  integer      :: I_m1m2

  integer,allocatable  :: list_nbnd1(:), list_nbnd2(:) 
  integer,allocatable  :: list_k_vec_loop(:)
  integer,allocatable  :: list_i_hw_min(:), list_i_hw_max(:) 

  logical      :: cycle_out

  character    ::  matdescra(6)


  real(dp)     :: time1, time2
  real(dp)     :: read_time1, read_time2, read_time
  real(dp)     :: c_time1, c_time2, c_time

  !--------------------------------------------------------------------------
  ! 1)
  ! The summed indices, {k n1 n2},  will be combined into a super-index knn.
  ! Only the combinations which can give a contribution are kept.
  !
  ! 2) Define a super-array
  !
  !	M(i_w , knn ) = 
  !		[f_{n1 k}-f_{n2 k+q}]* delta(hw+ E_{n1 k} -E_{n2 k+q})
  !
  !     Clearly, this array is sparse because of the narrowness of the delta!
  !
  ! 3) Define a projector array
  !
  !    b(knn, Iproj) = 
  !		sum_{m1m2} projector[{m1m2},Iproj] * U_k[n1,m1]*U^dagger_{k+q}[m2,n2]
  !
  ! 4) The contributions to chi0 that we seek to compute are of the form
  !	chi0[ i_hw, Iproj1, Iproj2] =
  !			sum_{knn} M(i_w , knn )  b(knn, Iproj1) b^*(knn, Iproj2) 
  !
  !	which, for a given Iproj1, Iproj2 combination, is a spare matrix
  !	multiplied by a vector.
  !--------------------------------------------------------------------------



  write(142,5) '#==========================================================='
  write(142,5) '#                          new block                        '
  write(142,5) '#-----------------------------------------------------------'
  
  !-------------------------------------------------------------------------
  ! First, probe the fk_m_fkq array to see how large the (n1n2) dimension 
  ! must be
  !-------------------------------------------------------------------------

  write(142,5) '#                                                           '
  write(142,5) '#  Tabulating bands and frequency indices ...               '
  call get_timing(time1)


  ! find out exactly how many non-zero terms there are
  k_nn_max = 0
  do k_vec_loop = 1,nk_vec
     do I_band = 1,num_wann*num_wann

          nbnd1  = modulo(I_band-1,num_wann)+1
          nbnd2  = (I_band-nbnd1)/num_wann+1

          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle
          k_nn_max = k_nn_max + 1

     end do
     
  end do


  ! next, allocate an array that will hold the data, and tabulate
  ! the band indices
  allocate(list_nbnd1(k_nn_max))
  allocate(list_nbnd2(k_nn_max))

  allocate(list_k_vec_loop(k_nn_max))


  allocate(list_i_hw_min(k_nn_max))
  allocate(list_i_hw_max(k_nn_max))

  list_k_vec_loop  = 0

  list_nbnd1    = 0
  list_nbnd2    = 0
  list_i_hw_min = 0
  list_i_hw_max = 0

  call get_timing(time2)
  write(142,10) '#         time :',time2-time1,' seconds '
  write(142, 5) '#-----------------------------------------------------------'

  !-------------------------------------------------------------------------
  ! Although it is deeply displeasing, we must write the frequency stuff
  ! to file or compute it twice because we don't know the size of the 
  ! arrays!; once we know how many non-zero elements the data array will
  ! have, we can read the data back in.
  !-------------------------------------------------------------------------

  matrix_data_nnz = 0

  write(142,5) '#                                                           '
  write(142,5) '#  computing number of non-zero entries in sparse matrix... '

  call get_timing(time1)

  d_i_hw_avg  = 0

  k_nn = 0
  do k_vec_loop = 1,nk_vec
     do I_band = 1,num_wann*num_wann

          nbnd1  = modulo(I_band-1,num_wann)+1
          nbnd2  = (I_band-nbnd1)/num_wann+1

          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle
          k_nn = k_nn + 1
  
          list_k_vec_loop(k_nn)  = k_vec_loop

  	  list_nbnd1(k_nn) =  nbnd1
          list_nbnd2(k_nn) =  nbnd2


          ek    =  eig_ks (nbnd1,k_vec_loop)
          ekq   =  eig_ksq(nbnd2,k_vec_loop)
          de    =  ek-ekq

          f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 

          delta = widths(nbnd1,nbnd2,k_vec_loop)

  	  call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)

          ! if cycle_out, there is no fequency contribution. Do not
          ! update indices, they are already set to zero 
	  if (cycle_out) cycle

  	  list_i_hw_min(k_nn) = i_hw_min
  	  list_i_hw_max(k_nn) = i_hw_max

	  d_i_hw_avg = d_i_hw_avg + &
			dble(i_hw_max-i_hw_min+1)/dble(k_nn_max)


          matrix_data_nnz = matrix_data_nnz+i_hw_max-i_hw_min+1

     end do
  end do
  call get_timing(time2)
  write(142, 5) '#         dimensions:'
  write(142, 7) '#      	-    nk_vec     = ',nk_vec
  write(142, 7) '#      	-    k_nn_max   = ',k_nn_max
  write(142, 5) '#'
  write(142, 7) '#            -   matrix_data_nnz   = ',matrix_data_nnz 
  write(142, 7) '#     	      - <i_hw_max-i_hw_min> = ',nint(d_i_hw_avg)

  write(142,10) '#         time :',time2-time1,' seconds '
  write(142, 5) '#-----------------------------------------------------------'

  !-------------------------------------------------------------------------
  ! We can now fill in the rows and columns of the matrix!
  !-------------------------------------------------------------------------

  allocate(matrix_data_rows(matrix_data_nnz))
  allocate(matrix_data_columns(matrix_data_nnz))
  allocate(matrix_data_coo_format_values(matrix_data_nnz))


  ! initialize arrays
  matrix_data_coo_format_values = cmplx_0

  write(142,5) '#                                                           '
  write(142,5) '#  Building frequency dependent data array ...              '


  loop_nnz = 0

  call get_timing(time1)
  do k_nn = 1, k_nn_max

          k_vec_loop  = list_k_vec_loop(k_nn)
        
  	  nbnd1 = list_nbnd1(k_nn)
          nbnd2 = list_nbnd2(k_nn)

          ek    =  eig_ks (nbnd1,k_vec_loop)
          ekq   =  eig_ksq(nbnd2,k_vec_loop)
          de    =  ek-ekq

          delta = widths(nbnd1,nbnd2,k_vec_loop)

  	  call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)

          ! if cycle_out, there is no frequency contribution
	  if (cycle_out) cycle

          f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 


          do i_hw = i_hw_min, i_hw_max
          	loop_nnz = loop_nnz + 1 

                matrix_data_coo_format_values(loop_nnz ) = f_m_f*d_chi0(i_hw)

                ! columns are knn, rows are fequencies
  		matrix_data_rows(loop_nnz)      = i_hw
  		matrix_data_columns(loop_nnz)   = k_nn

          end do

  end do ! k_nn
  call get_timing(time2)
  write(142,10) '#         time :',time2-time1,' seconds '
  write(142, 5) '#-----------------------------------------------------------'

  !--------------------------------------------------------------------------
  ! Next, compute the projector functions to project onto the 
  ! projected Wannier basis.
  !--------------------------------------------------------------------------
  write(142,5) '#                                                           '
  write(142,5) '#  Building projector arrays ...                            '

  call get_timing(time1)

  allocate(b_projectors(basis_size_intw,k_nn_max))

  b_projectors(:,:) = cmplx_0


  !$omp parallel default(none)                                   &
  !$omp shared(k_nn_max,list_nbnd1,list_nbnd2, u_ks, u_ksq)      &
  !$omp shared(cmplx_0,num_wann,list_i_hw_min, list_k_vec_loop)  &
  !$omp shared(list_m1_intw, list_m2_intw, b_projectors)         &
  !$omp shared( basis_size_intw, projectors_intw)                &
  !$omp private( k_nn ,i_hw_min, nbnd1, nbnd2, k_vec_loop )      &
  !$omp private( work_u1, work_u2, Iproj1, b_projector )         &
  !$omp private( work_V, I_m1m2 ,m1, m2 )           

  !$omp do 
  do k_nn = 1, k_nn_max

        i_hw_min = list_i_hw_min(k_nn)  

        ! don't bother to calculate if there is no contribution!
	if (i_hw_min == 0 ) cycle

  	nbnd1 = list_nbnd1(k_nn)
        nbnd2 = list_nbnd2(k_nn)

        k_vec_loop = list_k_vec_loop(k_nn)

	work_u1(:) =       u_ks (nbnd1,:,k_vec_loop)
	work_u2(:) = conjg(u_ksq(nbnd2,:,k_vec_loop))


        do Iproj1 = 1, basis_size_intw

           b_projector = cmplx_0
	   
           work_V (:) = projectors_intw(:,Iproj1)


           do I_m1m2 = 1, num_wann**2   

              m1 = list_m1_intw(I_m1m2)
              m2 = list_m2_intw(I_m1m2)

	      b_projector = b_projector   +        & 
			work_V(I_m1m2)*work_u1(m1)*work_u2(m2)

           end do ! I_m1m2
        
  	   b_projectors(Iproj1,k_nn) = b_projector 


       end do ! Iproj1

  end do ! k_nn
  !$omp end do
  !$omp end parallel

  call get_timing(time2)
  write(142,10) '#         time :',time2-time1,' seconds '
  write(142, 5) '#-----------------------------------------------------------'


  !--------------------------------------------------------------------------
  ! Next, compute the contribution  to each combination {Iproj1,Iproj2}
  ! of the projected Wannier basis.
  !--------------------------------------------------------------------------

  write(142,5) '#                                                           '
  write(142,5) '#  read, update, write chi0 array to file...                '

  record_index_chi0 = 0


  matdescra(1) = 'G' ! general matrix
  matdescra(4) = 'F' ! fortran indexing


  allocate(vector(k_nn_max))

  call get_timing(time1)
  do Iproj2 = 1,basis_size_intw
        do Iproj1 = 1,basis_size_intw
  

           call get_timing(c_time1)

           partial_chi0(:) = cmplx_0

           !$omp parallel default(none)               &
           !$omp shared(b_projectors, Iproj1, Iproj2) &
           !$omp shared(vector,k_nn_max)              &
           !$omp private( k_nn )
           !$omp do
           do k_nn = 1, k_nn_max
              vector(k_nn) = b_projectors(Iproj1,k_nn)*conjg(b_projectors(Iproj2,k_nn))
           end do
           !$omp end do
           !$omp end parallel

	  !write(177,'(A,I2,I2)') 'Iproj1, Iproj2 = ',Iproj1,Iproj2
          !do k_nn = 1, k_nn_max
	  !	write(177,'(I4,2F12.8)') k_nn, vector(k_nn)
          ! end do


	   ! sparse matrix * vector subroutine. See MKL manual
           ! for details.

        ! I *think* this one DOES NOT run in parallel with openmp
!	   call mkl_zcoomv('N',   	& ! "normal", as in "not transposed"
!			   n_hw,        & ! number of rows
!			   k_nn_max,    & ! number of columns
!		           cmplx_1,     & ! prefactor of operation: should be 1
!			   matdescra,   & ! descriptor for the matrix
!         matrix_data_coo_format_values, & ! matrix components, in COO format
!                    matrix_data_rows,   & ! non-zero row indices, in COO format
!  		 matrix_data_columns,   & ! non-zero column indices, in COO format
!	             matrix_data_nnz,   & ! dimension of 3 arrays above
!                              vector,   & ! "vector" to which matrix is applied
!			     cmplx_0,   & ! just a number, zero
!                          partial_chi0)   ! result of matrix * vector product!
!


        ! I *think* this one runs in parallel with openmp
        ! this is the simplified interface
	call mkl_zcoogemv(  'N',        & ! "normal", as in "not transposed"
                           n_hw,        & ! number of rows
         matrix_data_coo_format_values, & ! matrix components, in COO format
                    matrix_data_rows,   & ! non-zero row indices, in COO format
  		 matrix_data_columns,   & ! non-zero column indices, in COO format
	             matrix_data_nnz,   & ! dimension of 3 arrays above
                              vector,   & ! "vector" to which matrix is applied
                          partial_chi0)   ! result of matrix * vector product!



           call get_timing(c_time2)

	   c_time = c_time2-c_time1


           ! read-update-write to bin file!

           call get_timing(read_time1)
           record_index_chi0 = record_index_chi0 + 1 
	   read(chi0_io_unit ,rec=record_index_chi0) chi0_tmp(:)

           chi0_tmp(:) = chi0_tmp(:) + partial_chi0(:)

	   write(chi0_io_unit ,rec=record_index_chi0) chi0_tmp(:)
           call get_timing(read_time2)
           read_time = read_time2-read_time1

           write(142,15) '# I1,I2 = ',Iproj1,Iproj2,                  &
			 ': c time = ',c_time,' sec, read_time = ',read_time,' sec'
        end do ! Iproj1 
      end do !Iproj2
      call get_timing(time2)

      write(142,10) '#         time :',time2-time1,' seconds '
      write(142, 5) '#==========================================================='

      close(matrix_data_io_unit)


  5   format(A)
  7   format(A,I8)
  10  format(A,F8.4,A)
  15  format(A,2I4,A,E10.4,A,E10.4,A)

  end subroutine update_chi0_intw_sparse_matrix_algorithm

  subroutine compute_dchi0_cold_smearing_KK &
	(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
  !------------------------------------------------------------------------
  ! Compute the contribution to chi0 using Marzari-Vanderbilt width;
  ! only the imaginary part is calculated and the real part will be 
  ! obtained from Kramers-Kronig.
  !
  ! The cold smearing function, by Marzari and Vanderbilt, is given by
  !
  !	delta(x) = (2-sqrt(2)x)/sqrt(pi) exp(-(x-1/sqrt(2))^2)
  !
  !------------------------------------------------------------------------
  implicit none

  
  ! delta(x) < 10^-16*delta_max if x < x_minus or x > x_plus
  ! the parameters were obtained numerically.

  real(dp),parameter  :: x_minus = -5.524654195658513
  real(dp),parameter  :: x_plus  =  6.920282881677849

  ! input variables
  real(dp)     :: delta  ! width of gaussian
  real(dp)     :: de     ! e_k-e_k+q

  ! output variables
  complex(dp)  :: d_chi0(n_hw)
  logical      :: cycle_out
  integer      :: i_hw_min, i_hw_max


  ! computation variables
  real(dp)     :: e_y2, y, y2
  real(dp)     :: sqrt2, one_sqrt2
  complex(dp)  :: i_sqrt_pi_on_delta
  integer      :: i_hw
  real(dp)     :: const_hw1 , const_hw2, const_hw3 

  ! only compute the contribution to the delta-function
  ! where the function doesn't virtually vanish!
  const_hw1 = -de-hw_min
  const_hw2 = delta*x_minus 
  const_hw3 = delta*x_plus

  sqrt2     = sqrt(2.0_dp)
  one_sqrt2 = ONE/sqrt2

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


  i_sqrt_pi_on_delta = -cmplx_i*sqrt_pi/delta

  ! only compute the contributions to d_chi0 which
  ! are larger than machine-precision
  do i_hw = i_hw_min, i_hw_max

     y    = (list_hbar_omega(i_hw) + de)/delta
     y2   = (y-one_sqrt2)**2
     e_y2 = dexp(-y2)

     d_chi0(i_hw) = e_y2*(2.0_dp-sqrt2*y)*i_sqrt_pi_on_delta
  end do

  end subroutine compute_dchi0_cold_smearing_KK


  subroutine update_chi0_intw(nk_vec_max,nk_vec,           &
			eig_ks,eig_ksq,u_ks,u_ksq,         &
			fk_m_fkq, widths)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to chi0_intw, the
  ! array of coefficients in the MLWF basis.
  !
  ! The details of exactly how the contribution to chi0_intw is calculated
  ! (eg, only imaginary part vs. real and imaginary; gaussian vs. lorentzian
  !  broadening,etc...) is left to another subroutine.
  !
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max
  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)
  complex(dp)  :: u_ks   (num_wann,num_wann,nk_vec_max)
  complex(dp)  :: u_ksq  (num_wann,num_wann,nk_vec_max)

  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)
  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)

  ! computation variables
  real(dp)     :: delta,  de

  integer      :: k_vec_loop
  integer      :: nbnd1, nbnd2
  real(dp)     :: ek, ekq, f_m_f
  complex(dp)  :: numerator
  complex(dp)  :: d_chi0(n_hw)

  integer      :: i_hw_min, i_hw_max

  complex(dp)  :: b_projectors(basis_size_intw)

  integer      :: I_band
  logical      :: cycle_out

  integer      :: number_of_calls

  real(dp)     :: b_time, add_time, chi_time
  real(dp)     :: time1, time2

  b_time  =  0.0_dp
  add_time=  0.0_dp
  chi_time=  0.0_dp

  number_of_calls = 0
  do k_vec_loop = 1,nk_vec
     do I_band = 1,num_wann*num_wann

          nbnd1  = modulo(I_band-1,num_wann)+1
          nbnd2  = (I_band-nbnd1)/num_wann+1

          if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle


          ek    =  eig_ks (nbnd1,k_vec_loop)
          ekq   =  eig_ksq(nbnd2,k_vec_loop)
          de    =  ek-ekq

          f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 

          delta = widths(nbnd1,nbnd2,k_vec_loop)


          ! compute the contribution to chi0_intw. Let the subroutine
          ! decide what method to use!

         ! call get_timing(time1) 
  	  call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)
 	 ! call get_timing(time2) 
	 ! chi_time = chi_time+time2-time1 

	  if (cycle_out) cycle

         ! call get_timing(time1) 
          call build_b_projectors(k_vec_loop,nk_vec_max, nbnd1,nbnd2, &
						u_ks, u_ksq, b_projectors)
 	 ! call get_timing(time2) 
	 ! b_time = b_time+time2-time1 


         ! call get_timing(time1) 
  	  call add_d_chi0_intw(d_chi0,f_m_f, i_hw_min, i_hw_max, b_projectors)
 	 ! call get_timing(time2) 
	 ! add_time = add_time+time2-time1 
	 ! number_of_calls = number_of_calls  + 1 


     end do  ! I_band
  end do  ! kvec_loop

  !write(142, 5) '#==========================================================='
  !write(142, 7) '# number of calls:',number_of_calls
  !write(142,10) '# compute d_chi0_intw  time :',chi_time,' seconds '
  !write(142,10) '# compute b_projectors time :',b_time,' seconds '
  !write(142,10) '# add d_chi0_intw      time :',add_time,' seconds '
  !write(142, 5) '#==========================================================='

  5   format(A)
  7   format(A,I8)
  10  format(A,F8.4,A)
  15  format(A,2I4,A,E10.4,A,E10.4,A)


  end subroutine update_chi0_intw

  subroutine build_b_projectors(k_vec_loop,nk_vec_max, nbnd1,nbnd2, &
			u_ks, u_ksq, b_projectors)
  !------------------------------------------------------------------------
  ! This subroutine builds the projectors necessary to extract the 
  ! chi0_intw coefficients.
  !
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: k_vec_loop, nk_vec_max
  integer      :: nbnd1, nbnd2
  complex(dp)  :: u_ks (num_wann,num_wann,nk_vec_max)
  complex(dp)  :: u_ksq(num_wann,num_wann,nk_vec_max)

  ! output variables
  complex(dp)  :: b_projectors(basis_size_intw)

  ! local variables
  complex(dp)  :: work_u(num_wann**2)
  integer      :: Iproj, I_m1m2
  integer      :: m1, m2



  do I_m1m2 = 1, num_wann**2   

     m1 = list_m1_intw(I_m1m2)
     m2 = list_m2_intw(I_m1m2)

     work_u(I_m1m2) = u_ks (nbnd1,m1,k_vec_loop)*conjg(u_ksq(nbnd2,m2,k_vec_loop))

  end do ! I_m1m2


  ! use MKL routine to compute matrix-vector product
  ! computes y = alpha*A^T*x+beta*y
  call zgemv('T',                &      ! matrix should be transposed
              num_wann**2,       &      ! number of rows in matrix
              basis_size_intw,   &      ! number of columns in matrix
 	      cmplx_1,           &      ! prefactor in equation
	      projectors_intw,   &      ! matrix A
              num_wann**2,       &      ! first dimension of A
	      work_u,            &      ! vector x
	      1,                 &      ! increment of x
	      cmplx_0,           &      ! beta
	      b_projectors,      &      ! y array
	      1)                        ! increment of y



  end subroutine build_b_projectors


  subroutine add_d_chi0_intw ( d_chi0, f_m_f , i_hw_min, i_hw_max, b_projectors)
  !------------------------------------------------------------------------
  ! This subroutine adds the contribution of d_chi0 to chi0.
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  complex(dp)  :: d_chi0(n_hw)
  integer      :: i_hw_min, i_hw_max
  real(dp)     :: f_m_f

  complex(dp)  :: b_projectors(basis_size_intw)


  ! local variables

  complex(dp)  :: b_projectors_coefficient

  integer      :: Iproj1, Iproj2, II

  integer      :: i_hw
  complex(dp)  :: numerator


  real(dp)     :: time1, time2

!  call get_timing(time1) 

  ! combine the Iproj1 and Iproj2 loop

  ! OMP IS SLOWER THAN SERIAL!
  !!$omp parallel default(none)                               &
  !!$omp shared(basis_size_intw,b_projectors, f_m_f,d_chi0)   &
  !!$omp shared( chi0_intw)                                   &
  !!$omp private(II,Iproj1,Iproj2, b_projectors_coefficient ) 

  !!$omp do 
  do II = 1, basis_size_intw**2

     Iproj1 = modulo(II-1,basis_size_intw)+1
     Iproj2 = (II-Iproj1)/basis_size_intw+1


     b_projectors_coefficient    =          &
         b_projectors(Iproj1)*conjg(b_projectors(Iproj2))*f_m_f


     chi0_intw(:,Iproj1,Iproj2)  =            &
         chi0_intw(:,Iproj1,Iproj2)+b_projectors_coefficient*d_chi0(:)


  end do

  !!$omp end do
  !!$omp end parallel


  end subroutine add_d_chi0_intw

  subroutine normalize_chi0_intw()
  !------------------------------------------------------------------------
  ! This subroutine normalizes the response function.
  !------------------------------------------------------------------------
  implicit none

  ! computation variables
  real(dp)       :: volume,  rho_0, spin_sum_factor


  volume = alat**3*abs(                                &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))

  if (nspin == 1) then
	spin_sum_factor = 2.0_dp
  else 
	spin_sum_factor = 1.0_dp
  end if


  rho_0     = 1.0_dp/volume

  
  !--------------------------------------------------------------------
  ! normalize with respect to sum on k-points!
  chi0_intw(:,:,:) = chi0_intw(:,:,:)*rho_0/dble(nk1s*nk2s*nk3s)*spin_sum_factor 



  end subroutine normalize_chi0_intw

  subroutine apply_kramers_kronig_to_chi0_intw()
  !------------------------------------------------------------------------
  ! This subroutine computes the "real part" of chi0 using the
  ! "imaginary part", using the Kramers-Kronig (KK) FFT algorithm. 
  ! Note that since the KK equations are linear, they can be applied to
  ! complex functions, namely the "real part" need not be real, nor does
  ! the "imaginary part" need be imaginary.
  !------------------------------------------------------------------------
  implicit none

  integer :: chi0_io_unit
  real(dp)       :: wmin, wmax

  complex(dp),allocatable  :: list_Fi_w(:)
  complex(dp),allocatable  :: list_Fr_w(:)

  logical  :: double_data
  integer  :: iw_loop
  integer  :: Mfft_input
  integer  :: Iproj1, Iproj2
  integer  :: i_hw

  complex(dp)  :: partial_chi0(n_hw)
  complex(dp)  :: chi0_tmp(n_hw)



  ! set up the Kramers-Kronig scheme

  wmax = hw_max

  if (hw_min == -hw_max) then
     ! use the data as is
     wmin =  hw_min
     double_data = .False.
     Mfft_input  = n_hw-1
  else if( abs(hw_min) < eps_10) then
     ! symmetrize below w = 0
     wmin = -hw_max
     double_data = .True.
     Mfft_input  = 2*(n_hw-1)
  else
     ! something is wrong. This should be tested for 
     ! at the begining of the execution of the program...
     write(*,*) '**************************************'
     write(*,*) '*  Problems with frequencies:        *'
     write(*,*) '*  cannot perform Kramers-Kronig     *'
     write(*,*) '*     program stops.                 *'
     write(*,*) '**************************************'
     stop
  end if

 
  ! set up the KK scheme 
  call setup_kramers_kronig_arrays(wmin,wmax,Mfft_input)
  call setup_FFTW_kramers_kronig_plans()
  call get_W_trapezoid()
  call get_alpha_trapezoid()
  call get_beta_trapezoid()


  ! allocate the working arrays for the KK scheme
  allocate(list_Fi_w(Mfft+1))
  allocate(list_Fr_w(Mfft+1))

  ! loop on all components of the response function, which
  ! is assumed to only contain the "imaginary part" 

  do Iproj2 = 1,basis_size_intw
        do Iproj1 = 1,basis_size_intw

	   partial_chi0(:) = chi0_intw(:,Iproj1,Iproj2)

           ! dump the "imaginary part" of chi0 in the working array

           if ( double_data == .False. ) then
		  list_Fi_w(:) =  partial_chi0(:) 
	   else 
                  !**************************************
                  ! I'm not sure what is below is valid.
                  ! It comes from an earlier version of
                  ! the algorithm. For now, it is best
                  ! not to use double_data...
                  !**************************************
		  list_Fi_w(Mfft/2+1:Mfft+1) =  partial_chi0(:) 
				
                  if ( magnon == .True.) then
		     ! if this is a magnon calculation, zero pad
                     ! on the negative part of the axis.
                  	do iw_loop =1, Mfft/2
		     		list_Fi_w(iw_loop)  =  cmplx_0
                  	end do
                  else
                  	do iw_loop =1, Mfft/2
		     		list_Fi_w(iw_loop)  =  -list_Fi_w(Mfft+2-iw_loop)
                  	end do
                  end if


           end if

           ! apply the KK-FFT scheme
           call apply_kramers_kronig_FFT(list_Fi_w,list_Fr_w)


           ! add the "real part" of chi0 to the "imaginary part"
           if ( double_data == .False. ) then

               partial_chi0(:) = partial_chi0(:) + list_Fr_w(:)

           else
               partial_chi0(:) = partial_chi0(:) + list_Fr_w(Mfft/2+1:Mfft+1)
           end if


           ! write back to array
	   chi0_intw(:,Iproj1,Iproj2) = partial_chi0(:)

    end do !Iproj2
  end do !Iproj1



  ! clean up
  deallocate(list_Fi_w)
  deallocate(list_Fr_w)

  call deallocate_kramers_kronig_arrays()

  end subroutine apply_kramers_kronig_to_chi0_intw

  subroutine output_chi0_intw_to_netcdf_2()
  !------------------------------------------------------------------------
  ! This subroutine creates and initializes the netcdf file which will
  ! contain chi0 in nectdf format, in the Wannier basis,
  ! along with relevant information to interpret the data.
  !------------------------------------------------------------------------
  use netcdf
  implicit none



  ! computation variables
  character(256) :: filename

  ! netcdf variables
  integer :: nc_id        
  ! netcdf file id

  integer :: Iproj1, Iproj2
   
  integer :: var_hw_id, var_re_chi0_id,var_im_chi0_id, var_q_id
  ! variables id: frequency, response function, q vector

  integer :: dim_hw_id, dim_chi0_id
  integer :: dim_space_id
  integer :: dim_Iproj1_id, dim_Iproj2_id
  ! dimension id: frequency, response function, 
  !               3D space dimension, Wannier basis indices dimensions.

  integer :: dim_hw, dim_chi0,  dim_space, dim_Iproj
  ! dimensions: frequency, response function, 
  !	        dimension of 3D space (3!), wannier basis dimensions

  character (*), parameter :: hw_name       = "frequency"
  character (*), parameter :: q_name        = "q_vector"
  character (*), parameter :: re_chi0_name  = "real_chi_KS"
  character (*), parameter :: im_chi0_name  = "imaginary_chi_KS"
  character (*), parameter :: space_name    = "space_coordinates"
  character (*), parameter :: Iproj1_name   = "I_projected_1"
  character (*), parameter :: Iproj2_name   = "I_projected_2"

  ! Define the units
  character (*), parameter :: UNITS       = "units"

  character (*), parameter :: hw_UNITS    = "eV"
  character (*), parameter :: chi0_UNITS  = "1/(eV*a_0^3)"
  character (*), parameter :: space_UNITS = "crystal_basis"

  character(256) :: name


  character(256) :: time_stamp, hostname
  integer        :: time_values(8)



  !------------------------------------------------------------------
  ! open and prepare a netcdf file for the response function.
  !------------------------------------------------------------------

  ! prepare the netcdf file, writing meta-data

  ! create file
  filename = trim(chi0_file_prefix)//'.nc'

  call check_netcdf(nf90_create(filename, NF90_CLOBBER, nc_id))

  !------------------------------------------------------------------
  ! write meta-data to allow user to know what is in the file
  !------------------------------------------------------------------

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "title",  &
        		"Kohn-Sham susceptibility function in Wannier basis") )

  call check_netcdf( nf90_put_att(nc_id, NF90_GLOBAL, "source",  &
        "program INTW by Bruno Rousseau and Asier Eiguren (EHU-UPV / CSIC / DIPC)") )

  call date_and_time(VALUES=time_values)

  write(time_stamp,'(A,I2,A,I2,A,I4)') &
  	'Produced on dd/mm/yyyy :',time_values(3),'/',time_values(2),'/',time_values(1)

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "production_date",  time_stamp))

  call check_netcdf(nf90_put_att(nc_id, NF90_GLOBAL, "note",  &
	"This data is not meant to be fully self describing: you should already"//&
	" know what system this file describes!"))

  !------------------------------
  ! define the dimensions 
  !------------------------------

  ! create the "3D-space" dimension, 
  dim_space = 3
  call check_netcdf( nf90_def_dim(nc_id,space_name, dim_space, dim_space_id) )

  ! create the "frequency" dimension, as well as the frequency variable
  dim_hw = n_hw
  call check_netcdf( nf90_def_dim(nc_id,hw_name, dim_hw, dim_hw_id) )

  dim_Iproj = basis_size_intw
  ! create the "wannier bands" dimensions
  call check_netcdf( nf90_def_dim(nc_id,Iproj1_name, dim_Iproj, dim_Iproj1_id) )
  call check_netcdf( nf90_def_dim(nc_id,Iproj2_name, dim_Iproj, dim_Iproj2_id) )

  !------------------------------
  ! define the variables
  !------------------------------
  ! create the "q-vector" array
  call check_netcdf( nf90_def_var(nc_id,q_name, NF90_DOUBLE, dim_space_id, var_q_id) )

  ! create the "frequency" array
  call check_netcdf( nf90_def_var(nc_id,hw_name, NF90_DOUBLE, dim_hw_id, var_hw_id) )

  ! create the real and imaginary part of the response function 

  call check_netcdf( nf90_def_var(nc_id,re_chi0_name, NF90_DOUBLE,            &
     (/dim_hw_id,dim_Iproj1_id,dim_Iproj2_id /), var_re_chi0_id) )

  call check_netcdf( nf90_def_var(nc_id,im_chi0_name, NF90_DOUBLE,            &
     (/dim_hw_id,dim_Iproj1_id,dim_Iproj2_id /), var_im_chi0_id) )

  !------------------------------------------------
  ! Assign unit attributes to the variables
  !-------------------------------------------------

  ! frequency
  call check_netcdf( nf90_put_att(nc_id, var_hw_id, UNITS, hw_UNITS) )

  ! q vector
  call check_netcdf( nf90_put_att(nc_id, var_q_id, UNITS, space_UNITS) )

  ! response function
  call check_netcdf( nf90_put_att(nc_id, var_im_chi0_id, UNITS, chi0_UNITS) )
  call check_netcdf( nf90_put_att(nc_id, var_re_chi0_id, UNITS, chi0_UNITS) )

  !--------------------------------------------------------------------
  ! End define mode. This tells netCDF we are done defining meta-data.
  !--------------------------------------------------------------------
  call check_netcdf( nf90_enddef(nc_id) )

  !--------------------------------------------------------------------
  ! Write the data to the file!
  !--------------------------------------------------------------------

  ! write the q-vector
  call check_netcdf( nf90_put_var(nc_id, var_q_id,(/qpt1,qpt2,qpt3/)) )

  ! write the frequencies
  call check_netcdf( nf90_put_var(nc_id, var_hw_id, list_hbar_omega ))

  ! write the response function to file


  call check_netcdf( nf90_put_var(nc_id,var_re_chi0_id, &
					real(chi0_intw(:,:,:))) ) 

  call check_netcdf( nf90_put_var(nc_id,var_im_chi0_id, &
					aimag(chi0_intw(:,:,:)))) 

  ! close netcdf file
  call check_netcdf( nf90_close(nc_id) )


  end subroutine output_chi0_intw_to_netcdf_2

  subroutine update_chi0_intw_omp(nk_vec_max,nk_vec,       &
			eig_ks,eig_ksq,u_ks,u_ksq,         &
			fk_m_fkq, widths)
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
  use omp_lib
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max
  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)
  complex(dp)  :: u_ks   (num_wann,num_wann,nk_vec_max)
  complex(dp)  :: u_ksq  (num_wann,num_wann,nk_vec_max)

  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)
  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)

  ! computation variables
  real(dp)     :: delta,  de

  integer      :: k_vec_loop

  integer      :: knn
  integer      :: nbnd1, nbnd2
  real(dp)     :: ek, ekq, f_m_f
  complex(dp)  :: numerator
  complex(dp)  :: d_chi0(n_hw)

  integer      :: i_hw, i_hw_min, i_hw_max

  complex(dp)  :: b_projectors(basis_size_intw)

  complex(dp)  :: f_m_f_bb
  integer      :: num_wann_2,basis_size_intw_2

  integer      :: I_band
  integer      :: II, Iproj1, Iproj2
  logical      :: cycle_out


  ! timing variable
  real(dp)     :: time1, time2, time3, time4, time5, time6
  real(dp)     :: time_indices, time_fetching, time_dchi0, time_proj
  real(dp)     :: time_uploc, time_upglo

  integer      :: thread_id


  complex(dp)  :: chi0_local(n_hw,basis_size_intw,basis_size_intw)


  ! output timing to different log file!
  integer      :: io_unit
  character(*),parameter :: log_file = 'update_chi0_timing.log'

  io_unit = find_free_unit()

  open(unit=io_unit,file=log_file,status='unknown',access='append')

  write(io_unit,5) '#==================================================='
  write(io_unit,5) '# New timing block:'
  write(io_unit,5) '#==================================================='
  write(io_unit,5) '#         what            node       timing (sec)'
  write(io_unit,5) '#---------------------------------------------------'



  ! some convenience variables, to speed up the calculations
  num_wann_2        = num_wann**2
  basis_size_intw_2 = basis_size_intw**2

  ! combine the loops in k, n1 and n2

  !$omp parallel default(none)                                           &
  !$omp shared(cmplx_0,nk_vec,num_wann,fk_m_fkq, eig_ks, eig_ksq )       &
  !$omp shared(widths,nk_vec_max,u_ks, u_ksq, chi0_intw,n_hw,ZERO)       &
  !$omp shared(basis_size_intw,basis_size_intw_2,num_wann_2,io_unit)     &
  !$omp private(chi0_local,knn,I_band,nbnd1, nbnd2,ek, ekq)              &
  !$omp private(de, f_m_f, delta,d_chi0,cycle_out, i_hw_min, i_hw_max)   &
  !$omp private(time1,time2,time3,time4,time5,time6)                     &
  !$omp private(time_indices, time_fetching, time_dchi0, time_proj   )   &
  !$omp private(time_uploc, time_upglo,thread_id )                       &
  !$omp private(k_vec_loop,b_projectors,II, Iproj1, Iproj2, i_hw,f_m_f_bb )

  chi0_local(:,:,:) = cmplx_0

  time_indices  = ZERO
  time_fetching = ZERO
  time_dchi0    = ZERO
  time_proj     = ZERO
  time_uploc    = ZERO
  time_upglo    = ZERO

  ! First, build the chi0_local arrays, on each thread.
  !$omp do 
  
  do knn = 1,nk_vec*num_wann_2

     ! extract k_vec_loop, nbnd1, nbnd2

     call get_timing(time1)

     I_band      = modulo(knn-1,num_wann_2)+1
     k_vec_loop  = (knn-I_band)/num_wann_2+1

     nbnd1       = modulo(I_band-1,num_wann)+1
     nbnd2       = (I_band-nbnd1)/num_wann+1

     call get_timing(time2)

     time_indices = time_indices +time2-time1 

     if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle


     ek    =  eig_ks (nbnd1,k_vec_loop)
     ekq   =  eig_ksq(nbnd2,k_vec_loop)
     de    =  ek-ekq

     f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 

     delta = widths(nbnd1,nbnd2,k_vec_loop)

     call get_timing(time3)

     time_fetching = time_fetching+time3-time2

     ! compute the contribution to chi0_intw. Let the subroutine
     ! decide what method to use!

     call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)


     call get_timing(time4)

     time_dchi0 = time_dchi0+time4-time3


     if (cycle_out) cycle

     call build_b_projectors(k_vec_loop,nk_vec_max, nbnd1,nbnd2, &
						u_ks, u_ksq, b_projectors)

     call get_timing(time5)

     time_proj  = time_proj+time5-time4

      do II = 1, basis_size_intw_2

     	Iproj1 = modulo(II-1,basis_size_intw)+1
     	Iproj2 = (II-Iproj1)/basis_size_intw+1

 	! convenience variable
	f_m_f_bb = f_m_f*b_projectors(Iproj1)*conjg(b_projectors(Iproj2))

        do i_hw = i_hw_min, i_hw_max
     	   chi0_local(i_hw,Iproj1,Iproj2)  =                        &
         	chi0_local(i_hw,Iproj1,Iproj2)+d_chi0(i_hw)*f_m_f_bb 

!		performing these products inside the inner loop is wasteful!
!         	chi0_local(i_hw,Iproj1,Iproj2)+d_chi0(i_hw)*f_m_f*  &
!         	b_projectors(Iproj1)*conjg(b_projectors(Iproj2))
        end do !i_hw

      end do !II

     call get_timing(time6)

     time_uploc = time_uploc+time6-time5


  end do  ! knn
  !$omp end do


  !------------------------------------------------------------
  ! I think a barrier is not necessary.
  ! In fact, avoiding the barrier may be useful,
  ! that way different processors will reach the reduction
  ! area at different moments and may not hit a lock!
  !-----------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!$omp barrier 



  ! Next, dump chi0_local into chi0_intw;
  ! each thread should dump all its components here, so don't
  ! parallelize over the loop variables! We must still be in the
  ! "parallel" environment, however, or else chi0_local becomes
  ! undefined

  !---------------------------------------------------------
  ! Is it better to use critical or atomic???
  !---------------------------------------------------------

  call get_timing(time1)
  !$omp critical

  do Iproj2 = 1, basis_size_intw
     do Iproj1 = 1, basis_size_intw
        do i_hw = 1, n_hw
           !!!!!!!$omp atomic
     	   chi0_intw(i_hw,Iproj1,Iproj2)  =                 &
     	   		chi0_intw(i_hw,Iproj1,Iproj2)  +    &
     				chi0_local(i_hw,Iproj1,Iproj2) 
  	end do ! i_hw
     end do ! Iproj1
  end do ! Iproj2

  !$omp end critical
  call get_timing(time2)
  time_upglo = time2-time1 

  ! output timing! 
  thread_id = omp_get_thread_num()

  write(io_unit,10) trim('getting indices  '),thread_id, time_indices
  write(io_unit,10) trim('fetching data    '),thread_id, time_fetching 
  write(io_unit,10) trim('computing dchi0  '),thread_id, time_dchi0 
  write(io_unit,10) trim('computing b-proj '),thread_id, time_proj  
  write(io_unit,10) trim('local update     '),thread_id, time_uploc 
  write(io_unit,10) trim('global update    '),thread_id, time_upglo 
  !$omp end parallel


  close(io_unit)

  5   format(A)
  10  format(A20,5X,I4,5X,G14.4E3)

  end subroutine update_chi0_intw_omp

  subroutine apply_kramers_kronig_to_chi0_intw_2()
  !------------------------------------------------------------------------
  ! This subroutine computes the "real part" of chi0 using the
  ! "imaginary part", using the Kramers-Kronig (KK) relation 
  ! and the FFT algorithm to perform a convolution.
  ! Note that since the KK equations are linear, they can be applied to
  ! complex functions, namely the "real part" need not be real, nor does
  ! the "imaginary part" need be imaginary.
  !------------------------------------------------------------------------
  implicit none



  complex(dp),allocatable  :: list_Fi_w(:)
  complex(dp),allocatable  :: list_Fr_w(:)

  integer  :: Iproj1, Iproj2
  integer  :: i_hw

  complex(dp)  :: partial_chi0(n_hw)
  complex(dp)  :: chi0_tmp(n_hw)


 ! Define the size of the convolution, which is expressed
 ! through a global variable
  convolution_size = 2*n_hw

  ! allocate the working arrays for the KK scheme
  allocate(list_Fi_w(convolution_size))
  allocate(list_Fr_w(convolution_size))


  ! Set up the various arrays and plans for FFTW
  call setup_FFTW_kramers_kronig_plans_convolution()


  ! build the Kernel, in Fourier space
  call build_convolution_kernel()


  ! initialize the array
  list_Fi_w(:) = cmplx_0

  ! loop on all components of the response function, which
  ! is assumed to only contain the "imaginary part" 

  do Iproj2 = 1,basis_size_intw
        do Iproj1 = 1,basis_size_intw

	   partial_chi0(:) = chi0_intw(:,Iproj1,Iproj2)

           ! dump the "imaginary part" of chi0 in the working array

	   list_Fi_w(1:n_hw) =  partial_chi0(:) 

           ! apply the KK-FFT scheme
           call apply_kramers_kronig_FFT_convolution      &
					(list_Fi_w,list_Fr_w)


           ! add the "real part" of chi0 to the "imaginary part"
           partial_chi0(:) = partial_chi0(:) + list_Fr_w(1:n_hw)


           ! write back to array
	   chi0_intw(:,Iproj1,Iproj2) = partial_chi0(:)

    end do !Iproj2
  end do !Iproj1



  ! clean up
  deallocate(list_Fi_w)
  deallocate(list_Fr_w)
  call deallocate_kramers_kronig_convolution()

  end subroutine apply_kramers_kronig_to_chi0_intw_2

  subroutine update_chi0_intw_no_matrix_elements_omp(nk_vec_max,nk_vec,       &
			eig_ks,eig_ksq, fk_m_fkq, widths)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to chi0_intw, the
  ! array of coefficients in the MLWF basis.
  !
  ! This subroutine neglects the matrix elements entirely, assuming 
  ! the occupation factors already account for spin indices.
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
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)

  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)
  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)

  ! computation variables
  real(dp)     :: delta,  de

  integer      :: k_vec_loop

  integer      :: knn
  integer      :: nbnd1, nbnd2
  real(dp)     :: ek, ekq, f_m_f
  complex(dp)  :: numerator
  complex(dp)  :: d_chi0(n_hw)

  integer      :: i_hw, i_hw_min, i_hw_max


  integer      :: I_band
  logical      :: cycle_out


  real(dp)     :: time1, time2

  complex(dp)  :: chi0_local(n_hw,basis_size_intw,basis_size_intw)

  ! combine the loops in k, n1 and n2

  !$omp parallel default(none)                                            &
  !$omp shared(cmplx_0,nk_vec,num_wann,fk_m_fkq, eig_ks, eig_ksq )        &
  !$omp shared(widths,nk_vec_max,chi0_intw,n_hw)                          &
  !$omp private(chi0_local,knn,I_band,nbnd1,nbnd2,ek,ekq,k_vec_loop,i_hw) &
  !$omp private(de, f_m_f, delta,d_chi0,cycle_out, i_hw_min, i_hw_max) 


  chi0_local(:,:,:) = cmplx_0
  ! First, build the chi0_local arrays, on each thread.
  !$omp do 
  
  do knn = 1,nk_vec*num_wann**2

     ! extract k_vec_loop, nbnd1, nbnd2

     I_band      = modulo(knn-1,num_wann**2)+1
     k_vec_loop  = (knn-I_band)/num_wann**2+1

     nbnd1       = modulo(I_band-1,num_wann)+1
     nbnd2       = (I_band-nbnd1)/num_wann+1

     if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle

     if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) < 0) then
	write(*,*) '* Negative occupation factor! How is this possible??? *'
     end if

     ek    =  eig_ks (nbnd1,k_vec_loop)
     ekq   =  eig_ksq(nbnd2,k_vec_loop)
     de    =  ek-ekq

     f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 

     delta = widths(nbnd1,nbnd2,k_vec_loop)


     ! compute the contribution to chi0_intw. Let the subroutine
     ! decide what method to use!

     call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)

     if (cycle_out) cycle

     do i_hw = i_hw_min, i_hw_max
     	   chi0_local(i_hw,1,1)  = chi0_local(i_hw,1,1)+d_chi0(i_hw)*f_m_f
     end do !i_hw


  end do  ! knn
  !$omp end do

  !$omp barrier
  ! Next, dump chi0_local into chi0_intw;
  ! each thread should dump all its components here, so don't
  ! parallelize over the loop variables! We must still be in the
  ! "parallel" environment, however, or else chi0_local becomes
  ! undefined
  do i_hw = 1, n_hw
     !$omp atomic
     chi0_intw(i_hw,1,1)  = chi0_intw(i_hw,1,1)  + chi0_local(i_hw,1,1) 
  end do ! i_hw

  !$omp end parallel


  end subroutine update_chi0_intw_no_matrix_elements_omp

  subroutine update_chi0_intw_omp_2(nk_vec_max,nk_vec,       &
			  eig_ks,eig_ksq,u_ks,u_ksq,         &
			  fk_m_fkq, widths)
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
  ! We will try to squeeze more efficiency out of this subroutine, as
  ! this is really where all the time is spent!
  !------------------------------------------------------------------------
  use omp_lib
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max
  real(dp)     :: eig_ks (num_wann,nk_vec_max)
  real(dp)     :: eig_ksq(num_wann,nk_vec_max)
  complex(dp)  :: u_ks   (num_wann,num_wann,nk_vec_max)
  complex(dp)  :: u_ksq  (num_wann,num_wann,nk_vec_max)

  integer      :: fk_m_fkq(num_wann,num_wann,nk_vec_max)
  real(dp)     :: widths(num_wann,num_wann,nk_vec_max)

  ! computation variables
  real(dp)     :: delta,  de

  integer      :: k_vec_loop

  integer      :: knn
  integer      :: nbnd1, nbnd2
  real(dp)     :: ek, ekq, f_m_f
  complex(dp)  :: numerator
  complex(dp)  :: d_chi0(n_hw), d_chi0_matrix(n_hw,1)

  integer      :: i_hw_min, i_hw_max

  complex(dp)  :: b_projectors(basis_size_intw)

  complex(dp)  :: f_m_f_bb(basis_size_intw**2,1)

  integer      :: num_wann_2,basis_size_intw_2

  integer      :: I_band
  integer      :: II, Iproj1, Iproj2
  logical      :: cycle_out


  ! timing variable
  real(dp)     :: time1, time2, time3, time4, time5, time6
  real(dp)     :: time_indices, time_fetching, time_dchi0, time_proj
  real(dp)     :: time_uploc, time_upglo

  integer      :: thread_id


  ! combine the basis dimensions, for faster matrix operations!
  complex(dp)  :: chi0_local(n_hw,basis_size_intw**2) 

  integer      :: lda, ldb, ldc

  ! output timing to different log file!
  integer      :: io_unit
  character(*),parameter :: log_file = 'update_chi0_timing.log'

  io_unit = find_free_unit()

  open(unit=io_unit,file=log_file,status='unknown',access='append')

  write(io_unit,5) '#==================================================='
  write(io_unit,5) '# New timing block:'
  write(io_unit,5) '#==================================================='
  write(io_unit,5) '#         what            node       timing (sec)'
  write(io_unit,5) '#---------------------------------------------------'



  ! some convenience variables, to speed up the calculations
  num_wann_2        = num_wann**2
  basis_size_intw_2 = basis_size_intw**2

  chi0_local(:,:) = cmplx_0

  ! combine the loops in k, n1 and n2

  !$omp parallel default(none)                                           &
  !$omp shared(cmplx_0,nk_vec,num_wann,fk_m_fkq, eig_ks, eig_ksq )       &
  !$omp shared(widths,nk_vec_max,u_ks, u_ksq, chi0_intw,n_hw,ZERO,cmplx_1)&
  !$omp shared(basis_size_intw,basis_size_intw_2,num_wann_2,io_unit)     &
  !$omp firstprivate(chi0_local)                                         &
  !$omp private(knn,I_band,nbnd1, nbnd2,ek, ekq)                         &
  !$omp private(de, f_m_f, delta,d_chi0,cycle_out, i_hw_min, i_hw_max)   &
  !$omp private(d_chi0_matrix,time1,time2,time3,time4,time5,time6)       &
  !$omp private(time_indices, time_fetching, time_dchi0, time_proj)      &
  !$omp private(time_uploc, time_upglo,thread_id,lda,ldb,ldc )           &
  !$omp private(k_vec_loop,b_projectors,II, Iproj1, Iproj2, f_m_f_bb)    


  thread_id = omp_get_thread_num()


  time_indices  = ZERO
  time_fetching = ZERO
  time_dchi0    = ZERO
  time_proj     = ZERO
  time_uploc    = ZERO
  time_upglo    = ZERO

  ! First, build the chi0_local arrays, on each thread.
  !$omp do 
  do knn = 1,nk_vec*num_wann_2

     ! extract k_vec_loop, nbnd1, nbnd2

     call get_timing(time1)

     I_band      = modulo(knn-1,num_wann_2)+1
     k_vec_loop  = (knn-I_band)/num_wann_2+1

     nbnd1       = modulo(I_band-1,num_wann)+1
     nbnd2       = (I_band-nbnd1)/num_wann+1

     call get_timing(time2)

     time_indices = time_indices +time2-time1 

     if (fk_m_fkq(nbnd1,nbnd2,k_vec_loop) == 0) cycle


     ek    =  eig_ks (nbnd1,k_vec_loop)
     ekq   =  eig_ksq(nbnd2,k_vec_loop)
     de    =  ek-ekq

     f_m_f = dble(fk_m_fkq(nbnd1,nbnd2,k_vec_loop)) 

     delta = widths(nbnd1,nbnd2,k_vec_loop)

     call get_timing(time3)

     time_fetching = time_fetching+time3-time2

     ! compute the contribution to chi0_intw. Let the subroutine
     ! decide what method to use!

     call compute_dchi0(d_chi0, cycle_out,i_hw_min, i_hw_max,de,delta)

     d_chi0_matrix(:,1) = d_chi0(:)

     call get_timing(time4)

     time_dchi0 = time_dchi0+time4-time3


     if (cycle_out) cycle

     call build_b_projectors(k_vec_loop,nk_vec_max, nbnd1,nbnd2, &
						u_ks, u_ksq, b_projectors)


     ! perform the update of chi_local as a matrix-matrix operation!

     ! build one of the "matrices"
     do II = 1, basis_size_intw_2
     	Iproj1 = modulo(II-1,basis_size_intw)+1
     	Iproj2 = (II-Iproj1)/basis_size_intw+1

 	! convenience variable
	f_m_f_bb(II,1) = f_m_f*b_projectors(Iproj1)*conjg(b_projectors(Iproj2))
	!f_m_f_bb(II) = f_m_f*b_projectors(Iproj1)*conjg(b_projectors(Iproj2))
     end do

     call get_timing(time5)
     time_proj  = time_proj+time5-time4



     ! Use a BLAS routine to update chi0_local!

     ! operation performed,
     !		C = alpha*op(A)*op(B) + beta*C,
     ! where 
     !		alpha = 1
     !		beta  = 1
     !		 C    = chi0_local              dimension m x n
     !		op(A) = d_chi0,                 dimension m x k
     !		op(B) = transpose(f_m_f_bb),    dimension k x n
     !
     !	The dimensions are 
     !			m = n_hw
     !			n = basis_size_intw_2
     !			k = 1

     !zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)


!----------------------------------------------------------------------------------
! Issues with zgemm in an openMP parallel region
!----------------------------------------------------------------------------------
!
! Using zgemm to update chi0_local is SIGNIFICANTLY FASTER (factor 4 speed up!) 
! than updating explicitely with loops. It is thus worth while to use this 
! BLAS routine.
!
! HOWEVER, I have encountered strange bugs in this situation. I want to 
! briefly document what I have done in the hope that it will be useful in the
! future. Note that these notes are written with the ifort 11.1 compiler in 
! mind, and its corresponding MKL and omp versions. Also, I am not an expert
! and tests were not conducted rigorously, so comments below should be taken
! with a grain of salt.
!
! OBSERVATION:
! Calling the zgemm BLAS routine from MKL within an omp parallel region 
! CRASHES the code (ie segfault or glibc error), the SECOND TIME the subroutine 
! is called, IF arrays passed to the zgemm subroutine are sliced, ie
! passed as "A(imin:imax,:)". 
!
! Doing this would be useful if i_hw_max-i_hw_min+1 < n_hw.
!
! If arrays are not sliced, the code runs just fine.
!
! It would be interesting to be able to slice the arrays, but I won't for now
! so that the code will run. My current case of interest involves
! i_hw_max-i_hw_min+1 = n_hw, so I will feel no performance degradation:
! if you are using this subroutine in a different context, you might want to
! explore this bug a bit further!
!				Bruno Rousseau, post-doc, DIPC, November 2011.
!---------------------------------------------------------------------------------

! This is an example of a call that crashes the code!
!     call zgemm(	'N',                       & ! first array not transposed
!			'T',                       & ! second array is transposed 
!	                lda,                       & ! m dimension
!	  	        ldb,                       & ! n dimension
!		          1,                       & ! k dimension
!		    cmplx_1,                       & ! alpha
!	d_chi0_matrix(i_hw_min:i_hw_max,:),        & ! the A "matrix"
!	                lda,                       & ! first dimension of A
!		   f_m_f_bb,                       & ! the B "matrix"
!			ldb,                       & ! first dimension of B
!		    cmplx_1,                       & ! beta
!	chi0_local(i_hw_min:i_hw_max,:),           & ! C matrix
!			 ldc)                        ! first dimension of C matrix



! the version below appears to work properly
     lda = n_hw
     ldb = basis_size_intw_2
     ldc = lda

     call zgemm(	'N',                       & ! first array not transposed
			'T',                       & ! second array is transposed 
	                lda,                       & ! m dimension
	  	        ldb,                       & ! n dimension
		          1,                       & ! k dimension
		    cmplx_1,                       & ! alpha
	      d_chi0_matrix,                       & ! the A "matrix"
	                lda,                       & ! first dimension of A
		   f_m_f_bb,                       & ! the B "matrix"
			ldb,                       & ! first dimension of B
		    cmplx_1,                       & ! beta
		 chi0_local,                       & ! C matrix
			 ldc)                        ! first dimension of C matrix


!---------------------------------------------------
! If zgemm cannot be used, an explicit algorithm 
! can be used instead (much slower!). Note that
! some testing with matmul is also much slower
! than zgemm.
!---------------------------------------------------
!     do II = 1, basis_size_intw_2
!         chi0_local(i_hw_min:i_hw_max,II) =        &
!		chi0_local(i_hw_min:i_hw_max,II) + &
!                    d_chi0(i_hw_min:i_hw_max)*f_m_f_bb(II)
!     end do 

     ! matmul is very slow!
     !chi0_local = chi0_local + matmul(d_chi0_matrix,transpose(f_m_f_bb))


     call get_timing(time6)

     time_uploc = time_uploc+time6-time5


  end do  ! knn
  !$omp end do



  ! Next, dump chi0_local into chi0_intw;
  ! each thread should dump all its components here, so don't
  ! parallelize over the loop variables! We must still be in the
  ! "parallel" environment, however, or else chi0_local becomes
  ! undefined


  call get_timing(time1)
  !$omp critical
  do II = 1, basis_size_intw_2
	Iproj1 = modulo(II-1,basis_size_intw)+1
     	Iproj2 = (II-Iproj1)/basis_size_intw+1

     	chi0_intw(:,Iproj1,Iproj2)  = chi0_intw(:,Iproj1,Iproj2) + chi0_local(:,II) 
  end do ! II

  !$omp end critical
  call get_timing(time2)
  time_upglo = time2-time1 

  ! output timing! 

  write(io_unit,10) trim('getting indices  '),thread_id, time_indices
  write(io_unit,10) trim('fetching data    '),thread_id, time_fetching 
  write(io_unit,10) trim('computing dchi0  '),thread_id, time_dchi0 
  write(io_unit,10) trim('computing b-proj '),thread_id, time_proj  
  write(io_unit,10) trim('local update     '),thread_id, time_uploc 
  write(io_unit,10) trim('global update    '),thread_id, time_upglo 


  !$omp end parallel


  close(io_unit)

  5   format(A)
  10  format(A20,5X,I4,5X,F14.4)

  end subroutine update_chi0_intw_omp_2


!--------------------------------------------------------------------------------
!
end module intw_chi_0
!
!--------------------------------------------------------------------------------
