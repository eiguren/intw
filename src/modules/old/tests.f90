!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
module intw_tests
!----------------------------------------------------------------------------!
!
!       This module contains subroutines which performs TESTS, insuring 
!       that the tests subroutines behave as expected.  
!       
!       Testing is absolutely essential to code development. Having a
!       strong test battery builds confidence in the validity of the code.
!----------------------------------------------------------------------------!

use intw_useful_constants
use intw_input_parameters
use intw_utility
use intw_reading
use intw_w90
use intw_symmetries
use intw_matrix_elements
use intw_fft
use intw_plane_wave_setup

contains


  subroutine find_k_in_centered_1BZ_and_G(kpoint,nk1,nk2,nk3, &
				nG_shell_max_first_neighbors,kpt_cBZ,G_cBZ)
  !----------------------------------------------------------------------------!
  !     Given a kpoint(3) in crystal coordinates, this subroutine 
  !     generates:
  !            - kpt_cBZ(3),   which is equivalent to kpoint, but
  !			       which minimizes the cartesian distance to Gamma.
  !			       In other words, this point belongs to the traditional,
  !			       Gamma-centered 1BZ.
  !            - G(3),       G(i) an integer, the G vector which links kpoint and
  !                          kpt_cBZ.
  !
  !	Define
  !		kpoint = k_1BZ + G1
  !		       = (k_cBZ + G2) + G1
  !                    = k_cBZ + G,          with G = G1+G2
  !----------------------------------------------------------------------------!
  implicit none
  
  ! input
  integer        :: nk1, nk2, nk3                ! the mesh parameters
  integer        :: nG_shell_max_first_neighbors ! G vectors in 1st coordination shell
  real(dp)       :: kpoint(3)                    ! the kpoint of interest 
  
  ! output
  real(dp)       :: kpt_cBZ(3)   ! the k point in the 1BZ centered at Gamma
  integer        :: G_cBZ(3)     ! the translation vector

  ! internal variables
  integer        :: i, j, k   
  integer        :: G1(3), G2(3)
  real(dp)       :: kpt_in_1BZ(3)
  real(dp)       :: kpt_cryst(3), kpt_cart(3) 
  real(dp)       :: dist2, min_dist2

  ! transfer the kpoint to the 1BZ with boundaries between zero and one.
  call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpt_in_1BZ,G1)

  ! find the G vector which minimizes the distance between k_cBZ and Gamma


  ! initialize with the first G vector, which is Gamma
  kpt_cryst(:) = kpt_in_1BZ(:)-gvec(:,1)
  kpt_cart(:)  = kpt_cryst(1)*bg(:,1)+kpt_cryst(2)*bg(:,2)+kpt_cryst(3)*bg(:,3)
  min_dist2    = sum(kpt_cart(:)*kpt_cart(:))

  G2(:)        = gvec(:,1) ! this is Gamma
  kpt_cBZ(:)   = kpt_cryst(:)


  ! loop on first shell, looking for a G vector that yields a smaller distance
  do i =2,nG_shell_max_first_neighbors

	kpt_cryst(:) = kpt_in_1BZ(:)-gvec(:,i)
        kpt_cart(:)  = kpt_cryst(1)*bg(:,1)+kpt_cryst(2)*bg(:,2)+kpt_cryst(3)*bg(:,3)
	dist2        = sum(kpt_cart(:)*kpt_cart(:))

	if (dist2 < min_dist2) then
  		G2(:)        =  gvec(:,i) ! this is Gamma
  		kpt_cBZ(:)   =  kpt_cryst(:)
		min_dist2    =  dist2
        end if
  end do

  G_cBZ(:) = G1(:) + G2(:)

  end subroutine find_k_in_centered_1BZ_and_G




  subroutine test_symmetry_axis_angle()
  !-------------------------------------------------------------------------------!
  ! This subroutine tests the angle and the rotational axis related to the sym mat! 
  !-------------------------------------------------------------------------------!

  real(dp) :: axis_old(3), axis_new(3)
  real(dp) :: angle_old,   angle_new

  real(dp) :: daxis, daxis_plus, daxis_minus 
  integer :: isym , io_unit

  integer :: sym(3,3)

  logical :: test_angle, test_axis

  character(*),parameter :: log_file = 'test_symmetry_axis_angle.test'
  io_unit = find_free_unit()

  open(unit=io_unit,file=log_file,status='unknown')

  write(io_unit,'(a)') '#---------------------------------------------------'
  write(io_unit,'(a)') '# This test checks that the angle and axis of       '
  write(io_unit,'(a)') '# rotation for every symmetry operation is correct. '  
  write(io_unit,'(a)') '# It also checks the consistency of the older code  '  
  write(io_unit,'(a)') '# and the newer code.                               '  
  write(io_unit,'(a)') '# Rotation axis is in cartesian "a" coordinates.    '  
  write(io_unit,'(a)') '#---------------------------------------------------'


  test_angle = .true.
  test_axis  = .true.

  do isym=1,nsym

     sym = s(:,:,isym)

     call rotaxis(sym,axis_old,angle_old)
     call rotaxis_crystal(sym,axis_new,angle_new)

     write(io_unit, 3)'Symmetry Operation = ',isym
     write(io_unit, 5) 'version           Rotation axis               angle (on Pi)'
     write(io_unit,10) 'old :',axis_old(:),angle_old/pi
     write(io_unit,10) 'new :',axis_new(:),angle_new/pi
     write(io_unit, 5) ''


     if ( abs(angle_new-angle_old) < eps_8) then
     	write(io_unit, 5) 'Delta angle = 0 : PASSED'
     else
     	write(io_unit, 5) 'Delta angle = 0 : FAILED'
	test_angle = .false.
     end if

     if      ( abs(angle_old) < eps_8) then
     	write(io_unit, 5) ' angle = 0 : axis doesnt matter'

     else if ( abs(angle_old-pi) < eps_8) then
     	write(io_unit, 5) ' angle = pi : axis defined up to a sign'
        daxis_plus = dot_product(axis_old(:)+axis_new(:),axis_old(:)+axis_new(:))
        daxis_plus = sqrt(daxis_plus )

        daxis_minus= dot_product(axis_old(:)-axis_new(:),axis_old(:)-axis_new(:))
        daxis_minus= sqrt(daxis_minus)
     
        if ( daxis_minus < eps_8 .or. daxis_plus < eps_8) then
     	   write(io_unit, 5) 'axis_old +/- axis_new = 0 : PASSED'
        else
     	   write(io_unit, 5) 'axis_old +/- axis_new /= 0 : FAILED'
           test_axis  = .false.
        end if
     else 
     	write(io_unit, 5) ' angle /= 0, pi'
        daxis = dot_product(axis_old(:)-axis_new(:),axis_old(:)-axis_new(:))

        if ( daxis < eps_8 ) then
     	   write(io_unit, 5) 'axis_old - axis_new = 0 : PASSED'
        else
     	   write(io_unit, 5) 'axis_old - axis_new /= 0 : FAILED'
           test_axis  = .false.
        end if

     end if
        	

     write(io_unit, 5) '-------------------------------------------------'
  enddo 


  write(io_unit, 5) '================================================='

  if (.not. test_angle .or. .not. test_axis)  then
     	write(io_unit,5)  '|                  TEST FAILED!                =='
  else
     	write(io_unit,5)  '|                  TEST PASSED                 =='
  end if 

  write(io_unit, 5) '================================================='

  close(io_unit)

  3  format(A,I2)
  5  format(A)
  10 format(A,3F12.6,F12.6)

  end subroutine test_symmetry_axis_angle  

  subroutine test_switch_indices(nk_1,nk_2,nk_3)
  !----------------------------------------------------------------------------!
  ! This subroutine tests the subroutine "switch_indices", in intw_utility
  !----------------------------------------------------------------------------!
  implicit none

  integer        :: i,       j,       k     
  integer        :: ip,      jp,     kp     

  integer        :: nk_1,  nk_2,  nk_3

  integer        :: switch ,  io_unit

  integer        :: ikpt, n

  character(*),parameter :: log_file = 'test_switch_indices.test'
  
  io_unit = find_free_unit()
  open(unit=io_unit,file=log_file,status='unknown')


  write(io_unit,'(a)') '# This is a test of the switch_indices routine' 
  write(io_unit,'(a)') '# parameters:'
  write(io_unit,'(a,I4)') '# nk1 = ',nk_1
  write(io_unit,'(a,I4)') '# nk2 = ',nk_2
  write(io_unit,'(a,I4)') '# nk3 = ',nk_3
  write(io_unit,'(a)') ''
  
  write(io_unit,'(a)') '# First test :: (i,j,k)  --->   ikpt'
  write(io_unit,'(a)') '# second test::    n     ---> (i,j,k)'
  write(io_unit,'(a)') '#----------------------------------------------------------------'
  write(io_unit,'(a)') '#  i   j   k  ikpt   ||  ikpt-n   ||  n   i   j   k'
  write(io_unit,'(a)') '#----------------------------------------------------------------'
  n = 0
  do i=1,nk_1
    do j=1,nk_2
      do k=1,nk_3
        n  =  n + 1
        
        switch = +1
        call switch_indices(nk_1,nk_2,nk_3,ikpt,i,j,k,switch)
        
        switch = -1
        call switch_indices(nk_1,nk_2,nk_3,n,ip,jp,kp,switch)

        write(io_unit,'(4i4,8x,i4,7x,4i4)') i,j,k,ikpt,ikpt-n,n,ip,jp,kp

        if ( ikpt/= n) write(io_unit,'(a)') '!!! TEST FAILURE!!!' 
        if ( i  /= ip) write(io_unit,'(a)') '!!! TEST FAILURE!!!' 
        if ( j  /= jp) write(io_unit,'(a)') '!!! TEST FAILURE!!!' 
        if ( k  /= kp) write(io_unit,'(a)') '!!! TEST FAILURE!!!' 

      end do
    end do
  end do

  close(io_unit)
  end subroutine test_switch_indices
  
  subroutine echo_W90()
  !----------------------------------------------------------------------------!
  ! This subroutine simply ehoes the parameters read from the W90 environment, 
  ! using subroutine allocate_and_read_W90, in module intw_w90.
  !----------------------------------------------------------------------------!
  implicit none

  character(*),parameter :: log_file = 'echo_W90.test'
  integer                ::  io_unit

  io_unit = find_free_unit()
  open(unit=io_unit,file=log_file,form='formatted',status='unknown')

      write(io_unit,'(a)') "========================================="
      write(io_unit,'(a)') "=echo of W90 environment, from $seed.chk="
      write(io_unit,'(a)') "========================================="
      write(io_unit,'(a,I4)') "             num_bands       = ",num_bands
      write(io_unit,'(a,I4)') "             num_wann        = ",num_wann
      write(io_unit,'(a,I4)') "             nkpoints_QE     = ",nkpoints_QE
      write(io_unit,'(a,L4)') "   allocated(u_matrix)?      = ",allocated(u_matrix)
      write(io_unit,'(a,L4)') "   allocated(u_matrix_opt)?  = ",allocated(u_matrix_opt)
      write(io_unit,'(a,3I5)')"       shape(u_matrix)       = ",shape(u_matrix)
  if ( allocated(u_matrix_opt) ) then
      write(io_unit,'(a,3I5)')"   shape(u_matrix_opt)       = ",shape(u_matrix_opt)
  end if
      write(io_unit,'(a,a)') ""
      write(io_unit,'(a,a)') "========================================="

  close(io_unit)
  end subroutine echo_W90

  
  subroutine manipulate_u_matrix_opt(ikpt)
  !----------------------------------------------------------------------------!
  !	This subroutine will manipulate the u_matrix_opt object, for testing
  !	purposes. This matrix is read in intw_w90.
  !----------------------------------------------------------------------------!

  implicit none

  integer        :: i, j, k, io_unit, ikpt
  
  complex(dp)    :: uud(num_bands,num_bands),udu(num_wann,num_wann)


  do i=1,num_wann
     do j=1,num_wann
        udu(i,j) =  (0.0_dp,0.0_dp)
        do k = 1,num_bands  
           udu(i,j) =  udu(i,j) +          & 
           conjg(u_matrix_opt(k,i,ikpt))*u_matrix_opt(k,j,ikpt)
        end do
     end do
  end do
     
  do i=1,num_bands
     do j=1,num_bands
        uud(i,j) =  (0.0_dp,0.0_dp)
        do k = 1,num_wann
           uud(i,j) =  uud(i,j) +          & 
           u_matrix_opt(i,k,ikpt)*conjg(u_matrix_opt(j,k,ikpt))
        end do
     end do
  end do

  io_unit = find_free_unit()
  open(unit=io_unit,file=trim('u_matrix_opt')//'.dat',form='formatted',status='unknown')
  write(io_unit,'(a,I5)') "#  ikpt = ",ikpt
  write(io_unit,'(a)') "#------------------------------------------"

  write(io_unit,*) "#Re[u^dagger u]"
  write(io_unit,*) "#--------------"
  do i=1,num_wann
     write(io_unit,'(200F6.3)')  (real(udu(i,j)), j=1,num_wann)
  end do

  write(io_unit,*) "#Im[u^dagger u]"
  write(io_unit,*) "#--------------"
  do i=1,num_wann
     write(io_unit,'(200F6.3)')  (aimag(udu(i,j)), j=1,num_wann)
  end do

  write(io_unit,*) "#Re[u u^dagger]"
  write(io_unit,*) "#--------------"
  do i=1,num_bands
     write(io_unit,'(200F6.3)')  (real(uud(i,j)), j=1,num_bands)
  end do

  write(io_unit,*) "#Im[u u^dagger]"
  write(io_unit,*) "#--------------"
  do i=1,num_bands
     write(io_unit,'(200F6.3)')  (aimag(uud(i,j)), j=1,num_bands)
  end do

  close(io_unit)


  end subroutine manipulate_u_matrix_opt




  subroutine test_k_symmetry(kmesh,nkpoints_QE,kpoints_QE)
  !----------------------------------------------------------------------------!
  ! This subroutine tests the code which finds the inequivalent k points in the
  ! IBZ.  
  !----------------------------------------------------------------------------!
  implicit none

  integer        :: io_unit
 
  integer        :: nkpoints_QE

  real(dp)       :: kpoints_QE(3,nkpoints_QE)       ! k vectors in the QE folder 

  real(dp)       :: kmesh(3,nk1*nk2*nk3) ! k vectors on the mesh, in

  real(dp)       :: kQE(3), kMP(3), dk(3),  k_rot(3)
  real(dp)       :: norm 
  real(dp)       :: TR_sign


  integer        ::  ikpt, i_folder
  integer        ::  G(3)

  integer        ::  i_sym, inv
  integer        ::  i

  logical        ::  test 

  character(*),parameter :: log_file = 'test_k_symmetry.test'

  test  = .true.

  io_unit = find_free_unit()
  open(unit=io_unit,file=log_file,status='unknown')

  write(io_unit,*) '#=======================================================#'
  write(io_unit,*) '#              Test of the k symmetry code              #'
  write(io_unit,*) '#=======================================================#'
  write(io_unit,*) '#   Definitions:'
  write(io_unit,*) '#        -  kpoints_QE == kQE : k-point read in a QE folder'
  write(io_unit,*) '#'
  write(io_unit,*) '#        -     kmesh   == kMP : k-point on a MP mesh, in '
  write(io_unit,*) '#                               canonical order.         '
  write(io_unit,*) '#        -    i_folder        : index of folder where kQE'
  write(io_unit,*) '#                               was read.                '
  write(io_unit,*) '#'
  write(io_unit,*) '#=======================================================#'

  if (full_mesh ) then
        write(io_unit,*) '# A full mesh is available in the QE folders'
  else 
        write(io_unit,*) '# Only an IBZ is available in the QE folders'
  end if

  do ikpt = 1, nk1*nk2* nk3

     kMP = kmesh(:,ikpt)
     write(io_unit,*) '================================================================'
     write(io_unit,100) 'ikpt        = ',ikpt

     write(io_unit,75) 'kmesh = (',kMP(1),',', kMP(2),',',kMP(3),')'

     if (full_mesh) then
        i_folder = QE_folder_nosym(ikpt)
        G        = nosym_G(:,ikpt)
        kQE      =  kpoints_QE(:,i_folder)

        write(io_unit,*) ' '
        write(io_unit,*) '------ no symmetry ------'
        write(io_unit,*) ' '
        write(io_unit,60) 'i_folder   = ',i_folder
        write(io_unit,75) 'kpoint_QE  = (',kQE(1),',',kQE(2),',',kQE(3),')'
        write(io_unit,350)'nosym_G    = (',G(1),',',G(2),',', G(3),')'
        write(io_unit,*) ' '

        dk   = kQE-kMP-dble(G)
      
        norm = sqrt(dot_product(dk,dk))

        write(io_unit,300) '|kQE - kMP-G | = ',norm
        if ( norm < eps_8 ) then
           write(io_unit,310) 'PASSED'
        else
           write(io_unit,310) 'FAILED!'
           test = .false.
        end if 

     end if

     write(io_unit,*) ' '
     write(io_unit,*) '------ symmetry ------'
     write(io_unit,*) ' '

     i_folder = QE_folder_sym(ikpt)
     G        = sym_G(:,ikpt)
     kQE      = kpoints_QE(:,i_folder)
     i_sym    = symlink(ikpt,1)

     write(io_unit,60) 'i_folder  = ',i_folder
     write(io_unit,75) 'kpoint_QE = (',kQE(1),',',kQE(2),',',kQE(3),')'
     write(io_unit,350)'sym_G     = (',G(1),',',G(2),',', G(3),')'
     write(io_unit,60) 'i_sym     = ',i_sym

     if ( symlink(ikpt,2) == 1) then
          TR_sign =  -1.0_dp
          write(io_unit,210) 'TR        =      true'
     else 
          TR_sign =   1.0_dp
          write(io_unit,210) 'TR        =      false'
     end if
     
     write(io_unit,*) ' '
     write(io_unit,75) 'ftau      = (',ftau(1,i_sym),',',ftau(2,i_sym),',',ftau(3,i_sym),')'
     write(io_unit,*) ' '

     ! perform rotation
     ! REMEMBER the convention in crystal coordinates:
     !  k_rot(i) = sum_j s(i,j) k(j)

     inv   = inverse_indices(i_sym) 
     k_rot = matmul(s(:,:,inv),kQE)

     k_rot =  TR_sign*k_rot

     write(io_unit,75) 'T * S * kQE = (',  &
         k_rot(1),',', k_rot(2),',',k_rot(3),')'

     dk    = kMP-k_rot+TR_sign*dble(G)
     norm  = sqrt(dot_product(dk,dk))

     write(io_unit,300) '| [kMP + T*G] - T * S * kQE | = ', norm

     if ( norm <  eps_8 ) then
        write(io_unit,310) 'PASSED'
     else 
        write(io_unit,310) 'FAILED' 
        test = .false.
     end if

  write(io_unit,*) '----------------------------------------------------------------'
  end do
  write(io_unit,*)    '========================================'
  if (test) then 
     write(io_unit,*) '========  TEST IS PASSED     ==========='
  else 
     write(io_unit,*) '========  TEST IS FAILED     ==========='
  end if
  write(io_unit,*)    '========================================'

  close(io_unit)


  50  format(a,I6)
  60  format(10x,a,I6)
  75  format(10x,a,F10.4,a,F10.4,a,F10.4,a)
  100 format(a,I4,a,F10.4,a,F10.4,a,F10.4,a)
  200 format(a,I4,a)
  210 format(10x,a)
  300 format(20x,a,F10.4)
  310 format(20x,a)
  350 format(10x,a,I4,a,I4,a,I4,a)
  400 format(a,F10.4)
  450 format(a,E10.4,a)

  end subroutine test_k_symmetry

  subroutine test_point_group_symmetries()
  !----------------------------------------------------------------------------!
  ! This subroutine tests the symmetry operations read from QE.
  !
  !     The conventions are:
  !
  !             r' = R * ( r + f)
  !
  !     cartesian units:                        crystal units
  !     -----------------                       --------------
  !     r'_mu = {R}_{mu,nu} r_nu            X'_i = sum_{j}S(R)_{ji} X_j
  !
  !     where S(R)_{ij} =  a_i * R^{-1} * b_j
  !----------------------------------------------------------------------------!
  use intw_reading,   only: at, bg,  nsym, s, ftau, ityp, tau, nat,    &
                            ntyp, ityp, atom_labels
  use intw_useful_constants
  
  implicit none

  integer        ::  io_unit

  integer        ::  isym,  jsym, iat, jat
  integer        ::  jat_cart, jat_cryst

  integer        ::  i, j, l,    mu,  nu

  integer        ::  inv


  integer        ::  inverse_indices(nsym) 

  real(dp)       ::  tau_cryst(3,nat)

  real(dp)       ::  rot_cart(3,3,nsym)
  real(dp)       ::  S_mat(3,3,nsym)

  real(dp)       ::  S_mat_transpose(3,3)
  real(dp)       ::  A_mat_transpose(3,3)
  real(dp)       ::  B_mat_transpose(3,3)

  real(dp)       ::  ftau_cart(3,nsym)

  real(dp)       ::  rot1(3,3), rot2(3,3), prod(3,3)

  real(dp)       ::  norm
  real(dp)       ::  r1_cart(3), r2_cart(3), dr(3), norm_cart
  real(dp)       ::  r1_cryst(3), r2_cryst(3),  norm_cryst

  real(dp)       ::  n_vec(3)
  real(dp)       ::  eps

  logical        ::  found_cart, found_cryst, found
  logical        ::  pass ,  pass_same
  character(*),parameter :: log_file = 'test_point_group_symmetries.test'
  ! Convert ionic positions to crystal coordinates


  eps       = eps_8

  pass_same = .true.

  tau_cryst =  zero 

  io_unit = find_free_unit() 
  open(unit=io_unit,file=log_file,status='unknown')


  ! Basic output and test of the A and B matrices
  write(io_unit,10) '!============================================================='
  write(io_unit,10) '!             test_point_group_symmetries.test                '
  write(io_unit,10) '!           ------------------------------------              '
  write(io_unit,10) '!  This file contains various tests pertaining to the use of  '
  write(io_unit,10) '!  rotational symmetry. Quantum Espresso implements rotations '
  write(io_unit,10) '!  in sometimes confusing ways, so it is vital to check that  '
  write(io_unit,10) '!  INTW uses symmetries exactly in the right way.             '
  write(io_unit,10) '!============================================================='


  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) '!        Testing the A and B matrices            '
  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) ''
  write(io_unit,10) '! A = 1/a (  |   |   |  ),  B = a/2pi (  |   |   |  )  '
  write(io_unit,10) '!         ( a_1 a_2 a_3 )             ( b_1 b_2 b_3 )  '
  write(io_unit,10) '!         (  |   |   |  )             (  |   |   |  )  '




  write(io_unit,12) ' A =  ( ',at(1,1),',',at(1,2),',',at(1,3),')'
  write(io_unit,12) '      ( ',at(2,1),',',at(2,2),',',at(2,3),')'
  write(io_unit,12) '      ( ',at(3,1),',',at(3,2),',',at(3,3),')'
  write(io_unit,10) ''
  write(io_unit,12) ' B =  ( ',bg(1,1),',',bg(1,2),',',bg(1,3),')'
  write(io_unit,12) '      ( ',bg(2,1),',',bg(2,2),',',bg(2,3),')'
  write(io_unit,12) '      ( ',bg(3,1),',',bg(3,2),',',bg(3,3),')'

  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) '!        Testing the rotation matrices           '
  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) '!   define R a rotation, then                    '
  write(io_unit,10) '!         S(R) = A^T * R^T * B                   ' 
  write(io_unit,10) '!           R  = A   * S^T * B^T                 '
  write(io_unit,10) ''

  A_mat_transpose(:,:) = transpose(at)
  B_mat_transpose(:,:) = transpose(bg)

  do isym =1, nsym
  	S_mat(:,:,isym) = dble(s(:,:,isym)) 
  end do

  do isym =1, nsym

  	S_mat_transpose(:,:)= transpose(S_mat(:,:,isym))

	rot_cart(:,:,isym)  =    matmul(at(:,:), &
				 matmul(S_mat_transpose(:,:),B_mat_transpose(:,:)))


        ftau_cart(:,isym)   =    matmul(at(:,:),ftau(:,isym))

        write(io_unit,14) 'isym = ',isym
        write(io_unit,10) ''

        write(io_unit,16) 'cryst. s(',isym,') = (',       &
			             s(1,1,isym),',',     &
                                     s(1,2,isym),',',     &
                                     s(1,3,isym),')',     &
                          '  f = (',ftau(1,isym),')'

        write(io_unit,18) '(',                            &
			             s(2,1,isym),',',     &
                                     s(2,2,isym),',',     &
                                     s(2,3,isym),')',     &
                          '      (',ftau(2,isym),')'
        write(io_unit,18) '(',                            &
			             s(3,1,isym),',',     &
                                     s(3,2,isym),',',     &
                                     s(3,3,isym),')',     &
                          '      (',ftau(3,isym),')'
        write(io_unit,10) ''
        write(io_unit,20) 'cart.  R(',isym,') = (',       &
			      rot_cart(1,1,isym),',',     &
                              rot_cart(1,2,isym),',',     &
                              rot_cart(1,3,isym),')',     &
                          '  f = (',ftau_cart(1,isym),')'
        write(io_unit,22) '(',       &
			      rot_cart(2,1,isym),',',     &
                              rot_cart(2,2,isym),',',     &
                              rot_cart(2,3,isym),')',     &
                          '      (',ftau_cart(2,isym),')'
        write(io_unit,22) '(',       &
			      rot_cart(3,1,isym),',',     &
                              rot_cart(3,2,isym),',',     &
                              rot_cart(3,3,isym),')',     &
                          '      (',ftau_cart(3,isym),')'
        write(io_unit,10) ''
        write(io_unit,10) ''



  end do

  ! Find the index of the inverse matrix
  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) '!        Find the inverse of each symmetry       '
  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) ''

  write(io_unit,10) '!      isym      inverse isym'

  inverse_indices(:) = 0

  do isym =1, nsym

      rot1(:,:) = rot_cart(:,:,isym)

      found = .false.

      do jsym = 1, nsym
         rot2(:,:) = rot_cart(:,:,jsym)
         ! set the product to minus 1
         prod      = zero
         prod(1,1) = -1.0_dp
         prod(2,2) = -1.0_dp
         prod(3,3) = -1.0_dp

	 prod(:,:) = prod(:,:) + matmul(rot1,rot2)

         norm      = zero 

         do i=1,3
           do j=1,3
             norm = norm + prod(i,j)**2
           end do
         end do
         
         norm = sqrt(norm)

         if (norm < eps_8) then 
              inverse_indices(isym) = jsym
              found = .true.
              exit
         end if

      end do 

      if (found) then
        write(io_unit,24) isym,inverse_indices(isym) 
      else
        write(io_unit,26)         &
        'ERROR: no inverse matrix found for operation isym = ',isym
      end if 

  end do 


  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) '!        Testing that the point group operations '
  write(io_unit,10) '!        leave the ionic positions invariant.    '
  write(io_unit,10) '!------------------------------------------------'
  write(io_unit,10) '!'
  write(io_unit,10) '!  By looking at the QE code, it appears that    '
  write(io_unit,10) '!  the point group operations are defined as     '
  write(io_unit,10) '!'
  write(io_unit,10) '!     tau_j = R*(tau_i)-f (cartesian coordinates)'
  write(io_unit,10) '!'
  write(io_unit,10) '!     where tau_i, tau_j are ionic positions and '
  write(io_unit,10) '!     { R | f } is a point group operation.      '
  write(io_unit,10) '!  Note that this is NOT what the QE paper says! '
  write(io_unit,10) '!'
  write(io_unit,10) '! Ionic positions (CARTESIAN alat units)'
  write(io_unit,10) '! --------------------------------------'

  do iat = 1, nat
     write(io_unit,28) ' ion ',iat,' (',atom_labels(ityp(iat)),')',         &
                       ' tau = (',tau(1,iat),',',tau(2,iat),',',tau(3,iat),')'
  end do
  write(io_unit,10) ''

  ! find crystal coordinates

  write(io_unit,10) '! Ionic positions (CRYSTAL units, inside unit cube)'
  write(io_unit,10) '! -------------------------------------------------'
  do iat = 1, nat

    do i =1, 3 ! crystal index
        tau_cryst(i ,iat) = dot_product(bg(:,i),tau(:,iat))

      ! impose that the positions be in a box [0,1] x [0,1] x [0,1] 
        tau_cryst(i,iat) = tau_cryst(i,iat) - floor(tau_cryst(i,iat) )
    end do

     write(io_unit,28) ' ion ',iat,' (',atom_labels(ityp(iat)),')',         &
     ' tau = (',tau_cryst(1,iat),',',tau_cryst(2,iat),',',tau_cryst(3,iat),')'

  end do
  


  write(io_unit,10) ' -------------------------------------------------'
  write(io_unit,10) ' ! ion #    under operation #   goes to ion #     '
  write(io_unit,10) ' !                            (cart)    (cryst)   '
  write(io_unit,10) ' -------------------------------------------------'
  
  pass = .true.
  ! impose every symmetry, on every atom, in CARTESIAN coordinates

  do iat = 1, nat
     write(io_unit,30) iat

     r1_cart(:)  = tau(:,iat)
     r1_cryst(:) = tau_cryst(:,iat)

     do isym =1, nsym
        inv   = inverse_indices(isym) 


! CHECK
        !r2_cart(:) = matmul(rot_cart(:,:,inv),r1_cart(:)+ftau_cart(:,isym))

        ! This form is consistent with the crystal basis below
        r2_cart(:) = matmul(rot_cart(:,:,isym),r1_cart(:))-ftau_cart(:,isym)


	S_mat_transpose(:,:) = transpose(S_mat(:,:,isym))

        ! I think this is the correct way, consistent with Quantum Espresso.
        ! See the QE subroutines checkallsym, checksym and eqvect.
        r2_cryst(:) = matmul(S_mat_transpose, r1_cryst(:))-ftau(:,isym)

        ! bring back r2_cryst insize unit box
	r2_cryst(1) = r2_cryst(1) - floor(r2_cryst(1))
	r2_cryst(2) = r2_cryst(2) - floor(r2_cryst(2))
	r2_cryst(3) = r2_cryst(3) - floor(r2_cryst(3))

        found_cart  = .false.
        found_cryst = .false.

        ! find the equivalent atom by checking that, 
        ! for one atom in the basis, the difference between
        ! r2 and the position of the atom is a lattice vector.

        ! cartesian
        do jat = 1,nat
            dr(:) = r2_cart(:) - tau(:,jat)
            ! check that this is a lattice vector

            ! the matrix bg contains the reciprocal basis as columns.
            ! Are these coefficients integers?
            n_vec(1) = dot_product(bg(:,1),dr(:))
            n_vec(2) = dot_product(bg(:,2),dr(:))
            n_vec(3) = dot_product(bg(:,3),dr(:))

            n_vec(1) = n_vec(1) - nint(n_vec(1)) 
            n_vec(2) = n_vec(2) - nint(n_vec(2)) 
            n_vec(3) = n_vec(3) - nint(n_vec(3)) 

            norm     = sqrt(dot_product(n_vec,n_vec))

            if (norm < eps) then
                found_cart = .true.
		jat_cart   = jat
                exit
            end if

        end do

        ! crystal
        do jat = 1,nat

	    norm = sqrt(dot_product(r2_cryst(:)-tau_cryst(:,jat), &
	                            r2_cryst(:)-tau_cryst(:,jat)))
            if (norm < eps) then
                found_cryst= .true.
		jat_cryst  = jat
                exit
            end if
        end do


        if (.not. found_cart) then
                pass = .false.
                write(io_unit,'(a,I3,a,I3,a)')    &
                'FAILED to find image of ion ',iat,' under symmetry ',isym,' in CARTESIAN coordinates'
		jat_cart = 0
        end if
        if (.not. found_cryst) then
                pass = .false.
                write(io_unit,'(a,I3,a,I3,a)')    &
                'FAILED to find image of ion ',iat,' under symmetry ',isym,' in CRYSTAL coordinates'
		jat_cryst = 0
        end if


        write(io_unit,32) isym,jat_cart,jat_cryst

	if (jat_cart /= jat_cryst) then
		pass      = .false.
		pass_same = .false.
	end if
     end do 

  end do 

  if ( .not. pass_same ) then
  	write(io_unit,'(a)')  'Some atoms are not taken to the same image '
  	write(io_unit,'(a)')  'in crystal and cartesian coordinates!'
  end if
  write(io_unit,'(a)')  '============================='
  if ( pass ) then
    write(io_unit,'(a)')'=====    TEST PASSED    ====='
  else
    write(io_unit,'(a)')'=====    TEST FAILED    ====='
  end if
  write(io_unit,'(a)')  '============================='


  close(io_unit) 


10 format(A)
12 format(5X,A,F8.4,A,F8.4,A,F8.4,A)
14 format(5X,A,I4)
16 format(5X,A,I4,A,I8,A,I8,a,I8,A,5x,A,F8.4,A)  
18 format(22X,A,I8,A,I8,a,I8,a,5x,A,F8.4,A)  
20 format(5X,A,I4,A,F8.4,A,F8.4,a,F8.4,A,5x,A,F8.4,A)  
22 format(22X,A,F8.4,A,F8.4,a,F8.4,a,5x,A,F8.4,A)  
24 format(I8,6X,I8)
26 format(A,I8)
28 format(8x,A,I2,A,A,A,8x,A,F8.4,A,F8.4,A,F8.4,A)
30 format(5x,I4)
32 format(16X,I4,10X,I4,5X,I4)

  end subroutine test_point_group_symmetries


  subroutine test_matrix_elements_orthogonality( nk_1,nk_2,nk_3,method)
  !----------------------------------------------------------------------------!
  ! This subroutine tests the computation of plane wave matrix elements;
  ! it makes sure that < psi_n1 k | psi_n2 k> = delta_{n1n2} 
  ! mehtod is either "FFT" or "CONVOLUTION"
  !----------------------------------------------------------------------------!
  use intw_matrix_elements, only : get_plane_wave_matrix_element_convolution, &
                                   get_plane_wave_matrix_element_FFT
  use intw_reading

  implicit none

  character(*)   ::  method
  character(256) ::  filename 

  integer        ::  io_unit
  integer        ::  nk_1,nk_2,nk_3 

  integer        ::  ikpt
  integer        ::  nb1,  nb2 

  integer        ::  G(3)

  integer        :: list_iG(nG_max)
  complex(dp)    :: wfc(nG_max,num_bands)
!  complex(dp)    :: wfc(nG_max,nbands)
  real(dp)       :: QE_eig(num_bands)
!  real(dp)       :: QE_eig(nbands)

  real(dp)       :: eps 

  complex(dp)    :: pw_mat_el(num_bands,num_bands), test
!  complex(dp)    :: pw_mat_el(nbands,nbands), test
  
  logical        :: passed, global_passed

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio


  nbands_loc=num_bands

  eps = eps_8

  G = 0
  if (trim(method) .eq. 'FFT' ) then
        filename  = 'orthogonality.FFT.test'
  else if (trim(method) .eq. 'CONVOLUTION' ) then
        filename  = 'orthogonality.conv.test'
  else
        write(*,*) 'ERROR: this method does not exist'
        return
  end if
  
  global_passed = .true.

  io_unit = find_free_unit()
  open(unit=io_unit,file=trim(filename),status='unknown')

  write(io_unit,*) '###########################################'
  write(io_unit,*) '##  Test of matrix element orthogonality ##'
  write(io_unit,*) '##      small number = 1e-8              ##'
  write(io_unit,*) '##  the test is if                       ##'
  write(io_unit,*) '##  <psi_n1k| psi_n2k> = delta_{n1,n2}   ##'
  write(io_unit,*) '##  COMPUTATION DONE WITH '//trim(method)// '##'
  write(io_unit,*) '###########################################'
  write(io_unit,*) ''
  

  !loop on all points
  do ikpt = 1, nk_1*nk_2*nk_3
    write(*,*) 'ikpt = ',ikpt
    write(io_unit,*) '#----------------------------------------------------------------'
    write(io_unit,'(a,I4)') '# ikpt       = ',ikpt

    call get_K_folder_data(ikpt,list_iG,wfc,QE_eig)

    if (method .eq. 'CONVOLUTION') then 
            call get_plane_wave_matrix_element_convolution      &
                        (G,list_iG,list_iG, wfc,wfc,pw_mat_el)
    else 
            call get_plane_wave_matrix_element_FFT              &
                        (G,list_iG,list_iG, wfc,wfc,pw_mat_el)
    end if 

    passed        = .true.
    do nb1 = 1,nbands_loc
!    do nb1 = 1,nbands
      do nb2 = 1,nbands_loc
!      do nb2 = 1,nbands
           if ( nb1 .eq. nb2) then
              test = pw_mat_el(nb1,nb2)-cmplx_1
           else
              test = pw_mat_el(nb1,nb2)
           end if
           
           if ( abs(real(test)) > eps .or. abs(aimag(test)) > eps ) then
                passed        = .false.
                global_passed = .false.
           end if
      end do
    end do
    
    if (passed) then
        write(io_unit,*) '     test PASSED' 
    else 
        write(io_unit,*) '     test FAILED' 
    end if

    write(io_unit,*) '#----------------------------------------------------------------'

  end do

  write(io_unit,*)    '========================================'
  if (global_passed) then 
     write(io_unit,*) '========  TEST IS PASSED     ==========='
  else 
     write(io_unit,*) '========  TEST IS FAILED     ==========='
  end if
  write(io_unit,*)    '========================================'


  close(unit=io_unit)

  end subroutine test_matrix_elements_orthogonality


  subroutine test_CONVOLUTION ( G,list_iG_1,list_iG_2, wfc_1,wfc_2)
  !----------------------------------------------------------------------------!
  ! This subroutine uses the convolution algorithm to compute matrix elements,
  ! for testing purposes. Only the bands 1, 1 will be considered
  !
  ! NOTE: This test was instrumental in finding a HUGE bug in the original
  !       version of "get_plane_wave_matrix_element_convolution". Once again,
  !       testing is an ESSENTIAL part of software design.
  !----------------------------------------------------------------------------!

  use intw_useful_constants
  use intw_reading
  use intw_utility
  use intw_fft, only: nl, g_fft_map

  implicit none


  integer        :: G(3),    G1(3),  G2(3),  Gprime(3)

  integer        :: iG1,     iG2

  integer        :: nb1,  nb2

  integer        :: i,  j

  integer        :: n1, n2, n3 

  integer        :: io_unit

  integer        :: list_iG_1(nG_max),           list_iG_2(nG_max)

  complex(dp)    :: pw_el_mat

  complex(dp)    :: wfc_1(nG_max,num_bands),       wfc_2(nG_max,num_bands)
!  complex(dp)    :: wfc_1(nG_max,nbands),       wfc_2(nG_max,nbands)

  logical        :: found

  pw_el_mat  = cmplx_0

  io_unit = find_free_unit()
  open(unit=io_unit,file=('conv.test'),status='unknown')

  write(io_unit,*) '##########################################################'
  write(io_unit,*) '# Find the product < n1=1 k1 | e^{-i (G+k2-k1) r | n2=1 k2>   ' 
  write(io_unit,*) '# using the convolution algorithm.'
  write(io_unit,*) '##########################################################'
  write(io_unit,*) ''
  write(io_unit,*) ' G = ',G


  do  i = 1, nG_max
        write(io_unit,*) '#-----------------------------------------------------' 
        write(io_unit,'(a,I6)') ' i = ',i

        ! list_iG is zero-padded at the end. 
        ! When this point is reached, the computation is over.

        iG1 = list_iG_1(i) 

        if (iG1 == 0 ) then  
           write(io_unit,'(8x,a,I6,a)') 'iG1 = ',iG1,', ZERO-PADDED REGION' 
        else

           write(io_unit,'(8x,a,I6)') 'iG1 = ',iG1
           write(io_unit,*) ''

           G1  = gvec(:,iG1) ! make sure you find G1 by using its proper index, iG1 !

           write(io_unit,'(12x,a,I4,a,I4,a,I4,a)') '        G1  = (',G1(1),',',G1(2),',',G1(3),')'

           Gprime = G1+G 

           write(io_unit,'(12x,a,I4,a,I4,a,I4,a)') ' G+G1 = G2  = (',Gprime(1),',',Gprime(2),',',Gprime(3),')'

           ! find if Gprime is in the domain of the second wavefunction
           found = .false.
           do j = 1, nG_max

               iG2 = list_iG_2(j)

               G2  = gvec(:,iG2)
               
               ! if G2 = Gprime, add contribution to product 
               if ( (G2(1) == Gprime(1)) .and.  &
                    (G2(2) == Gprime(2)) .and.  &
                    (G2(3) == Gprime(3)) )   then  
                        
                    found = .true. 
                    write(io_unit,'(a,I6)') ' j = ',j
                    write(io_unit,'(8x,a,I6)') 'iG2 = ',iG2

                        pw_el_mat = pw_el_mat + CONJG(wfc_1(i,1))*wfc_2(j,1)

               else if (.not.found .and. j == nG_max) then 

                    write(io_unit,'(8x,a)') 'G2 not present in list_G2 '

               end if

           end do 

           write(io_unit,'(a,2F18.12)') ' mat_el  = ',pw_el_mat

        end if
  end do

  write(io_unit,*) '#-----------------------------------------------------' 
  write(io_unit,*) ' FINAL VALUE'

  write(io_unit,'(2F18.12)') pw_el_mat

  close(io_unit)

  end subroutine test_CONVOLUTION


  subroutine test_slow_FT ( G,list_iG_1,list_iG_2, wfc_1,wfc_2)
  !----------------------------------------------------------------------------!
  ! This subroutine computes the matrix product using Fourier transforms,
  ! but by a slow, direct computation algorithm instead of using FFT. 
  !----------------------------------------------------------------------------!

  use intw_useful_constants
  use intw_reading
  use intw_utility

  implicit none


  integer        :: G(3),    G1(3),  G2(3),  Gprime(3)

  integer        :: iG1,     iG2

  integer        :: nb1,     nb2,    nb

  integer        :: ir,  jr, kr,     i


  integer        :: io_unit

  integer        :: list_iG_1(nG_max),           list_iG_2(nG_max)

  complex(dp)    :: wfc_1(nG_max,num_bands),       wfc_2(nG_max,num_bands)
!  complex(dp)    :: wfc_1(nG_max,nbands),       wfc_2(nG_max,nbands)

  complex(dp)    :: wfc_1_r(nr1,nr2,nr3,num_bands),       wfc_2_r(nr1,nr2,nr3,num_bands)
!  complex(dp)    :: wfc_1_r(nr1,nr2,nr3,nbands),       wfc_2_r(nr1,nr2,nr3,nbands)

  complex(dp)    :: prod_r
  complex(dp)    :: phase1, phase2,  phase

  complex(dp)    :: mat_el(num_bands,num_bands)
!  complex(dp)    :: mat_el(nbands,nbands) 

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio


  nbands_loc=num_bands

  ! Transform to real space
  wfc_1_r  =  cmplx_0
  wfc_2_r  =  cmplx_0 

  do ir=1, nr1
    do jr=1, nr2
      do kr=1, nr3

            do i = 1, nG_max

                iG1 = list_iG_1(i)
                if (iG1 == 0 ) exit                

                G1  = gvec(:,iG1)

                phase1 = exp( tpi*cmplx_i   *  (       &
                     cmplx_1* G1(1) * (ir-1)/nr1  +  &
                     cmplx_1* G1(2) * (jr-1)/nr2  +  &
                     cmplx_1* G1(3) * (kr-1)/nr3     ))

               do nb = 1, nbands_loc
!               do nb = 1, nbands
                 wfc_1_r(ir,jr,kr,nb)  =  wfc_1_r(ir,jr,kr,nb)    + phase1*wfc_1(i,nb) 
               end do

            end do

            do i = 1, nG_max

                iG2 = list_iG_2(i)
                if (iG2 == 0 ) exit                

                G2  = gvec(:,iG2)

                phase2 = exp( tpi*cmplx_i   *  (       &
                                cmplx_1* G2(1) * (ir-1)/nr1  +  &
                                cmplx_1* G2(2) * (jr-1)/nr2  +  &
                                cmplx_1* G2(3) * (kr-1)/nr3     ))
                

               do nb = 1, nbands_loc
!               do nb = 1, nbands
                 wfc_2_r(ir,jr,kr,nb)  =  wfc_2_r(ir,jr,kr,nb)    + phase2*wfc_2(i,nb) 
               end do

            end do

      end do
    end do
  end do


  ! product and inv FT
  mat_el = cmplx_0

  do ir=1, nr1
    do jr=1, nr2
      do kr=1, nr3

        phase = exp( -tpi*cmplx_i*  (       &
                     cmplx_1* G(1) * (ir-1)/nr1  +  &
                     cmplx_1* G(2) * (jr-1)/nr2  +  &
                     cmplx_1* G(3) * (kr-1)/nr3     )) 

        do nb1 = 1, nbands_loc
!        do nb1 = 1, nbands
          do nb2 = 1, nbands_loc
!          do nb2 = 1, nbands

            prod_r  =  conjg(wfc_1_r(ir,jr,kr,nb1))*wfc_2_r(ir,jr,kr,nb2)

            mat_el(nb1,nb2) = mat_el(nb1,nb2) + phase*prod_r 

          end do
        end do

      end do
    end do
  end do

  mat_el(:,:) = mat_el(:,:)/(nr1*nr2*nr3)
  
  io_unit = find_free_unit()
  open(unit=io_unit,file=('slow_FT.test'),status='unknown')
  write(io_unit,*) '##########################################################'
  write(io_unit,*) '# Find the product < n1 k1 | e^{-i (G+k2-k1) r | n2 k2>   ' 
  write(io_unit,*) '# using fourier transforms, but not FFT.                  ' 
  write(io_unit,*) '##########################################################'
  write(io_unit,*) ''

  do nb2 = 1,nbands_loc 
!  do nb2 = 1,nbands
    do nb1 = 1,nbands_loc 
!    do nb1 = 1,nbands
      write(io_unit,'(2F18.12)') mat_el(nb1,nb2)
    end do
  end do


  close(io_unit)

  end subroutine test_slow_FT


  subroutine test_get_psi_k()
  !----------------------------------------------------------------------------!
  ! This subroutine compares the wfc computed directly from QE and wfc
  ! obtained by rotation from the appropriate k_irr. The subroutine thus builds
  ! The matrix 
  !             U_{mn} = < psi_m (QE) | psi_n (irr rotated) >
  !
  ! This matrix should be unitary!
  !----------------------------------------------------------------------------!

  implicit none

  integer        :: iG
  integer        :: n_bnd, m_bnd, p
  integer        :: sub_block
  integer        :: imin, imax
  integer        :: ikpt, ipol
  integer        :: G(3)

  integer        :: io_unit

  integer        :: list_iG_1(nG_max), list_iG_2(nG_max)

  complex(dp)    :: wfc_k_1(nG_max,num_bands,nspin), wfc_k_2(nG_max,num_bands,nspin)
!  complex(dp)    :: wfc_k_1(nG_max,nbands,nspin), wfc_k_2(nG_max,nbands,nspin)
  complex(dp)    :: wfc_work(nG_max,num_bands,nspin)
!  complex(dp)    :: wfc_work(nG_max,nbands,nspin)

  integer        :: permutations(num_bands)
!  integer        :: permutations(nbands)

  real(dp)       :: eps

  real(dp)       :: QE_eig_1(num_bands),  QE_eig_2(num_bands), d_eig(num_bands)
!  real(dp)       :: QE_eig_1(nbands),  QE_eig_2(nbands), d_eig(nbands)

  real(dp)       :: norm

  real(dp)       :: ftau_sym(3)

  real(dp)       :: test_z

  complex(dp)    :: U(num_bands,num_bands), clean_U(num_bands,num_bands)
!  complex(dp)    :: U(nbands,nbands), clean_U(nbands,nbands)
  integer        :: number_of_off_diagonal_elements

  complex(dp),allocatable    :: sub_U(:,:), sub_U_d(:,:), sub_id(:,:)

  complex(dp)    :: pw_mat_el(num_bands,num_bands,nspin,nspin)
!  complex(dp)    :: pw_mat_el(nbands,nbands,nspin,nspin)

  complex(dp)    :: identity(num_bands,num_bands)
!  complex(dp)    :: identity(nbands,nbands)

  logical        :: pass_G, pass_u, warning_u
  logical        :: global_pass_G, global_pass_u, global_pass_e, global_warning_u

  logical        :: use_IBZ_true, use_IBZ_false

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio


  nbands_loc=num_bands

  eps = eps_6

  global_pass_G    = .true.
  global_pass_u    = .true.
  global_pass_e    = .true.
  global_warning_u = .false.

  G   = 0

  use_IBZ_true    = .true.
  use_IBZ_false   = .false.

  ! open a test file
  io_unit = find_free_unit()
  open(unit=io_unit,file=trim('get_psi_k.test'))

  write(io_unit,10) '#=========================================================='
  write(io_unit,10) '# compare wfc_k from file and rotated from k_irr' 
  write(io_unit,10) '#=========================================================='
  write(io_unit,10) '#                                                          '
  write(io_unit,10) '# This test computes U_{mn} = < psi_m(QE) | psi_n(rot) >   '
  write(io_unit,10) '# where psi_m(QE) is read directly from file, and          '   
  write(io_unit,10) '# psi_n(rot) is obtained from the rotation of an irr       '   
  write(io_unit,10) '# k point.                                                 '
  write(io_unit,10) '#                                                          '
  write(io_unit,10) '# The test constitutes in checking that                    '
  write(io_unit,10) '#       - list_G is the same for both  wavefunctions       '
  write(io_unit,10) '#       - U is unitary and leaves eigenvalues unchanged    '
  write(io_unit,10) '#                                                          '
  write(io_unit,10) '# however:                                                 '
  write(io_unit,10) '# - Quantum Espresso cuts off the number of bands in an    '
  write(io_unit,10) '#   arbitrary, user defined way. That is to say it does    '
  write(io_unit,10) '#   not consider symmetry and may throw away degenerate    '
  write(io_unit,10) '#   states at certain k-points. Thus, if the test above    '
  write(io_unit,10) '#   fails, only a *warning* will be issued if the fail-    '
  write(io_unit,10) '#   ure is related to the last eigenvalue.                 '
  write(io_unit,10) '#                                                          '
  write(io_unit,10) '#=========================================================='

  ! loop on all k points in the coarse mesh
  do ikpt = 1,nk1*nk2*nk3
        ! get the wavefunctions 
        call get_psi_k(ikpt,use_IBZ_true ,list_iG_1,wfc_k_1,QE_eig_1)
        call get_psi_k(ikpt,use_IBZ_false,list_iG_2,wfc_k_2,QE_eig_2)

        ! sort the eigenvalues and corresponding wave functions

        call HPSORT_real(nbands_loc,QE_eig_1,permutations)
!        call HPSORT_real(nbands,QE_eig_1,permutations)
        do n_bnd = 1,nbands_loc
!        do n_bnd = 1,nbands
                p = permutations(n_bnd)
                do ipol = 1, nspin
                   wfc_work(:,n_bnd,ipol) = wfc_k_1(:,p,ipol)
		end do
        end do
        wfc_k_1  = wfc_work

        call HPSORT_real(nbands_loc,QE_eig_2,permutations)
!        call HPSORT_real(nbands,QE_eig_2,permutations)
        do n_bnd = 1,nbands_loc
!        do n_bnd = 1,nbands
                p = permutations(n_bnd)
                do ipol = 1, nspin
                   wfc_work(:,n_bnd,ipol) = wfc_k_2(:,p,ipol)
                end do
        end do
        wfc_k_2  = wfc_work


        call get_plane_wave_matrix_element_convolution      &
                        (G,list_iG_1,list_iG_2, wfc_k_1,wfc_k_2,pw_mat_el)
  
        U(:,:)  = cmplx_0

	do ipol = 1,nspin
	   U(:,:) = U(:,:) + pw_mat_el(:,:,ipol,ipol) 
        end do

        clean_U = U       

        number_of_off_diagonal_elements  = nbands_loc*nbands_loc
!        number_of_off_diagonal_elements  = nbands*nbands

        write(io_unit,10) '----------------------------------------------------'
        write(io_unit,20) 'ikpt = ',ikpt,', i_folder =',QE_folder_sym(ikpt)
        write(io_unit,10) ''
        if (symlink(ikpt,2) == 1) then
                write(io_unit,10) ' * TR = True'
        else
                write(io_unit,10) ' * TR = False'
        end if
        ! test G
        pass_G     = .true.
        do iG=1,nG_max
                if ( list_iG_1(iG)- list_iG_2(iG) /= 0) then
                        pass_G         = .false.
                        global_pass_G  = .false.
                end if
        end do

        if ( pass_G ) then 
           write(io_unit,10) '  * all G vectors are the same: PASSED'
        else
           write(io_unit,10) '  * different G vectors! : FAILED'
           write(io_unit,10) '  |  iG (rot)   |   iG(QE)     |  difference | '
           
           do iG=1,nG_max
                write(io_unit,15) list_iG_1(iG), list_iG_2(iG),  &
                                  list_iG_1(iG)-list_iG_2(iG)

                if (list_iG_1(iG) == 0 .and. list_iG_2(iG) == 0 ) exit
           end do

        end if

        write(io_unit,10) ''

        ! check eigenvalues

        d_eig     = QE_eig_1-QE_eig_2
        norm      = sqrt(dot_product(d_eig,d_eig))
        if (norm < eps) then
                write (io_unit,11) ' * | eigenvalues_ikpt - egenvalues_ifolder | = ',norm,' eV :: PASSED'
        else
                write (io_unit,11) ' * | eigenvalues_ikpt - egenvalues_ifolder | = ',norm,' eV :: FAILED'
                global_pass_e  = .false.
        end if
       
        ! Check U 

        identity = cmplx_0
        do n_bnd = 1, nbands_loc
!        do n_bnd = 1, nbands
                identity(n_bnd,n_bnd) = cmplx_1
        end do 

        pass_u     = .true.
        warning_u  = .false.


        imin       = 1
        imax       = imin

        sub_block  = 0

        write(io_unit,10) '   sub-block        (imin, imax) '

        do  
           ! find the block 
           do while (  abs(QE_eig_1(imin)-QE_eig_1(imax)) < eps ) 
                imax = imax +1 
                if ( imax == nbands_loc+1) exit
!                if ( imax == nbands+1) exit
           end do
           imax     = imax - 1 


           sub_block = sub_block + 1

           write(io_unit,21) sub_block,imin,imax


           allocate(sub_U  (imax-imin+1,imax-imin+1))
           allocate(sub_U_d(imax-imin+1,imax-imin+1))
           allocate(sub_id(imax-imin+1,imax-imin+1))

           sub_U   = U(imin:imax,imin:imax)
           sub_U_d = conjg(transpose(sub_U))
           sub_id  = identity(imin:imax,imin:imax)

           test_z  = sqrt(sum(abs((matmul(sub_U,sub_U_d)-sub_id)**2)))

!           test_z    = sqrt(sum(abs(matmul(                        &
!                                        U(imin:imax,imin:imax),    &
!                        conjg(transpose(U(imin:imax,imin:imax))))  &
!                           -     identity(imin:imax,imin:imax))**2))
           deallocate(sub_U  )
           deallocate(sub_U_d)
           deallocate(sub_id )


           if (abs(test_z) > eps ) then
              if (  abs(QE_eig_1(imax)-QE_eig_1(nbands_loc)) > eps )  then
!              if (  abs(QE_eig_1(imax)-QE_eig_1(nbands)) > eps )  then
                  pass_u        = .false.
                  global_pass_u = .false.
                  write(io_unit,12) '  * sub-block ',sub_block,       &
                                    ' of U matrix is NOT unitary!: FAILED'
               else
                  warning_u        = .true.
                  global_warning_u = .true.
                  write(io_unit,12) '  * sub-block ',sub_block,   &
                                ' of U matrix is NOT unitary.'//  &
                                ' However, non unitarity is in the last block: WARNING'
               end if
           end if

           clean_U(imin:imax,imin:imax) = cmplx_0

           number_of_off_diagonal_elements  = number_of_off_diagonal_elements  &
					      -  (imax-imin+1)**2

           imin      = imax+1
           imax      = imin
           if ( imax == nbands_loc+1) exit
!           if ( imax == nbands+1) exit
        end do

        ! check off-diagonal blocks; compute the average error

        test_z  = sqrt(sum(abs(matmul(clean_U,conjg(transpose(clean_U))))))/ &
			sqrt(dble(number_of_off_diagonal_elements)) 
                                        
        if (abs(test_z) > eps ) then
                pass_u        = .false.
                global_pass_u = .false.
                write(io_unit,10) '  * off diagonal U matrix is NOT zero! : FAILED'
        end if

        if ( pass_u .and. .not. warning_u) then 
           write(io_unit,10) '  * U matrix is unitary: PASSED'
        else if ( pass_u .and. warning_u) then 
           write(io_unit,10) '  * U matrix is not unitary in last block: WARNING'
        else
           write(io_unit,10) '  * U matrix is NOT unitary!: FAILED'
        end if

        write(io_unit,10) ' eigenvalues = '
        write(io_unit,35)  QE_eig_1

        write(io_unit,10) ' U ='
        do n_bnd = 1,nbands_loc
!        do n_bnd = 1,nbands
               write(io_unit,30)  U(n_bnd,:)
        end do

        write(io_unit,10) ''
       
  end do 

  write(io_unit,10) '#=========================================================='
  write(io_unit,10) '# SUMMARY:'
  if ( global_pass_G) then
        write(io_unit,10) '# * ALL G VECTORS THE SAME       : PASSED' 
  else
        write(io_unit,10) '# * SOME G VECTORS DIFFERENT!    : FAILED!' 
  end if

  if ( global_pass_e) then
        write(io_unit,10) '# * ALL EIGENVALUES THE SAME     : PASSED' 
  else
        write(io_unit,10) '# * SOME EIGENVALUES DIFFERENT!  : FAILED!' 
  end if

  if ( global_pass_u .and. .not. global_warning_u) then
        write(io_unit,10) '# * ALL U MATRICES UNITARY       : PASSED' 
  else if ( global_pass_u .and. global_warning_u) then
        write(io_unit,10) '# * SOME U MATRICES NOT UNITARY IN LAST BLOCK: WARNING' 
  else
        write(io_unit,10) '# * SOME U MATRICES NOT UNITARY! : FAILED!' 
  end if

  write(io_unit,10) '#=========================================================='

  close(io_unit)

  10 format(A)
  11 format(A,E10.4,A)
  12 format(A,I4,A)
  13 format(A,E10.4)
  15 format(10x3I10)
  20 format(A,I4,A,I4)
  21 format(I6,14X,I4,I4)
  30 format(100(2F8.4,5X))
  35 format(100F8.4)

  end subroutine test_get_psi_k

  subroutine test_generate_nl()
  !----------------------------------------------------------------------------!
  ! This subroutine tests the generate_nl() subroutine 
  !----------------------------------------------------------------------------!
  implicit none

  integer        ::  io_unit
  integer        ::  iG 
  integer        ::  n1, n2, n3 
  integer        ::  switch 
  character(*),parameter :: log_file = 'test_generate_nl.test'

  logical        ::  passed

  switch  = -1  ! singlet-to-triplet index

  io_unit = find_free_unit()
  call create_or_append_to_file(io_unit,log_file)

  write(io_unit,10) "========================================="
  write(io_unit,10) "= Test of generate_nl() subroutine      ="
  write(io_unit,10) "========================================="
  write(io_unit,10) "" 
  write(io_unit,10) "-----------------------------------------"
  write(io_unit,10) "---          Definitions              ---"
  write(io_unit,10) "-----------------------------------------"
  write(io_unit,15) "1) G is in crystal units ==> G = G(1) b1 + G(2) b2 + G(3) b3"
  write(io_unit,15) "" 
  write(io_unit,15) "2) ngm is the total number of G-vectors in gvec."
  write(io_unit,16)  "iG is the index of G in this list; 1 <= iG <= ngm"
  write(io_unit,15) "" 
  write(io_unit,15) "3) G_fft is the periodic image of G in the FFT box."
  write(io_unit,15) "" 
  write(io_unit,15) "4) n1,n2,n3 is the triplet index of G_fft, and iG_fft "
  write(io_unit,16) "the corresponding scalar index."
  write(io_unit,15) "" 
  write(io_unit,15) "5) nr1,nr2,nr3 represent the 3D-FFT mesh."
  write(io_unit,15) "" 
  write(io_unit,15) "Equalities: "
  write(io_unit,15) "------------"
  write(io_unit,15) ""
  write(io_unit,15) ""
  write(io_unit,20) "(nr1, nr2, nr3) = (",nr1,",",nr2,",",nr3,")"
  write(io_unit,25) "nr1 * nr2 * nr3 = ",nr1*nr2*nr3
  write(io_unit,25) "ngm             = ",ngm
  write(io_unit,15) "-> iG_fft  = nl(iG)"
  write(io_unit,15) "-> g_fft_map(n1,n2,n3) =   iG"
  write(io_unit,15) "test if iG = g_fft_map(n1,n2,n3)..."
  write(io_unit,10) ""
  
  !loop on all G vectors 

  passed = .true.

  do iG = 1, ngm 
    call switch_indices(nr1,nr2,nr3,nl(iG),n1,n2,n3,switch)
    if ( .not. g_fft_map(n1,n2,n3) .eq.  iG) then 

        write(io_unit,10) "#-----------------------------------"
        write(io_unit,30) "G = (",  &
                gvec(1,iG)," ,", gvec(2,iG)," ,", gvec(3,iG),")"

        write(io_unit,35) "iG = ",iG," iG_fft  = ",nl(iG), &
                          " n1,n2,n3   =  ", n1,n2,n3

        write(io_unit,40) "g_fft_map = ",g_fft_map(n1,n2,n3),' FAIL !'

        passed = .false.
    end if
  end do

  if ( passed ) then      
    write(io_unit,10) "TEST IS PASSED"
  else 
    write(io_unit,10) "TEST IS FAILED!"
  end if

  write(io_unit,10) "========================================="

  close(unit=io_unit)

  10 format(A) 
  15 format(4X,A) 
  16 format(6X,A) 
  20 format(4X,A,I4,A,I4,A,I4,A)
  25 format(4X,A,I8) 
  30 format(4X,A,I4,A,I4,A,I4,A) 
  35 format(4X,A,I4,A,I4,A,3I4) 
  40 format(4X,A,I4,A)

  end subroutine test_generate_nl

  subroutine test_G_shells(number_G_shells,nG_shell_max)
  !----------------------------------------------------------------------------!
  !     This subroutine checks that the G vectors belong to shells
  !----------------------------------------------------------------------------!

  implicit none
  
  character(*),parameter :: log_file = 'test_G_Shells.test'

  integer      :: number_G_shells, nG_shell_max

  integer      :: nshell
  integer      :: iG

  integer   :: io_unit

  integer   :: iflag, nvec

  real(dp)  :: G_cart(3)

  real(dp)  :: norm2, old_norm2


  io_unit   = find_free_unit()
  call create_or_append_to_file(io_unit,log_file)

  write(io_unit,10) "========================================="
  write(io_unit,10) "= Outputing of G-vectors in shells      ="
  write(io_unit,10) "========================================="
  write(io_unit,20) "number_G_shells   = ",number_G_shells
  write(io_unit,20) "nG_shell_max      = ",nG_shell_max

  write(io_unit,10) "-----------------------------------------"
  write(io_unit,10) "---         shells of G vectors       ---"
  write(io_unit,10) "-----------------------------------------"

  write(io_unit,10) 'shell # 1'
  write(io_unit,30) dble(gvec(:,1))

  nshell       = 1 

  nvec   = 1
  iflag  = 1

  old_norm2 = ZERO

  do iG = 2,ngm

        G_cart = dble(gvec(:,iG))
        call cryst_to_cart (nvec, G_cart, bg, iflag)

        norm2   = G_cart(1)**2+G_cart(2)**2+G_cart(3)**2

        if ( abs(norm2-old_norm2) < eps_8) then
                write(io_unit,30) G_cart 
        else
                nshell =  nshell + 1  
                if (nshell == number_G_shells+1) exit
                write(io_unit,10) '--------------------'
                write(io_unit,15) 'shell #',nshell
                write(io_unit,30) G_cart 
                
                old_norm2 = norm2 

        end if
  end do
  write(io_unit,10) "========================================="


 close(io_unit)

 10 format(A) 
 15 format(A,I4) 
 20 format(4X,A,I4) 
 30 format (4X,'G_cart = (', F8.4,',',F8.4,',',F8.4,')')

 end subroutine test_G_shells

  subroutine empty_lattice_test_eigenvalues(eig_int_ps,list_ps,nk_vec_max,nk_vec)
  !----------------------------------------------------------------------------!
  !     This subroutine computes the empty lattice (ie free electron) 
  !     eigenvalues, for the free electron gas corresponding to the 
  !     system's volume. This is mostly for testing.
  !     The algorithm does not assume that the p-vectors stored in list_ps
  !     are in the 1BZ.
  !----------------------------------------------------------------------------!

  implicit none

  integer     :: nk_vec_max ! the size of the kpoint array
  integer     :: nk_vec     ! the number of meaningful kpoints in the array

  real(dp)    :: eig_int_ps(num_wann,nk_vec_max)
  real(dp)    :: list_ps   (3,nk_vec_max)

  real(dp)    :: list_ks   (3,nk_vec_max)

  integer     :: ikpt,  nbnd
  integer     :: G(3)
  real(dp)    :: p(3)


  real(dp)    :: energy_conversion_to_eV

  energy_conversion_to_eV = Ha_in_eV*2.0_dp*pi**2/alat**2


  ! insure that the vectors are in the 1BZ
  list_ks = modulo(list_ps,ONE)

  eig_int_ps = ZERO

  do nbnd = 1, num_wann
        G = gvec(:,nbnd)
        do ikpt = 1, nk_vec

                ! note that the b vectors are stored as columns in the array
                p = matmul(list_ks(:,ikpt)+dble(G),bg)

                eig_int_ps(nbnd,ikpt) = dot_product(p,p)*energy_conversion_to_eV 
        end do
  end do 

 end subroutine empty_lattice_test_eigenvalues

  subroutine empty_lattice_test_matrix_elements(iG,pw_matrix_elements,       &
                                                        nk_vec,nk_vec_max)
  !----------------------------------------------------------------------------!
  !             THIS SUBROUTINE IS OBSOLETE
  !
  !     This subroutine computes the empty lattice (ie free electron) 
  !     matrix elements, for the free electron gas corresponding to the 
  !     system's volume. This is mostly for testing.
  !----------------------------------------------------------------------------!

  implicit none

  integer     :: iG         ! The index of the external G vector

  complex(dp) :: pw_matrix_elements(num_wann,num_wann,nk_vec_max)

  integer     :: nk_vec_max ! the size of the kpoint array
  integer     :: nk_vec     ! the number of meaningful kpoints in the array


  integer     :: ikpt
  integer     :: nbnd_1, nbnd_2
  integer     :: G(3), G_1(3), G_2(3), dG(3)

  pw_matrix_elements = cmplx_0

  G = gvec(:,iG)

  do nbnd_1 = 1, num_wann
     G_1 = gvec(:,nbnd_1)
     do nbnd_2 = 1, num_wann
        G_2 = gvec(:,nbnd_2)
        
        dG  = G+G_1-G_2

        if (dot_product(dG,dG) == 0) then
                pw_matrix_elements(nbnd_1,nbnd_2,:) = cmplx_1
        end if

     end do
  end do 


 end subroutine empty_lattice_test_matrix_elements


  subroutine empty_lattice_test_fix_interpolated_matrix_elements          &
                (nk_vec_max,nk_vec,iG,list_ks,list_ksq,pw_matrix_elements)
  !------------------------------------------------------------------
  ! 
  ! This subroutine SETS (instead of interpolates) the matrix elements 
  ! for all ks belonging to a given sub-block, for a given G, for
  ! the empty lattice test.
  !
  !------------------------------------------------------------------
  implicit none

  ! input variables
  integer       ::  nk_vec_max, nk_vec, iG

  real(dp)      :: list_ks(3,nk_vec_max), list_ksq(3,nk_vec_max)

  ! output variable
  complex(dp)   ::  pw_matrix_elements(num_wann,num_wann,nk_vec_max)

  ! local variables
  integer       ::  nbnd_1, nbnd_2
  integer       ::  ikpt
  integer       ::  G(3), G_1(3), G_2(3), G_n1(3), G_n2(3), dG(3)

  real(dp)      ::  time1, time2
  ! initialize


!  call get_timing(time1)
  ! Don't initialize the whole array, it takes time!!!
  ! set to zero only the bit which will actually be used.

  pw_matrix_elements(:,:,1:nk_vec)  =  cmplx_0

  G     = gvec(:,iG)

   
   do nbnd_1 =1, num_wann
      G_n1 = gvec(:,nbnd_1)

      do nbnd_2 =1, num_wann
         G_n2 = gvec(:,nbnd_2)

         do ikpt = 1, nk_vec

                G_1 = int(list_ks (:,ikpt) - modulo(list_ks (:,ikpt),ONE))
                G_2 = int(list_ksq(:,ikpt) - modulo(list_ksq(:,ikpt),ONE))

                dG  = G+G_2-G_1+G_n1-G_n2

                if (dot_product(dG,dG) == 0) then
                        pw_matrix_elements(nbnd_1,nbnd_2,ikpt) = cmplx_1
                end if

         end do

      end do
   end do

!  call get_timing(time2)
!  write(*,*) 'time in routine:',time2-time1,' seconds'

  end subroutine empty_lattice_test_fix_interpolated_matrix_elements          

  subroutine test_wfc_orthogonality()
  !----------------------------------------------------------------------------!
  ! This subroutine extracts the wfc computed directly from QE 
  ! and checks that the matrix 
  !             U_{mn} = < psi_m | psi_n >
  ! is the identity matrix for every k point in the QE folders.
  !
  !----------------------------------------------------------------------------!

  implicit none

  integer        :: iG
  integer        :: n_bnd, m_bnd
  integer        :: i_folder
  integer        :: G(3)

  integer        :: io_unit

  integer        :: list_iG(nG_max)
  integer        :: ipol

  complex(dp)    :: wfc_k(nG_max,num_bands,nspin)
!  complex(dp)    :: wfc_k(nG_max,nbands,nspin)

  real(dp)       :: eps

  real(dp)       :: QE_eig(num_bands)
!  real(dp)       :: QE_eig(nbands)

  real(dp)       :: norm

  complex(dp)    :: pw_mat_el(num_bands,num_bands,nspin,nspin)
!  complex(dp)    :: pw_mat_el(nbands,nbands,nspin,nspin)
  complex(dp)    :: U(num_bands,num_bands), identity(num_bands,num_bands)
!  complex(dp)    :: U(nbands,nbands), identity(nbands,nbands)

  logical        :: pass_u, global_pass_u

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio


  nbands_loc=num_bands

  ! initialize
  eps = eps_8
  global_pass_u = .true.
  G   = 0

  identity = cmplx_0
  do n_bnd = 1, nbands_loc
!  do n_bnd = 1, nbands
     identity(n_bnd,n_bnd) = cmplx_1
  end do 

  ! open a test file
  io_unit = find_free_unit()
  open(unit=io_unit,file=trim('orthogonality.test'))

  write(io_unit,10) '#=========================================================='
  write(io_unit,10) '# compute the inner product at all k points     ' 
  write(io_unit,10) '#=========================================================='
  write(io_unit,10) '#                                                          '
  write(io_unit,10) '# This test computes U_{mn} = < psi_m(QE) | psi_n(QE) >    '
  write(io_unit,10) '# where psi_m,n(QE) is read directly from file, regardless '   
  write(io_unit,10) '# of what the folders contain.                             '   
  write(io_unit,10) '#                                                          '
  write(io_unit,10) '# The test constitutes in checking that U is the identity  '
  write(io_unit,10) '# matrix.                                                  '
  write(io_unit,10) '#=========================================================='

  ! loop on all QE folders
  do i_folder = 1, nkpoints_QE

        ! get the wavefunctions 
        call get_K_folder_data(i_folder,list_iG,wfc_k,QE_eig)


        call get_plane_wave_matrix_element_convolution      &
                        (G,list_iG,list_iG, wfc_k,wfc_k,pw_mat_el)

	U(:,:) = cmplx_0	

	do ipol = 1,nspin
		U(:,:) = U(:,:) + pw_mat_el(:,:,ipol,ipol)
	end do

        write(io_unit,10) '----------------------------------------------------'
        write(io_unit,14) 'i_folder =', i_folder
        write(io_unit,10) ''

        ! Check U 
        pass_u     = .true.
        
        norm = sqrt(sum((real(U(:,:)-identity(:,:))**2+  &
                        aimag(U(:,:)-identity(:,:))**2)))


        if (norm > eps ) then
                pass_u        = .false.
                global_pass_u = .false.
                write(io_unit,10) '  * U matrix is NOT identity! : FAILED'

                write(io_unit,10) ' U ='
                do n_bnd = 1,nbands_loc
!                do n_bnd = 1,nbands
                     write(io_unit,30)  U(n_bnd,:)
                end do
                write(io_unit,10) ''
        end if
       
  end do 

  write(io_unit,10) '#=========================================================='
  write(io_unit,10) '# SUMMARY:'
  if ( global_pass_u ) then
        write(io_unit,10) '# * ALL U MATRICES IDENTITY      : PASSED' 
  else
        write(io_unit,10) '# * SOME U MATRICES NOT IDENTITY! : FAILED!' 
  end if

  write(io_unit,10) '#=========================================================='

  close(io_unit)

  10 format(A)
  11 format(A,E10.4,A)
  12 format(A,I4,A)
  13 format(A,E10.4)
  14 format(A,I6)
  15 format(10x3I10)
  20 format(A,I4,A,I4,A,I4)
  30 format(100(2F8.4,5X))
  35 format(100F8.4)

  end subroutine test_wfc_orthogonality

  subroutine update_test_density(test_density,nk_vec_max,nk_vec,       &
                         eig_ks,  pw_matrix_elements)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to test_density.
  !
  !------------------------------------------------------------------------
  implicit none

  ! input variables
  integer      :: nk_vec, nk_vec_max
  real(dp)     :: eig_ks (num_wann,nk_vec_max)


  complex(dp)  :: test_density(nG_shell_max,nG_shell_max,nspin,nspin)

  complex(dp)  :: pw_matrix_elements(num_wann,num_wann,nspin,nspin,   &
                                           nk_vec_max,nG_shell_max)

  complex(dp)  :: m1(nspin,nspin), m2(nspin,nspin)

  ! computation variables
  integer      :: iG_1, iG_2

  integer      :: ipol1, ipol2, ipol3, ipol4

  integer      :: k_vec_loop
  integer      :: nbnd1, nbnd2
  real(dp)     :: ek
  integer      :: I_band



  do k_vec_loop = 1,nk_vec

     do I_band = 1,num_wann*num_wann

          nbnd1  = modulo(I_band-1,num_wann)+1
          nbnd2  = (I_band-nbnd1)/num_wann+1

          ek    =  eig_ks (nbnd1,k_vec_loop)

	  if (ek  > chemical_potential ) cycle

          do iG_1 = 1, nG_shell_max
             do iG_2 = 1, nG_shell_max

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
               !---------------------------------------------------- 

               m1(:,:)=pw_matrix_elements(nbnd1,nbnd2,:,:,k_vec_loop,iG_1)

	       ! take transpose to account for comment above
               m2(:,:)=transpose(conjg(           &
			pw_matrix_elements(nbnd1,nbnd2,:,:,k_vec_loop,iG_2)))
             


               do ipol1 = 1,nspin
                do ipol4 = 1,nspin
                
                   do ipol2 = 1,nspin
                    do ipol3 = 1,nspin


                     test_density(iG_1,iG_2,ipol1,ipol4) =              &
                     		test_density(iG_1,iG_2,ipol1,ipol4) +   &
				  m1(ipol1,ipol2)*m2(ipol3,ipol4)

                    enddo
                   enddo

                enddo
               enddo


             end do 
          end do 

     end do  ! I_band

  end do  ! kvec_loop



11 format(A)
12 format(A,ES12.3,A)


  end subroutine update_test_density

  subroutine output_test_density(test_density)
  !------------------------------------------------------------------------
  ! This subroutine adds a sub-block contribution to chi0
  !------------------------------------------------------------------------
  implicit none


  complex(dp)  :: test_density(nG_shell_max,nG_shell_max,nspin,nspin)


  real(dp)     :: volume,  rho_0, prefactor, spin_sum_factor
  integer      :: io_unit
  integer      :: iG_1,  iG_2, iG

  integer      :: i1,i2,i3
  integer      :: ipol1, ipol2
  integer      :: G(3)


  character(256) :: filename

  ! normalize at the very end, to avoid truncation error

  ! In the paramagnetic case, all spin states are degenerate
  ! and there is a factor of 2 coming from a sum on spin. 
  ! This factor should not be present in the case of a non-paramagnetic 
  ! system.

  ! the sum on ipol2, ipol3 in update_test_density gives the
  ! trace of the nspin x nspin identity matrix,
  ! which is equal to 2 for nspin = 2. This is why
  ! we must divide by 2 if nspin =2; on the other hand,
  ! if nspin ==1, we are looking for the charge density, so there
  ! is a sum on spin to be performed, hence a factor of 2.


  if (nspin == 1) then
	spin_sum_factor = 2.0_dp
  else 
	spin_sum_factor = 0.5_dp
  end if


  volume = alat**3*abs(                                &
           at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))

  rho_0     = 1.0_dp/volume


  test_density(:,:,:,:) = test_density(:,:,:,:)/dble(nk1s*nk2s*nk3s) &
                          *rho_0*spin_sum_factor 


  io_unit  = find_free_unit()

  filename = trim('test_density.test')

  open(unit = io_unit, file = filename)

  write(io_unit,5) '#=========================================================='
  write(io_unit,5) '#   This file contains the computed density   '
  write(io_unit,5) '#   using the interpolated wavefunctions and  '
  write(io_unit,5) '#   eigenvalues.                              '
  write(io_unit,5) '#                                             '
  write(io_unit,5) '#   The q vector, in crystal coordinates:     '
  write(io_unit,5) '#                                             '
  write(io_unit,10)'#         q(1)  = ', qpt1
  write(io_unit,10)'#         q(2)  = ', qpt2
  write(io_unit,10)'#         q(3)  = ', qpt3
  write(io_unit,5) '#                                             '
  write(io_unit,5) '#                                             '

  if ( nspin == 1) then
  	write(io_unit,5) '#  1st   column: iG_1                         '
	write(io_unit,5) '#  2nd   column: iG_2                         '
	write(io_unit,5) '#  3rd   column: iG (index of G = G1-G2)      '
	write(io_unit,5) '#  4rd   column: Re[n_{G}] (1/a0^3)           '
	write(io_unit,5) '#  5th   column: Im[n_{G}] (1/a0^3)           '
	write(io_unit,5) '#                                             '                     
	write(io_unit,5) '#=========================================================='
	write(io_unit,5) '#  iG_1  iG_2  iG     Re[nG]          Im[nG]'
	write(io_unit,5) '#=========================================================='
  else
  	write(io_unit,5) '#  1st   column: iG_1                         '
  	write(io_unit,5) '#  2nd   column: iG_2                         '
	write(io_unit,5) '#  3rd   column: iG (index of G = G1-G2)      '
  	write(io_unit,5) '#  4rd   column: spin combination             '
  	write(io_unit,5) '#  5rd   column: Re[n_{G}] (1/a0^3)           '
  	write(io_unit,5) '#  6th   column: Im[n_{G}] (1/a0^3)           '
  	write(io_unit,5) '#                                             '                     
  	write(io_unit,5) '#=========================================================='
  	write(io_unit,5) '#  iG_1  iG_2  iG   spin     Re[nG]          Im[nG]'
  	write(io_unit,5) '#=========================================================='
  end if

  do iG_1 = 1,nG_shell_max
     do iG_2 = 1,nG_shell_max

	G  = gvec(:,iG_1)-gvec(:,iG_2)

        i1 = modulo(G(1), nr1)+1
        i2 = modulo(G(2), nr2)+1
        i3 = modulo(G(3), nr3)+1

	iG = g_fft_map(i1,i2,i3) 


        if (nspin == 2) then
          ! write to file
	  do ipol1 = 1, nspin
	     do ipol2 = 1, nspin
        	  write(io_unit,20) iG_1,iG_2,iG, ipol1, ipol2,          &
			  real(test_density(iG_1,iG_2,ipol1,ipol2)),     &
                         aimag(test_density(iG_1,iG_2,ipol1,ipol2))
             end do
          end do

          write(io_unit,5) ''

        else
        	  write(io_unit,22) iG_1,iG_2,iG,                  &
			  real(test_density(iG_1,iG_2,1,1)),       &
                         aimag(test_density(iG_1,iG_2,1,1))


        end if
     end do !iG_2
  end do !iG_1

  close(io_unit)


  5   format(A)
  10  format(A,3G18.8E3)
  15  format(A,I6)
  20  format(1X,I5,1X,I5,1X,I5,2X,2I1,2F16.8)
  22  format(1X,I5,1X,I5,1X,I5,4X,2F16.8)
 

  end subroutine output_test_density

!----------------------------------------------------------------------------!
!
!
end module intw_tests
!
!
!----------------------------------------------------------------------------!
