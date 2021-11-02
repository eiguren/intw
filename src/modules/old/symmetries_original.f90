!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
!

module intw_symmetries
  !
  !----------------------------------------------------------------------------!
  !
  !       The subroutines in this module handle the space group symmetries 
  !       and their effect on wave-functions and matrix elements. They 
  !       also allow the identification of an irreducible set of k-vectors in
  !       the 1BZ, forming the IBZ.
  !
  !       It is advantageous to perform symmetry operations in crystal coordinates;
  !       indeed, in these coordinates, the zone is simply a cube and MP meshes
  !       are nice and regular. However, the crystal coordinate system
  !       introduces minor complications in the application of point group 
  !       symmetry. A full description is given in algorithms.pdf.
  !
  !----------------------------------------------------------------------------!

  use intw_useful_constants 
  use intw_utility      ! provides useful functions!
  use intw_reading      
  use intw_input_parameters, only : TR_symmetry !Only this, not equiv and symlink!.

  use intw_fft

!Peio
!Number of Bloch Original states for the wannierization
  use w90_parameters, only: num_bands
!Peio

  !
  implicit none
  !
  save
  !
  integer,allocatable :: symlink(:,:)
  integer,allocatable :: sym_G(:,:)
  integer,allocatable :: nosym_G(:,:)

  integer,allocatable :: QE_folder_sym(:), QE_folder_nosym(:) 

  integer,allocatable :: inverse_indices(:)
  integer             :: identity_matrix_index

  ! this variable will contain the 2x2 matrices which rotate
  ! spinors.
  complex(dp),allocatable :: spin_symmetry_matrices(:,:,:)

  ! logical variables defining what is present in the QE folders
  logical             :: full_mesh, IBZ

  !Asier && Idoia 17 07 2014 
  !Identy of atoms under symmetry operation 
  integer, allocatable  :: rtau_index(:,:) 
  real(dp), allocatable :: tau_cryst(:,:)
  real(dp), allocatable :: rtau(:,:,:)
  real(dp), allocatable :: rtau_cryst(:,:,:)

  integer               :: symtable(48,48)

contains 

  subroutine allocate_symmetry_related_k(nk1,nk2,nk3,nsym)
    !------------------------------------------------------------------
    !
    !     nk1, nk2, nk3 are the MP coefficient of the mesh, and nsym is
    !     the total number of point group operations.
    !
    !     Define kpt to be a vector in the 1BZ, and ikpt to be its
    !     singlet index. Let kpt_irr be a symmetry equivalent k-point to kpt.
    !     Finally, let i_folder be an index which represents the QE folders.
    !
    !     This subroutine allocates the arrays
    ! 
    !     - QE_folder_sym   : this array indicates in which QE folder "i_folder"
    !                         the kpoint "kpt_irr", which is symmetry equivalent to
    !                         "kpt", can be found.
    !
    !     - symlink         : What rotation-like operation must be performed
    !                         on "kpt_irr" to obtain "kpt".
    !
    !     - sym_G           : What G translation must be applied to the rotated 
    !                         "kpt_irr" to obtain "kpt".
    !
    !     In equations:
    !             R        = symlink(ikpt)
    !             i_folder = QE_folder_sym(ikpt)
    !             G        = sym_G(ikpt)
    !
    !             =========>  R*kpt_irr = kpt+G
    !
    !     Also, in the case that a full mesh is present in the QE folders
    !     the following arrays will be allocated
    !     - QE_folder_nosym : this array indicates in which QE folder "i_folder"
    !                         the kpoint "kpt_1", which is translation equivalent to
    !                         "kpt", can be found.
    !
    !     - nosym_G         : What G translation must be applied to 
    !                         "kpt_1" to obtain "kpt".
    !
    !     In equations:
    !             i_folder = QE_folder_nosym(ikpt)
    !             G        = nosym_G(ikpt)
    !
    !             =========>  kpt_1 = kpt+G
    !------------------------------------------------------------------

    integer  :: nk1, nk2, nk3 , nsym
    integer  :: nkmesh

    nkmesh = nk1*nk2*nk3

    allocate(inverse_indices(nsym))
    allocate(QE_folder_sym(nkmesh))
    allocate(sym_G(3,nkmesh) )
    allocate(symlink(nkmesh,2) )

    allocate(QE_folder_nosym(nkmesh))

    allocate(nosym_G(3,nkmesh) )

  end subroutine allocate_symmetry_related_k

  subroutine deallocate_spin_symmetry_matrices()
    !------------------------------------------------------------------
    ! This subroutine deallocates the array spin_symmetry_matrices
    !------------------------------------------------------------------
    deallocate(spin_symmetry_matrices)

  end subroutine deallocate_spin_symmetry_matrices

  subroutine deallocate_symmetry_related_k()
    !------------------------------------------------------------------
    ! This subroutine deallocates the arrays equiv and symlink.
    !------------------------------------------------------------------

    if (allocated(symlink)) then
       deallocate(symlink)
    end if

    if (allocated(sym_G)) then
       deallocate(sym_G)
    end if

    if (allocated(QE_folder_sym)) then
       deallocate(QE_folder_sym)
    end if
    if (allocated(QE_folder_nosym)) then
       deallocate(QE_folder_nosym)
    end if
    if (allocated(inverse_indices)) then
       deallocate(inverse_indices)
    end if
    if (allocated(nosym_G)) then
       deallocate(nosym_G)
    end if



  end subroutine deallocate_symmetry_related_k

  subroutine find_inverse_symmetry_matrices_indices()
    !------------------------------------------------------------------
    ! For every point group operation isym, this subroutine finds 
    ! the index of the inverse point group operation. 
    !  
    ! This will be useful because, by convention, in crystal coordinates,
    !
    !     R_alpha * k ===>   \sum_j s(R^{-1}_alpha)_{ij} k_j
    !
    ! (This extra layer of complication is QE's fault, not mine) 
    !------------------------------------------------------------------
    use intw_useful_constants 
    use intw_reading,               only : s,  nsym, at, bg
    use intw_utility 

    implicit none

    real(dp)       ::  rot_cart(3,3,nsym)

    real(dp)       ::  rot1(3,3), rot2(3,3), prod(3,3)

    real(dp)       ::  norm 


    integer        ::  isym,  jsym, ns

    integer        ::  i, j,  l,   mu,  nu
    logical        ::  found

    ! Find the rotation matrices in cartesian coordinates
    rot_cart  =  zero

    do isym =1, nsym

       do mu=1,3
          do nu=1,3
             do i=1,3
                do j=1,3

                   rot_cart(mu,nu,isym)  =    rot_cart(mu,nu,isym)     &
                        +  bg(nu,i)*s(i,j,isym)*at(mu,j)
                end do
             end do
          end do
       end do

    end do ! isym

    ! Find the inverse of every symmetry matrix 
    inverse_indices = 0

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

          norm      = zero 

          do i=1,3
             do j=1,3
                do l = 1,3 
                   prod(i,j)  =  prod(i,j) + rot1(i,l)*rot2(l,j)
                end do
                norm = norm + prod(i,j)**2
             end do
          end do

          norm = sqrt(norm)

          if (norm < eps_8) then 
             inverse_indices(isym) = jsym
             found = .true.
             exit
          end if

       end do  ! jsym

       if (.not. found) then
          write(*,'(a,I3)')         &
               'ERROR: no inverse matrix found for operation isym = ',isym
          stop
       end if

    end do  ! isym

    do ns = 1, nsym
       if (s(1,1,ns) == 1 .and. s(1,2,ns) == 0 .and. s(1,3,ns) == 0 .and. &
            s(2,1,ns) == 0 .and. s(2,2,ns) == 1 .and. s(2,3,ns) == 0 .and. &
            s(3,1,ns) == 0 .and. s(3,2,ns) == 0 .and. s(3,3,ns) == 1) then

          identity_matrix_index = ns
          exit
       end if
    end do

  end subroutine find_inverse_symmetry_matrices_indices

  subroutine allocate_and_build_spin_symmetry_matrices(nsym)
    !------------------------------------------------------------------
    !
    ! This subroutine builds, once and for all, all the 2x2 spin 
    ! rotation matrices needed by the program and tabulates them in the
    ! array spin_symmetry_matrices. 
    !
    ! Some elementary testing is also implemented.
    !
    !------------------------------------------------------------------
    use  intw_useful_constants

    implicit none

    !input
    integer  :: nsym

    ! local variables
    integer  :: isym

    ! symmetry matrix, in crystal coordinates
    integer        :: sym(3,3)

    ! the 2x2 spin rotation matrix
    complex(dp)    :: S_u(2,2) 

    ! axis and angle of a given rotation matrix
    real(dp)       :: axis(3), angle

       !-Pauli matrices
        complex(kind=dp) :: I(2,2), sig_x(2,2), sig_y(2,2), sig_z(2,2)

    I2   (1,1)= cmplx_1 ;    I2(1,2)= cmplx_0 ;    I2(2,1)= cmplx_0  ;    I2(2,2)= cmplx_1
    sig_x(1,1)= cmplx_0 ; sig_x(1,2)= cmplx_1 ; sig_x(2,1)= cmplx_1  ; sig_x(2,2)= cmplx_0
    sig_y(1,1)= cmplx_0 ; sig_y(1,2)=-cmplx_i ; sig_y(2,1)= cmplx_i  ; sig_y(2,2)= cmplx_0
    sig_z(1,1)= cmplx_1 ; sig_z(1,2)= cmplx_0 ; sig_z(2,1)= cmplx_0  ; sig_z(2,2)=-cmplx_1


    allocate(spin_symmetry_matrices(2,2,nsym))

    ! build the spin symmetry matrices for all symmetry operations
    do isym = 1, nsym


       sym  = s(:,:,isym)

       call rotaxis_crystal(sym,axis,angle)

       ! debugging test
       ! angle = -angle

       S_u(:,:) = dcos(angle/2.d0) * I2(:,:)           - &
            cmplx_i * dsin(angle/2.d0)* (        &
            axis(1)*sig_x(:,:) + &
            axis(2)*sig_y(:,:) + &
            axis(3)*sig_z(:,:)   )


       spin_symmetry_matrices(:,:,isym)  = S_u(:,:) 

    end do !isym


  end subroutine allocate_and_build_spin_symmetry_matrices


  subroutine rotaxis_crystal(sym,axis,angle)
    !------------------------------------------------------------------
    !
    ! This subroutine finds the axis and the angle of rotation 
    ! for a given point group operation sym. Note that "sym" is in 
    ! crystal coordinates, and that it might be the composition of
    ! a rotation and the inversion.
    !
    ! Define a point group operation in cartesian coorinates S^{cart}
    ! such that
    ! 
    !		r'_{alpha} = \sum_{beta=1}^3 S^{cart}_{alpha beta} r_{beta} 
    !
    ! In crystal coordinates, and following QE's conventions,  this
    ! leads to
    !            x'_{i}  = \sum_{j=1}^3  S^{cryst}_{ji} x_{j}
    !
    !	where             
    !			S^{cryst}_{ji} = 1/2pi b_i * S^{cart} * a_j
    !	which becomes
    !			S^{cryst} = A * [S^{cart}^T] * B
    !
    !			with A = [ - a_1 - ]     B = [  |   |   |  ]
    !			         [ - a_2 - ]     B = [ b_1 b_2 b_3 ]
    !			         [ - a_3 - ]     B = [  |   |   |  ]
    !
    !	S^{cryst} is the matrix read from Quantum Espresso.
    !
    !	-	The angle is given by 1+2cos(angle) = Tr[S], which 
    !		is independent of the basis.
    !
    !	-	The angle is given by 1+2cos(angle) = Tr[S], which 
    !		is independent of the basis.
    !
    !------------------------------------------------------------------

    implicit none

    !input
    integer        :: sym(3,3)          ! symmetry matrix, in crystal coordinates

    ! output 
    real(dp)       :: axis(3), angle    ! axis and angle of a given rotation matrix

    ! local variables
    integer        :: determinant, trace

    real(dp)       :: two_sin_angle

    integer        :: Rot_cryst(3,3)   ! rotation part of the operation, in crystal
    ! coordinates

    real(dp)       :: Rot_cart(3,3)      ! rotation matrix, in cartesian coordinates

    integer        :: mu, nu, i, j

    real(dp)       :: x_tmp
    real(dp)       :: norm_axis(3), sign_axis(3)
    logical        :: vanish(3)

    ! First, find the determinant
    determinant = sym(1,1) * sym(2,2) * sym(3,3) +  &
         sym(1,2) * sym(2,3) * sym(3,1) +  &
         sym(1,3) * sym(2,1) * sym(3,2) -  &
         sym(1,3) * sym(2,2) * sym(3,1) -  &
         sym(1,1) * sym(2,3) * sym(3,2) -  &
         sym(1,2) * sym(2,1) * sym(3,3)

    ! Put the rotation part of the symmetry matrix in the Rotation array, 
    ! multiplied by the appropriate coefficient to account for inversion

    if ( determinant == 1) then
       Rot_cryst(:,:) = sym(:,:)
    else if ( determinant == -1) then
       Rot_cryst(:,:) = -sym(:,:)
    else
       write(*,*) '************************************************'
       write(*,*) '** ERROR: The determinant of the rotation     **'
       write(*,*) '**        matrix is not +/- 1.                **'
       write(*,*) '**        review code.                        **'
       write(*,*) '**          program stops.                    **'
       write(*,*) '************************************************'
       stop
    end if

    ! compute the rotation matrix in cartesian coordinates
    Rot_cart(:,:)  =   ZERO

    do mu=1,3
       do nu=1,3
          do i=1,3
             do j=1,3

                Rot_cart(mu,nu)  =    Rot_cart(mu,nu)     &
                     +  bg(nu,i)*Rot_cryst(i,j)*at(mu,j)
             end do
          end do
       end do
    end do

    ! Extract the rotation angle from the trace of the matrix
    trace = Rot_cryst(1,1) +  Rot_cryst(2,2) +  Rot_cryst(3,3) 

    ! there are only 5 possibilities in a crystal;
    ! tabulating insures there is no problem with picking the right quadrant.

    if     (trace == -1 ) then
       angle         =  pi
       two_sin_angle =  ZERO

    else if (trace == 0 ) then
       angle         =  2.0_dp*pi/3.0_dp
       two_sin_angle =  sqrt(3.0_dp)

    else if (trace == 1 ) then
       angle         =  pi/2.0_dp
       two_sin_angle =  2.0_dp

    else if (trace == 2 ) then
       angle         =  pi/3.0_dp
       two_sin_angle =  sqrt(3.0_dp)

    else if (trace == 3 ) then
       angle         =  ZERO
       two_sin_angle =  ZERO
    else
       write(*,*) '************************************************'
       write(*,*) '** ERROR: The trace of the rotation matrix    **'
       write(*,*) '**        is not in [-1,0,1,2,3]. review code.**'
       write(*,*) '**              program stops.                **'
       write(*,*) '************************************************'
       stop
    end if

    ! build the axis array
    if     (trace == -1 ) then
       ! This is the complicated case. Since angle = pi, 
       ! the cartesian rotation matrix is symmetric. 
       ! A bit of cleverness is required to extract the axis vector.

       ! First, find the norms of the coordinates 
       x_tmp = Rot_cart(1,1)+ONE 
       if (x_tmp > eps_8) then
          norm_axis(1) = sqrt(x_tmp)/sqrt(2.0_dp)
          vanish   (1) = .false.
       else
          norm_axis(1) = ZERO
          vanish(1)    = .true.
       end if

       x_tmp = Rot_cart(2,2)+ONE 
       if (x_tmp > eps_8) then
          norm_axis(2) = sqrt(x_tmp)/sqrt(2.0_dp)
          vanish   (2) = .false.
       else
          norm_axis(2) = ZERO
          vanish   (2) = .true.
       end if

       x_tmp = Rot_cart(3,3)+ONE 
       if (x_tmp > eps_8) then
          norm_axis(3) = sqrt(x_tmp)/sqrt(2.0_dp)
          vanish   (3) = .false.
       else
          norm_axis(3) = ZERO
          vanish   (3) = .true.
       end if


       if (.not. vanish(1) .and. .not. vanish(2) .and. .not. vanish(3) ) then
          ! if no component vanishes, arbitrarily set the sign of
          ! n3 to be positive.
          sign_axis(3) = ONE

          sign_axis(2) = Rot_cart(2,3)/(2.0_dp*norm_axis(2)*norm_axis(3))
          sign_axis(1) = Rot_cart(1,3)/(2.0_dp*norm_axis(1)*norm_axis(3))

       else if (.not. vanish(1) .and. .not. vanish(2) .and. vanish(3) ) then
          ! if one component vanishes, arbitrarily set the sign of the largest index
          ! component to be positive.
          sign_axis(3) = ZERO
          sign_axis(2) = ONE
          sign_axis(1) = Rot_cart(1,2)/(2.0_dp*norm_axis(1)*norm_axis(2))

       else if (.not. vanish(1) .and. vanish(2) .and. .not. vanish(3) ) then
          sign_axis(3) = ONE
          sign_axis(2) = ZERO
          sign_axis(1) = Rot_cart(1,3)/(2.0_dp*norm_axis(1)*norm_axis(3))

       else if (vanish(1) .and. .not. vanish(2) .and. .not. vanish(3) ) then
          sign_axis(3) = ONE
          sign_axis(1) = ZERO
          sign_axis(2) = Rot_cart(2,3)/(2.0_dp*norm_axis(2)*norm_axis(3))

       else if (vanish(1) .and. vanish(2) .and. .not. vanish(3) ) then
          ! if two components vanish, arbitrarily set the sign of the non
          ! vanishing component to be positive.
          sign_axis(1) = ZERO
          sign_axis(2) = ZERO
          sign_axis(3) = ONE

       else if (vanish(1) .and. .not. vanish(2) .and. vanish(3) ) then
          sign_axis(1) = ZERO
          sign_axis(2) = ONE
          sign_axis(3) = ZERO
       else if (.not. vanish(1) .and. vanish(2) .and. vanish(3) ) then
          sign_axis(1) = ONE 
          sign_axis(2) = ZERO
          sign_axis(3) = ZERO
       end if

       axis(1) = norm_axis(1)*sign_axis(1)
       axis(2) = norm_axis(2)*sign_axis(2)
       axis(3) = norm_axis(3)*sign_axis(3)

    else if (trace == 0 .or. trace ==1 .or. trace ==2 ) then
       ! For these cases, sin(alpha) is not zero
       ! we can extract the axis from the off diagonal elements
       ! of the rotation matrix.
       axis(1) = ( Rot_cart(3,2)-Rot_cart(2,3) )/two_sin_angle 
       axis(2) = ( Rot_cart(1,3)-Rot_cart(3,1) )/two_sin_angle 
       axis(3) = ( Rot_cart(2,1)-Rot_cart(1,2) )/two_sin_angle 

    else if (trace == 3 ) then
       ! This is the simplest case: it corresponds to
       ! the identity operation, and the direction of the axis
       ! is irrelevant.
       axis	  =  ZERO

    end if


  end subroutine rotaxis_crystal


  subroutine set_symmetry_relations(nk1, nk2, nk3, nkpoints_QE, kpoints_QE,kmesh,       &
       k_points_consistent, QE_folder_sym, sym_G, symlink)
    !--------------------------------------------------------------------------!
    !     This subroutine fills in the symmetry arrays allocated by
    !     the subroutine "allocate_symmetry_related_k", and tests if the 
    !     contents of the QE_folders are consistent with the MP mesh.
    !
    !     If the QE folders contain a full MP mesh, this subroutine will
    !     also tabulate the relationship between the k-points in the folders
    !     and the canonical k-points. This will be useful for testing.
    !--------------------------------------------------------------------------
    implicit none

    !input
    integer,intent(in)  :: nk1, nk2, nk3, nkpoints_QE                ! Number of k points found in the QE folders
    real(dp),intent(in)  :: kpoints_QE (3,nkpoints_QE) ! The kpoints read from the QE folders
    real(dp),intent(in)  :: kmesh      (3,nk1*nk2*nk3) ! The kmesh points, in canonical order 

    logical,intent(out)   :: k_points_consistent        ! test if the kpoints are consistent
    ! with the parameters of the input file

    !output
    !     This subroutine will fill in the global arrays
    !             - QE_folder_sym(nk1*nk2*nk3)
    !             - sym_G      (3,nk1*nk2*nk3)
    !             - symlink      (nk1*nk2*nk3,2)
    !
    !     IF a full mesh is potentially present, it will also fill the arrays
    !             - QE_folder_nosym(nk1*nk2*nk3)
    !             - nosym_G      (3,nk1*nk2*nk3)


    !Asier && Idoia 21 07 2014  
    integer, intent(out) :: QE_folder_sym(nk1*nk2*nk3)
    integer, intent(out) :: sym_G      (3,nk1*nk2*nk3)
    integer, intent(out) :: symlink      (nk1*nk2*nk3,2)



    !local variables

    real(dp):: kpt(3) , rotated_kpt(3), kpt_in_1BZ(3)
    real(dp):: kpt_on_mesh(3), d_kpt(3) 

    logical :: kpoint_is_found_sym(nk1*nk2*nk3), kpoint_is_found_nosym(nk1*nk2*nk3) 

    logical :: possible_full_mesh

    integer :: nkmesh
    integer :: ikpt, switch, i,  j,  k, i_folder

    integer :: isym

    integer :: G(3)


    ! initialize arrays to a negative number: if there is a bug
    ! in the code, it will thus look for inexistent folders. It
    ! is better to crash than to produce wrong results!
    QE_folder_sym(:)    = -4
    QE_folder_nosym(:)  = -4


    k_points_consistent      = .true.

    kpoint_is_found_sym(:)   = .false.
    kpoint_is_found_nosym(:) = .false.
    nkmesh = nk1*nk2*nk3

    if ( nkpoints_QE  == nkmesh ) then
       possible_full_mesh  = .true.
    else
       possible_full_mesh  = .false.
    end if

    switch = 1 ! triplet to singlet

    do i_folder = 1, nkpoints_QE

       ! coordinates of this k-point, which need not lie in the 1BZ
       kpt = kpoints_QE(:,i_folder)

       ! extract the triplet coordinates, the kpt in the 1BZ and
       ! the G vector
       call find_k_1BZ_and_G(kpt,nk1,nk2,nk3,i,j,k,kpt_in_1BZ,G)

       ! test that this triplet index indeed produces the k-point.
       ! This tests that the kpoint is indeed on a mesh consistent
       ! with the input file.
       kpt_on_mesh(1) = dble(i-1)/dble(nk1) 
       kpt_on_mesh(2) = dble(j-1)/dble(nk2) 
       kpt_on_mesh(3) = dble(k-1)/dble(nk3) 

       ! make sure this is identical to the original k-point
       d_kpt(:) = kpt_in_1BZ(:) - kpt_on_mesh(:)

       if (sqrt(dot_product(d_kpt,d_kpt)) > eps_8) then
          k_points_consistent = .false.     
          write(*,*) ' consistency FAILURE '
          write(*,'(A,3F8.4)') ' kpt        = ',kpt
          write(*,'(A,3F8.4)') ' kpt_in_1BZ = ',kpt_in_1BZ
          write(*,'(A,3F8.4)') ' kpt_on_mesh= ',kpt_on_mesh
          write(*,'(A,3I4)')   ' i, j , k   = ',i,j,k
          exit
       end if

       ! if the k-points are not consistent with the input parameters,
       ! an error will be thrown in the main code. Assume consistency
       ! from this point on.

       ! if there is the possibility of a full mesh being present, 
       ! find the correspondence
       if (possible_full_mesh) then

          call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,switch)

          kpoint_is_found_nosym(ikpt) = .true.
          QE_folder_nosym(ikpt)       = i_folder
          nosym_G(:,ikpt)             = G

       end if !possible_full_mesh

       ! loop on all point group operations 
       do isym = 1,nsym

          ! rotate the k-point.
          rotated_kpt =  matmul(dble(s(:,:,isym)),kpt)
          !  There is an added layer of complexity introduced by the fact
          !  that we are using crystal coordinates. The convention for the
          !  action of the crystal coordinate point group matrices is:
          !
          !  R_mu * k^{cart}  => sum_{j} s(R^{-1}_mu)_{ij} k^{cryst}_j
          !
          !  Thus, it is actually the index of the INVERSE which must be used.
          !  This may seem like a trivial change, but it affects the phase
          !  factor which must be introduced when a wavefunction is rotated,
          !  in the case of a non-symmorphic group.
          !  Find the corresponding k-point in the canonical 1BZ.

          ! extract the triplet coordinates, the kpt in the 1BZ and
          ! the G vector
          call find_k_1BZ_and_G(rotated_kpt,nk1,nk2,nk3,i,j,k,kpt_in_1BZ,G)


          ! Tabulate this point as found, but only if allowed to do so!

!          if (nspin == 1 .or. (nspin ==2 .and. (.not. can_use_TR(isym)))) then
             ! find its singlet coordinate
             call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,+1)

             ! if this point hasn't been found before, well, it's found now!
             if ( .not. kpoint_is_found_sym(ikpt) ) then
              	kpoint_is_found_sym(ikpt) = .true.
              	QE_folder_sym(ikpt)       = i_folder
              	sym_G(:,ikpt)             = G
              	symlink(ikpt,1)           = inverse_indices(isym)
              	symlink(ikpt,2)           = 0

             end if
!          end if

        ! --Asier&&Idoia 13 01 2015: Try first everything except TR 
        end do !

        do isym = 1,nsym
        ! -- end Asier&&Idoia 13 01 2015

          rotated_kpt =  matmul(dble(s(:,:,isym)),kpt)

          ! repeat with TR symmetry, if allowed (or required!) 
          if (can_use_TR(isym)) then
             rotated_kpt =  -rotated_kpt

             call find_k_1BZ_and_G(rotated_kpt,nk1,nk2,nk3,i,j,k,kpt_in_1BZ,G)

             ! find its singlet coordinate
             call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,switch)

             if ( .not. kpoint_is_found_sym(ikpt) ) then

                kpoint_is_found_sym(ikpt) = .true.
                QE_folder_sym(ikpt)       = i_folder
                sym_G(:,ikpt)             = -G
                symlink(ikpt,1)           = inverse_indices(isym)
                symlink(ikpt,2)           = 1

                ! CAREFUL! It is now -G which enters the sym_G array
                ! See the rotation code for details, as well as the
                ! symmetries document.

             end if ! not found
          end if ! TR

       enddo !isym
    enddo !i_folder


    ! test if all the kmesh points were found, and store the result
    ! in these GLOBAL logical variables


    IBZ         = all( kpoint_is_found_sym )
    full_mesh   = all( kpoint_is_found_nosym )

  end subroutine set_symmetry_relations
  
  subroutine find_the_irreducible_k_set_and_equiv2 (nspt, k_entire, k_irr, nk_irr, equiv, symlink, sym_G)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the a general k list 
    !------------------------------------------------------------------

    use intw_reading, only : s,  nsym
    use intw_input_parameters, only : TR_symmetry !Only this, not equiv and symlink!.

    use intw_utility ! provides useful functions!

    implicit none

    !input

    integer  :: nspt

    real(dp) :: k_entire (3,nspt), k_aux(1:3,48)

    ! output
    integer  :: nk_irr
    real(dp) :: k_irr  (3, nspt) 
    integer :: equiv(nspt) ! which is the equivalent point
    integer :: sym_G(1:3, nspt), symlink(nspt,1:2)

    !local variables
    real(dp):: k_rot(3), k_1BZ(3), dk(3), dist1, dist2

    integer :: nkpt             ! The total number of points 
    integer :: i,  j,   k       ! triplet indices 
    integer :: is, js,  ks      ! triplet indices  obtained by symmetry

    integer :: G(3)

    integer :: ns

    integer :: ikpt, ikpts      ! singlet index, singlet index obtained by symmetry

    integer :: switch, ii, jj, kk, nsp

    logical :: found(nspt) 
    real(dp), parameter :: eps=10E-7

    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    switch = 1  ! triplet-to-singlet coordinate

    found  = .false.
    nk_irr = 0        

    do i=1, nspt
       if (.not. found(i)) then
          found(i) = .true.
          nk_irr = nk_irr + 1  
          k_irr(1:3,nk_irr) = k_entire(1:3,i)
          equiv(i) = nk_irr 

          do ns=1,nsym
             k_rot = matmul(dble(s(:,:,ns)), k_irr(:,nk_irr))

             do j=1, nspt

                if (((abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)-1.0_dp)<eps).or. &
                     (abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)       )<eps)).and.&
                     ((abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)-1.0_dp)<eps).or. &
                     (abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)       )<eps)).and.&
                     ((abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)-1.0_dp)<eps).or. &
                     (abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)       )<eps)) ) then 

                   equiv(j)   = nk_irr 
                   found(j)   =.true.
                   symlink(j,1) = ns
                   symlink(j,2) = 0 
                   sym_G(:,j) = nint(k_entire(:,j) - k_rot(:))
                   cycle 
                endif

                if (TR_symmetry) then

                   if (((abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)-1.0_dp)<eps).or. &
                        (abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)       )<eps)).and.&
                        ((abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)-1.0_dp)<eps).or. &
                        (abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)       )<eps)).and.&
                        ((abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)-1.0_dp)<eps).or. &
                        (abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)       )<eps)) ) then
                      equiv(j) = nk_irr 
                      found(j) =.true.
                      symlink(j,1) = ns
                      symlink(j,2) = 1 
                      sym_G(:,j) = nint(k_entire(:,j) -(- k_rot(:)))

                   endif
                endif

             end do !j 

          enddo ! ns 

       end if !found

    enddo !i 


    do ikpt=1, nk_irr

       k_rot(:) = matmul(bg,k_irr(:,ikpt))

       dist1 = sum (k_rot(:)**2)

       do ii=-1,1
          do jj=-1,1
             do kk=-1,1

                k_rot(:) = matmul(bg,k_irr(:,ikpt))  + matmul(bg, dble((/ii,jj,kk/)))
                dist2 = sum (k_rot(:)**2) 

                if (dist2<dist1) then 
                   k_irr(:,ikpt) = matmul (transpose(at), k_rot)
                   dist1 = dist2
                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine find_the_irreducible_k_set_and_equiv2

  subroutine find_the_irreducible_k_set_and_equiv (nspt, k_entire, k_irr, nk_irr, equiv)!, symlink, sym_G)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the a general k list 
    !------------------------------------------------------------------

    use intw_reading, only : s,  nsym
    use intw_input_parameters, only : TR_symmetry !Only this, not equiv and symlink!.

    use intw_utility ! provides useful functions!

    implicit none

    !input

    integer  :: nspt

    real(dp) :: k_entire (3,nspt), k_aux(1:3,48)

    ! output
    integer  :: nk_irr
    real(dp) :: k_irr  (3, nspt) 
    integer :: equiv(nspt) ! which is the equivalent point
    integer :: sym_G(1:3, nspt), symlink(nspt,1:2)

    !local variables
    real(dp):: k_rot(3), k_1BZ(3), dk(3), dist1, dist2

    integer :: nkpt             ! The total number of points 
    integer :: i,  j,   k       ! triplet indices 
    integer :: is, js,  ks      ! triplet indices  obtained by symmetry

    integer :: G(3)

    integer :: ns

    integer :: ikpt, ikpts      ! singlet index, singlet index obtained by symmetry

    integer :: switch, ii, jj, kk, nsp

    logical :: found(nspt) 
    real(dp), parameter :: eps=10E-7

    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    switch = 1  ! triplet-to-singlet coordinate

    found  = .false.
    nk_irr = 0        

    do i=1, nspt
       if (.not. found(i)) then
          found(i) = .true.
          nk_irr = nk_irr + 1  
          k_irr(1:3,nk_irr) = k_entire(1:3,i)
          equiv(i) = nk_irr 

          do ns=1,nsym
             k_rot = matmul(dble(s(:,:,ns)), k_irr(:,nk_irr))

             do j=1, nspt

                if (((abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)-1.0_dp)<eps).or. &
                     (abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)       )<eps)).and.&
                     ((abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)-1.0_dp)<eps).or. &
                     (abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)       )<eps)).and.&
                     ((abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)-1.0_dp)<eps).or. &
                     (abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)       )<eps)) ) then 

                   equiv(j)   = nk_irr 
                   found(j)   =.true.
                   symlink(j,1) = ns
                   symlink(j,2) = 0 
                   sym_G(:,j) = nint(k_entire(:,j) - k_rot(:))
                   cycle 
                endif

                if (TR_symmetry) then

                   if (((abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)-1.0_dp)<eps).or. &
                        (abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)       )<eps)).and.&
                        ((abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)-1.0_dp)<eps).or. &
                        (abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)       )<eps)).and.&
                        ((abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)-1.0_dp)<eps).or. &
                        (abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)       )<eps)) ) then
                      equiv(j) = nk_irr 
                      found(j) =.true.
                      symlink(j,1) = ns
                      symlink(j,2) = 1 
                      sym_G(:,j) = nint(k_entire(:,j) -(- k_rot(:)))

                   endif
                endif

             end do !j 

          enddo ! ns 

       end if !found

    enddo !i 


    do ikpt=1, nk_irr

       k_rot(:) = matmul(bg,k_irr(:,ikpt))

       dist1 = sum (k_rot(:)**2)

       do ii=-1,1
          do jj=-1,1
             do kk=-1,1

                k_rot(:) = matmul(bg,k_irr(:,ikpt))  + matmul(bg, dble((/ii,jj,kk/)))
                dist2 = sum (k_rot(:)**2) 

                if (dist2<dist1) then 
                   k_irr(:,ikpt) = matmul (transpose(at), k_rot)
                   dist1 = dist2
                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine find_the_irreducible_k_set_and_equiv

  subroutine find_the_irreducible_k_set (kmesh, kpoints_irr,nk_irr)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the canonical 
    ! 1BZ mesh.
    !
    !------------------------------------------------------------------

    use intw_reading, only : s,  nsym
    use intw_input_parameters, only : TR_symmetry !Only this, not equiv and symlink!.

    use intw_utility ! provides useful functions!

    implicit none

    !input

    real(dp) :: kpoints_irr  (3,nk1* nk2* nk3) 
    ! The irreducible kpoints in crystal triplet coordinates.
    ! The size of the array is nk1* nk2* nk3 instead of nk_irr;
    ! it is supposed that we still do not know the value of nk_irr

    real(dp) :: kmesh (3,nk1* nk2* nk3)
    ! The full kmesh, in canonical order

    ! output
    integer  :: nk_irr 
    ! how many points were found?

    !local variables
    real(dp):: k_rot(3), k_1BZ(3), dk(3)

    integer :: nkpt             ! The total number of points 
    integer :: i,  j,   k       ! triplet indices 
    integer :: is, js,  ks      ! triplet indices  obtained by symmetry

    integer :: G(3)

    integer :: ns

    integer :: ikpt, ikpts      ! singlet index, singlet index obtained by symmetry

    integer :: switch

    logical :: found(nk1*nk2*nk3) 


    nkpt   =  nk1*nk2*nk3


    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    switch = 1  ! triplet-to-singlet coordinate

    found  = .false.
    nk_irr = 0        


    ! loop on the whole mesh, in the appropriate order
    do i=1,nk1
       do j=1,nk2
          do k=1,nk3
             ! find scalar index of point (i,j,k) 
             call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,switch)
             ! operate on this point only if it has not already been found!
             if (.not. found(ikpt)) then

                ! it's found now. This point is part of the IBZ.
                found(ikpt) = .true.  

                nk_irr = nk_irr + 1        

                kpoints_irr(1,nk_irr) = dble(i-1)/nk1
                kpoints_irr(2,nk_irr) = dble(j-1)/nk2
                kpoints_irr(3,nk_irr) = dble(k-1)/nk3

                ! loop on all symmetry operations
                do ns=1,nsym    
                   !perform matrix product
                   ! CAREFUL! since the matrix is in crystal coordinates, 
                   ! and it acts in reciprocal space, the convention is :
                   !          k_rot(i) = sum_j s(i,j)*k(j)

                   k_rot = matmul(dble(s(:,:,ns)), kpoints_irr(:,nk_irr))

                   ! find what point in the 1BZ this corresponds to
                   call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)

                   ! check that k_1BZ+G = k_rot. If not, k_rot isn't on the mesh,
                   ! and the algorithm in "find_k_1BZ_and_G" cannot be trusted.
                   dk = k_rot-(k_1BZ+dble(G))

                   if (sqrt(dot_product(dk,dk)) < eps_8) then
                      ! what is the scalar index 
                      call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,switch)

                      if (.not. found(ikpts)) found(ikpts) = .true.

                   end if ! dk

                   ! Repeat, with Time-Reversal symmetry if present
                   if (TR_symmetry) then

                      k_rot = -k_rot

                      ! find what point in the 1BZ this corresponds to
                      call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)

                      dk = k_rot-(k_1BZ+dble(G))

                      if (sqrt(dot_product(dk,dk)) < eps_8) then
                         ! what is the scalar index 
                         call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,switch)

                         if (.not. found(ikpts)) found(ikpts) = .true.
                      end if ! dk

                   end if !TR_symmetry


                end do ! ns

             end if ! found(ikpt)     

          end do ! k
       end do ! j
    end do ! i 
  end subroutine find_the_irreducible_k_set

  subroutine find_entire_nice_BZ (nk1, nk2, nk3, nspt, ksvec) !(kmesh, kpoints_irr,nk_irr)

    !------------------------------------------------------------------
    ! input: private nk1 nk2 nk3, output full BZ nk>nk1*nk2*nk3 for tria_diag by AE.&IGdG.  
    !------------------------------------------------------------------

    use intw_reading, only : s,  nsym
    use intw_input_parameters, only : TR_symmetry !Only this, not equiv and symlink!.

    use intw_utility ! provides useful functions!

    implicit none

    !input

    integer, intent (in) :: nk1, nk2, nk3

    ! output
    real(dp), intent(out) :: ksvec(3,2*nk1*nk2*nk3)
    integer, intent(out)  :: nspt

    !local variables
    real(dp):: k_rot(3), k_1BZ(3), dk(3)

    integer :: nkpt             ! The total number of points 
    integer :: i,  j,   k, ii,jj,kk       ! triplet indices 
    integer :: is, js,  ks      ! triplet indices  obtained by symmetry

    integer :: G(3)

    integer :: ns

    integer :: ikpt, ikpts      ! singlet index, singlet index obtained by symmetry

    integer :: switch

    logical :: found(nk1*nk2*nk3) 

    real(dp) :: k_aux(3,48)

    real(dp),parameter  :: eps=1d-6

    integer :: nsp, iss, nk_irr
    real(dp) :: kpoints_irr  (3,nk1*nk2*nk3), dist1, dist2


    nkpt   =  nk1*nk2*nk3


    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    switch = 1  ! triplet-to-singlet coordinate

    found  = .false.
    nk_irr = 0        


    ! loop on the whole mesh, in the appropriate order
    do i=1,nk1
       do j=1,nk2
          do k=1,nk3
             ! find scalar index of point (i,j,k) 
             call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,switch)
             ! operate on this point only if it has not already been found!
             if (.not. found(ikpt)) then

                ! it's found now. This point is part of the IBZ.
                found(ikpt) = .true.  

                nk_irr = nk_irr + 1        

                kpoints_irr(1,nk_irr) = dble(i-1)/nk1
                kpoints_irr(2,nk_irr) = dble(j-1)/nk2
                kpoints_irr(3,nk_irr) = dble(k-1)/nk3

                ! loop on all symmetry operations
                do ns=1,nsym    
                   !perform matrix product
                   ! CAREFUL! since the matrix is in crystal coordinates, 
                   ! and it acts in reciprocal space, the convention is :
                   !          k_rot(i) = sum_j s(i,j)*k(j)

                   k_rot = matmul(dble(s(:,:,ns)), kpoints_irr(:,nk_irr))

                   ! find what point in the 1BZ this corresponds to
                   call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)

                   ! check that k_1BZ+G = k_rot. If not, k_rot isn't on the mesh,
                   ! and the algorithm in "find_k_1BZ_and_G" cannot be trusted.
                   dk = k_rot-(k_1BZ+dble(G))

                   if (sqrt(dot_product(dk,dk)) < eps_8) then
                      ! what is the scalar index 
                      call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,switch)

                      if (.not. found(ikpts)) found(ikpts) = .true.

                   end if ! dk

                   ! Repeat, with Time-Reversal symmetry if present
                   if (TR_symmetry) then

                      k_rot = -k_rot

                      ! find what point in the 1BZ this corresponds to
                      call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)

                      dk = k_rot-(k_1BZ+dble(G))

                      if (sqrt(dot_product(dk,dk)) < eps_8) then
                         ! what is the scalar index 
                         call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,switch)

                         if (.not. found(ikpts)) found(ikpts) = .true.
                      end if ! dk

                   end if !TR_symmetry


                end do ! ns

             end if ! found(ikpt)     

          end do ! k
       end do ! j
    end do ! i

    do ikpt=1, nk_irr

       k_rot(:) = matmul(bg,kpoints_irr(:,ikpt))

       dist1 = sum (k_rot(:)**2)

       do ii=-1,1
          do jj=-1,1
             do kk=-1,1

                k_rot(:) = matmul(bg,kpoints_irr(:,ikpt))  + matmul(bg, dble((/ii,jj,kk/)))
                dist2 = sum (k_rot(:)**2) 

                if (dist2<dist1) then 
                   kpoints_irr(:,ikpt) = matmul (transpose(at), k_rot)
                   dist1 = dist2
                endif

             enddo
          enddo
       enddo
    enddo


    nspt=0
    do ikpt=1, nk_irr
       !-print*, "ikpt", ikpt
       nsp =1
       nspt=nspt+1
       k_aux(:, nsp) =  kpoints_irr(:,ikpt)
       ksvec(:,nspt) =  k_aux(:, nsp)
       s_l: do is=1,nsym
          k_rot = matmul(dble(s(:,:,is)), kpoints_irr(:,ikpt))
          iss_l: do iss=1,nsp
             if (sum(abs(k_aux(:, iss)-k_rot(:)))<eps)  cycle s_l
          enddo iss_l
          nsp=nsp+1
          k_aux(:, nsp) =  k_rot(:)
          nspt=nspt+1
          ksvec(:,nspt) =  k_aux(:, nsp)
       enddo s_l
    enddo

    do ikpt=1, nspt
       write(123,"(100f12.6)") ksvec(:,ikpt)
    enddo

  end subroutine find_entire_nice_BZ


  subroutine irr_kp_grid_to_full (nk_1, nk_2, nk_3, nk_irr, k_irr, &
       equiv_l, symlink_l)
    !-------------------------------------------------------------------------
    ! Given a irreducible k set, this subroutine finds the connection 
    ! to the full BZ.
    !
    ! input  : nk1 nk2 nk3 nk_irr (n. of irr k points) and k_irr(3,nk_irr)
    ! output : equiv_l, symlink_l  
    !
    !-------------------------------------------------------------------------

    use intw_reading, only : s, nsym
    use intw_input_parameters, only : TR_symmetry !Only this, not equiv and symlink!.

    implicit none

    !input
    integer :: nk_1, nk_2, nk_3

    integer :: nk_irr                  ! N. of irreducible k points

    integer :: k_irr(3,nk_irr)         ! triplet indices of the irreducible k-points

    !output
    integer :: equiv_l  (nk_1* nk_2* nk_3)   ! For a given k point in the entire BZ, 
    ! gives the symmetry related
    ! k point index in the irreducible k set

    integer :: symlink_l(nk_1* nk_2* nk_3,2) ! Symmetry link .


    !local
    integer :: nkpt, ikpts
    integer :: n2n3, n1n2, n1n3

    integer :: il

    integer :: irr

    integer :: i,  j,  k
    integer :: is, js, ks

    integer :: aux_sym(3)

    integer :: ns 
    integer :: switch

    logical :: done(nk_1* nk_2* nk_3)


    switch = 1 ! tripet - to - singlet
    nkpt   =  nk_1*nk_2*nk_3

    n1n3   =  nk_1*nk_3
    n2n3   =  nk_2*nk_3
    n1n2   =  nk_1*nk_2

    done(:)=  .false.


    ! CAREFUL!!! 
    !
    ! There is great potential for confusion here. 
    ! The kpoints are well labeled by either a triplet or a singlet
    ! index, and one can switch between them using switch_indices.
    ! HOWEVER, when we consider a set of irreducible points, these
    ! points may be listed sequentially, but the index of these
    ! points have NOTHING to do with the singlet index!
    !
    ! an example to make this clear: Consider an 8 x 8 x 8 mesh.
    ! there are 512 k-points in the zone, and (say) 29 irreducible
    ! k-points. If one performs a nscf QE computation to obtain the 
    ! wavefunctions of only these 29 k-points, the QE data folders
    ! will be labeled K0001, K0002, .., K00029. These numbers (1..29)
    ! however, are NOT the singlet indices of the points!
    !
    ! below, "k_irr" is the triplet index of a given irreducible k-point
    ! "irr" is a sequential index and ikpt is the singlet index 
    ! corresponding to the triplet index.
    !
    ! This subroutine will thus return the array
    ! equiv(ikpt)   = irr  
    !
    ! with irr representing the sequential number. This is useful when
    ! starting a computation with only the ireducible QE data.

    do irr = 1,nkpt

       i = k_irr(1,irr) 
       j = k_irr(2,irr) 
       k = k_irr(3,irr) 


       do ns  = 1,nsym

          do il= 1,3
             aux_sym(il) =    n2n3*s(il,1,ns) * (i-1)  & 
                  +  n1n3*s(il,2,ns) * (j-1)  & 
                  +  n1n2*s(il,3,ns) * (k-1)  
          enddo
          ! CAREFUL !!! "modulo(a,p)" is always positive for p positive; 
          !              "mod(a,p)" can be negative
          is = modulo ( aux_sym(1), nkpt )/n2n3 + 1
          js = modulo ( aux_sym(2), nkpt )/n1n3 + 1
          ks = modulo ( aux_sym(3), nkpt )/n1n2 + 1

          call switch_indices(nk_1,nk_2,nk_3,ikpts,is,js,ks,switch)

          if (.not.done(ikpts)) then
             done      (ikpts  ) = .true.
             equiv_l   (ikpts  ) = irr
             !         symlink_l (ikpts,1) = ns
             symlink_l(ikpts,1)  =  inverse_indices(ns)
             symlink_l (ikpts,2) = 0
          endif

          ! If time reversal symmetry is present, use it!
          if (TR_symmetry) then
             is = -is  
             js = -js  
             ks = -ks  
             call switch_indices(nk_1,nk_2,nk_3,ikpts,is,js,ks,switch)

             if (.not.done(ikpts)) then
                done      (ikpts  ) = .true.
                equiv_l   (ikpts  ) = irr
                !             symlink_l (ikpts,1) = ns
                symlink_l(ikpts,1) =  inverse_indices(ns)
                symlink_l (ikpts,2) = 1
             endif

          end if

       enddo ! ns
    enddo ! irr


  end subroutine irr_kp_grid_to_full

  subroutine calculate_star_r(v, vstar, nstar, symop)
    real(dp), intent(in) :: v(3)
    real(dp), intent(out) ::vstar(3,48)
    integer,intent(out) :: nstar, symop(48)
    integer :: isym, i
    real(dp) :: vrot(3)


    nstar=1
    vstar(1:3,nstar) = v(1:3)
    symop(1)=1

    do isym=1, nsym
       vrot(:) = matmul(dble(s(:,:,isym)), v(:))

       do i=1,nstar
          if  ( sum(abs(vrot(:)- vstar(1:3,i)))<10E-5  ) then
             goto 1987 
          end if
       enddo

       nstar = nstar + 1
       vstar(1:3,nstar) = vrot(1:3)
       symop(nstar)=isym

1987   continue
    enddo


  end subroutine calculate_star_r

  subroutine calculate_star(v, vstar, nstar, symop)
    integer, intent(in) :: v(3)
    integer, intent(out) ::vstar(3,48), symop(48)
    integer,intent(out) :: nstar
    integer :: isym, i
    real(dp) :: vrot(3)


    nstar=1
    vstar(1:3,nstar) = v(1:3)
    symop(1)=1

    do isym=1, nsym
       vrot(:) = nint(matmul(dble(s(:,:,isym)), v(:)))

       do i=1,nstar
          if  ( sum(abs(vrot(:)- vstar(1:3,i)))<10E-5  ) then
             goto 1984 
          end if
       enddo

       nstar = nstar + 1
       vstar(1:3,nstar) = vrot(1:3)
       symop(nstar)=isym

1984   continue
    enddo


  end subroutine calculate_star

  subroutine echo_symmetry_1BZ(nk_1,nk_2,nk_3,nk_irr,equiv,symlink)
    !--------------------------------------------------------------------------------
    !
    ! simple writing routine, for testing purposes
    !
    !--------------------------------------------------------------------------------
    implicit none

    integer  :: nk_1, nk_2, nk_3, nkr, nk_irr
    integer  :: io_unit
    integer  :: i

    integer  :: equiv  (nk_1* nk_2* nk_3)  
    integer  :: symlink(nk_1* nk_2* nk_3,2)

    nkr = nk_1*nk_2*nk_3

    io_unit = find_free_unit()
    open(io_unit,file='equivalence_of_k_points.test')
    write(io_unit,'(a,i4)')'          nk_irr  = ',nk_irr
    write(io_unit,'(a)')'----------------------------'

    write(io_unit,'(a)')'          Nk              Equivalent      Symmetry           Time reversal. '
    write(io_unit,'(a)')'                          in reduced l.   connection.        not used:0, used :1  '
    write(io_unit,'(a)')'-----------------------------------------------------------------------------------------------'

    do i=1,nkr
       write(io_unit,100) i,  equiv(i), symlink(i,:)
    enddo

    close(io_unit)

100 format(i10,12x,i10,12x,i10,12x,i10)
  end subroutine echo_symmetry_1BZ

  subroutine rot_atoms(nat,nsym, tau  )
    implicit none
    integer,  intent(in)  :: nat, nsym
    real(dp), intent(in)  :: tau (3,nat)
    !out global
!    integer,  intent  (out) :: rtau_index(nat, nsym)
!    real(dp), intent (out) :: rtau(3,nsym,nat)
!    real(dp), intent (out) :: tau_cryst(3,nat)
!    real(dp), intent (out) :: rtau_cryst(3,nsym,nat)

    integer              :: i,j,k,h
    real(dp)             :: vec(3),zero(3),rat(3)
    logical              :: eqatom
    integer              :: isym
    integer              :: nr(3)
    integer              :: a_index, b_index,  na, nb     
  
    integer :: ipol, jpol, kpol, lpol

    real (dp), parameter :: epsat= 1E-3

    real(dp) :: s_cart(3,3, nsym)

    do isym=1,nsym
     do ipol = 1, 3
     do jpol = 1, 3
        s_cart (ipol, jpol, isym) = 0.d0
        do kpol = 1, 3
           do lpol = 1, 3
              s_cart (ipol, jpol,isym) = s_cart (ipol, jpol, isym) + at (ipol, kpol) * &
                   s (lpol, kpol, isym)  * bg (jpol, lpol)
           enddo
        enddo
     enddo
    enddo
    enddo

    nr = (/nr1,nr2,nr3/)

    rtau_index = -11

    do i = 1, nat
       tau_cryst(:,i) = matmul(ainv(at),  tau (:, i))
    enddo

    do i=1,nat
       do isym=1,nsym
          do ipol=1,3

             rtau_cryst(ipol,isym,i) = s(1,ipol,isym) * tau_cryst(1,i) &
                              + s(2,ipol,isym) * tau_cryst(2,i) &
                              + s(3,ipol,isym) * tau_cryst(3,i)

             rtau_cryst(ipol,isym,i) = rtau_cryst (ipol,isym,i) - dble(ftau(ipol,isym))

          end do

       enddo
    enddo

    do i=1,nat
       do isym=1,nsym
          do j=1,nat
             if (eqvect (rtau_cryst(:,isym,i), tau_cryst(:,j), (/0.d0,0.0d0,0.0d0/) ) )  rtau_index(i, isym)=j
          enddo
      enddo
    enddo
 
    do isym=1,nsym
       do na = 1, nat
          a_index = rtau_index(na,isym)
          do  h=1,3
             rtau(h,isym,na) = s_cart(1,h, isym) * tau(1,na) &
                  + s_cart(2,h, isym) * tau(2,na) &
                  + s_cart(3,h, isym) * tau(3,na)

             rtau(h,isym,na)=rtau(h,isym,na)-tau(h,a_index)
          end do
       enddo
    enddo

    do i=1,nat
       do isym=1,nsym
          if (rtau_index(i,isym).eq.0) then
             write(*,*)'ERROR in rot_at: At least one atom does not map properly under sym. op.',isym,'atom:',i
          endif
       enddo
    enddo

  end subroutine rot_atoms

  subroutine rotate_wfc (wfc_k_irr,list_iG_k_irr,wfc_k, &
       list_iG_k, i_sym, sym_inverse_l, ftau_l, G_sym_l)
    !------------------------------------------------------------------
    ! This subroutine takes in the periodic part of a wavefunction psi_{nk} and 
    ! returns the periodic part of psi_{n Rk}, where Rk is the rotated k-vector.
    !
    ! The wavefunctions have the form 
    !			psi_{nk}(r) = e^{ikr}/sqrt{V} u_{nk}(r)
    !			u_{nk}(r)   = \sum_G e^{iGr} u_{nk}(G).
    !
    !
    !  a crystal rotation-like symmetry can be expressed as
    !		S = { R | f }
    !	where R is a rotation and f a fractional translation.
    !
    ! Note that symmetry is implemented in a very confusing way in Quantum Espresso.

    ! On the one hand, section A.4 of the Quantum Espresso reference paper,
    ! 		"Quantum Espresso: a modular and open-source software project for 
    !  			quantum simulations of materials",
    !     suggests the convention:
    !           r' =  { R | f } r = R * ( r + f )  !! NOTE: this is NOT the usual textbook definition!
    !
    ! HOWEVER: poking around in the code suggests that the convention actually
    !	     used in the code is 
    !           r' =  { R | f } r = R * r - f 
    !
    !	(This only matters in non-symmorphic systems)
    ! In what follows, the second convention will be used.
    !
    ! assumptions:  
    !			-  i_sym is the index of the symmetry operation
    !			-  sym_inverse_l is a point group operation:
    !                        it is the INVERSE of the actual operation; this is the
    !                        appropriate operator to act on k, in crystal coordinates
    !			-  k_irr is a k-point in the IBZ
    !			-  sym_inverse_l * k_irr  = k + G_sym, with k in the 1BZ
    !
    !	applying the point group operation yields
    !
    !           u_{nk}(sym_l*G+G_sym) =  e^{i R*G*tau} u_{nk_irr}(G)
    !
    !
    !------------------------------------------------------------------
    use intw_useful_constants 
    use intw_reading
    use intw_utility 
    use intw_fft,  only: find_iG
    !
    implicit none
    !
    ! input::
    complex(dp)    ::   wfc_k_irr(nG_max,num_bands,nspin)  ! wfc at point k_irr  in the IBZ
!    complex(dp)    ::   wfc_k_irr(nG_max,nbands,nspin)  ! wfc at point k_irr  in the IBZ
    integer        ::   list_iG_k_irr(nG_max)     ! G vector indices for k_irr 
    integer        ::   i_sym                     ! index of the symmetry operation
    integer        ::   sym_inverse_l (3,3)       ! inverse point group operation, 
    ! the one acting on k (in crystal coordinates)
    integer        ::   G_sym_l(3)                ! G vector such that  R*k + G_sym_l = sym_l * k_irr

    real(dp)       ::   ftau_l(3)                 ! fractional translation associated with point group operation
    ! output::
    complex(dp)    ::   wfc_k(nG_max,num_bands,nspin)      ! rotated wfc at point k, in the 1BZ
!    complex(dp)    ::   wfc_k(nG_max,nbands,nspin)      ! rotated wfc at point k, in the 1BZ
    integer        ::   list_iG_k(nG_max)         ! G vector indices for k, sorted


    ! local
    integer        ::  p_i,   i, j,     alpha
    integer        ::  iG_k_irr,   iG_k

    integer        ::  G_k_irr(3)              ! a vector for k in the IBZ
    integer        ::  RG_k_irr(3)             ! ( symmetry operation )* G_k
    integer        ::  G_k(3)                  ! a vector for Rk, the point in the 1BZ

    integer        ::   permutations(nG_max)   ! index permutation which orders list_G_k 

    integer        ::  nb, ipol, jpol 

    integer        ::  nG         ! counter on the number of G vectors in the array

    complex(dp)    ::  phases(nG_max)

    real(dp)       :: axis(3), angle

    complex(dp)    :: S_u(2,2) 

!Peio
!The variable we use instead of num_bands
    integer        :: nbands_loc
!Peio

    
    nbands_loc = num_bands

    phases    =  cmplx_0
    wfc_k     =  cmplx_0

    list_iG_k          =  0
    nG                 =  0
    permutations       =  0

    ! loop on all G_k_irr, the coefficients of the wavefunction at the IBZ k point


    do i=1,nG_max

       iG_k_irr = list_iG_k_irr(i)

       if (iG_k_irr == 0) exit  ! the index array is zero-padded at the end. 
       nG = nG+1 ! only increment if iG_k_irr /= 0!


       G_k_irr = gvec(:,iG_k_irr)


       do alpha=1,3

          ! REMEMBER the CONVENTION: k_rot(i) = sum_j s(i,j) k(j) 
          RG_k_irr(alpha)  =  sym_inverse_l(alpha,1)*G_k_irr(1)    &
               + sym_inverse_l(alpha,2)*G_k_irr(2)    &
               + sym_inverse_l(alpha,3)*G_k_irr(3)    


       end do


       phases(nG) =  exp(  cmplx_I*tpi*  ( dble(RG_k_irr(1)) * ftau_l(1)  &
            + dble(RG_k_irr(2)) * ftau_l(2)  &
            + dble(RG_k_irr(3)) * ftau_l(3)  ))

       G_k(:) = RG_k_irr(:)+G_sym_l(:)


       call find_iG(G_k,iG_k)

       list_iG_k(nG) = iG_k

    end do

    ! There is no guarantee that the indices in list_iG_k will be sorted in ascending
    ! order! This is not an absolute necessity, but it would be nice and consistent for
    ! the indices to be sorted. 
    ! Sort the indices using a canned heap sort subroutine.

    call HPSORT(nG,list_iG_k,permutations)

    ! To understand how this works, consider an example: 
    !
    !            i      f(i)        iG(i)   permutation(i)
    !            ---------------------------------------               
    !            1      0.1         4            2   
    !            2      0.2         1            4
    !            3      0.3         3            3
    !            4      0.4         2            1
    !
    !            j   sort(iG)(j)    sort(f)(j)
    !            ------------------------------------
    !            1      1               0.2  
    !            2      2               0.4
    !            3      3               0.3
    !            4      4               0.1
    !
    !             ===> sort(f) (j)  =   f( permutation(j) )
    ! 

    ! finally, populate the rotated wave function

    do i= 1, nG

       p_i = permutations(i)
       ! compute the wfc element
       do nb = 1,nbands_loc
!       do nb=1,nbands
          do ipol=1,nspin
             wfc_k(i,nb,ipol) =  phases(p_i)*wfc_k_irr(p_i,nb,ipol)
          enddo
       end do

    end do

    if (nspin==2) then

       S_u(:,:)  = spin_symmetry_matrices(:,:,i_sym)

       wfc_k_irr = wfc_k
       wfc_k     = cmplx_0

       do i= 1, nG
          do nb = 1,nbands_loc
!          do nb=1,nbands
             do ipol=1,nspin
                do jpol=1,nspin
                   wfc_k(i,nb,ipol) = wfc_k(i,nb,ipol) + S_u(ipol,jpol) * wfc_k_irr(i,nb,jpol)
                enddo
             end do
          end do
       end do

    endif ! If non-collinear

  end subroutine rotate_wfc

  subroutine rotate_wfc_test (wfc_k_irr,list_iG_k_irr,wfc_k, &
       list_iG_k, i_sym, sym_inverse_l, ftau_l, G_sym_l)
    !------------------------------------------------------------------
    ! This subroutine takes in the periodic part of a wavefunction psi_{nk} and 
    ! returns the periodic part of psi_{n Rk}, where Rk is the rotated k-vector.
    !
    ! The wavefunctions have the form 
    !			psi_{nk}(r) = e^{ikr}/sqrt{V} u_{nk}(r)
    !			u_{nk}(r)   = \sum_G e^{iGr} u_{nk}(G).
    !
    !
    !  a crystal rotation-like symmetry can be expressed as
    !		S = { R | f }
    !	where R is a rotation and f a fractional translation.
    !
    ! Note that symmetry is implemented in a very confusing way in Quantum Espresso.

    ! On the one hand, section A.4 of the Quantum Espresso reference paper,
    ! 		"Quantum Espresso: a modular and open-source software project for 
    !  			quantum simulations of materials",
    !     suggests the convention:
    !           r' =  { R | f } r = R * ( r + f )  !! NOTE: this is NOT the usual textbook definition!
    !
    ! HOWEVER: poking around in the code suggests that the convention actually
    !	     used in the code is 
    !           r' =  { R | f } r = R * r - f 
    !
    !	(This only matters in non-symmorphic systems)
    ! In what follows, the second convention will be used.
    !
    ! assumptions:  
    !			-  i_sym is the index of the symmetry operation
    !			-  sym_inverse_l is a point group operation:
    !                        it is the INVERSE of the actual operation; this is the
    !                        appropriate operator to act on k, in crystal coordinates
    !			-  k_irr is a k-point in the IBZ
    !			-  sym_inverse_l * k_irr  = k + G_sym, with k in the 1BZ
    !
    !	applying the point group operation yields
    !
    !           u_{nk}(sym_l*G+G_sym) =  e^{i R*G*tau} u_{nk_irr}(G)
    !
    !
    !------------------------------------------------------------------
    use intw_useful_constants 
    use intw_reading
    use intw_utility 
    use intw_fft,  only: find_iG
    !
    implicit none
    !
    ! input::
    complex(dp), intent(in)    ::   wfc_k_irr(nG_max,num_bands,nspin)  ! wfc at point k_irr  in the IBZ
!    complex(dp), intent(in)    ::   wfc_k_irr(nG_max,nbands,nspin)  ! wfc at point k_irr  in the IBZ
    integer, intent(in)        ::   list_iG_k_irr(nG_max)     ! G vector indices for k_irr 
    integer, intent(in)        ::   i_sym                     ! index of the symmetry operation
    integer, intent(in)       ::   sym_inverse_l (3,3)       ! inverse point group operation, 
    ! the one acting on k (in crystal coordinates)
    integer, intent(in)        ::   G_sym_l(3)                ! G vector such that  R*k + G_sym_l = sym_l * k_irr

    real(dp), intent(in)       ::   ftau_l(3)                 ! fractional translation associated with point group operation
    ! output::
    complex(dp), intent(out)    ::   wfc_k(nG_max,num_bands,nspin)      ! rotated wfc at point k, in the 1BZ
!    complex(dp), intent(out)    ::   wfc_k(nG_max,nbands,nspin)      ! rotated wfc at point k, in the 1BZ
    integer, intent(out)        ::   list_iG_k(nG_max)         ! G vector indices for k, sorted


    ! local
    complex(dp)                ::   wfc_k_aux(nG_max,num_bands,nspin)
!    complex(dp)                ::   wfc_k_aux(nG_max,nbands,nspin)
    integer        ::  p_i,   i, j,     alpha
    integer        ::  iG_k_irr,   iG_k

    integer        ::  G_k_irr(3)              ! a vector for k in the IBZ
    integer        ::  RG_k_irr(3)             ! ( symmetry operation )* G_k
    integer        ::  G_k(3)                  ! a vector for Rk, the point in the 1BZ

    integer        ::   permutations(nG_max)   ! index permutation which orders list_G_k 

    integer        ::  nb, ipol, jpol 

    integer        ::  nG         ! counter on the number of G vectors in the array

    complex(dp)    ::  phases(nG_max)

    real(dp)       :: axis(3), angle

    complex(dp)    :: S_u(2,2) 

!Peio
!The variable we use instead of num_bands
    integer        :: nbands_loc
    integer        :: ibnd
!Peio


    nbands_loc=num_bands

    phases    =  cmplx_0
    wfc_k     =  cmplx_0

    list_iG_k          =  0
    nG                 =  0
    permutations       =  0

    ! loop on all G_k_irr, the coefficients of the wavefunction at the IBZ k point

!Peio
!do i=1,nG_max
!   write(1111,'(a,i8,a,i8)') 'iG= ',i,' igk(iG)= ',list_iG_k_irr(i)
!enddo
!do ibnd=1,nbands_loc
!   do ipol=1,nspin
!      write(1112,'(a,i4,a,i4)')'ibnd= ',ibnd,' ipol= ',ipol
!      write(1112,*)'    ir   igk(ir)   wfc(ir)'
!      do i=1,nG_max
!         write(1112,'(i8,i8,f12.6)') i,list_iG_k_irr(i),abs(wfc_k_irr(i,ibnd,ipol))
!      enddo
!   enddo
!enddo
!stop
!Peio

    do i=1,nG_max

       iG_k_irr = list_iG_k_irr(i)

       if (iG_k_irr == 0) exit  ! the index array is zero-padded at the end. 
       nG = nG+1 ! only increment if iG_k_irr /= 0!


       G_k_irr = gvec(:,iG_k_irr)


       do alpha=1,3

          ! REMEMBER the CONVENTION: k_rot(i) = sum_j s(i,j) k(j) 
          RG_k_irr(alpha)  =  sym_inverse_l(alpha,1)*G_k_irr(1)    &
               + sym_inverse_l(alpha,2)*G_k_irr(2)    &
               + sym_inverse_l(alpha,3)*G_k_irr(3)    


       end do

       phases(nG) =  exp(  cmplx_I*tpi*  ( dble(RG_k_irr(1)) * ftau_l(1)  &
            + dble(RG_k_irr(2)) * ftau_l(2)  &
            + dble(RG_k_irr(3)) * ftau_l(3)  ))

       G_k(:) = RG_k_irr(:)+G_sym_l(:)

       call find_iG(G_k,iG_k)

       list_iG_k(nG) = iG_k

    end do

!Peio
!do i=1,nG_max
!   write(1111,'(a,i8,a,i8,a,i8)') 'iG= ',i,' igk(iG)= ',list_iG_k(i),' perm(iG)= ',permutations(i)
!enddo
!Peio

    call HPSORT(nG,list_iG_k,permutations)

!Peio
!do i=1,nG_max
!   write(1112,'(a,i8,a,i8,a,i8)') 'iG= ',i,' igk(iG)= ',list_iG_k(i),' perm(iG)= ',permutations(i)
!enddo
!Peio

    do i= 1, nG_max
       ! compute the wfc element
       p_i=permutations(i)
       do nb = 1,nbands_loc
!       do nb = 1,nbands
          do ipol=1,nspin
             if (p_i==0) then
                wfc_k(i,nb,ipol) =  cmplx_0
             else
                wfc_k(i,nb,ipol) =  phases(p_i)*wfc_k_irr(p_i,nb,ipol)
             endif 
          enddo
       end do

    end do

    if (nspin==2) then

       S_u(:,:)  = spin_symmetry_matrices(:,:,i_sym)

       wfc_k_aux = wfc_k
       wfc_k     = cmplx_0

       do i= 1, nG
          do nb = 1,nbands_loc
!          do nb = 1,nbands
             do ipol=1,nspin
                do jpol=1,nspin
                   wfc_k(i,nb,ipol) = wfc_k(i,nb,ipol) + S_u(ipol,jpol) * wfc_k_aux(i,nb,jpol)
                enddo
             end do
          end do
       end do

    endif ! If non-collinear
  end subroutine rotate_wfc_test

  subroutine get_psi_general_k(kpoint,use_IBZ,list_iG,wfc_k,QE_eig)

    !----------------------------------------------------------------------------!
    !     This subroutine is a driver which performs all the necessary operations
    !     to obtain the periodic part of the wavefunction at the specified k point.
    !
    !     Asier&&Idoia 24 06 2014 Compared to "get_psi_k", here the input 
    !     is a general k(3) (dp, crystal) which may be far from original kmesh.
    !
    !     It will extract the appropriate information from the QE folder, and
    !     rotate the wfc if necessary.
    !
    !     If use_IBZ is set to TRUE, the subroutine will obtain the wavefunction
    !     by rotation. If it is False, the subroutine will try to fetch
    !     the wavefunction from the QE directory; if this is not possible,
    !     it will throw an error.
    !----------------------------------------------------------------------------!

    implicit none

    real(dp), intent(in) :: kpoint(3)
    logical, intent(in)        :: use_IBZ


    complex(dp), intent(out)   :: wfc_k(nG_max,num_bands,nspin)
!    complex(dp), intent(out)   :: wfc_k(nG_max,nbands,nspin)
    integer, intent(out)       :: list_iG(nG_max)
    real(dp), intent(out)      :: QE_eig(num_bands)
!    real(dp), intent(out)      :: QE_eig(nbands)

    complex(dp)    :: wfc_k_irr(nG_max,num_bands,nspin)
!    complex(dp)    :: wfc_k_irr(nG_max,nbands,nspin)
    integer        :: list_iG_irr(nG_max)


    integer        :: ikpt, ikpt_irr, i_folder

    integer        :: i_sym, TR
    integer        :: G_sym(3), G_plus(3)
    integer        :: sym(3,3)


    real(dp)       :: ftau_sym(3)
    integer        :: ipol, i_1bz, j_1bz, k_1bz, i_, j_, k_

    real(dp)           :: kpoint_1BZ(3)

    complex(dp)   :: wfc_k_1(nG_max,num_bands,nspin)
!    complex(dp)   :: wfc_k_1(nG_max,nbands,nspin)

!Peio
!The variable we use instead of num_bands
    integer :: nbands_loc
!Peio


    nbands_loc=num_bands

    if (.not. use_IBZ .and. .not. full_mesh) then
       ! the parameters are not consistent
       write(*,*) '*****************************************************'
       write(*,*) '*                 subroutine get_psi_k              *'
       write(*,*) '*               --------------------------          *'
       write(*,*) '*      A wavefunction obtained without rotation     *'
       write(*,*) '*      was requested, but the QE folders do not     *'
       write(*,*) '*      contain a full mesh! Review code...          *'
       write(*,*) '*        Program continues at your own risks!       *'
       write(*,*) '*****************************************************'
       !stop
    end if

    call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i_1bz,j_1bz, k_1bz, kpoint_1bz, G_plus)

    call switch_indices (nk1,nk2,nk3,ikpt,i_1bz,j_1bz,k_1bz, +1)

    if (.not. use_IBZ ) then
       ! do not use the IBZ; just fetch the wavefunction from the QE folders

       i_folder = QE_folder_nosym(ikpt) 
       G_sym    = nosym_G(:,ikpt) + G_plus(:) !Asier&&Idoia 24 06 2014 
       ftau_sym = ZERO
       sym      = s       (:,:,identity_matrix_index)

       call get_K_folder_data(i_folder,list_iG_irr,wfc_k_irr,QE_eig)


       call rotate_wfc_test (wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
            identity_matrix_index, sym, ftau_sym, G_sym)

    else

       ! Use the IBZ and symmetry!

       ! identify the right folder
       i_folder = QE_folder_sym(ikpt) 

       call get_K_folder_data(i_folder,list_iG_irr,wfc_k_irr,QE_eig)

       ! The symmetry which takes k_irr to k 
       i_sym    = symlink(ikpt,1)
       TR       = symlink(ikpt,2)
       G_sym    =  sym_G(:,ikpt) !+ G_plus(:) !Asier&&Idoia 24 06 2014
       ftau_sym = ftau(:,i_sym)
       sym      = s(:,:,inverse_indices(i_sym))

       ! debugging test
       !        call rotate_wfc (wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
       !			inverse_indices(i_sym), sym, ftau_sym, G_sym)

       call rotate_wfc_test (wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
            i_sym, sym, ftau_sym, G_sym)

       ! If time-reversal is present, the wavefunction currently stored 
       ! in wfc_k is actually for (-k). Complex conjugation must now
       ! be applied to recover wfc_k.
       if ( TR == 1)  then
          call apply_TR_to_wfc (wfc_k,list_iG)
       end if

    end if

    list_iG_irr=list_iG

    call wfc_by_expigr (kpoint, nbands_loc, npol, ng_max, list_iG_irr, list_iG, wfc_k, G_plus)
!    call wfc_by_expigr (kpoint, nbands, npol, ng_max, list_iG_irr, list_iG, wfc_k, G_plus)
    !(list_iG_irr, list_iG, G_plus)

  end subroutine get_psi_general_k

  subroutine get_psi_k (ikpt,use_IBZ,list_iG,wfc_k,QE_eig)
    !----------------------------------------------------------------------------!
    !     This subroutine is a driver which performs all the necessary operations
    !     to obtain the periodic part of the wavefunction at the specified k point.
    !
    !     It will extract the appropriate information from the QE folder, and
    !     rotate the wfc if necessary.
    !
    !     If use_IBZ is set to TRUE, the subroutine will obtain the wavefunction
    !     by rotation. If it is False, the subroutine will try to fetch
    !     the wavefunction from the QE directory; if this is not possible,
    !     it will throw an error.
    !----------------------------------------------------------------------------!

    implicit none

    logical        :: use_IBZ
    integer        :: ikpt, ikpt_irr, i_folder

    integer        :: i_sym, TR
    integer        :: G_sym(3)
    integer        :: sym(3,3)

    complex(dp)    :: wfc_k(nG_max,num_bands,nspin)
!    complex(dp)    :: wfc_k(nG_max,nbands,nspin)
    complex(dp)    :: wfc_k_irr(nG_max,num_bands,nspin)
!    complex(dp)    :: wfc_k_irr(nG_max,nbands,nspin)

    integer        :: list_iG(nG_max)
    integer        :: list_iG_irr(nG_max)

    real(dp)       :: QE_eig(num_bands)
!    real(dp)       :: QE_eig(nbands)
    real(dp)       :: ftau_sym(3)

    integer        :: ipol


    if (.not. use_IBZ .and. .not. full_mesh) then
       ! the parameters are not consistent
       write(*,*) '*****************************************************'
       write(*,*) '*                 subroutine get_psi_k              *'
       write(*,*) '*               --------------------------          *'
       write(*,*) '*      A wavefunction obtained without rotation     *'
       write(*,*) '*      was requested, but the QE folders do not     *'
       write(*,*) '*      contain a full mesh! Review code...          *'
       write(*,*) '*        Program continues at your own risks!       *'
       write(*,*) '*****************************************************'
       !stop
    end if

    if (.not. use_IBZ ) then
       ! do not use the IBZ; just fetch the wavefunction from the QE folders

       i_folder = QE_folder_nosym(ikpt) 
       G_sym    = nosym_G(:,ikpt)
       ftau_sym = ZERO
       sym      = s       (:,:,identity_matrix_index)

       call get_K_folder_data(i_folder,list_iG_irr,wfc_k_irr,QE_eig)


       call rotate_wfc (wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
            identity_matrix_index, sym, ftau_sym, G_sym)

    else
       ! Use the IBZ and symmetry!

       ! identify the right folder
       i_folder = QE_folder_sym(ikpt) 

       call get_K_folder_data(i_folder,list_iG_irr,wfc_k_irr,QE_eig)

       ! The symmetry which takes k_irr to k 
       i_sym    = symlink(ikpt,1)
       TR       = symlink(ikpt,2)
       G_sym    = sym_G(:,ikpt)
       ftau_sym = ftau(:,i_sym)
       sym      = s(:,:,inverse_indices(i_sym))

       ! debugging test
       !        call rotate_wfc (wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
       !			inverse_indices(i_sym), sym, ftau_sym, G_sym)

       call rotate_wfc (wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
            i_sym, sym, ftau_sym, G_sym)

       ! If time-reversal is present, the wavefunction currently stored 
       ! in wfc_k is actually for (-k). Complex conjugation must now
       ! be applied to recover wfc_k.
       if ( TR == 1)  then
          call apply_TR_to_wfc (wfc_k,list_iG)
       end if

    end if

  end subroutine get_psi_k

  subroutine intw_check_mesh(nk_1,nk_2,nk_3,kmesh,k_irr,nk_irr,kpoints_QE)
    !----------------------------------------------------------------------------!
    !     This subroutine compares the mesh points kmesh to the global array
    !     kpoints (which comes from the QE folders) in order to determine if 
    !             1) kpoints = kmesh
    !             2) kpoints = irreducible zone 
    !             3) none of the above
    !----------------------------------------------------------------------------!
    implicit none

    integer   :: nk_1,nk_2,nk_3, nkmesh 
    real(dp)  :: kmesh(3,nk_1*nk_2*nk_3) 

    real(dp)  :: kpoints_QE(3,nkpoints_QE)

    integer   :: nk_irr 
    integer   :: k_irr(3,nk_irr) 

    real(dp)  :: kpt(3)
    integer   :: ikpt 
    real(dp)  :: norm

    nkmesh    = nk_1*nk_2*nk_3

    full_mesh = .true. 
    IBZ       = .true. 

    ! First, establish if a full mesh is present
    if (nkpoints_QE /= nkmesh) then
       full_mesh = .false.
    else
       do ikpt = 1, nkmesh

          norm =  sqrt( (kpoints_QE(1,ikpt)-kmesh(1,ikpt))**2  &
               + (kpoints_QE(2,ikpt)-kmesh(2,ikpt))**2  &
               + (kpoints_QE(3,ikpt)-kmesh(3,ikpt))**2 )

          if (norm > eps_8) then
             full_mesh = .false.
             exit
          end if

       end do

    end if

    ! if a full mesh is NOT present, is the IBZ present? 
    if (full_mesh .or. (nkpoints_QE /= nk_irr) ) then
       IBZ = .false.
    else
       do ikpt = 1,nk_irr
          kpt(1) = dble(k_irr(1,ikpt)-1)/nk_1
          kpt(2) = dble(k_irr(2,ikpt)-1)/nk_2
          kpt(3) = dble(k_irr(3,ikpt)-1)/nk_3

          norm   =  sqrt( (kpoints_QE(1,ikpt)-kpt(1))**2  &
               + (kpoints_QE(2,ikpt)-kpt(2))**2  &
               + (kpoints_QE(3,ikpt)-kpt(3))**2 )

          if (norm > eps_8) then
             IBZ = .false.
             exit
          end if

       end do
    end if

  end subroutine intw_check_mesh

  subroutine get_G_shells(number_G_shells,nG_shell_max)
    !----------------------------------------------------------------------------!
    !     This subroutine computes the index nG_shell_max such that 
    !     gvec(1:nG_shell_max) represents the first number_G_shells shells. 
    !
    !     It is assumed that gvec is already sorted according to the length
    !     of the G vectors.
    !----------------------------------------------------------------------------!
    use intw_reading,   only: gvec, ngm, bg

    implicit none

    integer   :: number_G_shells, nG_shell_max

    integer   :: nshell
    integer   :: iG

    integer   :: iflag, nvec

    real(dp)  :: G_cart(3)

    real(dp)  :: norm2, old_norm2

    nshell       = 1 
    nG_shell_max = 1

    nvec   = 1
    iflag  = 1

    old_norm2 = ZERO

    do iG = 2,ngm

       G_cart = dble(gvec(:,iG))
       call cryst_to_cart (nvec, G_cart, bg, iflag)

       norm2   = G_cart(1)**2+G_cart(2)**2+G_cart(3)**2

       if ( abs(norm2-old_norm2) < eps_8) then
          nG_shell_max  = nG_shell_max + 1 
       else
          nshell =  nshell + 1  
          if (nshell == number_G_shells+1) return
          nG_shell_max  = nG_shell_max + 1 
          old_norm2 = norm2 
       end if


    end do

  end subroutine get_G_shells

  subroutine apply_TR_to_wfc (wfc,list_iG)
    !------------------------------------------------------------------
    ! This subroutine takes in a wavefunction wfc = u_{n-k} and 
    ! returns, in the same array, wfc^* = u_{n-k}^*, which is equal to
    ! u_{nk} if time-reversal symmetry applies. 
    !
    ! Actually, -k is in the 1BZ, but k may not be. Define
    !                     k = k_I+G_TR
    !
    !     The wavefunctions are represented as
    !             u_{-k}(r) = \sum_{G} C_{-k}(G) e^{iGr} ===> wfc(iG) = C_{-k}(G)
    !
    !             u_{k }(r) = \sum_{G} C_{k}(G) e^{iGr}   
    !
    !             u_{k_I}(r) = \sum_{G} C_{k_I}(G) e^{iGr}   
    !
    !           ===>  C_{k_I} (G_TR-G) = C^*_{-k}(G)
    !
    !------------------------------------------------------------------
    use intw_fft,  only: find_iG
    !
    implicit none
    !
    ! input/ouput::
    complex(dp)    ::   wfc(nG_max,num_bands,nspin)        
!    complex(dp)    ::   wfc(nG_max,nbands,nspin)
    integer        ::   list_iG(nG_max)     

    ! local::
    complex(dp)    ::   wfc_tmp(nG_max,num_bands,nspin)
!    complex(dp)    ::   wfc_tmp(nG_max,nbands,nspin)
    integer        ::   list_iG_tmp(nG_max)              

    integer        ::   iG, i_minus_G, i
    integer        ::   G(3), minus_G(3)

    integer        ::  nG         ! counter on the number of G vectors in the array
    integer        ::  p_i
    integer        ::  permutations(nG_max)   ! index permutation which orders list_G 


    ! initialize
    nG           =  0
    permutations =  0
    list_iG_tmp  =  0
    wfc_tmp      = cmplx_0

    ! loop on all G
    do i=1,nG_max
       ! find G
       iG = list_iG(i)

       if (iG == 0) exit  ! the index array is zero-padded at the end. 

       nG = nG+1 ! only increment if iG /= 0!

       G  = gvec(:,iG)

       minus_G  = -G

       ! find the index of -G
       call find_iG(minus_G,i_minus_G)
       list_iG_tmp(nG) = i_minus_G

       ! conjugate the wavefunction
       wfc_tmp(nG,:,:)   = conjg(wfc(i,:,:))

    end do

    ! There is no guarantee that the indices in list_iG_k will be sorted in ascending
    ! order! This is not an absolute necessity, but it would be nice and consistent for
    ! the indices to be sorted. 
    ! Sort the indices using a canned heap sort subroutine.

    call HPSORT(nG,list_iG_tmp,permutations)


    ! To understand how this works, consider an example: 
    !
    !            i      f(i)        iG(i)   permutation(i)
    !            ---------------------------------------               
    !            1      0.1         4            2   
    !            2      0.2         1            4
    !            3      0.3         3            3
    !            4      0.4         2            1
    !
    !            j   sort(iG)(j)    sort(f)(j)
    !            ------------------------------------
    !            1      1               0.2  
    !            2      2               0.4
    !            3      3               0.3
    !            4      4               0.1
    !
    !             ===> sort(f) (j)  =   f( permutation(j) )
    ! 

    ! list_iG_tmp is now properly sorted, and can be dumped in the input/output variable
    list_iG = list_iG_tmp

    ! finally, populate the conjugated wave function

    do i= 1, nG
       p_i = permutations(i)
       ! compute the wfc element
       if (nspin==1) then
          wfc(i,:,:) =  wfc_tmp(p_i,:,:)
       elseif (nspin==2) then
          wfc(i,:,1) = - wfc_tmp(p_i,:,2)
          wfc(i,:,2) =   wfc_tmp(p_i,:,1)
       endif
    end do

  end subroutine apply_TR_to_wfc

  subroutine find_size_of_irreducible_k_set(nk_1, nk_2, nk_3, nk_irr)
    !------------------------------------------------------------------
    ! This subroutine finds the number of k-points in the IBZ, based
    ! on the symmetry operations to be used.
    !------------------------------------------------------------------
    implicit none

    !input
    integer  :: nk_1, nk_2, nk_3     ! The input k division

    !output
    integer  :: nk_irr               ! N. of irreducible k points found 
    ! for the nk_1 nk_2 nk_3 division.

    !local variables
    integer :: nkpt                         ! The total number of points 
    integer :: n2n3, n1n2, n1n3             ! useful combinations 
    integer :: i,  j,   k
    integer :: is, js , ks

    integer :: ns

    integer :: ikpt, ikpts

    integer :: aux_sym(3)     

    integer :: i_ID
    integer :: switch
    integer :: il

    logical :: found(nk_1*nk_2*nk_3) 

    nkpt   =  nk_1*nk_2*nk_3

    n1n3   =  nk_1*nk_3
    n2n3   =  nk_2*nk_3
    n1n2   =  nk_1*nk_2

    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    do ns = 1, nsym
       if (s(1,1,ns) == 1 .and. s(1,2,ns) == 0 .and. s(1,3,ns) == 0 .and. &
            s(2,1,ns) == 0 .and. s(2,2,ns) == 1 .and. s(2,3,ns) == 0 .and. &
            s(3,1,ns) == 0 .and. s(3,2,ns) == 0 .and. s(3,3,ns) == 1) then
          i_ID = ns 
       end if
    end do

    switch = 1  ! triplet-to-singlet coordinate

    ! initialize
    found  = .false.
    nk_irr = 0        

    ! loop on the whole mesh
    do i=1,nk_1
       do j=1,nk_2
          do k=1,nk_3

             ! find scalar index of point (i,j,k) 
             call switch_indices(nk_1,nk_2,nk_3,ikpt,i,j,k,switch)
             ! operate on this point only if it has not already been found!
             if (.not. found(ikpt)) then
                ! it's found now. This point is part of the IBZ.
                found(ikpt) = .true.  

                nk_irr = nk_irr + 1        

                ! find the images of this point under all symmetry operations
                do ns=1,nsym    
                   !perform matrix product
                   ! CAREFUL! since the matrix is in crystal coordinates, 
                   ! and it acts in reciprocal space, the convention is :
                   !          k_rot(i) = sum_j s(i,j)*k(j)
                   do il = 1,3

                      aux_sym(il) =     n2n3*s(il,1,ns) * (i-1)  & 
                           +  n1n3*s(il,2,ns) * (j-1)  & 
                           +  n1n2*s(il,3,ns) * (k-1)  
                   end do
                   ! find what point in the 1BZ this corresponds to

                   ! CAREFUL !!! "modulo" is always positive, "mod" can be negative
                   is = modulo ( aux_sym(1), nkpt )/n2n3 + 1
                   js = modulo ( aux_sym(2), nkpt )/n1n3 + 1
                   ks = modulo ( aux_sym(3), nkpt )/n1n2 + 1


                   ! find its singlet coordinate 
                   call switch_indices(nk_1,nk_2,nk_3,ikpts,is,js,ks,switch)

                   ! if its not already accounted for, mark it as found
                   if (.not. found(ikpts)) then
                      found    (ikpts)   = .true. 
                   end if

                   ! Repeat, with Time-Reversal symmetry if present
                   if (TR_symmetry) then

                      aux_sym(:) =  -aux_sym(:)

                      ! CAREFUL !!! "modulo" is always positive, "mod" can be negative
                      is = modulo ( aux_sym(1), nkpt )/n2n3 + 1
                      js = modulo ( aux_sym(2), nkpt )/n1n3 + 1
                      ks = modulo ( aux_sym(3), nkpt )/n1n2 + 1

                      call switch_indices(nk_1,nk_2,nk_3,ikpts,is,js,ks,switch)

                      if (.not. found(ikpts)) then
                         found    (ikpts)   = .true.
                      end if ! found(ikpts)     

                   end if !TR_symmetry

                end do ! ns

             end if ! found(ikpt)     

          end do ! k
       end do ! j
    end do ! i 

  end subroutine find_size_of_irreducible_k_set

  subroutine sort_G_vectors_in_symmetry_shells(nG_shell_max,q_direction,shells_table)
    !----------------------------------------------------------------------------!
    !	This subroutine determines what set of point group operations leave
    !	q_direction invariant. Using this restricted set of symmetries,
    !     the subroutine groups together G vectors belonging to a same 
    !	rotational shell, namely in groups such that all rotations of 
    !	the G vector are present.
    !		
    !     gvec(1:nG_shell_max) represents the first number_G_shells shells. 
    !
    !     It is assumed that gvec is already sorted according to the length
    !     of the G vectors.
    !----------------------------------------------------------------------------!
    use intw_reading,   only : s, nsym, can_use_TR,gvec

    implicit none

    ! input 
    integer   ::  nG_shell_max
    real(dp)  ::  q_direction(3)
    ! output 
    integer   ::  shells_table(nG_shell_max,nsym)

    ! local
    integer   :: iG, iG2
    integer   :: isym

    real(dp)  :: rot_q(3)
    integer   :: restricted_nsym

    integer   :: ishell,iG_shell

    logical   :: G_vector_found(nG_shell_max)

    integer   :: G(3), G2(3), RG(3)
    real(dp)  :: dG(3)

    real(dp)  :: norm

    real(dp)  :: local_rotations(3,3,nsym)
    logical   :: local_TR(nsym)



    ! first, find which symmetries can be used 

    restricted_nsym = 0

    do isym = 1, nsym
       rot_q(:) = matmul(s(:,:,isym),q_direction)

       if (can_use_TR(isym)) rot_q(:) = -rot_q(:)


       norm  = sqrt(dot_product(rot_q-q_direction,rot_q-q_direction))
       if (norm < eps_8) then
          ! add the symmetry to the  set
          restricted_nsym = restricted_nsym +1 

          local_rotations(:,:,restricted_nsym ) = s(:,:,isym)
          local_TR           (restricted_nsym ) = can_use_TR(isym)

       end if

    end do


    ! second, populate the table which contains the G vectors in shells!


    G_vector_found(:) = .false.

    ! This variable will contain the G vector indices, where each row
    ! will represent a rotation shell.
    shells_table(:,:) = 0

    ishell   = 0
    iG_shell = 0

    do iG = 1, nG_shell_max

       ! If the G vector has already been found, skip the rest and move
       ! on to the next index
       if ( G_vector_found(iG) ) cycle

       ! Start a new shell
       ishell   = ishell+1
       iG_shell = 1

       shells_table(ishell,iG_shell) = iG

       G(:) = gvec(:,iG)

       G_vector_found(iG) = .true. ! it is found now!

       ! apply all symmetries to this G vector!
       do isym = 1, restricted_nsym 
          RG(:) = matmul(local_rotations(:,:,isym),G)

          if (local_TR(isym)) RG(:) = -RG(:)

          ! look for this vector in the G vectors not yet found
          do iG2 = 1,nG_shell_max

             if (  G_vector_found(iG2) ) cycle

             G2(:) = gvec(:,iG2)

             dG(:) = ONE*(RG(:) - G2(:))
             norm  = sqrt(dot_product(dG,dG))

             if (norm < eps_8) then
                ! add the new vector to the shell
                G_vector_found(iG2) = .true. ! it is found now!

                iG_shell = iG_shell+1

                shells_table(ishell,iG_shell) = iG2
             end if

          end do ! iG2


       end do ! isym

    end do

  end subroutine sort_G_vectors_in_symmetry_shells

  subroutine output_sorted_G_vectors_in_symmetry_shells(nG_shell_max,q_direction,shells_table)
    !----------------------------------------------------------------------------!
    !	This subroutine writes the G-vector shells to an ASCII file
    !	for testing and checking purposes.
    !----------------------------------------------------------------------------!
    use intw_reading,   only : nsym, gvec, bg

    implicit none

    ! input 
    integer   ::  nG_shell_max
    real(dp)  ::  q_direction(3)
    ! output 
    integer   ::  shells_table(nG_shell_max,nsym)

    ! local
    integer   :: io_unit

    integer   :: ishell, iG_shell, iG1

    integer   :: q(3), G(3)


    real(dp)  ::  q_cart(3), G_cart(3)

    q_cart(:) = matmul(bg,dble(q_direction))

    io_unit = find_free_unit()
    open(unit = io_unit, file ='Symmetry_sorted_G_shells.dat')

    write(io_unit,5) '#=========================================================='
    write(io_unit,5) '#   This file contains the G-vectors sorted in symmetry shells.'
    write(io_unit,5) '#   These shells are obtained by applying all symmetry operations'
    write(io_unit,5) '#   which leaves the following q-vector invariant.  '
    write(io_unit,5) '#                                             '
    write(io_unit,5) '#   The q vector, in crystal coordinates:     '
    write(io_unit,5) '#                                             '
    write(io_unit,10)'#         q(1)  = ', q_direction(1)
    write(io_unit,10)'#         q(2)  = ', q_direction(2)
    write(io_unit,10)'#         q(3)  = ', q_direction(3)
    write(io_unit,5) '#                                             '
    write(io_unit,5) '#   which corresponds to, in cartesian 2pi/a coordinates: '
    write(io_unit,5) '#                                             '
    write(io_unit,10)'#         q_cart(1)  = ', q_cart(1)
    write(io_unit,10)'#         q_cart(2)  = ', q_cart(2)
    write(io_unit,10)'#         q_cart(3)  = ', q_cart(3)
    write(io_unit,5) '#                                             '
    write(io_unit,5) '#=========================================================='
    write(io_unit,5) '# shell number      iG     iG1 iG2 iG3   G vector (2pi/a)'
    write(io_unit,5) '#=========================================================='



    do ishell =1, nG_shell_max

       iG1 = shells_table(ishell,1)

       if (iG1 == 0) exit
       write(io_unit,5) '                                             '

       do iG_shell = 1, nsym
          iG1 = shells_table(ishell,iG_shell)
          if (iG1 == 0) exit 

          G(:)        = gvec(:,iG1)

          G_cart(:) = matmul(bg,dble(G))

          write(io_unit,15) ishell,iG1,G(:), G_cart

       end do

    end do

    close(unit = io_unit)

5   format(A)
10  format(A,3G18.8E3)
15  format(I8,8X,I6,5X,3I4,1x,3F8.2)


  end subroutine output_sorted_G_vectors_in_symmetry_shells

logical function eqvect (x, y, f)
  !-----------------------------------------------------------------------
  !
  implicit none
  real(dp) :: x (3), y (3), f (3)
  ! input: input vector
  ! input: second input vector
  ! input: fractionary translation
  real(dp) :: accep
  ! acceptance parameter
  parameter (accep = 1.0d-4)
  !
  eqvect = abs( x(1)-y(1)-f(1) - nint(x(1)-y(1)-f(1)) ).lt.accep .and. &
           abs( x(2)-y(2)-f(2) - nint(x(2)-y(2)-f(2)) ).lt.accep .and. &
           abs( x(3)-y(3)-f(3) - nint(x(3)-y(3)-f(3)) ).lt.accep
  return
end function eqvect

!-----------------------------------------------------------------------
subroutine multable (nsym, s, table)
  !-----------------------------------------------------------------------
  !
  !  sets up the multiplication table for a group represented by 3x3
  !  integer matrices and checks that {s} is a group indeed:
  !
  !  table(n,m) = index( s(n)*s(m) )
  !
  implicit none
  !
  !    here the dummy variables
  !
  integer :: nsym, s (3, 3, 48), table (48, 48)
  ! input: the number of symmetry of the
  ! input: the symmetry matrices
  ! output: the multiplication table
  !
  !  and here the local variables
  !
  integer :: irot, jrot, krot, ipol, jpol, kpol, ss (3, 3)
  ! \
  !   counter on rotations
  ! /
  ! \
  !   counters on polarizations
  ! /
  ! buffer multiplication matrix

  logical :: found, smn
  ! if true the table has been set
  ! used to check symmetries
  do irot = 1, nsym
     do jrot = 1, nsym
        !
        do ipol = 1, 3
           ! sets up th
           do jpol = 1, 3
              ! product
              ss (ipol, jpol) = 0
              ! matrix
              do kpol = 1, 3
                 !
                 ss (ipol, jpol) = ss (ipol, jpol) + s (ipol, kpol, jrot) * s ( &
                      kpol, jpol, irot)
                 ! ss=s(j)*s(
                 !
              enddo
              !
           enddo
           !
        enddo
        !
        !     here checks that the input matrices really form a group
        !     and sets the multiplication table
        !
        found = .false.
        do krot = 1, nsym
           smn = .true.
           do ipol = 1, 3
              do jpol = 1, 3
                 smn = smn.and. (s (ipol, jpol, krot) .eq.ss (ipol, jpol) )
              enddo
           enddo
           if (smn) then
              if (found) write(*,*)'something wrong in multable:1'
              found = .true.
              table (jrot, irot) = krot
           endif
        enddo

     enddo
     if (.not.found)     write(*,*)'something wrong in multable:2' 

  enddo
  return
end subroutine multable
  !----------------------------------------------------------------------------!
  !
  !
end module intw_symmetries
!
!
!----------------------------------------------------------------------------!
