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

  use kinds, only: dp

  implicit none
  !
  ! variables
  public :: symlink, sym_G, nosym_G, QE_folder_sym, QE_folder_nosym, inverse_indices, &
            identity_matrix_index, spin_symmetry_matrices, full_mesh, IBZ, rtau_index, &
            tau_cryst, rtau, rtau_cryst, symtable
  !
  ! subroutines
  public :: allocate_symmetry_related_k, deallocate_spin_symmetry_matrices, &
            deallocate_symmetry_related_k, find_inverse_symmetry_matrices_indices, &
            allocate_and_build_spin_symmetry_matrices, rotaxis_crystal, &
            set_symmetry_relations, find_the_irreducible_k_set_and_equiv2, &
            find_the_irreducible_k_set_and_equiv, find_the_irreducible_k_set, &
            find_entire_nice_BZ, irr_kp_grid_to_full, calculate_star_r, &
            calculate_star, echo_symmetry_1BZ, rot_atoms, rotate_wfc_test, &
            get_psi_general_k, intw_check_mesh, get_G_shells, apply_TR_to_wfc, &
            find_size_of_irreducible_k_set, find_the_irreducible_k_set_irr, &
            sort_G_vectors_in_symmetry_shells, output_sorted_G_vectors_in_symmetry_shells, &
            multable
  !
  ! functions
  public :: eqvect

  private
  !
  save
  !
  integer, allocatable :: symlink(:,:)
  integer, allocatable :: sym_G(:,:)
  integer, allocatable :: nosym_G(:,:)

  integer, allocatable :: QE_folder_sym(:), QE_folder_nosym(:)

  integer, allocatable :: inverse_indices(:)
  integer :: identity_matrix_index

  ! this variable will contain the 2x2 matrices which rotate
  ! spinors.
  complex(kind=dp), allocatable :: spin_symmetry_matrices(:,:,:)

  ! logical variables defining what is present in the QE folders
  logical :: full_mesh, IBZ

  !Asier && Idoia 17 07 2014
  !Identy of atoms under symmetry operation
  integer, allocatable :: rtau_index(:,:)
  real(kind=dp), allocatable :: tau_cryst(:,:)
  real(kind=dp), allocatable :: rtau(:,:,:)
  real(kind=dp), allocatable :: rtau_cryst(:,:,:)

  integer :: symtable(48,48)

contains
  !
  !-----------------------------------------------------------
  subroutine allocate_symmetry_related_k(nk1,nk2,nk3,nsym)
    !-----------------------------------------------------------
    !
    !------------------------------------------------------------------
    !   nk1, nk2, nk3 are the MP coefficient of the mesh, and nsym is
    !   the total number of point group operations.
    !
    !   Define kpt to be a vector in the 1BZ, and ikpt to be its
    !   singlet index. Let kpt_irr be a symmetry equivalent k-point to kpt.
    !   Finally, let i_folder be an index which represents the QE folders.
    !
    !   This subroutine allocates the arrays
    !
    !   - QE_folder_sym   : this array indicates in which QE folder "i_folder"
    !                       the kpoint "kpt_irr", which is symmetry equivalent to
    !                       "kpt", can be found.
    !
    !   - symlink         : What rotation-like operation must be performed
    !                       on "kpt_irr" to obtain "kpt".
    !
    !   - sym_G           : What G translation must be applied to the rotated
    !                       "kpt_irr" to obtain "kpt".
    !
    !   In equations:
    !           R        = symlink(ikpt)
    !           i_folder = QE_folder_sym(ikpt)
    !           G        = sym_G(ikpt)
    !
    !           =========>  R*kpt_irr = kpt+G
    !
    !   Also, in the case that a full mesh is present in the QE folders
    !   the following arrays will be allocated
    !   - QE_folder_nosym : this array indicates in which QE folder "i_folder"
    !                       the kpoint "kpt_1", which is translation equivalent to
    !                       "kpt", can be found.
    !
    !   - nosym_G         : What G translation must be applied to
    !                       "kpt_1" to obtain "kpt".
    !
    !   In equations:
    !           i_folder = QE_folder_nosym(ikpt)
    !           G        = nosym_G(ikpt)
    !
    !           =========>  kpt_1 = kpt+G
    !------------------------------------------------------------------

    implicit none

    !I/O variables

    integer,intent(in) :: nk1, nk2, nk3, nsym

    !local variables

    integer  :: nkmesh

    nkmesh = nk1*nk2*nk3
    !
    allocate(QE_folder_sym(nkmesh))
    allocate(sym_G(3,nkmesh))
    allocate(symlink(nkmesh,2))
    allocate(QE_folder_nosym(nkmesh))
    allocate(nosym_G(3,nkmesh))
    !
    return

  end subroutine allocate_symmetry_related_k
  !--------------------------------------------------
  !**************************************************
  !--------------------------------------------------
  subroutine deallocate_spin_symmetry_matrices()
    !--------------------------------------------------
    !
    !------------------------------------------------------------------
    ! This subroutine deallocates the array spin_symmetry_matrices
    !------------------------------------------------------------------

    deallocate(spin_symmetry_matrices)

  end subroutine deallocate_spin_symmetry_matrices
  !----------------------------------------------------
  !****************************************************
  !----------------------------------------------------
  subroutine deallocate_symmetry_related_k()
    !----------------------------------------------------
    !
    !------------------------------------------------------------------
    ! This subroutine deallocates the arrays equiv and symlink.
    !------------------------------------------------------------------

    if (allocated(symlink)) deallocate(symlink)
    !
    if (allocated(sym_G)) deallocate(sym_G)
    !
    if (allocated(QE_folder_sym)) deallocate(QE_folder_sym)
    !
    if (allocated(QE_folder_nosym)) deallocate(QE_folder_nosym)
    !
    if (allocated(inverse_indices)) deallocate(inverse_indices)
    !
    if (allocated(nosym_G)) deallocate(nosym_G)
    !
    return

  end subroutine deallocate_symmetry_related_k
  !-------------------------------------------------------
  !*******************************************************
  !-------------------------------------------------------
  subroutine find_inverse_symmetry_matrices_indices()
    !-------------------------------------------------------
    !
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

    use intw_useful_constants, only: ZERO, eps_8, eps_5
    use intw_reading, only: nsym, s, at, bg

    implicit none

    real(kind=dp) :: rot_cart(3,3,nsym), rot1(3,3), rot2(3,3), prod(3,3)
    real(kind=dp) :: norm
    integer :: isym, jsym, ns, i, j, l, mu, nu
    logical :: found

    ! Find the rotation matrices in cartesian coordinates
    !

    allocate(inverse_indices(nsym))

    rot_cart = zero
    !
    do isym=1,nsym
      !
      do mu=1,3
        do nu=1,3
          do i=1,3
            do j=1,3
              !
              rot_cart(mu,nu,isym) = rot_cart(mu,nu,isym) + bg(nu,i)*s(i,j,isym)*at(mu,j)
              !
            enddo !j
          enddo !i
        enddo !nu
      enddo !mu
      !
    enddo !isym
    !
    ! Find the inverse of every symmetry matrix
    !
    inverse_indices = 0
    !
    do isym=1,nsym
      !
      rot1(:,:) = rot_cart(:,:,isym)
      !
      found = .false.
      !
      do jsym=1,nsym
        !
        rot2(:,:) = rot_cart(:,:,jsym)
        !
        ! set the product to minus 1
        !
        prod = zero
        prod(1,1) = -1.0_dp
        prod(2,2) = -1.0_dp
        prod(3,3) = -1.0_dp
        !
        norm = zero
        !
        do i=1,3
          do j=1,3
            do l=1,3
              !
              prod(i,j) = prod(i,j) + rot1(i,l)*rot2(l,j)
              !
            enddo !l
            !
            norm = norm + prod(i,j)**2
            !
          enddo !j
        enddo !i
        !
        norm = sqrt(norm)
        !
        ! IGG -- Sim hexagonalean bestela ez du topatzen
        ! if (norm<eps_8*100 ) then
        if (norm<eps_5) then
          !End IGG
          !
          inverse_indices(isym) = jsym
          found = .true.
          exit
          !
        endif
        !
      enddo !jsym
      !
      if (.not.found) then
        !
        write(*,'(a,I3)') &
              'ERROR: no inverse matrix found for operation isym = ',isym
        stop
        !
      endif
      !
    enddo !isym
    !
    do ns=1,nsym
      !
      if (s(1,1,ns)==1 .and. s(1,2,ns)==0 .and. s(1,3,ns)==0 .and. &
          s(2,1,ns)==0 .and. s(2,2,ns)==1 .and. s(2,3,ns)==0 .and. &
          s(3,1,ns)==0 .and. s(3,2,ns)==0 .and. s(3,3,ns)==1) then
        !
        identity_matrix_index = ns
        !
        exit
        !
      endif
    enddo !ns
    !
    return

  end subroutine find_inverse_symmetry_matrices_indices
  !-------------------------------------------------------------
  !*************************************************************
  !-------------------------------------------------------------
  subroutine allocate_and_build_spin_symmetry_matrices(nsym)
    !-------------------------------------------------------------
    !
    !------------------------------------------------------------------
    !
    ! This subroutine builds, once and for all, all the 2x2 spin
    ! rotation matrices needed by the program and tabulates them in the
    ! array spin_symmetry_matrices.
    !
    ! Some elementary testing is also implemented.
    !
    !------------------------------------------------------------------
    use intw_useful_constants, only: cmplx_0, cmplx_1, cmplx_i, i2, sig_x, sig_y, sig_z, eps_5, pi
    use intw_reading, only: s, at, bg

    implicit none

    !I/O variables

    integer, intent(in) :: nsym

    !local variables

    integer :: isym
    integer :: sym(3,3)
    real(kind=dp) :: sym_cart(3,3) !symmetry matrix in crystal coordinates
    complex(kind=dp) :: S_u(2,2) !the 2x2 spin rotation matrix
    real(kind=dp) :: axis(3), angle, norm !axis and angle of a given rotation matrix
    real(kind=dp) :: axis_cryst(3), det, I3(3,3)
    integer :: is,js
    real(kind=dp) :: bgtrans(3,3)
    integer :: i,j,k,h

       I2(1,1) = cmplx_1 ;    I2(1,2) = cmplx_0 ;    I2(2,1) = cmplx_0  ;    I2(2,2) = cmplx_1
    sig_x(1,1) = cmplx_0 ; sig_x(1,2) = cmplx_1 ; sig_x(2,1) = cmplx_1  ; sig_x(2,2) = cmplx_0
    sig_y(1,1) = cmplx_0 ; sig_y(1,2) =-cmplx_i ; sig_y(2,1) = cmplx_i  ; sig_y(2,2) = cmplx_0
    sig_z(1,1) = cmplx_1 ; sig_z(1,2) = cmplx_0 ; sig_z(2,1) = cmplx_0  ; sig_z(2,2) =-cmplx_1

    I3(1,:) = (/1.0_DP,0.0_DP,0.0_DP/)
    I3(2,:) = (/0.0_DP,1.0_DP,0.0_DP/)
    I3(3,:) = (/0.0_DP,0.0_DP,1.0_DP/)

    !
    allocate(spin_symmetry_matrices(2,2,nsym))
    !
    ! build the spin symmetry matrices for all symmetry operations
    !
    bgtrans = transpose(bg)

    do isym=1,nsym
      !
      sym = s(:,:,isym)
      !
      !ASIER rotaxis MUST be replaced by another procedure
      !diagonalizin the sym_cart matriz
      !       call rotaxis_crystal(sym,axis,angle)
      !call rotaxis(sym,axis,angle)

      det = sym(1,1)*sym(2,2)*sym(3,3) +  &
            sym(1,2)*sym(2,3)*sym(3,1) +  &
            sym(1,3)*sym(2,1)*sym(3,2) -  &
            sym(1,3)*sym(2,2)*sym(3,1) -  &
            sym(1,1)*sym(2,3)*sym(3,2) -  &
            sym(1,2)*sym(2,1)*sym(3,3)

      !the inversion part does not affect the spin
      !we must remove this, otherwise the dimension of the
      !null space of sym-I is more than 1.
      if (det<0) sym = -sym


      sym_cart = 0.0_dp
      do i=1,3
        do j=1,3
          do k=1,3
            do h=1,3
              sym_cart(i,j) = sym_cart(i,j) + at(i,k) * sym(h,k) * bg(j,h)
            enddo
          enddo
        enddo
      enddo


      !Identitity element
      if (sum(abs(sym_cart-I3))<eps_5) then
        axis = (/1.0_dp,0.0_dp,0.0_dp/)
        angle = 0.0
        S_u = I2
        spin_symmetry_matrices(:,:,isym) = S_u(:,:)
      else
        !Other elements
        !
        call compute_rotation_axis(sym_cart,axis,angle)
        !
        S_u(:,:) = cos(angle/2.d0)*I2(:,:)-cmplx_i*sin(angle/2.d0)*(  axis(1)*sig_x(:,:) &
                                                                    + axis(2)*sig_y(:,:) &
                                                                    + axis(3)*sig_z(:,:) )
        spin_symmetry_matrices(:,:,isym) = S_u(:,:)

      end if

      !
    enddo !isym

    !
    return

  end subroutine allocate_and_build_spin_symmetry_matrices

  subroutine compute_rotation_axis(A,axis,angle)
    !ASIER
    !This subrutine is new 15/03/2022
    !determines the null space of (R-I)
    !which is the rotation axis
    use kinds, only : dp
    implicit none

    real(kind=dp), intent(in)  :: a(:,:)
    real(kind=dp), intent(out) :: axis(size(a,1))
    real(kind=dp), intent(out) :: angle

    !local variables
    real(kind=dp) ::  u(size(a,1),size(a,1))
    real(kind=dp) :: vt(size(a,1),size(a,1))
    real(kind=dp) ::  s(size(a,1))
    real(kind=dp) :: ai(size(a,1),size(a,1))

    !3X3 identity matrix
    real(kind=dp) :: i3(size(a,1),size(a,1))

    integer, parameter :: lwmax = 1000
    integer :: info, lwork
    real(kind=dp), allocatable :: work(:)
    integer :: m = 3
    integer :: n = 3
    integer :: ind, k, i
    real(kind=dp) :: u1(3), u2(3), u3(3), kosinu, sinu

    i3(1,:) = (/1.0_DP,0.0_DP,0.0_DP/)
    i3(2,:) = (/0.0_DP,1.0_DP,0.0_DP/)
    i3(3,:) = (/0.0_DP,0.0_DP,1.0_DP/)

    if (n/=m) then
      write(*,*) "rotation matrix must be 3x3"
      write(*,*)m,n
      stop
    end if

    AI = A - I3

    lwork = -1
    allocate(work(1))
    call dgesvd( 'A', 'A', N, N, AI, N, S, U, N, VT, N, &
                WORK, LWORK, INFO )
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    call dgesvd( 'A', 'A', N, N, AI, N, S, U, N, VT, N, &
                WORK, LWORK, INFO )

    !Here we check the dimension of null space of R-I
    k = 0
    do i=1,3
      if (abs(s(i))<1.e-6) then
        ind = i
        k = k + 1
      end if
    end do

    if (k>1) then
      write(*,*) "something wrong in svd part of rotaxis"
      write(*,*) "the dimension of the kernel of r-i must be = 1 ", k
    end if

    ! We construct an (oriented) frame as follows
    ! u1 is the rotation axis
    u1 = vt(ind,:)!/sqrt(sum(vt(ind,:)*vt(ind,:)))

    !Choose any vector of VT, we know it is orthogonal to u1 because VT is
    !an orthogonal matrix
    do i=1,3
      if (i/=ind) then
        u2(:) = vt(i,:)
      end if
    end do
    !3th vector by cross product to ensure orientation.
    !u1xu2
    u3(1) = u1(2) * u2(3) - u1(3) * u2(2)
    u3(2) = u1(3) * u2(1) - u1(1) * u2(3)
    u3(3) = u1(1) * u2(2) - u1(2) * u2(1)

    ! If we apply the rotation to u2: R*u2 =   Cos[angle] * u2 +  Sin[angle] *u3

    ! u2.R*u2 =  Cos[angle]
    ! u3.R*u2 =  Sin[angle]

    kosinu = sum(u2*matmul(A,u2))
    sinu   = sum(u3*matmul(A,u2))

    axis  = u1
    angle = atan2(sinu,kosinu)

  end subroutine compute_rotation_axis

  !*************************************************************
  !-------------------------------------------------------------
  subroutine rotaxis_crystal(sym,axis,angle)
    !-------------------------------------------------------------
    ! ASIER: DO NOT USE THIS 15/03/2022 (BECAUSE WRONG)
    ! USE NEW COMPUTE_ROTATION_AXIS OR OLDER ROTAXIS
    !------------------------------------------------------------------
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
    !------------------------------------------------------------------
    use intw_useful_constants, only: ZERO, ONE, pi, eps_8
    use intw_reading, only: at, bg

    implicit none

    !I/O variables
    integer, intent(in) :: sym(3,3) !symmetry matrix in crystal coordinates
    real(kind=dp), intent(out) :: axis(3), angle !axis and angle of a given rotation matrix

    !local variables

    integer :: determinant, trace, mu, nu, i, j
    real(kind=dp) :: two_sin_angle
    integer :: Rot_cryst(3,3) ! rotation part of the operation, in crystal coordinates
    real(kind=dp) :: Rot_cart(3,3) ! rotation matrix, in cartesian coordinates
    real(kind=dp) :: x_tmp, norm_axis(3), sign_axis(3)
    logical :: vanish(3)

    ! First, find the determinant
    !
    determinant = sym(1,1)*sym(2,2)*sym(3,3) +  &
                  sym(1,2)*sym(2,3)*sym(3,1) +  &
                  sym(1,3)*sym(2,1)*sym(3,2) -  &
                  sym(1,3)*sym(2,2)*sym(3,1) -  &
                  sym(1,1)*sym(2,3)*sym(3,2) -  &
                  sym(1,2)*sym(2,1)*sym(3,3)
    !
    ! Put the rotation part of the symmetry matrix in the Rotation array,
    ! multiplied by the appropriate coefficient to account for inversion
    !
    if (determinant==1) then
      !
      Rot_cryst(:,:) = sym(:,:)
      !
    elseif (determinant==-1) then
      !
      Rot_cryst(:,:) = -sym(:,:)
      !
    else
      !
      write(*,*) '************************************************'
      write(*,*) '** ERROR: The determinant of the rotation     **'
      write(*,*) '**        matrix is not +/- 1.                **'
      write(*,*) '**        review code.                        **'
      write(*,*) '**          program stops.                    **'
      write(*,*) '************************************************'
      !
      stop
      !
    endif
    !
    ! compute the rotation matrix in cartesian coordinates
    !
    Rot_cart(:,:) = ZERO
    !
    do mu=1,3
      do nu=1,3
        do i=1,3
          do j=1,3
            !
            Rot_cart(mu,nu) = Rot_cart(mu,nu) + bg(nu,i)*Rot_cryst(i,j)*at(mu,j)
            !
          enddo !j
        end do !i
      enddo !nu
    enddo !mu
    !
    ! Extract the rotation angle from the trace of the matrix
    !
    trace = Rot_cryst(1,1) + Rot_cryst(2,2) + Rot_cryst(3,3)
    !
    ! there are only 5 possibilities in a crystal;
    ! tabulating insures there is no problem with picking the right quadrant.
    !
    if (trace==-1) then
      !
      angle = pi
      two_sin_angle = ZERO
      !
    elseif (trace==0) then
      !
      angle = 2.0_dp*pi/3.0_dp
      two_sin_angle = sqrt(3.0_dp)
      !
    elseif (trace==1) then
      !
      angle = pi/2.0_dp
      two_sin_angle = 2.0_dp
      !
    elseif (trace==2) then
      !
      angle = pi/3.0_dp
      two_sin_angle = sqrt(3.0_dp)
      !
    elseif (trace==3) then
      !
      angle = ZERO
      two_sin_angle = ZERO
      !
    else
      !
      write(*,*) '************************************************'
      write(*,*) '** ERROR: The trace of the rotation matrix    **'
      write(*,*) '**        is not in [-1,0,1,2,3]. review code.**'
      write(*,*) '**              program stops.                **'
      write(*,*) '************************************************'
      !
      stop
      !
    endif
    !
    ! build the axis array
    !
    if (trace==-1) then
      !
      ! This is the complicated case. Since angle = pi,
      ! the cartesian rotation matrix is symmetric.
      ! A bit of cleverness is required to extract the axis vector.
      !
      ! First, find the norms of the coordinates
      !
      x_tmp = Rot_cart(1,1) + ONE
      !
      if (x_tmp > eps_8) then
        !
        norm_axis(1) = sqrt(x_tmp)/sqrt(2.0_dp)
        vanish(1) = .false.
        !
      else
        !
        norm_axis(1) = ZERO
        vanish(1) = .true.
        !
      endif
      !
      x_tmp = Rot_cart(2,2) + ONE
      !
      if (x_tmp > eps_8) then
        !
        norm_axis(2) = sqrt(x_tmp)/sqrt(2.0_dp)
        vanish(2) = .false.
        !
      else
        !
        norm_axis(2) = ZERO
        vanish(2) = .true.
        !
      endif
      !
      x_tmp = Rot_cart(3,3) + ONE
      !
      if (x_tmp > eps_8) then
        !
        norm_axis(3) = sqrt(x_tmp)/sqrt(2.0_dp)
        vanish(3) = .false.
        !
      else
        !
        norm_axis(3) = ZERO
        vanish(3) = .true.
        !
      endif
      !
      !
      if (.not.vanish(1) .and. .not.vanish(2) .and. .not.vanish(3)) then
        !
        ! if no component vanishes, arbitrarily set the sign of
        ! n3 to be positive.
        !
        sign_axis(3) = ONE
        sign_axis(2) = Rot_cart(2,3)/(2.0_dp*norm_axis(2)*norm_axis(3))
        sign_axis(1) = Rot_cart(1,3)/(2.0_dp*norm_axis(1)*norm_axis(3))
        !
      elseif (.not.vanish(1) .and. .not.vanish(2) .and. vanish(3)) then
        !
        ! if one component vanishes, arbitrarily set the sign of the largest index
        ! component to be positive.
        !
        sign_axis(3) = ZERO
        sign_axis(2) = ONE
        sign_axis(1) = Rot_cart(1,2)/(2.0_dp*norm_axis(1)*norm_axis(2))
        !
      elseif (.not.vanish(1) .and. vanish(2) .and. .not.vanish(3)) then
        !
        sign_axis(3) = ONE
        sign_axis(2) = ZERO
        sign_axis(1) = Rot_cart(1,3)/(2.0_dp*norm_axis(1)*norm_axis(3))
        !
      elseif (vanish(1) .and. .not.vanish(2) .and. .not.vanish(3)) then
        !
        sign_axis(3) = ONE
        sign_axis(1) = ZERO
        sign_axis(2) = Rot_cart(2,3)/(2.0_dp*norm_axis(2)*norm_axis(3))
        !
      elseif (vanish(1) .and. vanish(2) .and. .not.vanish(3)) then
        !
        ! if two components vanish, arbitrarily set the sign of the non
        ! vanishing component to be positive.
        !
        sign_axis(1) = ZERO
        sign_axis(2) = ZERO
        sign_axis(3) = ONE
        !
      elseif (vanish(1) .and. .not.vanish(2) .and. vanish(3)) then
        !
        sign_axis(1) = ZERO
        sign_axis(2) = ONE
        sign_axis(3) = ZERO
        !
      elseif (.not.vanish(1) .and. vanish(2) .and. vanish(3)) then
        !
        sign_axis(1) = ONE
        sign_axis(2) = ZERO
        sign_axis(3) = ZERO
        !
      endif
      !
      axis(1) = norm_axis(1)*sign_axis(1)
      axis(2) = norm_axis(2)*sign_axis(2)
      axis(3) = norm_axis(3)*sign_axis(3)
      !
    elseif (trace==0 .or. trace==1 .or. trace==2) then
      !
      ! For these cases, sin(alpha) is not zero
      ! we can extract the axis from the off diagonal elements
      ! of the rotation matrix.
      !
      axis(1) = ( Rot_cart(3,2)-Rot_cart(2,3) )/two_sin_angle
      axis(2) = ( Rot_cart(1,3)-Rot_cart(3,1) )/two_sin_angle
      axis(3) = ( Rot_cart(2,1)-Rot_cart(1,2) )/two_sin_angle
      !
  elseif (trace==3) then
      !
      ! This is the simplest case: it corresponds to
      ! the identity operation, and the direction of the axis
      ! is irrelevant.
      !
      axis = ZERO
      !
    endif
    !
    return

  end subroutine rotaxis_crystal
  !-------------------------------------------------------------------------------------
  !*************************************************************************************
  !-------------------------------------------------------------------------------------
  subroutine set_symmetry_relations(nk1_,nk2_,nk3_,nkpoints_QE_,kpoints_QE_,kmesh_, k_points_consistent_, &
                                    QE_folder_nosym_, QE_folder_sym_, nosym_G_, sym_G_, &
                                    symlink_, full_mesh_, IBZ_)
    !-------------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------!
    !     This subroutine fills in the symmetry arrays allocated by
    !     the subroutine "allocate_symmetry_related_k", and tests if the
    !     contents of the QE_folders are consistent with the MP mesh.
    !
    !     If the QE folders contain a full MP mesh, this subroutine will
    !     also tabulate the relationship between the k-points in the folders
    !     and the canonical k-points. This will be useful for testing.
    !--------------------------------------------------------------------------
    use intw_useful_constants, only: eps_8
    use intw_utility, only: find_k_1BZ_and_G, switch_indices
    use intw_reading, only: nsym, s, can_use_TR

    implicit none

    !I/O variables

    integer, intent(in) :: nk1_, nk2_, nk3_, nkpoints_QE_  ! Number of k points found in the QE folders
    real(kind=dp), intent(in) :: kpoints_QE_(3,nkpoints_QE_)  ! The kpoints read from the QE folders
    real(kind=dp), intent(in) :: kmesh_(3,nk1_*nk2_*nk3_)  ! The kmesh points, in canonical order
    logical, intent(out) :: k_points_consistent_  ! test if the kpoints are consistent with parameters of input
    integer, intent(out) :: QE_folder_nosym_(nk1_*nk2_*nk3_), QE_folder_sym_(nk1_*nk2_*nk3_)
    integer, intent(out) :: nosym_G_(3,nk1_*nk2_*nk3_), sym_G_(3,nk1_*nk2_*nk3_)
    integer, intent(out) :: symlink_(nk1_*nk2_*nk3_,2)
    logical, intent(out) :: full_mesh_, IBZ_

    !output
    !     This subroutine will fill in the global arrays
    !             - QE_folder_sym(nk1*nk2*nk3)
    !             - sym_G(3,nk1*nk2*nk3)
    !             - symlink(nk1*nk2*nk3,2)
    !
    !     IF a full mesh is potentially present, it will also fill the arrays
    !             - QE_folder_nosym(nk1*nk2*nk3)
    !             - nosym_G(3,nk1*nk2*nk3)

    !local variables

    real(kind=dp) :: kpt(3), rotated_kpt(3), kpt_in_1BZ(3), kpt_on_mesh(3), d_kpt(3)
    logical :: kpoint_is_found_sym(nk1_*nk2_*nk3_), kpoint_is_found_nosym(nk1_*nk2_*nk3_)
    logical :: possible_full_mesh
    integer :: nkmesh, G(3)
    integer :: ikpt, switch, i, j, k, i_folder, isym

    ! initialize arrays to a negative number: if there is a bug
    ! in the code, it will thus look for inexistent folders. It
    ! is better to crash than to produce wrong results!
    !
    QE_folder_sym_(:)   = -4
    QE_folder_nosym_(:) = -4
    !
    k_points_consistent_ = .true.
    !
    nosym_G_ = -4
    sym_G_ = -4
    !
    kpoint_is_found_sym(:)   = .false.
    kpoint_is_found_nosym(:) = .false.
    nkmesh = nk1_*nk2_*nk3_
    !
    if (nkpoints_QE_==nkmesh) then
      !
      possible_full_mesh = .true.
      !
    else
      !
      possible_full_mesh = .false.
      !
    endif
    !
    switch = 1 ! triplet to singlet
    !
    do i_folder=1,nkpoints_QE_
      !
      ! coordinates of this k-point, which need not lie in the 1BZ
      !
      kpt = kpoints_QE_(:,i_folder)
      !
      ! extract the triplet coordinates, the kpt in the 1BZ and
      ! the G vector
      !
      call find_k_1BZ_and_G(kpt,nk1_,nk2_,nk3_,i,j,k,kpt_in_1BZ,G)
      !
      ! test that this triplet index indeed produces the k-point.
      ! This tests that the kpoint is indeed on a mesh consistent
      ! with the input file.
      !
      kpt_on_mesh(1) = dble(i-1)/dble(nk1_)
      kpt_on_mesh(2) = dble(j-1)/dble(nk2_)
      kpt_on_mesh(3) = dble(k-1)/dble(nk3_)
      !
      ! make sure this is identical to the original k-point
      !
      d_kpt = kpt_in_1BZ - kpt_on_mesh
      !
      if (sqrt(dot_product(d_kpt,d_kpt)) > eps_8) then
        !
        k_points_consistent_ = .false.
        !
        write(*,*) ' consistency FAILURE '
        write(*,'(A,3F8.4)') ' kpt        = ',kpt
        write(*,'(A,3F8.4)') ' kpt_in_1BZ = ',kpt_in_1BZ
        write(*,'(A,3F8.4)') ' kpt_on_mesh= ',kpt_on_mesh
        write(*,'(A,3I4)')   ' i, j , k   = ',i,j,k
        !
        exit
        !
      endif
      !
      ! if the k-points are not consistent with the input parameters,
      ! an error will be thrown in the main code. Assume consistency
      ! from this point on.
      !
      ! if there is the possibility of a full mesh being present,
      ! find the correspondence
      !
      if (possible_full_mesh) then
        !
        ! if (sqrt(G(1)**2.0+G(2)**2+G(3)**2) < eps_8) then
        !   !
          call switch_indices(nk1_,nk2_,nk3_,ikpt,i,j,k,switch)
          !
          kpoint_is_found_nosym(ikpt) = .true.
          QE_folder_nosym_(ikpt)       = i_folder
          nosym_G_(:,ikpt)             = G
        ! endif
        !
      endif !possible_full_mesh
      !
      ! loop on all point group operations
      !
      do isym=1,nsym
        !
        ! rotate the k-point.
        !
        rotated_kpt = matmul(dble(s(:,:,isym)),kpt)
        !
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
        !
        ! extract the triplet coordinates, the kpt in the 1BZ and
        ! the G vector
        !
        call find_k_1BZ_and_G(rotated_kpt,nk1_,nk2_,nk3_,i,j,k,kpt_in_1BZ,G)
        !
        ! Tabulate this point as found, but only if allowed to do so!
        !
        ! find its singlet coordinate
        !
        call switch_indices(nk1_,nk2_,nk3_,ikpt,i,j,k,+1)
        !
        ! if this point hasn't been found before, well, it's found now!
        !
        if (.not.kpoint_is_found_sym(ikpt)) then
          !
          kpoint_is_found_sym(ikpt) = .true.
          QE_folder_sym_(ikpt)       = i_folder
          sym_G_(:,ikpt)             = G
          !ASIER 09/03/20222
          !we are taking the inverse of the inverse twice all over the code.
          !symlink(ikpt,1)           = inverse_index(isym)
          symlink_(ikpt,1)           = isym
          symlink_(ikpt,2)           = 0
          !
        endif
        !
      enddo ! isym
      !
      do isym=1,nsym
        !
        rotated_kpt = matmul(dble(s(:,:,isym)),kpt)
        !
        ! repeat with TR symmetry, if allowed (or required!)
        !
        if (can_use_TR(isym)) then
          !
          rotated_kpt = -rotated_kpt
          !
          call find_k_1BZ_and_G(rotated_kpt,nk1_,nk2_,nk3_,i,j,k,kpt_in_1BZ,G)
          !
          ! find its singlet coordinate
          !
          call switch_indices(nk1_,nk2_,nk3_,ikpt,i,j,k,switch)
          !
          if (.not.kpoint_is_found_sym(ikpt)) then
            !
            kpoint_is_found_sym(ikpt) = .true.
            QE_folder_sym_(ikpt)       = i_folder
            sym_G_(:,ikpt)             = -G
            !ASIER, same as above
            symlink_(ikpt,1)           = isym
            !symlink(ikpt,1)           = inverse_indices(isym)
            symlink_(ikpt,2)           = 1
            !
            ! CAREFUL! It is now -G which enters the sym_G array
            ! See the rotation code for details, as well as the
            ! symmetries document.
            !
          endif ! not found
          !
        endif ! TR
        !
      enddo !isym
      !
    enddo !i_folder
    !
    ! test if all the kmesh points were found, and store the result
    ! in these GLOBAL logical variables
    !
    IBZ_       = all(kpoint_is_found_sym)
    full_mesh_ = all(kpoint_is_found_nosym)
    !
    return

  end subroutine set_symmetry_relations
  !-----------------------------------------------

  subroutine find_the_irreducible_k_set_and_equiv2(nspt, k_entire, k_irr, nk_irr, equiv, symlink, sym_G)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the a general k list
    !------------------------------------------------------------------
    use intw_input_parameters, only: TR_symmetry
    use intw_reading, only: nsym, s, at, bg

    implicit none

    !input

    integer :: nspt

    real(kind=dp) :: k_entire(3,nspt)

    ! output
    integer :: nk_irr
    real(kind=dp) :: k_irr(3,nspt)
    integer :: equiv(nspt) ! which is the equivalent point
    integer :: sym_G(1:3,nspt), symlink(nspt,1:2)

    !local variables
    real(kind=dp) :: k_rot(3), dist1, dist2

    integer :: i, j ! triplet indices

    integer :: ns

    integer :: ikpt ! singlet index, singlet index obtained by symmetry

    integer :: switch, ii, jj, kk

    logical :: found(nspt)
    real(kind=dp), parameter :: eps = 10E-7

    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    switch = 1 ! triplet-to-singlet coordinate

    found  = .false.
    nk_irr = 0

    do i=1,nspt
      if (.not. found(i)) then
        found(i) = .true.
        nk_irr = nk_irr + 1
        k_irr(1:3,nk_irr) = k_entire(1:3,i)
        equiv(i) = nk_irr

        do ns=1,nsym

          k_rot = matmul(dble(s(:,:,ns)), k_irr(:,nk_irr))

          do j=1,nspt

            if (((abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)-1.0_dp)<eps) .or. &
                 (abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)       )<eps)) .and.&
                ((abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)-1.0_dp)<eps) .or. &
                 (abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)       )<eps)) .and.&
                ((abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)-1.0_dp)<eps) .or. &
                 (abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)       )<eps))) then

              equiv(j) = nk_irr
              found(j) = .true.
              symlink(j,1) = ns
              symlink(j,2) = 0
              sym_G(:,j) = nint(k_entire(:,j) - k_rot(:))
              cycle
            endif

            if (TR_symmetry) then

              if (((abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)-1.0_dp)<eps) .or. &
                   (abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)       )<eps)) .and.&
                  ((abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)-1.0_dp)<eps) .or. &
                   (abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)       )<eps)) .and.&
                  ((abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)-1.0_dp)<eps) .or. &
                   (abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)       )<eps))) then
                equiv(j) = nk_irr
                found(j) = .true.
                symlink(j,1) = ns
                symlink(j,2) = 1
                sym_G(:,j) = nint(k_entire(:,j) -(- k_rot(:)))

                endif
            endif

          end do !j

        enddo ! ns

      end if !found

    enddo !i


    do ikpt=1,nk_irr

      k_rot(:) = matmul(bg,k_irr(:,ikpt))

      dist1 = sum(k_rot(:)**2)

      do ii=-1,1
        do jj=-1,1
          do kk=-1,1

            k_rot(:) = matmul(bg,k_irr(:,ikpt)) + matmul(bg, dble((/ii,jj,kk/)))
            dist2 = sum(k_rot(:)**2)

            if (dist2<dist1) then
              k_irr(:,ikpt) = matmul(transpose(at), k_rot)
              dist1 = dist2
            endif

          enddo
        enddo
      enddo
    enddo

  end subroutine find_the_irreducible_k_set_and_equiv2

  subroutine find_the_irreducible_k_set_and_equiv(nspt, k_entire, k_irr, nk_irr, equiv)!, symlink, sym_G)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the a general k list
    !------------------------------------------------------------------
    use intw_input_parameters, only: TR_symmetry
    use intw_reading, only: nsym, s, at, bg

    implicit none

    !input

    integer  :: nspt

    real(kind=dp) :: k_entire(3,nspt)

    ! output
    integer  :: nk_irr
    real(kind=dp) :: k_irr(3,nspt)
    integer :: equiv(nspt) ! which is the equivalent point
    integer :: sym_G(1:3,nspt), symlink(nspt,1:2)

    !local variables
    real(kind=dp):: k_rot(3), dist1, dist2

    integer :: i, j ! triplet indices

    integer :: ns

    integer :: ikpt ! singlet index, singlet index obtained by symmetry

    integer :: switch, ii, jj, kk

    logical :: found(nspt)
    real(kind=dp), parameter :: eps = 10E-7

    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    switch = 1 ! triplet-to-singlet coordinate

    found  = .false.
    nk_irr = 0

    do i=1,nspt
      if (.not. found(i)) then
        found(i) = .true.
        nk_irr = nk_irr + 1
        k_irr(1:3,nk_irr) = k_entire(1:3,i)
        equiv(i) = nk_irr

        do ns=1,nsym
          k_rot = matmul(dble(s(:,:,ns)), k_irr(:,nk_irr))

          do j=1, nspt

            if (((abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)-1.0_dp)<eps) .or. &
                 (abs(modulo(k_entire(1,j)-k_rot(1),1.0_dp)       )<eps)) .and.&
                ((abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)-1.0_dp)<eps) .or. &
                 (abs(modulo(k_entire(2,j)-k_rot(2),1.0_dp)       )<eps)) .and.&
                ((abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)-1.0_dp)<eps) .or. &
                 (abs(modulo(k_entire(3,j)-k_rot(3),1.0_dp)       )<eps))) then

                equiv(j) = nk_irr
                found(j) = .true.
                symlink(j,1) = ns
                symlink(j,2) = 0
                sym_G(:,j) = nint(k_entire(:,j) - k_rot(:))
                cycle
            endif

            if (TR_symmetry) then

              if (((abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)-1.0_dp)<eps) .or. &
                   (abs(modulo(k_entire(1,j)+k_rot(1),1.0_dp)       )<eps)) .and.&
                  ((abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)-1.0_dp)<eps) .or. &
                   (abs(modulo(k_entire(2,j)+k_rot(2),1.0_dp)       )<eps)) .and.&
                  ((abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)-1.0_dp)<eps) .or. &
                   (abs(modulo(k_entire(3,j)+k_rot(3),1.0_dp)       )<eps))) then
                equiv(j) = nk_irr
                found(j) = .true.
                symlink(j,1) = ns
                symlink(j,2) = 1
                sym_G(:,j) = nint(k_entire(:,j) -(- k_rot(:)))

              endif
            endif

          end do !j

        enddo ! ns

      end if !found

    enddo !i


    do ikpt=1,nk_irr

      k_rot(:) = matmul(bg,k_irr(:,ikpt))

      dist1 = sum(k_rot(:)**2)

      do ii=-1,1
        do jj=-1,1
          do kk=-1,1

            k_rot(:) = matmul(bg,k_irr(:,ikpt)) + matmul(bg, dble((/ii,jj,kk/)))
            dist2 = sum(k_rot(:)**2)

            if (dist2<dist1) then
                k_irr(:,ikpt) = matmul(transpose(at), k_rot)
                dist1 = dist2
            endif

          enddo
        enddo
      enddo
    enddo

  end subroutine find_the_irreducible_k_set_and_equiv

  subroutine find_the_irreducible_k_set(nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the canonical
    ! 1BZ mesh.
    !
    !------------------------------------------------------------------
    use intw_input_parameters, only: TR_symmetry
    use intw_useful_constants, only: eps_8
    use intw_utility, only: find_k_1BZ_and_G, switch_indices
    use intw_reading, only: nsym, s

    implicit none

    !input
    integer, intent(in) :: nk1, nk2, nk3
    ! The input k division
    real(kind=dp), intent(in) :: kmesh(3,nk1*nk2*nk3)
    ! The full kmesh, in canonical order

    !output
    real(kind=dp), intent(out) :: kpoints_irr(3,nk1*nk2*nk3)
    ! The irreducible kpoints in crystal triplet coordinates.
    ! The size of the array is nk1* nk2* nk3 instead of nk_irr;
    ! it is supposed that we still do not know the value of nk_irr

    integer :: nk_irr
    ! How many points were found?

    !local variables
    real(kind=dp) :: k_rot(3), k_1BZ(3), dk(3)
    integer :: nkpt ! The total number of points
    integer :: i, j, k ! triplet indices
    integer :: is, js, ks ! triplet indices  obtained by symmetry
    integer :: G(3)
    integer :: isym
    integer :: ikpt, ikpts ! singlet index, singlet index obtained by symmetry
    logical :: found(nk1*nk2*nk3)

    nkpt = nk1*nk2*nk3
    found = .false.
    nk_irr = 0
    !
    ! loop on the whole mesh, in the appropriate order (very important!)
    do i=1,nk1
      do j=1,nk2
        do k=1,nk3
          !
          ! find scalar index of point (i,j,k)
          call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,+1)
          !
          ! operate on this point only if it has not already been found!
          if (.not. found(ikpt)) then
            !
            ! it's found now. This point is part of the IBZ.
            found(ikpt) = .true.
            !
            nk_irr = nk_irr + 1
            !
            kpoints_irr(1,nk_irr) = dble(i-1)/nk1
            kpoints_irr(2,nk_irr) = dble(j-1)/nk2
            kpoints_irr(3,nk_irr) = dble(k-1)/nk3
            !
            ! loop on all symmetry operations
            do isym=1,nsym
              !
              !perform matrix product
              ! CAREFUL! since the matrix is in crystal coordinates,
              ! and it acts in reciprocal space, the convention is :
              !          k_rot(i) = sum_j s(i,j)*k(j)
              !
              k_rot = matmul(dble(s(:,:,isym)), kpoints_irr(:,nk_irr))
              !
              ! find what point in the 1BZ this corresponds to
              call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)
              !
              ! check that k_1BZ+G = k_rot. If not, k_rot isn't on the mesh,
              ! and the algorithm in "find_k_1BZ_and_G" cannot be trusted.
              dk = k_rot - (k_1BZ + dble(G))
              if (sqrt(dot_product(dk,dk))<eps_8) then
                !
                ! what is the scalar index
                call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,+1)
                !
                if (.not. found(ikpts)) found(ikpts) = .true.
                !
              endif ! dk
              !
              ! Repeat, with Time-Reversal symmetry if present
              if (TR_symmetry) then
                !
                k_rot = -k_rot
                !
                ! find what point in the 1BZ this corresponds to
                call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)
                !
                ! we check again the value of dk, so if k_1BZ+G = k_rot
                dk = k_rot - (k_1BZ + dble(G))
                if (sqrt(dot_product(dk,dk))<eps_8) then
                  !
                  ! what is the scalar index
                  call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,+1)
                  !
                  if (.not.found(ikpts)) found(ikpts) = .true.
                  !
                endif ! dk
                !
              endif !TR_symmetry
              !
            enddo ! isym
            !
          endif ! found(ikpt)
          !
        enddo ! k
      enddo ! j
    enddo ! i

  end subroutine find_the_irreducible_k_set

  subroutine find_entire_nice_BZ(nk1, nk2, nk3, nspt, ksvec) !(kmesh, kpoints_irr,nk_irr)

    !------------------------------------------------------------------
    ! input: private nk1 nk2 nk3, output full BZ nk>nk1*nk2*nk3 for tria_diag by AE.&IGdG.
    !------------------------------------------------------------------
    use intw_input_parameters, only: TR_symmetry
    use intw_useful_constants, only: eps_8
    use intw_utility, only: find_k_1BZ_and_G, switch_indices
    use intw_reading, only: nsym, s, at, bg

    implicit none

    !input

    integer, intent(in) :: nk1, nk2, nk3

    ! output
    real(kind=dp), intent(out) :: ksvec(3,2*nk1*nk2*nk3)
    integer, intent(out) :: nspt

    !local variables
    real(kind=dp) :: k_rot(3), k_1BZ(3), dk(3)

    integer :: nkpt ! The total number of points
    integer :: i, j, k, ii, jj, kk ! triplet indices
    integer :: is, js, ks ! triplet indices  obtained by symmetry

    integer :: G(3)

    integer :: ns

    integer :: ikpt, ikpts ! singlet index, singlet index obtained by symmetry

    integer :: switch

    logical :: found(nk1*nk2*nk3)

    real(kind=dp) :: k_aux(3,48)

    real(kind=dp),parameter  :: eps = 1d-6

    integer :: nsp, iss, nk_irr
    real(kind=dp) :: kpoints_irr(3,nk1*nk2*nk3), dist1, dist2


    nkpt = nk1*nk2*nk3


    ! Find which symmetry operation is the identity
    ! most likely always the first element, but let's be sure

    switch = 1 ! triplet-to-singlet coordinate

    found = .false.
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
              dk = k_rot - (k_1BZ + dble(G))

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

                dk = k_rot - (k_1BZ + dble(G))

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

    do ikpt=1,nk_irr

      k_rot(:) = matmul(bg,kpoints_irr(:,ikpt))

      dist1 = sum(k_rot(:)**2)

      do ii=-1,1
        do jj=-1,1
          do kk=-1,1

            k_rot(:) = matmul(bg,kpoints_irr(:,ikpt)) + matmul(bg, dble((/ii,jj,kk/)))
            dist2 = sum(k_rot(:)**2)

            if (dist2<dist1) then
                kpoints_irr(:,ikpt) = matmul(transpose(at), k_rot)
                dist1 = dist2
            endif

          enddo
        enddo
      enddo
    enddo


    nspt = 0
    do ikpt=1,nk_irr
      nsp = 1
      nspt = nspt + 1
      k_aux(:,nsp) = kpoints_irr(:,ikpt)
      ksvec(:,nspt) = k_aux(:,nsp)
      s_l: do is=1,nsym
        k_rot = matmul(dble(s(:,:,is)), kpoints_irr(:,ikpt))
        iss_l: do iss=1,nsp
            if (sum(abs(k_aux(:,iss)-k_rot(:)))<eps) cycle s_l
        enddo iss_l
        nsp = nsp + 1
        k_aux(:,nsp) = k_rot(:)
        nspt = nspt + 1
        ksvec(:,nspt) = k_aux(:,nsp)
      enddo s_l
    enddo

    do ikpt=1, nspt
      write(123,"(100f12.6)") ksvec(:,ikpt)
    enddo

  end subroutine find_entire_nice_BZ


  subroutine irr_kp_grid_to_full(nk_1, nk_2, nk_3, nk_irr, k_irr, &
                                 equiv_l, symlink_l)
    !-------------------------------------------------------------------------
    ! Given a irreducible k set, this subroutine finds the connection
    ! to the full BZ.
    !
    ! input  : nk1 nk2 nk3 nk_irr (n. of irr k points) and k_irr(3,nk_irr)
    ! output : equiv_l, symlink_l
    !
    !-------------------------------------------------------------------------
    use intw_input_parameters, only: TR_symmetry
    use intw_utility, only: switch_indices
    use intw_reading, only: nsym, s

    implicit none

    !input
    integer :: nk_1, nk_2, nk_3

    integer :: nk_irr                  ! N. of irreducible k points

    integer :: k_irr(3,nk_irr)         ! triplet indices of the irreducible k-points

    !output
    integer :: equiv_l(nk_1* nk_2* nk_3)   ! For a given k point in the entire BZ,
    ! gives the symmetry related
    ! k point index in the irreducible k set

    integer :: symlink_l(nk_1* nk_2* nk_3,2) ! Symmetry link .


    !local
    integer :: nkpt, ikpts
    integer :: n2n3, n1n2, n1n3

    integer :: il

    integer :: irr

    integer :: i, j, k
    integer :: is, js, ks

    integer :: aux_sym(3)

    integer :: ns
    integer :: switch

    logical :: done(nk_1*nk_2*nk_3)


    switch = 1 ! tripet - to - singlet
    nkpt = nk_1*nk_2*nk_3

    n1n3 = nk_1*nk_3
    n2n3 = nk_2*nk_3
    n1n2 = nk_1*nk_2

    done(:) = .false.


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

    do irr=1,nkpt

      i = k_irr(1,irr)
      j = k_irr(2,irr)
      k = k_irr(3,irr)

      do ns=1,nsym

        do il=1,3
          aux_sym(il) = n2n3*s(il,1,ns) * (i-1) + &
                        n1n3*s(il,2,ns) * (j-1) + &
                        n1n2*s(il,3,ns) * (k-1)
        enddo
        ! CAREFUL !!! "modulo(a,p)" is always positive for p positive;
        !              "mod(a,p)" can be negative
        is = modulo( aux_sym(1), nkpt )/n2n3 + 1
        js = modulo( aux_sym(2), nkpt )/n1n3 + 1
        ks = modulo( aux_sym(3), nkpt )/n1n2 + 1

        call switch_indices(nk_1,nk_2,nk_3,ikpts,is,js,ks,switch)

        if (.not.done(ikpts)) then
          done(ikpts) = .true.
          equiv_l(ikpts) = irr
          !         symlink_l (ikpts,1) = ns
          symlink_l(ikpts,1)  =  inverse_indices(ns)
          symlink_l(ikpts,2) = 0
        endif

        ! If time reversal symmetry is present, use it!
        if (TR_symmetry) then
          is = -is
          js = -js
          ks = -ks
          call switch_indices(nk_1,nk_2,nk_3,ikpts,is,js,ks,switch)

          if (.not.done(ikpts)) then
            done(ikpts) = .true.
            equiv_l(ikpts) = irr
            !             symlink_l (ikpts,1) = ns
            symlink_l(ikpts,1) =  inverse_indices(ns)
            symlink_l(ikpts,2) = 1
          endif

        end if

      enddo ! ns
    enddo ! irr


  end subroutine irr_kp_grid_to_full

  subroutine calculate_star_r(v, vstar, nstar, symop)

    use intw_reading, only: nsym, s

    implicit none

    real(kind=dp), intent(in) :: v(3)
    real(kind=dp), intent(out) :: vstar(3,48)
    integer, intent(out) :: nstar, symop(48)
    integer :: isym, i
    real(kind=dp) :: vrot(3)


    nstar = 1
    vstar(1:3,nstar) = v(1:3)
    symop(1) = 1

    do isym=1,nsym
      vrot(:) = matmul(dble(s(:,:,isym)), v(:))

      do i=1,nstar
        if ( sum(abs(vrot(:)- vstar(1:3,i)))<10E-5 ) then
          goto 1987
        end if
      enddo

      nstar = nstar + 1
      vstar(1:3,nstar) = vrot(1:3)
      symop(nstar)=isym

1987  continue
    enddo


  end subroutine calculate_star_r

  subroutine calculate_star(v, vstar, nstar, symop)

    use intw_reading, only: nsym, s

    implicit none

    integer, intent(in) :: v(3)
    integer, intent(out) :: vstar(3,48), symop(48)
    integer, intent(out) :: nstar
    integer :: isym, i
    real(kind=dp) :: vrot(3)


    nstar = 1
    vstar(1:3,nstar) = v(1:3)
    symop(1) = 1

    do isym=1,nsym
      vrot(:) = nint(matmul(dble(s(:,:,isym)), v(:)))

      do i=1,nstar
        if ( sum(abs(vrot(:)-vstar(1:3,i)))<10E-5 ) then
          goto 1984
        end if
      enddo

      nstar = nstar + 1
      vstar(1:3,nstar) = vrot(1:3)
      symop(nstar) = isym

1984  continue
    enddo


  end subroutine calculate_star

  subroutine echo_symmetry_1BZ(nk_1,nk_2,nk_3,nk_irr,equiv,symlink)
    !--------------------------------------------------------------------------------
    !
    ! simple writing routine, for testing purposes
    !
    !--------------------------------------------------------------------------------
    use intw_utility, only: find_free_unit

    implicit none

    integer :: nk_1, nk_2, nk_3, nkr, nk_irr
    integer :: io_unit
    integer :: i

    integer :: equiv(nk_1*nk_2*nk_3)
    integer :: symlink(nk_1*nk_2*nk_3,2)

    nkr = nk_1*nk_2*nk_3

    io_unit = find_free_unit()
    open(io_unit,file='equivalence_of_k_points.test')
    write(io_unit,'(a,i4)')'          nk_irr  = ', nk_irr
    write(io_unit,'(a)')'----------------------------'

    write(io_unit,'(a)')'          Nk              Equivalent      Symmetry           Time reversal. '
    write(io_unit,'(a)')'                          in reduced l.   connection.        not used:0, used :1  '
    write(io_unit,'(a)')'-----------------------------------------------------------------------------------------------'

    do i=1,nkr
      write(io_unit,100) i, equiv(i), symlink(i,:)
    enddo

    close(io_unit)

100 format(i10,12x,i10,12x,i10,12x,i10)
  end subroutine echo_symmetry_1BZ

  subroutine rot_atoms(nat,nsym,tau)
    use intw_utility, only: ainv
    use intw_reading, only: nr1, nr2, nr3
    use intw_reading, only: s, ftau, at, bg

    implicit none
    integer, intent(in) :: nat, nsym
    real(kind=dp), intent(in) :: tau(3,nat)
    !out global
    !    integer, intent(out) :: rtau_index(nat,nsym)
    !    real(kind=dp), intent(out) :: rtau(3,nsym,nat)
    !    real(kind=dp), intent(out) :: tau_cryst(3,nat)
    !    real(kind=dp), intent(out) :: rtau_cryst(3,nsym,nat)

    integer :: i, j, h
    integer :: isym
    integer :: nr(3)
    integer :: a_index, na

    integer :: ipol, jpol, kpol, lpol

    real(kind=dp), parameter :: epsat = 1E-3

    real(kind=dp) :: s_cart(3,3,nsym)

    do isym=1,nsym
      do ipol=1,3
        do jpol=1,3
          s_cart(ipol,jpol,isym) = 0.d0
          do kpol=1,3
            do lpol=1,3
              s_cart(ipol,jpol,isym) = s_cart(ipol,jpol,isym) + at(ipol,kpol) * &
                                       s(lpol,kpol,isym) * bg(jpol,lpol)
            enddo
          enddo
        enddo
      enddo
    enddo

    nr = (/nr1,nr2,nr3/)

    rtau_index = -11

    do i=1,nat
      tau_cryst(:,i) = matmul(ainv(at), tau(:,i))
    enddo

    do i=1,nat
      do isym=1,nsym
        do ipol=1,3

          rtau_cryst(ipol,isym,i) = s(1,ipol,isym) * tau_cryst(1,i) + &
                                    s(2,ipol,isym) * tau_cryst(2,i) + &
                                    s(3,ipol,isym) * tau_cryst(3,i)

          rtau_cryst(ipol,isym,i) = rtau_cryst (ipol,isym,i) - dble(ftau(ipol,isym))

        end do

       enddo
    enddo

    do i=1,nat
      do isym=1,nsym
        do j=1,nat
          if (eqvect(rtau_cryst(:,isym,i), tau_cryst(:,j), (/0.d0,0.0d0,0.0d0/))) rtau_index(i,isym) = j
        enddo
      enddo
    enddo

    do isym=1,nsym
      do na=1,nat
        a_index = rtau_index(na,isym)
        do h=1,3
          rtau(h,isym,na) = s_cart(1,h, isym) * tau(1,na) + &
                            s_cart(2,h, isym) * tau(2,na) + &
                            s_cart(3,h, isym) * tau(3,na)

          rtau(h,isym,na) = rtau(h,isym,na) - tau(h,a_index)
        end do
       enddo
    enddo

    do i=1,nat
      do isym=1,nsym
        if (rtau_index(i,isym).eq.0) then
          write(*,*)'ERROR in rot_at: At least one atom does not map properly under sym. op.', isym, 'atom:', i
        endif
      enddo
    enddo

  end subroutine rot_atoms
  !-----------------------------------------------------------------------------
  !*****************************************************************************
  !-----------------------------------------------------------------------------
  subroutine rotate_wfc_test(wfc_k_irr,list_iG,wfc_k,list_iG_k, &
                             i_sym,sym,ftau,G_sym)
    !-----------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------------------------------------
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
    !			-  sym is a point group operation:
    !                        it is the INVERSE of the actual operation; this is the
    !                        appropriate operator to act on k, in crystal coordinates
    !			-  k_irr is a k-point in the IBZ
    !			-  sym * k_irr  = k + G_sym, with k in the 1BZ
    !
    !	applying the point group operation yields
    !
    !           u_{nk}(sym_l*G+G_sym) =  e^{i R*G*tau} u_{nk_irr}(G)
    !--------------------------------------------------------------------------------------------------------
    use intw_fft, only: find_iG
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i
    use intw_utility, only: HPSORT, cmplx_ainv_2
    use intw_reading, only: nG_max, gvec, nspin, num_bands_intw

    implicit none

    !I/O variables

    integer, intent(in) :: i_sym ! index of the symmetry operation
    integer, intent(in) :: G_sym(3) ! G vector such that  R*k + G_sym = sym_l * k_irr
    integer, intent(in) :: sym(3,3) ! inverse point group operation the one acting on k (cryst.coord.)
    integer, intent(in) :: list_iG(nG_max) ! G vector indices for k_irr
    real(kind=dp), intent(in) :: ftau(3) ! fractional translation associated with point group operation
    complex(kind=dp), intent(inout) :: wfc_k_irr(nG_max,num_bands_intw,nspin) ! wfc at point k_irr in the IBZ
    integer, intent(out) :: list_iG_k(nG_max) ! G vector indices for k, sorted
    complex(kind=dp), intent(out) :: wfc_k(nG_max,num_bands_intw,nspin) ! rotated wfc at point k in the 1BZ

    !local variables

    complex(kind=dp) :: wfc_k_aux(nG_max,num_bands_intw,nspin)
    integer :: p_i, i, alpha, iGk, iG_k
    integer :: Gk(3) ! a vector for k in the IBZ
    integer :: RGk(3) ! ( symmetry operation )* G_k
    integer :: G_k(3) ! a vector for Rk, the point in the 1BZ
    integer :: permutations(nG_max) ! index permutation which orders list_G_k
    integer :: ibnd, is, js
    integer :: nG  ! counter on the number of G vectors in the array
    complex(kind=dp) :: phases(nG_max)

    phases(:) = cmplx_0
    wfc_k(:,:,:) = cmplx_0
    !
    list_iG_k(:) = 0
    nG = 0
    permutations(:) = 0
    !
    ! loop on all Gk, the coefficients of the wavefunction at the IBZ k point
    !
    do i=1,nG_max
      !
      iGk = list_iG(i)
      !
      if (iGk==0) exit  ! the index array is zero-padded at the end.
      nG = nG + 1 ! only increment if iGk /= 0!
      !
      Gk = gvec(:,iGk)
      !
      RGk = matmul(sym,Gk)
      !
      !ASIER: Gk and NOT RGk!!!
      ! - sign is well checked below.
      phases(nG) = exp(- cmplx_I*tpi*sum(dble(Gk) * ftau))

      G_k(:) = RGk(:) + G_sym(:)
      call find_iG(G_k,iG_k)
      !
      list_iG_k(nG) = iG_k
      !
    enddo
    !
    call HPSORT(nG,list_iG_k,permutations)
    !
    do i=1,nG
      !
      ! compute the wfc element
      !
      p_i = permutations(i)
      !
      do ibnd=1,num_bands_intw
        do is=1,nspin
          !
          wfc_k(i,ibnd,is) = wfc_k_irr(p_i,ibnd,is) * phases(p_i)
          !
        enddo !is
      enddo !ibnd
      !
    enddo !i

    !
    if (nspin==2) then
      !
      !print*, "spin matrix--------------------", i_sym
      !do is=1,nspin
      !  write(*,"(100f6.2)") (spin_symmetry_matrices(is,js,i_sym), js=1,nspin)
      !end do
      !
      wfc_k_aux = wfc_k
      wfc_k     = cmplx_0
      !
      do i=1,nG
        do ibnd=1,num_bands_intw
          do is=1,nspin
            do js=1,nspin
              !
              ! JLB, MBR 29/06/2023
              !wfc_k(i,ibnd,is) = wfc_k(i,ibnd,is) + spin_symmetry_matrices(is,js,i_sym) * wfc_k_aux(i,ibnd,js)
              wfc_k(i,ibnd,is) = wfc_k(i,ibnd,is) + spin_symmetry_matrices(is,js,inverse_indices(i_sym)) * wfc_k_aux(i,ibnd,js)
              !
            enddo !js
          enddo !is
        enddo !ibnd
      enddo !i
      !
    endif ! If non-collinear
    !
    return

  end subroutine rotate_wfc_test
  !---------------------------------------------------------------------
  !*********************************************************************
  !---------------------------------------------------------------------
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
    use intw_input_parameters, only: nk1, nk2, nk3
    use intw_reading, only: get_K_folder_data
    use intw_reading, only: s, ftau, nG_max, npol, nspin, num_bands_intw
    use intw_useful_constants, only: ZERO
    use intw_fft, only: wfc_by_expigr
    use intw_utility, only: find_k_1BZ_and_G, switch_indices

    implicit none

    real(kind=dp), intent(in) :: kpoint(3)
    logical, intent(in) :: use_IBZ


    complex(kind=dp), intent(out) :: wfc_k(nG_max,num_bands_intw,nspin)
    ! complex(kind=dp), intent(out) :: wfc_k(nG_max,nbands,nspin)
    integer, intent(out) :: list_iG(nG_max)
    real(kind=dp), intent(out) :: QE_eig(num_bands_intw)
    ! real(kind=dp), intent(out) :: QE_eig(nbands)

    complex(kind=dp) :: wfc_k_irr(nG_max,num_bands_intw,nspin)
    ! complex(kind=dp) :: wfc_k_irr(nG_max,nbands,nspin)
    integer :: list_iG_irr(nG_max)


    integer :: ikpt, i_folder

    integer :: i_sym, TR
    integer :: G_sym(3), G_plus(3)
    integer :: sym(3,3)


    real(kind=dp) :: ftau_sym(3)
    integer :: i_1bz, j_1bz, k_1bz

    real(kind=dp) :: kpoint_1BZ(3)


    !!Peio
    !!The variable we use instead of num_bands
    !integer :: nbands_loc
    !!Peio
    !nbands_loc = num_bands_intw

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

    call switch_indices(nk1,nk2,nk3,ikpt,i_1bz,j_1bz,k_1bz, +1)

    if (.not.use_IBZ) then
      ! do not use the IBZ; just fetch the wavefunction from the QE folders
      i_folder = QE_folder_nosym(ikpt)
      G_sym    = nosym_G(:,ikpt) + G_plus(:) !Asier&&Idoia 24 06 2014
      ftau_sym = ZERO
      sym      = s(:,:,identity_matrix_index)
      call get_K_folder_data(i_folder,list_iG_irr,wfc_k_irr,QE_eig)

      call rotate_wfc_test(wfc_k_irr,list_iG_irr,wfc_k, list_iG,         &
                           identity_matrix_index, sym, ftau_sym, G_sym)
    else

      ! Use the IBZ and symmetry!

      ! identify the right folder
      i_folder = QE_folder_sym(ikpt)

      call get_K_folder_data(i_folder,list_iG_irr,wfc_k_irr,QE_eig)

      ! The symmetry which takes k_irr to k
      i_sym    = symlink(ikpt,1)
      TR       = symlink(ikpt,2)
      G_sym    = sym_G(:,ikpt) !+ G_plus(:) !Asier&&Idoia 24 06 2014
      ftau_sym = ftau(:,i_sym)
      sym      = s(:,:,inverse_indices(i_sym))

      ! debugging test
      ! call rotate_wfc(wfc_k_irr,list_iG_irr,wfc_k, list_iG, &
      !                  inverse_indices(i_sym), sym, ftau_sym, G_sym)

      call rotate_wfc_test(wfc_k_irr,list_iG_irr,wfc_k, list_iG, &
                           i_sym, sym, ftau_sym, G_sym)

      ! If time-reversal is present, the wavefunction currently stored
      ! in wfc_k is actually for (-k). Complex conjugation must now
      ! be applied to recover wfc_k.
      if ( TR == 1)  then
        call apply_TR_to_wfc(wfc_k,list_iG)
      end if

    end if

    list_iG_irr = list_iG

    call wfc_by_expigr(kpoint, num_bands_intw, npol, ng_max, list_iG_irr, list_iG, wfc_k, G_plus)

    !    call wfc_by_expigr(kpoint, nbands, npol, ng_max, list_iG_irr, list_iG, wfc_k, G_plus)
    !(list_iG_irr, list_iG, G_plus)

  end subroutine get_psi_general_k
  !---------------------------------------------------------------------------
  !***************************************************************************
  !---------------------------------------------------------------------------
  subroutine intw_check_mesh(nk_1,nk_2,nk_3,kmesh,k_irr,nk_irr,kpoints_QE)
    !----------------------------------------------------------------------------!
    !     This subroutine compares the mesh points kmesh to the global array
    !     kpoints (which comes from the QE folders) in order to determine if
    !             1) kpoints = kmesh
    !             2) kpoints = irreducible zone
    !             3) none of the above
    !----------------------------------------------------------------------------!
    use intw_useful_constants, only: eps_8
    use intw_reading, only: nkpoints_QE

    implicit none

    integer :: nk_1,nk_2,nk_3, nkmesh
    real(kind=dp) :: kmesh(3,nk_1*nk_2*nk_3)

    real(kind=dp), intent(in) :: kpoints_QE(3,nkpoints_QE)

    integer :: nk_irr
    integer :: k_irr(3,nk_irr)

    real(kind=dp) :: kpt(3)
    integer :: ikpt
    real(kind=dp) :: norm

    nkmesh = nk_1*nk_2*nk_3

    full_mesh = .true.
    IBZ       = .true.

    ! First, establish if a full mesh is present
    if (nkpoints_QE /= nkmesh) then
      full_mesh = .false.
    else
       do ikpt=1,nkmesh

        norm = sqrt( (kpoints_QE(1,ikpt)-kmesh(1,ikpt))**2 + &
                     (kpoints_QE(2,ikpt)-kmesh(2,ikpt))**2 + &
                     (kpoints_QE(3,ikpt)-kmesh(3,ikpt))**2 )

        if (norm > eps_8) then
          full_mesh = .false.
          exit
        end if

       end do

    end if

    ! if a full mesh is NOT present, is the IBZ present?
    if (full_mesh .or. (nkpoints_QE /= nk_irr)) then
      IBZ = .false.
    else
      do ikpt=1,nk_irr
        kpt(1) = dble(k_irr(1,ikpt)-1)/nk_1
        kpt(2) = dble(k_irr(2,ikpt)-1)/nk_2
        kpt(3) = dble(k_irr(3,ikpt)-1)/nk_3

        norm = sqrt( (kpoints_QE(1,ikpt)-kpt(1))**2 + &
                     (kpoints_QE(2,ikpt)-kpt(2))**2 + &
                     (kpoints_QE(3,ikpt)-kpt(3))**2 )

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
    use intw_useful_constants, only: ZERO, eps_8
    use intw_utility, only: cryst_to_cart
    use intw_reading, only: bg, ngm, gvec

    implicit none

    integer :: number_G_shells, nG_shell_max

    integer :: nshell
    integer :: iG

    integer :: iflag, nvec

    real(kind=dp) :: G_cart(3)

    real(kind=dp) :: norm2, old_norm2

    nshell       = 1
    nG_shell_max = 1

    nvec = 1
    iflag = 1

    old_norm2 = ZERO

    do iG=2,ngm

      G_cart = dble(gvec(:,iG))
      call cryst_to_cart(nvec, G_cart, bg, iflag)

      norm2 = G_cart(1)**2 + G_cart(2)**2 + G_cart(3)**2

      if (abs(norm2-old_norm2) < eps_8) then
        nG_shell_max = nG_shell_max + 1
      else
        nshell = nshell + 1
        if (nshell == number_G_shells + 1) return
        nG_shell_max = nG_shell_max + 1
        old_norm2 = norm2
      end if

    end do

  end subroutine get_G_shells
  !--------------------------------------------
  !********************************************
  !--------------------------------------------
  subroutine apply_TR_to_wfc(wfc,list_iG)
    !--------------------------------------------
    !
    !-----------------------------------------------------------------------------
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
    !-----------------------------------------------------------------------------
    use intw_fft, only: find_iG
    use intw_useful_constants, only: cmplx_0
    use intw_utility, only: HPSORT
    use intw_reading, only: nG_max, gvec, nspin, num_bands_intw

    implicit none

    !I/O variables

    integer, intent(inout) :: list_iG(nG_max)
    complex(kind=dp), intent(inout) :: wfc(nG_max,num_bands_intw,nspin)

    !local variables

    integer :: iG, i_minus_G, i, p_i
    integer :: G(3), minus_G(3)
    integer :: nG ! counter on the number of G vectors in the array
    integer :: permutations(nG_max) ! index permutation which orders list_G
    integer :: list_iG_tmp(nG_max)
    complex(kind=dp) :: wfc_tmp(nG_max,num_bands_intw,nspin)

    ! Initialize the different variables
    !
    nG = 0
    permutations = 0
    list_iG_tmp = 0
    wfc_tmp = cmplx_0
    !
    ! loop on all G
    !
    do i=1,nG_max
      !
      iG = list_iG(i)
      !
      ! We work with G that contribute in the wfc
      !
      if (iG == 0) exit  ! the index array is zero-padded at the end.
      !
      nG = nG + 1
      G = gvec(:,iG)
      !
      minus_G = -G
      !
      ! find the index of -G
      !
      call find_iG(minus_G,i_minus_G)
      !
      list_iG_tmp(nG) = i_minus_G
      !
      ! conjugate the wavefunction
      !
      wfc_tmp(nG,:,:) = conjg(wfc(i,:,:))
      !
    enddo
    !
    ! There is no guarantee that the indices in list_iG_k will be sorted in ascending
    ! order! This is not an absolute necessity, but it would be nice and consistent for
    ! the indices to be sorted.
    ! Sort the indices using a canned heap sort subroutine.
    !
    call HPSORT(nG,list_iG_tmp,permutations)
    !
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
    !
    ! list_iG_tmp is now properly sorted, and can be dumped in the input/output variable
    !
    list_iG = list_iG_tmp
    !
    ! finally, populate the conjugated wave function
    !
    do i=1,nG
      !
      p_i = permutations(i)
      !
      ! compute the wfc element
      !
      if (nspin==1) then
        !
        wfc(i,:,:) = wfc_tmp(p_i,:,:)
        !
      elseif (nspin==2) then
        !
        wfc(i,:,1) = -wfc_tmp(p_i,:,2)
        wfc(i,:,2) = wfc_tmp(p_i,:,1)
        !
      endif !nspin
      !
    enddo !i

  end subroutine apply_TR_to_wfc
  !----------------------------------------------------------------------
  !**********************************************************************
  !----------------------------------------------------------------------
  subroutine find_size_of_irreducible_k_set(nk1,nk2,nk3,kmesh,nk_irr)
    !------------------------------------------------------------------
    ! This subroutine finds the number of k-points in the IBZ, based
    ! on the symmetry operations to be used.
    !------------------------------------------------------------------
    use intw_input_parameters, only: TR_symmetry
    use intw_useful_constants, only: eps_8
    use intw_utility, only: find_k_1BZ_and_G, switch_indices
    use intw_reading, only: nsym, s

    implicit none

    !input
    integer, intent(in) :: nk1, nk2, nk3
    ! The input k division
    real(kind=dp), intent(in) :: kmesh(3,nk1*nk2*nk3)
    ! The full kmesh, in canonical order

    !output
    integer, intent(out) :: nk_irr
    ! N. of irreducible k points found for the nk_1 nk_2 nk_3 division.

    !local variables
    real(kind=dp) :: k_rot(3), k_1BZ(3), dk(3), k_irr(3)
    integer :: nkpt ! The total number of points
    integer :: i, j, k ! triplet indices
    integer :: is, js, ks ! triplet indices  obtained by symmetry
    integer :: G(3)
    integer :: isym
    integer :: ikpt, ikpts ! singlet index, singlet index obtained by symmetry
    logical :: found(nk1*nk2*nk3)

    nkpt = nk1*nk2*nk3
    found = .false.
    nk_irr = 0
    !
    ! loop on the whole mesh, in the appropriate order (very important!)
    do i=1,nk1
      do j=1,nk2
        do k=1,nk3
          !
          ! find scalar index of point (i,j,k)
          call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,+1)
          !
          ! operate on this point only if it has not already been found!
          if (.not. found(ikpt)) then
            !
            ! it's found now. This point is part of the IBZ.
            found(ikpt) = .true.
            !
            nk_irr = nk_irr + 1
            !
            k_irr(1) = dble(i-1)/nk1
            k_irr(2) = dble(j-1)/nk2
            k_irr(3) = dble(k-1)/nk3
            !
            ! loop on all symmetry operations
            do isym=1,nsym
              !
              !perform matrix product
              ! CAREFUL! since the matrix is in crystal coordinates,
              ! and it acts in reciprocal space, the convention is :
              !          k_rot(i) = sum_j s(i,j)*k(j)
              !
              k_rot(:) = matmul(dble(s(:,:,isym)),k_irr(:))
              !
              ! find what point in the 1BZ this corresponds to
              call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)
              !
              ! check that k_1BZ+G = k_rot. If not, k_rot isn't on the mesh,
              ! and the algorithm in "find_k_1BZ_and_G" cannot be trusted.
              dk = k_rot - (k_1BZ + dble(G))
              if (sqrt(dot_product(dk,dk))<eps_8) then
                !
                ! what is the scalar index
                call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,+1)
                !
                if (.not. found(ikpts)) found(ikpts) = .true.
                !
              endif ! dk
              !
              ! Repeat, with Time-Reversal symmetry if present
              if (TR_symmetry) then
                !
                k_rot = -k_rot
                !
                ! find what point in the 1BZ this corresponds to
                call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)
                !
                ! we check again the value of dk, so if k_1BZ+G = k_rot
                dk = k_rot - (k_1BZ + dble(G))
                if (sqrt(dot_product(dk,dk))<eps_8) then
                  !
                  ! what is the scalar index
                  call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,+1)
                  !
                  if (.not.found(ikpts)) found(ikpts) = .true.
                  !
                endif ! dk
                !
              endif !TR_symmetry
              !
            enddo ! isym
            !
          endif ! found(ikpt)
          !
        enddo ! k
      enddo ! j
    enddo ! i

  end subroutine find_size_of_irreducible_k_set

  subroutine find_the_irreducible_k_set_irr(nk1,nk2,nk3,kmesh,nk_irr,list_ik_irr,kpoints_irr,weight_irr)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the canonical
    ! 1BZ mesh.
    !
    !------------------------------------------------------------------
    use intw_input_parameters, only: TR_symmetry
    use intw_useful_constants, only: eps_8
    use intw_utility, only: find_k_1BZ_and_G, switch_indices
    use intw_reading, only: nsym, s

    implicit none

    !input
    integer, intent(in) :: nk1, nk2, nk3, nk_irr
    ! The input k division
    real(kind=dp), intent(in) :: kmesh(3,nk1*nk2*nk3)
    ! The full kmesh, in canonical order

    !output
    integer, intent(out) :: list_ik_irr(nk_irr)
    real(kind=dp), intent(out) :: kpoints_irr(3,nk_irr)
    ! The irreducible kpoints in crystal triplet coordinates.
    ! The size of the array is nk1* nk2* nk3 instead of nk_irr;
    ! it is supposed that we still do not know the value of nk_irr
    real(kind=dp), intent(out) :: weight_irr(nk_irr)

    !local variables
    real(kind=dp) :: k_rot(3), k_1BZ(3), dk(3)
    integer :: nkpt ! The total number of points
    integer :: i, j, k ! triplet indices
    integer :: is, js, ks ! triplet indices  obtained by symmetry
    integer :: G(3)
    integer :: isym, ik_irr
    integer :: ikpt, ikpts ! singlet index, singlet index obtained by symmetry
    logical :: found(nk1*nk2*nk3)

    nkpt = nk1*nk2*nk3
    found = .false.
    ik_irr = 0
    weight_irr(:) = 0.d0
    !
    ! loop on the whole mesh, in the appropriate order (very important!)
    do i=1,nk1
      do j=1,nk2
        do k=1,nk3
          !
          ! find scalar index of point (i,j,k)
          call switch_indices(nk1,nk2,nk3,ikpt,i,j,k,+1)
          !
          ! operate on this point only if it has not already been found!
          if (.not. found(ikpt)) then
            !
            ! it's found now. This point is part of the IBZ.
            found(ikpt) = .true.
            !
            ik_irr = ik_irr + 1
            !
            weight_irr(ik_irr) = weight_irr(ik_irr) + 1
            !
            list_ik_irr(ik_irr) = ikpt
            !
            kpoints_irr(:,ik_irr) = kmesh(:,ikpt)
            !                kpoints_irr(1,ik_irr)=dble(i-1)/nk1
            !                kpoints_irr(2,ik_irr)=dble(j-1)/nk2
            !                kpoints_irr(3,ik_irr)=dble(k-1)/nk3
            !
            ! loop on all symmetry operations
            do isym=1,nsym
              !
              !perform matrix product
              ! CAREFUL! since the matrix is in crystal coordinates,
              ! and it acts in reciprocal space, the convention is :
              !          k_rot(i) = sum_j s(i,j)*k(j)
              !
              k_rot = matmul(dble(s(:,:,isym)),kpoints_irr(:,ik_irr))
              !
              ! find what point in the 1BZ this corresponds to
              call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)
              !
              ! check that k_1BZ+G = k_rot. If not, k_rot isn't on the mesh,
              ! and the algorithm in "find_k_1BZ_and_G" cannot be trusted.
              dk = k_rot - (k_1BZ + dble(G))
              if (sqrt(dot_product(dk,dk))<eps_8) then
                !
                ! what is the scalar index
                call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,+1)
                !
                if (.not. found(ikpts)) then
                    found(ikpts) = .true.
                    weight_irr(ik_irr) = weight_irr(ik_irr) + 1
                endif
                !
              endif ! dk
              !
              ! Repeat, with Time-Reversal symmetry if present
              if (TR_symmetry) then
                !
                k_rot = -k_rot
                !
                ! find what point in the 1BZ this corresponds to
                call find_k_1BZ_and_G(k_rot,nk1,nk2,nk3,is,js,ks,k_1BZ,G)
                !
                ! we check again the value of dk, so if k_1BZ+G = k_rot
                dk = k_rot - (k_1BZ + dble(G))
                if (sqrt(dot_product(dk,dk))<eps_8) then
                  !
                  ! what is the scalar index
                  call switch_indices(nk1,nk2,nk3,ikpts,is,js,ks,+1)
                  !
                  if (.not.found(ikpts)) then
                    found(ikpts) = .true.
                    weight_irr(ik_irr) = weight_irr(ik_irr) + 1
                  endif
                  !
                endif ! dk
                !
              endif !TR_symmetry
              !
            enddo ! isym
            !
          endif ! found(ikpt)
          !
        enddo ! k
      enddo ! j
    enddo ! i

    weight_irr(:) = weight_irr(:)/nkpt

  end subroutine find_the_irreducible_k_set_irr

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
    use intw_useful_constants, only: eps_8, ONE
    use intw_reading, only: nsym, s, can_use_TR, gvec

    implicit none

    ! input
    integer :: nG_shell_max
    real(kind=dp) :: q_direction(3)
    ! output
    integer :: shells_table(nG_shell_max,nsym)

    ! local
    integer :: iG, iG2
    integer :: isym

    real(kind=dp) :: rot_q(3)
    integer :: restricted_nsym

    integer :: ishell, iG_shell

    logical :: G_vector_found(nG_shell_max)

    integer :: G(3), G2(3), RG(3)
    real(kind=dp) :: dG(3)

    real(kind=dp) :: norm

    real(kind=dp) :: local_rotations(3,3,nsym)
    logical :: local_TR(nsym)


    ! first, find which symmetries can be used

    restricted_nsym = 0

    do isym=1,nsym
      rot_q(:) = matmul(s(:,:,isym),q_direction)

      if (can_use_TR(isym)) rot_q(:) = -rot_q(:)


      norm = sqrt(dot_product(rot_q-q_direction,rot_q-q_direction))
      if (norm < eps_8) then
        ! add the symmetry to the  set
        restricted_nsym = restricted_nsym + 1

        local_rotations(:,:,restricted_nsym) = s(:,:,isym)
        local_TR(restricted_nsym) = can_use_TR(isym)

      end if

    end do


    ! second, populate the table which contains the G vectors in shells!


    G_vector_found(:) = .false.

    ! This variable will contain the G vector indices, where each row
    ! will represent a rotation shell.
    shells_table(:,:) = 0

    ishell = 0
    iG_shell = 0

    do iG=1,nG_shell_max

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
      do isym=1,restricted_nsym
        RG(:) = matmul(local_rotations(:,:,isym),G)

        if (local_TR(isym)) RG(:) = -RG(:)

        ! look for this vector in the G vectors not yet found
        do iG2=1,nG_shell_max

          if ( G_vector_found(iG2) ) cycle

          G2(:) = gvec(:,iG2)

          dG(:) = ONE*(RG(:) - G2(:))
          norm = sqrt(dot_product(dG,dG))

          if (norm < eps_8) then
            ! add the new vector to the shell
            G_vector_found(iG2) = .true. ! it is found now!

            iG_shell = iG_shell + 1

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
    use intw_utility, only: find_free_unit
    use intw_reading, only: nsym, bg, gvec

    implicit none

    ! input
    integer :: nG_shell_max
    real(kind=dp) :: q_direction(3)
    ! output
    integer :: shells_table(nG_shell_max,nsym)

    ! local
    integer :: io_unit

    integer :: ishell, iG_shell, iG1

    integer :: G(3)


    real(kind=dp) :: q_cart(3), G_cart(3)

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



    do ishell=1,nG_shell_max

      iG1 = shells_table(ishell,1)

      if (iG1 == 0) exit
      write(io_unit,5) '                                             '

      do iG_shell=1,nsym
        iG1 = shells_table(ishell,iG_shell)
        if (iG1 == 0) exit

        G(:) = gvec(:,iG1)

        G_cart(:) = matmul(bg,dble(G))

        write(io_unit,15) ishell,iG1,G(:), G_cart

      end do

    end do

    close(unit = io_unit)

5   format(A)
10  format(A,3G18.8E3)
15  format(I8,8X,I6,5X,3I4,1x,3F8.2)


  end subroutine output_sorted_G_vectors_in_symmetry_shells

  logical function eqvect(x,y,f)
    !-----------------------------------------------------------------------
    !
    use kinds, only: dp

    implicit none
    real(kind=dp) :: x(3), y(3), f(3)
    ! input: input vector
    ! input: second input vector
    ! input: fractionary translation
    real(kind=dp), parameter :: accep = 1.0d-4
    ! acceptance parameter
    !
    eqvect = abs( x(1)-y(1)-f(1) - nint(x(1)-y(1)-f(1)) ) .lt. accep .and. &
             abs( x(2)-y(2)-f(2) - nint(x(2)-y(2)-f(2)) ) .lt. accep .and. &
             abs( x(3)-y(3)-f(3) - nint(x(3)-y(3)-f(3)) ) .lt. accep
    return
  end function eqvect

  !----------------------------------------
  subroutine multable(nsym,s,table)
    !----------------------------------------
    !
    !--------------------------------------------------------------------
    !  sets up the multiplication table for a group represented by 3x3
    !  integer matrices and checks that {s} is a group indeed:
    !
    !  table(n,m) = index( s(n)*s(m) )
    !--------------------------------------------------------------------

    implicit none

    !I/O variables

    integer, intent(in) :: nsym, s(3,3,48) ! number of symmetry+symm matrices
    integer, intent(out) :: table(48,48) ! the multiplication table

    !local variables

    integer :: irot, jrot, krot, ipol, jpol, kpol, ss(3,3)
    ! counter on rotations
    ! counters on polarizations
    ! buffer multiplication matrix
    logical :: found, smn
    ! if true the table has been set used to check symmetries

    do irot=1,nsym
      do jrot=1,nsym
        !
        do ipol=1,3
          do jpol=1,3
            !
            ss(ipol,jpol) = 0
            !
            do kpol=1,3
              !
              ss(ipol,jpol) = ss(ipol,jpol) + &
                              s(ipol,kpol,jrot)*s(kpol,jpol,irot)
              !
            enddo !kpol
            !
          enddo !jpol
        enddo !ipol
        !
        ! here checks that the input matrices really form a group
        ! and sets the multiplication table
        !
        found = .false.
        !
        do krot=1,nsym
          !
          smn = .true.
          !
          do ipol=1,3
            do jpol=1,3
              !
              smn = smn .and. (s(ipol,jpol,krot).eq.ss(ipol,jpol))
              !
            enddo !jpol
          enddo !ipol
          !
          if (smn) then
            !
            if (found) write(*,*)'something wrong in multable:1'
            !
            found = .true.
            table(jrot,irot) = krot
            !
          endif
          !
        enddo !krot
      enddo !jrot
      !
      if (.not.found) write(*,*)'something wrong in multable:2'
      !
    enddo !irot
    !
    return

  end subroutine multable
  !----------------------------------------------------------------------------!
end module intw_symmetries
!----------------------------------------------------------------------------!
