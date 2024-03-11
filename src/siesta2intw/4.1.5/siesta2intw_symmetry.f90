module siesta2intw_symmetry

  use precision, only: dp
  use spglib_f08, only: spg_get_dataset, SpglibDataset

  implicit none

  public :: compute_symmetry, print_symmetry_data

  public :: nsym, s, ftau, t_rev

  private

  ! symmetry variables
  integer :: nsym
  integer, allocatable, dimension(:,:,:) :: s
  real(kind=dp), allocatable, dimension(:,:) :: ftau
  integer, allocatable, dimension(:) :: t_rev


contains

  subroutine compute_symmetry(use_sym)
    !
    ! Compute the symmetry of the system using spglib
    !

    ! variables
    use precision, only: dp
    use siesta_geom, only: na_u, ucell, xa, isa
    use m_spin, only: spin
    use parallel, only: ionode
    use siesta2intw_io, only: stdout
    ! functions and subroutines
    use siesta2intw_utils, only: atomic_magnetic_moment
    !
    implicit none
    !
    logical, intent(in) :: use_sym
    !
    real(kind=dp), dimension(3, na_u) :: xa_cryst
    integer, dimension(na_u) :: ityp
    real(kind=dp), dimension(3, 3) :: bg
    !
    type(SpglibDataset) :: spg_sym
    !
    real(kind=dp) :: Svec(3, na_u) ! Spin magnetic moment vector for each atom in unit cell
    real(kind=dp) :: Stot(na_u) ! Spin magnetic moment magnitude for each atom in unit cell

    logical :: magnetic_sym
    integer :: ia, isym, jsym


    call reclat( transpose(ucell), bg, -1 )
    !
    do ia = 1, na_u
      xa_cryst(:, ia) = matmul(bg, xa(:, ia))
      ityp(ia) = isa(ia)
    enddo
    !
    ! Check if there are any atoms with non-zero spin magnetic moments.
    if ( spin%none ) then
      !
      magnetic_sym = .false.
      !
    else
      !
      call atomic_magnetic_moment(Svec, Stot)
      !
      if  ( any(abs(Stot) > 1.e-6) ) then
        if (ionode) write(stdout, *) "WARNING: Net magnetic moment detected: Symmetries will not be used!"
        magnetic_sym = .true.
      else
        if (ionode) write(stdout, *) "WARNING: No net magnetic moment detected: Symmetries will be used!"
        magnetic_sym = .false.
      endif
      !
    endif
    !
    ! Compute symmetries
    if (use_sym .and. (.not. magnetic_sym)) then
      !
      spg_sym = spg_get_dataset(transpose(ucell), xa_cryst, ityp, na_u, 1.d-5)
      !
      nsym = spg_sym%n_operations
      allocate(s(3, 3, nsym))
      allocate(ftau(3, nsym))
      allocate(t_rev(nsym))
      !
      do isym = 1, nsym
        s(:, :, isym) = spg_sym%rotations(:, :, isym)
        ftau(:, isym) = -spg_sym%translations(:, isym)
        t_rev (isym) = 0 ! TODO: I should compute the magnetic space group symmetries when needed and check whether time reversal is enabled or not
      enddo
      !
    else
      !
      nsym = 1
      allocate(s(3, 3, nsym))
      allocate(ftau(3, nsym))
      allocate(t_rev(nsym))
      s(1, :, 1) = (/1, 0, 0/)
      s(2, :, 1) = (/0, 1, 0/)
      s(3, :, 1) = (/0, 0, 1/)
      ftau = 0.0_dp
      t_rev = 0
      !
    endif

  end subroutine compute_symmetry


  subroutine print_symmetry_data()
    !
    ! Print the symmetry data of the system in a nice format.
    !

    ! variables
    use siesta2intw_io, only: stdout

    implicit none

    integer :: isym


    write(stdout, *) "- Symmetry operations:"

    write(stdout, *) " nsym"
    write(stdout, *) nsym
    !
    23 format("  - [", i2, ", ", i2, " ,", i2, "]")
    24 format("  translation: [", f8.5, ", ", f8.5, ", ", f8.5, " ]")
    do isym = 1, nsym
      write(stdout, *) " s", isym
      write(stdout, 23) s(1, :, isym)
      write(stdout, 23) s(2, :, isym)
      write(stdout, 23) s(3, :, isym)
      write(stdout, 24) ftau(:, isym)
    enddo

  end subroutine print_symmetry_data

end module siesta2intw_symmetry
