!
! Copyright (C) 2024 INTW group
!
! This file is part of INTW.
!
! INTW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! INTW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
PROGRAM create_intw_q_dirs
  !-----------------------------------------------------------------------
  !
  USE Kinds
  USE mp_global, ONLY: mp_startup

  USE io_files, ONLY: prefix_QE => prefix
  USE environment, ONLY: environment_start, environment_end
  USE cell_base, ONLY: bg

  IMPLICIT NONE

  EXTERNAL :: errore, read_file

  ! I/O
  CHARACTER(len=256) :: prefix = " "
  CHARACTER(len=256) :: outdir = "./"
  CHARACTER(len=256) :: qlist_file = "qlist.txt"
  CHARACTER(len=256) :: phdir = "./"
  INTEGER :: nq1, nq2, nq3
  CHARACTER(len=256) :: reference_file

  real(kind=dp), allocatable, dimension(:,:) :: qirr_cryst
  real(kind=dp), dimension(3) :: qirr_cart
  integer :: nq_irr
  character(len=256) :: q_dir
  integer :: iq
  integer :: iounit, iostat, strlen
  integer, external :: find_free_unit


  NAMELIST / inputpp / prefix, outdir, phdir, nq1, nq2, nq3, reference_file


  !
  ! Initialize MPI and environment
#if defined(__MPI)
  CALL mp_startup( )
#endif

  CALL environment_start( "create_intw_q_dirs" )

  !
  ! Read input file
  READ (5, inputpp, iostat=iostat)
  if (iostat /= 0) call errore( "create_intw_q_dirs", "ERROR: create_intw_q_dirs: error reading inputpp", 1 )

  strlen = len_trim(outdir)
  if ( outdir(strlen:strlen+1) .ne. "/" ) outdir(strlen+1:strlen+2) = "/"
  strlen = len_trim(phdir)
  if ( phdir(strlen:strlen+1) .ne. "/" ) phdir(strlen+1:strlen+2) = "/"

  !
  ! Read QE data
  prefix_QE = prefix
  CALL read_file()

  !
  ! Find q-mesh and irreducible q-points
  allocate(qirr_cryst(3,nq1*nq2*nq3))
  call find_the_irreducible_k_set(nq1, nq2, nq3, qirr_cryst, nq_irr)

  !
  ! Create qlist_file and qq directories
  iounit = find_free_unit()
  open(unit=iounit, file=trim(outdir)//trim(qlist_file), status="replace", action="write", iostat=iostat)
  if ( iostat /= 0 ) call errore( "create_intw_q_dirs", "ERROR: create_intw_q_dirs: Error opening qlist_file", 1 )
  !
  do iq=1,nq_irr
    !
    call write_tag("qq", iq, q_dir)
    !
    ! Create directory
    call execute_command_line("mkdir -p "//trim(outdir)//trim(phdir)//trim(q_dir))
    !
    ! Write input file
    qirr_cart = matmul(bg, qirr_cryst(:, iq))
    call write_ph_in(trim(outdir)//trim(phdir)//trim(reference_file), &
                     trim(outdir)//trim(phdir)//trim(q_dir)//"/"//trim(reference_file), qirr_cart)
    !
    write(iounit, "(i3,3f18.10)") iq, qirr_cart
    !
  enddo
  !
  close(iounit)

  !
  ! End job
  call environment_end( "create_intw_q_dirs" )

contains

  subroutine find_the_irreducible_k_set(nk1, nk2, nk3, kirr_cryst, nkirr)
    !------------------------------------------------------------------
    ! This subroutine finds the irreducible k set for the canonical
    ! 1BZ mesh.
    !------------------------------------------------------------------
    USE symm_base, ONLY: nsym, s, t_rev, time_reversal


    implicit none

    ! I/O

    integer, intent(in) :: nk1, nk2, nk3
    real(kind=dp), intent(out) :: kirr_cryst(3,nk1*nk2*nk3)

    integer, intent(out) :: nkirr
    ! How many points were found?

    !local variables
    real(kind=dp) :: k_rot(3), k_1BZ(3), dk(3)
    integer :: G(3)
    integer :: i, j, k ! triple indices
    integer :: is, js, ks ! triple indices obtained by symmetry
    integer :: isym
    logical :: found(nk1,nk2,nk3)
    real(kind=dp), parameter :: eps_8 = 1.0e-8_dp



    kirr_cryst = 0.0_dp
    found = .false.
    nkirr = 0
    !
    ! loop on the whole mesh
    do i=1,nk1
      do j=1,nk2
        do k=1,nk3
          !
          ! operate on this point only if it has not already been found!
          if (found(i,j,k)) cycle
          !
          ! it's found now. This point is part of the IBZ.
          found(i,j,k) = .true.
          !
          nkirr = nkirr + 1
          !
          kirr_cryst(1,nkirr) = dble(i-1)/nk1
          kirr_cryst(2,nkirr) = dble(j-1)/nk2
          kirr_cryst(3,nkirr) = dble(k-1)/nk3
          !
          ! loop on all symmetry operations
          do isym=1,nsym
            !
            !perform matrix product
            ! CAREFUL! since the matrix is in crystal coordinates,
            ! and it acts in reciprocal space, the convention is :
            !          k_rot(i) = sum_j s(i,j)*k(j)
            !
            k_rot = matmul(dble(s(:,:,isym)), kirr_cryst(:,nkirr))
            !
            ! Time-Reversal symmetry
            if (t_rev(isym)==1) k_rot = -k_rot
            !
            ! find what point in the 1BZ this corresponds to
            call find_k_1BZ_and_G(k_rot, nk1, nk2, nk3, is, js, ks, k_1BZ, G)
            !
            ! we check again the value of dk, so if k_1BZ+G = k_rot
            dk = k_rot - (k_1BZ + dble(G))
            if ( norm2(dk) < eps_8 ) found(is,js,ks) = .true.
            !
            ! q -> -q
            if (time_reversal) then
              !
              k_rot = -k_rot
              !
              ! find what point in the 1BZ this corresponds to
              call find_k_1BZ_and_G(k_rot, nk1, nk2, nk3, is, js, ks, k_1BZ, G)
              !
              ! we check again the value of dk, so if k_1BZ+G = k_rot
              dk = k_rot - (k_1BZ + dble(G))
              if ( norm2(dk) < eps_8 ) found(is,js,ks) = .true.
              !
            endif
            !
          enddo ! isym
          !
        enddo ! k
      enddo ! j
    enddo ! i

  end subroutine find_the_irreducible_k_set


  subroutine find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpt_in_1BZ,G)
    !----------------------------------------------------------------------------!
    ! Given a kpoint(3) in crystal coordinates, this subroutine
    ! generates:
    !            - k_1BZ(3),   with 0 <= k_1BZ(i) < 1
    !            - G(3),       G(i) an integer
    !            - i,j,k       the triple coordinates of the point in the mesh.
    !
    ! It it EXTREMELY important to have a subroutine which performs this
    ! task in a defensive way. The fortran internal functions nint, floor,
    ! and modulo are DANGEROUS. IT IS CRUCIAL TO DO THIS RIGHT ONCE AND FOR
    ! ALL.
    !
    ! The basic relationship is
    !    kpoint = k_1BZ + G
    !
    ! Assume
    !    kpoint(l) =  G(l) + (i_l-1)/nk_l+ epsilon ,
    ! with
    !    l    = 1, 2, 3
    !    i_l  = i, j ,k
    !    nk_l = nk1, nk2, nk3
    !    G(l) an integer
    ! epsilon numerical noise
    !----------------------------------------------------------------------------!
    implicit none

    ! input
    integer, intent(in) :: nk1, nk2, nk3             ! the mesh parameters
    real(dp), intent(in) :: kpoint(3)                ! the kpoint of interest

    ! output
    integer, intent(out) :: i, j, k                  ! the triple coordinates of k_1BZ
    real(dp), intent(out) :: kpt_in_1BZ(3)           ! the k point in the 1BZ
    integer, intent(out) :: G(3)                     ! the translation vector

    ! internal variables
    integer :: nG_im1, nG_jm1, nG_km1
    integer :: im1, jm1, km1


    ! this step kills epsilon
    nG_im1 = nint(kpoint(1)*dble(nk1))
    nG_jm1 = nint(kpoint(2)*dble(nk2))
    nG_km1 = nint(kpoint(3)*dble(nk3))

    ! this step gets rid of G
    im1 = modulo(nG_im1,nk1)
    jm1 = modulo(nG_jm1,nk2)
    km1 = modulo(nG_km1,nk3)

    ! G can be extracted. This division must be exact.
    G(1) = (nG_im1 - im1)/nk1
    G(2) = (nG_jm1 - jm1)/nk2
    G(3) = (nG_km1 - km1)/nk3

    ! finally we have the triple coordinates
    i = im1 + 1
    j = jm1 + 1
    k = km1 + 1

    ! compute the k point in the 1BZ
    kpt_in_1BZ(1) = dble(i-1)/dble(nk1)
    kpt_in_1BZ(2) = dble(j-1)/dble(nk2)
    kpt_in_1BZ(3) = dble(k-1)/dble(nk3)

  end subroutine find_k_1BZ_and_G


  subroutine write_ph_in(reference_file, filename, qpoint_cart)
    !
    !
    !
    use, intrinsic :: iso_fortran_env, only: iostat_end ! end of file
    !
    implicit none
    !
    ! I/O
    !
    character(len=*), intent(in) :: reference_file, filename
    real(kind=dp), intent(in) :: qpoint_cart(3)
    !
    ! Local
    !
    character(len=256) :: line
    integer :: iounit_reference, iounit
    logical :: fildvscf_found = .false.
    integer :: ios
    integer, external :: find_free_unit


    !
    iounit_reference = find_free_unit()
    open(unit=iounit_reference, file=trim(outdir)//trim(reference_file), status="old", action="read", iostat=ios)
    if ( ios /= 0 ) stop "write_fdf: Error opening file outdir/reference_file"
    !
    iounit = find_free_unit()
    open(unit=iounit, file=trim(outdir)//trim(filename), status="replace", action="write", iostat=ios)
    if ( ios /= 0 ) stop "write_fdf: Error opening file outdir/filename"
    !
    do
      !
      ! Read each line of the input file
      read(unit=iounit_reference, fmt="(a)", iostat=ios) line
      !
      ! Some checks
      if (word_in_string("trans", line)) then
        ! If .true. the phonons are computed.
        ! If trans .and. epsil are .true. effective charges are calculated.
        ! By default it is true
        if (word_in_string("false", line)) call errore( &
          "create_q_dirs", "write_ph_in: trans = .true. must be used to compute phonons", 1 )
      endif
      !
      if (word_in_string("ldisp", line)) then
        ! If .true. the run calculates phonons for a grid of
        ! q-points specified by nq1, nq2, nq3 - for direct
        ! calculation of the entire phonon dispersion.
        ! By default it is false
        if (word_in_string("true", line)) call errore( &
          "create_q_dirs", "write_ph_in: ldisp = .false. must be used to compute phonons with create_q_dirs", 1 )
      endif
      !
      if (word_in_string("fildvscf", line)) then
        ! If it is not specified the induced potentials will not be written
        fildvscf_found = .true.
      endif
      !
      ! Exit the do loop if it is the last line
      if ( ios == iostat_end ) exit
      !
      ! Write the lane
      write(iounit,"(a)") trim(line)
      !
    end do
    !
    close(iounit_reference)
    !
    ! Check if the induced potentials will be written
    if (.not.fildvscf_found) call errore( "create_q_dirs","write_ph_in: fildvscf must be used to save the induced potentials", 1 )
    !
    ! Write q point
    write(iounit,"(3f15.8)") qpoint_cart
    !
    close(iounit)

  end subroutine write_ph_in


  function word_in_string(word,string)
    !
    ! Check if a word is present in a string.
    !
    implicit none
    !
    character(*), intent(in) :: word, string
    logical :: word_in_string
    !
    character(len=len(word)) :: tmp_word
    character(len=len(string)) :: tmp_string
    integer :: len_word, len_string


    tmp_word = string_tolower(adjustl(word))
    len_word = len(trim(tmp_word))
    tmp_string = string_tolower(adjustl(string))
    len_string = len(trim(tmp_string))
    !
    if ( index(tmp_string,tmp_word) == 0 ) then
      word_in_string = .false.
    else
      word_in_string = .true.
    endif

  end function word_in_string


  function string_tolower(string)
    !
    ! Return the string in lowercase.
    !
    implicit none
    !
    character(*) :: string
    character(len=len(string)) :: string_tolower
    !
    integer :: i, k, length


    length = len(string)
    string_tolower = string
    do i = 1,len(string)
      k = iachar(string(i:i))
      if ( k >= iachar("A") .and. k <= iachar("Z") ) then
        k = k + iachar("a") - iachar("A")
        string_tolower(i:i) = achar(k)
      endif
    enddo

  end function string_tolower


  subroutine write_tag(string,i,tag)
    !-----------------------------------------------
    ! This subroutine creates a character string of
    ! the form "string"integer, where the integer
    ! will be immediately after the end of "string",
    ! without blank spaces.
    !-----------------------------------------------
    implicit none

    integer :: i
    character(*) :: string
    character(256) :: integer_part, tag


    if (i < 10) then
      write(integer_part,100) i
    elseif (i < 100 ) then
      write(integer_part,200) i
    elseif (i < 1000 ) then
      write(integer_part,300) i
    elseif (i < 10000 ) then
      write(integer_part,400) i
    elseif (i < 100000 ) then
      write(integer_part,500) i
    end if

    tag = trim(string)//trim(integer_part)

    100 format(I1)
    200 format(I2)
    300 format(I3)
    400 format(I4)
    500 format(I5)

  end subroutine write_tag

end PROGRAM create_intw_q_dirs
