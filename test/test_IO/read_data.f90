integer function read_data_test() result(r)

  use kinds, only: dp
  use intw_test_module, only: assert
  use intw_reading, only: read_parameters_data_file_xml
  use intw_reading, only: lspinorb, nsym, atom_pfile, atom_labels, at, bg, alat, &
                          nat, ntyp, ityp, tau, amass, nr1, nr2, nr3, ngm, &
                          noncolin, spinorb_mag, ecutwfc, ecutrho
  use intw_input_parameters, only: mesh_dir, prefix
  use intw_useful_constants, only: twopi

  implicit none

  real(kind=dp), parameter :: prec = 1e-6

  integer :: ngm_test
  real(kind=dp) :: alat_test = 7.200000000
  real(kind=dp), dimension(3,3) :: at_test
  real(kind=dp), dimension(3,3) :: bg_test
  real(kind=dp) :: ecutwfc_test = 80.0000000000000
  real(kind=dp) :: ecutrho_test = 320.000000000000
  integer :: nr1_test = 30, nr2_test = 30, nr3_test = 30
  logical :: noncolin_test = .false.
  logical :: lspinorb_test = .false.
  logical ::spinorb_mag_test = .false.
  integer, parameter :: nat_test=1
  integer, parameter :: ntyp_test=1
  real(kind=dp), dimension(ntyp_test) :: amass_test = (/ 63.54600 /)
  character(len=3), dimension(ntyp_test)  :: atom_labels_test = (/ "Cu" /)
  character(len=80), dimension(ntyp_test) :: atom_pfile_test = (/ "Cu.RRKJ3.UPF" /)
  integer,  dimension(nat_test) :: ityp_test
  real(kind=dp), dimension(3,nat_test) :: tau_test
  integer :: nsym_test = 48

  r = 0

  mesh_dir = "./qe/"
  prefix = "cu"

  call read_parameters_data_file_xml()

  if (.not.assert(alat,alat_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(ecutwfc,ecutwfc_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(ecutrho,ecutrho_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(nr1,nr1_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(nr2,nr2_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(nr3,nr3_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(noncolin,noncolin_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(lspinorb,lspinorb_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(spinorb_mag,spinorb_mag_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(nat,nat_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(ntyp,ntyp_test,prec)) then
    r = 1
    return
  endif

  if (.not.assert(nsym,nsym_test,prec)) then
    r = 1
    return
  endif


end function

! interface assert
!
!   logical function assert_integer(a,b,prec) result(r)
!     use kinds, only: dp
!     integer, intent(in) :: a, b
!     real(kind=dp), intent(in) :: prec
!   end function
!
!   logical function assert_integer_vec(a,b,prec) result(r)
!     use kinds, only: dp
!     integer, intent(in), dimension(:) :: a, b
!     real(kind=dp), intent(in) :: prec
!   end function
!
!   logical function assert_real(a,b,prec) result(r)
!     use kinds, only: dp
!     real(kind=dp), intent(in) :: a, b
!     real(kind=dp), intent(in) :: prec
!   end function
!
!   logical function assert_real_vec(a,b,prec) result(r)
!     use kinds, only: dp
!     real(kind=dp), intent(in), dimension(:) :: a, b
!     real(kind=dp), intent(in) :: prec
!   end function
!
! end interface assert

! logical function assert_logical(a,b,prec) result(r)
!
!   use kinds, only: dp
!
!   logical, intent(in) :: a, b
!   real(kind=dp), intent(in) :: prec
!
!   r = .true.
!   if (a /= b) then
!     r = .false.
!   endif
!
! end function
!
!
! logical function assert_integer(a,b,prec) result(r)
!
!   use kinds, only: dp
!
!   integer, intent(in) :: a, b
!   real(kind=dp), intent(in) :: prec
!
!   r = .true.
!   if (a /= b) then
!     r = .false.
!   endif
!
! end function
!
! logical function assert_integer_vec(a,b,prec) result(r)
!
!   use kinds, only: dp
!
!   integer, intent(in), dimension(:) :: a, b
!   real(kind=dp), intent(in) :: prec
!
!   if (size(a) /= size(b)) then
!     r = .false.
!     write(6,*) "WARNING: assert_integer_vec: different size vectors"
!     return
!   endif
!
!   r = .true.
!   do i=1,size(a)
!     if (a(i) /= b(i)) then
!       r = .false.
!       return
!     endif
!   enddo
!
! end function
!
! logical function assert_real(a,b,prec) result(r)
!
!   use kinds, only: dp
!
!   real(kind=dp), intent(in) :: a, b
!   real(kind=dp), intent(in) :: prec
!
!   r = .true.
!   if (abs(a-b) > prec) then
!     r = .false.
!   endif
!
! end function
!
! logical function assert_real_vec(a,b,prec) result(r)
!
!   use kinds, only: dp
!
!   real(kind=dp), intent(in), dimension(:) :: a, b
!   real(kind=dp), intent(in) :: prec
!
!   if (size(a) /= size(b)) then
!     r = .false.
!     write(6,*) "WARNING: assert_real_vec: different size vectors"
!     return
!   endif
!
!   r = .true.
!   do i=1,size(a)
!     if (abs(a(i)-b(i)) > prec) then
!       r = .false.
!       return
!     endif
!   enddo
!
! end function
