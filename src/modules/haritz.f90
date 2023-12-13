module intw_haritz
  !
  ! This module contains general purpose variables and subroutines
  ! that I use in many programs.
  !
  use kinds, only: dp
  !
  implicit none
  !
  save
  !
  ! Declare public variables
  ! Conversion parameters
  public :: eV_to_Ha, eV_to_Ry, Ha_to_eV, Ha_to_Ry, Ry_to_eV, Ry_to_Ha
  ! I/O variables
  public :: stdout, stdin, CR, LF
  !
  private
  !
  ! Energy unit conversion parameters
  real(kind=dp), parameter :: eV_to_Ha = 1.0_dp/27.211383860484776_dp
  real(kind=dp), parameter :: eV_to_Ry = 2.0_dp/27.211383860484776_dp
  real(kind=dp), parameter :: Ha_to_eV = 27.211383860484776_dp
  real(kind=dp), parameter :: Ha_to_Ry = 2.0_dp
  real(kind=dp), parameter :: Ry_to_eV = 27.211383860484776_dp/2.0_dp
  real(kind=dp), parameter :: Ry_to_Ha = 0.5_dp
  !
  ! I/O variables
  integer         , parameter :: stdout=6    ! standard output unit
  integer         , parameter :: stdin=5     ! standard input unit
  character(len=1), parameter :: CR=CHAR(13) ! carriage return character
  character(len=1), parameter :: LF=CHAR(10) ! carriage return character
  !
end module intw_haritz
