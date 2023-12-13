module intw_useful_constants
  !       The name is "useful_constants" to avoid clashes with W90 "constants.f90"
  !
  use kinds, only: dp
  !
  implicit none
  !
  ! I/O variables
  public :: stdout, stdin
  !
  ! Energy conversion factors
  public :: eV_to_Ha, eV_to_Ry, Ha_to_eV, Ha_to_Ry, Ry_to_eV, Ry_to_Ha
  !
  public :: ZERO, ONE, e2, sqrt2, pi, tpi, fpi, sqrt_pi, &
            eps_2, eps_5, eps_6, eps_7, eps_8, eps_10, eps_14, &
            bohr, boltzmann, cmplx_i, cmplx_1, cmplx_0, &
            direct_io_factor, direct_io_factor_cmplx, double_complex, double_real, &
            I2, sig_x, sig_y, sig_z
  !
  private
  !
  save
  !
  !
  !NOTE: haritz 26/11/2021
  ! All this variables should be parameters, but if they are defined as parameters
  ! OpenMP gives an error in some points
  real(kind=dp), parameter :: ZERO       = 0.0_dp
  real(kind=dp), parameter :: ONE        = 1.0_dp
  real(kind=dp), parameter :: e2         = 2.0_dp
  real(kind=dp), parameter :: sqrt2      = sqrt(2.0_dp)
  real(kind=dp), parameter :: pi         = dacos(-1.0_dp)
  real(kind=dp), parameter :: tpi        = 2.0_dp*pi
  real(kind=dp), parameter :: fpi        = 4.0_dp*pi
  real(kind=dp), parameter :: sqrt_pi    = sqrt(pi)
  real(kind=dp), parameter :: eps_2      = 1.0e-2_dp
  real(kind=dp), parameter :: eps_5      = 1.0e-5_dp
  real(kind=dp), parameter :: eps_6      = 1.0e-6_dp
  real(kind=dp), parameter :: eps_7      = 1.0e-7_dp
  real(kind=dp), parameter :: eps_8      = 1.0e-8_dp
  real(kind=dp), parameter :: eps_10     = 1.0e-10_dp
  real(kind=dp), parameter :: eps_14     = 1.0e-14_dp
  real(kind=dp), parameter :: bohr       = 0.52917720859_dp
  real(kind=dp), parameter :: boltzmann  = 8.6173324*0.00001_dp
  complex(kind=dp), parameter :: cmplx_i = (0.0_dp, 1.0_dp)
  complex(kind=dp), parameter :: cmplx_1 = (1.0_dp, 0.0_dp)
  complex(kind=dp), parameter :: cmplx_0 = (0.0_dp, 0.0_dp)
  !
  ! Energy conversion factors
  real(kind=dp), parameter :: Ha_to_eV = 27.211383860484776_dp
  real(kind=dp), parameter :: eV_to_Ha = 1.0_dp/Ha_to_eV
  real(kind=dp), parameter :: Ha_to_Ry = 2.0_dp
  real(kind=dp), parameter :: Ry_to_Ha = 1.0_dp/Ha_to_Ry
  real(kind=dp), parameter :: Ry_to_eV = Ry_to_Ha*Ha_to_eV
  real(kind=dp), parameter :: eV_to_Ry = 1.0_dp/Ry_to_eV
  !
  ! I/O variables
  integer, parameter :: stdout = 6 ! standard output unit
  integer, parameter :: stdin = 5 ! standard input unit
  !
  !-Pauli matrices
  complex(kind=dp) :: I2(2,2), sig_x(2,2), sig_y(2,2), sig_z(2,2)

  !TODO: haritz 26/11/2021
  ! As the following parameters could be compiler dependent they should be
  ! removed from the code and use inquire to obtain the record legth
  ! this is for a real double
  integer, parameter :: direct_io_factor = 2
  ! this is for a complex double
  integer, parameter :: direct_io_factor_cmplx = 4
  ! what is the size in bytes of a complex double
  integer, parameter :: double_complex = 16
  ! what is the size in bytes of a complex real
  integer, parameter :: double_real = 8

end module  intw_useful_constants
