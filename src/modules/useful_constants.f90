module intw_useful_constants
  !       The name is "useful_constants" to avoid clashes with W90 "constants.f90"
  !
  use kinds, only: dp
  !
  implicit none
  !
  public :: ZERO, ONE, sqrt2, e2, cmplx_i, cmplx_1, cmplx_0, pi, sqrt_pi, &
            tpi, fpi, eps_2, eps_5, eps_6, eps_7, eps_8, &
            eps_10, eps_14, bohr, boltzmann, direct_io_factor, &
            direct_io_factor_cmplx, double_complex, double_real, &
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
  real(dp), parameter :: ZERO      = 0.0_dp
  real(dp), parameter :: ONE       = 1.0_dp
  real(dp), parameter :: e2        = 2.0_dp
  real(dp) :: pi        = dacos(-1.0_dp)
  real(dp), parameter :: sqrt_pi   = sqrt(dacos(-1.0_dp))
  real(dp), parameter :: sqrt2     = sqrt(2.0_dp)
  real(dp) :: tpi       = 2.0_dp*dacos(-1.0_dp)
  real(dp), parameter :: fpi       = 4.0_dp*dacos(-1.0_dp)
  real(dp), parameter :: eps_2     = 1.0e-2_dp
  real(dp) :: eps_5     = 1.0e-5_dp
  real(dp) :: eps_6     = 1.0e-6_dp
  real(dp), parameter :: eps_7     = 1.0e-7_dp
  real(dp), parameter :: eps_8     = 1.0e-8_dp
  real(dp), parameter :: eps_10    = 1.0e-10_dp
  real(dp), parameter :: eps_14    = 1.0e-14_dp
  real(dp), parameter :: bohr      = 0.52917720859_dp
  real(dp), parameter :: boltzmann = 8.6173324*0.00001_dp
  !
  complex(dp) :: cmplx_i = (0.0_dp,1.0_dp)
  complex(dp) :: cmplx_1 = (1.0_dp,0.0_dp)
  complex(dp) :: cmplx_0 = (0.0_dp,0.0_dp)
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
