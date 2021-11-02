module intw_useful_constants 
!       The name is "useful_constants" to avoid clashes with W90 "constants.f90"
!
use kinds,     only : dp
  !
  implicit none
  !
!haritz
  public :: ZERO, ONE,sqrt2,e2, cmplx_i, cmplx_1, cmplx_0, pi, sqrt_pi, &
            twopi, tpi, fpi, eps_2, eps_5, eps_6, eps_7, eps_8, &
            eps_10, eps_14, bohr, Ha_in_eV, boltzmann, direct_io_factor, &
            direct_io_factor_cmplx, double_complex, double_real, &
            mpi_rank, mpi_nproc, mpi_err, I2, sig_x, sig_y, sig_z

  private
!haritz
  !
  save
        !

        real(dp)    :: ZERO        = 0.0_dp
        real(dp)    :: ONE         = 1.0_dp
        real(dp)    :: e2         = 2.0_dp

        complex(dp) :: cmplx_i     = (0.0_dp,1.0_dp)
        complex(dp) :: cmplx_1     = (1.0_dp,0.0_dp)
        complex(dp) :: cmplx_0     = (0.0_dp,0.0_dp)
        real(dp)    :: pi          =  dacos(-1.0_dp)
        real(dp)    :: sqrt_pi     =  sqrt(dacos(-1.0_dp))
        real(dp)    :: sqrt2     =  sqrt(2.0_dp)
        real(dp)    :: twopi       =  2_dp*dacos(-1.0_dp)
        real(dp)    :: tpi         =  2.0_dp*dacos(-1.0_dp)
        real(dp)    :: fpi         =  4.0_dp*dacos(-1.0_dp)
        real(dp)    :: eps_2       =  1.0e-2_dp
        real(dp)    :: eps_5       =  1.0e-5_dp
        real(dp)    :: eps_6       =  1.0e-6_dp
        real(dp)    :: eps_7       =  1.0e-7_dp
        real(dp)    :: eps_8       =  1.0e-8_dp
        real(dp)    :: eps_10      =  1.0e-10_dp
        real(dp)    :: eps_14      =  1.0e-14_dp
        real(dp)    :: bohr        =  0.52917720859_dp 
        !
        real(dp)    :: Ha_in_eV    = 27.211383860484776_dp
        !
        real(dp)    :: boltzmann   =  8.6173324*0.00001_dp

        real(dp)    :: AUTOEV = 27.211383860484776_dp

        !careful with this: it might be system / compiler dependent...


        ! the variables direct_io_* are appropriate to 
        ! specify the length of a record when the memory UNIT used by
        ! fortran to write unformatted files is a "long word", which
        ! is 4 bytes (I guess a "word" would be 1 byte=8 bit ?).

        ! this is for a real double
        integer,parameter     :: direct_io_factor = 2 
        ! this is for a complex double
        integer,parameter     :: direct_io_factor_cmplx = 4 

        ! what is the size in bytes of a complex double
        integer,parameter     :: double_complex = 16

        ! what is the size in bytes of a complex real
        integer,parameter     :: double_real    = 8

        ! MPI rank and procs
        integer  :: mpi_rank, mpi_nproc, mpi_err

        !-Pauli matrices
        complex(kind=dp) :: I2(2,2), sig_x(2,2), sig_y(2,2), sig_z(2,2)


end module  intw_useful_constants
