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
module intw_useful_constants
  !
  use kinds, only: dp
  !
  implicit none
  !
  ! I/O variables
  public :: stdout, stdin
  !
  ! Useful numbers
  public :: ZERO, ONE, TWO, sqrt2, pi, tpi, fpi, sqrt_pi, &
            eps_2, eps_5, eps_6, eps_7, eps_8, eps_10, eps_14, &
            cmplx_i, cmplx_1, cmplx_0
  !
  ! Physical constants
  public :: bohr, boltzmann
  !
  ! Energy conversion factors
  public :: eV_to_Ha, eV_to_Ry, Ha_to_eV, Ha_to_Ry, Ry_to_eV, Ry_to_Ha
  !
  ! Mass conversion factors
  public :: pmass
  !
  ! Pauli matrices
  public :: I2, sig_x, sig_y, sig_z
  !
  private
  !
  save
  !
  !
  ! I/O variables
  integer, parameter :: stdout = 6 ! standard output unit
  integer, parameter :: stdin = 5 ! standard input unit
  !
  ! Useful numbers
  real(kind=dp), parameter :: ZERO       = 0.0_dp
  real(kind=dp), parameter :: ONE        = 1.0_dp
  real(kind=dp), parameter :: TWO        = 2.0_dp
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
  complex(kind=dp), parameter :: cmplx_i = (0.0_dp, 1.0_dp)
  complex(kind=dp), parameter :: cmplx_1 = (1.0_dp, 0.0_dp)
  complex(kind=dp), parameter :: cmplx_0 = (0.0_dp, 0.0_dp)
  !
  ! Physical constants
  real(kind=dp), parameter :: bohr       = 0.52917720859_dp
  real(kind=dp), parameter :: boltzmann  = 8.6173324*0.00001_dp
  !
  ! Energy conversion factors
  real(kind=dp), parameter :: Ha_to_eV = 27.211383860484776_dp
  real(kind=dp), parameter :: eV_to_Ha = 1.0_dp/Ha_to_eV
  real(kind=dp), parameter :: Ha_to_Ry = 2.0_dp
  real(kind=dp), parameter :: Ry_to_Ha = 1.0_dp/Ha_to_Ry
  real(kind=dp), parameter :: Ry_to_eV = Ry_to_Ha*Ha_to_eV
  real(kind=dp), parameter :: eV_to_Ry = 1.0_dp/Ry_to_eV
  !
  ! Mass conversion factors
  real(kind=dp), parameter :: pmass = 1822.88848426_dp ! 1 Dalton in atomic units
  !
  ! Pauli matrices
  complex(kind=dp), dimension(2,2), parameter ::    I2 = transpose(reshape((/ cmplx_1, cmplx_0, &
                                                                              cmplx_0, cmplx_1 /), (/2, 2/)))
  complex(kind=dp), dimension(2,2), parameter :: sig_x = transpose(reshape((/ cmplx_0, cmplx_1, &
                                                                              cmplx_1, cmplx_0 /), (/2, 2/)))
  complex(kind=dp), dimension(2,2), parameter :: sig_y = transpose(reshape((/ cmplx_0,-cmplx_i, &
                                                                              cmplx_i, cmplx_0 /), (/2, 2/)))
  complex(kind=dp), dimension(2,2), parameter :: sig_z = transpose(reshape((/ cmplx_1, cmplx_0, &
                                                                              cmplx_0,-cmplx_1 /), (/2, 2/)))

end module  intw_useful_constants
