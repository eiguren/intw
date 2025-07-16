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
module intw_fft

  !--------------------------------------------------------------!
  ! The subroutines in this module handle most of the g -> r and !
  ! r -> g tasks using FFT.                                      !
  !--------------------------------------------------------------!

  use kinds, only: dp
  !
  implicit none
  !
  ! variables
  public :: nl, g_fft_map, gvec_cart, gg, phase, strf
  !
  ! subroutines
  public :: allocate_fft, deallocate_fft, &
            generate_nl, wfc_by_expigr, wfc_from_g_to_r, wfc_from_r_to_g, &
            func_from_g_to_r, r_function_by_exp_igr, func_from_r_to_g, &
            find_iG, coarse_to_smooth
  !
  private
  !


  ! Correspondence between the nr1 nr2 nr3 FFT grid with the G list (calculations).
  integer, allocatable :: nl(:)
  integer, allocatable :: g_fft_map(:,:,:)
  ! The G vectors in cartesian coordinates
  real(dp), allocatable :: gvec_cart(:,:)
  ! The module squared of the G vectors
  real(dp), allocatable :: gg(:)

  ! The complex phases for each G and atomic position: exp( -2 pi i tau(1:nat).G )
  complex(dp), allocatable :: phase(:,:)
  ! The structure factor
  complex(dp), allocatable :: strf(:,:)
  !
  save


contains

  subroutine allocate_fft()
    !--------------------------------------------------------
    ! This subroutine simply allocates the arrays needed
    ! by the fft algorithms.
    !--------------------------------------------------------
    use intw_reading, only: nG, nr1, nr2, nr3, nat, ntyp

    implicit none

    allocate (nl(nG), g_fft_map( nr1, nr2, nr3) )
    allocate (phase(nG,nat) )
    allocate (gvec_cart(3,nG))
    allocate (gg(nG))
    allocate (strf(nG, ntyp))

  end subroutine allocate_fft

  subroutine deallocate_fft()
    !--------------------------------------------------------
    ! This subroutine simply deallocates the arrays needed
    ! by the fft algorithms.
    !--------------------------------------------------------
    deallocate (nl, g_fft_map)
    deallocate (phase)
    deallocate (gvec_cart)
    deallocate (gg)
    deallocate (strf)

  end subroutine deallocate_fft

  subroutine generate_nl()
    !------------------------------------------------------------------------
    !  This subroutine generates information important for the 3D-FFT
    !  algorithm.
    !
    !  A bit of theory:
    !  ---------------
    !
    !  Consider a mesh in real space described by nr1 x nr2 x nr3, such that
    !            r(i,j,k) =  (i-1)/nr1 a1 +  (j-1)/nr2 a2 + (k-1)/nr3 a3
    !
    !  with a1, a2, a3 the basis vectors of the real space lattice, and
    !                1 <= i <= nr1, 1 <= j <= nr2, 1 <= k <= nr3.
    !
    !  A reciprocal space lattice can then be defined:
    !
    !  G_fft(I,J,K)  =  (I-1) b1 +  (J-1) b2 + (K-1) b3
    !  with b1, b2, b3 the basis vectors of the reciprocal space lattice.
    !             As usual, a_i * b_j = 2 pi delta_{ij}
    !
    !  So that the phase factor is given by:
    !
    !     EXP [ i G_fft(I,J,K) * r(i,j,k) ] =
    !
    !             EXP[ 2 pi i (  (I-1)(i-1)   + (J-1)(j-1)   + (K-1)(k-1) )]
    !                [        (  ---------      ----------     ---------  )]
    !                [        (     nr1            nr2             nr3    )]
    !
    !     Then, clearly, the only independent phase factors correspond to
    !                1 <= I <= nr1, 1 <= J <= nr2, 1 <= K <= nr3.
    !     Any G vector outside of this box is periodically equivalent to
    !     a G vector inside the box.
    !
    !  So, what does this subroutine do?
    !  ---------------------------------
    !
    !  The G vectors usually used are not defined in the box
    !                1 <= I <= nr1, 1 <= J <= nr2, 1 <= K <= nr3,
    !  but rather they are symmetrically distributed around the origin.
    !  This makes sense from a physical point of view, which doesn't care
    !  about FFT meshes.
    !
    !  THIS SUBROUTINE:
    !
    !             - finds the relationship between the usually defined
    !               G vectors, stored in gvec in crystal coordinates, and
    !               the FFT G vector mesh.
    !
    !             - It creates an array
    !                  G_fft_map(n1,n2,n3) = index of G vector in gvec which
    !                                        corresponds to G_fft(n1,n2,n3)
    !
    !             - It creates an array nl, which returns the scalar
    !               index "nl" corresponding to the triple index (n1,n2,n3)
    !               of the G vector gvec(ig).
    !
    !             - It computes the phases
    !                    phase(ig,ia) = Exp ( -2 pi i gvec(ig) * tau(ia) )
    !
    !------------------------------------------------------------------------
    use intw_reading, only: nG, gvec, bg, nat, tau, ntyp, ityp, nr1, nr2, nr3, tpiba2
    use intw_utility, only: triple_to_joint_index_r, cryst_to_cart
    use intw_useful_constants, only: tpi, cmplx_i, cmplx_0

    implicit none

    !local

    logical :: assigned(nr1,nr2,nr3)
    integer :: n1, n2, n3
    integer :: ig, na, nt


    g_fft_map(:,:,:) = 0
    nl(:) = 0
    assigned(:,:,:) = .false.
    ! loop on all G vectors in the global array gvec
    do ig = 1, nG
      ! find the triple index corresponding to the G_fft mesh
      ! NOTE: the function "modulo" always returns a positive number in FORTRAN90
      !       the function "mod" is more dangerous.

      n1 = modulo(gvec(1,ig), nr1)+1

      n2 = modulo(gvec(2,ig), nr2)+1

      n3 = modulo(gvec(3,ig), nr3)+1

      if (.not. assigned(n1,n2,n3) ) then

        assigned(n1,n2,n3) = .true.
        g_fft_map(n1,n2,n3) = ig

        ! compute the scalar index corresponding to n1,n2,n3 and
        ! assign it to nl(ig)

        call triple_to_joint_index_r(nr1,nr2,nr3,nl(ig),n1,n2,n3)

      else
        write(*,*) 'ERROR in generate_nl. FFT mesh too small?'
        write(*,*) '    More than one G-vector in the gvec array are being'
        write(*,*) '    assigned to the same FFT triple (n1,n2,n3);       '
        write(*,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
        write(*,*) '    too small.                                        '

        stop
      endif

    end do ! ig

    ! Obtain the G vectors in cartesian coordinates
    gvec_cart(1:3,1:nG) = gvec(1:3,1:nG)
    call cryst_to_cart (nG, gvec_cart, bg, 1)

    ! Calculate module squared of the G vectors
    do ig=1, nG
      gg(ig) = tpiba2*dot_product(gvec_cart(:,ig), gvec_cart(:,ig))
    enddo

    ! Compute the phases
    do na=1,nat
      do ig=1,nG
        ! the tau vectors are in cartesian, alat units.
        phase(ig,na) = exp ( -tpi*cmplx_i*dot_product(gvec_cart(:,ig), tau(:,na)) )
      enddo
    enddo

    ! Compute structure factor
    strf(:,:) = cmplx_0
    do nt = 1, ntyp
      do na = 1, nat
        if ( ityp(na) == nt ) then
          strf(:,nt) = strf(:,nt) + phase(:,na)
        endif
      enddo
    enddo

  end subroutine generate_nl


  subroutine wfc_by_expigr(num_bands, nspin, G, list_iG, wfc)

    use intw_reading, only: gvec, nGk_max
    use intw_utility, only: hpsort_integer
    use intw_useful_constants, only: cmplx_0

    implicit none

    !I/O variables

    integer, intent(in) :: num_bands, nspin
    integer, intent(in) :: G(3) ! G vector such that k_out = k_in + G
    integer, intent(inout) :: list_iG(nGk_max) ! On input, G vector indices for k, sorted
                                              ! On output, G vector indices for k + G, sorted
    complex(dp), intent(inout) :: wfc(nGk_max,num_bands,nspin) ! On input, wave function components for k
                                                              ! On output, wave function components for k + G

    !local variables

    integer :: list_iG_k_irr(nGk_max)
    complex(dp) :: wfc_k_irr(nGk_max,num_bands,nspin)
    integer :: p_i, i, iG_k_irr, iG_k
    integer :: G_k(3) ! a vector for Rk, the point in the 1BZ
    integer :: permutations(nGk_max) ! index permutation which orders list_G_k
    integer :: nb, ispin, nGk


    !Initialization
    !
    list_iG_k_irr = list_iG
    list_iG = 0
    !
    ! loop on all G_k_irr, the coefficients of the wavefunction at the IBZ k point
    !
    nGk = 0
    do i = 1, nGk_max
      !
      iG_k_irr = list_iG_k_irr(i)
      !
      if (iG_k_irr == 0) exit ! the index array is zero-padded at the end.
      !
      nGk = nGk + 1
      !
      G_k(:) = gvec(:,iG_k_irr) - G(:) ! minus, zeren horrela da konbentzioa exp(-igr) (testatuta dago intw2wan).
      !
      call find_iG(G_k, iG_k)
      !
      list_iG(nGk) = iG_k
      !
    enddo
    !
    call hpsort_integer(nGk, list_iG, permutations)
    wfc_k_irr = wfc
    wfc = cmplx_0
    !
    do i = 1, nGk
      !
      p_i = permutations(i)
      !
      ! compute the wfc element
      !
      do nb = 1, num_bands
        do ispin = 1, nspin
          !
          wfc(i,nb,ispin) = wfc_k_irr(p_i,nb,ispin)
          !
        enddo
      enddo
      !
    enddo

  end subroutine wfc_by_expigr


  subroutine wfc_from_g_to_r (list_iG,wfc_g, wfc_r)
    !--------------------------------------------------------
    !  This subroutine is a driver which uses the 3D-FFT
    !  code to transform a wavefunction in G space to
    !  a wavefunction in r space.
    !
    !     in ::        wfc_g(nGk_max)        : the u_{nk}(G) coefficients
    !                  list_iG               : the indices of the G vectors used
    !                                          in the wave function
    !
    !     out::  wfc_r(nr1*nr2*nr3)          : The periodic part of the wave functions
    !                                          in real space, with the space index
    !                                          represented by a scalar.
    !--------------------------------------------------------
    use intw_reading, only: nGk_max, nr1, nr2, nr3
    use intw_useful_constants, only: cmplx_0

    implicit none

    external :: cfftnd

    integer :: i, iG
    integer :: list_iG(nGk_max)

    complex(dp), intent(in) :: wfc_g(nGk_max)
    complex(dp), intent(out) :: wfc_r(nr1*nr2*nr3)



    ! initialize work array
    wfc_r(:) = cmplx_0

    ! put wfc_g in wfc_r
    do i=1,nGk_max
      ! identify the G vector by its index, as stored in list_iG
      iG = list_iG(i)

      if (iG == 0) exit
      ! use nl to identify which G_fft vector G corresponds to,
      ! and assign the value of the wave function in the aux array
      wfc_r(nl(iG)) = wfc_g(i)
    enddo

    ! perform fourier transform in place wfc_g(G) -> wfc_r(r)
    ! CONVENTION BY ASIER
     call cfftnd(3,(/nr1,nr2,nr3/),1,wfc_r) !
                               ! this convention reproduces
                               ! the results of pw2wannier EXACTLY

  end subroutine wfc_from_g_to_r

  subroutine wfc_from_r_to_g (list_iG,wfc_r, wfc_g)
    !--------------------------------------------------------
    !  This subroutine is a driver which uses the 3D-FFT
    !  code to transform a wavefunction in G space to
    !  a wavefunction in r space.
    !
    !     in ::        wfc_g(nGk_max)        : the u_{nk}(G) coefficients
    !                  list_iG               : the indices of the G vectors used
    !                                          in the wave function
    !
    !     out::  wfc_r(nr1*nr2*nr3)          : The periodic part of the wave functions
    !                                          in real space, with the space index
    !                                          represented by a scalar.
    !--------------------------------------------------------
    use intw_reading, only: nGk_max, nr1, nr2, nr3
    use intw_useful_constants, only: cmplx_0

    implicit none

    external :: cfftnd

    integer :: i, iG
    integer :: list_iG(nGk_max)

    complex(dp), intent(out) :: wfc_g(nGk_max)
    complex(dp), intent(in) :: wfc_r(nr1*nr2*nr3)

    complex(dp) :: aux(nr1*nr2*nr3)


    ! initialize work array
    aux = wfc_r
    wfc_g(:) = cmplx_0

    call cfftnd(3,(/nr1,nr2,nr3/),-1,aux) !
                               ! this convention reproduces
                               ! the results of pw2wannier EXACTLY

  do i=1,nGk_max
    ! identify the G vector by its index, as stored in list_iG
    iG = list_iG(i)
    if (iG == 0) exit
    ! use nl to identify which G_fft vector G corresponds to,
    ! and assign the value of the wave function in the aux array
    wfc_g(i) = aux(nl(ig))
  enddo

  end subroutine wfc_from_r_to_g

  subroutine func_from_g_to_r (nfunc, fg, fr)
    !--------------------------------------------------------
    !  This subroutine is a driver which uses the 3D-FFT
    !  code to go from f(G) to f(r).
    !--------------------------------------------------------
    use intw_reading, only: nr1, nr2, nr3, nG
    use intw_useful_constants, only: cmplx_0

    implicit none

    integer, intent(in) :: nfunc
    complex(dp), intent(in) :: fg(nG,nfunc)
    complex(dp), intent(out) :: fr(nr1*nr2*nr3,nfunc)

    integer :: mode, ig
    complex(dp) :: aux(nr1*nr2*nr3)

    external :: cfftnd


    do mode=1,nfunc
      ! initialize work array
      aux(:) = cmplx_0

      ! put fg in aux
      do ig=1,nG
        aux(nl(ig)) = fg(ig, mode)
      enddo

      ! perform fourier transform in place aux(fg) -> aux(fr)
      call cfftnd(3, (/nr1,nr2,nr3/), 1, aux) ! This is the right convention

      ! put aux in fr
      fr(:,mode) = aux(:)

    enddo

  end subroutine func_from_g_to_r


  subroutine r_function_by_exp_igr (g_cryst, nfunc, nr1, nr2, nr3, fr, fr_exp_igr)

    use intw_utility, only: joint_to_triple_index_g
    use intw_useful_constants, only: tpi, cmplx_i

    implicit none

    integer,intent(in) :: g_cryst(3), nfunc, nr1, nr2, nr3
    complex(dp),intent(in) :: fr(nr1*nr2*nr3,nfunc)
    complex(dp),intent(out) :: fr_exp_igr(nr1*nr2*nr3,nfunc)

    integer :: mode, i, j, k, ir
    real(dp) :: gr


    do mode=1,nfunc

      do ir=1,nr1*nr2*nr3

        call joint_to_triple_index_g(nr1,nr2,nr3,ir,i,j,k)

        gr = tpi*(g_cryst(1)*(i-1)/nr1 + g_cryst(1)*(j-1)/nr2 + g_cryst(1)*(k-1)/nr3)

        fr_exp_igr(ir,mode) = fr(ir,mode) * exp( cmplx_i * gr )

      enddo
    enddo

  end subroutine r_function_by_exp_igr


  subroutine func_from_r_to_g (nfunc, fr, fg)
    !--------------------------------------------------------
    !  This subroutine is a driver which uses the 3D-FFT
    !  code to go from f(r) to f(G).
    !--------------------------------------------------------
    use intw_reading, only: nr1, nr2, nr3, nG
    use intw_useful_constants, only: cmplx_0

    implicit none

    integer, intent(in) :: nfunc
    complex(dp), intent(in) :: fr(nr1*nr2*nr3,nfunc)
    complex(dp), intent(out) :: fg(nG,nfunc)

    integer :: mode, ig, ir
    complex(dp) :: aux(nr1*nr2*nr3)

    external :: cfftnd


    do mode=1,nfunc

      aux(:) = cmplx_0

      do ir=1,nr1*nr2*nr3
        aux(ir) = fr(ir,mode)
      enddo

      call cfftnd(3, (/nr1,nr2,nr3/), -1, aux) ! this is the right convention

      do ig=1,nG
        fg(ig,mode) = aux(nl(ig))
      enddo
    enddo

  end subroutine func_from_r_to_g

  subroutine find_iG(G,iG)
    !----------------------------------------------------------------------------!
    !     Given a G vector in crystal coordinates, this subroutine
    !     finds the index iG to which this G vector corresponds to in the
    !     global list of G vectors. In other words:
    !
    !     G = gvec(iG)
    !
    !     If G is not found, an error is thrown.
    !----------------------------------------------------------------------------!
    use intw_reading, only: nr1, nr2, nr3

    implicit none

    integer :: G(3)
    integer :: iG, n1, n2, n3

    ! Find the FFT G vector corresponding to G
    n1 = modulo(G(1), nr1)+1
    n2 = modulo(G(2), nr2)+1
    n3 = modulo(G(3), nr3)+1

    ! use the tabulated values to find iG
    iG = g_fft_map(n1,n2,n3)

  end subroutine find_iG


  subroutine coarse_to_smooth(n1,n2,n3,FR_coarse, n1s,n2s,n3s,FR_smooth)
    !-------------------------------------------
    ! This subroutine puts the coarse force
    ! constants in the smooth force constant
    ! array
    !-------------------------------------------
    use intw_useful_constants, only: cmplx_0
    implicit none

    ! input/output variables
    complex(dp) :: FR_coarse(n1,n2,n3)
    complex(dp) :: FR_smooth(n1s,n2s,n3s)

    integer :: n1, n2, n3, n1s, n2s, n3s

    ! local variables
    integer :: i1, i2, i3, i1s, i2s, i3s

    ! local variables

    FR_smooth(:,:,:) = cmplx_0

    do i1 = 1, n1
      if (i1 <= n1/2 ) then
        i1s = i1
      else
        i1s = n1s-n1+i1
      end if

      do i2 = 1, n2
        if (i2 <= n2/2 ) then
          i2s = i2
        else
          i2s = n2s-n2+i2
        end if

        do i3 = 1, n3
          if (i3 <= n3/2 ) then
            i3s = i3
          else
            i3s = n3s-n3+i3
          end if

          FR_smooth(i1s,i2s,i3s) = FR_coarse(i1,i2,i3)

        end do
      end do
    end do

 end subroutine coarse_to_smooth


!----------------------------------------------------------------------------!
!
!
  end module intw_fft
!
!
!----------------------------------------------------------------------------!
