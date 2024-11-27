!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
!
module intw_fft
  !
  !----------------------------------------------------------------------------!
  !
  !       The subroutines in this module handle most of the g -> r and
  !       r -> g tasks using FFT.
  !
  !----------------------------------------------------------------------------!
  !
  use kinds, only: dp
  !
  implicit none

  !
  ! variables
  public :: nl, g_fft_map, ig1, ig2, ig3, gvec_cart, gg, phase, eigts1, &
            eigts2, eigts3, strf, mill
  !
  ! subroutines
  public :: allocate_fft, deallocate_fft, &
            generate_nl, wfc_by_expigr, wfc_from_g_to_r, wfc_from_r_to_g, &
            func_from_g_to_r, r_function_by_exp_igr, func_from_r_to_g, &
            find_iG, coarse_to_smooth
  !
  private
  !


  ! correspondence between the nr1 nr2 nr3 FFT grid with the G list (calculations).
  integer, allocatable :: nl(:)
  integer, allocatable :: g_fft_map(:,:,:)
  integer, allocatable :: ig1(:), ig2(:), ig3(:)

  real(dp), allocatable :: gvec_cart(:,:), gg(:)

  !The complex phases for each G and atomic position: exp( -2 pi i tau(1:nat).G )
  complex(dp), allocatable :: phase(:,:)

  complex(dp), allocatable :: eigts1(:,:), eigts2(:,:), eigts3(:,:), strf(:,:)
  integer, allocatable :: mill(:,:)
  !
  save


contains

  subroutine allocate_fft()
    !--------------------------------------------------------
    ! This subroutine simply allocates the arrays needed
    ! by the fft algorithms.
    !--------------------------------------------------------
    use intw_reading, only: ngm, nr1, nr2, nr3, nat, ntyp

    implicit none

    allocate (nl(ngm), g_fft_map( nr1, nr2, nr3) )
    allocate (ig1(ngm), ig2(ngm), ig3(ngm))
    allocate (phase(ngm,nat) )
    allocate (gvec_cart(3,ngm))
    allocate (gg(ngm))
    allocate (eigts1(-nr1:nr1, nat))
    allocate (eigts2(-nr2:nr2, nat))
    allocate (eigts3(-nr3:nr3, nat))
    allocate (strf(ngm, ntyp))

    allocate (mill(3,ngm))

  end subroutine allocate_fft

  subroutine deallocate_fft()
    !--------------------------------------------------------
    ! This subroutine simply deallocates the arrays needed
    ! by the fft algorithms.
    !--------------------------------------------------------
    deallocate (nl, g_fft_map )
    deallocate (ig1, ig2, ig3)
    deallocate (phase )
    deallocate (gvec_cart)
    deallocate (gg)

    deallocate (eigts1)
    deallocate (eigts2)
    deallocate (eigts3)
    deallocate (strf  )

    deallocate (mill)

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
    !             - It creates an array nl(ng), which returns the scalar
    !               index "nl" corresponding to the triple index (n1,n2,n3)
    !               of the G vector gvec(ng).
    !
    !             - It computes the phases
    !                    phase(ng,na) = Exp ( -2 pi i gvec(ng) * tau(na) )
    !
    !------------------------------------------------------------------------
    use intw_reading, only: ngm, gvec, bg, nat, tau, ntyp, ityp, nr1, nr2, nr3
    use intw_utility, only: triple_to_joint_index_r, cryst_to_cart
    use intw_useful_constants, only: tpi, cmplx_i

    implicit none

    !local

    integer :: n1, n2, n3
    integer :: ng, na

    logical :: assigned(nr1,nr2,nr3)

    integer :: nt, ipol

    real(dp) :: arg, bgtau(3)
    !output :: nl and ig1 ig2 ig3 and phase

    g_fft_map(:,:,:) = 0
    nl(:) = 0
    assigned = .false.



    ! loop on all G vectors in the global array gvec
    do ng = 1, ngm
      ! find the triple index corresponding to the G_fft mesh
      ! NOTE: the function "modulo" always returns a positive number in FORTRAN90
      !       the function "mod" is more dangerous.

      !n1 = nint (sum(g (:, ng) * at (:, 1))) + 1
      !mill (1,ng) = n1 - 1

      n1 = modulo(gvec(1,ng), nr1)+1 !modulo(gvec(1,ng), nr1)+1
      ig1 (ng) = n1 - 1

      n2 = modulo(gvec(2,ng), nr2)+1
      ig2 (ng) = n2 - 1

      n3 = modulo(gvec(3,ng), nr3)+1
      ig3 (ng) = n3 - 1


      mill (1, ng) = gvec(1,ng) !ig1 (ng)
      mill (2, ng) = gvec(2,ng) !ig2 (ng)
      mill (3, ng) = gvec(3,ng) !ig3 (ng)

      if (.not. assigned(n1,n2,n3) ) then

        assigned(n1,n2,n3) = .true.
        g_fft_map(n1,n2,n3) = ng

        ! compute the scalar index corresponding to n1,n2,n3 and
        ! assign it to nl(ng)

        call triple_to_joint_index_r(nr1,nr2,nr3,nl(ng),n1,n2,n3)

      else
        write(*,*) 'ERROR in generate_nl. FFT mesh too small?'
        write(*,*) '    More than one G-vector in the gvec array are being'
        write(*,*) '    assigned to the same FFT triple (n1,n2,n3);       '
        write(*,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
        write(*,*) '    too small.                                        '

        stop
      endif
    end do

    !--------------------------------------------------
    ! I think the code that follows is not useful
    ! and should be removed (Bruno)
    !--------------------------------------------------

    ! Obtain the G vectors in cartesian coordinates.
    gvec_cart(1:3,1:ngm) = gvec(1:3,1:ngm)
    call cryst_to_cart (ngm, gvec_cart, bg, 1)

    do ng=1, ngm
      gg(ng) = (gvec_cart(1,ng)**2 + gvec_cart(2,ng)**2 + gvec_cart(3,ng)**2 )
    enddo

    ! Compute the phases
    do na=1,nat
      do ng=1,ngm
        ! the tau vectors are in cartesian, alat units.

        phase(ng,na) = Exp ( -tpi*cmplx_i* &
                          (  gvec_cart(1,ng)*tau(1,na) &
                           + gvec_cart(2,ng)*tau(2,na) &
                           + gvec_cart(3,ng)*tau(3,na) ) )
      enddo
    enddo

    strf(:,:) = (0.d0,0.d0)
    do nt = 1, ntyp
      do na = 1, nat
        if (ityp (na) .eq.nt) then
          do ng = 1, ngm
            arg = ( gvec_cart (1, ng) * tau (1, na) + &
                    gvec_cart (2, ng) * tau (2, na) + &
                    gvec_cart (3, ng) * tau (3, na) ) * tpi
            strf (ng, nt) = strf (ng, nt) + CMPLX(cos (arg), -sin (arg),kind=DP)
          enddo
        endif
      enddo
    enddo

    do na = 1, nat
      do ipol = 1, 3
        bgtau (ipol) = bg (1, ipol) * tau (1, na) + &
                       bg (2, ipol) * tau (2, na) + &
                       bg (3, ipol) * tau (3, na)
      enddo
      do n1 = - nr1, nr1
        arg = tpi * n1 * bgtau (1)
        eigts1 (n1, na) = CMPLX(cos (arg), - sin (arg),kind=DP)
      enddo
      do n2 = - nr2, nr2
        arg = tpi * n2 * bgtau (2)
        eigts2 (n2, na) = CMPLX(cos (arg), - sin (arg),kind=DP)
      enddo
      do n3 = - nr3, nr3
        arg = tpi * n3 * bgtau (3)
        eigts3 (n3, na) = CMPLX(cos (arg), - sin (arg),kind=DP)
      enddo
    enddo

  end subroutine generate_nl


  subroutine wfc_by_expigr(num_bands, nspin, G, list_iG, wfc)

    use intw_reading, only: gvec, nG_max
    use intw_utility, only: hpsort_integer
    use intw_useful_constants, only: cmplx_0

    implicit none

    !I/O variables

    integer, intent(in) :: num_bands, nspin
    integer, intent(in) :: G(3) ! G vector such that k_out = k_in + G
    integer, intent(inout) :: list_iG(nG_max) ! On input, G vector indices for k, sorted
                                              ! On output, G vector indices for k + G, sorted
    complex(dp), intent(inout) :: wfc(ng_max,num_bands,nspin) ! On input, wave function components for k
                                                              ! On output, wave function components for k + G

    !local variables

    integer :: list_iG_k_irr(nG_max)
    complex(dp) :: wfc_k_irr(ng_max,num_bands,nspin)
    integer :: p_i, i, iG_k_irr, iG_k
    integer :: G_k(3) ! a vector for Rk, the point in the 1BZ
    integer :: permutations(nG_max) ! index permutation which orders list_G_k
    integer :: nb, ispin, nG


    !Initialization
    !
    list_iG_k_irr = list_iG
    list_iG = 0
    !
    ! loop on all G_k_irr, the coefficients of the wavefunction at the IBZ k point
    !
    nG = 0
    do i = 1, nG_max
      !
      iG_k_irr = list_iG_k_irr(i)
      !
      if (iG_k_irr == 0) exit ! the index array is zero-padded at the end.
      !
      nG = nG+1
      !
      G_k(:) = gvec(:,iG_k_irr) - G(:) ! minus, zeren horrela da konbentzioa exp(-igr) (testatuta dago intw2wan).
      !
      call find_iG(G_k, iG_k)
      !
      list_iG(nG) = iG_k
      !
    enddo
    !
    call hpsort_integer(nG, list_iG, permutations)
    wfc_k_irr = wfc
    wfc = cmplx_0
    !
    do i = 1, nG
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
    !     in ::        wfc_g(nG_max)         : the u_{nk}(G) coefficients
    !                  list_iG               : the indices of the G vectors used
    !                                          in the wave function
    !
    !     out::  wfc_r(nr1*nr2*nr3)          : The periodic part of the wave functions
    !                                          in real space, with the space index
    !                                          represented by a scalar.
    !--------------------------------------------------------
    use intw_reading, only: nG_max, nr1, nr2, nr3
    use intw_useful_constants, only: cmplx_0

    implicit none

    external :: cfftnd

    integer :: i, iG
    integer :: list_iG(nG_max)

    complex(dp), intent(in) :: wfc_g(nG_max)
    complex(dp), intent(out) :: wfc_r(nr1*nr2*nr3)



    ! initialize work array
    wfc_r(:) = cmplx_0

    ! put wfc_g in wfc_r
    do i=1,nG_max
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
    !     in ::        wfc_g(nG_max)         : the u_{nk}(G) coefficients
    !                  list_iG               : the indices of the G vectors used
    !                                          in the wave function
    !
    !     out::  wfc_r(nr1*nr2*nr3)          : The periodic part of the wave functions
    !                                          in real space, with the space index
    !                                          represented by a scalar.
    !--------------------------------------------------------
    use intw_reading, only: nG_max, nr1, nr2, nr3
    use intw_useful_constants, only: cmplx_0

    implicit none

    external :: cfftnd

    integer :: i, iG
    integer :: list_iG(nG_max)

    complex(dp), intent(out) :: wfc_g(nG_max)
    complex(dp), intent(in) :: wfc_r(nr1*nr2*nr3)

    complex(dp) :: aux(nr1*nr2*nr3)


    ! initialize work array
    aux = wfc_r
    wfc_g(:) = cmplx_0

    call cfftnd(3,(/nr1,nr2,nr3/),-1,aux) !
                               ! this convention reproduces
                               ! the results of pw2wannier EXACTLY

  do i=1,nG_max
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
    use intw_reading, only: nr1, nr2, nr3, nspin, ngm
    use intw_useful_constants, only: cmplx_0

    implicit none

    integer, intent(in) :: nfunc
    complex(dp), intent(in) :: fg(ngm,nfunc)
    complex(dp), intent(out) :: fr(nr1*nr2*nr3,nfunc)

    integer :: mode, ig
    complex(dp) :: aux(nr1*nr2*nr3)

    external :: cfftnd


    do mode=1,nfunc
      ! initialize work array
      aux(:) = cmplx_0

      ! put fg in aux
      do ig=1,ngm
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
    use intw_reading, only: nr1, nr2, nr3, ngm
    use intw_useful_constants, only: cmplx_0

    implicit none

    integer, intent(in) :: nfunc
    complex(dp), intent(in) :: fr(nr1*nr2*nr3,nfunc)
    complex(dp), intent(out) :: fg(ngm,nfunc)

    integer :: mode, ig, ir
    complex(dp) :: aux(nr1*nr2*nr3)

    external :: cfftnd


    do mode=1,nfunc

      aux(:) = cmplx_0

      do ir=1,nr1*nr2*nr3
        aux(ir) = fr(ir,mode)
      enddo

      call cfftnd(3, (/nr1,nr2,nr3/), -1, aux) ! this is the right convention

      do ig=1,ngm
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
