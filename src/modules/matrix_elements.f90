!----------------------------------------------------------------------------!
!	intw project.
!----------------------------------------------------------------------------!
module intw_matrix_elements

  use kinds, only: dp

  implicit none
  !
  ! subroutines
  public :: get_plane_wave_matrix_element_convolution_map, &
            get_plane_wave_matrix_element_FFT, &
            compute_index_interpolation_mesh, &
            get_spin_component, &
            write_matrix_elements, &
            wfc_G_from_1D_to_3D
  !
  private
  !

contains

  subroutine get_plane_wave_matrix_element_convolution_map(G, list_iG_1, ngk1, list_iG_2, ngk2, wfc_1, wfc_2, pw_mat_el)
    !--------------------------------------------------------------------------------------------------
    ! Given two wavefunctions wfc_1 and wfc_2, this subroutine computes
    !
    !              < wfc_1 | e^{-i G r} | wfc_2 >
    !
    ! where wfc is the periodic part of the wavefunction:
    !              u(r) = sum_G' e^{i G' r} wfc(G');
    !
    ! which leads to
    !              < wfc_1 | e^{-i G r} | wfc_2 > = sum_G1 wfc_1(G1)^* wfc_2(G1+G).
    !
    ! The computation is done over all bands.
    !
    ! The G-vectors are referenced by their indices in list_iG1, list_iG2,
    ! which refer to the global list gvec(3,ngm), which should be defined
    ! BEFORE using this subroutine.
    !--------------------------------------------------------------------------------

    use intw_useful_constants, only: zero, one, cmplx_0
    use intw_reading, only: nG_max, gvec, nspin, num_bands_intw
    use intw_fft, only : find_iG

    implicit none

    ! I/O variables

    integer, intent(in) :: G(3), ngk1, ngk2
    integer, intent(in) :: list_iG_1(nG_max), list_iG_2(nG_max)
    complex(dp), intent(in) :: wfc_1(nG_max,num_bands_intw,nspin), wfc_2(nG_max,num_bands_intw,nspin)
    complex(dp), intent(out) :: pw_mat_el(num_bands_intw,num_bands_intw,nspin,nspin)

    ! local variables

    integer :: jd(1)
    integer :: i, ibnd, jbnd, is, js, iG_1
    complex(dp) :: pw_mat_el_local(num_bands_intw,num_bands_intw,nspin,nspin)


    pw_mat_el = cmplx_0

    !$omp parallel default(none) &
    !$omp shared(list_iG_1,list_iG_2,ngk1,ngk2) &
    !$omp shared(num_bands_intw,nspin,wfc_1,wfc_2, pw_mat_el,nG_max) &
    !$omp shared(cmplx_0) &
    !$omp private(i,jd,ibnd,jbnd,is,js) &
    !$omp private( pw_mat_el_local)
    !
    ! First, build the pw_mat_el_local arrays, on each thread.
    !
    pw_mat_el_local = cmplx_0
    !
    !$omp do
    do i=1,nGk1
      !
      call find_iG(gvec(:,list_iG_1(i)) + G, iG_1)
      !
#if __GFORTRAN__ & __GNUC__ < 9
      jd(1) = 0
      do j=1,nG_max
        if (list_iG_2(j)==iG_1) then
          jd(1) = j
          exit
        endif
      enddo
#else
      jd = findloc(list_iG_2, iG_1)
#endif

      if (jd(1)==0) cycle
      !
      do js=1,nspin
        do is=1,nspin
          do jbnd=1,num_bands_intw
            do ibnd=1,num_bands_intw
              !
              pw_mat_el_local(ibnd,jbnd,is,js) = pw_mat_el_local(ibnd,jbnd,is,js) &
                                                 + conjg(wfc_1(i,ibnd,is)) * wfc_2(jd(1),jbnd,js)
              !
            enddo !ibnd
          enddo !jbnd
        enddo !is
      enddo !js
      !
    enddo ! i loop
    !$omp end do
    !
    !$omp barrier
    ! Next, dump pw_mat_el_local into pw_mat_el;
    ! each thread should dump all its components here, so don't
    ! parallelize over the loop variables! We must still be in the
    ! "parallel" environment, however, or else pw_mat_el_local becomes
    ! undefined
    !
    do js=1,nspin
      do is=1,nspin
        do jbnd=1,num_bands_intw
          do ibnd=1,num_bands_intw
            !
            ! make sure the update is atomic!
            !$omp atomic
            !
            pw_mat_el(ibnd,jbnd,is,js) = pw_mat_el(ibnd,jbnd,is,js) &
                                         + pw_mat_el_local(ibnd,jbnd,is,js)
            !
          enddo !ibnd
        enddo !jbnd
      enddo !is
    enddo !js
    !
    !$omp end parallel

  end subroutine get_plane_wave_matrix_element_convolution_map


  subroutine get_plane_wave_matrix_element_FFT(G, list_iG_1, list_iG_2, wfc_1, wfc_2, pw_mat_el)
    !------------------------------------------------------------------------
    ! Given two wavefunctions wfc_1 and wfc_2, this subroutine computes
    !
    !              < wfc_1 | e^{-i G r} | wfc_2 >
    !
    ! where wfc is the periodic part of the wavefunction:
    !              u(r) = sum_G' e^{i G' r} wfc(G');
    !
    ! which leads to
    !              < wfc_1 | e^{-i G r} | wfc_2 > = sum_G1 wfc_1(G1)^* wfc_2(G1+G)
    !
    ! The computation is done over all bands.
    !
    ! The G-vectors are referenced by their indices in list_iG1, list_iG2,
    ! which refer to the global list gvec(3,ngm), which should be defined
    ! BEFORE using this subroutine.
    !
    ! The matrix element is obtained using the FFT routines.
    !-------------------------------------------------------------------------

    use intw_fft, only: wfc_from_g_to_r, nl, find_iG
    use intw_reading, only: nr1, nr2, nr3, nG_max, nspin, num_bands_intw
    use intw_useful_constants, only: zero, one, cmplx_0

    implicit none

    external :: cfftnd

    !I/O variables

    integer, intent(in) :: G(3)
    integer, intent(in) :: list_iG_1(nG_max), list_iG_2(nG_max)
    complex(dp), intent(in) :: wfc_1(nG_max,num_bands_intw,nspin), wfc_2(nG_max,num_bands_intw,nspin)
    complex(dp), intent(inout) :: pw_mat_el(num_bands_intw,num_bands_intw,nspin,nspin)

    !local variables

    integer :: iG, iG_fft, ir
    integer :: ibnd, jbnd, is, js
    complex(dp) :: wfc_r1(nr1*nr2*nr3), wfc_r2(nr1*nr2*nr3)
    complex(dp) :: uu(nr1*nr2*nr3)


    pw_mat_el = cmplx_0
    !
    ! find the index of G in the global list gvec
    !
    call find_iG(G,iG)
    !
    ! find its scalar FFT index
    !
    iG_fft = nl(iG)
    !
    !loop on spin (diagonal is, js)
    !
    do is=1,nspin
      do js=1,nspin
        !
        ! loop on bands
        !
        do ibnd=1,num_bands_intw
          !
          ! fourier- transform the 1st wavefunction
          !
          call wfc_from_g_to_r(list_iG_1, wfc_1(:,ibnd,is), wfc_r1)
          !
          do jbnd=1,num_bands_intw
            !
            ! fourier- transform the 2nd wavefunction
            !
            call wfc_from_g_to_r(list_iG_2, wfc_2(:,jbnd,js), wfc_r2)
            !
            ! compute the product in real space
            !
            do ir=1,nr1*nr2*nr3
              !
              uu(ir) = conjg(wfc_r1(ir))*wfc_r2(ir)
              !
            enddo !ir
            !
            ! FFT to G in place
            ! call cfftnd(3,nr,1,uu) ! CONVENTION BY ASIER
            !
            call cfftnd(3, (/nr1,nr2,nr3/), -1, uu) ! OPPOSITE CONVENTION
                                                    ! this convention reproduces
                                                    ! the results of pw2wannier EXACTLY
            !
            ! extract the desired component
            !
            pw_mat_el(ibnd,jbnd,is,js) = pw_mat_el(ibnd,jbnd,is,js) + uu(iG_fft) ! Sum is over the spin.
            !
          enddo !jbnd
        enddo !ibnd
      enddo !js
    enddo !is

  end subroutine get_plane_wave_matrix_element_FFT


  subroutine compute_index_interpolation_mesh(iqpt, list_ikpt1, list_G1, list_ikpt2, list_G2)
    !----------------------------------------------------------------------------!
    ! This routine will be useful in computing matrix elements of the form:
    !              < psi{n1 k} | e^{-i (G+q) r} | psi_{n2 k+q} >.
    !
    ! Assume:
    !              k   = k1 + G1
    !              k+q = k2 + G2
    !
    ! where k1,k2 are in the mesh describing the 1BZ.
    !
    ! This subroutine, given the index of a qpoint, computes G1, G2, k1, k2
    ! for all k in a mesh appropriate for cubic interpolation.
    !--------------------------------------------------------------------------------

    use intw_input_parameters, only: nk1, nk2, nk3
    use intw_useful_constants, only: zero, one
    use intw_utility, only: triple_to_joint_index_g, joint_to_triple_index_g

    implicit none

    ! triple indices
    integer :: i_k,  j_k,  k_k  ! triple indices for k
    integer :: i_k1, j_k1, k_k1 ! triple indices for k1
    integer :: i_k2, j_k2, k_k2 ! triple indices for k2
    integer :: i_kq, j_kq, k_kq ! triple indices for k+q
    integer :: i_q,  j_q,  k_q  ! triple indices for q

    ! scalar indices
    integer :: ikpt1 ! scalar index for k1
    integer :: ikpt2 ! scalar index for k2
    integer :: iqpt  ! scalar index for q

    integer :: icm ! i cubic mesh: loop index over the extended mesh

    integer :: G1(3)
    integer :: G2(3)

    integer :: list_G1(3,(nk1+3)*(nk2+3)*(nk3+3))
    integer :: list_G2(3,(nk1+3)*(nk2+3)*(nk3+3))
    integer :: list_ikpt1((nk1+3)*(nk2+3)*(nk3+3))
    integer :: list_ikpt2((nk1+3)*(nk2+3)*(nk3+3))


    ! find the triple of indices which corresponds to iqpt
    call joint_to_triple_index_g(nk1,nk2,nk3,iqpt,i_q,j_q,k_q)

    ! loop over the mesh points for a cubic interpolation,
    ! including the extra layers.
    ! the third index, k, loops fastest!

    ! loop over the extended mesh, indentifying G1,G2, ikpt1 ikpt2
    icm = 0

    do i_k=0,nk1+2

      G1(1) = floor(1.0_dp*(i_k-1)/nk1)
      i_k1  = i_k-nk1*G1(1)
      i_kq  = i_k+i_q
      G2(1) = floor(1.0_dp*(i_kq-1)/nk1)
      i_k2  = i_kq-nk1*G2(1)


      do j_k=0,nk2+2
        G1(2) = floor(1.0_dp*(j_k-1)/nk2)
        j_k1  = j_k-nk2*G1(2)
        j_kq  = j_k+j_q
        G2(2) = floor(1.0_dp*(j_kq-1)/nk2)
        j_k2  = j_kq-nk2*G2(2)

        do k_k=0,nk3+2
          G1(3) = floor(1.0_dp*(k_k-1)/nk3)
          k_k1  = k_k-nk3*G1(3)
          k_kq  = k_k+k_q
          G2(3) = floor(1.0_dp*(k_kq-1)/nk3)
          k_k2  = k_kq-nk3*G2(3)

          ! find the indices
          call triple_to_joint_index_g(nk1,nk2,nk3,ikpt1,i_k1,j_k1,k_k1)
          call triple_to_joint_index_g(nk1,nk2,nk3,ikpt2,i_k2,j_k2,k_k2)

          icm = icm + 1 ! increment the mesh index

          ! save
          list_ikpt1(icm) = ikpt1
          list_ikpt2(icm) = ikpt2

          list_G1(:,icm) = G1
          list_G2(:,icm) = G2

        end do
      end do
    end do

  end subroutine compute_index_interpolation_mesh


  subroutine get_spin_component(wfc, spin)
    !-------------------------------------------------------------------------------
    ! Given wavefunction wfc, this subroutine computes
    !
    !              < wfc | sigma_alpha | wfc >
    !
    ! where wfc is the periodic part of the wavefunction:
    !              u(r) = sum_G e^{i G r} wfc(G);
    !
    ! The computation is done over all bands.
    !--------------------------------------------------------------------------------
    use intw_useful_constants, only: zero, one, cmplx_0, sig_x, sig_y, sig_z
    use intw_reading, only: nG_max, nspin, num_bands_intw

    implicit none

    !I/O variables

    complex(dp), intent(in) :: wfc(nG_max,num_bands_intw,nspin)
    complex(dp), intent(out) :: spin(num_bands_intw,3)

    !local variables

    integer :: ibnd, is, js, iG

    spin = cmplx_0
    !
    do iG=1,nG_max
      !
      do ibnd=1,num_bands_intw
        do is=1,nspin
          do js=1,nspin
            !
            spin(ibnd,1) = spin(ibnd,1) + 0.5d0*conjg(wfc(iG,ibnd,is))*sig_x(is,js)*wfc(iG,ibnd,js)
            !
            spin(ibnd,2) = spin(ibnd,2) + 0.5d0*conjg(wfc(iG,ibnd,is))*sig_y(is,js)*wfc(iG,ibnd,js)
            !
            spin(ibnd,3) = spin(ibnd,3) + 0.5d0*conjg(wfc(iG,ibnd,is))*sig_z(is,js)*wfc(iG,ibnd,js)
            !
          enddo !js
        enddo !is
      enddo !ibnd
      !
    enddo !iG

  end subroutine get_spin_component


  subroutine write_matrix_elements(filename, matrix_elements, num_bands, nb1, nb2)
    !------------------------------------------------------------------
    ! This subroutine writes out the matrix elements in a file
    !------------------------------------------------------------------
    use intw_utility, only: find_free_unit
    use intw_input_parameters, only: nk1, nk2, nk3
    use intw_useful_constants, only: zero, one

    implicit none

    complex(kind=dp), intent(in) :: matrix_elements(num_bands,num_bands,0:nk1+2,0:nk2+2,0:nk3+2)
    integer, intent(in) :: num_bands, nb1, nb2

    integer :: io_unit
    integer :: i, j, k
    real(kind=dp) :: x, y, z
    character(*) :: filename


    io_unit = find_free_unit()
    open(unit=io_unit, file=trim(filename), status='unknown')

    do i=1,nk1+1
      x = 1.0_dp*(i-1)/nk1

      do j=1,nk2+1
        y = 1.0_dp*(j-1)/nk2

        do k=1,nk3+1
          z = 1.0_dp*(k-1)/nk3

          write(io_unit,"(5f12.8)") x, y, z, matrix_elements(nb1,nb2,i,j,k)

        end do
      end do
    end do

    close(io_unit)

  end subroutine write_matrix_elements


  subroutine wfc_G_from_1D_to_3D(list_iG, wfc_G_1D, wfc_G_3D)
    !--------------------------------------------------------
    ! This subroutine puts a wavefunction, which is indexed
    ! by a scalar iG index, into a 3D array where the G
    ! vector is identified by a triple index.
    !--------------------------------------------------------

    use intw_reading, only: nr1, nr2, nr3, nG_max, nspin, num_bands_intw
    use intw_useful_constants, only: zero, one, cmplx_0
    use intw_utility, only: joint_to_triple_index_g
    use intw_fft, only: nl

    implicit none

    ! input
    integer, intent(in) :: list_iG(nG_max)
    complex(dp), intent(in) :: wfc_G_1D(nG_max,num_bands_intw,nspin)

    ! output
    complex(dp), intent(out) :: wfc_G_3D(nr1,nr2,nr3,num_bands_intw,nspin)

    ! computation variables
    integer :: i, iG

    integer :: i_joint
    integer :: n1, n2, n3


    ! initialize output array
    wfc_G_3D = cmplx_0

    ! loop on all G vectors in the array

    do i = 1,nG_max
      ! identify the G vector by its index, as stored in list_iG
      iG = list_iG(i)

      if (iG == 0) exit

      ! extract the scalar FFT index of this G vector
      i_joint = nl(iG)
      ! compute the triple index corresponding to iG
      call joint_to_triple_index_g(nr1,nr2,nr3,i_joint,n1,n2,n3)

      ! dump 1D wavefunction in 3D array

      ! careful! the wavefunction is indexed by i, not iG!
      wfc_G_3D(n1,n2,n3,:,:) = wfc_G_1D(i,:,:)
    enddo

  end subroutine wfc_G_from_1D_to_3D

end module intw_matrix_elements