program proba_pp

  ! Module variables
  use kinds, only: dp
  !
  use intw_useful_constants, only: cmplx_0, zero
  !
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3, &
                                   ph_dir, dvscf_dir, dvscf_name, &
                                   nq1, nq2, nq3, nqirr
  !
  use intw_pseudo
  !
  use intw_haritz, only: gamma_only, nlm
  !
  use intw_utility, only: find_free_unit, switch_indices_zyx, cryst_to_cart
  !
  use intw_eqv, only: vlocq
  use intw_ph, only: eigqts
  !
  use intw_haritz, only: get_gamma_only_haritz
  !
  USE intw_reading
  USE intw_uspp_param
  USE intw_uspp
  !
  !
  implicit none
  !
  ! dv related variables
  real(kind=dp), allocatable :: dv_local(:) ! local part of the variation of the PP in real space (nr1*nr2*nr3)

  ! Loop variables
  integer :: imode
  integer :: ir

  ! I/O variables
  integer :: dv_read_unit
  integer :: dv_write_unit
  integer :: rl
  integer :: ios
  character(len=256) :: dv_read_file
  character(len=256) :: dv_write_file

  ! xsf format varibales
  logical :: xsf=.false.
  integer :: ia, id, i, j, k, kind, jind
  real(kind=dp), dimension(3) :: x0, e1, e2
  real(kind=dp), allocatable, dimension(:,:,:) :: rho
  character(len=256) :: xsf_file
  !
  !
  !
  !================================================================================
  !       Specify input parameters
  !================================================================================
  !
  mesh_dir = '/home/haritz/Kalkuluak/Probak/Fe-MgO-Ag100/bat/'
  prefix = 'Fe-MgO-Ag100'
  nk1 = 1
  nk2 = 1
  nk3 = 1
  ph_dir='PH/'
  dvscf_dir='./'
  dvscf_name='fe-mgo-ag100.dvscf_q1'
  nq1=1
  nq2=1
  nq3=1
  nqirr=1
  !
  !
  !================================================================================
  !       read the parameters from the SCF QE calculation
  !================================================================================
  !
  call read_parameters_data_file_xml()
  !
  call get_gamma_only_haritz()
  !
  ! allocate(gvec(3,ngm)) ! it is deallocated in the bottom of this main code.
  ! call get_gvec()
  ! !
  ! call allocate_fft()
  ! !
  ! call generate_nl()
  !
  print*, "Reading PPs"
  call read_all_pseudo ()
  if (.not.lspinorb) call average_pp(ntyp)
  !
  call allocate_nlpot1
  !
  print*, ""
  !
  do ia=1,ntyp
    print*, upf(ia)%psd
    print*, upf(ia)%lloc
    print*, upf(ia)%nbeta
    print*, upf(ia)%lll
    ! do i=1,upf(ia)%nbeta
    !     print*, ( upf(ia)%dion(i,j) , j=1,upf(ia)%nbeta )
    ! enddo
    ! print*, upf(ia)%
  enddo
  !
  print*, ""
  !
  print*, lmaxkb
  print*, nhm

contains

  subroutine generate_nl()
    !------------------------------------------------------------------------
    !  This subroutine generates information important for the 3D-FFT
    !  algorithm.
    !------------------------------------------------------------------------
    ! use intw_fft, only: gvec_cart, eigts1, eigts2, eigts3, nl, strf
    use intw_fft, only: gvec_cart, nl, strf
    use intw_useful_constants, only: tpi
    !
    implicit none
    !
    !local
    integer  :: n1, n2, n3
    integer  :: ng, na
    !
    integer  :: switch
    !
    logical  :: assigned(nr1,nr2,nr3)
    !
    integer :: nt
    !
    real(dp) :: arg
    !output :: nl and ig1 ig2 ig3 and phase
    nl(:)            = 0
    assigned         = .false.
    !
    !
    switch    =   1   ! triplet - to - singlet index
    !
    ! loop on all G vectors in the global array gvec
    do ng=1,ngm
      !
      n1 = modulo( gvec(1,ng), nr1 ) + 1
      n2 = modulo( gvec(2,ng), nr2 ) + 1
      n3 = modulo( gvec(3,ng), nr3 ) + 1
      !
      if ( .not. assigned(n1,n2,n3) ) then
        !
        assigned(n1,n2,n3) = .true.
        !
        ! compute the scalar index corresponding to n1,n2,n3 and
        ! assign it to nl(ng)
        call switch_indices_zyx(nr1,nr2,nr3,nl(ng),n1,n2,n3,switch)
        !
      else
        !
        write(*,*) 'ERROR in generate_nl. FFT mesh too small?'
        write(*,*) '    More than one G-vector in the gvec array are being'
        write(*,*) '    assigned to the same FFT triplet (n1,n2,n3);      '
        write(*,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
        write(*,*) '    too small.                                        '
        !
        stop
        !
      endif
      !
      if ( gamma_only ) then
        !
        ! Negative G
        !
        ! Gamma
        if ( gvec(1,ng)==-gvec(1,ng) .and. gvec(2,ng)==-gvec(2,ng) .and. gvec(3,ng)==-gvec(3,ng) ) assigned(n1,n2,n3)=.false.
        !
        n1 = modulo( -gvec(1,ng), nr1 ) + 1
        n2 = modulo( -gvec(2,ng), nr2 ) + 1
        n3 = modulo( -gvec(3,ng), nr3 ) + 1
        !
        if ( .not. assigned(n1,n2,n3) ) then
          !
          assigned(n1,n2,n3) = .true.
          !
          ! compute the index corresponding to n1,n2,n3 and
          ! assign it to nl(ng)
          call switch_indices_zyx(nr1,nr2,nr3,nlm(ng),n1,n2,n3,switch)
          !
        else
          !
          write(*,*) 'ERROR in generate_nl. FFT mesh too small?'
          write(*,*) '    More than one G-vector in the gvec array are being'
          write(*,*) '    assigned to the same FFT triplet (n1,n2,n3);      '
          write(*,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
          write(*,*) '    too small.                                        '
          stop
          !
        endif
        !
      endif
      !
    end do
    !
    ! Obtain the G vectors in cartesian coordinates.
    gvec_cart(1:3,1:ngm) = gvec(1:3,1:ngm)
    call cryst_to_cart (ngm, gvec_cart , bg, 1)
    !
    strf(:,:) = (0.d0,0.d0)
    do nt = 1, ntyp
       do na = 1, nat
          if (ityp (na) .eq.nt) then
             do ng = 1, ngm
                arg = (gvec_cart (1, ng) * tau (1, na) + gvec_cart (2, ng) * tau (2, na) &
                     + gvec_cart (3, ng) * tau (3, na) ) * tpi
                strf (ng, nt) = strf (ng, nt) + CMPLX(cos (arg), -sin (arg),kind=DP)
             enddo
          endif
       enddo
    enddo
    !
  end subroutine generate_nl



  subroutine allocate_fft()
    !--------------------------------------------------------
    ! This subroutine simply allocates the arrays needed
    ! by the fft algorithms.
    !--------------------------------------------------------
    !
    use intw_fft, only: gvec_cart, nl, strf
    !
    allocate (nl(ngm))
    if (gamma_only) allocate(nlm(ngm))
    allocate (gvec_cart(3,ngm))
    allocate (strf(ngm, ntyp))
    !
  end subroutine allocate_fft

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





end program proba_pp
