program add_pp_local
  !
  ! Program to add the derivative of the local part of the PP
  ! to the induced variation of the potential.
  !

  ! Module variables
  use kinds, only: dp
  !
  use intw_useful_constants, only: cmplx_0, zero
  !
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3, &
                                   ph_dir, dvscf_dir, dvscf_name, &
                                   nq1, nq2, nq3, nqirr
  !
  use intw_reading, only: lspinorb, ntyp, ngm, gvec, nr1, nr2, nr3, nat, ityp, &
                          atom_labels, tau, at, alat, bg
  !
  use intw_pseudo, only: read_all_pseudo, average_pp
  !
  use intw_haritz, only: gamma_only, nlm
  !
  ! Module subroutines and functions
  use intw_reading, only: read_parameters_data_file_xml, get_gvec
  !
  use intw_utility, only: find_free_unit, switch_indices_zyx, cryst_to_cart
  !
  use intw_eqv, only: vlocq
  use intw_ph, only: eigqts
  !
  use intw_haritz, only: get_gamma_only_haritz
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
  mesh_dir = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/ONE_FE/SPIN/'
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
  allocate(gvec(3,ngm)) ! it is deallocated in the bottom of this main code.
  call get_gvec()
  !
  call allocate_fft()
  !
  call generate_nl()


  ! allocate (strf(ngm, ntyp))
  ! allocate (mill(3,ngm))
  !
  call read_all_pseudo ()
  if (.not.lspinorb) call average_pp(ntyp)
  !
  call allocate_nlpot1
  !
  ! call allocate_phq
  allocate (eigqts(nat))
  allocate (vlocq ( ngm , ntyp))

  !
  call init_us_1
  !
  call allocate_nlpot2 ( (/ 0.0_dp, 0.0_dp, 0.0_dp /) )
  !
  call phq_init( (/ 0.0_dp, 0.0_dp, 0.0_dp /) )
  !
  print*, 'Enter on calculate_local_part_dv'
  allocate(dv_local(nr1*nr2*nr3))
  if (xsf) allocate(rho(2,nr1,nr2))
  !
  !
  !
  ! Record length of the variable to read
  inquire(iolength=rl) dv_local
  !
  ! Open the files
  dv_read_file= trim(mesh_dir) // trim(ph_dir) // trim(dvscf_dir) // trim(dvscf_name)
  dv_write_file= trim(mesh_dir) // trim(ph_dir) // trim(dvscf_dir) // trim(dvscf_name) // trim("_PP")
  !
  dv_read_unit=find_free_unit()
  dv_write_unit=find_free_unit()
  !
  open( unit=dv_read_unit , file=trim(dv_read_file) , iostat=ios, form='unformatted', &
        status='old', action='read' , access='direct', recl=rl )
  open( unit=dv_write_unit, file=trim(dv_write_file), iostat=ios, form='unformatted', &
        status='new', action='write', access='direct', recl=rl )

  do imode=1,3*nat
    !
    ! Read the variation of the potential
    read(unit=dv_read_unit,rec=imode,iostat=ios) dv_local
    !
    ! Add the local part of the PP
    call calculate_local_part_dv( (/ 0.0_dp, 0.0_dp, 0.0_dp /), imode, dv_local)
    !
    ! Save the variation of the potential
    write(unit=dv_write_unit,rec=imode,iostat=ios) dv_local
    !
    ! Write the variation of the potential in xsf format if needed
    if (xsf) then
      !
      ia=(imode-1)/3+1
      id=modulo(imode-1,3)+1
      print*, imode, ia, id
      !
      if (imode<10) write(xsf_file,'(a17,i1)') 'fe.xsf_ionic_mode', imode
      if (imode>9) write(xsf_file,'(a17,i2)') 'fe.xsf_ionic_mode', imode
      open(unit=123, file=xsf_file, status="replace", action="write")
      !
      call xsf_struct (alat, at, nat, tau, atom_labels, ityp, 123)
      !
      x0=0.0
      x0(3) = tau(3,ia)
      e1 = at(:,1)
      e2 = at(:,2)
      if ( tau(3,ia) < 0.0 ) then
        kind = nint( tau(3,ia)/at(3,3)*nr3 + 1 ) + nr3
      else
        kind = nint( tau(3,ia)/at(3,3)*nr3 + 1 )
      endif
      if ( tau(2,ia) < 0.0 ) then
        jind = nint( tau(2,ia)/at(2,2)*nr2 + 1 ) + nr2
      else
        jind = nint( tau(2,ia)/at(2,2)*nr2 + 1 )
      endif
      kind=nr3/2
      print*, jind, kind
      !
      do ir=1,nr1*nr2*nr3
        !
        call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)
        if ( k==kind ) then
          rho(1,i,j) = real ( dv_local(ir) )
          rho(2,i,j) = 0.0
          if ( j==jind ) then
            write(999,*) real(dv_local(ir))
          endif
        endif
        !
      enddo ! ir
      !
      call xsf_datagrid_2d (rho, nr1, nr2, 1.0_dp, 1.0_dp, x0, e1, e2, alat, 123)
      !
      close(unit=123)
      !
    endif ! xsf
    !
  enddo
  !
  ! Close dv_unit
  close( unit =  dv_read_unit )
  close( unit = dv_write_unit )



contains



  subroutine calculate_local_part_dv(qpoint, mode, dvq_local )
    !
    !======================================================================
    ! We have dV_scf as input and we add to it the derivative of the PP   !
    !======================================================================
    !
    USE intw_reading , ONLY :  nr1, nr2, nr3, ngm, tau
    USE intw_fft, ONLY : nl
    use intw_reading, only: tpiba
    use intw_useful_constants, only: cmplx_i, tpi
    use intw_fft, only: gvec_cart

    !
    implicit none
    !
    !I/O variables
    !
    real(dp),intent(in) :: qpoint(1:3) ! crystal coord.
    integer,intent(in) :: mode
    real(dp),intent(inout) :: dvq_local(nr1*nr2*nr3) ! spin idependentea da baina koherentzia mantenduko dugu.
    !
    !local variables
    !
    integer :: na, ipol, ig, nt, ir
    complex(dp) :: aux (nr1*nr2*nr3), fact, gtau
    real(dp) :: qcart(3) ! qpoint in cart.
    real(dp) :: arg
    !
    qcart=matmul(bg,qpoint)
    !
    !
    na=(mode-1)/3+1
    ipol=modulo(mode-1,3)+1
    nt=ityp (na)
    !
    aux(:)=cmplx_0
    !
    fact=-tpiba*cmplx_i*eigqts(na)
    !
    do ig=1,ngm
      !
      arg = tpi * ( gvec_cart(1,ig)*tau(1,na) + gvec_cart(2,ig)*tau(2,na) + gvec_cart(3,ig)*tau(3,na) ) ! G*tau
      gtau=exp(-cmplx_i*arg)
      !
      aux(nl(ig))=aux(nl(ig))+vlocq(ig,nt)*(qcart(ipol)+gvec_cart(ipol,ig))*gtau*fact
      !
    enddo !ig
    !
    if (gamma_only) then
      do ig = 1, ngm
        aux(nlm(ig)) = CONJG(aux(nl(ig)))
      enddo
    end if
    !
    call cfftnd(3,(/nr1,nr2,nr3/),1,aux)
    !
    do ir=1,nr1*nr2*nr3
      !
      dvq_local(ir)=dvq_local(ir)+real(aux(ir))
      !
    enddo !ir
    !
    return
    !
  end subroutine calculate_local_part_dv

  subroutine calculate_local_part_dv_2(qpoint, nat, mode, dvq_local )
    !
    !======================================================================
    ! We have dV_scf as input and we add to it the derivative of the PP   !
    !======================================================================
    !
    USE intw_reading , ONLY :  nr1, nr2, nr3, ngm, tau
    USE intw_fft, ONLY : nl
    use intw_reading, only: tpiba
    use intw_useful_constants, only: cmplx_i, tpi
    use intw_fft, only: gvec_cart, strf
    !
    implicit none
    !
    !I/O variables
    !
    real(dp),intent(in) :: qpoint(1:3) ! crystal coord.
    integer,intent(in) :: nat, mode
    complex(dp),intent(inout) :: dvq_local(nr1*nr2*nr3) ! spin idependentea da baina koherentzia mantenduko dugu.
    !
    !local variables
    !
    integer :: na, ipol, ig, nt, ir, ia, it
    complex(dp) :: aux1 (nr1*nr2*nr3), aux2 (nr1*nr2*nr3), fact
    real(dp) :: qcart(3) ! qpoint in cart.
    real(dp) :: arg, dx
    !
    !
    dx = 0.0001_dp * alat
    !
    qcart=matmul(bg,qpoint)
    !
    !
    na=(mode-1)/3+1
    ipol=modulo(mode-1,3)+1
    nt=ityp (na)
    !
    aux1(:)=cmplx_0
    aux2(:)=cmplx_0
    aux1(:)=cmplx_0
    aux2(:)=cmplx_0
    !
    !
    ! Calculate the structure factor with the displaced atom
    strf(:,:) = (0.d0,0.d0)
    tau(ipol,na) = tau(ipol,na) + dx/alat
    do it = 1, ntyp
      do ia = 1, nat
        if (ityp (ia) .eq. it) then
          do ig = 1, ngm
            !
            arg = (gvec_cart (1, ig) * tau (1, ia) + gvec_cart (2, ig) * tau (2, ia) &
                  + gvec_cart (3, ig) * tau (3, ia) ) * tpi
            !
            strf (ig, nt) = strf (ig, nt) + CMPLX(cos (arg), -sin (arg),kind=DP)
            !
          enddo
        endif
      enddo
    enddo
    tau(ipol,na) = tau(ipol,na) - dx/alat
    !
    !
    ! Calculate the ionic potential with the displaced atom
    do ig=1,ngm
      !
      ! aux1(nl(ig))=aux1(nl(ig))+vlocq(ig,nt)*strf(ig,nt)
      aux1(nl(ig))=vlocq(ig,nt)*strf(ig,nt)
      !
    enddo !ig
    !
    if (gamma_only) then
      do ig = 1, ngm
        aux1(nlm(ig)) = CONJG(aux1(nl(ig)))
      enddo
    end if
    !
    !
    ! Calculate the structure factor of the unit cell
    strf(:,:) = (0.d0,0.d0)
    do it = 1, ntyp
      do ia = 1, nat
        if (ityp (ia) .eq. it) then
          do ig = 1, ngm
            !
            arg = (gvec_cart (1, ig) * tau (1, ia) + gvec_cart (2, ig) * tau (2, ia) &
                  + gvec_cart (3, ig) * tau (3, ia) ) * tpi
            !
            strf (ig, nt) = strf (ig, nt) + CMPLX(cos (arg), -sin (arg),kind=DP)
            !
          enddo
        endif
      enddo
    enddo
    !
    !
    ! Calculate the ionic potential of the unit cell
    do ig=1,ngm
      !
      aux2(nl(ig))=aux2(nl(ig))+vlocq(ig,nt)*strf(ig,nt)
      !
    enddo !ig
    !
    if (gamma_only) then
      do ig = 1, ngm
        aux2(nlm(ig)) = CONJG(aux2(nl(ig)))
      enddo
    end if
    !
    !
    ! FFT to obtain the potential in real space
    call cfftnd(3,(/nr1,nr2,nr3/),1,aux1)
    call cfftnd(3,(/nr1,nr2,nr3/),1,aux2)
    !
    !
    ! Use finit differences to obtain the derivative of the ionic potential
    do ir=1,nr1*nr2*nr3
      !
      dvq_local(ir) = ( (aux1(ir)) - (aux2(ir)) ) / ( 1.0_dp * dx )
      !
    enddo !ir
    !
    return
    !
  end subroutine calculate_local_part_dv_2


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

  subroutine read_dv(mode, dv)
    !
    ! This subroutine reads the variation of the potential induced by
    ! the cartesian displacement # mode
    !
    implicit none
    !
    ! I/O variables
    integer, intent(in) :: mode
    real(kind=dp), intent(out) :: dv(:)
    !
    ! local variables
    integer :: read_unit, rl, ios
    character(len=256) :: dv_file
    !
    !
    ! Initialize the output variable
    dv = ZERO
    !
    ! Complete path to the file
    dv_file= trim(mesh_dir) // trim(ph_dir) // trim(dvscf_dir) // trim(dvscf_name)
    !
    ! Record length of the variable to read
    inquire(iolength=rl) dv
    !
    ! Open the file
    read_unit=find_free_unit()
    open( unit=read_unit, file=trim(dv_file), iostat=ios, form='unformatted', &
          status='old', action='read', access='direct', recl=rl )
    !
    ! Read the variation of the potential
    read(unit=read_unit,rec=mode,iostat=ios) dv
    !
    ! Close the file
    close(unit=read_unit)
    !
  end subroutine read_dv



end program add_pp_local
