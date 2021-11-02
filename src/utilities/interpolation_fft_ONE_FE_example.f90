program proba
  !
  use kinds, only: dp
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3, &
                                   ph_dir, dvscf_dir, dvscf_name, &
                                   calc_epmat, nq1, nq2, nq3, nqirr
  use intw_reading, only: bg, at, alat, tau, atom_labels, ityp, tpiba2, &
                          read_parameters_data_file_xml
  use intw_fft_interp
  use intw_utility, only: switch_indices_zyx
  use intw_useful_constants, only: pi, cmplx_0, ZERO
  !
  implicit none
  !
  integer :: nr1=10, nr2=10, nr3=10, &
             nr1s=15, nr2s=15, nr3s=15
  !
  complex(kind=dp), allocatable :: dv_sc_r(:)
  real(kind=dp), allocatable :: dv_sc_rs(:)
  real(kind=dp) :: x, y, z
  integer :: rl_dv_sc_r, rl_dv_sc_rs
  integer :: ia, id, imode, nat, ia_hand, imode_hand
  integer :: ios

  ! xsf format varibales
  logical :: xsf=.true.
  integer :: ja, jd, ir, i, j, k, kind
  real(kind=dp), dimension(3) :: x0, e1, e2
  real(kind=dp), allocatable, dimension(:,:,:) :: rho
  character(len=256) :: xsf_file

  mesh_dir = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/ONE_FE/SPIN/'
  prefix = 'Fe-MgO-Ag100'
  nk1 = 1
  nk2 = 1
  nk3 = 1
  calc_epmat=.true.
  ph_dir='PH/'
  dvscf_dir='./'
  dvscf_name='fe-mgo-ag100.dvscf_q1_PP'
  nq1=1
  nq2=1
  nq3=1
  nqirr=1

  call read_parameters_data_file_xml()


  bg(:,1)=(/1.0,0.0,0.0/)
  bg(:,2)=(/0.0,1.0,0.0/)
  bg(:,3)=(/0.0,0.0,1.0/)
  tpiba2 = 0.0_dp


  nr1=128
  nr2=128
  nr3=480
  nr1s=180
  nr2s=180
  nr3s=648
  nat=273


  allocate( dv_sc_r ( nr1* nr2* nr3) )
  allocate( dv_sc_rs(nr1s*nr2s*nr3s) )
  dv_sc_r = cmplx_0
  dv_sc_rs = ZERO
  !
  inquire(iolength=rl_dv_sc_r) dv_sc_r(:)
  inquire(iolength=rl_dv_sc_rs) dv_sc_rs(:)
  print*, rl_dv_sc_r, rl_dv_sc_rs
  !
  ! Open the files to read the dv
  open ( unit = 1234, file = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/M_ELEMENTS/mgo-ag100_R441.dvscf_q1', iostat = ios, &
         form = 'unformatted', status = 'old', access = 'direct', &
         action = 'read', recl = rl_dv_sc_r )
  open ( unit = 1235, file = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/M_ELEMENTS/fe-mgo-ag100.dvscf_hand_only', iostat = ios, &
         form = 'unformatted', status = 'old', access = 'direct', &
         action = 'read', recl = rl_dv_sc_rs )
  !
  !
  ! The dv calculated by hand are ordered in this way:
  ! 1.: atm = 273 ! Fe_1
  ! 3.: atm = 172 ! O_1
  ! 4.: atm = 186 ! O_2
  do ia=1,nat
    if ( .not. (ia==1) ) cycle ! Nahi dudan atomoa
    do id=1,3
      !
      if ( .not. (id==3) ) cycle ! Nahi dudan norabidea
      !
      imode = id + 3*(ia-1)
      !
      print*, imode, "/", 3 + 3*(nat-1)
      !
      if (ia==172) then
        ! O_1
        ia_hand=3
        imode_hand = id + 3*(ia_hand-1)
        !
        ! read the dv (this is already on the fine grid)
        read(unit=1235, rec=imode_hand, iostat=ios) dv_sc_rs
        !
      elseif (ia==186) then
        ! O_2
        ia_hand=4
        imode_hand = id + 3*(ia_hand-1)
        !
        ! read the dv (this is already on the fine grid)
        read(unit=1235, rec=imode_hand, iostat=ios) dv_sc_rs
        !
      elseif (ia==273) then
        ! Fe_1
        ia_hand=1
        imode_hand = id + 3*(ia_hand-1)
        !
        ! read the dv (this is already on the fine grid)
        read(unit=1235, rec=imode_hand, iostat=ios) dv_sc_rs
        !
      else
        !
        ! read dv on the coarse grid
        read(unit=1234, rec=imode, iostat=ios) dv_sc_r
        !
        ! interpolate to the fine grid using Fourier interpolation
        call fft_interp_3d_real(nr1, nr2, nr3, nr1s, nr2s, nr3s, huge(1.0_dp), real(dv_sc_r), dv_sc_rs)
        !
      endif
      !
      ! Write the variation of the potential in xsf format if needed
      if (xsf) then
        !
        ja=(imode-1)/3+1
        jd=modulo(imode-1,3)+1
        print*, imode, ja, jd
        !
        if ( imode> 0 .and. imode<  10 ) write(xsf_file,'(a18,i1)') 'fe.xsf_scf_nr_mode', imode
        if ( imode> 9 .and. imode< 100 ) write(xsf_file,'(a18,i2)') 'fe.xsf_scf_nr_mode', imode
        if ( imode>99 .and. imode<1000 ) write(xsf_file,'(a18,i3)') 'fe.xsf_scf_nr_mode', imode
        open(unit=123, file=trim(xsf_file), status="replace", action="write")
        !
        call xsf_struct (alat, at, nat, tau, atom_labels, ityp, 123)
        !
        allocate(rho(2,nr1,nr2))
        x0=0.0
        x0(3) = tau(3,ja)
        e1 = at(:,1)
        e2 = at(:,2)
        if ( tau(3,ja) < 0.0 ) then
          kind = nint( tau(3,ja)/at(3,3)*nr3 + 1 ) + nr3
        else
          kind = nint( tau(3,ja)/at(3,3)*nr3 + 1 )
        endif
        print*, 'kind=', kind
        !
        do ir=1,nr1*nr2*nr3
          !
          call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)
          if ( k==kind ) then
            rho(1,i,j) = real ( dv_sc_r(ir) )
            rho(2,i,j) = 0.0
          endif
          !
        enddo ! ir
        !
        call xsf_datagrid_2d (rho, nr1, nr2, 1.0_dp, 1.0_dp, x0, e1, e2, alat, 123)
        deallocate(rho)
        !
        close(unit=123)
        !
        !
        !
        if ( imode> 0 .and. imode<  10 ) write(xsf_file,'(a18,i1)') 'fe.xsf_scf_nrs_mode', imode
        if ( imode> 9 .and. imode< 100 ) write(xsf_file,'(a18,i2)') 'fe.xsf_scf_nrs_mode', imode
        if ( imode>99 .and. imode<1000 ) write(xsf_file,'(a18,i3)') 'fe.xsf_scf_nrs_mode', imode
        open(unit=1234, file=trim(xsf_file), status="replace", action="write")
        !
        call xsf_struct (alat, at, nat, tau, atom_labels, ityp, 1234)
        !
        allocate(rho(2,nr1s,nr2s))
        x0=0.0
        x0(3) = tau(3,ja)
        e1 = at(:,1)
        e2 = at(:,2)
        if ( tau(3,ja) < 0.0 ) then
          kind = nint( tau(3,ja)/at(3,3)*nr3s + 1 ) + nr3s
        else
          kind = nint( tau(3,ja)/at(3,3)*nr3s + 1 )
        endif
        print*, 'kind=', kind
        !
        do ir=1,nr1s*nr2s*nr3s
          !
          call switch_indices_zyx(nr1s,nr2s,nr3s,ir,i,j,k,-1)
          if ( k==kind ) then
            rho(1,i,j) = real ( dv_sc_rs(ir) )
            rho(2,i,j) = 0.0
          endif
          !
        enddo ! ir
        !
        call xsf_datagrid_2d (rho, nr1s, nr2s, 1.0_dp, 1.0_dp, x0, e1, e2, alat, 1234)
        deallocate(rho)
        !
        close(unit=1234)
        !
      endif ! xsf
      !
    enddo
    !
  enddo
  !
  deallocate(dv_sc_r, dv_sc_rs)
  !
  close(1234)
  close(1235)
  close(4321)
  !
  !



contains

  function f (xx,yy,zz)

    real(dp), intent(in) :: xx,yy,zz
    real(dp)             :: f

    f = -sin(2*pi*xx)*exp(-xx**2)*sin(2*pi*yy)*exp(-2*yy**2)*sin(4*pi*zz)*exp(-zz**2)

  end function f


end program proba
