program dv_hand
  !
  ! This program reads the data of the supercell frozen-phonon calculation to
  ! calculate the force-constant matrix (and then obtain the phonon frequencies)
  ! and also the potential change due to the phonon.
  !
  use kinds, only: dp
  use intw_input_parameters, only: mesh_dir, prefix, ph_dir, dvscf_dir, &
                                   dvscf_name, data_dir, data_name, &
                                   nk1, nk2, nk3, fc_mat, &
                                   qlist, nq1, nq2, nq3, nqirr
  use intw_reading, only: at, bg, alat, nat, ntyp, tau, amass, nr1, nr2, nr3, &
                          nkpoints_QE, ityp, atom_labels
  use intw_reading, only: write_tag, read_parameters_data_file_xml, &
                          read_kpoints_data_file_xml
  use intw_ph, only: read_ph_information_xml, readfc, mat_inv_four_t

  use intw_utility, only: find_free_unit, generate_kmesh, switch_indices_zyx

  use intw_fft, only: wfc_by_expigr

  use intw_symmetries, only: find_the_irreducible_k_set

  use intw_useful_constants
  use iotk_module
  USE fft_base,         ONLY : dfftp
  USE ions_base,        ONLY : zv, atm
  !
  implicit none
  !
  real(kind=dp), parameter :: pmass=1822.88848426_dp, aumev=  27211.396_dp
  !
  !
  complex(kind=dp), allocatable, dimension(:,:,:,:,:,:,:) :: frc_0 ! fc matrix at Gamma
  complex(kind=dp), allocatable, dimension(:,:) :: dyn_q ! dynamical matrix
  real(kind=dp), allocatable, dimension(:) :: w2 ! frequencies
  !
  real(kind=dp), allocatable, dimension(:,:) :: f_p, f_m
  real(kind=dp), allocatable, dimension(:) :: V_total_p, V_total_m
  real(kind=dp), allocatable, dimension(:) :: V_ionic_p, V_ionic_m
  real(kind=dp), allocatable, dimension(:,:) :: delta_V_scf
  real(kind=dp), allocatable, dimension(:,:) :: delta_V_loc
  !
  character(len=1024) :: scf_file_p, scf_file_m
  character(len=256) :: dummy1, dummy2
  character(len=1), dimension(3) :: disp
  real(kind=dp) :: dx = 0.001
  !
  integer :: ia, ja, id, jd, it, imode, i, j, k, ir
  !
  integer, parameter :: stdout=6
  ! character(len=1) :: CR=CHAR(13)
  !
  logical :: full_zone
  integer :: nkmesh, nqmesh, nk_irr, nq_irr
  real(kind=dp), allocatable, dimension(:,:) :: kmesh, kpoints_irr, kpoints_QE
  !
  integer :: ia_original, ja_original
  real(kind=dp), dimension(3) :: R_l, R_m, R
  integer :: m1, m2, m3
  !
  integer :: plot_num, ibrav
  real(kind=dp) :: ecutwfc, dual, gcutm
  real(kind=dp), dimension(6) :: celldm
  character (len=75) :: title
  !
  real(kind=dp), dimension(3) :: x0, e1, e2
  real(kind=dp), allocatable, dimension(:,:,:) :: rho
  character(len=256) :: xsf_file
  integer :: ounit
  integer :: kind
  !
  !
  real(kind=dp), parameter :: b2a=0.529177249_dp
  complex(kind=dp), allocatable, dimension(:,:) :: z



  20 format(A)
  !
#ifdef DEBUG
  2100 FORMAT('DEBUG: ',a15,x,i4)
  2200 FORMAT('DEBUG: ',a15,x,3i4)
  2300 FORMAT('DEBUG: ',a15,x,f16.8)
  2400 FORMAT('DEBUG: ',a15,x,3f16.8)
#endif
  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(stdout,20) '======================================================'
  write(stdout,20) '|                 program dv_hand.x                  |'
  write(stdout,20) '|        ------------------------------------        |'
  write(stdout,20) '======================================================'
  !
  !================================================================================
  !
  ! Directory and files of the SC calculation
  mesh_dir = '/home/haritz/Kalkuluak/Probak/grafenoa/SC/0.001/SCF/v0/'
  prefix = 'g'
  nk1 = 1
  nk2 = 1
  nk3 = 1
  fc_mat='g.fc'
  data_dir='data_dir/'
  data_name='data-file.'
  dvscf_dir='dvscf_dir/'
  dvscf_name='g.dvscf_q'
  !
  !
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)
  !
  !
  !
  write(stdout,20) 'Reading the QE calculation data of the SC'
  !
  call read_parameters_data_file_xml()
  !
  allocate(kpoints_QE(3,nkpoints_QE))
  call read_kpoints_data_file_xml(kpoints_QE)
  !
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Calculate the FC matrix'
  write(stdout,20) ''
  !
  ! Calculate the FC matrix
  allocate( frc_0(1,1,1,3,3,nat,nat) )
  allocate(f_m(3,nat),f_p(3,nat))
  frc_0 = cmplx_0
  f_m = zero
  f_p = zero
  !
  disp(1) = 'x'
  disp(2) = 'y'
  disp(3) = 'z'
  dx = 0.001_dp * alat
  !
  mesh_dir = '/home/haritz/Kalkuluak/Probak/grafenoa/SC/0.001/SCF/'
  !
  !
  !
  do ia=1,nat ! for each atom in the unit cell
    do id=1,3 ! for each direction
      !
      ! Read forces for + and - displacaments
      write(dummy1,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'+'
      scf_file_p = trim(mesh_dir)//'/'//trim(dummy1)//'/'//trim(dummy1)//'.scf.out'
      write(dummy2,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'-'
      scf_file_m = trim(mesh_dir)//'/'//trim(dummy2)//'/'//trim(dummy2)//'.scf.out'
      call read_f(scf_file_p,nat,f_p)
      call read_f(scf_file_m,nat,f_m)
      !
      do ja=1,nat
        do jd=1,3
          !
          ! Calculate the FC matrix with the derivative of force
          frc_0(1,1,1,id,jd,ia,ja) = -( f_p(jd,ja) - f_m(jd,ja) ) / ( 2.0_dp * dx )
          !
        enddo
      enddo
      !
    enddo
  enddo
  !
  ! Dealocate unnecessary arrays
  deallocate(f_p,f_m)
  !
  !
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Calculate the dynamical matrix'
  write(stdout,20) ''
  !
  allocate( dyn_q(3*nat,3*nat), w2(3*nat) )
  allocate(z(3,nat))
  !
  dyn_q=cmplx_0
  w2=zero
  !
  open(unit=112, file='vib.xyz')
  !
  !
  call mat_inv_four_t( (/ 0.0_dp, 0.0_dp, 0.0_dp /) ,1,1,1,3*nat,frc_0,dyn_q(:,:))
  !
  ! diagonalize the dynamical matrix to obtain eigenmodes
  call diagonalize_cmat( 3*nat, dyn_q, w2 )
  !
  do imode=1,3*nat
    !
    ! Normalize the eigenvector
    dyn_q(:,imode) = dyn_q(:,imode) / sqrt(sum(dyn_q(:,imode)**2.0_dp))
    !
    ! Avoid imaginary frequencies
    if (w2(imode)<zero) then
      w2(imode)=-sqrt(abs(w2(imode)))
    else
      w2(imode)=sqrt(w2(imode))
    endif
    !
    write(*,'(a,i3,2(f12.6,x,a3))') 'mode', imode, w2(imode), 'meV', w2(imode)*0.24180, 'THz'
    !
    ! save polarization
    write(112,*) nat
    write(112,'(a6,i4,a5,3f8.4,a5,1f10.6,a5)') 'Mode #',imode, ' q = ', (/ 0.0_dp, 0.0_dp, 0.0_dp /), ' f = ', w2(imode),' meV.'
    do ia = 1,nat
      do id=1,3
        z(id,ia) = dyn_q(3*(ia-1)+id,imode)
      enddo
      write(112,'(a4,x,3(2f10.5))') atom_labels(ityp(ia)), ( tau(id,ia)/nq1*alat*b2a ,id=1,3 ), (real(z(id,ia)), id=1,3)
    enddo
    !
  enddo ! imode
  !
  !
  close(112)
  deallocate(dyn_q,w2,z)
  !
  !
  ! Calculate the potential change dV
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Calculate the potential change'
  write(stdout,20) ''
  !

  allocate( V_total_p(nr1*nr2*nr3), V_total_m(nr1*nr2*nr3) )
  allocate( V_ionic_p(nr1*nr2*nr3), V_ionic_m(nr1*nr2*nr3) )
  allocate( delta_V_scf(nr1*nr2*nr3,3*nat) )
  allocate( delta_V_loc(nr1*nr2*nr3,3*nat) )
  allocate( rho(2,nr1,nr2) )
  rho = 0.0_dp
  !
  do ia=1,nat ! for each atom in the unit cell
    do id=1,3 ! for each direction
      !
      imode = id + 3*(ia-1)
      !
      ! Read total potential for + and - displacaments
      write(dummy1,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'+'
      scf_file_p = trim(mesh_dir)//'/'//trim(dummy1)//'/'//trim(dummy1)//'_total_pot'
      write(dummy2,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'-'
      scf_file_m = trim(mesh_dir)//'/'//trim(dummy2)//'/'//trim(dummy2)//'_total_pot'
      !
      ! read everything (requires that all variables that are read are allocated
      ! with the correct dimensions!)
      CALL plot_io (trim(scf_file_p), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_total_p, - 1)
      CALL plot_io (trim(scf_file_m), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_total_m, - 1)
      !
      ! Read ionic potential for + and - displacaments
      write(dummy1,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'+'
      scf_file_p = trim(mesh_dir)//'/'//trim(dummy1)//'/'//trim(dummy1)//'_ionic_pot'
      write(dummy2,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'-'
      scf_file_m = trim(mesh_dir)//'/'//trim(dummy2)//'/'//trim(dummy2)//'_ionic_pot'
      !
      ! read everything (requires that all variables that are read are allocated
      ! with the correct dimensions!)
      CALL plot_io (trim(scf_file_p), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_ionic_p, - 1)
      CALL plot_io (trim(scf_file_m), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_ionic_m, - 1)

      delta_V_scf(:,imode) = ( (V_total_p - V_ionic_p) - (V_total_m - V_ionic_m) ) / ( 2.0_dp * dx )
      delta_V_loc(:,imode) = ( V_ionic_p - V_ionic_m ) / ( 2.0_dp * dx )
      !
      ! save scf part
      if (imode<10) write(xsf_file,'(a23,i1)') 'graphene_xsf_dvscf_mode', imode
      if (imode>9)  write(xsf_file,'(a23,i2)') 'graphene_xsf_dvscf_mode', imode
      ounit = find_free_unit()
      open(unit=ounit, file=xsf_file, status="replace", action="write")
      call xsf_struct (alat, at, nat, tau, atom_labels, ityp, ounit)
      x0=0.0
      x0(3) = tau(3,ia)
      e1 = at(:,1)
      e2 = at(:,2)
      if ( tau(3,ia) < 0.0 ) then
        kind = nint( tau(3,ia)/at(3,3)*nr3 + 1 ) + nr3
      else
        kind = nint( tau(3,ia)/at(3,3)*nr3 + 1 )
      endif
      print*, 'kind', kind

      do ir=1,nr1*nr2*nr3
        call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)
        if (k==kind) then
          rho(1,i,j) = real ( delta_V_scf(ir, imode) )
        endif
      enddo
      call xsf_datagrid_2d (rho, nr1, nr2, 1.0_dp, 1.0_dp, x0, e1, e2, alat, ounit)
      close(unit=ounit)
      !
      ! save local part of the pp
      if (imode<10) write(xsf_file,'(a23,i1)') 'graphene_xsf_dvloc_mode', imode
      if (imode>9)  write(xsf_file,'(a23,i2)') 'graphene_xsf_dvloc_mode', imode
      ounit = find_free_unit()
      open(unit=ounit, file=xsf_file, status="replace", action="write")
      call xsf_struct (alat, at, nat, tau, atom_labels, ityp, ounit)
      x0=0.0
      x0(3) = tau(3,ia)
      e1 = at(:,1)
      e2 = at(:,2)
      if ( tau(3,ia) < 0.0 ) then
        kind = nint( tau(3,ia)/at(3,3)*nr3 + 1 ) + nr3
      else
        kind = nint( tau(3,ia)/at(3,3)*nr3 + 1 )
      endif
      print*, 'kind', kind

      do ir=1,nr1*nr2*nr3
        call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)
        if (k==kind) then
          rho(1,i,j) = real ( delta_V_loc(ir, imode) )
        endif
      enddo
      call xsf_datagrid_2d (rho, nr1, nr2, 1.0_dp, 1.0_dp, x0, e1, e2, alat, ounit)
      close(unit=ounit)


      !
    enddo
  enddo





contains

  subroutine read_f(scf_file,nat,f)
    !
    ! subroutine to read the forces from a scf file
    !
    use, intrinsic :: iso_fortran_env, only: iostat_end   ! end of file
    !
    implicit none
    !
    ! input variables
    character(len=1024), intent(in) :: scf_file
    integer, intent(in) :: nat
    real(kind=dp), dimension(3,nat), intent(out) :: f
    !
    ! local variables
    character(256) :: line
    integer :: stat
    integer :: iounit, ios, i, dummy_i1, dummy_i2
    !
    !
    iounit = find_free_unit()
    ! open(unit=iounit, file=trim(scf_file), iostat=ios, status="old", action="read")
    open(unit=iounit, file=trim(scf_file), status="old", action="read")
    ios=0
    if ( ios /= 0 ) then
      write(stdout,'(a)') "read_f: Error opening scf_file "//trim(scf_file)
      stop
    else
      write(stdout,'(a)') "read_f: Reading forces from scf_file "//trim(scf_file)
    endif
    !
    f=0.0_dp
    !
    do
      !
      read(unit=iounit, fmt='(a)', iostat=stat) line
      !
      if (trim(adjustl(line))=='Forces acting on atoms (Ry/au):') then
        !
        read(unit=iounit,fmt=*) ! read white line
        !
        do i=1,nat
          read(unit=iounit,fmt=9035) dummy_i1, dummy_i2, f(:,i)
          if (dummy_i1/=i) write(stdout,'(a)') 'read_f: WARNING, atoms are not properly ordered on scf file'
        enddo
        !
        close(unit=iounit)
        return
        !
      endif
      !
      if ( stat == iostat_end ) then
        write(stdout,*) 'Error reading file '//trim(scf_file)//'. End of file before reading forces.'
        close(unit=iounit)
        return
      end if
      !
    end do
    !
    close(unit=iounit)
    !
    9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
    !
  end subroutine read_f

  subroutine diagonalize_cmat(n,a,w)
    !
    ! subroutine to diagonalize a complex matrix
    !
    implicit none
    !
    ! input variables
    integer, intent(in)  :: n
    complex(dp),intent(inout) :: a(n,n) ! matrix to diagonalize on input, eigenvector on output
    real(dp),intent(out) :: w(n) ! eigenvalues

    complex(dp) :: a_pack(n*(n+1)/2)

    integer :: i,j, nfound

    complex(dp) :: cwork(2*n)
    real   (dp) :: rwork(7*n)
    integer     :: iwork(5*n), ifail(n), info


    do j=1,n
       do i=1,j
          a_pack(i+((j-1)*j)/2)=a(i,j)
       enddo
    enddo

    ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
    call ZHPEVX('V', 'A', 'U', n, a_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp,  &
         nfound, w, a, n,  cwork,  rwork,       iwork,          ifail, info)

  end subroutine diagonalize_cmat


end program dv_hand
