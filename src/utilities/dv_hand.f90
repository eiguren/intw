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
  ! Parameters of the simple Unit Cell
  integer :: nat_uc!, ntyp_uc, nr1_uc, nr2_uc, nr3_uc, nkpoints_QE_uc
  ! integer, allocatable, dimension(:) :: ityp_uc
  ! real(kind=dp) :: alat_uc
  ! real(kind=dp), dimension(3,3) :: at_uc, bg_uc
  ! real(kind=dp), allocatable, dimension(:) :: amass_uc
  ! real(kind=dp), allocatable, dimension(:,:) :: tau_uc
  ! character(len=3), allocatable, dimension(:)  :: atom_labels_uc
  !
  ! Parameters of the Super Cell
  integer :: nat_sc, ntyp_sc, nr1_sc, nr2_sc, nr3_sc, nkpoints_QE_sc
  integer, allocatable, dimension(:) :: ityp_sc
  real(kind=dp) :: alat_sc
  real(kind=dp), dimension(3,3) :: at_sc, bg_sc
  real(kind=dp), allocatable, dimension(:) :: amass_sc
  real(kind=dp), allocatable, dimension(:,:) :: tau_sc
  character(len=3), allocatable, dimension(:)  :: atom_labels_sc
  !
  !
  complex(kind=dp), allocatable, dimension(:,:,:,:,:,:,:) :: frc_0 ! fc matrix at Gamma
  complex(kind=dp), allocatable, dimension(:,:) :: dyn_q ! dynamical matrix
  real(kind=dp), allocatable, dimension(:) :: w2 ! frequencies
  !
  real(kind=dp), allocatable, dimension(:,:) :: f_p, f_m
  real(kind=dp), allocatable, dimension(:) :: V_total_p, V_total_m, &
                                              V_ionic_p, V_ionic_m
  real(kind=dp), allocatable, dimension(:,:) :: delta_V
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
  real(kind=dp), allocatable, dimension(:,:) :: kmesh, qmesh, &
                                                kpoints_irr, qpoints_irr, &
                                                kpoints_QE_sc
  !
  !
  !
  !
  !
  integer :: ia_original, ja_original
  real(kind=dp), dimension(3) :: R_l, R_m, R
  integer :: m1, m2, m3

  integer :: plot_num, ibrav
  real(kind=dp) :: ecutwfc, dual, gcutm
  real(kind=dp), dimension(6) :: celldm
  character (len=75) :: title

  real(kind=dp), dimension(3) :: x0, e1, e2
  real(kind=dp), allocatable, dimension(:,:,:) :: rho
  character(len=256) :: xsf_file
  integer :: ounit
  !
  !
  logical :: sc_from_uc=.false.

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
  mesh_dir = 'SCF/v0/'
  prefix = 'g'
  nk1 = 12
  nk2 = 12
  nk3 =  1
  !
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)
  !
  !
  !
  write(stdout,20) 'Reading the QE calculation data of the SC'
  !
  ! Read from data-file:
  ! -alat (Bohr)
  ! -at (in units of alat), a1=at(:,1)
  ! -volume0 (Bohr)
  ! -bg (2pi/alat)
  ! -ecutwfc (Hartree)
  ! -ecutrho (Hartree)
  ! -nr1, nr2, nr3 (FFT_GRID)
  ! -noncolin
  ! -lspinorb
  ! -spinorb_mag
  !   if (noncolin) then
  !     nspin = 2
  !   else
  !     nspin = 1
  !   end if
  ! -nat
  ! -ntyp
  ! allocate(atom_labels(ntyp))
  ! allocate(atom_pfile(ntyp))
  ! allocate(ityp(nat))
  ! allocate(tau(3,nat)) (alat)
  ! allocate(amass(ntyp)) (amu)
  ! -nsym
  ! allocate(s(3,3,nsym))
  ! allocate(ftau(3,nsym))
  ! allocate(can_use_TR(nsym))
  ! -nkpoints_QE
  ! -nG_max (MAX_NUMBER_OF_GK-VECTORS)
  ! -nbands

  call read_parameters_data_file_xml()
  !
  alat_sc = alat
  at_sc = at ! a1=at(:,1) (alat)
  bg_sc = bg ! (2pi/alat)
  nr1_sc=nr1
  nr2_sc=nr2
  nr3_sc=nr3
  nat_sc=nat
  ntyp_sc=ntyp
  allocate(atom_labels_sc(ntyp_sc))
  atom_labels_sc=atom_labels
  allocate(ityp_sc(nat_sc))
  ityp_sc=ityp
  allocate(tau_sc(3,nat_sc)) ! alat_sc
  tau_sc=tau
  allocate(amass_sc(ntyp_sc))
  amass_sc=amass
  nkpoints_QE_sc=nkpoints_QE

  allocate(kpoints_QE_sc(3,nkpoints_QE_sc))
  call read_kpoints_data_file_xml(kpoints_QE_sc)
  !
  allocate(kpoints_irr(3,nk1*nk2*nk3))
  call find_the_irreducible_k_set(nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
  !
  !
  if ( nkmesh == nkpoints_QE_sc ) then
    full_zone = .true.
    write(stdout,20) 'The QE calculation contains all k points of the mesh'
  else
    full_zone = .false.
    if ( nk_irr == nkpoints_QE_sc ) then
      write(stdout,20) 'The QE calculation contains all irreaducible k points of the mesh'
    else
      write(stdout,20) '*****************************************************'
      write(stdout,20) '*      The number of kpoints present in the QE      *'
      write(stdout,20) '*      folders are not consistent with a full       *'
      write(stdout,20) '*      Brillouin Zone or an irreducible Brillouin   *'
      write(stdout,20) '*      zone! Review your input...                   *'
      write(stdout,20) '*                   Program stops.                  *'
      write(stdout,20) '*****************************************************'
      write(stdout,20) '* debug information:                                *'
      write(stdout,21) 'nkpoints_QE = ',nkpoints_QE_sc
      write(stdout,21) 'nkmesh      = ',nkmesh
      write(stdout,21) 'nk_irr      = ',nk_irr
      write(stdout,20) '*****************************************************'
      21 FORMAT('*',a20,i4,T53,'*')
      stop
    end if
  end if
  !
  !
!! BEGIN_DEBUG
!! Print information readed from QE
#ifdef DEBUG
  if (.false.) then
    write(stdout,2100)     'nat_sc:', nat_sc
    write(stdout,2100)     'ntyp_sc:', ntyp_sc
    write(stdout,2300)     'alat_sc:', alat_sc
    write(stdout,20) ''
    do id=1,3
      write(stdout,2400)   'at_sc:', at_sc(:,id)
    enddo
    write(stdout,20) ''
    do id=1,3
      write(stdout,2400)   'bg_sc:', bg_sc(:,id)
    enddo
    write(stdout,2200)     'nr1_sc, nr2_sc, nr3_sc:', nr1_sc, nr2_sc, nr3_sc
    write(stdout,20) ''
    do ia=1,nat
      write(stdout,2400)   'tau_sc:', tau_sc(:,ia)
    enddo
    write(stdout,20) ''
    do it=1,ntyp
      write(stdout,2400)   'amass_sc:', amass_sc
    enddo
  endif
#endif
!! END_DEBUG


  !
  !
  mesh_dir='./'
  ph_dir='PH/221/'
  fc_mat='g.fc'
  data_dir='data_dir/'
  data_name='data-file.'
  dvscf_dir='dvscf_dir/'
  dvscf_name='g.dvscf_q'
  qlist='irr-qpoints.txt'
  !
  ! If the SC can be created by a repetition of a simpler UC
  if (sc_from_uc) then
    !
    nq1   = 2
    nq2   = 2
    nq3   = 1
    nqirr = 2
    nat_uc = 2
    !
    ! read phonon related data:
    ! -nmodes = 3*nat
    ! -allocate(q_irr(3,nqirr)) (2pi/alat)
    ! -allocate (u_irr(nmodes,nmodes,nqirr))

    nat = nat_uc
    call read_ph_information_xml() ! To read q_irr
    !
    !
  else
    !
    nq1   = 1
    nq2   = 1
    nq3   = 1
    nqirr = 1
    nat_uc = nat_sc
    !
  endif ! sc_from_uc
  !
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !
  allocate(qpoints_irr(3,nq1*nq2*nq3))
  call find_the_irreducible_k_set(nq1,nq2,nq3,qmesh,qpoints_irr,nq_irr)


!! BEGIN_DEBUG
!! Print information about phonons readed from QE
#ifdef DEBUG
  if (.false.) then
    write(stdout,*) 'sc_from_uc', sc_from_uc
    do i=1,nqmesh
      write(stdout,'(i2,3f10.5)') i, qmesh(:,i)
    enddo
    write(stdout,2100) 'nqirr:', nq_irr
    do i=1,nq_irr
      write(stdout,'(i2,3f10.5)') i, qpoints_irr(:,i)
    enddo
    write(stdout,2200)     'nq1, nq2, nq3:', nq1, nq2, nq3
    write(stdout,2100)     'nqirr:', nqirr
  endif
#endif
!! END_DEBUG


  !
  !
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Calculate the FC matrix'
  write(stdout,20) ''
  !
  ! Calculate the FC matrix
  allocate( frc_0(nq1,nq2,nq3,3,3,nat_uc,nat_uc) )
  allocate(f_m(3,nat_sc),f_p(3,nat_sc))
  frc_0 = cmplx_0
  f_m = zero
  f_p = zero
  !
  disp(1) = 'x'
  disp(2) = 'y'
  disp(3) = 'z'
  dx = 0.0001_dp * alat_sc
  !
  mesh_dir = 'SCF/'
  !
  ! put tau in units of the graphene cell parameter
  if ( nq1==nq2 ) then
    tau_sc = tau_sc * nq1
  else
    write(stdout,20) 'nq1 and nq2 must be equal'
    stop
  endif
  !
#ifdef DEBUG
  if (.false.) then
    do ia=1,nat_sc
      write(stdout,2400)   'tau_sc:', tau_sc(:,ia)
    enddo
  endif
#endif
  !
  do ia=1,nat_sc ! for each atom in the unit cell
    do id=1,3 ! for each direction
      !
      ! Read forces for + and - displacaments
      write(dummy1,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'+'
      scf_file_p = trim(mesh_dir)//'/'//trim(dummy1)//'/'//trim(dummy1)//'.scf.out'
      write(dummy2,'(a,i1,a,a)') trim(prefix), ia, trim(adjustl(disp(id))),'-'
      scf_file_m = trim(mesh_dir)//'/'//trim(dummy2)//'/'//trim(dummy2)//'.scf.out'
      call read_f(scf_file_p,nat_sc,f_p)
      call read_f(scf_file_m,nat_sc,f_m)
      !
      ! equivalent atom of the original cell
      ia_original = mod(ia-1,nat_uc) + 1
      R_m(:) = (tau_sc(:,ia) - tau_sc(:,ia_original))
      !
      do ja=1,nat_sc
        do jd=1,3
          !
          ! equivalent atom of the original cell
          ja_original = mod(ja-1,nat_uc) + 1
          R_l(:) = (tau(:,ja) - tau(:,ja_original))
          !
          ! R = R_l - R_m
          R = R_l - R_m
          R = matmul(transpose(bg),R)
          !
          ! find the FC related to C_ij(R)
          m1 = mod(nint(R(1))+1,nq1)
          if(m1.le.0) m1=m1+nq1
          m2 = mod(nint(R(2))+1,nq2)
          if(m2.le.0) m2=m2+nq2
          m3 = mod(nint(R(3))+1,nq3)
          if(m3.le.0) m3=m3+nq3
          !
          ! Calculate the FC matrix with the derivative of force
          frc_0(m1,m2,m3,id,jd,ia_original,ja_original) = -( f_p(jd,ja) - f_m(jd,ja) ) / ( 2.0_dp * dx )
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
  ! calculate the dynamical matrix
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Calculate the dynamical matrix'
  write(stdout,20) ''
  !
  allocate( dyn_q(3*nat_uc,3*nat_uc), w2(3*nat_uc) )
  allocate(z(3,nat_sc))
  !
  dyn_q=cmplx_0
  w2=zero
  !
  open(unit=112, file='vib.xyz')
  !
  do i=1,nqmesh
    !
    call mat_inv_four_t(qmesh(:,i),nq1,nq2,nq3,3*nat_uc,frc_0,dyn_q(:,:))
    !
    ! diagonalize the dynamical matrix to obtain eigenmodes
    call diagonalize_cmat( 3*nat_uc, dyn_q, w2 )
    !
    write(stdout,'(a,3f12.6)') 'q=',  qmesh(:,i)
    !
    do imode=1,3*nat_uc
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
      write(stdout,'(a,i3,2(f12.6,x,a3))') 'mode', imode, w2(imode), 'meV', w2(imode)*0.24180, 'THz'
      !
      ! save polarization
      write(112,*) nat_sc
      write(112,'(a6,i4,a5,3f8.4,a5,1f10.6,a5)') 'Mode #',imode, ' q = ', qmesh(:,i), ' f = ', w2(imode),' meV.'
      do ia = 1,nat_sc
        ia_original = mod(ia-1,nat_uc) + 1
        R(:) = (tau_sc(:,ia) - tau_sc(:,ia_original))
        R = matmul(transpose(bg),R)
        print*, twopi*cmplx_i*( qmesh(1,i)*R(1) + qmesh(2,i)*R(2) + qmesh(3,i)*R(3) )
        do id=1,3
          z(id,ia) = dyn_q(3*(ia_original-1)+id,imode)*sqrt(real(nat_uc)/nat_sc)
          z(id,ia) = z(id,ia)*exp( twopi*cmplx_i*( qmesh(1,i)*R(1) + qmesh(2,i)*R(2) + qmesh(3,i)*R(3) ) )
        enddo
        write(112,'(a4,x,3(2f10.5))') atom_labels_sc(ityp_sc(ia)), ( tau_sc(id,ia)/nq1*alat_sc*b2a ,id=1,3 ), (real(z(id,ia)), id=1,3)
      enddo
      !
    enddo ! imode
    !
  enddo ! iq
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
  allocate( delta_V(nr1*nr2*nr3,3*nat_sc) )
  allocate( rho(2,nr1,nr2) )

  do ia=1,nat_sc ! for each atom in the unit cell
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
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_sc, ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_total_p, - 1)
      CALL plot_io (trim(scf_file_m), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_sc, ntyp, ibrav, celldm, at,&
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
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_sc, ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_ionic_p, - 1)
      CALL plot_io (trim(scf_file_m), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
            dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_sc, ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_ionic_m, - 1)

      delta_V(:,imode) = ( (V_total_p - V_ionic_p) - (V_total_m - V_ionic_m) ) / ( 2.0_dp * dx )
      !
      if (imode<10) write(xsf_file,'(a17,i1)') 'fe.xsf_dvscf_mode', imode
      if (imode>9) write(xsf_file,'(a17,i2)') 'fe.xsf_dvscf_mode', imode
      ounit = find_free_unit()
      open(unit=ounit, file=xsf_file, status="replace", action="write")
      call xsf_struct (alat, at, nat_sc, tau, atom_labels, ityp, ounit)
      x0=0.0
      x0(3) = 2.5/2.0
      e1 = at(:,1)
      e2 = at(:,2)
      do ir=1,nr1*nr2*nr3
        call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)
        if (k==nr3/2) then
          rho(1,i,j) = real ( delta_V(ir, imode) )
          rho(2,i,j) = 0.0
        endif
      enddo
      ! print*, maxval(abs(rho(1,:,:)))/ maxval(abs(rho(2,:,:)))
      call xsf_datagrid_2d (rho, nr1, nr2, 1.0_dp, 1.0_dp, x0, e1, e2, alat, ounit)
      close(unit=ounit)

      !
    enddo
  enddo



  ! This subroutine does this:
  !   read (iunplot, '(a)') title
  !   read (iunplot, * ) nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
  !   read (iunplot, * ) ibrav, celldm
  !   if (ibrav == 0) then
  !     do i = 1,3
  !       read ( iunplot, * ) ( at(ipol,i),ipol=1,3 )
  !     enddo
  !   endif
  !   read (iunplot, * ) gcutm, dual, ecut, plot_num
  !   read (iunplot, '(i4,3x,a2,3x,f5.2)') &
  !         (ndum, atm(nt), zv(nt), nt=1, ntyp)
  !   read (iunplot, *) (ndum,  (tau (ipol, na), ipol = 1, 3), &
  !         ityp(na), na = 1, nat)
  !   read (iunplot, * ) (plot (ir), ir = 1, nr1x * nr2x * nr3)






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
