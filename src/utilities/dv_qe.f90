program dv_hand
  !
  ! This program reads the data of the DFPT calculation of QE to calculate the
  ! potential change for each atom of the supercell.
  !
  use kinds, only: dp
  use intw_input_parameters, only: mesh_dir, prefix, ph_dir, dvscf_dir, &
                                   dvscf_name, data_dir, data_name, &
                                   nk1, nk2, nk3, fc_mat, &
                                   qlist, nq1, nq2, nq3, nqirr
  use intw_reading, only: at, bg, alat, nat, ntyp, tau, amass, nr1, nr2, nr3, &
                          nkpoints_QE, nsym, npol, atom_labels, ityp
  use intw_reading, only: write_tag, read_parameters_data_file_xml, &
                          read_kpoints_data_file_xml
  use intw_ph, only: frc, q_irr_cryst, QE_folder_sym_q, sym_G_q, symlink_q, &
                     read_allq_dvr, read_ph_information_xml, readfc, &
                     mat_inv_four_t, dvq_local, get_dv, q_irr

  use intw_utility, only: find_free_unit, generate_kmesh, switch_indices_zyx

  use intw_fft, only: wfc_by_expigr

  use intw_symmetries, only: find_the_irreducible_k_set, &
                             find_inverse_symmetry_matrices_indices, &
                             allocate_symmetry_related_k, &
                             allocate_and_build_spin_symmetry_matrices, &
                             set_symmetry_relations, rot_atoms, &
                             rtau_index, rtau, tau_cryst, rtau_cryst

  use intw_useful_constants
  use iotk_module
  !
  implicit none
  !
  real(kind=dp), parameter :: pmass=1822.88848426_dp, aumev=  27211.396_dp
  !
  complex(kind=dp), allocatable, dimension(:,:,:) :: dyn_q ! dynamical matrix
  real(kind=dp), allocatable, dimension(:,:) :: w2 ! frequencies
  !
  character(len=256) :: fc_file_name, xsf_file
  !
  integer :: ia, it, id, i, j, k, irr1, irr2, irr3, ia_sc, idir, imode_sc
  !
  integer, parameter :: stdout=6
  ! character(len=1) :: CR=CHAR(13)

  !
  logical :: full_zone
  integer :: nkmesh, nqmesh, nk_irr, nq_irr
  real(kind=dp), allocatable, dimension(:,:) :: kmesh, qmesh, &
                                                kpoints_irr, qpoints_irr, &
                                                kpoints_QE
  !
  complex(kind=dp), allocatable, dimension(:,:,:,:) :: dv_sc
  real(kind=dp), allocatable, dimension(:,:,:) :: rho
  !
  ! Supercell variables
  integer :: nr1_sc, nr2_sc, nr3_sc
  integer :: ir1, ir2, ir3, ir, iq, imode
  !
  real(kind=dp), dimension(3) :: qpoint, x0, e1, e2
  !
  logical                  :: q_points_consistent
  !
  integer                  :: nr(3)
  !
  real(kind=dp) :: qrr, qr, m
  !
  integer :: ounit
  !
  ! animation
  real(kind=dp) :: t, dt
  integer :: N, times, ix, iy
  real(kind=dp) :: xx, yy ,zz
  real(kind=dp), dimension(3) :: fi

  ! wigner_seitz
  integer :: nrr_k
  integer :: nrpts
  integer, allocatable :: ndegen(:), irvec(:,:)


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
  write(stdout,20) '|                  program dv_qe.x                   |'
  write(stdout,20) '|        ------------------------------------        |'
  write(stdout,20) '======================================================'
  !
  !================================================================================
  !
  ! Directory and files
  mesh_dir = 'PH/'
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
  write(*,20) 'Reading the QE calculation data'
  !
  ! Read data-file
  call read_parameters_data_file_xml()
  !
  allocate(kpoints_QE(3,nkpoints_QE))
  call read_kpoints_data_file_xml(kpoints_QE)
  !
  allocate(kpoints_irr(3,nk1*nk2*nk3))
  call find_the_irreducible_k_set(nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
  !
  !
  if ( nkmesh == nkpoints_QE ) then
    full_zone = .true.
    write(stdout,20) 'The QE calculation contains all k points of the mesh'
  else
    full_zone = .false.
    if ( nk_irr == nkpoints_QE ) then
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
      write(stdout,21) 'nkpoints_QE = ',nkpoints_QE
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
    write(*,2100)     'nat:', nat
    write(*,2100)     'ntyp:', ntyp
    write(*,2300)     'alat:', alat
    write(*,20) ''
    do id=1,3
      write(*,2400)   'at:', at(:,id)
    enddo
    write(*,20) ''
    do id=1,3
      write(*,2400)   'bg:', bg(:,id)
    enddo
    write(*,2200)     'nr1, nr2, nr3:', nr1, nr2, nr3
    write(*,20) ''
    do ia=1,nat
      write(*,2400)   'tau:', tau(:,ia)
    enddo
    write(*,20) ''
    do it=1,ntyp
      write(*,2400)   'amass:', amass
    enddo
  endif
#endif
!! END_DEBUG

  !
  !
  !
  !
  ph_dir='221/'
  fc_mat='g.fc'
  data_dir='data_dir/'
  data_name='data-file.'
  dvscf_dir='dvscf_dir/'
  dvscf_name='g.dvscf_q'
  qlist='irr-qpoints.txt'
  nq1   = 2
  nq2   = 2
  nq3   = 1
  nqirr = 2

  !
  ! read phonon related data
  call read_ph_information_xml()
  !
  !
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !
  allocate(qpoints_irr(3,nqmesh))
  qpoints_irr = 0.0_dp
  call find_the_irreducible_k_set(nq1,nq2,nq3,qmesh,qpoints_irr,nq_irr)
  !
  !


!! BEGIN_DEBUG
!! Print information about phonons readed from QE
#ifdef DEBUG
  if (.false.) then
    do i=1,nqmesh
      write(*,'(i2,3f10.5)') i, qmesh(:,i)
    enddo
    write(*,2100) 'nqirr:', nq_irr
    do i=1,nq_irr
      write(*,'(i2,3f10.5)') i, qpoints_irr(:,i)
    enddo
    write(*,2200)     'nq1, nq2, nq3:', nq1, nq2, nq3
    write(*,2100)     'nqirr:', nqirr
  endif
#endif
!! END_DEBUG


  !
  !
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Read the FC matrix'
  write(stdout,20) ''
  !
  fc_file_name=trim(trim(mesh_dir)//"/"//trim(ph_dir)//"/"//trim(fc_mat))
  !
  call readfc(fc_file_name,nq1,nq2,nq3,nat,alat,at,ntyp,amass)
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
  allocate( dyn_q(3*nat,3*nat,nq_irr), w2(3*nat,nq_irr) )
  dyn_q=cmplx_0
  w2=zero
  !
  !
  !
  do iq=1,nq_irr
    !
    call mat_inv_four_t(qpoints_irr(:,iq),nq1,nq2,nq3,3*nat,frc,dyn_q(:,:,iq))
    !
    ! diagonalize the dynamical matrix to obtain eigenmodes
    call diagonalize_cmat( 3*nat, dyn_q(:,:,iq), w2(:,iq) )
    !
    write(*,'(a,3f12.6)') 'q=',  qpoints_irr(:,iq)
    !
    do imode=1,3*nat
      !
      ! Normalize the eigenvector
      dyn_q(:,imode,iq) = dyn_q(:,imode,iq) / sqrt(sum(dyn_q(:,imode,iq)**2.0_dp))
      !
      ! Avoid imaginary frequencies
      if (w2(imode,iq)<zero) then
        w2(imode,iq)=-sqrt(abs(w2(imode,iq)))
      else
        w2(imode,iq)=sqrt(w2(imode,iq))
      endif
      !
      write(*,'(a,i3,2(f12.6,x,a3))') 'mode', imode, w2(imode,iq), 'meV', w2(imode,iq)*0.24180, 'THz'
      do k=1,nat
        write(*,'(a,i3,3(2f10.5))') 'atom', k, ( dyn_q(3*(k-1)+i,imode,iq), i=1, 3)
      enddo
      !
    enddo
    !
    !
  enddo
  !
  !


!! BEGIN_DEBUG
!! Save the polarization of the modes on a file to visualize with Jmol
#ifdef DEBUG
  !
  if (.false.) then
    !
    open(unit=111, file='mode_polarizations.xyz')
    qpoint(1)=0.098
    qpoint(2)=0.0
    qpoint(3)=0.0
    iq=1
    call mat_inv_four_t(qpoint,nq1,nq2,nq3,3*nat,frc,dyn_q(:,:,iq))
    call diagonalize_cmat( 3*nat, dyn_q(:,:,iq), w2(:,iq) )
    do imode=1,3*nat
      dyn_q(:,imode,iq) = dyn_q(:,imode,iq) / sqrt(sum(dyn_q(:,imode,iq)**2.0_dp))
    enddo
    imode=1
    write(*,'(a,3f15.6)') 'q=', qpoint
    write(*,'(a,i3,f15.6)') 'mode, freq (THz)=', 1, w2(imode,iq)*0.24180
    do k=1,nat
      write(*,'(a,i3,3(2f10.5))') 'polarization', k, ( dyn_q(3*(k-1)+i,imode,iq), i=1, 3)
    enddo
    9000 format(a,3 ('(',f15.6,',',f15.6,')' ) )
    9100 format (a,1x,3 ('(',f10.6,',',f10.6,1x,')',5x))

    N=80
    times=2
    dt=1.0_dp/(N*times-1)*times
    do i=1,N*times
      t=dt*(i-1)
      write(111,*) nat*10*10
      write(111,*) t
      do ix=1,20
        do iy=1,5
          do ia=1,nat
            do id=1,3
              fi(id) = 0.4*real( dyn_q(3*(ia-1)+id,imode,iq) &
                        * exp( cmplx_i*2.0_dp*pi*2.4*t ) &
                        * exp( cmplx_i * 2.0_dp*pi * ( qpoint(1)*ix + qpoint(2)*iy ) ) )
            enddo
            xx = ( tau(1,ia) + ix*at(1,1) + iy*at(1,2) )*alat*bohr + fi(1)
            yy = ( tau(2,ia) + ix*at(2,1) + iy*at(2,2) )*alat*bohr + fi(2)
            zz =   tau(3,ia)                            *alat*bohr + fi(3)
            write(111,*) atom_labels(ityp(ia)), xx, yy, zz
          end do
        enddo
      enddo
    end do

    close(111)
    !
    ! stop
    !
  endif
  !
#endif
!! END_DEBUG

  !
  ! Calculate the potential change dV
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Read the potential change'
  write(stdout,20) ''
  !
  call read_allq_dvr(nqirr,3*nat)
  !
  !
  ! dvscf_cart(nr1*nr2*nr3,nqirr,3*nat,1:npol,1:npol) daukeu ia
  !
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Calculate the potential change in the SC'
  write(stdout,20) ''
  !
  !
  allocate(rtau_index(nat,nsym))
  allocate(rtau(3,nsym,nat))
  allocate(tau_cryst(3,nat))
  allocate(rtau_cryst(3,nsym,nat))
  !
  nr=(/nr1,nr2,nr3/)
  !
  call rot_atoms(nat,nsym,tau)
  !
  !================================================================================
  !      allocate the symmetry arrays
  !      CAREFUL! the subroutine needs to know the global value of "full_mesh",
  !      so it is crucial that this allocation occurs AFTER setting full_mesh
  !================================================================================
  !
  call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)
  !
  !================================================================================
  !     Compute the indices of the inverse rotation matrices.
  !     This will be useful later in the exectution.
  !================================================================================
  !
  call find_inverse_symmetry_matrices_indices()
  !
  call allocate_and_build_spin_symmetry_matrices(nsym)
  !
  allocate(QE_folder_sym_q(nqmesh))
  allocate(q_irr_cryst(3,nqirr))
  allocate(sym_G_q(3,nqmesh))
  allocate(symlink_q(nqmesh,2))
  !
  do iq=1,nqirr
    q_irr_cryst(:,iq)=matmul(transpose(at),q_irr(:,iq))
    print*, iq, q_irr_cryst(:,iq)
  enddo
  !
  !================================================================================
  ! We need the symmtry relations between irreducible/full in order to obtain
  ! the induced potential for any q point.
  !================================================================================
  !
  call set_symmetry_relations(nq1,nq2,nq3,nqirr,q_irr_cryst,qmesh, &
              q_points_consistent,QE_folder_sym_q,sym_G_q,symlink_q)
  !
  !
  nr1_sc = nr1 * nq1
  nr2_sc = nr2 * nq2
  nr3_sc = nr3 * nq3
  !
  nrr_k=0
  call hamiltonian_wigner_seitz(.true.,nq1,nq2,nq3)
  print*, nrpts
  allocate(ndegen(nrpts), irvec(3,nrpts))
  call hamiltonian_wigner_seitz(.false.,nq1,nq2,nq3)
  do i=1,nrpts
    print*, irvec(:,i), ndegen(i)
  enddo
  !
  !
  allocate( dv_sc(nr1_sc,nr2_sc,nr3_sc, 3*nat*nq1*nq2*nq3) )
  allocate( rho(2,nr1_sc,nr2_sc) )
  allocate( dvq_local(nr1*nr2*nr3,3*nat,npol,npol) )
  !
  dv_sc = cmplx_0
  rho = cmplx_0
  !
  !


!! BEGIN_DEBUG
!! Calculate delta_v (the periodic part) for G/2 and -G/2
#ifdef DEBUG
  !
  if (.false.) then
    !
    ! delta_v_q
    qpoint = (/ +0.5, +0.5, 0.0 /)
    print*, qpoint
    call get_dv(1,qpoint,3*nat,dvq_local)
    do i=1,nr1*nr2*nr3
      write(123,'(2f12.7)') dvq_local(i,1,1,1)
    enddo
    !
    ! delta_v_-q
    qpoint = (/ -0.5, -0.5, 0.0 /)
    print*, qpoint
    call get_dv(1,qpoint,3*nat,dvq_local)
    do i=1,nr1*nr2*nr3
      write(124,'(2f12.7)') dvq_local(i,1,1,1)
    enddo
    !
    ! delta_v_-q * e^iqr
    do ir3=1,nr3_sc
      do ir2=1,nr2_sc
        do ir1=1,nr1_sc
          i=mod(ir1-1,nr1)+1
          j=mod(ir2-1,nr2)+1
          k=mod(ir3-1,nr3)+1
          call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,+1)
          qr=tpi*(qpoint(1)*real(ir1-1,dp)/nr1 + qpoint(2)*real(ir2-1,dp)/nr2 + qpoint(3)*real(ir3-1,dp)/nr3)
          write(125,'(2f12.7)') dvq_local(ir, 1,1,1) * exp(cmplx_i*qr)
        enddo
      enddo
    enddo
    !
    stop
    !
  endif
  !
#endif
!! END_DEBUG





!! BEGIN_DEBUG
!! Calculate delta_v_q * e^iqr for each q and save on a file
#ifdef DEBUG
  !
  if (.false.) then
    do iq=1,nrpts
      !
      qpoint(1) = real(irvec(1,iq),dp)/nq1
      qpoint(2) = real(irvec(2,iq),dp)/nq2
      qpoint(3) = real(irvec(3,iq),dp)/nq3
      m=1.0_dp/(ndegen(iq)*nqmesh)
      !
      write(*,'(a,i3,3f12.6,2x,f12.6)') 'q=', iq, qpoint(:),m
      !
      dvq_local=cmplx_0
      call get_dv(iq,qpoint,3*nat,dvq_local)
      !
      do ir3=1,nr3
        do ir2=1,nr2
          do ir1=1,nr1
            !
            i=mod(ir1-1,nr1)+1
            j=mod(ir2-1,nr2)+1
            k=mod(ir3-1,nr3)+1
            call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,+1)
            qr=tpi*(qpoint(1)*real(ir1-1,dp)/nr1 + qpoint(2)*real(ir2-1,dp)/nr2 + qpoint(3)*real(ir3-1,dp)/nr3)
            write(1000+iq,'(2f12.7)') dvq_local(ir, 1,1,1) * exp(cmplx_i*qr)
            !
          enddo
        enddo
      enddo
      !
    enddo
    ! stop
  endif

#endif
!! END_DEBUG



  ! for each cell of the supercell
  do irr1=1,nq1
    do irr2=1,nq2
      do irr3=1,nq3
        !
        print*, 'R=', irr1, irr2, irr3
        !
        ! Sum in q
        do iq=1,nqmesh
          qpoint = qmesh(:,iq)
          m=1.0_dp/(nqmesh)
        ! do iq=1,nrpts
        !   qpoint(1) = real(irvec(1,iq),dp)/nq1
        !   qpoint(2) = real(irvec(2,iq),dp)/nq2
        !   qpoint(3) = real(irvec(3,iq),dp)/nq3
        !   m=1.0_dp/(ndegen(iq)*nqmesh)
          !
          write(*,'(a,i3,3f12.6,2x,f12.6)') 'q=', iq, qpoint(:),m
          !
          qrr = tpi*(qpoint(1)*(irr1-1) + qpoint(2)*(irr2-1) + qpoint(3)*(irr3-1))
          !
          dvq_local=cmplx_0
          call get_dv(iq,qpoint,3*nat,dvq_local)
          !
          do ir3=1,nr3_sc
            do ir2=1,nr2_sc
              do ir1=1,nr1_sc
                !
                ! Calculate the equivalent r point in the Unit Cell
                i=mod(ir1-1,nr1)+1
                j=mod(ir2-1,nr2)+1
                k=mod(ir3-1,nr3)+1
                !
                call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,+1)
                !
                ! Calculate the phase
                qr=tpi*(qpoint(1)*real(ir1-1,dp)/nr1 + qpoint(2)*real(ir2-1,dp)/nr2 + qpoint(3)*real(ir3-1,dp)/nr3)
                !
                do imode=1,3*nat
                  !
                  ia = (imode-1)/3 + 1
                  idir = mod(imode-1,3)+1
                  ia_sc = ia + nat*(nq2*(irr1-1)+(irr2-1)+nq1*nq2*(irr3-1))
                  ! The reason of this orden (irr2 first and then irr1) is that when I
                  ! create the suppercell increase the Y direction first, and then
                  ! the X direction
                  imode_sc = idir + 3*(ia_sc-1)
                  !
                  ! if (iq==1.or.iq==3.or.iq==8.or.iq==12.or.iq==17.or.iq==19) then
                  !   dv_sc(ir1,ir2,ir3, imode_sc) = &
                  !         dv_sc(ir1,ir2,ir3, imode_sc) + &
                  !         m* ( real(dvq_local(ir, imode,1,1) * exp(cmplx_i*qr))* exp(-cmplx_i*qrr) )! + &
                  ! else
                  dv_sc(ir1,ir2,ir3, imode_sc) = &
                        dv_sc(ir1,ir2,ir3, imode_sc) + &
                        m* real( dvq_local(ir, imode,1,1) * exp(cmplx_i*qr)* exp(-cmplx_i*qrr) )! + &
                        ! m*conjg( dvq_local(ir, imode,1,1) * exp(cmplx_i*qr)* exp(-cmplx_i*qrr) )
                  ! endif
                  !
                enddo !imode
                !
              enddo !ir1
            enddo !ir2
          enddo !ir3
          !
        enddo !iq
        !
        do imode=1,3*nat
          ia = (imode-1)/3 + 1
          idir = mod(imode-1,3)+1
          ia_sc = ia + nat*( (irr2-1) + nq2*(irr1-1) + nq1*nq2*(irr3-1) )
          imode_sc = idir + 3*(ia_sc-1)
          ounit = find_free_unit()
          if (imode_sc<10) write(xsf_file,'(a13,3i1,a5,i1)') 'fe.xsf_real_R', irr1, irr2, irr3, '_mode', imode_sc
          if (imode_sc>9) write(xsf_file,'(a13,3i1,a5,i2)') 'fe.xsf_real_R', irr1, irr2, irr3, '_mode', imode_sc
          open(unit=ounit, file=xsf_file, status="replace", action="write")
          call xsf_struct (alat, at, nat, tau, atom_labels, ityp, ounit)
          x0=0.0
          x0(3) = 2.5
          e1 = at(:,1)*nq1
          e2 = at(:,2)*nq2
          do ir=1,nr1_sc*nr2_sc*nr3_sc
            call switch_indices_zyx(nr1_sc,nr2_sc,nr3_sc,ir,i,j,k,-1)
            if (k==nr3_sc/2) then
              rho(1,i,j) = real ( dv_sc(i,j,k, imode_sc) )
              rho(2,i,j) = aimag( dv_sc(i,j,k, imode_sc) )
            endif
          enddo
          ! print*, maxval(abs(rho(1,:,:)))/ maxval(abs(rho(2,:,:)))
          call xsf_datagrid_2d (rho, nr1_sc, nr2_sc, 1.0_dp, 1.0_dp, x0, e1, e2, alat, ounit)
          close(unit=ounit)
        enddo
        do imode=1,3*nat
          ia = (imode-1)/3 + 1
          idir = mod(imode-1,3)+1
          ia_sc = ia + nat*(nq2*(irr1-1)+(irr2-1)+nq1*nq2*(irr3-1))
          imode_sc = idir + 3*(ia_sc-1)
          ounit = find_free_unit()
          if (imode_sc<10) write(xsf_file,'(a13,3i1,a5,i1)') 'fe.xsf_imag_R', irr1, irr2, irr3, '_mode', imode_sc
          if (imode_sc>9) write(xsf_file,'(a13,3i1,a5,i2)') 'fe.xsf_imag_R', irr1, irr2, irr3, '_mode', imode_sc
          open(unit=ounit, file=xsf_file, status="replace", action="write")
          call xsf_struct (alat, at, nat, tau, atom_labels, ityp, ounit)
          x0=0.0
          x0(3) = 2.5
          e1 = at(:,1)*nq1
          e2 = at(:,2)*nq2
          do ir=1,nr1_sc*nr2_sc*nr3_sc
            call switch_indices_zyx(nr1_sc,nr2_sc,nr3_sc,ir,i,j,k,-1)
            if (k==nr3_sc/2) then
              rho(2,i,j) = real ( dv_sc(i,j,k, imode_sc) )
              rho(1,i,j) = aimag( dv_sc(i,j,k, imode_sc) )
            endif
          enddo
          call xsf_datagrid_2d (rho, nr1_sc, nr2_sc, 1.0_dp, 1.0_dp, x0, e1, e2, alat, ounit)
          close(unit=ounit)
        enddo
        !
      enddo
    enddo
  enddo

  ! fasea = dv_sc(3,1,1,1)/abs(dv_sc(3,1,1,1))
  ! do irr1=1,nr1_sc
  !   print'(a1,2f15.10,a1,x,a1,2f15.10,a1)', '(',dv_sc(irr1,1,1,1),')','(', dv_sc(irr1,1,1,1)/fasea, ')'
  ! enddo

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

  subroutine hamiltonian_wigner_seitz(count_pts,nq1,nq2,nq3)
      !================================================================================!
      ! Calculates a grid of points that fall inside of (and eventually on the         !
      ! surface of) the Wigner-Seitz supercell centered on the origin of the B         !
      ! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3  !
      !================================================================================!

      use w90_constants,  only : eps7,eps8
      use w90_io,         only : io_error,io_stopwatch,stdout
      use w90_parameters, only : iprint,mp_grid,real_metric,timing_level

      ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
      !                   irvec(1:3,irpt) in the basis of the lattice vectors
      ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
      ! nrpts             number of Wigner-Seitz grid points

      implicit none

      logical, intent(in) :: count_pts
      integer, intent(in) :: nq1, nq2, nq3

      integer       :: ndiff (3)
      real(kind=dp) :: dist(125),tot,dist_min
      integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j

      if (timing_level>1) call io_stopwatch('hamiltonian: wigner_seitz',1)

      ! The Wannier functions live in a supercell of the real space unit cell
      ! this supercell is mp_grid unit cells long in each direction
      !
      ! We loop over grid points r on a unit cell that is 8 times larger than this
      ! primitive supercell.
      !
      ! One of these points is in the W-S cell if it is closer to R=0 than any of the
      ! other points, R (where R are the translation vectors of the supercell)

      ! In the end nrpts contains the total number of grid
      ! points that have been found in the Wigner-Seitz cell

      nrpts = 0
      do n1 = -nq1 , nq1
        do n2 = -nq2, nq2
          do n3 = -nq3, nq3
            ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=1, or icnt=14
            icnt = 0
            do i1 = -2, 2
              do i2 = -2, 2
                do i3 = -2, 2
                  icnt = icnt + 1
                  ! Calculate distance squared |r-R|^2
                  ndiff(1) = n1 - i1 * nq1
                  ndiff(2) = n2 - i2 * nq2
                  ndiff(3) = n3 - i3 * nq3
                  dist(icnt) = 0.0_dp
                  do i = 1, 3
                    do j = 1, 3
                      dist(icnt) = dist(icnt) + real(ndiff(i),dp) * at(i,j) &
                            * real(ndiff(j),dp)
                    enddo
                  enddo
                enddo
              enddo

              ! AAM: On first pass, we reference unallocated variables (ndegen,irvec)

            enddo
            dist_min=minval(dist)
            if (abs(dist(63) - dist_min ) .lt. eps7 ) then
              nrpts = nrpts + 1
              if(.not. count_pts) then
                ndegen(nrpts)=0
                do i=1,125
                  if (abs (dist (i) - dist_min) .lt. eps7 ) ndegen(nrpts)=ndegen(nrpts)+1
                end do
                irvec(1, nrpts) = n1
                irvec(2, nrpts) = n2
                irvec(3, nrpts) = n3
              endif
            end if

            !n3
          enddo
          !n2
        enddo
        !n1
      enddo
      !
      if(count_pts) return


      if(iprint>=3) then
        write(stdout,'(1x,i4,a,/)') nrpts,  ' lattice points in Wigner-Seitz supercell:'
        do i=1,nrpts
          write(stdout,'(4x,a,3(i3,1x),a,i2)') '  vector ', irvec(1,i),irvec(2,i),&
            irvec(3,i),'  degeneracy: ', ndegen(i)
        enddo
      endif
      ! Check the "sum rule"
      tot = 0.0_dp
      do i = 1, nrpts
        tot = tot + 1.0_dp/real(ndegen(i),dp)
      enddo
      if (abs (tot - real(nq1 * nq2 * nq3,dp) ) > eps8) then
        print*, tot
        write(*,*) 'ERROR in hamiltonian_wigner_seitz: error in finding Wigner-Seitz points'
      endif

      if (timing_level>1) call io_stopwatch('hamiltonian: wigner_seitz',2)

      return

  end subroutine hamiltonian_wigner_seitz


end program dv_hand
