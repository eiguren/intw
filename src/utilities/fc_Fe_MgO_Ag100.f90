program fc
  !
  ! This program reads the 4x4x1 q-mesh FC matrix of the clean MgO/Ag(100) surface
  ! and comletes the 1x1x1 q-mesh FC matrix of the Fe/MgO/Ag(100) 4x4x1 supercell.
  !
  ! Import modules
  use kinds, only: dp
  !
  use intw_input_parameters, only: mesh_dir, prefix, ph_dir, dvscf_dir, &
                                   dvscf_name, data_dir, data_name, &
                                   qlist, fc_mat, nq1, nq2, nq3, nqirr
  !
  use intw_reading, only: at, bg, alat, nat, ntyp, tau, amass, nr1, nr2, nr3, &
                          atom_labels, atom_pfile, ityp, s, ftau, can_use_TR, &
                          read_parameters_data_file_xml, nsym
  !
  use intw_ph, only: frc, read_ph_information_xml, readfc, mat_inv_four_t
  !
  use intw_useful_constants
  !
  use intw_utility, only: ainv, generate_kmesh
  !
  use intw_symmetries, only: rtau_index, rot_atoms, tau_cryst, rtau_cryst, rtau
  !
  !
  implicit none
  !
  !
  ! Declare variables
  !
  !
  ! MgO-Ag(100) surface
  character(len=256)                                      :: mesh_dir_mgo
  character(len=256)                                      :: prefix_mgo
  character(len=256)                                      :: ph_dir_mgo
  character(len=256)                                      :: dvscf_dir_mgo
  character(len=256)                                      :: dvscf_name_mgo
  character(len=256)                                      :: data_dir_mgo
  character(len=256)                                      :: data_name_mgo
  character(len=256)                                      :: qlist_mgo
  character(len=256)                                      :: fc_mat_mgo
  integer                                                 :: nat_mgo
  integer                                                 :: ntyp_mgo
  real(kind=dp)                                           :: alat_mgo
  real(kind=dp), dimension(3,3)                           :: at_mgo
  real(kind=dp), dimension(3,3)                           :: bg_mgo
  integer                                                 :: nr1_mgo
  integer                                                 :: nr2_mgo
  integer                                                 :: nr3_mgo
  real(kind=dp), allocatable, dimension(:,:)              :: tau_mgo
  real(kind=dp), allocatable, dimension(:)                :: amass_mgo
  integer, allocatable, dimension(:)                      :: ityp_mgo
  character(len=3), allocatable, dimension(:)             :: atom_labels_mgo
  character(len=256)                                      :: fc_file_name_mgo
  integer                                                 :: nq1_mgo
  integer                                                 :: nq2_mgo
  integer                                                 :: nq3_mgo
  integer                                                 :: nqirr_mgo
  integer                                                 :: nqmesh_mgo
  complex(kind=dp), allocatable, dimension(:,:,:,:,:,:,:) :: frc_mgo
  !
  ! Fe-MgO-Ag(100) surface
  character(len=256)                             :: mesh_dir_fe
  character(len=256)                             :: prefix_fe
  character(len=256)                             :: ph_dir_fe
  character(len=256)                             :: dvscf_dir_fe
  character(len=256)                             :: dvscf_name_fe
  character(len=256)                             :: data_dir_fe
  character(len=256)                             :: data_name_fe
  character(len=256)                             :: qlist_fe
  character(len=256)                             :: fc_mat_fe
  integer                                        :: nat_fe
  integer                                        :: ntyp_fe
  real(kind=dp)                                  :: alat_fe
  real(kind=dp), dimension(3,3)                  :: at_fe
  real(kind=dp), dimension(3,3)                  :: bg_fe
  integer                                        :: nr1_fe
  integer                                        :: nr2_fe
  integer                                        :: nr3_fe
  real(kind=dp), allocatable, dimension(:,:)     :: tau_fe
  real(kind=dp), allocatable, dimension(:)       :: amass_fe
  integer, allocatable, dimension(:)             :: ityp_fe
  character(len=3), allocatable, dimension(:)    :: atom_labels_fe
  integer                                        :: nq1_fe
  integer                                        :: nq2_fe
  integer                                        :: nq3_fe
  integer                                        :: nqirr_fe
  integer                                        :: nqmesh_fe
  real(kind=dp), allocatable, dimension(:,:,:,:) :: frc_fe
  !
  ! loop variables
  integer :: iq ! loop in q points
  integer :: ia ! loop in atoms
  integer :: ja ! loop in atoms
  integer :: id ! loop in directions (x, y, z)
  integer :: jd ! loop in directions (x, y, z)
  integer :: imode ! loop in modes
  integer :: jmode ! loop in modes
  integer :: irr1 ! loop in unit cell
  integer :: irr2 ! loop in unit cell
  integer :: irr3 ! loop in unit cell
  !
  ! Unit cell <--> Supercell
  integer :: ia_original ! index of the corresponding atom on the original unit cell
  integer :: ja_original ! index of the corresponding atom on the original unit cell
  integer :: ia_sc ! index of the corresponding atom on the supercell
  integer :: ja_sc ! index of the corresponding atom on the supercell
  integer, dimension(3) :: R_m ! tau(ia) - tau(ia_original)
  integer, dimension(3) :: R_l ! tau(ja) - tau(ja_original)
  integer, dimension(3) :: R ! R = R_l - R_m
  integer :: m1 ! index of the corresponding unit cell
  integer :: m2 ! index of the corresponding unit cell
  integer :: m3 ! index of the corresponding unit cell
  !
  ! Dynamical matrix
  real(kind=dp), parameter :: pmass = 1822.88848426_dp
  real(kind=dp), parameter :: aumev = 27211.396_dp
  complex(kind=dp), allocatable, dimension(:,:,:) :: dyn_q ! dynamical matrix for each q point
  complex(kind=dp), allocatable, dimension(:,:) :: dyn ! dynamical matrix for GAMMA
  real(kind=dp), allocatable, dimension(:) :: w2 ! frequencies
  !
  ! Forces and derivatives
  real(kind=dp), allocatable, dimension(:,:) :: f_0, f_fe_p, f_fe_m
  real(kind=dp) :: dx
  integer :: atm, dir
  character(len=256) :: scf_fe
  !
  ! Phonon DOS variables
  integer :: ie, Ne=1000
  real(kind=dp) :: Emin=-10.0, Emax=100.0, de, sigma=0.008, E, DOS
  !
  ! Sum Rule variables
  real(kind=dp), allocatable, dimension(:,:,:) :: zeu
  real(kind=dp) :: tmp
  !
  ! Polarizations
  complex(kind=dp), allocatable, dimension(:,:) :: z
  !
  ! Fourier transform varibales
  real(kind=dp), allocatable, dimension(:,:) :: qmesh ! complete q mesh in crystal coordinates
  real(kind=dp), dimension(3) :: qpoint ! single q point in crystal coordinates
  real(kind=dp) :: qr ! q.R dot product
  real(kind=dp) :: w ! weight of the q point
  !
  ! Other variables
  character(len=20) :: dummy
  !
  !
  20 format(A)
  !
  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(*,20) '====================================================='
  write(*,20) '|                   program fc.x                    |'
  write(*,20) '|        ---------------------------------          |'
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  !       Define parameters
  !================================================================================
  !
  ! MgO-Ag(100) directory and files
  ! mesh_dir_mgo='/home/haritz/Kalkuluak/Ag/Fe-MgO-Ag100/1ML/ph/MgO-Ag100/'
  mesh_dir_mgo='/home/haritz/tmp/MgO-Ag100/v0/'
  prefix_mgo='mgo-ag100'
  ! ph_dir_mgo='ph/'
  ph_dir_mgo='../ph/'
  fc_mat_mgo='mgo-ag100_441.fc'
  data_dir_mgo='data_dir/'
  data_name_mgo='data-file.'
  ! dvscf_dir_mgo='dvscf_dir/'
  ! dvscf_name_mgo='mgo-ag100.dvscf_q'
  qlist_mgo='qirred441.dat'
  nq1_mgo=4
  nq2_mgo=4
  nq3_mgo=1
  nqirr_mgo=6
  nqmesh_mgo=nq1_mgo*nq2_mgo*nq3_mgo
  !
  ! Fe-MgO-Ag(100) directory and files
  ! mesh_dir_fe='/home/haritz/Kalkuluak/Ag/Fe-MgO-Ag100/1ML/ph/Fe-MgO-Ag100/smearing_0.01/'
  mesh_dir_fe='/home/haritz/tmp/Fe-MgO-Ag100/colinear/v0/'
  prefix_fe='Fe-MgO-Ag100'
  ! ph_dir_fe='ph/'
  ph_dir_fe='../ph/'
  fc_mat_fe='fe-mgo-ag100_441.fc'
  ! data_dir_fe='data_dir/'
  ! data_name_fe='data-file.'
  ! dvscf_dir_fe='dvscf_dir/'
  ! dvscf_name_fe='mgo-ag100.dvscf_q'
  ! qlist_fe='qirred441.dat'
  nq1_fe=1
  nq2_fe=1
  nq3_fe=1
  nqirr_fe=1
  nqmesh_fe=nq1_fe*nq2_fe*nq3_fe
  !
  !
  !================================================================================
  !    Read the QE MgO-Ag(100) surface calculation data
  !================================================================================
  !
  mesh_dir=mesh_dir_mgo
  prefix=prefix_mgo
  ph_dir=ph_dir_mgo
  fc_mat=fc_mat_mgo
  data_dir=data_dir_mgo
  data_name=data_name_mgo
  ! dvscf_dir=dvscf_dir_mgo
  ! dvscf_name=dvscf_name_mgo
  qlist=qlist_mgo
  nq1=nq1_mgo
  nq2=nq2_mgo
  nq3=nq3_mgo
  nqirr=nqirr_mgo
  !
  call read_parameters_data_file_xml()
  !
  ! read clean MgO-Ag(100) surface phonon related data
  call read_ph_information_xml()
  !
  ! save readed parameters of the MgO-Ag(100) surface
  nat_mgo=nat
  ntyp_mgo=ntyp
  alat_mgo=alat
  at_mgo=at
  bg_mgo=bg
  nr1_mgo=nr1
  nr2_mgo=nr2
  nr3_mgo=nr3
  allocate(tau_mgo(3,nat_mgo))
  allocate(amass_mgo(ntyp_mgo))
  allocate(ityp_mgo(nat_mgo))
  allocate(atom_labels_mgo(ntyp_mgo))
  tau_mgo = tau
  amass_mgo = amass
  ityp_mgo = ityp
  atom_labels_mgo = atom_labels
  deallocate(tau)
  deallocate(amass)
  deallocate(ityp)
  deallocate(atom_labels)
  deallocate(atom_pfile)
  deallocate(s)
  deallocate(ftau)
  deallocate(can_use_TR)
  !

  !! BEGIN_DEBUG
  !! Print some variables readed from QE files
  if (.false.) then
    print*,'DEBUG:','nat_mgo', nat_mgo
    print*,'DEBUG:','ntyp_mgo', ntyp_mgo
    print*,'DEBUG:','alat_mgo', alat_mgo
    do id=1,3
      write(*,'(2a,x,3f16.8)') 'DEBUG:','at_fe', at_mgo(:,id)
    enddo
    print*,'DEBUG:','bg_mgo', bg_mgo
    print*,'DEBUG:','nr1_mgo, nr2_mgo, nr3_mgo', nr1_mgo, nr2_mgo, nr3_mgo
    do ia=1,nat_mgo
      print'(a,a,x,a3,i3,3f15.7)','DEBUG:','label ityp tau_mgo', atom_labels_mgo(ityp_mgo(ia)), ityp_mgo(ia), tau_mgo(:,ia)*alat
    enddo
    do ia=1,nat_mgo
      print'(a,a,x,i3,x,a,i3,f10.5)','DEBUG:','nat label ityp amass_mgo', ia, atom_labels_mgo(ityp_mgo(ia)), ityp_mgo(ia), amass_mgo(ityp_mgo(ia))
    enddo
    print*,'DEBUG:','nq1_mgo, nq2_mgo, nq3_mgo', nq1_mgo, nq2_mgo, nq3_mgo
    print*,'DEBUG:','nqirr_mgo', nqirr_mgo
  endif
  !! END_DEBUG
  !
  !
  !================================================================================
  !   Read the force constant matrix of the clean MgO-Ag(100) surface
  !================================================================================
  !
  fc_file_name_mgo=trim(trim(mesh_dir_mgo)//"/"//trim(ph_dir_mgo)//"/"//trim(fc_mat_mgo))
  !
  allocate(amass(ntyp))
  call readfc(fc_file_name_mgo,nq1,nq2,nq3,nat,alat,at,ntyp,amass)
  deallocate(amass)
  !
  ! save the force constants of MgO-Ag(100) surface
  allocate( frc_mgo(nq1,nq2,nq3,3,3,nat,nat) )
  frc_mgo = frc
  deallocate(frc)
  !
  !
  !================================================================================
  !    Read the QE Fe-MgO-Ag(100) surface calculation data from V0
  !================================================================================
  !
  mesh_dir=mesh_dir_fe
  prefix=prefix_fe
  ph_dir=ph_dir_fe
  fc_mat=fc_mat_fe
  data_dir=data_dir_fe
  data_name=data_name_fe
  ! dvscf_dir=dvscf_dir_fe
  ! dvscf_name=dvscf_name_fe
  qlist=qlist_fe
  !
  call read_parameters_data_file_xml()
  !
  ! save readed parameters of the Fe-MgO-Ag(100) surface
  !
  nat_fe=nat
  ntyp_fe=ntyp
  alat_fe=alat
  at_fe=at
  bg_fe=bg
  nr1_fe=nr1
  nr2_fe=nr2
  nr3_fe=nr3
  allocate(tau_fe(3,nat_fe))
  allocate(amass_fe(ntyp_fe))
  allocate(ityp_fe(nat_fe))
  allocate(atom_labels_fe(ntyp_fe))
  tau_fe=tau
  amass_fe=amass
  ityp_fe = ityp
  atom_labels_fe = atom_labels
  deallocate(tau)
  deallocate(amass)
  deallocate(ityp)
  deallocate(atom_labels)
  deallocate(atom_pfile)
  !

  !! BEGIN_DEBUG
  !! Print some variables readed from QE files
  if (.false.) then
    print*,'DEBUG:','nat_fe', nat_fe
    print*,'DEBUG:','ntyp_fe', ntyp_fe
    print*,'DEBUG:','alat_fe', alat_fe
    do id=1,3
      write(*,'(2a,x,3f16.8)') 'DEBUG:','at_fe', at_fe(:,id)
    enddo
    print*,'DEBUG:','bg_fe', bg_fe
    print*,'DEBUG:','nr1_fe, nr2_fe, nr3_fe', nr1_fe, nr2_fe, nr3_fe
    do ia=1,nat_fe
      print'(a,a,x,a3,i3,3f15.7)','DEBUG:','label ityp tau_fe', atom_labels_fe(ityp_fe(ia)), ityp_fe(ia), tau_fe(:,ia)*alat
    enddo
    do ia=1,nat_fe
      print'(a,a,x,i3,x,a,i3,f10.5)','DEBUG:','nat label ityp amass_fe', ia, atom_labels_fe(ityp_fe(ia)), ityp_fe(ia), amass_fe(ityp_fe(ia))
    enddo
    print*,'DEBUG:','nq1_fe, nq2_fe, nq3_fe', nq1_fe, nq2_fe, nq3_fe
    print*,'DEBUG:','nqirr_fe', nqirr_fe
  endif
  !! END_DEBUG

  !
  !
  !================================================================================
  !    Calculate the dynamical matrix of the clean surface from the FC matrix
  !================================================================================
  !
  ! Create the q mesh
  allocate(qmesh(3,nqmesh_mgo))
  call generate_kmesh(qmesh,nq1_mgo,nq2_mgo,nq3_mgo)
  !
  ! Calculate the dynamical matrix for each point of the q mesh
  allocate(dyn_q(3*nat_mgo,3*nat_mgo,nqmesh_mgo))
  allocate(tau(3,nat_mgo),amass(ntyp_mgo),ityp(nat_mgo))
  tau = tau_mgo
  amass=amass_mgo
  ityp=ityp_mgo
  bg=bg_mgo
  nat=nat_mgo
  do iq=1,nqmesh_mgo
    print*, qmesh(:,iq)
    call mat_inv_four_t(qmesh(:,iq),nq1_mgo,nq2_mgo,nq3_mgo,3*nat_mgo,frc_mgo,dyn_q(:,:,iq))
  enddo
  deallocate(tau, amass, ityp)
  !
  !
  !================================================================================
  !    Calculate the FC matrix of the clean surface from the dynamical matrix
  !================================================================================
  !
  ! calculate the force constant
  frc_mgo=cmplx_0
  ! for each cell of the supercell
  do irr1=1,nq1_mgo
    do irr2=1,nq2_mgo
      do irr3=1,nq3_mgo
        !
        print*, 'R=', irr1, irr2, irr3
        !
        ! Sum in q
        do iq=1,nqmesh_mgo
          qpoint = qmesh(:,iq)
          w=1.0_dp/(nqmesh_mgo)
          write(*,'(a,i3,3f12.6,2x,f12.6)') 'q=', iq, qpoint(:)
          !
          qr = tpi*(qpoint(1)*(irr1-1) + qpoint(2)*(irr2-1) + qpoint(3)*(irr3-1))
          !
          do imode=1,3*nat_mgo
            !
            ia = (imode-1)/3 + 1
            id = mod(imode-1,3) + 1
            ia_sc = ia + nat_mgo*(nq2_mgo*(irr1-1)+(irr2-1)+nq1_mgo*nq2_mgo*(irr3-1))
            !
            do jmode=1,3*nat_mgo
              !
              ja = (jmode-1)/3 + 1
              jd = mod(jmode-1,3) + 1
              ja_sc = ja + nat_mgo*(nq2_mgo*(irr1-1)+(irr2-1)+nq1_mgo*nq2_mgo*(irr3-1))
              !
              frc_mgo(irr1,irr2,irr3,id,jd,ia,ja) = &
                  frc_mgo(irr1,irr2,irr3,id,jd,ia,ja) + &
                  w*dyn_q(imode,jmode,iq) * exp(-cmplx_i*qr) * &
                  sqrt(amass_mgo(ityp_mgo(ia))*amass_mgo(ityp_mgo(ja))) * &
                  (pmass/2.d0) / & !Up to this in Ry.
                  (aumev/2.d0)**2 !Now in meV.
              !
            enddo
          enddo
          !
        enddo ! iq
        !
      enddo
    enddo
  enddo
  !
  deallocate(qmesh, dyn_q)
  !
  !
  !================================================================================
  !    Create the FC matrix of the clean Fe/MgO/Ag(100) supercell surface
  !================================================================================
  !
  nat=nat_fe
  allocate(frc_fe(3,3,nat,nat))
  frc_fe=zero
  !
  ! put tau_fe in units of alat_mgo
  allocate(tau(3,nat))
  tau = tau_fe(:,1:nat) * alat_fe/alat_mgo
  !
  do ia=1,nat_mgo*nqmesh_mgo
    !
    ! equivalent atom of the original cell
    ia_original = mod(ia-1,nat_mgo) + 1
    R_m(:) = nint(tau(:,ia) - tau(:,ia_original))
    !
    do ja=1,nat_mgo*nqmesh_mgo
      !
      ! equivalent atom of the original cell
      ja_original = mod(ja-1,nat_mgo) + 1
      R_l(:) = nint(tau(:,ja) - tau(:,ja_original))
      !
      ! R = R_l - R_m
      R = R_l - R_m
      if ( sqrt((((tau_fe(3,ia_original)-tau_fe(3,ja_original))*alat_fe)**2)) > 30.0_dp ) cycle
      !
      m1 = mod(R(1)+1,nq1_mgo)
      if(m1.le.0) m1=m1+nq1_mgo
      m2 = mod(R(2)+1,nq2_mgo)
      if(m2.le.0) m2=m2+nq2_mgo
      m3 = mod(R(3)+1,nq3_mgo)
      if(m3.le.0) m3=m3+nq3_mgo
      !
      ! complete the FC matrix of the supercell
      frc_fe(:,:,ia,ja) = frc_mgo(m1,m2,m3,:,:,ia_original,ja_original)
      !
    end do
  end do
  !
  deallocate(frc_mgo, tau)
  !
  !
  !
  !================================================================================
  !    Simple Acoustic Sum Rule on force constants in real space
  !================================================================================
  ! !
  ! allocate(zeu(3,3,nat))
  ! zeu =0.d0
  ! !
  ! do id=1,3
  !   do jd=1,3
  !     tmp=0.0d0
  !     do ia=1,nat
  !       tmp = tmp + zeu(id,jd,ia)
  !     end do
  !     do ia=1,nat
  !       zeu(id,jd,ia) = zeu(id,jd,ia) - tmp/nat
  !     end do
  !   end do
  ! end do
  ! !
  ! do id=1,3
  !   do jd=1,3
  !     do ia=1,nat
  !       tmp=0.0d0
  !       do ja=1,nat
  !         tmp=tmp+frc_fe(id,jd,ia,ja)
  !       end do
  !       frc_fe(id,jd,ia,ia) = frc_fe(id,jd,ia,ia) - tmp
  !     end do
  !   end do
  ! end do
  ! !
  ! deallocate(zeu)
  ! !
  ! !
  !================================================================================
  !   Save DOS and polarizations
  !================================================================================
  !
  allocate( dyn(3*nat,3*nat), w2(3*nat) )
  dyn=cmplx_0
  w2=zero
  !
  do ia=1,nat
    do ja=1,nat
      do id=1,3
        do jd=1,3
          dyn((ia-1)*3+id, (ja-1)*3+jd)= dyn((ia-1)*3+id, (ja-1)*3+jd) + &
                          frc_fe(id, jd, ia, ja) / &
                          sqrt(amass_fe(ityp_fe(ia))*amass_fe(ityp_fe(ja))) / &
                          (pmass/2.d0) * & !Up to this in Ry.
                          (aumev/2.d0)**2 !Now in meV.
          !
        end do
      end do
    end do
  end do
  !
  dyn = 0.5d0 * (dyn + conjg(transpose(dyn)))
  !
  ! diagonalize the dynamical matrix to obtain eigenmodes
  call diagonalize_cmat( 3*nat, dyn, w2 )
  !
  do imode=1,3*nat
    write(dummy,'(1f20.10)') w2(imode)
    read(dummy,*) w2(imode)
    if (w2(imode)<zero) then
      w2(imode)=-sqrt(abs(w2(imode)))
    else
      w2(imode)=sqrt(w2(imode))
    endif
  enddo
  !
  !
  ! calculate the phonon DOS
  ! w2(1)=0.0_dp
  ! w2(2)=0.0_dp
  ! w2(3)=0.0_dp
  !
  Emin = 0.0
  Emax = 80.0
  de=(Emax-Emin)/(Ne-1)
  sigma=0.7
  open(unit=990, file='dos0.dat')
  !
  do ie=1,Ne
    E = Emin + de*(ie-1)
    DOS = 0.0_dp
    do imode=1,3*nat
      ! if (w2(imode)<0.4_dp) then
      !   w2(imode)=-10.0
      ! endif
      DOS = DOS &
            !+ 1.0_dp/(sigma*sqrt(2.0_dp*3.1416_dp))*exp(-0.5_dp*((+w2(imode)-E)/sigma)**2.0) &
            !- 1.0_dp/(sigma*sqrt(2.0_dp*3.1416_dp))*exp(-0.5_dp*((-w2(imode)-E)/sigma)**2.0)
            + 1.0_dp/( sigma*3.1416_dp* (1+( ( w2(imode)-E ) / sigma )**2 ) ) &
            - 1.0_dp/( sigma*3.1416_dp* (1+( (-w2(imode)-E ) / sigma )**2 ) )
    enddo
    write(990,*) E, DOS/(3*nat)
  enddo
  !
  close(990)
  !
  ! print polarizations
  allocate(z(3,nat))
  open(unit=112, file='vib0.xyz')
  do imode=1,3*nat
    write(112,*) nat
    write(112,'(a6,i4,a5,1f10.6,a5)') 'Mode #',imode,' f = ', w2(imode),' meV.'
    do ia = 1,nat
      do id=1,3
        z(id,ia) = dyn(3*(ia-1)+id,imode)
      enddo
      write(112,'(a4,x,2(3f15.10))') atom_labels_fe(ityp_fe(ia)), ( tau_fe(id,ia)*alat_fe*bohr ,id=1,3 ), (real(z(id,ia)), id=1,3)
    enddo
  end do
  close(112)
  deallocate(dyn,w2,z)
  !
  !
  !================================================================================
  !   Read equilibrium forces
  !================================================================================
  !
  allocate(f_0(3,nat_fe),f_fe_p(3,nat_fe),f_fe_m(3,nat_fe))
  f_0=zero
  f_fe_p=zero
  f_fe_m=zero
  !
  scf_fe=trim(mesh_dir_fe)//'fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_0)
  !


  !! BEGIN_DEBUG
  !! Print the forces of the relaxed structure
  if (.false.) then
    do ia=1,nat_fe
      write(*,'(i,x,3f14.8)') ia, f_0(:,ia)
    enddo
  endif
  !! END_DEBUG

  !! BEGIN_DEBUG
  !! Check the if the derivative is the same using the + - displacements or
  !! using only the + displacement
  if (.false.) then
    !
    atm = 273 ! number of the moved atom
    dir = 1 ! direction in which the atom has been moved
    dx = 0.0001*alat ! displacement of the atom
    !
    scf_fe=trim(mesh_dir_fe)//'Fe_1_x/mgo-ag100.scf.out'
    call read_f(scf_fe,nat_fe,f_fe_p)
    scf_fe=trim(mesh_dir_fe)//'Fe_1_-x/mgo-ag100.scf.out'
    call read_f(scf_fe,nat_fe,f_fe_m)
    !
    do ia=1,nat_fe
      write(1001,*) atm,ia, -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
      write(1001,*) atm,ia, -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
      write(1001,*) atm,ia, -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
      !
      write(1001,*) ia,atm, -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
      write(1001,*) ia,atm, -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
      write(1001,*) ia,atm, -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
    enddo
    !
    atm = 273 ! number of the moved atom
    dir = 1 ! direction in which the atom has been moved
    dx = 0.0001*alat ! displacement of the atom
    !
    scf_fe=trim(mesh_dir_fe)//'Fe_1_x/mgo-ag100.scf.out'
    call read_f(scf_fe,nat_fe,f_fe_p)
    scf_fe=trim(mesh_dir_fe)//'V0/mgo-ag100.scf.out'
    call read_f(scf_fe,nat_fe,f_fe_m)
    !
    do ia=1,nat_fe
      write(1002,*) atm,ia, -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( dx )
      write(1002,*) atm,ia, -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( dx )
      write(1002,*) atm,ia, -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( dx )
      !
      write(1002,*) ia,atm, -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( dx )
      write(1002,*) ia,atm, -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( dx )
      write(1002,*) ia,atm, -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( dx )
    enddo
    !
  endif
  !! END_DEBUG


  !
  !
  !================================================================================
  !   Add iron atoms to the FC matrix
  !================================================================================
  !
  ! read scf file of Fe displaced on x direction
  atm = 273 ! number of the moved atom
  dir = 1 ! direction in which the atom has been moved
  dx = 0.001*alat ! displacement of the atom
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'273x+/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_p)
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'273x-/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_m)
  !
  ! add this to the FC matrix
  do ia=1,nat
    ! if ( sqrt(((tau_fe(3,atm)-tau_fe(3,ia))**2.0))*alat_fe > 30.0_dp ) cycle
    frc_fe(dir,1,atm,ia) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,2,atm,ia) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,3,atm,ia) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
    !
    frc_fe(1,dir,ia,atm) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(2,dir,ia,atm) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(3,dir,ia,atm) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
  enddo
  !
  ! read scf file of Fe_1 displaced on y direction
  atm = 273 ! number of the moved atom
  dir = 2 ! direction in which the atom has been moved
  dx = 0.001*alat ! displacement of the atom
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'273y+/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_p)
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'273y-/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_m)
  !
  ! add this to the FC matrix
  do ia=1,nat
    ! if ( sqrt(((tau_fe(3,atm)-tau_fe(3,ia))**2.0))*alat_fe > 30.0_dp ) cycle
    frc_fe(dir,1,atm,ia) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,2,atm,ia) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,3,atm,ia) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
    !
    frc_fe(1,dir,ia,atm) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(2,dir,ia,atm) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(3,dir,ia,atm) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
  enddo
  !
  ! read scf file of Fe_1 displaced on z direction
  atm = 273 ! number of the moved atom
  dir = 3 ! direction in which the atom has been moved
  dx = 0.001*alat ! displacement of the atom
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'273z+/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_p)
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'273z-/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_m)
  !
  ! add this to the FC matrix
  do ia=1,nat
    ! if ( sqrt(((tau_fe(3,atm)-tau_fe(3,ia))**2.0))*alat_fe > 30.0_dp ) cycle
    frc_fe(dir,1,atm,ia) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,2,atm,ia) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,3,atm,ia) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
    !
    frc_fe(1,dir,ia,atm) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(2,dir,ia,atm) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(3,dir,ia,atm) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
  enddo
  !
  !
  print*, "Fe(273)-O(186)"
  print*, frc_fe(1,1,273,186), frc_fe(1,2,273,186), frc_fe(1,3,273,186)
  print*, "Fe(273)-Fe(273)"
  print*, frc_fe(1,1,273,273), frc_fe(1,2,273,273), frc_fe(1,3,273,273)
  !
  !
  !
  !================================================================================
  !   Replace oxigen atoms on the FC matrix
  !================================================================================
  !
  !
  ! read scf file of O displaced on x direction
  atm = 186 ! number of the moved atom
  dir = 1 ! direction in which the atom has been moved
  dx = 0.001*alat ! displacement of the atom
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'186x+/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_p)
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'186x-/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_m)
  !
  ! add this to the FC matrix
  do ia=1,nat_fe
    if ( sqrt(((tau_fe(3,atm)-tau_fe(3,ia))**2.0))*alat_fe > 30.0_dp ) cycle
    frc_fe(dir,1,atm,ia) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,2,atm,ia) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,3,atm,ia) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
    !
    frc_fe(1,dir,ia,atm) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(2,dir,ia,atm) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(3,dir,ia,atm) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
  enddo
  !
  ! read scf file of O_2 displaced on y direction
  atm = 186 ! number of the moved atom
  dir = 2 ! direction in which the atom has been moved
  dx = 0.001*alat ! displacement of the atom
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'186y+/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_p)
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'186y-/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_m)
  !
  ! add this to the FC matrix
  do ia=1,nat_fe
    if ( sqrt(((tau_fe(3,atm)-tau_fe(3,ia))**2.0))*alat_fe > 30.0_dp ) cycle
    frc_fe(dir,1,atm,ia) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,2,atm,ia) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,3,atm,ia) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
    !
    frc_fe(1,dir,ia,atm) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(2,dir,ia,atm) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(3,dir,ia,atm) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
  enddo
  !
  ! read scf file of O_2 displaced on z direction
  atm = 186 ! number of the moved atom
  dir = 3 ! direction in which the atom has been moved
  dx = 0.001*alat ! displacement of the atom
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'186z+/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_p)
  scf_fe=trim(mesh_dir_fe)//trim(ph_dir_fe)//'186z-/fe-mgo-ag100.scf.out'
  call read_f(scf_fe,nat_fe,f_fe_m)
  !
  ! add this to the FC matrix
  do ia=1,nat_fe
    if ( sqrt(((tau_fe(3,atm)-tau_fe(3,ia))**2.0))*alat_fe > 30.0_dp ) cycle
    frc_fe(dir,1,atm,ia) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,2,atm,ia) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(dir,3,atm,ia) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
    !
    frc_fe(1,dir,ia,atm) = -( f_fe_p(1,ia) - f_fe_m(1,ia) ) / ( 2.0_dp*dx )
    frc_fe(2,dir,ia,atm) = -( f_fe_p(2,ia) - f_fe_m(2,ia) ) / ( 2.0_dp*dx )
    frc_fe(3,dir,ia,atm) = -( f_fe_p(3,ia) - f_fe_m(3,ia) ) / ( 2.0_dp*dx )
  enddo
  !
  !
  deallocate(f_0,f_fe_p,f_fe_m)
  !
  !
  print*, "Fe(186)-O(186)"
  print*, frc_fe(1,1,186,186), frc_fe(1,2,186,186), frc_fe(1,3,186,186)
  print*, "Fe(186)-Fe(273)"
  print*, frc_fe(1,1,186,273), frc_fe(1,2,186,273), frc_fe(1,3,186,273)
  !
  !
  !================================================================================
  !   Apply symmetries
  !================================================================================
  !
  ! allocate(tau_cryst(3,nat), rtau_cryst(3,nsym,nat), rtau(3,nsym,nat), rtau_index(nat,nsym))
  ! !
  ! call rot_atoms(nat,nsym, tau_fe)
  ! !
  ! ! Apply symmetries to the FC matrix
  ! call apply_all_sym()
  !
  !
  !================================================================================
  !    Simple Acoustic Sum Rule on force constants in real space
  !================================================================================
  !
  allocate(zeu(3,3,nat))
  zeu =0.d0
  !
  do id=1,3
    do jd=1,3
      tmp=0.0d0
      do ia=1,nat
        tmp = tmp + zeu(id,jd,ia)
      end do
      do ia=1,nat
        zeu(id,jd,ia) = zeu(id,jd,ia) - tmp/nat
      end do
    end do
  end do
  !
  do id=1,3
    do jd=1,3
      do ia=1,nat
        tmp=0.0d0
        do ja=1,nat
          tmp=tmp+frc_fe(id,jd,ia,ja)
        end do
        frc_fe(id,jd,ia,ia) = frc_fe(id,jd,ia,ia) - tmp
      end do
    end do
  end do
  !
  deallocate(zeu)
  !
  !
  !================================================================================
  !    Calculate the dynamical matrix at Gamma
  !================================================================================
  !
  allocate( dyn(3*nat,3*nat), w2(3*nat) )
  dyn=cmplx_0
  w2=zero
  !
  do ia=1,nat
    do ja=1,nat
      do id=1,3
        do jd=1,3
          dyn((ia-1)*3+id, (ja-1)*3+jd)= dyn((ia-1)*3+id, (ja-1)*3+jd) + &
                          frc_fe(id, jd, ia, ja) / &
                          sqrt(amass_fe(ityp_fe(ia))*amass_fe(ityp_fe(ja))) / &
                          (pmass/2.d0) * & !Up to this in Ry.
                          (aumev/2.d0)**2 !Now in meV.
          !
        end do
      end do
    end do
  end do
  !
  dyn = 0.5d0 * (dyn + conjg(transpose(dyn)))
  !
  ! diagonalize the dynamical matrix to obtain eigenmodes
  call diagonalize_cmat( 3*nat, dyn, w2 )
  !
  ! round the eigenvalues to 10th decimal (to avoid small but negative eigenvalues)
  do imode=1,3*nat
    write(dummy,'(1f20.10)') w2(imode)
    read(dummy,*) w2(imode)
    if (w2(imode)<zero) then
      w2(imode)=-sqrt(abs(w2(imode)))
    else
      w2(imode)=sqrt(w2(imode))
    endif
  enddo
  !
  ! calculate the phonon DOS
  ! w2(1)=0.0_dp
  ! w2(2)=0.0_dp
  ! w2(3)=0.0_dp
  !
  Emin = 0.0
  Emax = 80.0
  de=(Emax-Emin)/(Ne-1)
  sigma=0.7
  open(unit=991,file='dos.dat')
  !
  do ie=1,Ne
    E = Emin + de*(ie-1)
    DOS = 0.0_dp
    do imode=1,3*nat
      ! if (w2(imode)<0.4_dp) then
      !   w2(imode)=-10.0
      ! endif
      DOS = DOS &
            !+ 1.0_dp/(sigma*sqrt(2.0_dp*3.1416_dp))*exp(-0.5_dp*((+w2(imode)-E)/sigma)**2.0) &
            !- 1.0_dp/(sigma*sqrt(2.0_dp*3.1416_dp))*exp(-0.5_dp*((-w2(imode)-E)/sigma)**2.0)
            + 1.0_dp/( sigma*3.1416_dp* (1+( ( w2(imode)-E ) / sigma )**2 ) ) &
            - 1.0_dp/( sigma*3.1416_dp* (1+( (-w2(imode)-E ) / sigma )**2 ) )
    enddo
    write(991,*) E, DOS/(3*nat)
  enddo
  close(991)
  !
  ! print polarizations
  allocate(z(3,nat))
  open(unit=112, file='vib.xyz')
  do imode=1,3*nat
    write(112,*) nat
    write(112,'(a6,i4,a5,1f10.6,a5)') 'Mode #',imode,' f = ', w2(imode),' meV.'
    do ia = 1,nat
      do id=1,3
        z(id,ia) = dyn(3*(ia-1)+id,imode)
      enddo
      write(112,'(a4,x,2(3f15.10))') atom_labels_fe(ityp_fe(ia)), ( tau_fe(id,ia)*alat_fe*bohr ,id=1,3 ), (real(z(id,ia)), id=1,3)
    enddo
  end do
  close(112)
  !
  9000 format(a,3 ('(',f15.6,',',f15.6,')' ) )
  9100 format (2a,1x,3 ('(',f10.6,',',f10.6,1x,')',5x))
  !
  deallocate(z)
  !

  !
  deallocate(dyn)




contains

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


  subroutine read_f(scf_file,nat,f)
    !
    ! subroutine to read the forces from a scf file
    !
    use, intrinsic :: iso_fortran_env, only: iostat_end   ! end of file
    !
    use intw_utility, only: find_free_unit
    !
    implicit none
    !
    ! input variables
    character(len=256), intent(in) :: scf_file
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
    open(unit=iounit, file=trim(scf_file), iostat=ios, status="old", action="read")
    if ( ios /= 0 ) then
      write(*,'(a)') "read_f: Error opening scf_file "//trim(scf_file)
      stop
    else
      write(*,'(a)') "read_f: Reading forces from scf_file "//trim(scf_file)
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
          if (dummy_i1/=i) write(*,'(a)') 'read_f: WARNING, atoms are not properly ordered on scf file'
        enddo
        !
        close(unit=iounit)
        return
        !
      endif
      !
      if ( stat == iostat_end ) then
        write(*,*) 'Error reading file '//trim(scf_file)//'. End of file before reading forces.'
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
  !
  !
  subroutine apply_sym(isym)
    !
    implicit none
    !
    integer, intent(in) :: isym
    integer :: Sia, Sja
    real(kind=dp), dimension(3) :: R, R_sym
    real(kind=dp), dimension(3,3) :: tmp, SFC
    !
    !
    do ia=1,nat_fe
      do ja=1,nat_fe
        !
        Sia = rtau_index(ia,isym)
        Sja = rtau_index(ja,isym)
        R(:) = tau_fe(:,ia) - tau_fe(:,ja)
        R_sym(:) = tau_fe(:,Sia) - tau_fe(:,Sja)
        !
        print*, ''
        print*, 'R', R
        print*, 'SR', matmul(s(:,:,isym),R(:))
        print*, 'R_sym', R_sym
        !
        if ( sqrt( sum( (matmul(s(:,:,isym),R(:)) - R_sym)**2 ) ) < 0.00001_dp ) then
          !
          if ( maxval( matmul( ainv(real(s(:,:,isym),dp)), matmul( frc_fe(:,:,Sia,Sja),s(:,:,isym)) ) - frc_fe(:,:,ia,ja) ) > 0.0001_dp ) then
            !
            print*, ia, ja, Sia, Sja
            print*, 'FC(:,:,ia,ja):'
            print'(3f20.15)', frc_fe(:,:,ia,ja)
            ! print'(3f20.15)', frc_fe(1,:,ia,ja)
            ! print'(3f20.15)', frc_fe(2,:,ia,ja)
            ! print'(3f20.15)', frc_fe(3,:,ia,ja)
            print*, 'S^-1 FC(:,:,Sia,Sja) S:'
            SFC(:,:) = matmul( ainv(real(s(:,:,isym),dp)), matmul( frc_fe(:,:,Sia,Sja),s(:,:,isym)) )
            print'(3f20.15)', SFC(:,:)
            read(*,*)
            !
            tmp(:,:) = ( frc_fe(:,:,ia,ja) + SFC(:,:) ) / 2.0_dp
            frc_fe(:,:,ia,ja) = tmp(:,:)
            frc_fe(:,:,Sia,Sja) = matmul( s(:,:,isym), matmul( frc_fe(:,:,Sia,Sja),ainv(real(s(:,:,isym),dp) ) ) )
          endif
          !
        else
          !
          ! print*, 'Do not apply'
          !
        endif
        !
      enddo
    enddo
    !
  end subroutine apply_sym
  !
  !
  subroutine apply_all_sym()
    !
    implicit none
    !
    integer :: isym
    integer :: Sia, Sja, deg
    real(kind=dp), dimension(3) :: R, R_sym
    real(kind=dp), dimension(3,3) :: tmp, SFC
    !
    !
    do ia=1,nat_fe
      do ja=1,nat_fe
        !
        deg=0
        tmp=0.0_dp
        do isym=1,nsym
          !
          Sia = rtau_index(ia,isym)
          Sja = rtau_index(ja,isym)
          R(:) = tau_fe(:,ia) - tau_fe(:,ja)
          R_sym(:) = tau_fe(:,Sia) - tau_fe(:,Sja)
          ! print*, ia, ja, Sia, Sja
          ! print'(3f10.6)', R
          ! print'(3f10.6)', R_sym
          !
          if ( sqrt( sum( (matmul(s(:,:,isym),R(:)) - R_sym)**2 ) ) < 0.00001_dp ) then
            !
            deg=deg+1
            ! print*, 'Apply', deg
            !
            SFC(:,:) = matmul( ainv(real(s(:,:,isym),dp)), matmul( frc_fe(:,:,Sia,Sja),s(:,:,isym)) )
            !
            tmp(:,:) = tmp(:,:) + SFC(:,:)
            !
          else
            !
            ! print*, 'Do not apply', deg
            !
          endif
          !
          ! read(*,*)
        enddo ! isym
        !
        do isym=1,nsym
          !
          Sia = rtau_index(ia,isym)
          Sja = rtau_index(ja,isym)
          R(:) = tau_fe(:,ia) - tau_fe(:,ja)
          R_sym(:) = tau_fe(:,Sia) - tau_fe(:,Sja)
          !
          ! print*, ia, ja, Sia, Sja
          ! print'(3f10.6)', R
          ! print'(3f10.6)', R_sym
          if ( sqrt( sum( (matmul(s(:,:,isym),R(:)) - R_sym)**2 ) ) < 0.00001_dp ) then
            ! print*, 'Apply', deg
            ! print*, 'FC(:,:,Sia,Sja)'
            ! print'(3f18.10)', frc_fe(:,:,Sia,Sja)
            ! print*, 'tmp(:,:)'
            ! print'(3f18.10)', matmul( (real(s(:,:,isym),dp)), matmul( tmp(:,:)/deg,ainv(real(s(:,:,isym),dp))) )
            frc_fe(:,:,Sia,Sja) = matmul( (real(s(:,:,isym),dp)), matmul( tmp(:,:)/deg,ainv(real(s(:,:,isym),dp))) )
            ! frc_fe(:,:,Sia,Sja) = tmp(:,:)/deg
            ! read(*,*)
          endif
          !
        enddo
        !
      enddo
    enddo
    !
  end subroutine apply_all_sym



  ! subroutine read_input_file()
  !   !
  !   ! This subroutine sets default values for variables of input namelist
  !   ! and then reads them from input file if there are specified
  !   !
  !   implicit none
  !   !
  !   character(len=256) :: input_file_name
  !   ! Namelists
  !   !
  !   ! input variables namelist
  !   NAMELIST / input_fc_fe / mesh_dir, prefix, ph_dir, dvscf_dir, dvscf_name, &
  !                            ep_mat_file, nat_move, ia_move_to_ia, calc_epmat, &
  !                            nbands_initial, nbands_final
  !   !
  !   ! Set default values
  !   mesh_dir = '/home/haritz/Kalkuluak/Probak/Fe-O/K/'
  !   prefix = 'fe-o'
  !   nk1 = 1
  !   nk2 = 1
  !   nk3 = 1
  !   calc_epmat=.true.
  !   ph_dir='./'
  !   dvscf_dir='./dvscf_dir/'
  !   dvscf_name='fe-o.dvscf_q1'
  !   ep_mat_file=trim("ep_mat_is_1_iband_1274_VL.dat")
  !   nq1=1
  !   nq2=1
  !   nq3=1
  !   nqirr=1
  !   nat_move = 0
  !   ia_move_to_ia = 0
  !   nbands_initial = 1
  !   nbands_final = 15
  !   !
  !   INQUIRE(stdin, NAME=input_file_name)
  !   if (input_file_name(1:4)=="/dev") then
  !     ! there is no input file: use default values
  !     return
  !   else
  !     read(stdin,input_ep_melements)
  !   endif
  !   !
  !   ! Reopen terminal as input file
  !   ! close(unit=stdin)
  !   ! open(unit=stdin,file='/dev/tty')
  !   !
  !   !
  !   num_bands = nbands_final - nbands_initial + 1
  !   !
  !   !
  ! end subroutine read_input_file



end program fc
