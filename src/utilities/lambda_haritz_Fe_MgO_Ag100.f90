program lambda_prog
  !
  !
  !
  use kinds, only: dp
  !
  use intw_input_parameters, only: mesh_dir, prefix, ph_dir, dvscf_dir, &
                                   dvscf_name, data_dir, data_name, &
                                   qlist, fc_mat, nq1, nq2, nq3, nqirr
  !
  use intw_reading, only: alat, nat, tau, nbands, amass, ityp, atom_labels, &
                          read_parameters_data_file_xml
  !
  use intw_reading, only: nat
  !
  use intw_useful_constants
  !
  use intw_utility, only:  find_free_unit
  !
  use intw_haritz
  !
  implicit none
  !
  character(len=256) :: polarizations_file
  character(len=256) :: ep_mat_el_file
  character(len=256) :: ep_mat_dir
  !
  ! Eigenmodes energy and polarization
  real(kind=dp), allocatable, dimension(:) :: w
  real(kind=dp), allocatable, dimension(:,:) :: polarization
  !
  ! Matrix element
  complex(kind=dp), allocatable :: ep_mat_el_ia(:,:)
  complex(kind=dp), allocatable :: ep_mat_el_nu(:,:)
  complex(kind=dp), allocatable :: ep_mat_el_nl(:,:)
  real(kind=dp), allocatable :: ep_mat_el_vl(:,:)
  !
  ! Loop variables
  integer :: iband
  integer :: imode
  integer :: jmode
  !
  !
  integer :: ia
  integer :: alpha
  real(kind=dp), dimension(3) :: tau_re
  !
  !
  real(kind=dp), parameter :: pmass = 1822.88848426_dp
  real(kind=dp), parameter :: aumev = 27211.396_dp
  !
  !
  real(kind=dp) :: lambda, g, sigma
  real(kind=dp), allocatable, dimension(:) :: e
  !
  ! Input parameters
  integer :: initial = 1274 ! or 1280

  integer :: nl_unit, nl_rl, ios, im, ib, tot, haundiak
  logical :: nonl
  !
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
  ! Fe-MgO-Ag(100) directory and files
  ! general parameters of intw
  mesh_dir = "/home/haritz/Kalkuluak/Ag/Fe-MgO-Ag100/1ML/ph/Fe-MgO-Ag100/spin/smearing_0.01/"
  mesh_dir = "/home/haritz/Kalkuluak/Probak/Fe-MgO-Ag100/bat/"
  prefix = "Fe-MgO-Ag100"
  ph_dir = "ph/"
  data_dir = "data_dir/"
  data_name = "data-file."
  ! my own parameters
  ep_mat_dir = "matdir/"
  polarizations_file = "vib.xyz"
  ep_mat_el_file = "ep_mat_is_1_iband_1274_VL.dat"
  !
  !
  !================================================================================
  !    Read the QE Fe-MgO-Ag(100) surface calculation data from V0
  !================================================================================
  !
  call read_parameters_data_file_xml()
  !
  !
  !================================================================================
  !    Rearrange the atoms so that the last ones are the irons
  !================================================================================
  !
  !
  ! tau_re = tau(:,189) ! save the iron coordinates
  ! do ia=190,nat
  !   tau(:,ia-1) = tau(:,ia)
  !   ityp(ia-1) = ityp(ia)
  ! enddo
  ! tau(:,nat) = tau_re
  ! ityp(nat) = 4
  ! !
  ! tau_re = tau(:,171) ! save the iron coordinates
  ! do ia=172,nat-1
  !   tau(:,ia-1) = tau(:,ia)
  !   ityp(ia-1) = ityp(ia)
  ! enddo
  ! tau(:,nat-1) = tau_re
  ! ityp(nat-1) = 4
  !
  !
  !================================================================================
  !    Read eigenmode energies and polarizations
  !================================================================================
  !
  ! Allocate the variables
  allocate( w(3*nat) )
  allocate( polarization(3*nat,3*nat) )
  w=ZERO
  polarization=ZERO
  !
  ! Read them from file
  call read_eigenmodes(polarizations_file)
  !
  !
  !================================================================================
  !    Read matrix elements
  !================================================================================
  !
  ! Allocate the variables
  allocate( ep_mat_el_vl(nbands,3*nat) ) ! ep_mat_el_ia( n, i alpha )
  allocate( ep_mat_el_ia(nbands,3*nat) ) ! ep_mat_el_ia( n, i alpha )
  allocate( ep_mat_el_nl(50,3*nat) ) ! ep_mat_el_ia( n, i alpha )
  ep_mat_el_vl=cmplx_0
  ep_mat_el_ia=cmplx_0
  ep_mat_el_nl=cmplx_0
  !
  ! Read local part of matrix elemnts from file
  call read_ep_mat_elements(ep_mat_el_file)


  ep_mat_el_vl(1274,:) = 0.0
  ep_mat_el_vl(1:1229,:) = 0.0
  ep_mat_el_vl(1328:nbands,:) = 0.0

  print*, "  1-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,795) = 0.0
  print*, "  2-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,797) = 0.0
  print*, "  3-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,796) = 0.0
  print*, "  4-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,798) = 0.0
  print*, "  5-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,799) = 0.0
  print*, "  6-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,800) = 0.0
  print*, "  7-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,792) = 0.0
  print*, "  8-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,793) = 0.0
  print*, "  9-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,794) = 0.0
  print*, " 10-", maxloc(abs(ep_mat_el_vl)), maxval(abs(ep_mat_el_vl))
  ep_mat_el_vl(1257,802) = 0.0

  print*, ""


  ep_mat_el_ia = ep_mat_el_vl
  !
  ! Read non local part
  nonl=.true.
  if (nonl) then
    ep_mat_el_file="ep_mat_is_1_iband_1274_jband_1229-1278_NL.dat_1"
    !
    nl_unit = find_free_unit()
    inquire(iolength=nl_rl) ep_mat_el_nl(:,1)
    open( unit=nl_unit, file=trim(mesh_dir)//trim(ph_dir)//trim(ep_mat_dir)//trim(ep_mat_el_file), iostat=ios, &
          form='unformatted', status='old', access='direct', recl=nl_rl )
    if ( ios /= 0 ) stop "read_g_elements nl: Error opening file"
    !
    do im=1,3*nat
      read( unit=nl_unit, rec=im ) ep_mat_el_nl(:,im)
    enddo
    !
    close(unit=nl_unit)
    !
    ib = 0
    do iband=1229,1278
      ib = ib + 1
      ! do im=1,3*nat
      !   if (abs(ep_mat_el_nl(ib,im))>1) print*, iband, im, ep_mat_el_ia(iband,im), ep_mat_el_ia(iband,im) + real(ep_mat_el_nl(ib,im))
      ! enddo
      ep_mat_el_ia(iband,:) = ep_mat_el_vl(iband,:) + real(ep_mat_el_nl(ib,:))
    enddo
    !
    !
    ep_mat_el_file="ep_mat_is_1_iband_1274_jband_1279-1328_NL.dat_1"
    !
    nl_unit = find_free_unit()
    inquire(iolength=nl_rl) ep_mat_el_nl(:,1)
    open( unit=nl_unit, file=trim(mesh_dir)//trim(ph_dir)//trim(ep_mat_dir)//trim(ep_mat_el_file), iostat=ios, &
          form='unformatted', status='old', access='direct', recl=nl_rl )
    if ( ios /= 0 ) stop "read_g_elements nl: Error opening file"
    !
    do im=1,3*nat
      read( unit=nl_unit, rec=im ) ep_mat_el_nl(:,im)
    enddo
    !
    close(unit=nl_unit)
    !
    ib = 0
    do iband=1279,1328
      ib = ib + 1
      ! do im=1,3*nat
      !   if (abs(ep_mat_el_nl(ib,im))>1) print*, iband, im, ep_mat_el_ia(iband,im), ep_mat_el_ia(iband,im) + real(ep_mat_el_nl(ib,im))
      ! enddo
      ep_mat_el_ia(iband,:) = ep_mat_el_vl(iband,:) + real(ep_mat_el_nl(ib,:))
    enddo
  endif

  !
  tot=0
  haundiak=0
  do im=1,3*nat
    do ib=1229,1328
      tot = tot + 1
      if (abs(ep_mat_el_ia(ib,im))>1) haundiak = haundiak + 1
    enddo
  enddo
  print*, tot, haundiak, real(haundiak,dp)/tot*100
  print*, ""


  print*, "  1-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(4,818) = 0.0
  print*, "  2-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(3,817) = 0.0
  print*, "  3-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(4,817) = 0.0
  print*, "  4-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(3,818) = 0.0
  print*, "  5-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(50,817) = 0.0
  print*, "  6-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(50,818) = 0.0
  print*, "  7-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(17,311) = 0.0
  print*, "  8-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(17,719) = 0.0
  print*, "  9-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(16,565) = 0.0
  print*, " 10-", 1229-1+maxloc(abs(ep_mat_el_nl)), maxval(abs(ep_mat_el_nl))
  ep_mat_el_nl(16,463) = 0.0

  print*, ""

  stop
  !
  !
  !================================================================================
  !   g_m(n,nu) = sum_{ia,alpha} epsil({ia,alpha}, nu) * g_m(n,{ia,alpha}) / sqrt( 2*M(ia)*w(nu) )
  !================================================================================
  !
  ! Allocate variables
  allocate( ep_mat_el_nu(nbands,3*nat) ) ! ep_mat_el_nu( n, i alpha )
  ep_mat_el_nu=cmplx_0
  !
  ! do iband=1,nbands
  do iband=1229,1328
    write(unit=stdout,fmt="(a1,i4,a1,i4)",advance="no") CR, iband, "/", nbands
    flush(stdout)
    write(unit=stdout,fmt="(a1)", advance="no") LF
    do imode=1,3*nat
      !
      do jmode=1,3*nat ! sum_{ia,alpha}
        !
        ia = (jmode-1)/3 + 1
        alpha = mod(jmode-1,3) + 1
        !
        ep_mat_el_nu(iband,imode) = ep_mat_el_nu(iband,imode) + &
                                    ep_mat_el_ia(iband,jmode) * Ry_to_eV * polarization(jmode,imode) / &
                                    sqrt( 2.0_dp * amass(ityp(ia))*pmass * w(imode)*eV_to_Ha )
        !
      enddo
      !
    enddo
  enddo
  !
  ! ep_mat_el_nu(1:1229,:) = cmplx_0
  ! ep_mat_el_nu(1328:nbands,:) = cmplx_0
  !
  write(unit=1234,fmt="(f20.16)") abs(ep_mat_el_nu)
  !
  print*, "1-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
  ep_mat_el_nu(1274,9) = 0.0
  print*, "2-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
  ep_mat_el_nu(1274,4) = 0.0
  print*, "3-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
  ep_mat_el_nu(1274,10) = 0.0
  print*, "4-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
  ep_mat_el_nu(1274,6) = 0.0
  print*, "5-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
  ep_mat_el_nu(1274,5) = 0.0
  !
  ep_mat_el_nu(1274,:) = 0.0
  print*, ""
  !
  if (nonl) then
    print*, "1-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1282,23) = 0.0
    print*, "2-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1232,23) = 0.0
    print*, "3-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1257,59) = 0.0
    print*, "4-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1257,17) = 0.0
    print*, "5-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1277,59) = 0.0
  else
    print*, "1-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1257,59) = 0.0
    print*, "2-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1257,17) = 0.0
    print*, "3-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1277,59) = 0.0
    print*, "4-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1257,16) = 0.0
    print*, "5-", maxloc(abs(ep_mat_el_nu)), maxval(abs(ep_mat_el_nu))
    ep_mat_el_nu(1257,129) = 0.0
  endif

  stop

  !
  !
  !================================================================================
  !   Calculate lambda parameter
  !================================================================================
  !
  stop
  !
  sigma = 0.0001_dp
  !
  lambda = 0.0_dp
  !
  do imode=1,3*nat
    !
    if ( abs(w(imode)) < 0.0001_dp ) then
      ! ep_mat_el_nu(:,imode) = 0.0_dp
      cycle
    endif
    !
    do iband=1,nbands
      g = 1.0_dp/(sigma*sqrt(2.0_dp*3.1416_dp)) * exp( -0.5_dp*((e(initial)-e(iband)+w(imode))/sigma)**2 )
      lambda = lambda + ep_mat_el_nu(iband,imode)**2/w(imode) * g
      g = 1.0_dp/(sigma*sqrt(2.0_dp*3.1416_dp)) * exp( -0.5_dp*((e(initial)-e(iband)-w(imode))/sigma)**2 )
      lambda = lambda + ep_mat_el_nu(iband,imode)**2/w(imode) * g
    enddo
    !
  enddo
  !
  !
contains

  subroutine read_ep_mat_elements(file_name)
    !
    ! This subroutine reads matrix elements form a file
    !
    ! Import modules
    use intw_utility, only:  find_free_unit
    !
    implicit none
    !
    ! Input variables
    character(len=256), intent(in) :: file_name
    !
    ! Local variables
    !
    ! loop variables
    integer :: im
    ! I/O variables
    integer :: record_length
    integer :: file_unit
    integer :: ios
    !
    !
    ! Open the file
    file_unit = find_free_unit()
    inquire(iolength=record_length) ep_mat_el_vl(:,1)
    open( unit=file_unit, file=trim(mesh_dir)//trim(ph_dir)//trim(ep_mat_dir)//trim(file_name), iostat=ios, &
          form='unformatted', status='old', access='direct', recl=record_length )
    if ( ios /= 0 ) stop "read_g_elements: Error opening file"
    !
    do im=1,3*nat
      read( unit=file_unit, rec=im ) ep_mat_el_vl(:,im)
    enddo
    !
    ! Colse the file
    close(unit=file_unit)
    !
    !
  end subroutine read_ep_mat_elements



  subroutine read_eigenmodes(file_name)
    !
    ! This subroutine reads the energies and polarizations of all vibrational
    ! eigenmodes of the sistem from a Jmol suited .xyz file.
    !
    ! Import modules
    use intw_utility, only:  find_free_unit
    !
    implicit none
    !
    ! Input variables
    character(len=256), intent(in) :: file_name
    !
    ! Local variables
    !
    ! to read the file
    integer :: nat_file
    integer :: nmode_file
    real(kind=dp), dimension(3) :: tau_file
    character(len=256) :: dummy_file
    ! loop variables
    integer :: im
    integer :: ia
    integer :: id
    ! I/O variables
    integer :: file_unit
    integer :: ios
    ! to normalize the polarization vector
    real(kind=dp) :: norm
    !
    !
    ! Open the file
    file_unit = find_free_unit()
    ! open(unit=file_unit, file=trim(mesh_dir)//trim(ph_dir)//trim(file_name), iostat=ios)
    open(unit=file_unit, file=trim(file_name), iostat=ios)
    if ( ios /= 0 ) stop "read_eigenmodes: Error opening file"
    !
    ! For each mode
    do im=1,3*nat
      !
      ! Read the number of atoms
      read(file_unit,"(i)") nat_file
      ! if (nat_file /= nat) stop "read_eigenmodes: The number of atoms does not match"
      !
      ! Read the mode number and energy
      ! Mode #   1 f =  -0.290132 meV.
      read(file_unit,"(a6,i4,a5,1f10.6,a5)") dummy_file(1:6), nmode_file, dummy_file(7:11), w(im),dummy_file(12:16)
      w(im) = w(im)/1000.0_dp ! from meV to eV
      ! print*, "haritz: nmode_file, w: ", nmode_file, w(im)
      !
      ! Read polarizations
      do ia=1,nat
        !
        read(file_unit,"(a4,x,2(3f15.10))") dummy_file(1:5), ( tau_file(id) ,id=1,3 ), ( polarization(3*(ia-1)+id,im), id=1,3 )
        !
        tau_file = tau_file/bohr/alat ! angstrom to bohr to alat
        ! print"(i4,a4,3f20.16)", ia, trim(dummy_file(1:5)), tau_file
        ! print"(i4,a4,3f20.16)", ia, trim(atom_labels(ityp(ia))), tau(:,ia)
        ! if ( sqrt( sum( ( tau_file(:)-tau(:,ia) )**2 ) ) > 1.0E-2 )  stop "read_eigenmodes: The coordinates of the atoms do not match"
        !
      enddo
      !
      norm = sqrt( sum( ( polarization(:,im) )**2 ) )
      polarization(:,im) = polarization(:,im) / norm
      !
    enddo
    !
    ! Close the file
    close(unit=file_unit)
    !
    !
  end subroutine read_eigenmodes

end program lambda_prog
