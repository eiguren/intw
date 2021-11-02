program dv_hand
  !
  ! This program reads the data of the supercell frozen-phonon calculation to
  ! calculate the potential change due to the phonon.
  !
  use kinds, only: dp
  !
  use intw_input_parameters, only: mesh_dir, prefix
  !
  use intw_reading, only: at, alat, nat, ntyp, tau, nr1, nr2, nr3, ityp, nspin
  !
  use intw_reading, only: read_parameters_data_file_xml
  !
  use intw_utility, only: find_free_unit
  !
  use intw_haritz, only: get_gamma_only_haritz, get_spin_haritz, stdout
  !
  use fft_base, only: dfftp
  !
  use ions_base, only: zv, atm
  !
  implicit none
  !
  ! Potentials and change of potential (in real space)
  real(kind=dp), allocatable, dimension(:) :: V_total_p, V_total_m
  real(kind=dp), allocatable, dimension(:) :: V_ionic_p, V_ionic_m
  real(kind=dp), allocatable, dimension(:,:) :: delta_V
  !
  ! Loop variables
  integer :: ia, ia_move, id, imode, is
  !
  ! Variables needed to read the potential with plot_io
  character(len=1024) :: scf_file_p, scf_file_m
  character(len=256) :: dummy1, dummy2
  character(len=1), dimension(3) :: disp
  integer :: plot_num, ibrav
  real(kind=dp) :: ecutwfc, dual, gcutm
  real(kind=dp), dimension(6) :: celldm
  character (len=75) :: title
  !
  ! Output file variables
  integer :: ounit, ios, rl_delta_V
  !
  ! Input parameters
  real(kind=dp) :: dx
  integer :: nat_move ! The number of atoms moved in the calculation
  integer, dimension(100) :: ia_move_to_ia


  20 format(A)
  !
!#ifdef DEBUG
  2100 FORMAT('DEBUG: ',a15,x,i4)
  2200 FORMAT('DEBUG: ',a15,x,3i4)
  2300 FORMAT('DEBUG: ',a15,x,f16.8)
  2400 FORMAT('DEBUG: ',a15,x,3f16.8)
!#endif
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
  !
  !================================================================================
  !       Specify input parameters
  !================================================================================
  !
  call read_input_file()
  !
  write(stdout,20) 'Reading the QE calculation data of v0'
  !
  !================================================================================
  !       read the parameters from the SCF QE calculation
  !================================================================================
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
  !
  ! Remember that this subroutine reads nspin=2 only if noncolin=.true. !!!!!!!
  !
  call read_parameters_data_file_xml()
  !
  call get_gamma_only_haritz()
  !
  call get_spin_haritz()
  !
  !
  ! Calculate the potential change dV
  !
  disp(1) = 'x'
  disp(2) = 'y'
  disp(3) = 'z'
  !
  write(stdout,20) ''
  write(stdout,20) '======================================================'
  write(stdout,20) ''
  write(stdout,20) 'Calculate the potential change'
  write(stdout,20) ''
  !

  allocate( V_total_p(nr1*nr2*nr3), V_total_m(nr1*nr2*nr3) )
  allocate( V_ionic_p(nr1*nr2*nr3), V_ionic_m(nr1*nr2*nr3) )
  allocate( delta_V(nr1*nr2*nr3,nspin) )
  !
  inquire(iolength=rl_delta_V) delta_V(:,:)
  print*, rl_delta_V
  !
  ! Open the file where the dv will be written
  ounit = find_free_unit()
  open ( unit=ounit, file = 'fe-mgo-ag100.dvscf_hand_only', iostat = ios, &
         form = 'unformatted', status = 'new', access = 'direct', &
         action = 'write', recl = rl_delta_V )
  !
  !
  do ia_move=1,nat_move ! for each atom moved in the unit cell
    !
    ! identify the atom of the unit cell
    ia = ia_move_to_ia(ia_move)
    !
    do id=1,3 ! for each direction
      !
      ! identify the vibrational mode
      imode = id + 3*(ia_move-1)
      !
      do is=1,nspin ! for each spin
        !
        !
        write(dummy1,"(i3,a1,a1)") ia, disp(id), "+"
        write(dummy2,"(i3,a1,a1)") ia, disp(id), "-"
        !
        ! Read total potential for + and - displacaments
        if (nspin==1) then
          scf_file_p = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy1))//'/'//'total_pot'
          scf_file_m = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy2))//'/'//'total_pot'
        elseif (nspin==2) then
          if (is==1) then
            scf_file_p = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy1))//'/'//'total_pot_s1'
            scf_file_m = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy2))//'/'//'total_pot_s1'
          elseif (is==2) then
            scf_file_p = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy1))//'/'//'total_pot_s2'
            scf_file_m = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy2))//'/'//'total_pot_s2'
          endif
        endif
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
        scf_file_p = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy1))//'/'//'ionic_pot'
        scf_file_m = trim(adjustl(mesh_dir))//'/../'//trim(adjustl(dummy2))//'/'//'ionic_pot'
        !
        ! read everything (requires that all variables that are read are allocated
        ! with the correct dimensions!)
        CALL plot_io (trim(scf_file_p), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
              dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
              gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_ionic_p, - 1)
        CALL plot_io (trim(scf_file_m), title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
              dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
              gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, V_ionic_m, - 1)

        delta_V(:,is) = ( (V_total_p - V_ionic_p) - (V_total_m - V_ionic_m) ) / ( 2.0_dp * dx*alat )
        !
        !
        !
      enddo ! is
      !
      write(unit=ounit, rec=imode, iostat= ios ) delta_V(:,:)
      !
    enddo ! id
    !
  enddo ! ia_move
  !
  !
  close(ounit)
  !
  !


contains

  subroutine read_input_file()
    !
    ! This subroutine sets default values for variables of input namelist
    ! and then reads them from input file if there are specified
    !
    use intw_haritz, only: stdin
    !
    ! Variables to read
    use intw_input_parameters, only: mesh_dir, prefix
    !
    implicit none
    !
    character(len=256) :: input_file_name
    ! Namelists
    !
    ! input variables namelist
    NAMELIST / input_dv_hand / mesh_dir, prefix, dx, nat_move, ia_move_to_ia
    !
    ! Set default values
    mesh_dir = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/RELAX/Fe1/ph/v0/'
    prefix = 'Fe-MgO-Ag100'
    nat_move = 0
    ia_move_to_ia = 0
    dx = 0.001_dp
    !
    INQUIRE(stdin, NAME=input_file_name)
    if (input_file_name(1:4)=="/dev") then
      ! there is no input file: use default values
      return
    else
      read(stdin,input_dv_hand)
    endif
    !
    ! Reopen terminal as input file
    ! close(unit=stdin)
    ! open(unit=stdin,file='/dev/tty')
    !
    !
  end subroutine read_input_file


end program dv_hand
