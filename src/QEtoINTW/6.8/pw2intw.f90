!-----------------------------------------------------------------------
PROGRAM pw2intw
  !-----------------------------------------------------------------------
  !
  USE Kinds
  USE mp_global, ONLY: mp_startup

  USE io_files, ONLY: prefix_QE => prefix
  USE environment, ONLY: environment_start, environment_end
  USE basis, ONLY: starting_wfc

  IMPLICIT NONE

  EXTERNAL :: errore, read_file, wfcinit, cryst_to_cart, davcio, system

  ! I/O
  CHARACTER(len=256) :: prefix = " "
  CHARACTER(len=256) :: mesh_dir = "./"
  CHARACTER(len=256) :: data_dir = "./"
  CHARACTER(len=256) :: rho_dir = "./"
  LOGICAL :: phonons = .false.
  CHARACTER(len=256) :: qlist_file = "qlist.txt"
  CHARACTER(len=256) :: ph_dir = "./"
  CHARACTER(len=256) :: dvscf_dir = "./"
  integer :: nqirr = 0
  logical :: dynxml = .false.
  CHARACTER(len=256) :: fildyn = "matdyn"
  CHARACTER(len=256) :: fildvscf = "dvscf"

  CHARACTER(len=256) :: intwdir
  INTEGER :: iostat, strlen

  NAMELIST / inputpp / prefix, mesh_dir, phonons, data_dir, dvscf_dir, rho_dir, qlist_file, nqirr, dynxml, ph_dir, fildyn, fildvscf


  !
  ! Initialize MPI and environment
#if defined(__MPI)
  CALL mp_startup( )
#endif

  CALL environment_start( "pw2intw" )

  !
  ! Read input file
  READ (5, inputpp, iostat=iostat)
  if (iostat /= 0) call errore( "pw2intw", "ERROR: pw2intw: error reading inputpp", iostat )

  strlen = len_trim(mesh_dir)
  if ( mesh_dir(strlen:strlen+1) .ne. "/" ) mesh_dir(strlen+1:strlen+2) = "/"
  strlen = len_trim(data_dir)
  if ( data_dir(strlen:strlen+1) .ne. "/" ) data_dir(strlen+1:strlen+2) = "/"
  strlen = len_trim(dvscf_dir)
  if ( dvscf_dir(strlen:strlen+1) .ne. "/" ) dvscf_dir(strlen+1:strlen+2) = "/"
  strlen = len_trim(rho_dir)
  if ( rho_dir(strlen:strlen+1) .ne. "/" ) rho_dir(strlen+1:strlen+2) = "/"
  strlen = len_trim(ph_dir)
  if ( ph_dir(strlen:strlen+1) .ne. "/" ) ph_dir(strlen+1:strlen+2) = "/"

  !
  ! Create intwdir
  intwdir = trim(mesh_dir)//trim(prefix)//".save.intw/"

  call system("mkdir -p "//intwdir)

  !
  ! Read QE data
  prefix_QE = prefix
  starting_wfc = "file"
  CALL read_file()
  CALL wfcinit()

  !
  ! Write all data in intw format
  call write_pp_intw()

  call write_crystal_info_and_bands ()

  call scf_v_and_rho()

  call write_FFT_information ()

  call write_wfc()

  if (phonons) call write_phonon_info()

  !
  ! End job
  call environment_end( "pw2intw" )

contains


  SUBROUTINE write_pp_intw()

    USE ions_base, ONLY: nsp
    USE uspp_param, ONLY: upf

    implicit none

    integer :: io_unit, ios
    integer :: is, ir, nb, nb1
    integer, external :: find_free_unit
    character(len=256) :: pp_file, datafile

    100 format(i1"-KBPP.txt")
    200 format(i2"-KBPP.txt")

    do is=1,nsp

      if (is<9) then
        write(pp_file,100) is
      else if ((is>9).and.(is<19) ) then
        write(pp_file,200) is
      else
        print*, "ERROR: The num. of species is bigger than 19"
      end if
      io_unit= find_free_unit()
      datafile = trim(intwdir)//trim(pp_file)
      open(unit=io_unit, file=datafile, status="unknown", iostat= ios)
      if (ios /= 0) call errore( "pw2intw", "write_pp_intw: error opening pp_file", iostat )

      write(unit=io_unit,fmt="(a)")"ATOM LABEL"
      write(unit=io_unit,fmt=*) upf(is)%psd

      !TODO: Hau zuzendu, batzuetan ez du idazten
      write(unit=io_unit,fmt="(a)")"IS RELAT."
      write(unit=io_unit,fmt=*) upf(is)%has_so

      write(unit=io_unit,fmt="(a)")"HAS SO COUPLING"
      write(unit=io_unit,fmt=*) upf(is)%has_so

      write(unit=io_unit,fmt="(a)")"NLCC"
      write(unit=io_unit,fmt=*) upf(is)%nlcc

      write(unit=io_unit,fmt="(a)")"Z VALENCE"
      write(unit=io_unit,fmt=*) upf(is)%zp

      write(unit=io_unit,fmt="(a)")"L LOC"
      write(unit=io_unit,fmt=*) upf(is)%lloc

      write(unit=io_unit,fmt="(a)")"LMAX"
      write(unit=io_unit,fmt=*) upf(is)%lmax

      write(unit=io_unit,fmt="(a)")"NBETA"
      write(unit=io_unit,fmt=*) upf(is)%nbeta

      write(unit=io_unit,fmt="(a)")"NUMBER of GRID POINTS FOR EACH BETA"
      write(unit=io_unit,fmt=*) upf(is)%kbeta(1:upf(is)%nbeta)

      write(unit=io_unit,fmt="(a)")"L FOR EACH BETA"
      write(unit=io_unit,fmt="(100i4)") upf(is)%lll(1:upf(is)%nbeta)

      write(unit=io_unit,fmt="(a)")"J FOR EACH BETA (IF SO COUPLING)"
      if (upf(is)%has_so) then
        write(unit=io_unit,fmt="(100f6.2)") upf(is)%jjj(1:upf(is)%nbeta)
      else
        write(unit=io_unit,fmt=*)"NOTHING"
      end if
      write(unit=io_unit,fmt="(a)")"N OF MESH POINTS"
      write(unit=io_unit,fmt=*) upf(is)%mesh

      write(unit=io_unit,fmt="(a)")"RC_LOC"
      write(unit=io_unit,fmt=*) upf(is)%rcloc

      write(unit=io_unit,fmt="(a)")"KB D{ij}"
      do nb=1,upf(is)%nbeta
      write(unit=io_unit,fmt="(100es18.8)") (upf(is)%dion(nb,nb1),nb1=1,upf(is)%nbeta)
      end do

      write(unit=io_unit,fmt="(a)")"    R                 RAB              VLOC"//&
                                  "              BETA1              BETA2             .."
      do ir=1,upf(is)%mesh
        write(unit=io_unit,fmt="(100es18.8)") &
      upf(is)%r(ir), upf(is)%rab(ir), upf(is)%vloc(ir),(upf(is)%beta(ir,nb), nb=1,upf(is)%nbeta)
      end do

    end do

  end SUBROUTINE write_pp_intw


  SUBROUTINE scf_v_and_rho()
    !ASIER: 22/02/2022
    !Here we asume that:
    ! nspin=1 (no magnetism at all)
    ! nspin=2 (collinear)
    ! nspin=4 (non-collinear)
    ! It is clear that we must fix this, because
    ! in many places we asume that npol=2 is non-collinear.

    USE lsda_mod, ONLY: nspin
    USE scf, ONLY: rho, v

    implicit none

    integer :: io_unit, ios
    integer :: ispin
    integer, external :: find_free_unit

    character(len=256) :: datafile
    integer :: record_length

    io_unit = find_free_unit()
    datafile = trim(intwdir)//"scf_vr.dat"
    inquire(iolength=record_length) v%of_r(:,1)
    open(unit=io_unit, file=datafile, iostat=ios, &
        status="unknown", action="write", form="unformatted", access="direct", recl=record_length)
    if (ios /= 0) call errore( "pw2intw", "scf_v_and_rho: error opening scf_vr.dat", ios )
    !
    do ispin=1, nspin
      write (unit=io_unit, rec=ispin) v%of_r(:,ispin)
    end do
    !
    close(unit=io_unit)

    !
    io_unit = find_free_unit()
    datafile = trim(intwdir)//"scf_rhor.dat"
    inquire(iolength=record_length) rho%of_r(:,1)
    open(unit=io_unit, file=datafile, iostat=ios, &
        status="unknown", action="write", form="unformatted", access="direct", recl=record_length)
    if (ios /= 0) call errore( "pw2intw", "scf_v_and_rho: error opening scf_rhor.dat", ios )
    !
    do ispin=1, nspin
      write (unit=io_unit,rec=ispin) rho%of_r(:,ispin)
    end do
    !
    close(unit=io_unit)

  END SUBROUTINE scf_v_and_rho


  subroutine write_phonon_info()

    implicit none

    !
    ! Irreducible patterns
    call write_paterns()

    !
    ! Induced potentials
    call write_dV()

    !
    ! Dynamical matrices
    if (dynxml) call errore( "pw2intw", "write_phonon_info: dynxml not implemented yet", 1 )
    call write_dyn()

    !
    ! Induced charge density
    call write_drho()

  end subroutine write_phonon_info


  subroutine write_paterns()

    USE xmltools, ONLY: xml_openfile, xml_closefile, &
                        xmlr_opentag, xmlr_readtag, xmlr_closetag, i2c
    USE ions_base, ONLY: nat

    implicit none

    ! loop variables
    integer :: imode, jmode, iq
    ! irreducible patterns variables
    complex(kind=dp) :: u_irr(1:3*nat,1:3*nat)
    integer :: nirr, ipert, imode0, irr
    integer, dimension(48) :: npert=-1
    ! i/o variables
    character(len=256) :: datafile, q_dir
    integer :: io_unit_write, io_unit_read, ios
    logical :: existitu
    integer, external :: find_free_unit


    ! Open file for writing
    io_unit_write = find_free_unit()
    datafile = trim(intwdir)//"/irrq_patterns.dat"
    open(unit=io_unit_write, file=datafile, status="replace", form="formatted", iostat=ios)
    if ( ios /= 0 ) call errore( "pw2intw", "write_phonon_info: error opening irrq_patterns.dat", ios )
    !
    do iq=1,nqirr
      !
      call write_tag("qq", iq, q_dir)
      !
      ! Read irreducible patterns
      datafile = trim(ph_dir)//trim(q_dir)//"/_ph0/"//trim(prefix)//".phsave/patterns.1.xml"
      inquire(file=datafile, exist=existitu)
      if (.not. existitu) call errore( "pw2intw", "write_phonon_info: patterns.1.xml not found", 1 )
      !
      io_unit_read = xml_openfile( datafile )
      !
      CALL xmlr_opentag( "IRREPS_INFO" )
      CALL xmlr_readtag( "NUMBER_IRR_REP", nirr )
      imode0 = 0
      DO irr = 1, nirr
        !
        CALL xmlr_opentag( "REPRESENTION."//i2c(irr) )
        CALL xmlr_readtag( "NUMBER_OF_PERTURBATIONS", npert(irr) )
        DO ipert = 1, npert(irr)
          imode = imode0 + ipert
          CALL xmlr_opentag( "PERTURBATION."// i2c(ipert) )
          CALL xmlr_readtag( "DISPLACEMENT_PATTERN", u_irr(:,imode) )
          CALL xmlr_closetag( )
        ENDDO
        imode0 = imode0 + npert(irr)
        CALL xmlr_closetag( )
        !
        END DO !irr
        !
      CALL xmlr_closetag("IRREPS_INFO")
      CALL xml_closefile ( )
      !
      ! Write irreducible patterns
      write(io_unit_write,"(a,i4)") "q", iq
      do imode=1,3*nat
        write(io_unit_write,"(10000(a,f16.10,a,f16.10,a))") &
          ("(", real(u_irr(jmode,imode),dp), ",", aimag(u_irr(jmode,imode)), ") ", jmode=1,3*nat)
      enddo
      !
    enddo !iq
    !
    close(unit=io_unit_write)

  end subroutine write_paterns

  subroutine write_dV()

    USE ions_base, ONLY: nat
    USE fft_base, ONLY: dfftp
    USE spin_orb, ONLY: domag
    USE noncollin_module, ONLY: noncolin
    USE lsda_mod, ONLY: lsda

    implicit none

    ! loop variables
    integer :: imode, iq
    ! induced potential variables
    complex(kind=dp), allocatable :: dvq(:, :)
    ! i/o variables
    character(len=256) :: datafile, q_dir, dv_file
    integer :: rl, io_unit_write, io_unit_read, ios
    integer, external :: find_free_unit


    ! Allocate dvq array
    if (lsda) then
      allocate(dvq(dfftp%nr1*dfftp%nr2*dfftp%nr3, 2))
    else
      if (noncolin .and. domag) then
        allocate(dvq(dfftp%nr1*dfftp%nr2*dfftp%nr3, 4))
      else
        allocate(dvq(dfftp%nr1*dfftp%nr2*dfftp%nr3, 1))
      endif
    endif
    dvq = (0.0_dp, 0.0_dp)
    inquire(iolength=rl) dvq(:, :)
    !
    do iq=1,nqirr
      !
      call write_tag("qq", iq, q_dir)
      !
      ! Open file for reading
      io_unit_read = find_free_unit()
      datafile = trim(ph_dir)//trim(q_dir)//"/_ph0/"//trim(prefix)//"."//trim(fildvscf)//"1"
      open( unit=io_unit_read, file=trim(datafile), iostat=ios, &
            form="unformatted", status="old", action="read", access="direct", recl=rl )
      if ( ios /= 0 ) call errore( "pw2intw", "write_phonon_info: error opening dvscf file to read", ios )
      !
      ! Open file for writing
      io_unit_write = find_free_unit()
      call write_tag(trim(prefix)//".dvscf_q", iq, dv_file)
      datafile = trim(intwdir)//trim(dv_file)
      open( unit=io_unit_write, file=trim(datafile), iostat=ios, &
            form="unformatted", status="replace", action="write", access="direct", recl=rl )
      if ( ios /= 0 ) call errore( "pw2intw", "write_phonon_info: error opening dvscf file to write", ios )
      !
      do imode=1,3*nat
        !
        read(unit=io_unit_read, rec=imode, iostat=ios) dvq(:, :)
        !
        if (lsda) then
          ! Transform the colinear potential (rho_up, rho_down) into a 2x2 non-colinear potential (rho, m_x, m_y, m_z) for compatibility with intw
          write(unit=io_unit_write, rec=2*(imode-1)+1, iostat=ios) (dvq(:, 1)+dvq(:, 2))/2.0_dp, 0.0_dp*dvq(:, 1)
          write(unit=io_unit_write, rec=2*(imode-1)+2, iostat=ios) 0.0_dp*dvq(:, 2), (dvq(:, 1)-dvq(:, 2))/2.0_dp
        else
          write(unit=io_unit_write, rec=imode, iostat=ios) dvq(:, :)
        endif
        !
      enddo
      !
      close(unit=io_unit_read)
      close(unit=io_unit_write)
      !
    enddo

  end subroutine write_dV


  subroutine write_dyn()

    USE cell_base, ONLY: at
    USE ions_base, ONLY: nat, nsp, tau

    implicit none

    ! loop variables
    integer :: iq
    ! i/o variables
    character(len=256) :: datafile, q_dir, dyn_file
    integer :: io_unit_write, io_unit_read, ios
    logical :: existitu
    integer, external :: find_free_unit
    character(len=100) :: comentario
    integer :: ntyp_, nat_, ibrav_
    integer :: ia, ja, it, i, j
    real(dp) :: celldm_(6), at_(3,3), fracpos(3)
    real(dp) :: qpoint_cart(3), qpoint_cryst(3)
    real(dp) ::  dynq_re(3,3), dynq_im(3,3)
    ! TODO add to useful_constants
    real(dp), parameter :: Ry2Hartree = 0.5_dp


    do iq=1,nqirr
      !
      call write_tag("qq", iq, q_dir)
      !
      !!!!!!!!!!!!!!! JUST COPY DYN FILES !!!!!!!!!!!!!!!
      datafile = trim(ph_dir)//trim(q_dir)//"/"//trim(fildyn)
      inquire(file=datafile, exist=existitu)
      if (.not. existitu) call errore( "pw2intw", "write_phonon_info: fildyn not found: check fildyn input variable", 1 )
      !
      call write_tag(trim(prefix)//".dyn_q", iq, dyn_file)
      call system("cp " // &
                  trim(datafile)//" " // &
                  trim(intwdir)//trim(dyn_file) )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !!!! READ AND WRITE ONLY IRREDUCIBLE Q DYN MAT !!!!
      !
      ! Open file for reading
      io_unit_read = find_free_unit()
      datafile = trim(ph_dir)//trim(q_dir)//"/"//trim(fildyn)
      open( unit=io_unit_read, file=trim(datafile), iostat=ios, status="old", action="read" )
      if ( ios /= 0 ) call errore( "pw2intw", "write_phonon_info: error opening fildyn file to read", ios )
      !
      ! Open file for writing
      io_unit_write = find_free_unit()
      call write_tag(trim(prefix)//".dyn_q", iq, dyn_file)
      datafile = trim(intwdir)//trim(dyn_file)//"_sym"
      open( unit=io_unit_write, file=trim(datafile), iostat=ios, status="replace", action="write" )
      if ( ios /= 0 ) call errore( "pw2intw", "write_phonon_info: error opening fildyn file to write", ios )
      !
      !
      ! Read QE dyn file header: Information about lattice and atoms (not used but checked for
      ! consistency in the parameters)
      !
      ! Dummy lines
      read(io_unit_read,"(a)") comentario
      read(io_unit_read,"(a)") comentario
      ! Number of atomic species, number of atoms...
      read(io_unit_read,*) ntyp_, nat_, ibrav_, celldm_
      if ( ntyp_ /= nsp .or. nat_ /= nat ) call errore( "pw2intw", "write_dyn: Error reading parameters in dyn file", ios )
      !
      ! Cell parameters
      if ( ibrav_ == 0 ) then
        read(io_unit_read,"(a)") comentario ! symm_type not used
        read(io_unit_read,*) ((at_(i,j),i=1,3),j=1,3)
        if ( any( abs(at_-at)>1.d-8 ) ) call errore( "pw2intw", "write_dyn: Error reading parameters in dyn file", ios )
      end if
      !
      ! Atomic species
      do it=1,nsp
        read(io_unit_read,"(a)") comentario !i, atom, amass, not used
      end do
      !
      ! Atomic positions
      do ia = 1,nat
        read(io_unit_read,*) i,j, fracpos(:)
        if ( any( abs(fracpos-tau(:,ia))>1.d-8 ) ) call errore( "pw2intw", "write_dyn: Error reading parameters in dyn file", ios )
      end do
      !
      ! Dummy line
      read(io_unit_read,"(/,a,/)") comentario
      !
      !
      ! Read irreducible q point
      read(io_unit_read,"(11x, 3(f14.9), / )") qpoint_cart
      !
      ! Transform from cartesian to crystalline
      qpoint_cryst = qpoint_cart
      call cryst_to_cart (1, qpoint_cryst, at, -1)
      !
      ! Write irreducible q point
      write(io_unit_write,"(a11,3f14.9,a2)") "q_cryst = ( ",  qpoint_cryst, ")"
      !
      !
      ! Read dynamical matrix of the irreducible q point ONLY and write it
      do ia = 1,nat
        do ja = 1,nat
          !
          read(io_unit_read,*) i, j ! not used
          !
          ! cartesian 3x3 block of this atom pair in dynq matrix (without the mass factor)
          read(io_unit_read,*) ((dynq_re(i,j), dynq_im(i,j), j=1,3), i=1,3) ! in Ry/Bohr^2
          write(io_unit_write, "(3(a,f16.10,a,f16.10,a))") (("(", dynq_re(i,j)*Ry2Hartree, &
                                                             ",", dynq_im(i,j)*Ry2Hartree, ") ", j=1,3), i=1,3) ! in a.u.
          !
        end do ! ja
      end do ! ia
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    enddo ! iq

  end subroutine write_dyn

  subroutine write_drho()

    implicit none

    ! loop variables
    integer :: iq
    ! i/o variables
    character(len=256) :: datafile, q_dir, rho_file
    logical :: existitu


    do iq=1,nqirr
      !
      call write_tag("qq", iq, q_dir)
      !
      datafile = trim(ph_dir)//trim(q_dir)//"/_ph0/"//trim(prefix)//".rho1"
      inquire(file=datafile, exist=existitu)
      if (existitu) then
        call write_tag(trim(prefix)//".rho_q", iq, rho_file)
        call system("cp " // &
                    trim(datafile)//" " // &
                    trim(intwdir)//trim(rho_file) )
      end if
      !
    enddo

  end subroutine write_drho


  SUBROUTINE write_fft_information()

    USE cell_base, ONLY: at
    USE gvect, ONLY: ngm, g
    USE wvfct, ONLY: npwx
    USE klist, ONLY: nks, igk_k, ngk

    IMPLICIT NONE

    INTEGER, EXTERNAL :: find_free_unit

    integer :: io_unit, ios
    integer :: ig
    CHARACTER(LEN=256) :: datafile
    integer :: ik, i

    io_unit = find_free_unit()
    datafile = trim(intwdir)//"gvectors.dat"
    open(unit=io_unit, file=datafile, iostat=ios, status="unknown", action="write", form="unformatted")
    if (ios /= 0) call errore( "pw2intw", "write_fft_information: error opening gvectors.dat", ios )
    !
    write(unit=io_unit) ngm
    do ig=1,ngm
       write(unit=io_unit) ( nint(sum(at(:,i)*g(:,ig))), i=1,3 )
    end do
    !
    close(unit=io_unit)

    io_unit = find_free_unit()
    datafile = trim(intwdir)//"iGlist.dat"
    open(unit=io_unit, file=datafile, iostat=ios, status="unknown", action="write", form="unformatted")
    if (ios /= 0) call errore( "pw2intw", "write_fft_information: error opening iGlist.dat", ios )
    !
    write(unit=io_unit) npwx
    do ik=1,nks
       write(unit=io_unit) ngk(ik)
       write(unit=io_unit) igk_k(1:ngk(ik),ik)
    end do
    !
    close(unit=io_unit)

  end SUBROUTINE  write_fft_information


  subroutine write_wfc()

    USE wvfct, ONLY: nbnd, npwx
    USE wavefunctions, ONLY: evc
    USE io_files, ONLY: nwordwfc, iunwfc
    USE klist, ONLY: nks
    USE klist, ONLY: igk_k, ngk
    USE noncollin_module, ONLY: noncolin
    USE wvfct, ONLY: nbnd, et
    USE constants, ONLY: rytoev
    USE lsda_mod, ONLY: lsda

    implicit none

    integer :: ik, ibnd, is, npol
    character(256) :: wfc_file, datafile
    integer :: io_unit, ios
    integer, external :: find_free_unit
    !
    complex(dp), allocatable :: evc_down(:, :) !JLB

    100 format("wfc"I5.5".dat")

    io_unit = find_free_unit()

    npol = 1
    if (noncolin) npol=2

    ! JLB
    ! Special case: Collinear spin-polarized calculation
    if (lsda) then
      !
      ! The number of k-points is artificially doubled in QE for this case,
      ! and the spin-up and down bands separated
      allocate(evc_down(npwx*npol,nbnd))
      write(*,*) nbnd, nks, nks/2
      do ik=1,nks/2

        write(wfc_file,100) ik
        datafile = trim(intwdir)//trim(wfc_file)
        open(unit=io_unit, file=datafile, status="unknown", action="write", form="unformatted", iostat=ios)
        if (ios /= 0) call errore( "pw2intw", "write_wfc: error opening wfc_file", ios )
        !
        ! Read spin-up wf contained in ik
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
        ! Read spin-down wf contained in nks/2+ik
        CALL davcio (evc_down, 2*nwordwfc, iunwfc, nks/2+ik, -1 )
        ! Write both up and down wfcs as spinors
        write(unit=io_unit) ngk(ik)
        write(unit=io_unit) igk_k(1:ngk(ik), ik)
        write(unit=io_unit) ( et(ibnd,ik)*rytoev, et(ibnd,nks/2+ik)*rytoev, ibnd=1,nbnd )
        do ibnd=1, nbnd
          write(unit=io_unit) evc(1:ngk(ik), ibnd), 0.0*evc(1:ngk(ik), ibnd)
          write(unit=io_unit) 0.d0*evc_down(1:ngk(ik), ibnd), evc_down(1:ngk(ik), ibnd)
        enddo
        write(*,"(4ES18.6)") evc(1, nbnd), evc(npwx+1, nbnd)
        write(*,"(4ES18.6)") evc_down(1, nbnd), evc_down(npwx+1, nbnd)
        !
        close(unit=io_unit)

      end do !ik
      !
      deallocate(evc_down)
      !
    !end JLB
    else
      !
      do ik=1,nks

        write(wfc_file,100) ik
        datafile = trim(intwdir)//trim(wfc_file)
        open(unit=io_unit, file=datafile, status="unknown", action="write", form="unformatted", iostat=ios)
        if (ios /= 0) call errore( "pw2intw", "write_wfc: error opening wfc_file", ios )
        !
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
        write(unit=io_unit) ngk(ik)
        write(unit=io_unit) igk_k(1:ngk(ik), ik)
        write(unit=io_unit) ( et(ibnd,ik)*rytoev, ibnd=1,nbnd )
        do ibnd=1,nbnd
          ! MBR
          !write(unit=io_unit)evc(1:npol*ngk(ik), ibnd)
          write(unit=io_unit) evc(1:ngk(ik), ibnd), evc(npwx+1:npwx+ngk(ik), ibnd)
        enddo
        !
        close(unit=io_unit)

      end do !ik
      !
    end if

  end subroutine write_wfc


  SUBROUTINE write_crystal_info_and_bands

    USE wvfct, ONLY: nbnd
    USE constants, ONLY: rytoev
    USE klist, ONLY: xk, nkstot
    USE cell_base, ONLY: at, bg, alat
    USE spin_orb, ONLY: lspinorb, domag
    USE symm_base, ONLY: nsym, s, ft, t_rev, ft
    USE fft_base, ONLY: dfftp
    USE gvecw, ONLY: ecutwfc
    USE gvect, ONLY: ecutrho
    USE spin_orb, ONLY: lspinorb
    USE ions_base, ONLY: nat, ityp, nsp, tau, atm, amass
    USE io_files, ONLY: psfile
    USE wvfct, ONLY: npwx
    USE noncollin_module, ONLY: noncolin
    USE control_flags, ONLY: gamma_only
    USE lsda_mod, ONLY: lsda

    IMPLICIT NONE
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    INTEGER ik, i, j

    INTEGER :: isym
    CHARACTER(LEN=256) :: datafile
    integer :: io_unit, ios

    ! JLB new variables to support collinear spin-polarized calculation
    LOGICAL :: noncolin_intw, domag_intw, lsda_intw
    INTEGER :: nkstot_intw, nbnd_intw

    io_unit = find_free_unit()
    datafile = trim(intwdir)//"crystal.dat"
    open(unit=io_unit, file=datafile, status="unknown", action="write", form="formatted", iostat=ios)
    if (ios /= 0) call errore( "pw2intw", "write_crystal_info_and_bands: error opening crystal.dat", ios )

    ! JLB: adapt spin-polarized calculation from QE to spinor in INTW
    ! (see also subroutine write_wfc() below)
    ! Here I define new variables to leave the ones read from QE modules unchanged
    if (lsda) then
      noncolin_intw = .true.
      domag_intw = .true. ! not sure about this one
      lsda_intw = .false.
      nkstot_intw = nkstot / 2
      nbnd_intw = nbnd * 2
    else
      ! Just leave as in QE
      noncolin_intw = noncolin
      domag_intw = domag
      lsda_intw = lsda
      nkstot_intw = nkstot
      nbnd_intw = nbnd
    end if
    ! end JLB

    write(unit=io_unit,fmt=*)"ALAT"
    write(unit=io_unit,fmt="(f16.10)")alat

    write(unit=io_unit,fmt=*)"AT"
    do i=1,3
       write(unit=io_unit, fmt="(3f12.6)")(at(i,j), j=1,3)
    enddo

    write(unit=io_unit,fmt=*)"BG"
    do i=1,3
       write(unit=io_unit, fmt="(3f12.6)")(bg(i,j), j=1,3)
    enddo

    write(unit=io_unit,fmt=*)"FFT GRID"
    write(unit=io_unit,fmt=*)dfftp%nr1,dfftp%nr2, dfftp%nr3

    write(unit=io_unit,fmt=*)"ECUTWFC"
    write(unit=io_unit,fmt=*)ecutwfc

    write(unit=io_unit,fmt=*)"ECUTRHO"
    write(unit=io_unit,fmt=*)ecutrho

    write(unit=io_unit,fmt=*)"LSDA"
    write(unit=io_unit,fmt=*) lsda_intw !lsda

    write(unit=io_unit,fmt=*)"NONCOLIN"
    write(unit=io_unit,fmt=*)noncolin_intw !noncolin

    write(unit=io_unit,fmt=*)"LSPINORB"
    write(unit=io_unit,fmt=*)lspinorb

    write(unit=io_unit,fmt=*)"SPINORB_MAG (KONTUZ, hau aldatu da)"
    write(unit=io_unit,fmt=*)domag_intw !domag

    write(unit=io_unit,fmt=*)"NAT"
    write(unit=io_unit,fmt=*)nat

    write(unit=io_unit,fmt=*)"NTYP"
    write(unit=io_unit,fmt=*)nsp

    write(unit=io_unit,fmt=*)"ATOM_LABELS, MASS AND PP_FILE (1:NTYP)"
    do i=1,nsp
       write(unit=io_unit,fmt="(a,f16.5,x,a)")atm(i), amass(i), trim(psfile(i))
    enddo

    write(unit=io_unit,fmt=*)"POSITIONS (1:NAT)"
    do i=1,nat
       write(unit=io_unit,fmt="(a,i4,3f16.8)")atm(ityp(i)),ityp(i), tau(:,i)
    end do

    write(unit=io_unit,fmt=*)"NSYM"
    write(unit=io_unit,fmt=*)nsym

    do isym=1,nsym
       write(unit=io_unit,fmt="(i8)") isym
       do i=1, 3
          write(unit=io_unit, fmt="(3i8)") (s(i,j,isym), j=1,3)
       enddo
       write(unit=io_unit,fmt="(100f16.10)")  &
            real(ft(1,isym),dp), real(ft(2,isym),dp), real(ft(3,isym),dp)
       write(unit=io_unit,fmt="(48i3)")t_rev (isym)
    enddo

    write(unit=io_unit,fmt=*)"NKS"
    write(unit=io_unit,fmt=*)nkstot_intw !nkstot

    write(unit=io_unit,fmt=*)"GAMMA ONLY"
    write(unit=io_unit,fmt=*) gamma_only

    write(unit=io_unit,fmt=*)"NGMAX"
    write(unit=io_unit,fmt=*)npwx

    write(unit=io_unit,fmt=*)"NBAND"
    write(unit=io_unit,fmt=*)nbnd_intw !nbnd

    close(unit=io_unit)


    io_unit = find_free_unit()
    datafile = trim(intwdir)//"kpoints.dat"
    open(unit=io_unit, file=datafile, status="unknown", action="write", form="formatted", iostat=ios)
    if (ios /= 0) call errore( "pw2intw", "write_crystal_info_and_bands: error openeing kpoints.dat", ios )
    !
    DO ik=1,nkstot_intw !nkstot
      write(unit=io_unit,fmt="(100f16.10)")( xk(i,ik), i = 1, 3 )
      !write(unit=io_unit,fmt="(100f16.10)") ( et(ibnd,ik)*rytoev, ibnd=1,nbnd )
    ENDDO
    !
    close(unit=io_unit)

  END SUBROUTINE write_crystal_info_and_bands


  subroutine write_tag(string,i,tag)
    !-----------------------------------------------
    ! This subroutine creates a character string of
    ! the form "string"integer, where the integer
    ! will be immediately after the end of "string",
    ! without blank spaces.
    !-----------------------------------------------
    implicit none

    integer :: i
    character(*) :: string
    character(256) :: integer_part, tag


    if (i < 10) then
      write(integer_part,100) i
    elseif (i < 100 ) then
      write(integer_part,200) i
    elseif (i < 1000 ) then
      write(integer_part,300) i
    elseif (i < 10000 ) then
      write(integer_part,400) i
    elseif (i < 100000 ) then
      write(integer_part,500) i
    end if

    tag = trim(string)//trim(integer_part)

    100 format(I1)
    200 format(I2)
    300 format(I3)
    400 format(I4)
    500 format(I5)

  end subroutine write_tag

end PROGRAM pw2intw
