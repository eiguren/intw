module pw2intw_write

  use kinds

  use pw2intw_io

  implicit none

  ! functions and subroutines
  public :: copy_nscf_files, write_pp_intw, scf_v_and_rho, &
            write_phonon_info, write_paterns, write_dV, write_dyn, &
            write_drho, write_fft_information, write_wfc, write_crystal_info

  ! variables
  public :: intwdir

  private

  external :: errore

  character(len=256) :: intwdir


contains

  subroutine copy_nscf_files()
    !
    ! Copy scf schema and charge_density for later QE nscf jobs
    !
    implicit none

    CHARACTER(len=256) :: datafile


    datafile = trim(outdir)//trim(prefix)//".save/charge-density.dat"
    call execute_command_line("cp "//trim(datafile)//" "//trim(intwdir) )

    datafile = trim(outdir)//trim(prefix)//".save/data-file-schema.xml"
    call execute_command_line("cp "//trim(datafile)//" "//trim(intwdir) )

  end subroutine copy_nscf_files


  subroutine write_pp_intw()

    use ions_base, only: nsp
    use uspp_param, only: upf

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
      if (ios /= 0) call errore( "pw2intw", "write_pp_intw: error opening pp_file", ios )

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

      close(io_unit)

    end do

  end subroutine write_pp_intw


  subroutine scf_v_and_rho()
    !ASIER: 22/02/2022
    !Here we asume that:
    ! nspin=1 (no magnetism at all)
    ! nspin=2 (collinear)
    ! nspin=4 (non-collinear)
    ! It is clear that we must fix this, because
    ! in many places we asume that npol=2 is non-collinear.

    use lsda_mod, only: nspin
    use scf, only: rho, v

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

  end subroutine scf_v_and_rho

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

    use xmltools, only: xml_openfile, xml_closefile, &
                        xmlr_opentag, xmlr_readtag, xmlr_closetag, i2c
    use ions_base, only: nat

    use pw2intw_utils, only: write_tag

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
      datafile = trim(phdir)//trim(q_dir)//"/_ph0/"//trim(prefix)//".phsave/patterns.1.xml"
      inquire(file=datafile, exist=existitu)
      if (.not. existitu) call errore( "pw2intw", "write_phonon_info: patterns.1.xml not found", 1 )
      !
      io_unit_read = xml_openfile( datafile )
      !
      call xmlr_opentag( "IRREPS_INFO" )
      call xmlr_readtag( "NUMBER_IRR_REP", nirr )
      imode0 = 0
      do irr = 1, nirr
        !
        call xmlr_opentag( "REPRESENTION."//i2c(irr) )
        call xmlr_readtag( "NUMBER_OF_PERTURBATIONS", npert(irr) )
        do ipert = 1, npert(irr)
          imode = imode0 + ipert
          call xmlr_opentag( "PERTURBATION."// i2c(ipert) )
          call xmlr_readtag( "DISPLACEMENT_PATTERN", u_irr(:,imode) )
          call xmlr_closetag( )
        enddo
        imode0 = imode0 + npert(irr)
        call xmlr_closetag( )
        !
      enddo !irr
      !
      call xmlr_closetag("IRREPS_INFO")
      call xml_closefile ( )
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

    use ions_base, only: nat
    use fft_base, only: dfftp
    use spin_orb, only: domag
    use noncollin_module, only: noncolin
    use lsda_mod, only: lsda

    use pw2intw_utils, only: write_tag

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
      datafile = trim(phdir)//trim(q_dir)//"/_ph0/"//trim(prefix)//"."//trim(fildvscf)//"1"
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

    use cell_base, only: at
    use ions_base, only: nat, nsp, tau

    use pw2intw_utils, only: write_tag

    implicit none

    external :: cryst_to_cart

    ! loop variables
    integer :: iq
    ! i/o variables
    character(len=256) :: datafile, q_dir, dyn_file
    integer :: io_unit_write, io_unit_read, ios
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
      !!!! READ AND WRITE ONLY IRREDUCIBLE Q DYN MAT !!!!
      !
      ! Open file for reading
      io_unit_read = find_free_unit()
      datafile = trim(phdir)//trim(q_dir)//"/"//trim(fildyn)
      open( unit=io_unit_read, file=trim(datafile), iostat=ios, status="old", action="read" )
      if ( ios /= 0 ) call errore( "pw2intw", "write_phonon_info: error opening fildyn file to read", ios )
      !
      ! Open file for writing
      io_unit_write = find_free_unit()
      call write_tag(trim(prefix)//".dyn_q", iq, dyn_file)
      datafile = trim(intwdir)//trim(dyn_file)
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
      !
      close(unit=io_unit_read)
      close(unit=io_unit_write)
      !
    enddo ! iq

  end subroutine write_dyn

  subroutine write_drho()

    use pw2intw_utils, only: write_tag

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
      datafile = trim(phdir)//trim(q_dir)//"/_ph0/"//trim(prefix)//".rho1"
      inquire(file=datafile, exist=existitu)
      if (existitu) then
        call write_tag(trim(prefix)//".rho_q", iq, rho_file)
        call execute_command_line("cp " // &
                                  trim(datafile)//" " // &
                                  trim(intwdir)//trim(rho_file) )
      end if
      !
    enddo

  end subroutine write_drho


  subroutine write_fft_information()

    use cell_base, only: at
    use gvect, only: ngm, g
    use wvfct, only: npwx
    use klist, only: nks, igk_k, ngk

    implicit none

    integer, external :: find_free_unit

    integer :: io_unit, ios
    integer :: ig
    character(len=256) :: datafile
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

  end subroutine  write_fft_information


  subroutine write_wfc()

    use wvfct, only: nbnd, npwx
    use wavefunctions, only: evc
    use io_files, only: nwordwfc, iunwfc
    use klist, only: nks
    use klist, only: igk_k, ngk
    use noncollin_module, only: noncolin
    use wvfct, only: nbnd, et
    use constants, only: rytoev
    use lsda_mod, only: lsda

    implicit none

    external :: davcio

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
        call davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
        ! Read spin-down wf contained in nks/2+ik
        call davcio (evc_down, 2*nwordwfc, iunwfc, nks/2+ik, -1 )
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
        call davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
        write(unit=io_unit) ngk(ik)
        write(unit=io_unit) igk_k(1:ngk(ik), ik)
        write(unit=io_unit) ( et(ibnd,ik)*rytoev, ibnd=1,nbnd )
        do ibnd=1,nbnd
          write(unit=io_unit) ( evc((is-1)*npwx+1:(is-1)*npwx+ngk(ik), ibnd), is=1,npol )
        enddo
        !
        close(unit=io_unit)

      end do !ik
      !
    end if

  end subroutine write_wfc


  subroutine write_crystal_info

    use wvfct, only: nbnd
    use constants, only: rytoev
    use klist, only: xk, nkstot
    use cell_base, only: at, bg, alat
    use spin_orb, only: lspinorb, domag
    use symm_base, only: nsym, s, ft, t_rev, ft
    use fft_base, only: dfftp
    use gvecw, only: ecutwfc
    use gvect, only: ecutrho, ngm
    use spin_orb, only: lspinorb
    use ions_base, only: nat, ityp, nsp, tau, atm, amass
    use io_files, only: psfile
    use wvfct, only: npwx
    use noncollin_module, only: noncolin
    use lsda_mod, only: lsda

    implicit none

    integer, external :: find_free_unit

    integer ik, i, j, isym
    character(len=256) :: datafile
    integer :: io_unit, ios

    ! JLB new variables to support collinear spin-polarized calculation
    LOGICAL :: lspin, lmag
    integer :: nkstot_intw, nbnd_intw


    ! Write crystal.dat file

    io_unit = find_free_unit()
    datafile = trim(intwdir)//"crystal.dat"
    open(unit=io_unit, file=datafile, status="unknown", action="write", form="formatted", iostat=ios)
    if (ios /= 0) call errore( "pw2intw", "write_crystal_info: error opening crystal.dat", ios )

    ! JLB: adapt spin-polarized calculation from QE to spinor in INTW
    ! (see also subroutine write_wfc() below)
    ! Here I define new variables to leave the ones read from QE modules unchanged
    lspin = lsda .or. noncolin
    if (lsda) then
      lmag = .true.
      nkstot_intw = nkstot / 2
      nbnd_intw = nbnd * 2
    else
      ! Just leave as in QE
      lmag = domag
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

    write(unit=io_unit,fmt=*)"NTYP"
    write(unit=io_unit,fmt=*)nsp

    write(unit=io_unit,fmt=*)"ATOM_LABELS, MASS AND PP_FILE (1:NTYP)"
    do i=1,nsp
      write(unit=io_unit,fmt="(a,f16.5,x,a)")atm(i), amass(i), trim(psfile(i))
    enddo

    write(unit=io_unit,fmt=*)"NAT"
    write(unit=io_unit,fmt=*)nat

    write(unit=io_unit,fmt=*)"POSITIONS (1:NAT)"
    do i=1,nat
      write(unit=io_unit,fmt="(a,i4,3f16.8)")atm(ityp(i)),ityp(i), tau(:,i)
    end do

    write(unit=io_unit,fmt=*)"NSYM"
    write(unit=io_unit,fmt=*)nsym

    do isym=1,nsym
      write(unit=io_unit,fmt=*)"SYM"
      write(unit=io_unit,fmt="(i8)")isym
      do i=1, 3
          write(unit=io_unit, fmt="(3i8)") (s(i,j,isym), j=1,3)
      enddo
      write(unit=io_unit,fmt="(100f16.10)")  &
            real(ft(1,isym),dp), real(ft(2,isym),dp), real(ft(3,isym),dp)
      write(unit=io_unit,fmt="(48i3)")t_rev (isym)
    enddo

    write(unit=io_unit,fmt=*)"LSPIN"
    write(unit=io_unit,fmt=*) lspin

    write(unit=io_unit,fmt=*)"LSPINORB"
    write(unit=io_unit,fmt=*)lspinorb

    write(unit=io_unit,fmt=*)"LMAG"
    write(unit=io_unit,fmt=*)lmag

    write(unit=io_unit,fmt=*)"NKS"
    write(unit=io_unit,fmt=*)nkstot_intw !nkstot

    write(unit=io_unit,fmt=*)"NBAND"
    write(unit=io_unit,fmt=*)nbnd_intw !nbnd

    write(unit=io_unit,fmt=*)"FFT GRID"
    write(unit=io_unit,fmt=*)dfftp%nr1,dfftp%nr2, dfftp%nr3

    write(unit=io_unit,fmt=*)"ECUTWFC"
    write(unit=io_unit,fmt=*)ecutwfc

    write(unit=io_unit,fmt=*)"ECUTRHO"
    write(unit=io_unit,fmt=*)ecutrho

    write(unit=io_unit,fmt=*)"NG"
    write(unit=io_unit,fmt=*)ngm

    write(unit=io_unit,fmt=*)"NGK_MAX"
    write(unit=io_unit,fmt=*)npwx

    close(unit=io_unit)

    ! Write kpoints.dat file

    io_unit = find_free_unit()
    datafile = trim(intwdir)//"kpoints.dat"
    open(unit=io_unit, file=datafile, status="unknown", action="write", form="formatted", iostat=ios)
    if (ios /= 0) call errore( "pw2intw", "write_crystal_info: error openeing kpoints.dat", ios )
    !
    do ik=1,nkstot_intw !nkstot
      write(unit=io_unit,fmt="(100f16.10)")( xk(i,ik), i = 1, 3 )
      !write(unit=io_unit,fmt="(100f16.10)") ( et(ibnd,ik)*rytoev, ibnd=1,nbnd )
    enddo
    !
    close(unit=io_unit)

  end subroutine write_crystal_info

end module pw2intw_write