module siesta2intw_pp
  !
  !     Create a the intw pseudopotential file with the KB projectors written in
  !     a X.ion file created by siesta.
  !
  use precision, only: dp

  implicit none

  public :: intwPPs, write_PP_files

  private




  integer, parameter :: max_nkb = 14 ! the maximum number of KBs
  integer, parameter :: max_l = 4 ! maximum orbital angular momentum
  integer, parameter :: siesta_nrpts = 500 ! grid size for radial functions in the ion.xml file
  integer, parameter :: intw_max_nrpts = 3500 ! maximum grid size for radial functions in the intw PP file

  !
  ! ion.xml file variables

  ! This type contains only PP related info (not PAOs)
  type :: ion_xml_pp
    !
    character(len=2) :: symbol ! symbol of the atomic specie
    character(len=20) :: label ! label used by siesta
    real(kind=dp) :: Z ! the ionic charge
    real(kind=dp) :: Zval ! the ionic pseudo charge
    integer :: lmax ! the maximum angular momentum of the KB projectors
    logical :: lj_projs
    ! Header
    character(len=2) :: psf_name ! Atomic specie
    character(len=2) :: psf_icorr ! XC functional
    character(len=3) :: psf_irel ! non-relativisitic or relativistic PP
    character(len=4) :: psf_nicore ! NLCC
    real(kind=dp), dimension(0:max_l) :: psf_ocup ! orbital occupation
    real(kind=dp), dimension(0:max_l) :: psf_rcs ! cutoff radius
    character(len=2), dimension(0:max_l) :: psf_labels ! label of the orbitals
    !
    ! real(kind=dp), dimension(siesta_nrpts) :: vna, r_vna
    ! real(kind=dp), dimension(siesta_nrpts) :: chlocal, r_chlocal
    !
    logical :: nlcc ! true if nlcc pseudopotential
    real(kind=dp), dimension(siesta_nrpts) :: core, r_core ! NLCC charge
    !
    real(kind=dp) :: rcloc ! cut-off radius for local PP
    real(kind=dp), dimension(siesta_nrpts) :: reduced_vlocal, r_reduced_vlocal ! local PP
    !
    integer :: nkb ! number of KB projectors
    integer, dimension(max_nkb) :: ns ! main quantum number of each KB projector
    integer, dimension(max_nkb) :: ls ! orbital angular momentum of each KB projector
    real(kind=dp), dimension(max_nkb) :: js ! total angular momentum of each KB projector (only rel=2)
    real(kind=dp), dimension(max_nkb) :: rcut ! cut-off radius of each KB  projector
    real(kind=dp), dimension(max_nkb) :: ekb ! reference energy of each KB projetor
    real(kind=dp), dimension(siesta_nrpts,max_nkb) :: kb, r_kb ! KB projectors
    !
  end type ion_xml_pp


  type :: intw_pp
    !
    character(len=20) :: label ! label
    real(kind=dp) :: Z ! the ionic charge
    real(kind=dp) :: Zval ! the ionic pseudo charge
    integer :: lmax ! the maximum angular momentum of the KB projectors
    integer :: rel ! 0 nonrelativistic calculation, 1 scalar relativistic calculation, 2 calculation with the full dirac equation
    ! Radial grid
    integer :: mesh ! the actual number of mesh points
    real(kind=dp), dimension(intw_max_nrpts) :: r   ! the radial mesh
    real(kind=dp), dimension(intw_max_nrpts) :: rab ! d r(x) / d x where x is the linear grid
    !
    integer :: lloc ! l component of the local PP
    real(kind=dp) :: rcloc ! cut-off radius for local PP
    real(kind=dp), dimension(intw_max_nrpts) :: vpsloc ! local PP
    !
    logical :: nlcc ! true if nlcc pseudopotential
    real(kind=dp), dimension(intw_max_nrpts) :: rhoc ! NLCC core charge
    !
    integer :: nbeta ! the number of KB projectors
    integer, dimension(max_nkb) :: ls ! orbital angular momentum of each KB projector
    real(kind=dp), dimension(max_nkb) :: js ! total angular momentum of each KB projector (only rel=2)
    integer, dimension(max_nkb) :: ikk ! cutoff radius index for each KB projector
    real(kind=dp), dimension(intw_max_nrpts,max_nkb) :: betas ! KB projectrors
    real(kind=dp), dimension(max_nkb,max_nkb) :: bmat ! DIJ
  end type intw_pp


  type(ion_xml_pp), allocatable, dimension(:) :: siestaPPs

  type(intw_pp), allocatable, dimension(:) :: intwPPs




contains

  subroutine write_PP_files()
    !
    !
    !

    ! variables
    use atm_types, only: nspecies, species
    use siesta2intw_io, only: stdout, outdir, prefix

    implicit none

    integer :: nt
    character(len=256) :: ion_xml_filename, pp_filename
    character(len=2) :: tag


    allocate(siestaPPs(nspecies))
    allocate(intwPPs(nspecies))


    do nt = 1, nspecies

      ! Open read the ion.xml file
      ion_xml_filename = trim(species(nt)%label)//".ion.xml"
      call read_pp_from_ion_xml(ion_xml_filename, siestaPPs(nt))
      !
      ! Check the labels
      if (trim(species(nt)%label) /= trim(siestaPPs(nt)%label)) stop "ERROR: write_PP_files: wrong label"
      !
      ! Transform the PP
      call convert_pp(siestaPPs(nt), intwPPs(nt))
      !
      ! Save PP in intw format
      if (nt<9) then
        write(tag,"(i1)") nt
      else if ((nt>9).and.(nt<19) ) then
        write(tag,"(i2)") nt
      else
        write(stdout,*) "ERROR: write_PP_files: The num. of species is bigger than 19"
      end if
      pp_filename = trim(outdir)//trim(prefix)//".save.intw/"//trim(tag)//"-KBPP.txt"
      !
      call write_pp_intw(intwPPs(nt), pp_filename)

    enddo

  end subroutine write_PP_files


  subroutine read_pp_from_ion_xml(filename, siestaPP)
    !
    ! Read all the needed information from the label.ion.xml file
    !

    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit

    implicit none
    !
    character(len=*), intent(in) :: filename
    type(ion_xml_pp), intent(out) :: siestaPP
    !
    integer :: iounit, iostat
    !
    character(len=256) :: line

    ! Open the ion.xml file
    iounit = find_free_unit()
    open(unit=iounit, file=trim(filename), status="old", action="read", iostat=iostat)
    if (iostat /= 0) stop "ERROR: read_pp_from_ion_xml: ion.xml file could not be opened"
    !
    !
    siestaPP%nlcc = .false.
    siestaPP%lj_projs = .false.
    !
    ! Read the .ion.xml file
    do
      read(iounit,"(a)") line
      if (line(1:8) == "<symbol>")read(line(9:),*) siestaPP%symbol
      if (line(1:7) == "<label>") read(line(8:),*) siestaPP%label
      if (line(1:3) == "<z>") read(line(4:),*) siestaPP%Z
      if (line(1:9) == "<valence>") read(line(10:),*) siestaPP%Zval
      ! if (line(1:6) == "<mass>")
      ! if (line(1:13) == "<self_energy>") read(line(14:),*) etots
      ! if (line(1:12) == "<lmax_basis>")
      ! if (line(1:10) == "<norbs_nl>")
      if (line(1:12) == "<lmax_projs>") read(line(13:),*) siestaPP%lmax
      if (line(1:11) == "<nprojs_nl>") read(line(12:),*) siestaPP%nkb
      if (line(1:10) == "<lj_projs>") read(line(11:),*) siestaPP%lj_projs
      if (trim(line) == "<preamble>") call read_preamble()
      if (trim(line) == "<paos>") call read_paos()
      if (trim(line) == "<kbs>") call read_kbs()
      if (trim(line) == "<vna>") call read_vna()
      if (trim(line) == "<chlocal>") call read_chlocal()
      if (trim(line) == "<reduced_vlocal>") call read_reduced_vlocal()
      if (trim(line) == "<core>") call read_core()
      if (trim(line) == "</ion>") exit
    enddo
    !
    ! A safety check
    if (siestaPP%psf_nicore == 'fcec') stop 'read_pp_from_ion_xml: nicore=fcec case not tested! We must check it.'
    !
    !
    contains
    !
    !
    subroutine read_preamble()
      ! Read the preamble
      implicit none
      integer :: ihead
      character(len=256), dimension(100) :: header

      character(len=256) :: dummyc
      integer :: l
      logical :: data


      data = .false.
      ihead = 0
      do
        read(iounit,"(a)") line
        if (trim(line) == "</pseudopotential_header>") data = .false.
        if (data) ihead = ihead + 1
        if (data) header(ihead) = trim(line)
        if (trim(line) == "<pseudopotential_header>") data = .true.
        if (trim(line) == "</preamble>") exit
      enddo
      !
      ! pp_header
      8005 format(1x,a2,1x,a2,1x,a3,1x,a4)
      8000 format(1x,i2)
      8015 format(1x,2i3,i5,4g20.12)
      8030 format(4(g20.12))
      8040 format(1x,a)
      8091 format(a2,f5.2,a4,f5.2)
      !
      ! Read XC functional, relativistic type and nlcc type from the header
      read(header(1),8005) siestaPP%psf_name, siestaPP%psf_icorr, siestaPP%psf_irel, siestaPP%psf_nicore
      !
      ! Read the labels and ocupations from the header
      line = header(3)
      do l=0,siestaPP%lmax
        read(line(17*l+2:),8091) siestaPP%psf_labels(l), siestaPP%psf_ocup(l), dummyc, siestaPP%psf_rcs(l)
      end do
    end subroutine read_preamble
    !
    !
    subroutine read_paos()
      ! Read the PAOs
      implicit none

      ! logical :: data
      ! integer :: ir, npts
      ! integer :: z ! the index of the orbital for the orbital channel
      ! integer :: ispol
      ! npao = 0
      ! ir = 0
      ! data = .false.
      do
        read(iounit,"(a)") line
        ! if (trim(line) == "<orbital") npao = npao + 1
        ! if (line(1:3) == " z=") read(line(5:29),*) z
        ! if (line(1:7) == " ispol=") read(line(9:33),*) ispol
        ! if (line(1:12) == " ref_energy=") read(line(14:31),*) epseu
        ! if (line(1:6) == "<npts>") read(line(7:),*) npts
        ! if (line(1:8) == "<cutoff>") read(line(9:),*) rcut(npao)
        ! if (npts.ne.siesta_nrpts) stop "STOP: npts .ne. 500"
        ! if (line == "</data>") data = .false.
        ! if (line == "</data>") ir = 0
        ! if (line == "</data>".and.(.not.(z == 1 .and. ispol == 0))) npao = npao - 1
        ! if (data) ir = ir + 1
        ! if (data) read(line,*) r_pao(ir,npao), pao(ir,npao)
        ! if (line == "<data>" .and. z == 1 .and. ispol == 0) data = .true.
        if (trim(line) == "</paos>") exit
      enddo
    end subroutine read_paos
    !
    !
    subroutine read_kbs()
      ! Read the KBs
      implicit none

      integer :: ikb
      integer :: ir, npts
      logical :: data

      npts = 500
      ikb = 0
      ir = 0
      data = .false.
      do
        read(iounit,"(a)") line
        if (trim(line) == "<projector") ikb = ikb + 1
        if (trim(line) == "</projector") exit
        if (line(1:3) == " l=") read(line(5:29),*) siestaPP%ls(ikb)
        if (line(1:3) == " j=") then
          read(line(5:22),*) siestaPP%js(ikb) ! this is present only if lj_projs is true
          if ( siestaPP%psf_irel == "nrl" ) then
            if ( abs( siestaPP%js(ikb) - siestaPP%ls(ikb) - 0.5_dp ) < 0.0001_dp .and. siestaPP%ls(ikb) > 0 ) ikb = ikb - 1
          endif
        endif
        if (line(1:3) == " n=") read(line(5:29),*) siestaPP%ns(ikb) ! this indicates the n within the projectors with same l
        if (line(1:12) == " ref_energy=") read(line(14:31),*) siestaPP%ekb(ikb)
        if (line(1:6) == "<npts>") read(line(7:),*) npts
        if (line(1:8) == "<cutoff>") read(line(9:),*) siestaPP%rcut(ikb)
        if (npts.ne.siesta_nrpts) stop "STOP: npts .ne. 500"
        if (line == "</data>") data = .false.
        if (line == "</data>") ir = 0
        if (data) ir=ir+1
        if (data) read(line,*) siestaPP%r_kb(ir,ikb), siestaPP%kb(ir,ikb)
        if (line == "<data>") data = .true.
        if (trim(line) == "</kbs>") exit
      enddo
      if (ikb/=siestaPP%nkb) stop "ERROR: read_kbs: wrong number of KBs"
    end subroutine read_kbs
    !
    !
    subroutine read_vna()
      ! Read vna
      implicit none

      ! integer :: ir, npts
      ! logical :: data

      ! ir = 0
      ! data = .false.
      do
        read(iounit,"(a)") line
        ! if (line(1:6) == "<npts>") read(line(7:),*) npts
        ! ! if (line(1:8) == "<cutoff>") read(line(9:),*) rc
        ! if (npts.ne.siesta_nrpts) stop "STOP: npts .ne. 500"
        ! if (line == "</data>") data = .false.
        ! if (line == "</data>") ir = 0
        ! if (data) ir=ir+1
        ! if (data) read(line,*) r_vna(ir), vna(ir)
        ! if (line == "<data>") data = .true.
        if (trim(line) == "</vna>") exit
      enddo
    end subroutine read_vna
    !
    !
    subroutine read_chlocal()
      ! Read chlocal
      implicit none

      ! integer :: ir, npts
      ! logical :: data

      ! ir = 0
      ! data = .false.
      do
        read(iounit,"(a)") line
        ! if (line(1:6) == "<npts>") read(line(7:),*) npts
        ! ! if (line(1:8) == "<cutoff>") read(line(9:),*) rc
        ! if (npts.ne.siesta_nrpts) stop "STOP: npts .ne. 500"
        ! if (line == "</data>") data = .false.
        ! if (line == "</data>") ir = 0
        ! if (data) ir=ir+1
        ! if (data) read(line,*) r_chlocal(ir), chlocal(ir)
        ! if (line == "<data>") data = .true.
        if (trim(line) == "</chlocal>") exit
      enddo
    end subroutine read_chlocal
    !
    !
    subroutine read_reduced_vlocal()
      ! Read reduced_vlocal
      implicit none

      integer :: ir, npts
      logical :: data

      npts = 500
      ir = 0
      data = .false.
      do
        read(iounit,"(a)") line
        if (line(1:6) == "<npts>") read(line(7:),*) npts
        if (line(1:8) == "<cutoff>") read(line(9:),*) siestaPP%rcloc
        if (npts.ne.siesta_nrpts) stop "STOP: npts .ne. 500"
        if (line == "</data>") data = .false.
        if (line == "</data>") ir = 0
        if (data) ir=ir+1
        if (data) read(line,*) siestaPP%r_reduced_vlocal(ir), siestaPP%reduced_vlocal(ir)
        if (line == "<data>") data = .true.
        if (trim(line) == "</reduced_vlocal>") exit
      enddo
    end subroutine read_reduced_vlocal
    !
    !
    subroutine read_core()
      ! Read core
      implicit none

      integer :: ir, npts
      logical :: data

      ! Check if the PP has NLCC
      if (siestaPP%psf_nicore == 'nc  ') stop "ERROR: read_core: PP doesn't have NLCC"
      npts = 500
      ir = 0
      data = .false.
      do
        read(iounit,"(a)") line
        if (line(1:6) == "<npts>") read(line(7:),*) npts
        ! if (line(1:8) == "<cutoff>") read(line(9:),*) rc
        if (npts.ne.siesta_nrpts) stop "STOP: npts .ne. 500"
        if (line == "</data>") data = .false.
        if (line == "</data>") ir = 0
        if (data) ir=ir+1
        if (data) read(line,*) siestaPP%r_core(ir), siestaPP%core(ir)
        if (line == "<data>") data = .true.
        if (trim(line) == "</core>") exit
      enddo
    end subroutine read_core
    !
    !
  end subroutine read_pp_from_ion_xml


  subroutine convert_pp(siestaPP, intwPP)
    !
    ! Convert PP to upf format
    !

    ! variables
    use units, only: pi
    use siesta2intw_io, only: stdout

    implicit none

    type(ion_xml_pp), intent(in) :: siestaPP
    type(intw_pp), intent(out) :: intwPP

    real(kind=dp), parameter :: fpi = 4.0_dp*pi

    integer :: ir, ikb

    ! Radial grid variables
    real(kind=dp), parameter :: xmin = -7.0_dp, dx = 0.0125_dp, rmax = 100.0_dp ! default UPF values


    !
    !
    ! pp_info
    intwPP%label = siestaPP%label
    intwPP%Z = siestaPP%Z
    intwPP%Zval = siestaPP%Zval
    intwPP%lmax = siestaPP%lmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set the type of relativistic pseudo-potential
    intwPP%rel = -1
    if ( siestaPP%psf_irel == "nrl" ) intwPP%rel = 0
    if ( siestaPP%psf_irel == "rel" .and. (.not. siestaPP%lj_projs) ) intwPP%rel = 1
    if ( siestaPP%psf_irel == "rel" .and. siestaPP%lj_projs ) intwPP%rel = 2
    ! 11-01-2022:
    ! In QE each projector is duplicated, one for j=l-1/2 and the other for j=l+1/2. Additionally, a new block is added at the end of
    ! the UPF file to specify the l and j of each projector.
    ! In the on-site approximation of SIESTA, even if the calculation includes SO interaction the projetors are the same as if it is
    ! not an SO calculation. The SO components of the Hamiltonian are added separately. In the off-site approach, the projectors are
    ! duplicated as in QE. If the PSF is rel but the calculation is not SO I think that just the Down Pseudopotential is used on the
    ! calculation. Thus, it is clear that:
    !  1 - If the PP is not rel and the calculation is not SO: rel = 0 (no relativistic)
    !  2 - If the PP is rel and the calculation is SO: rel = 2 (full relativistic)
    !  3 - If the PP is rel and the calculation is not SO: rel = 1 (scalar relativistic)
    !      In this case, only the Down Pseudopotential of the PSF file is used, whic is, indeed, the scalar relativistic potential:
    !        Vdown = Vion_l = 1/(2l+1) [ (l+1)V_{l+1/2} + lV_{l-1/2} ]
    !        Vup = Vso_l = 2/(2l+1) [ V_{l+1/2} - V_{l-1/2} ]
    ! But what happens in the other situation:
    !  4 - If the PP is not rel but the calculation is SO?
    !      In this second case, the PSF does not contain any relativistic effect for sure, but the projectors are duplicated. Thus,
    !      the duplicated projetors are exactly equal and just one of them could be used to obtain the exactly same results. Maybe, in
    !      this case, the best would be to remove the duplicated projectors and keep just one, and use rel = 0 (no relativistic).
    !      Anyway, I am not sure if the Dij's or ekb's must be modified. Furthermore, I am not sure if keeping all the
    !      projetors and using rel = 0 (no relativistic) is correct.
    ! 12-01-2022:
    ! I found that in the 4th case I should use rel = 0 and take just one of the duplicated projectors letting the Dij's as they are.
    ! Another option is to use rel = 0 and keep both projectors, but then, for the duplicated projectors Dij's must be divided by 2.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !

    ! Set the XC functional
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! More functionals should be added, and check this ones
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if (siestaPP%psf_icorr == "pb") call set_dft_from_name("PBE")
    ! if (siestaPP%psf_icorr == "ca") call set_dft_from_name("CA")
    !
    !
    ! pp_mesh
    intwPP%mesh = (log(rmax*intwPP%Z)-xmin)/dx+1
    intwPP%mesh = 2*(intwPP%mesh/2) + 1 ! mesh must be odd for simpson integration.
    if(intwPP%mesh+1 > intw_max_nrpts) stop 'ERROR: convert_pp: intw_max_nrpts is too small'
    !
    intwPP%r = 0.0_dp
    intwPP%rab = 0.0_dp
    do ir=1,intwPP%mesh
      intwPP%r(ir) = exp(xmin+dble(ir-1)*dx)/intwPP%Z
      intwPP%rab(ir) = intwPP%r(ir)*dx
    end do
    !
    !
    ! pp_local
    intwPP%lloc = -1 ! Siesta uses its own method to calculate the local part
    intwPP%rcloc = siestaPP%rcloc
    !
    intwPP%vpsloc = 0.0_dp
    call dosplineint( siestaPP%r_reduced_vlocal(1:siesta_nrpts), siestaPP%reduced_vlocal(1:siesta_nrpts), &
                      intwPP%r(1:intwPP%mesh),                   intwPP%vpsloc(1:intwPP%mesh)               )
    do ir=1,intwPP%mesh
      if (intwPP%r(ir) .lt. intwPP%rcloc) then
        intwPP%vpsloc(ir) = ( intwPP%vpsloc(ir) - 2*intwPP%Zval )/intwPP%r(ir)
      else
        intwPP%vpsloc(ir) = - 2*intwPP%Zval/intwPP%r(ir)
      endif
    enddo
    !
    !
    ! pp_nlcc
    intwPP%nlcc = siestaPP%psf_nicore /= 'nc  '
    !
    intwPP%rhoc = 0.0_dp
    if (intwPP%nlcc) then
      call dosplineint( siestaPP%r_core(1:siesta_nrpts), siestaPP%core(1:siesta_nrpts), &
                        intwPP%r(1:intwPP%mesh),           intwPP%rhoc(1:intwPP%mesh)       )
      do ir=1,intwPP%mesh
        if (intwPP%r(ir) .lt. siestaPP%r_core(siesta_nrpts)) then
          intwPP%rhoc(ir) = intwPP%rhoc(ir)*fpi*intwPP%r(ir)**2
        else
          intwPP%rhoc(ir) = 0.0_dp
        endif
      enddo
    end if
    !
    !
    ! pp_nonlocal
    intwPP%nbeta = siestaPP%nkb
    intwPP%ls = siestaPP%ls
    intwPP%js = siestaPP%js
    !
    !
    ! pp_beta
    intwPP%betas = 0.0_dp
    do ikb=1,intwPP%nbeta
      !
      call dosplineint( siestaPP%r_kb(1:siesta_nrpts, ikb), siestaPP%kb(1:siesta_nrpts, ikb),  &
                        intwPP%r(1:intwPP%mesh),              intwPP%betas(1:intwPP%mesh, ikb) )
      do ir=1,intwPP%mesh
        if (intwPP%r(ir) .lt. siestaPP%rcut(ikb)) then
          intwPP%ikk(ikb) = ir ! cutoff_radius_index
          intwPP%betas(ir,ikb) = intwPP%betas(ir,ikb)*(intwPP%r(ir)**(intwPP%ls(ikb)+1))
        else
          intwPP%betas(ir,ikb) = 0.0_dp
        endif
      enddo
      !
    enddo
    !
    !
    ! pp_dij
    intwPP%bmat = 0.0_dp
    do ikb=1,intwPP%nbeta
      intwPP%bmat(ikb,ikb) = siestaPP%ekb(ikb)
    enddo

  end subroutine convert_pp


  subroutine write_pp_intw(intwPP, filename)
    !
    !
    !

    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit

    implicit none

    type(intw_pp), intent(in) :: intwPP
    character(len=*), intent(in) :: filename

    integer :: iounit, iostat
    integer :: ir, ikb, jkb


    iounit = find_free_unit()
    open(unit=iounit, file=filename, status="replace", iostat=iostat)
    if (iostat /= 0) stop "ERROR: write_pp_intw: Error opening KBPP.txt"


    write(unit=iounit,fmt="(a)")"ATOM LABEL"
    write(unit=iounit,fmt=*) trim(intwPP%label)

    write(unit=iounit,fmt="(a)")"IS RELAT."
    write(unit=iounit,fmt=*) intwPP%rel == 1 ! whether the PSF PP is relativistic or not

    write(unit=iounit,fmt="(a)")"HAS SO COUPLING"
    write(unit=iounit,fmt=*) intwPP%rel == 2 ! whether the KB PP has SO porjectors or not

    write(unit=iounit,fmt="(a)")"NLCC"
    write(unit=iounit,fmt=*) intwPP%nlcc

    write(unit=iounit,fmt="(a)")"Z VALENCE"
    write(unit=iounit,fmt=*) intwPP%Zval

    write(unit=iounit,fmt="(a)")"L LOC"
    write(unit=iounit,fmt=*) intwPP%lloc

    write(unit=iounit,fmt="(a)")"L MAX"
    write(unit=iounit,fmt=*) intwPP%lmax

    write(unit=iounit,fmt="(a)")"NBETA"
    write(unit=iounit,fmt=*) intwPP%nbeta

    write(unit=iounit,fmt="(a)")"NUMBER of GRID POINTS FOR EACH BETA"
    write(unit=iounit,fmt=*) (intwPP%ikk(ikb), ikb=1,intwPP%nbeta)

    write(unit=iounit,fmt="(a)")"L FOR EACH BETA"
    write(unit=iounit,fmt="(100i4)") intwPP%ls(1:intwPP%nbeta)

    write(unit=iounit,fmt="(a)")"J FOR EACH BETA (IF SO COUPLING)"
    if (intwPP%rel == 2) then
      write(unit=iounit,fmt="(100f6.2)") intwPP%js(1:intwPP%nbeta)
    else
      write(unit=iounit,fmt=*)"NOTHING"
    end if
    write(unit=iounit,fmt="(a)")"N OF MESH POINTS"
    write(unit=iounit,fmt=*) intwPP%mesh

    write(unit=iounit,fmt="(a)")"RC_LOC"
    write(unit=iounit,fmt=*) intwPP%rcloc

    write(unit=iounit,fmt="(a)")"KB DIJ"
    do ikb=1,intwPP%nbeta
      write(unit=iounit,fmt="(100es18.8)") (intwPP%bmat(ikb,jkb), jkb=1,intwPP%nbeta)
    end do

    write(unit=iounit,fmt="(6a18)") "R", "RAB", "VLOC", "BETA1", "BETA2", "..."
    do ir=1,intwPP%mesh
      write(unit=iounit,fmt="(100es18.8)") &
        intwPP%r(ir), intwPP%rab(ir), intwPP%vpsloc(ir), (intwPP%betas(ir,ikb), ikb=1,intwPP%nbeta)
    end do


  end subroutine write_pp_intw


  SUBROUTINE dosplineint( old_mesh, old_vec, new_mesh, new_vec)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL (DP), INTENT(IN)  :: old_mesh(:), new_mesh(:)
    REAL (DP), INTENT(IN)  :: old_vec(:)
    REAL (DP), INTENT(OUT) :: new_vec(:)
    !
    REAL (DP), ALLOCATABLE :: d2y(:)
    INTEGER                :: i
    INTEGER                :: old_dim, new_dim
    !
    !
    old_dim = SIZE( old_vec )
    new_dim = SIZE( new_vec )
    !
    ! IF ( old_dim /= SIZE( old_mesh ) ) CALL errore( 'dosplineint', 'dimensions of old_mesh and old_vec do not match', 1 )
    IF ( old_dim /= SIZE( old_mesh ) ) stop 'ERROR: dosplineint: dimensions of old_mesh and old_vec do not match'
    !
    ! IF ( new_dim /= SIZE( new_mesh ) ) CALL errore( 'dosplineint', 'dimensions of new_mesh and new_vec do not match', 1 )
    IF ( new_dim /= SIZE( new_mesh ) ) stop 'ERROR: dosplineint: dimensions of new_mesh and new_vec do not match'
    !
    ALLOCATE( d2y( old_dim ) )
    !
    d2y = 0
    !
    CALL spline( old_mesh , old_vec(:), 0.0_DP, 0.0_DP, d2y  )
    !
    DO i = 1, new_dim
        !
        new_vec(i) = splint( old_mesh, old_vec(:), d2y, new_mesh(i) )
        !
    END DO
    !
    DEALLOCATE( d2y )
    !
  END SUBROUTINE dosplineint


  SUBROUTINE spline( xdata, ydata, startu, startd, d2y )
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: xdata(:), ydata(:), startu, startd
    REAL(DP), INTENT(OUT) :: d2y(:)
    !
    INTEGER               :: i, k, ydim
    REAL(DP)              :: p, sig
    REAL(DP), ALLOCATABLE :: u(:)
    !
    !
    ydim = SIZE( ydata )
    !
    ALLOCATE( u( ydim ) )
    !
    u(1)   = startu
    d2y(1) = startd
    !
    DO  i = 2, ydim - 1
        !
        sig    = ( xdata(i) - xdata(i-1) ) / ( xdata(i+1) - xdata(i-1) )
        p      = sig * d2y(i- 1) + 2.0_DP
        d2y(i) = ( sig - 1.0_DP ) / p
        u(i)   = ( 6.0_DP * ( ( ydata(i+1) - ydata(i) ) / &
                  ( xdata(i+1) - xdata(i) ) - ( ydata(i) - ydata(i-1) ) / &
                  ( xdata(i) - xdata(i-1) ) ) / &
                  ( xdata(i+1) - xdata(i-1) ) - sig * u(i-1) ) / p
        !
    END DO
    !
    d2y(ydim) = 0
    !
    DO  k = ydim - 1, 1, -1
        !
        d2y(k) = d2y(k) * d2y(k+1) + u(k)
        !
    END DO
    !
    DEALLOCATE( u )
    !
  END SUBROUTINE spline


  FUNCTION splint( xdata, ydata, d2y, x )
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xdata(:), ydata(:), d2y(:)
    REAL(DP), INTENT(IN) :: x
    !
    REAL(DP) :: splint
    INTEGER  :: khi, klo, xdim
    REAL(DP) :: a, b, h
    !
    !
    xdim = SIZE( xdata )
    !
    klo = 1
    khi = xdim
    !
    klo = MAX( MIN( locate( xdata, x ), ( xdim - 1 ) ), 1 )
    !
    khi = klo + 1
    !
    h = xdata(khi) - xdata(klo)
    !
    a = ( xdata(khi) - x ) / h
    b = ( x - xdata(klo) ) / h
    !
    splint = a * ydata(klo) + b * ydata(khi) + &
              ( ( a**3 - a ) * d2y(klo) + ( b**3 - b ) * d2y(khi) ) * &
              ( h**2 ) / 6.0_DP

  END FUNCTION splint


  FUNCTION locate( xx, x )
    !-------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xx(:)
    REAL(DP), INTENT(IN) :: x
    !
    INTEGER :: locate
    INTEGER :: n, jl, jm, ju
    LOGICAL :: ascnd
    !
    !
    n     = SIZE( xx )
    ascnd = ( xx(n) >= xx(1) )
    jl    = 0
    ju    = n + 1
    !
    main_loop: DO
      !
      IF ( ( ju - jl ) <= 1 ) EXIT main_loop
      !
      jm = ( ju + jl ) / 2
      !
      IF ( ascnd .EQV. ( x >= xx(jm) ) ) THEN
          !
          jl = jm
          !
      ELSE
          !
          ju = jm
          !
      END IF
      !
    END DO main_loop
    !
    IF ( x  ==  xx(1) ) THEN
      !
      locate = 1
      !
    ELSE IF ( x  ==  xx(n) ) THEN
      !
      locate = n - 1
      !
    ELSE
      !
      locate = jl
      !
    END IF
    !
  END FUNCTION locate


end module siesta2intw_pp