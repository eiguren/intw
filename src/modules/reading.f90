!----------------------------------------------------------------------------!
!     intw project.
!
!----------------------------------------------------------------------------!
!
module intw_reading
  !
  !----------------------------------------------------------------------------!
  !       The subroutines in this module read the QE data stored in xml format.
  !
  !       IN THIS MODULE:
  !            read_parameters_data_file_xml
  !          read_kpoints_data_file_xml
  !               get_G_data
  !          get_K_folder_data
  !          write_tag
  !            write_parameters_data_file_xml
  !----------------------------------------------------------------------------!
  !
  use kinds, only: dp
  !
  implicit none
  !
  save
  !
  !
  ! variables
  public :: nG_max, nbands, lspinorb, nspin, npol, nkpoints_QE, kpoints_QE, &
            s, can_use_TR, ftau, nrot, nsym, atom_pfile, atom_labels, &
            at, bg, alat, volume0, nat, ntyp, ityp, tau, amass, nr1, nr2, nr3, &
            ngm, gvec, QE_ngk, QE_igk, lsda, noncolin, spinorb_mag, ecutwfc, ecutrho, &
            tpiba, tpiba2, dual, gcutm, gamma_only, &
            num_bands_intw, num_wann_intw, num_exclude_bands_intw, band_excluded_intw
  !
  ! subroutines
  public :: read_parameters_data_file_xml, &
            read_kpoints_data_file_xml, &
            get_ngm, &
            get_ngmax, &
            get_gvec, &
            get_K_folder_data, get_K_folder_data_with_nG,&
            write_tag, &
            deallocate_reading_variables, &
            set_num_bands, &
            scan_file_to
  !
  private
  !

  !----------------------------------------------------------------------------!
  ! Variables strongly inspired by QE variables
  ! (but not necessarily exactly equivalent!)
  !
  ! NOTE: All variables will be read from .xml files produced by QE.
  !       The order of the atoms in the cell will thus be taken directly
  !     as the order in the .xml file.
  !----------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!
    ! Bands variables
  !!!!!!!!!!!!!!!!!!!

  integer :: nG_max
  ! the maximum number of G vectors for any  k point

  integer :: nbands
  ! the number of bands computed in the QE calculation

  !! JLB: These are not read from QE but I think it's cleanest to put them here.
  !!      They are set in the subroutine "set_num_bands" in this module.
  !!      To be discussed.
  integer :: num_bands_intw
  ! Number of bands forming the subspace from which Wannier functions are to be obtained
  ! Might be different from nbands if exclude_bands is used

  integer :: num_wann_intw
  ! Number of Wannier functions, or equivalently bands/KS states to be interpolated
  ! Might be different from num_bands_intw if disentanglement is used

  integer :: num_exclude_bands_intw
  ! Number of bands to be excluded

  logical, allocatable :: band_excluded_intw(:)
  ! Array determining the bands to be excluded

  logical :: lspinorb
  ! if true, spin-orbit non-collinear calculation.

  integer :: nspin
  ! spinor dimension. If nspin=2 non-collinear spin.
  ! warning: only lspinor is read from file.
  ! if lspinor=T > nspin=2, else nspin=1.

  integer :: npol
  ! how many spin polarizations are there in arrays
  ! that are not wavefunctions.
  ! we have npol = nspin, except for magnon = .true.,
  ! in which case nspin = 2 and npol = 1.

  integer :: nkpoints_QE
  !ASIER
  real(kind=dp), allocatable :: kpoints_QE(:,:)
  ! the number of kpoints in the QE folders

  !!!!!!!!!!!!!!!!!!!!!!
    ! Symmetry variables
  !!!!!!!!!!!!!!!!!!!!!!

  integer, allocatable :: s(:,:,:)
  ! symmetry matrices, in crystal coordinates,
  ! CAREFUL!!! since the symmetry matrices are expressed
  ! in the crystal coordinates, they may not *look* like
  ! rotation matrices: convert to cartesian coordinates
  ! and they'll have a sensible form.

  logical, allocatable :: can_use_TR(:)
  ! This array indicates if TR can (should) be used
  ! with a given point group operation. TR is optional
  ! if the ground state is not magnetically polarized, but
  ! can become mandatory; for a magnetic ground state,
  ! certain rotations can only be permitted in conjunction
  ! with time reversal symmetry (ie: the rotation flips B, and
  ! TR flips it back). This array will take note of
  ! which point group symmetries require TR.

  real(dp), allocatable :: ftau(:,:)
  ! fractional translations, in crystal coordinates.

  integer :: nrot
  ! number of bravais lattice symmetries

  integer :: nsym
  ! number of crystal symmetries

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! system / configuration variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  character(len=80), allocatable :: atom_pfile(:)
  ! atom_pfile( j )  = name of pseudopotential file for
  ! the j-th atomic type (or species)

  character(len=3), allocatable :: atom_labels(:)
  ! atom_label( j )  = name of the j-th atomic type (or species)

  real(dp) :: at(3,3)
  ! The a1 a2 and a3 lattice vectors

  real(dp) :: bg(3,3)
  ! The reciprocal lattice vectors

  real(dp) :: alat
  ! The lattice parameter

  real(dp) :: volume0

  integer :: nat
  ! the number of atoms in the cell

  integer :: ntyp
  ! the number types of atoms in the cell

  integer, allocatable :: ityp(:)
  ! ityp( i ) = type of the i-th atom

  real(dp), allocatable :: tau(:,:)
  ! tau( 1:3, i ) = position of the i-th atom (in cartesian, alat units)

  real(dp), allocatable :: amass(:)
  ! amass(1:nat)  = atomic masses

  integer :: nr1, nr2, nr3
  ! The FFT  grid

  integer :: ngm
  ! The true max num. of g vectors

  integer, allocatable :: gvec(:,:)
  ! The global g vectors (to which indices refer).

  integer, allocatable :: QE_ngk(:)
  ! nG for all ks

  integer, allocatable :: QE_igk(:,:)
  ! iG_list for all ks

  logical :: noncolin
  ! if a noncolinear calculation noncolin=T

  logical :: spinorb_mag
  ! if spinorb_mag=T It is a situation with TR broken and spinor calculation.
  ! if spinorb_mag=F TR sym. is present!

  real(dp) :: ecutwfc
  ! cutoff energy for wfc  (Hartree)

  real(dp) :: ecutrho
  ! cutoff energy for wfc (Hartree)

  real(kind=dp):: tpiba != twopi/alat
  real(kind=dp) :: tpiba2 != tpiba**2

  real(dp) :: dual
  real(dp) :: gcutm

  logical :: gamma_only

  logical :: lsda

contains

  subroutine read_parameters_data_file_xml()
    !------------------------------------------------------------------
    ! This subroutine reads in atomic positions and composition,
    ! as well as symmetry data from the "data-file.xml" file,
    ! using the iotk library. See document QE/iotk/doc/manual.txt
    ! for details regarding the iotk subroutines.
    !------------------------------------------------------------------
    use intw_input_parameters, only: mesh_dir, prefix, TR_symmetry
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: twopi

    implicit none

    integer :: io_unit ! input/output file unit
    integer :: ii ! a looping variable

    character(256) :: datafile ! full path of the data-file.xml file in the .xml file
    ! note: mesh_dir and prefix are defined in the input_parameter module.


    integer :: max_nsym ! maximum possible value of nsym. Relevant
    ! variable for the magnetic case.

    integer :: i_sym

    integer, dimension(48) :: trev

    logical :: datafile_exists ! for checking if files in .save QE folder exist

    integer, allocatable, dimension(:,:,:) :: ss
    real(kind=dp), allocatable, dimension(:,:) :: fftau
    character(len=256) :: dummy

    integer :: i, j
    ! First, read the value of the parameter ngm
    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"gvectors.dat")

    inquire(file=datafile, exist=datafile_exists)

    if (.not.datafile_exists) then
      write(*,20) '*****************************************************'
      write(*,20) '*   ERROR                                           *'
      write(*,20) '*   The file gvectors.dat cannot be found,          *'
      write(*,20) '*   check if the folder .save.intw                  *'
      write(*,20) '*   exists or if the prefix in intw input file      *'
      write(*,20) '*   has been correctly chosen.                      *'
      write(*,20) '*   STOPPING                                        *'
      write(*,20) '*****************************************************'
      stop
    end if

    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="unformatted")
    read(unit=io_unit)ngm
    close(unit=io_unit)

    !read the data related to the crystal structure
    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"crystal.dat")
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="formatted")
    !ALAT
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)alat
    !AT
    read(unit=io_unit,fmt=*)dummy
    do i=1,3
      read(unit=io_unit,fmt=*)(at (i,j),j=1,3)
    enddo
    volume0 = alat**3*abs( at(1,1)*(at(2,2)*at(3,3)-at(2,3)*at(3,2)) + &
                           at(1,2)*(at(2,3)*at(3,1)-at(2,1)*at(3,3)) + &
                           at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)) )

    !BG
    read(unit=io_unit,fmt=*)dummy
    do i=1,3
      read(unit=io_unit,fmt=*)(bg (i,j),j=1,3)
    enddo
    !FFT grid nr1 nr2 nr3
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)nr1, nr2, nr3

    !ECUT_WFCS
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)ecutwfc
    !ECUT_POT
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)ecutrho

    tpiba = twopi/alat
    tpiba2 = tpiba*tpiba
    dual = ecutrho/ecutwfc
    gcutm = dual * ecutwfc / tpiba2

    !LSDA
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)lsda
    !NONCOLIN
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)noncolin
    !LSO
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)lspinorb
    !DOMAG
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)spinorb_mag
    if (noncolin) then
       nspin = 2
    else
      nspin = 1
    end if

    !NAT
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)nat

    !NTYP
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)ntyp
    allocate(atom_labels(ntyp))
    allocate(atom_pfile(ntyp))
    allocate(ityp(nat))
    allocate(tau(3,nat))
    allocate(amass(ntyp))
    !ATOM LABELS and PP files
    read(unit=io_unit,fmt=*)dummy
    do i=1,ntyp
      read(unit=io_unit,fmt=*)atom_labels(i), amass(i), atom_pfile(i)
    end do

    !POSITIONS
    read(unit=io_unit,fmt=*)dummy
    do i=1,nat
      read(unit=io_unit,fmt=*)dummy,ityp(i), tau(:,i)
    end do

    !NSYM
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)nsym
    allocate(s(3,3,nsym))
    allocate(ftau(3,nsym))
    allocate(can_use_TR(nsym))
    !
    do ii=1,nsym
      read(unit=io_unit,fmt=*) i_sym
      do i=1, 3
        read(unit=io_unit, fmt=*) (s(i,j,ii), j=1,3)
      enddo
      read(unit=io_unit, fmt=*) ( ftau(j,ii), j=1,3)
      read(unit=io_unit, fmt=*)  trev(ii)
    end do !ii


    if (nspin == 1 .or. nspin ==2.and. TR_symmetry) then
      ! This is the easy case. Simply read in the symmetries

      do ii=1,nsym
        if (nspin == 1) then
          can_use_TR(ii) = TR_symmetry
        else if (nspin ==2) then

          if ( trev(ii) == 1) then
            can_use_TR(ii) = .true.
          else
            can_use_TR(ii) = .false.
            if (.not. spinorb_mag) can_use_TR(ii) = .true.
          end if
        end if !nspin 1 or 2
      end do
    else if (nspin == 2 .and. .not. TR_symmetry) then
      ! This case is slightly more complicated. QE by default seems to
      ! assume time reversal symmetry is sometimes ok even in the presence
      ! of a B field. This is only true if B is colinear! I don't think it is
      ! true in the non-colinear case. Thus, rotations combined with TR must
      ! be suppressed!

      max_nsym = nsym

      nsym = 0

      do ii=1,max_nsym
        if ( trev(ii) == 0) nsym = nsym + 1
      end do

      ! next, populate the symmetry array
      allocate(ss(3,3,max_nsym))
      allocate(fftau(3,max_nsym))

      can_use_TR(:) = .false.
      ss(:,:,1:max_nsym) = s(:,:,1:max_nsym)
      fftau(1:3,1:max_nsym) = ftau(1:3,1:max_nsym)

      s = 0
      ftau = 0.0_dp

      i_sym = 0
      do ii=1,max_nsym
        if ( trev(ii) == 0) then
          i_sym = i_sym + 1
          s(:,:,i_sym) = ss(:,:,ii)
          ftau(:,i_sym) = fftau(:,ii)
        end if
      end do

    end if !nspin

    !-KONTUZ
    npol=nspin
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)nkpoints_QE

    !GAMMA_ONLY
    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)gamma_only
    if (gamma_only) stop "ERROR: intw does not support gamma_only calculations yet"

    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)nG_max


    read(unit=io_unit,fmt=*)dummy
    read(unit=io_unit,fmt=*)nbands


    close(unit=io_unit)
20  format(A)
    return

  end subroutine read_parameters_data_file_xml


  subroutine read_kpoints_data_file_xml(kpoints_cryst)
    !----------------------------------------------------------------------
    !       This subroutine reads in the kpoints from the "data-file.xml"
    !       file, using the iotk library.
    !
    !       The data read is :
    !               kpoint(3,nk)    : the kpoints
    !
    !       NOTE: The tasks performed by this subroutine are NOT performed by
    !             the previous subroutine, "read_parameters_data_file_xml"
    !             because, for some obscure reason, FORTRAN does not allow an
    !             allocatable array to be allocated in a subroutine.
    !          (I no longer think this is true... but it doesn't matter).
    !
    !      The k-points are in cartesian 2pi/a units: they will be converted
    !       to the more convenient crystal units by using the basis vectors
    !       a1 a2 a3, where ai * bj = 2pi delta_{ij}
    !------------------------------------------------------------------------
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit

    implicit none

    integer :: io_unit
    character(256) :: datafile

    real(dp) :: k(3), kpoints_cryst(3,nkpoints_QE)

    integer :: i, j

    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"kpoints.dat")
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="formatted")

    do i=1,nkpoints_QE
      read(unit=io_unit,fmt=*)k(1:3)
      do j=1,3
        kpoints_cryst(j,i) = sum(at(:,j)*k(:)) ! lehen at/alat eginda zegoen ...
      enddo
    end do

    close(unit=io_unit)

  end subroutine read_kpoints_data_file_xml

  subroutine get_ngm ()

    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit

    implicit none

    integer :: io_unit

    !output:
    ! 'ngm', the real max num. of g vectors in the global list.
    !local :
    character(256) :: datafile

    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"gvectors.dat")
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="unformatted")
    read(unit=io_unit)ngm
    close(unit=io_unit)

  end subroutine get_ngm

  subroutine get_ngmax ()

    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit

    implicit none

    integer :: io_unit

    !output:
    ! 'ngm', the real max num. of g vectors in the global list.
    !local :
    character(256) :: datafile

    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"iGlist.dat")
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="unformatted")
    read(unit=io_unit)nG_max
    close(unit=io_unit)

  end subroutine get_ngmax


  subroutine get_gvec ()

    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit

    implicit none

    integer :: io_unit
    !input
    ! output
    ! gvec  : is the real max num. of g vectors in the global list.

    !local
    character(256) :: datafile
    integer :: ig, ik

    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"gvectors.dat")
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="unformatted")
    read(unit=io_unit)ngm
    !allocate(gvec(1:2,ngm))

    do ig=1,ngm
      read(unit=io_unit)gvec(1:3,ig)
    end do

    close(unit=io_unit)

    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"iGlist.dat")
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="unformatted")
    read(unit=io_unit)nG_max

    allocate(QE_igk(nG_max,nkpoints_QE))
    allocate(QE_ngk(nkpoints_QE))
    QE_igk = 0
    do ik=1,nkpoints_QE
      read(unit=io_unit)QE_ngk(ik)
      read(unit=io_unit)QE_igk(1:QE_ngk(ik),ik)
    end do

    close(unit=io_unit)

  end subroutine get_gvec

  subroutine get_K_folder_data(ik,list_iG,wfc,QE_eig)
    !-------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! For the kpoint labeled by ik, this subroutine reads in all the
    ! wavefunctions for bands 1,.., nbands and stores them in the array
    ! wfc(nG_max_k,nbands). It reads the G vectors index array list_iG,
    ! which refers to the global list of G vectors gvecs, from the folder ./K%ik
    ! using the iotk library. It also reads the QE eigenvalues.
    !
    ! Do not be confused! G means "reciprocal lattice vector".
    !                     K means "point in the 1BZ"
    !
    ! NOTE: In fortran, "the first index varies fastest". I take this to mean
    ! that arrays are stored column-wise, namely, for a matrix
    !
    ! M = [ m_{11}   m_{12}   m_{13}]
    ! [ m_{21}   m_{22}   m_{23}]
    ! [ m_{31}   m_{32}   m_{33}]
    !
    ! m_{21} is closer in memory to m_{11} than m_{12}, and in fact
    ! M is stored as
    ! M ~ [ m_{11} m_{21} m_{31} m_{12} m_{22} m_{32} m_{13} m_{23} m_{33}]
    ! Thus, it makes GOOD SENSE to put the G index first, as this is
    ! the index that will be used to perform inner products.
    !
    ! JLB: Excluded bands are taken care for at this point.
    !      Only relevant "num_bands_intw" wave functions are stored on output.
    !      This subroutine is only called in symmetries.f90,
    !      where only the wfc-s within "num_bands_intw" are rotated.
    !
    !------------------------------------------------------------------------
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: cmplx_0

    implicit none

    !I/O variables

    integer,intent(in) :: ik
    integer,intent(out) :: list_iG(nG_max)
    real(dp),intent(out) :: QE_eig(nbands)
    complex(dp),intent(out) :: wfc(nG_max,num_bands_intw,nspin)

     integer :: nG

    !logical variables

    character(256) :: wfc_file, datafile
    integer :: io_unit,is,i,n_yes,ibnd
    integer :: nexclude
    real(dp), parameter :: ha_to_ev = 27.211383860484784
    complex(dp) :: wfc_all(nG_max,nbands,nspin)
    real(dp) :: QE_eig_all(nbands)



    !band_excluded(:) = .false.
    !!
    !do i=1,num_exclude_bands_intw
    !   !
    !   nexclude = exclude_bands(i)
    !   band_excluded(nexclude) = .true.
    !   !
    !enddo
    !
    ! initialize the arrays to zero (zero will be broadcasted)
    !
    list_iG(:) = 0
    wfc_all(:,:,:) = cmplx_0
    wfc    (:,:,:) = cmplx_0
    !
    write(wfc_file,100) ik
    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//trim(wfc_file))
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="unformatted")
    read(unit=io_unit)nG !

    read(unit=io_unit)list_iG(1:nG)
    read(unit=io_unit)QE_eig_all(1:nbands)

    do ibnd=1,nbands
      read(unit=io_unit) (wfc_all((is-1)*nG+1:is*nG,ibnd,is),is=1,nspin)
    enddo

    n_yes = 0
    !
    do ibnd=1,nbands
      !
      if (band_excluded_intw(ibnd)) cycle
      !
      n_yes=n_yes+1
      wfc(:,n_yes,:) = wfc_all(:,ibnd,:)
      QE_eig(n_yes) = QE_eig_all(ibnd)
      !
      !endif
      !
    enddo
    close(io_unit)

100 format('wfc'I5.5'.dat')!
    return

  end subroutine get_K_folder_data

  subroutine get_K_folder_data_with_nG(ik,list_iG,wfc,QE_eig,nG)
    !-------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! For the kpoint labeled by ik, this subroutine reads in all the
    ! wavefunctions for bands 1,.., nbands and stores them in the array
    ! wfc(nG_max_k,nbands). It reads the G vectors index array list_iG,
    ! which refers to the global list of G vectors gvecs, from the folder ./K%ik
    ! using the iotk library. It also reads the QE eigenvalues.
    !
    ! Do not be confused! G means "reciprocal lattice vector".
    !                     K means "point in the 1BZ"
    !
    ! NOTE: In fortran, "the first index varies fastest". I take this to mean
    ! that arrays are stored column-wise, namely, for a matrix
    !
    ! M = [ m_{11}   m_{12}   m_{13}]
    ! [ m_{21}   m_{22}   m_{23}]
    ! [ m_{31}   m_{32}   m_{33}]
    !
    ! m_{21} is closer in memory to m_{11} than m_{12}, and in fact
    ! M is stored as
    ! M ~ [ m_{11} m_{21} m_{31} m_{12} m_{22} m_{32} m_{13} m_{23} m_{33}]
    ! Thus, it makes GOOD SENSE to put the G index first, as this is
    ! the index that will be used to perform inner products.
    !
    ! JLB 07/2023: Shouldn't output wfc have dimensions num_bands_intw instead of nbands?
    !
    !------------------------------------------------------------------------
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: cmplx_0

    implicit none

    !I/O variables

    integer,intent(in) :: ik
    integer,intent(out) :: list_iG(nG_max)
    !real(dp),intent(out) :: QE_eig(nbands)
    real(dp),intent(out) :: QE_eig(num_bands_intw)
    !complex(dp),intent(out) :: wfc(nG_max,nbands,nspin)
    complex(dp),intent(out) :: wfc(nG_max,num_bands_intw,nspin)
    integer, intent(out) :: nG

    !logical variables

    character(256) :: wfc_file, datafile
    integer :: io_unit,is,i,n_yes,ibnd
    integer :: nexclude
    real(dp), parameter :: ha_to_ev = 27.211383860484784
    complex(dp) :: wfc_all(nG_max*nspin,nbands)
    real(dp) :: QE_eig_all(nbands)
    integer :: iG

    !band_excluded(:) = .false.
    !!
    !do i=1,num_exclude_bands
    !   !
    !   nexclude = exclude_bands(i)
    !   band_excluded(nexclude) = .true.
    !   !
    !enddo
    !
    ! initialize the arrays to zero (zero will be broadcasted)
    !
    list_iG(:) = 0
    wfc_all(:,:  ) = cmplx_0
    wfc    (:,:,:) = cmplx_0
    !
    write(wfc_file,100) ik
    datafile = trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//trim(wfc_file))
    io_unit = find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="read",form="unformatted")
    read(unit=io_unit)nG !
    read(unit=io_unit)list_iG(1:nG)
    read(unit=io_unit)QE_eig_all(1:nbands)

    do ibnd=1,nbands
            read(unit=io_unit) (wfc_all((is-1)*nG+1:is*nG,ibnd),is=1,nspin)
    enddo

    n_yes = 0
    !
    do ibnd=1,nbands
      !
      if (band_excluded_intw(ibnd)) cycle
      !
      n_yes=n_yes+1
      do is=1,nspin
        do iG=1,nG
          wfc(iG,n_yes,is) = wfc_all((is-1)*nG+iG,ibnd)
        end do !iG
      end do!is
      QE_eig(n_yes) = QE_eig_all(ibnd)
      !
    enddo
    close(io_unit)

100 format('wfc'I5.5'.dat')!
    return

  end subroutine get_K_folder_data_with_nG


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

    return

  end subroutine write_tag


  subroutine deallocate_reading_variables()

    deallocate(atom_labels)
    deallocate(atom_pfile)
    deallocate(amass)
    deallocate(ityp)
    deallocate(tau)
    deallocate(s)
    deallocate(ftau)

    ! JLB: Not sure whether this should go somewhere else
    if (allocated(band_excluded_intw)) deallocate(band_excluded_intw)

  end subroutine deallocate_reading_variables


  !----------------------------------------------
  subroutine scan_file_to (nnkp_unit,keyword)
    !----------------------------------------------
    !
    !----------------------------------------------------------------------!
    ! This subroutine reads a file all the way to the line
    !            begin $keyword
    ! This is useful when extracting parameters from the ascii file $seed.nnkp.
    ! The subroutine is heavily inspired by the subroutine
    !     QE/PP/pw2pwannier.f90{scan_file_to}
    !
    ! JLB 07/2023: Moved this subroutine here from intw2wannier.f90
    !-----------------------------------------------------------------------!

      implicit none

      character(len=*)   :: keyword
      character(len=256) :: word1, word2
      logical            :: found, test
      integer            :: ios,   nnkp_unit

      found=.false.
      !
      do
         !
         read(nnkp_unit,*,iostat=ios) word1, word2
         !
         test=(trim(word1).eq.'begin').and.(trim(word2).eq.keyword)
         !
         if (test) exit
         !
         if (ios.ne.0) then
            !
            write (*,*) keyword," data-block missing "
            stop
            !
         endif
         !
      enddo
      !
      return

      end subroutine scan_file_to
    !------------------------------------------


  subroutine set_num_bands()
    !----------------------------------------------------------------------
    !       This subroutine sets the sizes of the different band subspaces.
    !       Very important for all subsequent calculations,
    !       should be called at the beginning of all utilities.
    !
    !       For the moment it detects whether a .nnkp file exists,
    !       and in that case it reads them from there.
    !       If there's no such file, the number of bands is set to nbands
    !       (the number of bands in the QE calculation).
    !       JLB: I think it would be useful to add the option to set them on input.
    !
    !       WARNING: This subroutine must always be called after "read_parameters_data_file_xml"
    !                so that nbands is set.
    !
    !------------------------------------------------------------------------

    use intw_utility, only: find_free_unit
    use intw_input_parameters, only: prefix

    implicit none

    character(len=256) :: nnkp_filename
    logical :: have_nnkp
    integer :: nnkp_unit, i, nn

    ! Check if .nnkp file exists
    nnkp_filename = trim(prefix)//trim('.nnkp')
    inquire (file=nnkp_filename,exist=have_nnkp)

    ! Write info
    if(have_nnkp) then
      write(*,'(A)') '|       - .nnkp file found                          |'
      write(*,'(A)') '|         Setting the number of bands from .nnkp    |'
      write(*,'(A)') '|           ---------------------------------       |'
    else
      write(*,'(A)') '|       - .nnkp file NOT found                      |'
      write(*,'(A)') '|         Setting the number of bands from QE       |'
      write(*,'(A)') '|           ---------------------------------       |'
    end if

    ! Set number of bands

    ! from .nnkp
    if(have_nnkp) then
      !
      nnkp_unit=find_free_unit()
      open(unit=nnkp_unit,file=nnkp_filename,status='old')
      !
      ! Number of wannier functions (after disentanglement)
      ! must be the same as number of projections
      if (noncolin) then
        call scan_file_to (nnkp_unit,'spinor_projections')
        read(nnkp_unit,*) num_wann_intw
      else
        call scan_file_to (nnkp_unit,'projections')
        read(nnkp_unit,*) num_wann_intw
      end if
      !
      ! Excluded bands, if any
      call scan_file_to (nnkp_unit,'exclude_bands ')
      read(nnkp_unit,*) num_exclude_bands_intw
      allocate(band_excluded_intw(nbands))
      band_excluded_intw(:)=.false.
      do i=1, num_exclude_bands_intw
         read(nnkp_unit,*) nn
         band_excluded_intw(nn)=.true.
      enddo
      !
      ! Number of bands (before disentanglement)
      num_bands_intw = nbands - num_exclude_bands_intw
      !
    ! all bands from QE
    else
      !
      allocate(band_excluded_intw(nbands))
      band_excluded_intw(:) = .false.
      num_bands_intw = nbands
      num_wann_intw = nbands
      !
    end if

  end subroutine


end module intw_reading
!
!
!----------------------------------------------------------------------------!
