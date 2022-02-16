module intw_pseudo
  !
  use kinds, only: dp
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! variables
  public :: INTWPSEUDO, upf, vlocq
  public :: nqxq, nqx, dq, qrad, tab, spline_ps, tab_d2y, npsx, nh, nhm, &
            nbetam, lmaxkb, lmaxx, nkb, indv, nhtol, nhtolm, ijtoh, vkb, vkqb, &
            dvan, dvan_so, nhtoj, beta
  !
  ! subroutines
  public :: read_all_pseudo
  !
  private
  !
  TYPE INTWPSEUDO
    CHARACTER(LEN=2) :: psd=' '      ! Element label
    CHARACTER(LEN=20) :: typ=' '     ! Pseudo type ( NC or US or PAW)
    CHARACTER(len=6) :: rel=' '      ! relativistic: {no|scalar|full}
    LOGICAL :: nlcc                  ! Non linear core corrections
    REAL(DP) :: zp                   ! z valence
    REAL(DP) :: ecutwfc              ! suggested cut-off for wfc
    !
    INTEGER :: lmax                  ! maximum l component in beta
    INTEGER :: nbeta                 ! number of projectors
    INTEGER, DIMENSION(:), ALLOCATABLE :: kbeta   ! kbeta(nbeta): number of grid points used for betas
                                    ! this defines the cutoff radius for each of them.
    !
    INTEGER :: kkbeta
    INTEGER, DIMENSION(:), ALLOCATABLE :: lll     ! lll(nbeta) l of each projector
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beta   ! beta(mesh,nbeta) projectors
    !
    INTEGER :: mesh                  ! number of points in the radial mesh
    REAL(DP), DIMENSION(:), ALLOCATABLE :: r      ! r(mesh)  radial grid
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rab    ! rab(mesh) dr(x)/dx (x=linear grid)
    INTEGER :: lloc                 ! L of channel used to generate local potential
    REAL(DP) :: rcloc               ! vloc = v_ae for r > rcloc
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vloc    ! vloc(mesh) local atomic potential
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dion  ! dion(nbeta,nbeta) atomic D_{mu,nu}

    LOGICAL :: has_so             ! if .true. includes spin-orbit
    REAL(DP), DIMENSION(:), ALLOCATABLE :: jjj   ! jjj(nbeta) j=l+1/2 or l-1/2 of beta
  END TYPE INTWPSEUDO

  TYPE (INTWPSEUDO), DIMENSION(:), ALLOCATABLE :: UPF

  REAL (DP), ALLOCATABLE :: vlocq(:,:)

  LOGICAL :: spline_ps=.true.

  !Former US in QE
  INTEGER :: &
       nqxq,             &! size of interpolation table
       nqx                ! number of interpolation points
  REAL(DP), PARAMETER:: &
       dq = 0.01D0           ! space between points in the pseudopotential tab.
  REAL(DP), ALLOCATABLE :: &
       qrad(:,:,:,:),         &! radial FT of Q functions
       tab(:,:,:)           ! interpolation table for PPs

  REAL(DP), ALLOCATABLE :: &
       tab_d2y(:,:,:)            ! for cubic splines
  !Former USPP_PARAM in QE (INTW)
  !TYPE (INTWPSEUDO),  ALLOCATABLE, TARGET :: upf(:)

  INTEGER, PARAMETER :: npsx=10

  INTEGER :: &
       nh(npsx),             &! number of beta functions per atomic type
       nhm,                  &! max number of different beta functions per atom
       nbetam                 ! max number of beta functions
  INTEGER :: &
       lmaxkb                ! max angular momentum
  !
  !Former USPP
  INTEGER, PARAMETER :: &
       lmaxx  = 3      ! max non local angular momentum (l=0 to lmaxx)
  !
  INTEGER :: nkb        ! total number of beta functions, with struct.fact.
  !
  INTEGER, ALLOCATABLE ::&
       indv(:,:),        &! indes linking  atomic beta's to beta's in the solid
       nhtol(:,:),       &! correspondence n <-> angular momentum l
       nhtolm(:,:),      &! correspondence n <-> combined lm index for (l,m)
       ijtoh(:,:,:)       ! correspondence beta indexes ih,jh -> composite index ijh
  !
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       vkb(:,:), vkqb(:,:)             ! all beta functions in reciprocal space
  REAL(DP), ALLOCATABLE :: &
       becsum(:,:,:)           ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  REAL(DP), ALLOCATABLE :: &
       dvan(:,:,:),           &! the D functions of the solid
       nhtoj(:,:)              ! correspondence n <-> total angular momentum
  !
  COMPLEX(DP), ALLOCATABLE :: & ! variables for spin-orbit/noncolinear case:
       dvan_so(:,:,:,:),      & ! D_{nm}
       DKB(:,:,:,:)             !
  !
  REAL(DP), ALLOCATABLE :: &
       beta(:,:,:)           ! beta functions for CP (without struct.factor)

contains

  !---------------------------------------------------------------------
  subroutine read_all_pseudo ()
    !
    use intw_utility, only: find_free_unit
    USE intw_reading, only: ntyp
    use intw_input_parameters, only: mesh_dir, prefix
    !
    IMPLICIT NONE
    !
    !     Local variables
    !
    INTEGER :: ios, nr, is, nb, ir, nb1
    INTEGER :: io_unit, ierr
    CHARACTER(256) :: file_pseudo
    CHARACTER(256) :: dum
    CHARACTER(1)   :: tag1
    CHARACTER(2)   :: tag2


    ALLOCATE( UPF(ntyp) )
    !
    ierr = 0
    do is=1,ntyp

      io_unit = find_free_unit()

      if ((is>0).and.(is<9)) then
        write(tag1,"(i1)")is
        write(*,*)"here ...", is, tag1
        file_pseudo=trim(trim(adjustl(mesh_dir)))//trim(prefix)//".save.intw/"//tag1//"-KBPP.txt"
      else if ((is>9).and.(is<19) ) then
        write(tag2,"(i2)")is
        file_pseudo=trim(trim(adjustl(mesh_dir)))//trim(prefix)//".save.intw/"//tag2//"-KBPP.txt"
      else
        print*, "ERROR: The num. of species is bigger than 19 (or <0)"
      end if

      OPEN(UNIT=io_unit,FILE=file_pseudo,STATUS='old',FORM='formatted', IOSTAT=ios)
      write(*,*)"Opening PP of specie num ", is, "ios ", ios

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%psd

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%rel

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%has_so

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%nlcc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%zp

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%ecutwfc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%lloc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%lmax

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%nbeta

      allocate( UPF(is)%kbeta( UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%kbeta
      UPF(is)%kkbeta = MAXVAL( UPF(is)%kbeta(:) )

      allocate( UPF(is)%lll( UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%lll

      allocate( UPF(is)%jjj( UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      if (UPF(is)%has_so) then
        read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%jjj
      else
        read(unit=io_unit,fmt=*,iostat=ierr) dum
      end if

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%mesh

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%rcloc

      allocate( UPF(is)%dion( UPF(is)%nbeta,UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      do nb=1,UPF(is)%nbeta
        read(unit=io_unit,fmt=*,iostat=ierr) (UPF(is)%dion(nb,nb1),nb1=1,UPF(is)%nbeta)
      end do

      nr=UPF(is)%mesh
      nb=UPF(is)%nbeta
      allocate( UPF(is)%r(nr), UPF(is)%rab(nr), UPF(is)%vloc(nr) )
      allocate( UPF(is)%beta(nr,nb) )

      read(unit=io_unit,fmt="(a)",iostat=ierr)dum
      do ir=1,nr
        read(unit=io_unit,fmt=*,iostat=ierr) &
        UPF(is)%r(ir), UPF(is)%rab(ir), UPF(is)%vloc(ir), (UPF(is)%beta(ir,nb), nb=1,UPF(is)%nbeta)
      end do

      IF (ierr /= 0) then
         write(unit=*,fmt=*)"ERROR reading PP, ", file_pseudo
         stop
      END IF
      !
      CLOSE(io_unit)

    END DO !is
    !
  END subroutine read_all_pseudo

end module intw_pseudo
