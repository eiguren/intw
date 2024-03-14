module intw_pseudo

  use kinds, only: dp

  implicit none

  save

  ! variables
  public :: nt_max, intwpseudo, upf

  ! subroutines
  public :: read_all_pseudo

  private


  integer, parameter :: nt_max = 99 ! Max number of different atomic types
  ! NOTE: Haritz 14/03/2024:
  ! In QE the maximum number of atomic types is 10, however, I didn't find any limit in SIESTA.
  ! In intw we just set a limit of 2 digits due to the PP filename format (see read_all_pseudo below).


  type intwpseudo
    character(len=2)                           :: psd = ' ' ! Element label
    character(len=6)                           :: rel = ' ' ! relativistic: {no|scalar|full}
    logical                                    :: nlcc      ! Non linear core corrections
    real(kind=dp)                              :: zp        ! z valence
    real(kind=dp)                              :: ecutwfc   ! Suggested cut-off for wfc
    logical                                    :: has_so    ! If .true. includes spin-orbit
    !
    integer                                    :: mesh      ! Number of points in the radial mesh
    real(kind=dp), dimension(:),   allocatable :: r         ! r(mesh)  radial grid
    real(kind=dp), dimension(:),   allocatable :: rab       ! rab(mesh) dr(i)/di
    !
    real(kind=dp)                              :: rcloc     ! vloc = v_ae for r > rcloc
    integer                                    :: lloc      ! l of channel used to generate local potential
    real(kind=dp), dimension(:),   allocatable :: vloc      ! vloc(mesh) local atomic potential
    !
    integer                                    :: nbeta     ! Number of beta projectors
    integer                                    :: lmax      ! Max l component in beta
    integer,       dimension(:),   allocatable :: kbeta     ! kbeta(nbeta): number of grid points used for betas
                                                            ! this defines the cutoff radius for each of them.
    integer                                    :: kkbeta    ! Max number of grid points used for betas
    integer,       dimension(:),   allocatable :: lll       ! lll(nbeta) l of each projector
    real(kind=dp), dimension(:),   allocatable :: jjj       ! jjj(nbeta) j=l+1/2 or l-1/2 of beta
    real(kind=dp), dimension(:,:), allocatable :: beta      ! beta(mesh,nbeta) projectors
    real(kind=dp), dimension(:,:), allocatable :: dion      ! dion(nbeta,nbeta) atomic D_{mu,nu}
  end type intwpseudo


  type(intwpseudo), dimension(:), allocatable :: upf



contains

  !---------------------------------------------------------------------
  subroutine read_all_pseudo()
    !
    use intw_utility, only: find_free_unit
    use intw_reading, only: ntyp
    use intw_input_parameters, only: mesh_dir, prefix
    !
    implicit none
    !
    !     Local variables
    !
    integer :: ios, nr, nt, nb, ir, nb1
    integer :: io_unit, ierr
    character(256) :: file_pseudo
    character(256) :: dum
    character(1)   :: tag1
    character(2)   :: tag2


    allocate(upf(ntyp))
    !
    ierr = 0
    do nt = 1, ntyp

      io_unit = find_free_unit()

      if ( nt>=1 .and. nt<=9 ) then
        write(tag1,"(i1)") nt
        write(*,20)"|       - Reading:   "//tag1//"-KBPP.txt"//" ..                  |"
        file_pseudo=trim(mesh_dir)//trim(prefix)//".save.intw/"//tag1//"-KBPP.txt"
      else if ( nt>=10 .and. nt<=nt_max ) then
        write(tag2,"(i2)") nt
        file_pseudo=trim(mesh_dir)//trim(prefix)//".save.intw/"//tag2//"-KBPP.txt"
      else
        stop "ERROR: read_all_pseudo: ntyp > nt_max"
      end if

      open(unit=io_unit ,file=file_pseudo, status='old', form='formatted', iostat=ios)

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%psd
      write(*,20)"|                 .. for the specie "//trim(upf(nt)%psd)//"              |"

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%rel

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%has_so

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%nlcc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%zp

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%ecutwfc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%lloc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%lmax

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%nbeta

      allocate(upf(nt)%kbeta(upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%kbeta
      upf(nt)%kkbeta = maxval( upf(nt)%kbeta(:) )

      allocate(upf(nt)%lll(upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%lll

      allocate(upf(nt)%jjj(upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      if (upf(nt)%has_so) then
        read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%jjj
      else
        read(unit=io_unit,fmt=*,iostat=ierr) dum
      end if

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%mesh

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%rcloc

      allocate(upf(nt)%dion( upf(nt)%nbeta,upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      do nb = 1, upf(nt)%nbeta
        read(unit=io_unit,fmt=*,iostat=ierr) (upf(nt)%dion(nb,nb1),nb1=1,upf(nt)%nbeta)
      end do

      nr=upf(nt)%mesh
      nb=upf(nt)%nbeta
      allocate(upf(nt)%r(nr), upf(nt)%rab(nr), upf(nt)%vloc(nr))
      allocate(upf(nt)%beta(nr,nb))

      read(unit=io_unit,fmt="(a)",iostat=ierr)dum
      do ir = 1, nr
        read(unit=io_unit,fmt=*,iostat=ierr) &
        upf(nt)%r(ir), upf(nt)%rab(ir), upf(nt)%vloc(ir), (upf(nt)%beta(ir,nb), nb=1,upf(nt)%nbeta)
      end do

      if (ierr /= 0) then
         write(unit=*,fmt=*)"ERROR reading PP, ", file_pseudo
         stop
      end if
      !
      close(io_unit)

    end do !nt
    !

20 format(A)
30 format(A,F8.2,6X,A)

  end subroutine read_all_pseudo

end module intw_pseudo
