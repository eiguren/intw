!----------------------------------------------------------------------------!
!  intw project.
!
!  This module is strongly inspired from the module "input_parameters"
!  found in the QE distribution.
!
!----------------------------------------------------------------------------!
!
module intw_input_parameters
  !
  !----------------------------------------------------------------------------!
  !
  !  this module contains the definitions of all input parameters for intw
  !  as well as subroutines to read and test the input.
  !
  !----------------------------------------------------------------------------!
  !
  use kinds, only: dp
  !
  implicit none
  !
  save
  !
  ! variables
  public :: mesh_dir, prefix, nk1, nk2, nk3, TR_symmetry, chemical_potential
  public :: intw2W_fullzone, intw2W_method, compute_mmn, compute_amn
  public :: ph_dir, qlist, fc_mat, dvscf_dir, dvscf_name, ep_mat_file, &
            nq1, nq2, nq3, nqirr, calc_epmat

  ! subroutines
  public :: read_input
  !
  private

  !----------------------------------------------------------------------------!
  !  &input namelist nnput parameters
  !----------------------------------------------------------------------------!
  !
  character(len=256) :: mesh_dir = 'unassigned'
  ! the directory where the coarse mesh QE calculations are stored.

  character(len=256) :: prefix = 'unassigned'
  ! the QE prefix for the computations

  integer :: nk1 = 0, nk2 = 0, nk3 = 0
  ! Monkhorst-Pack mesh indices for the coarse mesh

  real(dp) :: chemical_potential = 0.0_dp
  ! The  value which determines the occupation factors, in eV.

  logical :: TR_symmetry
  !If TR symmetry is present TR_symmetry=.true.

  logical :: compute_mmn = .true.
  ! If  True, the code produces the $prefix.mmn  and $prefix.eig files.

  logical :: compute_amn = .true.
  ! If  True, the code produces the $prefix.amn file.

  logical :: intw2W_fullzone = .False.
  ! If  True, the code wil assume that a full zone QE calculation
  ! has been performed and that wavefunctions for every k-point
  ! are available. This is mostly for testing and directly comparing
  ! the results of intw2W90 and pw2wannier.

  character(256) :: intw2W_method = 'CONVOLUTION'
  ! What method should be used to compute matrix elements;
  ! CONVOLUTION or FFT?

  logical :: calc_epmat = .false.
  !Indicates if we want to calculate ep matrix elements

  character(256) :: ph_dir = './'
  character(256) :: dvscf_dir = './'
  character(256) :: dvscf_name = 'dvscf_q'
  character(256) :: qlist = 'qlist.txt'
  character(256) :: fc_mat = '--.fc'

  character(256) :: ep_mat_file = "ep_mat.dat"

  integer :: nq1 = -1, nq2 = -1, nq3 = -1, nqirr = -1


  NAMELIST / input / mesh_dir, prefix, nk1, nk2, nk3, &
                     TR_symmetry, &
                     chemical_potential

  NAMELIST / intw2W / intw2W_fullzone, intw2W_method, compute_mmn, &
                      compute_amn

  NAMELIST / ph / ph_dir, qlist, fc_mat, dvscf_dir, dvscf_name, ep_mat_file, &
                  nq1, nq2, nq3, nqirr, calc_epmat



  ! ----------------------------------------------------------------------
  !  END namelist
  ! ----------------------------------------------------------------------

contains

  subroutine read_input(read_status)
    !------------------------------------------------------------------
    ! This subroutine reads outdir and prefix from standard input
    !------------------------------------------------------------------
    use intw_useful_constants, only: eps_10

    implicit none

    integer :: ios, strlen
    logical :: read_status

    read_status = .false.


    write(*,20) '|       - Reading standard input file...            |'
    write(*,20) '|         namelist             ios                  |'
    write(*,20) '|         --------             ----                 |'


    read( 5, input, iostat = ios)
    write(*,22) '|           &input             ',ios,'                   |'

    read( 5, intw2W, iostat = ios)
    write(*,22) '|           &intw2W            ',ios,'                   |'

    read( 5, ph, iostat = ios)
    write(*,22) '|           &ph                ',ios,'                   |'

    write(*,20) '====================================================='


    !  Test the various read parameters

    if ( mesh_dir .eq. 'unassigned' ) then
      read_status = .true.
      write(*,*) 'MISSING mesh_dir'
    end if

    if ( prefix .eq. 'unassigned' ) then
      read_status = .true.
      write(*,*) 'MISSING prefix'
    end if

    if ( nk1 .eq. 0 ) then
      read_status = .true.
      write(*,*) 'MISSING nk1'
    end if

    if ( nk2 .eq. 0 ) then
      read_status = .true.
      write(*,*) 'MISSING nk2'
    end if

    if ( nk3 .eq. 0 ) then
      read_status = .true.
      write(*,*) 'MISSING nk3'
    end if


    if ( read_status ) then
      write(*,*) "PROBLEM!: the input should be of the form:"
      write(*,*) "&input"
      write(*,*) "             mesh_dir        = 'directory'"
      write(*,*) "             prefix          = 'prefix'"
      write(*,*) "             nk1             = integer "
      write(*,*) "             nk2             = integer "
      write(*,*) "             nk3             = integer "
      write(*,*) "             TR_symmetry     = T or F  "
      write(*,*) "/"
      write(*,*) "&intw2W"
      write(*,*) "             intw2W_fullzone = T or F"
      write(*,*) "             intw2W_method   = CONVOLUTION or FFT"
      write(*,*) "             compute_amn     = T or F"
      write(*,*) "             compute_mmn     = T or F"
      write(*,*) "/"
      write(*,*) "&ph"
      write(*,*) "             ph_dir = 'directory'"
      write(*,*) "             qlist = 'file'"
      write(*,*) "             fc_mat = 'file'"
      write(*,*) "             dvscf_dir = 'directory'"
      write(*,*) "             dvscf_name = 'file'"
      write(*,*) "             ep_mat_file = 'file'"
      write(*,*) "             nq1 = integer"
      write(*,*) "             nq2 = integer"
      write(*,*) "             nq3 = integer"
      write(*,*) "             nqirr = integer"
      write(*,*) "             calc_epmat = T or F"
      write(*,*) "/"
    end if

    strlen = len_trim(mesh_dir)
    if ( mesh_dir(strlen:strlen+1) .ne. "/" ) mesh_dir(strlen+1:strlen+2) = "/"
    strlen = len_trim(ph_dir)
    if ( ph_dir(strlen:strlen+1) .ne. "/" ) ph_dir(strlen+1:strlen+2) = "/"
    strlen = len_trim(dvscf_dir)
    if ( dvscf_dir(strlen:strlen+1) .ne. "/" ) dvscf_dir(strlen+1:strlen+2) = "/"

    return

20 format(A)
22 format(A,I2,A)

  end subroutine read_input


END MODULE intw_input_parameters
