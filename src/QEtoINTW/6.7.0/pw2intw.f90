!-----------------------------------------------------------------------
PROGRAM pw2intw
  !-----------------------------------------------------------------------
  !
  USE Kinds
  USE mp_global,  ONLY : mp_startup
  USE mp_pools,   ONLY : npool
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE mp_images,  ONLY : intra_image_comm

  USE cell_base,  ONLY : at, bg
  USE lsda_mod,   ONLY : nspin, isk
  USE klist,      ONLY : nkstot
  USE io_files,   ONLY : prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE control_flags,    ONLY : gamma_only
  USE environment,ONLY : environment_start, environment_end
!  USE wannier
  USE basis,                ONLY : natomwfc, starting_wfc
  USE uspp_param, ONLY : upf
  USE lsda_mod, ONLY : lsda
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios
  CHARACTER(len=256) :: mesh_dir

  LOGICAL:: phonons=.false.
  CHARACTER(len=256) :: data_dir=""
  CHARACTER(len=256) :: dvscf_dir=""
  CHARACTER(len=256) :: rho_dir=""
  CHARACTER(len=256) :: qlist_file="qlist.txt", ph_dir=""
  integer :: nqirr=0
  logical :: dynxml=.false.
  CHARACTER(len=256) :: intwdir
  INTEGER :: kunittmp, strlen, file_exists

  NAMELIST / inputpp / prefix, mesh_dir, phonons, data_dir, dvscf_dir, rho_dir, qlist_file, nqirr, dynxml, ph_dir

  starting_wfc="file"
  prefix=" "
  mesh_dir="./"


#if defined(MPI)
  CALL mp_startup ( )
#endif

  CALL environment_start ( 'pw2intw' )

  CALL start_clock( 'init_pw2intw' )


  READ (5, inputpp, iostat=ios)

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

  intwdir=trim(trim(mesh_dir)//trim(prefix)//".save.intw")

  call system("mkdir -p "//intwdir)

  WRITE(*,*)
  WRITE(*,*) ' Reading nscf_save data'
  CALL read_file()
  WRITE(*,*)

  call write_pp_intw()

  CALL wfcinit()

  call write_crystal_info_and_bands ()

  call scf_v_and_rho()

  call write_FFT_information ()

  call write_wfc()

  if (phonons) call write_phonon_info()
!
contains



  SUBROUTINE write_pp_intw()
  USE xmltools
  USE ions_base,        ONLY : nsp
  integer :: is,io_unit,ir,nb, nb1
  integer, external :: find_free_unit
  character(len=1) :: tag1
  character(len=2) :: tag2

  do is=1,nsp
    io_unit= find_free_unit()
    if (is<9) then
     write(tag1,"(i1)")is
     open(unit=io_unit,file=trim(trim(mesh_dir))//trim(prefix)//".save.intw/"//tag1//"-KBPP.txt",status="unknown")
    else if ((is>9).and.(is<19) ) then
     write(tag2,"(i2)")is
     open(unit=io_unit,file=trim(trim(mesh_dir))//trim(prefix)//".save.intw/"//tag2//"-KBPP.txt",status="unknown")
    else
     print*, "ERROR: The num. of species is bigger than 19"
    end if

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

    write(unit=io_unit,fmt="(a)")"ECUT WFC"
    write(unit=io_unit,fmt=*) upf(is)%ecutwfc

    write(unit=io_unit,fmt="(a)")"L LOC"
    write(unit=io_unit,fmt=*) upf(is)%lloc

    write(unit=io_unit,fmt="(a)")"LMAX"
    write(unit=io_unit,fmt=*) upf(is)%lmax

    write(unit=io_unit,fmt="(a)")"NBETA"
    write(unit=io_unit,fmt=*) upf(is)%nbeta

    write(unit=io_unit,fmt="(a)")"NUMBER of GRID POINTS FOR EACH BETA"
    write(unit=io_unit,fmt=*) upf(is)%kbeta

    write(unit=io_unit,fmt="(a)")"L FOR EACH BETA"
    write(unit=io_unit,fmt="(100i4)") upf(is)%lll

    write(unit=io_unit,fmt="(a)")"J FOR EACH BETA (IF SO COUPLING)"
    if (upf(is)%has_so) then
      write(unit=io_unit,fmt="(100f6.2)") upf(is)%jjj
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

  USE lsda_mod,   ONLY : nspin
  USE scf, ONLY : rho,v
  implicit none
  integer ::io_unit, ispin
  integer, external :: find_free_unit

  character(len=256) :: datafile
  integer :: record_length

  io_unit=find_free_unit()
  datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"scf_vr.dat")
  inquire(iolength=record_length) v%of_r(:,1)
  open(unit=io_unit,file=datafile,status="unknown", action="write",form="unformatted",access='direct',recl=record_length)

  do ispin=1, nspin
   write (unit=io_unit,rec=ispin) v%of_r(:,ispin)
  end do
  close(unit=io_unit)

  io_unit=find_free_unit()
  datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"scf_rhor.dat")
  inquire(iolength=record_length) rho%of_r(:,1)
  open(unit=io_unit,file=datafile,status="unknown", action="write",form="unformatted",access='direct',recl=record_length)

  do ispin=1, nspin
   write (unit=io_unit,rec=ispin) rho%of_r(:,ispin)
  end do
  close(unit=io_unit)

  END SUBROUTINE scf_v_and_rho


!
  SUBROUTINE write_phonon_info()
  USE xmltools
         USE ions_base,        ONLY : nat
  logical :: existitu
  integer :: iq
  character(len=256) :: numq
  character(len=1) :: nn
  character(len=4) :: num
!
  integer :: io_unit, nirr, ipert, imode0, irr, imode, jmode
  integer, dimension(48) :: npert=-1
  complex(kind=dp) :: u_irr(1:3*nat,1:3*nat,nqirr)
  integer, external :: find_free_unit
  character(len=256) :: datafile, type_string, dv_file, rho_file
!
  inquire(FILE=trim(qlist_file), EXIST=existitu)

!
  do iq=1,nqirr

       write(unit=numq,fmt=*)iq

       datafile=trim(trim(adjustl(ph_dir))//"/"//"qq"//trim(adjustl(numq))//"/_ph0/"&
          //trim(prefix)//".phsave/"//"patterns.1.xml")
       io_unit = xml_openfile( datafile )

       CALL xmlr_opentag( "IRREPS_INFO" )
       CALL xmlr_readtag( "NUMBER_IRR_REP", nirr )
       imode0 = 0
         DO irr = 1, nirr
           CALL xmlr_opentag( "REPRESENTION."//i2c(irr) )
           CALL xmlr_readtag( "NUMBER_OF_PERTURBATIONS", npert(irr) )
             DO ipert = 1, npert(irr)
             imode = imode0 + ipert
             CALL xmlr_opentag( "PERTURBATION."// i2c(ipert) )
             CALL xmlr_readtag( "DISPLACEMENT_PATTERN", u_irr(:,imode,iq) )
             CALL xmlr_closetag( )
             ENDDO
        imode0 = imode0 + npert(irr)
        CALL xmlr_closetag( )

        END DO !irr
!
       CALL xmlr_closetag("IRREPS_INFO")
       CALL xml_closefile ( )
!
    enddo !iq
!
  io_unit= find_free_unit()
  datafile=trim(trim(adjustl(intwdir)))//"/irrq_patterns.dat"

  open(unit=io_unit,file=datafile,status="replace",form="formatted")

  do iq=1, nqirr
    write(unit=io_unit,fmt="(a,i4)")"q", iq
    write(unit=numq,fmt=*)iq
!
    do imode=1,3*nat
      write(unit=io_unit,fmt="(10000(a,f16.10,a,f16.10,a))") ( "(", real(u_irr(jmode,imode,iq),dp), ",", aimag(u_irr(jmode,imode,iq)), ")",jmode=1,3*nat)
    enddo
!
  enddo
  close(unit=io_unit)
!
  do iq=1, nqirr
    write(unit=numq,fmt=*)iq
    dv_file=trim(trim(adjustl(ph_dir)))//"/"//"qq"//trim(adjustl(numq))//"/_ph0/"//trim(prefix)//".dvscf1"
    call system("cp "//dv_file//"   "//&
    trim(adjustl(intwdir))//"/"//trim(prefix)//".dvscf_q"//trim(adjustl(numq)) )

    rho_file=trim(trim(adjustl(ph_dir)))//"/"//"qq"//trim(adjustl(numq))//"/_ph0/"//trim(prefix)//".rho1"
    INQUIRE(FILE=rho_file, EXIST=existitu)
    if (existitu) then
      call system("cp "//rho_file//"   "//&
        trim(adjustl(intwdir))//"/"//trim(prefix)//".rho_q"//trim(adjustl(numq)) )
    end if
  enddo
!
  do iq=1, nqirr
    write(unit=numq,fmt=*)iq
    if (dynxml) then
    !call system("cp "//trim(trim(adjustl(ph_dir)))//"/"//"qq"//trim(adjustl(numq))//"/"//trim(prefix)//".dyn.xml ./"//trim(prefix)//".dyn"//&
    !        trim(adjustl(numq))//".xml")
    else
    call system("cp "//trim(trim(adjustl(ph_dir)))//"/"//"qq"//trim(adjustl(numq))//"/"//trim(prefix)//".dyn "//&
    trim(adjustl(intwdir))//"/"//trim(prefix)//".dyn_q"//trim(adjustl(numq)) )
    end if
  enddo
!
  end SUBROUTINE write_phonon_info

  SUBROUTINE write_fft_information ()
    USE gvect,     ONLY : ngm, g
    USE wvfct,     ONLY : npwx
    USE klist,     ONLY : nks, xk, igk_k, ngk

    IMPLICIT NONE
    INTEGER, EXTERNAL :: find_free_unit

    integer :: ig, io_unit
    CHARACTER(LEN=256) :: datafile
    integer :: ik,i

    io_unit= find_free_unit()
    datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"gvectors.dat")
    open(unit=io_unit,file=datafile,status="unknown", action="write",form="unformatted")

    write(unit=io_unit)ngm
    do ig=1,ngm
       write(unit=io_unit)  (nint(sum(at(:,i)*g(:,ig))), i=1,3)
    end do
    close(unit=io_unit)

    io_unit= find_free_unit()
    datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"iGlist.dat")
    open(unit=io_unit,file=datafile,status="unknown", action="write",form="unformatted")

    write(unit=io_unit)npwx
    do ik=1,nks
       write(unit=io_unit) ngk(ik)
       write(unit=io_unit) igk_k(1:ngk(ik),ik)
    end do
    close(unit=io_unit)

  end SUBROUTINE  write_fft_information

  subroutine write_wfc()
    USE wvfct, ONLY : nbnd, npwx
    USE wavefunctions, ONLY : evc
    USE io_files,        ONLY : nwordwfc, iunwfc
    USE klist, ONLY : nks
    USE klist,           ONLY : nkstot, xk, igk_k, ngk
    USE noncollin_module, ONLY : noncolin
    USE wvfct, ONLY : nbnd, et
    USE constants, ONLY: rytoev
    use intw_utility, only: find_free_unit


    integer :: ik, ibnd, is, npol, ig
    character(256) ::  wfc_file, datafile
    integer :: io_unit

    npol=1
    io_unit = find_free_unit()

    if (noncolin) npol=2

    do ik=1,nks

       write(wfc_file,100) ik
       datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//trim(wfc_file))
       open(unit=io_unit,file=datafile,status="unknown", action="write",form="unformatted")
       CALL davcio (evc, 2*nwordwfc, iunwfc, ik, -1 )
       write(unit=io_unit) ngk(ik)
       write(unit=io_unit) igk_k(1:ngk(ik), ik)
       write(unit=io_unit) ( et(ibnd,ik)*rytoev, ibnd=1,nbnd )

       do ibnd=1, nbnd
          write(unit=io_unit)evc(1:npol*ngk(ik), ibnd)
       enddo
       close(unit=io_unit)

    end do !ik
100 format('wfc'I5.5'.dat')

  end subroutine write_wfc

  SUBROUTINE write_crystal_info_and_bands
    USE wvfct, ONLY : nbnd, et
    USE klist, ONLY : nks
    USE constants, ONLY: rytoev
    USE klist,                ONLY : xk, ngk, nks, nkstot
    USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
    USE spin_orb,             ONLY : lspinorb, domag
    USE symm_base,            ONLY : nrot, nsym, invsym, s, ft, irt, &
         t_rev, sname, time_reversal, no_t_rev, ft
    !USE symme,            ONLY : nrot, nsym, invsym, s, ft, irt, &
    !     t_rev, sname, time_reversal, no_t_rev, ftau
    USE fft_base,         ONLY : dfftp
    USE gvecw,            ONLY : ecutwfc
    USE gvect,     ONLY : ecutrho
    USE spin_orb,             ONLY : lspinorb
    USE ions_base,        ONLY : nat, ityp, nsp, tau, atm, amass, tau_format
    USE io_files,  ONLY : psfile
    USE wvfct,     ONLY : npwx
    USE ener,           ONLY : ef

    IMPLICIT NONE
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    INTEGER ik, ibnd, ibnd1, ikevc, i, j

    REAL(DP), allocatable, dimension(:,:) :: eigval
    INTEGER :: isym
    CHARACTER(LEN=256) :: datafile
    integer :: io_unit

    io_unit= find_free_unit()
    datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"crystal.dat")
    open(unit=io_unit,file=datafile,status="unknown", action="write",form="formatted")

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

    write(unit=io_unit,fmt=*)"NONCOLIN"
    write(unit=io_unit,fmt=*)noncolin
    write(unit=io_unit,fmt=*)"LSPINORB"
    write(unit=io_unit,fmt=*)lspinorb
    write(unit=io_unit,fmt=*)"SPINORB_MAG (KONTUZ, hau aldatu da)"
    write(unit=io_unit,fmt=*)domag
    write(unit=io_unit,fmt=*)"NAT"
    write(unit=io_unit,fmt=*)nat

    write(unit=io_unit,fmt=*)"GAMMA ONLY"
    write(unit=io_unit,fmt=*) gamma_only

    write(unit=io_unit,fmt=*)"EFERMI"
    write(unit=io_unit,fmt=*) ef

    write(unit=io_unit,fmt=*)"LSDA"
    write(unit=io_unit,fmt=*) lsda

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
       write(unit=io_unit,fmt="(i8)"), isym
       do i=1, 3
          write(unit=io_unit, fmt="(3i8)"), (s(i,j,isym), j=1,3)
       enddo
       write(unit=io_unit,fmt="(100f16.10)")  &
            real(ft(1,isym),dp), real(ft(2,isym),dp), real(ft(3,isym),dp)
       write(unit=io_unit,fmt="(48i3)")t_rev (isym)
    enddo

    write(unit=io_unit,fmt=*)"NKS"
    write(unit=io_unit,fmt=*)nkstot

    write(unit=io_unit,fmt=*)"NGMAX"
    write(unit=io_unit,fmt=*)npwx

    write(unit=io_unit,fmt=*)"NBAND"
    write(unit=io_unit,fmt=*)nbnd
    close(unit=io_unit)


    io_unit= find_free_unit()
    datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"kpoints.dat")
    open(unit=io_unit,file=datafile,status="unknown", action="write",form="formatted")

    DO ik=1, nkstot
       write(unit=io_unit,fmt="(100f16.10)")( xk(i,ik), i = 1, 3 )
       !write(unit=io_unit,fmt="(100f16.10)") ( et(ibnd,ik)*rytoev, ibnd=1,nbnd )
    ENDDO
    close(unit=io_unit)

    RETURN
  END SUBROUTINE write_crystal_info_and_bands

    subroutine write_tag(string,i,tag)
    !-----------------------------------------------
    ! This subroutine creates a character string of
    ! the form "string"integer, where the integer
    ! will be immediately after the end of "string",
    ! without blank spaces.
    !-----------------------------------------------
    implicit none

    integer         ::     i
    character(*)    ::     string
    character(256)  ::     integer_part,    tag


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

    tag     =       trim(string)//trim(integer_part)

100 format(I1)
200 format(I2)
300 format(I3)
400 format(I4)
500 format(I5)

    return

  end subroutine write_tag



end PROGRAM pw2intw
