! Plot the screened potential of a phonon
program surface_dv

  use kinds, only: dp
  use intw_input_parameters
  use intw_reading
  use intw_ph

  !================================================================================
  !       Declare the variables
  !================================================================================

  implicit none


  integer                  :: q_index_irr, imode, jmode
  integer                  :: i, j, k
  integer                  :: ir, iq
  logical                  :: read_status
  complex(kind=dp), allocatable, dimension(:,:,:) :: pot
  real(kind=dp) :: max_z, dz, surface_z

  CHARACTER(len=256) :: fc_file_name



! fc
  INTEGER ibrav, nq1_fc, nq2_fc, nq3_fc, nat_fc, ntype_fc
  REAL(kind=dp) :: alat_fc
  REAL(kind=dp), dimension(3,3) :: at_fc(3,3)
  complex(kind=dp), allocatable, dimension(:,:,:) :: dyn_q
  real(kind=dp), allocatable, dimension(:,:) :: w2
  real(kind=dp), dimension(3) :: q_point
  complex(kind=dp), allocatable, dimension(:,:,:) :: dv_modes
! fc

  CHARACTER(len=8) :: file_name
  CHARACTER(len=1024) :: files, cmd


  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(*,20) '====================================================='
  write(*,20) '|               program surface_dv                  |'
  write(*,20) '|        ---------------------------------          |'
  write(*,20) '====================================================='
  write(*,20) '|    waiting for input file...                      |'
  !
  !================================================================================
  !       read the input file
  !       Read in the necessary information from standard input
  !================================================================================
  !
  call read_input(read_status)
  !
  if (read_status) then
     !
     stop
     !
  endif
  !
  !================================================================================
  !       read the parameters from the SCF QE calculation, in
  !       particular, read in the symmetry matrices!!!
  !================================================================================
  !
  call read_parameters_data_file_xml()
  !
  allocate( pot(nr1,nr2,nr3) )
  !


  ! DEBUG
  ! print symmetries
  if (.false.) then
  write(*,'(a,i3)') 'haritz: nsym=', nsym
  do i=1,nsym
    do j=1,3
      write(*,'(a,i2,a,3i5,2x,a,f6.3,a)') 's(' , i, ')', s(:,j,i), '(', ftau(j,i), ')'
    enddo
    write(*,20) ''
  enddo
  endif
  ! DEBUG



  ! DEBUG
  ! print the name of the files to be readed
  if (.false.) then
  write(*,20) 'read_parameters_data_file_xml() reads:'
  write(*,20) ''
  write(*,20) '  datafile from:', 'mesh_dir/prefix.save/gvectors', trim(trim(mesh_dir)//trim(prefix)//".save/"//"gvectors.dat")
  write(*,20) ''
  write(*,20) 'datafile from:', 'mesh_dir/prefix.save/data-file.xml', trim(mesh_dir)//trim(prefix)//".save/"//trim("data-file.xml")
  write(*,20) ''
  !
  write(*,20) 'read_ph_information_xml() reads:'
  write(*,20) ''
  write(*,20) '  qlist from:', 'mesh_dir/ph_dir/qlist', trim(trim(mesh_dir)//trim(ph_dir)//trim(qlist))
  write(*,20) ''
  write(*,20) '  datafile (read polarizations) for each qirred from:', 'mesh_dir/ph_dir/data_dir/data_nameNum.xml', &
              trim(mesh_dir)//"/"//trim(ph_dir)//"/"//trim(data_dir)//"/"//trim(trim(data_name)//adjustl('num'))//".xml"
  write(*,20) ''
  !
  write(*,20) 'read_allq_dvr(nqirr,nmodes) reads:'
  write(*,20) ''
  write(*,20) '  dv_name from', 'mesh_dir/ph_dir/dvscf_dir/dvscf_nameNum', trim(trim(mesh_dir)//trim(ph_dir)//trim(dvscf_dir)//trim(dvscf_name)//trim('num'))
  write(*,20) ''
  write(*,20) ''
  write(*,20) ''
  endif
  ! DEBUG



  !
  ! read phonon related data
  !
  call read_ph_information_xml()


  ! DEBUG
  ! print irreducible q points
  if (.false.) then
  write(*,'(a,i3)') 'haritz: nqirr=', nqirr
  do i=1,nqirr
    write(*,'(a,i2,a,3f15.8)') 'q(', i, ')', q_irr(:,i)
  enddo
stop
  endif
  !
  ! print rotated irreducible q points
  if (.false.) then
  write(*,20) 'haritz: rotated q_irr'
  do i=1,nqirr
    do j=1,nsym
      write(*,'(2(a,i2,a),3f10.4)') 'q_irr(', i, ')', 'sym(', j, ')', ( mod(1.0_dp+s(1,k,j)*q_irr(1,i)+s(2,k,j)*q_irr(2,i)+s(3,k,j)*q_irr(3,i),1.0_dp) , k=1,3 )
    enddo
  enddo
  endif
  ! DEBUG


  !
  !================================================================================
  !       Read the force constant matrix from the QE directory
  !================================================================================
  !
  fc_file_name=trim(trim(mesh_dir)//"/"//trim(ph_dir)//"/"//trim(fc_mat))
  !
!  if (.not.lspinorb) then
!     call readfc(fc_file_name,nq1_,nq2_,nq3_,nat,alat,at_frc,ntyp,amass)
!  else
!     call read_fc_from_XML ()
!  endif
  !
  nq1_fc=nq1
  nq2_fc=nq2
  nq3_fc=nq3
  nat_fc=nat
  alat_fc=alat
  at_fc=at
  ntype_fc=ntyp
  !
  call readfc(fc_file_name,nq1_fc,nq2_fc,nq3_fc,nat_fc,alat_fc,at_fc,ntype_fc,amass)
  !
  !================================================================================
  ! Build the phonon qmesh corresponding to the parameters in the input file
  !================================================================================
  !
  nqmesh=nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  !
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !


  ! DEBUG
  !print qmesh
  if (.false.) then
  do i=1,nqmesh
    write(*,'(a,3f10.4)') 'haritz: qmesh(i)', qmesh(:,i)
  enddo
  endif
  ! DEBUG


  !================================================================================
  ! Allocate dynamical matrix and phonon energy array
  !================================================================================
  !
  allocate(dyn_q(3*nat,3*nat,nqirr),w2(3*nat,nqirr))
  !
  dyn_q(:,:,:) = cmplx_0
  w2(:,:) = 0.0_dp
  !
  do iq=3,3
    !
    !-matrize dinamikoa kalkulatu q puntu jakin batentzat.
    !
    q_point = q_irr(:,iq)
    call mat_inv_four_t(q_point,nq1_fc,nq2_fc,nq3_fc,3*nat_fc,frc,dyn_q(:,:,iq))
    !
    !diagonalizatu goian kalkulatutako matrize dinamikoa: polarizazio bektore eta energiak lortu ditugu ia.
    !
    call diagonalize_cmat( 3*nat, dyn_q(:,:,iq), w2(:,iq) )

    ! DEBUG
    ! print frequencies
    if (.false.) print'(a,3f10.5)', 'haritz: q=',q_point(:)
    if (.false.) print*, 'haritz: w=',sqrt(w2(:,iq))
    ! print polarization vectors
    if (.false.) then
    do j=1, nmodes
      dyn_q(:,j,iq) = dyn_q(:,j,iq) / sqrt(sum(dyn_q(:,j,iq)**2.0_dp))
    enddo
    do j=1, nmodes
      do k=1,nat
        write(*,'(a,i,3(2f10.5))') 'haritz:', j, ( dyn_q(3*(k-1)+i,j,iq), i=1, 3)
      enddo
    enddo
    endif
    ! DEBUG
    !
  enddo


  ! DEBUG
  ! print some data
  if (.false.) then
  write(*,20) ''
  write(*,'(3(a,i3,x))') 'nat:', nat, ', nqirr:',nqirr, ', nmodes:', nmodes
  write(*,'(a,3i4)') 'nr1, nr2, nr3:', nr1, nr2, nr3
  write(*,20) ''
  write(*,'(a,f10.6)') 'alat:', alat
  do i=1,3
    write(*,'(a1,i1,a1,3f15.10)') 'a',i, ':', at(:,i)
  enddo
  write(*,20) ''
  do i=1,nat
    write(*,'(a,i2,a,3f15.10)') 'tau(', i, '):', tau(:,i)
  enddo
  write(*,20) ''
  endif
  ! DEBUG



  ! get the size of the unit cell along different axes
!  max_x = sqrt( sum( at(:,1)**2.0 ) ) * alat
!  max_y = sqrt( sum( at(:,2)**2.0 ) ) * alat
!  max_z = sqrt( sum( at(:,3)**2.0 ) ) * alat
  max_z = at(3,3)*alat
  !
  dz = max_z/(nr3-1)
  !
  surface_z = tau(3,nat)*alat
!  surface_z=0.0
  !
  ! read dvscf potential
  !
  nqirr=3 ! read only 3 irreducibles due to the small RAM of my computer
  call read_allq_dvr(nqirr,nmodes)
  !
  ! rotate dv_cart to normal modes
  !
  allocate( dv_modes(nr1*nr2*nr3, nqirr, nmodes) )
  !
  dv_modes(:,:,:) = cmplx_0
  !
  do iq=3,3!1,nqirr
    do imode=1,nmodes
      do jmode=1,nmodes
        dv_modes(:,iq,imode)=dv_modes(:,iq,imode)+(dyn_q(jmode,imode,iq))*dvscf_cart(:,iq,jmode,1,1)
      enddo
    enddo
  enddo
  !
  ! print dvsf for mode imode and qirr q_index_irr
  !
  q_index_irr=3
  !
  files = 'fort.100'
  !
  do imode=1,nmodes
    !
    ! create a list with the files with the potential of each mode
    write(file_name,'(a5,i3)') 'fort.',100+imode
    files = trim(files)//' '//trim(file_name)
    !
    do ir=1,nr1*nr2*nr3
      !
      call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)

      ! DEBUG
      ! save dv_modes along z axis for different x,y values (x=16, y=1:nr2)
      if (.false.) then
      if ( imode==23 .or. imode==24 ) then
        if (i==16)   write(2000+j,*) k, abs(dv_modes(ir,q_index_irr,imode))
        if (i==16 .and. k==nr3 )  write(2000+j,*) ''
      endif
      endif
      ! DEBUG

      pot(i,j,k) = dvscf_cart(ir,q_index_irr,imode,1,1)
!      pot(i,j,k) = dv_modes(ir,q_index_irr,imode)
      !
    enddo
    !


    ! DEBUG
    ! save dv_modes on xy plane for z=100
    if (.false.) then
    do i=1,nr1
      do j=1,nr2
        write(1000+imode,'(2i4,2f18.14)') i, j, pot(i,j,100)
      enddo
    enddo
    endif
    ! DEBUG

    !
    ! save dv_modes along z axis taking the average of the xy plane
    do k=1,nr3
      !
      write(100+imode,"(f18.10)") abs(sum( (pot(:,:,k)) )/nr1/nr2)
      !
    enddo
    !
  enddo
  !
  !
  ! save z coordinates
  do k=1,nr3
    !
    write(100,"(f18.10)") dz*(k-1)-surface_z
    !
  enddo
  !
  !
  cmd = "paste -d ' ' "//trim(files)//' > dv_modes.dat'
  write(*,20) trim(cmd)
  call system(trim(cmd))
  cmd = 'rm '//trim(files)
  write(*,20) trim(cmd)
  call system(trim(cmd))
  !
  ! save data to create a gif
  !
  if (.false.) then
    !
    do i=1000,1200
      !
      do k=1,nr3
        !
        write(i,"(i4,2es18.10)") k, real( sum( pot(:,:,k) )/nr1/nr2 * exp(cmplx_i*i*0.1) )
        !
      enddo
      !
    enddo
    !
    open(unit=11,file="plotme.gp",status="replace",action="write")
    !
    write(11,20) "set term gif animate"
    write(11,20) "set out 'proba.gif'"
    write(11,20) "set xrange [1:384]"
    write(11,20) "set yrange [-3.0:3.0]"
    !
    do i=1000,1200
      !
      write(11,"(a,i4,a)") "p 'fort.",i,"' w l"
      !
    end do
    !
    close(11)
    !
  endif
  !


print*, 'end'

20 format(A)


contains


  subroutine diagonalize_cmat (n,a,w)
    !----------------------------------------
    integer, intent(in)  :: n
    complex(dp),intent(inout) :: a(n,n)
    real(dp),intent(out) :: w(n)

    complex(dp) :: a_pack(n*(n+1)/2)

    integer :: i,j, nfound

    complex(dp) :: cwork(2*n)
    real   (dp) :: rwork(7*n)
    integer     :: iwork(5*n), ifail(n), info


    do j=1,n
       do i=1,j
          a_pack(i+((j-1)*j)/2)=a(i,j)
       enddo
    enddo

    ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
    call ZHPEVX('V', 'A', 'U', n, a_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp,  &
         nfound, w, a, n,  cwork,  rwork,       iwork,          ifail, info)

  end subroutine diagonalize_cmat


end program surface_dv
