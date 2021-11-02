module tria_const
  integer, parameter :: i8=selected_int_kind(7)

  PUBLIC :: I8
 

end module tria_const



module intw_tria_diag

use intw_useful_constants

use  tria_const

public :: cross, scalar, norma

integer(KIND=I8), allocatable :: ndif_p_list(:), fsh_index(:,:)
integer(KIND=I8), allocatable :: ndif_list(:)



contains 

subroutine tria_diag(tag_in, Emin, Emax, fsh_sparse, alat, bcell, nbnd, bnd_index, nmax, nfsh_bnd, nfsh_bnd_bare )


  implicit none
  
  character(len=10), intent(in) :: tag_in

  logical :: fsh_sparse
  integer :: nmax
!  integer,intent (in)                      :: nfsh, nbnd
  integer,intent (in)                       ::  nbnd
  integer,dimension(:),intent (inout)    ::  nfsh_bnd
  integer,dimension(:),intent (inout)    ::  nfsh_bnd_bare

  real(kind=dp), intent(in)                :: alat
  real(kind=dp), intent(in),dimension(3,3) :: bcell
   integer, dimension(nbnd), intent(in)        :: bnd_index

  integer(kind=I8) :: n,ntri, ndif_p_tot
  real(kind=dp), allocatable,dimension(:) :: x , y , z, veloc, veloc_d !<-

  integer(kind=I8), allocatable, dimension(:) :: dif_vertex_index !<-
  integer(kind=I8), allocatable, dimension(:,:) :: triangle_node !<- ! n/3 n. vertex. ikusi/konprobatu betetzen dela
  integer(kind=I8), allocatable, dimension(:,:) :: near_vert_tri
  integer(kind=I8), allocatable, dimension(:,:,:) :: near_vert_vert
  integer(kind=I8), allocatable, dimension(:) :: p_index,pi_index 

  real(kind=dp), allocatable, dimension(:,:)  :: near_vert_distance
  real(kind=dp), allocatable, dimension(:,:)  :: near_vert_cot
  real(kind=dp), allocatable, dimension(:  )  :: near_vert_area
  real(kind=dp), allocatable, dimension(:  )  :: near_vert_angle


  integer(kind=I8), allocatable, dimension(:) :: n_of_tri_vert, n_of_vert_vert
  integer(kind=I8)  :: ndif_p
  !local
  integer(kind=I8) :: i, j, k, id, jd, kd, it, jt, kt, kd1, kd2, ii, jj, idd, jdd
  real(kind=dp), allocatable, dimension(:)     :: xd, yd, zd !<-
  real(kind=dp), allocatable, dimension(:,:)     :: vd !<-

  integer(kind=I8) :: ndif
  integer(kind=I8), allocatable, dimension(:) :: aux

  integer(kind=I8), parameter :: nt_max=20
  integer(kind=I8), dimension(3) :: t1, t2

  integer(kind=I8), dimension (nt_max,nt_max) :: mat
  logical :: errepika
  real(kind=dp), dimension(1:3) :: vi, vii, vjj
  real(kind=dp) :: tau1, tau2, alpha, beta, &
       total_area, area1, area2,ctn1,ctn2, cosi
  real(kind=dp) :: ang
  integer(kind=I8) :: summat
  integer  :: ibnd
  !diag
  integer(kind=I8) :: l,inf
  real(kind=dp), allocatable ::  a(:,:),a_bare(:,:),eig(:),eig_bare(:), work(:), b(:,:), b_bare(:,:)
   
  real(kind=dp), parameter :: eps=10e-8

  real(kind=dp), allocatable :: dos_ef (:)

  real(kind=dp) ::  bz_vol, v_mean,eskala
  character(len=8) :: num 
  integer :: ibn

  !-full fermi surface-space variables/parameters
  real(kind=dp), allocatable :: energy_list(:,:),el(:)
  !integer(KIND=I8), allocatable :: ndif_p_list(:), ndif_list(:), fsh_index(:,:)
  integer :: num_fsh, indx(3)

  !asier TEST
   
  real(kind=dp), allocatable :: V1(:), V2(:)
  real(KIND=DP)  :: PROD
  ! FEAST
  integer :: nnz
  integer, allocatable :: isa(:), jsa(:),isb(:), jsb(:), rowind_a(:),colind_a(:), rowind_b(:),colind_b(:)
  double precision , allocatable :: sa(:), sb(:), sb_bare(:), sac(:), sbc(:) 
  double precision  :: row_sum

  integer :: job(8), info
  integer :: ifail,ij,nd
  integer :: nr,ibase1,ibase2,locat
  double precision   :: rfail
  !
  integer,dimension(64) :: feastparam 
  double precision  :: epsout
  integer :: loop
  character(len=1) :: UPLO='F'
!!!!!!!!!!!!!!!!! Matrix declaration variable
  integer :: M0,M
  double precision  :: Emin,Emax
  double precision ,dimension(:,:),allocatable :: Xv, Xv_bare ! eigenvectors
  double precision ,dimension(:),allocatable :: E,E_bare,res ! eigenvalue+residual
  
  integer :: ind(nt_max)
  double precision :: row_values(nt_max),  ind_r (nt_max)
!
  real(dp),allocatable :: eulerc(:)
  real(dp) :: normali, normali_b

  write(*,*)"Allocating, ndif_p_list"
  if (.not.allocated(ndif_p_list)) &
        allocate (ndif_p_list(nbnd))
  write(*,*)"Allocating ndif_list(nbnd)"
  if (.not.allocated(ndif_list)) &
        allocate (ndif_list(nbnd))
  write(*,*)"Allocating dos_ef(nbnd)"
  if (.not.allocated(dos_ef)) &
        allocate (dos_ef(nbnd))
  write(*,*)"Allocating eulerc (nbnd)"
  if (.not.allocated(eulerc )) &
        allocate (eulerc (nbnd)) 

  eulerc(:)=0.0_dp

  dos_ef(:)=0.0_dp
 
  !bcell=transpose(bcell)

  do ibn=1,nbnd 
     ibnd = bnd_INDEX(ibn)

     num = writenum(ibnd)

     write(*,*)"bnd_INDEX(ibn)", ibn, bnd_INDEX(ibn),"num", trim(num)
     open(unit=14,file=TRIM("TRI_"//TRIM(tag_in)//".DAT_"//NUM),status="unknown")
 
     open(unit=15,file=TRIM("A_MAT_"//TRIM(tag_in)//".DAT_"//NUM),status="unknown")

     read(14,*)ntri

     write(*,*)"ntri", ntri
     write(*,*)"Allocating x (ntri), y(ntri) , z(ntri), veloc(Ntri))"
     allocate(x (ntri), y(ntri) , z(ntri), veloc(Ntri)) !n
     write(*,*)"Allocating xd(ntri), yd(ntri), zd(ntri), veloc_d(ntri)"
     allocate(xd(ntri), yd(ntri), zd(ntri), veloc_d(ntri)) !n
     write(*,*)"Allocating triangle_node(3,ntri/3) "
     allocate(triangle_node(3,ntri/3) ) !N
     write(*,*)"Allocating dif_vertex_index (ntri)"
     allocate(dif_vertex_index (ntri)) !n

     k=0

     do j=1,ntri
        k=k+1
        read(14,*) x(k), y(k), z(k), veloc(k)
     enddo

     close(unit=14)

      total_area=0.0
      v_mean = 0.0
      do i=1,ntri/3
         j = (i-1)*3 + 1
         id= (i-1)*3 + 2
         idd =(i-1)*3 + 3
         vi  = (/ x(id ), y(id ), z(id )/) -  (/ x(j), y(j), z(j)/) 
         vi  = matmul(bcell,vi)
         vii = (/ x(idd), y(idd), z(idd)/) -  (/ x(j), y(j), z(j)/) 
         vii  = matmul(bcell,vii)

         v_mean = v_mean + (veloc(i)+ veloc(iD)+veloc(iDD))/3*abs(area_vec(vi,vii))
         total_area = total_area + abs(area_vec(vi,vii))
      enddo

     WRITE(*,*)"Total_area", total_area, "mean v", v_mean/total_area
     eskala=v_mean/total_area

     ndif=1
     dif_vertex_index(1)=1
     xd(1) = x(1) ; yd(1) = y(1) ; zd(1) = z(1) ; veloc_d(1)=veloc(1)

     do i=2,ntri
        do j=1,ndif
           if (abs(xd(j)-x(i))+(abs(yd(j)-y(i))+(abs(zd(j)-z(i))))<10E-7) then
              dif_vertex_index(i) = j
              goto 11
           endif
        enddo
        ndif=ndif+1
        xd(ndif) = x(i) ; yd(ndif) = y(i) ; zd(ndif) = z(i) ; veloc_d(ndif) = veloc(i)
        dif_vertex_index(i) = ndif
11      continue
     enddo

     ndif_list(ibn) = ndif 
     write(*,*)"Allocating p_index(ndif),pi_index(ndif)"
     allocate(p_index(ndif),pi_index(ndif) ) 

     ndif_p=0
     p_index(:) = 0
     do i=1,ndif
        if (( 0 > xd(i) > 1).OR.( 0 > yd(i) > 1).OR.( 0 > zd(i) > 1)) THEN
           write(*,*)"ERROR, x y z beyond the defined range"
           stop 
        else
           if ((abs(xd(i)-1.0_dp)>eps).and.(abs(yd(i)-1.0_dp)>eps).and.(abs(zd(i)-1.0_dp)>eps)) then
              ndif_p=ndif_p+1
              p_index(ndif_p) = i
           endif
        endif

     enddo

     ndif_p_list(ibn) = ndif_p 

     pi_index=0
     do i=1,ndif
        do j=1,ndif_p

           if (((abs(modulo(xd(p_index(j))-xd(i),one)-1.0_dp)<eps).or.(abs(modulo(xd(p_index(j))-xd(i),one))<eps)).and.&
               ((abs(modulo(yd(p_index(j))-yd(i),one)-1.0_dp)<eps).or.(abs(modulo(yd(p_index(j))-yd(i),one))<eps)).and.&
               ((abs(modulo(zd(p_index(j))-zd(i),one)-1.0_dp)<eps).or.(abs(modulo(zd(p_index(j))-zd(i),one))<eps)) )  then

              pi_index (i) = j
           endif
        enddo

        if (pi_index (i)==0) then
           write(*,*)"errorea pi_index" ; stop
        endif

     enddo

     open(unit=123,file=TRIM("TRI_VERT_INDX_"//TRIM(tag_in)//".DAT_"//NUM),status="unknown")
     do i=1,ntri,3
        write(123,*) dif_vertex_index(i),dif_vertex_index(i+1),dif_vertex_index(i+2)
     enddo
     close(unit=123)

     open(unit=123,file=TRIM("VERT_COORD_"//TRIM(tag_in)//".DAT_"//NUM),status="unknown")
     write(123,*)ndif
     do i=1,ndif
        vi=(/xd(i), yd(i),zd(i)/)
        vi=matmul(bcell,vi)
        write(123,"(100F18.12)") vi, veloc_d(i)
     enddo
     close(unit=123)

     write(*,*)"Allocating vd(3,ndif))"
     allocate(vd(3,ndif))

     vd(:,:)=0.0
     vd(1,1:ndif)=xd(1:ndif)
     vd(2,1:ndif)=yd(1:ndif)
     vd(3,1:ndif)=zd(1:ndif)

     if (ntri/3*3/=ntri) then
        write(*,*)"Something strange happened, bukatuko dut.", n, n/3*3
        stop
     endif

     k=0
     do i=1,ntri/3 !triangles
        k=k+1
        triangle_node(1,i) =  dif_vertex_index(k)
        k=k+1
        triangle_node(2,i) =  dif_vertex_index(k)
        k=k+1
        triangle_node(3,i) =  dif_vertex_index(k)
     enddo

     !aurkitu vertize eta triangeluen konexioa. hau, da erpin bakoitza
     !zein triangelukin lotua dagoen

     !do i=1,ntri/3 
     !  write(123,*) pi_index(triangle_node(1,i)), pi_index(triangle_node(2,i)),pi_index(triangle_node(3,i))
     !enddo

     write(*,*)"Allocating near_vert_tri(nt_max,ndif),n_of_tri_vert(ndif)"
     allocate(near_vert_tri(nt_max,ndif),n_of_tri_vert(ndif))
     write(*,*)"Allocating near_vert_vert(nt_max,3,ndif),n_of_vert_vert(ndif)"
     allocate(near_vert_vert(nt_max,3,ndif),n_of_vert_vert(ndif))
     write(*,*)"Allocating ear_vert_distance(nt_max,ndif)"
     allocate(near_vert_distance(nt_max,ndif))
     write(*,*)"Allocating near_vert_cot     (nt_max,ndif)"
     allocate(near_vert_cot     (nt_max,ndif))
     write(*,*)"Allocating near_vert_area    (  ndif)"
     allocate(near_vert_area    (  ndif))
     write(*,*)"Allocating near_vert_angle    (  ndif)"
     allocate(near_vert_angle    (  ndif))
     write(*,*)"Allocating aux(nt_max)"
     allocate(aux(nt_max))

     near_vert_tri(:,:)=0
     n_of_tri_vert(:)=0
     near_vert_vert(:,:,:)=0
     n_of_vert_vert(:)=0

     do i=1,ntri/3 !triangeluetan zehar
        do j=1,3
           jd  = triangle_node(j,i)
           jdd = pi_index(jd) !triangely bakoitzeko erpinak ezberdinen listan        

           n_of_tri_vert(jdd) = n_of_tri_vert(jdd)  + 1 
           if (n_of_tri_vert(jdd)<=nt_max) then
              near_vert_tri(n_of_tri_vert(jdd),jdd) = i ! jd erpin bakoitza zein triangelurekin lotuta 
           else
              write(*,*)"errorea", jd,jdd,n_of_tri_vert(jdd) ; STOP
           end if                                    ! konturatu ez daudela ordenatuta oraindik.
           ! ez daude elkarren
           ! ondoan
        enddo !j
     enddo !i

     !
     id_l: do id=1,ndif_p
        k=1
        j_l: do j=1,n_of_tri_vert(id)
           it=near_vert_tri(j,id)
           i_l: do i=1,n_of_tri_vert(id)
              jt= near_vert_tri(i,id) 
              t1(1:3)= pi_index(triangle_node(1:3, it))
              t2(1:3)= pi_index(triangle_node(1:3, jt))

              mat(i,j) = overlap(t1,t2)!overlap((/pi_index(t1(1)),pi_index(t1(2)),pi_index(t1(3))/),&
              !        (/pi_index(t2(1)),pi_index(t2(2)),pi_index(t2(3))/))


              summat=sum(mat)

              if ( overlap(t1,t2)==2) then

                 call ij_angle_dist(t1,t2,id,jd,kd1,kd2)

                 errepika=.false. 
                 do ii=1,k-1
                    if (near_vert_vert(ii,1,id)==jd) errepika=.true.
                 enddo

                 if ((.not.errepika).and.(k<=n_of_tri_vert(id)+1)) then
                    near_vert_vert(k,1,id) = jd
                    near_vert_vert(k,2,id) = kd1
                    near_vert_vert(k,3,id) = kd2
                    k=k+1

                 else if (k==n_of_tri_vert(id)+1) then
                    !exit i_l

                 end if
              endif


           enddo i_l
        enddo j_l

        !zuloak detektatu -baldin badaude-
        if (sum(mat(1:n_of_tri_vert(id),1:n_of_tri_vert(id)))/=&
               (4+n_of_tri_vert(id))*n_of_tri_vert(id)) then
           write(*,*)"ERROR. The overlap matrix does "
           write(*,*)"not have the correct structure"
           write(*,*)vd(:,p_index(id))
           do i=1,n_of_tri_vert(id)
              write(*,"(100I4)")(mat(i,j),j=1,n_of_tri_vert(id))
           enddo
           WRITE(*,*)


           do i=1,n_of_tri_vert(id)
              write(*,"(100I4)")P_INDEX(near_vert_vert(i,:,id))!(mat(i,j),j=1,n_of_tri_vert(id))
              write(*,"(100F12.6)") vd(1:3,p_index(near_vert_vert(i,1,id))) 
              write(*,"(100F12.6)") vd(1:3,p_index(near_vert_vert(i,2,id)))
              write(*,"(100F12.6)") vd(1:3,p_index(near_vert_vert(i,3,id))) 

           end do
           stop
        else if (minval(near_vert_vert(1:n_of_tri_vert(id),:,id))==0) then
           write(*,*)"errorea, near_vert_vert:", id, near_vert_vert(1:n_of_tri_vert(id),:,id)
           stop
        endif


     enddo id_l

     if (.not.fsh_sparse) then
      write(*,*)"NOT SPARSE. allocating matrix dim ", ndif_P," X ", ndif_P
      write(*,*)"Allocating a(ndif_p,ndif_p),a_bare(ndif_p,ndif_p),b(ndif_p,ndif_p), b_bare(ndif_p,ndif_p) ,&
                    eig(ndif_p),eig_bare(ndif_p),work(ndif_p*(3+ndif_p/2))" 
      allocate(a(ndif_p,ndif_p),a_bare(ndif_p,ndif_p),b(ndif_p,ndif_p), b_bare(ndif_p,ndif_p) ,&
                    eig(ndif_p),eig_bare(ndif_p),work(ndif_p*(3+ndif_p/2)))
      a(:,:)=0.0_dp
     endif

     total_area = 0.0_dp
     near_vert_distance=0.0_dp

     nnz=0 
     do id=1,ndif_p
        PRINT*, "id=", id, n_of_tri_vert(id)
        near_vert_area     (id) = 0.0
        near_vert_angle    (id) = 0.0       

        idd= p_index(id)

        do j=1,n_of_tri_vert(id)
           !IF (ID==500) THEN
           WRITE(*,*) id, idd,p_index(near_vert_vert(j,1,id)) 

           !ENDIF
           vi(1:3) =  vd(1:3,p_index(near_vert_vert(j,1,id))) -  vd(1:3,idd)
           call mod05(vi)
           vi=matmul(bcell,vi)

        end do
        ang=0.0_dp

        row_sum=0.0_dp
        do j=1,n_of_tri_vert(id)   
           !
           vi(1:3) =  vd(1:3,p_index(near_vert_vert(j,1,id))) -  vd(1:3,idd)
           vii(1:3) = vd(1:3,p_index(near_vert_vert(j,2,id))) -  vd(1:3,idd)
           call mod05(vi)
           call mod05(vii)
           call mod05_2(vi,vii)
           vi =matmul(bcell,vi)
           vii=matmul(bcell,vii)
           !

           tau1  = datan2(1.0_dp,ctn(vi,vii))
           !
           vi(1:3) =  vd(1:3,p_index(near_vert_vert(j,1,id))) -  vd(1:3,idd)
           vii(1:3) = vd(1:3,p_index(near_vert_vert(j,3,id))) -  vd(1:3,idd)    
           call mod05(vi)
           call mod05(vii)
           call mod05_2(vi,vii)
           vi =matmul(bcell,vi)
           vii=matmul(bcell,vii)

           tau2 =  datan2(1.0_dp,ctn(vi,vii))
           !          !
           vi(1:3)  = vd(1:3,idd)                     - vd(1:3,p_index(near_vert_vert(j,3,id)))
           vii(1:3) = vd(1:3,p_index(near_vert_vert(j,1,id)))  - vd(1:3,p_index(near_vert_vert(j,3,id)))
           call mod05(vi)
           call mod05(vii)
           call mod05_2(vi,vii)

           vi=matmul(bcell,vi)
           vii=matmul(bcell,vii)

           ctn2 = ctn(vi,vii)
           beta = datan2(1.0_dp,ctn2)
           !          
           !          !
           vi(1:3)  = vd(1:3,idd)                     - vd(1:3,p_index(near_vert_vert(j,2,id)))
           vii(1:3) = vd(1:3,p_index(near_vert_vert(j,1,id)))  - vd(1:3,p_index(near_vert_vert(j,2,id)))
           call mod05(vi)
           call mod05(vii)
           call mod05_2(vi,vii)

           vi=matmul(bcell,vi)
           vii=matmul(bcell,vii)

           !         
           ctn1 = ctn(vi,vii) 
           alpha = datan2(1.0_dp,ctn1)

           !
           !          !-bARIZENTRUAREN AZALERA
           vi (1:3)  = (vd(1:3,p_index(near_vert_vert(j,1,id))) - vd(1:3,idd))                            
           !vii(1:3)  = (vd(1:3,p_index(near_vert_vert(j,2,id))) - vd(1:3,p_index(near_vert_vert(j,1,id))))
           !vii(1:3)  = (vd(1:3,p_index(near_vert_vert(j,3,id))) + vd(1:3,p_index(near_vert_vert(j,1,id))) - 2* vd(1:3,idd) )

           vii(1:3)  = (vd(1:3,p_index(near_vert_vert(j,2,id))) -  vd(1:3,idd) )
           vjj(1:3)  = (vd(1:3,p_index(near_vert_vert(j,1,id))) -  vd(1:3,idd) )

         
           call mod05(vi)
           call mod05(vii)
           call mod05(vjj)
           vii = vii + vjj

           !call mod05_2(vi,vii)

           vi=matmul(bcell,vi)/2.0_DP
           vii=matmul(bcell,vii)/3.0_DP

           area1 = ABS(area_vec(vi,vii))
           !         
           vi (1:3)  = (vd(1:3,p_index(near_vert_vert(j,1,id))) - vd(1:3,idd))                            
           !vii(1:3)  = (vd(1:3,p_index(near_vert_vert(j,3,id))) - vd(1:3,p_index(near_vert_vert(j,1,id))))
           !vii(1:3)  = (vd(1:3,p_index(near_vert_vert(j,3,id))) + vd(1:3,p_index(near_vert_vert(j,1,id))) - 2* vd(1:3,idd) )
           
           vii(1:3)  = (vd(1:3,p_index(near_vert_vert(j,3,id))) -  vd(1:3,idd) )
           vjj(1:3)  = (vd(1:3,p_index(near_vert_vert(j,1,id))) -  vd(1:3,idd) )

           call mod05(vi)
           call mod05(vii)
           call mod05(vjj)
           vii = vii + vjj

           !call mod05_2(vi,vii)

           vi =matmul(bcell,vi )/2.0_DP
           vii=matmul(bcell,vii)/3.0_DP


           area2 = ABS(area_vec(vi,vii))

           near_vert_area     (id)  = near_vert_area    (id)   + (area1 + area2) 
           near_vert_angle    (id)  = near_vert_angle    (id)  + (tau1  + tau2)/(4*pi)


           nnz = nnz+1
           if(.not.fsh_sparse) then
             a(id,near_vert_vert(j,1,id)) = 0.5_dp *( ctn1 + ctn2)
           else       
            row_values(j) = 0.5_dp *( ctn1 + ctn2)
            row_sum = row_sum + 0.5_dp *( ctn1 + ctn2)
            ind_r(j) = real(near_vert_vert(j,1,id),dp)
           endif
        enddo
  
      eulerc(ibn) =  eulerc(ibn) + (1 -  near_vert_angle    (id) )  

           nnz = nnz+1
           if(fsh_sparse) then
            ind_r(n_of_tri_vert(id) +1 ) = real(id,dp)
            vi =matmul(bcell, vd(:,idd)-0.5)
            row_values(n_of_tri_vert(id) +1) = -row_sum !+ (vi(3))*0.01 

           
            ind=0
            call hp_sort_r (n_of_tri_vert(id) +1,  ind_r(1:(n_of_tri_vert(id) +1))  , ind)
            do j=1,n_of_tri_vert(id) +1
              write(unit=15,fmt="(2i8,2x,se18.10)") id, nint(ind_r(j)), row_values(ind(j))
            enddo
           end if  
 

     end do !id MAIN

     
    
     write(*,*)"Number of non zero elements ", nnz
     if(fsh_sparse) then

     endfile (unit=15) !
     rewind (unit=15)  !prestatu irakurtzeko 
  
     write(*,*)"Allocating isa(ndif_p+1)"
     allocate(isa(ndif_p+1) )
     write(*,*)"Allocating jsa(nnz) "
     allocate(jsa(nnz) )
     write(*,*)"Allocating sa  (nnz) "
     allocate(sa  (nnz) )
     write(*,*)"Allocating sac (nnz) "
     allocate(sac (nnz) )
     write(*,*)"Allocating rowind_a(nnz)"
     allocate(rowind_a(nnz) )
     write(*,*)"Allocating colind_a(nnz) )"
     allocate(colind_a(nnz) )

     write(*,*)"Allocating isb(ndif_p+1)"
     allocate(isb(ndif_p+1) )
     write(*,*)"Allocating jsb(ndif_p)"
     allocate(jsb(ndif_p) )
     write(*,*)"Allocating sb (ndif_p) "
     allocate(sb (ndif_p) )
     write(*,*)"Allocating sb_bare (ndif_p)"
     allocate(sb_bare (ndif_p) )
     write(*,*)"Allocating sbc (nnz) "
     allocate(sbc (nnz) )
     write(*,*)"Allocating rowind_b(ndif_p)"
     allocate(rowind_b(ndif_p) )
     write(*,*)"Allocating colind_b(ndif_p)"
     allocate(colind_b(ndif_p) )
     
     
     do i=1,nnz
       read(unit=15,fmt=*) rowind_a(i), colind_a(i), sac(i)
     enddo
 
     do i=1,ndif_p
       rowind_b(i) = i 
       colind_b(i)=  i
       sb(i) = near_vert_area  (i)/ veloc_d(p_index(i)) 
       sb_bare(i) = near_vert_area  (i)
     enddo

  ifail=0
  rfail=0.d0
  info = 0
  job(2)=1
  job(3)=1
  !job(4)=2

  job(5)=nnz
  job(1)=1
  job(6)=0
  call mkl_dcsrcoo (job, ndif_p, sa, jsa,isa,nnz   , sac, rowind_a,colind_a,info)

  job(5)=ndif_p
  job(1)=1
  job(6)=0
  call mkl_dcsrcoo (job, ndif_p, sb, jsb,isb,ndif_p, sb, rowind_b,colind_b,info)

  call mkl_dcsrcoo (job, ndif_p, sb_bare, jsb,isb,ndif_p, sb_bare, rowind_b,colind_b,info)


!. FEAST
  call feastinit(feastparam)
  feastparam(1)=1
  
  write(*,*)"Allocating Xv(ndif_p, nfsh_bnd(ibn))"
  allocate(Xv(ndif_p, nfsh_bnd(ibn)))
  write(*,*)"Allocating Xv_bare(ndif_p, nfsh_bnd(ibn))"
  allocate(Xv_bare(ndif_p, nfsh_bnd(ibn)))
  write(*,*)"Allocating e(nfsh_bnd(ibn)), res(nfsh_bnd(ibn))"
  allocate(e(nfsh_bnd(ibn)), res(nfsh_bnd(ibn)))
  write(*,*)"Allocating e_bare(nfsh_bnd(ibn))"
  allocate(e_bare(nfsh_bnd(ibn)))

  !Emin= -0.01_dp
  !Emax= 1000.0_dp
  M0=nfsh_bnd(ibn)
  sa = -sa
  call dfeast_scsrgv('F',ndif_p, sa,isa,jsa, sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,E,Xv,M,res,info)
  nfsh_bnd(ibn)=M0


  print*,  "eskala", eskala

 M0=nfsh_bnd_bare(ibn)
 call dfeast_scsrgv('F',ndif_p, sa,isa,jsa, sb_bare,isb,jsb,feastparam,epsout,loop,Emin,Emax/eskala,M0,E_bare,Xv_bare,M,res,info)
 nfsh_bnd_bare(ibn)=M0

  else
     b(:,:)=0.0_dp
     b_bare =0.0_dp
     do i=1,ndif_P
        b(i,i) = near_vert_area  (i)/ veloc_d(p_index(i)) ! /total_area
        b_bare(i,i) = near_vert_area  (i)
     enddo

     do i=1,ndif_P
        do j=1,ndif_P
           if (abs(a(i,j)-a(j,i))>10E-2) then
              write(*,*)"kontuz, ez da simetrikoa..", i,j, a(i,j), a(j,i)
           endif
        enddo
     enddo

     a(:,:)=-a(:,:)
     do i=1,ndif_p
        a(i,i) = - sum (a(i,:))
     enddo
     a_bare=a

     write(*,*)"Regular diagonalization ..."

     l=ndif_p*(3+ndif_p/2)
     do i=1,ndif
        id=pi_index(i)
        vi=(/xd(i),yd(i),zd(i)/)
     enddo
     
     CALL dsygv( 1, 'V', 'U', NDIF_P, a, NDIF_P, b, NDIF_P, EIG, work, l, inf )
     CALL dsygv( 1, 'V', 'U', NDIF_P, a_bare, NDIF_P, b_bare, NDIF_P, EIG_bare, work, l, inf )
  
  end if !sparse


     open(unit=100,file=TRIM("FSH_"//TRIM(tag_in)//".DAT_"//NUM),status="unknown")
     open(unit=102,file=TRIM("FSH_ENERGY_"//TRIM(tag_in)//".DAT_"//NUM),status="unknown")
     open(unit=200,file=TRIM("FSH_"//TRIM(tag_in)//"_bare.DAT_"//NUM),status="unknown")
     open(unit=202,file=TRIM("FSH_ENERGY_"//TRIM(tag_in)//"_bare.DAT_"//NUM),status="unknown")

     if (fsh_sparse) then
     do i=1,min(ndif_p,nfsh_bnd(ibn))
        write(102,*)i,e(i)
     enddo
     do i=1,min(ndif_p,nfsh_bnd_bare(ibn))
        write(202,*)i,e_bare(i)
     enddo


     normali  =0.0_dp
     normali_b=0.0_dp

     do i=1, ndif_p
        normali   = normali   +  near_vert_area  ( i )/ veloc_d(p_index(i))
        normali_b = normali_b +  near_vert_area  ( i )
     enddo
     normali  =sqrt(normali  )
     normali_b=sqrt(normali_b)


     do j=1,min(ndif_p,nfsh_bnd(ibn))
         !if (xd(p_index(maxval( Xv(1:ndif_p,j)))).lt.0) Xv(:,j) = -Xv(:,j)
     enddo
 
     do i=1,ndif
        vi=(/xd(i),yd(i),zd(i)/)
        vi=matmul(bcell,vi)
        write(100,"(1000(x,se18.12))")vi, 1/veloc_d(i), near_vert_area  (pi_index(i)), (1-near_vert_angle  (pi_index(i))),&
                     ( normali *Xv(pi_index(i),j), j=1,min(ndif_p,nfsh_bnd(ibn)) )
        write(200,"(1000(x,se18.12))")vi, 1/veloc_d(i), near_vert_area  (pi_index(i)), (1-near_vert_angle  (pi_index(i))),&
                    ( normali_b*Xv_bare(pi_index(i),j), j=1,min(ndif_p,nfsh_bnd_bare(ibn)) )
     enddo

     else ! not sparse
     
     do i=1,min(ndif_p,nfsh_bnd(ibn))
        write(102,*)i,eig(i)
        write(202,*)i,eig_bare(i)
     enddo

     normali  =0.0_dp
     normali_b=0.0_dp
     do i=1, ndif_p
        normali   = normali   +  near_vert_area  ( i )/ veloc_d(p_index(i))
        normali_b = normali_b +  near_vert_area  ( i )
     enddo
     write(2468,*) 'regular', ibn, normali, normali_b
     normali  =sqrt(normali  )
     normali_b=sqrt(normali_b)
     
     do i=1,ndif
        vi=(/xd(i),yd(i),zd(i)/)
        vi=matmul(bcell,vi)
        write(100,"(1000F18.12)")vi, 1/veloc_d(i), near_vert_area  (pi_index(i)), &
                near_vert_angle  (pi_index(i)), (normali*a(pi_index(i),j), j=1,min(ndif_p,nfsh_bnd(ibn)) )
        write(200,"(1000F18.12)")vi, 1/veloc_d(i), near_vert_area  (pi_index(i)), &
                near_vert_angle  (pi_index(i)), (normali_b*a_bare(pi_index(i),j), j=1,min(ndif_p,nfsh_bnd(ibn)) )
     enddo

     endif

     total_area=0.0_dp
     do id=1,ndif_p
        total_area = total_area +  near_vert_area    (id)
     enddo


     print*, "Euler C    ", eulerc(ibn)
     print*, "Total area ", total_area

     print*,"Ndif       ", ndif
     print*,"ndif_p     ", ndif_p

     dos_ef(ibn)=0.0d0
     do i=1,ndif_P
      dos_ef(ibn) = dos_ef(ibn) + near_vert_area    (I)/ veloc_d(p_index(i)) 
     enddo

     print*, "dos_ef_bare", dos_ef(ibn)

     bz_vol = scalar(bcell (:,1), cross( bcell (:,2),bcell (:,3) ) ) 

     dos_ef(ibn)=dos_ef(ibn) *  2 *  (2*PI)**2 / bz_vol / (2*PI)**3 
                             !SPINT  D^2K         ! vbz
     print*, "dos_ef ", DOS_EF(ibn) , " (st/eV/Spin)"

     open(unit=500, file=TRIM("DOS_"//TRIM(tag_in)//".DAT_"//num), status='unknown')
     write(500, "(10F18.12)") DOS_EF(ibn) 
     close (unit=500)

     open(unit=600, file=TRIM("AREA_"//TRIM(tag_in)//".DAT_"//num), status='unknown')
     write(600, "(10F18.12)") total_area
     close (unit=600)

     open(unit=700, file=TRIM("EULER_C_"//TRIM(tag_in)//".DAT_"//num), status='unknown')
     write(700, "(10F18.12)") eulerc(ibn)
     close (unit=700)

     

     close(unit=100)
     close(unit=200)
     close(unit=102)
     close(unit=202)
     close(unit=15)

     write(*,*)"...bukatuta..", num
     write(*,*)"...dos_ef total ", sum(dos_ef), sum(dos_ef)*bz_vol, sum(dos_ef)*bz_vol*pi, sum(dos_ef)*pi

     !nORM TEST
     !ALLOCATE(V1(ndif_p), V2(ndif_p))
     !DO II=1,10
     ! DO JJ=1,10 
     !  V1(:) = a(1:ndif_p,II)
     !  V2(:) = a(1:ndif_p,JJ)
     !  PROD=0.0
     !  DO I=1, ndif_p 
     !  PROD=PROD + near_vert_area  (i) /veloc_d(p_index(i)) * V1(I) * V2(I)
     !  ENDDO 
     !  PRINT*,I,J,"PROD", PROD
     ! ENDDO
     !ENDDO
     !DEALLOCATE(V1,V2)

!deallocate
     deallocate(x , y , z, veloc) !n
     deallocate(xd, yd, zd, veloc_d) !n
     deallocate(triangle_node ) !N
     deallocate(dif_vertex_index ) !n
     deallocate(p_index,pi_index ) 
     deallocate(vd)
     deallocate(near_vert_tri,n_of_tri_vert)
     deallocate(near_vert_vert,n_of_vert_vert)
     deallocate(near_vert_distance)
     deallocate(near_vert_cot    )
     deallocate(near_vert_area  )
     deallocate(near_vert_angle  )
     deallocate(aux)

     if (.not.fsh_sparse) then
      deallocate(a,b,eig,eig_bare,work,  b_bare,a_bare)
     else

     deallocate(isa )
     deallocate(jsa )
     deallocate(sa  )
     deallocate(sac )
     deallocate(rowind_a )
     deallocate(colind_a )
     deallocate(isb )
     deallocate(jsb )
     deallocate(sb  )
     deallocate(sb_bare)
     deallocate(sbc )
     deallocate(rowind_b )
     deallocate(colind_b )
     deallocate(Xv )
     deallocate(Xv_bare )
     deallocate(e, res) 
     deallocate(e_bare) 

     endif
  enddo !bands

 !integratzeko orduan, kontutan hartu beharreko espazioa banda guztien fermi gainazalen batura da.

 !guk, banda bakoitzari dagokion aspiespazioak, bakoitza bere aldetik diagonalizatu dugu, baina hori posiblea 
 !izan da solik elkarrekintza matrizea, bloke diagonal eran agertzen delako (gainazal ezberdinak ez dira gurutzatzen). 

 !hala ere, esan bezala dimentsio osoa
 
 !                 sum_ibn = ndin_p(ibnd) dimentsiokoa da.

 !orduan, energiak ordenatzeko orduan ere batura espazioa da kontutan hartu beharrekoa.

 !hurrengo lerroetan inplementatzen dugu hori.


 WRITE(*,*)"ndif_p_list :", ndif_p_list 
 WRITE(*,*)"ndif_list :", ndif_list 
 WRITE(*,*)"nfs_bnd :", nfsh_bnd
 WRITE(*,*)"nfs_bnd_bare :", nfsh_bnd_bare

   !deallocate(ndif_p_list )
   !deallocate(ndif_list )
!   deallocate(nfsh_bnd)
!   deallocate(nfsh_bnd_bare)
   deallocate(dos_ef)
   deallocate(eulerc)


!  allocate(vd(3,sum(ndif_list)))
!  allocate(energy_list(NBND,maxval(ndif_p_list) )) ! dimentsiorik handienarekin alokatu banda guztien f. gainazala.
!  allocate(el(nbnd))
!  allocate(fsh_index(nbnd,sum(ndif_p_list)) )
!  deallocate(vd)
!  deallocate(el)

end subroutine tria_diag

!--------------------------------------------
  function writenum(i)
  integer,intent(in) :: i
  character(len=4) ::  writenum

     if (1<i<=9) then
        write(writenum,"(i1)")i
     else if (10<=i<=99) then 
        write(writenum,"(i2)")i
     else if (100<=i<=999) then
        write(writenum,"(i3)")i
     else
        write(writenum,*)i
     endif

  end function writenum
!--------------------------------------------
  function numlines_file(iu)
  integer, intent(in) :: iu
  integer ::  numlines_file
  integer :: i,ios
  character(len=1) :: dumm
  integer, parameter :: mxl=100000000

  numlines_file=0
  do 
   read(iu,*,iostat=ios) dumm
   if (ios/=0) exit
   if (numlines_file>mxl) then
    write(*,*)"max number of lines reached:: ERROR numlines_file"   
   endif
   numlines_file=numlines_file+1
  enddo
  rewind(iu)
  end function numlines_file
!--------------------------------------------
  function overlap(t1,t2)
    integer(kind=I8), dimension(3),intent(in) :: t1,t2
    integer(kind=I8) :: overlap 
    integer(kind=I8) :: i,j     

    if ((size(t1)>3).or.(size(t2)>3)) then
       write(unit=*,fmt=*)"errorea. bukatutko dut. overlap", size(t1), size(t2)
       stop
    endif

    overlap=0
    do i=1, size(t1)
       do j=1, size(t2) 
          if (t1(i)==t2(j)) then
             overlap=overlap+1
          endif
       enddo
    enddo

  end function overlap
!--------------------------------------------
  subroutine ij_angle_dist(t1,t2,id,jd,kd1,kd2)
    integer(kind=I8), dimension(3),intent(in) :: t1,t2
    integer(kind=I8), intent(in)  :: id ! erpin komuna.
    integer(kind=I8), intent(out) :: jd, kd1, kd2 ! beste erpinak.

    integer(kind=I8) :: overlap 
    integer(kind=I8) :: i,j,k     

    overlap=0
    do i=1, size(t1)
       if (t1(i)==id) then
          overlap=overlap+1
       endif
    enddo

    do i=1, size(t2)
       if (t2(i)==id) then
          overlap=overlap+1
       endif
    enddo

    if (overlap/=2) then 
       write(*,*) "errorea:id ez da erpin banatu bat"
       write(*,*)"t1", t1
       write(*,*)"t2", t2
       !stop
    endif

    i_l: do i=1, size(t1)
       j_l: do j=1, size(t2)
          if (t1(i)==t2(j).and.(t1(i)/=id)) then
             jd=t1(i)
             exit i_l
          endif
       enddo j_l
    enddo i_l

    do i=1, size(t1)
       if ((t1(i)/=id).and.(t1(i)/=jd)) then
          kd1= t1(i)
          exit
       endif
    enddo

    do i=1, size(t2)
       if ((t2(i)/=id).and.(t2(i)/=jd)) then
          kd2= t2(i)
          exit
       endif
    enddo

  end subroutine ij_angle_dist
!--------------------------------------------
  function ctn(v1,v2)
    real(kind=dp), dimension(3),intent(in) :: v1, v2
    real(kind=dp) :: cosi,ctn

    cosi= scalar(v1,v2)/dsqrt(scalar(v1,v1)*scalar(v2,v2))

    ctn = cosi/dsqrt(1.0_dp-cosi**2) 

  end function ctn
!--------------------------------------------
  function scalar(v1,v2)
    real(kind=dp), dimension(3),intent(in) :: v1, v2

    integer(kind=I8) :: i
    real(kind=dp) :: scalar

    scalar=0.0_dp
    do i=1,3
       scalar = scalar + v1(i)*v2(i)
    enddo

      end function SCALAR
 
    function norma(v)
    real(kind=dp), dimension(3),intent(in) :: v

    integer(kind=I8) :: i
    real(kind=dp) :: norma

     NORMA=0.0_dp
    do i=1,3
     NORMA = norma + v(i)*v(i)
    enddo

    NORMA=SQRT(NORMA)

  end function NORMA
!--------------------------------------------
  function area_vec(v1,v2)
    real(kind=dp), dimension(3),intent(in) :: v1, v2
    real(kind=dp) :: area_vec

    area_vec = (v1(2)*v2(3)-v1(3)*v2(2))**2 + &
         (v1(3)*v2(1)-v1(1)*v2(3))**2 + &
         (v1(1)*v2(2)-v1(2)*v2(1))**2   

    area_vec = sqrt(area_vec) * 0.5_dp


  end function area_vec
!--------------------------------------------
  subroutine mod05(x)
    real(kind=dp), dimension(3),intent(inout) :: x
    integer(kind=I8) :: i 

    do i=1, 3
       if (x(i)     >  0.5_dp) then
          x(i)=x(i)-1.0_dp
       else if (x(i)<=-0.5_dp) then
          x(i)=x(i)+1.0_dp
       endif
    enddo
  end subroutine mod05
!--------------------------------------------
  subroutine mod05_2(x1,x2)
    real(kind=dp), dimension(3),intent(inout) :: x1,x2
    integer(kind=I8) :: i 

    do i=1, 3
       if (x2(i)-x1(i)     >  0.5_dp) then
          x2(i)=x2(i)-1.0_dp
       else if (x2(i)-x1(i)<=-0.5_dp) then
          x1(i)=x1(i)-1.0_dp
       endif
       if (abs(x2(i)-x1(i))>0.5) then 
          write(*,*)"errorea mod05_2"
          WRITE(*,*)"x1", x1
          WRITE(*,*)"x2", x2
          WRITE(*,*)"x2-x1", x2-x1
       endif
    enddo
  end subroutine mod05_2
!--------------------------------------------
  function cross(a, b)
  real(kind=dp), dimension(3)  :: CROSS
  real(kind=dp), dimension(3), intent(in) :: A, B

  CROSS(1) = A(2) * B(3) - A(3) * B(2)
  CROSS(2) = A(3) * B(1) - A(1) * B(3)
  CROSS(3) = A(1) * B(2) - A(2) * B(1)
  end function CROSS

  subroutine hp_sort_r (n, ra, ind)  
    !Subroutine directly copied from QE
    !---------------------------------------------------------------------
    ! sort an array ra(1:n) into ascending order using heapsort algorithm.
    ! n is input, ra is replaced on output by its sorted rearrangement.
    ! create an index table (ind) by making an exchange in the index array
    ! whenever an exchange is made on the sorted data array (ra).
    ! in case of equal values in the data array (ra) the values in the
    ! index array (ind) are used to order the entries.
    ! if on input ind(1)  = 0 then indices are initialized in the routine,
    ! if on input ind(1) != 0 then indices are assumed to have been
    !                initialized before entering the routine and these
    !                indices are carried around during the sorting process
    !
    ! no work space needed !
    ! free us from machine-dependent sorting-routines !
    !
    ! adapted from Numerical Recipes pg. 329 (new edition)
    !
    use kinds, only : DP
    implicit none  
    !-input/output variables
    integer :: n  
    integer :: ind (*)  
    real(DP) :: ra (*)  
    !-local variables
    integer :: i, ir, j, l, iind  
    real(DP) :: rra  
    ! initialize index array
    if (ind (1) .eq.0) then  
       do i = 1, n  
          ind (i) = i  
       enddo
    endif
    ! nothing to order
    if (n.lt.2) return  
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1  
    ir = n  
10  continue  
    ! still in hiring phase
    if (l.gt.1) then  
       l = l - 1  
       rra = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    else  
       ! clear a space at the end of the array
       rra = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       if (ir.eq.1) then  
          ! the least competent worker at all !
          ra (1) = rra  
          !
          ind (1) = iind  
          return  
       endif
    endif
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    do while (j.le.ir)  
       if (j.lt.ir) then  
          ! compare to better underling
          if (ra (j) .lt.ra (j + 1) ) then  
             j = j + 1  
          elseif (ra (j) .eq.ra (j + 1) ) then  
             if (ind (j) .lt.ind (j + 1) ) j = j + 1  
          endif
       endif
       ! demote rra
       if (rra.lt.ra (j) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       elseif (rra.eq.ra (j) ) then  
          ! demote rra
          if (iind.lt.ind (j) ) then  
             ra (i) = ra (j)  
             ind (i) = ind (j)  
             i = j  
             j = j + j  
          else  
             ! set j to terminate do-while loop
             j = ir + 1  
          endif
          ! this is the right place for rra
       else  
          ! set j to terminate do-while loop
          j = ir + 1  
       endif
    enddo
    ra (i) = rra  
    ind (i) = iind  
    goto 10  
    !
  end subroutine hp_sort_r


end module intw_tria_diag


