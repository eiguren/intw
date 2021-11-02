module tria_const
  integer, parameter :: i8=selected_int_kind(7)

  PUBLIC :: I8
 

end module tria_const



module intw_tria_diag

use intw_useful_constants

use  tria_const

public :: cross, scalar, norma

integer(KIND=I8), allocatable, public :: ndif_p_list(:), fsh_index(:,:)
integer(KIND=I8), allocatable, private :: ndif_list(:)


contains 

subroutine tria_diag(alat, bcell, nbnd, bnd_index, nmax, nfsh )


  implicit none

  integer :: nmax
  integer,intent (in)                      :: nfsh, nbnd

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

  integer(kind=I8), parameter :: nt_max=14
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
  real(kind=dp), allocatable ::  a(:,:),eig(:), work(:), b(:,:)
   
  real(kind=dp), parameter :: eps=10e-8

  real(kind=dp), allocatable :: dos_ef (:)

  real(kind=dp) ::  bz_vol, v_mean
  character(len=8) :: num 
  integer :: ibn

  !-full fermi surface-space variables/parameters
  real(kind=dp), allocatable :: energy_list(:,:),el(:)
  !integer(KIND=I8), allocatable :: ndif_p_list(:), ndif_list(:), fsh_index(:,:)
  integer :: num_fsh, indx(3)

  !asier TEST
   
  real(kind=dp), allocatable :: V1(:), V2(:)
  real(KIND=DP)  :: PROD

  

  allocate (ndif_p_list(nbnd))
  allocate (ndif_list(nbnd))
  allocate (dos_ef(nbnd))

  dos_ef(:)=0.0
 
  !bcell=transpose(bcell)

  do ibn=1,nbnd 
     ibnd = bnd_INDEX(ibn)

     num = writenum(ibnd)

     write(*,*)"bnd_INDEX(ibn)", ibn, bnd_INDEX(ibn),"num", trim(num)
     open(unit=14,file=TRIM("TRI_01.DAT_"//NUM),status="unknown")

     read(14,*)ntri

     write(*,*)"ntri", ntri

     allocate(x (ntri), y(ntri) , z(ntri), veloc(Ntri)) !n
     allocate(xd(ntri), yd(ntri), zd(ntri), veloc_d(ntri)) !n

     allocate(triangle_node(3,ntri/3) ) !N
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

     !open(unit=123,file=TRIM("TRI_VERT_INDX_01.DAT_"//NUM),status="unknown")
     !do i=1,ntri,3
     !   write(123,*) dif_vertex_index(i),dif_vertex_index(i+1),dif_vertex_index(i+2)
     !enddo
     !close(unit=123)

     !open(unit=123,file=TRIM("VERT_COORD_01.DAT_"//NUM),status="unknown")
     !write(123,*)ndif
     !do i=1,ndif
     !   vi=(/xd(i), yd(i),zd(i)/)
     !   vi=matmul(bcell,vi)
     !   write(123,"(100F18.12)") vi, veloc_d(i)
     !enddo
     !close(unit=123)

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

     allocate(near_vert_tri(nt_max,ndif),n_of_tri_vert(ndif))
     allocate(near_vert_vert(nt_max,3,ndif),n_of_vert_vert(ndif))
     allocate(near_vert_distance(nt_max,ndif))
     allocate(near_vert_cot     (nt_max,ndif))
     allocate(near_vert_area    (  ndif))
     allocate(near_vert_angle    (  ndif))
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

     write(*,*)"allocating", ndif_P," X ", ndif_P, " matrix."
     allocate(a(ndif_p,ndif_p),b(ndif_p,ndif_p),eig(ndif_p),work(ndif_p*(3+ndif_p/2)))
     a(:,:)=0.0_dp

     total_area = 0.0_dp
     near_vert_distance=0.0_dp

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
           near_vert_angle    (id)  = near_vert_angle    (id)  + (tau1  + tau2)/(2*pi)

           a(id,near_vert_vert(j,1,id)) = 0.5_dp *( ctn1 + ctn2)

        enddo

     end do !id MAIN

     total_area=0.0_dp
     do id=1,ndif_p
        total_area = total_area +  near_vert_area    (id)
     enddo

     print*,"Total area ", total_area

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
      
     b(:,:)=0.0_dp
     do i=1,ndif_P
        b(i,i) = near_vert_area  (i)/ veloc_d(p_index(i)) ! /total_area
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

     write(*,*)"diagonali...."

     l=ndif_p*(3+ndif_p/2)
     do i=1,ndif
        id=pi_index(i)
        vi=(/xd(i),yd(i),zd(i)/)
     enddo
     
     CALL dsygv( 1, 'V', 'U', NDIF_P, a, NDIF_P, b, NDIF_P, EIG, work, l, inf )

     open(unit=100,file=TRIM("FSH_01.DAT_"//NUM),status="unknown")
     open(unit=102,file=TRIM("FSH_ENERGY.DAT_"//NUM),status="unknown")

     do i=1,ndif_p
        write(102,*)i,eig(i)
     enddo
 
     !e_fsh(ibn,1:min(nfsh,ndif_p)) = eig(1:min(nfsh,ndif_p))

     do i=1,ndif
        id=pi_index(i)
 
        !fsh(i, ibn,1:min(nfsh,ndif_p)) =  a(pi_index(i), 1:min(nfsh,ndif_p) ) 

        vi=(/xd(i),yd(i),zd(i)/)
        vi=matmul(bcell,vi)
        write(100,"(1000F18.12)")vi, 1/veloc_d(i), near_vert_area  (pi_index(i)), (a(pi_index(i),j),j=1,min(ndif_p,nfsh))
     enddo

     close(unit=100)
     close(unit=102)

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
     deallocate(a,b,eig,work) 

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


!  allocate(vd(3,sum(ndif_list)))
!  allocate(energy_list(NBND,maxval(ndif_p_list) )) ! dimentsiorik handienarekin alokatu banda guztien f. gainazala.
!  allocate(el(nbnd))
!  allocate(fsh_index(nbnd,sum(ndif_p_list)) )
!
!  energy_list=10**100
!
!  do ibn=1,nbnd 
!
!     ibnd = bnd_INDEX(ibn)
!
!     num = writenum(ibnd)
!
!     open(unit=1000+ibn,file=TRIM("FSH_ENERGY.DAT_"//NUM),status="unknown")
!       do i=1,ndif_p_list(ibn)
!         read(unit=1000+ibn,fmt=*)j,energy_list(ibn,i)
!       enddo 
!     close(unit=1000+ibn)
!  enddo
!
!  open(unit=1004,file="FULLSPACE_FSH_LIST.DAT",status="unknown")
!
!
!  fsh_index(:,1) = 1
!  do ibn=1,nbnd
!     el(ibn) = energy_list(ibn,1)
!     write(1004,"(3(2X,I4),4X,E12.6)")ibn , bnd_INDEX(ibn), 1,  energy_list(ibn,1)
!  enddo
! 
!  num_fsh=1
!
!  !-gainazal ezberdinen konbinazioak emandako fsh moduak ordenatu energian. 
!  do i=1, sum (ndif_p_list)-nbnd
!
!    !write(1004,"(3(2X,I4),4X,E12.6)")num_fsh, bnd_INDEX(minloc(el)), fsh_index(minloc(el), num_fsh), minval(el)
!
!    num_fsh=num_fsh+1
!    do ibn=1,nbnd
!     el(ibn) = energy_list(ibn,fsh_index(ibn,num_fsh-1)+1)
!    enddo
!    
!    fsh_index(:, num_fsh)          = fsh_index(:, num_fsh-1)
!    fsh_index(minloc(el), num_fsh) = fsh_index(minloc(el), num_fsh-1)+1
!
!    write(1004,"(3(2X,I4),4X,E12.6)")num_fsh+ibn, bnd_INDEX(minloc(el)), fsh_index(minloc(el), num_fsh), minval(el) !fsh_index(:, num_fsh)
!  enddo
!
!  close(unit=1004)
!
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
  integer, parameter :: mxl=100000

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

end module intw_tria_diag


