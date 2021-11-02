
   

MODULE WS 

PUBLIC :: WignerSeitz, gett3, det3x3d, vecsize

REAL*8,PARAMETER   :: MAXERR=10d-10, MAXERR2=10D-8
INTEGER, PARAMETER :: MAX_POLY=100   , &
                            MAX_VERT=400   , &
                            MAX_PL_VER=100 , &
                            MAX_NORD = 100


CONTAINS


      SUBROUTINE WignerSeitz(vec)
      implicit none

      REAL*8 vec(3,3)


!     this SUBROUTINE is not efficient but is very simple !!!
      REAL*8 :: vertex(3,MAX_VERT), & !VERTICES OF VORONOI POLIHEDRA
      point(4,26),       & !equation of plane: ax+by+cz=d; point = (a,b,c,d)
      vect1(3),          & !arbitrary vector 1
      vect2(3),          & !arbitrary vector 2
      vnml(3),            & ! nml=vect1 x vect2
      t(3)

      INTEGER    ::  vindex(3,MAX_VERT), & !information of vertex's indices
      indexlist(MAX_POLY)              , & !list of all-different indices
      iverplane(MAX_PL_VER,MAX_POLY)   , & !indices of vertices in plane/polygon
      iverplane2(MAX_PL_VER,MAX_POLY)  , & !copy of iverplane
      nverplane(MAX_POLY)              , & !number of vertices in plane/polygon
      icc(3)
      LOGICAL ::  activ(MAX_PL_VER,MAX_POLY)  , &  !flag for iverplane(,)
      newvertex, lnew, lplane, equal12, equal13, equal23, equalvert 
      INTEGER i, j, k, ic, jc, kc, icc1, lc, lc1, ideb, ii, jj, kk, &
      iverIC, iverJC, iverKC, l, nind, nplane, nver, nvpl, nvpl2, npoi, iver 
      REAL*8 :: det, detx, dety, detz, di, dj, dk, dist1, dist2, &
      px, py, pz, vsize
      


!     NPOI.........number of point
!     NVER.........number of vertices
!     NIND.........number of different indices and planes that polihedra 
!                  is made of
!     NVPL.........number of vertex in polygon/plane

!     ===========================================
!           |11   21   31|           |x1  x2  x3|
!     VEC = |12   22   32|  => VEC = |y1  y2  y3|
!           |13   23   33|           |z1  z2  z3|
!     ===========================================

!     DEBUG
!      print *,'DEBUG> ', word
!      print *,'DEBUG> ', vec
!     /DEBUG

      do i=1,MAX_POLY
      do j=1,MAX_PL_VER
      activ(MAX_PL_VER,MAX_POLY) = .false.
      enddo
      enddo


!     make all neighbour points of origin (0,0,0); there are 26 such points
      NPOI=0
      do i = -1, 1
      do j = -1, 1
      do k = -1,1
      if(i.ne.0 .or. j.ne.0 .or. k.ne.0)then !disregard point(0,0,0)
      di=dble(i)
      dj=dble(j)
      dk=dble(k)
      NPOI=NPOI+1
!     equation of plane: ax+by+cz=d, point = (a,b,c,d)
      point(1,NPOI)=di*vec(1,1) + dj*vec(2,1) + dk*vec(3,1)
      point(2,NPOI)=di*vec(1,2) + dj*vec(2,2) + dk*vec(3,2)
      point(3,NPOI)=di*vec(1,3) + dj*vec(2,3) + dk*vec(3,3)
      point(4,NPOI)=0.5d0 * & 
      ( point(1,NPOI)*point(1,NPOI) + point(2,NPOI)*point(2,NPOI) + point(3,NPOI)*point(3,NPOI))
      endif
      enddo
      enddo
      enddo

!     find vertices as intersections of three planes; 
!     planes are made at the middle of origin and points
      NVER=0
      do i=1,NPOI-2
      do j=i+1,NPOI-1
      do k=j+1,NPOI
      det=det3x3d(                                 &
                      point(1,i), point(1,j), point(1,k),     &
                      point(2,i), point(2,j), point(2,k),     &
                      point(3,i), point(3,j), point(3,k))
      if(abs(det).gt.MAXERR)then   
      detx=det3x3d(                             &
                         point(4,i), point(4,j), point(4,k),  &
                         point(2,i), point(2,j), point(2,k),  &
                         point(3,i), point(3,j), point(3,k))
  
      dety=det3x3d(                             &
                         point(1,i), point(1,j), point(1,k),  &
                         point(4,i), point(4,j), point(4,k),  &
                         point(3,i), point(3,j), point(3,k))

      detz=det3x3d(                             &
                         point(1,i), point(1,j), point(1,k),  &
                         point(2,i), point(2,j), point(2,k),  &
                         point(4,i), point(4,j), point(4,k))

      px=detx/det
      py=dety/det
      pz=detz/det

!     if vertex is closer to origin than to any other point, 
!     than we have new vertex
      l=0
      newvertex=.true.
      do while (newvertex .and. l.lt.NPOI)
      l=l+1
      dist1=dsqrt(px*px + py*py + pz*pz)
      dist2=dsqrt(                           &
                            (px-point(1,l))*(px-point(1,l)) + &
                            (py-point(2,l))*(py-point(2,l)) + &
                            (pz-point(3,l))*(pz-point(3,l)))
      if(dist1.gt.(dist2 + MAXERR)) newvertex=.false.
      enddo
      
      if(newvertex) then
      NVER=NVER+1
      vertex(1,NVER) = px
      vertex(2,NVER) = py
      vertex(3,NVER) = pz
      vindex(1,NVER) = i
      vindex(2,NVER) = j
      vindex(3,NVER) = k
      !write(111,"(100f12.6)")vertex(1:3,NVER)
      endif
      endif
      enddo
      enddo
      enddo

!     --- --- DEBUG_BEGIN --- ---
      !PRINT *,'ALL VERICES OF VORONOI: NVER=',NVER
      do i=1,NVER
      !WRITE(*,*) (vertex(j,i),j=1,3),'; --> ',                &
      !          (vindex(j,i),j=1,3)
      
      !write(111,"(100f12.6)")vertex(1:3,i)
      enddo
      !WRITE(*,*) ' '
!     --- --- DEBUG_END --- ---

!     number of planes that a polyhedra is made of, is the number of different
!     indices, that appear in vindex(,). Each index represents one plane.

!     search for all different indices
      NIND=0
      do i=1,NVER
      do ii=1,3
      lnew=.true.
      j=1
      do while(lnew .and. j.le.NIND)                  
      if(indexlist(j).eq.vindex(ii,i)) lnew=.false.
      j=j+1
      enddo
      if(lnew)then
      NIND=NIND+1
      indexlist(NIND)=vindex(ii,i)
      endif
      enddo         
      enddo

!     --- --- DEBUG_BEGIN --- ---
      !WRITE(*,*)'INDEX LIST:'
      !do i=1,NIND
      !WRITE(*,*) i,':  ',indexlist(i)
      !enddo
      !WRITE(*,*)' '
!     --- --- DEBUG_END --- ---


!     so we have NIND different planes; if vertex has index M, that means
!     that it is a member of Mth plane; make plane data structure
!
!     some vertices that belong to one plane may be identical or colinear, 
!     get rid of that vertices
      NPLANE=0
      do i=1,NIND
      NVPL=0                 !number of vertex in polygon/plane
      lplane=.true.
      do j=1,NVER
      do jj=1,3
      if(vindex(jj,j).eq.indexlist(i)) then 
!     vertex 'j' belongs to plane 'i'
      NVPL=NVPL+1
      iverplane(NVPL,i)=j
      activ(NVPL,i)=.true.
      nverplane(i)=NVPL
      endif
      enddo            
      enddo
!     if NVPL<3 -> this is not a plane; disregard
      if(NVPL.lt.3) lplane=.false.
!     get rid of colinear and identical points
      if(lplane)then

      do ic=1,NVPL-2
      if(activ(ic,i)) then

      do jc=ic+1,NVPL-1   
      if(activ(jc,i))then

      do kc=jc+1,NVPL
      if(activ(kc,i))then
      
      iverIC=iverplane(ic,i)
      iverJC=iverplane(jc,i)
      iverKC=iverplane(kc,i)
!     get rid of identical points
      equal12 =                            &
             abs(vertex(1,iverIC)-vertex(1,iverJC)).lt.MAXERR2 .and. &
             abs(vertex(2,iverIC)-vertex(2,iverJC)).lt.MAXERR2 .and. &
             abs(vertex(3,iverIC)-vertex(3,iverJC)).lt.MAXERR2

      equal13 =                            &
             abs(vertex(1,iverIC)-vertex(1,iverKC)).lt.MAXERR2 .and. &
             abs(vertex(2,iverIC)-vertex(2,iverKC)).lt.MAXERR2 .and. &
             abs(vertex(3,iverIC)-vertex(3,iverKC)).lt.MAXERR2

      equal23 =                            &
             abs(vertex(1,iverJC)-vertex(1,iverKC)).lt.MAXERR2 .and. &
             abs(vertex(2,iverJC)-vertex(2,iverKC)).lt.MAXERR2 .and. &
             abs(vertex(3,iverJC)-vertex(3,iverKC)).lt.MAXERR2
      
      equalvert=.true.
      if(equal12.and.equal13)then 
!     all three points are identical
!     disregard with lower indices: ic < jc < kc
      activ(ic,i)=.false.
      activ(jc,i)=.false.
      elseif(equal12)then
      activ(ic,i)=.false.
      elseif(equal13)then
      activ(ic,i)=.false.
      elseif(equal23)then
      activ(jc,i)=.false.
      else
      equalvert=.false.
      endif

!     --- --- DEBUG_BEGIN --- ---
!                              PRINT *,'i=',i,';  equalvert=',equalvert
!     --- --- DEBUG_BEGIN --- ---

      
      if(.not.equalvert)then
      do lc=1,3
      vect1(lc) = vertex(lc,iverJC)- vertex(lc,iverIC) 
      vect2(lc) = vertex(lc,iverKC)- vertex(lc,iverIC) 
      enddo

      CALL VecProduct(vect1,vect2,vnml)
      vsize = VecSize(vnml)

!     --- --- DEBUG_BEGIN --- ---
!                              PRINT *,'SIZE>',
!                                       (vertex(lc,iverplane(ic,i))
!                                    ,vertex(lc,iverplane(jc,i))
!                                    ,vertex(lc,iverplane(kc,i)),lc=1,3)
!                              PRINT *,'SIZE::',vsize,'MAXERR2::',MAXERR2
!     --- --- DEBUG_END --- ---
      if(abs(vsize).lt.MAXERR2)then
!     colinear points; disregard middle point
      icc(1)=ic
      icc(2)=jc
      icc(3)=kc
      t(1)=0.0d0
      t(2)=1.0d0
      t(3)=gett3(vertex(1,iverIC),vertex(1,iverJC),vertex(1,iverKC))
      do lc=2,3
      icc1=icc(lc)
      do lc1=lc-1,1,-1
      if(t(lc1).lt.t(icc1))  goto 101
      icc(lc1+1)=icc(lc1)
      enddo
      lc1=0
 101                                   icc(lc1+1)=icc1
      enddo

!     --- --- DEBUG_BEGIN --- ---
!                                    PRINT *,'i=',i,                 &
!                                          ';SORTING ORDER OF VERT:'
!                                    PRINT *,'ic=',ic,';  ',         &
!           (vertex(iver,iverplane(icc(1),i)),iver=1,3)
!                                    PRINT *,'jc=',jc,';  ',         &
!           (vertex(iver,iverplane(icc(2),i)),iver=1,3)
!                                    PRINT *,'kc=',kc,';  ',         &
!           (vertex(iver,iverplane(icc(3),i)),iver=1,3)

!     --- --- DEBUG_END --- ---
      activ(icc(2),i)=.false. !middle point disregarded
      endif !if(abs(vsize).lt.MAXERR2)then
      endif !if(.not.equalvert)then
      endif
      enddo   !kc
      endif
      enddo         !jc
      endif
      enddo               !ic
      
      if(lplane)then
      NVPL2=0
      do ic=1,NVPL
!     --- --- DEBUG_BEGIN --- ---
!                  PRINT *,'ic=',ic,';  activ(ic,i)=',activ(ic,i)
!                  PRINT *,(vertex(lc,iverplane(ic,i)),lc=1,3)
!     --- --- DEGUB_END --- ---
      if(activ(ic,i))then
      NVPL2=NVPL2+1
      iverplane2(NVPL2,nplane+1)=iverplane(ic,i)
      nverplane(nplane+1)=NVPL2
      endif
      enddo
      endif

      endif
     
      open (unit=123,file="polyhedra.dat", status="unknown")
 
      if(lplane .and. NVPL2.gt.2) then
      NPLANE=NPLANE+1
!     --- --- DEBUG_BEGIN --- ---
      !WRITE(*,*)'PLANE N.:',i,'; NVPL2=',NVPL2,'; NVPL',NVPL,  &
      !             ';  INDEX LIST/VERTEX LIST'
      write(123,*) nplane, NVPL2 
      do ideb=1,NVPL2
       write(123,"(3(x,f18.15))")  (vertex(j,iverplane2(ideb,nplane)),j=1,3)
      enddo
      !WRITE(*,*)'NUMBER OF PLANES: ',nplane
!     --- --- DEBUG_END --- ---

      NVPL=0
!     now we have nverplane(i) vertices in plane i; make Convex Hull
!     for that plane
      !CALL ConvexHull(nplane,iverplane2,NVPL2,vertex)
      endif
      enddo

      close(unit =123)

      RETURN
      END SUBROUTINE WignerSeitz

!     ====================================
!     equation of line: x = x1 + t(x2-x1)
!     function return a t3 for third point
!     ====================================
      FUNCTION GETT3(ver1,ver2,ver3)      
      IMPLICIT none
      REAL*8 ver1(3), ver2(3), ver3(3)  !three vertices
      REAL*8 dx, dy, dz, sum,GETT3
      
      dx = ver2(1)-ver1(1) 
      dy = ver2(2)-ver1(2) 
      dz = ver2(3)-ver1(3)
      sum = dx+dy+dz
      if(abs(sum).lt.MAXERR) then
      if(abs(dx).gt.MAXERR) then
      sum=dx
      elseif(abs(dy).gt.MAXERR) then
      sum=dy
      elseif(abs(dy).gt.MAXERR) then
      sum=dz
      else
      gett3=0.0d0  !three points overlap
      return
      endif
      endif

      gett3 = (ver3(1)+ver3(2)+ver3(3) - ver1(1)+ver1(2)+ver1(3)) / sum
      return
      END FUNCTION GETT3

      FUNCTION DET3X3D(x11,x12,x13,x21,x22,x23,x31,x32,x33)
      REAL*8 :: x11,x12,x13,x21,x22,x23,x31,x32,x33, DET3X3D
      
      DET3X3D=x11*x22*x33 + x12*x23*x31 + x13*x21*x32- &
             x31*x22*x13 - x32*x23*x11 - x33*x21*x12
      RETURN
      END FUNCTION DET3X3D

      SUBROUTINE VecProduct(vec1,vec2,resvec)
      REAL*8 :: vec1(3),vec2(3),resvec(3)
      
      resvec(1) = vec1(2)*vec2(3) - vec2(2)*vec1(3)
      resvec(2) = vec1(3)*vec2(1) - vec2(3)*vec1(1)
      resvec(3) = vec1(1)*vec2(2) - vec2(1)*vec1(2)
      RETURN
      END SUBROUTINE VecProduct

!     ============================
      FUNCTION VecSize(vec)
      REAL*8 :: vec(3), VecSize
      
      VecSize = dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
      RETURN
      END FUNCTION VecSize

END MODULE WS
