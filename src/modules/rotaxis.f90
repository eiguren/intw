      subroutine rotaxis(symm,v,ang)
      use intw_reading, only: at,bg
      implicit none
      integer :: symm(3,3)
      real(8) :: A(3,3), P(3,3), v(3),vst(1,3),vs(3,1),ang,det, Su(3,3),S(3,3)
      real(8) :: I3(3,3), Inv(3,3),diff, eps
      complex(8) :: eig(3)
      integer :: i,j, nu,mu
      logical :: inv_used
      external :: det

      a=0.d0

      do mu=1,3
      do nu=1,3
        do i=1,3
          do j=1,3
              a(mu,nu)  =    a(mu,nu)     &
                                      +  at(mu,j)*symm(i,j)*bg(nu,i)
          end do
        end do
      end do
      end do

      eps=0.00001d0

      I3(:,:)=0.d0
      I3(1,1)=1.d0
      I3(2,2)=1.d0
      I3(3,3)=1.d0
      Inv=I3*(-1.d0)

      call diagreal(A,eig,ang,v)

      !write(*,*)'ANG',ANG

      vst(1,:)=v(:) ; vs(:,1)=v(:)

      Su(:,:)=0.d0
      Su(1,2)=-v(3) ; Su(2,1)= v(3)
      Su(1,3)= v(2) ; Su(3,1)=-v(2)
      Su(2,3)=-v(1) ; Su(3,2)= v(1)

      !Reconstruct the symmetry matrix. with the unit vector u and the angle=ang.
      S=dcos(ang)*I3+(1.d0-dcos(ang))*matmul(vs,vst)+ Su*sin(ang)

      if (det(a).lt.0) S=-S

      diff=0.d0
      do i=1,3
       do j=1,3
        diff=diff+abs(A(i,j)-S(i,j))
       enddo
      enddo

      if (diff.gt.eps) then
       !If v does not work we try -v
       v=-v

       vst(1,:)=v(:) ; vs(:,1)=v(:)

       Su(:,:)=0.d0

       Su(1,2)=-v(3) ; Su(2,1)= v(3)
       Su(1,3)= v(2) ; Su(3,1)=-v(2)
       Su(2,3)=-v(1) ; Su(3,2)= v(1)

       S=dcos(ang)*I3+(1.d0-dcos(ang))*matmul(vs,vst)+ Su*sin(ang)

       if (det(a).lt.0) S=-S

       diff=0.d0
       do i=1,3
        do j=1,3
         diff=diff+abs(A(i,j)-S(i,j))
        enddo
       enddo
      endif

      if (diff.gt.0.01) then
        write(*,*)'ERROR rotaxis.f90: I am not able to reproduce'
        write(*,*)'the S matrix.', diff

        do i=1,3
         write(*,'(100f12.6)')(A(i,j),j=1,3)
        enddo
        write(*,*)
        do i=1,3
         write(*,'(100f12.6)')(S(i,j),j=1,3)
        enddo

        stop
      endif

      inv_used=.false.
      if (abs(det(a)+1.0).lt.0.01) inv_used=.true.

      return
      end

      Subroutine diagreal(A,eig,ang,axis)
      implicit none
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=3,NOUT=3)
      INTEGER          NB, NMAX
      PARAMETER        (NB=64,NMAX=3)
      INTEGER          LDA, LDVS, LWORK
      PARAMETER        (LDA=NMAX,LDVS=NMAX,LWORK=(2+NB)*NMAX)
!     .. Local Scalars ..
      INTEGER          I, IFAIL, INFO, J, LWKOPT, N, SDIM
!     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), VS(LDVS,NMAX), WI(NMAX), &
                       WORK(LWORK), WR(NMAX)
      LOGICAL          BWORK(NMAX)
!     .. External Functions ..
      LOGICAL          SELECT
      EXTERNAL         SELECT
!     .. External Subroutines ..
      EXTERNAL         DGEES
!     .. Executable Statements ..
      COMPLEX(8)       EIG(3)
      REAL(8)          det,Tr, Ar(3,3),ang, I3(3,3), Inv(3,3), axis(3),deter

      I3(:,:)= 0.d0
      I3(1,1)=1.d0
      I3(2,2)=1.d0
      I3(3,3)=1.d0
      Inv=I3*(-1.d0)

      Ar(:,:)=A(:,:)

      deter=det(A)
      if (det(A).lt.0.) Ar=matmul(Inv,A)

      tr=Ar(1,1)+Ar(2,2)+Ar(3,3)
      ang = acos((tr-1.d0)/2.d0)!/(2.*acos(-1.d0))*360.

      N=3
!
         CALL DGEES('Vectors (Schur)','Sort',SELECT,N,Ar,LDA,SDIM,WR,WI, &
                   VS,LDVS,WORK,LWORK,BWORK,INFO)
         LWKOPT = WORK(1)

      axis(:)=Vs(:,1)

      do i=1,3
       eig(i)=cmplx(wr(i),wi(i))
      enddo

      ang=abs(imag(log(eig(2))))

      END

      LOGICAL FUNCTION SELECT(AR,AI)
!
      DOUBLE PRECISION        AI, AR
!     .. Local Scalars ..
      LOGICAL                 D
      COMPLEX(8)              C

      C=cmplx(ar,ai)
!     .. Executable Statements ..
      IF (abs(c-(1.d0,0.d0)).lt.0.0001) THEN
         D = .TRUE.
      ELSE
         D = .FALSE.
      END IF
!
      SELECT = D
!
      RETURN
      END

      function det(A)
      real(8) :: A(3,3), det

      det = A(1,1) *A(2,2) *A(3,3) +  A(1,2) *A(2,3) *A(3,1) + A(1,3) *A(2,1) *A(3,2) &
          - A(1,3) *A(2,2) *A(3,1) -  A(1,1) *A(2,3) *A(3,2) - A(1,2) *A(2,1) *A(3,3)

      end
