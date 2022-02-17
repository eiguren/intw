
SUBROUTINe multiply_psi_by_KB( k_point, q_point, npol, nbands,ng_max, list_iGk, list_iGkq,psi, dvnl_psi)

  USE intw_useful_constants, only: pi,tpi,fpi, cmplx_0, cmplx_i
  USE intw_reading, only : tpiba, alat,  bg, nat 
  USE intw_reading,    ONLY : ntyp, volume0, ityp
  USE intw_fft, ONLY : gvec_cart
  USE kinds, only : dp 
  USE intw_pseudo, ONLY : upf, lmaxkb, nhm, nh
  USE intw_pseudo,       ONLY : nkb, vkqb, vkb, nhtol, nhtolm, indv,  dvan
  USE intw_reading, ONLY : noncolin
  implicit none

  integer, intent(in)               :: npol, nG_max, list_iGk(nG_max), list_iGkq(nG_max), nbands
  real(dp), intent(in)              :: k_point(3), q_point(3)

  complex(dp), intent(inout)          :: psi(nG_max, nbands,npol), dvnl_psi(nG_max,nbands, npol,3*nat)

  complex(dp)                       :: projec_d1(nat,nhm,3,npol), projec_d2(nat,nhm,3,npol)
  integer                           :: spol, ipol, jpol, na, ih, ig
  real(dp)                          :: k_(3), q_(3)
  integer :: nt, mode, ikb, iGk, ibnd, iGkq

  k_(:)=matmul(bg, k_point)
  q_(:)=matmul(bg, q_point)



  do ibnd=1,nbands

     projec_d1(:,:,:,:) = cmplx_0
     projec_d2(:,:,:,:) = cmplx_0 


     ! Asier: KB potentziala hurrengo eran emanik dago: sum_l |b(l)> <b(l)| 
     !               beraz deribatuak bi gai dauzka:
     !               sum_l d|b(l,r)> <b(l,r)| + |b(l,r)> d<b(l,r)| ~ Fourier ~
     !               sum_l i(k+G)|b(l,G)> <b(l,G)| + |b(l,G)> <b(l,G)|
     !               

     do spol = 1,npol !spin
        do ipol = 1, 3  !cart coord.
           ikb=0
           do nt=1,ntyp
              do na = 1, nat !
                 if (ityp(na) == nt) then 

                 do ih=1,nh(ityp(na))
                    ikb=ikb+1

                    do ig=1,nG_max
                       iGk= list_iGk(iG)
                       if (iGk==0) exit
                       projec_d1(na,ih,ipol,spol) =   &
                       projec_d1(na,ih,ipol,spol) -   ( tpiba * cmplx_i) * conjg(vkb(iG,ikb)) * psi(iG, ibnd, spol)

                       projec_d2(na,ih,ipol,spol) =   &
                       projec_d2(na,ih,ipol,spol) +   ( tpiba * cmplx_i )* conjg(vkb(iG,ikb)) * psi(iG, ibnd, spol)* & 
                            ( k_(ipol)  + gvec_cart(ipol,iGk) ) 
                             
                    enddo !ig

                 enddo !ih
                 end if
              enddo !na
           enddo ! nt
        enddo !ipol
     enddo !spol

     do spol = 1,npol !spin
        do ipol = 1, 3  !cart coord.
           ikb=0 
           do nt=1,ntyp
           do na = 1, nat !
              if (ityp(na) == nt) then 
              mode=(na-1)*3+ipol

              do ih=1,nh(ityp(na))
                 ikb=ikb+1
                 do ig=1,nG_max
                    iGk= list_iGkq(iG)
                    if (iGk==0) exit                               
                    dvnl_psi(iG, ibnd,spol,mode) = dvnl_psi(iG,ibnd,spol,mode)&
                         + dvan (ih, ih,ityp(na) ) * projec_d2(na,ih,ipol,spol) * vkqb(iG,ikb)
                    dvnl_psi(iG, ibnd,spol,mode) = dvnl_psi(iG,ibnd,spol,mode)&
                         + dvan (ih, ih, ityp(na)) * projec_d1(na,ih,ipol,spol) * vkqb(iG,ikb)* &
                         (k_(ipol) +  q_(ipol) + gvec_cart(ipol,iGk) )
                 enddo !ig
              enddo!ih
              end if
            enddo! na
           enddo !nt
        enddo !i pol
     enddo !spol : This refers to the spinor components. 

  enddo !ibnd

  return
end SUBROUTINe multiply_psi_by_KB
