
SUBROUTINe multiply_psi_by_vKB( k_point, npol, nbands, ng_max, list_iGk, psi, dvnl_psi)
  !INTW project: KB projection by wave functions.
  ! 

  USE intw_useful_constants, only: pi,tpi,fpi, cmplx_0, cmplx_i
  USE intw_reading, only : tpiba, alat,  bg, nat, lspinorb
  USE intw_reading,    ONLY : ntyp, volume0, ityp
  USE intw_fft, ONLY : gvec_cart
  USE kinds, only : dp
  USE intw_pseudo, ONLY : upf, lmaxkb, nhm, nh
  USE intw_pseudo,       ONLY : nkb, vkb, DKB
  USE intw_reading, ONLY : noncolin
  implicit none

  integer, intent(in)               :: npol, nG_max, list_iGk(nG_max), nbands
  real(dp), intent(in)              :: k_point(3)
  complex(dp), intent(inout)        :: psi(nG_max, nbands,npol), dvnl_psi(nG_max,nbands, npol)

  complex(dp)                       :: projec_d(nhm,3,npol)
  integer                           :: spol1, spol2 , ipol, jpol, na, ih, ig
  real(dp)                          :: k_(3)
  integer :: nt, mode, ikb, iGk, ibnd

  k_(:)=matmul(bg, k_point)

  do ibnd=1,nbands

     projec_d(:,:,:,:) = cmplx_0

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

                       do spol1 = 1,npol !spin

                          ! DKB(nhm,nhm,nspin,nspin,ntyp)
                          projec_d(ih,ipol,spol1) = projec_d(ih,ipol,spol1) +   conjg(vkb(iG,ikb)) * psi(iG, ibnd, spol1)

                       enddo !spol

                    enddo !ig

                 enddo !ih
              end if
           enddo !na
        enddo ! nt
     enddo !ipol

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

                       if (upf(nt)%has_so.and.lspinorb) then

                          do spol1 = 1,npol !spin
                             do spol2 = 1,npol !spin
                                dvnl_psi(iG, ibnd,spol1,mode) = dvnl_psi(iG,ibnd,spol1,mode)&
                                     + DKB (ih, ih,spol1,spol2,ityp(na) ) * projec_d(ih,ipol,spol2) * vkb(iG,ikb)
                             end do !spol2
                          end do !spol1

                       else ! no SO

                          do spol1 = 1,npol !spin

                             dvnl_psi(iG, ibnd,spol1,mode) = dvnl_psi(iG,ibnd,spol1,mode)&
                                  + DKB (ih, ih,spol1,spol1,ityp(na) )* projec_d(na,ih,ipol,spol1) * vkb(iG,ikb)

                          end do
                       end if

                    enddo !ig
                 enddo!ih
              end if
           enddo! na
        enddo !nt
     enddo !i pol

  enddo !ibnd

  return
end SUBROUTINe multiply_psi_by_vKB
