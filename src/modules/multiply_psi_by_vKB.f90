
SUBROUTINe multiply_psi_by_vKB( k_point, npol, nbands, ng_max, list_iGk, psi, dvnl_psi)
  !INTW project: KB projection by wave functions.
  !

  USE kinds, only : dp
  USE intw_useful_constants, only: pi, tpi, fpi, cmplx_0, cmplx_i
  USE intw_reading, only : bg, nat, lspinorb, ityp
  USE intw_pseudo, ONLY : upf, nhm, nh, vkb, DKB
  implicit none

  integer, intent(in)               :: npol, nG_max, list_iGk(nG_max), nbands
  real(dp), intent(in)              :: k_point(3)
  complex(dp), intent(inout)        :: psi(nG_max, nbands,npol), dvnl_psi(nG_max,nbands, npol)

  complex(dp)                       :: projec_d(nat,nhm,npol)
  integer                           :: spol1, spol2 , na, ih, ig
  real(dp)                          :: k_(3)
  integer :: nt, ikb, iGk, ibnd

  k_(:)=matmul(bg, k_point)

  do ibnd=1,nbands

     projec_d(:,:,:) = cmplx_0

     ikb=0
     do na = 1, nat !
        do ih=1,nh(ityp(na))
           ikb=ikb+1

           do ig=1,nG_max
              iGk= list_iGk(iG)
              if (iGk==0) exit

              do spol1 = 1,npol !spin

                 ! DKB(nhm,nhm,nspin,nspin,ntyp)
                 projec_d(na,ih,spol1) = projec_d(na,ih,spol1) +   conjg(vkb(iG,ikb)) * psi(iG, ibnd, spol1)

              enddo !spol

           enddo !ig

        enddo !ih
     enddo !na

     ikb=0
     do na = 1, nat !
        do ih=1,nh(ityp(na))
           ikb=ikb+1
           do ig=1,nG_max
              iGk= list_iGk(iG)
              if (iGk==0) exit

              if (upf(nt)%has_so.and.lspinorb) then

                 do spol1 = 1,npol !spin
                    do spol2 = 1,npol !spin
                       dvnl_psi(iG, ibnd,spol1) = dvnl_psi(iG,ibnd,spol1)&
                            + DKB (ih, ih,spol1,spol2,ityp(na) ) * projec_d(na,ih,spol2) * vkb(iG,ikb)
                    end do !spol2
                 end do !spol1

              else ! no SO

                 do spol1 = 1,npol !spin

                    dvnl_psi(iG, ibnd,spol1) = dvnl_psi(iG,ibnd,spol1)&
                         + DKB (ih, ih,spol1,spol1,ityp(na) )* projec_d(na,ih,spol1) * vkb(iG,ikb)

                 end do
              end if

           enddo !ig
        enddo!ih
     enddo! na

  enddo !ibnd

  return
end SUBROUTINe multiply_psi_by_vKB
