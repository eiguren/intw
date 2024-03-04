
SUBROUTINe multiply_psi_by_dvKB( k_point, q_point, nspin, nbands,ng_max, list_iGk, list_iGkq,psi, dvnl_psi)

  USE intw_useful_constants, only: pi,tpi,fpi, cmplx_0, cmplx_i
  USE intw_reading, only : tpiba, alat, bg, nat, lspinorb
  USE intw_reading,    ONLY : ntyp, volume0, ityp
  USE intw_fft, ONLY : gvec_cart
  USE kinds, only : dp
  USE intw_pseudo, ONLY : upf, lmaxkb, nhm, nh
  USE intw_pseudo,       ONLY : nkb, vkqb, vkb, DKB
  USE intw_reading, ONLY : noncolin
  implicit none

  integer, intent(in)        :: nspin, nG_max, list_iGk(nG_max), list_iGkq(nG_max), nbands
  real(dp), intent(in)       :: k_point(3), q_point(3)

  complex(dp), intent(inout) :: psi(nG_max, nbands,nspin), dvnl_psi(nG_max,nbands, nspin, nspin, 3*nat)

  complex(kind=dp)           :: Dij(nkb,nkb,nspin,nspin)
  complex(dp)                :: projec_1(nkb,3,nspin), projec_2(nkb,3,nspin)
  complex(dp)                :: Dij_projec_1(nkb,3,nspin,nspin), Dij_projec_2(nkb,3,nspin,nspin)

  real(dp)                   :: k_(3), q_(3)
  integer                    :: nt, ntj, na, naj, ih, jh, ikb, jkb, ispin, jspin, ipol, jpol, ig
  integer                    :: mode, ibnd, iGk, iGkq


  k_(:)=matmul(bg, k_point)
  q_(:)=matmul(bg, q_point)

  ! build D matrix with ikb index
  Dij = cmplx_0
  ikb = 0
  do nt=1,ntyp
    do na=1,nat
      !
      if (ityp(na) == nt) then
        !
        do ih=1,nh(ityp(na))
          ikb = ikb + 1
          jkb = 0
          do ntj=1,ntyp
            do naj=1,nat
              if (ityp(naj) == ntj) then
                do jh=1,nh(ityp(naj))
                  jkb = jkb + 1
                  if (na==naj) then
                    if (lspinorb) then
                      Dij(ikb,jkb,:,:) = DKB(ih,jh,:,:,ityp(na))
                    else
                      do ispin=1,nspin
                        Dij(ikb,jkb,ispin,ispin) = DKB(ih,jh,1,1,ityp(na))
                      enddo
                    endif
                  endif
                enddo !jh
              endif
            enddo !naj
          enddo !ntj
        enddo !ih
        !
      endif
      !
    enddo !na
  enddo !nt




  do ibnd=1,nbands

    projec_1 = cmplx_0
    projec_2 = cmplx_0

    ! Asier: KB potentziala hurrengo eran emanik dago: sum_l |b(l)> <b(l)|
    !               beraz deribatuak bi gai dauzka:
    !               sum_l d|b(l,r)> <b(l,r)| + |b(l,r)> d<b(l,r)| ~ Fourier ~
    !               sum_l i(k+G)|b(l,G)> <b(l,G)| + |b(l,G)> <b(l,G)|
    !

    do ipol=1,3  !cart coord.
      !
      ikb = 0
      do nt=1,ntyp
        do na=1, nat
          !
          if (ityp(na) == nt) then
            !
            do ih=1,nh(ityp(na))
              !
              ikb=ikb+1
              !
              do ig=1,nG_max
              !
                iGk = list_iGk(iG)
                if (iGk==0) exit
                !
                do ispin=1,nspin !spin

                  projec_1(ikb,ipol,ispin) = projec_1(ikb,ipol,ispin) + &
                    conjg(vkb(iG,ikb)) * psi(iG,ibnd,ispin)

                  projec_2(ikb,ipol,ispin) =  projec_2(ikb,ipol,ispin) + &
                    conjg(vkb(iG,ikb)) * psi(iG,ibnd,ispin) * &
                    tpiba * cmplx_i * ( k_(ipol) + gvec_cart(ipol,iGk) )

                enddo !ispin
                !
              enddo !ig
              !
            enddo !ih
            !
          end if
          !
        enddo !na
      enddo ! nt
      !
    enddo !ipol


    ! multiplay the projections <\beta_j|\psi_n> by the matrix Dij
    Dij_projec_1 = cmplx_0
    Dij_projec_2 = cmplx_0
    do ipol=1,3
      !
      if (lspinorb) then
        do ispin=1,nspin
          do jspin=1,nspin
            Dij_projec_1(:,ipol,ispin,jspin) = matmul( Dij(:,:,ispin,jspin), projec_1(:,ipol,jspin) )
            Dij_projec_2(:,ipol,ispin,jspin) = matmul( Dij(:,:,ispin,jspin), projec_2(:,ipol,jspin) )
          enddo !jspin
        enddo !ispin
      else
        do ispin=1,nspin
          Dij_projec_1(:,ipol,ispin,ispin) = matmul( Dij(:,:,ispin,ispin), projec_1(:,ipol,ispin) )
          Dij_projec_2(:,ipol,ispin,ispin) = matmul( Dij(:,:,ispin,ispin), projec_2(:,ipol,ispin) )
        enddo !ispin
      endif
      !
    enddo !ipol


    do ipol=1,3  !cart coord.
      !
      ikb = 0
      do nt=1,ntyp
        do na =1,nat
          !
          if (ityp(na) == nt) then
            !
            mode = (na-1)*3 + ipol
            !
            do ih=1,nh(ityp(na))
              !
              ikb = ikb + 1
              !
              do ig=1,nG_max
                !
                iGkq = list_iGkq(iG)
                if (iGkq==0) exit
                !
                if (upf(nt)%has_so.and.lspinorb) then
                  !
                  do ispin=1,nspin !spin
                    do jspin=1,nspin !spin

                      dvnl_psi(iG,ibnd,ispin,jspin,mode) = dvnl_psi(iG,ibnd,ispin,jspin,mode) + &
                           Dij_projec_2(ikb,ipol,ispin,jspin) * vkqb(iG,ikb)
                      dvnl_psi(iG,ibnd,ispin,jspin,mode) = dvnl_psi(iG,ibnd,ispin,jspin,mode) - &
                           Dij_projec_1(ikb,ipol,ispin,jspin) * vkqb(iG,ikb) * &
                           ( tpiba * cmplx_i )*( k_(ipol) +  q_(ipol) + gvec_cart(ipol,iGkq) )

                    end do !jspin
                  end do !ispin
                  !
                else ! no SO
                  !
                  do ispin = 1,nspin !spin

                     dvnl_psi(iG,ibnd,ispin,ispin,mode) = dvnl_psi(iG,ibnd,ispin,ispin,mode) + &
                          Dij_projec_2(ikb,ipol,ispin,ispin) * vkqb(iG,ikb)
                     dvnl_psi(iG,ibnd,ispin,ispin,mode) = dvnl_psi(iG,ibnd,ispin,ispin,mode) - &
                          Dij_projec_1(ikb,ipol,ispin,ispin) * vkqb(iG,ikb)* &
                          ( tpiba * cmplx_i )*( k_(ipol) +  q_(ipol) + gvec_cart(ipol,iGkq) )

                  end do !ispin
                  !
                end if

              enddo !ig
              !
            enddo !ih
            !
          end if
          !
        enddo !na
      enddo !nt
      !
    enddo !ipol

  enddo !ibnd

  return
end SUBROUTINe multiply_psi_by_dvKB
