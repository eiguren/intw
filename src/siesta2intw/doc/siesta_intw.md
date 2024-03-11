# SIESTA-intw interface.

`siesta-master` direktorioan`Src/write_intw.f90` bezala sartuta dagoen moduloaren azalpen labur bat da hau. Siestako exekutatuko den direktorian hurrengo fitxategia dagoela suposatzen du `s2intw.in`, izen horrekin eta itxura hau izan behar du.

```bash
&inputpp
prefix="FePt_soc"
nsym=1
symfile="sym.dat"
nk1=3, nk2=3, nk3=3,
mesh_dir="./"
/
```

`nk1,nk2,nk3` hurrengo batean kristal batean kalkulatzeko puntu irreduzibleak kalkulatzeko sartuta dago. Kalkulatzen du baina siestako inputean sartzea zen ideia, oraindik ez dagoena eginda. Oraingoz, eskuz sartu behar dira kalkulatu nahi diren puntuak. Adibidez,



```bash 
DM.UseSaveDM true
Denchar.PlotWaveFunctions true
WaveFuncKPointsScale ReciprocalLatticeVectors
%block WaveFuncKPoints
0.000000    0.000000    0.000000
0.000000    0.000000    0.125000
0.000000    0.000000    0.25000
%endblock WaveFuncKPoints
```

blokea gehituz siestako input fitxategian. Goikoak hiru puntu kalkulatzen ditu baina guk puntu irreduzibleak kalkulatu ondoren siestari horiek kalkulatzeko esan geniezaioke.

Siestako oinarriko elementuak $\phi_{i,l,m}(r)$ sare erradial batean emanda daude eta programa honek egiten duena da $\Psi_{k,n}(r)=\sum_{R,j} e^{ikR} c_j \phi_{j,l,m}(r)$ bezala datorren uhin bat transformatu uhin lauetako oinarrira. Horretarako,
$$
e^{i\mathbf{k}\mathbf{r}} = 4 \pi \sum_{l=0}^{\infty}\sum_{m=-l}^{l} i^{l} j_l(k r) Y_{l,m}(\hat{k})Y^*_{l,m}(\hat{r})
$$
erabiliko dugu. Izatez, oinarriko elementuak $lm$ bikote bat finkatuta izanik, uhin lau bat berarekin proiektatzean (espazio errealean) integral erradial bat geratzen da eta bietako harmoniko esferiko bat desagertzen da, ortogonala izateagatik beste guztiekin.
$$
\tilde{\phi}_{i,l,m}(k)=\int_0^\infty dr \ r^2 j_{l}(kr) \phi_{i,l,m}(r).
$$
Izatez, $r$ integrala bakarrik oinarriko elementuentzat definituta dagoen $r_c$ cutoff bateraino bakarrik kalkulatzen da, oraingoz `Simpon erregela` zikin bat erabiliz.  Izatez, $k+G$ bekore batentzak kalkulatu behar izango da goikoa $k$ bektore bat lotuta daukan uhin batentzat, kodea horrela pentsatuz eginda dago.



```fortran
!-------------------------------------------------------------------------------    
  FUNCTION g_of_phi(is,io,gg)
    !-------------------------------------------------------------------------------
    !-Runge-Kutta4 (=Simpson) integrazioa \int^rc_0 r**2 * bessel(l,k*r) * phi(is,io,r) dr 
    implicit none
    integer, intent(in) :: is, io ! specie eta orbital indizea
    real(kind=dp), intent(in)    :: gg     ! g+k gektorearen moduloa
    real(kind=dp) :: g_of_phi
    !-lokal
    integer :: ir
    integer,parameter   :: nr=500 ! integratzeko 500 puntu [0,rc]. Ez dakit merezi duen orokortzea 
    real(kind=dp)       :: r, k1, k2, k4, phi,dphidr,rh,h
    integer :: l

    l =  species(is)%orb_l(io)

    h = rcut(is,io)/nr

    g_of_phi= 0.0_dp
    r       = 0.0_dp
    ! ez da efizientea baina bai garbia (puntu erdiak 2 aldiz kalkulatze dira)
    do ir=1, nr
       r =(ir-1)*h ! (nr-1)*rc/nr=rc-h
       rh=r
       call rphiatm(is,io,rh,phi,dphidr)
       k1 = phi * rh**2 * sphb(l,rh*gg)
       rh = r+0.5_dp*h 
       call rphiatm(is,io,rh,phi,dphidr)
       k2 = phi * rh**2 * sphb(l,rh*gg)
       rh = r+h  
       call rphiatm(is,io,rh,phi,dphidr)
       k4 = phi * rh**2 * sphb(l,rh*gg)
       r = r + h
       g_of_phi =  g_of_phi + h * (k1 + 4*k2 +k4) / 6.0_dp    ! 
    end do ! ir

```

Goikoak $k$ sare batentzat kalkulatzen du Fourier transformatua. Goikoa elementu guztientzat egin eta gorde hurrengo subrutinarekin,

```fortran 
 !-------------------------------------------------------------------------------
  SUBROUTINE build_basis_g_table(nq)
    !-------------------------------------------------------------------------------
    USE siesta_options,  only: g2cut
    use atm_types, only: species_info, species, nspecies
    implicit none 
    integer, intent(in):: nq
    integer :: is, io, iq, l
    real(kind=dp) :: gg(nq),q,f

    !ASIER: Hemen m ezberdinetako Fourier Tr. berdinak dira
    !       eta dena den kalkulatzen da errazagoa izateagatik.
    phi_g_table=0.0_dp

    do is=1,nspecies

       do io=1, species(is)%norbs
          !print*, " -- io --",io,species(is)%orb_l(io),species(is)%orb_m(io)
          !print*, "rc", rcut(is,io)
          l =  species(is)%orb_l(io)
          do iq=1,nq
             gg_list(iq)=(iq-1)*sqrt(g2cut)/nq
             phi_g_table(is,io,iq) = g_of_phi(is,io,gg_list(iq))
          enddo
          call spline_dp ( gg_list, phi_g_table(is,io,:), nq, phi_g_y2(is,io,:) ) 
       end do
    end do

  END SUBROUTINE build_basis_g_table

```

Goian `phi_g_table` aldagaian gordentzen da lista jakin bateko $q$ puntuen Fourier osagaiak eta spline interpolazioa egiteko bigarren deribatuak ere gordetzen dira`phi_g_y2`, oso azkar beste edozein $k+G$ balioentzat kalkulatu ahal izateko hurrengoarekin,

```fortran
  !-------------------------------------------------------------------------------    
  FUNCTION spline_g_of_phi(is,io,q)
    !-------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: is, io ! specie eta orbital indizea
    real(kind=dp), intent(in)    :: q
    real(kind=dp) :: spline_g_of_phi
    !
    integer :: nq

    nq=size(gg_list) 
    if (q> maxval(gg_list)) then 
       spline_g_of_phi=0.0_dp
       return
    end if
    call splint_dp (gg_list, phi_g_table(is,io,:),phi_g_y2(is,io,:),&
          nq, q, spline_g_of_phi ) 

  end FUNCTION spline_g_of_phi

```

###  

Beraz, orain artekoarekin printzipioz oinarriko elementuen $\int d^r e^{i (\mathbf{(k+G)r})} \phi_{l,m}(\mathbf{r})$ integralak eginda daude **specie**bakoitzerako eta atomo bakoitza zein espzieri dagokion jakinin, $e^{ikR_a}$ fase batekin biderkatu behar da Fourier garapenean desplazatutako atomo bati egokitzeko.

```fortran
 !-------------------------------------------------------------------------------
  subroutine write_wfc()
    !-------------------------------------------------------------------------------

    USE precision
    use siesta_geom
    USE siesta_options,  only: g2cut
    USE units, only : pi
    USE writewave, only: nwk, wfk
    USE s2intw
    use m_spin, only: NonCol, SpOrb, nspin
    use atomlist, only: indxuo, iaorb, no_u
    use atmparams,       only: lmaxd
    use spher_harm

    implicit none
    integer :: ik, ibnd, jbnd, is, npol, ig
    character(256) ::  wfc_file, datafile
    integer :: io_unit
    real(dp) :: q 
    integer :: iat, io, ia, iq
    real(dp) :: ain(3,3), kpoint(3), gv(3),kr, kpg, kg_unit(1:3)
    integer :: ipol 
    complex(dp) :: phase
    integer :: l,m, ind

    real(kind=dp), allocatable :: rylm(:), grylm(:,:) , orbital_gk(:,:,:,:)
    complex(kind=dp) :: ss 

    npol=1
    if (nspin>=4) npol=2 !ASIER KONTUZ!!!

    allocate(evc(npol*npwx, nbnd, nwk))
    evc=(0.0_dp,0.0_dp)
 
   allocate(rylm( 0:(lmaxd+1)**2-1))
    allocate(grylm(1:3, 0:(lmaxd+1)**2-1))

    datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"kpoints.dat")
    io_unit=find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="write",form="formatted")
    do ik=1,nwk
     write(io_unit,"(3f18.10)") kpoint_wfsx(1:3,ik) / (2*pi/alat)!matmul(ainv(bg*2*pi/alat),kpoint_wfsx(1:3,ik))
    end do
    close(unit=io_unit) 
 
    !lehenik transformatu esferikoa egin espezie bakoitzarentzat. Gero atomoa mugitzean xa(:,ia)
    !fase bat besterik ez du aldatuko
    allocate(orbital_gk(nspecies, maxval(species(:)%norbs ) ,npwx,nwk))

    !Behean specie bakoitzarentzat orbitalentzat Fourier egin jatorrian egongo balira bezala
    do ik=1,nwk
       print*,"ik = ..", ik, ngk(ik)
       kpoint(1:3) = kpoint_wfsx(1:3,ik) 
       do iG=1,ngk(ik)
          kpg = sqrt(sum( (  kpoint(:)+ gvec_cart(:,list_iG(ig,ik)) )**2 )) ! |k+G| moduloa.
          kg_unit=  (kpoint(:) + gvec_cart(:,list_iG(ig,ik)))/kpg           ! (k+G)/|k+G| unitate bektorea 
          call rlylm( lmaxd, kg_unit , rylm, grylm )                        ! lamx bateraiko kalkulatzen ditu Ylm guztiak
                                                                            ! do l=0,lmax
                                                                            !  do m=-l,l
                do is=1,nspecies                                            ! indize bateratua: l**2 + l + m + 1 
                  do io=1,species(is)%norbs
                   l=species(is)%orb_l(io)    !orbital horren l
                   m=species(is)%orb_m(io)    !orbital horren m
                   ind = l**2 + l + m + 1     ! hau indize bateratua lm  
                   phase = (0.0_dp,1.0_dp)**l ! i^l faktorea exp(kr) = 4pi \sum_{lm} i^l j_l(kr) Y_lm(^k) Y_lm(^r)  
                   orbital_gk(is,io,iG,ik) = 4*pi * phase *  spline_g_of_phi(is,io, kpg) * rylm(ind) 
                end do !io
              end do !is
       end do !iG
    end do !ik

    ! Orbitalen Fourier aprobetxatu uhien Fourier T. kalkulatzeko
    do ik=1,nwk
       kpoint(1:3) = kpoint_wfsx(1:3,ik) ! 2*pi/alat biderkatuta atala dauka
       do iG=1,ngk(ik)
          do ibnd=1, nbnd 
             do ipol=1,npol
                do io=1,no_u
                   ia = iaorb(io) ! zein atomori dagokion orbital hori
                   is = isa(ia)   ! zein den atomo speciea 
                   kr = sum(kpoint(:)*xa(:,ia)) ! xa-k alat biderkatuta dauka alat, eta kpoit 2pi/alat
                   phase = cmplx(cos(kr),-sin(kr))
                   evc((ipol-1)*ngk(ik)+ig, ibnd,ik) = phase * orbital_gk(is,io,iG,ik) * wf(io,ipol,ibnd, ik)  
                end do !io
             end do !ipol
          end do !ibnd
       end do !iG
    end do !ik 
    !ain=ainv(at)
    !do ia=1,na_u 
    !  print"(a,i4,3f12.6)","ia", ia, matmul(ain,xa(:,ia))/alat
    !enddo

    do ik=1,nwk
       write(wfc_file,100) ik
       datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//trim(wfc_file))
       io_unit=find_free_unit()
       open(unit=io_unit,file=datafile,status="unknown", action="write",form="unformatted")
       write(unit=io_unit) ngk(ik) 
       write(unit=io_unit) list_iG(1:ngk(ik), ik)
       write(unit=io_unit) ( et(ibnd,ik), ibnd=1,nbnd ) 
       do ibnd=1,nbnd
         write(unit=io_unit) ((evc((ipol-1)*ngk(ik)+ig, ibnd,ik), ig=1,ngk(ik)), ipol=1,npol )
       end do
       close(unit=io_unit) 
    end do !ik

    !ortogonltasun testa 
    do ik=1,1!nwk
         do ibnd=1, 4!nbnd
          do jbnd=1, 4!nbnd
               ss=(0.0,0.0)
              do ig=1, ngk(ik)
               ss=ss+ conjg(evc(ig, ibnd,ik)) * evc(ig, jbnd,ik) 
              end do
              write(*,"(a,2i3,2f14.8)")"ibnd jbnd", ibnd, jbnd, ss
          end do
        end do
   end do !ik
stop

  deallocate(orbital_gk)
100 format('wfc'I5.5'.dat')

  end subroutine write_wfc

```

```fortran
!------------------------------------------------------------------------------- 
  recursive function sphb(n,x) result(sphb_)
    !-------------------------------------------------------------------------------
    implicit none
    ! ASIER:: Errekursiboki definitu dut (kontuz)
    ! behar bada siestan badago beste zerbait
    real(kind=dp), intent(in) :: x
    integer, intent(in) :: n
    real(kind=dp) :: sphb_
    if (n==0) then 
       if (abs(x)<epsilon(x)) then
          sphb_ = 1.0_dp
       else
          sphb_ = sin(x)/x
       end if
    else if (n==1) then
       if (abs(x)<epsilon(x)) then
          sphb_ = 0.0_dp
       else
          sphb_ = sin(x)/x**2 - cos(x)/x
       end if
    else if (n>1) then
       if (abs(x)<epsilon(x)) then
          sphb_ = 0.0_dp
       else
          sphb_ = - sphb(n-2,x) + (2.0_dp* n - 1.0_dp ) * sphb(n-1,x)/x 
       end if

    end if

  end function sphb

```

###    `write_intw_file` subrutina.

Honek, intw-k irakurriko dituen aldagai garrantzitsuenak idazten ditu eta beste subrutina batzuei deitu, hala nola fft sarea kalkulatzeko etb. Oso zikin dago oraindik.

```fortran
  !-------------------------------------------------------------------------------
  subroutine write_intw_file()
  !-------------------------------------------------------------------------------
    use FDF
    use datuak
    use s2intw

    use siesta_geom
    use atomlist, only: indxuo, iaorb, lastkb, lasto, datm, no_l, &
         iphkb, no_u, no_s, iza, iphorb, rmaxo, indxua
    USE siesta_options,  only: g2cut
    use units,     only : Ang, eV, deg, pi
    use m_spin, only: NonCol, SpOrb, nspin
    use atm_types
    use atmfuncs
    USE writewave, only: nwk, wfk

    implicit none

    integer :: ix, ia, i, j
    integer :: io_unit
    character(len=256) :: datafile

    real(kind=dp) :: b(3)
    integer :: nri(3)

    real(dp):: gg
    integer :: is,io,iq

    call read_input()

    io_unit= find_free_unit()
    datafile=trim(trim(mesh_dir)//trim(prefix)//&
         ".save.intw/"//"crystal.dat")
    open(unit=io_unit,file=trim(datafile),&
         status="unknown", action="write",form="formatted")

    alat = fdf_get( 'LatticeConstant', 0.0_dp, 'Bohr' ) 
    write(unit=io_unit,fmt=*)"ALAT"
    write(unit=io_unit,fmt="(f16.10)")alat
    write(unit=io_unit,fmt=*)"AT"
    at=ucell/alat
    do i=1,3
       write(unit=io_unit, fmt="(3f12.6)")(at(i,j), j=1,3)
    enddo
    call reclat( at, bg, -1 )
    write(unit=io_unit,fmt=*)"BG"
    do i=1,3
       write(unit=io_unit, fmt="(3f12.6)")(bg(i,j), j=1,3)
    enddo
    !g2cut = fdf_get('MeshCutoff',300._dp,'Ry')
    do i=1,3
       b = at(:,i)*alat
       nri(i)= nint(sqrt(G2cut)* sqrt(sum( b*b ))/pi) + 1  !
       !bakoitia bada bikoiti bihurtu
       if (nri(i)/2*2/=nri(i)) nri(i) = nri(i) +1
    end do
    write(unit=io_unit,fmt=*)"FFT GRID"
    write(unit=io_unit,fmt=*)nri(:)
    write(unit=io_unit,fmt=*)"ECUTWFC"
    write(unit=io_unit,fmt="(f12.6)")real(nint(g2cut/4))
    write(unit=io_unit,fmt=*)"ECUTRHO"
    write(unit=io_unit,fmt="(f12.6)")g2cut
    write(unit=io_unit,fmt=*)"NONCOLIN(kontuz)"
    !ASIER kontuz hemen. Ikusi nola irakurtzen diren siestako wfc 
    ! hor ezartzen dira noncoll deitutako beste balio batzuk
    noncol=.false.
    if (nspin>=4) noncol=.true.
    write(unit=io_unit,fmt=*) NonCol
    write(unit=io_unit,fmt=*)"LSPINORB(kontuz)"
    write(unit=io_unit,fmt=*) SpOrb
    write(unit=io_unit,fmt=*)"SPINORB_MAG (KONTUZ, hau aldatu da)"
    if (nspin==8) then 
       write(unit=io_unit,fmt=*) .true.
    else 
       write(unit=io_unit,fmt=*) .false. 
    end if

    write(unit=io_unit,fmt=*)"NAT"
    write(unit=io_unit,fmt=*)na_u
    write(unit=io_unit,fmt=*)"NTYP"
    write(unit=io_unit,fmt=*)nspecies
    write(unit=io_unit,fmt=*)"ATOM_LABELS, MASS AND PP_FILE (1:NTYP)"
    do i=1,nspecies
       write(unit=io_unit,fmt="(a,f16.5,x,a)")trim(species(i)%symbol), species(i)%mass, trim(species(i)%label)
    enddo

    write(unit=io_unit,fmt=*)"POSITIONS (1:NAT)"
    do i=1,na_u
       write(unit=io_unit,fmt="(a,i4,3f16.8)") trim(species(isa(i))%symbol), isa(i), xa(:,i)
    end do

    write(unit=io_unit,fmt=*)"NSYM"
    write(unit=io_unit,fmt=*)nsym
    do isym=1,nsym 
       write(unit=io_unit,fmt="(i8)"), isym 
       do i=1, 3
          write(unit=io_unit, fmt="(3i8)"), (s(i,j,isym), j=1,3) 
       enddo
       write(unit=io_unit,fmt="(100f16.10)")  &
            real(ftau(1:3,isym),dp)/nri(1:3)
       write(unit=io_unit,fmt="(48i3)")t_rev (isym)  
    end do
    print*, "nspin", nspin

    write(unit=io_unit,fmt=*)"NKS"
    write(unit=io_unit,fmt=*)nwk
    !FFT sarea definitu.
    call  write_fft_information (alat, at, bg, g2cut, nri ) 
    !siestako uhinak irakurri. Oraingoz nk, nbnd guztiak allocatzen dira
    !baina aldatu beharko da behar bada RAM-ik ez badago
    call read_siesta_wfcs()

    write(unit=io_unit,fmt=*)"NGMAX"
    write(unit=io_unit,fmt=*)npwx
   
    write(unit=io_unit,fmt=*)"NBAND"
    ! siestak ez du beti banda kopuru berdina erabiltzen k guztietarako
    ! beraz, minimoa erabiliko dugu.
    nbnd=minval(number_of_wfns(:))
    write(unit=io_unit,fmt=*)nbnd
    close(unit=io_unit)

    ! hemen oinarriaren g taula allocate bat egin
    call allocate_basis_g_table(nq,g2cut)
    ! hemen k eskalar batzuentzat oinarriko elementuen Fourier egiten 
    ! da eta spline interpolatzeko bigarren deribatuak gorde.
    call build_basis_g_table(nq)
    !hemen bukaerako transformazioa egiten da.
    call write_wfc()

  end  subroutine write_intw_file
```





