      subroutine check_band_states(zb,ev,pnt_,gb,nbnd,nbas,ac)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: nbnd,nbas
      real*8, intent(in) :: pnt_(3),gb(3,nplwbas)
      complex*16, intent(in) :: zb(nfun,nbnd),ev(nbasmpw,nbnd)
      real*8, intent(out) :: ac(3,nbnd)
      integer :: nbndt,igb,isort,iatom,lm,ii,l,m,k,k1,ibnd,n,ipp,ih,is1,
     &           limj,ll,ll0,ist,i,j,mj,lmj,mjs,lmm,i11,lget,li
      real*8 :: pi4om,pi2,pi2a,gki(3),gi,sc,cs,sn,fr,sq2,fii,ffr,
     &          ffi,drr,di,bdr,bdi
      complex*16 :: cc,ccd,cc1,ccd1,ccc
	real*8, allocatable :: ylm(:),bes(:),besd(:)
	complex*16, allocatable :: mtl(:,:,:,:),itl(:,:,:,:),itl0(:,:,:,:)
	allocate(ylm((maxylm+1)**2),bes(0:maxylm+1),besd(0:maxylm+1))
	allocate(itl(nrel*(maxylm+1)**2,nsort,2,nbnd))
	allocate(itl0(nrel*(maxylm+1)**2,nsort,2,nbasmpw))
	nbndt=nbas/nrel
      pi4om=4.d0*pi/sqrt(amega)
      pi2=2.d0*pi
      pi2a=pi2/par
      sq2=sqrt(2.d0)
      itl0=(0.d0,0.d0)
      do igb=1,nbndt
        gki=(pnt_(:)+gb(:,igb))*pi2a
        gi=sqrt(dot_product(gki,gki))
        if(gi.lt.1.d-4) then
          ylm(1)=1.d0/sqrt(4.d0*pi)
          if(maxb.gt.0) call dcopy((maxylm+1)**2-1,0.d0,0,ylm(2),1)
        else
          call sphharm(gki(1)/gi,gki(2)/gi,gki(3)/gi,maxylm,ylm)
        endif
        do isort=1,nsort
          iatom=iat_1(isort)
          call bessr(gi,smt(isort),bes,besd,maxylm)
          sc=par*dot_product(gki,tau(:,iatom))
          cs=cos(sc)
          sn=sin(sc)
          ccc=dcmplx(cs,sn)
          if(irel<=1) then
            lm=0
            do l=0,maxylm
              cc=pi4om*bes(l)*(0.d0,1.d0)**l
              ccd=pi4om*besd(l)*(0.d0,1.d0)**l
              do m=-l,l
                lm=lm+1
                cc1=cc*ylm(lm)
                ccd1=ccd*ylm(lm)
                itl0(lm,isort,1,igb)=cc1*ccc
                itl0(lm,isort,2,igb)=ccd1*ccc
              enddo   !!! over m
            enddo   !!! over l
          else if(irel==2) then
            do ipp=1,2     !! over s
              ih=(ipp-1)*nbndt+igb
              if(ipp.eq.1) is1=-1
              if(ipp.eq.2) is1=1
              limj=0
              do l=0,maxylm
                ll=l*(l+2)
                ll0=l*(l+1)+1
                ist=-1
                if(l.eq.0) ist=1
                do i=ist,1,2        !! over i
                  ii=0
                  if(i.eq.-1) ii=1
                  if(i.eq.1) ii=2
                  j=l+l+i
                  do mj=-j,j,2        !! =2*mj in (III-61-1)
                    limj=limj+1
                    lmj=ll+(mj+3)/2
                    mjs=(mj-is1)/2     !! =mj-s in (III-61-1)
                    if(iabs(mjs)>l) cycle
                    lm=ll0+mjs
                    lmm=ll0-mjs
                    if(mjs.lt.0) then
                      fr=ckg2(ii,ipp,lmj)*ylm(lmm)/sq2
                      fii=ckg2(ii,ipp,lmj)*ylm(lm)/sq2
                    else if(mjs.gt.0) then
                      i11=(-1)**mjs
                      fr=ckg2(ii,ipp,lmj)*ylm(lm)*i11/sq2
                      fii=-ckg2(ii,ipp,lmj)*ylm(lmm)*i11/sq2
                    else if(mjs.eq.0) then
                      fr=ckg2(ii,ipp,lmj)*ylm(lm)
                      fii=0.d0
                    endif
                    ffr=fr*rr(l)-fii*ri(l)
                    ffi=fr*ri(l)+fii*rr(l)
                    if(rel_interst) then
                      drr=0.0d0
                      di=0.0d0
                      bdr=0.0d0
                      bdi=0.0d0
                    else
                      drr=pi4om*bes(l)*ffr
                      di=pi4om*bes(l)*ffi
                      bdr=pi4om*besd(l)*ffr
                      bdi=pi4om*besd(l)*ffi
                    endif
                    itl0(limj,isort,1,ih)=dcmplx(drr*cs-di*sn,
     &                                           drr*sn+di*cs)
                    itl0(limj,isort,2,ih)=dcmplx(bdr*cs-bdi*sn,
     &                                           bdr*sn+bdi*cs)
                  enddo   !!! over mj
                enddo   !!! over i
              enddo   !!! over l
            enddo   !!! over ipp
          endif
        enddo   !!! over isort
      enddo  !!! over igb
      n=nrel*(maxylm+1)**2*nsort*2
	call zgemm('n','n',n,nbnd,nbas,(1.d0,0.d0),itl0,n,ev,nbasmpw,
     &	         (0.d0,0.d0),itl,n)
	deallocate (itl0,bes,besd,ylm)
	allocate(mtl(nrel*(maxylm+1)**2,nsort,2,nbnd))
	mtl=(0.d0,0.d0)
	do isort=1,nsort
	  iatom=iat_1(isort)
	  do lm=1,nrel*(lmb(isort)+1)**2
	    if(irel<=1) l=lget(lm)
	    if(irel==2) call getlimj(lm,l,i,mj,li,0)
	    ii=n_e_apw(l,isort)  !! augmentation energy
          k=io_lem(iatom)+indbasa(1,ii,lm,isort)-1
          k1=io_lem(iatom)+indbasa(2,ii,lm,isort)-1
          do ibnd=1,nbnd
            mtl(lm,isort,1,ibnd)=zb(k,ibnd)
            mtl(lm,isort,2,ibnd)=zb(k1,ibnd)
          enddo
	  enddo
	enddo
	do ibnd=1,nbnd
	  do isort=1,nsort
	    do lm=1,nrel*(maxylm+1)**2
	      ac(1,ibnd)=ac(1,ibnd)+abs(mtl(lm,isort,1,ibnd)
     &	                             -itl(lm,isort,1,ibnd))
	      ac(2,ibnd)=ac(2,ibnd)+abs(mtl(lm,isort,2,ibnd)
     &	                             -itl(lm,isort,2,ibnd))
          enddo
        enddo
      enddo
	deallocate(mtl,itl)
      end
