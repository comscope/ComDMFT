      subroutine gety(gb,yr,yi,pnt_,nbas,ind_y)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer :: ind_y(2*limlb,natom),nbas,igb,ii,isort,lm,n,l,
     &	         iatom,ind,ind1,ipp,ih,is1,limj,ll,ll0,i,j,mj,lmj,
     &           mjs,lmm,i11,iii,nbndt,m,ist
      real*8 :: gb(3,nbas/nrel),gki(3),yr(max_bs,nbas),yi(max_bs,nbas),
     &	        pnt_(3),pi4om,pi2,pi2a,g2i,gi,ss,ssd,ssr,ssi,sdr,sdi,sc,
     &          co,si,sq2,c2,constb,c05,eni,c2i,anormi,cony,conb,fr,fii,
     &          ffr,ffi,cory,coiy,corb,coib,cc,ccc,bdr,bdi,drr,di,sinus,
     &          sr1,si1,sdr1,sdi1
	real*8, allocatable :: ylm(:),bes(:),besd(:)
	allocate (ylm(limlb),bes(0:maxb+1),besd(0:maxb+1))
	nbndt=nbas/nrel
      pi4om=4.d0*pi/sqrt(amega)
      pi2=2.d0*pi
      pi2a=pi2/par
      if(irel.ne.2) then
        do igb=1,nbas
          g2i=0.d0
          do ii=1,3
            gki(ii)=(pnt_(ii)+gb(ii,igb))*pi2a
            g2i=g2i+gki(ii)**2
          enddo
          gi=sqrt(g2i)
          if(gi.lt.1.d-4) then
            ylm(1)=1.d0/sqrt(4.d0*pi)
            if(maxb.gt.0) call dcopy((maxb+1)**2-1,0.d0,0,ylm(2),1)
          else
            call sphharm(gki(1)/gi,gki(2)/gi,gki(3)/gi,maxb,ylm)
          endif
          do isort=1,nsort
            call bessr(gi,smt(isort),bes,besd,lmb(isort))
            lm=0
            n=(lmb(isort)+1)**2
            do l=0,lmb(isort)
              ss=pi4om*bes(l)
              ssd=pi4om*besd(l)
              ssr=ss*rr(l)
              ssi=ss*ri(l)
              sdr=ssd*rr(l)
              sdi=ssd*ri(l)
              do m=-l,l
                lm=lm+1
                sr1=ssr*ylm(lm)
                si1=ssi*ylm(lm)
                sdr1=sdr*ylm(lm)
                sdi1=sdi*ylm(lm)
                do iatom=1,natom
                  if(is(iatom).eq.isort) then
                    ind=ind_y(lm,iatom)
                    sc=par*dot_product(gki,tau(:,iatom))
                    co=cos(sc)
                    si=sin(sc)
                    yr(ind,igb)=sr1*co-si1*si
                    yi(ind,igb)=sr1*si+si1*co
                    ind1=ind_y(n+lm,iatom)
                    yr(ind1,igb)=sdr1*co-sdi1*si
                    yi(ind1,igb)=sdr1*si+sdi1*co
                  endif    !!! over iatom=isort
                enddo    !!! over iatom=isort
              enddo   !!! over m
            enddo   !!! over l
          enddo   !!! over isort
        enddo  !!! over igb
      else if(irel.eq.2) then
        sq2=sqrt(2.d0)
        c2=clight*clight
        constb=-pi4om*c2
        c05=0.5d0*clight
        do igb=1,nbndt
          gki=pi2a*(pnt_+gb(:,igb))
          g2i=dot_product(gki,gki)
          gi=sqrt(g2i)
          eni=c05*(sqrt(c2+(2.d0*gi)**2)-clight)
          c2i=c2+eni
          anormi=c2i/sqrt(c2i*c2i+(clight*gi)**2)
          cony=pi4om*anormi
          conb=constb*gi*anormi/c2i
          if(gi.lt.1.d-4) then
            ylm(1)=1.d0/sqrt(4.d0*pi)
            if(maxb.gt.0) call dcopy((maxb+1)**2-1,0.d0,0,ylm(2),1)
          else
            call sphharm(gki(1)/gi,gki(2)/gi,gki(3)/gi,maxb,ylm)
          endif
          do isort=1,nsort
            call bessr(gi,smt(isort),bes,besd,lmb(isort)+1)
            n=2*(lmb(isort)+1)**2
            do ipp=1,2     !! over s
              ih=(ipp-1)*nbndt+igb
              if(ipp.eq.1) is1=-1
              if(ipp.eq.2) is1=1
              limj=0
              do l=0,lmb(isort)
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
                      cory=cony*ffr
                      coiy=cony*ffi
                      iii=(-1)**((i-1)/2)
                      corb=conb*ffr*iii
                      coib=conb*ffi*iii
                      cc=bes(l)
                      ccc=bes(l+i)
                      bdr=ccc*corb
                      bdi=ccc*coib
                      drr=cc*cory
                      di=cc*coiy
                    else
                      drr=pi4om*bes(l)*ffr
                      di=pi4om*bes(l)*ffi
                      bdr=pi4om*besd(l)*ffr
                      bdi=pi4om*besd(l)*ffi
                    endif
                    do iatom=1,natom
                      if(is(iatom).eq.isort) then
                        sc=par*dot_product(gki,tau(:,iatom))
                        co=cos(sc)
                        sinus=sin(sc)
                        ind=ind_y(limj,iatom)
                        yr(ind,ih)=drr*co-di*sinus
                        yi(ind,ih)=drr*sinus+di*co
                        ind1=ind_y(n+limj,iatom)
                        yr(ind1,ih)=bdr*co-bdi*sinus
                        yi(ind1,ih)=bdr*sinus+bdi*co
                      endif    !!! over iatom=isort
                    enddo    !!! over iatom=isort
                  enddo   !!! over mj
                enddo   !!! over i
              enddo   !!! over l
            enddo   !!! over ipp
          enddo   !!! over isort
        enddo   !!! over igb
      endif     !!! over irel
	deallocate (ylm,bes,besd)
      end
