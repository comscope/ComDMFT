      subroutine setdiag(gb,b,a,pnt_,nbas,ndim,yr,yi,ind_y,ar0,
     &	                 ai0,br0,bi0,ar1,ai1,br1,bi1,ispin,indg)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: nbas,indg(nbas/nrel)
      integer :: ind_y(2*limlb,natom),ndim,ispin,
     &	         iatom,isort,llmm,ll,ll1,mm,ie1,jbas,lm,l,l1,
     &           m,ie,ibas,li1,ist1,i1,jj1,mj1,li,jj,mj,jatom,
     &           jsort,n,ind,ie0,iba,ia,ib,ic,jgb0,igb0,jp,jh,ipp,ih,
     &           nbndt,ist,i,ind1,jgb,igb,ii,ispin0,ig0,ig,iatom0
      real*8 :: gb(3,nbas/nrel),yr(max_bs,nbas),yi(max_bs,nbas),
     &	        ar0(nbas,nlcloc),ai0(nbas,nlcloc),
     &          br0(nbas,nlcloc),bi0(nbas,nlcloc),ar1(nbas,nbas),
     &          ai1(nbas,nbas),br1(nbas,nbas),bi1(nbas,nbas),
     &          pnt_(3),cr(3),gki(3),gkj(3),arr,brr,smtcon,drr,di,bbr,
     &          bbi,aar,aai,ffh,pi2,pi2a,c1,gj1,gj2,gj3,gj,gj22,vhv0,
     &          vhs0,shv0,shs0,vv,vs,ss,gi1,gi2,gi3,gi,x,fac,faci,sp,
     &          sum,c05,c2,enj,c2j,anormj,const,cs,g2i,scal,amult,sr,si,
     &          cons,consi,am,abskap,sumi,con,g2j,cs1
      complex*16 :: a(ndim,ndim),b(ndim,ndim)
	real*8, allocatable :: p(:),pp(:),fc(:),gvv(:,:),gss(:,:),
     &	                     gvhv(:,:),gshs(:,:),g1(:,:,:),g2(:,:,:),
     &                       ppr(:),ppi(:),bes(:),besd(:)
	nbndt=nbas/nrel
	ar0=0.d0
	ai0=0.d0
	br0=0.d0
	bi0=0.d0
	allocate(p(0:maxb),gvv(nrel*maxb+1,nsort),
     &	       gss(nrel*maxb+1,nsort),gvhv(nrel*maxb+1,nsort),
     &         gshs(nrel*maxb+1,nsort),g1(nrel*maxb+1,nsort,nbndt),
     &         g2(nrel*maxb+1,nsort,nbndt),fc(maxb+1),pp(0:maxb),
     &         ppr(nrel*maxb+1),ppi(nrel*maxb+1),bes(0:maxb+1),
     &         besd(0:maxb+1))
c --------- ME LO(lo)-LO(lo) (spherical) between local functions -------
      do iatom=1,natom
        isort=is(iatom)
        if(irel.ne.2) then
	    ispin0=ispin
	    iatom0=iat_1(isort)
	    do ig0=1,ngroup
	      if(iatom0==ip(iatom,ig0)) then
	        ig=ig0
	        exit
	      endif
	    enddo
	    llmm=0
          do ll=0,lmb(isort)
            ll1=ll+1
            do mm=-ll,ll
	        llmm=llmm+1
              do ie1=1,ntle(ll,isort)
                if(augm(ie1,ll,isort)=='LOC') then
	            jbas=indbas(ie1,llmm,iatom)  !! the position in basis
	            lm=0
                  do l=0,lmb(isort)
                    l1=l+1
                    do m=-l,l
	                lm=lm+1
                      do ie=1,ntle(l,isort)
                        if(augm(ie,l,isort)=='LOC') then
                          ibas=indbas(ie,lm,iatom)
                          if(ll.eq.l.and.mm.eq.m) then
                            arr=0.5d0*(ffhmt(1,1,ie,ie1,l1,isort,ispin0)
     &				                 +ffhmt(1,1,ie1,ie,l1,isort,ispin0))
                            brr=ffsmt(1,1,ie,ie1,l1,isort,ispin0)
                            a(ibas,jbas)=dcmplx(arr,0.d0)
                            b(ibas,jbas)=dcmplx(brr,0.d0)
                          endif
                        endif
                      enddo
                    enddo
                  enddo
                endif
              enddo   !! over ie1
            enddo  !! over mm
          enddo
        else if(irel.eq.2) then
          li1=0
	    llmm=0
          do ll=0,lmb(isort)
            ist1=-1
            if(ll.eq.0) ist1=1
            do i1=ist1,1,2
              li1=li1+1
              jj1=ll+ll+i1
              do mj1=-jj1,jj1,2
	          llmm=llmm+1
                do ie1=1,ntle(ll,isort)
                  if(augm(ie1,ll,isort)=='LOC') then
                    jbas=indbas(ie1,llmm,iatom)
                    li=0
	              lm=0
                    do l=0,lmb(isort)
                      ist=-1
                      if(l.eq.0) ist=1
                      do i=ist,1,2
                        li=li+1
                        jj=l+l+i
                        do mj=-jj,jj,2
	                    lm=lm+1
                          do ie=1,ntle(l,isort)
                            if(augm(ie,l,isort)=='LOC') then
                              ibas=indbas(ie,lm,iatom)
                              if(li1.eq.li.and.mj.eq.mj1) then
                                arr=0.5d0*
     &                           (ffhmt(1,1,ie,ie1,li,isort,ispin)
     &  	 	                 +ffhmt(1,1,ie1,ie,li,isort,ispin))
                                brr=ffsmt(1,1,ie,ie1,li,isort,ispin)
                                a(ibas,jbas)=dcmplx(arr,0.d0)
                                b(ibas,jbas)=dcmplx(brr,0.d0)
                              endif
                            endif
                          enddo
                        enddo
                      enddo
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        endif
      enddo   !!! over iatom
c ----------- PW-LO(lo) spherical matrix elements ------------------------
      do jatom=1,natom
        jsort=is(jatom)
	  smtcon=0.5d0*smt(jsort)**2
	  if(irel.ne.2) then
 	    ispin0=ispin
	    iatom0=iat_1(jsort)
	    do ig0=1,ngroup
	      if(iatom0==ip(jatom,ig0)) then
	        ig=ig0
	        exit
	      endif
	    enddo
          n=(lmb(jsort)+1)**2
          do l=0,lmb(jsort)
            ie=n_e_apw(l,jsort)
            if(ie==0) cycle
            l1=l+1
            do m=-l,l
              lm=l*l1+m+1
              ind=ind_y(lm,jatom)
              ind1=ind_y(n+lm,jatom)
              do ie0=1,ntle(l,jsort)  !! over local orbitals
                if(augm(ie0,l,jsort)=='LOC') then
                  jbas=indbas(ie0,lm,jatom)
                  do ibas=nlcloc+1,ndim
                    iba=ibas-nlcloc
	              drr=yr(ind,iba)
	              di=-yi(ind,iba)
                    bbr=drr*ffsmt(1,1,ie,ie0,l1,jsort,ispin0)
                    bbi=di*ffsmt(1,1,ie,ie0,l1,jsort,ispin0)
                    aar=0.5d0*drr*(ffhmt(1,1,ie,ie0,l1,jsort,ispin0)
     &			                  +ffhmt(1,1,ie0,ie,l1,jsort,ispin0))
                    aai=0.5d0*di*(ffhmt(1,1,ie,ie0,l1,jsort,ispin0)
     &			 	             +ffhmt(1,1,ie0,ie,l1,jsort,ispin0))
c --------- Standard LAPW calculation for slope-typefunction --------------------
	              drr=yr(ind1,iba)
	              di=-yi(ind1,iba)
                    ffh=0.5d0*(ffhmt(2,1,ie,ie0,l1,jsort,ispin0)
     &                        +ffhmt(1,2,ie0,ie,l1,jsort,ispin0))
                    aar=aar+drr*ffh
                    bbr=bbr+drr*ffsmt(2,1,ie,ie0,l1,jsort,ispin0)
                    aai=aai+di*ffh
                    bbi=bbi+di*ffsmt(2,1,ie,ie0,l1,jsort,ispin0)
                    ar0(iba,jbas)=aar
                    br0(iba,jbas)=bbr
                    ai0(iba,jbas)=aai
                    bi0(iba,jbas)=bbi
                  enddo
                endif
              enddo   !!! over ie0
            enddo   !!! over m
          enddo   !!! over l
        else if(irel.eq.2) then
          n=2*(lmb(jsort)+1)**2
          li=0
          do l=0,lmb(jsort)
            ie=n_e_apw(l,jsort)
            if(ie==0) cycle
            ist=-1
            if(l.eq.0) ist=1
            do i=ist,1,2
              li=li+1
              jj=l+l+i
              do mj=-jj,jj,2
                call getlimj(lm,l,i,mj,li,1)
                ind=ind_y(lm,jatom)
                ind1=ind_y(n+lm,jatom)
                do ie0=1,ntle(l,jsort)  !! over local orbitals
                  if(augm(ie0,l,jsort)=='LOC') then
                    jbas=indbas(ie0,lm,jatom)
                    do ibas=nlcloc+1,ndim
                      iba=ibas-nlcloc
	                drr=yr(ind,iba)
	                di=-yi(ind,iba)
                      bbr=drr*ffsmt(1,1,ie,ie0,li,jsort,ispin)
                      bbi=di*ffsmt(1,1,ie,ie0,li,jsort,ispin)
                      aar=0.5d0*drr*(ffhmt(1,1,ie,ie0,li,jsort,ispin)
     &			                    +ffhmt(1,1,ie0,ie,li,jsort,ispin))
                      aai=0.5d0*di*(ffhmt(1,1,ie,ie0,li,jsort,ispin)
     &			                   +ffhmt(1,1,ie0,ie,li,jsort,ispin))
c --------- Standard LAPW calculation for slope-typefunction --------------------
	                drr=yr(ind1,iba)
	                di=-yi(ind1,iba)
                      ffh=0.5d0*(ffhmt(2,1,ie,ie0,li,jsort,ispin)
     &                          +ffhmt(1,2,ie0,ie,li,jsort,ispin))
                      aar=aar+drr*ffh
                      bbr=bbr+drr*ffsmt(2,1,ie,ie0,li,jsort,ispin)
                      aai=aai+di*ffh
                      bbi=bbi+di*ffsmt(2,1,ie,ie0,li,jsort,ispin)
                      ar0(iba,jbas)=aar
                      br0(iba,jbas)=bbr
                      ai0(iba,jbas)=aai
                      bi0(iba,jbas)=bbi
                    enddo
                  endif
                enddo   !!! over ie0
              enddo   !!! over mj
            enddo   !!! over i
          enddo   !!! over l
        endif
      enddo   !!! over jatom
c ----------- PW-PW spherical+interstitial matrix elements ------------------------
      pi2=2.d0*pi
      pi2a=pi2/par
      c1=amegaint/amega
      if(irel.ne.2) then
        do l1=1,maxb+1
          fc(l1)=2.d0*sqrt(pi/amega)
        enddo
        do jgb=1,nbndt
          gj1=pnt_(1)+gb(1,jgb)
          gj2=pnt_(2)+gb(2,jgb)
          gj3=pnt_(3)+gb(3,jgb)
          gj=sqrt(gj1*gj1+gj2*gj2+gj3*gj3)
          do isort=1,nsort
            call BESSR(pi2a*gj,smt(isort),BES,BESD,lmb(isort))
c ------ See V-59-3,4 --------------------------------------------
            do l1=1,lmb(isort)+1
              g1(l1,isort,jgb)=fc(l1)*bes(l1-1)
              g2(l1,isort,jgb)=fc(l1)*besd(l1-1)
            enddo
          enddo
        enddo
        do jgb=1,nbndt
          jgb0=indg(jgb)
          gj1=pnt_(1)+gb(1,jgb)
          gj2=pnt_(2)+gb(2,jgb)
          gj3=pnt_(3)+gb(3,jgb)
          gj=gj1*gj1+gj2*gj2+gj3*gj3
          gj22=pi2a**2*gj
          gj=sqrt(gj)
          do isort=1,nsort
	      smtcon=0.5d0*smt(isort)**2
            do l1=1,lmb(isort)+1
              ie=n_e_apw(l1-1,isort)
              if(ie==0) then
                gshs(l1,isort)=0.d0
                gvv(l1,isort)=0.d0
                gvhv(l1,isort)=0.d0
                gss(l1,isort)=0.d0
                cycle
              endif
              vhv0=ffhmt(1,1,ie,ie,l1,isort,ispin)
c -------- Adding correction VII-1-5 ----------------------------------
              vhs0=0.5d0*(ffhmt(1,2,ie,ie,l1,isort,ispin)
     &                   +ffhmt(2,1,ie,ie,l1,isort,ispin))
	        shv0=vhs0
              shs0=ffhmt(2,2,ie,ie,l1,isort,ispin)
              vv=ffsmt(1,1,ie,ie,l1,isort,ispin)
              vs=ffsmt(1,2,ie,ie,l1,isort,ispin)
              ss=ffsmt(2,2,ie,ie,l1,isort,ispin)
              gshs(l1,isort)=g1(l1,isort,jgb)*shv0
     &                      +g2(l1,isort,jgb)*shs0
              gvv(l1,isort)=g1(l1,isort,jgb)*vv
     &                     +g2(l1,isort,jgb)*vs
              gss(l1,isort)=g1(l1,isort,jgb)*vs
     &                     +g2(l1,isort,jgb)*ss
              gvhv(l1,isort)=g1(l1,isort,jgb)*vhv0
     &		                +g2(l1,isort,jgb)*vhs0
c ------ See V-59-5 --------------------------------------------
            enddo
          enddo
          do igb=1,nbndt
            gi1=pnt_(1)+gb(1,igb)
            gi2=pnt_(2)+gb(2,igb)
            gi3=pnt_(3)+gb(3,igb)
            gi=sqrt(gi1*gi1+gi2*gi2+gi3*gi3)
            igb0=indg(igb)
            ia=igbs(1,jgb0)-igbs(1,igb0)
            ib=igbs(2,jgb0)-igbs(2,igb0)
            ic=igbs(3,jgb0)-igbs(3,igb0)
C     *************************************************************
            p(0)=1.d0
            if(maxb.gt.0.and.gi.gt.1.d-8.and.gj.gt.1.d-8) then
              x=(gj1*gi1+gj2*gi2+gj3*gi3)/gi/gj
              p(1)=x
              do l=1,maxb-1
                p(l+1)=(dfloat(l+l+1)*x*p(l)-dfloat(l)*p(l-1))
     &			                                      /dfloat(l+1)
              enddo
            else
              do l=1,maxb
                p(l)=0.d0
              enddo
            endif
            do l=0,maxb
              p(l)=p(l)*(l+l+1.d0)
            enddo
            do isort=1,nsort
              fac=0.d0
              if(complex_ro) faci=0.d0
              do iatom=1,natom
                if(is(iatom).eq.isort) then
                  sp=(gj1-gi1)*tau(1,iatom)+(gj2-gi2)*tau(2,iatom)+
     &               (gj3-gi3)*tau(3,iatom)
                  sp=pi2*sp
                  fac=fac+cos(sp)
                  if(complex_ro) faci=faci+sin(sp)
                endif
              enddo
C     ***     O_mt       *********************************************
              sum=0.d0
              do l1=1,lmb(isort)+1
                sum=sum+p(l1-1)*(g1(l1,isort,igb)*gvv(l1,isort)+
     &                           g2(l1,isort,igb)*gss(l1,isort))
              enddo
              br1(igb,jgb)=br1(igb,jgb)+fac*sum
              if(complex_ro) bi1(igb,jgb)=bi1(igb,jgb)+faci*sum
C     ***     H_mt       *********************************************
              sum=0.d0
              do l1=1,lmb(isort)+1
                sum=sum+p(l1-1)*(g1(l1,isort,igb)*gvhv(l1,isort)+
     &                           g2(l1,isort,igb)*gshs(l1,isort))
              enddo
              ar1(igb,jgb)=ar1(igb,jgb)+fac*sum
              if(complex_ro) ai1(igb,jgb)=ai1(igb,jgb)+faci*sum
            enddo            !!! over isort
          enddo           !!! over igb
        enddo               !!! over jgb
      else if(irel.eq.2) then
        c05=0.5d0*clight
        c2=clight*clight
        fc(1)=2.d0*sqrt(pi/amega)
        do jgb=1,nbndt
          gj1=pnt_(1)+gb(1,jgb)
          gj2=pnt_(2)+gb(2,jgb)
          gj3=pnt_(3)+gb(3,jgb)
          g2j=gj1*gj1+gj2*gj2+gj3*gj3
          gj=sqrt(g2j)
          enj=c05*(sqrt(c2+4.d0*g2j*pi2a**2)-clight)
          c2j=c2+enj
          anormj=c2j/dsqrt(c2j*c2j+c2*g2j*pi2a**2)
          do isort=1,nsort
            const=fc(1)*anormj*c2*pi2a*gj/c2j
            call BESSR(pi2a*gj,smt(isort),BES,BESD,lmb(isort)+1)
            li=0
            do l=0,lmb(isort)
              ist=-1
              if(l.eq.0) ist=1
              if(rel_interst) then
                cs=fc(1)*bes(l)*anormj
                do ii=ist,1,2
                  li=li+1
                  g1(li,isort,jgb)=cs
                  g2(li,isort,jgb)=const*bes(l+ii)
                  if(ii.eq.1) g2(li,isort,jgb)=-g2(li,isort,jgb)
                enddo
              else if(.not.rel_interst) then
                cs=fc(1)*bes(l)
                cs1=fc(1)*besd(l)
                do ii=ist,1,2
                  li=li+1
                  g1(li,isort,jgb)=cs
                  g2(li,isort,jgb)=cs1
                enddo
              endif
            enddo
          enddo
        enddo
        do jgb=1,nbndt
          jgb0=indg(jgb)
          gkj=pi2a*(pnt_+gb(:,jgb))
          g2j=dot_product(gkj,gkj)
          gj=sqrt(g2j)
          do isort=1,nsort
	      smtcon=0.5d0*smt(isort)**2
            do li=1,2*lmb(isort)+1
	        l=li/2
              ie=n_e_apw(l,isort)
              vhv0=ffhmt(1,1,ie,ie,li,isort,ispin)
c -------- Adding correction VII-2-9 ----------------------------------
              vhs0=0.5d0*(ffhmt(1,2,ie,ie,li,isort,ispin)
     &                   +ffhmt(2,1,ie,ie,li,isort,ispin))
	        shv0=vhs0
              shs0=ffhmt(2,2,ie,ie,li,isort,ispin)
              vv=ffsmt(1,1,ie,ie,li,isort,ispin)
              vs=ffsmt(1,2,ie,ie,li,isort,ispin)
              ss=ffsmt(2,2,ie,ie,li,isort,ispin)
              gshs(li,isort)=g1(li,isort,jgb)*shv0
     &                      +g2(li,isort,jgb)*shs0
              gvv(li,isort)=g1(li,isort,jgb)*vv+g2(li,isort,jgb)*vs
              gss(li,isort)=g1(li,isort,jgb)*vs+g2(li,isort,jgb)*ss
              gvhv(li,isort)=g1(li,isort,jgb)*vhv0+g2(li,isort,jgb)*vhs0
            enddo
          enddo
          do igb=1,nbndt
            gki=pi2a*(pnt_+gb(:,igb))
            g2i=dot_product(gki,gki)
            gi=sqrt(g2i)
            igb0=indg(igb)
            ia=igbs(1,jgb0)-igbs(1,igb0)
            ib=igbs(2,jgb0)-igbs(2,igb0)
            ic=igbs(3,jgb0)-igbs(3,igb0)
            scal=dot_product(gki,gkj)
            call cross(cr,gki,gkj)
            p(0)=1.d0
            pp(0)=0.d0
            amult=gi*gj
            if(maxb.gt.0.and.amult>1.d-8) then
              x=scal/gi/gj
              p(1)=x
              pp(1)=1.d0
              do l=1,maxb-1
                p(l+1)=(dfloat(l+l+1)*x*p(l)-dfloat(l)*p(l-1))
     &			                                    /dfloat(l+1)
                pp(l+1)=dfloat(l+1)*p(l)+x*pp(l)
              enddo
            else
              do l=1,maxb
                p(l)=0.d0
                pp(l)=0.d0
              enddo
            endif
            do jp=1,2
              jh=(jp-1)*nbndt+jgb
              do ipp=1,2
                ih=(ipp-1)*nbndt+igb
C     ******************************************************************
                sr=0.d0
                si=0.d0
                if(amult<1.d-8) then
                  sr=0.d0
                  si=0.d0
                else
                  if(ipp.eq.1.and.jp.eq.1) then
                    sr=0.d0
                    si=-cr(3)/amult
                  else if(ipp.eq.1.and.jp.eq.2) then
                    sr=-cr(2)/amult
                    si=cr(1)/amult
                  else if(ipp.eq.2.and.jp.eq.1) then
                    sr=cr(2)/amult
                    si=cr(1)/amult
                  else if(ipp.eq.2.and.jp.eq.2) then
                    sr=0.d0
                    si=cr(3)/amult
                  endif
                endif
                li=0
                do l=0,maxb
                  ist=-1
                  if(l.eq.0) ist=1
                  cons=pp(l)*sr
                  consi=pp(l)*si
                  do ii=ist,1,2
                    li=li+1
                    am=1.d0
                    if(ii.eq.1) am=-am
                    ppr(li)=cons*am
                    ppi(li)=consi*am
                    if(ipp.eq.jp) then
                      abskap=dfloat(ii)/2.d0
                      abskap=abs(abskap*(abskap+l+l+1)+0.25d0)
                      ppr(li)=ppr(li)+abskap*p(l)
                    endif
                  enddo   !!! over ii
                enddo   !!! over l
                do isort=1,nsort
                  fac=0.d0
                  if(complex_ro) faci=0.d0
                  do iatom=1,natom
                    if(is(iatom).eq.isort) then
                      sp=(gkj(1)-gki(1))*tau(1,iatom)+
     &                   (gkj(2)-gki(2))*tau(2,iatom)+
     &                   (gkj(3)-gki(3))*tau(3,iatom)
                      sp=par*sp
                      fac=fac+cos(sp)
                      if(complex_ro) faci=faci+sin(sp)
                    endif
                  enddo
C     ***     O_mt    **   *******************************************
                  sum=0.d0
                  sumi=0.d0
                  do li=1,2*lmb(isort)+1
                    con=g1(li,isort,igb)*gvv(li,isort)+
     &                  g2(li,isort,igb)*gss(li,isort)
                    sum=sum+ppr(li)*con
                    sumi=sumi+ppi(li)*con
                  enddo
                  br1(ih,jh)=br1(ih,jh)+fac*sum
                  bi1(ih,jh)=bi1(ih,jh)+fac*sumi
                  if(complex_ro) then
                    br1(ih,jh)=br1(ih,jh)-faci*sumi
                    bi1(ih,jh)=bi1(ih,jh)+faci*sum
                  endif
C     ***     H_mt    **   *******************************************
                  sum=0.d0
                  sumi=0.d0
                  do li=1,2*lmb(isort)+1
                    con=g1(li,isort,igb)*gvhv(li,isort)+
     &                  g2(li,isort,igb)*gshs(li,isort)
                    sum=sum+ppr(li)*con
                    sumi=sumi+ppi(li)*con
                  enddo
                  ar1(ih,jh)=ar1(ih,jh)+fac*sum
                  ai1(ih,jh)=ai1(ih,jh)+fac*sumi
                  if(complex_ro) then
                    ar1(ih,jh)=ar1(ih,jh)-faci*sumi
                    ai1(ih,jh)=ai1(ih,jh)+faci*sum
                  endif
                enddo  !!! over isort
              enddo  !!! over ip
            enddo  !!! over jp
          enddo  !!! over igb
        enddo  !!! over jgb
      endif  !!! over irel
	deallocate(p,gvv,gss,gvhv,gshs,g1,g2,fc,pp,ppr,ppi,bes,besd)
      end
