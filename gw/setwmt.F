      subroutine setwmt(gb,pnt_,nbas,ar1,ai1,br1,bi1,indg,ispin,key)
c     key = 0 Only Hamiltonian matrix
c     key = 1 Both
c     key = 2 Only Overlap matrix
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: nbas,indg(nbas/nrel),ispin,key
      real*8, intent(in) :: gb(3,nbas/nrel),pnt_(3)
      real*8, intent(out) :: ar1(nbas,nbas),ai1(nbas,nbas),
     &	                     br1(nbas,nbas),bi1(nbas,nbas)
	integer :: jgb,jgb0,igb,igb0,ia,ib,ic,ind,jp,jh,ipp,ih,nbndt,iv
      real*8 :: cr(3),gki(3),gkj(3),pi2a,c1,gj1,conbr(3),conbi(3),
     &	        gj2,gj3,gj,gj22,gi1,gi2,gi3,gi,gij22,gi22,c05,c2,g2j,
     &          enj,c2j,anormj,g2i,eni,c2i,anormi,const,const1,const2,
     &          scal,ann,coni,conr,amult,pr,pii,rel0r,rel0i,relr,reli,
     &          enij,xr0,xi0,pi2,brij,biij,ddot,hr1,hr2,hr3,hi1,hi2,hi3
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(key<=1) then
	  ai1=0.d0
	  ar1=0.d0
	endif
      if(key>=1) then
	  br1=0.d0
	  bi1=0.d0
	endif
      pi2=2.d0*pi
      pi2a=pi2/par
      c1=amegaint/amega
	nbndt=nbas/nrel
      if(irel.ne.2.or..not.rel_interst) then
        do jgb=1,nbndt
          jgb0=indg(jgb)
          gj1=pnt_(1)+gb(1,jgb)
          gj2=pnt_(2)+gb(2,jgb)
          gj3=pnt_(3)+gb(3,jgb)
          gj=gj1*gj1+gj2*gj2+gj3*gj3
          gj22=pi2a**2*gj
          gj=sqrt(gj)
          do igb=jgb,nbndt
            gi1=pnt_(1)+gb(1,igb)
            gi2=pnt_(2)+gb(2,igb)
            gi3=pnt_(3)+gb(3,igb)
            gi=gi1*gi1+gi2*gi2+gi3*gi3
	      gij22=gj22
            gi22=pi2a**2*gi
	      gij22=0.5d0*(gi22+gj22)
            gi=sqrt(gi)
            igb0=indg(igb)
            ia=igbs(1,jgb0)-igbs(1,igb0)
            ib=igbs(2,jgb0)-igbs(2,igb0)
            ic=igbs(3,jgb0)-igbs(3,igb0)
            ind=indplw(ia,ib,ic)
C     ***     O_int   *********************************************
            if(key>=1) then
              if(igb.eq.jgb) then
                br1(igb,jgb)=c1
                if(complex_ro) bi1(igb,jgb)=0.d0
              else
                br1(igb,jgb)=sovr(ind)
                if(complex_ro) bi1(igb,jgb)=sovi(ind)
              endif
            endif
C     ***     H_int   *********************************************
            if(key<=1) then
              ar1(igb,jgb)=gij22*br1(igb,jgb)+war(ind,ispin)
              if(complex_ro) ai1(igb,jgb)=gij22*bi1(igb,jgb)
     &                                   +wari(ind,ispin)
            endif
C     * For IREL==2 but no spin-orbit in inrestitial ************
            if(irel==2) then
              if(key<=1) then
                ar1(nbndt+igb,nbndt+jgb)=ar1(igb,jgb)
                ai1(nbndt+igb,nbndt+jgb)=ai1(igb,jgb)
                if(magn==2) then
                  ar1(igb,nbndt+jgb)=ar1(igb,nbndt+jgb)
     &                              +warb(1,ind)-warbi(2,ind)
                  ai1(igb,nbndt+jgb)=ai1(igb,nbndt+jgb)
     &                              +warbi(1,ind)+warb(2,ind)
                  ar1(nbndt+igb,jgb)=ar1(nbndt+igb,jgb)
     &                              +warb(1,ind)+warbi(2,ind)
                  ai1(nbndt+igb,jgb)=ai1(nbndt+igb,jgb)
     &                              +warbi(1,ind)-warb(2,ind)
                  ar1(nbndt+igb,nbndt+jgb)=ar1(nbndt+igb,nbndt+jgb)
     &                                    +warb(3,ind)
                  ai1(nbndt+igb,nbndt+jgb)=ai1(nbndt+igb,nbndt+jgb)
     &                                    +warbi(3,ind)
                  ar1(igb,jgb)=ar1(igb,jgb)-warb(3,ind)
                  ai1(igb,jgb)=ai1(igb,jgb)-warbi(3,ind)
                endif
              endif
              if(key>=1) then
                br1(nbndt+igb,nbndt+jgb)=br1(igb,jgb)
                bi1(nbndt+igb,nbndt+jgb)=bi1(igb,jgb)
              endif
            endif
c --- We add external magnetic field here -------------------------
            if(irel<=1) then
              if(ispin==1) then
                ar1(igb,jgb)=ar1(igb,jgb)+b_extval*sovr(ind)
                if(complex_ro) ai1(igb,jgb)=ai1(igb,jgb)
     &                                     +b_extval*sovi(ind)
              endif
              if(ispin==nspin) then
                ar1(igb,jgb)=ar1(igb,jgb)-b_extval*sovr(ind)
                if(complex_ro) ai1(igb,jgb)=ai1(igb,jgb)
     &                                     -b_extval*sovi(ind)
              endif
            else if(irel==2) then
              ar1(igb,jgb)=ar1(igb,jgb)+b_extval*sovr(ind)
              ar1(nbndt+igb,nbndt+jgb)=ar1(nbndt+igb,nbndt+jgb)
     &                                -b_extval*sovr(ind)
              if(complex_ro) then
                ai1(igb,jgb)=ai1(igb,jgb)+b_extval*sovi(ind)
                ai1(nbndt+igb,nbndt+jgb)=ai1(nbndt+igb,nbndt+jgb)
     &                                  -b_extval*sovi(ind)
              endif
            endif
          enddo           !!! over igb
        enddo               !!! over jgb
      else if(irel.eq.2) then
        c05=0.5d0*clight
        c2=clight*clight
        do jgb=1,nbndt
          jgb0=indg(jgb)
          gkj=pi2a*(pnt_+gb(:,jgb))
          g2j=dot_product(gkj,gkj)
          gj=sqrt(g2j)
          enj=c05*(sqrt(c2+4.d0*g2j)-clight)
          c2j=c2+enj
          anormj=c2j/dsqrt(c2j*c2j+c2*g2j)
          do igb=1,nbndt
            gki=pi2a*(pnt_+gb(:,igb))
            g2i=dot_product(gki,gki)
            gi=sqrt(g2i)
            eni=c05*(sqrt(c2+4.d0*g2i)-clight)
            c2i=c2+eni
            anormi=c2i/sqrt(c2i*c2i+c2*g2i)
            igb0=indg(igb)
            ia=igbs(1,jgb0)-igbs(1,igb0)
            ib=igbs(2,jgb0)-igbs(2,igb0)
            ic=igbs(3,jgb0)-igbs(3,igb0)
            ind=indplw(ia,ib,ic)
            biij=0.0d0
            if(igb.eq.jgb) then
              brij=c1
              if(complex_ro) biij=0.d0
            else
              brij=sovr(ind)
              if(complex_ro) biij=sovi(ind)
            endif
            const=c2/c2i/c2j
            const1=clight/c2i
            const2=clight/c2j
            scal=ddot(3,gki,1,gkj,1)
            call cross(cr,gki,gkj)
            ann=anormi*anormj
            coni=0.d0
            if(complex_ro) coni=ann*biij
            conr=ann*brij
            do iv=1,3
              conbr(iv)=0.d0
              conbi(iv)=0.d0
            enddo
            if(magn.eq.2) then
              do iv=1,3
                if(inv.lt.2) conbr(iv)=ann*warb(iv,ind)
                if(inv.ne.1) conbi(iv)=ann*warbi(iv,ind)
              enddo
            endif
            amult=gi*gj
            do jp=1,2
              jh=(jp-1)*nbndt+jgb
              do ipp=1,2
                pr=0.0d0
                pii=0.0d0
                ih=(ipp-1)*nbndt+igb
c     &&&&& Begin  IV-40-24 &&&&&&&&&&&&&&&&&&&
                if(ipp.eq.1.and.jp.eq.1) then
                  pr=scal
                  pii=-cr(3)
                else if(ipp.eq.1.and.jp.eq.2) then
                  pr=-cr(2)
                  pii=cr(1)
                else if(ipp.eq.2.and.jp.eq.1) then
                  pr=cr(2)
                  pii=cr(1)
                else if(ipp.eq.2.and.jp.eq.2) then
                  pr=scal
                  pii=cr(3)
                endif
                rel0r=const*pr
                rel0i=const*pii
c     &&&&& End  IV-40-24 &&&&&&&&&&&&&&&&&&&
                if(ipp.eq.jp) then
                  relr=ann*(1.d0+rel0r)
                else
                  relr=ann*rel0r
                endif
                reli=ann*rel0i
c     &&&&& End  IV-40-25 &&&&&&&&&&&&&&&&&&&
C     ***     O_int   *********************************************
                if(key>=1) then
                  if(.not.complex_ro) then
                    br1(ih,jh)=relr*brij
                    bi1(ih,jh)=reli*brij
                  else
                    br1(ih,jh)=relr*brij-reli*biij
                    bi1(ih,jh)=relr*biij+reli*brij
                  endif
                endif
c ----------- Symmetrization -------------------------------
                enij=0.5d0*(eni+enj)	            
                xr0=enij*rel0r
                xi0=enij*rel0i
                if(ipp.eq.jp) then
                  xr0=xr0+0.5d0*c2*(g2i/c2i+g2j/c2j)
                endif
C     ***     H_int   *********************************************
                if(key<=1) then
                  ar1(ih,jh)=conr*xr0-coni*xi0
                  ai1(ih,jh)=conr*xi0+coni*xr0
                  if(.not.complex_ro) then
                    ar1(ih,jh)=ar1(ih,jh)+relr*war(ind,ispin)
                    ai1(ih,jh)=ai1(ih,jh)+reli*war(ind,ispin)
                  else
                    ar1(ih,jh)=ar1(ih,jh)+relr*war(ind,ispin)
     &                                   -reli*wari(ind,ispin)
                    ai1(ih,jh)=ai1(ih,jh)+relr*wari(ind,ispin)
     &                                   +reli*war(ind,ispin)
                  endif
                endif
                if(magn.eq.2) then
                  hr1=0.0d0
                  hr2=0.0d0
                  hr3=0.0d0
                  hi1=0.0d0
                  hi2=0.0d0
                  hi3=0.0d0
                  if(ipp.eq.1.and.jp.eq.1) then
                    hr1=-gki(1)*gkj(3)-gkj(1)*gki(3)
                    hr1=-const*hr1
                    hi1=cr(1)*const
                    hr2=-gki(2)*gkj(3)-gkj(2)*gki(3)
                    hr2=-const*hr2
                    hi2=cr(2)*const
                    hr3=-2.d0*gki(3)*gkj(3)+scal
                    hr3=-1.d0-const*hr3
                    hi3=cr(3)*const
                  else if(ipp.eq.1.and.jp.eq.2) then
                    hr1=gki(1)*gkj(1)+gkj(1)*gki(1)-scal
                    hr1=1.d0-const*hr1
                    hi1=gki(1)*gkj(2)+gkj(1)*gki(2)
                    hi1=-const*hi1
                    hr2=gki(2)*gkj(1)+gkj(2)*gki(1)
                    hr2=-const*hr2
                    hi2=gki(2)*gkj(2)+gkj(2)*gki(2)-scal
                    hi2=1.d0-const*hi2
                    hr3=gki(3)*gkj(1)+gkj(3)*gki(1)
                    hi3=gki(3)*gkj(2)+gkj(3)*gki(2)
                    hr3=-hr3*const
                    hi3=-hi3*const
                  else if(ipp.eq.2.and.jp.eq.1) then
                    hr1=gki(1)*gkj(1)+gkj(1)*gki(1)-scal
                    hr1=1.d0-const*hr1
                    hi1=-gki(1)*gkj(2)-gkj(1)*gki(2)
                    hi1=-const*hi1
                    hr2=gki(2)*gkj(1)+gkj(2)*gki(1)
                    hr2=-const*hr2
                    hi2=-gki(2)*gkj(2)-gkj(2)*gki(2)+scal
                    hi2=-1.d0-const*hi2
                    hr3=gki(3)*gkj(1)+gkj(3)*gki(1)
                    hi3=-gki(3)*gkj(2)-gkj(3)*gki(2)
                    hr3=-hr3*const
                    hi3=-hi3*const
                  else if(ipp.eq.2.and.jp.eq.2) then
                    hr1=gki(1)*gkj(3)+gkj(1)*gki(3)
                    hr1=-const*hr1
                    hi1=cr(1)*const
                    hr2=gki(2)*gkj(3)+gkj(2)*gki(3)
                    hr2=-const*hr2
                    hi2=cr(2)*const
                    hr3=2.d0*gki(3)*gkj(3)-scal
                    hr3=1.d0-const*hr3
                    hi3=cr(3)*const
                  endif
                  ar1(ih,jh)=ar1(ih,jh)+conbr(1)*hr1-conbi(1)*hi1
                  ai1(ih,jh)=ai1(ih,jh)+conbr(1)*hi1+conbi(1)*hr1
                  ar1(ih,jh)=ar1(ih,jh)+conbr(2)*hr2-conbi(2)*hi2
                  ai1(ih,jh)=ai1(ih,jh)+conbr(2)*hi2+conbi(2)*hr2
                  ar1(ih,jh)=ar1(ih,jh)+conbr(3)*hr3-conbi(3)*hi3
                  ai1(ih,jh)=ai1(ih,jh)+conbr(3)*hi3+conbi(3)*hr3
                endif   !!! over magn=2
C     ******************************************************************
              enddo  !!! over ip
            enddo  !!! over jp
          enddo  !!! over igb
        enddo  !!! over jgb
      endif  !!! over irel
      end
