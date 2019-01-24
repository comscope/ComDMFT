      subroutine overlap_lapw(gb,pnt_,nbas,ov,indg)
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: nbas,indg(nbas/nrel)
      real*8, intent(in) :: gb(3,nbas/nrel),pnt_(3)
      complex*16, intent(out) :: ov(nbas,nbas)
	integer :: jgb,jgb0,igb,igb0,ia,ib,ic,ind,jp,jh,ipp,ih,nbndt
      real*8 :: cr(3),gki(3),gkj(3),pi2a,c1,gj1,
     &	        gj2,gj3,gj,gj22,gi1,gi2,gi3,gi,gij22,gi22,c05,c2,g2j,
     &          enj,c2j,anormj,g2i,eni,c2i,anormi,const,const1,const2,
     &          scal,ann,coni,conr,amult,pr,pii,rel0r,rel0i,relr,reli,
     &          pi2,brij,biij,ddot,br,bi
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ov=(0.d0,0.d0)
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
          do igb=1,nbndt
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
            if(igb.eq.jgb) then
              ov(igb,jgb)=c1
            else
              if(.not.complex_ro) ov(igb,jgb)=sovr(ind)
              if(complex_ro) ov(igb,jgb)=dcmplx(sovr(ind),sovi(ind))
            endif
C     * For IREL==2 but no spin-orbit in inrestitial ************
            if(irel==2) then
              ov(nbndt+igb,nbndt+jgb)=ov(igb,jgb)
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
            brij=c1
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
            amult=gi*gj
            do jp=1,2
              jh=(jp-1)*nbndt+jgb
              do ipp=1,2
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
                if(.not.complex_ro) then
                  br=relr*brij
                  bi=reli*brij
                else
                  br=relr*brij-reli*biij
                  bi=relr*biij+reli*brij
                endif
                ov(ih,jh)=dcmplx(br,bi)
              enddo  !!! over ip
            enddo  !!! over jp
          enddo  !!! over igb
        enddo  !!! over jgb
      endif  !!! over irel
      end
