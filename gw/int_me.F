      subroutine int_me(a,pt,ispin,ndim,ndim0,nbas,indg,evbn,w,wi)
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin,ndim,ndim0,nbas,indg(nbndtm)
	real*8, intent(in) :: pt(3),w(nplwdif,nspin),wi(nplwdif,nspin)
	complex*16, intent(in) :: evbn(nbasmpw,nbndf)
	real*8, intent(inout) :: a(ndim0,ndim0)
	integer :: nbndt,jgb,jgb0,igb,igb0,ia,ib,ic,ind,jp,jh,ipp,ih,iv
	real*8 :: gkj(3),conbr(3),conbi(3),cr(3),gki(3),pi2,pi2a,c1,c05,
     &	        c2,g2j,gj,enj,c2j,anormj,g2i,gi,eni,c2i,anormi,const,
     &          const1,const2,scal,ann,amult,pr,pii,rel0r,rel0i,
     &          relr,reli,arr,aii,hr1,hi1,hr2,hi2,hr3,hi3
      complex*16 :: bs
	complex*16, allocatable :: ai(:,:),tmp(:,:),tmp1(:,:)
	nbndt=nbas/nrel
	allocate(ai(nbas,nbas),tmp(nbas,ndim))
	allocate(tmp1(ndim,ndim))
	pi2=2.d0*pi
      pi2a=pi2/par
      c1=amegaint/amega
      ai=(0.d0,0.d0)
      if(irel.ne.2.or..not.rel_interst) then
        do jgb=1,nbndt
          jgb0=indg(jgb)
          do igb=jgb,nbndt
            igb0=indg(igb)
            ia=igbs(1,jgb0)-igbs(1,igb0)
            ib=igbs(2,jgb0)-igbs(2,igb0)
            ic=igbs(3,jgb0)-igbs(3,igb0)
            ind=indplw(ia,ib,ic)
	      if(.not.complex_ro) then
	        ai(igb,jgb)=dcmplx(w(ind,ispin),0.d0)
	        bs=b_extval*sovr(ind)
            else if(complex_ro) then
              ai(igb,jgb)=dcmplx(w(ind,ispin),wi(ind,ispin))
	        bs=b_extval*dcmplx(sovr(ind),sovi(ind))
	      endif
C     * Here IREL==2 but no spin-orbit in inrestitial ************
            if(irel==2) ai(nbndt+igb,nbndt+jgb)=ai(igb,jgb)
c --- We add external magnetic field here -------------------------
            if(irel<=1) then
c              if(ispin==1) ai(igb,jgb)=ai(igb,jgb)+bs
c              if(ispin==nspin) ai(igb,jgb)=ai(igb,jgb)-bs
            else if(irel==2) then
              ai(igb,jgb)=ai(igb,jgb)+bs
              ai(nbndt+igb,nbndt+jgb)=ai(nbndt+igb,nbndt+jgb)-bs
            endif
C     *************************************************************
          enddo           !!! over igb
        enddo               !!! over jgb
      else if(irel.eq.2) then
        c05=0.5d0*clight
        c2=clight*clight
        do jgb=1,nbndt
          jgb0=indg(jgb)
          gkj=pi2a*(pt+gbs(:,jgb0))
          g2j=dot_product(gkj,gkj)
          gj=sqrt(g2j)
          enj=c05*(sqrt(c2+4.d0*g2j)-clight)
          c2j=c2+enj
          anormj=c2j/dsqrt(c2j*c2j+c2*g2j)
          do igb=1,nbndt
            igb0=indg(igb)
            gki=pi2a*(pt+gbs(:,igb0))
            g2i=dot_product(gki,gki)
            gi=sqrt(g2i)
            eni=c05*(sqrt(c2+4.d0*g2i)-clight)
            c2i=c2+eni
            anormi=c2i/sqrt(c2i*c2i+c2*g2i)
            ia=igbs(1,jgb0)-igbs(1,igb0)
            ib=igbs(2,jgb0)-igbs(2,igb0)
            ic=igbs(3,jgb0)-igbs(3,igb0)
            ind=indplw(ia,ib,ic)
            const=c2/c2i/c2j
            const1=clight/c2i
            const2=clight/c2j
            scal=dot_product(gki,gkj)
            call cross(cr,gki,gkj)
            ann=anormi*anormj
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
                ih=(ipp-1)*nbndt+igb
c     &&&&& Begin  IV-40-24 &&&&&&&&&&&&&&&&&&&
                pr=0.0d0
                pii=0.0d0
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
                if(inv.ne.0) then
                  ai(ih,jh)=dcmplx(relr,reli)*w(ind,ispin)
                else
                  arr=relr*w(ind,1)-reli*wi(ind,1)
                  aii=relr*wi(ind,1)+reli*w(ind,1)
                  ai(ih,jh)=dcmplx(arr,aii)
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
                  arr=conbr(1)*hr1-conbi(1)*hi1
                  aii=conbr(1)*hi1+conbi(1)*hr1
                  arr=arr+conbr(2)*hr2-conbi(2)*hi2
                  aii=aii+conbr(2)*hi2+conbi(2)*hr2
                  arr=arr+conbr(3)*hr3-conbi(3)*hi3
                  aii=aii+conbr(3)*hi3+conbi(3)*hr3
                  ai(ih,jh)=ai(ih,jh)+dcmplx(arr,aii)
                endif   !!! over magn=2
              enddo  !!! over ipp
            enddo  !!! over jp
          enddo  !!! over igb
        enddo  !!! over jgb
      endif  !!! over irel
C     ******************************************************************
      call zhemm('l','l',nbas,ndim,(1.d0,0.d0),ai,nbas,evbn,nbasmpw,
     &           (0.d0,0.d0),tmp,nbas)
      call zgemm('c','n',ndim,ndim,nbas,(1.d0,0.d0),evbn,nbasmpw,tmp,
     &           nbas,(0.d0,0.d0),tmp1,ndim)
	call pack_hermit(tmp1,a,ndim,ndim,ndim0,1.d0,1.d0)
	deallocate(ai,tmp,tmp1)
      end
